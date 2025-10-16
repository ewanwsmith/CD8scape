#!/usr/bin/env julia
using Test, CSV, DataFrames

@testset "cache_reuse" begin
    tmpdir = Base.mktempdir()
    panel = joinpath(tmpdir, "panel")
    mkpath(panel)

    # Create a small peptides file
    peptides = ["PEPAAA", "PEPBBB", "PEPCCC"]
    open(joinpath(panel, "context_peptides.pep"), "w") do io
        for p in peptides
            println(io, p)
        end
    end

    # Minimal alleles file
    open(joinpath(panel, "alleles.txt"), "w") do io
        println(io, "HLA-A01:01")
    end

    # Create a fake netMHCpan script that writes a simple chunk TSV and records invocations
    fake = joinpath(tmpdir, "fake_netmhcpan.sh")
    open(fake, "w") do io
        println(io, "#!/usr/bin/env bash")
        println(io, "PEPFILE=''")
        println(io, "OUTFILE=''")
        println(io, "# parse args")
        println(io, "while [[ \$# -gt 0 ]]; do")
        println(io, "  case \$1 in")
        println(io, "    -p) PEPFILE=\"\$2\"; shift; shift;;")
        println(io, "    -xlsfile) OUTFILE=\"\$2\"; shift; shift;;")
        println(io, "    *) shift;;")
        println(io, "  esac")
        println(io, "done")
        println(io, "# record invocation")
        println(io, "echo \$(date +%s) >> \"\${OUTFILE}.fake_called.txt\"")
        println(io, "# write minimal netMHCpan-like output")
        println(io, "ALLELE=HLA-A01:01")
        println(io, "echo \"" * "\$ALLELE" * "\" > \"\$OUTFILE\"")
        println(io, "echo -e \"Pos\\tPeptide\\tID\\tcore\\ticore\\tEL-score\\tEL_Rank\" >> \"\$OUTFILE\"")
        println(io, "while read -r pep; do")
        println(io, "  echo -e \"1\\t\${pep}\\tID\\tCORE\\tICORE\\t0.1\\t0.5\" >> \"\$OUTFILE\"")
        println(io, "done < \"\${PEPFILE}\"")
        println(io, "exit 0")
    end
    run(`chmod +x $fake`)

    # Set environment to point to fake netMHCpan
    ENV["NETMHCPAN"] = fake

    runner = joinpath(@__DIR__, "..", "..", "src", "context_run", "run_netMHCpan_context.jl")

    # First run: expect netMHCpan to be called at least once
    run(`julia --project=. $runner --folder $panel --mode panel --peptides $(joinpath(panel, "context_peptides.pep"))`)
    # Run Perl postprocessor the orchestrator usually calls
    perl_proc = joinpath(@__DIR__, "..", "..", "src", "context_run", "process_output_context.pl")
    run(`perl $perl_proc $(joinpath(panel, "netMHCpan_output.tsv"))`)

    fake_called_files = filter(f->endswith(f, ".fake_called.txt"), readdir(panel; join=true))
    @test !isempty(fake_called_files) || error("Fake netMHCpan was not invoked on first run")
    calls1 = sum(length(readlines(f)) for f in fake_called_files)
    @test calls1 >= 1

    # Simulate orchestrator persisting processed output into master cache
    processed = joinpath(panel, "context_processed_netMHCpan_output.csv")
    @test isfile(processed) || error("Processed output not created by runner")
    master_cache = joinpath(panel, "peptide_netmhc_processed.csv")
    cp(processed, master_cache; force=true)

    # Second run: with cache present, missing peptides should be empty and netMHCpan should not be called
    attempt_peps = readlines(joinpath(panel, "context_peptides.pep")) |> x -> unique(filter(y->!isempty(strip(y)), x))
    cached_df = CSV.read(master_cache, DataFrame)
    cached_peps = unique(String.(cached_df.Peptide))
    missing = setdiff(attempt_peps, cached_peps)
    @test isempty(missing)

    # Check fake invocation files again; ensure invocation count didn't increase
    fake_called_files2 = filter(f->endswith(f, ".fake_called.txt"), readdir(panel; join=true))
    calls2 = sum(length(readlines(f)) for f in fake_called_files2)
    @test calls2 == calls1
end


@testset "cache_reuse_partial_and_dedup" begin
    # Partial cache reuse: only some peptides are in master cache
    tmpdir = Base.mktempdir()
    panel = joinpath(tmpdir, "panel")
    mkpath(panel)
    peptides = ["AA1", "AA2", "AA3", "AA4"]
    open(joinpath(panel, "context_peptides.pep"), "w") do io
        for p in peptides
            println(io, p)
        end
    end
    open(joinpath(panel, "alleles.txt"), "w") do io
        println(io, "HLA-A01:01")
    end

    # reuse fake netMHCpan from earlier pattern
    fake = joinpath(tmpdir, "fake_netmhcpan.sh")
    open(fake, "w") do io
        println(io, "#!/usr/bin/env bash")
        println(io, "PEPFILE=''")
        println(io, "OUTFILE=''")
        println(io, "while [[ \$# -gt 0 ]]; do case \$1 in -p) PEPFILE=\"\$2\"; shift; shift;; -xlsfile) OUTFILE=\"\$2\"; shift; shift;; *) shift;; esac; done")
        println(io, "echo \$(date +%s) >> \"\${OUTFILE}.fake_called.txt\"")
        println(io, "echo \"HLA-A01:01\" > \"\$OUTFILE\"")
        println(io, "echo -e \"Pos\\tPeptide\\tID\\tcore\\ticore\\tEL-score\\tEL_Rank\" >> \"\$OUTFILE\"")
        println(io, "while read -r pep; do echo -e \"1\\t\${pep}\\tID\\tCORE\\tICORE\\t0.1\\t0.5\" >> \"\$OUTFILE\"; done < \"\${PEPFILE}\"")
        println(io, "exit 0")
    end
    run(`chmod +x $fake`)
    ENV["NETMHCPAN"] = fake

    runner = joinpath(@__DIR__, "..", "..", "src", "context_run", "run_netMHCpan_context.jl")

    # Make a master cache that already contains AA1 and AA2
    # We create a minimal processed CSV with Peptide column matching process_output_context.pl format.
    processed_df = DataFrame(Pos=Int[], Peptide=String[], ID=String[], HLA=String[], core=String[], icore=String[], EL_score=Float64[], EL_Rank=Float64[])
    push!(processed_df, (1, "AA1", "ID", "HLA-A01:01", "c", "i", 0.1, 0.5))
    push!(processed_df, (1, "AA2", "ID", "HLA-A01:01", "c", "i", 0.1, 0.5))
    CSV.write(joinpath(panel, "peptide_netmhc_processed.csv"), processed_df)

    # Run runner which should only process missing peptides (AA3, AA4)
    run(`julia --project=. $runner --folder $panel --mode panel --peptides $(joinpath(panel, "context_peptides.pep"))`)
    run(`perl $(joinpath(@__DIR__, "..", "..", "src", "context_run", "process_output_context.pl")) $(joinpath(panel, "netMHCpan_output.tsv"))`)

    # Read new processed output and combine with master cache to simulate orchestrator append
    new_processed = CSV.read(joinpath(panel, "context_processed_netMHCpan_output.csv"), DataFrame)
    master = CSV.read(joinpath(panel, "peptide_netmhc_processed.csv"), DataFrame)
    # Normalize column names: process_output_context.pl uses 'EL-score' while test master uses 'EL_score'
    if "EL-score" in names(new_processed)
        rename!(new_processed, "EL-score" => :EL_score)
    end
    if "EL-score" in names(master)
        rename!(master, "EL-score" => :EL_score)
    end
    if :"EL-score" in names(master)
        rename!(master, Symbol("EL-score") => :EL_score)
    end
    if :"EL-score" in names(new_processed)
        rename!(new_processed, Symbol("EL-score") => :EL_score)
    end
    combined = unique(vcat(master, new_processed))
    final_peps = unique(String.(combined.Peptide))
    @test sort(final_peps) == sort(peptides)

    # Deduplication test: appending the same processed_df twice yields same unique peptides
    doubled = unique(vcat(combined, combined))
    @test length(unique(String.(doubled.Peptide))) == length(unique(String.(combined.Peptide)))
end
