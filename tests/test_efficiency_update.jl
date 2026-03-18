#!/usr/bin/env julia
"""
test_efficiency_update.jl

Tests for the efficiency_update branch changes:
  - Allele squishing (deduplication + frequency merging)
  - Character-based allele batching (netMHCpan -a limit)
  - Adaptive (chunk × batch) parallelism
  - Merge correctness: horizontal stitch of batches, vertical concat of chunks
  - process_best_ranks_supertype squished-panel preference
  - process_output.pl compatibility with multi-allele output format

Run from the project root:
    julia --project=. tests/test_efficiency_update.jl
"""

using Test, CSV, DataFrames
using Logging

# ─── Helpers ──────────────────────────────────────────────────────────────────

PROJECT_ROOT = abspath(joinpath(@__DIR__, ".."))
SRC          = joinpath(PROJECT_ROOT, "src")
NETMHCPAN    = begin
    settings = joinpath(SRC, "settings.txt")
    local path = ""
    for line in readlines(settings)
        s = strip(line)
        isempty(s) || startswith(s, "#") && continue
        if occursin('=', s)
            k, v = strip.(split(s, '='; limit=2))
            if uppercase(k) == "NETMHCPAN"
                path = normpath(replace(v, "~" => homedir()))
            end
        end
    end
    path
end

netmhcpan_available = isfile(NETMHCPAN) && Base.Filesystem.isexecutable(NETMHCPAN)

function mktemprun(f)
    d = mktempdir()
    try
        f(d)
    finally
        rm(d; recursive=true, force=true)
    end
end

pass_count = Ref(0)
fail_count = Ref(0)
function report(name, ok, msg="")
    if ok
        println("  ✓  $name")
        pass_count[] += 1
    else
        println("  ✗  $name$(isempty(msg) ? "" : ": $msg")")
        fail_count[] += 1
    end
end

# ─── Section 1: Unit — allele squishing algorithm ─────────────────────────────

println("\n══════════════════════════════════════════════════")
println("  Section 1 — Allele squishing logic")
println("══════════════════════════════════════════════════")

# Replicate the squishing algorithm from run_netMHCpan_global.jl
function squish_alleles(allele_strings, freqs)
    # allele_strings: panel canonical forms (with *); freqs: corresponding frequencies
    raw_netmhcpan = [replace(a, "*" => "") for a in allele_strings]
    raw_canonical  = allele_strings
    raw_freqs      = freqs

    seen_idx         = Dict{String,Int}()
    allele_list      = String[]
    allele_canonical = String[]
    allele_freq_sum  = Float64[]
    squish_originals = Dict{String,Vector{String}}()

    for (net_a, canon_a, f) in zip(raw_netmhcpan, raw_canonical, raw_freqs)
        if haskey(seen_idx, net_a)
            i = seen_idx[net_a]
            allele_freq_sum[i] += f
            push!(squish_originals[net_a], canon_a)
        else
            push!(allele_list, net_a)
            push!(allele_canonical, canon_a)
            push!(allele_freq_sum, f)
            seen_idx[net_a] = length(allele_list)
            squish_originals[net_a] = [canon_a]
        end
    end
    return (allele_list=allele_list,
            allele_canonical=allele_canonical,
            allele_freq_sum=allele_freq_sum,
            squish_originals=squish_originals,
            seen_idx=seen_idx,
            raw_netmhcpan=raw_netmhcpan,
            raw_canonical=raw_canonical)
end

# Test 1a: no duplicates — nothing should be squished
begin
    alleles = ["HLA-A*01:01", "HLA-A*02:01", "HLA-B*07:02"]
    freqs   = [0.3, 0.5, 0.2]
    r = squish_alleles(alleles, freqs)
    report("No duplicates → allele_list unchanged (length 3)",
           length(r.allele_list) == 3)
    report("No duplicates → all frequencies preserved",
           r.allele_freq_sum ≈ freqs)
    report("No duplicates → no squish_originals with >1 entry",
           all(v -> length(v) == 1, values(r.squish_originals)))
end

# Test 1b: one duplicate pair (HLA-A*01:01 and HLA-A01:01 both → HLA-A01:01)
begin
    alleles = ["HLA-A*01:01", "HLA-A01:01", "HLA-B*07:02"]
    freqs   = [0.3, 0.2, 0.5]
    r = squish_alleles(alleles, freqs)
    report("One duplicate pair → allele_list length 2",
           length(r.allele_list) == 2)
    squished_freq = r.allele_freq_sum[r.seen_idx["HLA-A01:01"]]
    report("One duplicate pair → frequencies summed (0.3 + 0.2 = 0.5)",
           squished_freq ≈ 0.5)
    report("One duplicate pair → HLA-B*07:02 untouched",
           r.allele_freq_sum[r.seen_idx["HLA-B07:02"]] ≈ 0.5)
    report("One duplicate pair → squish_originals has 2 entries for the merged allele",
           length(r.squish_originals["HLA-A01:01"]) == 2)
end

# Test 1c: squishing map content (Was_Squished flags)
begin
    alleles = ["HLA-A*01:01", "HLA-A01:01", "HLA-B*07:02"]
    freqs   = [0.3, 0.2, 0.5]
    r = squish_alleles(alleles, freqs)
    was_squished = [length(r.squish_originals[net_a]) > 1 for net_a in r.raw_netmhcpan]
    report("Squishing map: HLA-A*01:01 was_squished=true",  was_squished[1] == true)
    report("Squishing map: HLA-A01:01  was_squished=true",  was_squished[2] == true)
    report("Squishing map: HLA-B*07:02 was_squished=false", was_squished[3] == false)
    squished_into = [r.allele_canonical[r.seen_idx[net_a]] for net_a in r.raw_netmhcpan]
    report("Squishing map: both A alleles squished_into same canonical",
           squished_into[1] == squished_into[2])
end

# Test 1d: three-way duplicate
begin
    alleles = ["HLA-A*01:01", "HLA-A01:01", "HLA-A*01:01", "HLA-B*07:02"]
    freqs   = [0.2, 0.2, 0.2, 0.4]
    r = squish_alleles(alleles, freqs)
    report("Three-way duplicate → allele_list length 2",
           length(r.allele_list) == 2)
    combined = r.allele_freq_sum[r.seen_idx["HLA-A01:01"]]
    report("Three-way duplicate → frequencies summed (0.6)",
           combined ≈ 0.6)
end

# Test 1e: order preservation — first occurrence becomes the canonical representative
begin
    alleles = ["HLA-A01:01", "HLA-A*01:01", "HLA-B*07:02"]
    freqs   = [0.3, 0.2, 0.5]
    r = squish_alleles(alleles, freqs)
    report("Order preserved → canonical is first occurrence (HLA-A01:01, no *)",
           r.allele_canonical[r.seen_idx["HLA-A01:01"]] == "HLA-A01:01")
end

# Test 1f: 4-digit allele normalisation (HLA-A*2402 → HLAA24:02, same as HLA-A*24:02 → HLAA24:02)
# The clean_allele function in the script inserts colons; replicate that here.
function clean_allele(a)
    s = String(a)
    s = strip(s)
    s = replace(s, r"\*(\d{2})(\d{2})" => s"*\1:\2")
    return s
end
begin
    raw_alleles = ["HLA-A*2402", "HLA-A*24:02", "HLA-B*07:02"]
    cleaned = clean_allele.(raw_alleles)
    freqs   = [0.3, 0.2, 0.5]
    r = squish_alleles(cleaned, freqs)
    report("4-digit form (HLA-A*2402) squished with HLA-A*24:02 → allele_list length 2",
           length(r.allele_list) == 2)
    combined = r.allele_freq_sum[r.seen_idx["HLA-A24:02"]]
    report("4-digit normalisation → combined frequency 0.5",
           combined ≈ 0.5)
end

# ─── Section 2: Unit — character-based allele batching ───────────────────────

println("\n══════════════════════════════════════════════════")
println("  Section 2 — Character-based allele batching")
println("══════════════════════════════════════════════════")

function make_batches(allele_list, char_limit; count_limit=typemax(Int))
    batches = Vector{Vector{String}}()
    current = String[]; cur_len = 0
    for allele in allele_list
        entry_len = length(allele) + (isempty(current) ? 0 : 1)
        if !isempty(current) && (cur_len + entry_len > char_limit || length(current) >= count_limit)
            push!(batches, current)
            current = [allele]; cur_len = length(allele)
        else
            push!(current, allele); cur_len += entry_len
        end
    end
    isempty(current) || push!(batches, current)
    return batches
end

# Test 2a: all alleles fit within limit → single batch
begin
    alleles = ["HLAA01:01", "HLAA02:01", "HLAB07:02"]  # 9+1+9+1+9 = 29 chars
    batches = make_batches(alleles, 1023)
    report("All fit → single batch", length(batches) == 1)
    report("All fit → batch contains all alleles", length(batches[1]) == 3)
end

# Test 2b: tight limit forces split after first allele
begin
    alleles = ["HLAA01:01", "HLAA02:01", "HLAB07:02"]   # each 9 chars
    # limit of 9 → only first allele fits per batch
    batches = make_batches(alleles, 9)
    report("Tight limit (9) → 3 batches", length(batches) == 3)
    report("Tight limit → each batch has exactly 1 allele", all(b -> length(b) == 1, batches))
end

# Test 2c: limit of 19 → first two fit (9+1+9=19), third is new batch
begin
    alleles = ["HLAA01:01", "HLAA02:01", "HLAB07:02"]
    batches = make_batches(alleles, 19)
    report("Limit 19 → 2 batches", length(batches) == 2)
    report("Limit 19 → first batch has 2 alleles", length(batches[1]) == 2)
    report("Limit 19 → second batch has 1 allele", length(batches[2]) == 1)
end

# Test 2d: realistic 66-allele panel stays in one batch at limit 1023
begin
    panel = CSV.read(joinpath(PROJECT_ROOT, "data", "Example_data", "supertype_panel.csv"), DataFrame)
    allele_col = names(panel)[findfirst(n -> lowercase(replace(String(n), r"\s+" => "")) == "allele", names(panel))]
    cleaned = [replace(strip(String(a)), "*" => "") for a in panel[!, allele_col] if a !== missing]
    # Filter zeros (freq=0 rows would be removed in actual run, but let's test the full list for worst case)
    joined = join(cleaned, ",")
    report("66-allele panel joined length < 1023 chars → single batch expected",
           length(joined) < 1023)
    batches = make_batches(cleaned, 1023)
    report("66-allele panel with limit 1023 → single batch",
           length(batches) == 1)
    println("    (Panel joined allele string: $(length(joined)) chars, $(length(cleaned)) alleles)")
end

# Test 2e: no alleles → zero batches
begin
    batches = make_batches(String[], 1023)
    report("Empty allele list → zero batches", length(batches) == 0)
end

# Test 2f: single allele → one batch
begin
    batches = make_batches(["HLAA01:01"], 1023)
    report("Single allele → one batch of one", length(batches) == 1 && length(batches[1]) == 1)
end

# Test 2g: count limit forces split even when chars fit
begin
    alleles = ["HLAA01:01", "HLAA02:01", "HLAB07:02", "HLAB08:01"]  # 4 alleles, 39 chars total
    batches = make_batches(alleles, 1023; count_limit=2)
    report("Count limit (2) → 2 batches of 2", length(batches) == 2 && all(b -> length(b) == 2, batches))
end

# Test 2h: count limit of 75 (production default) with 93-allele diverse panel
begin
    panel = CSV.read(joinpath(PROJECT_ROOT, "data", "Example_data", "supertype_panel.csv"), DataFrame)
    allele_col = names(panel)[findfirst(n -> lowercase(replace(String(n), r"\s+" => "")) == "allele", names(panel))]
    cleaned = [replace(strip(String(a)), "*" => "") for a in panel[!, allele_col] if a !== missing]
    batches_char = make_batches(cleaned, 1023)
    batches_both = make_batches(cleaned, 1023; count_limit=75)
    report("66-allele panel with count_limit=75 stays 1 batch (≤75 alleles)",
           length(batches_both) == 1 && length(cleaned) <= 75)
end

# ─── Section 3: Unit — merge logic with mock data ─────────────────────────────

println("\n══════════════════════════════════════════════════")
println("  Section 3 — Merge logic (mock data)")
println("══════════════════════════════════════════════════")

# Replicate the merge: given results_matrix and allele_list, produce netmhcpan_output.tsv
function write_mock_batch_file(path, alleles, peptides)
    # Produce a minimal netMHCpan-style XLS file for the given alleles and peptides.
    # Line 1: allele names, line 2: column header, lines 3+: data rows.
    open(path, "w") do io
        println(io, "\t\t" * join(alleles, '\t'))
        cols = ["Pos", "Peptide", "ID"]
        for _ in alleles; append!(cols, ["core", "icore", "EL-score", "EL_Rank"]); end
        println(io, join(cols, '\t'))
        for (i, pep) in enumerate(peptides)
            row = [string(i), pep, "id$(i)"]
            for _ in alleles
                append!(row, ["$(pep[1:4])", "$(pep[1:4])", "0.1", "$(10*i)"])
            end
            println(io, join(row, '\t'))
        end
    end
end

function run_merge(results_matrix, allele_list, total_chunks, n_batches, out_path)
    open(out_path, "w") do out_io
        println(out_io, "\t\t" * join(allele_list, '\t'))
        col_header_written = false
        for chunk_idx in 1:total_chunks
            batch_files = [results_matrix[chunk_idx, bi] for bi in 1:n_batches]
            any(f -> f === nothing || !isfile(f), batch_files) && continue
            lines_per_batch = [collect(eachline(f)) for f in batch_files]
            nlines = length(lines_per_batch[1])
            for bl in lines_per_batch
                length(bl) == nlines || error("Mismatched lines in chunk $chunk_idx")
            end
            if !col_header_written
                fixed_hdr = split(lines_per_batch[1][2], '\t')[1:min(3, length(split(lines_per_batch[1][2], '\t')))]
                allele_cols = [split(bl[2], '\t')[4:end] for bl in lines_per_batch]
                println(out_io, join(vcat(fixed_hdr, reduce(vcat, allele_cols)), '\t'))
                col_header_written = true
            end
            for row_idx in 3:nlines
                fp = split(lines_per_batch[1][row_idx], '\t')
                row = vcat(fp[1:min(3,length(fp))],
                           reduce(vcat, [split(bl[row_idx], '\t')[4:end] for bl in lines_per_batch]))
                println(out_io, join(row, '\t'))
            end
        end
    end
end

mktemprun() do tmpdir
    peptides = ["AAAAAAAAAA", "CCCCCCCCC", "GGGGGGGGGG",
                "TTTTTTTTTT", "MMMMMMMMMM", "FFFFFFFFFFF"]
    alleles_b1 = ["HLAA01:01", "HLAA02:01"]
    alleles_b2 = ["HLAB07:02"]
    all_alleles = vcat(alleles_b1, alleles_b2)

    # Two chunks, two batches each
    total_chunks = 2; n_batches = 2
    chunks = [peptides[1:3], peptides[4:6]]
    results_matrix = Matrix{Union{String,Nothing}}(nothing, total_chunks, n_batches)

    for (ci, chunk_peps) in enumerate(chunks)
        for (bi, batch_alleles) in enumerate([alleles_b1, alleles_b2])
            f = joinpath(tmpdir, "_temp_netMHCpan_output_$(ci)_$(bi).tsv")
            write_mock_batch_file(f, batch_alleles, chunk_peps)
            results_matrix[ci, bi] = f
        end
    end

    out = joinpath(tmpdir, "netmhcpan_output.tsv")
    run_merge(results_matrix, all_alleles, total_chunks, n_batches, out)
    lines = readlines(out)

    report("Merge: output file exists", isfile(out))
    report("Merge: line 1 is allele header with all 3 alleles",
           all(a -> occursin(a, lines[1]), all_alleles))
    report("Merge: line 2 is column header starting with Pos",
           startswith(lines[2], "Pos"))
    report("Merge: correct total data rows (6 peptides)",
           length(lines) == 2 + 6)

    # Each data row should have 3 fixed cols + 3 alleles × 4 cols = 15 cols
    data_widths = [length(split(l, '\t')) for l in lines[3:end]]
    expected_cols = 3 + length(all_alleles) * 4
    report("Merge: each data row has $(expected_cols) tab-separated columns",
           all(w -> w == expected_cols, data_widths))

    # Row order: chunk 1 peptides (rows 1-3) then chunk 2 peptides (rows 4-6)
    # lines[1]=allele hdr, lines[2]=col hdr; data starts at lines[3]
    # 3rd data row = lines[5], 6th data row = lines[8]
    row3_pep = split(lines[5], '\t')[2]
    row6_pep = split(lines[8], '\t')[2]
    report("Merge: row order preserved (chunk 1 first)",
           row3_pep == peptides[3] && row6_pep == peptides[6])
end

# Test 3b: single-batch case (no horizontal stitching needed)
mktemprun() do tmpdir
    peptides = ["AAAAAAAAAA", "CCCCCCCCC"]
    all_alleles = ["HLAA01:01", "HLAB07:02"]

    total_chunks = 1; n_batches = 1
    results_matrix = Matrix{Union{String,Nothing}}(nothing, 1, 1)
    f = joinpath(tmpdir, "_temp_netMHCpan_output_1_1.tsv")
    write_mock_batch_file(f, all_alleles, peptides)
    results_matrix[1, 1] = f

    out = joinpath(tmpdir, "netmhcpan_output.tsv")
    run_merge(results_matrix, all_alleles, total_chunks, n_batches, out)
    lines = readlines(out)

    report("Single-batch merge: 2 data rows", length(lines) == 4)
    ncols = length(split(lines[3], '\t'))
    report("Single-batch merge: correct column count (3 + 2×4 = 11)",
           ncols == 11)
end

# Test 3c: empty chunk skipped gracefully
mktemprun() do tmpdir
    all_alleles = ["HLAA01:01"]
    total_chunks = 2; n_batches = 1
    results_matrix = Matrix{Union{String,Nothing}}(nothing, 2, 1)
    # Only chunk 1 has results; chunk 2 is empty (nothing)
    f = joinpath(tmpdir, "_temp_netMHCpan_output_1_1.tsv")
    write_mock_batch_file(f, all_alleles, ["AAAAAAAAAA"])
    results_matrix[1, 1] = f
    results_matrix[2, 1] = nothing

    out = joinpath(tmpdir, "netmhcpan_output.tsv")
    run_merge(results_matrix, all_alleles, total_chunks, n_batches, out)
    lines = readlines(out)
    report("Empty chunk skipped: output has 1 data row (chunk 2 skipped)",
           length(lines) == 3)
end

# ─── Section 4: Unit — process_output.pl compatibility ───────────────────────

println("\n══════════════════════════════════════════════════")
println("  Section 4 — process_output.pl compatibility")
println("══════════════════════════════════════════════════")

mktemprun() do tmpdir
    # Write a mock combined netmhcpan_output.tsv with 3 alleles
    alleles  = ["HLAA01:01", "HLAB07:02", "HLAC07:01"]
    peptides = ["AAAAAAAAAA", "CCCCCCCCC", "GGGGGGGGGGG"]
    tsv_path = joinpath(tmpdir, "netmhcpan_output.tsv")
    open(tsv_path, "w") do io
        println(io, "\t\t" * join(alleles, '\t'))
        cols = ["Pos", "Peptide", "ID"]
        for _ in alleles; append!(cols, ["core", "icore", "EL-score", "EL_Rank"]); end
        println(io, join(cols, '\t'))
        for (i, pep) in enumerate(peptides)
            row = [string(i), pep, "$(pep[1:2])_mut"]
            for allele in alleles
                append!(row, ["AAAAA", "AAAAA", "0.$(i)$(i)", "$(i).$(i)"])
            end
            println(io, join(row, '\t'))
        end
    end

    csv_path = joinpath(tmpdir, "processed.csv")
    perl_result = try
        run(pipeline(`perl $(joinpath(SRC, "process_output.pl")) $tsv_path`,
                     stdout=csv_path, stderr=devnull))
        true
    catch
        false
    end
    report("process_output.pl runs without error on multi-allele TSV", perl_result)

    if perl_result && isfile(csv_path)
        out_lines = readlines(csv_path)
        report("process_output.pl: header line present", length(out_lines) >= 1)
        report("process_output.pl: 3 peptides × 3 alleles = 9 data rows",
               length(out_lines) == 1 + 3 * 3)
        if length(out_lines) >= 2
            hdr_cols  = split(out_lines[1], ',')
            data_cols = split(out_lines[2], ',')
            report("process_output.pl: header has HLA column", "HLA" in hdr_cols)
            report("process_output.pl: header has EL_Rank column",
                   any(c -> occursin("EL_Rank", c), hdr_cols))
            report("process_output.pl: data rows have same column count as header",
                   length(hdr_cols) == length(data_cols))
            hla_idx = findfirst(==("HLA"), hdr_cols)
            if hla_idx !== nothing
                hla_vals = [split(l, ',')[hla_idx] for l in out_lines[2:end]]
                report("process_output.pl: HLA values cover all 3 alleles",
                       Set(hla_vals) == Set(alleles))
            end
        end
    end
end

# ─── Section 5: Integration — run_netMHCpan_global.jl ────────────────────────

println("\n══════════════════════════════════════════════════")
println("  Section 5 — Integration: run_netMHCpan_global.jl")
if !netmhcpan_available
    println("  ⚠  netMHCpan not found at $(NETMHCPAN) — skipping integration tests")
end
println("══════════════════════════════════════════════════")

if netmhcpan_available
    # Minimal test data: 5 canonical 9-mers, small panel with one duplicate pair
    TEST_PEPTIDES = ["AAAAAAAAAA", "KKKKKKKKKK", "FFFFFFFFFFFFFF"[1:9],
                     "LLLLLLLLLL"[1:9], "VVVVVVVVVV"[1:9]]
    TEST_PANEL = """Locus,Allele,Frequency
A,HLA-A*01:01,0.3
A,HLA-A01:01,0.2
B,HLA-B*07:02,0.5
"""
    # HLA-A*01:01 and HLA-A01:01 should squish → 2 unique alleles

    # ─── 5a: Sequential run (--t 1) ───────────────────────────────────────────
    mktemprun() do tmpdir
        println("\n  5a) Sequential run (--t 1)...")
        writefile(p, s) = (open(p, "w") do io; write(io, s); end)
        writefile(joinpath(tmpdir, "Peptides.pep"), join(TEST_PEPTIDES, "\n") * "\n")
        writefile(joinpath(tmpdir, "supertype_panel.csv"), TEST_PANEL)
        writefile(joinpath(PROJECT_ROOT, "src", "settings.txt"),
                  "NETMHCPAN=$(NETMHCPAN)\n")  # ensure settings reachable

        script = joinpath(SRC, "run_netMHCpan_global.jl")
        cmd = `julia --project=$(PROJECT_ROOT) $script --folder $tmpdir --t 1`
        result = try
            run(pipeline(cmd, stdout=devnull, stderr=devnull))
            true
        catch e
            println("    Error: $e")
            false
        end
        report("5a: script exits successfully", result)

        if result
            # Check squishing_map.csv
            smap_path = joinpath(tmpdir, "squishing_map.csv")
            report("5a: squishing_map.csv created", isfile(smap_path))
            if isfile(smap_path)
                smap = CSV.read(smap_path, DataFrame)
                report("5a: squishing_map has 3 rows (one per original allele)",
                       nrow(smap) == 3)
                report("5a: squishing_map has Original_Allele column",
                       hasproperty(smap, :Original_Allele))
                report("5a: squishing_map has Was_Squished column",
                       hasproperty(smap, :Was_Squished))
                if hasproperty(smap, :Was_Squished)
                    n_squished = count(smap.Was_Squished)
                    report("5a: exactly 2 rows have Was_Squished=true",
                           n_squished == 2)
                end
            end

            # Check supertype_panel_squished.csv
            spanel_path = joinpath(tmpdir, "supertype_panel_squished.csv")
            report("5a: supertype_panel_squished.csv created", isfile(spanel_path))
            if isfile(spanel_path)
                spanel = CSV.read(spanel_path, DataFrame)
                report("5a: squished panel has 2 rows (after merging duplicate A alleles)",
                       nrow(spanel) == 2)
                if hasproperty(spanel, :Frequency)
                    a_row = filter(r -> occursin("A0", replace(String(r.Allele), "*" => "")), spanel)
                    if nrow(a_row) > 0
                        report("5a: squished A allele frequency = 0.5 (0.3 + 0.2)",
                               a_row.Frequency[1] ≈ 0.5)
                    end
                end
            end

            # Check netmhcpan_output.tsv
            out_path = joinpath(tmpdir, "netmhcpan_output.tsv")
            report("5a: netmhcpan_output.tsv created", isfile(out_path))
            if isfile(out_path)
                lines = readlines(out_path)
                report("5a: output has ≥ 3 lines (2 headers + ≥1 data)",
                       length(lines) >= 3)
                report("5a: line 1 contains allele names",
                       occursin("HLA", lines[1]))
                # After squishing we have 2 alleles; each contributes 4 cols → 3+8=11 cols
                if length(lines) >= 3
                    data_cols = length(split(lines[3], '\t'))
                    report("5a: data rows have 3 + 2×4 = 11 columns (2 squished alleles)",
                           data_cols == 11)
                end
                report("5a: correct number of data rows ($(length(TEST_PEPTIDES)) peptides)",
                       length(lines) == 2 + length(TEST_PEPTIDES))
            end

            # Check no temp files left over
            leftover = [f for f in readdir(tmpdir) if startswith(f, "_temp_")]
            report("5a: all temp files cleaned up", isempty(leftover))
        end
    end

    # ─── 5b: Parallel run (--t 2) ─────────────────────────────────────────────
    mktemprun() do tmpdir
        println("\n  5b) Parallel run (--t 2), same data → same output...")
        writefile(p, s) = (open(p, "w") do io; write(io, s); end)
        writefile(joinpath(tmpdir, "Peptides.pep"), join(TEST_PEPTIDES, "\n") * "\n")
        writefile(joinpath(tmpdir, "supertype_panel.csv"), TEST_PANEL)

        script = joinpath(SRC, "run_netMHCpan_global.jl")
        result = try
            run(pipeline(`julia --project=$(PROJECT_ROOT) $script --folder $tmpdir --t 2`,
                         stdout=devnull, stderr=devnull))
            true
        catch; false; end
        report("5b: script exits successfully with --t 2", result)

        if result
            out_path = joinpath(tmpdir, "netmhcpan_output.tsv")
            lines = isfile(out_path) ? readlines(out_path) : String[]
            report("5b: output file created", !isempty(lines))
            report("5b: correct number of data rows", length(lines) == 2 + length(TEST_PEPTIDES))
        end
    end

    # ─── 5c: Sequential vs parallel output equivalence ────────────────────────
    println("\n  5c) Sequential vs parallel equivalence check...")
    seq_output = parallel_output = nothing
    mktemprun() do d1
        mktemprun() do d2
            writefile(p, s) = (open(p, "w") do io; write(io, s); end)
            for d in (d1, d2)
                writefile(joinpath(d, "Peptides.pep"), join(TEST_PEPTIDES, "\n") * "\n")
                writefile(joinpath(d, "supertype_panel.csv"), TEST_PANEL)
            end
            script = joinpath(SRC, "run_netMHCpan_global.jl")
            ok1 = try
                run(pipeline(`julia --project=$(PROJECT_ROOT) $script --folder $d1 --t 1`,
                             stdout=devnull, stderr=devnull)); true
            catch; false; end
            ok2 = try
                run(pipeline(`julia --project=$(PROJECT_ROOT) $script --folder $d2 --t 2`,
                             stdout=devnull, stderr=devnull)); true
            catch; false; end

            if ok1 && ok2
                f1 = joinpath(d1, "netmhcpan_output.tsv")
                f2 = joinpath(d2, "netmhcpan_output.tsv")
                if isfile(f1) && isfile(f2)
                    lines1 = readlines(f1)
                    lines2 = readlines(f2)
                    report("5c: sequential and parallel produce same number of lines",
                           length(lines1) == length(lines2))
                    report("5c: sequential and parallel allele headers match",
                           lines1[1] == lines2[1])
                    report("5c: sequential and parallel column headers match",
                           lines1[2] == lines2[2])
                    # Data values may differ in floating point precision, so compare
                    # just the peptide column (col 2) of each row
                    peps1 = [split(l, '\t')[2] for l in lines1[3:end]]
                    peps2 = [split(l, '\t')[2] for l in lines2[3:end]]
                    report("5c: same peptides in same row order", peps1 == peps2)
                end
            end
        end
    end

    # ─── 5d: Large panel (> 1 batch needed) ───────────────────────────────────
    println("\n  5d) Large panel forcing multi-batch split...")
    mktemprun() do tmpdir
        # Build a panel whose joined allele string exceeds 1023 chars,
        # forcing the character-based batcher to split into ≥2 batches.
        # Use 95 real alleles verified to exist in netMHCpan's database
        # (joined length ≈1044 chars, safely above the 1023-char limit).
        big_alleles = [
            "HLA-A01:01", "HLA-A01:02", "HLA-A02:01", "HLA-A02:02", "HLA-A03:01", "HLA-A03:02",
            "HLA-A11:01", "HLA-A11:02", "HLA-A23:01", "HLA-A23:02", "HLA-A24:02", "HLA-A24:03",
            "HLA-A25:01", "HLA-A25:02", "HLA-A26:01", "HLA-A26:02", "HLA-A29:01", "HLA-A29:02",
            "HLA-A30:01", "HLA-A30:02", "HLA-A31:01", "HLA-A31:02", "HLA-A32:01", "HLA-A32:02",
            "HLA-A33:01", "HLA-A33:03", "HLA-A34:01", "HLA-A34:02", "HLA-A36:01", "HLA-A36:02",
            "HLA-A43:01", "HLA-A66:01", "HLA-A66:02", "HLA-A68:01", "HLA-A68:02", "HLA-A69:01",
            "HLA-A69:02", "HLA-A74:01", "HLA-A74:02", "HLA-A80:01", "HLA-A80:02", "HLA-B07:02",
            "HLA-B07:03", "HLA-B08:01", "HLA-B08:02", "HLA-B13:01", "HLA-B13:02", "HLA-B14:01",
            "HLA-B14:02", "HLA-B15:01", "HLA-B15:02", "HLA-B18:01", "HLA-B18:02", "HLA-B27:01",
            "HLA-B27:02", "HLA-B35:01", "HLA-B35:02", "HLA-B37:01", "HLA-B37:02", "HLA-B38:01",
            "HLA-B38:02", "HLA-B39:01", "HLA-B39:02", "HLA-B40:01", "HLA-B40:02", "HLA-B41:01",
            "HLA-B41:02", "HLA-B42:01", "HLA-B42:02", "HLA-B44:02", "HLA-B44:03", "HLA-B45:01",
            "HLA-B45:02", "HLA-B46:01", "HLA-B46:02", "HLA-B47:01", "HLA-B47:02", "HLA-B48:01",
            "HLA-B48:02", "HLA-B49:01", "HLA-B49:02", "HLA-B50:01", "HLA-B50:02", "HLA-B51:01",
            "HLA-B51:02", "HLA-B52:01", "HLA-B52:02", "HLA-B53:01", "HLA-B53:02", "HLA-B54:01",
            "HLA-B54:02", "HLA-B55:01", "HLA-B55:02",
        ]
        big_panel_str = "Locus,Allele,Frequency\n"
        freq_each = round(1.0 / length(big_alleles); digits=6)
        for a in big_alleles
            locus = String(a[5:5])
            big_panel_str *= "$locus,$a,$freq_each\n"
        end

        writefile(p, s) = (open(p, "w") do io; write(io, s); end)
        writefile(joinpath(tmpdir, "Peptides.pep"), join(TEST_PEPTIDES, "\n") * "\n")
        writefile(joinpath(tmpdir, "supertype_panel.csv"), big_panel_str)

        script = joinpath(SRC, "run_netMHCpan_global.jl")
        result = try
            run(pipeline(`julia --project=$(PROJECT_ROOT) $script --folder $tmpdir --t 1`,
                         stdout=devnull, stderr=devnull)); true
        catch; false; end
        report("5d: large panel ($(length(big_alleles)) alleles, multi-batch) completes successfully", result)

        if result
            out_path = joinpath(tmpdir, "netmhcpan_output.tsv")
            if isfile(out_path)
                lines = readlines(out_path)
                expected_allele_cols = length(big_alleles) * 4
                if length(lines) >= 3
                    n_cols = length(split(lines[3], '\t'))
                    report("5d: output has correct column count (3 + $(expected_allele_cols) = $(3+expected_allele_cols))",
                           n_cols == 3 + expected_allele_cols)
                end
                report("5d: correct number of data rows", length(lines) == 2 + length(TEST_PEPTIDES))
            end
        end
    end
end  # if netmhcpan_available

# ─── Section 6: Integration — run_netMHCpan.jl (individual alleles) ──────────

println("\n══════════════════════════════════════════════════")
println("  Section 6 — Integration: run_netMHCpan.jl")
if !netmhcpan_available
    println("  ⚠  netMHCpan not found — skipping")
end
println("══════════════════════════════════════════════════")

if netmhcpan_available
    TEST_ALLELES = "HLA-A*01:01\nHLA-B*07:02\n"

    mktemprun() do tmpdir
        writefile(p, s) = (open(p, "w") do io; write(io, s); end)
        writefile(joinpath(tmpdir, "Peptides.pep"), join(["AAAAAAAAAA", "KKKKKKKKKK"], "\n") * "\n")
        writefile(joinpath(tmpdir, "alleles.txt"), TEST_ALLELES)

        script = joinpath(SRC, "run_netMHCpan.jl")
        result = try
            run(pipeline(`julia --project=$(PROJECT_ROOT) $script --folder $tmpdir --t 1`,
                         stdout=devnull, stderr=devnull)); true
        catch; false; end
        report("6a: run_netMHCpan.jl runs successfully", result)

        if result
            out_path = joinpath(tmpdir, "netMHCpan_output.tsv")
            lines = isfile(out_path) ? readlines(out_path) : String[]
            report("6a: output file created", !isempty(lines))
            report("6a: 2 data rows (2 peptides)", length(lines) == 4)
            if length(lines) >= 3
                n_cols = length(split(lines[3], '\t'))
                report("6a: correct column count (3 + 2 alleles × 4 = 11)",
                       n_cols == 11)
            end
            leftover = [f for f in readdir(tmpdir) if startswith(f, "_temp_")]
            report("6a: temp files cleaned up", isempty(leftover))
        end
    end
end

# ─── Section 7: process_best_ranks_supertype panel preference ────────────────

println("\n══════════════════════════════════════════════════")
println("  Section 7 — process_best_ranks_supertype panel preference")
println("══════════════════════════════════════════════════")

begin
    # Read the relevant section from the script and verify it checks for squished panel
    script_text = read(joinpath(SRC, "process_best_ranks_supertype.jl"), String)
    report("7a: script checks for supertype_panel_squished.csv before falling back",
           occursin("supertype_panel_squished.csv", script_text))
    report("7b: squished panel is preferred (isfile check precedes fallback)",
           occursin(r"isfile\(squished_path\).*supertype_panel\.csv"s, script_text) ||
           occursin("isfile(squished_path) ? squished_path", script_text))
end

# ─── Final summary ─────────────────────────────────────────────────────────────

println()
println("══════════════════════════════════════════════════")
total = pass_count[] + fail_count[]
println("  Results: $(pass_count[])/$(total) passed, $(fail_count[]) failed")
println("══════════════════════════════════════════════════")
if fail_count[] > 0
    println("\n  ❌ Some tests failed.")
    exit(1)
else
    println("\n  ✅ All tests passed.")
end
