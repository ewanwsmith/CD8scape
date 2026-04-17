#!/usr/bin/env julia
"""
process_best_ranks.jl

Processes best ranks for ancestral and derived peptides for each locus.
Computes harmonic mean best rank (HMBR) across all panel alleles and the
resulting log2 fold change (derived / ancestral) per variant.

Usage:
    julia process_best_ranks.jl <folder_path> [--suffix <name>] [--latest|--no-latest] [--per-allele]

Arguments:
    <folder_path>   Path to the folder containing processed_peptides.csv.

Options:
    --suffix <name>     Suffix appended to input/output filenames.
    --latest            Use the most recently modified input file when ambiguous (default).
    --no-latest         Error on ambiguous input files instead.
    --per-allele        Perform the log2 fold change calculation per allele in the provided
                        genome, writing a separate per_allele_best_ranks.csv file containing one
                        row per (Locus, Mutation, allele) with columns:
                          Frame              - protein / region label
                          MHC                - the allele identifier
                          ELBR_A            - best ancestral EL rank for that allele
                          ELBR_D            - best derived EL rank for that allele
                          foldchange_BR      - EL_Rank_D / EL_Rank_A
                          log2_foldchange_BR - log2(foldchange_BR)
                        Only alleles where ancestral EL rank ≤ 2% are included.
"""

using DataFrames, CSV, Statistics, StatsBase
include("path_utils.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    # Parse positional folder and optional flags
    if length(ARGS) < 1
        println("Usage: julia process_best_ranks.jl <folder_path> [--suffix <name>] [--latest]")
        exit(1)
    end

    folder_path = ARGS[1]
    # Helper to parse flags from the remainder of ARGS to avoid soft-scope warnings
    function _parse_suffix_latest(argv::Vector{String})
        sfx = ""
        lat = true
        me = false
        i = 1
        while i <= length(argv)
            a = argv[i]
            if a == "--suffix"
                if i + 1 <= length(argv) && !startswith(argv[i+1], "--")
                    i += 1
                    sfx = argv[i]
                end
            elseif a == "--latest"
                lat = true
            elseif a == "--no-latest"
                lat = false
            elseif a == "--per-allele"
                me = true
            end
            i += 1
        end
        return sfx, lat, me
    end
    suffix, latest, per_allele = _parse_suffix_latest(ARGS[2:end])

    input_file = resolve_read(joinpath(folder_path, "processed_peptides.csv"); suffix=suffix, latest=latest)

    # -------------------------------------------------------------------------
    # Memory-safe, fast streaming read
    # -------------------------------------------------------------------------
    # The processed_peptides.csv file can grow to hundreds of millions of rows
    # for supertype (panel) runs over large segments (e.g. ~530M rows for HA
    # simulated on the 2078-allele human panel). Loading the whole file into a
    # DataFrame requires 20–30 GB of RAM and OOMs.
    #
    # Instead, we stream the CSV row-by-row with `CSV.Rows` and accumulate the
    # running minimum `EL_Rank` per `(Locus, MHC, Mutation)` group for the
    # ancestral (`_A`) and derived (`_D`) peptide types simultaneously. Peak
    # memory is proportional to the number of unique groups (~millions) rather
    # than the number of rows (~hundreds of millions).
    #
    # Performance notes — critical to avoid the ~20× slowdown of naive
    # `CSV.Rows`:
    #   * `types=...` lets CSV.jl parse Locus/EL_Rank as Int/Float64 natively,
    #     skipping per-row String→number allocations.
    #   * `reusebuffer=true` reuses the row object across iterations.
    #   * `parse_label` is a hand-written parser (no regex, no Regex.replace
    #     per row) that recovers `(mutation, frame)` from `Peptide_label`.
    #   * MHC / Mutation / Frame are interned through a shared String pool
    #     so the dict's string memory stays O(unique values), not O(groups).
    #
    # Dict value layout: (best_EL_Rank::Float64, frame::String, sequence::String)
    println("Reading input file: $input_file")
    println("Streaming CSV to compute per-(Locus,MHC,Mutation) best ranks...")

    GroupKey  = Tuple{Int,String,String}            # (Locus, MHC, Mutation)
    GroupVal  = Tuple{Float64,String,String}        # (best_EL_Rank, Frame, Sequence)
    best_A_dict = Dict{GroupKey, GroupVal}()
    best_D_dict = Dict{GroupKey, GroupVal}()

    # String interning pool — many rows share the same MHC/Mutation/Frame.
    # Storing a single canonical `String` per distinct value caps memory.
    string_pool = Dict{String,String}()
    @inline intern(s::AbstractString) = begin
        k = String(s)
        get!(string_pool, k, k)
    end

    # Fast, allocation-light parser for `Peptide_label`. Expected shape:
    #     "<mutation>_<...>_<frame>_<digits>_<A|C|D|V>"
    # Returns (mutation, frame) as `SubString` views into `label` when
    # possible, avoiding an intermediate `String` allocation.
    @inline function parse_label(label::AbstractString)
        L = lastindex(label)
        # Strip trailing "_[A|C|D|V]"
        if L >= 2 && label[prevind(label, L)] == '_'
            c = label[L]
            if c == 'A' || c == 'C' || c == 'D' || c == 'V'
                L = prevind(label, L, 2)
            end
        end
        # Strip trailing "_<digits>"
        j = L
        seen_digit = false
        while j > 0 && isdigit(label[j])
            seen_digit = true
            j = prevind(label, j)
        end
        if seen_digit && j > 0 && label[j] == '_'
            L = prevind(label, j)
        end
        # mutation = prefix before first '_'
        p1 = findnext(==('_'), label, firstindex(label))
        mut = p1 === nothing || p1 > L ? SubString(label, firstindex(label), L) :
                                         SubString(label, firstindex(label), prevind(label, p1))
        # frame = suffix after last '_' in [1:L]
        p2 = findprev(==('_'), label, L)
        frm = p2 === nothing ? SubString(label, firstindex(label), L) :
                               SubString(label, nextind(label, p2), L)
        return mut, frm
    end

    # Specify types so CSV.jl parses natively (big speedup vs default lazy
    # strings-only). We still tolerate missing/NA via `missingstring=...`.
    col_types = Dict(
        :Locus         => Int,
        :EL_Rank       => Float64,
        :Peptide_label => String,
        :MHC           => String,
        :Peptide       => String,
    )

    row_count = 0
    processed_count = 0
    try
        for row in CSV.Rows(input_file;
                            reusebuffer   = true,
                            types         = col_types,
                            missingstring = ["", "NA"])
            row_count += 1

            # --- Classify peptide type by label suffix (ancestral / derived) ----
            plabel = row.Peptide_label
            plabel === missing && continue
            is_A = endswith(plabel, "_A")
            is_D = endswith(plabel, "_D")
            (is_A || is_D) || continue

            # --- Usable numeric fields? -----------------------------------------
            er = row.EL_Rank
            er === missing && continue
            (isnan(er) || er <= 0) && continue
            loc = row.Locus
            loc === missing && continue

            # --- Derive mutation + frame (no regex, no String allocs here) ------
            mut_sv, frm_sv = parse_label(plabel)
            mutation    = intern(mut_sv)
            frame_label = intern(frm_sv)
            mhc         = intern(row.MHC)

            # Peptide sequences are highly variable → don't intern (pool bloat).
            sequence = String(row.Peptide)

            key    = (loc, mhc, mutation)
            target = is_A ? best_A_dict : best_D_dict
            existing = get(target, key, nothing)
            if existing === nothing || er < existing[1]
                target[key] = (er, frame_label, sequence)
            end

            processed_count += 1
            if row_count % 10_000_000 == 0
                println("  $(row_count ÷ 1_000_000)M rows scanned "
                        * "($(processed_count) retained, "
                        * "$(length(best_A_dict))+$(length(best_D_dict)) groups)")
            end
        end
    catch e
        println("Error streaming input file: $e")
        exit(1)
    end
    println("Successfully processed $(row_count) rows "
            * "($(processed_count) passed filters, "
            * "$(length(best_A_dict)) ancestral groups, "
            * "$(length(best_D_dict)) derived groups)")

    # -------------------------------------------------------------------------
    # Convert accumulated dicts into DataFrames matching the downstream schema.
    # Memory here is bounded by the number of unique groups (already reduced).
    # -------------------------------------------------------------------------
    function dict_to_df(d::Dict, peptide_type::String)
        n = length(d)
        Loci  = Vector{Int}(undef, n);       MHCs      = Vector{String}(undef, n)
        Muts  = Vector{String}(undef, n);    Ranks     = Vector{Float64}(undef, n)
        Types = fill(peptide_type, n);       Frames    = Vector{String}(undef, n)
        Seqs  = Vector{String}(undef, n)
        i = 1
        for ((loc, mhc, mut), (rank, frm, seq)) in d
            Loci[i]   = loc;   MHCs[i]   = mhc;   Muts[i]  = mut
            Ranks[i]  = rank;  Frames[i] = frm;   Seqs[i]  = seq
            i += 1
        end
        return DataFrame(
            Locus        = Loci,
            MHC          = MHCs,
            Mutation     = Muts,
            Best_EL_Rank = Ranks,
            Peptide_Type = Types,
            Frame        = Frames,
            Sequence     = Seqs,
        )
    end

    println("Materialising best ranks for ancestral peptides (_A)...")
    best_C = dict_to_df(best_A_dict, "A")
    println("Found best ranks for ancestral peptides: $(nrow(best_C)) entries")

    println("Materialising best ranks for derived peptides (_D)...")
    best_V = dict_to_df(best_D_dict, "D")
    println("Found best ranks for derived peptides: $(nrow(best_V)) entries")

    # Free the accumulator dicts early; everything we need now lives in the
    # small-ish best_C / best_V DataFrames. Help the GC reclaim the pools.
    empty!(best_A_dict); empty!(best_D_dict); empty!(string_pool)
    GC.gc()

# Combine ancestral and derived results
    best_ranks = vcat(best_C, best_V)

# --- Map Locus to protein Description from frames.csv ---
    # Context:
    # Previously, `Frame` was parsed from `Peptide_label` by taking the
    # last underscore-delimited token, which truncated protein names.
    # We now replace `Frame` using frames.csv by matching each row's
    # `Locus` to the Region ranges and taking the corresponding full
    # `Description`. If no match is found, we keep the original value.
    frames_file = resolve_read(joinpath(folder_path, "frames.csv"); suffix="", latest=latest)
    if isfile(frames_file)
        frames_df = CSV.read(frames_file, DataFrame)

        # Parse a Region string like "77,496" or multiple coords "77,496;606,980"
        # into a vector of (start, end) tuples
        function region_to_segments(region_str::AbstractString)::Vector{Tuple{Int,Int}}
            rs = strip(String(region_str))
            isempty(rs) && return Tuple{Int,Int}[]
            segments = Tuple{Int,Int}[]
            for part in split(rs, ';')
                se = split(part, ',')
                if length(se) == 2
                    start_pos = try
                        parse(Int, se[1])
                    catch
                        continue
                    end
                    end_pos = try
                        parse(Int, se[2])
                    catch
                        continue
                    end
                    push!(segments, (start_pos, end_pos))
                end
            end
            return segments
        end

        # Build a lookup table of regions → description
        regions_lookup = [(region_to_segments(fr.Region), String(fr.Description)) for fr in eachrow(frames_df)]

        # Replace Frame (derived from Peptide_label) with full Description matched by Locus.
        # If multiple frames overlap a locus, the first match in frames.csv is used.
        if !isempty(best_ranks)
            best_ranks[!, :Frame] = [
                begin
                    loc = r.Locus
                    # find first frames.csv entry whose any segment contains locus
                    desc = nothing
                    for (segments, d) in regions_lookup
                        for (s,e) in segments
                            if loc >= s && loc <= e
                                desc = d
                                break
                            end
                        end
                        if desc !== nothing
                            break
                        end
                    end
                    desc === nothing ? String(r.Frame) : desc
                end for r in eachrow(best_ranks)
            ]
        end
    else
        println("Warning: frames.csv not found in $folder_path. Protein labels will not be mapped.")
    end

# Reorder: put Frame, Locus, Mutation first
    best_ranks = select(best_ranks, :Frame, :Locus, :Mutation, Not([:Frame, :Locus, :Mutation]))

# Compute one Frame per (Locus, Mutation): prefer current mapped, otherwise first seen
    description_roots = DataFrame(Locus=Int[], Mutation=String[], Frame=String[])
    if !isempty(best_ranks)
        description_roots = combine(groupby(best_ranks, [:Locus, :Mutation])) do sdf
            (; Locus = sdf.Locus[1], Mutation = sdf.Mutation[1], Frame = sdf.Frame[1])
        end
    end

# Save best_ranks.csv
    best_ranks_file = resolve_write(joinpath(folder_path, "best_ranks.csv"); suffix=suffix)
    CSV.write(best_ranks_file, best_ranks)
    println("Saved best ranks to $best_ranks_file")

# Pivot best_ranks to have separate columns for HMBR_A and HMBR_D
    println("Calculating harmonic mean best ranks (HMBR) per Frame/Locus/Mutation...")
    if !isempty(best_ranks)
        # Compute harmonic mean while ignoring NaN/missing/zero values
        function safe_harmmean(xs)
            vals = [float(x) for x in xs if !(ismissing(x) || (x isa AbstractFloat && isnan(x)) || !(x > 0))]
            isempty(vals) && return NaN
            return harmmean(vals)
        end

        pivot_df = unstack(
            combine(groupby(best_ranks, [:Frame, :Locus, :Mutation, :Peptide_Type]),
                    :Best_EL_Rank => safe_harmmean => :HMBR),
            :Peptide_Type, :HMBR)

    # Only rename columns if they exist
    colnames = names(pivot_df)
    rename_pairs = Pair{Symbol,Symbol}[]
    if "A" in colnames
        push!(rename_pairs, Symbol("A") => :HMBR_A)
    end
    if "D" in colnames
        push!(rename_pairs, Symbol("D") => :HMBR_D)
    end
    if !isempty(rename_pairs)
        rename!(pivot_df, rename_pairs...)
    end

    # Warn if no derived data is present
    if !("HMBR_D" in names(pivot_df))
        println("Warning: No derived peptides found. Skipping fold change calculations and HMBR_D output.")
    end

    # Identify and report missing/NaN values before fold change calculation (per Locus, Mutation)
    missing_msgs = Set{Tuple{Int,String,String}}()  # (Locus, Mutation, which)
    for row in eachrow(pivot_df)
        locus = row.Locus
        change = get(row, :Mutation, "?")
        hA = get(row, :HMBR_A, NaN)
        hD = get(row, :HMBR_D, NaN)
        if (ismissing(hA) || (hA isa AbstractFloat && isnan(hA)))
            push!(missing_msgs, (locus, String(change), "ancestral"))
        end
        if (ismissing(hD) || (hD isa AbstractFloat && isnan(hD)))
            push!(missing_msgs, (locus, String(change), "derived"))
        end
    end
    for (locus, change, which) in missing_msgs
        println("Fold change could not be calculated for locus $(locus) (change $(change)) due to missing $(which) rank.")
    end

    # Filter out loci where both HMBR_A and HMBR_D are greater than 2 (non-binding)
    before_filter = nrow(pivot_df)
    pivot_df = filter(row -> begin
            hA = get(row, :HMBR_A, NaN)
            hD = get(row, :HMBR_D, NaN)
            # keep only rows where both HMBR values are finite numbers
            !(ismissing(hA) || (hA isa AbstractFloat && isnan(hA)) || ismissing(hD) || (hD isa AbstractFloat && isnan(hD))) && !(hA > 2 && hD > 2)
        end, pivot_df)
    removed_count = before_filter - nrow(pivot_df)
    println("Removed $removed_count loci where both ancestral and derived states were predicted to be non-binding (HMBR > 2)")

    # Calculate fold change (Derived (_D) / Ancestral (_A)) for all valid rows
    if ("HMBR_A" in names(pivot_df)) && ("HMBR_D" in names(pivot_df))
        pivot_df.foldchange_HMBR = pivot_df.HMBR_D ./ pivot_df.HMBR_A
        pivot_df.log2_foldchange_HMBR = log2.(pivot_df.foldchange_HMBR)
    else
        println("DEBUG: Cannot calculate foldchange_HMBR or log2_foldchange_HMBR due to missing columns.")
    end

    # Frame is already present from grouping; no join needed

    # Keep Description exactly as mapped; ensure no trailing underscores or numeric suffixes
    # No further cleanup needed for Frame

    # Calculate fold change and log2 after all joins/cleaning
    if ("HMBR_A" in names(pivot_df)) && ("HMBR_D" in names(pivot_df))
        pivot_df.foldchange_HMBR = pivot_df.HMBR_D ./ pivot_df.HMBR_A
        pivot_df.log2_foldchange_HMBR = log2.(pivot_df.foldchange_HMBR)
    end

    # Compute per-allele log2 fold change if requested
    if per_allele
        best_A_pa = filter(row -> row.Peptide_Type == "A", best_ranks)
        best_D_pa = filter(row -> row.Peptide_Type == "D", best_ranks)
        if !isempty(best_A_pa) && !isempty(best_D_pa)
            per_allele_df = innerjoin(
                select(best_A_pa, :Frame, :Locus, :MHC, :Mutation, :Best_EL_Rank => :ELBR_A),
                select(best_D_pa, :Locus, :MHC, :Mutation, :Best_EL_Rank => :ELBR_D),
                on = [:Locus, :MHC, :Mutation]
            )
            filter!(row -> row.ELBR_A <= 2.0, per_allele_df)
            if !isempty(per_allele_df)
                per_allele_df[!, :foldchange_BR] = per_allele_df.ELBR_D ./ per_allele_df.ELBR_A
                per_allele_df[!, :log2_foldchange_BR] = log2.(per_allele_df.foldchange_BR)
                select!(per_allele_df, :Frame, :Locus, :Mutation, :MHC, :ELBR_A, :ELBR_D, :foldchange_BR, :log2_foldchange_BR)
                per_allele_file = resolve_write(joinpath(folder_path, "per_allele_best_ranks.csv"); suffix=suffix)
                CSV.write(per_allele_file, per_allele_df)
                println("Saved per-allele escape log2 fold changes to $per_allele_file")
            else
                println("Warning: No alleles with ancestral EL_Rank ≤ 2 found. per_allele_best_ranks.csv will not be written.")
            end
        else
            println("Warning: Missing ancestral or derived peptide data. per_allele_best_ranks.csv will not be written.")
        end
    end

    # Prepare final output columns (use already computed values)
    output_cols = [:Frame, :Locus, :Mutation, :HMBR_A, :HMBR_D, :foldchange_HMBR, :log2_foldchange_HMBR]
    pivot_df = select(pivot_df, output_cols...)

    # Save harmonic mean results with fold change
    harmonic_mean_file = resolve_write(joinpath(folder_path, "harmonic_mean_best_ranks.csv"); suffix=suffix)
    CSV.write(harmonic_mean_file, pivot_df)
    println("Saved harmonic mean best ranks to $harmonic_mean_file")

    # Minimal test: write Locus and log2_foldchange_HMBR to a separate CSV for debugging
    if :log2_foldchange_HMBR in names(pivot_df)
        minimal_test_file = resolve_write(joinpath(folder_path, "harmonic_mean_best_ranks_log2_test.csv"); suffix=suffix)
        minimal_df = select(pivot_df, :Locus, :log2_foldchange_HMBR)
        CSV.write(minimal_test_file, minimal_df)
        println("Minimal test CSV written to $minimal_test_file")
    else
    end
    else
        println("No valid best rank data available. Skipping harmonic mean calculations.")
    end
end