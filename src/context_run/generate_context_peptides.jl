#!/usr/bin/env julia

"""
This script generates context peptides from genomic frames data.

Purpose:
- Reads a CSV file containing genomic frames.
- Samples random loci within these frames.
- Introduces random nucleotide mutations.
- Translates DNA sequences to amino acids.
- Generates peptides centered on mutated amino acid loci.
- Filters peptides containing stop codons or identical consensus/variant sequences.
- Outputs labeled peptides to CSV and .pep files.

Usage:
    julia generate_context_peptides.jl --folder <path_to_folder> [--n_loci <number_of_loci>] [--seed <random_seed>]

Arguments:
    --folder   Path to the folder containing 'frames.csv'.
    --n_loci   Number of loci to simulate (default: 1000).
    --seed     Random number seed (default: 1320).
"""

using ArgParse
using CSV, DataFrames
using Random
using FilePathsBase
# load variant expansion helper
include(joinpath("..", "variant_expansion.jl"))
using .VariantExpansion

# Define codon to amino acid translation dictionary
const CODON_DICT = Dict(
    "ATA" => "I", "ATC" => "I", "ATT" => "I", "ATG" => "M",
    "ACA" => "T", "ACC" => "T", "ACG" => "T", "ACT" => "T",
    "AAC" => "N", "AAT" => "N", "AAA" => "K", "AAG" => "K",
    "AGC" => "S", "AGT" => "S", "AGA" => "R", "AGG" => "R",
    "CTA" => "L", "CTC" => "L", "CTG" => "L", "CTT" => "L",
    "CCA" => "P", "CCC" => "P", "CCG" => "P", "CCT" => "P",
    "CAC" => "H", "CAT" => "H", "CAA" => "Q", "CAG" => "Q",
    "CGA" => "R", "CGC" => "R", "CGG" => "R", "CGT" => "R",
    "GTA" => "V", "GTC" => "V", "GTG" => "V", "GTT" => "V",
    "GCA" => "A", "GCC" => "A", "GCG" => "A", "GCT" => "A",
    "GAC" => "D", "GAT" => "D", "GAA" => "E", "GAG" => "E",
    "GGA" => "G", "GGC" => "G", "GGG" => "G", "GGT" => "G",
    "TCA" => "S", "TCC" => "S", "TCG" => "S", "TCT" => "S",
    "TTC" => "F", "TTT" => "F", "TTA" => "L", "TTG" => "L",
    "TAC" => "Y", "TAT" => "Y", "TAA" => "*", "TAG" => "*",
    "TGC" => "C", "TGT" => "C", "TGA" => "*", "TGG" => "W"
)


# Function to translate DNA sequence to protein sequence
function translate_dna_to_protein(dna_sequence::String)::String
    join([get(CODON_DICT, uppercase(dna_sequence[i:i+2]), "?") 
          for i in 1:3:length(dna_sequence)-2], "")
end

function generate_peptides(sequence::String, aa_locus::Int, lengths::Vector{Int})::Vector{String}
    peptides = String[]
    for len in lengths
        for i in 1:(length(sequence) - len + 1)
            start_pos = i
            end_pos = i + len - 1
            # Include peptide if amino acid locus is within the peptide boundaries
            if start_pos ≤ aa_locus ≤ end_pos
                push!(peptides, sequence[start_pos:end_pos])
            end
        end
    end
    return peptides
end

# Function to print help message for command-line usage
function print_help()
    println("""
    Usage:
        julia generate_context_peptides.jl --folder <path_to_folder> [--n_loci <number_of_loci>]

    Arguments:
        --folder   Path to the folder containing 'frames.csv'.
        --n_loci   Number of loci to simulate (default: 1000).
        --seed     Random number seed (default: 1320).
    """)
end

s = ArgParseSettings()
@add_arg_table! s begin
    "--folder"
        help = "Path to folder containing frames.csv"
        arg_type = String
    "--n_loci"
        help = "Number of loci to simulate"
        arg_type = Int
        default = 1000
    "--seed"
        help = "Random number seed"
        arg_type = Int
        default = 1320
end

parsed_args = parse_args(s)
folder = parsed_args["folder"]
n_loci = parsed_args["n_loci"]
seed = parsed_args["seed"]
frames_path = joinpath(folder, "frames.csv")

# Load the CSV containing genomic frames data
frames = CSV.read(frames_path, DataFrame; normalizenames=true)

# Define output directory based on frames_path
outdir = dirname(frames_path)
outdir = normpath(outdir)

# Extract Start and End positions from the Region column
frames.Start = map(r -> parse(Int, split(split(r, ";")[1], ",")[1]), frames.Region)
frames.End   = map(r -> parse(Int, split(split(r, ";")[end], ",")[end]), frames.Region)

# Set random seed for reproducibility
Random.seed!(seed)

# Generate random loci uniformly across all unique loci in all frames

# Build a list of all possible loci in all frames
loci_set = Set{Tuple{Int,Int}}() # (frame_idx, locus)
for i in 1:nrow(frames)
    for pos in frames.Start[i]:frames.End[i]
        push!(loci_set, (i, pos))
    end
end
all_loci = collect(loci_set)

# Define possible nucleotides for mutation (as Chars to match collected sequence elements)
nucleotides = ['A', 'C', 'G', 'T']

function main()
    # Collect n_loci non-synonymous random mutations
    context_rows = DataFrame()
    tries = 0
    while nrow(context_rows) < n_loci
        tries += 1
        # Randomly select a locus
        (frame_idx, locus) = rand(all_loci)
        frame = frames[frame_idx, :]
        relative_locus = locus - frame.Start
        consensus_seq = frame.Consensus_sequence
        if relative_locus < 1 || relative_locus > length(consensus_seq)
            continue
        end
    consensus_base = consensus_seq[relative_locus]
        # Determine candidate alt bases that produce a non-synonymous change
        candidate_alts = setdiff(nucleotides, [consensus_base])
    non_syn_alts = Char[]
        cons_aa = translate_dna_to_protein(consensus_seq)
        aa_locus_tmp = floor(Int, relative_locus / 3) + 1
        for alt in candidate_alts
            seq_tmp = collect(consensus_seq)
            seq_tmp[relative_locus] = alt
            vseq_tmp = join(seq_tmp)
            vprot_tmp = translate_dna_to_protein(vseq_tmp)
            if aa_locus_tmp >= 1 && aa_locus_tmp <= min(length(cons_aa), length(vprot_tmp))
                if cons_aa[aa_locus_tmp] != vprot_tmp[aa_locus_tmp]
                    push!(non_syn_alts, alt)
                end
            end
        end
        # If no non-synonymous alt exists, skip this locus
        if isempty(non_syn_alts)
            continue
        end
        # RNG chooses one alternative base for this locus
    variant_base = rand(non_syn_alts)
        # Edit consensus sequence to introduce chosen variant
        seq = collect(consensus_seq)
        seq[relative_locus] = variant_base
        variant_seq = join(seq)
        # Translate both sequences
        consensus_aa = translate_dna_to_protein(consensus_seq)
        variant_aa = translate_dna_to_protein(variant_seq)
        aa_locus = floor(Int, relative_locus / 3) + 1
        # Only keep if non-synonymous
        if aa_locus < 1 || aa_locus > min(length(consensus_aa), length(variant_aa))
            continue
        end
        if consensus_aa[aa_locus] == variant_aa[aa_locus]
            continue
        end
        # Build row
        new_row = DataFrame(
            Locus = locus,
            Relative_Locus = relative_locus,
            AA_Locus = aa_locus,
            Consensus_sequence = consensus_seq,
            Variant_sequence = variant_seq,
            Consensus_AA_sequence = consensus_aa,
            Variant_AA_sequence = variant_aa,
            Description = frame.Description
        )
        context_rows = vcat(context_rows, new_row)
    end
    # ...existing code...

    # ...existing code...

    # The rest of the script that uses context_rows should be moved inside this function.
    # For example, replace all uses of 'joined' with 'context_rows' below.

    joined = context_rows

    # Translate consensus and variant DNA sequences to amino acid sequences
    joined.Consensus_AA_sequence = [translate_dna_to_protein(seq) for seq in joined.Consensus_sequence]
    joined.Variant_AA_sequence = [translate_dna_to_protein(seq) for seq in joined.Variant_sequence]

    # Function to generate peptides containing the amino acid locus for given lengths
    substr_lengths = [8, 9, 10, 11]
    joined.AA_Locus = floor.(Int, joined.Relative_Locus ./ 3) .+ 1

    # Prepare DataFrame to hold generated peptides and labels
    peptides_df = DataFrame(
        Locus = Int[],
        Relative_Locus = Int[],
        AA_Locus = Int[],
        Consensus_Peptide = String[],
        Variant_Peptide = String[],
        Peptide_label = String[]
    )

    # saturation log path
    saturation_log = joinpath(outdir, "saturation_log.csv")
    skipped_saturated = 0

    # Generate peptides and labels for each row in joined data
    for row in eachrow(joined)
        # For context-run we may allow multiple alternative bases per locus as long as each
        # resulting amino-acid substitution is distinct. Build candidate alt bases by
        # testing all possible nucleotide changes and keep only those that are non-synonymous
        # and produce distinct amino-acid substitutions.
    ref_base = string(row.Consensus_sequence[row.Relative_Locus])
    alt_base = string(row.Variant_sequence[row.Relative_Locus])
        site = Site(row.Relative_Locus, ref_base, [alt_base], nothing)

        combos = expand_sites([site]; max_alts_per_site=1, max_combinations=10000, sample_on_saturate=true, saturation_log=saturation_log, stop_on_saturation=true)
        if isempty(combos)
            # saturated and policy says stop -> skip this window
            skipped_saturated += 1
            continue
        end

        # Map combos (each is a Vector{String} of length 1) to peptides. Each combo choice is either ref_base or alt_base.
        choices = Set([c[1] for c in combos])
        cps = generate_peptides(row.Consensus_AA_sequence, row.AA_Locus, substr_lengths)
        vps = generate_peptides(row.Variant_AA_sequence, row.AA_Locus, substr_lengths)

        # Only proceed if consensus and variant peptide counts match
        if length(cps) == length(vps)
            consensus_aa = (1 ≤ row.AA_Locus ≤ length(row.Consensus_AA_sequence)) ? row.Consensus_AA_sequence[row.AA_Locus] : "?"
            variant_aa = (1 ≤ row.AA_Locus ≤ length(row.Variant_AA_sequence)) ? row.Variant_AA_sequence[row.AA_Locus] : "?"

            change_label = "$(consensus_aa)$(row.AA_Locus)$(variant_aa)"
            base = "$(change_label)_$(replace(String(row.Description), " " => "_"))"

            # If both choices are present, include both peptides (dedupe later)
            for (i, (c, v)) in enumerate(zip(cps, vps))
                if ref_base in choices
                    label = "$(base)_$(i)"
                    push!(peptides_df, (row.Locus, row.Relative_Locus, row.AA_Locus, c, v, label))
                end
                if alt_base in choices
                    label = "$(base)_$(i)"
                    push!(peptides_df, (row.Locus, row.Relative_Locus, row.AA_Locus, c, v, label))
                end
            end
        end
    end

    # Filter out peptides that contain stop codons or where consensus and variant peptides are identical
    initial_rows = nrow(peptides_df)
    dropped_df = filter(row -> row.Consensus_Peptide == row.Variant_Peptide ||
                                occursin('*', row.Consensus_Peptide) ||
                                occursin('*', row.Variant_Peptide), peptides_df)
    removed_loci = length(unique(dropped_df.Locus))

    filtered_df = filter(row -> row.Consensus_Peptide != row.Variant_Peptide &&
                                !occursin('*', row.Consensus_Peptide) &&
                                !occursin('*', row.Variant_Peptide), peptides_df)

    final_rows = nrow(filtered_df)
    removed_rows = initial_rows - final_rows
    # ...existing code...

    # Flatten DataFrame to separate consensus and variant peptides for output
    flattened_df = DataFrame(
        Locus = Int[],
        Peptide = String[],
        Peptide_label = String[]
    )

    for row in eachrow(filtered_df)
        # Add consensus peptide with label suffix '_C'
        push!(flattened_df, (
            row.Locus,
            row.Consensus_Peptide,
            "$(row.Peptide_label)_C"
        ))

        # Add variant peptide with label suffix '_V'
        push!(flattened_df, (
            row.Locus,
            row.Variant_Peptide,
            "$(row.Peptide_label)_V"
        ))
    end

    # Deduplicate peptides by Locus, Peptide, and State (C or V)
    flattened_df[!, :State] = last.(split.(flattened_df.Peptide_label, "_"))
    dedup_df = unique(flattened_df, [:Locus, :Peptide, :State])
    select!(dedup_df, Not(:State))

    # Save deduplicated peptides with labels to CSV file
    CSV.write(joinpath(outdir, "context_peptides_labels.csv"), dedup_df)
    # ...existing code...

    # Save peptides only (no headers) to .pep file
    open(joinpath(outdir, "context_peptides.pep"), "w") do io
        for row in eachrow(dedup_df)
            println(io, row.Peptide)
        end
    end
    # ...existing code...
end

# Call main only after all definitions
main()