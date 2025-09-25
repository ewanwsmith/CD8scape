
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

# Function to generate peptides containing the amino acid locus for given lengths
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
    julia generate_context_peptides.jl --folder <path_to_folder> [--n_loci <number_of_loci>]

Arguments:
    --folder   Path to the folder containing 'frames.csv'.
    --n_loci   Number of loci to simulate (default: 1000).
"""

using ArgParse
using CSV, DataFrames
using Random

using FilePathsBase

# Function to print help message for command-line usage
function print_help()
    println("""
    Usage:
        julia generate_context_peptides.jl --folder <path_to_folder> [--n_loci <number_of_loci>]

    Arguments:
        --folder   Path to the folder containing 'frames.csv'.
        --n_loci   Number of loci to simulate (default: 1000).
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
end

parsed_args = parse_args(s)
folder = parsed_args["folder"]
n_loci = parsed_args["n_loci"]
frames_path = joinpath(folder, "frames.csv")

# Load the CSV containing genomic frames data
frames = CSV.read(frames_path, DataFrame; normalizenames=true)

# Define output directory based on frames_path
outdir = dirname(frames_path)
outdir = normpath(outdir)

# Extract Start and End positions from the Region column
frames.Start = map(r -> parse(Int, split(split(r, ";")[1], ",")[1]), frames.Region)
frames.End   = map(r -> parse(Int, split(split(r, ";")[end], ",")[end]), frames.Region)

# Read and set random seed for reproducibility
seed_path = joinpath(@__DIR__, "random_seed.txt")
seed = parse(Int, strip(read(seed_path, String)))
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

# Define possible nucleotides for mutation
nucleotides = ["A", "C", "G", "T"]

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
        consensus_base = string(consensus_seq[relative_locus])
        # Randomly select a mutated base different from consensus
    variant_base = rand(setdiff(nucleotides, [consensus_base]))[1] # ensure Char
    # Edit consensus sequence to introduce variant
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
    println("Generated $n_loci non-synonymous random mutations in $tries attempts.")

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

    # Generate peptides and labels for each row in joined data
    for row in eachrow(joined)
        cps = generate_peptides(row.Consensus_AA_sequence, row.AA_Locus, substr_lengths)
        vps = generate_peptides(row.Variant_AA_sequence, row.AA_Locus, substr_lengths)

        # Only proceed if consensus and variant peptide counts match
        if length(cps) == length(vps)
            # Extract amino acids at the mutation locus for labeling
            consensus_aa = (1 ≤ row.AA_Locus ≤ length(row.Consensus_AA_sequence)) ? row.Consensus_AA_sequence[row.AA_Locus] : "?"
            variant_aa = (1 ≤ row.AA_Locus ≤ length(row.Variant_AA_sequence)) ? row.Variant_AA_sequence[row.AA_Locus] : "?"

            change_label = "$(consensus_aa)$(row.AA_Locus)$(variant_aa)"
            base = "$(change_label)_$(replace(String(row.Description), " " => "_"))"

            # Label and store each peptide pair
            for (i, (c, v)) in enumerate(zip(cps, vps))
                label = "$(base)_$(i)"
                push!(peptides_df, (row.Locus, row.Relative_Locus, row.AA_Locus, c, v, label))
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
    println("Removed $removed_rows peptides from $removed_loci loci due to stop codons.")

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
    println("context_peptides_labels.csv file has been written.")

    # Save peptides only (no headers) to .pep file
    open(joinpath(outdir, "context_peptides.pep"), "w") do io
        for row in eachrow(dedup_df)
            println(io, row.Peptide)
        end
    end
    println("context_peptides.pep file has been written.")
end