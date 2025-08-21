#!/usr/bin/env julia
using ArgParse
using CSV, DataFrames
using Random

using FilePathsBase
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

# Load the CSV
frames = CSV.read(frames_path, DataFrame; normalizenames=true)

# Define output directory based on frames_path
outdir = joinpath(dirname(frames_path), "..")
outdir = normpath(outdir)

# Extract Start and End from Region
frames.Start = map(r -> parse(Int, split(split(r, ";")[1], ",")[1]), frames.Region)
frames.End   = map(r -> parse(Int, split(split(r, ";")[end], ",")[end]), frames.Region)

seed_path = joinpath(@__DIR__, "random_seed.txt")
seed = parse(Int, strip(read(seed_path, String)))
Random.seed!(seed)

# Generate random loci uniformly across all unique loci in all frames
loci_set = Set{Int}()
for i in 1:nrow(frames)
    for pos in frames.Start[i]:frames.End[i]
        push!(loci_set, pos)
    end
end
all_loci = collect(loci_set)
sampled_loci = rand(all_loci, n_loci)

# Create loci DataFrame
random_loci_df = DataFrame(Locus = sampled_loci)

# Cross join with frames and filter to those where the locus falls inside frame
joined = crossjoin(random_loci_df, frames)
filter!(row -> row.Start ≤ row.Locus ≤ row.End, joined)

# Compute Relative_Locus
joined.Relative_Locus = joined.Locus .- joined.Start

# Pull base from Consensus_sequence
joined.Pulled_Base = [
    (1 ≤ rl ≤ length(seq) ? string(seq[rl]) : missing)
    for (rl, seq) in zip(joined.Relative_Locus, joined.Consensus_sequence)
]

# Define possible nucleotides
nucleotides = ["A", "C", "G", "T"]

# Generate random mutations that are different from the consensus base
joined.Mutated_Base = [
    ismissing(base) ? missing :
    rand(setdiff(nucleotides, [base]))
    for base in joined.Pulled_Base
]

# Rename for compatibility
joined.Variant = joined.Mutated_Base

# Edit consensus sequence with variant
joined.Variant_sequence = similar(joined.Consensus_sequence)
for row in eachrow(joined)
    rl = row.Relative_Locus
    consensus = row.Consensus_sequence
    variant_nt = row.Variant

    if !ismissing(rl) && 1 ≤ rl ≤ length(consensus)
        seq = collect(consensus)
        seq[rl] = variant_nt[1]
        row.Variant_sequence = join(seq)
    else
        row.Variant_sequence = consensus
    end
end

# Translation dictionary
if !isdefined(Main, :CODON_DICT)
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
end

# Translation function
function translate_dna_to_protein(dna_sequence::String)::String
    join([get(CODON_DICT, uppercase(dna_sequence[i:i+2]), "?") 
          for i in 1:3:length(dna_sequence)-2], "")
end

# Translate consensus and variant sequences
joined.Consensus_AA_sequence = [translate_dna_to_protein(seq) for seq in joined.Consensus_sequence]
joined.Variant_AA_sequence = [translate_dna_to_protein(seq) for seq in joined.Variant_sequence]


# Generate peptides centered on AA_Locus
function generate_peptides(sequence::String, aa_locus::Int, lengths::Vector{Int})::Vector{String}
    peptides = String[]
    for len in lengths
        for i in 1:(length(sequence) - len + 1)
            start_pos = i
            end_pos = i + len - 1
            if start_pos ≤ aa_locus ≤ end_pos
                push!(peptides, sequence[start_pos:end_pos])
            end
        end
    end
    return peptides
end

# Peptide generation parameters
substr_lengths = [8, 9, 10, 11]
joined.AA_Locus = floor.(Int, joined.Relative_Locus ./ 3) .+ 1

# Prepare output DataFrame
peptides_df = DataFrame(
    Locus = Int[],
    Relative_Locus = Int[],
    AA_Locus = Int[],
    Consensus_Peptide = String[],
    Variant_Peptide = String[],
    Peptide_label = String[]
)

# Generate peptides and labels
for row in eachrow(joined)
    cps = generate_peptides(row.Consensus_AA_sequence, row.AA_Locus, substr_lengths)
    vps = generate_peptides(row.Variant_AA_sequence, row.AA_Locus, substr_lengths)

    if length(cps) == length(vps)
        consensus_aa = (1 ≤ row.AA_Locus ≤ length(row.Consensus_AA_sequence)) ? row.Consensus_AA_sequence[row.AA_Locus] : "?"
        variant_aa = (1 ≤ row.AA_Locus ≤ length(row.Variant_AA_sequence)) ? row.Variant_AA_sequence[row.AA_Locus] : "?"

        change_label = "$(consensus_aa)$(row.AA_Locus)$(variant_aa)"
        base = "$(change_label)_$(replace(String(row.Description), " " => "_"))"

        for (i, (c, v)) in enumerate(zip(cps, vps))
            label = "$(base)_$(i)"
            push!(peptides_df, (row.Locus, row.Relative_Locus, row.AA_Locus, c, v, label))
        end
    end
end



# Filter out peptides with stop codons or identical sequences
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

# Flatten to separate consensus and variant peptides
flattened_df = DataFrame(
    Locus = Int[],
    Peptide = String[],
    Peptide_label = String[]
)

for row in eachrow(filtered_df)
    push!(flattened_df, (
        row.Locus,
        row.Consensus_Peptide,
        "$(row.Peptide_label)_C"
    ))

    push!(flattened_df, (
        row.Locus,
        row.Variant_Peptide,
        "$(row.Peptide_label)_V"
    ))
end

# Deduplicate peptides
flattened_df[!, :State] = last.(split.(flattened_df.Peptide_label, "_"))
dedup_df = unique(flattened_df, [:Locus, :Peptide, :State])
select!(dedup_df, Not(:State))

# Save to CSV
CSV.write(joinpath(outdir, "context_peptides_labels.csv"), dedup_df)
println("context_peptides_labels.csv file has been written.")

# Save to .pep (no headers)
open(joinpath(outdir, "context_Peptides.pep"), "w") do io
    for row in eachrow(dedup_df)
        println(io, row.Peptide)
    end
end
println("context_Peptides.pep file has been written.")

exit()