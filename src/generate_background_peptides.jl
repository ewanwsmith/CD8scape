#!/usr/bin/env julia

using CSV
using DataFrames

function print_help()
    println("Usage: julia generate_background_peptides.jl <data_folder>")
    println("")
    println("This script processes frames.csv in the specified data_folder and generates:")
    println("  1. background_peptides_labels.csv: CSV with peptide and label columns")
    println("  2. Background_Peptides.pep: plain text file of peptides (no headers)")
end

function translate_dna_to_protein(seq::String)::String
    codon_table = Dict(
        "ATA"=>"I", "ATC"=>"I", "ATT"=>"I", "ATG"=>"M",
        "ACA"=>"T", "ACC"=>"T", "ACG"=>"T", "ACT"=>"T",
        "AAC"=>"N", "AAT"=>"N", "AAA"=>"K", "AAG"=>"K",
        "AGC"=>"S", "AGT"=>"S", "AGA"=>"R", "AGG"=>"R",
        "CTA"=>"L", "CTC"=>"L", "CTG"=>"L", "CTT"=>"L",
        "CCA"=>"P", "CCC"=>"P", "CCG"=>"P", "CCT"=>"P",
        "CAC"=>"H", "CAT"=>"H", "CAA"=>"Q", "CAG"=>"Q",
        "CGA"=>"R", "CGC"=>"R", "CGG"=>"R", "CGT"=>"R",
        "GTA"=>"V", "GTC"=>"V", "GTG"=>"V", "GTT"=>"V",
        "GCA"=>"A", "GCC"=>"A", "GCG"=>"A", "GCT"=>"A",
        "GAC"=>"D", "GAT"=>"D", "GAA"=>"E", "GAG"=>"E",
        "GGA"=>"G", "GGC"=>"G", "GGG"=>"G", "GGT"=>"G",
        "TCA"=>"S", "TCC"=>"S", "TCG"=>"S", "TCT"=>"S",
        "TTC"=>"F", "TTT"=>"F", "TTA"=>"L", "TTG"=>"L",
        "TAC"=>"Y", "TAT"=>"Y", "TAA"=>"*", "TAG"=>"*",
        "TGC"=>"C", "TGT"=>"C", "TGA"=>"*", "TGG"=>"W"
    )

    aa_seq = IOBuffer()
    for i in 1:3:length(seq)-2
        codon = uppercase(seq[i:i+2])
        aa = get(codon_table, codon, "?")
        print(aa_seq, aa)
    end
    return String(take!(aa_seq))
end

function generate_background_peptides(aa_sequence::String, region::String, description::String)
    peptide_lengths = [8, 9, 10, 11]
    peptides = String[]
    labels = String[]

    region_parts = split(region, ";")
    start_nt, _ = parse.(Int, split(region_parts[1], ","))
    _, end_nt = parse.(Int, split(region_parts[end], ","))

    description_clean = replace(description, " " => "_")

    for len in peptide_lengths
        for i in 1:(length(aa_sequence) - len + 1)
            peptide = aa_sequence[i:(i+len-1)]
            peptide_start_nt = start_nt + (i - 1) * 3
            peptide_end_nt = start_nt + (i + len - 2) * 3
            push!(peptides, peptide)
            label = "$(peptide_start_nt)-$(peptide_end_nt)_$(description_clean)_A"
            push!(labels, label)
        end
    end

    return peptides, labels
end

function create_background_peptide_dataframe(frames::DataFrame)
    peps, labels = String[], String[]
    for row in eachrow(frames)
        aa_seq = if :Consensus_AA_sequence in propertynames(row)
            row.Consensus_AA_sequence
        elseif :Consensus_sequence in propertynames(row)
            translate_dna_to_protein(String(row.Consensus_sequence))
        else
            continue
        end

        if !ismissing(aa_seq) && aa_seq != "missing"
            p, l = generate_background_peptides(
                String(aa_seq),
                String(row.Region),
                String(row.Description)
            )
            append!(peps, p); append!(labels, l)
        end
    end
    return DataFrame(Peptide = peps, Peptide_label = labels)
end

function write_peptides_file_no_headers(df::DataFrame, folder_path::String)
    file_path = joinpath(folder_path, "Background_Peptides.pep")
    open(file_path, "w") do io
        for row in eachrow(df)
            println(io, row.Peptide)
        end
    end
    println("Background_Peptides.pep written to: $file_path")
end

if length(ARGS) == 1 && ARGS[1] == "--help"
    print_help()
    exit(0)
elseif length(ARGS) < 1
    println("Error: No data folder provided.")
    println("Run `julia generate_background_peptides.jl --help` for usage.")
    exit(1)
end

data_folder = ARGS[1]
frames_path = joinpath(data_folder, "frames.csv")
frames = CSV.read(frames_path, DataFrame)

# If Consensus_AA_sequence is missing, derive it from Consensus_sequence
if :Consensus_AA_sequence ∉ names(frames) && :Consensus_sequence ∈ names(frames)
    frames.Consensus_AA_sequence = [translate_dna_to_protein(seq) for seq in frames.Consensus_sequence]
end

bg_df = create_background_peptide_dataframe(frames)
csv_path = joinpath(data_folder, "background_peptides_labels.csv")
CSV.write(csv_path, bg_df)
println("background_peptides_labels.csv written to: $csv_path")

write_peptides_file_no_headers(bg_df, data_folder)