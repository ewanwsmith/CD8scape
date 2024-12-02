using DataFrames

function read_locus_file(filepath)
    # Initialize arrays to store the data
    loci = Int[]
    original_bases = String[]
    variant_bases = String[]

    # Open and read the file line by line
    open(filepath, "r") do f
        for line in eachline(f)
            tokens = split(line)
            if length(tokens) >= 3
                # Parse the first three columns
                locus = parse(Int, tokens[1])
                original_base = tokens[2]
                variant_base = tokens[3]

                # Append the parsed data to the arrays
                push!(loci, locus)
                push!(original_bases, original_base)
                push!(variant_bases, variant_base)
            else
                @warn "Line with insufficient data: $line"
            end
        end
    end

    # Create a DataFrame from the collected data
    df = DataFrame(
        locus = loci,
        original_base = original_bases,
        variant_base = variant_bases
    )
    return df
end

# Example usage:
loci = read_locus_file("/Users/e.smith.5/Documents/PhD/CD8scape/data/RSV_example/F/single_locus_trajectories10_F_Fusion_protein.out")

using DataFrames

function read_haps_file(filepath)
    # Initialize an array to store the haplotypes
    haps = String[]
    
    # Define a regular expression pattern to match lines starting with haplotype sequences
    pattern = r"^[ACGT]+"
    
    # Open and read the file line by line
    open(filepath, "r") do f
        for line in eachline(f)
            line = strip(line)
            if isempty(line)
                continue  # Skip empty lines
            end
            if occursin(pattern, line)
                # Split the line into tokens
                tokens = split(line)
                # The first token is the haplotype
                hap = tokens[1]
                # Verify that the rest of the tokens are numbers
                if all(x -> tryparse(Float64, x) !== nothing, tokens[2:end])
                    push!(haps, hap)
                else
                    @warn "Skipping line with non-numeric data: $line"
                end
            else
                # Line doesn't start with a haplotype sequence; skip it
                @warn "Skipping line: $line"
            end
        end
    end
    
    # Create a DataFrame with the collected haplotypes
    df = DataFrame(haps = haps)
    return df
end

# Example usage:
haps = read_haps_file("/Users/e.smith.5/Documents/PhD/CD8scape/data/RSV_example/F/Inference_14_0.out")

using DataFrames

function construct_bases_dataframe(loci_df::DataFrame, haps_df::DataFrame)
    # Step 1: Extract the unique loci and sort them
    loci = unique(loci_df.locus)
    sort!(loci)
    num_loci = length(loci)

    # Step 2: Initialize a DataFrame with columns for each locus
    # We'll also include a column for haplotype identifiers
    bases = DataFrame()

    # Add the haplotype strings as an identifier column
    bases[!, :haps] = haps_df.haps

    # Step 3: Check that each haplotype string has the same length as the number of loci
    for (i, hap_string) in enumerate(haps_df.haps)
        hap_length = length(hap_string)
        if hap_length != num_loci
            error("Length of haplotype string (length $hap_length) does not match number of loci ($num_loci) for haplotype at row $(i).")
        end
    end

    # Step 4: Add columns for each locus, filling in the bases from each haplotype
    # The position in the haplotype string corresponds to the sorted loci positions
    for (idx, locus) in enumerate(loci)
        # Column name based on the locus
        column_name = Symbol("locus_$(locus)")
        # Extract the base at this position for each haplotype
        bases[!, column_name] = [hap[idx] for hap in haps_df.haps]
    end

    return bases
end

# Example usage:
bases = construct_bases_dataframe(loci, haps)