using DataFrames, CSV

# Load the CSV file into a DataFrame
input_file = "input_file.csv"
df = CSV.read(input_file, DataFrame)

# Columns for which net values are needed
columns_of_interest = [
    "netmhcpan_ba percentile",
    "immunogenicity score",
    "proteasome score",
    "tap score",
    "mhc score",
    "processing score",
    "processing total score"
]

# Ensure column names are symbols for easier processing
rename!(df, Symbol)

# Extract all unique `seq #` values
seq_numbers = unique(df.:seq)
seq_1 = df[df.:seq .== 1, :]  # Rows where `seq #` is 1

# Prepare a DataFrame to store results
results = DataFrame()

for seq_val in seq_numbers
    if seq_val == 1  # Skip the `seq #` of 1 for pairwise comparisons
        continue
    end

    # Rows with current `seq #`
    seq_other = df[df.:seq .== seq_val, :]

    # Pair rows based on peptide (assumes peptide column can match pairs)
    pairs = innerjoin(seq_1, seq_other, on=:peptide, suffix=("_seq1", "_seq_other"))

    # Calculate net differences for the columns of interest
    for col in columns_of_interest
        pairs[Symbol("net_", col)] = pairs[Symbol(col, "_seq_other")] - pairs[Symbol(col, "_seq1")]
    end

    # Append to results
    append!(results, pairs)
end

# Save the results to a new CSV file in the same directory as the input
output_file = "net_scores.csv"
CSV.write(output_file, results)