using CSV, DataFrames

function process_csv(file_path::String)
    # Read the first row as a list (raw content of the file)
    first_row = CSV.File(file_path, header=false) |> first |> collect

    # Load the CSV file starting from the second row
    df = CSV.File(file_path, header=2) |> DataFrame

    # Column renaming logic
    column_counts = Dict{String, Int}()
    new_column_names = Vector{String}(undef, length(names(df)))

    for (i, col_name) in enumerate(names(df))
        if haskey(column_counts, col_name)
            # Increment the count for duplicates and assign a suffix
            column_counts[col_name] += 1
            new_column_names[i] = "$(col_name)_$(column_counts[col_name])"
        else
            # First occurrence: Check if it will have duplicates
            if count(x -> x == col_name, names(df)) > 1
                column_counts[col_name] = 1
                new_column_names[i] = "$(col_name)_1"
            else
                column_counts[col_name] = 1
                new_column_names[i] = col_name
            end
        end
    end

    # Apply the new column names to the DataFrame
    rename!(df, Dict(zip(names(df), new_column_names)))

    return first_row, df
end

# Test the function
file_path = "/Users/e.smith.5/Documents/PhD/CD8scape/data/RSV_example/RSV.csv"
first_row, df = process_csv(file_path)

println("First Row:")
println(first_row)

println("Processed Headers:")
println(names(df))

println("Processed DataFrame:")
display(df)