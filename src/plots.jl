#!/usr/bin/env julia

using DataFrames, CSV, RCall

# Load dataset
function load_data(folder::String)
    path = joinpath(folder, "harmonic_mean_best_ranks.csv")
    df = CSV.read(path, DataFrame)
    # Remove rows with missing data in HMBR_C or HMBR_V
    dropmissing!(df, [:HMBR_C, :HMBR_V])
    return df
end

# Save plot
function save_plot(df::DataFrame, output_path::String)
    @rput df output_path
    R"""
    library(ggplot2)

    p <- ggplot(df, aes(x = HMBR_C, y = HMBR_V, color = as.factor(Locus))) +
        geom_point(size = 3) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
        labs(x = "HMBR_C", y = "HMBR_V", color = "Locus") +
        theme(legend.position = "bottom")

    ggsave(output_path, plot = p, device = "jpeg")
    """
end

# Main function
function main()
    input_path = ARGS[1]  # Get folder path from command line arguments
    df = load_data(input_path)
    output_path = joinpath(input_path, "output_plot.jpeg")
    save_plot(df, output_path)
end

main()