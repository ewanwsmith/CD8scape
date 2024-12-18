#!/usr/bin/env julia

using CSV
using DataFrames
using RCall

"""
    read_net_scores(folder::String) -> DataFrame

Reads the `net_scores.csv` file from the specified folder into a DataFrame.

# Arguments
- `folder::String`: Path to the folder containing `net_scores.csv`.

# Returns
- A `DataFrame` containing the data from `net_scores.csv`.

# Errors
Throws an error if the file is not found.
"""
function read_net_scores(folder::String)::DataFrame
    file_path = joinpath(folder, "net_scores.csv")
    if isfile(file_path)
        return CSV.read(file_path, DataFrame)
    else
        error("File `net_scores.csv` not found in folder: $folder")
    end
end

# Get folder path from command line argument
if length(ARGS) == 0
    println("Usage: julia script_name.jl <folder_path>")
    exit(1)
end

folder_path = ARGS[1]
println("Extracting plot data from: ", folder_path)

# Read scores data
scores = read_net_scores(folder_path)

# Suppress RCall messages temporarily
redirect_stdout(devnull) do
    redirect_stderr(devnull) do
        # Pass the scores DataFrame and folder path to the R environment
        @rput scores
        @rput folder_path

        R"""
        # Load required R libraries
        suppressMessages({
            library(ggplot2)
            library(viridis)
            library(forcats)
            library(ggpubr)
        })

        options(warn = -1)  # Suppress warnings

        # Scatter Plot: EL Rank Comparison
        rank_scatterplot = ggplot(scores, aes(x = EL_rank_V, y = EL_rank_C, color = as.factor(Locus))) +
            geom_point() +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
            labs(
                title = "EL Ranks by MHC", 
                x = "Variant EL Rank", 
                y = "Consensus EL Rank", 
                color = "Locus"
            ) +
            scale_color_viridis(discrete = TRUE) +
            theme_minimal() +
            facet_wrap(~ MHC)

        # Boxplot: Net EL Rank with Point Size = 1 and Horizontal Line
        rank_boxplot = ggplot(scores, aes(x = fct_rev(as.factor(Locus)), y = EL_rank_net, color = as.factor(Locus))) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.2, size = 1) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
            ylim(-50, 50) +
            labs(
                title = "Net EL Rank by MHC", 
                x = "Locus", 
                y = "Net EL Rank", 
                color = "Locus"
            ) +
            scale_color_viridis(discrete = TRUE) +
            theme_minimal() +
            coord_flip() +
            facet_wrap(~ MHC)

        # Arrange Scatter Plot and Boxplot Vertically
        rank_plots = ggarrange(
            rank_scatterplot, 
            rank_boxplot, 
            ncol = 1, nrow = 2, common.legend = TRUE, legend = "right"
        )

        # Save the combined plot as a JPEG file
        output_path <- file.path(folder_path, "net_scores.jpeg")
        ggsave(output_path, plot = rank_plots, width = 10, height = 12, units = "in", dpi = 300)

        options(warn = 0)  # Restore warning settings
        """
    end
end