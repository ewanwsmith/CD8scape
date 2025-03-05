using Pkg

# Activate the environment in the current directory silently
Pkg.activate(".", io = devnull)

# Install RCall first, silently
Pkg.add("RCall"; io = devnull)

# Load RCall after installation
using RCall

# List of other Julia packages to install silently
packages = [
    "DataFrames",
    "CSV",
    "FilePathsBase",
    "ArgParse",
    "CodecZlib",
    "Statistics",
    "StatsBase",
]

# Add remaining Julia packages without output
for pkg in packages
    Pkg.add(pkg; io = devnull)
end

println("Julia dependencies installed.")

# Silently install R packages using RCall
R"""
options(warn = -1)  # Suppress warnings
suppressMessages({
    if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "https://cran.r-project.org")
    if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis", repos = "https://cran.r-project.org")
    if (!requireNamespace("forcats", quietly = TRUE)) install.packages("forcats", repos = "https://cran.r-project.org")
    if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr", repos = "https://cran.r-project.org")
})
options(warn = 0)  # Restore warnings
"""

println("R dependencies installed.")