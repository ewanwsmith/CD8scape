using Pkg

# Activate the environment in the current directory without output
Pkg.activate(".", io=devnull)

# List of packages to add
packages = [
    "DataFrames",
    "CSV",
    "FilePathsBase",
]

# Add each package without output
for pkg in packages
    Pkg.add(pkg; io=devnull)
end