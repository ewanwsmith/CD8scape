using Pkg

# Activate the environment in the current directory
Pkg.activate(".")

# List of packages to add
packages = [
    "DataFrames",
    "CSV",
    # Add other necessary packages
]

# Add each package
for pkg in packages
    Pkg.add(pkg)
end