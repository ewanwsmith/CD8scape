
#!/usr/bin/env julia
include("env.jl")
using CSV, DataFrames

function parse_arguments()
    args = Dict()
    for (i, arg) in enumerate(ARGS)
        if arg in ["--folder", "-f"]
            args["folder"] = ARGS[i + 1]
        end
    end
    if !haskey(args, "folder")
        error("Usage: ./run_netMHCpan.jl --folder /path/to/data")
    end
    return args
end

function find_settings_file()
    current_dir = pwd()
    while current_dir != "/"
        candidate = joinpath(current_dir, "src", "settings.txt")
        if isfile(candidate)
            return candidate
        end
        current_dir = dirname(current_dir)
    end
    error("settings.txt not found in any 'src' directory above current working directory.")
end

function get_netMHCpan_path()
    settings_file = find_settings_file()

    if !isfile(settings_file)
        error("settings.txt not found at $settings_file. Please provide the file.")
    end

    for line in readlines(settings_file)
        line = strip(line)
        if !isempty(line) && !startswith(line, "#")
            path = replace(line, r"['\"]" => "") |> normpath
            println("NetMHCpan path -> ", path)

            if isdir(path)
                return path
            else
                error("ERROR: Path '$path' does not exist.")
            end
        end
    end

    error("ERROR: No valid netMHCpan path found in settings.txt.")
end

function main()
    args = parse_arguments()
    folder_path = args["folder"]

    netMHCpan_path = get_netMHCpan_path()
    println("Using NetMHCpan from: ", netMHCpan_path)

    # File paths
    function find_representatives_file()
        current_dir = pwd()
        while current_dir != "/"
            candidate = joinpath(current_dir, "src", "supertype_panel.csv")
            if isfile(candidate)
                return candidate
            end
            current_dir = dirname(current_dir)
        end
        error("supertype_panel.csv not found in any 'src' directory above current working directory.")
    end
    representatives_file = find_representatives_file()
    peptides_file = joinpath(folder_path, "Peptides.pep")
    xlsfile_path = joinpath(folder_path, "netMHCpan_output.tsv")

    # Check if required files exist
    if !isfile(peptides_file)
        error("File $peptides_file not found. Please ensure it exists in the specified folder.")
    end

    # Test write permissions
    test_file = joinpath(folder_path, "test_write.txt")
    try
        open(test_file, "w") do io
            write(io, "test")
        end
        rm(test_file)
    catch e
        error("ERROR: Cannot write to output folder. Check permissions!")
    end

    # Read alleles from the supertype_panel.csv file
    df = CSV.read(representatives_file, DataFrame)
    df = filter(row -> row.Frequency != 0, df)
    if !(:Allele in names(df))
        error("ERROR: Column 'Allele' not found in $(representatives_file).")
    end
    allele_list = [replace(String(a), "*" => "") for a in collect(skipmissing(df.Allele)) if !isempty(String(a))]
    alleles = join(allele_list, ",")
    cmd = Cmd([joinpath(netMHCpan_path, "netMHCpan"), "-p", peptides_file, "-xls", "-a", alleles, "-xlsfile", xlsfile_path])
    println("Running NetMHCpan with command:")
    println(cmd)
    
    try
        run(cmd)

        if isfile(xlsfile_path)
            println("NetMHCpan output found at $xlsfile_path")
        else
            error("ERROR: NetMHCpan did not create the expected output file.")
        end
        println("NetMHCpan completed successfully.")
    catch e
        println("Error running NetMHCpan: ", e)
        exit(1)
    end
end

main()
