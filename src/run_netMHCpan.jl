#!/usr/bin/env julia

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

function get_netMHCpan_path()
    settings_file = "/Users/e.smith.5/Documents/PhD/CD8scape/src/settings.txt"

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
    alleles_file = joinpath(folder_path, "alleles.txt")
    peptides_file = joinpath(folder_path, "Peptides.pep")
    xlsfile_path = joinpath(folder_path, "netMHCpan_output.tsv")

    # Check if required files exist
    for file in [alleles_file, peptides_file]
        if !isfile(file)
            error("File $file not found. Please ensure it exists in the specified folder.")
        end
    end

    # Test write permissions
    test_file = joinpath(folder_path, "test_write.txt")
    try
        open(test_file, "w") do io
            write(io, "test")
        end
        println("DEBUG: Write test succeeded. Removing test file...")
        rm(test_file)
    catch e
        error("ERROR: Cannot write to output folder. Check permissions!")
    end

    # Read alleles
    allele_list = open(alleles_file) do file
        [split(line, r"\s+")[1] for line in readlines(file) if !isempty(line)]
    end
    alleles = join(allele_list, ",")

    cmd = `$netMHCpan_path/netMHCpan -p $peptides_file -xls -a $alleles -xlsfile $xlsfile_path`
    
    println("DEBUG: Running NetMHCpan with command:")
    println(cmd)
    
    try
        run(cmd)

        println("DEBUG: Searching for .tsv files in output folder...")
        run(`find $folder_path -name "*.tsv"`)

        println("DEBUG: Checking if output file was created: $xlsfile_path")
        if isfile(xlsfile_path)
            println("DEBUG: NetMHCpan output found at $xlsfile_path")
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