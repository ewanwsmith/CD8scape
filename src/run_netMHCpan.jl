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

function get_netMHCpan_executable()
    settings_file = find_settings_file()

    if !isfile(settings_file)
        error("settings.txt not found at $settings_file. Please provide the file.")
    end

    for line in readlines(settings_file)
        line = strip(line)
        if !isempty(line) && !startswith(line, "#")
            raw = replace(line, r"['\"]" => "") |> normpath
            println("NetMHCpan path -> ", raw)

            # Case 1: settings.txt contains a directory path; expect the executable inside it
            if isdir(raw)
                exe = joinpath(raw, "netMHCpan")
                if isfile(exe)
                    return exe
                else
                    error("ERROR: Directory '$raw' does not contain 'netMHCpan' executable.")
                end
            end

            # Case 2: settings.txt contains a direct path to the executable
            if isfile(raw)
                return raw
            end

            error("ERROR: Path '$raw' does not exist.")
        end
    end

    error("ERROR: No valid netMHCpan path found in settings.txt.")
end

function main()
    args = parse_arguments()
    folder_path = args["folder"]

    netMHCpan_exe = get_netMHCpan_executable()
    println("Using NetMHCpan executable: ", netMHCpan_exe)

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
        rm(test_file)
    catch e
        error("ERROR: Cannot write to output folder. Check permissions!")
    end

    # Read alleles
    allele_list = open(alleles_file) do file
        [replace(split(line, r"\s+")[1], "*" => "") for line in readlines(file) if !isempty(line)]
    end
    alleles = join(allele_list, ",")

    cmd = `$netMHCpan_exe -p $peptides_file -xls -a $alleles -xlsfile $xlsfile_path`
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
