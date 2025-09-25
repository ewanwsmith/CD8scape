function normalize_allele(a)
    # Remove '*' if present, but keep ':'
    replace(a, "*" => "")
end

function to_colon_format(a)
    # Converts compact allele format (e.g., HLA-A3201) to colon format (e.g., HLA-A32:01)
    m = match(r"(HLA-[A-Z]+)(\d{2})(\d{2,3})", a)
    if m !== nothing
        return string(m.captures[1], m.captures[2], ":", m.captures[3])
    else
        return a
    end
end
#!/usr/bin/env julia
"""
run_netMHCpan.jl

Runs NetMHCpan for all peptides using provided alleles.

Usage:
    julia run_netMHCpan.jl --folder <folder_path>

Arguments:
    --folder   Path to the folder containing input files.
"""

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
        s = strip(line)
        if isempty(s) || startswith(s, "#")
            continue
        end
        if occursin('=', s)
            k, v = strip.(split(s, '=', limit=2))
            if uppercase(k) == "NETMHCPAN"
                p = normpath(replace(v, "~" => homedir()))
                println("NetMHCpan path (settings) -> ", p)
                if isfile(p)
                    # After locating NetMHCpan, copy MHC_pseudo.dat into CD8scape
                    netMHCpan_dir = dirname(p)
                    pseudo_src = joinpath(netMHCpan_dir, "data", "MHC_pseudo.dat")
                    pseudo_dst = joinpath(dirname(settings_file), "MHC_pseudo.dat")
                    if isfile(pseudo_src)
                        cp(pseudo_src, pseudo_dst; force=true)
                        println("Copied MHC_pseudo.dat from $pseudo_src to $pseudo_dst")
                    else
                        println("Warning: Could not find MHC_pseudo.dat at $pseudo_src")
                    end
                    return p
                else
                    error("ERROR: Path '$p' does not exist.")
                end
            end
        else
            p = normpath(replace(s, "~" => homedir()))
            println("NetMHCpan path (settings) -> ", p)
            if isfile(p)
                # After locating NetMHCpan, copy MHC_pseudo.dat into CD8scape
                netMHCpan_dir = dirname(p)
                pseudo_src = joinpath(netMHCpan_dir, "data", "MHC_pseudo.dat")
                pseudo_dst = joinpath(dirname(settings_file), "..", "MHC_pseudo.dat")
                if isfile(pseudo_src)
                    cp(pseudo_src, pseudo_dst; force=true)
                    println("Copied MHC_pseudo.dat from $pseudo_src to $pseudo_dst")
                else
                    println("Warning: Could not find MHC_pseudo.dat at $pseudo_src")
                end
                return p
            else
                error("ERROR: Path '$p' does not exist.")
            end
        end
    end

    error("ERROR: No valid netMHCpan path found in settings.txt.")
end

function main()
    args = parse_arguments()
    folder_path = args["folder"]

    settings_file = find_settings_file()
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

    # Read alleles robustly (handle single-line or multi-line)
    allele_list = open(alleles_file) do file
        lines = [strip(line) for line in readlines(file) if !isempty(strip(line))]
        if length(lines) == 1 && occursin(",", lines[1])
            # Single comma-separated line
            [strip(a) for a in split(lines[1], ",")]
        else
            # One allele per line
            [strip(split(line, r"\s+")[1]) for line in lines]
        end
    end

    # Define pseudo_file before using it
    pseudo_file = joinpath(dirname(settings_file), "MHC_pseudo.dat")
    if !isfile(pseudo_file)
        error("Could not find MHC_pseudo.dat at $pseudo_file")
    end
    # Build a map from normalized allele (no star, no colon) to the format found in pseudo list
    allele_format_map = Dict{String,String}()
    for line in readlines(pseudo_file)
        s = strip(line)
        if isempty(s) || startswith(s, "#")
            continue
        end
        allele = split(s)[1]
        # Normalize: remove star and colon for lookup
        key = replace(replace(allele, "*" => ""), ":" => "")
        allele_format_map[key] = allele
    end
    valid_alleles = Set(values(allele_format_map))

    normalized_valid_alleles = Set{String}(normalize_allele(a) for a in valid_alleles)
    colon_valid_alleles = Set{String}(to_colon_format(normalize_allele(a)) for a in valid_alleles)
    filtered_alleles = String[]
    skipped_alleles = String[]
    for allele in allele_list
        # Lookup normalized allele in map
        key = replace(replace(allele, "*" => ""), ":" => "")
        if haskey(allele_format_map, key)
            push!(filtered_alleles, allele_format_map[key])
        else
            push!(skipped_alleles, allele)
        end
    end
    if !isempty(skipped_alleles)
        println("Warning: The following alleles were not found in MHC_pseudo.dat and will be skipped:")
        println(join(skipped_alleles, ", "))
    end
    # Only chunk supported alleles
    allele_list = filtered_alleles

    # ...existing code...

    # At the end, print how many alleles were skipped
    println("NetMHCpan run complete.")
    println("Number of alleles skipped (not found in MHC_pseudo.dat): ", length(skipped_alleles))

    # Chunk alleles by total character length (≤1024 chars for -a argument)
    max_chars = 1024
    allele_chunks = Vector{Vector{String}}()
    current_chunk = String[]
    current_len = 0
    for allele in allele_list
        # +1 for comma separator (except first allele)
        add_len = length(allele) + (isempty(current_chunk) ? 0 : 1)
        if current_len + add_len > max_chars
            push!(allele_chunks, current_chunk)
            current_chunk = String[]
            current_len = 0
        end
        push!(current_chunk, allele)
        current_len += add_len
    end
    if !isempty(current_chunk)
        push!(allele_chunks, current_chunk)
    end
    temp_files = String[]

    for (i, chunk) in enumerate(allele_chunks)
        if isempty(chunk)
            println("Skipping empty chunk $(i)")
            continue
        end
            # Convert alleles in chunk to colon format for NetMHCpan
        # For each allele, submit both compact and colon formats if either is present in MHC_pseudo.dat
        submitted_alleles = String[]
        for a in chunk
            # a is already in the correct format from filtered_alleles
            if !(a in submitted_alleles)
                push!(submitted_alleles, a)
            end
        end
        println("Chunk $(i): $(length(submitted_alleles)) alleles submitted")
        println("Allele string for chunk $(i): ", join(submitted_alleles, ","))
        alleles = join(submitted_alleles, ",")
        chunk_out = joinpath(folder_path, "netMHCpan_output_chunk$(i).tsv")
        cmd = Cmd([netMHCpan_exe, "-p", peptides_file, "-xls", "-a", alleles, "-xlsfile", chunk_out])
        println("Running NetMHCpan chunk $(i) with command:")
        println(cmd)
        try
            run(cmd)
            if isfile(chunk_out)
                println("NetMHCpan output found at $chunk_out")
                push!(temp_files, chunk_out)
            else
                println("ERROR: NetMHCpan did not create the expected output file for chunk $(i). Allele string: ", alleles)
                error("ERROR: NetMHCpan did not create the expected output file for chunk $(i).")
            end
            println("NetMHCpan chunk $(i) completed successfully.")
        catch e
            println("Error running NetMHCpan chunk $(i): ", e)
            println("Allele string for failed chunk $(i): ", alleles)
            exit(1)
        end
    end

    # Join chunk outputs
    println("Joining NetMHCpan outputs...")
    open(xlsfile_path, "w") do out_io
        for (i, temp_file) in enumerate(temp_files)
            open(temp_file, "r") do in_io
                for (j, line) in enumerate(eachline(in_io))
                    # Write header only for first chunk
                    if i == 1 || j > 1
                        write(out_io, line, "\n")
                    end
                end
            end
            rm(temp_file)
        end
    end
    println("NetMHCpan outputs joined to $xlsfile_path")
end

main()
