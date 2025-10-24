using Test

@testset "no DIAG prints in source files" begin
    src_dir = joinpath(@__DIR__, "..", "..", "src")
    # collect .jl source files recursively under src/ using walkdir for compatibility
    files = String[]
    for (root, _dirs, names) in walkdir(src_dir)
        for name in names
            if endswith(name, ".jl")
                push!(files, joinpath(root, name))
            end
        end
    end
    offenders = String[]
    for f in files
        try
            txt = read(f, String)
            if occursin("[DIAG]", txt)
                push!(offenders, f)
            end
        catch e
            @warn "Could not read file during DIAG scan" file=f error=e
        end
    end

    if !isempty(offenders)
        @info "Found [DIAG] diagnostics in files: $(join(offenders, ", "))"
    end
    @test isempty(offenders)
end
