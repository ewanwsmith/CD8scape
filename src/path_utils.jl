"""
path_utils.jl

Shared helpers for handling suffixed filenames and latest-file discovery.

Convention: insert suffix before the extension using an underscore,
e.g., name_suffix.ext. When `suffix == ""`, paths are left unchanged.

Discovery rules (for reads when no suffix provided):
- If exactly one candidate exists for base name, use it.
- If multiple candidates and `latest=true`, pick the most recent by mtime.
- If multiple candidates and `latest=false`, throw with a helpful message.
"""

import Base.Filesystem: isfile

"""
Insert `_<suffix>` before the extension of `path` when `suffix` is non-empty.
Returns `path` unchanged if `suffix == ""`.
"""
function with_suffix(path::AbstractString, suffix::AbstractString="")::String
    if isempty(suffix)
        return String(path)
    end
    dir = dirname(path)
    base = basename(path)
    name, ext = splitext(base)
    return joinpath(dir, string(name, "_", suffix, ext))
end

"""
Return all candidate files in `dir` that match base name `name` and extension
`ext` (e.g., frames + .csv). Candidates include the exact base (no suffix)
and any files of the form `name_<anything><ext>`.
"""
function _candidates_for(dir::AbstractString, name::AbstractString, ext::AbstractString)::Vector{String}
    paths = String[]
    # Exact base file
    base = joinpath(dir, string(name, ext))
    if isfile(base)
        push!(paths, base)
    end
    # Suffixed variants
    for f in readdir(dir)
        # skip directories
        full = joinpath(dir, f)
        if !isfile(full)
            continue
        end
        # must end with ext and start with name followed by "_"
        if endswith(f, ext) && startswith(f, string(name, "_"))
            push!(paths, full)
        end
    end
    # Deduplicate while preserving order
    return unique(paths)
end

"""
Discover a file for the given `base_path` (e.g., "/.../frames.csv").
When `latest=true` and multiple candidates exist, choose the one with
the most recent modification time. Otherwise, error on ambiguity.
"""
function discover_path(base_path::AbstractString; latest::Bool=false)::String
    dir = dirname(base_path)
    base = basename(base_path)
    name, ext = splitext(base)
    candidates = _candidates_for(dir, name, ext)
    if isempty(candidates)
        error("No file found for base $(base) in $(dir). Expected either $(base) or $(name)_<suffix>$(ext)")
    elseif length(candidates) == 1
        return candidates[1]
    else
        if latest
            # Pick by most recent mtime
            mtimes = map(p -> stat(p).mtime, candidates)
            # argmax over mtimes
            idx = argmax(mtimes)
            return candidates[idx]
        else
            # Build message listing candidates
            msg = "Multiple candidates found for $(base) in $(dir):\n" * join(candidates, "\n") * "\nProvide a suffix or set latest=true."
            error(msg)
        end
    end
end

"""
Resolve a path for reading given `base_path`, optional `suffix` and `latest`.
- If `suffix != ""`, return `with_suffix(base_path, suffix)` (must exist).
- Else, return `discover_path(base_path; latest=latest)`.
"""
function resolve_read(base_path::AbstractString; suffix::AbstractString="", latest::Bool=false)::String
    if !isempty(suffix)
        p = with_suffix(base_path, suffix)
        if !isfile(p)
            error("Expected file not found: $(p).")
        end
        return p
    else
        return discover_path(base_path; latest=latest)
    end
end

"""
Resolve a path for writing given `base_path` and optional `suffix`.
This does not validate existence; it only builds the path.
"""
function resolve_write(base_path::AbstractString; suffix::AbstractString="")::String
    return with_suffix(base_path, suffix)
end
