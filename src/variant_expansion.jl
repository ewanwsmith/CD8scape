module VariantExpansion

export Site, expand_sites

using Random
using Dates

struct Site
    pos::Int
    ref::String
    alts::Vector{String}
    af::Union{Nothing, Vector{Float64}}
end

"""
expand_sites(sites; max_alts_per_site=3, max_combinations=10000, rng=Random.GLOBAL_RNG, sample_on_saturate=true)
Return a Vector of combinations. Each combination is Vector{String} with chosen allele for each site (REF or one of ALTs).
- If the full product <= max_combinations: returns full cartesian expansion.
- Otherwise, if sample_on_saturate==true: returns up to max_combinations sampled combos (probabilities by AF if available).
- Else returns deterministic reduced set using only top-1 ALT per site.
"""
function expand_sites(sites::Vector{Site};
    max_alts_per_site::Int=3,
    max_combinations::Union{Nothing,Int}=nothing,
    rng::AbstractRNG=Random.GLOBAL_RNG,
    sample_on_saturate::Bool=true,
    saturation_log::Union{Nothing,String}=nothing,
    stop_on_saturation::Bool=false)

    # Allow overriding default max_combinations via ENV var CD8S_MAX_VARIANT_COMBINATIONS
    if max_combinations === nothing
        max_combinations = try
            parse(Int, get(ENV, "CD8S_MAX_VARIANT_COMBINATIONS", "10000"))
        catch
            10000
        end
    end

    # Build options per site (include REF as the first option)
    opts = Vector{Vector{String}}(undef, length(sites))
    probs = Vector{Union{Nothing,Vector{Float64}}}(undef, length(sites))
    for (i, s) in enumerate(sites)
        alts_keep = s.alts[1:min(max_alts_per_site, length(s.alts))]
        opts[i] = [s.ref; alts_keep]
        if s.af !== nothing
            afs = s.af[1:min(max_alts_per_site, length(s.alts))]
            ref_af = max(0.0, 1.0 - sum(afs))
            p = [ref_af; afs]
            ps = p ./ sum(p)
            probs[i] = ps
        else
            probs[i] = nothing
        end
    end

    # compute product
    prod_count = prod(length(o) for o in opts)
    if prod_count <= max_combinations
        # full enumeration
        combos = Vector{Vector{String}}()
        idxs = [1 for _ in opts]
        while true
            push!(combos, [opts[i][idxs[i]] for i in 1:length(opts)])
            # increment
            k = length(opts)
            carry = true
            while k >= 1 && carry
                idxs[k] += 1
                if idxs[k] > length(opts[k])
                    idxs[k] = 1
                    k -= 1
                else
                    carry = false
                end
            end
            if carry
                break
            end
        end
        return combos
    else
        # log saturation if requested
        try
            if saturation_log !== nothing
                strategy = stop_on_saturation ? "stop" : (sample_on_saturate ? "sampling" : "fallback")
                sites_pos = join([string(s.pos) for s in sites], ";")
                # Ensure directory exists
                try mkpath(dirname(saturation_log)) catch end
                header = "timestamp,n_sites,product,threshold,strategy,site_positions"
                line = join([
                    Dates.format(Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
                    string(length(sites)),
                    string(prod_count),
                    string(max_combinations),
                    strategy,
                    '"'*sites_pos*'"'
                ], ",")
                # If file doesn't exist, write a creation comment and header first
                if !isfile(saturation_log)
                    open(saturation_log, "w") do io
                        println(io, "# Saturation log created: " * Dates.format(Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
                        println(io, header)
                        println(io, line)
                    end
                else
                    open(saturation_log, "a") do io
                        println(io, "# entry: " * Dates.format(Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
                        println(io, line)
                    end
                end
            end
        catch err
            # ignore logging errors but keep a hint in stderr
            try println(stderr, "saturation_log error: ", err) catch end
        end

        if stop_on_saturation
            # Caller requested that we stop (return no combinations) when saturated.
            return Vector{Vector{String}}()
        elseif sample_on_saturate
            combos_set = Set{Vector{String}}()
            attempts = 0
            max_attempts = max_combinations * 10
            while length(combos_set) < max_combinations && attempts < max_attempts
                choice = Vector{String}(undef, length(opts))
                for i in 1:length(opts)
                    if probs[i] !== nothing
                        p = probs[i]
                        # simple weighted sampling without external packages
                        r = rand(rng)
                        cs = cumsum(p)
                        idx = searchsortedfirst(cs, r)
                        if idx < 1
                            idx = 1
                        end
                        choice[i] = opts[i][idx]
                    else
                        choice[i] = rand(rng, opts[i])
                    end
                end
                push!(combos_set, choice)
                attempts += 1
            end
            return collect(combos_set)
        else
            # deterministic fallback: keep only REF and the top ALT per site
            combos = Vector{Vector{String}}()
            push!(combos, [opts[i][1] for i in 1:length(opts)]) # all REF
            push!(combos, [opts[i][min(2, length(opts[i]))] for i in 1:length(opts)])
            return combos
        end
    end
end

end # module
