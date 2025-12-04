# utils.jl

function compute_mass(c::Vector{Float64}, dx::Float64)
    sum(c) * dx
end

function compute_half_life(c0::Vector{Float64}, results::Array{Float64,2},
                           dx::Float64, tspan)
    tvec = collect(tspan)
    m0 = compute_mass(c0, dx)

    for (j, t) in enumerate(tvec)
        if compute_mass(results[:, j], dx) <= 0.5 * m0
            return t
        end
    end

    return nothing
end
