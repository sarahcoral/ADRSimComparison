# 1D linear
function upfdm_1d_linear(c0::Vector{Float64}, D::Float64, dx::Float64, dt::Float64,
                     tspan, k_func::Function)

    tvec = collect(tspan)
    n = length(c0)
    c = copy(c0)
    results = zeros(n, length(tvec))

    α = D * dt / dx^2  # diffusion number

    for (j, t) in enumerate(tvec)
        results[:, j] .= c
        k = k_func(t)

        new_c = similar(c)

        for i in 2:n-1
            num = c[i] + α * (c[i+1] - 2c[i] + c[i-1])
            denom = 1 + dt * k
            new_c[i] = num / denom
        end

        # Boundary conditions (Dirichlet)
        new_c[1] = c[1]
        new_c[end] = c[end]

        # Enforce positivity
        @inbounds for i in 1:n
            new_c[i] = max(new_c[i], 0.0)
        end

        c .= new_c
    end

    return results
end

# 1D non-linear
function upfdm_1d_nonlinear(c0:: Vector{Float64})
    # Placeholder for 1D non-linear FDM implementation
    return nothing
end

# 2D linear
function upfdm_2d_linear(c0:: Matrix{Float64})
    # Placeholder for 2D linear FDM implementation
    return nothing
end

# 2D non-linear
function upfdm_2d_nonlinear(c0:: Matrix{Float64})
    # Placeholder for 2D non-linear FDM implementation
    return nothing
end