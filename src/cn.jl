# 1D linear
function fdm_1d_linear(c0::Vector{Float64}, x::Vector{Float64}, tspan::Vector{Float64}, vspan::Vector{Float64}, v_func::Function,  d::Float64, k_func::Function)
    # c0 is a vector of initial concentrations
    # x is a vector of spatial points (allows non-uniform grids)
    # tspan is a vector of time points (allows non-uniform time steps)
    # vspan is a vector of advection velocities at each spatial point
    # d is the diffusion coefficient (assumed constant here)
    # k_func is a function k(t) returning reaction rate at time t

    n = length(c0) # number of spatial points
    c = copy(c0) # initial concentrations
    results = zeros(n, length(tspan))

    for (j, t) in enumerate(tspan)
        results[:, j] .= c
        dt = j == 1 ? tspan[1] : tspan[j] - tspan[j-1] # time step

        #reaction term
        k = k_func(t)
        # Explicit Euler decay (no diffusion yet)

        for i in 2:n-1
            # Advection term (upwind scheme)
            vspan_i = v_func(x[i],t)
            advective_flux = vspan_i * (c[i] - c[i-1]) / (x[i] - x[i-1])
            # Diffusion term (central difference)
            diffusive_flux = -d * (c[i+1] - 2c[i] + c[i-1]) / ((x[i+1] - x[i]) * (x[i] - x[i-1]))
            # Reaction term
            reaction_term = k * c[i]
            # Update concentration
            c[i] += dt * (advective_flux + diffusive_flux + reaction_term)
        end

        # Keep boundary values
        c[1] = c[1]
        c[end] = c[end]
    end

    return results
end

# 1D non-linear
function fdm_1d_nonlinear(c0:: Vector{Float64})
    # Placeholder for 1D non-linear FDM implementation
    return nothing
end

# 2D linear
function fdm_2d_linear(c0:: Matrix{Float64})
    # Placeholder for 2D linear FDM implementation
    return nothing
end

# 2D non-linear
function fdm_2d_nonlinear(c0:: Matrix{Float64})
    # Placeholder for 2D non-linear FDM implementation
    return nothing
end