# 1D linear
function upfdm_1d_linear(Ua, D, f, x, t, u0, ua, ub)
    Nx = length(x)
    Nt = length(t)
    dx = x[2] - x[1]
    dt = t[2] - t[1]

    λ = dt/dx
    β = dt/dx^2

    u = zeros(Float64, Nt, Nx)
    u[1,:] = u0

    for n in 1:Nt-1
        # boundary conditions
        u[n,1]   = ua(t[n])
        u[n,end] = ub(t[n])

        for i in 2:Nx-1
            u[n+1,i] = (u[n,i] + (λ+β)*u[n,i-1] + β*u[n,i+1]) /
                        (1 + dt + λ + 2β)
        end

        # right boundary: Neumann ux(b,t) = -u(b,t)
        u[n+1,end] = u[n+1,end-1] - dx*u[n+1,end]
    end

    return u
end


# 1D non-linear
function upfdm_1d_nonlinear(D, r, x, t, u0)
    Nx = length(x)
    Nt = length(t)
    dx = x[2] - x[1]
    dt = t[2] - t[1]

    λ = D*dt/dx^2

    u = zeros(Float64, Nt, Nx)
    u[1,:] = u0

    for n in 1:Nt-1
        u[n,1] = 0
        u[n,end] = 0

        for i in 2:Nx-1
            num = (1 + r*dt)*u[n,i] + λ*(u[n,i+1] - u[n,i])
            den = 1 + 2λ + r*dt*u[n,i]
            u[n+1,i] = num/den
        end
    end

    return u
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