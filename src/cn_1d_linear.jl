using SparseArrays, LinearAlgebra
function cn_1d_linear(
    v, D, R, xspan, tspan, Nx, Nt, u0;
    bc_left = 0.0, bc_right = 0.0
)

    # Discretization
    x0, xL = xspan
    t0, tF = tspan

    dx = (xL - x0) / (Nx - 1)
    dt = (tF - t0) / (Nt - 1)

    x = range(x0, xL, length=Nx)
    U = zeros(Float64, Nx, Nt)
    U[:, 1] = u0.(x)

    # Finite difference coefficients
    α = D / dx^2
    β = v / (2dx)

    # Crank–Nicolson matrices
    # Left (implicit part): (I - dt/2 * A)
    # Right (explicit part): (I + dt/2 * A)
    main = zeros(Float64, Nx)
    upper = zeros(Float64, Nx-1)
    lower = zeros(Float64, Nx-1)

    # Fill A matrix coefficients
    for i in 2:Nx-1
        lower[i-1] =  α - β
        main[i]     = -2α + R
        upper[i]   =  α + β
    end

    # Sparse system matrices
    A = spdiagm(-1 => lower, 0 => main, 1 => upper)

    L = I - dt/2 * A
    Rmat = I + dt/2 * A

    # Time stepping
    for n in 1:Nt-1
        rhs = Rmat * U[:, n]

        # Apply Dirichlet boundary conditions
        rhs[1]  = bc_left
        rhs[end] = bc_right

        # Solve linear system
        U[:, n+1] = L \ rhs

        # Enforce boundary conditions directly
        U[1, n+1] = bc_left
        U[end, n+1] = bc_right
    end

    return U
end
