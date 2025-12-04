using SparseArrays, LinearAlgebra
function cn_2d_linear(
    vx, vy, Dx, Dy, R, xspan, yspan, tspan, Nx, Ny, Nt, u0;
    bc_left = 0.0, bc_right = 0.0, bc_top = 0.0, bc_bottom = 0.0
)
    # Discretization
    x0, xL = xspan
    y0, yL = yspan
    t0, tF = tspan

    dx = (xL - x0) / (Nx - 1)
    dy = (yL - y0) / (Ny - 1)
    dt = (tF - t0) / (Nt - 1)

    x = range(x0, xL, length=Nx)
    y = range(y0, yL, length=Ny)
    tvals = range(t0, tF, length=Nt)

    # Solution storage
    U = zeros(Float64, Nx, Ny, Nt)
    #U[:, :, 1] = u0.(x, y)
    # Initialize the solution at t=0 using a 2D grid evaluation of the initial condition
    for i in 1:Nx
        for j in 1:Ny
            U[i, j, 1] = u0(x[i], y[j])
        end
    end

    # Coefficients for diffusion and advection
    αx = Dx / dx^2
    βx = vx / (2 * dx)
    αy = Dy / dy^2
    βy = vy / (2 * dy)

    # Crank-Nicolson matrices for x and y directions
    # x-direction (diffusion + advection)
    lower_x = zeros(Float64, Nx-1)
    main_x  = zeros(Float64, Nx)
    upper_x = zeros(Float64, Nx-1)

    for i in 2:Nx-1
        lower_x[i-1] = αx - βx
        main_x[i]    = -2*αx
        upper_x[i]   = αx + βx
    end

    A_x = spdiagm(-1 => lower_x, 0 => main_x, 1 => upper_x)

    # y-direction (diffusion + advection)
    lower_y = zeros(Float64, Ny-1)
    main_y  = zeros(Float64, Ny)
    upper_y = zeros(Float64, Ny-1)

    for i in 2:Ny-1
        lower_y[i-1] = αy - βy
        main_y[i]    = -2*αy
        upper_y[i]   = αy + βy
    end

    A_y = spdiagm(-1 => lower_y, 0 => main_y, 1 => upper_y)

    # Crank-Nicolson operator matrices
    L_x = I - dt / 2 * A_x     # implicit operator for x direction
    R_x = I + dt / 2 * A_x     # explicit operator for x direction

    L_y = I - dt / 2 * A_y     # implicit operator for y direction
    R_y = I + dt / 2 * A_y     # explicit operator for y direction

    # Time stepping
    # Precompute LU factorizations for speed
    Fx = lu(L_x)      # L_x = I - dt/2 * A_x   for half step
    Fy = lu(L_y)      # L_y = I - dt   * A_y   for full-step (we scale below)

    # Time stepping with Strang splitting
    for n in 1:Nt-1
        un = U[:,:,n]

        # --- 1. First half x-step ---
        # Solve   L_x * u_half = R_x * un
        u_half = similar(un)
        for j in 1:Ny
            rhs = R_x * un[:,j]     # rhs :: Vector{Float64} (Nx)
            u_half[:,j] = Fx \ rhs  # solve Nx-system per y-slice
        end

        # --- 2. Full y-step ---
        # Solve   L_y * u_star = R_y * u_half
        u_star = similar(un)
        for i in 1:Nx
            rhs = R_y * u_half[i,:]    # rhs :: Vector{Float64} (Ny)
            u_star[i,:] = Fy \ rhs
        end

        # --- 3. Reaction step (explicit Euler or exponential) ---
        # u_react = exp(R*dt) * u_star   (exact)
        # or use explicit: u_react = u_star + dt * R * u_star
        u_react = u_star .* (1 + dt*R)    # stable for constant R

        # --- 4. Second half x-step ---
        unew = similar(un)
        for j in 1:Ny
            rhs = R_x * u_react[:,j]
            unew[:,j] = Fx \ rhs
        end

        # Apply Dirichlet BCs
        unew[1,:]  .= bc_left
        unew[end,:].= bc_right
        unew[:,1]  .= bc_bottom
        unew[:,end].= bc_top

        U[:,:,n+1] = unew
    end

    return U
end