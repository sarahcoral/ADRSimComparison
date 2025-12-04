using .MicroplasticsDegradation
using DataFrames, CSV

println("Running comparisons: CN, UPFDM")

# Linear ADR

# Parameters
v = 1.0      # advection speed
D = 0.1      # diffusion coefficient
R = -0.2     # reaction rate

xspan = (0.0, 1.0)
tspan = (0.0, 0.5)

Nx = 101
Nt = 200

# Initial condition
u0(x) = exp(-100*(x-0.3)^2)

#U = fdm_1d_linear(v, D, R, xspan, tspan, Nx, Nt, u0)

println("1D linear ADR simulation complete.")

# Parameters
vx = 1.0     # advection speed in x-direction
vy = 1.0     # advection speed in y-direction
Dx = 0.1     # diffusion coefficient in x-direction
Dy = 0.1     # diffusion coefficient in y-direction
R = -0.2     # reaction coefficient (linear)

# Spatial and temporal domains
xspan = (0.0, 1.0)
yspan = (0.0, 1.0)
tspan = (0.0, 0.5)

# Grid size
Nx = 101
Ny = 101
Nt = 200

# Initial condition (Gaussian bump)
u0(x, y) = exp(-100 * ((x - 0.5)^2 + (y - 0.5)^2))
# Solve the PDE
U = cn_2d_linear(vx, vy, Dx, Dy, R, xspan, yspan, tspan, Nx, Ny, Nt, u0)

# U is a 3D array with shape (Nx, Ny, Nt), holding concentration at each grid point over time

println("2D linear ADR simulation complete.")

X, Y, Z = size(U)

df = DataFrame(
    i = repeat(1:X, outer=Y*Z),
    j = repeat(repeat(1:Y, inner=X), outer=Z),
    k = repeat(1:Z, inner=X*Y),
    concentration = vec(U)
)
CSV.write("cn_2d_linear_results.csv", df)