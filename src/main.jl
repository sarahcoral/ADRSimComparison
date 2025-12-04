using .MicroplasticsDegradation

println("Running comparisons: CN, UPFDM, EUPFDM, SixOrdFDM")

# domain
n = 100 # number of spatial points (for evenly spaced grid)
x = range(0, 1, length=n) # spatial grid
dx = x[2] - x[1] 
c0 = ones(n)

# time
tspan = 0:0.1:10

# diffusion
D = 0.0001

# velocity
v = 0.01
vspan = fill(v, n)
v_func = (x,t) -> v # constant velocity

# reaction rate model
k_func = t -> degradation_rate(:HDPE, 20.0, 0.0)

# FDM
fdm = solve_fdm(c0, x, tspan, vspan, v_func, D, k_func)
println("FDM half-life: ", compute_half_life(c0, fdm, dx, tspan))
