using .MicroplasticsDegradation

# 1D domain
n = 100
x = range(0.0, 1.0, length=n)
dx = x[2] - x[1]
c0 = ones(n)

# Time
tspan = 0.0:0.1:10.0

# Example degradation function (HDPE, T=20, S=0)
k_func = t -> degradation_rate(:HDPE, 20.0, 0.0)

# Solve FDM
fdm_results = solve_fdm(c0, dx, 0.1, tspan, k_func)
half_life_fdm = compute_half_life(c0, fdm_results, dx, tspan)
println("FDM half-life: $half_life_fdm")

# Solve FEM (placeholder, uses FDM for now)
fem_results = solve_fem(x, c0, tspan, k_func)
half_life_fem = compute_half_life(c0, fem_results, dx, tspan)
println("FEM half-life: $half_life_fem")

# Particle tracking
particles = [Particle(xi, 1.0) for xi in x]
for t in tspan
    evolve_particles!(particles, 0.1, k_func)
end
total_mass = sum(p.m for p in particles)
println("Particle total mass: $total_mass")
