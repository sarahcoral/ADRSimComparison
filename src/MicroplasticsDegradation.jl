module MicroplasticsDegradation

using DifferentialEquations, LinearAlgebra, Plots

include("models.jl")
include("fdm.jl")
include("fem.jl")
include("particle.jl")
include("utils.jl")

export Particle, degradation_rate, solve_fdm, solve_fem, evolve_particles!, compute_mass, compute_half_life

end