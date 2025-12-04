module MicroplasticsDegradation

include("cn.jl")
include("upfdm.jl")
include("eupfdm.jl")
include("utils.jl")
include("models.jl")

export degradation_rate,
       fdm_1d_linear,
       fdm_1d_nonlinear,
       fdm_2d_linear,
       fdm_2d_nonlinear,
       upfdm_1d_linear,
       upfdm_1d_nonlinear,
       upfdm_2d_linear,
       upfdm_2d_nonlinear,
       eupfdm_1d_linear,
       eupfdm_1d_nonlinear,
       eupfdm_2d_linear,
       eupfdm_2d_nonlinear,
       compute_mass,
       compute_half_life

end # module