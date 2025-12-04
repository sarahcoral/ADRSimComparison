module MicroplasticsDegradation

include("cn_1d_linear.jl")
include("cn_1d_nonlinear.jl")
include("cn_2d_linear.jl")
include("cn_2d_nonlinear.jl")
include("up_1d_linear.jl")
include("up_1d_nonlinear.jl")
include("up_2d_linear.jl")
include("up_2d_nonlinear.jl")
include("utils.jl")
include("models.jl")
include("test.jl")

export degradation_rate,
       cn_1d_linear,
       cn_1d_nonlinear,
       cn_2d_linear,
       cn_2d_nonlinear,
       up_1d_linear,
       up_1d_nonlinear,
       up_2d_linear,
       up_2d_nonlinear,
       test,
       compute_mass,
       compute_half_life

end # module