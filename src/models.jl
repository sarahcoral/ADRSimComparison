# models.jl

"""
    degradation_rate(material, T, S)

Compute a simple degradation rate:
- material: :HDPE, :LDPE, :PLA
- T: temperature (°C)
- S: sun exposure factor (0 = none, 1 = high)
"""
function degradation_rate(material::Symbol, T::Float64, S::Float64)
    base = material == :HDPE ? 1e-6 :
           material == :LDPE ? 3e-6 :
           material == :PLA  ? 1e-5 : 1e-6

    Ea = 5000.0        # activation energy (example)
    R = 8.314
    kT = base * exp(-Ea / (R * (T + 273.15)))

    return kT * (1 + S)
end 