# --------------------------------
# ADVECTION
# --------------------------------
# velocity (cm/min)
vx = 0.1
vy = 0

# --------------------------------
# DIFFUSION
# --------------------------------
# diffusion coefficient
D = 0.3

# --------------------------------
# PHOTODEGRADATION
# --------------------------------
# ligth sensitivity coefficient
λ = 0.05
# light intensity decline rate
γ = 0.2
# ambient light intensity
I0 = 1.0

# --------------------------------
# GRID PARAMETERS
# --------------------------------
# spatial domain (cm)
x0 = 0.0
xL = 10.0
y0 = 0.0
yL = 10.0
# temporal domain (m)
t0 = 0.0
tF = 5.0
# grid size
Nx = 101
Ny = 101
Nt = 200

# --------------------------------
# INITIAL CONDITIONS
# --------------------------------
# initial concentration function
u0(x, y) = exp(-0.1 * ((x - 5.0)^2 + (y - 5.0)^2))
