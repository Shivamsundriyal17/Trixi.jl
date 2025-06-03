using Trixi
using OrdinaryDiffEq
using LinearAlgebra
using Plots

# switch to turn damping on/off
damping = true

# initial condition with non-zero moments in x_2 component to start with a bended beam
function initial_condition(x, t, ::Damped_Full_Hyperbolic)
    return SVector{12}([zeros(10); 1.0; 0.0])
end

# length of the beam
L = 4.0

# flexibility, mass matrix, and damping matrix describing the beam's materials
flexMat = inv(diagm([1.0e3, 1.0e3, 1.0e3, 500.0, 500.0, 500.0]))
massMat = diagm([1.0, 1.0, 1.0, 20.0, 10.0, 10.0])
dampingMat = inv(diagm([1.0e3, 1.0e3, 1.0e3, 10.0, 10.0, 10.0]))
# no damping if the switch is set to false
if !damping
    dampingMat = 0.0 * dampingMat
end

# zero initial curvature
k_init = [0.0 0.0 0.0]

# zero external forces/moments
f_ext_func(x, t) = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0]

# set up hyperbolic equations
equations_hyperbolic = Damped_Full_Hyperbolic(flexMat, massMat, dampingMat, k_init, f_ext_func)

# set up parabolic equations based on hyperbolic equations
equations_parabolic = Damped_Full_Parabolic(equations_hyperbolic,
                                            dampingMat,
                                            flexMat,
                                            massMat,
                                            equations_hyperbolic.Gamma_inv,
                                            equations_hyperbolic.Psi,
                                            k_init)

# divergence flux at boundaries
@inline function boundary_condition_damped(flux_inner, 
                                           u_inner,
                                           normal::Int64,
                                           direction::Int64,
                                           x,
                                           t,
                                           operator_type::Trixi.Divergence,
                                           equations_parabolic::Damped_Full_Parabolic)
    if direction == 1
        return [flux_inner[1:6]; zeros(6)]
    elseif direction == 2
        return [zeros(6); flux_inner[7:12]]
    else
        error("Invalid direction specified. Use 1 for left boundary and 2 for right boundary.")
    end
end

# gradient flux at boundaries
@inline function boundary_condition_damped(flux_inner,
                                           u_inner,
                                           normal::Any,
                                           direction::Int64,
                                           x,
                                           t,
                                           operator_type::Trixi.Gradient,
                                           equations_parabolic::Damped_Full_Parabolic)
    if direction == 1
        return [zeros(6); u_inner[7:12]]
    elseif direction == 2
        return [u_inner[1:6]; u_inner[7:12]]
    else
        error("Invalid direction specified. Use 1 for left boundary and 2 for right boundary.")
    end
end

boundary_conditions_hyperbolic = boundary_conditions_ib     
boundary_conditions_parabolic = boundary_condition_damped

# define Discontinuous Galerkin Spectral Element Method (DGSEM) solver
solver = DGSEM(polydeg = 3,
               surface_flux = flux_upwind)

# set coordinates for the mesh
coordinates_min = (0.0,)
coordinates_max = (L,)

# create the mesh using TreeMesh
mesh = TreeMesh(coordinates_min,
                coordinates_max,
                initial_refinement_level = 3,
                n_cells_max = 10_000,
                periodicity = false)

# set up semidiscretization object with hyperbolic and parabolic equations
semi = SemidiscretizationHyperbolicParabolic(mesh,
                                             (equations_hyperbolic, equations_parabolic),
                                             initial_condition,
                                             solver,
                                             boundary_conditions = (boundary_conditions_hyperbolic,
                                                                    boundary_conditions_parabolic))
# time span for the simulation
tspan = (0.0, 10.0)

# create ode object for DGSEM semi discretization
ode = semidiscretize(semi, tspan)

# callbacks for the simulation
stepsize_callback = StepsizeCallback(cfl=0.5)
alive_callback = AliveCallback(alive_interval = 1000)
callbacks = CallbackSet(stepsize_callback, alive_callback)

# solve the ODE with specified options
sol = solve(ode,
            CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0,                                             # dummy timestep to satisfy interface (will be overwritten by cfl callback)
            save_everystep = true,
            callback = callbacks)

# spatial nodes
x = semi.cache.elements._node_coordinates
Nx = length(x)

# time nodes
t = sol.t
Nt = length(t)

# reshape solution vector into 3D array (variable, node, time)
U = Array{Float64}(undef, 12, Nx, Nt)
for j in 1:Nt
    U[:, :, j] = reshape_solution(sol.u[j])
end

# compute strains for every time step
γ = Array{Float64}(undef, 3, Nx, Nt)
κ = Array{Float64}(undef, 3, Nx, Nt)
for j in 1:Nt
  strains = flexMat * U[7:12, :, j]

  γ[:, :, j] = strains[1:3, :]
  κ[:, :, j] = strains[4:6, :]
end

# determine position line of beam
S0 = diagm([1.0; 1.0; 1.0])
r = Array{Float64}(undef, 3, Nx, Nt)
# helper to represent zero curvature (k_init) as a function of x because
# calc_deformation_trapez requires this
k_func(x) = [0.0 0.0 0.0]
for j in 1:Nt
    # compute position line from strains
    r[:, :, j] = calc_deformation_trapez(S0, x, k_func, κ[:, :, j], γ[:, :, j])
end

# plot animation in x1-x3-plane
anim = @animate for n in 1:50:Nt
  plot(r[1, :, n], r[3, :, n],
       xlims = (-0.1, 4.5),
       ylims = (-0.02, 0.02),
       xlabel = "[m]",
       ylabel = "[m]",
       legend = false,
       linewidth = 3,
       title = string("t = ", round(t[n], digits = 2), "s"))
end
if damping
    prefix = "damped"
else
    prefix = "undamped"
end
gif(anim, prefix*"_swinging_beam.gif", fps=100)
