using Trixi
using OrdinaryDiffEq
using LinearAlgebra
using Plots

# initial condition with non-zero moments in x_2 component to start with a bended beam
function initial_condition(x, t, ::IntrinsicBeamEquation)
    return SVector{12}([zeros(4); 1.0; zeros(7)])
end

# domain (beam length)
L = 4.0

# flexibility, mass matrix
flexMat = inv(diagm([1.0e3, 1.0e3, 1.0e3, 500.0, 500.0, 500.0]))
massMat = diagm([1.0, 1.0, 1.0, 20.0, 10.0, 10.0])

# initial curvature, external forces/moments
k_func(x) = [0.0; 0.0; 0.0]
f_ext_func(x,t) = [0.0; 0.0; 0.0]
m_ext_func(x,t) = [0.0; 0.0; 0.0]

equation = IntrinsicBeamEquation(flexMat, massMat, k_func, f_ext_func, m_ext_func)

solver = DGSEM(polydeg=3, surface_flux=flux_upwind)

coordinates_min = (0.0, )
coordinates_max = (L, )
# create the mesh using TreeMesh
mesh = TreeMesh(coordinates_min,
                coordinates_max,
                initial_refinement_level = 3,
                n_cells_max = 10_000,
                periodicity = false)

source_term = source_term_intrinsic_beam_diag
boundary_conditions = boundary_conditions_ib

semi = SemidiscretizationHyperbolic(mesh,
                                    equation,
                                    initial_condition,
                                    solver,
                                    boundary_conditions=boundary_conditions,
                                    source_terms=source_term)

# simulation time 
T = 10.0
tspan = (0.0, T)
cfl = 0.5

ode = semidiscretize(semi, tspan)

stepsize_callback = StepsizeCallback(cfl=cfl)
alive_callback = AliveCallback(alive_interval = 1000)
callbacks = CallbackSet(stepsize_callback, alive_callback)

sol = solve(ode,
            CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0,                                             # dummy timestep to satisfy interface (will be overwritten by cfl callback)
            save_everystep=true,
            callback=callbacks)

# spatial nodes
x = semi.cache.elements._node_coordinates
Nx = length(x)

# time nodes
t = sol.t
Nt = length(t)

# reshape solution vector into 3D array (variable, node, time)
U = Array{Float64}(undef, 12, Nx, Nt)
for j in 1:Nt
    U[:,:,j] = reshape_solution(sol.u[j])
end

# apply constitutive laws to determine strains
γ = Array{Float64}(undef, 3, Nx, Nt)
κ = Array{Float64}(undef, 3, Nx, Nt)
for j in 1:Nt
    γ[:,:,j], κ[:,:,j], _, _ = constitutive_laws(U[:,:,j], flexMat, massMat)
end

# set initital value for recursion to determine rotation matrix
S0 = diagm([1.0; 1.0; 1.0])

# determine position line of beam
r = Array{Float64}(undef, 3, Nx, Nt)
for j in 1:Nt
    r[:,:,j] = calc_deformation_trapez(S0, x, k_func, κ[:,:,j], γ[:,:,j])
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
gif(anim, "undamped_swinging_beam.gif", fps=100)