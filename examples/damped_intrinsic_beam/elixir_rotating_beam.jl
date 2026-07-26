using LinearAlgebra: Diagonal
using OrdinaryDiffEqLowStorageRK
using Trixi

# Full nonlinear spin-up from rest; no steady-branch reduction is imposed.
beam_length = 4.0
flexibility_matrix = Diagonal(1.0 ./ [1.0e3, 1.0e3, 1.0e3,
                                  500.0, 500.0, 500.0])
mass_matrix = Diagonal([1.0, 1.0, 1.0, 20.0, 10.0, 10.0])
damping_matrix = Diagonal(1.0 ./ [1.0e3, 1.0e3, 1.0e3,
                              10.0, 10.0, 10.0])

terminal_angular_speed = 2.8523
ramp_duration = 1.0

function initial_condition(x, t, equations::DampedIntrinsicBeamEquations1D)
    return zero(SVector{12, eltype(x)})
end

function root_velocity(x, t, equations)
    ramp = min(max(t / ramp_duration, zero(t)), one(t))
    return SVector(0.0, 0.0, 0.0, 0.0, 0.0,
                   terminal_angular_speed * ramp)
end

equations_hyperbolic = DampedIntrinsicBeamEquations1D(;
                                                      mass_matrix,
                                                      flexibility_matrix,
                                                      damping_matrix,
                                                      initial_curvature = zeros(3))
equations_parabolic = DampedIntrinsicBeamDiffusion1D(equations_hyperbolic)

boundary_condition = BoundaryConditionDampedIntrinsicBeam(;
                                                          left_velocity = root_velocity)
boundary_conditions = (boundary_condition, boundary_condition)

polydeg = 3
solver = DGSEM(polydeg = polydeg, surface_flux = flux_upwind)
solver_parabolic = ViscousFormulationLocalDG()
initial_refinement_level = 3
mesh = TreeMesh((0.0,), (beam_length,);
                initial_refinement_level,
                n_cells_max = 10_000,
                periodicity = false)

semi = SemidiscretizationHyperbolicParabolic(mesh,
                                             (equations_hyperbolic,
                                              equations_parabolic),
                                             initial_condition, solver;
                                             solver_parabolic,
                                             boundary_conditions)

tspan = (0.0, 40.0)
ode = semidiscretize(semi, tspan)

# The paper protocol uses k=3, eight cells, alternating auxiliary traces, and
# the conservative CFL 0.01. The first dt is a dummy value replaced by the
# callback before the first step.
cfl = 0.01
stepsize_callback = StepsizeCallback(cfl = cfl)
save_times = range(first(tspan), last(tspan); length = 401)
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0,
            adaptive = false,
            callback = stepsize_callback,
            saveat = save_times,
            ode_default_options()...)

@assert all(isfinite, sol.u[end])
energy_history = [Trixi.integrate(entropy, state, semi; normalize = false)
                  for state in sol.u]

# Analytic terminal straight branch used only for post-processing/comparison.
function steady_rotating_solution(x)
    a = terminal_angular_speed *
        sqrt(mass_matrix[2, 2] * flexibility_matrix[1, 1])
    b = sqrt(flexibility_matrix[1, 1] / mass_matrix[2, 2])
    coefficient = inv(flexibility_matrix[1, 1] *
                      cos(a * beam_length))
    f1 = coefficient * cos(a * x) - inv(flexibility_matrix[1, 1])
    v2 = coefficient * b * sin(a * x)
    return SVector(0.0, v2, 0.0, 0.0, 0.0, terminal_angular_speed,
                   f1, 0.0, 0.0, 0.0, 0.0, 0.0)
end
