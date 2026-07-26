using LinearAlgebra: Diagonal
using OrdinaryDiffEqLowStorageRK
using Trixi

beam_length = 1.0
flexibility_matrix = Diagonal(1.0 ./ [1.0e3, 1.0e3, 1.0e3,
                                  500.0, 500.0, 500.0])
mass_matrix = Diagonal([1.0, 1.0, 1.0, 20.0, 10.0, 10.0])
damping_matrix = 2.0 *
                 Diagonal(1.0 ./ [1.0e3, 1.0e3, 1.0e3,
                              10.0, 10.0, 10.0])

equations_hyperbolic = DampedIntrinsicBeamEquations1D(;
                                                      mass_matrix,
                                                      flexibility_matrix,
                                                      damping_matrix,
                                                      initial_curvature = zeros(3))
equations_parabolic = DampedIntrinsicBeamDiffusion1D(equations_hyperbolic)

function initial_condition(x, t, equations::DampedIntrinsicBeamEquations1D)
    axial_resultant = 0.5 < x[1] < 0.75 ? 3000.0 : 0.0
    return SVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   axial_resultant, 0.0, 0.0, 0.0, 0.0, 0.0)
end

boundary_condition = BoundaryConditionDampedIntrinsicBeam()
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

tspan = (0.0, 10.0)
ode = semidiscretize(semi, tspan)

cfl = 0.01
stepsize_callback = StepsizeCallback(cfl = cfl)
save_times = sort!(unique!(vcat(collect(range(0.0, 0.1; length = 101)),
                                collect(range(0.1, 1.0; length = 91)),
                                collect(range(1.0, last(tspan); length = 181)))))
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0,
            adaptive = false,
            callback = stepsize_callback,
            saveat = save_times,
            ode_default_options()...)

@assert all(isfinite, sol.u[end])
energy_history = [Trixi.integrate(entropy, state, semi; normalize = false)
                  for state in sol.u]
maximum_axial_resultant = [maximum(abs, @view(state[7:12:length(state)]))
                           for state in sol.u]
maximum_axial_strain = flexibility_matrix[1, 1] .* maximum_axial_resultant
