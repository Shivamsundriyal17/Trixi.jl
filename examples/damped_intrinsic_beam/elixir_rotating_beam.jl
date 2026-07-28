using LinearAlgebra: Diagonal, dot
using OrdinaryDiffEqLowStorageRK
using Trixi

# Full nonlinear spin-up from rest; no steady-branch reduction is imposed.
case_name = get(ENV, "ROTATING_CASE", "baseline")
damping_multiplier = parse(Float64,
                           get(ENV, "ROTATING_DAMPING_MULTIPLIER", "1.0"))
steady_initial_condition = lowercase(get(ENV, "ROTATING_STEADY_INITIAL", "false")) in ("1",
                                                                                       "true",
                                                                                       "yes")
record_online_ledger = lowercase(get(ENV, "ROTATING_RECORD_LEDGER", "true")) in ("1",
                                                                                 "true",
                                                                                 "yes")

beam_length = 4.0
flexibility_matrix = Diagonal(1.0 ./ [1.0e3, 1.0e3, 1.0e3,
                                  500.0, 500.0, 500.0])
mass_matrix = Diagonal([1.0, 1.0, 1.0, 20.0, 10.0, 10.0])
damping_matrix = damping_multiplier *
                 Diagonal(1.0 ./ [1.0e3, 1.0e3, 1.0e3,
                              10.0, 10.0, 10.0])

terminal_angular_speed = 2.8523
ramp_duration = 1.0

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

function initial_condition(x, t, equations::DampedIntrinsicBeamEquations1D)
    if steady_initial_condition
        return SVector{12, eltype(x)}(steady_rotating_solution(x[1]))
    else
        return zero(SVector{12, eltype(x)})
    end
end

function root_velocity(x, t, equations)
    ramp = steady_initial_condition ?
           one(t) : min(max(t / ramp_duration, zero(t)), one(t))
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

t_end = parse(Float64,
              get(ENV, "ROTATING_T_END",
                  steady_initial_condition ? "10.0" : "40.0"))
tspan = (0.0, t_end)
ode = semidiscretize(semi, tspan)

# The online ledger accumulates every accepted time step. This avoids the
# quadrature error obtained by integrating only the sparsely saved plot states.
node_coordinates = semi.cache.elements.node_coordinates
element_order = sortperm([sum(node_coordinates[1, :, element])
                          for element in axes(node_coordinates, 3)])
quadrature_weights = solver.basis.weights
volume_jacobians = [abs(inv(semi.cache.elements.inverse_jacobian[element]))
                    for element in axes(node_coordinates, 3)]
capacity_abs_flux = equations_hyperbolic.capacity_matrix *
                    equations_hyperbolic.propagation_matrix_abs
damping_operator_inverse = iszero(damping_multiplier) ?
                           zeros(6, 6) :
                           inv(Matrix(equations_hyperbolic.damping_operator))
ledger_parabolic_workspace = similar(ode.u0)

function beam_quadrature_sum(function_)
    value = 0.0
    for element in axes(node_coordinates, 3), i in eachindex(quadrature_weights)
        value += volume_jacobians[element] * quadrature_weights[i] *
                 function_(i, element)
    end
    return value
end

function rotating_ledger_terms(state, t)
    Trixi.rhs_parabolic!(ledger_parabolic_workspace, state, semi, t)
    state_array = reshape(state, 12, length(quadrature_weights),
                          length(volume_jacobians))
    gradients = semi.cache_parabolic.viscous_container.gradients

    material_dissipation = if iszero(damping_multiplier)
        0.0
    else
        beam_quadrature_sum() do i, element
            state_node = SVector{12}(state_array[:, i, element])
            gradient_node = SVector{12}(gradients[:, i, element])
            damping_resultant = intrinsic_beam_damping_resultant(state_node,
                                                                 gradient_node,
                                                                 equations_parabolic)
            dot(damping_resultant,
                damping_operator_inverse * damping_resultant)
        end
    end

    jump_dissipation = 0.0
    for index in 1:(length(element_order) - 1)
        left_element = element_order[index]
        right_element = element_order[index + 1]
        jump = SVector{12}(state_array[:, end, left_element] -
                           state_array[:, 1, right_element])
        jump_dissipation += 0.5 *
                            dot(jump, capacity_abs_flux * jump)
    end

    left_element = first(element_order)
    right_element = last(element_order)
    left_state = SVector{12}(state_array[:, 1, left_element])
    right_state = SVector{12}(state_array[:, end, right_element])
    left_velocity = SVector{6}(left_state[1:6])
    left_resultant = SVector{6}(left_state[7:12])
    right_resultant = SVector{6}(right_state[7:12])
    left_gradient = SVector{12}(gradients[:, 1, left_element])
    left_damping_resultant = intrinsic_beam_damping_resultant(left_state,
                                                              left_gradient,
                                                              equations_parabolic)
    prescribed_velocity = root_velocity(SVector(0.0), t,
                                        equations_hyperbolic)

    left_boundary_dissipation = dot(left_velocity,
                                    equations_hyperbolic.left_impedance *
                                    left_velocity)
    right_boundary_dissipation = dot(right_resultant,
                                     equations_hyperbolic.right_impedance *
                                     right_resultant)
    physical_root_power = -dot(prescribed_velocity,
                               left_resultant + left_damping_resultant)
    sat_data_power = dot(prescribed_velocity,
                         equations_hyperbolic.left_impedance * left_velocity)

    return SVector(material_dissipation, jump_dissipation,
                   left_boundary_dissipation, right_boundary_dissipation,
                   physical_root_power, sat_data_power)
end

# The paper protocol uses k=3, eight cells, alternating auxiliary traces, and
# the conservative CFL 0.01. The first dt is a dummy value replaced by the
# callback before the first step.
cfl = 0.01
stepsize_callback = StepsizeCallback(cfl = cfl)
save_count = parse(Int, get(ENV, "ROTATING_SAVE_COUNT", "401"))
save_times = range(first(tspan), last(tspan); length = save_count)

ledger_integrals = zeros(6)
ledger_previous_terms = zeros(6)
ledger_previous_time = Ref(first(tspan))

function initialize_ledger_callback(callback, state, t, integrator)
    ledger_previous_terms .= rotating_ledger_terms(state, t)
    ledger_previous_time[] = t
    return nothing
end

ledger_condition(state, t, integrator) = true

function update_ledger!(integrator)
    current_terms = rotating_ledger_terms(integrator.u, integrator.t)
    step = integrator.t - ledger_previous_time[]
    ledger_integrals .+= 0.5 * step .* (ledger_previous_terms .+
                                        current_terms)
    ledger_previous_terms .= current_terms
    ledger_previous_time[] = integrator.t
    return nothing
end

ledger_callback = DiscreteCallback(ledger_condition, update_ledger!;
                                   initialize = initialize_ledger_callback,
                                   save_positions = (false, false))
callbacks = record_online_ledger ?
            CallbackSet(stepsize_callback, ledger_callback) :
            stepsize_callback

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0,
            adaptive = false,
            callback = callbacks,
            saveat = save_times,
            ode_default_options()...)

@assert all(isfinite, sol.u[end])
energy_history = [Trixi.integrate(entropy, state, semi; normalize = false)
                  for state in sol.u]
steady_energy = beam_quadrature_sum() do i, element
    entropy(steady_rotating_solution(node_coordinates[1, i, element]),
            equations_hyperbolic)
end
ledger_cumulative_residual = if record_online_ledger
    energy_history[end] - energy_history[1] +
    sum(ledger_integrals[1:4]) - sum(ledger_integrals[5:6])
else
    NaN
end
