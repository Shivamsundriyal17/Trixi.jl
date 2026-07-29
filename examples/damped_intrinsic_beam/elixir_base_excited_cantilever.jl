using LinearAlgebra: Diagonal, dot
using OrdinaryDiffEqRosenbrock
using Printf
using SparseArrays
using Trixi

# Experimentally validated in-vacuo cantilever of
# Farokhi, Xia, and Erturk, Nonlinear Dynamics 107 (2022), 457--475.
#
# The calculation is nondimensionalized with
#   x = s / L,  tau = t / T,  T = L^2 * sqrt(rho*A / (E*I)).
# Thus, the beam length, transverse mass per unit length, and weak-axis
# bending stiffness are all one. The prescribed root velocity represents the
# measured base motion in an inertial frame. Gravity is added below after
# reconstructing the material angle from the intrinsic curvature.

const FAROKHI_E = 200.0e9
const FAROKHI_RHO = 7800.0
const FAROKHI_H = 0.0762e-3
const FAROKHI_B = 9.0e-3
const FAROKHI_L = 81.5e-3
const FAROKHI_G = 9.81
const FAROKHI_NU = 0.30
const FAROKHI_SHEAR_CORRECTION = 5.0 / 6.0
const FAROKHI_REFERENCE_OMEGA1 = 1.875104068711961^2

const FAROKHI_A = FAROKHI_B * FAROKHI_H
const FAROKHI_I_WEAK = FAROKHI_B * FAROKHI_H^3 / 12.0
const FAROKHI_I_STRONG = FAROKHI_H * FAROKHI_B^3 / 12.0
const FAROKHI_T = FAROKHI_L^2 *
                   sqrt(FAROKHI_RHO * FAROKHI_A /
                        (FAROKHI_E * FAROKHI_I_WEAK))
const FAROKHI_GAMMA = FAROKHI_RHO * FAROKHI_A * FAROKHI_G *
                       FAROKHI_L^3 / (FAROKHI_E * FAROKHI_I_WEAK)
const FAROKHI_CHI_WEAK = FAROKHI_I_WEAK /
                          (FAROKHI_A * FAROKHI_L^2)
const FAROKHI_CHI_STRONG = FAROKHI_I_STRONG /
                            (FAROKHI_A * FAROKHI_L^2)
const FAROKHI_CHI_POLAR = FAROKHI_CHI_WEAK + FAROKHI_CHI_STRONG
const FAROKHI_AXIAL_STIFFNESS = FAROKHI_A * FAROKHI_L^2 /
                                 FAROKHI_I_WEAK
const FAROKHI_SHEAR_STIFFNESS = FAROKHI_SHEAR_CORRECTION /
                                 (2.0 * (1.0 + FAROKHI_NU)) *
                                 FAROKHI_AXIAL_STIFFNESS
const FAROKHI_STRONG_BENDING_STIFFNESS = FAROKHI_I_STRONG /
                                          FAROKHI_I_WEAK
const FAROKHI_ETA_D = 0.0037

acceleration_rms_g = parse(Float64,
                           get(ENV, "CANTILEVER_ACCELERATION_RMS_G", "0.2"))
frequency_hz = parse(Float64,
                     get(ENV, "CANTILEVER_FREQUENCY_HZ", "9.2"))
omega = 2.0 * pi * frequency_hz * FAROKHI_T
acceleration = sqrt(2.0) * acceleration_rms_g * FAROKHI_G *
               FAROKHI_T^2 / FAROKHI_L

# The physical ratios are the default. The multiplier supports a transparent
# slender-limit sensitivity study without changing the weak-axis bending
# stiffness. A value below one softens only axial/shear constraints.
constraint_multiplier = parse(Float64,
                              get(ENV, "CANTILEVER_CONSTRAINT_MULTIPLIER",
                                  "1.0"))
axial_stiffness = constraint_multiplier * FAROKHI_AXIAL_STIFFNESS
shear_stiffness = constraint_multiplier * FAROKHI_SHEAR_STIFFNESS

# Components are [v1,v2,v3,omega1,omega2,omega3] and
# [f1,f2,f3,m1,m2,m3]. The experiment is invariant in the (1,2) plane.
mass_matrix = Diagonal([1.0, 1.0, 1.0,
                        FAROKHI_CHI_POLAR, FAROKHI_CHI_STRONG,
                        FAROKHI_CHI_WEAK])
torsional_stiffness = 4.0 / (2.0 * (1.0 + FAROKHI_NU))
stiffnesses = [axial_stiffness, shear_stiffness, shear_stiffness,
               torsional_stiffness, FAROKHI_STRONG_BENDING_STIFFNESS, 1.0]
flexibility_matrix = Diagonal(1.0 ./ stiffnesses)

# Farokhi et al. apply Kelvin--Voigt damping only to weak-axis bending.
# Since r_tau = C_tau*C^{-1}*kappa_dot in the Trixi equation and the
# nondimensional weak-axis bending stiffness is one, C_tau[6,6] = eta_d.
damping_matrix = Diagonal([0.0, 0.0, 0.0, 0.0, 0.0, FAROKHI_ETA_D])

function initial_condition(x, t, equations::DampedIntrinsicBeamEquations1D)
    # Exact straight equilibrium under the initially vertical dead load:
    # f1_x - gamma = 0, f1(1) = 0.
    f1 = -FAROKHI_GAMMA * (1.0 - x[1])
    return SVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   f1, 0.0, 0.0, 0.0, 0.0, 0.0)
end

ramp_cycles = parse(Float64, get(ENV, "CANTILEVER_RAMP_CYCLES", "8"))
period = 2.0 * pi / omega
ramp_duration = ramp_cycles * period

@inline function smooth_ramp(t, ramp_duration)
    if t <= 0.0
        return 0.0
    elseif t >= ramp_duration
        return 1.0
    else
        return sinpi(0.5 * t / ramp_duration)^2
    end
end

mutable struct HarmonicBaseVelocity
    omega::Float64
    acceleration::Float64
    ramp_duration::Float64
end

function (base::HarmonicBaseVelocity)(x, t, equations)
    # Harmonic phase is chosen so that the computation starts from rest.
    # The full-amplitude acceleration has peak value `acceleration`.
    v2 = smooth_ramp(t, base.ramp_duration) *
         base.acceleration / base.omega * sin(base.omega * t)
    return SVector(0.0, v2, 0.0, 0.0, 0.0, 0.0)
end

root_velocity = HarmonicBaseVelocity(omega, acceleration, ramp_duration)

equations_hyperbolic = DampedIntrinsicBeamEquations1D(;
                                                      mass_matrix,
                                                      flexibility_matrix,
                                                      damping_matrix,
                                                      initial_curvature = zeros(3))
equations_parabolic = DampedIntrinsicBeamDiffusion1D(equations_hyperbolic)
boundary_condition = BoundaryConditionDampedIntrinsicBeam(;
                                                          left_velocity = root_velocity)
boundary_conditions = (boundary_condition, boundary_condition)

polydeg = parse(Int, get(ENV, "CANTILEVER_POLYDEG", "4"))
refinement_level = parse(Int,
                         get(ENV, "CANTILEVER_REFINEMENT_LEVEL", "1"))
solver = DGSEM(polydeg = polydeg, surface_flux = flux_upwind)
solver_parabolic = ViscousFormulationLocalDG()
mesh = TreeMesh((0.0,), (1.0,);
                initial_refinement_level = refinement_level,
                n_cells_max = 10_000,
                periodicity = false)
semi = SemidiscretizationHyperbolicParabolic(mesh,
                                             (equations_hyperbolic,
                                              equations_parabolic),
                                             initial_condition, solver;
                                             solver_parabolic,
                                             boundary_conditions)

# Polynomial antiderivative matrix on the reference element. If `values`
# contains a degree-k nodal polynomial, Q*values contains its integrals from
# -1 to each Lobatto node. This is used for the nonlocal angle reconstruction.
lobatto_nodes = collect(solver.basis.nodes)
nnodes = length(lobatto_nodes)
vandermonde = [lobatto_nodes[i]^(power - 1)
               for i in 1:nnodes, power in 1:nnodes]
integrated_monomials = [(lobatto_nodes[i]^power - (-1.0)^power) / power
                        for i in 1:nnodes, power in 1:nnodes]
integration_matrix = integrated_monomials / vandermonde

node_coordinates = semi.cache.elements.node_coordinates
nelements = size(node_coordinates, 3)
element_order = sortperm([sum(node_coordinates[1, :, element])
                          for element in 1:nelements])
volume_jacobians = [abs(inv(semi.cache.elements.inverse_jacobian[element]))
                    for element in 1:nelements]
angle_cache = zeros(nnodes, nelements)

struct CantileverRHSParameters{Semi, VectorT, MatrixT, ElementOrder,
                               Jacobians, Weights}
    semi::Semi
    parabolic_derivative::VectorT
    angle_cache::MatrixT
    integration_matrix::MatrixT
    element_order::ElementOrder
    volume_jacobians::Jacobians
    quadrature_weights::Weights
    nnodes::Int
    nelements::Int
    bending_flexibility::Float64
    gamma::Float64
end

function reconstruct_angle!(angles, state, parameters::CantileverRHSParameters)
    state_array = reshape(state, 12, parameters.nnodes,
                          parameters.nelements)
    fill!(angles, 0.0)
    left_angle = 0.0
    for element in parameters.element_order
        jacobian = parameters.volume_jacobians[element]
        for i in 1:parameters.nnodes
            integral = 0.0
            for j in 1:parameters.nnodes
                # Total curvature is kappa_0 + C_m*m. Here kappa_0=0.
                integral += parameters.integration_matrix[i, j] *
                            parameters.bending_flexibility *
                            state_array[12, j, element]
            end
            angles[i, element] = left_angle + jacobian * integral
        end
        full_integral = 0.0
        for j in 1:parameters.nnodes
            full_integral += parameters.quadrature_weights[j] *
                             parameters.bending_flexibility *
                             state_array[12, j, element]
        end
        left_angle += jacobian * full_integral
    end
    return angles
end

function add_dead_gravity!(derivative, state,
                           parameters::CantileverRHSParameters)
    reconstruct_angle!(parameters.angle_cache, state, parameters)
    derivative_array = reshape(derivative, 12, parameters.nnodes,
                               parameters.nelements)
    for element in 1:parameters.nelements, i in 1:parameters.nnodes
        angle = parameters.angle_cache[i, element]
        # R(angle)' * (-gamma*E1) in the material frame.
        derivative_array[1, i, element] -= parameters.gamma * cos(angle)
        derivative_array[2, i, element] += parameters.gamma * sin(angle)
    end
    return derivative
end

split_problem = semidiscretize(semi, (0.0, period))
parabolic_derivative = similar(split_problem.u0)
rhs_parameters = CantileverRHSParameters(semi, parabolic_derivative,
                                         angle_cache, integration_matrix,
                                         element_order, volume_jacobians,
                                         collect(solver.basis.weights),
                                         nnodes, nelements,
                                         flexibility_matrix[6, 6],
                                         FAROKHI_GAMMA)

function cantilever_rhs!(derivative, state,
                         parameters::CantileverRHSParameters, t)
    Trixi.rhs!(derivative, state, parameters.semi, t)
    Trixi.rhs_parabolic!(parameters.parabolic_derivative, state,
                         parameters.semi, t)
    derivative .+= parameters.parabolic_derivative
    add_dead_gravity!(derivative, state, parameters)
    return nothing
end

settling_cycles = parse(Int, get(ENV, "CANTILEVER_SETTLING_CYCLES", "50"))
measurement_cycles = parse(Int,
                           get(ENV, "CANTILEVER_MEASUREMENT_CYCLES", "2"))
samples_per_cycle = parse(Int,
                          get(ENV, "CANTILEVER_SAMPLES_PER_CYCLE", "96"))
total_cycles = ceil(Int, ramp_cycles) + settling_cycles + measurement_cycles
t_end = total_cycles * period
measurement_start = t_end - measurement_cycles * period
save_times = range(measurement_start, t_end;
                   length = measurement_cycles * samples_per_cycle + 1)
function numerical_jacobian_sparsity(rhs!, initial_state, parameters)
    state = copy(initial_state)
    # Activate all algebraic couplings that vanish at the straight equilibrium.
    for index in eachindex(state)
        state[index] += 0.03 * sin(0.731 * index)
    end
    base_derivative = similar(state)
    perturbed_derivative = similar(state)
    rhs!(base_derivative, state, parameters, 0.37 * period)

    rows = Int[]
    columns = Int[]
    for column in eachindex(state)
        original_value = state[column]
        step = sqrt(eps(Float64)) * max(1.0, abs(original_value))
        state[column] = original_value + step
        rhs!(perturbed_derivative, state, parameters, 0.37 * period)
        state[column] = original_value

        column_scale = maximum(abs(perturbed_derivative[row] -
                                   base_derivative[row])
                               for row in eachindex(state))
        threshold = max(1.0e-11 * column_scale, 1.0e-13)
        for row in eachindex(state)
            difference = abs(perturbed_derivative[row] -
                             base_derivative[row])
            if difference > threshold || row == column
                push!(rows, row)
                push!(columns, column)
            end
        end
    end
    return sparse(rows, columns, ones(length(rows)),
                  length(state), length(state))
end

use_sparse_jacobian = lowercase(get(ENV, "CANTILEVER_SPARSE_JACOBIAN",
                                    "true")) in ("1", "true", "yes")
jacobian_prototype = use_sparse_jacobian ?
                     numerical_jacobian_sparsity(cantilever_rhs!,
                                                 split_problem.u0,
                                                 rhs_parameters) :
                     nothing
ode_function = ODEFunction(cantilever_rhs!;
                           jac_prototype = jacobian_prototype)
problem = ODEProblem(ode_function, copy(split_problem.u0), (0.0, t_end),
                     rhs_parameters)

relative_tolerance = parse(Float64,
                           get(ENV, "CANTILEVER_RELTOL", "1e-4"))
absolute_tolerance = parse(Float64,
                           get(ENV, "CANTILEVER_ABSTOL", "1e-6"))
maximum_step = parse(Float64,
                     get(ENV, "CANTILEVER_DTMAX",
                         string(period / 32.0)))

# Finite differences are used because the Trixi caches are specialized to
# Float64. The sparsity pattern is detected once at a generic perturbed state;
# sparse coloring then avoids one residual evaluation per degree of freedom.
linear_solver = use_sparse_jacobian ?
                OrdinaryDiffEqRosenbrock.LinearSolve.KLUFactorization() :
                nothing
algorithm = Rodas5P(;
                    autodiff = OrdinaryDiffEqRosenbrock.AutoFiniteDiff(),
                    linsolve = linear_solver)
solution = solve(problem, algorithm;
                 reltol = relative_tolerance,
                 abstol = absolute_tolerance,
                 dtmax = maximum_step,
                 saveat = save_times,
                 save_start = false,
                 save_everystep = false,
                 maxiters = 10^7)

@assert SciMLBase.successful_retcode(solution)
@assert all(state -> all(isfinite, state), solution.u)

function reconstruct_tip(state)
    state_array = reshape(state, 12, nnodes, nelements)
    reconstruct_angle!(angle_cache, state, rhs_parameters)
    vertical = 0.0
    transverse = 0.0
    for element in element_order
        jacobian = volume_jacobians[element]
        for i in 1:nnodes
            f1 = state_array[7, i, element]
            f2 = state_array[8, i, element]
            gamma1 = 1.0 + flexibility_matrix[1, 1] * f1
            gamma2 = flexibility_matrix[2, 2] * f2
            angle = angle_cache[i, element]
            weight = jacobian * solver.basis.weights[i]
            vertical += weight *
                        (cos(angle) * gamma1 - sin(angle) * gamma2)
            transverse += weight *
                          (sin(angle) * gamma1 + cos(angle) * gamma2)
        end
    end
    tip_angle = angle_cache[end, last(element_order)]
    return SVector(vertical - 1.0, transverse, tip_angle)
end

capacity_abs_flux = equations_hyperbolic.capacity_matrix *
                    equations_hyperbolic.propagation_matrix_abs

function beam_quadrature_sum(function_)
    value = 0.0
    for element in axes(node_coordinates, 3), i in eachindex(solver.basis.weights)
        value += volume_jacobians[element] * solver.basis.weights[i] *
                 function_(i, element)
    end
    return value
end

function gravitational_potential(state)
    state_array = reshape(state, 12, nnodes, nelements)
    reconstruct_angle!(angle_cache, state, rhs_parameters)
    vertical_positions = zeros(nnodes, nelements)
    left_vertical = 0.0

    for element in element_order
        jacobian = volume_jacobians[element]
        vertical_derivative = zeros(nnodes)
        for i in 1:nnodes
            f1 = state_array[7, i, element]
            f2 = state_array[8, i, element]
            gamma1 = 1.0 + flexibility_matrix[1, 1] * f1
            gamma2 = flexibility_matrix[2, 2] * f2
            angle = angle_cache[i, element]
            vertical_derivative[i] =
                cos(angle) * gamma1 - sin(angle) * gamma2
        end
        vertical_positions[:, element] .=
            left_vertical .+
            jacobian .* (integration_matrix * vertical_derivative)
        left_vertical += jacobian *
                         dot(solver.basis.weights, vertical_derivative)
    end

    return FAROKHI_GAMMA * beam_quadrature_sum() do i, element
        vertical_positions[i, element]
    end
end

function cantilever_total_energy(state)
    intrinsic_energy = Trixi.integrate(entropy, state, semi;
                                       normalize = false)
    return intrinsic_energy + gravitational_potential(state)
end

function cantilever_ledger_terms(state, t)
    # Only the auxiliary gradient is needed for the Kelvin--Voigt resultant.
    # Reusing Trixi's LDG gradient makes this diagnostic consistent with the
    # semidiscrete damping operator.
    state_parabolic = Trixi.wrap_array(state, mesh, equations_parabolic,
                                       solver, semi.cache_parabolic)
    viscous_container = semi.cache_parabolic.viscous_container
    Trixi.transform_variables!(viscous_container.u_transformed,
                               state_parabolic, mesh, equations_parabolic,
                               solver, solver_parabolic, semi.cache,
                               semi.cache_parabolic)
    Trixi.calc_gradient!(viscous_container.gradients,
                         viscous_container.u_transformed, t, mesh,
                         equations_parabolic,
                         semi.boundary_conditions_parabolic, solver,
                         solver_parabolic, semi.cache,
                         semi.cache_parabolic)

    state_array = reshape(state, 12, nnodes, nelements)
    gradients = viscous_container.gradients
    material_dissipation = beam_quadrature_sum() do i, element
        state_node = SVector{12}(state_array[:, i, element])
        gradient_node = SVector{12}(gradients[:, i, element])
        damping_resultant =
            intrinsic_beam_damping_resultant(state_node, gradient_node,
                                              equations_parabolic)
        damping_resultant[6]^2 / FAROKHI_ETA_D
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
    left_damping_resultant =
        intrinsic_beam_damping_resultant(left_state, left_gradient,
                                          equations_parabolic)
    prescribed_velocity = root_velocity(SVector(0.0), t,
                                        equations_hyperbolic)

    left_boundary_dissipation =
        dot(left_velocity,
            equations_hyperbolic.left_impedance * left_velocity)
    right_boundary_dissipation =
        dot(right_resultant,
            equations_hyperbolic.right_impedance * right_resultant)
    physical_root_power =
        -dot(prescribed_velocity, left_resultant + left_damping_resultant)
    sat_data_power =
        dot(prescribed_velocity,
            equations_hyperbolic.left_impedance * left_velocity)

    return SVector(material_dissipation, jump_dissipation,
                   left_boundary_dissipation, right_boundary_dissipation,
                   physical_root_power, sat_data_power)
end

function cantilever_cycle_ledger(states, times)
    length(states) == length(times) ||
        throw(DimensionMismatch("states and times must have equal lengths"))
    length(states) >= 2 ||
        throw(ArgumentError("at least two cycle samples are required"))

    integrated_terms = zeros(6)
    previous_terms = cantilever_ledger_terms(first(states), first(times))
    for index in 2:length(states)
        current_terms = cantilever_ledger_terms(states[index], times[index])
        step = times[index] - times[index - 1]
        integrated_terms .+=
            0.5 * step .* (previous_terms .+ current_terms)
        previous_terms = current_terms
    end

    energy_change = cantilever_total_energy(last(states)) -
                    cantilever_total_energy(first(states))
    residual = energy_change + sum(integrated_terms[1:4]) -
               sum(integrated_terms[5:6])
    scale = max(abs(energy_change),
                sum(abs, integrated_terms[1:4]),
                sum(abs, integrated_terms[5:6]),
                eps(Float64))

    return (;
            total_energy_change = energy_change,
            material_dissipation = integrated_terms[1],
            jump_dissipation = integrated_terms[2],
            left_boundary_dissipation = integrated_terms[3],
            right_boundary_dissipation = integrated_terms[4],
            physical_root_work = integrated_terms[5],
            sat_data_work = integrated_terms[6],
            ledger_residual = residual,
            relative_ledger_residual = abs(residual) / scale)
end

tip_history = [reconstruct_tip(state) for state in solution.u]
last_cycle = (length(tip_history) - samples_per_cycle):length(tip_history)
previous_cycle = ((length(tip_history) - 2 * samples_per_cycle):
                  (length(tip_history) - samples_per_cycle))
last_tip = tip_history[last_cycle]
previous_tip = tip_history[previous_cycle]

transverse_peak = maximum(abs(value[2]) for value in last_tip)
longitudinal_minimum = minimum(value[1] for value in last_tip)
rotation_peak = maximum(abs(value[3]) for value in last_tip)
periodicity_error = maximum(maximum(abs, last_tip[i] - previous_tip[i])
                            for i in eachindex(last_tip))

result = (;
          acceleration_rms_g,
          frequency_hz,
          omega,
          frequency_over_unloaded_omega1 =
              omega / FAROKHI_REFERENCE_OMEGA1,
          acceleration,
          constraint_multiplier,
          polydeg,
          refinement_level,
          cells = 2^refinement_level,
          dofs = length(split_problem.u0),
          jacobian_nonzeros = use_sparse_jacobian ?
                              nnz(jacobian_prototype) :
                              length(split_problem.u0)^2,
          transverse_peak,
          longitudinal_minimum,
          rotation_peak,
          periodicity_error,
          accepted_steps = solution.destats.naccept,
          rejected_steps = solution.destats.nreject,
          rhs_evaluations = solution.destats.nf,
          retcode = string(solution.retcode))

if abspath(PROGRAM_FILE) == @__FILE__
    @printf("Farokhi cantilever: %.1fg RMS, f = %.6f Hz, Omega = %.8f\n",
            acceleration_rms_g, frequency_hz, omega)
    @printf("  k = %d, cells = %d, dofs = %d, constraint multiplier = %.3e\n",
            polydeg, 2^refinement_level, length(split_problem.u0),
            constraint_multiplier)
    @printf("  max|w_tip| = %.10e, min(u_tip) = %.10e, max|psi_tip| = %.10e\n",
            transverse_peak, longitudinal_minimum, rotation_peak)
    @printf("  cycle mismatch = %.3e, accepted/rejected = %d/%d, retcode = %s\n",
            periodicity_error, solution.destats.naccept,
            solution.destats.nreject, string(solution.retcode))
end
