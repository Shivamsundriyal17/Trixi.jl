if !isdefined(@__MODULE__, :BeamRunHelpers)
    Base.include(@__MODULE__, joinpath(@__DIR__, "beam_run_helpers.jl"))
end
using .BeamRunHelpers

using LinearAlgebra: dot
using Printf
using Serialization: deserialize

example_directory = @__DIR__
checkpoint_path = isempty(ARGS) ?
                  joinpath(example_directory, "results",
                           "farokhi_05g_k3n4_refined_up_sweep_" *
                           "checkpoint_point_007.jls") :
                  abspath(ARGS[1])
output_path = length(ARGS) < 2 ?
              joinpath(example_directory, "results",
                       "base_excited_cantilever_gravity_work_audit.csv") :
              abspath(ARGS[2])

checkpoint = deserialize(checkpoint_path)
checkpoint.format_version in (3, 4, 5) ||
    error("expected a format-3, format-4, or format-5 cantilever checkpoint")
isempty(checkpoint.rows) && error("checkpoint has no accepted response point")
accepted_row = last(checkpoint.rows)

ENV["CANTILEVER_ACCELERATION_RMS_G"] = string(checkpoint.acceleration_rms_g)
ENV["CANTILEVER_FREQUENCY_HZ"] = string(accepted_row.frequency_hz)
ENV["CANTILEVER_POLYDEG"] = string(checkpoint.polydeg)
ENV["CANTILEVER_REFINEMENT_LEVEL"] = string(checkpoint.refinement_level)
ENV["CANTILEVER_CONSTRAINT_MULTIPLIER"] = string(checkpoint.constraint_multiplier)
checkpoint_gravity(checkpoint) == 1.0 ||
    error("this gravity work audit requires a unit-gravity checkpoint")
ENV["CANTILEVER_GRAVITY_MULTIPLIER"] = string(checkpoint_gravity(checkpoint))
ENV["CANTILEVER_RAMP_CYCLES"] = "0"
ENV["CANTILEVER_SETUP_ONLY"] = "true"
ENV["CANTILEVER_RELTOL"] = "1e-5"
ENV["CANTILEVER_ABSTOL"] = "1e-7"

Base.include(@__MODULE__,
             joinpath(example_directory,
                      "elixir_base_excited_cantilever.jl"))

# Release a large-deformation state from the periodic upper branch with the
# root fixed. The resulting motion is a genuine damped transient, not a
# cycle-integrated balance in which the potential change is nearly zero.
root_velocity.acceleration = 0.0
root_velocity.ramp_duration = 0.0
root_velocity.omega = accepted_row.omega
local_period = 2.0 * pi / root_velocity.omega
audit_end = 0.0625 * local_period
audit_times = range(0.0, audit_end; length = 5)
audit_problem = remake(problem;
                       u0 = copy(checkpoint.continuation_state),
                       tspan = (0.0, audit_end))
audit_solution = solve(audit_problem, algorithm;
                       reltol = relative_tolerance,
                       abstol = absolute_tolerance,
                       dtmax = local_period / 64.0,
                       saveat = audit_times,
                       save_start = true,
                       save_everystep = false,
                       maxiters = 10^7)
require_complete_solution(audit_solution, last(audit_problem.tspan))

function gravity_power(state)
    gravity_derivative = zeros(eltype(state), length(state))
    add_dead_gravity!(gravity_derivative, state, rhs_parameters)
    state_array = reshape(state, 12, nnodes, nelements)
    derivative_array = reshape(gravity_derivative, 12, nnodes, nelements)
    return beam_quadrature_sum() do i, element
        velocity = state_array[1:6, i, element]
        gravity_acceleration = derivative_array[1:6, i, element]
        dot(velocity, mass_matrix * gravity_acceleration)
    end
end

# Differentiate the discrete reconstruction itself. The angle rate is the
# antiderivative of the curvature rate, while the vertical centreline velocity
# follows by differentiating R(angle)*(e1 + C*f) at every Lobatto node.
function potential_directional_derivative(state, state_derivative)
    state_array = reshape(state, 12, nnodes, nelements)
    derivative_array = reshape(state_derivative, 12, nnodes, nelements)
    angles = similar(angle_cache)
    angle_rates = similar(angle_cache)
    reconstruct_angle!(angles, state, rhs_parameters)
    reconstruct_angle!(angle_rates, state_derivative, rhs_parameters)
    vertical_velocity = zeros(nnodes, nelements)
    left_vertical_velocity = 0.0

    for element in element_order
        jacobian = volume_jacobians[element]
        vertical_derivative_rate = zeros(nnodes)
        for i in 1:nnodes
            force1 = state_array[7, i, element]
            force2 = state_array[8, i, element]
            force1_rate = derivative_array[7, i, element]
            force2_rate = derivative_array[8, i, element]
            gamma1 = 1.0 + flexibility_matrix[1, 1] * force1
            gamma2 = flexibility_matrix[2, 2] * force2
            gamma1_rate = flexibility_matrix[1, 1] * force1_rate
            gamma2_rate = flexibility_matrix[2, 2] * force2_rate
            angle = angles[i, element]
            angle_rate = angle_rates[i, element]
            vertical_derivative_rate[i] = cos(angle) * gamma1_rate -
                                          sin(angle) * gamma2_rate -
                                          angle_rate *
                                          (sin(angle) * gamma1 + cos(angle) * gamma2)
        end
        vertical_velocity[:, element] .= left_vertical_velocity .+
                                         jacobian .*
                                         (integration_matrix * vertical_derivative_rate)
        left_vertical_velocity += jacobian * dot(solver.basis.weights,
                                      vertical_derivative_rate)
    end

    return cantilever_gamma * beam_quadrature_sum() do i, element
        vertical_velocity[i, element]
    end
end

mkpath(dirname(output_path))
maximum_absolute_residual = Ref(0.0)
maximum_relative_residual = Ref(0.0)
open(output_path, "w") do io
    println(io,
            "time,forcing_period_fraction,gravitational_potential," *
            "potential_time_derivative,gravity_power,residual," *
            "relative_residual")
    for (time, state) in zip(audit_solution.t, audit_solution.u)
        state_derivative = similar(state)
        cantilever_rhs!(state_derivative, state, rhs_parameters, time)
        potential_rate = potential_directional_derivative(state, state_derivative)
        power = gravity_power(state)
        residual = potential_rate + power
        scale = max(abs(potential_rate), abs(power), eps(Float64))
        relative_residual = abs(residual) / scale
        maximum_absolute_residual[] = max(maximum_absolute_residual[],
                                          abs(residual))
        maximum_relative_residual[] = max(maximum_relative_residual[],
                                          relative_residual)
        @printf(io, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                time, time/local_period, gravitational_potential(state),
                potential_rate, power, residual, relative_residual)
    end
end

@printf("Wrote %s\n", output_path)
@printf("maximum absolute gravity-work residual = %.6e\n",
        maximum_absolute_residual[])
@printf("maximum relative gravity-work residual = %.6e\n",
        maximum_relative_residual[])
maximum_relative_residual[] <= 1.0e-3 ||
    error("gravity work--potential audit exceeds the 1e-3 tolerance")
