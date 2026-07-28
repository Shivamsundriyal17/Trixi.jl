using LinearAlgebra: norm
using Serialization: serialize

Base.include(@__MODULE__, joinpath(@__DIR__, "elixir_rotating_beam.jl"))

output_file = haskey(ENV, "ROTATING_OUTPUT_FILE") ?
              abspath(ENV["ROTATING_OUTPUT_FILE"]) :
              isempty(ARGS) ?
              joinpath(@__DIR__, "results", "rotating_beam.jls") :
              abspath(ARGS[1])
mkpath(dirname(output_file))

function rotating_field_metrics(state)
    state_array = reshape(state, 12, length(quadrature_weights),
                          length(volume_jacobians))
    f1_error_squared = 0.0
    v2_error_squared = 0.0
    f1_reference_squared = 0.0
    v2_reference_squared = 0.0
    f1_error_maximum = 0.0
    v2_error_maximum = 0.0
    f1_reference_maximum = 0.0
    v2_reference_maximum = 0.0

    for element in axes(node_coordinates, 3), i in eachindex(quadrature_weights)
        x = node_coordinates[1, i, element]
        reference = steady_rotating_solution(x)
        weight = volume_jacobians[element] * quadrature_weights[i]
        f1_error = state_array[7, i, element] - reference[7]
        v2_error = state_array[2, i, element] - reference[2]
        f1_error_squared += weight * f1_error^2
        v2_error_squared += weight * v2_error^2
        f1_reference_squared += weight * reference[7]^2
        v2_reference_squared += weight * reference[2]^2
        f1_error_maximum = max(f1_error_maximum, abs(f1_error))
        v2_error_maximum = max(v2_error_maximum, abs(v2_error))
        f1_reference_maximum = max(f1_reference_maximum, abs(reference[7]))
        v2_reference_maximum = max(v2_reference_maximum, abs(reference[2]))
    end

    return (f1_l2 = sqrt(f1_error_squared),
            f1_relative_l2 = sqrt(f1_error_squared /
                                  f1_reference_squared),
            f1_linf = f1_error_maximum,
            f1_relative_linf = f1_error_maximum /
                               f1_reference_maximum,
            v2_l2 = sqrt(v2_error_squared),
            v2_relative_l2 = sqrt(v2_error_squared /
                                  v2_reference_squared),
            v2_linf = v2_error_maximum,
            v2_relative_linf = v2_error_maximum /
                               v2_reference_maximum)
end

function instantaneous_ledger(state, t)
    hyperbolic_rhs = similar(state)
    parabolic_rhs = similar(state)
    Trixi.rhs!(hyperbolic_rhs, state, semi, t)
    Trixi.rhs_parabolic!(parabolic_rhs, state, semi, t)
    state_array = reshape(state, 12, length(quadrature_weights),
                          length(volume_jacobians))
    rhs_array = reshape(hyperbolic_rhs + parabolic_rhs, 12,
                        length(quadrature_weights),
                        length(volume_jacobians))
    energy_rate = beam_quadrature_sum() do i, element
        state_node = SVector{12}(state_array[:, i, element])
        rhs_node = SVector{12}(rhs_array[:, i, element])
        dot(equations_hyperbolic.capacity_matrix * state_node,
            rhs_node)
    end
    terms = rotating_ledger_terms(state, t)
    residual = energy_rate + sum(terms[1:4]) - sum(terms[5:6])
    scale = max(abs(energy_rate), maximum(abs, terms), 1.0)
    return (energy_rate,
            material_dissipation = terms[1],
            jump_dissipation = terms[2],
            left_boundary_dissipation = terms[3],
            right_boundary_dissipation = terms[4],
            physical_root_power = terms[5],
            sat_data_power = terms[6],
            residual,
            relative_residual = abs(residual) / scale)
end

field_history = [rotating_field_metrics(state) for state in sol.u]
ledger_history = [instantaneous_ledger(state, time)
                  for (state, time) in zip(sol.u, sol.t)]
final_metrics = last(field_history)
late_start = max(first(tspan), last(tspan) - 5)
late_indices = findall(time -> time >= late_start, sol.t)
late_f1_relative_linf = maximum(field_history[index].f1_relative_linf
                                for index in late_indices)
late_v2_relative_linf = maximum(field_history[index].v2_relative_linf
                                for index in late_indices)
maximum_ledger_relative_residual = maximum(item.relative_residual
                                           for item in ledger_history)
ledger_integral(index) = record_online_ledger ? ledger_integrals[index] : NaN

final_state = reshape(sol.u[end], 12, length(quadrature_weights),
                      length(volume_jacobians))
final_beam_length = beam_quadrature_sum() do i, element
    strain_curvature = flexibility_matrix *
                       final_state[7:12, i, element]
    norm(SVector(1 + strain_curvature[1],
                 strain_curvature[2], strain_curvature[3]))
end
steady_beam_length = beam_quadrature_sum() do i, element
    reference = steady_rotating_solution(node_coordinates[1, i, element])
    strain_curvature = flexibility_matrix * reference[7:12]
    norm(SVector(1 + strain_curvature[1],
                 strain_curvature[2], strain_curvature[3]))
end

diagnostics = (final_metrics,
               energy_final = energy_history[end],
               energy_steady = steady_energy,
               energy_relative_error = abs(energy_history[end] -
                                           steady_energy) / steady_energy,
               final_beam_length,
               steady_beam_length,
               beam_length_relative_error = abs(final_beam_length -
                                                steady_beam_length) /
                                            steady_beam_length,
               late_interval = (late_start, last(tspan)),
               late_f1_relative_linf,
               late_v2_relative_linf,
               maximum_ledger_relative_residual,
               cumulative_material_dissipation = ledger_integral(1),
               cumulative_jump_dissipation = ledger_integral(2),
               cumulative_left_boundary_dissipation = ledger_integral(3),
               cumulative_right_boundary_dissipation = ledger_integral(4),
               cumulative_physical_root_work = ledger_integral(5),
               cumulative_sat_data_work = ledger_integral(6),
               cumulative_ledger_residual = ledger_cumulative_residual)

node_coordinates_flat = vec(copy(semi.cache.elements.node_coordinates))
state_shape = (12, size(semi.cache.elements.node_coordinates, 2),
               size(semi.cache.elements.node_coordinates, 3))
states = [Array(reshape(state, state_shape)) for state in sol.u]
steady_states = [collect(steady_rotating_solution(x))
                 for x in node_coordinates_flat]

data = (times = collect(sol.t),
        node_coordinates = node_coordinates_flat,
        states,
        flexibility_matrix = Matrix(flexibility_matrix),
        initial_curvature = collect(equations_hyperbolic.initial_curvature),
        energy_history,
        f1_relative_linf_history = getproperty.(field_history,
                                                :f1_relative_linf),
        v2_relative_linf_history = getproperty.(field_history,
                                                :v2_relative_linf),
        terminal_angular_speed,
        ramp_duration,
        steady_states,
        case_name,
        damping_multiplier,
        steady_initial_condition,
        polydeg,
        ncells = 2^initial_refinement_level,
        cfl,
        diagnostics)
serialize(output_file, data)
println("wrote ", output_file)
