if !isdefined(@__MODULE__, :BeamRunHelpers)
    Base.include(@__MODULE__, joinpath(@__DIR__, "beam_run_helpers.jl"))
end
using .BeamRunHelpers

using Printf
using Serialization: deserialize

example_directory = @__DIR__
checkpoint_path = isempty(ARGS) ?
                  joinpath(example_directory, "results",
                           "farokhi_05g_k3n4_refined_up_sweep_" *
                           "checkpoint_point_007.jls") :
                  abspath(ARGS[1])
output_directory = length(ARGS) < 2 ?
                   joinpath(example_directory, "results",
                            "base_excited_cantilever_extreme_cycle") :
                   abspath(ARGS[2])
frames_per_cycle = parse(Int,
                         get(ENV, "CANTILEVER_CYCLE_FRAMES", "96"))
points_per_element = parse(Int,
                           get(ENV, "CANTILEVER_GEOMETRY_POINTS_PER_ELEMENT",
                               "64"))
frames_per_cycle >= 8 ||
    throw(ArgumentError("CANTILEVER_CYCLE_FRAMES must be at least 8"))
points_per_element >= 4 ||
    throw(ArgumentError("CANTILEVER_GEOMETRY_POINTS_PER_ELEMENT must be at least 4"))

checkpoint = deserialize(checkpoint_path)
checkpoint.format_version in (4, 5) ||
    error("expected a format-4 or format-5 cantilever checkpoint")
length(checkpoint.rows) >= 1 ||
    error("checkpoint does not contain an accepted frequency")
accepted_row = last(checkpoint.rows)

for (key, value) in cantilever_checkpoint_environment(checkpoint, accepted_row)
    ENV[key] = value
end
ENV["CANTILEVER_RAMP_CYCLES"] = "0"
ENV["CANTILEVER_SETUP_ONLY"] = "true"

Base.include(@__MODULE__,
             joinpath(example_directory,
                      "elixir_base_excited_cantilever.jl"))

root_velocity.omega = accepted_row.omega
root_velocity.ramp_duration = 0.0
cycle_period = 2.0 * pi / root_velocity.omega
save_times_cycle = range(0.0, cycle_period;
                         length = frames_per_cycle + 1)
cycle_problem = remake(problem;
                       u0 = copy(checkpoint.continuation_state),
                       tspan = (0.0, cycle_period))
cycle_solution = solve(cycle_problem, algorithm;
                       reltol = relative_tolerance,
                       abstol = absolute_tolerance,
                       dtmax = maximum_step,
                       saveat = save_times_cycle,
                       save_start = true,
                       save_everystep = false,
                       maxiters = 10^7)
require_complete_solution(cycle_solution, last(cycle_problem.tspan))

function barycentric_weights(nodes)
    weights = ones(length(nodes))
    for i in eachindex(nodes), j in eachindex(nodes)
        if i != j
            weights[i] /= nodes[i] - nodes[j]
        end
    end
    return weights
end

function polynomial_value(nodes, weights, values, coordinate)
    for i in eachindex(nodes)
        if isapprox(coordinate, nodes[i]; atol = 8 * eps(Float64),
                    rtol = 0.0)
            return values[i]
        end
    end
    numerator = 0.0
    denominator = 0.0
    for i in eachindex(nodes)
        term = weights[i] / (coordinate - nodes[i])
        numerator += term * values[i]
        denominator += term
    end
    return numerator / denominator
end

barycentric = barycentric_weights(lobatto_nodes)

function centerline_geometry(state)
    state_array = reshape(state, 12, nnodes, nelements)
    material_coordinates = Float64[]
    vertical_positions = Float64[]
    transverse_positions = Float64[]
    angles = Float64[]
    vertical = 0.0
    transverse = 0.0
    angle = 0.0

    for (ordered_index, element) in enumerate(element_order)
        jacobian = volume_jacobians[element]
        center = sum(node_coordinates[1, :, element]) / nnodes
        f1_values = @view state_array[7, :, element]
        f2_values = @view state_array[8, :, element]
        m3_values = @view state_array[12, :, element]

        function kinematics(coordinate, local_angle)
            f1 = polynomial_value(lobatto_nodes, barycentric, f1_values,
                                  coordinate)
            f2 = polynomial_value(lobatto_nodes, barycentric, f2_values,
                                  coordinate)
            m3 = polynomial_value(lobatto_nodes, barycentric, m3_values,
                                  coordinate)
            gamma1 = 1.0 + flexibility_matrix[1, 1] * f1
            gamma2 = flexibility_matrix[2, 2] * f2
            return (jacobian *
                    (cos(local_angle) * gamma1 -
                     sin(local_angle) * gamma2),
                    jacobian *
                    (sin(local_angle) * gamma1 +
                     cos(local_angle) * gamma2),
                    jacobian * flexibility_matrix[6, 6] * m3)
        end

        coordinates = range(-1.0, 1.0;
                            length = points_per_element + 1)
        if ordered_index == 1
            push!(material_coordinates, center - jacobian)
            push!(vertical_positions, vertical)
            push!(transverse_positions, transverse)
            push!(angles, angle)
        end

        for index in 1:points_per_element
            left = coordinates[index]
            right = coordinates[index + 1]
            step = right - left
            midpoint = 0.5 * (left + right)

            k1 = kinematics(left, angle)
            k2 = kinematics(midpoint, angle + 0.5 * step * k1[3])
            k3 = kinematics(midpoint, angle + 0.5 * step * k2[3])
            k4 = kinematics(right, angle + step * k3[3])
            vertical += step / 6.0 *
                        (k1[1] + 2.0 * k2[1] +
                         2.0 * k3[1] + k4[1])
            transverse += step / 6.0 *
                          (k1[2] + 2.0 * k2[2] +
                           2.0 * k3[2] + k4[2])
            angle += step / 6.0 *
                     (k1[3] + 2.0 * k2[3] +
                      2.0 * k3[3] + k4[3])

            push!(material_coordinates, center + jacobian * right)
            push!(vertical_positions, vertical)
            push!(transverse_positions, transverse)
            push!(angles, angle)
        end
    end

    return (material_coordinates, vertical_positions,
            transverse_positions, angles)
end

mkpath(output_directory)
geometry_path = joinpath(output_directory, "centerline_cycle.csv")
tip_path = joinpath(output_directory, "tip_cycle.csv")
ledger_path = joinpath(output_directory, "cycle_ledger.csv")

open(geometry_path, "w") do io
    println(io, "frame,phase,material_coordinate,vertical,transverse,angle")
    for (frame, (state, time)) in enumerate(zip(cycle_solution.u, cycle_solution.t))
        phase = time / cycle_period
        coordinates, vertical, transverse, angles = centerline_geometry(state)
        for index in eachindex(coordinates)
            @printf(io, "%d,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    frame-1, phase, coordinates[index],
                    vertical[index], transverse[index], angles[index])
        end
    end
end

tips = [reconstruct_tip(state) for state in cycle_solution.u]
open(tip_path, "w") do io
    println(io, "frame,phase,longitudinal,transverse,rotation")
    for (frame, (tip, time)) in enumerate(zip(tips, cycle_solution.t))
        @printf(io, "%d,%.17g,%.17g,%.17g,%.17g\n",
                frame-1, time/cycle_period,
                tip[1], tip[2], tip[3])
    end
end

ledger = cantilever_cycle_ledger(cycle_solution.u, cycle_solution.t)
net_boundary_dissipation = ledger.left_boundary_dissipation +
                           ledger.right_boundary_dissipation -
                           ledger.sat_data_work
compact_residual = ledger.total_energy_change +
                   ledger.material_dissipation +
                   ledger.jump_dissipation +
                   net_boundary_dissipation -
                   ledger.physical_root_work
open(ledger_path, "w") do io
    println(io,
            "normalized_frequency,periodicity_error,total_energy_change," *
            "physical_root_work,material_dissipation,jump_dissipation," *
            "net_boundary_dissipation,compact_ledger_residual," *
            "compact_relative_ledger_residual")
    periodicity_error = maximum(abs.(tips[end] - tips[1]))
    @printf(io, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            accepted_row.normalized_frequency, periodicity_error,
            ledger.total_energy_change, ledger.physical_root_work,
            ledger.material_dissipation, ledger.jump_dissipation,
            net_boundary_dissipation, compact_residual,
            abs(compact_residual)/abs(ledger.physical_root_work))
end

println("Wrote ", geometry_path)
println("Wrote ", tip_path)
println("Wrote ", ledger_path)
