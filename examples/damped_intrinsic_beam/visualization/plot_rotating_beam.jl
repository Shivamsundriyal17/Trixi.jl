using Plots
using Serialization: deserialize

include(joinpath(@__DIR__, "beam_geometry.jl"))

input_file = isempty(ARGS) ?
             joinpath(@__DIR__, "..", "results", "rotating_beam.jls") :
             abspath(ARGS[1])
output_directory = length(ARGS) < 2 ?
                   joinpath(@__DIR__, "..", "results") :
                   abspath(ARGS[2])
mkpath(output_directory)
data = deserialize(input_file)

snapshot_indices = unique(round.(Int,
                                 range(1, length(data.times); length = 10)))
colors = palette(:viridis, length(snapshot_indices))
body_plot = plot(aspect_ratio = :equal, xlabel = "x1", ylabel = "x2",
                 title = "Rotating beam: root-attached frame", legend = false)
inertial_plot = plot(aspect_ratio = :equal, xlabel = "x1", ylabel = "x2",
                     title = "Rotating beam: inertial frame", legend = false)

for (color_index, state_index) in enumerate(snapshot_indices)
    _, centerline = reconstruct_centerline(data.node_coordinates,
                                           data.states[state_index],
                                           data.flexibility_matrix,
                                           data.initial_curvature)
    plot!(body_plot, centerline[1, :], centerline[2, :],
          color = colors[color_index], linewidth = 2)

    angle = root_rotation_angle(data.times[state_index],
                                data.terminal_angular_speed,
                                data.ramp_duration)
    inertial_centerline = rotation_about_x3(angle) * centerline
    plot!(inertial_plot, inertial_centerline[1, :],
          inertial_centerline[2, :],
          color = colors[color_index], linewidth = 2)
end

steady_state = reduce(hcat, data.steady_states)
_, steady_centerline = reconstruct_centerline(data.node_coordinates,
                                              steady_state,
                                              data.flexibility_matrix,
                                              data.initial_curvature)
plot!(body_plot, steady_centerline[1, :], steady_centerline[2, :],
      color = :black, linestyle = :dash, linewidth = 2, label = "steady")

energy_plot = plot(data.times, data.energy_history,
                   xlabel = "time", ylabel = "discrete energy",
                   title = "Rotating beam energy", linewidth = 2,
                   legend = false)

savefig(body_plot, joinpath(output_directory, "rotating_body_frame.pdf"))
savefig(inertial_plot,
        joinpath(output_directory, "rotating_inertial_frame.pdf"))
savefig(energy_plot, joinpath(output_directory, "rotating_energy.pdf"))
println("wrote rotating-beam figures to ", output_directory)
