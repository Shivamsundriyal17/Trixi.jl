using Plots
using Serialization: deserialize

include(joinpath(@__DIR__, "beam_geometry.jl"))

input_file = isempty(ARGS) ?
             joinpath(@__DIR__, "..", "results", "nonsmooth_resultant.jls") :
             abspath(ARGS[1])
output_directory = length(ARGS) < 2 ?
                   joinpath(@__DIR__, "..", "results") :
                   abspath(ARGS[2])
mkpath(output_directory)
data = deserialize(input_file)

target_times = (0.0, 0.02, 0.05, 0.1, 1.0, data.times[end])
snapshot_indices = [argmin(abs.(data.times .- target_time))
                    for target_time in target_times]
colors = palette(:plasma, length(snapshot_indices))

strain_plot = plot(xlabel = "x", ylabel = "axial strain gamma1",
                   title = "Localized resultant relaxation")
configuration_plot = plot(aspect_ratio = :equal,
                          xlabel = "x1", ylabel = "x2",
                          title = "Beam centerline", legend = :outerright)

order = sortperm(data.node_coordinates; alg = Base.Sort.MergeSort)
for (color_index, state_index) in enumerate(snapshot_indices)
    state = reshape(data.states[state_index], 12, :)
    axial_strain = vec(transpose(data.flexibility_matrix[1, :]) *
                       state[7:12, :])
    label = "t=$(round(data.times[state_index], digits = 3))"
    plot!(strain_plot, data.node_coordinates[order], axial_strain[order],
          color = colors[color_index], linewidth = 2, label = label)

    _, centerline = reconstruct_centerline(data.node_coordinates,
                                           data.states[state_index],
                                           data.flexibility_matrix,
                                           data.initial_curvature)
    plot!(configuration_plot, centerline[1, :], centerline[2, :],
          color = colors[color_index], linewidth = 2, label = label)
end

diagnostics_plot = plot(data.times, data.maximum_axial_strain,
                        xlabel = "time", ylabel = "maximum absolute value",
                        label = "max |gamma1|", linewidth = 2)
plot!(diagnostics_plot, data.times, data.energy_history,
      label = "energy", linewidth = 2)

savefig(strain_plot,
        joinpath(output_directory, "nonsmooth_axial_strain.pdf"))
savefig(configuration_plot,
        joinpath(output_directory, "nonsmooth_centerline.pdf"))
savefig(diagnostics_plot,
        joinpath(output_directory, "nonsmooth_diagnostics.pdf"))
println("wrote nonsmooth-resultant figures to ", output_directory)
