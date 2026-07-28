using Plots
using Serialization: deserialize

include(joinpath(@__DIR__, "beam_geometry.jl"))

input_directory = isempty(ARGS) ?
                  joinpath(@__DIR__, "..", "results",
                           "rotating_beam_campaign") :
                  abspath(ARGS[1])
output_file = length(ARGS) < 2 ?
              joinpath(input_directory, "rotating_beam_validation.pdf") :
              abspath(ARGS[2])
mkpath(dirname(output_file))

baseline = deserialize(joinpath(input_directory, "baseline.jls"))
double_damping = deserialize(joinpath(input_directory,
                                      "double_damping.jls"))
undamped = deserialize(joinpath(input_directory, "undamped.jls"))

blue = RGB(0 / 255, 114 / 255, 178 / 255)
orange = RGB(213 / 255, 94 / 255, 0 / 255)
green = RGB(0 / 255, 158 / 255, 115 / 255)
magenta = RGB(204 / 255, 121 / 255, 167 / 255)
gray = RGB(0.25, 0.25, 0.25)

configuration_plot = plot(xlabel = "x₁", ylabel = "x₂",
                          title = "(a) Co-rotating centerline",
                          legend = :bottomleft, grid = :on)
tip_displacements = map(eachindex(baseline.times)) do index
    _, centerline = reconstruct_centerline(baseline.node_coordinates,
                                           baseline.states[index],
                                           baseline.flexibility_matrix,
                                           baseline.initial_curvature)
    centerline[2, end]
end
ramp_index = argmin(abs.(baseline.times .- baseline.ramp_duration))
rebound_search = (ramp_index + 1):(length(baseline.times) - 1)
rebound_index = findfirst(rebound_search) do index
    tip_displacements[index - 1] < tip_displacements[index] &&
        tip_displacements[index] >= tip_displacements[index + 1] &&
        tip_displacements[index] > 0
end
rebound_index = isnothing(rebound_index) ?
                argmax(tip_displacements) :
                rebound_search[rebound_index]
late_index = argmin(abs.(baseline.times .-
                         min(10.0, baseline.times[end])))
configuration_indices = unique((1, ramp_index, rebound_index,
                                late_index, length(baseline.times)))
configuration_colors = (gray, orange, magenta, green, blue)
for (index, color) in zip(configuration_indices, configuration_colors)
    time = baseline.times[index]
    _, centerline = reconstruct_centerline(baseline.node_coordinates,
                                           baseline.states[index],
                                           baseline.flexibility_matrix,
                                           baseline.initial_curvature)
    plot!(configuration_plot, centerline[1, :], centerline[2, :],
          color = color, linewidth = 2,
          label = "t=$(round(time, digits = 2))")
    scatter!(configuration_plot, [centerline[1, end]],
             [centerline[2, end]], color = color, markersize = 3,
             markerstrokewidth = 0, label = "")
end
steady_state = reduce(hcat, baseline.steady_states)
_, steady_centerline = reconstruct_centerline(baseline.node_coordinates,
                                              steady_state,
                                              baseline.flexibility_matrix,
                                              baseline.initial_curvature)
plot!(configuration_plot, steady_centerline[1, :],
      steady_centerline[2, :], color = :black, linestyle = :dash,
      linewidth = 2, label = "analytic steady")
scatter!(configuration_plot, [steady_centerline[1, end]],
         [steady_centerline[2, end]], color = :black, markersize = 3,
         marker = :star5, label = "")

energy_plot = plot(xlabel = "time", ylabel = "Eₕ / E⋆",
                   title = "(b) Normalized energy",
                   legend = :topright, grid = :on)
for (data, label, color, style) in ((baseline, "Cτ", blue, :solid),
                                    (double_damping, "2Cτ", orange, :dash),
                                    (undamped, "undamped", gray, :dot))
    plot!(energy_plot, data.times,
          data.energy_history ./ baseline.diagnostics.energy_steady,
          color = color, linestyle = style, linewidth = 2,
          label = label)
end
hline!(energy_plot, [1.0], color = :black, linestyle = :dashdot,
       linewidth = 1.5, label = "analytic steady")

error_indices = findall(time -> time >= baseline.ramp_duration,
                        baseline.times)
error_plot = plot(baseline.times[error_indices],
                  baseline.f1_relative_linf_history[error_indices],
                  yscale = :log10, color = blue, linewidth = 2,
                  xlabel = "time", ylabel = "relative L∞ error",
                  title = "(c) Relaxation to steady branch",
                  label = "f₁", legend = :topright, grid = :on)
plot!(error_plot, baseline.times[error_indices],
      baseline.v2_relative_linf_history[error_indices],
      color = orange, linewidth = 2, linestyle = :dash,
      label = "v₂")

order = sortperm(baseline.node_coordinates; alg = Base.Sort.MergeSort)
x = baseline.node_coordinates[order]
final_state = reshape(baseline.states[end], 12, :)[:, order]
reference_state = reduce(hcat, baseline.steady_states)[:, order]
f1_scale = maximum(abs, reference_state[7, :])
v2_scale = maximum(abs, reference_state[2, :])
profile_plot = plot(x, reference_state[7, :] ./ f1_scale,
                    color = blue, linestyle = :dash, linewidth = 2,
                    xlabel = "x", ylabel = "normalized field",
                    title = "(d) Final profiles at T=40",
                    label = "f₁ analytic", legend = :right, grid = :on)
scatter!(profile_plot, x, final_state[7, :] ./ f1_scale,
         color = blue, marker = :circle, markersize = 2.5,
         markerstrokewidth = 0, label = "f₁ Trixi")
plot!(profile_plot, x, reference_state[2, :] ./ v2_scale,
      color = orange, linestyle = :dash, linewidth = 2,
      label = "v₂ analytic")
scatter!(profile_plot, x, final_state[2, :] ./ v2_scale,
         color = orange, marker = :diamond, markersize = 2.5,
         markerstrokewidth = 0, label = "v₂ Trixi")

figure = plot(configuration_plot, energy_plot, error_plot, profile_plot;
              layout = (2, 2), size = (1120, 820),
              left_margin = 5Plots.mm, bottom_margin = 4Plots.mm,
              top_margin = 3Plots.mm, titlefontsize = 11,
              guidefontsize = 10, tickfontsize = 9,
              legendfontsize = 8)
savefig(figure, output_file)
println("wrote ", output_file)
