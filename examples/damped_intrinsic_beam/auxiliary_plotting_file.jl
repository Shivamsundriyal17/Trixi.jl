#This file just contains some random functions used to generate the plots in the thesis


using Trixi          # Framework for DG simulations
using OrdinaryDiffEq # Differential equation solvers
using LinearAlgebra  # Linear algebra operations
using Plots          # Plotting functionality
using Printf         # Formatted output
using Dates          # Date utilities
using LaTeXStrings   # LaTeX support for text and labels
using ColorSchemes   # Color schemes for plots
using LsqFit
using DataFrames
using JLD2

# elongation_rod_length_damp
# Generate a timestamp in the format "YYYYMMDD_HHMMSS"
timestamp_str = Dates.format(now(), "yyyyMMdd_HHmmss")
timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
# Load results from the clockwise rotation simulation
@load "elongation_rod_length_damp.jld2" x t U γ κ r r_fixed L energy
x1, t1, U1, γ1, κ1, r1, r_fixed1, L1, energy1 = x, t, U, γ, κ, r, r_fixed, L, energy
@load "rotation_papermatrix_dampingperfect.jld2" x t U γ κ r r_fixed L energy
# x1, t1, U1, γ1, κ1, r1, r_fixed1, L1, energy1 = x, t, U, γ, κ, r, r_fixed, L, energy

@load "rotation_papermatrix_nodamp.jld2" x t U γ κ r r_fixed L energy
max_value = 18976.0

tip_pos = r[2,end,:]
plot(t,tip_pos, label="Tip Position", xlabel="Time", ylabel="Tip Position", legend=:topright)
plot(t, L, 
         xlabel = "t (seconds)", 
         ylabel = "Beam Length ", 
         linewidth = 3,  # Line thickness
         color=RGBA(0, 0, 1, 1),                      # Line color (blue)
         grid = true,    # Add grid for better readability
         legend = false, # Disable legend if not needed
         xguidefontsize=20,                           # Font size for X-axis label
         yguidefontsize=20,                           # Font size for Y-axis label
         xtickfontsize=14,                            # Font size for X-axis ticks
         ytickfontsize=14,                            # Font size for Y-axis ticks
         gridalpha=0.2,                               # Set grid transparency
         framestyle=:box,                             # Box-style frame for clarity
        #  xlims = (0, 50),                              # Set X-axis limits
        #  ylims = (0, 5.2),                              # Set Y-axis limits
         size = (800, 600)) 
         vline!(2, color=:red, linestyle=:dash, linewidth=2)  # Red dashed vertical line at x=2

        zero_crossing_times = t[second_derivative_zeros]
        zero_crossing_length = L[second_derivative_zeros]
        scatter!(zero_crossing_times, zero_crossing_length, color=:blue, label="Inflection Point")
     
         savefig("C:\\Users\\sund_si\\Desktop\\Solution_Plots\\discrete_energy\\beam_lengthvstime_discontinuous.pdf")

# # Load results from the anti-clockwise rotation simulation
# @load "rotation_anti_clockwise_1.jld2" x t U γ κ r r_fixed L energy
# x2, t2, U2, γ2, κ2, r2, r_fixed2, L2, energy2 = x, t, U, γ, κ, r, r_fixed, L, energy
min_value = minimum(energy)  # Slightly below the minimum energy
max_value = maximum(energy) + 100  # Slightly above the maximum energy
min_value = 0
using Plots

# Define the plot
using Plots

first_15_zeros = 1:20
t_base = t_3
energy_base = energy_upw_3
zero_crossing_times = t[second_derivative_zeros[first_15_zeros]]
zero_crossing_energies = energy[second_derivative_zeros[first_15_zeros]]
p = plot(t_3_2, energy_upw_3_2,
         xlabel="t (seconds)",                 # Proper LaTeX-style label
     ylabel="Total Discrete Energy",          
     xguidefontsize=16,                           
     yguidefontsize=16,                           
     xtickfontsize=14,                            
     ytickfontsize=14,                            
     grid=true,                                   
     gridalpha=0.2,                               
     color=RGBA(0, 0, 1, 1),                      
     linewidth=2,                                 
     framestyle=:box,                             
     legend=:topright,                         
     legendfontsize=14,                           
     label="Increased Damping: \$3\\mathbf{C}_\\tau\$",  
     size=(800, 600)       
)        

plot!(t_3, energy_upw_3, 
     label="Baseline Damping: \$\\mathbf{C}_\\tau\$", 
     color=RGBA(1, 0, 0, 1), 
     linewidth=1.4)

# Save as a high-quality PDF
savefig("C:\\Users\\sund_si\\Desktop\\Solution_Plots\\discrete_energy\\discrete_energy_comparison.pdf")


# Define the plot
p = plot(t, energy,
         xlabel="t (seconds)",                 # Proper LaTeX-style label
         ylabel="Total Discrete Energy",          # Label for Y-axis with text
         xguidefontsize=20,                           # Font size for X-axis label
         yguidefontsize=20,                           # Font size for Y-axis label
         xtickfontsize=14,                            # Font size for X-axis ticks
         ytickfontsize=14,                            # Font size for Y-axis ticks
         grid=true,                                   # Add grid for better readability
         gridalpha=0.2,                               # Set grid transparency
         color=RGBA(0, 0, 1, 1),                      # Line color (blue)
         linewidth=3,                                 # Line thickness
         framestyle=:box,                             # Box-style frame for clarity
         legend=:topright,                         # Position legend at bottom right
         legendfontsize=14,                           # Font size for legend
         label="Discrete Energy",                     # Legend label for the plot
         size=(800, 600),
         #xlims=(0,5),
         #ylims=(0,26)
         #ylims=(0,1650)
         )         # Custom Y-axis limits
           # Adjust figure size for clarity
        scatter!(zero_crossing_times, zero_crossing_energies, color=:blue, label="Inflection Point")
     
# Save as a high-quality PDF
savefig("C:\\Users\\sund_si\\Desktop\\Solution_Plots\\discrete_energy\\discontinuous_beam_energy.pdf")

# # Plot the energies
# p = plot(t1, energy1, xlabel = L"t", xguidefontsize = 20,
#          xtickfontsize = 12,
#          ytickfontsize = 12,
#          ylabel = "Total discrete energy",
#          yguidefontsize = 16,
#          label = "Clockwise Rotation",
#          fmt = :pdf,
#          color = RGBA(0, 0, 1, 1),  # Blue color for the clockwise energy
#          linewidth = 2,
#          size = (700, 500))

# # Add the anti-clockwise rotation energy to the same plot
# plot!(t2, energy2, label = "Anti-Clockwise Rotation",
#       color = RGBA(1, 0, 0, 1),  # Red color for the anti-clockwise energy
#       linewidth = 2)
# Save the plot with the timestamp
#savefig("C:\\Users\\sund_si\\Desktop\\Solution_Plots\\discrete_energy\\total_discrete_energy_plot_$timestamp_str.png")

function find_local_maxima(data)
    maxima = Int[]
    for i in 2:(length(data)-1)  # Ensure correct range
        if data[i] > data[i-1] && data[i] > data[i+1]
            push!(maxima, i)
        end
    end
    return maxima
end

function find_local_minima(data)
    minima = Int[]  # Array to store indices of local minima
    for i in 2:(length(data)-1)  # Ensure correct range
        if data[i] < data[i-1] && data[i] < data[i+1]  # Check for local minimum
            push!(minima, i)
        end
    end
    return minima
end
minima_indices = find_local_minima(energy)
local_minima = energy[minima_indices]
local_minima_times = t[minima_indices]
# Use the function to find local maxima
maxima_indices = find_local_maxima(energy)
local_maxima = energy[maxima_indices]
local_maxima_times = t[maxima_indices]

    # Exponential decay model
    function exponential_decay(t, p)
        A, γ, c = p
        return @. A * exp(-γ * t) + c
    end

    # Initial parameter guesses
    A_guess = maximum(local_maxima)
    γ_guess = 0.1
    c_guess = minimum(local_maxima)
    
    p0 = [A_guess, γ_guess, c_guess]
    
    # Perform the fit
    fit = curve_fit(exponential_decay, local_maxima_times, local_maxima, p0)
    
    # Extract parameters
    A_fit = fit.param[1]
    γ_fit = fit.param[2]
    c_fit = fit.param[3]
    half_life = log(2) / γ_fit

    t_fit = range(minimum(local_maxima_times), stop=maximum(local_maxima_times), length=500)
fitted_curve = exponential_decay(t_fit, fit.param)

plot(t, energy, label="Energy", xlabel="Time", ylabel="Energy", lw=2, color=:blue)
scatter!(local_maxima_times, local_maxima, label="Local Maxima", color=:red)
plot!(t_fit, fitted_curve, label="Fitted Exponential Decay", lw=2, color=:green)
savefig("C:\\Users\\sund_si\\Desktop\\Solution_Plots\\discrete_energy\\energy_decay_plot_$(timestamp).png")

function find_second_derivative_zeros(data, times)
    first_derivative = diff(data) ./ diff(times)
    second_derivative = diff(first_derivative) ./ diff(times[1:end-1])
    zero_crossings = Int[]  # To store indices
    for i in 1:length(second_derivative)-1
        if second_derivative[i] * second_derivative[i+1] < 0
            push!(zero_crossings, i + 1)  # +1 to account for second_derivative indexing
        end
    end

    return zero_crossings
end

# Apply the function
second_derivative_zeros = find_second_derivative_zeros(energy, t)
zero_crossing_times = t[second_derivative_zeros]
zero_crossing_energies = energy[second_derivative_zeros]
zero_crossing_length = L[second_derivative_zeros]
# Assuming you already have the variables and functions defined:
second_derivative_zeros = find_second_derivative_zeros(energy, t)  # Find second derivative zeros
first_four_zeros = second_derivative_zeros[1:4]  # Get the first 4 indices of the zeros
first_four_zeros_time = t[first_four_zeros]  # Get the times of the first 4 zeros



# Define your intervals (for time)
using Plots

# Define your intervals (for time)
intervals = [
    (0.0, first_four_zeros_time[1]),      # Interval from 0 to first value
    (first_four_zeros_time[1], first_four_zeros_time[2]),  # Interval from first value to second value
    (first_four_zeros_time[2], first_four_zeros_time[3]),  # Interval from second value to third value
    (first_four_zeros_time[3], first_four_zeros_time[4])   # Interval from third value to fourth value
]
intervals = [
    ( first_four_zeros_time[1], first_four_zeros_time[2]),      # Interval from 0 to first value
    (first_four_zeros_time[2], first_four_zeros_time[3]),  # Interval from first value to second value
    (first_four_zeros_time[3], first_four_zeros_time[4]),  # Interval from second value to third value
    (first_four_zeros_time[4], first_four_zeros_time[5])   # Interval from third value to fourth value
]
num_intervals = 4                     # Number of subfigures
num_frames = 10                        # Number of frames to show per subplot
color_gradient = cgrad([:blue, :red], num_frames)  # Generate a gradient of colors between blue and red



# Define your intervals (for time)
# intervals = [
#     (0.0, first_four_zeros_time[1]),      # Interval from 0 to first value
#     (first_four_zeros_time[1], first_four_zeros_time[2]),  # Interval from first value to second value
#     (first_four_zeros_time[2], first_four_zeros_time[3]),  # Interval from second value to third value
#     (first_four_zeros_time[3], first_four_zeros_time[4])   # Interval from third value to fourth value
# ]

# num_frames = 10                        # Number of frames to show per interval
# color_gradient = cgrad([:blue, :red], num_frames)  # Generate a gradient of colors between blue and red

# # Loop over each time interval and create individual plots
# for (interval_idx, (start, stop)) in enumerate(intervals)
#     # Create num_frames equally spaced points within the interval [start, stop]
#     points_in_interval = range(start, stop, length=num_frames)

#     # Find the closest points in t for the generated time points
#     closest_points = [t[argmin(abs.(t .- p))] for p in points_in_interval]

#     # Create a new plot for this interval
#     p = plot(size=(1000, 400), dpi=400, left_margin=5Plots.mm, bottom_margin=5Plots.mm)  # Adjust left margin

#     # Plot multiple time steps for the interval
#     for (i, closest_time) in enumerate(closest_points)
#         # Find the closest index in t to the generated point
#         idx = argmin(abs.(t .- closest_time))
        
#         # Get color and alpha
#         color = color_gradient[i]  # Use color from gradient
#         max_opacity = 0.8         # Maximum opacity
#         min_opacity = 0.3         # Minimum opacity for the last frame
#         alpha_plot = max_opacity - (max_opacity - min_opacity) * (i / (num_frames - 1))  # Smoothly decreasing alpha

#         # Plot closest points in the interval
#         plot!(p, r[1, :, idx], r[2, :, idx],  # Assuming r holds your data
#               label="t = $(round(t[idx], digits=2))",  # Label for each time step in the legend
#               color=color,  # Use color from gradient
#               alpha=alpha_plot,  # Gradually decreasing opacity
#               linewidth=1.5)  # Thinner line for scientific clarity
#     end

#     # Customize the plot
#     plot!(p,
#         xlims=(0, 4.9), ylims=(-1.0, 1.0),
#         xlabel = L"x_1", ylabel = L"x_2",
#         legend=:topleft,  # Move legend outside the plot area to the right
#         grid = true,  # Remove grid lines for cleaner visuals
#         tickfontsize=10,  # Smaller ticks for more compactness
#         guidefontsize=12)  # Slightly smaller font for axis labels and title

#     # Save the plot as an individual file
#     savefig(p, "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\beam_movingcam_$(interval_idx).png")
# end
# using Plots

# # Define your intervals (for time)
intervals = [
    (0.0, first_four_zeros_time[1]),      # Interval from 0 to first value
    (first_four_zeros_time[1], first_four_zeros_time[2]),  # Interval from first value to second value
    (first_four_zeros_time[2], first_four_zeros_time[3]),  # Interval from second value to third value
    (first_four_zeros_time[3], first_four_zeros_time[4])   # Interval from third value to fourth value
]

num_frames = 10                        # Number of frames to show per interval
color_gradient = cgrad([:blue, :red], num_frames)  # Generate a gradient of colors between blue and red

# Loop over each time interval and create individual plots
for (interval_idx, (start, stop)) in enumerate(intervals)
    # Create num_frames equally spaced points within the interval [start, stop]
    points_in_interval = range(start, stop, length=num_frames)

    # Find the closest points in t for the generated time points
    closest_points = [t[argmin(abs.(t .- p))] for p in points_in_interval]

    # Create a new plot for this interval with two subplots
    p = plot(layout=(1, 2), size=(1200, 400), dpi=400, left_margin=5Plots.mm, bottom_margin=5Plots.mm)

    # Plot for the moving camera perspective
    for (i, closest_time) in enumerate(closest_points)
        # Find the closest index in t to the generated point
        idx = argmin(abs.(t .- closest_time))
        
        # Get color and alpha
        color = color_gradient[i]  # Use color from gradient
        max_opacity = 0.8         # Maximum opacity
        min_opacity = 0.3         # Minimum opacity for the last frame
        alpha_plot = max_opacity - (max_opacity - min_opacity) * (i / (num_frames - 1))  # Smoothly decreasing alpha

        # Plot moving camera perspective
        plot!(p[1], r[1, :, idx], r[2, :, idx],  # Assuming r holds your data
              #label="t = $(round(t[idx], digits=2))",
              label =false,
              color=color,
              alpha=alpha_plot,
              linewidth=1.5)
    end

    # Customize moving camera subplot
    plot!(p[1],
        xlims=(0, 4.9), ylims=(-1.0, 1.0),
        xlabel = L"x_1", ylabel = L"x_2",
        grid=true)

    # Plot for the fixed camera perspective
    for (i, closest_time) in enumerate(closest_points)
        # Find the closest index in t to the generated point
        idx = argmin(abs.(t .- closest_time))
         # Get color and alpha
         color = color_gradient[i]  # Use color from gradient
         max_opacity = 0.8         # Maximum opacity
         min_opacity = 0.3         # Minimum opacity for the last frame
         alpha_plot = max_opacity - (max_opacity - min_opacity) * (i / (num_frames - 1))  # Smoothly decreasing alpha
        
        # Plot fixed camera perspective
        plot!(p[2], r_fixed[1, :, idx], r_fixed[2, :, idx],  # Assuming r_fixed holds your data
        label="t = $(round(t[idx], digits=2))",
        color=color,
        alpha=alpha_plot,
        linewidth=1.5,
        legendfontsize=14)  # Adjust this value to make the legend text larger
  
    end

    # Customize fixed camera subplot
    plot!(p[2],
        xlims=(-5, 5), ylims=(-5.0, 5.0),
        xlabel = L"x_1", ylabel = L"x_2",
        grid=true)

    # Add a shared legend
    plot!(p, legend=:outerright, tickfontsize=10, guidefontsize=12)

    # Save the plot as an individual file
    savefig(p, "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\beam_trial_interval_$(interval_idx).png")
end

start_index = second_derivative_zeros[1]  # Replace with your desired starting maximum index
end_index = second_derivative_zeros[40] 
# Time steps between the two local maxima
time_range_indices = start_index:end_index
total_steps = length(time_range_indices)

# Define animation step size (adjust skip as needed)
target_frames = 250  # Desired number of frames in the GIF
skip = max(1, div(total_steps, target_frames))

# Generate the animation
anim = @animate for n in 1:skip:total_steps
    current_index = time_range_indices[n]  # Current index for the beam data
    current_time = t[current_index]    # Current time step
    
    # Rotating Beam Plot
    rotating_beam_plot = plot(
        r[1, :, current_index],  # X-coordinates of the rotating beam
        r[2, :, current_index],  # Y-coordinates of the rotating beam
        xlims=(0, 5.0),           # X-axis limits
        ylims=(-3, 3.0),           # Y-axis limits
        linewidth=3,
        aspect_ratio=:equal,
        legend=false,
        title="Rotating Perspective at t = $(round(current_time, digits=2))"
    )
    
    # Fixed Perspective Beam Plot
    fixed_beam_plot = plot(
        r_fixed[1, :, current_index],  # X-coordinates of the fixed perspective
        r_fixed[2, :, current_index],  # Y-coordinates of the fixed perspective
        xlims=(-6, 6),                 # X-axis limits
        ylims=(-6, 6),                 # Y-axis limits
        linewidth=3,
        aspect_ratio=:equal,
        legend=false,
        title="Fixed Perspective at t = $(round(current_time, digits=2))"
    )
    
    # Energy Plot with Moving Marker
    energy_plot = plot(
        t, energy,                # Energy vs. time
        label="Energy",
        xlabel="Time",
        ylabel="Energy",
        linewidth=2,
        title="Energy vs. Time"
    )
    scatter!([current_time], [energy[current_index]], color=:red, label="Current Time")  # Moving marker
    # Highlight zero-crossings of the second derivative
    # zero_crossing_times = t[second_derivative_zeros]
    # zero_crossing_energies = energy[second_derivative_zeros]
    # scatter!(zero_crossing_times, zero_crossing_energies, color=:blue, label="Zero Crossings (2nd Derivative)")

    # Combine the three plots
    plot(rotating_beam_plot, fixed_beam_plot, energy_plot, layout = @layout([a b; c]), size=(1000, 800))
end

# Save the animation as a GIF
timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
gif(anim, "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\combined_animation_with_fixed_full_$(timestamp).gif", fps=30)

second_derivative_zeros = find_second_derivative_zeros(energy, t)  # Find second derivative zeros
first_four_zeros = second_derivative_zeros[1:4]  # Get the first 4 indices of the zeros
first_four_zeros_time = t[first_four_zeros]
# Define intervals
intervals = [
    (1, first_four_zeros[1]),      # Interval from 0 to first value
    (first_four_zeros[1], first_four_zeros[2]),  # Interval from first value to second value
    (first_four_zeros[2], first_four_zeros[3]),  # Interval from second value to third value
    (first_four_zeros[3], first_four_zeros[4])   # Interval from third value to fourth value
]

intervals = [
    (0.0, first_four_zeros_time[1]),      # Interval from 0 to first value
    (first_four_zeros_time[1], first_four_zeros_time[2]),  # Interval from first value to second value
    (first_four_zeros_time[2], first_four_zeros_time[3]),  # Interval from second value to third value
    (first_four_zeros_time[3], first_four_zeros_time[4])   # Interval from third value to fourth value
]
intervals = [
    ( first_four_zeros_time[1], first_four_zeros_time[2]),      # Interval from 0 to first value
    (first_four_zeros_time[2], first_four_zeros_time[3]),  # Interval from first value to second value
    (first_four_zeros_time[3], first_four_zeros_time[4]),  # Interval from second value to third value
    (first_four_zeros_time[4], first_four_zeros_time[5])   # Interval from third value to fourth value
]
num_frames = 10                        # Number of frames to show per interval
color_gradient = cgrad([:blue, :red], num_frames)  # Generate a gradient of colors between blue and red

# Loop over each time interval and create individual plots
for (interval_idx, (start, stop)) in enumerate(intervals)
    # Create num_frames equally spaced points within the interval [start, stop]
    points_in_interval = range(start, stop, length=num_frames)
    print(points_in_interval)
    # Find the closest points in t for the generated time points
    closest_points = points_in_interval
    #closest_points = [t[argmin(abs.(t .- p))] for p in points_in_interval]
    print(closest_points)
    # Create plot for the moving camera perspective
    p_moving = plot(size=(800, 600), dpi=400)
    
    # Plot for the moving camera perspective
    for (i, closest_time) in enumerate(closest_points)
        # Find the closest index in t to the generated point
        idx = argmin(abs.(t .- closest_time))
        
        # Get color and alpha
        color = color_gradient[i]  # Use color from gradient
        max_opacity = 0.9         # Maximum opacity
        min_opacity = 0.4         # Minimum opacity for the last frame
        alpha_plot = max_opacity - (max_opacity - min_opacity) * (i / (num_frames - 1))  # Smoothly decreasing alpha

        # Plot moving camera perspective
        plot!(p_moving, r[1, :, idx], r[2, :, idx],  # Assuming r holds your data
              label=false,
              color=color,
              alpha=alpha_plot,
              linewidth=1.5)
    end

    # Customize moving camera plot
    plot!(p_moving,
        xlims=(0, 4.5), ylims=(-2.5, 2.5),
        xlabel = L"x_1", ylabel = L"x_2",
        xlabelfontsize=14, ylabelfontsize=14,  # Increase the font size of the labels
        xtickfontsize=12, ytickfontsize=12,    # Increase the font size of tick labels
        grid=true)

    # Save the plot for moving camera as an individual file
    savefig(p_moving, "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\rotating_beam_interval_1_$(interval_idx)_moving.pdf")

    # Create plot for the fixed camera perspective
    p_fixed = plot(size=(800, 600), dpi=400)

    # Plot for the fixed camera perspective
    for (i, closest_time) in enumerate(closest_points)
        # Find the closest index in t to the generated point
        idx = argmin(abs.(t .- closest_time))
        
        # Get color and alpha
        color = color_gradient[i]  # Use color from gradient
        max_opacity = 0.9         # Maximum opacity
        min_opacity = 0.4         # Minimum opacity for the last frame
        alpha_plot = max_opacity - (max_opacity - min_opacity) * (i / (num_frames - 1))  # Smoothly decreasing alpha
        
        # Plot fixed camera perspective
        plot!(p_fixed, r_fixed[1, :, idx], r_fixed[2, :, idx],  # Assuming r_fixed holds your data
              label="t = $(round(t[idx], digits=2))",
              color=color,
              alpha=alpha_plot,
              linewidth=1.5,
              legendfontsize=18)  # Adjust this value to make the legend text larger
    end

    # Customize fixed camera plot
    plot!(p_fixed,legend=:outerright,
        xlims=(-5, 5), ylims=(-5.0, 5.0),
        xlabel = L"x_1", ylabel = L"x_2",
        xlabelfontsize=14, ylabelfontsize=14,  # Increase the font size of the labels
        xtickfontsize=12, ytickfontsize=12,    # Increase the font size of tick labels
        grid=true)

    # Save the plot for fixed camera as an individual file
    savefig(p_fixed, "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\beam_interval_fifth_$(interval_idx)_fixed.pdf")
end



#########
using Plots
using Colors

# Define your time and energy data

# Create a color gradient
color_gradient = cgrad([:blue, :red])  # Define the color gradient

# Generate colors corresponding to time points
colors = [color_gradient[ti / maximum(t)] for ti in t]

# Plot with gradient-based colors
p = scatter(t, energy, 
            color=colors,           # Map colors to time points
            markerstrokewidth=0.0, # Remove marker borders for clarity
            xlabel="t (seconds)", 
            ylabel="Total Discrete Energy",
            xguidefontsize=20,
            yguidefontsize=20,
            xtickfontsize=14,
            ytickfontsize=14,
            grid=true,
            gridalpha=0.5,
            framestyle=:box,
            legend=false,           # Disable legend for this plot
            size=(800, 600))

# Save the plot as a PDF
savefig("C:\\Users\\sund_si\\Desktop\\Solution_Plots\\discrete_energy\\total_discrete_energy_plot.pdf")

################
second_derivative_zeros = find_second_derivative_zeros(energy, t)  # Find second derivative zeros
first_four_zeros = second_derivative_zeros[4:8]  # Get the first 4 indices of the zeros
intervals = [
    ( first_four_zeros[1], first_four_zeros[2]),      # Interval from 0 to first value
    (first_four_zeros[2], first_four_zeros[3]),  # Interval from first value to second value
    (first_four_zeros[3], first_four_zeros[4]),  # Interval from second value to third value
    (first_four_zeros[4], first_four_zeros[5])   # Interval from third value to fourth value
]

intervals = [
    ( second_derivative_zeros[734], second_derivative_zeros[950]),      # Interval from 0 to first value
    (first_four_zeros[2], first_four_zeros[3]),  # Interval from first value to second value
    (first_four_zeros[3], first_four_zeros[4]),  # Interval from second value to third value
    (first_four_zeros[4], first_four_zeros[5])   # Interval from third value to fourth value
]
# Define intervals using indices directly
intervals = [
    (1, first_four_zeros[1]),      # Interval from 1 to the second zero index
    (first_four_zeros[1], first_four_zeros[2]),  # Interval from second to third zero index
    (first_four_zeros[2], first_four_zeros[3]),  # Interval from third to fourth zero index
    (first_four_zeros[3], first_four_zeros[4])   # Interval from fourth zero index to the end
]
num_frames = 10                        # Number of frames to show per interval
color_gradient = cgrad([:blue, :red], num_frames)  # Generate a gradient of colors between blue and red

# Loop over each time interval and create individual plots
for (interval_idx, (start_idx, stop_idx)) in enumerate(intervals)
    # Create num_frames equally spaced indices within the interval [start_idx, stop_idx]
    indices_in_interval = range(start_idx, stop_idx, length=num_frames)

    # Convert indices to integers
    indices_in_interval = round.(Int, indices_in_interval)

    # Create plot for the moving camera perspective
    p_moving = plot(size=(800, 600), dpi=400)
    
    # Plot for the moving camera perspective
    for (i, idx) in enumerate(indices_in_interval)
        # Get color and alpha
        color = color_gradient[i]  # Use color from gradient
        max_opacity = 0.9         # Maximum opacity
        min_opacity = 0.4         # Minimum opacity for the last frame
        alpha_plot = max_opacity - (max_opacity - min_opacity) * (i / (num_frames - 1))  # Smoothly decreasing alpha

        # Plot moving camera perspective
        plot!(p_moving, r[1, :, idx], r[2, :, idx],  # Assuming r holds your data
              label=false,
              color=color,
              alpha=alpha_plot,
              linewidth=3.0)
    end

    # Customize moving camera plot
    plot!(p_moving,
        xlims=(0, 5.0), ylims=(-2.7, 2.7),
        xlabel = L"x_1", ylabel = L"x_2",
        xlabelfontsize=22, ylabelfontsize=22,  # Increase the font size of the labels
        xtickfontsize=16, ytickfontsize=16,    # Increase the font size of tick labels
        grid=true)

    # Save the plot for moving camera as an individual file
    savefig(p_moving, "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\rotating_beam_plot\\rotating_beam_moving$(interval_idx)_3.pdf")

    # Create plot for the fixed camera perspective
    p_fixed = plot(size=(800, 600), dpi=400)

    # Plot for the fixed camera perspective
    for (i, idx) in enumerate(indices_in_interval)
        # Get color and alpha
        color = color_gradient[i]  # Use color from gradient
        max_opacity = 0.9         # Maximum opacity
        min_opacity = 0.4         # Minimum opacity for the last frame
        alpha_plot = max_opacity - (max_opacity - min_opacity) * (i / (num_frames - 1))  # Smoothly decreasing alpha
        
        # Plot fixed camera perspective
        plot!(p_fixed, r_fixed[1, :, idx], r_fixed[2, :, idx],  # Assuming r_fixed holds your data
              label="t = $(round(t[idx], digits=2))",
              color=color,
              alpha=alpha_plot,
              linewidth=2.5,
              legendfontsize=20)  # Adjust this value to make the legend text larger
    end

    # Customize fixed camera plot
    plot!(p_fixed,legend=:outerright,
        xlims=(-5, 5), ylims=(-5.0, 5.0),
        xlabel = L"x_1", ylabel = L"x_2",
        xlabelfontsize=22, ylabelfontsize=22,  # Increase the font size of the labels
        xtickfontsize=16, ytickfontsize=16,    # Increase the font size of tick labels
        grid=true)

    # Save the plot for fixed camera as an individual file
    savefig(p_fixed, "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\rotating_beam_plot\\rotating_beam_fixed$(interval_idx)_3.pdf")
end
