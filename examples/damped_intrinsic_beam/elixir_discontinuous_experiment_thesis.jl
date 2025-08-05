#Experiment with discontinuous initial condiiton from the thesis

using Trixi          # Trixi framework for DG simulations
using OrdinaryDiffEq # Differential equation solvers
using LinearAlgebra  # Linear algebra operations
using Plots          # Plotting functionality
using Printf         # Formatted output
using Dates           # Date utilities
using LaTeXStrings    # LaTeX support for text and labels
using JLD


flexMat = inv(diagm([10.0^3, 10.0^3, 10.0^3, 500.0, 500.0, 500.0]))
dampingMat =2*inv(diagm([10.0^3, 10.0^3, 10.0^3, 10.0, 10.0, 10.0]))
massMat = diagm([1.0, 1.0, 1.0, 20.0, 10.0, 10.0])

# Initial curvature and external forces/moments
k_init = [0.0 0.0 0.0]  # Initial curvature
f_ext_func(x, t) = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0]  # External forces/moments function

# Set up hyperbolic equations with damping
equations_hyperbolic = Damped_Full_Hyperbolic(flexMat, massMat, dampingMat, k_init,
                                              f_ext_func)

# Set up parabolic equations based on hyperbolic equations
equations_parabolic = Damped_Full_Parabolic(equations_hyperbolic,
                                            dampingMat, flexMat, massMat,
                                            equations_hyperbolic.Gamma_inv,
                                            equations_hyperbolic.Psi,
                                            k_init)

@inline function (boundary_condition_damped)(flux_inner, u_inner, normal::Int64,
                                             direction::Int64,
                                             x, t, operator_type::Trixi.Divergence,
                                             equations_parabolic::Damped_Full_Parabolic)
    if direction == 1
        return [flux_inner[1:6]; zeros(6)]

    elseif direction == 2
        return [zeros(6); flux_inner[7:12]]

    else
        error("Invalid direction specified. Use 1 for left boundary and 2 for right boundary.")
    end
end

@inline function (boundary_condition_damped)(flux_inner, u_inner, normal::Any,
                                             direction::Int64,
                                             x, t, operator_type::Trixi.Gradient,
                                             equations_parabolic::Damped_Full_Parabolic)
    u_boundary = zeros(12)
    if direction == 1
        y_B = zeros(6)  # Boundary condition vector for velocity (left boundary)
        u_boundary = [y_B[1:6]; u_inner[7:12]]

    elseif direction == 2
        s_B = zeros(6)  # Boundary condition vector for strain (right boundary)
        u_boundary = [u_inner[1:6]; u_inner[7:12]]

    else
        error("Invalid direction specified. Use 1 for left boundary and 2 for right boundary.")
    end

    return u_boundary
end

# Define boundary conditions
boundary_conditions_hyperbolic = boundary_conditions_ib
boundary_conditions_parabolic = boundary_condition_damped

# Define the Discontinuous Galerkin Spectral Element Method (DGSEM) solver
solver = DGSEM(polydeg = 4,                 # Polynomial degree for DGSEM
               surface_flux = flux_upwind)

# Set coordinates for the mesh
coordinates_min = (0.0,)  # Minimum coordinates (e.g., min(x), min(y))
coordinates_max = (1.0,)  # Maximum coordinates (e.g., max(x), max(y))

# Create the mesh using TreeMesh
mesh = TreeMesh(coordinates_min,                      # Minimum coordinates of the domain
                coordinates_max,                      # Maximum coordinates of the domain
                initial_refinement_level = 4,           # Initial level of mesh refinement
                n_cells_max = 10_000,                   # Maximum number of cells in the mesh
                periodicity = false)



function initial_condition_strain_discontinuity_full(x, t, ::Damped_Full_Hyperbolic)
    return SVector{12}([zeros(6); 3000.0 * (0.5 < x[1] < 0.75); zeros(5)])
end


# Example usage: select an initial condition
initial_condition = initial_condition_strain_discontinuity_full

# Define the flag Elongationrod
Elongationrod = false

# Set up semidiscretization with hyperbolic and parabolic equations
semi = SemidiscretizationHyperbolicParabolic(mesh,
                                             (equations_hyperbolic, equations_parabolic),
                                             initial_condition,
                                             solver;
                                             boundary_conditions = (boundary_conditions_hyperbolic,
                                                                    boundary_conditions_parabolic))
# Define time span for the simulation
tspan = (0.0, 5.0)

# Perform semidiscretization
ode = semidiscretize(semi, tspan)

# Define callbacks
stepsize_callback = StepsizeCallback()
alive_callback = AliveCallback(alive_interval = 1000)
callbacks = CallbackSet(alive_callback)

# Define time integration tolerance
time_int_tol = 1.0e-6

# Solve the ODE with specified options
sol = solve(ode,
            RDPK3SpFSAL49();
            abstol = time_int_tol,
            reltol = time_int_tol,
            ode_default_options()...,  # Ensure to pass all default options
            save_everystep = true,      # Save results at every step
            callback = callbacks)

###########################################################
# Postprocessing and Plotting
###########################################################


timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")

# spatial nodes
x = semi.cache.elements._node_coordinates
Nx = length(x)

# time nodes
t = sol.t
Nt = length(t)

# reshape solution vector into 3D array (variable, node, time)
U = Array{Float64}(undef, 12, Nx, Nt)
for j in 1:Nt
    U[:, :, j] = reshape_solution(sol.u[j])
end

flattened_sol_u = reduce(vcat, sol.u)
max_value = maximum(flattened_sol_u) + 0.5
min_value = minimum(flattened_sol_u) - 0.5

# Initialize array for strain measures
γ = Array{Float64}(undef, 3, Nx, Nt)
κ = Array{Float64}(undef, 3, Nx, Nt)
# Compute strain measures

for j in 1:Nt
    strain_result = flexMat * U[7:12, :, j]
    γ[:, :, j] = strain_result[1:3, :]
    κ[:, :, j] = strain_result[4:6, :]
end
global_min = minimum(γ[1, :, :])
global_max = maximum(γ[1, :, :])

function plot_strain(γ, step, x)
    p = plot(x, γ[1, :, step],
             #xlims = (0.0, 2.5),
             ylims = (global_min, global_max+1),
             legend = false,
             xlabel = "Spatial Coordinate (x)",
             ylabel = "Strain γ₁",
             xguidefontsize = 18,                # Font size for X-axis label
             yguidefontsize = 18,                # Font size for Y-axis label
             xtickfontsize = 12,                 # Font size for X-axis ticks
             ytickfontsize = 12,                 # Font size for Y-axis ticks
             linewidth = 2,
             color = :red,
             size = (700, 500),
             framestyle = :box)
    return p
end
S0 = diagm([1.0; 1.0; 1.0])
k_func(x) = [0.0; 0.0; 0.0]
# determine position line of beam
r = Array{Float64}(undef, 3, Nx, Nt)
for j in 1:Nt
    r[:, :, j] = calc_deformation_trapez(S0, x, k_func, κ[:, :, j], γ[:, :, j])
end

L = zeros(Nt)
for n in 1:Nt
    diffs = diff(r[:, :, n], dims = 2)  # Compute differences between consecutive points
    L[n] = sum(sqrt.(sum(diffs .^ 2, dims = 1)))  # Sum Euclidean distances
end
# Function to plot the beam deformation and node positions
function plot_deformation_with_nodes(r, step, x, γ)
    p = plot(r[1, :, step], r[3, :, step], line_z = γ[1, :, step],
             xlims = (0,3.0), ylims = (global_min, global_max),
             xguidefontsize = 10,                # Font size for X-axis label
             yguidefontsize = 10,                # Font size for Y-axis label
             xtickfontsize = 10,                 # Font size for X-axis ticks
             ytickfontsize = 10,                 # Font size for Y-axis ticks
             clims = (global_min, global_max),
             xlabel = L"x_1", ylabel = L"x_3",
             aspect_ratio = :equal, linewidth = 3,
             legend = false, cbar = true, cbar_title = L"γ_1",
             title = string("Deformation at t = ", round(t[step], digits = 2)))

    # Add smaller markers for the first and last nodes
    scatter!(p, [r[1, 1, step]], [r[3, 1, step]], label = "Start Node", color = :green,
             marker = :circle, markersize = 4)
    scatter!(p, [r[1, end, step]], [r[3, end, step]], label = "End Node", color = :blue,
             marker = :circle, markersize = 4)

    return p
end






target_frames = 200
x_min = minimum(r[1, :, :])
x_max = maximum(r[1, :, :])
x_limits = (x_min, x_max + 0.5)
plot_frequency = max(1, div(Nt, target_frames))

# Create the animation with two plots side by side
anim_force_x = @animate for n in 1:plot_frequency:Nt
    p1 = plot_deformation_with_nodes(r, n, x, γ)
    p2 = plot_strain(γ, n, x)
    plot(p1, p2, layout = @layout [a{0.5w} b{0.5w}])
end

#Save the combined animation as a GIF
gif(anim_force_x,
    "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\beam_deformation\\beam_$(timestamp).gif",
    fps = 10)

function plot_linearvelocities(pd, step, semi)
    # Extract the time step for the title
    time_step = round(sol.t[step], digits = 5)
    pd = PlotData1D(sol.u[step], semi)

    # Create individual plots with increased linewidth
    p1 = plot(pd["V1"], xguide = "Spatial coordinate", yguide = "V1",
              xlims = (0, 2), ylims = (min_value, max_value),
              title = "V1 Time t=$time_step", linewidth = 2,
              titlefont_size = 10,  # Adjust title font size
              xguidefontsize = 8,  # Adjust x-axis label font size
              yguidefontsize = 8)  # Adjust y-axis tick label font size

    p2 = plot(pd["V2"], xguide = "Spatial coordinate", yguide = "V2",
              xlims = (0, 2), ylims = (min_value, max_value),
              title = "V2", linewidth = 2,
              titlefont_size = 10,  # Adjust title font size
              xguidefontsize = 8,  # Adjust x-axis label font size
              yguidefontsize = 8)  # Adjust y-axis tick label font size

    p3 = plot(pd["V3"], xguide = "Spatial coordinate", yguide = "V3",
              xlims = (0, 2), ylims = (min_value, max_value),
              title = "V3", linewidth = 2,
              titlefont_size = 10,  # Adjust title font size
              xguidefontsize = 8,  # Adjust x-axis label font size
              yguidefontsize = 8)  # Adjust y-axis tick label font size

    # Combine the plots into a single figure with a 2x2 layout for compactness
    combined_plot = plot(p1, p2, p3)  # Adjust size and layout
    return combined_plot
end

function plot_linearstrains(pd, step, semi)
    # Extract the time step for the title
    time_step = round(sol.t[step], digits = 5)
    pd = PlotData1D(sol.u[step], semi)

    # Create individual plots with increased linewidth
    p1 = plot(pd["x2_1"], xguide = "Spatial coordinate", yguide = "Linear_Strain_1",
              xlims = (0, 2), ylims = (min_value, max_value),
              title = "x2_1 Time t=$time_step", linewidth = 2,
              titlefont_size = 10,  # Adjust title font size
              xguidefontsize = 8,  # Adjust x-axis label font size
              yguidefontsize = 8)  # Adjust y-axis tick label font size

    p2 = plot(pd["x2_2"], xguide = "Spatial coordinate", yguide = "Linear_Strain_2",
              xlims = (0, 2), ylims = (min_value, max_value),
              title = "x2_2", linewidth = 2,
              titlefont_size = 10,  # Adjust title font size
              xguidefontsize = 8,  # Adjust x-axis label font size
              yguidefontsize = 8)  # Adjust y-axis tick label font size

    p3 = plot(pd["x2_3"], xguide = "Spatial coordinate", yguide = "Linear_Strain_3",
              xlims = (0, 2), ylims = (min_value, max_value),
              title = "V3", linewidth = 2,
              titlefont_size = 10,  # Adjust title font size
              xguidefontsize = 8,  # Adjust x-axis label font size
              yguidefontsize = 8)  # Adjust y-axis tick label font size

    # Combine the plots into a single figure with a 2x2 layout for compactness
    combined_plot = plot(p1, p2, p3)  # Adjust size and layout
    return combined_plot
end

function plot_elongationrod(pd, step, semi)
    # Extract the time step for the title
    time_step = round(sol.t[step], digits = 5)
    pd = PlotData1D(sol.u[step], semi)

    # Create individual plots with increased linewidth
    p1 = plot(pd["V1"], xguide = "Spatial coordinate", yguide = "V1",
              xlims = (0, 2), ylims = (min_value, max_value),
              title = "V1 Time t=$time_step", linewidth = 2,
              titlefont_size = 10,  # Adjust title font size
              xguidefontsize = 8,  # Adjust x-axis label font size
              yguidefontsize = 8)  # Adjust y-axis tick label font size

    p2 = plot(pd["x2_1"], xguide = "Spatial coordinate", yguide = "Linear_Strain_1",
              xlims = (0, 2), ylims = (min_value, max_value),
              title = "x2_1 Time t=$time_step", linewidth = 2,
              titlefont_size = 10,  # Adjust title font size
              xguidefontsize = 8,  # Adjust x-axis label font size
              yguidefontsize = 8)  # Adjust y-axis tick label font size

    # Combine the plots into a single figure with a 2x2 layout for compactness
    combined_plot = plot(p1, p2)  # Adjust size and layout
    return combined_plot
end

# Calculate the dynamic skip value
total_steps = length(sol.u)
skip = max(1, div(total_steps, target_frames))
pd = PlotData1D(sol)

timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS")

# Conditional execution based on the Elongationrod flag
if Elongationrod
    filename_elongation_rod = "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\Linearstrains_plot_$(timestamp).gif"

    anim_sol_elongation_rod = @animate for step in 1:skip:total_steps
        plot_elongationrod(pd, step, semi)
    end

    gif(anim_sol_elongation_rod, filename_elongation_rod, fps = 10)
else
    filename_vel = "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\Linearvelocity_plot_$(timestamp).gif"
    filename_strain = "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\Linearstrains_plot_$(timestamp).gif"

    anim_sol_linear_velocity = @animate for step in 1:skip:total_steps
        plot_linearvelocities(pd, step, semi)
    end
    gif(anim_sol_linear_velocity, filename_vel, fps = 10)

    anim_sol_linear_strains = @animate for step in 1:skip:total_steps
        plot_linearstrains(pd, step, semi)
    end
    gif(anim_sol_linear_strains, filename_strain, fps = 10)
end

energy = Trixi.compute_energy(U, t, semi)

# plot discrete energy
p = plot(t, energy, xlabel = L"t", xguidefontsize = 20,
         xtickfontsize = 12,
         ytickfontsize = 12,
         ylabel = "Total discrete energy",
         yguidefontsize = 16,
         label = "Upwind flux",
         fmt = :pdf,
         color = RGBA(0, 0, 1, 1),
         linewidth = 2,
         size = (700, 500))
savefig("C:\\Users\\sund_si\\Desktop\\Solution_Plots\\discrete_energy\\total_discrete_energy_plot_$(timestamp).png")

##########

function plot_strains_in_interval(γ, interval_steps, x, num_frames)
    # Create a color gradient (e.g., blue to red) for the number of frames
    color_gradient = cgrad([:blue, :red], num_frames)

    # Maximum and minimum opacities for fading effect
    max_opacity = 0.9
    min_opacity = 0.3

    # Initialize the plot
    p = plot(size=(800, 600), dpi=400)

    # Loop through selected steps in the interval
    for (i, step) in enumerate(interval_steps)
        # Convert step to an integer index
        step_idx = Int(round(step))  # Ensures it can be used as an index

        # Calculate the color and opacity for this frame
        color = color_gradient[i]
        alpha_plot = max_opacity - (max_opacity - min_opacity) * (i / (num_frames - 1))
        time_step = round(sol.t[step_idx], digits = 5)

        # Add strain curve to the plot
        plot!(p, x, γ[1, :, step_idx],  # Assuming γ is a 3D array (strain data)
              #label = "t=$time_step",  # Optional legend label for each curve
              color = color,
              alpha = alpha_plot,
              linewidth = 3,
              legend=false
              )
    end

    # Customize the plot
    plot!(p,
          xlabel = "Spatial Coordinate (x)",
          ylabel = "Strain γ₁",
          ylims = (-0.5,0.5),
          #ylims = (-1,3.2),
          xlabelfontsize = 20,
          ylabelfontsize = 20,
          xtickfontsize = 16,
          ytickfontsize = 16,
          grid = true,
          alpha = 0.6,
          legend = false)  # Place legend outside the plot

    return p
end
## 273-521, 521-762, 762- 1115, 1115-1489, 1489-1836, 2557-2930 , 2930-3281, 3281-4005 down, 3281-3654 up
## 4 - 1115 , 2557- good
start_step = 53663
stop_step = 57257
interval_steps = range(start_step, stop_step, length=num_frames)

# Call the function
p = plot_strains_in_interval(γ, interval_steps, x, num_frames)

# Save the plot or display it
savefig(p, "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\discontinuousbeam\\strain_curves_6.pdf")
#display(p)

p = plot_beam_with_time_axis(r, interval_steps, x, num_frames)

savefig(p, "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\discontinuousbeam\\disc_beam_6.pdf")
#savefig(p, "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\discontinuousbeam\\disc_beam_5829s.pdf")
display(p)
# @save "elongation_rod_length_nodamp.jld2" x t U γ κ r r_fixed L energy dampingMat

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
# Assuming you already have the variables and functions defined

second_intervals = second_derivative_zeros[1:40]
function plot_beam_with_time_axis(r, interval_steps, x, num_frames)
    # Create a color gradient (e.g., blue to red) for the number of frames
    color_gradient = cgrad([:blue, :red], num_frames)

    # Maximum and minimum opacities for fading effect
    max_opacity = 0.9
    min_opacity = 0.3

    # Initialize the plot
    p = plot(size=(800, 600), dpi=400)

    # Vertical spacing between plots
    y_spacing = 0.5

    # Loop through selected steps in the interval
    for (i, step) in enumerate(interval_steps)
        # Convert step to an integer index
        step_idx = Int(round(step))  # Ensures it can be used as an index

        # Calculate the color and opacity for this frame
        color = color_gradient[i]
        alpha_plot = max_opacity - (max_opacity - min_opacity) * (i / (num_frames - 1))

        # Get the corresponding time for this step
        time_step = round(sol.t[step_idx], digits=5)

        # Calculate vertical offset for this plot
        y_offset = -y_spacing * i

        # Add the beam deformation curve to the plot with vertical offset
        plot!(p, r[1, :, step_idx], r[3, :, step_idx] .+ y_offset,  # Assuming r[3, :, :] is beam deformation along y-axis
              label = "t=$time_step",
              color = color,
              alpha = alpha_plot,
              linewidth = 3,
              legendfontsize = 20)
    end

    # Customize the overall plot
    plot!(p,
          xlabel = L"x_1",
          ylabel = "",       # Ensure no y-axis label
          yticks = false,    # Hide y-axis ticks
          #margin = 5Plots.mm,      # Adjust margins to prevent clipping
          xlims = (0, 1.8),
          xlabelfontsize = 20,
          xtickfontsize = 16,
          grid = true,
          alpha = 1.6,
          legend = :outerright,  # Place legend outside the plot
          ylims = (-y_spacing * (length(interval_steps) + 1), y_spacing))  # Adjust y-limits to fit all rows

    return p
end

start_step = 50000
stop_step = 60000
num_frames = 10
interval_steps = range(start_step, stop_step, length=num_frames)

# Call the function
p = plot_beam_with_time_axis(r, interval_steps, x, num_frames)

savefig(p, "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\discontinuousbeam\\disc_beam_1.pdf")
display(p)

function save_plots_for_intervals(γ, r, x, second_derivative_zeros, num_frames)
    # Iterate over adjacent values in second_derivative_zeros
    for i in 1:(length(second_derivative_zeros) - 1)
        # Define the interval using adjacent values
        start_step = second_derivative_zeros[i]
        stop_step = second_derivative_zeros[i + 1]

        # Generate interval steps for the current interval
        interval_steps = range(start_step, stop_step, length=num_frames)

        # Generate and save the strain plot
        p_strain = plot_strains_in_interval(γ, interval_steps, x, num_frames)
        savefig(p_strain, "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\discontinuousbeam\\strain_interval_$(i).pdf")

        # Generate and save the beam deformation plot
        p_beam = plot_beam_with_time_axis(r, interval_steps, x, num_frames)
        savefig(p_beam, "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\discontinuousbeam\\beam_interval_$(i).pdf")
    end
end

# Example usage
#save_plots_for_intervals(γ, r, x, second_intervals, num_frames)
