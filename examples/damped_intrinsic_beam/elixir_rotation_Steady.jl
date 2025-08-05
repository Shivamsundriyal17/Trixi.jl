#Rotational Beam Experiment from Thesis
#Uncomment lines in boundary conditions and velocity function based on wether starting from rest
#or steady state solution
using Trixi          # Trixi framework for DG simulations
using OrdinaryDiffEq # Differential equation solvers
using LinearAlgebra  # Linear algebra operations
using Plots          # Plotting functionality
using Printf         # Formatted output
using Dates          # Date utilities
using LaTeXStrings   # LaTeX support for text and labels
using ColorSchemes

# Define initial conditions
function initial_condition_zero_full(x, t, ::Damped_Full_Hyperbolic)
    return SVector{12}(zeros(12))
end

# Define flexibility, mass matrix, and damping matrix
flexMat = inv(diagm([1.0e3, 1.0e3, 1.0e3, 500.0, 500.0, 500.0]))
massMat = diagm([1.0, 1.0, 1.0, 20.0, 10.0, 10.0])
dampingMat = 0 * inv(diagm([1.0e3, 1.0e3, 1.0e3, 10.0, 10.0, 10.0]))

# Define Constant
L = 4.0
constAngSpeed = (pi / 2 - 1.21) / (sqrt(massMat[2, 2] * flexMat[1, 1]) * L)
k = (flexMat[1, 1] * cos(constAngSpeed * sqrt(massMat[2, 2] * flexMat[1, 1]) * L))^(-1)

# Steady State Solution\Initial Condition
function steady_state_initial_condition(x, t, ::Damped_Full_Hyperbolic)
    return SVector{12}([
                           0.0,  # v1
                           sqrt(flexMat[1, 1] / massMat[2, 2]) * k *
                           sin(constAngSpeed * sqrt(massMat[2, 2] * flexMat[1, 1]) * x[1]),  # v2
                           0.0, 0.0, 0.0,  # v3, omega1, omega2
                           constAngSpeed,  # omega3
                           k *
                           cos(constAngSpeed * sqrt(massMat[2, 2] * flexMat[1, 1]) * x[1]) -
                           1 / flexMat[1, 1],  # linear_Strain_1
                           0.0, 0.0, 0.0, 0.0, 0.0  # Remaining linear and angular strains
                       ])
end

# Initial curvature and external forces/moments
k_init = [0.0 0.0 0.0]
f_ext_func(x, t) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Set up hyperbolic and parabolic equations
equations_hyperbolic = Damped_Full_Hyperbolic(flexMat, massMat, dampingMat, k_init,
                                              f_ext_func)
equations_parabolic = Damped_Full_Parabolic(equations_hyperbolic,
                                            dampingMat, flexMat, massMat,
                                            equations_hyperbolic.Gamma_inv,
                                            equations_hyperbolic.Psi,
                                            k_init)

# Define velocity and strain functions

##Use this function when starting in the rest condition
# function velocities(t)
#     velocity_value = t <= 1.0 ? constAngSpeed * (t / 1.0) : constAngSpeed
#     return [zeros(5); velocity_value]
# end

#Use this when starting with the steady state solution
velocities(t) = [zeros(5); constAngSpeed]

strains(t) = zeros(6)

####################

# Define boundary conditions
boundary_conditions_hyperbolic = (u_inner, orientation, direction,
x, t, surface_flux_function,
equations) -> boundary_conditions_ib(u_inner, orientation, direction,
                                     x, t, surface_flux_function,
                                     velocities, strains,
                                     equations)

@inline function (boundary_condition_damped)(flux_inner, u_inner, normal::Int64,
                                             direction::Int64,
                                             x, t, operator_type::Trixi.Divergence,
                                             equations_parabolic::Damped_Full_Parabolic)
    if direction == 1
        return [flux_inner[1:6]; zeros(6)]

    elseif direction == 2
        return [zeros(6); flux_inner[7:12]]
    end
end

@inline function (boundary_condition_damped)(flux_inner, u_inner, normal::Any,
                                             direction::Int64,
                                             x, t, operator_type::Trixi.Gradient,
                                             equations_parabolic::Damped_Full_Parabolic)
    u_boundary = zeros(12)
    ##Use this function when starting in the rest condition
    # function velocities(t)
    #     velocity_value = t <= 1.0 ? constAngSpeed * (t / 1.0) : constAngSpeed
    #     return [zeros(5); velocity_value]
    # end
    #Use this when starting with the steady state solution
    velocities(t) = [zeros(5); constAngSpeed]
    if direction == 1
        y_B = zeros(6)  # Boundary condition vector for velocity (left boundary)
        #u_boundary = [y_B; u_inner[7:12]]
        u_boundary = [velocities(t); u_inner[7:12]]

    elseif direction == 2
        s_B = zeros(6)  # Boundary condition vector for strain (right boundary)
        u_boundary = [u_inner[1:6]; u_inner[7:12]]
    end

    return u_boundary
end

# Define Parabolic boundary conditions
boundary_conditions_parabolic = boundary_condition_damped

# Define the Discontinuous Galerkin Spectral Element Method (DGSEM) solver
solver = DGSEM(polydeg = 3,                 # Polynomial degree for DGSEM
               surface_flux = flux_upwind)

# Set coordinates for the mesh
coordinates_min = (0.0,)  # Minimum coordinates (e.g., min(x), min(y))
coordinates_max = (4.0,)  # Maximum coordinates (e.g., max(x), max(y))

# Create the mesh using TreeMesh
mesh = TreeMesh(coordinates_min,                      # Minimum coordinates of the domain
                coordinates_max,                      # Maximum coordinates of the domain
                initial_refinement_level = 3,           # Initial level of mesh refinement
                n_cells_max = 10_000,                   # Maximum number of cells in the mesh
                periodicity = false)

#initial_condition = initial_condition_zero_full
initial_condition = initial_condition_zero_full

# Set up semidiscretization with hyperbolic and parabolic equations
semi = SemidiscretizationHyperbolicParabolic(mesh,
                                             (equations_hyperbolic, equations_parabolic),
                                             initial_condition,
                                             solver;
                                             boundary_conditions = (boundary_conditions_hyperbolic,
                                                                    boundary_conditions_parabolic))
# Define time span for the simulation
tspan = (0.0, 40.0)

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

#################################################
# Postprocessing
#################################################

"""
    calc_rot_mat_y(U, t)

Calculates rotation matrix for rotation around y-axis by integrating the
associated angular velocities in U
"""
function calc_rot_mat_y_damped(U, t)
    Nt = length(t)
    R = Array{Float64}(undef, 3, 3, Nt)
    R[:, :, 1] = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    alpha = zeros(Nt)
    for j in 1:(Nt - 1)
        alpha[j + 1] = alpha[j] + (t[j + 1] - t[j]) * U[5, 1, j]
        R[:, :, j + 1] = [cos(alpha[j + 1]) 0 -sin(alpha[j + 1]); 0 1 0;
                          sin(alpha[j + 1]) 0 cos(alpha[j + 1])]
    end

    return R, alpha
end

"""
    calc_rot_mat_z(U, t)

    Calculates rotation matrix for rotation around z-axis by integrating the
        associated angular velocities in U
"""
function calc_rot_mat_z_damped(U, t)
    Nt = length(t)
    R = Array{Float64}(undef, 3, 3, Nt)
    R[:, :, 1] = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    alpha = zeros(Nt)
    for j in 1:(Nt - 1)
        alpha[j + 1] = alpha[j] + (t[j + 1] - t[j]) * U[6, 1, j]
        R[:, :, j + 1] = [cos(alpha[j + 1]) sin(alpha[j + 1]) 0;
                          -sin(alpha[j + 1]) cos(alpha[j + 1]) 0; 0 0 1]
    end

    return R, alpha
end

timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")

# spatial nodes
x = semi.cache.elements._node_coordinates
Nx = length(x)

# time nodes
t = sol.t
Nt = length(t)

num_frames = 200
step_size = max(1, Int(floor(Nt / num_frames)))

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

function plot_strain(γ, step, x)
    p = plot(x, γ[1, :, step],
             ylim = (minimum(γ[1, :, :]), maximum(γ[1, :, :])),
             legend = true,
             label = "Strain γ₁",
             title = "Strain at t=$(round(t[step], digits=5))",
             xlabel = "Spatial Coordinate (x)",
             ylabel = "Strain γ₁",
             linewidth = 2,
             color = :red)
    return p
end

S0 = diagm([1.0; 1.0; 1.0])
k_func(x) = [0.0; 0.0; 0.0]
S, alpha = calc_rot_mat_y_damped(U, t)
S, alpha = calc_rot_mat_z_damped(U, t)

# determine position line of beam
r = Array{Float64}(undef, 3, Nx, Nt)
r_fixed = Array{Float64}(undef, 3, Nx, Nt)
for j in 1:Nt
    r_fixed[:, :, j] = calc_deformation_trapez(S[:, :, j], x, k_func, κ[:, :, j],
                                               γ[:, :, j])
    r[:, :, j] = calc_deformation_trapez(S0, x, k_func, κ[:, :, j], γ[:, :, j])
end

# determine beam's approximate length at every time step
L = zeros(Nt)
for n in 1:Nt
    diffs = diff(r[:, :, n], dims = 2)  # Compute differences between consecutive points
    L[n] = sum(sqrt.(sum(diffs .^ 2, dims = 1)))  # Sum Euclidean distances
end

# anim = @animate for n in 1:step_size:Nt
#     plot(r[1, :, n], r[3, :, n], xlims = (-5, 5), ylims = (-2.5, 2.5),
#          aspect_ratio = :equal, legend = false,
#          title = string("t = ", round(t[n], digits = 2)))
# end

# gif(anim, "C:\\Users\\sund_si\\Desktop\\testrotationdamped$(timestamp).gif", fps = 100)
# Function to plot the beam deformation and node positions
# determine position line of beam

target_frames = 200
x_min = minimum(r[1, :, :])
x_max = maximum(r[1, :, :])
x_limits = (x_min, x_max + 0.5)
plot_frequency = max(1, div(Nt, target_frames))
# Calculate the dynamic skip value
total_steps = length(sol.u)
skip = max(1, div(total_steps, target_frames))
pd = PlotData1D(sol)

# plot rotating beam 1:skip:total_steps
anim = @animate for n in 1:skip:total_steps
    plot(r_fixed[1, :, n], r_fixed[2, :, n],
         xlims = (-6, 6),
         ylims = (-6, 6),
         linewidth = 3,
         aspect_ratio = :equal,
         legend = false,
         title = string("t = ", round(t[n], digits = 2)))
end
gif(anim, "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\moving_beam$(timestamp).gif",
    fps = 100)

anim = @animate for n in 1:skip:total_steps
    plot(r[1, :, n], r[2, :, n],
         xlims = (-6, 6),
         ylims = (-6, 6),
         linewidth = 3,
         aspect_ratio = :equal,
         legend = false,
         title = string("t = ", round(t[n], digits = 2)))
end
gif(anim, "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\moving_beamcam$(timestamp).gif",
    fps = 100)

# Function to plot the beam deformation and node positions
function plot_deformation_with_nodes(r, step, x, γ)
    p = plot(r[1, :, step], r[3, :, step], line_z = γ[1, :, step],
             xlims = (-3.0, 3.0), ylims = (-3.0, 3.0),
             clims = (min_value, max_value),
             xlabel = L"x_1", ylabel = L"x_3",
             aspect_ratio = :equal, linewidth = 3,
             legend = false, cbar = true,
             title = string("Deformation at t = ", round(t[step], digits = 2)))

    # Add smaller markers for the first and last nodes
    scatter!(p, [r[1, 1, step]], [r[3, 1, step]], label = "Start Node", color = :green,
             marker = :circle, markersize = 4)
    scatter!(p, [r[1, end, step]], [r[3, end, step]], label = "End Node", color = :blue,
             marker = :circle, markersize = 4)

    return p
end

# Function to plot the beam deformation and node positions
function plot_deformation_with_nodes_z(r, step, x, γ)
    p = plot(r[1, :, step], r[3, :, step], line_z = γ[3, :, step],
             xlims = (-3.0, 3.0), ylims = (-3.0, 3.0),
             clims = (global_min, global_max),
             xlabel = L"x_1", ylabel = L"x_3",
             aspect_ratio = :equal, linewidth = 3,
             legend = false, cbar = true,
             title = string("Deformation at t = ", round(t[step], digits = 2)))

    # Add smaller markers for the first and last nodes
    scatter!(p, [r[1, 1, step]], [r[3, 1, step]], label = "Start Node", color = :green,
             marker = :circle, markersize = 4)
    scatter!(p, [r[1, end, step]], [r[3, end, step]], label = "End Node", color = :blue,
             marker = :circle, markersize = 4)

    return p
end

function plot_linearvelocities(pd, step, semi)
    # Extract the time step for the title
    time_step = round(sol.t[step], digits = 5)
    pd = PlotData1D(sol.u[step], semi)

    # Create individual plots with increased linewidth
    p1 = plot(pd["V1"], xguide = "Spatial coordinate", yguide = "V1",
              xlims = (0, 2), ylims = (-0.25, 0.25),
              title = "V1 Time t=$time_step", linewidth = 2,
              titlefont_size = 10,  # Adjust title font size
              xguidefontsize = 8,  # Adjust x-axis label font size
              yguidefontsize = 8)  # Adjust y-axis tick label font size

    p2 = plot(pd["V2"], xguide = "Spatial coordinate", yguide = "V2",
              xlims = (0, 2), ylims = (-0.0, 5.0),
              title = "V2", linewidth = 2,
              titlefont_size = 10,  # Adjust title font size
              xguidefontsize = 8,  # Adjust x-axis label font size
              yguidefontsize = 8)  # Adjust y-axis tick label font size

    p3 = plot(pd["V3"], xguide = "Spatial coordinate", yguide = "V3",
              xlims = (0, 2),
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
              xlims = (0, 2), ylim = (minimum(γ[1, :, :]), maximum(γ[1, :, :])),
              title = "x2_1 Time t=$time_step", linewidth = 2,
              titlefont_size = 10,  # Adjust title font size
              xguidefontsize = 8,  # Adjust x-axis label font size
              yguidefontsize = 8)  # Adjust y-axis tick label font size

    p2 = plot(pd["x2_2"], xguide = "Spatial coordinate", yguide = "Linear_Strain_2",
              xlims = (0, 2), ylim = (minimum(γ[2, :, :]), maximum(γ[2, :, :])),
              title = "x2_2", linewidth = 2,
              titlefont_size = 10,  # Adjust title font size
              xguidefontsize = 8,  # Adjust x-axis label font size
              yguidefontsize = 8)  # Adjust y-axis tick label font size

    p3 = plot(pd["x2_3"], xguide = "Spatial coordinate", yguide = "Linear_Strain_3",
              xlims = (0, 2), ylim = (minimum(γ[3, :, :]), maximum(γ[3, :, :])),
              title = "x2_3", linewidth = 2,
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

# filename_vel = "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\Linearvelocity_plot_$(timestamp).gif"
# filename_strain = "C:\\Users\\sund_si\\Desktop\\Solution_Plots\\Linearstrains_plot_$(timestamp).gif"

# anim_sol_linear_velocity = @animate for step in 1:skip:total_steps
#     plot_linearvelocities(pd, step, semi)
# end
# gif(anim_sol_linear_velocity, filename_vel, fps=10)

# anim_sol_linear_strains = @animate for step in 1:skip:total_steps
#     plot_linearstrains(pd, step, semi)
# end
# gif(anim_sol_linear_strains, filename_strain, fps=10)

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

function find_second_derivative_zeros(data, times)
    first_derivative = diff(data) ./ diff(times)
    second_derivative = diff(first_derivative) ./ diff(times[1:(end - 1)])
    zero_crossings = Int[]  # To store indices
    for i in 1:(length(second_derivative) - 1)
        if second_derivative[i] * second_derivative[i + 1] < 0
            push!(zero_crossings, i + 1)  # +1 to account for second_derivative indexing
        end
    end

    return zero_crossings
end

function find_local_maxima(data)
    maxima = Int[]
    for i in 2:(length(data) - 1)  # Ensure correct range
        if data[i] > data[i - 1] && data[i] > data[i + 1]
            push!(maxima, i)
        end
    end
    return maxima
end

maxima_indices = find_local_maxima(energy)
local_maxima = energy[maxima_indices]
local_maxima_times = t[maxima_indices]

# Apply the function
second_derivative_zeros = find_second_derivative_zeros(energy, t)
zero_crossing_times = t[second_derivative_zeros]
zero_crossing_energies = energy[second_derivative_zeros]
max_value = 1650.8176079116051
#max_value = maximum(energy) + 100
p = plot(t, energy,
         xlabel = "t (seconds)",                 # Proper LaTeX-style label
         ylabel = "Total Discrete Energy",          # Label for Y-axis with text
         xguidefontsize = 20,                           # Font size for X-axis label
         yguidefontsize = 20,                           # Font size for Y-axis label
         xtickfontsize = 14,                            # Font size for X-axis ticks
         ytickfontsize = 14,                            # Font size for Y-axis ticks
         grid = true,                                   # Add grid for better readability
         gridalpha = 0.2,                               # Set grid transparency
         color = RGBA(0, 0, 1, 1),                      # Line color (blue)
         linewidth = 2,                                 # Line thickness
         framestyle = :box,                             # Box-style frame for clarity
         legend = :bottomright,                         # Position legend at bottom right
         legendfontsize = 14,                           # Font size for legend
         label = "Discrete Energy",                     # Legend label for the plot
         size = (800, 600), xlims = (0, 40), ylims = (0, max_value))         # Custom Y-axis limits
zero_crossing_times = t[second_derivative_zeros]
zero_crossing_energies = energy[second_derivative_zeros]
#scatter!(zero_crossing_times, zero_crossing_energies, color=:blue, label="Inflection Point")

# Save as a high-quality PDF
savefig("C:\\Users\\sund_si\\Desktop\\Solution_Plots\\discrete_energy\\total_discrete_energy_nodamp$(timestamp).pdf")