if !isdefined(@__MODULE__, :BeamRunHelpers)
    Base.include(@__MODULE__, joinpath(@__DIR__, "beam_run_helpers.jl"))
end
using .BeamRunHelpers

using OrdinaryDiffEqStabilizedRK
using Trixi

# A fully active physical Kelvin--Voigt MMS. Unlike elixir_mms_physical.jl,
# this case exercises the geometric and nonlinear terms in the lower
# compatibility equation and uses exact inhomogeneous split boundary data.
beam_length = 1.0
mms_lambda = parse(Float64, get(ENV, "RICH_MMS_LAMBDA", "2.0"))
mms_profile = get(ENV, "RICH_MMS_PROFILE", "exp1")
mms_profile in ("quartic", "exp1", "exp2", "poly7") ||
    throw(ArgumentError("unknown RICH_MMS_PROFILE='$mms_profile'; " *
                        "use quartic, exp1, exp2, or poly7"))
mms_growth_rate = mms_lambda - 1.0
mms_direction = SVector(3.0, 2.0, 3.0, 1.0, 1.0, 1.0)
mms_geometry_action = SVector(0.0, 1.0, -1.0, 0.0, 0.0, 0.0)
mms_axial_shift = SVector(1.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# Use exactly the same dense coefficient fixture as the primary MMS.
flexibility_matrix = [1.0 0.3 0.3 0.12 0.05 0.025;
                      0.3 2.0 0.47 0.2 0.12 0.05;
                      0.3 0.47 3.0 0.35 0.2 0.12;
                      0.12 0.2 0.35 4.0 0.34 0.2;
                      0.05 0.12 0.2 0.34 5.0 0.34;
                      0.025 0.05 0.12 0.2 0.34 6.0]
mass_matrix = [1.84147 0.25403 0.34207 0.11081 0.05841 0.02270;
               0.25403 2.84147 0.43782 0.22524 0.11081 0.05841;
               0.34207 0.43782 3.84147 0.32702 0.22524 0.11081;
               0.11081 0.22524 0.32702 4.84147 0.32161 0.22524;
               0.05841 0.11081 0.22524 0.32161 5.84147 0.32161;
               0.02270 0.05841 0.11081 0.22524 0.32161 6.84147]
damping_matrix = 0.01 * flexibility_matrix
flexibility_inverse = inv(flexibility_matrix)

if !isdefined(@__MODULE__, :RichMMSFunctions)
    Base.include(@__MODULE__, joinpath(@__DIR__, "rich_mms_functions.jl"))
end
using .RichMMSFunctions

mms_parameters = RichBeamMMS(mms_lambda, flexibility_inverse, mms_profile)
manufactured_solution = RichField(mms_parameters, :value)
manufactured_solution_t = RichField(mms_parameters, :time)
manufactured_solution_x = RichField(mms_parameters, :space)
manufactured_solution_xx = RichField(mms_parameters, :space2)
physical_external_force = RichForce(mms_parameters)

function initial_condition(x, t, equations::DampedIntrinsicBeamEquations1D)
    return RichField(equations.external_force.parameters, :value)(x[1], t)
end

equations_hyperbolic = DampedIntrinsicBeamEquations1D(;
                                                      mass_matrix,
                                                      flexibility_matrix,
                                                      damping_matrix,
                                                      initial_curvature = zeros(3),
                                                      external_force = physical_external_force)
equations_parabolic = DampedIntrinsicBeamDiffusion1D(equations_hyperbolic)

@inline mms_coordinate(x::Real) = x
@inline mms_coordinate(x) = first(x)

function exact_left_velocity(x, t, equations)
    u = RichField(equations.external_force.parameters, :value)(mms_coordinate(x), t)
    return SVector{6}(u[1:6])
end

function exact_right_resultant(x, t, equations)
    u = RichField(equations.external_force.parameters, :value)(mms_coordinate(x), t)
    return SVector{6}(u[7:12])
end

function exact_right_damping_resultant(x, t, equations)
    u_t = RichField(equations.external_force.parameters, :time)(mms_coordinate(x), t)
    return equations.damping_matrix * SVector{6}(u_t[7:12])
end

function verify_manufactured_solution()
    maximum_source_lower_residual = 0.0
    maximum_lower_residual = 0.0
    maximum_normalized_source_lower_residual = 0.0
    maximum_normalized_lower_residual = 0.0
    maximum_constitutive_residual = 0.0
    minimum_geometry_activity = Inf
    minimum_nonlinear_activity = Inf
    minimum_u1_component = Inf
    minimum_u2_component = Inf
    minimum_r_tau_component = Inf
    for x in range(0.0, beam_length; length = 9),
        t in range(0.0, 1.0; length = 9)

        source = manufactured_source_damped_intrinsic_beam(x, t,
                                                           equations_parabolic,
                                                           manufactured_solution,
                                                           manufactured_solution_t,
                                                           manufactured_solution_x,
                                                           manufactured_solution_xx)
        u = manufactured_solution(x, t)
        u_t = manufactured_solution_t(x, t)
        u_x = manufactured_solution_x(x, t)
        u1 = SVector{6}(u[1:6])
        u2 = SVector{6}(u[7:12])
        elastic_strain = equations_hyperbolic.flexibility_matrix * u2
        geometry_term = transpose(equations_hyperbolic.geometry_matrix) * u1
        nonlinear_term = transpose(intrinsic_beam_l1(u1)) * elastic_strain
        cancellation_scale = max(1.0,
                                 maximum(abs, SVector{6}(u_x[1:6])),
                                 maximum(abs, geometry_term),
                                 6.0 * maximum(abs, u1) *
                                 maximum(abs, elastic_strain))
        source_cancellation_scale = max(1.0,
                                        6.0 * maximum(abs, flexibility_inverse) *
                                        cancellation_scale)
        compatibility_residual = equations_hyperbolic.flexibility_matrix *
                                 SVector{6}(u_t[7:12]) -
                                 (SVector{6}(u_x[1:6]) - geometry_term +
                                  nonlinear_term)
        source_lower_residual = maximum(abs, source[7:12])
        normalized_source_lower_residual = source_lower_residual /
                                           source_cancellation_scale
        normalized_lower_residual = maximum(abs, compatibility_residual) /
                                    cancellation_scale
        maximum_source_lower_residual = max(maximum_source_lower_residual,
                                            source_lower_residual)
        maximum_lower_residual = max(maximum_lower_residual,
                                     maximum(abs, compatibility_residual))
        maximum_normalized_source_lower_residual = max(maximum_normalized_source_lower_residual,
                                                       normalized_source_lower_residual)
        maximum_normalized_lower_residual = max(maximum_normalized_lower_residual,
                                                normalized_lower_residual)
        minimum_geometry_activity = min(minimum_geometry_activity,
                                        maximum(abs, geometry_term))
        minimum_nonlinear_activity = min(minimum_nonlinear_activity,
                                         maximum(abs, nonlinear_term))
        @assert normalized_source_lower_residual < 100.0 * eps(Float64)
        @assert normalized_lower_residual < 100.0 * eps(Float64)
        @assert maximum(abs, geometry_term) > 0.0
        @assert maximum(abs, nonlinear_term) > 0.0
        @assert isapprox(nonlinear_term, mms_lambda * geometry_term;
                         atol = 5.0e-12, rtol = 5.0e-12)

        damping_resultant = intrinsic_beam_damping_resultant(u, u_x,
                                                             equations_parabolic)
        constitutive_resultant = damping_matrix * SVector{6}(u_t[7:12])
        maximum_constitutive_residual = max(maximum_constitutive_residual,
                                            maximum(abs,
                                                    damping_resultant -
                                                    constitutive_resultant))
        minimum_u1_component = min(minimum_u1_component, minimum(u1))
        minimum_u2_component = min(minimum_u2_component, minimum(u2))
        minimum_r_tau_component = min(minimum_r_tau_component,
                                      minimum(damping_resultant))
        @assert isapprox(damping_resultant, constitutive_resultant;
                         atol = 5.0e-12, rtol = 5.0e-12)
        @assert minimum(u1) > 0.0
        @assert minimum(u2) > 0.0
        @assert minimum(damping_resultant) > 0.0
    end
    return (; maximum_source_lower_residual, maximum_lower_residual,
            maximum_normalized_source_lower_residual,
            maximum_normalized_lower_residual,
            maximum_constitutive_residual,
            minimum_geometry_activity, minimum_nonlinear_activity,
            minimum_u1_component, minimum_u2_component,
            minimum_r_tau_component)
end
rich_mms_exact_audit = verify_manufactured_solution()

sigma = 1.0
flux_mms = RichFlux(sigma)

boundary_condition = BoundaryConditionDampedIntrinsicBeam(;
                                                          left_velocity = exact_left_velocity,
                                                          right_resultant = exact_right_resultant,
                                                          right_damping_resultant = exact_right_damping_resultant)
boundary_conditions = (boundary_condition, boundary_condition)

function verify_boundary_operators()
    maximum_left_hyperbolic_residual = 0.0
    maximum_right_hyperbolic_residual = 0.0
    maximum_left_gradient_residual = 0.0
    maximum_right_gradient_residual = 0.0
    maximum_left_divergence_residual = 0.0
    maximum_right_divergence_residual = 0.0
    maximum_split_target_residual = 0.0
    maximum_total_resultant_residual = 0.0
    orientation = 1
    for t in range(0.0, 1.0; length = 9),
        (x, direction) in ((0.0, 1), (beam_length, 2))

        coordinate = SVector(x)
        u = manufactured_solution(x, t)
        u_x = manufactured_solution_x(x, t)

        exact_hyperbolic_flux = flux_mms(u, u, orientation,
                                         equations_hyperbolic)
        boundary_hyperbolic_flux = boundary_condition(u, orientation,
                                                      direction, coordinate, t,
                                                      flux_mms,
                                                      equations_hyperbolic)
        hyperbolic_residual = maximum(abs, boundary_hyperbolic_flux -
                                           exact_hyperbolic_flux)
        if direction == 1
            maximum_left_hyperbolic_residual = max(maximum_left_hyperbolic_residual,
                                                   hyperbolic_residual)
        else
            maximum_right_hyperbolic_residual = max(maximum_right_hyperbolic_residual,
                                                    hyperbolic_residual)
        end
        @assert isapprox(boundary_hyperbolic_flux, exact_hyperbolic_flux;
                         atol = 5.0e-12, rtol = 5.0e-12)

        boundary_gradient_trace = boundary_condition(u, u, orientation,
                                                     direction, coordinate, t,
                                                     Trixi.Gradient(),
                                                     equations_parabolic)
        gradient_residual = maximum(abs, boundary_gradient_trace - u)
        if direction == 1
            maximum_left_gradient_residual = max(maximum_left_gradient_residual,
                                                 gradient_residual)
        else
            maximum_right_gradient_residual = max(maximum_right_gradient_residual,
                                                  gradient_residual)
        end
        @assert isapprox(boundary_gradient_trace, u;
                         atol = 5.0e-12, rtol = 5.0e-12)

        exact_parabolic_flux = flux(u, u_x, orientation,
                                    equations_parabolic)
        boundary_parabolic_flux = boundary_condition(exact_parabolic_flux,
                                                     nothing, orientation,
                                                     direction, coordinate, t,
                                                     Trixi.Divergence(),
                                                     equations_parabolic)
        divergence_residual = maximum(abs, boundary_parabolic_flux -
                                           exact_parabolic_flux)
        if direction == 1
            maximum_left_divergence_residual = max(maximum_left_divergence_residual,
                                                   divergence_residual)
        else
            maximum_right_divergence_residual = max(maximum_right_divergence_residual,
                                                    divergence_residual)
        end
        @assert isapprox(boundary_parabolic_flux, exact_parabolic_flux;
                         atol = 5.0e-12, rtol = 5.0e-12)
    end

    for t in range(0.0, 1.0; length = 9)
        coordinate = SVector(beam_length)
        u = manufactured_solution(beam_length, t)
        u2 = SVector{6}(u[7:12])
        r_tau = exact_right_damping_resultant(coordinate, t,
                                              equations_hyperbolic)
        split_target_residual = max(maximum(abs,
                                            exact_right_resultant(coordinate, t,
                                                                  equations_hyperbolic) -
                                            u2),
                                    maximum(abs,
                                            exact_right_damping_resultant(coordinate, t,
                                                                          equations_hyperbolic) -
                                            r_tau))
        maximum_split_target_residual = max(maximum_split_target_residual,
                                            split_target_residual)
        prescribed_total_resultant = exact_right_resultant(coordinate, t,
                                                           equations_hyperbolic) +
                                     r_tau
        maximum_total_resultant_residual = max(maximum_total_resultant_residual,
                                               maximum(abs,
                                                       prescribed_total_resultant -
                                                       (u2 + r_tau)))
        @assert isapprox(prescribed_total_resultant, u2 + r_tau;
                         atol = 5.0e-12, rtol = 5.0e-12)
    end
    return (; maximum_left_hyperbolic_residual,
            maximum_right_hyperbolic_residual,
            maximum_left_gradient_residual,
            maximum_right_gradient_residual,
            maximum_left_divergence_residual,
            maximum_right_divergence_residual,
            maximum_split_target_residual,
            maximum_total_resultant_residual)
end
rich_mms_boundary_audit = verify_boundary_operators()

polydeg = 3
solver = DGSEM(polydeg = polydeg, surface_flux = flux_mms)
solver_parabolic = ViscousFormulationLocalDG()
initial_refinement_level = 3
mesh = TreeMesh((0.0,), (beam_length,);
                initial_refinement_level,
                n_cells_max = 10_000,
                periodicity = false)

semi = SemidiscretizationHyperbolicParabolic(mesh,
                                             (equations_hyperbolic,
                                              equations_parabolic),
                                             initial_condition, solver;
                                             solver_parabolic,
                                             boundary_conditions)

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)
analysis_callback = AnalysisCallback(semi, interval = 0)

time_int_tol = 1.0e-12
sol = solve(ode, ROCK4();
            abstol = time_int_tol,
            reltol = time_int_tol,
            ode_default_options()...)

require_complete_solution(sol, last(tspan))
