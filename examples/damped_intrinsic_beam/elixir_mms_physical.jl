using OrdinaryDiffEqStabilizedRK
using Trixi

# This manufactured solution verifies the physical Kelvin--Voigt model:
# its forcing has the block structure [distributed force/moment; 0].
beam_length = 1.0
wavenumber = pi / (2 * beam_length)
mms_direction = SVector(1.0, 2.0, 3.0, 4.0, 0.0, 0.0)

# Dense symmetric positive-definite matrices exercise anisotropic coupling.
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
flexibility_times_direction = flexibility_matrix \ mms_direction

@inline spatial_profile(x) = sin(wavenumber * x)
@inline spatial_profile_x(x) = wavenumber * cos(wavenumber * x)
@inline spatial_profile_xx(x) = -wavenumber^2 * sin(wavenumber * x)
@inline spatial_profile_xxx(x) = -wavenumber^3 * cos(wavenumber * x)

function manufactured_solution(x, t)
    u1 = t * spatial_profile(x) * mms_direction
    u2 = 0.5 * t^2 * spatial_profile_x(x) * flexibility_times_direction
    return SVector{12}(u1..., u2...)
end

function manufactured_solution_t(x, t)
    u1_t = spatial_profile(x) * mms_direction
    u2_t = t * spatial_profile_x(x) * flexibility_times_direction
    return SVector{12}(u1_t..., u2_t...)
end

function manufactured_solution_x(x, t)
    u1_x = t * spatial_profile_x(x) * mms_direction
    u2_x = 0.5 * t^2 * spatial_profile_xx(x) * flexibility_times_direction
    return SVector{12}(u1_x..., u2_x...)
end

function manufactured_solution_xx(x, t)
    u1_xx = t * spatial_profile_xx(x) * mms_direction
    u2_xx = 0.5 * t^2 * spatial_profile_xxx(x) * flexibility_times_direction
    return SVector{12}(u1_xx..., u2_xx...)
end

function initial_condition(x, t, equations::DampedIntrinsicBeamEquations1D)
    return manufactured_solution(x[1], t)
end

function physical_external_force(x, t, equations)
    return manufactured_force_damped_intrinsic_beam(x, t, equations_parabolic,
                                                    manufactured_solution,
                                                    manufactured_solution_t,
                                                    manufactured_solution_x,
                                                    manufactured_solution_xx)
end

equations_hyperbolic = DampedIntrinsicBeamEquations1D(;
                                                      mass_matrix,
                                                      flexibility_matrix,
                                                      damping_matrix,
                                                      initial_curvature = zeros(3),
                                                      external_force = physical_external_force)
equations_parabolic = DampedIntrinsicBeamDiffusion1D(equations_hyperbolic)

function verify_manufactured_solution()
    for x in range(0.0, beam_length; length = 5), t in (0.0, 0.25, 1.0)
        source = manufactured_source_damped_intrinsic_beam(x, t,
                                                           equations_parabolic,
                                                           manufactured_solution,
                                                           manufactured_solution_t,
                                                           manufactured_solution_x,
                                                           manufactured_solution_xx)
        @assert maximum(abs, source[7:12]) < 2.0e-12

        u = manufactured_solution(x, t)
        u_x = manufactured_solution_x(x, t)
        damping_resultant = intrinsic_beam_damping_resultant(u, u_x,
                                                             equations_parabolic)
        constitutive_resultant = damping_matrix *
                                 manufactured_solution_t(x, t)[7:12]
        @assert isapprox(damping_resultant, constitutive_resultant;
                         atol = 2.0e-12, rtol = 2.0e-12)
    end

    zero_field = (x, t) -> zero(SVector{12, Float64})
    inconsistent_time_derivative = (x, t) -> SVector{12}(zeros(6)..., ones(6)...)
    rejected = false
    try
        manufactured_force_damped_intrinsic_beam(0.0, 0.0,
                                                 equations_parabolic,
                                                 zero_field,
                                                 inconsistent_time_derivative,
                                                 zero_field, zero_field)
    catch error
        rejected = error isa ArgumentError
    end
    @assert rejected
    return nothing
end
verify_manufactured_solution()

sigma = 1.0
@inline function flux_mms(u_ll, u_rr, orientation,
                          equations::DampedIntrinsicBeamEquations1D)
    central_flux = 0.5 * (flux(u_ll, orientation, equations) +
                    flux(u_rr, orientation, equations))
    dissipation = 0.5 * sigma * equations.propagation_matrix_abs *
                  (u_ll - u_rr)
    return central_flux + dissipation
end

boundary_condition = BoundaryConditionDampedIntrinsicBeam()
boundary_conditions = (boundary_condition, boundary_condition)

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
