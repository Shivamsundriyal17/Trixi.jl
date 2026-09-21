module TestDampedIntrinsicBeamRegressions

using Test
using LinearAlgebra
using Trixi
using SciMLBase: ReturnCode
using OrdinaryDiffEqStabilizedRK

include(joinpath(@__DIR__, "..", "examples", "damped_intrinsic_beam",
                 "beam_run_helpers.jl"))
using .BeamRunHelpers

function core_allocations(parabolic)
    u = SVector{12}(ntuple(i -> 0.01 * i, Val(12)))
    gradient = 0.1 * u
    x = SVector(0.2)
    intrinsic_beam_gradient_source(u, gradient, x, 0.0, parabolic)
    return @allocated intrinsic_beam_gradient_source(u, gradient, x, 0.0, parabolic)
end

function periodic_beam(damping; polydeg = 2, level = 2, periodic = true)
    mass = Matrix(Diagonal(1.0:6.0))
    flexibility = Matrix(Diagonal(2.0:7.0))
    equations = DampedIntrinsicBeamEquations1D(; mass_matrix = mass,
                                               flexibility_matrix = flexibility,
                                               damping_matrix = damping * flexibility,
                                               initial_curvature = [0.02, -0.03, 0.01])
    parabolic = DampedIntrinsicBeamDiffusion1D(equations)
    initial_condition(x, t, equations) = SVector{12}(ntuple(i -> 0.02 *
                                                                 sinpi(2 * x[1] + i / 7),
                                                            Val(12)))
    mesh = TreeMesh((0.0,), (1.0,); initial_refinement_level = level,
                    n_cells_max = 100, periodicity = periodic)
    solver = DGSEM(; polydeg, surface_flux = flux_upwind)
    semi = SemidiscretizationHyperbolicParabolic(mesh, (equations, parabolic),
                                                 initial_condition, solver;
                                                 solver_parabolic = ViscousFormulationLocalDG(),
                                                 boundary_conditions = periodic ?
                                                                       (boundary_condition_periodic,
                                                                        boundary_condition_periodic) :
                                                                       (BoundaryConditionDampedIntrinsicBeam(),
                                                                        BoundaryConditionDampedIntrinsicBeam()))
    return semidiscretize(semi, (0.0, 0.01))
end

function energy_residual(ode; periodic = true)
    semi = ode.p
    equations = semi.equations
    state = copy(ode.u0)
    hyperbolic_rhs = similar(state)
    parabolic_rhs = similar(state)
    Trixi.rhs!(hyperbolic_rhs, state, semi, 0.0)
    Trixi.rhs_parabolic!(parabolic_rhs, state, semi, 0.0)
    shape = (12, Trixi.nnodes(semi.solver), length(semi.cache.elements.inverse_jacobian))
    u = reshape(state, shape)
    du = reshape(hyperbolic_rhs + parabolic_rhs, shape)
    gradient = semi.cache_parabolic.viscous_container.gradients
    weights = semi.solver.basis.weights
    energy_rate = 0.0
    material = 0.0
    for element in axes(u, 3), i in axes(u, 2)
        value = beam_node(u, i, element)
        rate = beam_node(du, i, element)
        derivative = beam_node(gradient, i, element)
        u1 = SVector{6}(value[1:6])
        u2 = SVector{6}(value[7:12])
        argument = SVector{6}(derivative[1:6]) - transpose(equations.geometry_matrix) * u1 +
                   transpose(intrinsic_beam_l1(u1)) * (equations.flexibility_matrix * u2)
        weight = weights[i] / abs(semi.cache.elements.inverse_jacobian[element])
        energy_rate += weight * dot(value, equations.capacity_matrix * rate)
        material += weight * dot(argument, equations.damping_operator * argument)
    end
    order = sortperm(vec(semi.cache.elements.node_coordinates[1, 1, :]))
    jumps = 0.0
    for index in 1:(length(order) - (periodic ? 0 : 1))
        left = order[index]
        right = order[mod1(index + 1, length(order))]
        jump = beam_node(u, size(u, 2), left) - beam_node(u, 1, right)
        jumps += 0.5 * dot(jump,
                     equations.capacity_matrix * equations.propagation_matrix_abs * jump)
    end
    boundary = 0.0
    if !periodic
        left = beam_node(u, 1, first(order))
        right = beam_node(u, size(u, 2), last(order))
        velocity = SVector{6}(left[1:6])
        resultant = SVector{6}(right[7:12])
        boundary = dot(velocity, equations.left_impedance * velocity) +
                   dot(resultant, equations.right_impedance * resultant)
    end
    return energy_rate + material + jumps + boundary
end

# Carpenter-Kennedy 2N54 amplification, evaluated independently of the ODE solve.
function ck_amplification(z)
    a = (0.0, -567301805773 / 1357537059087, -2404267990393 / 2016746695238,
         -3550918686646 / 2091501179385, -1275806237668 / 842570457699)
    b = (1432997174477 / 9575080441755, 5161836677717 / 13612068292357,
         1720146321549 / 2090206949498, 3134564353537 / 4481467310338,
         2277821191437 / 14882151754819)
    value = one(z)
    residual = zero(z)
    for stage in 1:5
        residual = a[stage] * residual + z * value
        value += b[stage] * residual
    end
    return value
end

function linearized_operator(ode)
    state = zero(ode.u0)
    du = similar(state)
    dp = similar(state)
    matrix = zeros(length(state), length(state))
    step = 1.0e-6
    for column in eachindex(state)
        state[column] = step
        Trixi.rhs!(du, state, ode.p, 0.0)
        Trixi.rhs_parabolic!(dp, state, ode.p, 0.0)
        matrix[:, column] .= (du .+ dp) ./ step
        state[column] = 0
    end
    return matrix
end

@testset "Beam constructor and kernels" begin
    for T in (Float32, Float64)
        identity = Matrix{T}(I, 6, 6)
        equations = DampedIntrinsicBeamEquations1D(; mass_matrix = identity,
                                                   flexibility_matrix = identity,
                                                   damping_matrix = T(0.1) * identity)
        @test equations.maximum_wave_speed isa T
        @test eltype(equations.initial_curvature) == T
    end
    for T in (Float32, Float64)
        velocity = SVector{6, T}(1, 2, 3, 4, 5, 6)
        resultant = SVector{6, T}(2, -1, 3, -2, 5, 1)
        v, omega = SVector{3}(velocity[1:3]), SVector{3}(velocity[4:6])
        force, moment = SVector{3}(resultant[1:3]), SVector{3}(resultant[4:6])
        @test intrinsic_beam_l1(velocity) * resultant ==
              SVector{6}(cross(omega, force)...,
                         (cross(v, force) + cross(omega, moment))...)
        @test intrinsic_beam_l2(resultant) * velocity ==
              SVector{6}(cross(force, omega)...,
                         (cross(force, v) + cross(moment, omega))...)
    end
    identity = Matrix{Float64}(I, 6, 6)
    for scale in (1.0e-16, 1.0, 1.0e16)
        equations = DampedIntrinsicBeamEquations1D(; mass_matrix = scale * identity,
                                                   flexibility_matrix = identity,
                                                   damping_matrix = 0.1 * identity)
        @test equations.mass_matrix ≈ scale * identity
    end
    for invalid in (NaN, Inf, -1.0)
        bad = copy(identity)
        bad[1, 1] = invalid
        @test_throws ArgumentError DampedIntrinsicBeamEquations1D(; mass_matrix = bad,
                                                                  flexibility_matrix = identity,
                                                                  damping_matrix = identity)
    end
    for damping in (zeros(6, 6), Matrix(Diagonal([0.1, 0, 0.2, 0, 0, 0])))
        equations = DampedIntrinsicBeamEquations1D(; mass_matrix = identity,
                                                   flexibility_matrix = identity,
                                                   damping_matrix = damping)
        @test equations.damping_operator ≈ damping
        parabolic = DampedIntrinsicBeamDiffusion1D(equations)
        core_allocations(parabolic)
        @test core_allocations(parabolic) == 0
        @test abs(energy_residual(periodic_beam(damping))) < 1.0e-11
        @test abs(energy_residual(periodic_beam(damping; periodic = false);
                                  periodic = false)) < 1.0e-11
    end
end

@testset "Beam run integrity" begin
    good = (; retcode = ReturnCode.Success, t = [0.0, 1.0], u = [[0.0], [1.0]])
    @test require_complete_solution(good, 1.0) === good
    @test_throws ErrorException require_complete_solution(good, 2.0)
    @test_throws ErrorException require_complete_solution(merge(good,
                                                                (;
                                                                 retcode = ReturnCode.MaxIters)),
                                                          1.0)
    @test_throws ErrorException require_complete_solution(merge(good,
                                                                (; u = [[NaN], [1.0]])),
                                                          1.0)
    @test checkpoint_gravity((;)) == 1.0
    @test checkpoint_gravity((; gravity_multiplier = 0.0)) == 0.0
    checkpoint = (; acceleration_rms_g = 0.5, polydeg = 3, refinement_level = 2,
                  constraint_multiplier = 100.0, gravity_multiplier = 0.0,
                  run_configuration = (; reltol = 1e-5, abstol = 1e-7,
                                       steps_per_period = 64, use_sparse_jacobian = true))
    replay = cantilever_checkpoint_environment(checkpoint,
                                               (; frequency_hz = 9.0, omega = 2pi))
    @test replay["CANTILEVER_GRAVITY_MULTIPLIER"] == "0.0"
    @test parse(Float64, replay["CANTILEVER_RELTOL"]) == 1e-5
    @test parse(Float64, replay["CANTILEVER_DTMAX"]) == 1 / 64
    @test_throws ErrorException checkpoint_gravity((; gravity_multiplier = Inf))
    @test require_matching_configuration((; tolerance = 1e-12), (; tolerance = 1e-12)) ===
          nothing
    @test_throws ErrorException require_matching_configuration((; tolerance = 1e-12),
                                                               (; tolerance = 1e-8))
    configuration = (; final_time = 1.0, save_count = 2)
    provenance = (; source_identity = "fixture", julia_version = string(VERSION),
                  threads = 1)
    cached = (; run_configuration = configuration, provenance, completed = true,
              times = [0.0, 1.0], states = [[0.0], [1.0]])
    @test validate_rotating_result(cached, configuration, provenance) === cached
    @test_throws ErrorException validate_rotating_result(cached,
                                                         (; final_time = 40.0,
                                                          save_count = 401), provenance)
    @test_throws ErrorException validate_rotating_result(cached, configuration,
                                                         merge(provenance,
                                                               (;
                                                                source_identity = "changed")))
end

@testset "Beam explicit diffusion safeguard" begin
    for (degree, level, damping) in ((1, 1, 0.0), (2, 1, 0.1), (3, 2, 10.0))
        ode = periodic_beam(damping; polydeg = degree, level)
        dt = beam_explicit_timestep(ode, 0.01)
        @test 0 < dt <= StepsizeCallback(cfl = 0.01)(ode)
        eigenvalues = eigvals(linearized_operator(ode))
        @test maximum(abs, ck_amplification.(dt .* eigenvalues)) <= 1 + 1.0e-10
    end
    coarse = periodic_beam(10.0; level = 1)
    fine = periodic_beam(10.0; level = 2)
    @test beam_explicit_timestep(fine, 0.01) ≈ beam_explicit_timestep(coarse, 0.01) / 4
    @test_throws ArgumentError beam_explicit_timestep(fine, 0.01; diffusion_cfl = 1.0)
end

@testset "Non-beam parabolic regression" begin
    equations = LinearScalarAdvectionEquation1D(0.1)
    diffusion = LaplaceDiffusion1D(0.1, equations)
    initial(x, t, equations) = SVector(sinpi(2 * (x[1] - 0.1 * t)) *
                                       exp(-0.1 * (2pi)^2 * t))
    mesh = TreeMesh((0.0,), (1.0,); initial_refinement_level = 3,
                    n_cells_max = 100, periodicity = true)
    solver = DGSEM(polydeg = 3, surface_flux = flux_lax_friedrichs)
    semi = SemidiscretizationHyperbolicParabolic(mesh, (equations, diffusion), initial,
                                                 solver)
    ode = semidiscretize(semi, (0.0, 0.01))
    solution = solve(ode, ROCK4(); abstol = 1e-12, reltol = 1e-12,
                     save_everystep = false)
    require_complete_solution(solution, 0.01)
    errors = AnalysisCallback(semi, interval = 0)(solution)
    @test maximum(errors.l2) < 2e-4
    @test maximum(errors.linf) < 5e-4
    println("NONBEAM_ERRORS l2=", errors.l2, " linf=", errors.linf)
end

module RichRegression
using Trixi
end

function forcing_allocations(force, equations)
    force(SVector(0.2), 0.1, equations)
    return @allocated force(SVector(0.2), 0.1, equations)
end

@testset "Paper rich MMS" begin
    Trixi.trixi_include(RichRegression,
                        joinpath(@__DIR__, "..", "examples", "damped_intrinsic_beam",
                                 "elixir_mms_rich_physical.jl");
                        polydeg = 2, initial_refinement_level = 2, tspan = (0.0, 0.1),
                        time_int_tol = 1.0e-12)
    errors = Base.invokelatest(RichRegression.analysis_callback, RichRegression.sol)
    @test maximum(errors.l2)≈0.0008978808136296493 rtol=2.0e-7
    @test maximum(errors.linf)≈0.0029557718801171973 rtol=2.0e-7
    @test RichRegression.rich_mms_exact_audit.maximum_normalized_lower_residual <
          100 * eps()
    @test maximum(values(RichRegression.rich_mms_boundary_audit)) < 1e-10
    Base.invokelatest(forcing_allocations, RichRegression.physical_external_force,
                      RichRegression.equations_hyperbolic)
    @test Base.invokelatest(forcing_allocations, RichRegression.physical_external_force,
                            RichRegression.equations_hyperbolic) == 0
end

end # module
