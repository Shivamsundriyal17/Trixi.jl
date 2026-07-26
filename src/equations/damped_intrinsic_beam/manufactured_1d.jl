@inline _intrinsic_beam_coordinate(x::Real) = x
@inline _intrinsic_beam_coordinate(x) = first(x)

"""
    manufactured_source_damped_intrinsic_beam(
        x, t, equations, solution, solution_t, solution_x, solution_xx)

Return the full twelve-component capacity-form forcing required by arbitrary
manufactured fields for the damped intrinsic beam equations.

The last six components force the kinematic/constitutive block. Thus, a
nonzero lower block verifies a generalized mixed system, not the physical
Kelvin--Voigt model. Use
[`manufactured_force_damped_intrinsic_beam`](@ref) to enforce the physical
six-component forcing structure.
"""
function manufactured_source_damped_intrinsic_beam(x, t,
                                                   equations::DampedIntrinsicBeamDiffusion1D,
                                                   solution, solution_t,
                                                   solution_x, solution_xx)
    hyperbolic = equations.equations_hyperbolic
    coordinate = _intrinsic_beam_coordinate(x)
    u = solution(coordinate, t)
    u_t = solution_t(coordinate, t)
    u_x = solution_x(coordinate, t)
    u_xx = solution_xx(coordinate, t)

    u1 = _intrinsic_beam_block(u, 1)
    u2 = _intrinsic_beam_block(u, 7)
    u1_t = _intrinsic_beam_block(u_t, 1)
    u2_t = _intrinsic_beam_block(u_t, 7)
    u1_x = _intrinsic_beam_block(u_x, 1)
    u2_x = _intrinsic_beam_block(u_x, 7)
    u1_xx = _intrinsic_beam_block(u_xx, 1)

    l1 = intrinsic_beam_l1(u1)
    l2 = intrinsic_beam_l2(u2)
    l1_x = intrinsic_beam_l1(u1_x)
    elastic_strain = hyperbolic.flexibility_matrix * u2
    damping_resultant = hyperbolic.damping_operator *
                        (u1_x -
                         transpose(hyperbolic.geometry_matrix) * u1 +
                         transpose(l1) * elastic_strain)
    damping_resultant_x = hyperbolic.damping_operator *
                          (u1_xx -
                           transpose(hyperbolic.geometry_matrix) * u1_x +
                           transpose(l1_x) * elastic_strain +
                           transpose(l1) *
                           (hyperbolic.flexibility_matrix * u2_x))

    internal_upper = hyperbolic.geometry_matrix * u2 -
                     l1 * (hyperbolic.mass_matrix * u1) -
                     l2 * elastic_strain -
                     intrinsic_beam_l2(damping_resultant) * elastic_strain +
                     hyperbolic.geometry_matrix * damping_resultant
    internal_lower = -transpose(hyperbolic.geometry_matrix) * u1 +
                     transpose(l1) * elastic_strain

    residual_upper = hyperbolic.mass_matrix * u1_t -
                     u2_x - damping_resultant_x - internal_upper
    residual_lower = hyperbolic.flexibility_matrix * u2_t -
                     u1_x - internal_lower
    return SVector{12}(residual_upper..., residual_lower...)
end

"""
    manufactured_force_damped_intrinsic_beam(
        x, t, equations, solution, solution_t, solution_x, solution_xx;
        atol=nothing, rtol=nothing)

Return the six physical distributed forces and moments for a manufactured
solution. Throw an `ArgumentError` if the manufactured fields require a
nonzero forcing in the lower kinematic/constitutive block.
"""
function manufactured_force_damped_intrinsic_beam(x, t,
                                                  equations::DampedIntrinsicBeamDiffusion1D,
                                                  solution, solution_t,
                                                  solution_x, solution_xx;
                                                  atol = nothing, rtol = nothing)
    source = manufactured_source_damped_intrinsic_beam(x, t, equations,
                                                       solution, solution_t,
                                                       solution_x, solution_xx)
    ScalarT = typeof(real(zero(eltype(source))))
    absolute_tolerance = isnothing(atol) ? 100 * eps(ScalarT) : atol
    relative_tolerance = isnothing(rtol) ? 100 * eps(ScalarT) : rtol
    scale = max(norm(source, Inf), one(ScalarT))
    tolerance = absolute_tolerance + relative_tolerance * scale
    lower_norm = norm(SVector{6}(source[7:12]), Inf)

    lower_norm <= tolerance ||
        throw(ArgumentError("manufactured fields require a nonzero lower-block " *
                            "source (infinity norm $lower_norm); they are not " *
                            "a physical Kelvin-Voigt manufactured solution"))

    return SVector{6}(source[1:6])
end
