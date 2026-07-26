@doc raw"""
    DampedIntrinsicBeamDiffusion1D(equations_hyperbolic)

Kelvin--Voigt operator paired with
[`DampedIntrinsicBeamEquations1D`](@ref). It computes
```math
r_\tau = C_\tau C^{-1}
\left(u_1' - E^\top u_1 + L_1(u_1)^\top C u_2\right)
```
from Trixi's auxiliary gradient and contributes the capacity-scaled flux
``\Gamma^{-1}[r_\tau;0]`` and nonlinear source ``\Gamma^{-1}Q(u,r)``.
"""
struct DampedIntrinsicBeamDiffusion1D{EquationsHyperbolic} <:
       AbstractEquationsParabolic{1, 12, GradientVariablesConservative}
    equations_hyperbolic::EquationsHyperbolic
end

function DampedIntrinsicBeamDiffusion1D(equations_hyperbolic::DampedIntrinsicBeamEquations1D)
    return DampedIntrinsicBeamDiffusion1D{typeof(equations_hyperbolic)}(equations_hyperbolic)
end

function varnames(variable_mapping,
                  equations_parabolic::DampedIntrinsicBeamDiffusion1D)
    return varnames(variable_mapping, equations_parabolic.equations_hyperbolic)
end

@inline function _intrinsic_beam_block(vector, first_index)
    return SVector{6}(ntuple(index -> vector[first_index + index - 1], 6))
end

"""
    intrinsic_beam_damping_resultant(u, gradient,
                                     equations::DampedIntrinsicBeamDiffusion1D)

Evaluate the six-component Kelvin--Voigt auxiliary resultant ``r_\tau`` at one
point.
"""
@inline function intrinsic_beam_damping_resultant(u, gradient,
                                                  equations::DampedIntrinsicBeamDiffusion1D)
    hyperbolic = equations.equations_hyperbolic
    u1 = _intrinsic_beam_block(u, 1)
    u2 = _intrinsic_beam_block(u, 7)
    u1_x = _intrinsic_beam_block(gradient, 1)
    elastic_strain = hyperbolic.flexibility_matrix * u2
    constitutive_argument = u1_x -
                            transpose(hyperbolic.geometry_matrix) * u1 +
                            transpose(intrinsic_beam_l1(u1)) * elastic_strain
    return hyperbolic.damping_operator * constitutive_argument
end

@inline function flux(u, gradient, orientation::Integer,
                      equations::DampedIntrinsicBeamDiffusion1D)
    hyperbolic = equations.equations_hyperbolic
    damping_resultant = intrinsic_beam_damping_resultant(u, gradient, equations)
    upper_flux = hyperbolic.mass_inverse * damping_resultant
    ScalarT = eltype(upper_flux)
    return SVector{12}(ntuple(index -> index <= 6 ? upper_flux[index] :
                                       zero(ScalarT), 12))
end

"""
    intrinsic_beam_gradient_source(u, gradient, x, t,
                                   equations::DampedIntrinsicBeamDiffusion1D)

Evaluate the capacity-scaled nonlinear and external source at one point.
"""
@inline function intrinsic_beam_gradient_source(u, gradient, x, t,
                                                equations::DampedIntrinsicBeamDiffusion1D)
    hyperbolic = equations.equations_hyperbolic
    u1 = _intrinsic_beam_block(u, 1)
    u2 = _intrinsic_beam_block(u, 7)
    damping_resultant = intrinsic_beam_damping_resultant(u, gradient, equations)

    elastic_strain = hyperbolic.flexibility_matrix * u2
    external_force = hyperbolic.external_force(x, t, hyperbolic)
    length(external_force) == 6 ||
        throw(ArgumentError("external_force must return exactly six entries"))
    external_force_static = SVector{6}(external_force)

    upper_source = hyperbolic.geometry_matrix * u2 -
                   intrinsic_beam_l1(u1) * (hyperbolic.mass_matrix * u1) -
                   intrinsic_beam_l2(u2) * elastic_strain -
                   intrinsic_beam_l2(damping_resultant) * elastic_strain +
                   hyperbolic.geometry_matrix * damping_resultant +
                   external_force_static
    lower_source = -transpose(hyperbolic.geometry_matrix) * u1 +
                   transpose(intrinsic_beam_l1(u1)) * elastic_strain

    scaled_upper = hyperbolic.mass_inverse * upper_source
    scaled_lower = hyperbolic.flexibility_inverse * lower_source
    return SVector{12}(scaled_upper..., scaled_lower...)
end

function add_gradient_dependent_source_terms!(du, gradients, u, t,
                                              mesh::TreeMesh{1},
                                              equations::DampedIntrinsicBeamDiffusion1D,
                                              dg::DG, cache)
    @unpack node_coordinates = cache.elements

    @threaded for element in eachelement(dg, cache)
        for i in eachnode(dg)
            u_local = get_node_vars(u, equations, dg, i, element)
            gradient_local = get_node_vars(gradients, equations, dg, i, element)
            x_local = get_node_coords(node_coordinates, equations, dg, i, element)
            source_local = intrinsic_beam_gradient_source(u_local, gradient_local,
                                                          x_local, t, equations)
            add_to_node_vars!(du, source_local, equations, dg, i, element)
        end
    end

    return nothing
end

@inline function (boundary_condition::BoundaryConditionDampedIntrinsicBeam)(flux_inner,
                                                                            u_inner,
                                                                            normal,
                                                                            direction,
                                                                            x, t,
                                                                            ::Gradient,
                                                                            equations::DampedIntrinsicBeamDiffusion1D)
    if isodd(direction)
        hyperbolic = equations.equations_hyperbolic
        left_velocity = _intrinsic_beam_boundary_value(boundary_condition.left_velocity,
                                                       x, t, hyperbolic)
        inner_resultant = _intrinsic_beam_block(u_inner, 7)
        return SVector{12}(left_velocity..., inner_resultant...)
    else
        return u_inner
    end
end

@inline function (boundary_condition::BoundaryConditionDampedIntrinsicBeam)(flux_inner,
                                                                            u_inner,
                                                                            normal,
                                                                            direction,
                                                                            x, t,
                                                                            ::Divergence,
                                                                            equations::DampedIntrinsicBeamDiffusion1D)
    if isodd(direction)
        return flux_inner
    else
        hyperbolic = equations.equations_hyperbolic
        damping_resultant = _intrinsic_beam_boundary_value(boundary_condition.right_damping_resultant,
                                                           x, t, hyperbolic)
        upper_flux = hyperbolic.mass_inverse * damping_resultant
        ScalarT = eltype(upper_flux)
        return SVector{12}(ntuple(index -> index <= 6 ? upper_flux[index] :
                                           zero(ScalarT), 12))
    end
end
