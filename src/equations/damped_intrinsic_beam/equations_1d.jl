@inline function _zero_intrinsic_beam_data(x, t, equations)
    ScalarT = promote_type(eltype(x), typeof(t))
    return zero(SVector{6, ScalarT})
end

function _intrinsic_beam_matrix(matrix::AbstractMatrix{<:Real}, name;
                                positive_definite)
    size(matrix) == (6, 6) ||
        throw(ArgumentError("$name must be a 6-by-6 matrix"))

    ScalarT = float(eltype(matrix))
    dense_matrix = Matrix{ScalarT}(matrix)
    scale = max(LinearAlgebra.opnorm(dense_matrix, Inf), one(ScalarT))
    tolerance = 100 * eps(ScalarT) * scale
    isapprox(dense_matrix, transpose(dense_matrix);
             atol = tolerance, rtol = 100 * eps(ScalarT)) ||
        throw(ArgumentError("$name must be symmetric"))

    symmetric_matrix = LinearAlgebra.Symmetric(0.5 *
                                               (dense_matrix +
                                                transpose(dense_matrix)))
    eigenvalues = LinearAlgebra.eigvals(symmetric_matrix)
    if positive_definite
        minimum(eigenvalues) > tolerance ||
            throw(ArgumentError("$name must be positive definite"))
    else
        minimum(eigenvalues) >= -tolerance ||
            throw(ArgumentError("$name must be positive semidefinite"))
    end

    return SMatrix{6, 6, ScalarT}(symmetric_matrix)
end

function _intrinsic_beam_spd_sqrt(matrix::SMatrix{6, 6})
    decomposition = LinearAlgebra.eigen(LinearAlgebra.Symmetric(Matrix(matrix)))
    square_root = decomposition.vectors *
                  LinearAlgebra.Diagonal(sqrt.(decomposition.values)) *
                  transpose(decomposition.vectors)
    return SMatrix{6, 6}(square_root)
end

@doc raw"""
    DampedIntrinsicBeamEquations1D(; mass_matrix, flexibility_matrix,
                                    damping_matrix, initial_curvature,
                                    external_force)

Constant-coefficient first-order part of the geometrically exact intrinsic
beam equations with Kelvin--Voigt damping. The state is
```math
u = [u_1; u_2] = [v;\omega;f;m],
```
where `u_1` contains linear and angular velocities and `u_2` contains elastic
sectional force and moment resultants.

The matrices `mass_matrix` ``M`` and `flexibility_matrix` ``C`` must be
symmetric positive definite. The `damping_matrix` ``C_\tau`` may be symmetric
positive semidefinite, including zero, and must satisfy the energy
compatibility condition that ``C_\tau C^{-1}`` is symmetric positive
semidefinite. Coefficients are constant in space.

The callable `external_force(x, t, equations)` must return the six distributed
forces and moments. The default is zero.

Use this equation together with [`DampedIntrinsicBeamDiffusion1D`](@ref). The
constant capacity inverse is included in both operators, so no modification of
Trixi's generic right-hand side is required.
"""
struct DampedIntrinsicBeamEquations1D{RealT <: Real, ExternalForce} <:
       AbstractEquations{1, 12}
    mass_matrix::SMatrix{6, 6, RealT, 36}
    flexibility_matrix::SMatrix{6, 6, RealT, 36}
    damping_matrix::SMatrix{6, 6, RealT, 36}
    initial_curvature::SVector{3, RealT}
    external_force::ExternalForce
    geometry_matrix::SMatrix{6, 6, RealT, 36}
    mass_inverse::SMatrix{6, 6, RealT, 36}
    flexibility_inverse::SMatrix{6, 6, RealT, 36}
    damping_operator::SMatrix{6, 6, RealT, 36}
    capacity_matrix::SMatrix{12, 12, RealT, 144}
    capacity_inverse::SMatrix{12, 12, RealT, 144}
    propagation_matrix::SMatrix{12, 12, RealT, 144}
    propagation_matrix_abs::SMatrix{12, 12, RealT, 144}
    propagation_matrix_plus::SMatrix{12, 12, RealT, 144}
    propagation_matrix_minus::SMatrix{12, 12, RealT, 144}
    left_impedance::SMatrix{6, 6, RealT, 36}
    right_impedance::SMatrix{6, 6, RealT, 36}
    maximum_wave_speed::RealT
end

function DampedIntrinsicBeamEquations1D(;
                                        mass_matrix::AbstractMatrix{<:Real},
                                        flexibility_matrix::AbstractMatrix{<:Real},
                                        damping_matrix::AbstractMatrix{<:Real},
                                        initial_curvature = zeros(3),
                                        external_force = _zero_intrinsic_beam_data)
    length(initial_curvature) == 3 ||
        throw(ArgumentError("initial_curvature must contain exactly three entries"))

    RealT = promote_type(float(eltype(mass_matrix)),
                         float(eltype(flexibility_matrix)),
                         float(eltype(damping_matrix)),
                         float(eltype(initial_curvature)))
    mass = SMatrix{6, 6, RealT}(_intrinsic_beam_matrix(mass_matrix, "mass_matrix";
                                                       positive_definite = true))
    flexibility = SMatrix{6, 6, RealT}(_intrinsic_beam_matrix(flexibility_matrix,
                                                              "flexibility_matrix";
                                                              positive_definite = true))
    damping = SMatrix{6, 6, RealT}(_intrinsic_beam_matrix(damping_matrix,
                                                          "damping_matrix";
                                                          positive_definite = false))
    curvature = SVector{3, RealT}(initial_curvature)

    mass_inverse = inv(mass)
    flexibility_inverse = inv(flexibility)
    damping_operator_raw = damping * flexibility_inverse
    scale = max(LinearAlgebra.opnorm(damping_operator_raw, Inf), one(RealT))
    tolerance = 500 * eps(RealT) * scale
    isapprox(damping_operator_raw, transpose(damping_operator_raw);
             atol = tolerance, rtol = 500 * eps(RealT)) ||
        throw(ArgumentError("damping_matrix * inv(flexibility_matrix) must be symmetric"))
    damping_operator = SMatrix{6, 6, RealT}(0.5 *
                                            (damping_operator_raw +
                                             transpose(damping_operator_raw)))
    minimum(LinearAlgebra.eigvals(LinearAlgebra.Symmetric(Matrix(damping_operator)))) >=
    -tolerance ||
        throw(ArgumentError("the compatible damping operator must be positive semidefinite"))

    zero66 = zero(SMatrix{6, 6, RealT, 36})
    identity66 = SMatrix{6, 6, RealT}(I)
    capacity = SMatrix{12, 12, RealT}([mass zero66;
                                       zero66 flexibility])
    capacity_inverse = SMatrix{12, 12, RealT}([mass_inverse zero66;
                                               zero66 flexibility_inverse])
    pi_matrix = SMatrix{12, 12, RealT}([zero66 -identity66;
                                        -identity66 zero66])
    propagation = capacity_inverse * pi_matrix

    mass_sqrt = _intrinsic_beam_spd_sqrt(mass)
    flexibility_sqrt = _intrinsic_beam_spd_sqrt(flexibility)
    capacity_sqrt = SMatrix{12, 12, RealT}([mass_sqrt zero66;
                                            zero66 flexibility_sqrt])
    capacity_inverse_sqrt = inv(capacity_sqrt)

    symmetric_propagation = capacity_inverse_sqrt * pi_matrix *
                            capacity_inverse_sqrt
    decomposition = LinearAlgebra.eigen(LinearAlgebra.Symmetric(Matrix(symmetric_propagation)))
    symmetric_propagation_abs = decomposition.vectors *
                                LinearAlgebra.Diagonal(abs.(decomposition.values)) *
                                transpose(decomposition.vectors)
    propagation_abs = SMatrix{12, 12, RealT}(capacity_inverse_sqrt *
                                             symmetric_propagation_abs *
                                             capacity_sqrt)
    propagation_plus = 0.5 * (propagation + propagation_abs)
    propagation_minus = 0.5 * (propagation - propagation_abs)

    capacity_propagation_abs = capacity * propagation_abs
    left_impedance = SMatrix{6, 6, RealT}(capacity_propagation_abs[1:6, 1:6])
    right_impedance = SMatrix{6, 6, RealT}(capacity_propagation_abs[7:12, 7:12])
    maximum_wave_speed = maximum(abs, decomposition.values)

    geometry_matrix = intrinsic_beam_e(curvature)
    return DampedIntrinsicBeamEquations1D{RealT, typeof(external_force)}(mass,
                                                                         flexibility,
                                                                         damping,
                                                                         curvature,
                                                                         external_force,
                                                                         geometry_matrix,
                                                                         mass_inverse,
                                                                         flexibility_inverse,
                                                                         damping_operator,
                                                                         capacity,
                                                                         capacity_inverse,
                                                                         propagation,
                                                                         propagation_abs,
                                                                         propagation_plus,
                                                                         propagation_minus,
                                                                         left_impedance,
                                                                         right_impedance,
                                                                         maximum_wave_speed)
end

function varnames(::typeof(cons2cons), ::DampedIntrinsicBeamEquations1D)
    return ("v1", "v2", "v3", "omega1", "omega2", "omega3",
            "f1", "f2", "f3", "m1", "m2", "m3")
end

function varnames(::typeof(cons2prim), equations::DampedIntrinsicBeamEquations1D)
    varnames(cons2cons, equations)
end

@inline cons2prim(u, ::DampedIntrinsicBeamEquations1D) = u

@inline function cons2entropy(u, equations::DampedIntrinsicBeamEquations1D)
    return equations.capacity_matrix * u
end

@inline function entropy(u, equations::DampedIntrinsicBeamEquations1D)
    return 0.5f0 * dot(u, equations.capacity_matrix * u)
end

@inline function flux(u, orientation::Integer,
                      equations::DampedIntrinsicBeamEquations1D)
    return equations.propagation_matrix * u
end

"""
    flux_upwind(u_ll, u_rr, orientation, equations::DampedIntrinsicBeamEquations1D)

Return the exact characteristic upwind flux of the constant-coefficient
intrinsic-beam propagation operator.
"""
@inline function flux_upwind(u_ll, u_rr, orientation::Integer,
                             equations::DampedIntrinsicBeamEquations1D)
    return equations.propagation_matrix_plus * u_ll +
           equations.propagation_matrix_minus * u_rr
end

@inline have_constant_speed(::DampedIntrinsicBeamEquations1D) = True()

@inline function max_abs_speeds(equations::DampedIntrinsicBeamEquations1D)
    return equations.maximum_wave_speed
end

@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                     equations::DampedIntrinsicBeamEquations1D)
    return equations.maximum_wave_speed
end

struct BoundaryConditionDampedIntrinsicBeam{LeftVelocity, RightResultant,
                                            RightDampingResultant}
    left_velocity::LeftVelocity
    right_resultant::RightResultant
    right_damping_resultant::RightDampingResultant
end

"""
    BoundaryConditionDampedIntrinsicBeam(;
        left_velocity=_zero_intrinsic_beam_data,
        right_resultant=_zero_intrinsic_beam_data,
        right_damping_resultant=_zero_intrinsic_beam_data)

Create the split intrinsic-beam boundary closure. The three callables have the
signature `(x, t, equations)` and prescribe ``u_1`` at the left boundary,
``u_2`` at the right boundary, and the Kelvin--Voigt resultant ``r_\tau`` at
the right boundary, respectively.

For a prescribed total right resultant ``\beta=u_2+r_\tau``, the two right
targets must be chosen consistently. Homogeneous cantilever data use the
defaults.
"""
function BoundaryConditionDampedIntrinsicBeam(;
                                              left_velocity = _zero_intrinsic_beam_data,
                                              right_resultant = _zero_intrinsic_beam_data,
                                              right_damping_resultant = _zero_intrinsic_beam_data)
    return BoundaryConditionDampedIntrinsicBeam(left_velocity, right_resultant,
                                                right_damping_resultant)
end

@inline function _intrinsic_beam_boundary_value(function_, x, t, equations)
    value = function_(x, t, equations)
    length(value) == 6 ||
        throw(ArgumentError("intrinsic-beam boundary data must contain six entries"))
    return SVector{6}(value)
end

@inline function (boundary_condition::BoundaryConditionDampedIntrinsicBeam)(u_inner,
                                                                            orientation,
                                                                            direction,
                                                                            x, t,
                                                                            surface_flux_function,
                                                                            equations::DampedIntrinsicBeamEquations1D)
    if isodd(direction)
        left_velocity = _intrinsic_beam_boundary_value(boundary_condition.left_velocity,
                                                       x, t, equations)
        difference = SVector{6}(u_inner[1:6]) - left_velocity
        outer_resultant = SVector{6}(u_inner[7:12]) +
                          equations.left_impedance * difference
        u_outer = SVector{12}(left_velocity..., outer_resultant...)
        return surface_flux_function(u_outer, u_inner, orientation, equations)
    else
        right_resultant = _intrinsic_beam_boundary_value(boundary_condition.right_resultant,
                                                         x, t, equations)
        difference = SVector{6}(u_inner[7:12]) - right_resultant
        outer_velocity = SVector{6}(u_inner[1:6]) -
                         equations.right_impedance * difference
        u_outer = SVector{12}(outer_velocity..., right_resultant...)
        return surface_flux_function(u_inner, u_outer, orientation, equations)
    end
end
