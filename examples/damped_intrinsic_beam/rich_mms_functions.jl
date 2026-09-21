module RichMMSFunctions

using Trixi
using Trixi: SMatrix

export RichBeamMMS, RichField, RichForce, RichFlux

"""Typed parameters constructed after the elixir's overridable assignments."""
struct RichBeamMMS{T, Profile}
    lambda::T
    flexibility_inverse::SMatrix{6, 6, T, 36}
end

function RichBeamMMS(lambda, flexibility_inverse, profile)
    profile in ("quartic", "exp1", "exp2", "poly7") ||
        throw(ArgumentError("unknown rich MMS profile $profile"))
    return RichBeamMMS{typeof(lambda), Symbol(profile)}(lambda, flexibility_inverse)
end

@inline profile(::RichBeamMMS{T, :quartic}, x, t) where {T} = ((1 + 0.25 * (x + t))^4,
                                                               (1 + 0.25 * (x + t))^3,
                                                               0.75 *
                                                               (1 + 0.25 * (x + t))^2)
@inline profile(::RichBeamMMS{T, :exp1}, x, t) where {T} = (exp(x + t), exp(x + t),
                                                            exp(x + t))
@inline profile(::RichBeamMMS{T, :exp2}, x, t) where {T} = (exp(2 * (x + t)),
                                                            2 * exp(2 * (x + t)),
                                                            4 * exp(2 * (x + t)))
@inline profile(::RichBeamMMS{T, :poly7}, x, t) where {T} = ((1 + 0.5 * (x + t))^7,
                                                             3.5 * (1 + 0.5 * (x + t))^6,
                                                             10.5 * (1 + 0.5 * (x + t))^5)

struct RichField{Derivative, Parameters}
    parameters::Parameters
end
RichField(parameters, derivative) = RichField{derivative, typeof(parameters)}(parameters)

@inline function (field::RichField{Derivative})(x, t) where {Derivative}
    parameters = field.parameters
    value, first_derivative, second_derivative = profile(parameters, x, t)
    direction = SVector(3.0, 2.0, 3.0, 1.0, 1.0, 1.0)
    geometry_action = SVector(0.0, 1.0, -1.0, 0.0, 0.0, 0.0)
    axial_shift = SVector(1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    growth_rate = parameters.lambda - 1
    rotated_direction = direction + growth_rate * t * geometry_action
    if Derivative === :value
        upper = value * rotated_direction
        lower = parameters.flexibility_inverse * (upper - parameters.lambda * axial_shift)
    elseif Derivative === :time
        upper = growth_rate * value * geometry_action + first_derivative * rotated_direction
        lower = parameters.flexibility_inverse * upper
    elseif Derivative === :space
        upper = first_derivative * rotated_direction
        lower = parameters.flexibility_inverse * upper
    else
        upper = second_derivative * rotated_direction
        lower = parameters.flexibility_inverse * upper
    end
    return SVector{12}(upper..., lower...)
end

struct RichForce{Parameters}
    parameters::Parameters
end

struct RichFlux{T}
    sigma::T
end

@inline function (surface_flux::RichFlux)(u_ll, u_rr, orientation, equations)
    return 0.5 * (flux(u_ll, orientation, equations) + flux(u_rr, orientation, equations)) +
           0.5 * surface_flux.sigma * equations.propagation_matrix_abs * (u_ll - u_rr)
end

@inline function (force::RichForce)(x, t, equations)
    parameters = force.parameters
    return manufactured_force_damped_intrinsic_beam(x, t,
                                                    DampedIntrinsicBeamDiffusion1D(equations),
                                                    RichField(parameters, :value),
                                                    RichField(parameters, :time),
                                                    RichField(parameters, :space),
                                                    RichField(parameters, :space2))
end

end # module
