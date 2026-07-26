"""
    intrinsic_beam_skew(vector)

Return the skew-symmetric matrix associated with a three-component vector, i.e.,
the matrix satisfying `intrinsic_beam_skew(a) * b == cross(a, b)`.
"""
@inline function intrinsic_beam_skew(vector)
    length(vector) == 3 ||
        throw(ArgumentError("the vector must contain exactly three entries"))

    return @SMatrix [zero(vector[1]) -vector[3] vector[2];
                     vector[3] zero(vector[1]) -vector[1];
                     -vector[2] vector[1] zero(vector[1])]
end

"""
    intrinsic_beam_l1(velocity)

Return the six-by-six intrinsic-beam operator
```math
L_1([v;\\omega]) =
\\begin{bmatrix}
\\widetilde{\\omega} & 0 \\\\
\\widetilde{v} & \\widetilde{\\omega}
\\end{bmatrix}.
```
"""
@inline function intrinsic_beam_l1(velocity)
    length(velocity) == 6 ||
        throw(ArgumentError("the velocity must contain exactly six entries"))

    linear_velocity = SVector(velocity[1], velocity[2], velocity[3])
    angular_velocity = SVector(velocity[4], velocity[5], velocity[6])
    zero_block = zero(SMatrix{3, 3, eltype(linear_velocity), 9})

    return SMatrix{6, 6}([intrinsic_beam_skew(angular_velocity) zero_block;
                          intrinsic_beam_skew(linear_velocity) intrinsic_beam_skew(angular_velocity)])
end

"""
    intrinsic_beam_l2(resultant)

Return the six-by-six intrinsic-beam operator
```math
L_2([f;m]) =
\\begin{bmatrix}
0 & \\widetilde{f} \\\\
\\widetilde{f} & \\widetilde{m}
\\end{bmatrix}.
```
"""
@inline function intrinsic_beam_l2(resultant)
    length(resultant) == 6 ||
        throw(ArgumentError("the resultant must contain exactly six entries"))

    force = SVector(resultant[1], resultant[2], resultant[3])
    moment = SVector(resultant[4], resultant[5], resultant[6])
    zero_block = zero(SMatrix{3, 3, eltype(force), 9})

    return SMatrix{6, 6}([zero_block intrinsic_beam_skew(force);
                          intrinsic_beam_skew(force) intrinsic_beam_skew(moment)])
end

"""
    intrinsic_beam_e(initial_curvature)

Return the constant six-by-six intrinsic-beam geometry operator
```math
E(\\kappa_0) =
\\begin{bmatrix}
\\widetilde{\\kappa_0} & 0 \\\\
\\widetilde{e_1} & \\widetilde{\\kappa_0}
\\end{bmatrix}.
```
"""
@inline function intrinsic_beam_e(initial_curvature)
    length(initial_curvature) == 3 ||
        throw(ArgumentError("the initial curvature must contain exactly three entries"))

    curvature = SVector(initial_curvature[1], initial_curvature[2],
                        initial_curvature[3])
    ScalarT = eltype(curvature)
    e1 = SVector(one(ScalarT), zero(ScalarT), zero(ScalarT))
    zero_block = zero(SMatrix{3, 3, ScalarT, 9})

    return SMatrix{6, 6}([intrinsic_beam_skew(curvature) zero_block;
                          intrinsic_beam_skew(e1) intrinsic_beam_skew(curvature)])
end
