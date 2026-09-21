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

    v1, v2, v3, w1, w2, w3 = velocity
    z = zero(v1)
    # Scalar entries avoid the heap-allocated block concatenation of matrices.
    return @SMatrix [z -w3 w2 z z z;
                     w3 z -w1 z z z;
                     -w2 w1 z z z z;
                     z -v3 v2 z -w3 w2;
                     v3 z -v1 w3 z -w1;
                     -v2 v1 z -w2 w1 z]
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

    f1, f2, f3, m1, m2, m3 = resultant
    z = zero(f1)
    return @SMatrix [z z z z -f3 f2;
                     z z z f3 z -f1;
                     z z z -f2 f1 z;
                     z -f3 f2 z -m3 m2;
                     f3 z -f1 m3 z -m1;
                     -f2 f1 z -m2 m1 z]
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

    k1, k2, k3 = initial_curvature
    return intrinsic_beam_l1(SVector(one(k1), zero(k1), zero(k1), k1, k2, k3))
end
