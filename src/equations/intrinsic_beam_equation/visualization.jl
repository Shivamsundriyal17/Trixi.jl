

"""
    calc_deformation_trapez(C0, node_coordinates, k, κ, γ)

Calculates rotation matrix S^{-1}, where S contains the body attached basis
using S0 as initial value for the recursion formula. Afterwards, the posiiton
line of the deformed beam is determined.
"""
function calc_deformation_trapez(S0, node_coordinates, k, κ, γ)
    N = length(node_coordinates)
    II = diagm(ones(3))

    # determine cross product matrices
    k_tilde = Array{Float64}(undef, 3, 3, N)
    κ_tilde = Array{Float64}(undef, 3, 3, N)
    γ_tilde = Array{Float64}(undef, 3, 3, N)
    dx = Vector{Float64}(undef, N-1)
    for j in 1:N-1
        k_tilde[:,:,j] = tilde(k(node_coordinates[j]))
        κ_tilde[:,:,j] = tilde(κ[:,j])
        γ_tilde[:,:,j] = tilde(γ[:,j])
        dx[j] = node_coordinates[j+1] - node_coordinates[j]
    end
    k_tilde[:,:,end] = tilde(k(node_coordinates[end]))
    κ_tilde[:,:,end] = tilde(κ[:,end])
    γ_tilde[:,:,end] = tilde(γ[:,end])

    # determine S^{-1} using recursion formula
    S_inv = Array{Float64}(undef, 3, 3, N)
    S_inv[:,:,1] = S0
    for j in 1:N-1
        S_inv[:,:,j+1] = inv( II + 0.5*dx[j]*(κ_tilde[:,:,j+1] + k_tilde[:,:,j+1]) ) * ( II - 0.5*dx[j]*(κ_tilde[:,:,j] + k_tilde[:,:,j]) ) * S_inv[:,:,j]
    end

    # determine position line using recursion formula. note that S is orthogonal as rotation matrix, which is why S^{-1} = S^T
    r = Matrix{Float64}(undef, 3, length(node_coordinates))
    r[:,1] = [0.0; 0.0; 0.0]
    e1 = [1.0; 0.0; 0.0]
    for j in 1:N-1
        r[:,j+1] = r[:,j] + 0.5*dx[j] * ( S_inv[:,:,j]' *(γ[:,j]+e1) + S_inv[:,:,j+1]' *(γ[:,j+1]+e1) )
    end

    return r
end



"""
    calc_rot_mat_y(U, t)

Calculates rotation matrix for rotation around y-axis by integrating the
associated angular velocities in U
"""
function calc_rot_mat_y(U, t)
    Nt = length(t)
    R = Array{Float64}(undef, 3, 3, Nt)
    R[:,:,1] = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    alpha = zeros(Nt)
    for j in 1:Nt-1
        alpha[j+1] = alpha[j] + (t[j+1]-t[j]) * U[11,1,j]
        R[:,:,j+1] = [cos(alpha[j+1]) 0 -sin(alpha[j+1]); 0 1 0; sin(alpha[j+1]) 0 cos(alpha[j+1])]
    end

    return R, alpha
end


"""
    calc_rot_mat_z(U, t)

    Calculates rotation matrix for rotation around z-axis by integrating the
        associated angular velocities in U
"""
function calc_rot_mat_z(U, t)
    Nt = length(t)
    R = Array{Float64}(undef, 3, 3, Nt)
    R[:,:,1] = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    alpha = zeros(Nt)
    for j in 1:Nt-1
        alpha[j+1] = alpha[j] + (t[j+1]-t[j]) * U[12,1,j]
        R[:,:,j+1] = [cos(alpha[j+1]) sin(alpha[j+1]) 0; -sin(alpha[j+1]) cos(alpha[j+1]) 0; 0 0 1]
    end

    return R, alpha
end