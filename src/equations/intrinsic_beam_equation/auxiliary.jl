"""
    tilde(v::Vector{T}) where T<:Real

Computes cross product matrix, s.t.
    v x w = tilde(v) * w
"""
function tilde(v::Vector{T}) where T<:Real
    if length(v) != 3
        error("Vector must be of length three")
    end

    return [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
end


"""
    reshape_solution(u)

For a solution vector u of the intrinsic beam equation of length 12*(polydeg+1)*ncells at one timestep an 2D-array of size 
12 x (polydeg+1) is returned
"""
function reshape_solution(u)
    nnodes = Int(length(u) / 12)
    u_mat = Matrix{Float64}(undef, 12, nnodes)

    for j in 1:nnodes
        u_mat[:,j] = u[(j-1)*12+1:j*12] 
    end

    return u_mat
end


"""
    function generate_posdef_matrix(n::Integer)

Function that generates positive definite matrix
"""
function generate_posdef_matrix(n::Integer)
    A = round.( 2*rand(n,n) )
    A = A' * A + I(n)
end


"""
    integrate_trapez(t, f::Vector{T}) where T<:Real

Integrates a discrete scalar function f, by applying the trapeze rule.
The vector f should contain the values of function at nodes t
"""
function integrate_trapez(t, f::Vector{T}) where T<:Real
    N = length(t)
    
    Δt = Vector{Float64}(undef, N-1)
    for j in 1:N-1
        Δt[j] = t[j+1] - t[j]
    end

    integral = 0.0
    for j in 1:N-1
        integral += 0.5 * Δt[j] * (f[j+1] + f[j]) 
    end

    return integral
end


"""
    integrate_trapez(t, f::Matrix{T}) where T<:Real

Integrates vector valued function f, by applying the trapeze rule.
The matrix f should contain the values of the function at nodes t
"""
function integrate_trapez(t, f::Matrix{T}) where T<:Real
    N = length(t)
    
    if size(f, 1) != N
        f = f'
        transpose_result = true
    end

    Δt = Vector{Float64}(undef, N-1)
    for j in 1:N-1
        Δt[j] = t[j+1] - t[j]
    end

    integral = zeros(size(f,2))
    for j in 1:N-1
        integral += 0.5 * Δt[j] * (f[j+1,:] + f[j,:]) 
    end

    return integral
end
