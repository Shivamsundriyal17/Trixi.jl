function generate_L1_matrix(x::Vector{Float64})
    @assert length(x) == 6 "Input vector must have 6 entries."
    
    v = x[1:3]
    ω = x[4:6]
    L1 = zeros(6, 6) 
    # Fill in the blocks
    L1[1:3, 1:3] = tilde(ω)
    L1[4:6, 1:3] = tilde(v)
    L1[4:6, 4:6] = tilde(ω)
    
    return L1
end

function generate_L2_matrix(x::Vector{Float64})
    @assert length(x) == 6 "Input vector must have 6 entries."
    f = x[1:3]
    m = x[4:6]
    L2 = zeros(6, 6)
    L2[1:3, 4:6] = tilde(f)
    L2[4:6, 1:3] = tilde(f)
    L2[4:6, 4:6] = tilde(m)
    
    return L2
end

function generate_E_matrix(κ0::Matrix{Float64})
    @assert length(κ0) == 3 "Input vector κ0 must have 3 entries."
    
    e1 = [1.0, 0.0, 0.0]
    E = zeros(6, 6)
    E[1:3, 1:3] = tilde(κ0)
    E[4:6, 1:3] = tilde(e1)
    E[4:6, 4:6] = tilde(κ0)
    
    return E
end

function tilde(v::AbstractMatrix{T}) where T<:Real
    # Check that the matrix has exactly 1 row and 3 columns
    if size(v) != (1, 3)
        error("Matrix must be of size 1x3")
    end

    # Extract the single row vector from the matrix
    v_row = v[1, :]

    # Construct the skew-symmetric matrix
    return [0 -v_row[3] v_row[2]; 
            v_row[3] 0 -v_row[1]; 
            -v_row[2] v_row[1] 0]
end
