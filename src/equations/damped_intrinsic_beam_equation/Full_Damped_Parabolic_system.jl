struct Damped_Full_Parabolic{RealT<:Real} <: Trixi.AbstractEquationsParabolic{1, 12, GradientVariablesConservative}
    equations_hyperbolic::Damped_Full_Hyperbolic{RealT}    # E
    dampingMat::Matrix{RealT}                              # T1
    flexMat::Matrix{RealT}                                 # T2
    massMat::Matrix{RealT}                                 # T3
    gammainverse::Matrix{RealT}                            # T4
    Psi::Matrix{RealT}                                     # T5
    k::Matrix{RealT}                                       # T6
    Ematrix::Matrix{RealT}                                 
    Gamma::Matrix{RealT}
    Pi::Matrix{RealT}
    flexMatInverse::Matrix{RealT}                                 # T2


    # Constructor
    function Damped_Full_Parabolic{RealT}(equations_hyperbolic, dampingMat, flexMat, massMat, 
        gammainverse, Psi, k) where {RealT<:Real}
        flexMatInverse = inv(flexMat)
        Ematrix = generate_E_matrix(k)
        Gamma = inv(gammainverse)
        Pi = -[zeros(6,6) diagm(ones(6));  diagm(ones(6)) zeros(6,6)]

        new(equations_hyperbolic, dampingMat, flexMat, massMat, gammainverse, Psi, k, Ematrix, Gamma, Pi, flexMatInverse)
    end
end



# Constructor function
Damped_Full_Parabolic(equations_hyperbolic::Damped_Full_Hyperbolic{RealT},
                                  dampingMat::Matrix{RealT},
                                  flexMat::Matrix{RealT},
                                  massMat::Matrix{RealT},
                                  gammainverse::Matrix{RealT},
                                  Psi::Matrix{RealT},
                                  k::Matrix{RealT}) where {RealT<:Real} = 
                                  Damped_Full_Parabolic{RealT}(equations_hyperbolic,
                                      dampingMat,
                                      flexMat,
                                      massMat,
                                      gammainverse,
                                      Psi,
                                      k)

function varnames(variable_mapping, equations_parabolic::Damped_Full_Parabolic)
    varnames(variable_mapping, equations_parabolic.equations_hyperbolic)
end


function flux(u, gradients, orientation::Integer, equations_parabolic::Damped_Full_Parabolic) 
    @unpack dampingMat, flexMat, Ematrix, flexMat, flexMatInverse = equations_parabolic
    L1matrix = generate_L1_matrix(u[1:6])
    q = [dampingMat * flexMatInverse * ( gradients[1:6] - Ematrix'*u[1:6] + L1matrix'*flexMat*u[7:12] ); zeros(6)]
    return q
end


function manufactured_soln(u, x, t, equations::Damped_Full_Parabolic, solution, solution_dot, solution_prime, solution_dot_prime, solution_prime_prime)
    @unpack Pi, Gamma, flexMat, massMat, dampingMat, Ematrix, flexMatInverse = equations

    u_exact = solution(x[1], t)
    u_dot = solution_dot(x[1], t)
    u_prime = solution_prime(x[1],t)
    u_prime_prime = solution_prime_prime(x[1],t)

    L1matrix = generate_L1_matrix(u_exact[1:6])
    L2matrix = generate_L2_matrix(u_exact[7:12])
    L1_prime_matrix = generate_L1_matrix(u_prime[1:6])

    rhat = dampingMat * flexMatInverse * (u_prime[1:6] - Ematrix' * u_exact[1:6] + L1matrix' * flexMat * u_exact[7:12])
    rhat_prime = dampingMat * flexMatInverse * (u_prime_prime[1:6] - Ematrix' * u_prime[1:6] + L1_prime_matrix' * flexMat * u_exact[7:12] + L1matrix' * flexMat * u_prime[7:12])
    L2rmatrix = generate_L2_matrix(rhat)

    Bmatrix = [zeros(6, 6) -Ematrix; Ematrix' zeros(6, 6)]
    j1matrix = [L1matrix * massMat zeros(6, 6); zeros(6, 6) zeros(6, 6)]
    j2matrix = [zeros(6, 6) L2matrix * flexMat; zeros(6, 6) -L1matrix' * flexMat]
    j3matrix = [zeros(6, 6) L2rmatrix * flexMat; zeros(6, 6) zeros(6, 6)]
    j4matrix = [-Ematrix * rhat; zeros(6, 1)]

    part_1 = (Bmatrix + j1matrix + j2matrix + j3matrix) * u_exact
    part_2 = j4matrix - [equations.equations_hyperbolic.f_ext(x[1], t); zeros(6, 1)]

    Q = part_1 + part_2
    return (Gamma * u_dot + Pi * u_prime + Pi * [zeros(6, 1); rhat_prime] + Q)

end