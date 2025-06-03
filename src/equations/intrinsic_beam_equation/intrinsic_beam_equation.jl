"""
    constitutive_laws(U, flexMat, massMat)
    
Applies constitutive laws and determines γ, κ, P, H on the basis of U=(F,M,V,Ω)
"""
function constitutive_laws(U, flexMat, massMat)
    if size(U, 1) != 12
        error("Try to reshape solution into a (12 x n)-matrix!")
    end

    F = U[1:3,:]
    M = U[4:6,:]
    V = U[7:9,:]
    Ω = U[10:12,:]

    v = flexMat * [F; M]
    w = massMat * [V; Ω]

    (γ, κ) = (v[1:3,:], v[4:6,:])
    (P, H) = (w[1:3,:], w[4:6,:])

    return γ, κ, P, H
end


"""
    constitutive_laws_kappa(U, flexMat, massMat)
    
Applies constitutive laws and determines κ on the basis of U=(F,M,V,Ω)
"""
function constitutive_laws_kappa(U, flexMat, massMat)
    if size(U, 1) != 12
        error("Try to reshape solution into a (12 x n)-matrix!")
    end

    F = U[1:3,:]
    M = U[4:6,:]

    return flexMat[4:6, 1:6] * [F; M]
end


struct IntrinsicBeamEquation{RealT<:Real} <: Trixi.AbstractEquations{1, 12}
    A::Matrix{RealT}
    A_plus::Matrix{RealT}
    A_minus::Matrix{RealT}
    Psi::Matrix{RealT}
    Pi::Matrix{RealT}
    Gamma::Matrix{RealT}
    Gamma_inv::Matrix{RealT}
    flexMat::Matrix{RealT}
    massMat::Matrix{RealT}
    k
    f_ext
    m_ext

    function IntrinsicBeamEquation{RealT}(flexMat, massMat, k, f_ext, m_ext) where {RealT<:Real}
        A = -[zeros(6,6) inv(flexMat); inv(massMat) zeros(6,6)]
        
        if isdiag(massMat) && isdiag(flexMat)
            lambda = sqrt.(flexMat*massMat)    
            T_inv = [sqrt.(flexMat) sqrt.(massMat); sqrt.(flexMat) -sqrt.(massMat)]
            T = inv(T_inv)
            Psi = inv(flexMat*massMat)
            D = [-inv(lambda) zeros(6,6); zeros(6,6) inv(lambda)]
        else
            lambda = diagm(sqrt.(eigvals(inv(flexMat*massMat))))
            Psi = flexMat^(-1/2)*inv(massMat)*flexMat^(-1/2)
            U_inv = LinearAlgebra.eigvecs(Psi)
            U = inv(U_inv)
            T = 0.5 * [flexMat^(-1/2)*inv(U) flexMat^(-1/2)*inv(U); flexMat^(1/2)*inv(U)*lambda -flexMat^(1/2)*inv(U)*lambda]
            T_inv = [U*flexMat^(1/2) inv(lambda)*U*flexMat^(-1/2); U*flexMat^(1/2) -inv(lambda)*U*flexMat^(-1/2)]
            D = [-lambda zeros(6,6); zeros(6,6) lambda]
        end

        
        A_plus = 0.5 * ( A + T*abs.(D)*T_inv)
        A_minus = 0.5 * ( A - T*abs.(D)*T_inv)

        Pi = -[zeros(6,6) diagm(ones(6));  diagm(ones(6)) zeros(6,6)]

        Gamma = [flexMat zeros(6,6); zeros(6,6) massMat]
        Gamma_inv = inv(Gamma)

        new(A, A_plus, A_minus, Psi, Pi, Gamma, Gamma_inv, flexMat, massMat, k, f_ext, m_ext)
    end
end
# constructor
IntrinsicBeamEquation(flexMat::Matrix{RealT}, massMat::Matrix{RealT}, k, f_ext, m_ext) where {RealT<:Real} = IntrinsicBeamEquation{RealT}(flexMat, massMat, k, f_ext, m_ext)


varnames(::Any, ::IntrinsicBeamEquation) = ("F1", "F2", "F3", "M1", "M2", "M3", "V1", "V2", "V3", "Ω1", "Ω2", "Ω3", )


"""
    function flux(u, orientation::Integer, equations::IntrinsicBeamEquation)

Computes flux of capacity form
"""
@inline function flux(u, orientation::Integer, equations::IntrinsicBeamEquation) 
    @unpack Pi = equations

    return Pi * u
end


"""
    flux_upwind(u_ll, u_rr, orientation, equations::IntrinsicBeamEquation)

Computes upwind numerical flux for capacity form of intrinsic beam equation
"""
function flux_upwind(u_ll, u_rr, orientation, equations::IntrinsicBeamEquation)
    @unpack A_plus, A_minus, Gamma = equations

    return Gamma * ( A_plus*u_ll + A_minus*u_rr )
end


entropy(u, ::IntrinsicBeamEquation) = 0.5*u.^2
cons2prim(u,::IntrinsicBeamEquation) = u
cons2entropy(u, ::IntrinsicBeamEquation) = u


wavespeeds(equations::IntrinsicBeamEquation) = abs.(eigvals(equations.A))
max_abs_speeds(u, equations::IntrinsicBeamEquation) = wavespeeds(equations)
max_abs_speed_naive(u_ll, u_rr, orientation::Integer, equations::IntrinsicBeamEquation) = maximum(wavespeeds(equatios))


"""
    initial_condition_constant(x, t, ::IntrinsicBeamEquation)

Returns zero initial condition
"""
function initial_condition_constant(x, t, ::IntrinsicBeamEquation)
    return SVector{12}(zeros(12))
end


"""
    boundary_conditions_ib(u_inner, orientation, direction, x, t, surface_flux_function, equations::IntrinsicBeamEquation)

Sets uo the consistent outer values to impose boundary conditions for velocities at left boundary and for strains at right boundary
and evaluates the boundary flux assuming zero velocities and strains
"""
function boundary_conditions_ib(u_inner, orientation, direction,
                                x, t, surface_flux_function,
                                equations::IntrinsicBeamEquation)
    @unpack A_plus, A_minus, Gamma, Psi, flexMat = equations
    if direction == 1
        y_B = zeros(6)
        u_boundary = [u_inner[1:6] + flexMat^(-1/2)*Psi^(-1/2)*flexMat^(-1/2)*(u_inner[7:12]-y_B); y_B]
        flux = Gamma * (A_plus * u_boundary + A_minus * u_inner)
    elseif direction == 2
        s_B = zeros(6)
        u_boundary = [s_B; u_inner[7:12]-flexMat^(1/2)*Psi^(1/2)*flexMat^(1/2)*(u_inner[1:6]-s_B)]
        flux = Gamma * (A_plus * u_inner + A_minus * u_boundary)
    end
    return flux
end


"""
    boundary_conditions_ib(u_inner, orientation, direction, x, t, surface_flux_function, equations::IntrinsicBeamEquation)

Sets uo the consistent outer values to impose boundary conditions for velocities at left boundary and for strains at right boundary
and evaluates the boundary flux
"""
function boundary_conditions_ib(u_inner, orientation, direction,
                                x, t, surface_flux_function,
                                velocities, strains,
                                equations::IntrinsicBeamEquation)
    @unpack A_plus, A_minus, Gamma, Psi, flexMat = equations
    if direction == 1
        y_B = velocities(t)
        u_boundary = [u_inner[1:6] + flexMat^(-1/2)*Psi^(-1/2)*flexMat^(-1/2)*(u_inner[7:12]-y_B); y_B]
        flux = Gamma * (A_plus * u_boundary + A_minus * u_inner)
    elseif direction == 2
        s_B = strains(t)
        u_boundary = [s_B; u_inner[7:12]-flexMat^(1/2)*Psi^(1/2)*flexMat^(1/2)*(u_inner[1:6]-s_B)]
        flux = Gamma * (A_plus * u_inner + A_minus * u_boundary)
    end
    return flux
end

"""
    source_term_intrinsic_beam(u, x, t, equation)

Computes source term of capacity form
    Γ u_t + Π u_x = Q_cap(u)
"""
function source_term_intrinsic_beam(u, x, t, equations::IntrinsicBeamEquation)
    @unpack A, flexMat, massMat, k, f_ext, m_ext = equations

    # U = [Th; Xi; V; Ω]
    Th = u[1:3]
    Xi = u[4:6]
    V = u[7:9]
    Ω = u[10:12]

    k_ev = k(x)
    f_ext_ev = f_ext(x,t)
    m_ext_ev = m_ext(x,t)

    # flexMat = [F1 F2; F2^T F3]
    F1 = flexMat[1:3, 1:3]
    F2 = flexMat[1:3, 4:6]
    F3 = flexMat[4:6, 4:6]

    # massMat = [M1 M2; M2^T M3]
    M1 = massMat[1:3, 1:3]
    M2 = massMat[1:3, 4:6]
    M3 = massMat[4:6, 4:6]

    # entries in the first summand of J(u)
    J1 = F2'*Th + F3*Xi
    J2 = F1*Th + F2*Xi

    M1V = M1 * V
    M2Ω = M2 * Ω

    e1 = [1.0;0.0;0.0]

    Q = [
        cross(k_ev+J1,V) + cross(e1+J2,Ω) ;
        cross(k_ev+J1,Ω) ;
        cross(k_ev+J1,Th) - cross(Ω,M1V) - cross(Ω, M2Ω) + f_ext_ev;
        cross(e1+J2,Th) + cross(k_ev+J1, Xi) - cross(V,M1V) - cross(Ω,M2'*V) - cross(V,M2Ω) - cross(Ω, M3*Ω) + m_ext_ev
        ]

    
    return Q
end


"""
    source_term_intrinsic_beam_diag(u, x, t, equation)

Berechnet Quellterm der zur intrinsic beam equation aequivalenten Bilanzgleichung unter der Annahme,
dass Flexibilitaets- und Massenmatrix diagonal sind und
    Γ u_t + Π u_x = Q_cap(u)
"""
function source_term_intrinsic_beam_diag(u, x, t, equations)
    @unpack A, flexMat, massMat, k, f_ext, m_ext = equations

    # U = [Th; Xi; V; Ω]
    Th = u[1:3]
    Xi = u[4:6]
    V = u[7:9]
    Ω = u[10:12]

    k_ev = k(x)
    f_ext_ev = f_ext(x,t)
    m_ext_ev = m_ext(x,t)

    # flexMat = [F1 0; 0 F3]
    F1 = flexMat[1:3, 1:3]
    F3 = flexMat[4:6, 4:6]

    # massMat = [M1 0; 0 M3]
    M1 = massMat[1:3, 1:3]
    M3 = massMat[4:6, 4:6]

    # entries in the first summand of J(u)
    J1 = F3 * Xi
    J2 = F1 * Th

    M1V = M1 * V

    e1 = [1.0;0.0;0.0]

    Q = [
        cross(k_ev+J1,V) + cross(e1+J2,Ω) ;
        cross(k_ev+J1,Ω) ;
        cross(k_ev+J1,Th) - cross(Ω,M1V) + f_ext_ev;
        cross(e1+J2,Th) + cross(k_ev+J1, Xi) - cross(V,M1V) - cross(Ω, M3*Ω) + m_ext_ev
        ]
    
    return Q
end


"""
    man_sol(u, x, t, equations, solution, solution_dot, solution_prime)

Computes additional source term that emerges by inserting a function into capacity form
"""
function man_sol(u, x, t, equations, solution, solution_dot, solution_prime)
    @unpack Pi, Gamma, flexMat, massMat, k, f_ext, m_ext = equations

    u_exact = solution(x[1], t)
    u_dot = solution_dot(x[1], t)
    u_prime = solution_prime(x[1],t)

    # U = [Th; Xi; V; Ω]
    Th = u_exact[1:3]
    Xi = u_exact[4:6]
    V = u_exact[7:9]
    Ω = u_exact[10:12]

    k_ev = k(x)
    f_ext_ev = f_ext(x,t)
    m_ext_ev = m_ext(x,t)

    # flexMat = [F1 F2; F2^T F3]
    F1 = flexMat[1:3, 1:3]
    F2 = flexMat[1:3, 4:6]
    F3 = flexMat[4:6, 4:6]

    # massMat = [M1 M2; M2^T M3]
    M1 = massMat[1:3, 1:3]
    M2 = massMat[1:3, 4:6]
    M3 = massMat[4:6, 4:6]

    # entries in the first summand of J(u)
    J1 = F2'*Th + F3*Xi
    J2 = F1*Th + F2*Xi

    M1V = M1 * V
    M2Ω = M2 * Ω

    e1 = [1.0;0.0;0.0]

    Q = [
        cross(k_ev+J1,V) + cross(e1+J2,Ω) ;
        cross(k_ev+J1,Ω) ;
        cross(k_ev+J1,Th) - cross(Ω,M1V) - cross(Ω, M2Ω) + f_ext_ev;
        cross(e1+J2,Th) + cross(k_ev+J1, Xi) - cross(V,M1V) - cross(Ω,M2'*V) - cross(V,M2Ω) - cross(Ω, M3*Ω) + m_ext_ev
        ]
    
    return ( Gamma*u_dot + Pi*u_prime - Q )
end

"""
    compute_energy(U, t, semi)

Computes total mechanical energy at points in time t for given solution U
"""
function compute_energy(U, t, semi)
    @unpack Gamma, flexMat, massMat = semi.equations

    Nx = length(semi.cache.elements._node_coordinates)

    Nt = length(t)

    # initialize and compute integrand at every node, for every timestep
    integrand_energy = Matrix{Float64}(undef, Nx, Nt)
    for i in 1:Nx
        for j in 1:Nt
            integrand_energy[i,j] = U[:,i,j]' * Gamma * U[:,i,j]
        end
    end

    # determine jacobian and quadrature weights
    jacobian = 1 ./ semi.cache.elements.inverse_jacobian
    weights = repeat(semi.solver.basis.weights, length(jacobian), 1)

    # initialize and compute energy by applying LGL formula for integrand at every timestep and scale with jacobian
    # note that equidistand spatial discrtization is assumed currently
    energy = Vector{Float64}(undef, Nt)
    for j in 1:Nt
        energy[j] = (jacobian[1] * weights' * integrand_energy[:,j])[1]
    end

    return energy
end