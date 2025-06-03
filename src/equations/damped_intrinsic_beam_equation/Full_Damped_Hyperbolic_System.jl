using LinearAlgebra
struct Damped_Full_Hyperbolic{RealT <: Real} <: Trixi.AbstractEquations{1, 12}
    A::Matrix{RealT}
    A_plus::Matrix{RealT}
    A_minus::Matrix{RealT}
    Psi::Matrix{RealT}
    sqrtm_Psi::Matrix{RealT}           
    sqrtm_inv_Psi::Matrix{RealT}        
    Pi::Matrix{RealT}
    Gamma::Matrix{RealT}
    Gamma_inv::Matrix{RealT}
    flexMat::Matrix{RealT}
    flexMat_sqrt::Matrix{RealT}         
    flexMat_inv_sqrt::Matrix{RealT}     
    massMat::Matrix{RealT}
    dampingMat::Matrix{RealT}
    k::Matrix{RealT}
    f_ext

    function Damped_Full_Hyperbolic{RealT}(flexMat, massMat, dampingMat, k, f_ext) where {RealT<:Real}

        A = -[zeros(6, 6) inv(massMat); inv(flexMat) zeros(6, 6)]
        flexMat_sqrt = (flexMat)^(1/2)
        flexMat_inv_sqrt = (flexMat)^(-1/2)
        if isdiag(massMat) && isdiag(flexMat)
            lambda = sqrt.(flexMat * massMat)
            T = [flexMat_sqrt flexMat_sqrt; flexMat_inv_sqrt -flexMat_inv_sqrt]
            T_inv = inv(T)
            Psi = inv(flexMat * massMat)
            D = [-inv(lambda) zeros(6, 6); zeros(6, 6) inv(lambda)]
        else
            lambda = diagm(sqrt.(eigvals(inv(flexMat * massMat))))
            Psi = flexMat_inv_sqrt * inv(massMat) * flexMat_inv_sqrt
            U_inv = LinearAlgebra.eigvecs(Psi)
            U = inv(U_inv)
            T = [flexMat_sqrt * inv(U) flexMat_sqrt * inv(U);
                 flexMat_inv_sqrt * inv(U) * inv(lambda) -flexMat_inv_sqrt * inv(U) * inv(lambda)]
            T_inv = 0.5 * [U * flexMat_inv_sqrt (lambda) * U * flexMat_sqrt;
                           U * flexMat_inv_sqrt -(lambda) * U * flexMat_sqrt]
            D = [-lambda zeros(6, 6); zeros(6, 6) lambda]
        end
        sqrtm_Psi = Psi^(1/2)
        sqrtm_inv_Psi = Psi^(-1/2)
        A_plus = 0.5 * (A + T * abs.(D) * T_inv)
        A_minus = 0.5 * (A - T * abs.(D) * T_inv)

        Pi = -[zeros(6, 6) diagm(ones(6)); diagm(ones(6)) zeros(6, 6)]

        Gamma = [massMat zeros(6, 6); zeros(6, 6) flexMat]
        Gamma_inv = inv(Gamma)

        new(A, A_plus, A_minus, Psi, sqrtm_Psi, sqrtm_inv_Psi, Pi, Gamma, Gamma_inv, flexMat, flexMat_sqrt, flexMat_inv_sqrt, massMat, dampingMat, k, f_ext)
    end
end


Damped_Full_Hyperbolic(flexMat::Matrix{RealT}, massMat::Matrix{RealT}, dampingMat::Matrix{RealT} ,k::Matrix{RealT}, f_ext) where {RealT<:Real} = Damped_Full_Hyperbolic{RealT}(flexMat, massMat, dampingMat, k, f_ext)



#Return the variable names for the `Damped_Full_Hyperbolic` system.
Trixi.varnames(::typeof(cons2prim), ::Damped_Full_Hyperbolic)=("V1", "V2", "V3", "ω1", "ω2", "ω3", "x2_1", "x2_2", "x2_3", "x2_4", "x2_5", "x2_6")
Trixi.varnames(::typeof(cons2cons), ::Damped_Full_Hyperbolic)=("V1", "V2", "V3", "ω1", "ω2", "ω3", "x2_1", "x2_2", "x2_3", "x2_4", "x2_5", "x2_6")

entropy(u, ::Damped_Full_Hyperbolic) = 0.5*u.^2
cons2prim(u,::Damped_Full_Hyperbolic) = u
cons2entropy(u, ::Damped_Full_Hyperbolic) = u

"""
    function flux(u, orientation::Integer, equations::Damped_Full_Hyperbolic)

Computes flux of capacity form
"""
@inline function flux(u, orientation::Integer, equations::Damped_Full_Hyperbolic) 
    @unpack Pi = equations
    return Pi * u
end

"""
    flux_upwind(u_ll, u_rr, orientation, equations::Damped_Full_Hyperbolic)

Computes upwind numerical flux for capacity form of Damped_Full_Hyperbolic beam equation
"""
function flux_upwind(u_ll, u_rr, orientation, equations::Damped_Full_Hyperbolic)
    @unpack A_plus, A_minus, Gamma = equations

    return Gamma * ( A_plus*u_ll + A_minus*u_rr )
end

wavespeeds(equations::Damped_Full_Hyperbolic) = abs.(eigvals(equations.A))
max_abs_speeds(u, equations::Damped_Full_Hyperbolic) = wavespeeds(equations)
max_abs_speed_naive(u_ll, u_rr, orientation::Integer, equations::Damped_Full_Hyperbolic) = maximum(wavespeeds(equatios))


"""
    boundary_conditions_ib(u_inner, orientation, direction, x, t, surface_flux_function, equations::Damped_Full_Hyperbolic)

Sets uo the consistent outer values to impose boundary conditions for velocities at left boundary and for strains at right boundary
and evaluates the boundary flux assuming zero velocities and strains
"""
function boundary_conditions_ib(u_inner, orientation, direction, x, t, surface_flux_function, equations::Damped_Full_Hyperbolic)
    @unpack A_plus, A_minus, Gamma, sqrtm_Psi, sqrtm_inv_Psi, flexMat_sqrt, flexMat_inv_sqrt = equations  # Use precomputed fields

    if direction == 1
        y_B = zeros(6)
        u_boundary = [y_B; u_inner[7:12] + flexMat_inv_sqrt * sqrtm_inv_Psi * flexMat_inv_sqrt * (u_inner[1:6] - y_B)]
        flux = Gamma * (A_plus * u_boundary + A_minus * u_inner)
    elseif direction == 2
        s_B = zeros(6)
        u_boundary = [u_inner[1:6] - flexMat_sqrt * sqrtm_Psi * flexMat_sqrt * (u_inner[7:12] - s_B); s_B]
        flux = Gamma * (A_plus * u_inner + A_minus * u_boundary)
    end

    return flux
end


function boundary_conditions_ib(u_inner, orientation, direction,
                                x, t, surface_flux_function, velocities, strains,
                                equations::Damped_Full_Hyperbolic)
    @unpack A_plus, A_minus, Gamma, sqrtm_Psi, sqrtm_inv_Psi, flexMat_sqrt, flexMat_inv_sqrt = equations  # Use precomputed fields

    if direction == 1
        y_B = velocities(t)
        u_boundary = [y_B; u_inner[7:12] + flexMat_inv_sqrt * sqrtm_inv_Psi * flexMat_inv_sqrt * (u_inner[1:6] - y_B)]
        flux = Gamma * (A_plus * u_boundary + A_minus * u_inner)
    elseif direction == 2
        s_B = strains(t)
        u_boundary = [u_inner[1:6] - flexMat_sqrt * sqrtm_Psi * flexMat_sqrt * (u_inner[7:12] - s_B); s_B]
        flux = Gamma * (A_plus * u_inner + A_minus * u_boundary)
    end

    return flux
end

