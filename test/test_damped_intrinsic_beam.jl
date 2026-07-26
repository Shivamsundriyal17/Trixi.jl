module TestDampedIntrinsicBeam

using LinearAlgebra: Diagonal
using Test
using Trixi

include("test_trixi.jl")

@testset "Damped intrinsic beam" begin
    flexibility = Matrix(Diagonal(1.0:6.0))
    mass = Matrix(Diagonal(2.0:7.0))
    damping = 0.1 * flexibility
    equations = DampedIntrinsicBeamEquations1D(;
                                               mass_matrix = mass,
                                               flexibility_matrix = flexibility,
                                               damping_matrix = damping,
                                               initial_curvature = zeros(3))
    equations_parabolic = DampedIntrinsicBeamDiffusion1D(equations)

    @test equations.propagation_matrix_plus +
          equations.propagation_matrix_minus ≈ equations.propagation_matrix
    @test equations.capacity_matrix * equations.propagation_matrix_abs ≈
          transpose(equations.capacity_matrix *
                    equations.propagation_matrix_abs)

    u = SVector{12}(ntuple(index -> 0.03 * index, 12))
    gradient = SVector{12}(ntuple(index -> 0.02 - 0.01 * index, 12))
    @test flux_upwind(u, u, 1, equations) ≈ flux(u, 1, equations)
    @test entropy(u, equations) isa Number
    damping_resultant = intrinsic_beam_damping_resultant(u, gradient,
                                                         equations_parabolic)
    @test flux(u, gradient, 1, equations_parabolic)[1:6] ≈
          equations.mass_inverse * damping_resultant

    u1 = SVector{6}(u[1:6])
    u2 = SVector{6}(u[7:12])
    @test transpose(intrinsic_beam_l1(u1)) * u1 ≈ zeros(6)
    @test intrinsic_beam_l1(u1) * u2 ≈
          transpose(intrinsic_beam_l2(u2)) * u1

    zero66 = zeros(6, 6)
    geometry = equations.geometry_matrix
    l1 = intrinsic_beam_l1(u1)
    l2 = intrinsic_beam_l2(u2)
    l2_damping = intrinsic_beam_l2(damping_resultant)
    source_matrix = [zero66 geometry; -transpose(geometry) zero66] +
                    [-l1*equations.mass_matrix zero66; zero66 zero66] +
                    [zero66 -l2*equations.flexibility_matrix;
                     zero66 transpose(l1)*equations.flexibility_matrix] +
                    [zero66 -l2_damping*equations.flexibility_matrix;
                     zero66 zero66]
    source_reference = source_matrix * u +
                       [geometry * damping_resultant; zeros(6)]
    source_optimized = intrinsic_beam_gradient_source(u, gradient,
                                                      SVector(0.2), 0.3,
                                                      equations_parabolic)
    @test equations.capacity_matrix * source_optimized ≈ source_reference

    nonsymmetric_mass = copy(mass)
    nonsymmetric_mass[1, 2] = 1.0
    @test_throws ArgumentError DampedIntrinsicBeamEquations1D(;
                                                              mass_matrix = nonsymmetric_mass,
                                                              flexibility_matrix = flexibility,
                                                              damping_matrix = damping)

    incompatible_damping = copy(damping)
    incompatible_damping[1, 1] = 0.2
    incompatible_damping[1, 2] = incompatible_damping[2, 1] = 0.05
    @test_throws ArgumentError DampedIntrinsicBeamEquations1D(;
                                                              mass_matrix = mass,
                                                              flexibility_matrix = flexibility,
                                                              damping_matrix = incompatible_damping)

    zero_field = (x, t) -> zero(SVector{12, Float64})
    inconsistent_time_derivative = (x, t) -> SVector{12}(zeros(6)..., ones(6)...)
    @test_throws ArgumentError manufactured_force_damped_intrinsic_beam(0.0, 0.0,
                                                                        equations_parabolic,
                                                                        zero_field,
                                                                        inconsistent_time_derivative,
                                                                        zero_field,
                                                                        zero_field)
end

@trixi_testset "TreeMesh1D: elixir_mms_physical.jl" begin
    @test_trixi_include(joinpath(examples_dir(), "damped_intrinsic_beam",
                                 "elixir_mms_physical.jl"),
                        l2=[
                            4.234853472065755e-7,
                            1.0893933103973503e-6,
                            1.2460689375612526e-6,
                            1.5130404051980891e-6,
                            1.1274766408530033e-7,
                            1.7302794606238937e-7,
                            9.280790475106985e-7,
                            1.1244731465306087e-6,
                            2.9672936132746254e-6,
                            2.9528962220465824e-6,
                            1.3554102630853058e-7,
                            1.1784307340513368e-7
                        ],
                        linf=[
                            2.6087661337825807e-6,
                            6.611107070453315e-6,
                            7.402831811553767e-6,
                            9.48275078060945e-6,
                            5.112532626837479e-7,
                            8.432929719967199e-7,
                            2.774635358429589e-6,
                            5.672078576351991e-6,
                            7.940434284009479e-6,
                            8.423702675147693e-6,
                            3.068489687584053e-7,
                            6.147217452116827e-7
                        ])
end

@trixi_testset "TreeMesh1D: elixir_rotating_beam.jl" begin
    @test_trixi_include(joinpath(examples_dir(), "damped_intrinsic_beam",
                                 "elixir_rotating_beam.jl"),
                        tspan=(0.0, 1.0e-4),
                        save_times=[0.0, 1.0e-4])
end

@trixi_testset "TreeMesh1D: elixir_nonsmooth_resultant.jl" begin
    @test_trixi_include(joinpath(examples_dir(), "damped_intrinsic_beam",
                                 "elixir_nonsmooth_resultant.jl"),
                        tspan=(0.0, 1.0e-4),
                        save_times=[0.0, 1.0e-4])
end

end # module
