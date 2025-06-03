

# Linear scalar diffusion for use in linear scalar advection-diffusion problems
abstract type AbstractLaplaceDiffusion{NDIMS, NVARS} <:
              AbstractEquationsParabolic{NDIMS, NVARS, GradientVariablesConservative} end
include("laplace_diffusion_1d.jl")
include("laplace_diffusion_2d.jl")
include("laplace_diffusion_3d.jl")

include("laplace_diffusion_entropy_variables.jl")
include("laplace_diffusion_entropy_variables_1d.jl")
include("laplace_diffusion_entropy_variables_2d.jl")
include("laplace_diffusion_entropy_variables_3d.jl")

# Compressible Navier-Stokes equations
abstract type AbstractCompressibleNavierStokesDiffusion{NDIMS, NVARS, GradientVariables} <:
              AbstractEquationsParabolic{NDIMS, NVARS, GradientVariables} end
include("compressible_navier_stokes.jl")
include("compressible_navier_stokes_1d.jl")
include("compressible_navier_stokes_2d.jl")
include("compressible_navier_stokes_3d.jl")
