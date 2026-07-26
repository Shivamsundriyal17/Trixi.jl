# Damped intrinsic beam examples

These examples implement the constant-coefficient Kelvin--Voigt system used in
the accompanying DG analysis. They keep simulation, convergence analysis, and
visualization separate.

Instantiate the example environment once from the Trixi repository root:

```bash
julia --project=examples/damped_intrinsic_beam -e \
  'using Pkg; Pkg.instantiate()'
```

The committed manifest was generated with Julia `1.11.6` and records all
transitive versions. The local preference disables optional loop vectorization
for this example environment; it does not change the DG formulation.

On Julia 1.11, the current Trixi `main` branch may hit an unrelated manual
precompile assertion in its SciML callback declarations. If that occurs, add
`--compiled-modules=no` to the Julia commands below. This is slower but leaves
the numerical results unchanged.

Run the physical manufactured-solution case:

```bash
julia --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/elixir_mms_physical.jl
```

Run the complete convergence campaign:

```bash
julia --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/run_mms_convergence.jl \
  examples/damped_intrinsic_beam/results/mms_convergence.csv
```

The defaults are polynomial degrees `1,2,3`, cell counts
`4,8,16,32,64,128`, characteristic upwinding (`sigma=1`), both alternating
LDG and BR1 auxiliary traces, final time `T=1`, and absolute/relative ROCK4
tolerances `1e-12`. The example environment pins
`OrdinaryDiffEqStabilizedRK` to `1.4.0`. The CSV preamble records the Julia
version, Trixi commit and worktree state, command, and numerical settings.

For a smaller smoke run, override comma-separated settings through environment
variables:

```bash
MMS_POLYDEGS=2 MMS_REFINEMENT_LEVELS=2,3 MMS_AUXILIARY_FLUXES=alternating \
  julia --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/run_mms_convergence.jl
```

Available variables are `MMS_POLYDEGS`, `MMS_REFINEMENT_LEVELS`,
`MMS_SIGMAS`, `MMS_AUXILIARY_FLUXES`, and `MMS_TIME_TOL`. Refinement level
`r` means `2^r` uniform cells.

The ordinary Trixi `source_terms` callback belongs to the hyperbolic split and
is not reused by the Kelvin--Voigt operator. Physical distributed forces and
moments are passed as `external_force` when constructing
`DampedIntrinsicBeamEquations1D`.
