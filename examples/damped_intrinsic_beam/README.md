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

The committed full-run data and finest-grid EOC summary are in
[`reference/`](reference/).

The defaults are polynomial degrees `1,2,3`, cell counts
`4,8,16,32,64,128`, characteristic upwinding (`sigma=1`), both alternating
LDG and BR1 auxiliary traces, final time `T=1`, and absolute/relative ROCK4
tolerances `1e-12`. The example environment pins
`OrdinaryDiffEqStabilizedRK` to `1.4.0`. The CSV preamble records the Julia
version, Trixi commit and worktree state, command, and numerical settings.
The campaign also writes a sibling `_componentwise.csv` containing `L2` and
`Linf` errors and EOCs for all twelve primary fields.

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

Run the targeted finest-grid temporal-error check with:

```bash
JULIA_NUM_THREADS=1 julia --compiled-modules=no \
  --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/run_mms_temporal_sanity.jl
```

Its defaults are `k=3`, `N=128`, alternating LDG, and ROCK4 tolerances
`1e-8,1e-10,1e-12,1e-13`. The `MMS_TEMPORAL_*` environment variables
override these settings.

The two canonical nonsmooth/nonhomogeneous experiments are:

```bash
julia --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/elixir_rotating_beam.jl

julia --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/elixir_nonsmooth_resultant.jl
```

Both use characteristic upwinding, the alternating LDG traces
`hat(u1)=u1^-` and `hat(r_tau)=r_tau^+`, polynomial degree three, eight
uniform cells, and CFL `0.01`. The rotating case runs to `T=40`; the localized
resultant-jump case runs to `T=10`. They save solution states and scalar
diagnostics but do not load a plotting package or write figures.

Run the complete rotating-beam validation campaign with:

```bash
JULIA_NUM_THREADS=1 julia --compiled-modules=no \
  --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/run_rotating_beam_campaign.jl
```

The campaign runs the baseline and doubled-damping spin-ups, the singular
undamped comparison, and a steady-profile initialization check. It writes
serialized arrays and a `summary.csv` to
`results/rotating_beam_campaign/`. The baseline run accumulates material,
interface, boundary, physical-root, and SAT-data contributions at every
accepted time step. The saved-state analysis separately checks the
instantaneous semi-discrete ledger and errors against the analytic rotating
branch. The CSV retains every cumulative ledger channel, its absolute closure
residual, and a scale-normalized closure residual.

The reference campaign uses one Julia thread. Since this test has only eight
elements, additional threads add scheduling overhead without changing the
numerical protocol.

Use `ROTATING_CAMPAIGN_QUICK=true` for a four-case end-to-end smoke test.
Set `ROTATING_REUSE_RESULTS=true` to regenerate only the summary from existing
serialized results. Individual cases can be configured through
`ROTATING_DAMPING_MULTIPLIER`, `ROTATING_STEADY_INITIAL`,
`ROTATING_T_END`, `ROTATING_SAVE_COUNT`, `ROTATING_RECORD_LEDGER`,
`ROTATING_POLYDEG`, `ROTATING_REFINEMENT_LEVEL`, and `ROTATING_CFL`.

Run the inexpensive steady-initialized mesh check, which separates spatial
steady-branch drift from spin-up relaxation, with:

```bash
JULIA_NUM_THREADS=1 julia --compiled-modules=no \
  --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/run_rotating_steady_mesh_check.jl
```

The base-excited experimental cantilever benchmark has its own model mapping,
continuation protocols, quantitative assessment, and work-ledger discussion
in [`BASE_EXCITED_CANTILEVER.md`](BASE_EXCITED_CANTILEVER.md). For example,
run the 0.2g upper branch with:

```bash
examples/damped_intrinsic_beam/run_base_excited_cantilever_campaign.sh \
  02g_upper
```

The four baseline campaigns and their experimental comparisons are archived
under `reference/`; generated checkpoints and per-frequency states remain in
the ignored `results/` directory.

## Optional visualization

Exporters run the simulation environment and save plain Julia arrays:

```bash
julia --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/save_rotating_beam_data.jl

julia --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/save_nonsmooth_resultant_data.jl
```

Plotting has a separate environment:

```bash
julia --project=examples/damped_intrinsic_beam/visualization -e \
  'using Pkg; Pkg.instantiate()'

julia --project=examples/damped_intrinsic_beam/visualization \
  examples/damped_intrinsic_beam/visualization/plot_rotating_beam.jl

julia --project=examples/damped_intrinsic_beam/visualization \
  examples/damped_intrinsic_beam/visualization/plot_rotating_beam_validation.jl \
  examples/damped_intrinsic_beam/results/rotating_beam_campaign

julia --project=examples/damped_intrinsic_beam/visualization \
  examples/damped_intrinsic_beam/visualization/plot_nonsmooth_resultant.jl

julia --project=examples/damped_intrinsic_beam/visualization \
  examples/damped_intrinsic_beam/visualization/plot_mms_convergence.jl
```

Each exporter and plotter accepts an explicit input/output path through command
line arguments. Generated data and figures go to `results/` by default and are
ignored by Git.

The ordinary Trixi `source_terms` callback belongs to the hyperbolic split and
is not reused by the Kelvin--Voigt operator. Physical distributed forces and
moments are passed as `external_force` when constructing
`DampedIntrinsicBeamEquations1D`.
