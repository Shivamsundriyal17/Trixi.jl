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
  examples/damped_intrinsic_beam/elixir_mms_rich_physical.jl
```

Run the complete convergence campaign:

```bash
julia --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/run_mms_convergence.jl \
  examples/damped_intrinsic_beam/results/mms_convergence.csv
```

The committed [`reference/`](reference/) directory contains the compatible
exponential campaign generated from clean implementation commit `b5b1747b`
and promoted in reference-data commit `483e594b`. It includes the standard
campaign, tight cubic audit, reconstructed-resultant errors, and the
`lambda=1.5` coefficient guard.

The defaults are polynomial degrees `1,2,3`, cell counts
`4,8,16,32,64,128`, characteristic upwinding (`sigma=1`), both alternating
LDG and BR1 auxiliary traces, final time `T=1`, and absolute/relative ROCK4
tolerances `1e-12`. The example environment pins
`OrdinaryDiffEqStabilizedRK` to `1.4.0`. The CSV preamble records the Julia
version, Trixi commit and worktree state, command, and numerical settings.
The campaign also writes sibling `_componentwise.csv` and
`_viscous_componentwise.csv` files. The first contains `L2` and `Linf` errors
and EOCs for all twelve primary fields. The second reconstructs all six
Kelvin--Voigt resultant components from Trixi's actual discrete auxiliary
gradient and records the same norms and rates.

For a smaller smoke run, override comma-separated settings through environment
variables:

```bash
MMS_POLYDEGS=2 MMS_REFINEMENT_LEVELS=2,3 MMS_AUXILIARY_FLUXES=alternating \
  julia --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/run_mms_convergence.jl
```

Available variables are `MMS_POLYDEGS`, `MMS_REFINEMENT_LEVELS`,
`MMS_SIGMAS`, `MMS_AUXILIARY_FLUXES`, `MMS_TIME_TOL`, and `MMS_ELIXIR`.
Refinement level `r` means `2^r` uniform cells.

The paper-facing physical MMS activates all twelve state components, all six
Kelvin--Voigt resultant components, both lower compatibility terms, and exact
inhomogeneous split boundary data. Its default traveling profile is
`exp(x + t)`. Run a compact alternating-LDG audit with:

```bash
JULIA_NUM_THREADS=1 \
RICH_MMS_LAMBDA=2 \
MMS_POLYDEGS=1,2,3 MMS_REFINEMENT_LEVELS=2,3,4,5 \
MMS_AUXILIARY_FLUXES=alternating \
julia --compiled-modules=no \
  --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/run_mms_convergence.jl \
  examples/damped_intrinsic_beam/results/mms_rich_convergence.csv
```

Set `RICH_MMS_LAMBDA` to override the default value `2.0`; the independent
cancellation guard uses `1.5`. Set `RICH_MMS_PROFILE=quartic` to reproduce the
first rich-profile audit. The older homogeneous sine MMS remains available as
an independent regression with
`MMS_ELIXIR=examples/damped_intrinsic_beam/elixir_mms_physical.jl`.

For the full campaign, use refinement levels `2,3,4,5,6,7` and both
`alternating,br1` traces. With the default exponential profile, the alternating
`k=3` primary-state EOCs on the finest pair are approximately `3.92,3.96` for
`u1,u2`. A `1e-13` repeat at `N=64,128` gives resolved EOCs
`3.922,3.960,3.900` for `u1,u2,r_tau`. The quartic profile instead reaches the
temporal floor at the finest cubic mesh and is retained only as a supplementary
structural regression.

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

Once the `05g_refined_upper` campaign is complete, export and visualize its
final periodic orbit with:

```bash
julia --compiled-modules=no \
  --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/save_base_excited_cantilever_cycle.jl

python3 \
  examples/damped_intrinsic_beam/plot_base_excited_cantilever_cycle.py
```

The exporter writes centerline snapshots, tip histories, and a one-cycle
energy ledger. The plotter produces a static PNG/PDF and a supplementary GIF.
See [`BASE_EXCITED_CANTILEVER.md`](BASE_EXCITED_CANTILEVER.md) for the
refinement-diagnostic command and the exact archived state used.

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

## Pre-push validation and run integrity

Run the focused beam checks and a small matched MMS convergence/performance report:

```bash
JULIA_NUM_THREADS=1 julia --compiled-modules=no --check-bounds=yes \
  --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/run_pre_push_validation.jl
```

This runs constructor/precision checks, zero and rank-deficient damping energy
identities (periodic and split cantilever boundaries), explicit-step stability
checks, a non-beam diffusion regression, the paper's rich MMS, allocation checks,
and a degree-two MMS refinement on 4/8/16 cells at `T=0.1`. The CSV is written to
`results/pre_push_validation.csv`. Solve timings include compilation; the forcing
and ledger call timings are warmed measurements. Archived `reference/` data are
never regenerated by this command.

All example solves must succeed, reach the requested final time, and contain
finite saved states before exporting metrics. The explicit rotating and
nonsmooth examples now cap the wave-CFL step with the conservative estimate
`0.1*h_min^2/(norm(M^-1*H, Inf)*(p+1)^4)`. This checks the principal diffusion
restriction on the fixed mesh; it does not assert nonlinear stability for
arbitrary data. The original archived campaigns retain their recorded timestep
protocol and generating commits. New explicit runs can use smaller steps.

Rotating result caches carry their physical/numerical configuration and a hash
of executable sources and pinned environments. Reuse rejects incompatible or
legacy caches, including quick-mode results requested by a full campaign.
Regenerate rejected caches into `results/`. Cantilever checkpoints use format 5
and validate the complete saved configuration before resuming. Legacy formats
remain usable by the replay tools, but cannot be resumed without complete
provenance. Cycle replay restores the saved gravity multiplier and, for format
5, the integration tolerances, step cap, and Jacobian setting.

The rich MMS uses typed callable parameters created after `trixi_include`
overrides. Its forcing and boundary callbacks do not depend on mutable globals.
The coefficient constructor preserves Float32 inputs when curvature is omitted;
positive definiteness is checked without an absolute eigenvalue floor tied to
units. Beam algebra uses direct static-matrix entries and fixed-length tuples.

A separate integration smoke test generates/reuses the four short rotating
campaigns, rejects a quick/full cache mismatch, and replays a portable
zero-gravity legacy checkpoint:

```bash
JULIA_NUM_THREADS=1 julia --compiled-modules=no \
  --project=examples/damped_intrinsic_beam \
  test/test_damped_intrinsic_beam_campaigns.jl
```
