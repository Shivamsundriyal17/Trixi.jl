# Beam pre-push validation

The cleanup preserves the archived paper reference CSVs. It adds completion
checks, validated cache reuse and checkpoint resume, gravity-aware replay,
precision/scale-safe material construction, a conservative explicit diffusion
step cap, typed rich-MMS callbacks, and allocation-free pointwise beam algebra.
Existing uncommitted gravity work was retained and integrated.

## Environment and commands

Julia 1.11.6, one Julia thread, the pinned beam Manifest and LocalPreferences
(`loop_vectorization=false`). CPU only; CPU/CUDA parity is not applicable to
this implementation. The local DG CUDA skills target a different CMake solver
and cannot be used as Trixi validation.

The existing beam suite plus new regressions was run with:

```bash
JULIA_NUM_THREADS=1 julia --compiled-modules=no --check-bounds=yes \
  --project=/tmp/beam-prepush-test-env test/test_damped_intrinsic_beam.jl
```

The temporary test environment copies the beam Project/Manifest/LocalPreferences,
uses an absolute local Trixi source path, and adds TrixiTest 0.1.7 and
JuliaFormatter 1.0.60. The full beam suite passed **94 assertions**: 13 original
algebra checks, 27 old-MMS checks, two original example smoke checks, and 52 new
constructor, energy, completion/cache, explicit-step, scalar-diffusion, and rich-MMS
checks. The complete Trixi repository test matrix was not run.

The checked-in numerical/performance runner is:

```bash
JULIA_NUM_THREADS=1 julia --compiled-modules=no --check-bounds=yes \
  --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/run_pre_push_validation.jl
```

The portable campaign integration test can be run with:

```bash
JULIA_NUM_THREADS=1 julia --compiled-modules=no \
  --project=examples/damped_intrinsic_beam \
  test/test_damped_intrinsic_beam_campaigns.jl \
  examples/damped_intrinsic_beam/results/pre_push_campaigns
```

The fresh-process campaign integration run passed **11 assertions**, exit code 0:
rotating generation/reuse and quick/full mismatch rejection (3; 171.5 s),
format-5 cantilever resume and changed-tolerance rejection (5; 94.5 s), and
legacy zero-gravity checkpoint replay (3; 44.1 s). These test-set timings include
compilation but exclude initial package loading. The durable local log is
`results/pre_push_campaigns.txt`.

The two-cycle test exposed an omitted initial sample when measurement begins
at time zero. Both cantilever solve paths now save that initial sample, and
periodicity measurement requires at least two cycles. The integration test
passes with this fix. A separate earlier check also replayed the existing
legacy zero-gravity checkpoint
`gravity_sensitivity_05g_zero_k4n2_checkpoint_point_002.jls`; all three assertions
passed.

Additional checks: all six Python scripts parsed, the campaign shell passed
`bash -n`, changed Julia files were formatted with the repository's formatter
version, and `git diff --check` passed. An independent comparison against the
saved pre-edit algebra passed 120 exact Float32/Float64 matrix comparisons.

## Numerical results

Matched rich exponential MMS: degree 2, cells 4/8/16, final time 0.1, ROCK4
absolute/relative tolerance 1e-12. Values are the maximum componentwise norms.
These coarse-mesh EOCs are comparisons with the unchanged baseline, not a new
claim of asymptotic order. Full paper EOC campaigns were not regenerated.

| Cells | Baseline L2 | Baseline Linf |
|---:|---:|---:|
| 4 | 8.978808136296493e-4 | 2.9557718801171973e-3 |
| 8 | 1.6118473879984626e-4 | 5.820597618919621e-4 |
| 16 | 2.6333893081516855e-5 | 7.861423492894204e-5 |

Post-cleanup results from `results/pre_push_validation.csv`:

| Cells | L2 | Linf | L2 EOC | Linf EOC | Runtime (s, including JIT) |
|---:|---:|---:|---:|---:|---:|
| 4 | 8.978808136301955e-4 | 2.9557718801269672e-3 | — | — | 60.8541 |
| 8 | 1.6118473879977017e-4 | 5.820597618786394e-4 | 2.477809 | 2.344296 | 26.9163 |
| 16 | 2.6333893081154575e-5 | 7.861423493693565e-5 | 2.613722 | 2.888305 | 26.3503 |

Errors agree with the baseline to roundoff. Baseline solve runtimes were
68.0251, 27.2799, and 28.5857 seconds; these include JIT compilation.
The scalar non-beam regression gives L2
`9.3738639138675e-5` and Linf `3.164960186754426e-4` at T=0.01, degree 3,
eight cells, using the analytic periodic advection-diffusion solution.

The assembled energy tests cover periodic and homogeneous split cantilever
closures with zero and rank-deficient damping, tolerance 1e-11. The explicit
step tests check the Carpenter–Kennedy amplification factors of assembled
linearized beam operators for degrees 1/2/3, multiple mesh sizes, and damping
0/0.1/10; the refinement check verifies the h-squared diffusion restriction.
This is not a general nonlinear stability proof.

## Performance interpretation

Baseline warmed measurements (100 calls, minimum of three timing samples):
rich forcing 38,368 bytes and 12.054 microseconds/call; rotating ledger
225,136 bytes and 42.10218 microseconds/call. Core-source and rich-force
allocation tests require zero allocations in typed callers. The final runner
records the following post-cleanup measurements on an Intel Core i9-14900K,
using one Julia thread:

| Operation | Baseline bytes/call | Final bytes/call | Baseline µs/call | Final µs/call |
|---|---:|---:|---:|---:|
| Rich MMS forcing | 38,368 | 0 | 12.054 | 0.20153 |
| Rotating ledger (8 cells) | 225,136 | 3,408 | 42.10218 | 8.85821 |

The isolated L1/L2 constructors also dropped from 624 bytes each to zero.

MMS solve timings include JIT compilation and should not be read as end-to-end
speedup measurements. The smaller explicit timestep cap can increase the number
of explicit steps. No overall campaign speedup is claimed.

## Compatibility notes

New cantilever checkpoints use format 5 and reject resume with changed numerical
controls or source/environment identity. Legacy checkpoints remain replayable;
they cannot supply the missing provenance required for resume. Rotating caches
from older code must be regenerated before reuse. Full-duration rotating and
experimental cantilever sweeps were not rerun. The existing cantilever Jacobian
sparsity discovery remains a finite-difference heuristic.
