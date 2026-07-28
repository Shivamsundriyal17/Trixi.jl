# Reference convergence data

`mms_convergence.csv` is the complete physical MMS campaign generated from
clean Trixi commit `f4acef5a259a397aa22fcea93ac9c3e96d2ece7b` with Julia
`1.11.6`, one thread, and `OrdinaryDiffEqStabilizedRK 1.4.0`.

The campaign uses characteristic upwinding, polynomial degrees `1,2,3`,
uniform cell counts `4,8,16,32,64,128`, alternating LDG and BR1 auxiliary
traces, final time `T=1`, and ROCK4 absolute/relative tolerances `1e-12`.
There are 72 data rows: two state blocks for each of 36 spatial runs.

Finest-grid EOCs are:

| degree | auxiliary traces | `u1` | `u2` |
|---:|:---|---:|---:|
| 1 | alternating | 1.8214 | 1.9807 |
| 1 | BR1 | 1.9976 | 1.9824 |
| 2 | alternating | 2.9121 | 2.8890 |
| 2 | BR1 | 2.9903 | 2.9916 |
| 3 | alternating | 3.9303 | 3.9518 |
| 3 | BR1 | 3.9924 | 3.0652 |

Regenerate the file from the repository root with:

```bash
julia --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/run_mms_convergence.jl \
  examples/damped_intrinsic_beam/reference/mms_convergence.csv
```

`rotating_beam_validation.csv` is the quantitative rotating-beam campaign
used by the paper. Its metadata identifies the exact clean commit and Julia
thread count. `NaN` in a cumulative column means that online accumulation was
deliberately disabled for that comparison; instantaneous ledger defects were
still evaluated at every saved state.

Regenerate the rotating reference with:

```bash
JULIA_NUM_THREADS=1 julia --compiled-modules=no \
  --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/run_rotating_beam_campaign.jl \
  examples/damped_intrinsic_beam/results/rotating_beam_campaign \
  examples/damped_intrinsic_beam/reference/rotating_beam_validation.csv
```
