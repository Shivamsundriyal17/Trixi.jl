# Reference convergence data

`mms_convergence.csv` is the complete compatible exponential MMS campaign
generated from clean Trixi commit
`b5b1747b7a9f0a0fc34452d4222d35b0a26ac943` with Julia `1.11.6`, one thread,
and `OrdinaryDiffEqStabilizedRK 1.4.0`.  It uses
`lambda=2`, `q=(3,2,3,1,1,1)`, and `psi(x,t)=exp(x+t)`.  The exact fields
satisfy the physical lower compatibility equation without a manufactured
source while activating all twelve primary components, both geometric terms,
and all six Kelvin--Voigt resultants.

The campaign uses characteristic upwinding, polynomial degrees `1,2,3`,
uniform cell counts `4,8,16,32,64,128`, alternating LDG and BR1 auxiliary
traces, final time `T=1`, and ROCK4 absolute/relative tolerances `1e-12`.
There are 108 data rows: the two primary state blocks and reconstructed
`r_tau` block for each of 36 spatial runs.
`mms_convergence_componentwise.csv` records all twelve primary-component
`L2` and `Linf` errors for both trace closures, and
`mms_convergence_viscous_componentwise.csv` records the corresponding six
reconstructed-resultant errors.  `mms_temporal_sanity.csv` documents matched
`1e-12`/`1e-13` calculations on 64 and 128 cells.  In the tight cubic
alternating comparison, the primary-component `L2` EOCs range from `3.8063`
to `4.0693`, and the reconstructed-resultant EOCs range from `3.8101` to
`4.0227`.  Tightening the tolerance on 128 cells changes the displayed `u1`
and `u2` block errors by approximately `0.171%` and `0.0033%`, respectively;
the primary-state rates are therefore not limited by temporal error.

Finest-grid EOCs are:

| degree | auxiliary traces | `u1` | `u2` | reconstructed `r_tau` |
|---:|:---|---:|---:|---:|
| 1 | alternating | 1.8275 | 1.9751 | 1.6007 |
| 1 | BR1 | 2.0434 | 1.9783 | 1.0086 |
| 2 | alternating | 2.8771 | 2.8588 | 2.8031 |
| 2 | BR1 | 3.0185 | 3.0010 | 2.9232 |
| 3 | alternating | 3.9192 | 3.9603 | 3.5342 |
| 3 | BR1 | 4.0404 | 3.0883 | 2.9898 |

The `r_tau` entries in this table use the common `1e-12` campaign.  Its cubic
alternating rate is `3.8997` in the tighter `1e-13` comparison.  These
reconstructed-field rates are diagnostics; the paper does not assert an
`O(h^(k+1))` error theorem for `r_tau`, particularly for BR1.

Regenerate the file from the repository root with:

```bash
julia --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/run_mms_convergence.jl \
  examples/damped_intrinsic_beam/reference/mms_convergence.csv
```

The same driver writes the primary and viscous componentwise files through its
second and third output arguments.  The temporal reference combines matched
fine-grid runs with `MMS_TIME_TOL=1e-12` and `1e-13`.  The old homogeneous
`q_star` MMS and the compatible `lambda=3/2` coefficient guard remain
regression/reproducibility tests and are not the paper's displayed MMS.
The files `mms_convergence_tight_k3*.csv` contain the complete two-mesh tight
audit, including its computed EOCs.  The files `mms_lambda_1p5_guard*.csv`
record the independent coefficient guard on 16, 32, and 64 cells; it changes
the relation between the compatibility terms and guards specifically against
correlated coefficient errors that could cancel for `lambda=2`.

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

`rotating_steady_mesh_check.csv` records the steady-profile initialization
check on 4, 8, and 16 cells. It is regenerated with
`run_rotating_steady_mesh_check.jl`.
For the cubic discretization, the finest-pair `Linf` EOCs are `3.9868` for
`f1` and `3.7297` for `v2`; on 16 cells their relative `Linf` errors are
`6.15e-9` and `5.21e-10`. The maximum instantaneous ledger defect is at
roundoff on every mesh.

## Base-excited cantilever

The `farokhi_*_k4n2_*.csv` files are the raw degree-four/two-element
continuation campaigns for the 0.2g and 0.5g experimental comparison. The
up-sweep files intentionally retain post-fold or maximum-cycle rows; the
analysis script accepts only the explicitly documented frequency and
periodicity ranges. The two down-sweeps contain the strict lower stable
branches. `farokhi_02g_k3n4_refined_check.csv` is the independent
degree-three/four-element lower-branch check.
`farokhi_05g_k3n4_refined_up_sweep.csv` is the seven-point refined upper path
through the extreme 0.5g response.

`base_excited_cantilever_branch_comparison.csv` contains every accepted
numerical/experimental pair. `base_excited_cantilever_branch_summary.csv`
contains the branchwise errors, and
`base_excited_cantilever_energy_summary.csv` contains the final-cycle work
decomposition at each branch's maximum accepted numerical response. The raw
experimental markers are in
`farokhi_2022_experimental_frequency_response.csv`.
`base_excited_cantilever_refinement_comparison.csv` matches every refined
point to its baseline and experimental marker and records both response
changes and signed work fractions.

Regenerate the four baseline raw campaigns with:

```bash
for campaign in 02g_upper 02g_lower 05g_upper 05g_lower; do
  CANTILEVER_SWEEP_RESUME=false \
    examples/damped_intrinsic_beam/run_base_excited_cantilever_campaign.sh \
    "${campaign}"
done
```

Then regenerate the derived tables with:

```bash
python3 examples/damped_intrinsic_beam/analyze_base_excited_cantilever_sweeps.py \
  --results-directory examples/damped_intrinsic_beam/results \
  --output-directory examples/damped_intrinsic_beam/results
```

The exact frequency paths, cycle budgets, tolerances, and output names are
versioned in `run_base_excited_cantilever_campaign.sh`. Runs are restartable
both between frequencies and between additional settling blocks. See
`../BASE_EXCITED_CANTILEVER.md` for the normalization, marker extraction,
branch-selection rules, strict interpretation, and energy identity.

`base_excited_cantilever_validation.csv` and
`base_excited_cantilever_numerical_checks.csv` predate the full campaigns.
They are retained as legacy regression records and must not be used instead
of the branch tables for validation claims.
