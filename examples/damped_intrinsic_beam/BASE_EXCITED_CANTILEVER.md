# Base-excited experimental cantilever validation

This benchmark reproduces selected frequency-response measurements from
Farokhi, Xia, and Erturk, “Experimentally validated geometrically exact model
for extreme nonlinear motions of cantilevers,” *Nonlinear Dynamics* 107
(2022), 457–475
([DOI](https://doi.org/10.1007/s11071-021-07023-9)).

It is validation work only. No values from this directory have been copied
into the intrinsic-beam paper.

## Model mapping

The experiment is a vertical, initially straight 1095 spring-steel strip:

- \(E=200\) GPa, \(\rho=7800\) kg/m³;
- \(L=81.5\) mm, \(b=9\) mm, \(h=0.0762\) mm;
- horizontal base acceleration of 0.2g or 0.5g RMS;
- gravity and weak-axis Kelvin–Voigt damping \(\eta_d=0.0037\).

The elixir uses \(L=1\), transverse mass per length \(=1\), and weak-axis
bending stiffness \(=1\). The dimensional data reproduce the article's
\(\gamma=0.4280\), \(\chi=7.28\times10^{-8}\), and excitation amplitudes
\(a_z=0.1211\) and \(0.3027\).

The prescribed root velocity represents base excitation in an inertial frame.
Gravity is a dead load. Since its material components depend on cross-section
orientation, the benchmark reconstructs angle by integrating intrinsic
curvature during every residual evaluation and adds the rotated gravitational
source. The initial axial resultant is the straight compressed equilibrium;
its residual is \(4.23\times10^{-12}\).

The Trixi beam is shearable and extensible, while the reference model is
Euler–Bernoulli and inextensible. Physical axial and shear stiffnesses are
retained. At the validated large-amplitude points, maximum axial and shear
strains remain below \(7.1\times10^{-7}\) and \(1.7\times10^{-6}\),
respectively, so the intended slender limit is actually reached rather than
assumed.

## Running

Instantiate the example environment, then run a single fixed-frequency case:

```bash
julia --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/elixir_base_excited_cantilever.jl
```

The article's horizontal coordinate is excitation frequency normalized by the
system natural frequency. The campaign therefore normalizes by the first
eigenfrequency of the discrete, gravity-loaded straight equilibrium. The
unloaded Euler--Bernoulli value is retained in the output as a separate
diagnostic; it is not used for matching the experimental markers.

Run any one of the archived protocols by name:

```bash
examples/damped_intrinsic_beam/run_base_excited_cantilever_campaign.sh \
  02g_upper
```

The four baseline campaign names are `02g_upper`, `02g_lower`, `05g_upper`,
and `05g_lower`. The lower branches are deliberately swept from high to low
frequency. Two targeted refinement checks are named `02g_refined_lower` and
`05g_refined_upper`. A controlled gravity sensitivity run is named
`05g_zero_gravity`. The wrapper contains the exact frequency paths, cycle
budgets, tolerances, discretizations, and output names used for the archived
results.

Each accepted point writes its CSV row and continuation state atomically.
The current frequency is also checkpointed before every extra settling block,
so an interrupted expensive run repeats at most one block. Per-frequency state
archives make local fold refinement possible without retracing the entire
branch. Set `CANTILEVER_SWEEP_RESUME=false` to start a campaign from its
initial equilibrium instead of its checkpoint.

The CSV `accepted_steps` and `rejected_steps` fields describe the final solve
block from which the response and ledger are measured; they are not cumulative
work counters over all settling blocks.

The baseline uses two elements of degree four. This is the least expensive
tested discretization that agrees with four degree-three elements at the
checked response points. An isolated fixed-frequency start in the coexistence
region can converge to either attractor; branch direction and continuation
state are therefore part of the numerical protocol.

On the current Julia 1.11/SciML combination, this Trixi branch has an unrelated
manual-precompile assertion. Until that upstream compatibility issue is fixed,
add `--compiled-modules=no` to the commands above. The example project targets
Julia 1.10, but this benchmark has not yet been rerun there.

Important controls include:

- `CANTILEVER_SWEEP_NORMALIZED_FREQUENCIES`;
- `CANTILEVER_GRAVITY_MULTIPLIER`;
- `CANTILEVER_SWEEP_RELTOL` and `CANTILEVER_SWEEP_ABSTOL`;
- `CANTILEVER_SWEEP_PERIODICITY_TOLERANCE`;
- `CANTILEVER_POLYDEG` and `CANTILEVER_REFINEMENT_LEVEL`;
- `CANTILEVER_SWEEP_RESUME` and `CANTILEVER_SWEEP_ARCHIVE_STATES`.

## Gravity implementation audit

Gravity is not applied as constant material-frame components. At each residual
evaluation, the implementation reconstructs the current cross-section angle
from the intrinsic curvature and rotates the fixed spatial dead load into the
material frame. Two targeted checks exercise this path.

First run the named refined upper-branch campaign, then release its final
large-deformation state with the root fixed and compare gravitational power
with the directional derivative of the reconstructed gravitational potential:

```bash
examples/damped_intrinsic_beam/run_base_excited_cantilever_campaign.sh \
  05g_refined_upper
julia --compiled-modules=no \
  --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/audit_base_excited_cantilever_gravity.jl
```

For the recorded degree-three/four-element state, the largest relative
work--potential residual over the released transient is
\(6.35\times10^{-5}\). This is a semidiscrete compatibility audit rather than
an ODE finite-difference check: the potential is differentiated analytically
along the computed semidiscrete right-hand side.

The second check repeats the high-frequency part of the 0.5g up-sweep with
gravity disabled, preserving the discretization, continuation protocol, and
normalized frequency path:

```bash
examples/damped_intrinsic_beam/run_base_excited_cantilever_campaign.sh \
  05g_zero_gravity
python3 \
  examples/damped_intrinsic_beam/analyze_base_excited_cantilever_gravity.py
```

Gravity lowers the discrete linear frequency from 9.383789 Hz to 9.124964 Hz
(2.758%). Relative to each configuration's own first frequency, the last
accepted upper-branch point moves from \(f/f_1=1.030893\) without gravity to
\(1.041175\) with gravity. In absolute frequency it instead moves from
9.673681 Hz to 9.500685 Hz because the linear-frequency reduction dominates.
The comparison therefore supports both effects reported for this benchmark:
gravity lowers the resonance scale while modestly increasing the normalized
hardening range.

After all four baseline campaigns are present, reproduce the comparison tables
with:

```bash
python3 examples/damped_intrinsic_beam/analyze_base_excited_cantilever_sweeps.py
```

Generate a four-panel inspection plot from the archived comparison with:

```bash
python3 \
  examples/damped_intrinsic_beam/plot_base_excited_cantilever_comparison.py
```

## Extreme-cycle visualization

After completing the `05g_refined_upper` campaign, export one measured cycle
from its final accepted state:

```bash
julia --compiled-modules=no \
  --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/save_base_excited_cantilever_cycle.jl
```

The exporter reads the archived state at normalized frequency 1.04118,
reintegrates exactly one period at 97 uniformly spaced phases, and writes
plain centerline, tip-response, and cycle-ledger CSV files. It includes the
cantilever elixir in setup-only mode, so it does not perform a discarded
default-frequency solve first.

Generate the paper-ready static figure and supplementary GIF with:

```bash
python3 \
  examples/damped_intrinsic_beam/plot_base_excited_cantilever_cycle.py
```

Generate the spatial-refinement and work-partition diagnostic with:

```bash
python3 \
  examples/damped_intrinsic_beam/plot_base_excited_cantilever_refinement.py
```

All outputs go to the ignored `results/` directory. The static cycle figure
combines centerline snapshots, the tip orbit, and periodic tip histories. The
refinement figure is deliberately separate because a nearly converged
transverse response does not imply a grid-independent split between physical
and numerical dissipation.

## Experimental marker extraction

The article does not provide a raw frequency-response table. The archived
experimental CSV was extracted from the vector circle markers in Figures 4
and 7, not manually read from raster pixels:

```bash
python3 examples/damped_intrinsic_beam/extract_farokhi_2022_frequency_response.py \
  /path/to/Farokhi2022.pdf \
  /tmp/farokhi_2022_frequency_response.csv
```

The PDF is not redistributed. Figure 4(b) has displayed limits
\([-0.65,0.05]\), rather than tick-aligned limits \([-0.7,0]\). The extractor
therefore calibrates the longitudinal coordinate from the labeled 0 and -0.6
ticks. Treating the plot-box edges as those tick values introduces a spurious
\(-0.05L\) offset in every 0.2g longitudinal marker. Marker uncertainty still
includes plot-axis rounding, line thickness, and uncertainty in the original
image-processing experiment.

## Current strict assessment

The table reports RMSE normalized by the largest experimental magnitude on
each accepted branch. An upper marker is the larger response at a duplicated
frequency; a lower marker is the smaller one. Rows that exceed the campaign's
periodicity tolerance, including the post-fold jump in the 0.5g up-sweep, are
not counted.

| acceleration | branch | points | transverse NRMSE | longitudinal NRMSE | max periodicity defect |
|:---:|:---:|---:|---:|---:|---:|
| 0.2g | upper | 16 | 2.20% | 2.79% | \(2.00\times10^{-3}\) |
| 0.2g | lower | 14 | 17.53% | 33.57% | \(1.81\times10^{-4}\) |
| 0.5g | upper | 18 | 1.86% | 3.95% | \(4.71\times10^{-4}\) |
| 0.5g | lower | 14 | 9.09% | 17.15% | \(1.92\times10^{-4}\) |

The upper branches are genuinely strong validation results, not merely
selected-point agreement. At 0.2g and normalized frequency 1.01855, for
example, the numerical and experimental transverse amplitudes are 0.76290
and 0.76221; the longitudinal minima are -0.47702 and -0.47714. At 0.5g,
the numerical upper branch is within 0.2% transversely at normalized
frequency 1.04118. The longitudinal error there is 2.9%.

The lower branches show a systematic amplitude bias. At 0.2g and normalized
frequency 1.04140, the degree-four/two-element transverse amplitude is
0.18941, versus 0.13905 experimentally. Repeating that lower-branch point
with degree three on four elements gives 0.18925 and changes the longitudinal
magnitude by only 0.17%. Thus spatial under-resolution does not explain this
particular disagreement; damping-model, clamp, and specimen-model differences
are more plausible hypotheses. They are not identified parameters and should
not be tuned after seeing the validation data.

At 0.5g the numerical up-sweep remains on the upper branch through normalized
frequency 1.04118, then jumps before an accepted solution is obtained at
1.04323. The experimental upper markers continue to approximately 1.0453.
This brackets a discrepancy in the fold region, but it does not locate the
numerical saddle-node: smaller steps restarted from the archived 1.04118
state are required before making a bifurcation claim.

The complete degree-three/four-element upper-branch check covers seven
frequencies from 0.97742 to 1.04118. At the final point, refinement changes
the transverse response from 0.85383 to 0.85706 (0.38%), compared with the
experimental 0.85534. It changes the longitudinal minimum from -0.99826 to
-1.03726 (3.76%), compared with the experimental -0.96978. Thus the extreme
transverse plateau is spatially robust, but the shortening and rotation are
more sensitive and the refined shortening is about 7.0% too large.

## Cycle energy ledger

For the final measured cycle, the code checks

\[
 \Delta E + D_{\mathrm{KV}} + D_{\mathrm{jump}} + D_L + D_R
 - W_{\mathrm{root}} - W_{\mathrm{SAT}} = r.
\]

Here \(E\) includes intrinsic and gravitational potential energy,
\(D_{\mathrm{KV}}\) is physical Kelvin--Voigt dissipation, and
\(D_{\mathrm{jump}}\) is upwind interface dissipation. Since the imposed
root velocity enters through an SAT boundary state, the interpretable signed
boundary contribution is \(D_L+D_R-W_{\mathrm{SAT}}\); positive values are
dissipative and negative values inject energy. Reporting \(D_L\) or
\(W_{\mathrm{SAT}}\) alone exposes a large but artificial cancellation.

At the maximum accepted upper response of each campaign, the percentages of
physical root work are:

| acceleration/discretization | \(\Delta E\) | Kelvin--Voigt | interface | net boundary |
|:---|---:|---:|---:|---:|
| 0.2g, degree 4 / 2 cells | 0.58% | 94.90% | 1.89% | 2.63% |
| 0.5g, degree 4 / 2 cells | 0.013% | 77.76% | 8.17% | 14.05% |
| 0.5g, degree 3 / 4 cells | 0.024% | 91.74% | 3.91% | 4.32% |

The compact residual divided by physical root work is below
\(1.6\times10^{-6}\) for the baseline peak rows and
\(4.8\times10^{-6}\) for the refined 0.5g peak. Closure is therefore
excellent, but closure alone does not make the numerical losses physical.
At the 0.5g extreme response, refinement reduces interface plus net boundary
loss from 22.22% to 8.24% while changing transverse amplitude by only 0.38%.
The baseline work partition is therefore not grid-independent. The refined
ledger is a much better diagnostic, but a precise physical/numerical
dissipation split should not be claimed without another refinement level.

The old selected-point file is retained only as a legacy regression record.
Its rows used shorter settling runs, and its original 0.2g longitudinal
markers contained the Figure 4(b) axis-offset error. Use the branch comparison
and summary files for scientific conclusions.

Archived files:

- `reference/farokhi_2022_experimental_frequency_response.csv`;
- `reference/farokhi_02g_k4n2_sweep.csv`;
- `reference/farokhi_02g_k4n2_down_sweep.csv`;
- `reference/farokhi_05g_k4n2_up_sweep.csv`;
- `reference/farokhi_05g_k4n2_down_sweep.csv`;
- `reference/base_excited_cantilever_branch_comparison.csv`;
- `reference/base_excited_cantilever_branch_summary.csv`;
- `reference/base_excited_cantilever_energy_summary.csv`;
- `reference/base_excited_cantilever_refinement_comparison.csv`;
- `reference/farokhi_02g_k3n4_refined_check.csv`;
- `reference/farokhi_05g_k3n4_refined_up_sweep.csv`;
- `reference/base_excited_cantilever_validation.csv` (legacy);
- `reference/base_excited_cantilever_numerical_checks.csv` (legacy).

### Checkpoint compatibility after the pre-push cleanup

New sweeps write format-5 checkpoints with source/environment identity and the
physical, integration, sampling, and continuation settings. Resume rejects a
changed configuration; legacy checkpoints remain available for cycle replay.
The cycle exporter restores saved gravity, including zero gravity, and restores
format-5 integration settings. The unit-gravity work audit rejects a checkpoint
from a different gravity configuration. Direct sweep output now defaults to
`results/`, so running an exploratory sweep does not replace archived references.
