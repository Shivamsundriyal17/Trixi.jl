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

Trace a stable branch by increasing normalized frequency:

```bash
CANTILEVER_ACCELERATION_RMS_G=0.2 \
julia --project=examples/damped_intrinsic_beam \
  examples/damped_intrinsic_beam/run_base_excited_cantilever_sweep.jl
```

The sweep defaults to two elements of degree four. This is the least expensive
tested discretization that agrees with four degree-three elements. An isolated
start near the nonlinear fold can converge to the wrong stable branch; use the
sweep for frequency-response comparisons.

On the current Julia 1.11/SciML combination, this Trixi branch has an unrelated
manual-precompile assertion. Until that upstream compatibility issue is fixed,
add `--compiled-modules=no` to the commands above. The example project targets
Julia 1.10, but this benchmark has not yet been rerun there.

Important controls include:

- `CANTILEVER_SWEEP_NORMALIZED_FREQUENCIES`;
- `CANTILEVER_SWEEP_RELTOL` and `CANTILEVER_SWEEP_ABSTOL`;
- `CANTILEVER_SWEEP_PERIODICITY_TOLERANCE`;
- `CANTILEVER_POLYDEG` and `CANTILEVER_REFINEMENT_LEVEL`.

## Experimental marker extraction

The article does not provide a raw frequency-response table. The archived
experimental CSV was extracted from the vector circle markers in Figures 4
and 7, not manually read from raster pixels:

```bash
python3 examples/damped_intrinsic_beam/extract_farokhi_2022_frequency_response.py \
  /path/to/Farokhi2022.pdf \
  /tmp/farokhi_2022_frequency_response.csv
```

The PDF is not redistributed. Marker uncertainty still includes plot-axis
rounding, line thickness, and uncertainty in the original image-processing
experiment.

## Current strict assessment

The archived checks support the following conclusions:

- Two degree-three elements are not adequate: at \(f/f_1=1.0082\), refinement
  raises the transverse amplitude by 35.7% (the coarse value is 26.3% below
  the refined result).
- Two degree-four elements and four degree-three elements agree to 0.17% in
  transverse amplitude and 0.48% in longitudinal amplitude.
- Reducing the time tolerances and maximum step by a factor of four changes
  the reconstructed tip history by at most \(1.54\times10^{-7}\).
- At 0.2g and \(f/f_1=1.0206\), the converged upper branch gives
  \(w_\mathrm{tip}=0.7690\) and \(u_\mathrm{tip,min}=-0.4907\), versus
  digitized experimental values \(0.7864\) and \(-0.5771\). Transverse
  agreement is strong (2.2%); longitudinal shortening is underpredicted by
  about 15%.
- At 0.5g and the same normalized frequency, the converged-amplitude estimate
  gives \(w_\mathrm{tip}=0.8298\) and \(u_\mathrm{tip,min}=-0.7010\), versus
  interpolated experimental values \(0.8178\) and \(-0.6236\). Transverse
  agreement remains strong (1.5%); longitudinal shortening is overpredicted
  by about 12%.

These are good selected-point results, but not yet a complete publication
validation. Before using the experiment in a paper, run and archive the full
upper and lower stable branches at both acceleration levels, add a cycle work
ledger, and confirm at least one more point with the degree-three/four-element
discretization.

Archived files:

- `reference/farokhi_2022_experimental_frequency_response.csv`;
- `reference/base_excited_cantilever_validation.csv`;
- `reference/base_excited_cantilever_numerical_checks.csv`.
