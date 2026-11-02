# Device Fermi-level search

The device uses safeguarded Brent iteration for collinear, spin-locked and
SOC charge searches, both on the semicircular contour and imaginary axis.
Periodic lead Fermi levels use band filling separately; see
[BAND_FERMI.md](BAND_FERMI.md).

`device::SolveDeviceFermi` expands a bracket around the current device
shift, within the existing device bounds, and calls `numeric::BrentRoot`.
The initial half-width is max(0.25 eV, 2*FERMIACC); this is a bracketing
heuristic, not an accuracy target. The root uses FERMIACC and at most 5*MAX
iterations. The accepted point is evaluated again so the returned density
and populations correspond to the reported Fermi energy, even when Brent
returns a previously evaluated endpoint. Search variables are device
shifts; printed Fermi energies have the opposite sign.

The final charge residual must satisfy the existing device tolerance
CHARGEACC*(NCDEl+QEXCESS). This preserves the previous device convention;
1D lead band filling instead uses CHARGEACC as an absolute electron count
per chosen repeated cell. Missing brackets, nonfinite values, iteration
exhaustion and unacceptable final charge residuals stop explicitly.
Brent controls root-search error; it cannot fix unconverged density
quadrature or a failing device spectral-weight check.

FIXEFERMI bypasses the contour root search and now sets its success flag
explicitly. The fixed point is evaluated to obtain its density. No
hard-coded 0.1 eV fallback tolerance is used.

The extracted production wrappers and helper pass NVIDIA and Intel model
tests for nonlinear charge roots, both spin branches, SOC, imaginary-axis
dispatch, final density state, fixed Fermi, missing brackets, nonfinite
charge and residual failure.

## Device spectral-weight windows

The preliminary spectral-weight check first integrates from -1000 to
+1000 eV on the existing complex contour. A failed quadrature, nonfinite
result or mismatch of at least 0.1 electrons from the expected spectral
weight triggers one retry from -10000 to +10000 eV. If the second check
fails, stop explicitly. The 0.1-electron spectral-weight acceptance and
successive-quadrature criteria are unchanged. These temporary check windows
do not alter the density-contour bounds or device root bracket. The old
incremental energy-bound expansion is removed.

## Au10 production comparison (2026-11-01)

Both eight-thread runs used the same saved starting density, ETA=0.01,
NATPC1/2=1, NPC=3, FULLACC, CHARGEACC=FERMIACC=1e-4 and PACC=1e-7.
Both used the mandatory band-based lead solver; only the device root
implementation differed. Both terminated normally after six SCF cycles.

| Device root | Wall time (s) | Final device EF (eV) | Integration warnings | T(EF) |
|---|---:|---:|---:|---:|
| Original Muller/secant/bisection | 68.23 | -6.63166 | 1 | 1.9983 |
| Safeguarded Brent | 76.11 | -6.630973 | 0 | 1.9983 |

The final Brent charge residual was -7.35e-11 electrons, against the
preserved device tolerance of 0.0076 electrons. Its final shift bracket was
3.91e-6 eV, below FERMIACC. Typical searches used five Brent iterations,
eight charge evaluations including bracketing and the final density
refresh. These concurrent runs do not establish a controlled speed benefit.

The energies differ by 2.80e-4 Ha; finite transmission spectra have a
maximum difference of 0.5345 and RMS difference 0.02447 on the relative
energy grid. Matching T(EF) does not prove the full spectra converged.
The final density refresh changes the old stateful callback behavior,
and the old run also had an integration warning.

Au30 (NATPC1/2=3, NPC=4, ETA=0.01) found both lead Fermi levels promptly
but remained in the preliminary device spectral-weight check for 301 s.
The last completed 31-point quadrature estimate was 263.72075, target 264,
with device bounds [-128.84, 14217.65] eV. The diagnostic run was stopped
before any device Fermi search. Thus the Au30 device root is not yet
validated: this spectral-weight integration is a separate prerequisite
bottleneck, unaffected by replacing the root finder with Brent.

## Historical lead-search implementation

The following notes describe the superseded lead method, not current lead
behavior.



Included in master for 2.8.4 (implementation commit 55f3fe2).

## Implementation

numeric::BrentRoot maintains a sign-changing bracket, attempts inverse
quadratic/secant interpolation, and falls back to bisection when a proposal
fails the safeguards. It returns status, iterations, function calls, residual,
and bracket width. Invalid arguments, missing brackets, nonfinite callback
values, and iteration exhaustion are explicit failures.

OneDLead::FindFermi uses BrentRoot and reports lead-search call counts and wall
time. The full original bracket is retained, but the current lead code uses
FermiAcc as absolute energy accuracy in eV and ChargeAcc as absolute charge
accuracy in electrons per primitive cell. ChargeAcc governs successive-order
integration convergence and lead spectral-bound checks. The special zero-charge
shortcut is removed, so successive estimates must agree. Nonfinite/nonpositive
ChargeAcc is rejected. Per-lead warm-start bracketing is now implemented; see WARM_START.md.
There is no persistent charge-value cache.
OneDLead::TotCharge now stops
with status 1 on quadrature exhaustion instead of returning the failed estimate.
Other users of numeric::Bisec are unchanged.

## Tests

python3 tests/test_brent.py
FC=ifx python3 tests/test_brent.py

Both NVIDIA and Intel pass with bounds checking. Cases include increasing and
decreasing functions, reversed intervals, flat and steep roots, exact endpoint
roots, invalid tolerance/bracket, nonfinite evaluations, and iteration limits.
A smooth exponential example takes 10 charge-function calls versus bisection's
37 at 1e-10 tolerance. Production routines are extracted, not reimplemented.
The existing NVIDIA OneDLead parallel-charge regression also passes.

## Original coarse-tolerance Au10 comparison, 31 October 2026

Separate baseline/candidate builds and runs:
/home/palacios/ANT-runs/onedlead-brent/

Same NPC=3 input, identical saved density files, four Gaussian/OpenMP threads,
MKL_NUM_THREADS=1, LEADDOS disabled. Runs executed sequentially, baseline first.
The reference has the same charge-failure guard and instrumentation but calls
original Bisec. Both runs terminate normally in two SCF cycles without charge
quadrature or Dyson convergence warnings. Both spectra have 1001 finite samples
on the same -5..5 eV grid.

| Measurement | Bisection | Brent |
| --- | ---: | ---: |
| Initial lead charge evaluations, each | 16 | 13 |
| Final lead charge evaluations, each | 16 | 11 |
| Total lead search evaluations | 64 | 48 |
| Total lead Fermi search time (s) | 88.513004 | 49.889899 |
| Full calculation time (s) | 185.217014 | 146.663649 |
| Initial Fermi energy (eV), both leads | -6.86035 | -6.79657 |
| Final Fermi energy (eV), both leads | -6.81152 | -6.80294 |
| Transmission at zero relative energy | 1.9989 | 1.9988 |
| SCF energy (hartree) | -1356.43720656 | -1356.43698024 |

Lead searches took 43.6% less time; the full run took 20.8% less time. These are
single runs per method, not statistical performance estimates. Fewer function
calls alone do not fully predict runtime because quadrature orders vary with
trial energy.

## Coarse-tolerance accuracy limitation

The maximum absolute transmission difference over the grid is 1.3063 near
sharp spectral structure. The initial Fermi-energy difference is 0.06378 eV
and the final difference is 0.00858 eV, both within the historical 0.1 eV search
tolerance, but that tolerance is inadequate for proving matching spectra.
Brent's final charge residuals are about 0.00510 (initial leads) and 0.00217
(final leads), with final bracket widths 0.08993 and 0.03882 eV. Charge quadrature
still uses its historical 0.1 absolute successive-estimate threshold.

These results establish that Brent works and is faster in this case, not
scientific equivalence at stringent accuracy. For a stringent accuracy comparison, tighten and
validate the charge-integration and Fermi-energy tolerances together, then
repeat the identical-density comparison and check spectral convergence.


## Configured-tolerance repeat, 31 October 2026

All active 0.1 convergence thresholds in OneDLead are removed. FindFermi uses
FERMIACC; TotCharge and FindEBounds use CHARGEACC. This scope does not change
the independent device integration or its spectral-weight-check tolerances.
The serial charge regression reference now uses the same configured tolerance
rather than retaining its original 0.1 threshold. NVIDIA/Intel Brent and 1D-charge
regressions pass, with bounds checking.

Repeat directory: /home/palacios/ANT-runs/onedlead-brent-tight/.
Both runs use FERMIACC=1e-4 eV, CHARGEACC=1e-4 electrons, NPC=3, four threads,
identical starting saved densities, and the same remaining settings as above.
Bisection runs first, then Brent, without overlap. Both terminate normally in
two SCF cycles without integration/Dyson warnings. No quadrature limit was raised.

| Measurement | Bisection | Brent |
| --- | ---: | ---: |
| Initial Fermi-search evaluations, each lead | 23 | 15 |
| Final Fermi-search evaluations, each lead | 23 | 13 |
| Total Fermi-search evaluations | 92 | 56 |
| Summed Fermi-search time (s) | 271.404204 | 125.795183 |
| Full run time (s) | 391.244920 | 245.774639 |
| Initial Fermi energies (eV) | -6.80599 | -6.80597 |
| Final Fermi energies (eV) | -6.80695 | -6.80675 |
| SCF energy (hartree) | -1356.43701318 | -1356.43701308 |
| Transmission at zero relative energy | 1.9988 | 1.9988 |

Brent used 39.1% fewer Fermi-search evaluations, 53.7% less search time,
and 37.2% less full-run time. These remain single-run measurements, not
statistical estimates. Both spectra have 1001 finite samples on identical grids.
Maximum absolute transmission difference: 0.035 at +4.64 eV; RMS difference:
0.0021686. The previous coarse comparison had a maximum difference of 1.3063.

The initial lead Fermi values agree within 2e-5 eV. The final values differ by
2e-4 eV. A root-bracket tolerance bounds the root of the numerically evaluated
charge function, not integration error or differences between the final SCF
Hamiltonians. CHARGEACC=1e-4 can therefore still limit physical Fermi accuracy.
Brent reports charge residuals 1.56e-8 and -2.49e-9 for its approximate charge
function; these do not bound the underlying integration error. Tighter CHARGEACC
would be the next step if closer agreement at narrow spectral features is needed.

Merged into master for 2.8.4. Measurements above describe comparison builds;
rebuild to update installed executables.


## Accuracy-policy correction and warm-start follow-up

The 1e-4 values quoted for the configured-tolerance comparisons above are the
input targets FermiA/ChargeA. ANT.F90 scales effective FermiAcc/ChargeAcc by ten
in the first two lead-enabled SCF cycles unless FULLACC is set. The two-cycle
restart comparisons therefore used effective 1e-3 tolerances at those stages.
Removing the hard-coded lead thresholds did not disable that existing adaptive
policy. WARM_START.md documents the independent per-electrode histories, new
regressions and controlled full-bracket versus warm-bracket comparison.
