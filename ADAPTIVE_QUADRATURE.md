> Lead Fermi levels now use mandatory band filling. The lead-charge contour
> experiment below is historical; ADAPTQUAD remains available for device
> density contours. See [BAND_FERMI.md](BAND_FERMI.md).

# Adaptive contour quadrature (experimental, opt-in)

Add ADAPTQUAD on its own line in the .ant input. The default is off, leaving
the legacy integration rules unchanged. The current implementation covers
1D lead charge contours and collinear device density-matrix contours.
Spinor/SOC density contours, imaginary-axis integration, nonequilibrium
real-axis corrections and device spectral-weight checks retain their
existing integration rules.

The embedded Gauss 7 / Kronrod 15 rule starts from eight panels clustered
near the contour endpoints. It estimates local error and bisects the panel
with the largest normalized error. Existing panel integrals are retained;
only the two children require new Green-function evaluations. Each panel
batches 15 independent nodes across OpenMP threads, then accumulates in a
fixed order. Context is passed explicitly to avoid thread-dependent captured
callback state. Nested teams are suppressed.

Lead charge uses effective CHARGEACC (or its existing per-call override).
Device integration monitors every real density-matrix element with PACC,
plus an explicit Tr[P*S] component with half of min(PACC, CHARGEACC) per
spin. Cancellation in total charge alone cannot declare convergence.
The global error criterion sums per-panel maximum normalized component
error estimates. Estimates use the embedded difference, variation rescaling
and a floating-point roundoff floor. They are numerical estimators, not
rigorous mathematical bounds.

The limit is 128 panels (at most 3720 nodes with the current eight-panel
start). Nonfinite integrands, callback failures, resource exhaustion or
unresolved precision produce explicit failures, rather than accepting an
unconverged density. Panel storage scales with the matrix size; for a
660-orbital collinear device, retained matrix panels require about 446 MB,
plus point batches, Green-function workspaces and normal calculation storage.

Validation extracts production routines, not rewritten algorithms:
- test_adaptive_gk.py: narrow resonance near an endpoint, cancelling vector
  components, 1-8 thread reproducibility, nonfinite data, exhaustion and
  invalid inputs.
- test_adaptive_density.py: a rotated two-orbital Green function with an
  analytic Lorentzian integral, off-diagonal density entries, both spins,
  and exact thread reproducibility.
- test_onedlead_charge.py: real lead charge/Dyson routines with and without
  ADAPTQUAD on isolated, hopping, nonorthogonal and complex model leads.

NVIDIA and Intel bounds-checked tests pass. The opt-in should be benchmarked
against a converged reference before changing the default. Identical input
PACC values do not make the old charge-based estimator equivalent to the
new elementwise matrix-error check.

## First Au10 experiment, 1 November 2026

Same initial density, ETA=0.01, NATPC1/2=1, NPC=3, eight threads and input
tolerances as Au10-onset-coarse-eta001-20261101. ADAPTQUAD completed normally
in six SCF cycles, 298.59 s, with no integration warnings. The prior legacy
run took 66.65 s and six cycles but reported five failed density integrations.
The adaptive run used 9180 lead-charge nodes and 29460 device-density nodes
across all searches. This experiment ran alongside Au30; the wall times are
not a controlled performance benchmark. The stronger matrix criterion also
makes nominal PACC values inequivalent between methods.

Both output spectra/densities are finite. T(EF) is 1.9983 (adaptive) versus
1.9986 (legacy). Maximum density-coefficient difference is 2.0041e-4;
maximum transmission difference is 3.92792 and RMS difference 0.13481 on the
same 1001-point grid. Sharp off-Fermi structure is therefore not established
as converged by this comparison. A tighter reference is needed before
interpreting those differences physically or adopting this method by default.

Artifacts: /home/palacios/ANT-runs/Au10-adaptive-gk-20261101/.
Build: /home/palacios/ANT-runs/adaptive-build-20261101/.
The subsequent inversion-status guards were regression-tested and rebuilt;
the experiment above reported no inversion warnings.
