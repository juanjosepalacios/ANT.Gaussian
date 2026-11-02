> Historical method: 1D lead Fermi searches now use mandatory band filling.
> See [BAND_FERMI.md](BAND_FERMI.md).

# Per-electrode warm-start Fermi brackets

Included in master for 2.8.4 (implementation commit 55f3fe2).

## Behavior

Each lead number has an independent validity flag and previous Fermi energy.
The energy is saved on successful convergence, before AdjustFermi shifts that
lead's matrices. Lead 2 never borrows lead 1's solution, even if the example is
symmetric. Histories survive SCF electrode cleanup/reinitialization.

The first search for each electrode uses that electrode's full spectral
interval. Later searches start within +/-0.25 eV of its own previous energy,
clamped to its own current bounds. The half-width is a search heuristic, not
a convergence tolerance; it is at least twice the effective energy tolerance.
If the charge does not change sign across the interval, the bracket doubles
until it brackets a root or covers the full spectral interval. Failure to
bracket even on the full interval, or nonfinite charge values, is an explicit
error rather than an accepted result.

No charge values persist across SCF steps or between electrodes. Fresh endpoint
charges used to establish the current bracket are passed into BrentRoot to avoid
evaluating those same endpoints twice. Existing callers without optional endpoint
values keep their behavior. Counts include the actual bracketing and root-search
charge evaluations, and timers include both parts.

## Validation

python3 tests/test_warm_fermi.py
FC=ifx python3 tests/test_warm_fermi.py

NVIDIA and Intel pass with bounds checking. Tests extract the real FindFermi,
CountedLeadCharge, BracketRoot, BrentRoot and history declarations. Distinct lead
models have initial roots -7/+3 eV and electron targets 7/11, then change to
-6.97/+3.08 with reversed update order. They verify that cold lead 2 does not use
lead 1 history, each later search uses its own stored energy, a large lead-1 shift
expands to the full interval while preserving lead 2 history, and public shifted
EFermi fields are not mistaken for unshifted guesses. Additional cases cover
changed bounds, reuse of fixed endpoints, absent brackets, invalid inputs and
nonfinite charges. Existing Brent tests pass on both compilers, and the NVIDIA
parallel 1D-charge regression passes.

## Au10 timing comparison, 31 October 2026

Runs/builds: /home/palacios/ANT-runs/onedlead-brent-warm/.
Full-bracket Brent reference first, then warm-start Brent. Identical saved density,
NPC=3, four threads, no LEADDOS. The reference has the same code and diagnostics
but disables the warm-start branch. Both runs terminate normally in two SCF cycles
without charge/Dyson warnings and have 1001 finite transmission samples.

| Measurement | Full bracket | Per-lead warm bracket |
| --- | ---: | ---: |
| Initial search calls, each lead | 15 | 15 |
| Final search calls, each lead | 13 | 5 |
| Total Fermi-search charge calls | 56 | 40 |
| Total Fermi-search time (s) | 125.998082 | 110.817154 |
| Full run time (s) | 245.821200 | 230.689109 |
| Initial Fermi energies (eV) | -6.80597 | -6.80597 |
| Final Fermi energies (eV) | -6.80675 | -6.80672 |
| SCF energy (hartree) | -1356.43701308 | -1356.43701308 |

Final searches require 62% fewer calls, but wall savings are smaller because
quadrature cost varies with trial energy. Overall search time decreased by 12.0%
and full-run time by 6.2% in this single sequential comparison. No statistical
performance claim is made. Actual unequal-electrode physics was not benchmarked;
separate histories/electron targets are tested with deliberately unequal models.
Transmission at EF is 1.9988 for both runs; maximum absolute spectrum difference
is 0.0057 and RMS difference is 0.0003553 on identical energy grids.

## Accuracy-policy clarification

Input FERMIACC and CHARGEACC were both 1e-4. ANT.F90 adapts the effective FermiAcc
and ChargeAcc during SCF: in the first two lead-enabled cycles it multiplies the
input targets FermiA/ChargeA by ten unless FULLACC is requested. These restart
comparisons therefore used effective 1e-3 values during those stages. Earlier
reports describing the effective tolerance as 1e-4 stated the input rather than
the scaled runtime value. The warm-start patch preserves this existing policy.
Warm-search final brackets of 9.7338e-4 eV satisfy that effective 1e-3 tolerance.
Use FULLACC to request the input accuracy throughout SCF; that is a separate
accuracy/performance comparison, not enabled silently by this patch.
