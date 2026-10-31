# Dyson self-energy mixing (2.8.4)

The default L1DALPHA is now 0.75 (formerly 0.25). Set it in the .ant input:

```text
L1DALPHA = 0.75
SELFACC = 1.0d-6
```

Sigma_new = alpha F(Sigma_old) + (1-alpha) Sigma_old.
Allowed range: 0 < alpha <= 1. This differs from density mixing ALPHA.
SELFACC defaults to 1e-6 and is not adapted during SCF.

Convergence uses mean_ij abs(F(Sigma)-Sigma), before mixing, reset each
iteration. This corrects the mixing-dependent residual and false convergence
at alpha=1. Invalid inputs, nonfinite residuals, failed LAPACK inversions and
exhausted iterations produce explicit failures.

## Au10 benchmark, 31 October 2026

Identical saved-density restarts, NPC=3, eight threads, SELFACC=1e-6,
input FERMIACC=CHARGEACC=1e-4, no LEADDOS. Existing SCF adaptation made
the effective Fermi/charge tolerances 1e-3 in these two-cycle restarts.

| Mixing | Full runtime (s) | Speedup |
| --- | ---: | ---: |
| 0.25 | 165.1055 | 1.00 |
| 0.50 | 81.8291 | 2.02 |
| 0.70 | 58.1645 | 2.84 |
| 0.80 | 52.3773 | 3.15 |

All runs terminated normally in two SCF cycles, without warnings.
Final Fermi energies were -6.80696 eV for both leads; T(EF)=1.9988.
At 0.8 versus 0.25, total energy differed by 9e-8 hartree; maximum
transmission difference was 0.0076 and RMS difference 0.000326 over
the identical 1001-point grid.

A frozen-matrix sweep tested 2080 energy/direction cases with actual Au10
lead matrices. All cases converged for every tested mixing from 0.5 to 0.9.
Mean iterations were lowest at 0.7 (55.792); 0.75 gave 56.518,
versus 65.587 at 0.5. Trials used SELFACC=1e-6 with a 1e-9 reference.

The selected default 0.75 lies between the frozen-grid optimum near 0.7
and the fastest full run tested, 0.8. A full run at exactly 0.75 was not
performed. These are case-specific measurements, not a universal optimum.

## Regression tests and artifacts

Run python3 tests/test_dyson_mixing.py, optionally with FC=ifx.
Tests extract the production solver and cover analytic surface self-energy,
mixing-independent residuals, alpha=1 false convergence, invalid inputs,
iteration exhaustion and singular inversions.

Benchmark builds, logs, inputs and spectra are outside the repository:
 /home/palacios/ANT-runs/dyson-mixing-au10/
Its matrix-export helper is benchmark-only, not production code.
