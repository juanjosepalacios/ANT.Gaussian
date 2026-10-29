# TotCharge parallelization (2.8.3)

## Why the previous parallel loop failed

DDOS used device's module-level ispin as its spin-loop counter, while
gplus0 read the same shared variable to select HD and lead self-energies.
Independent quadrature points therefore changed one another's spin state.
Changing the reduction clause did not resolve this race.

DDOS now uses a local spin_index and passes it explicitly to gplus0.
The optional gplus0 argument preserves existing callers that use module
spin. Both routines are recursive so their local scratch is automatic.
The NVIDIA build uses -Mrecursive/-Mreentrant/-Mnosave; Intel uses automatic
locals. Do not build with static local scratch.

## Quadrature implementation

Generate nodes once per order, then distribute independent points with an
OpenMP static schedule on NVIDIA and Intel. Every point writes to its own
contribution slot. After joining the team, sum in original index order.
This avoids atomics, critical regions and reduction-order differences.
Nodes, tolerances, refinement and bounds checks remain unchanged.

Bethe-lattice and ghost lead paths are enabled. CompSelfEnergyBL and its
Dyson solver use call-local scratch and read electrode data in this path.
1D leads retain explicit saved scratch, so TotCharge falls back to serial.
Calls from an existing OpenMP region avoid creating a nested team.

## Regression tests

Run from this worktree:

    python3 tests/test_totcharge.py
    FC=ifx python3 tests/test_totcharge.py

Tests extract actual gplus0, DDOS and TotCharge into a two-orbital analytic
fixture, with an analytic inverse and zero lead self-energy. This isolates
spin selection, scheduling and accumulation, not full electrode physics.
Checks cover one/two spin channels, actual team sizes 1 through 8,
three repeats, unchanged module spin, optional-spin compatibility,
exact charge agreement and serial 1D fallback, with bounds checking.
Both NVIDIA pgf90 and Intel ifx passed. Full Gaussian builds/runs were NVIDIA;
full Intel Gaussian and magnetic-device validation were not performed.

## Real-data benchmark

28-atom Pt SOC saved-density restart, 252 basis functions, eight Gaussian
and OpenMP threads, MKL_NUM_THREADS=1, 4 GB Gaussian memory.
The serial reference is the same implementation with parallel_charge
forced false in a separate build. Only quadrature parallelism differs.
Run order: serial, parallel, parallel, serial.

| Mode | Repeat | TotCharge (s) | Full run (s) |
|---|---:|---:|---:|
| serial | 1 | 3.074434 | 52.265 |
| parallel | 1 | 0.434601 | 50.011 |
| parallel | 2 | 0.418965 | 50.300 |
| serial | 2 | 3.061789 | 52.390 |

Mean TotCharge: 3.0681 s serial versus 0.4268 s parallel (7.19x).
Mean full-run time: 52.327 s versus 50.155 s (4.2% reduction).
These are two runs per variant, not statistical confidence intervals.

All four runs converged in two cycles with energy -770.753121135 hartree.
Charge diagnostics agree at printed precision. All 601 transmission values
agree at output precision. The spectrum also matches the earlier fixed-build
eight-thread Pt result byte for byte.

Logs, binaries and JSON: /home/palacios/ANT-runs/totcharge/.
Serial build: /home/palacios/ANT-runs/totcharge-serial/.
Main checkout, installed executables and v2.8.2 tag remain unchanged.
Release v2.8.3 identifies this implementation and its recorded validation.
