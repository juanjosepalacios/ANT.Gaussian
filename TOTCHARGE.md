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
Device TotCharge still falls back to serial for 1D electrodes; the separate
OneDLead bulk-charge integration is now parallel (see below).
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


## OneDLead bulk charge: OpenMP fix

`OneDLead::TotCharge` now distributes contour points across an OpenMP team.
Each point calls `ChargePoint1D` with its own Green matrix; the new
`CompGreensFuncWorkspace` uses call-local self-energy matrices. The existing
`CompGreensFunc` wrapper preserves the module `GX` output for serial DOS calls.
Removed the explicit blanket SAVE from `CompSelfEnergy1D` and marked the
point-evaluation call chain recursive so local scratch is automatic.
`numeric::CInv` must also be built with automatic locals, as in the existing
NVIDIA `-Mrecursive/-Mreentrant/-Mnosave` and Intel recursive configurations.

Weighted contributions are saved per point and summed in original order.
Thread counts therefore produce identical results in the regression tests.
Compared with the original expression, storing a contribution can introduce
last-bit rounding differences (for example, by separating a fused multiply-add).
Contour, spin factor, 0.1 stopping tolerance and maximum order remain unchanged.
The previous estimate is initialized explicitly because Fortran need not
short-circuit the first convergence test. Calls inside an existing OpenMP
region use an inactive inner team. Lead initialization/selection must remain
serial: WhichLead, HX, SX and dimensions are shared read-only during integration.
This does not enable concurrent initialization of different electrodes.

Run `python3 tests/test_onedlead_charge.py` or
`FC=ifx python3 tests/test_onedlead_charge.py`. The test extracts actual charge,
Green-function, Dyson, hopping/overlap accessor, inversion and quadrature routines,
and links real LAPACK/BLAS. It checks one/two spin channels, nonzero hopping,
nonorthogonal overlaps, complex Hamiltonians, both contour branches, 127-point
concurrent evaluation, and threads 1 through 8 with three repeats. It verifies
agreement within 1e-12 with the archived original serial integration, exact
agreement across thread counts, analytic isolated-orbital Green functions,
unchanged legacy module scratch, actual team size, and nested-team fallback.
Both NVIDIA and Intel pass with bounds checking. The existing device-charge
regression also passes. ANT_TEST_LIBS can override the test link flags; libraries
must use the same 8-byte integer ABI as the compiler flags.

The full NVIDIA executable builds in `/home/palacios/ANT-runs/onedlead-openmp`.
The Au-Au-Au-Au 1D example was compared using the prior build at one thread and
the candidate at one/four threads. Both leads initialize with the same printed
Fermi energy (-1.58691 eV), identical quadrature refinement counts and byte-identical
Lead1DOS/Lead2DOS files. The observed lower-bound charge difference from the old
build is about 1.4e-16. All three calculations subsequently fail in the existing
device Fermi search (root not bracketed), so final transmission and full-run
speedup are not validated. H6-chain likewise develops nonfinite values after
lead initialization with the original build and was stopped. These downstream
issues were not changed as part of the parallelization fix.

Logs, baseline and candidate runs, and the separate executable are preserved
under `/home/palacios/ANT-runs/onedlead-openmp/`. Installed binaries are unchanged.


## Full-supercell charge region and DOS output follow-up

For two 1D electrodes, InitDevice now excludes NPC*NEmbed atoms at each end
from the reduced-device population region, matching the supercell on which
self-energies act. NCDEl and NCDAO1/NCDAO2 are computed from the same region.
The input NEmbed still specifies atoms per primitive cell; it is not changed.
A meaningful reduced-device region must remain between the excluded buffers.

Au10/lanl2dz at four threads and NPC=3 previously aborted with integrated
spectral weight 352.40420 versus expected 352. Excluding complete supercells
changes the reduced orbital range from 23..198 to 67..154: the check reports
176.00000 versus expected 176 and passes. Both initial lead Fermi energies
remain -7.37305 eV. The separate no-LEADDOS test then terminated normally.
Build and calculation: /home/palacios/ANT-runs/onedlead-buffer-check/.

The earlier NPC=2 unlimited baseline/parallel runs were stopped at the user's
request. Both failed 2047-point quadrature at -12.5 eV with a successive-estimate
difference about 2954.4. The baseline reported -8.98438 eV after another failed
integration; this is not a validated Fermi level. NPC=3 avoided those failures
in the tested initial lead searches. These runs do not establish speedup.

DOS output previously mixed the unshifted temporary HX with already shifted
lead hopping matrices and an additional -EFermi energy offset. When LEADDOS
is requested, shift temporary HX by -EFermi*SX for each spin and evaluate at
energy+i*gamma, consistently referencing all matrices to the lead Fermi level.
The temporary matrices are deallocated after lead initialization; this output
correction does not change the transport Hamiltonian.

The matching DOS snapshots from Au10 SCF cycle 26 are finite and identical
between leads at saved precision, with Fermi energies -6.66504 eV. They cover
-10..10 eV at 0.01 eV spacing. Raw nonmagnetic files contain one spin channel;
total DOS includes a factor of two. GAMMA=0.03 plus ETA=0.01 gives an effective
imaginary shift of 0.04 eV. The snapshot does not establish final SCF convergence
of the DOS-enabled run. Separate build/run: /home/palacios/ANT-runs/onedlead-dos-check/.

These tests do not resolve ignored LAPACK status codes, the Dyson error
accumulator, or propagation of failed lead quadrature to the Fermi search.
Those pre-existing issues remain separate work.
