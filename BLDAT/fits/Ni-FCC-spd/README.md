How the `FCC-spd` Ni Bethe-lattice set in `BLDAT/BL028.dat` was made: a minimal spd basis for
the CRENBS ECP, a CRYSTAL23 LSDA calculation of bulk fcc Ni with that basis, and a fit of the
Bethe-lattice parameters to its DOS, in the same way as the Fe `BCC-spd` and `FCC-spd` sets
(`../Fe-BCC-spd/README.md`, `../Fe-FCC-spd/README.md`). Use the set with `BLPAR1/2 = FCC-spd` and the
basis and ECP in `Ni-crenbs-spd.gbs`.

## Basis (`Ni-crenbs-spd.gbs`)

CRENBS s and d shells and ECP (Ar core, 10 valence electrons), the basis of `FCC-sd`
(`../Ni-FCC-sd/Ni-crenbs-sd.gbs`), plus the STO-3G Ni 4p shell (Basis Set Exchange) with all three
exponents multiplied by 0.60.

The scale factor minimizes the CRYSTAL23 LSDA (SVWN) total energy of ferromagnetic fcc Ni,
a = 3.524 A (`crystal/all.txt`; the first line is the sd basis):

| 4p scale | E (Ha) | Moment (muB) |
|---|---|---|
| none (sd) | -40.081566 | 0.699 |
| 0.4 | -40.213140 | 0.569 |
| 0.5 | -40.227679 | 0.547 |
| 0.55 | -40.231929 | 0.538 |
| 0.6 | -40.233032 | 0.533 |
| 0.65 | -40.231298 | 0.532 |
| 0.7 | -40.228280 | 0.535 |
| 0.8 | -40.217381 | 0.548 |
| 1.0 (STO-3G) | -40.195890 | 0.586 |
| 1.2 | -40.184247 | 0.621 |
| 1.5 | -40.168504 | 0.651 |
| 2.0 | -40.137532 | 0.669 |
| 3.0 | -40.097708 | 0.686 |

A parabola through the five points from 0.5 to 0.7 has its minimum at 0.601, so the computed 0.60
point is used (`crystal/s0.6/out`). The fcc lattice is the ground state of Ni, so the scale is
optimized on fcc (for Fe it was optimized on bcc).

`crystal/grid.txt` scales the two tighter and the most diffuse exponent separately (0.45, 0.6, 0.75
each). The lowest point, 0.6 / 0.45, is only 0.004 mHa below the uniform 0.6 scale, so a single scale
factor was kept.

## CRYSTAL23 (`crystal/s0.6/`)

The settings of `../Ni-FCC-sd` (`SHRINK 24 24`, `SMEAR 0.01`, `SPINLOCK 2 5`). The run takes 19 SCF
cycles and gives E = -40.233032 Ha, a moment of 0.533 muB (`crystal/moment`) and E_F = 0.20298 Ha.
`d3w` writes the s, p, d and total DOS of both spins from E_F - 0.58 to E_F + 0.57 Ha
(`DOSS.win.DAT`). E_F lies on a sharp peak at the top of the minority d band: DOS(E_F) is 2.23
states/eV for spin down and 0.11 for spin up.

## Fit (`fit/`)

The parameters are fitted on a 12-direction (fcc, <110>) Bethe lattice, orthogonal, without
overlaps, with the targets, weights and bounds of the Fe spd sets: the s, d and total DOS of both
spins from -10 to +4 eV around E_F and the CRYSTAL23 moment, 0.533 muB, 10 electrons; the p DOS is
not fitted. Hoppings are the same for both spins, the spin-up d on-site energy is fixed at 1.3812 eV
as in `FCC-sd` and the Fe sets, and ETA = 0.1 eV. `fitspd.py` uses the model of
`../Fe-BCC-spd/fit/bethe9.py`.

| Start | From | Cost | Moment (muB) | At a bound |
|---|---|---|---|---|
| A | Ni `FCC-sd` + the p parameters of Fe `FCC-spd` | 42.417 | 0.536 | `ddp`, `ddd` |
| B | Ni `Papacon-old` (BL028.dat, orthogonal spd), spin-averaged hoppings | 41.981 | 0.558 | `ddp` |
| C | Fe `FCC-spd` with spin splittings of 0.6 eV (d) and 0 (s, p) | 41.484 | 0.550 | `dds`, `ddp`, `ddd` |
| D | Ni `FCC-sd` + the p parameters of B | 42.056 | 0.518 | `ddd` |

(The `active bounds` line of the logs, from scipy, misses these; they were found by comparing the
parameters with the bounds.)

`FCC-spd` is `fit/best_B.json` converted to Rydberg. The four fits are almost equally good. B
reproduces the minority DOS at E_F best and has only one parameter at a bound, `ddp` = +0.05 eV; C
has a slightly lower cost but all three d-d hoppings at their bounds. Comparison (`cmpspd.py`; DOS
per spin in states/eV, Bethe / CRYSTAL23, relative L2 misfit from -10 to +4 eV):

| | DOS(E_F) up | down | Misfit up / down |
|---|---|---|---|
| s | 0.054 / 0.039 | 0.033 / 0.023 | 0.65 / 0.63 |
| p (not fitted) | 0.005 / 0.007 | 0.012 / 0.001 | 0.86 / 0.76 |
| d | 0.105 / 0.065 | 1.390 / 2.209 | 0.27 / 0.23 |
| total | 0.164 / 0.111 | 1.436 / 2.234 | 0.26 / 0.22 |

`dos_fit.png` compares B and C with CRYSTAL23. The width and position of the d bands of both spins
and the bottom of the s band are reproduced; the sharp CRYSTAL23 peaks, including the minority peak
at E_F, are not, so the minority DOS at E_F is about 35 % too low.

Fermi level (`efcheck.py`): the occupations below E_F are N_d = 4.63 / 4.64 (Bethe / CRYSTAL23) for
spin up and 4.03 / 4.11 for spin down. The Bethe lattice puts about 0.25 electrons per spin in p
and 0.40 in s, where the CRYSTAL23 Mulliken projection gives about 0 in p (its p projection is
partly negative) and 0.65 in s; s+p is 0.61 / 0.63 (up) and 0.66 / 0.70 (down). A rigid shift of
the Bethe DOS does not lower the d and total misfit (best shift 0.00 eV), so E_F of the fit and of
CRYSTAL23 are consistent. With ETA = 0.1 eV the fit has E_F = 2.900 eV on the scale of the on-site
energies.

## Check in ANT

ANT with `FCC-spd` on both electrodes of the Ni4-Ni-Ni-Ni4 geometry (now `examples/Ni-nanocontacts/Ni4-Ni-Ni-Ni4`) with
the basis of `Ni-crenbs-spd.gbs` (`LEADDOS`, links from `bin/`, 8 November 2026) finds 12 <110>
directions and a Bethe-lattice E_F of 2.818 eV (ETA -> 0, so slightly below the fit value). Its
`Lead1DOS.dat` agrees with `bethe9.py` at the same E_F and ETA = 0.001 eV to 2e-4 states/eV;
DOS(E_F) is 0.130 (up) and 1.565 (down) states/eV.

## Reproducing

- `crystal/scan.sh` runs the sd basis and the scale scan with `run.sh`, which writes
  `s<scale>/INPUT` (or `sd/INPUT`) from `base.d12` with `mk.py` and runs CRYSTAL23; `grid.sh` runs the
  two-scale grid with `mk2.py`. They call
  `/home/palacios/CRYSTAL23OMP/bin/Linux-ifx_i64_omp/v1.0.1/crystalOMP` with `NT` (default 4)
  threads; change the path. The sd run recorded in `../Ni-FCC-sd/crystal/` is `run.sh sd`.
- `crystal/s0.6/d3w` is the CRYSTAL23 properties input that writes `DOSS.DAT`; `DOSS.win.DAT` is a
  copy of it.
- `cd fit; python3 fitspd.py B` repeats the fit and writes `best_B.json` (needs numpy and scipy;
  about 6 s per call, 520 calls). Starts A and D read `../../Ni-FCC-sd/fit/best_sd.json`.
  `python3 cmpspd.py B` prints the comparison above and `python3 efcheck.py B` the Fermi-level
  check. `python3 plotspd.py B C` writes `dos_fit.npz`, and `plotspd.py --plot`, run with a Python
  that has matplotlib, draws `dos_fit.png` from it.

CRYSTAL23 1.0.1 (OpenMP, ifx) on 4 threads; fit with numpy 1.19 and scipy 1.5; 8 November 2026.
