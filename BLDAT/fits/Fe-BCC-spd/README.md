How the `BCC-spd` Fe Bethe-lattice set in `BLDAT/BL026.dat` was made: a minimal spd basis
for the CRENBS ECP, a CRYSTAL23 LSDA calculation of bulk bcc Fe with that basis, and a fit of
the Bethe-lattice parameters to its DOS. Use the set with `BLPAR1/2 = BCC-spd` and the basis in
`Fe-crenbs-spd.gbs`.

## Basis (`Fe-crenbs-spd.gbs`)

CRENBS s and d shells and ECP (Ar core, 8 valence electrons), the same as the sd basis of the
Fe-wire example, plus the STO-3G Fe 4p shell with all three exponents multiplied by 0.586.

The scale factor minimizes the CRYSTAL23 LSDA (SVWN) total energy of ferromagnetic bcc Fe,
a = 2.87 A (`crystal/all.txt`):

| 4p scale | E (Ha) | Moment (muB) |
|---|---|---|
| 0.4 | -22.267313 | 2.885 |
| 0.5 | -22.281741 | 2.831 |
| 0.55 | -22.284895 | 2.819 |
| 0.6 | -22.285115 | 2.815 |
| 0.65 | -22.283472 | 2.819 |
| 0.7 | -22.279027 | 2.830 |
| 0.8 | -22.265605 | 2.859 |
| 1.0 (STO-3G) | -22.234086 | 2.905 |
| 1.2 | -22.214907 | 2.919 |
| 1.5 | -22.192649 | 2.911 |
| 2.0 | -22.154022 | 2.884 |
| 3.0 | -22.111020 | 2.884 |

A parabola through the five points from 0.5 to 0.7 has its minimum at 0.586. At 0.586 the energy
is -22.285306 Ha and the moment 2.81 muB (`crystal/s0.586/out`).

`crystal/grid.txt` scales the two tighter and the most diffuse exponent separately (0.45, 0.586,
0.75 each). The lowest point, 0.586 / 0.45, is only 0.27 mHa below the uniform 0.586 scale, so a
single scale factor was kept.

## Fit (`fit/`)

The parameters are fitted on an 8-direction (bcc) Bethe lattice, orthogonal, without overlaps,
to the CRYSTAL23 DOS at the 0.586 scale (`crystal/s0.586/DOSS.win.DAT`). The fit matches the s,
d and total DOS of both spins from -10 to +4 eV around E_F, and the moment of 2.81 muB; the p DOS
is not fitted. Hoppings are the same for both spins, the d on-site energy of spin up is fixed to
that of `BCC-sd` (the energy reference), and ETA = 0.1 eV.

`bethe9.py` is a Python version of ANT's Bethe-lattice Green's function (`BetheLattice.F90`) for
s+p+d sets. `fitspd.py` runs a bounded least-squares fit from one of four starting points:

| Start | Cost | Moment (muB) | Note |
|---|---|---|---|
| A | 47.449 | 2.916 | from `BCC-sd` plus guessed p parameters |
| B | 47.445 | 2.917 | `dds` at its bound |
| C | 47.445 | 2.916 | `dds` at its bound |
| D | 46.983 | 2.804 | no parameter at a bound; used for `BCC-spd` |

`BCC-spd` is `fit/best_D.json` converted to Rydberg. `cmpspd.py` compares it with CRYSTAL23
(DOS per spin in states/eV, relative L2 misfit from -10 to +4 eV):

| | DOS(E_F) Bethe / CRYSTAL, up | down | Misfit up / down |
|---|---|---|---|
| s | 0.055 / 0.046 | 0.046 / 0.024 | 0.63 / 0.51 |
| p (not fitted) | 0.033 / 0.032 | 0.042 / 0.100 | 0.82 / 1.12 |
| d | 0.043 / 0.022 | 0.756 / 0.759 | 0.21 / 0.36 |
| total | 0.131 / 0.100 | 0.844 / 0.884 | 0.21 / 0.30 |

## Reproducing

- `crystal/run.sh <scale>` writes `s<scale>/INPUT` from `base.d12` with `mk.py` and runs
  CRYSTAL23; `run2.sh <s1> <s2>` does the same for two scales with `mk2.py`. Both call
  `/home/palacios/CRYSTAL23OMP/bin/Linux-ifx_i64_omp/v1.0.1/crystalOMP`; change the path.
- `crystal/s0.586/d3w` is the CRYSTAL23 properties input that writes `DOSS.win.DAT`.
- `crystal/dos.py` compares the sd and spd DOS; it reads the sd run from a directory that is not
  included.
- `python3 fit/fitspd.py D` repeats the fit and writes `best_D.json` (needs numpy and scipy;
  about 6 s per call, 700 calls). It reads the DOS from an absolute path on lian; change it to
  `../crystal/s0.586/DOSS.win.DAT`.

CRYSTAL23 1.0.1 (OpenMP, ifx), 6 November 2026.
