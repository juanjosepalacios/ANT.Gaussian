How the `FCC-spd` Fe Bethe-lattice set in `BLDAT/BL026.dat` was made: a CRYSTAL23 LSDA
calculation of bulk fcc Fe with the spd basis of `BCC-spd`, and a fit of the Bethe-lattice
parameters to its DOS. Use the set with `BLPAR1/2 = FCC-spd` and the basis in
`../Fe-BCC-spd/Fe-crenbs-spd.gbs`.

## CRYSTAL23 (`crystal/`)

Ferromagnetic fcc Fe, a = 3.515 A, SVWN, the same settings as the bcc spd run (`SHRINK 24 24`,
`SMEAR 0.01`, `SPINLOCK 4 5`). The basis is not re-optimized: the 4p scale stays at 0.586, the
value that minimizes the bcc energy (see `../Fe-BCC-spd/README.md`). `INPUT` differs from
`../Fe-BCC-spd/crystal/s0.586/INPUT` only in the space group (225) and lattice constant.

| | fcc (this set) | bcc (`BCC-spd`) |
|---|---|---|
| SCF cycles | 18 | 17 |
| E (Ha) | -22.282440 | -22.285306 |
| Moment (muB) | 2.797 | 2.81 |

fcc is 2.9 mHa (78 meV) above bcc. `d3w` is the properties input that writes the s, p, d and
total DOS of both spins (`DOSS.win.DAT`, E_F = 0.3832 Ha).

## Fit (`fit/`)

The parameters are fitted on a 12-direction (fcc, <110>) Bethe lattice, orthogonal, without
overlaps, with the same targets, weights and bounds as `BCC-spd`: the s, d and total DOS of both
spins from -10 to +4 eV around E_F and the CRYSTAL23 moment, 2.797 muB; the p DOS is not fitted.
Hoppings are the same for both spins, the spin-up d on-site energy is fixed at 1.3812 eV as in
`FCC-sd`, `BCC-sd` and `BCC-spd`, and ETA = 0.1 eV. `fitfcc.py` uses the Bethe-lattice model in
`../Fe-BCC-spd/fit/bethe9.py`.

| Start | From | Cost | Moment (muB) |
|---|---|---|---|
| A | `FCC-sd` + the p parameters of `BCC-spd` | 39.383 | 2.794 |
| B | `BCC-spd` with the d-d hoppings of `FCC-sd` | 52.512 | 2.719 |
| C | the starting point of the `BCC-spd` fit | 41.611 | 2.786 |
| D | `BCC-spd` | 39.412 | 2.797 |

`FCC-spd` is `fit/best_D.json` converted to Rydberg. A has a slightly lower cost, but D
reproduces the minority DOS at E_F, the part that matters most for transport, much better
(`cmpfcc.py`; DOS per spin in states/eV, Bethe / CRYSTAL23, relative L2 misfit from -10 to +4 eV):

| | A: DOS(E_F) up | down | D: DOS(E_F) up | down | D: misfit up / down |
|---|---|---|---|---|---|
| s | 0.062 / 0.040 | 0.024 / 0.017 | 0.054 / 0.040 | 0.022 / 0.017 | 0.60 / 0.43 |
| p (not fitted) | 0.018 / 0.031 | 0.044 / 0.059 | 0.029 / 0.031 | 0.045 / 0.059 | 1.28 / 0.73 |
| d | 0.053 / 0.024 | 0.603 / 0.823 | 0.063 / 0.024 | 0.788 / 0.823 | 0.25 / 0.34 |
| total | 0.134 / 0.096 | 0.671 / 0.899 | 0.147 / 0.096 | 0.855 / 0.899 | 0.24 / 0.31 |

`dos_fit.png` compares the d and total DOS of A and D with CRYSTAL23.

In D the three d-d hoppings sit at their bounds: dds = -0.10 eV, ddp = +0.05 eV and
ddd = -0.01 eV. The width and shape of the d band come instead from the s-d and p-d hoppings and
the e_g/t_2g on-site splitting. `FCC-sd` also has dds and ddd at their bounds. Neither fit
reproduces the sharp peaks of the CRYSTAL23 DOS.

## Reproducing

- Run CRYSTAL23 with `crystal/INPUT`, then `properties` with `crystal/d3w` (writes `DOSS.DAT`;
  `DOSS.win.DAT` is a copy of it).
- `cd fit; python3 fitfcc.py D` repeats the fit and writes `best_D.json` (needs numpy and
  scipy; about 6 s per call, 600 calls). `python3 cmpfcc.py D` prints the comparison above.

CRYSTAL23 1.0.1 (OpenMP, ifx) on 4 threads; fit with numpy 1.26 and scipy 1.13; 8 November 2026.
