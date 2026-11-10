How the `FCC-sd` Ni Bethe-lattice set in `BLDAT/BL028.dat` was made: a CRYSTAL23 LSDA
calculation of bulk fcc Ni with the CRENBS sd basis and ECP, and a fit of the Bethe-lattice
parameters to its DOS, in the same way as the Fe sets (`../Fe-BCC-spd/README.md`,
`../Fe-FCC-spd/README.md`). Use the set with `BLPAR1/2 = FCC-sd` and the basis and ECP in
`Ni-crenbs-sd.gbs`, the ones of `examples/Ni-nanocontacts`.

## Basis (`Ni-crenbs-sd.gbs`)

CRENBS s and d shells and ECP (Ar core, 10 valence electrons) from the Basis Set Exchange. The Fe
CRENBS entries there are identical to the Fe sd basis and ECP used for `BCC-sd` and `FCC-sd` of Fe.

## CRYSTAL23 (`crystal/`)

Ferromagnetic fcc Ni, a = 3.524 A (experiment), SVWN, the settings of the Fe runs (`SHRINK 24 24`,
`SMEAR 0.01`, `FMIXING 80`, `BROYDEN`, `TOLDEE 8`). The spin is locked for the first 5 cycles at
`SPINLOCK 2 5` (10 electrons need an even N_alpha - N_beta; Fe used 4).

| SCF cycles | E (Ha) | Moment (muB) | E_F (Ha) |
|---|---|---|---|
| 19 | -40.081566 | 0.699 | 0.20362 |

`d3w` writes the s, d and total DOS of both spins from E_F - 0.58 to E_F + 0.57 Ha
(`DOSS.win.DAT`, the window of the Fe runs). The integrated moment of this DOS is 0.73 muB; the fit
uses the Mulliken moment, 0.699 muB (`moment`), as for Fe. E_F lies on the upper edge of the
minority d band: DOS(E_F) is 1.52 states/eV for spin down and 0.10 for spin up.

## Fit (`fit/`)

The parameters are fitted on a 12-direction (fcc, <110>) Bethe lattice, orthogonal, without
overlaps, with the targets, weights and bounds of `BCC-sd`: the s (weight 3) and d DOS of both spins
from -10 to +4 eV around E_F and the moment (weight 5), 10 electrons. Hoppings are the same for both
spins, the spin-up d on-site energy is fixed at 1.3812 eV as in the Fe sets (ANT shifts all on-site
energies to put E_F at zero, so this only fixes the reference), and ETA = 0.1 eV. `fitsd.py` uses
the model of `../Fe-BCC-spd/fit/bethe9.py` with the p shell decoupled (no p hoppings, p on-site
energy at 200 eV), which is the s+d Green's function.

| Start | From | Cost | Moment (muB) | Note |
|---|---|---|---|---|
| A | Fe `FCC-sd` with d spin splittings of 0.6 eV and s of 0.2 eV | 89.843 | 0.769 | `dds` at its bound |
| B | Ni `Papacon-old` (BL028.dat, orthogonal), spin-averaged hoppings | 89.600 | 0.684 | `dds` at its bound |
| C | B with the s-s and s-d hoppings of Fe `FCC-sd` | 86.719 | 0.757 | no parameter at a bound; used for `FCC-sd` |
| D | A with the d-d hoppings of B | 90.081 | 0.698 | |

`FCC-sd` is `fit/best_C.json` (copied to `best_sd.json`, the start of the spd fit) converted to
Rydberg. C has the lowest cost and reproduces the minority DOS at E_F best (`cmpsd.py`; DOS per spin
in states/eV, Bethe / CRYSTAL23, relative L2 misfit from -10 to +4 eV):

| | DOS(E_F) up | down | Misfit up / down |
|---|---|---|---|
| s | 0.058 / 0.025 | 0.027 / 0.011 | 0.92 / 0.91 |
| d | 0.182 / 0.077 | 1.381 / 1.508 | 0.28 / 0.25 |
| total | 0.239 / 0.102 | 1.408 / 1.520 | 0.34 / 0.31 |

D reproduces the moment exactly but gives 1.10 states/eV for the minority DOS at E_F.
`dos_fit.png` compares C and D with CRYSTAL23. The d band width, its position relative to E_F and
the exchange splitting are reproduced; its peaks are not. The s misfit is large because the
CRYSTAL23 s DOS is small and dominated by the band-bottom peak at -8 eV and the ripples of the
Legendre expansion, which a Bethe lattice cannot reproduce.

Fermi level (`efcheck.py`): the occupations below E_F are N_s = 0.65 / 0.56 (Bethe / CRYSTAL23)
and N_d = 4.73 / 4.77 for spin up, and N_s = 0.56 / 0.57 and N_d = 4.06 / 4.04 for spin down. A rigid
shift of the Bethe DOS does not lower the d and total misfit (best shift 0.00 eV), so E_F of the
fit and of CRYSTAL23 are consistent. With ETA = 0.1 eV the fit has E_F = 3.550 eV on the scale
of the on-site energies.

## Check in ANT

ANT with `FCC-sd` on both electrodes of the Ni4-Ni-Ni-Ni4 geometry (now `examples/Ni-nanocontacts/Ni4-LANL1MB-Papacon`)
(`LEADDOS`, links from `bin/`, 8 November 2026) finds 12 <110> directions and a Bethe-lattice
E_F of 3.467 eV (ETA -> 0, so slightly below the fit value). Its `Lead1DOS.dat` agrees with
`bethe9.py` at the same E_F and ETA = 0.001 eV to 1e-4 states/eV; DOS(E_F) is 0.218 (up) and 1.569
(down) states/eV.

## Reproducing

- Run CRYSTAL23 with `crystal/INPUT`, then `properties` with `crystal/d3w` (writes `DOSS.DAT`;
  `DOSS.win.DAT` is a copy of it).
- `cd fit; python3 fitsd.py C` repeats the fit and writes `best_C.json` (needs numpy and scipy;
  about 4 s per call, 340 calls). `python3 cmpsd.py C` prints the comparison above and
  `python3 efcheck.py C` the Fermi-level check. `python3 plotsd.py C D` writes `dos_fit.npz`, and
  `plotsd.py --plot`, run with a Python that has matplotlib, draws `dos_fit.png` from it.

CRYSTAL23 1.0.1 (OpenMP, ifx) on 4 threads; fit with numpy 1.19 and scipy 1.5; 8 November 2026.
