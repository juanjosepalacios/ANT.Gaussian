Mn nanocontacts between Bethe-lattice electrodes, spin-polarized (UBPBE), with the book `Papacon`
set of `BLDAT/BL025.dat` (bcc paramagnetic Mn, the default electrode set).

## Jobs

| Job | Structure | Settings |
|---|---|---|
| `Mn9-Mn4-Mn-Mn-Mn4-Mn9` | 28 atoms, multiplicity 141 | `scf=(conver=4,maxcycle=200,dsymm)`; `ALPHA = 0.020`, `NPULAY = 4`, accuracies 1.0d-5/1.0d-8, `SL = 2.0d-2`, `MAX = 5` |
| `Mn9-Mn4-Mn-Mn-Mn4-Mn9.lan` | the same with LANL2DZ on atoms 10-19 | as above, `nodamp`; `PFIX` with `P.Mn9-Mn4-Mn-Mn-Mn4-Mn9.dat` and `NFIX = 10` (atoms 10-19) |
| `Mn14-Mn14` | 28 atoms, multiplicity 141 | `conver=6`; `ALPHA = 0.03`, accuracies 1.0d-5/1.0d-8, `MAX = 4`, `NSPINLOCK = 4`, `LDOS_BEG/END` 14/15 |
| `Mn9-Mn4-Mn` | half contact given as a z-matrix, multiplicity 71, geometry optimization, ghost second electrode | `opt(z-matrix)`; `ALPHA = 0.030`, `NPULAY = 3` |

Run `Mn9-Mn4-Mn-Mn-Mn4-Mn9` before the `.lan` job: it writes `P.Mn9-Mn4-Mn-Mn-Mn4-Mn9.dat`, which the
`.lan` job reads through `PFIX` (not stored here, 3.5 MB). The `.gjf` files substitute links from
`../bin`; change this and `%NProcShared` to your installation.

## Results

ANT.Gaussian 3.0.0 (commit 7709e49 sources), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads,
13 November 2026, from scratch. Outputs: `<job>.log`, `T.<job>.dat` (spin up and spin down, E - E_F
from -3 to 3 eV), `<job>.xyz`.

| Job | Basis functions | SCF energy (a.u.) | Cycles | Device E_F (eV) | T(E_F) up / down (e^2/h) | Unsuccessful integrations | Wall time |
|---|---|---|---|---|---|---|---|
| Mn9-Mn4-Mn-Mn-Mn4-Mn9 | 252 | -438.722578 | 22 | -5.5320 | 0.069 / 0.001 | 9 | 152 s |
| Mn9-Mn4-Mn-Mn-Mn4-Mn9.lan | 382 | -1321.088073 | 28 | -5.5398 | 0.061 / 0.046 | 7 | 401 s |
| Mn14-Mn14 | 252 | -438.562168 | 65 | -5.3747 | 0.070 / 0.002 | 157 | 240 s |

The Bethe-lattice Fermi energy is 9.976 eV. The three runs terminate normally.

## Notes

- There were no stored outputs before.
- `Mn9-Mn4-Mn` (geometry optimization) did not finish: after 101 SCF cycles of its first
  geometry ANT stopped with `Bisec/Error: Root is not bracketed`. No outputs are stored for it.
- The book set is for paramagnetic bcc Mn with one spin; ANT uses it for both spins, so the
  electrodes are not magnetic. The geometries are fcc-like pyramids, so the bcc parameters are only
  an approximation for these electrodes.
- Mn14-Mn14 reports 157 unsuccessful density integrations (point limit), and the transmissions are
  small; treat these numbers with care.
