Au9-Au4-Au-Au4-Au9 (27 atoms) between `Papacon` Bethe-lattice electrodes, LSDA.

`.ant`: `SL = 1.0d-3`, `OVERLAP = -1.0`. The route has `nosymm`. A run in this directory starts from
`P.Au9-Au4-Au-Au4-Au9.dat` if present.

## Outputs

`Au9-Au4-Au-Au4-Au9.log`, `T.Au9-Au4-Au-Au4-Au9.dat` (E - E_F from -3 to 3 eV),
`Au9-Au4-Au-Au4-Au9.xyz` and `P.Au9-Au4-Au-Au4-Au9.dat`, the density of this run.

## Results

ANT.Gaussian 3.0.0 (commit 7709e49 sources), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads, 13 November 2026, from scratch.

| Quantity | Value |
|---|---|
| Basis functions | 256 |
| SCF energy (a.u.) | -1041.555335 |
| SCF cycles | 16 |
| Device Fermi energy (eV) | -8.0647 |
| T(E_F) (e^2/h) | 2.335 |
| Wall time | 125 s |

The `.ant` had `PFIX` with `P.Au9-Au4-Au-Au4-Au9.dat0`, which is not in the repository; the `PFIX`
lines are removed.
