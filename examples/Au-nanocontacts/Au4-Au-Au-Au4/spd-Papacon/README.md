Au4-Au-Au-Au4 (10 atoms) between `Papacon` Bethe-lattice electrodes (`BLDAT/BL079.dat`), LSDA, with
a minimal spd basis and a 68-electron ECP on all atoms.

`.ant`: `ALPHA = 0.02`. A run in this directory starts from `P.Au4-Au-Au-Au4.dat` if present.

## Outputs

`Au4-Au-Au-Au4.log`, `T.Au4-Au-Au-Au4.dat` (E - E_F from -3 to 3 eV), `Au4-Au-Au-Au4.xyz` and
`P.Au4-Au-Au-Au4.dat`, the density of this run.

## Results

ANT.Gaussian 3.0.0 (commit 7709e49 sources), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads, 13 November 2026, from scratch.

| Quantity | Value |
|---|---|
| Basis functions | 90 |
| SCF energy (a.u.) | -348.018120 |
| SCF cycles | 12 |
| Device Fermi energy (eV) | -8.3298 |
| T(E_F) (e^2/h) | 1.261 |
| Wall time | 5 s |

The previous ANT 2.5.2 output, with the earlier Au `Papacon` set, had a Bethe-lattice Fermi energy of
2.709 eV, E = -347.701106 a.u. and T(E_F) = 0.862.
