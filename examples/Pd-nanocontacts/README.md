Pd nanocontact Pd9-Pd4-Pd-Pd-Pd4-Pd9 (28 atoms) between Bethe-lattice electrodes, LSDA.

## Settings

`scf=(conver=5,maxcycle=200,dsymm,nodamp)`; `.ant`: `ALPHA = 0.03`, `NPULAY = 4`, `FERMIACC` and
`CHARGEACC` 1.0d-5, `PACC = 1.0d-8`, `SL = 1.0d-2`, `LEADDOS`, `BLPAR1 = Papacon-old`.

Only `BLPAR1` is set, so the left electrode uses `Papacon-old` of `BLDAT/BL046.dat` and the right
one the default `Papacon`, which is now the book set (Papaconstantopoulos, Table 5.64). The two
electrodes are therefore different. The `.gjf` substitutes links from `../bin`; change this and
`%NProcShared` to your installation.

## Results

ANT.Gaussian 3.0.0 (commit 7709e49 sources), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads,
13 November 2026, from scratch. Outputs: `Pd9-Pd4-Pd-Pd-Pd4-Pd9.log`, `T.Pd9-Pd4-Pd-Pd-Pd4-Pd9.dat`
(E - E_F from -3 to 3 eV), `Pd9-Pd4-Pd-Pd-Pd4-Pd9.xyz`, `Lead1DOS.dat` and `Lead2DOS.dat`.

| Quantity | Value |
|---|---|
| Basis functions | 252 |
| SCF energy (a.u.) | -854.728053 |
| SCF cycles | 31 |
| Device Fermi energy (eV) | -7.8097 |
| T(E_F) (e^2/h) | 0.044 |
| Unsuccessful density integrations | 148 |
| Wall time | 110 s |

## Notes

- There were no stored outputs before.
- The result is not a sensible Pd contact: T(E_F) = 0.044 and 148 unsuccessful density
  integrations. The `Papacon-old` set of the left electrode has d on-site energies near 0.02 Ry, far
  from those of the book set on the right. For a symmetric contact with the book parameters, remove
  the `BLPAR1` line or set `BLPAR2` as well.
