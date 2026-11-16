Au4-Au-Au--Au4 (10 atoms) with one s function (and ECP) on the 8 outer atoms and LANL2DZ on the apex
atoms, between the 1-orbital `FCC-s` Bethe-lattice electrodes, LSDA.

`.ant`: `ALPHA = 0.02`, `BLPAR1 = BLPAR2 = FCC-s`. The `FCC-s` set (fitted for this s basis and ECP,
`BLDAT/fits/Au-FCC-s`) matches the s-only electrode planes.

## Outputs

`Au4-Au-Au--Au4.log`, `T.Au4-Au-Au--Au4.dat` (E - E_F from -3 to 3 eV), `Au4-Au-Au--Au4.xyz`.

## Results

ANT.Gaussian 3.0.0 (commit 7709e49 sources), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads, 13 November 2026, from scratch.

| Quantity | Value |
|---|---|
| Basis functions | 52 |
| SCF energy (a.u.) | -272.526621 |
| SCF cycles | 14 |
| Device Fermi energy (eV) | -2.5120 |
| T(E_F) (e^2/h) | 1.578 |
| Wall time | 1 s |

The Bethe-lattice Fermi energy of `FCC-s` is 0 eV.

## Notes

- The input had no `.ant` and the title `Au-nanocontact`. ANT takes the job name from the title, so
  it read no `.ant` and used the 9-orbital default `Papacon`, which stops with `NAOrbs read from BL
  data base not equal to NAOrbs of electrode plane`. The title is now the file name and the `.ant`
  selects `FCC-s`.
- With the earlier 1-orbital set `lsda-6s`, now removed, it gave E = -272.550088 a.u.,
  E_F = -2.4688 eV and T(E_F) = 1.553.
