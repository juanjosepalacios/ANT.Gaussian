Au nanocontacts between Bethe-lattice electrodes, LSDA, with the book `Papacon` set of
`BLDAT/BL079.dat` (the default electrode set).

## Jobs

| Job | Structure | `.ant` settings | Notes |
|---|---|---|---|
| `Au4-Au-Au-Au4` | 10 atoms | `ALPHA = 0.02` | starts from `P.Au4-Au-Au-Au4.dat` if present |
| `Au9-Au4-Au-Au4-Au9` | 27 atoms, one apex atom | `SL = 1.0d-3`, `OVERLAP = -1.0` | starts from `P.Au9-Au4-Au-Au4-Au9.dat` if present |
| `Au9-Au4-Au-Au-Au4-Au9` | 28 atoms | none (no `.ant`, defaults) | |
| `Au9-Au4-Au--Au4-Au9`, `Au9-Au4--Au-Au4-Au9` | 27 atoms with a gap on either side of the apex atom | `BIASVOLT = -3.0`, `SL = 1.0d-3`, `MULLIKEN` | |
| `Au4-Au` | 5 atoms | ghost electrodes | not rerun; no Bethe lattice |
| `Au4-Au-Au--Au4` | 10 atoms, one s function on the 8 outer atoms, LANL2DZ on the apex atoms | `ALPHA = 0.02`, `BLPAR1 = BLPAR2 = lsda-6s` | the 1-orbital `lsda-6s` electrodes match the s-only electrode planes |

The `.gjf` files substitute links from `../bin` or `../../bin`; change these and `%NProcShared` to
your installation. ANT takes the job name from the title line of the `.gjf`, which matches the file
names here.

## Results

ANT.Gaussian 3.0.0 (commit 7709e49 sources), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads,
13 November 2026, from scratch (the stored `P` files were not used). The outputs are
`<job>.log`, `T.<job>.dat` (E - E_F from -3 to 3 eV), `<job>.xyz`, and for the bias jobs
`Q.<job>.dat` (Mulliken charges).

| Job | Basis functions | SCF energy (a.u.) | Cycles | Device E_F (eV) | T(E_F) (e^2/h) | Wall time |
|---|---|---|---|---|---|---|
| Au4-Au-Au-Au4 | 90 | -348.018120 | 12 | -8.3298 | 1.261 | 5 s |
| Au9-Au4-Au-Au4-Au9 | 256 | -1041.555335 | 16 | -8.0647 | 2.335 | 125 s |
| Au9-Au4-Au-Au-Au4-Au9 | 252 | -975.625433 | 17 | -7.9370 | 1.060 | 34 s |
| Au9-Au4-Au--Au4-Au9 (-3 V) | 243 | -940.559212 | 27 | -8.0060 | 0.069 | 208 s |
| Au9-Au4--Au-Au4-Au9 (-3 V) | 243 | -940.555680 | 27 | -7.9672 | 0.171 | 215 s |
| Au4-Au-Au--Au4 (`lsda-6s`, 2 threads) | 52 | -272.550088 | 14 | -2.4688 | 1.553 | 2 s |

The Bethe-lattice Fermi energy is 6.931 eV with `Papacon` and -4.971 eV with `lsda-6s`. All runs terminate normally; the two bias
runs report 30 and 16 unsuccessful density integrations, the others none.

## Notes

- The previous outputs came from ANT 2.5.2 with the earlier Au `Papacon` set (book values with `edd`
  and `edt` swapped). For Au4-Au-Au-Au4 the stored log had a Bethe-lattice Fermi energy of 2.709 eV,
  E = -347.701106 a.u. and T(E_F) = 0.862; the differences also include all the changes to ANT
  since 2.5.2.
- `Au9-Au4-Au-Au4-Au9.ant` had `PFIX` with `P.Au9-Au4-Au-Au4-Au9.dat0`, which is not in the
  repository; the `PFIX` lines are removed.
- `Au4-Au-Au--Au4` had no `.ant` and the title `Au-nanocontact`. ANT takes the job name from the
  title, so it read no `.ant` and used the 9-orbital default `Papacon`, which stops with `NAOrbs read
  from BL data base not equal to NAOrbs of electrode plane`. The title is now the file name and the
  new `.ant` selects `lsda-6s`.
- `P.Au4-Au-Au-Au4.dat` and `P.Au9-Au4-Au-Au4-Au9.dat` are the densities of these runs: a new run in
  this directory starts from them.
