Au nanocontact Au9-Au4-Au-Au4-Au9 (27 atoms) between Bethe-lattice electrodes, LSDA: 9-atom and
4-atom layers on each side of a single apex atom.

Each subdirectory is one variation of the calculation for this structure, with its own inputs,
outputs and README. The files keep the name of the system, `Au9-Au4-Au-Au4-Au9`.

- [spd-Papacon/](spd-Papacon/README.md): `SL = 1.0d-3`, `OVERLAP = -1.0`, `nosymm`.
- [spd-Papacon-noant/](spd-Papacon-noant/README.md): the same input without `.ant` (ANT defaults)
  and without `nosymm`; it was the reference contact of
  [../../benzene/input_Electrodo_Au_3_benzene](../../benzene/input_Electrodo_Au_3_benzene/README.md).

Both use a minimal spd basis with a 68-electron ECP, and LANL2DZ on the apex atom (atom 14).

## Comparison

| Variation | E (a.u.) | Cycles | E_F (eV) | T(E_F) (e^2/h) |
|---|---|---|---|---|
| [spd-Papacon/](spd-Papacon/README.md) | -1041.555335 | 16 | -8.0647 | 2.335 |
| [spd-Papacon-noant/](spd-Papacon-noant/README.md) | -1041.523349 | 37 | -8.0632 | 2.333 |
