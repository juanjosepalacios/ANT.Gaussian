Small Au nanocontact Au4-Au-Au-Au4 (10 atoms) between Bethe-lattice electrodes, LSDA: two 4-atom
layers joined by two apex atoms. The gapped version is [../Au4-Au-Au--Au4](../Au4-Au-Au--Au4/README.md).

Each subdirectory is one variation of the calculation for this structure, with its own inputs,
outputs and README. The files keep the name of the system, `Au4-Au-Au-Au4`.

- [spd-Papacon/](spd-Papacon/README.md): minimal spd basis and ECP on all atoms, `Papacon` electrodes.
- [s-FCC-s/](s-FCC-s/README.md): one s function on the 8 layer atoms and LANL2DZ on the apex atoms,
  1-orbital `FCC-s` electrodes. Not run.

## Comparison

| Variation | Basis and electrodes | E_F (eV) | T(E_F) (e^2/h) |
|---|---|---|---|
| [spd-Papacon/](spd-Papacon/README.md) | spd, `Papacon` | -8.3298 | 1.261 |
| [s-FCC-s/](s-FCC-s/README.md) | s + LANL2DZ, `FCC-s` | not run | |
