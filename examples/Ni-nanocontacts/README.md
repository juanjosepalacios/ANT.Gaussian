Ni nanocontacts between Bethe-lattice electrodes, spin-polarized LSDA, in two geometries:

- Ni4-Ni-Ni-Ni4 (10 atoms): two square layers of 4 atoms at z = -3.065 and +3.065 A joined by two
  apex atoms on the axis at z = -1.30 and +1.30 A. Charge 0, starting multiplicity 15.
- Ni9-Ni4-Ni-Ni-Ni4-Ni9 (28 atoms): 9-atom outer layers at z = -4.75 and +4.75 A, 4-atom layers at
  -3.00 and +3.00 A and apex atoms at -1.25 and +1.25 A (the `doll` case uses a slightly more open
  version, see its README). Charge 0, starting multiplicity 31.

Each subdirectory is one variation of the calculation, with its own inputs, outputs and README:

- [Ni4-LANL1MB-Papacon/](Ni4-LANL1MB-Papacon/README.md): Ni4-Ni-Ni-Ni4 in a Zeeman field of 20 T
  along z, minimal spd basis LANL1MB with the Hay-Wadt large-core ECP on all atoms, Bethe-lattice
  electrodes `Papacon`.
- [Ni9-LANL1MB-Papacon/](Ni9-LANL1MB-Papacon/README.md): Ni9-Ni4-Ni-Ni-Ni4-Ni9 with parallel
  electrodes, LANL1MB basis, `Papacon` electrodes.
- [Ni9-spd-F/](Ni9-spd-F/README.md): Ni9-Ni4-Ni-Ni-Ni4-Ni9 with parallel electrodes, CRENBS s+d
  basis and ECP plus a 4p shell, earlier `Papacon` electrodes (removed from `BL028.dat`).
- [Ni9-spd-AF/](Ni9-spd-AF/README.md): as Ni9-spd-F with antiparallel electrodes and a domain wall in
  the constriction, started from the Ni9-spd-F density.
- [Ni9-doll/](Ni9-doll/README.md): Ni9-Ni4-Ni-Ni-Ni4-Ni9 with parallel electrodes, CRENBS sd basis
  and ECP on the outer layers and an all-electron basis on the central atoms 10-19, electrodes
  `Papacon-sd`.

`Papacon` is the orthogonal first-neighbour Slater-Koster fit of ferromagnetic fcc Ni from
Papaconstantopoulos, Handbook of the Band Structure of Elemental Solids, 2nd ed. (2015). The
Ni9-spd-F and Ni9-spd-AF results were computed with an earlier non-orthogonal set of that name and
cannot be reproduced with the current `BL028.dat`.

## Comparison

| Variation | Basis and electrodes | E_F (eV) | T(E_F) up + down (e^2/h) | G (G0) |
|---|---|---|---|---|
| [Ni4-LANL1MB-Papacon/](Ni4-LANL1MB-Papacon/README.md) | spd LANL1MB, `Papacon`, 20 T | -9.052 | 2.212 (spinor) | 1.11 |
| [Ni9-LANL1MB-Papacon/](Ni9-LANL1MB-Papacon/README.md) | spd LANL1MB, `Papacon` | -8.752 | 0.391 + 2.528 = 2.919 | 1.46 |
| [Ni9-spd-F/](Ni9-spd-F/README.md) | spd CRENBS, earlier `Papacon` | -4.364 | 0.837 + 2.864 = 3.700 | 1.85 |
| [Ni9-spd-AF/](Ni9-spd-AF/README.md) | spd CRENBS, earlier `Papacon`, antiparallel | -4.280 | 1.202 + 1.227 = 2.429 | 1.21 |
| [Ni9-doll/](Ni9-doll/README.md) | sd CRENBS + all-electron centre, `Papacon-sd` | -6.428 | 0.871 + 2.548 = 3.418 | 1.71 |

The device Fermi energies are on the energy scale of each electrode set, so they are comparable only
between runs with the same electrodes. G = T_total * e^2/h = (T_total / 2) G0.
