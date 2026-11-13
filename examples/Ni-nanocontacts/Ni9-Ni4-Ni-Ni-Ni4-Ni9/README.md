Ni nanocontact Ni9-Ni4-Ni-Ni-Ni4-Ni9 (28 atoms) between Bethe-lattice electrodes, spin-polarized
LSDA.

Two 9-atom outer layers at z = -4.75 and +4.75 A, two 4-atom layers at z = -3.00 and +3.00 A and two
apex atoms on the axis at z = -1.25 and +1.25 A (fcc stacking with a = 3.50 A). The charge is 0 and
the starting multiplicity is 31. The `doll-Papacon-sd` and `doll-FCC-spd` variations use the same
arrangement with a = 3.53 A (layers at +-4.83, +-3.065 and +-1.30 A) and starting multiplicity 37.

Each subdirectory is one variation of the calculation for this contact, with its own inputs, outputs
and README. The files keep the name of the system, `Ni9-Ni4-Ni-Ni-Ni4-Ni9`.

- [LANL1MB-Papacon-F/](LANL1MB-Papacon-F/README.md): parallel electrodes, minimal spd basis
  LANL1MB with the Hay-Wadt large-core ECP on all atoms, Bethe-lattice electrodes `Papacon`.
  Includes the converged density matrix used to start the AF case.
- [LANL1MB-Papacon-AF/](LANL1MB-Papacon-AF/README.md): as LANL1MB-Papacon-F with antiparallel
  electrodes and a domain wall in the constriction, started from its density.
- [CRENBS-spd-F/](CRENBS-spd-F/README.md): parallel electrodes, CRENBS s+d basis and ECP plus a 4p
  shell, earlier `Papacon` electrodes that have been removed from `BLDAT/BL028.dat`, so its results
  cannot be reproduced with the current file.
- [doll-Papacon-sd/](doll-Papacon-sd/README.md): parallel electrodes, CRENBS sd basis and ECP on
  the outer layers and an all-electron basis on the central atoms 10-19, electrodes `Papacon-sd`.
- [doll-FCC-spd/](doll-FCC-spd/README.md): as doll-Papacon-sd with the p-optimized CRENBS spd basis
  on the outer layers and electrodes `FCC-spd`.

## Comparison

| Variation | Basis and electrodes | E_F (eV) | T(E_F) up + down (e^2/h) | G (G0) |
|---|---|---|---|---|
| [LANL1MB-Papacon-F/](LANL1MB-Papacon-F/README.md) | spd LANL1MB, `Papacon` | -8.752 | 0.391 + 2.528 = 2.919 | 1.46 |
| [LANL1MB-Papacon-AF/](LANL1MB-Papacon-AF/README.md) | spd LANL1MB, `Papacon`, antiparallel | -8.759 | 0.622 + 0.622 = 1.244 | 0.62 |
| [CRENBS-spd-F/](CRENBS-spd-F/README.md) | spd CRENBS, earlier `Papacon` | -4.364 | 0.837 + 2.864 = 3.700 | 1.85 |
| [doll-Papacon-sd/](doll-Papacon-sd/README.md) | sd CRENBS + all-electron centre, `Papacon-sd` | -6.428 | 0.871 + 2.548 = 3.418 | 1.71 |
| [doll-FCC-spd/](doll-FCC-spd/README.md) | spd CRENBS + all-electron centre, `FCC-spd` | -8.06 | 1.022 + 2.663 = 3.685 | 1.84 |

The device Fermi energies are on the energy scale of each electrode set, so they are comparable only
between runs with the same electrodes. G = T_total * e^2/h = (T_total / 2) G0.
