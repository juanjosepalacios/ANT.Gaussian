Small Ni nanocontact Ni4-Ni-Ni-Ni4 (10 atoms) between Bethe-lattice electrodes, spin-polarized LSDA.

Two square layers of 4 atoms at z = -3.065 and +3.065 A, with the Ni-Ni distance of fcc Ni
(2.496 A, a = 3.530 A), are joined by two apex atoms on the axis at z = -1.30 and +1.30 A. The charge
is 0 and the starting multiplicity is 15.

Each subdirectory is one variation of the calculation for this contact, with its own inputs, outputs
and README. The files keep the name of the system, `Ni4-Ni-Ni-Ni4`.

- [LANL1MB-Papacon/](LANL1MB-Papacon/README.md): minimal spd basis LANL1MB with the Hay-Wadt
  large-core ECP on all atoms, Bethe-lattice electrodes `Papacon`, Zeeman field of 20 T along z.

## Comparison

| Variation | Basis and electrodes | E_F (eV) | T(E_F) total (e^2/h) | G (G0) |
|---|---|---|---|---|
| [LANL1MB-Papacon/](LANL1MB-Papacon/README.md) | spd LANL1MB, `Papacon`, 20 T | -9.052 | 2.212 (spinor) | 1.11 |
