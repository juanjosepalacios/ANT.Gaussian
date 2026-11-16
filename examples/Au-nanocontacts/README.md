Au nanocontacts between Bethe-lattice electrodes, LSDA. Unless noted, the electrodes are the book
`Papacon` set of `BLDAT/BL079.dat` (the default Au set) and the basis is a minimal spd set with a
68-electron ECP. Au chains between 1D leads are in [../Au-chain](../Au-chain/README.md).

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation (basis, electrodes, settings), with its own inputs, outputs and
README. The files keep the name of the system.

- [Au4-Au-Au-Au4/](Au4-Au-Au-Au4/README.md): 10 atoms, two 4-atom layers joined by two apex atoms
  - [spd-Papacon/](Au4-Au-Au-Au4/spd-Papacon/README.md)
  - [s-FCC-s/](Au4-Au-Au-Au4/s-FCC-s/README.md)
- [Au4-Au-Au--Au4/](Au4-Au-Au--Au4/README.md): the same 10 atoms with a gap between the apex atoms
  - [s-FCC-s/](Au4-Au-Au--Au4/s-FCC-s/README.md)
- [Au4-Au/](Au4-Au/README.md): 5 atoms, one 4-atom layer and an apex atom
  - [spd-ghost/](Au4-Au/spd-ghost/README.md)
- [Au9-Au4-Au-Au4-Au9/](Au9-Au4-Au-Au4-Au9/README.md): 27 atoms, one apex atom
  - [spd-Papacon/](Au9-Au4-Au-Au4-Au9/spd-Papacon/README.md)
  - [spd-Papacon-noant/](Au9-Au4-Au-Au4-Au9/spd-Papacon-noant/README.md)
- [Au9-Au4-Au-Au-Au4-Au9/](Au9-Au4-Au-Au-Au4-Au9/README.md): 28 atoms, two apex atoms
  - [spd-Papacon/](Au9-Au4-Au-Au-Au4-Au9/spd-Papacon/README.md)
- [Au9-Au4-Au--Au4-Au9/](Au9-Au4-Au--Au4-Au9/README.md) and
  [Au9-Au4--Au-Au4-Au9/](Au9-Au4--Au-Au4-Au9/README.md): 27 atoms with a gap on either side of the
  apex atom, at -3 V
  - [spd-Papacon-bias/](Au9-Au4-Au--Au4-Au9/spd-Papacon-bias/README.md)
  - [spd-Papacon-bias/](Au9-Au4--Au-Au4-Au9/spd-Papacon-bias/README.md)

The `.gjf` files substitute links from `../../../../bin`; change this and `%NProcShared` to your
installation. ANT takes the job name from the title line of the `.gjf`, which matches the file
names.

## Results

ANT.Gaussian 3.0.0 (commit 7709e49 sources), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads
(8 for spd-Papacon-noant), 13 November 2026, from scratch.

| System / variation | Basis functions | Device E_F (eV) | T(E_F) (e^2/h) |
|---|---|---|---|
| [Au4-Au-Au-Au4/spd-Papacon](Au4-Au-Au-Au4/spd-Papacon/README.md) | 90 | -8.3298 | 1.261 |
| [Au4-Au-Au--Au4/s-FCC-s](Au4-Au-Au--Au4/s-FCC-s/README.md) | 52 | -2.5120 | 1.578 |
| [Au9-Au4-Au-Au4-Au9/spd-Papacon](Au9-Au4-Au-Au4-Au9/spd-Papacon/README.md) | 256 | -8.0647 | 2.335 |
| [Au9-Au4-Au-Au4-Au9/spd-Papacon-noant](Au9-Au4-Au-Au4-Au9/spd-Papacon-noant/README.md) | 256 | -8.0632 | 2.333 |
| [Au9-Au4-Au-Au-Au4-Au9/spd-Papacon](Au9-Au4-Au-Au-Au4-Au9/spd-Papacon/README.md) | 252 | -7.9370 | 1.060 |
| [Au9-Au4-Au--Au4-Au9/spd-Papacon-bias](Au9-Au4-Au--Au4-Au9/spd-Papacon-bias/README.md) | 243 | -8.0060 | 0.069 |
| [Au9-Au4--Au-Au4-Au9/spd-Papacon-bias](Au9-Au4--Au-Au4-Au9/spd-Papacon-bias/README.md) | 243 | -7.9672 | 0.171 |

The Bethe-lattice Fermi energy is 6.931 eV with `Papacon` and 0 eV with `FCC-s`. The device Fermi
energies are on the energy scale of the electrodes.

The previous outputs came from ANT 2.5.2 with the earlier Au `Papacon` set (book values with `edd`
and `edt` swapped).
