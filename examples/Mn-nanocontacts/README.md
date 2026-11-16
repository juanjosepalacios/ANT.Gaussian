Mn nanocontacts between Bethe-lattice electrodes, spin-polarized (UBPBE), with the book `Papacon`
set of `BLDAT/BL025.dat` (bcc paramagnetic Mn, the default electrode set) and an spd basis with ECP.

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation (basis, electrodes, settings), with its own inputs, outputs and
README. The files keep the name of the system.

- [Mn9-Mn4-Mn-Mn-Mn4-Mn9/](Mn9-Mn4-Mn-Mn-Mn4-Mn9/README.md): 28 atoms, two apex atoms
  - [spd-Papacon/](Mn9-Mn4-Mn-Mn-Mn4-Mn9/spd-Papacon/README.md)
  - [spd-LANL2DZ-Papacon-PFIX/](Mn9-Mn4-Mn-Mn-Mn4-Mn9/spd-LANL2DZ-Papacon-PFIX/README.md)
- [Mn14-Mn14/](Mn14-Mn14/README.md): 28 atoms, two 14-atom halves
  - [spd-Papacon/](Mn14-Mn14/spd-Papacon/README.md)
- [Mn9-Mn4-Mn/](Mn9-Mn4-Mn/README.md): half contact, geometry optimization
  - [spd-ghost-opt/](Mn9-Mn4-Mn/spd-ghost-opt/README.md)

The `.gjf` files substitute links from `../../../../bin`; change this and `%NProcShared` to your
installation. ANT takes the job name from the title line of the `.gjf`, which matches the file
names.

## Results

ANT.Gaussian 3.0.0 (commit 7709e49 sources), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads, 13 November 2026, from scratch. Spin up / down.

| System / variation | Basis functions | Device E_F (eV) | T(E_F) up / down (e^2/h) |
|---|---|---|---|
| [Mn9-Mn4-Mn-Mn-Mn4-Mn9/spd-Papacon](Mn9-Mn4-Mn-Mn-Mn4-Mn9/spd-Papacon/README.md) | 252 | -5.5320 | 0.069 / 0.001 |
| [Mn9-Mn4-Mn-Mn-Mn4-Mn9/spd-LANL2DZ-Papacon-PFIX](Mn9-Mn4-Mn-Mn-Mn4-Mn9/spd-LANL2DZ-Papacon-PFIX/README.md) | 382 | -5.5398 | 0.061 / 0.046 |
| [Mn14-Mn14/spd-Papacon](Mn14-Mn14/spd-Papacon/README.md) | 252 | -5.3747 | 0.070 / 0.002 |

The Bethe-lattice Fermi energy is 9.976 eV. The book set is for paramagnetic bcc Mn with one spin;
ANT uses it for both spins, so the electrodes are not magnetic. The geometries are fcc-like pyramids,
so the bcc parameters are only an approximation for these electrodes.
