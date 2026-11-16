Au atomic chains between one-dimensional leads (`TYPE1 = TYPE2 = 1DLEAD`), LSDA. Au nanocontacts
between Bethe-lattice electrodes are in [../Au-nanocontacts](../Au-nanocontacts/README.md), and the H
chains that were here are in [../H-chain](../H-chain/README.md).

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation (basis, electrodes, settings), with its own inputs, outputs and
README. The files keep the name of the system.

- [Au5-chain/](Au5-chain/README.md): 5 atoms
  - [spd-1DLEAD/](Au5-chain/spd-1DLEAD/README.md)
- [Au6-chain/](Au6-chain/README.md): 6 atoms
  - [spd-1DLEAD/](Au6-chain/spd-1DLEAD/README.md)
- [Au7-chain/](Au7-chain/README.md): 7 atoms
  - [spd-1DLEAD/](Au7-chain/spd-1DLEAD/README.md)
- [Au8-chain/](Au8-chain/README.md): 8 atoms
  - [spd-1DLEAD/](Au8-chain/spd-1DLEAD/README.md)
- [Au9-chain/](Au9-chain/README.md): 9 atoms
  - [spd-1DLEAD/](Au9-chain/spd-1DLEAD/README.md)
- [Au10-chain/](Au10-chain/README.md): 10 atoms
  - [LANL2DZ-1DLEAD/](Au10-chain/LANL2DZ-1DLEAD/README.md)
- [Au-Au-Au-Au/](Au-Au-Au-Au/README.md): 4 atoms
  - [s-1DLEAD/](Au-Au-Au-Au/s-1DLEAD/README.md)
- [Au-Au-Au-Au-Au-Au/](Au-Au-Au-Au-Au-Au/README.md): 6 atoms
  - [s-1DLEAD/](Au-Au-Au-Au-Au-Au/s-1DLEAD/README.md)

The `.gjf` files substitute links from `../../../../bin`; change this and `%NProcShared` to your
installation. ANT takes the job name from the title line of the `.gjf`, which matches the file
names.

## Stored outputs

None of these examples has been rerun with ANT.Gaussian 3.0.0. The outputs are from older versions:

| System / variation | Version | Termination | T(E_F) (e^2/h) |
|---|---|---|---|
| [Au5-chain/spd-1DLEAD](Au5-chain/spd-1DLEAD/README.md) | 2.5.2 | normal | 0.786 |
| [Au6-chain/spd-1DLEAD](Au6-chain/spd-1DLEAD/README.md) | 2.5.2 | normal | 2.000 |
| [Au7-chain/spd-1DLEAD](Au7-chain/spd-1DLEAD/README.md) | 2.5.2 | error |  |
| [Au8-chain/spd-1DLEAD](Au8-chain/spd-1DLEAD/README.md) | 2.5.2 | normal | 2.000 |
| [Au9-chain/spd-1DLEAD](Au9-chain/spd-1DLEAD/README.md) | 2.7.0 | normal | 1.992 |
| [Au10-chain/LANL2DZ-1DLEAD](Au10-chain/LANL2DZ-1DLEAD/README.md) | 2.5.2 | error |  |
| [Au-Au-Au-Au/s-1DLEAD](Au-Au-Au-Au/s-1DLEAD/README.md) | 2.5.2 | normal | 1.373 |
| [Au-Au-Au-Au-Au-Au/s-1DLEAD](Au-Au-Au-Au-Au-Au/s-1DLEAD/README.md) | 2.5.2 | normal | 1.994 |
