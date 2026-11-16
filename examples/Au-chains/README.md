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
  - [LANL2DZ-1DLEAD/](Au8-chain/LANL2DZ-1DLEAD/README.md)
  - [s-1DLEAD-NATPC1-NPC2/](Au8-chain/s-1DLEAD-NATPC1-NPC2/README.md)
  - [s-1DLEAD-NATPC1-NPC3/](Au8-chain/s-1DLEAD-NATPC1-NPC3/README.md)
  - [s-1DLEAD-NATPC2-NPC1/](Au8-chain/s-1DLEAD-NATPC2-NPC1/README.md)
- [Au9-chain/](Au9-chain/README.md): 9 atoms
  - [spd-1DLEAD/](Au9-chain/spd-1DLEAD/README.md)
- [Au10-chain/](Au10-chain/README.md): 10 atoms
  - [LANL2DZ-1DLEAD/](Au10-chain/LANL2DZ-1DLEAD/README.md)

The `.gjf` files substitute links from `../../../../bin`; change this and `%NProcShared` to your
installation. ANT takes the job name from the title line of the `.gjf`, which matches the file
names.

## Results with ANT.Gaussian 3.0.0

Au8-chain, 16 November 2026: the minimal spd basis, LANL2DZ and one s orbital per atom (`SHF` ECP):

| System / variation | E_F (eV) | T(E_F) (e^2/h) |
|---|---|---|
| [Au8-chain/spd-1DLEAD](Au8-chain/spd-1DLEAD/README.md) | -6.2101 | 1.971 |
| [Au8-chain/LANL2DZ-1DLEAD](Au8-chain/LANL2DZ-1DLEAD/README.md) | -6.6671 | 1.967 |
| [Au8-chain/s-1DLEAD-NATPC1-NPC2](Au8-chain/s-1DLEAD-NATPC1-NPC2/README.md) | -2.8570 | 1.980 |
| [Au8-chain/s-1DLEAD-NATPC1-NPC3](Au8-chain/s-1DLEAD-NATPC1-NPC3/README.md) | -2.7779 | 1.981 |
| [Au8-chain/s-1DLEAD-NATPC2-NPC1](Au8-chain/s-1DLEAD-NATPC2-NPC1/README.md) | -2.9269 | 0.143 |

A 2-atom lead unit (`NATPC = 2`) opens a gap at E_F in a one-orbital chain; use a 1-atom unit.

## Stored outputs

The other examples have not been rerun with ANT.Gaussian 3.0.0. Their outputs are from older
versions:

| System / variation | Version | Termination | T(E_F) (e^2/h) |
|---|---|---|---|
| [Au5-chain/spd-1DLEAD](Au5-chain/spd-1DLEAD/README.md) | 2.5.2 | normal | 0.786 |
| [Au6-chain/spd-1DLEAD](Au6-chain/spd-1DLEAD/README.md) | 2.5.2 | normal | 2.000 |
| [Au7-chain/spd-1DLEAD](Au7-chain/spd-1DLEAD/README.md) | 2.5.2 | error |  |
| [Au9-chain/spd-1DLEAD](Au9-chain/spd-1DLEAD/README.md) | 2.7.0 | normal | 1.992 |
| [Au10-chain/LANL2DZ-1DLEAD](Au10-chain/LANL2DZ-1DLEAD/README.md) | 2.5.2 | error |  |
