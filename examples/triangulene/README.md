[4]Triangulene (C33H15, a triangular graphene flake with zigzag edges and spin S = 3/2 when
isolated and neutral) between the Al and Au pyramids of
[../benzene/Al-benzene-Al-axis](../benzene/Al-benzene-Al-axis/README.md), BLYP, spin-unrestricted.

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation, with its own inputs, outputs and README. The files keep the name of
the system.

- [Al-triangulene-Al/](Al-triangulene-Al/README.md): [4]triangulene with its plane containing the
  transport axis, 88 atoms, STO-3G on C and H
  - [sp-Papacon-UHF/](Al-triangulene-Al/sp-Papacon-UHF/README.md)
  - [sp-Papacon-UPLUS-UHF/](Al-triangulene-Al/sp-Papacon-UPLUS-UHF/README.md)
- [Au-triangulene-Au/](Au-triangulene-Au/README.md): the same structure with Au pyramids, LANL2DZ
  - [spd-Papacon-UHF/](Au-triangulene-Au/spd-Papacon-UHF/README.md)
  - [spd-Papacon-UPLUS-UHF/](Au-triangulene-Au/spd-Papacon-UPLUS-UHF/README.md)

The `.gjf` files substitute links from `../../../../bin`; change this, `%NProcShared` and `%mem` to
your installation. ANT takes the job name from the title line of the `.gjf`, which matches the file
names.
