C60 between the Al and Au pyramids of [../benzene/Al-benzene-Al](../benzene/Al-benzene-Al/README.md),
BLYP, spin-unrestricted, LANL2DZ basis with a CRENBS basis on the electrode planes that touch the
Bethe lattice.

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation, with its own inputs, outputs and README. The files keep the name of
the system.

- [Al-C60-Al/](Al-C60-Al/README.md): C60 with a hexagon facing each Al apex, 100 atoms
  - [sp-Papacon-UHF/](Al-C60-Al/sp-Papacon-UHF/README.md)
  - [sp-Papacon-UPLUS-UHF/](Al-C60-Al/sp-Papacon-UPLUS-UHF/README.md)
- [Au-C60-Au/](Au-C60-Au/README.md): the same structure with Au pyramids
  - [spd-Papacon-UHF/](Au-C60-Au/spd-Papacon-UHF/README.md)
  - [spd-Papacon-UPLUS-UHF/](Au-C60-Au/spd-Papacon-UPLUS-UHF/README.md)

The `.gjf` files substitute links from `../../../../bin`; change this, `%NProcShared` and `%mem` to
your installation. ANT takes the job name from the title line of the `.gjf`, which matches the file
names.
