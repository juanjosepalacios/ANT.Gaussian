Ni nanocontacts between Bethe-lattice electrodes, spin-polarized LSDA.

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation (basis, electrodes, magnetic configuration), with its own inputs,
outputs and README. The files keep the name of the system.

- [Ni4-Ni-Ni-Ni4/](Ni4-Ni-Ni-Ni4/README.md): 10 atoms, two 4-atom layers joined by two apex atoms
  - [LANL1MB-Papacon/](Ni4-Ni-Ni-Ni4/LANL1MB-Papacon/README.md)
- [Ni9-Ni4-Ni-Ni-Ni4-Ni9/](Ni9-Ni4-Ni-Ni-Ni4-Ni9/README.md): 28 atoms, 9-atom and 4-atom layers on
  each side of two apex atoms
  - [LANL1MB-Papacon-F/](Ni9-Ni4-Ni-Ni-Ni4-Ni9/LANL1MB-Papacon-F/README.md)
  - [LANL1MB-Papacon-AF/](Ni9-Ni4-Ni-Ni-Ni4-Ni9/LANL1MB-Papacon-AF/README.md)
  - [CRENBS-spd-F/](Ni9-Ni4-Ni-Ni-Ni4-Ni9/CRENBS-spd-F/README.md)
  - [doll-Papacon-sd/](Ni9-Ni4-Ni-Ni-Ni4-Ni9/doll-Papacon-sd/README.md)

`Papacon` is the orthogonal first-neighbour Slater-Koster fit of ferromagnetic fcc Ni from
Papaconstantopoulos, Handbook of the Band Structure of Elemental Solids, 2nd ed. (2015). The
CRENBS-spd-F results were computed with an earlier non-orthogonal set of that name and cannot be
reproduced with the current `BL028.dat`.

Each structure README compares its variations. The `.gjf` files substitute links from
`../../../../bin`.
