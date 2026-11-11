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

## Warning: Gaussian may reorient the contact

Unless the route line has `nosymm`, Gaussian moves the molecule to its standard orientation, and ANT
works in that frame, not in the input one. `nosymm` inside `scf=()` only switches off symmetry in the
SCF and does not prevent it. This matters for magnetic systems: the field of `ZM` (`Bx`, `By`, `Bz`)
and `ZMATOM`, the angles of `SPINROTATOM` and, with SOC, the direction of the moments relative to the
contact all refer to Gaussian's frame. The Bethe-lattice d on-site energies (`edd` on x^2-y^2 and
3z^2-r^2, `edt` on xy, xz and yz) are also set in that fixed frame, so the cubic axes of the
electrodes must lie along X, Y and Z.

- Ni4-Ni-Ni-Ni4/LANL1MB-Papacon and Ni9-Ni4-Ni-Ni-Ni4-Ni9/doll-Papacon-sd are reoriented:
  X = y_in, Y = x_in, Z = -z_in. The transport axis stays on Z, but reversed, so the Ni4 field
  `Bz = 20.0` points along -z_in.
- The LANL1MB-Papacon-F, LANL1MB-Papacon-AF and CRENBS-spd-F Ni9 cases give `nosymm` on the route
  line and keep the input frame.

In all cases the cubic axes stay along X, Y and Z, so the electrodes are correct. To keep the input
frame, give `nosymm` on the route line.

Each structure README compares its variations. The `.gjf` files substitute links from
`../../../../bin`.
