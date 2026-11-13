Al nanocontacts between Bethe-lattice electrodes, LSDA (spin-restricted, non-magnetic).

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation (basis, electrodes, bias), with its own inputs, outputs and README.
The files keep the name of the system.

- [Al9-Al4-Al-Al-Al4-Al9/](Al9-Al4-Al-Al-Al4-Al9/README.md): 28 atoms, 9-atom and 4-atom layers on
  each side of two apex atoms
  - [small-Papacon/](Al9-Al4-Al-Al-Al4-Al9/small-Papacon/README.md)
  - [dz-Papacon/](Al9-Al4-Al-Al-Al4-Al9/dz-Papacon/README.md)
  - [dzp-Papacon/](Al9-Al4-Al-Al-Al4-Al9/dzp-Papacon/README.md)
  - [dzp-Papacon-bias/](Al9-Al4-Al-Al-Al4-Al9/dzp-Papacon-bias/README.md)
  - [small-Papacon-nonortho/](Al9-Al4-Al-Al-Al4-Al9/small-Papacon-nonortho/README.md)

The electrodes use the Al sets of `BLDAT/BL013.dat` from Papaconstantopoulos, Handbook of the Band
Structure of Elemental Solids, 2nd ed. (2015), Table 7.10 (s and p first-neighbour terms of the
2-center fit): `Papacon`, the orthogonal fit and the `BLPAR` default, and `Papacon-nonortho`, the
non-orthogonal fit with overlaps.

## Warning: Gaussian reorients the contact

None of the cases has `nosymm` on the route line (`nosymm` inside `scf=()` does not prevent it), so
Gaussian moves the contact to its standard orientation, X = -y_in, Y = -x_in, Z = -z_in, and ANT
works in that frame. The Al sets have only s and p orbitals and the contact is non-magnetic, so the
results do not depend on it.

The structure README compares the variations. The `.gjf` files substitute links from
`../../../../bin`.
