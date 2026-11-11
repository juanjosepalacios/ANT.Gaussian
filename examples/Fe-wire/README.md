Fe nanowires, spin-polarized LSDA.

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation (basis, electrodes), with its own inputs, outputs and README. The files
keep the name of the system.

- [Fe64-bcc-wire/](Fe64-bcc-wire/README.md): 64 atoms, a bcc Fe column along z
  - [BCC-sd/](Fe64-bcc-wire/BCC-sd/README.md)
  - [BCC-spd/](Fe64-bcc-wire/BCC-spd/README.md)
  - [1DLEAD-spd/](Fe64-bcc-wire/1DLEAD-spd/README.md)

## Warning: keep `nosymm` on the route line

Unless the route line has `nosymm`, Gaussian moves the molecule to its standard orientation, and ANT
works in that frame, not in the input one. `nosymm` inside `scf=()` only switches off symmetry in the
SCF and does not prevent it. This matters for magnetic systems: the field of `ZM` (`Bx`, `By`, `Bz`)
and `ZMATOM`, the angles of `SPINROTATOM` and, with SOC, the direction of the moments relative to the
wire all refer to Gaussian's frame. The Bethe-lattice d on-site energies (`edd` on x^2-y^2 and
3z^2-r^2, `edt` on xy, xz and yz) are also set in that fixed frame, so the cubic axes of the
electrodes must lie along X, Y and Z.

All Fe64-bcc-wire inputs, the QC seeds included, give `nosymm` on the route line and keep the input
frame: the wire runs along z and the bcc cubic axes lie along x, y and z. Keep it there in new
inputs.

Each structure README compares its variations. The `.gjf` files substitute links from
`../../../../bin`.
