Examples with ghost electrodes (`TYPE = GHOST`): one side of the system is a Bethe lattice and the
other has no electrode. Al with an sp basis and ECP, STO-3G on C and H.

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation (basis, electrodes, settings), with its own inputs, outputs and
README. The files keep the name of the system.

- [Al9-Al4-C/](Al9-Al4-C/README.md): a C atom on an Al9-Al4 pyramid, 15 atoms
  - [sp-ghost-opt/](Al9-Al4-C/sp-ghost-opt/README.md)
- [Al9-Al4-H2/](Al9-Al4-H2/README.md): an H2 molecule on an Al9-Al4 pyramid
  - [sp-ghost/](Al9-Al4-H2/sp-ghost/README.md)

The `.gjf` files substitute links from `../../../../bin`; change this and `%NProcShared` to your
installation. ANT takes the job name from the title line of the `.gjf`, which matches the file
names.
