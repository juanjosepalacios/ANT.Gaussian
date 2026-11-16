Graphene examples.

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation (basis, electrodes, settings), with its own inputs, outputs and
README. The files keep the name of the system.

- [nanographene/](nanographene/README.md): 170-atom nanographene flake with a carbon Bethe-lattice
  electrode around its edge
  - [Papacon/](nanographene/Papacon/README.md): `SL = 0.1` and the local DOS of a central atom
  - [Papacon-SL0.01/](nanographene/Papacon-SL0.01/README.md): the original input, `SL = 1.0d-2`
- [nanographene2/](nanographene2/README.md): a smaller flake with ghost electrodes
  - [sp-ghost/](nanographene2/sp-ghost/README.md)
- [graphene-ribbon/](graphene-ribbon/README.md): a 304-atom graphene ribbon
  - [sp-graphene/](graphene-ribbon/sp-graphene/README.md)
- [hydrophene-ribbon/](hydrophene-ribbon/README.md): the ribbon with H atoms in place of C
  - [s-Papacon/](hydrophene-ribbon/s-Papacon/README.md)

Only `nanographene/Papacon` has been run with ANT.Gaussian 3.0.0; `graphene-ribbon` has a stored
output from ANT 2.3.4, and the others are inputs only.

The `.gjf` files substitute links from `../../../../bin`; change this and `%NProcShared` to your
installation. ANT takes the job name from the title line of the `.gjf`, which matches the file
names.
