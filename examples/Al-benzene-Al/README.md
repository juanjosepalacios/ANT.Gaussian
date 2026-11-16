Benzene between Al electrodes, BLYP. These inputs were in `examples/benzene`; the HTML manual
(`doc/manual_2_g03.html`) describes them.

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation (basis, electrodes, settings), with its own inputs, outputs and
README. The files keep the name of the system.

- [Al-benzene-Al/](Al-benzene-Al/README.md): benzene between two Al pyramids, 52 atoms
  - [sp-Papacon/](Al-benzene-Al/sp-Papacon/README.md)
  - [sp-Papacon-NEMBED10/](Al-benzene-Al/sp-Papacon-NEMBED10/README.md)

The `.gjf` files substitute links from `../../../../bin`; change this and `%NProcShared` to your
installation. ANT takes the job name from the title line of the `.gjf`, which matches the file
names.

`Al-benzene-Al-bethe.gjf` keeps the placeholder `/path/to/ANT.Gaussian/bin` of the manual.
