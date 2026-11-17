Hydrogen atomic chains between one-dimensional leads (`TYPE1 = TYPE2 = 1DLEAD`), LSDA/STO-3G. They
were in `examples/Au-chains`.

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation (basis, electrodes, settings), with its own inputs, outputs and
README. The files keep the name of the system.

- [H10-chain/](H10-chain/README.md): 10 atoms
  - [STO-3G-1DLEAD/](H10-chain/STO-3G-1DLEAD/README.md)

The `.gjf` files substitute links from `../../../../bin`; change this and `%NProcShared` to your
installation. ANT takes the job name from the title line of the `.gjf`, which matches the file
names.

## Stored outputs

| System / variation | Version | T(E_F) (e^2/h) |
|---|---|---|
| [H10-chain/STO-3G-1DLEAD](H10-chain/STO-3G-1DLEAD/README.md) | 2.5.2 | 2.000 |
