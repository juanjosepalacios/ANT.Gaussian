Ag atomic chains between one-dimensional leads, LSDA, with a minimal spd basis and ECP (36-electron
core, 11 valence electrons) on all atoms.

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation (basis, electrodes), with its own inputs, outputs and README. The files
keep the name of the system.

- [Ag12-chain-gap/](Ag12-chain-gap/README.md): 12 atoms on the z axis with 3 A spacing and a 6 A gap
  between atoms 6 and 7
  - [spd-1DLEAD/](Ag12-chain-gap/spd-1DLEAD/README.md)
- [Ag12-chain/](Ag12-chain/README.md): the same chain without the gap, 12 atoms with
  uniform 3 A spacing
  - [spd-1DLEAD/](Ag12-chain/spd-1DLEAD/README.md)

Both use periodic 1D leads built from the cluster (`TYPE1 = TYPE2 = 1DLEAD`) with a 3-atom repeated
unit (`NATPC1 = NATPC2 = 3`, `NPC = 1`), and run in two steps: a plain-Gaussian QC SCF to
`conver=3` that writes the checkpoint, then the ANT run from it (`guess=read`). Both steps use
`nosymm`, so Gaussian keeps the input orientation and atom order, which the 1D leads rely on. The
`.gjf` files substitute links from `../../../../bin`; change this and `%NProcShared` to your
installation.

## Results

ANT.Gaussian 3.0.0, Gaussian 09 D.01, NVIDIA HPC SDK build, 8 threads, 14 November 2026.

| System / variation | E_F (eV) | T(E_F) (e^2/h) | G (G0) |
|---|---|---|---|
| [Ag12-chain-gap/spd-1DLEAD](Ag12-chain-gap/spd-1DLEAD/README.md) | -6.0107 | 0.0649 | 0.032 |
| [Ag12-chain/spd-1DLEAD](Ag12-chain/spd-1DLEAD/README.md) | -5.1326 | 2.000 | 1.000 |

The perfect chain conducts through one open spin-degenerate channel; the 6 A gap reduces T(E_F) by
a factor of about 30. Both transmissions vanish below E - E_F = -1.05 eV and in two narrow windows
near -0.66 and +0.86 eV. These come from the extracted 1D leads, which are only approximately
periodic with two 3-atom units per side, not from the gap.
