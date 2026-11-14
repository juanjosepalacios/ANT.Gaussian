Ag nanocontacts between Bethe-lattice electrodes, LSDA, with a minimal spd basis and ECP
(36-electron core) on all atoms. Ag atomic chains between 1D leads are in
[../Ag-chain](../Ag-chain/README.md).

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation (basis, electrodes, SOC, bias), with its own inputs, outputs and
README. The files keep the name of the system.

- [Ag9-Ag4-Ag-Ag-Ag4-Ag9/](Ag9-Ag4-Ag-Ag-Ag4-Ag9/README.md): 28 atoms, 9-atom and 4-atom layers on
  each side of two apex atoms 3 A apart
  - [spd-Papacon/](Ag9-Ag4-Ag-Ag-Ag4-Ag9/spd-Papacon/README.md)
- [Ag16-Ag9-Ag4-Ag-Ag-Ag4-Ag9-Ag16/](Ag16-Ag9-Ag4-Ag-Ag-Ag4-Ag9-Ag16/README.md): 60 atoms, apex
  atoms 8 A apart (tunnel junction)
  - [spd-Papacon-SOC-bias/](Ag16-Ag9-Ag4-Ag-Ag-Ag4-Ag9-Ag16/spd-Papacon-SOC-bias/README.md)
- [Ag16-Ag9-Ag4-Ag-Ag-Ag4-Ag9-Ag16-nogap/](Ag16-Ag9-Ag4-Ag-Ag-Ag4-Ag9-Ag16-nogap/README.md): the
  same 60 atoms with the gap closed, apex atoms 3 A apart
  - [spd-Papacon-SOC-bias/](Ag16-Ag9-Ag4-Ag-Ag-Ag4-Ag9-Ag16-nogap/spd-Papacon-SOC-bias/README.md)

`Papacon` is the book set of `BLDAT/BL047.dat` (the default Ag electrode set), from
Papaconstantopoulos, Handbook of the Band Structure of Elemental Solids, 2nd ed. (2015).

The `.gjf` files substitute links from `../../../../bin`; change this and `%NProcShared` to your
installation. ANT takes the job name from the title line of the `.gjf`, which matches the file
names.

## Results

| System / variation | Version | E_F (eV) | T(E_F) (e^2/h) |
|---|---|---|---|
| [Ag9-Ag4-Ag-Ag-Ag4-Ag9/spd-Papacon](Ag9-Ag4-Ag-Ag-Ag4-Ag9/spd-Papacon/README.md) | 3.0.0 | -7.2933 | 1.169 |
| [Ag16-Ag9-Ag4-Ag-Ag-Ag4-Ag9-Ag16/spd-Papacon-SOC-bias](Ag16-Ag9-Ag4-Ag-Ag-Ag4-Ag9-Ag16/spd-Papacon-SOC-bias/README.md) | 3.0.0 | -6.9195 | 0.00026 |
| [Ag16-Ag9-Ag4-Ag-Ag-Ag4-Ag9-Ag16-nogap/spd-Papacon-SOC-bias](Ag16-Ag9-Ag4-Ag-Ag-Ag4-Ag9-Ag16-nogap/spd-Papacon-SOC-bias/README.md) | 3.0.0 | -7.1610 | 1.295 |

The device Fermi energies are on the energy scale of the electrodes.
