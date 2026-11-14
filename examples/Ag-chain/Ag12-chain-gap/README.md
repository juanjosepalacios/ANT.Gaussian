Ag chain Ag12-chain-gap (12 atoms on the z axis) between one-dimensional leads, LSDA.

Two 6-atom segments with 3 A spacing, z = 0 to 15 A and z = 21 to 36 A, separated by a 6 A gap
between atoms 6 and 7. Each side holds two copies of the 3-atom repeated unit of the 1D leads.
[Ag12-chain](../Ag12-chain/README.md) is the same chain without the gap. The charge is 0
and the multiplicity 1.

Each subdirectory is one variation of the calculation for this chain, with its own inputs, outputs
and README. The files keep the name of the system, `Ag12-chain-gap`.

- [spd-1DLEAD/](spd-1DLEAD/README.md): spd basis and ECP on all atoms, one-dimensional leads with a
  3-atom repeated unit (`NATPC1 = NATPC2 = 3`, `NPC = 1`).

## Comparison

| Variation | Basis and electrodes | E_F (eV) | T(E_F) (e^2/h) |
|---|---|---|---|
| [spd-1DLEAD/](spd-1DLEAD/README.md) | spd, `1DLEAD`, 3-atom unit | -6.0107 | 0.0649 |
