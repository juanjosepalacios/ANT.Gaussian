Ag chain Ag12-chain (12 atoms on the z axis, no gap) between one-dimensional leads, LSDA.

Twelve atoms with 3 A spacing from z = 0 to 33 A. It is
[Ag12-chain-gap](../Ag12-chain-gap/README.md) with the 6 A gap closed, as the reference for a
perfect chain. The charge is 0 and the multiplicity 1.

Each subdirectory is one variation of the calculation for this chain, with its own inputs, outputs
and README. The files keep the name of the system, `Ag12-chain`.

- [spd-1DLEAD/](spd-1DLEAD/README.md): spd basis and ECP on all atoms, one-dimensional leads with a
  3-atom repeated unit (`NATPC1 = NATPC2 = 3`, `NPC = 1`).

## Comparison

| Variation | Basis and electrodes | E_F (eV) | T(E_F) (e^2/h) |
|---|---|---|---|
| [spd-1DLEAD/](spd-1DLEAD/README.md) | spd, `1DLEAD`, 3-atom unit | -5.1326 | 2.000 |
