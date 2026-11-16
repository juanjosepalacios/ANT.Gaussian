Mn nanocontact Mn9-Mn4-Mn-Mn-Mn4-Mn9 (28 atoms, multiplicity 141) between Bethe-lattice electrodes,
UBPBE.

Each subdirectory is one variation of the calculation for this structure, with its own inputs,
outputs and README. The files keep the name of the system, `Mn9-Mn4-Mn-Mn-Mn4-Mn9`.

- [spd-Papacon/](spd-Papacon/README.md): spd basis and ECP on all atoms.
- [spd-LANL2DZ-Papacon-PFIX/](spd-LANL2DZ-Papacon-PFIX/README.md): LANL2DZ on atoms 10-19, which
  start from the density of spd-Papacon through `PFIX` (job `Mn9-Mn4-Mn-Mn-Mn4-Mn9.lan`). Run
  spd-Papacon first.

## Comparison

| Variation | E (a.u.) | E_F (eV) | T(E_F) up / down (e^2/h) |
|---|---|---|---|
| [spd-Papacon/](spd-Papacon/README.md) | -438.722578 | -5.5320 | 0.069 / 0.001 |
| [spd-LANL2DZ-Papacon-PFIX/](spd-LANL2DZ-Papacon-PFIX/README.md) | -1321.088073 | -5.5398 | 0.061 / 0.046 |
