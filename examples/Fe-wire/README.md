Fe nanowire (Fe64, bcc, 64 atoms), spin-polarized LSDA.

The wire is a bcc Fe column along z with a = 2.87 A: 13 square layers of 4 atoms at
(+-1.435, +-1.435) alternate with 12 single atoms on the axis, from z = -17.22 to +17.22 A. The charge
is 0 and the starting multiplicity is 163.

Each subdirectory is one variation of the calculation for this wire, with its own inputs, outputs and
README:

- [BCC-sd/](BCC-sd/README.md): CRENBS s+d basis and ECP, Bethe-lattice electrodes `BCC-sd`. Its README
  describes the two-step procedure (QC seed, then ANT) that the other variations also use.
- [BCC-spd/](BCC-spd/README.md): the same basis plus a scaled STO-3G 4p shell, Bethe-lattice
  electrodes `BCC-spd`.

## Comparison

| Basis and electrodes | E_F (eV) | T(E_F) up + down (e^2/h) | G (G0) |
|---|---|---|---|
| sd, Bethe `BCC-sd` ([BCC-sd/](BCC-sd/README.md)) | -2.114 | 2.801 + 1.658 = 4.458 | 2.23 |
| sd, 1D (`1DLEAD`, `NATPC1/2 = 5`, `NPC = 3`), restarted from a converged Bethe-lattice density (not included) | -2.181 | 2.679 + 1.572 = 4.251 | 2.13 |
| spd, Bethe `BCC-spd` ([BCC-spd/](BCC-spd/README.md)) | -5.583 | 2.893 + 2.440 = 5.333 | 2.67 |
