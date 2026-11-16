Benzene between two Au tips (32 atoms), LSDA, with `Papacon` Bethe-lattice electrodes. The Au
atoms that join the electrodes have one s, one p and one d shell (9 orbitals, with a 68-electron
ECP), which fits the 9-orbital `Papacon` set.

Each subdirectory is one variation of the calculation for this structure, with its own inputs,
outputs and README. The files keep the name of the system, `input_Electrodo_Au_3_benzene`.

- [spd-Papacon/](spd-Papacon/README.md): transport run with the molecular-gap keywords of ANT 2.5.x.
- [spd-Papacon-opt/](spd-Papacon-opt/README.md): geometry optimization with `LDOS` and Mulliken
  charges. This input was in `examples/benzene`.
