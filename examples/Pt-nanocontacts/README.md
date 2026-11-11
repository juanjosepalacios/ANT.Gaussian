Pt nanocontacts between Bethe-lattice electrodes, LSDA.

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation (basis, spin-orbit coupling, magnetic configuration), with its own
inputs, outputs and README. The files keep the name of the system.

- [Pt9-Pt4-Pt-Pt-Pt4-Pt9/](Pt9-Pt4-Pt-Pt-Pt4-Pt9/README.md): 28 atoms, 9-atom and 4-atom layers on
  each side of two apex atoms
  - Spin-orbit coupling:
    - [SOC-small/](Pt9-Pt4-Pt-Pt-Pt4-Pt9/SOC-small/README.md)
    - [SOC-large/](Pt9-Pt4-Pt-Pt-Pt4-Pt9/SOC-large/README.md)
  - Spin rotation test on a contact forced to be magnetic:
    - [mag/](Pt9-Pt4-Pt-Pt-Pt4-Pt9/mag/README.md)
    - [mag-rotz/](Pt9-Pt4-Pt-Pt-Pt4-Pt9/mag-rotz/README.md)
    - [mag-rotx/](Pt9-Pt4-Pt-Pt-Pt4-Pt9/mag-rotx/README.md)

All cases use the default Pt Bethe-lattice set `Papacon` in `BLDAT/BL078.dat`.

Each structure README compares its variations. The `.gjf` files substitute links from
`../../../../bin`.
