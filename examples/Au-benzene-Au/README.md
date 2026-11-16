Benzene between Au electrodes, LSDA. The Au contact without benzene is now
[../Au-nanocontacts/Au9-Au4-Au-Au4-Au9/spd-Papacon-noant](../Au-nanocontacts/Au9-Au4-Au-Au4-Au9/spd-Papacon-noant/README.md).

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation (basis, electrodes, settings), with its own inputs, outputs and
README. The files keep the name of the system.

- [input_Electrodo_Au_3_benzene/](input_Electrodo_Au_3_benzene/README.md): benzene between two Au
  tips, 32 atoms
  - [spd-Papacon/](input_Electrodo_Au_3_benzene/spd-Papacon/README.md)
  - [spd-Papacon-opt/](input_Electrodo_Au_3_benzene/spd-Papacon-opt/README.md)

The stored outputs come from ANT 2.5.x. spd-Papacon does not run with ANT.Gaussian 3.0.0
(removed keywords); spd-Papacon-opt has been rebuilt for 3.0.0, but its new run is not stored yet.

The `.gjf` files substitute links from `../../../../bin`; change this and `%NProcShared` to your
installation. ANT takes the job name from the title line of the `.gjf`, which matches the file
names.
