Benzene between Al and Au electrodes. The inputs of `Al-benzene-Al` and
`input_Electrodo_Au_3_benzene` were in `examples/benzene`; the HTML manual (`doc/manual_2_g03.html`)
describes `Al-benzene-Al`. The Au contact without benzene is
[../Au-nanocontacts/Au9-Au4-Au-Au4-Au9/spd-Papacon-noant](../Au-nanocontacts/Au9-Au4-Au-Au4-Au9/spd-Papacon-noant/README.md).

Each subdirectory is one atomic structure, named after it; inside it, each subdirectory is one
variation of the calculation (basis, electrodes, settings), with its own inputs, outputs and
README. The files keep the name of the system.

- [Al-benzene-Al/](Al-benzene-Al/README.md): benzene flat between two Al pyramids, 52 atoms, BLYP
  - [sp-Papacon/](Al-benzene-Al/sp-Papacon/README.md)
- [Al-benzene-Al-axis/](Al-benzene-Al-axis/README.md): the same benzene and Al pyramids, with the
  ring plane containing the transport axis and the tips moved apart to an 8.69 A apex gap
  - [sp-Papacon/](Al-benzene-Al-axis/sp-Papacon/README.md)
  - [sp-Papacon-UPLUS/](Al-benzene-Al-axis/sp-Papacon-UPLUS/README.md)
  - [sp-Papacon-UHF/](Al-benzene-Al-axis/sp-Papacon-UHF/README.md)
  - [sp-Papacon-UPLUS-UHF/](Al-benzene-Al-axis/sp-Papacon-UPLUS-UHF/README.md)
- [Au-benzene-Au-axis/](Au-benzene-Au-axis/README.md): `Al-benzene-Al-axis` with Au pyramids,
  spin-unrestricted BLYP, LANL2DZ
  - [spd-Papacon-UHF/](Au-benzene-Au-axis/spd-Papacon-UHF/README.md)
  - [spd-Papacon-UPLUS-UHF/](Au-benzene-Au-axis/spd-Papacon-UPLUS-UHF/README.md)
- [input_Electrodo_Au_3_benzene/](input_Electrodo_Au_3_benzene/README.md): benzene between two
  small Au tips, 32 atoms, LSDA. Its stored outputs come from ANT 2.5.x; `spd-Papacon` does not run
  with ANT.Gaussian 3.0.0 or later (removed keywords), `spd-Papacon-opt` has been rebuilt for 3.0.0
  but its new run is not stored yet.
  - [spd-Papacon/](input_Electrodo_Au_3_benzene/spd-Papacon/README.md)
  - [spd-Papacon-opt/](input_Electrodo_Au_3_benzene/spd-Papacon-opt/README.md)

The `.gjf` files substitute links from `../../../../bin`; change this, `%NProcShared` and `%mem` to
your installation. ANT takes the job name from the title line of the `.gjf`, which matches the file
names.
