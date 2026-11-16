Benzene between two Al pyramids (52 atoms), BLYP, with an sp basis and ECP on Al and STO-3G on C and
H, between `Papacon` Bethe-lattice electrodes (`BLDAT/BL026.dat`).

Each subdirectory is one variation of the calculation for this structure, with its own inputs,
outputs and README. The files keep the job name of their `.gjf` title, `Al-benzene-Al` or
`Al-benzene-Al-bethe`; the geometry is the same.

- [sp-Papacon/](sp-Papacon/README.md): ANT defaults except `ALPHA = 0.02`.
- [sp-Papacon-NEMBED10/](sp-Papacon-NEMBED10/README.md): the electrode settings written out, with
  `NEMBED1 = NEMBED2 = 10` and `NATOMEL1 = NATOMEL2 = 20`.
