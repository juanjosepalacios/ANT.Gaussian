Benzene between Au tips, `Papacon` electrodes, LSDA, with a geometry optimization
(`opt(maxcycle=100,loose)`, `nosymm`, `Scf=(conver=7,...)`).

`.ant`: `ALPHA = 0.015`, `SL = 0.01`, accuracies 1.0d-5/1.0d-7, `NPULAY = 4`, `LDOS_BEG = 14`,
`LDOS_END = 15`, `MULLIKEN`, `NSPINLOCK = 0`, `NCHANNELS = 6`, `BLPAR1 = BLPAR2 = Papacon`.

Basis: the contact Au atoms 1-6 and 27-32 have the CRENBS basis and ECP (9 orbitals, matching
`Papacon`), the inner Au atoms 7-10 and 23-26 LANL2DZ. `IOp(1/7=10000)` sets the RMS force threshold
to 0.01 au.

This input was in `examples/benzene`. It was rebuilt for 3.0.0 (see `CHANGELOG.md`); its new run has
not been stored yet. The stored `input_Electrodo_Au_3_benzene.log` is the old ANT 2.5.0 run, which
stopped at once with `NAOrbs read from BL data base not equal to NAOrbs of electrode plane` (1-orbital
`lsda-6s` electrodes against a full spd Au basis). `Q.input_Electrodo_Au_3_benzene.dat` and the
`.xyz` also come from an older run.
