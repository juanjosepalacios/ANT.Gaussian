Au4-Au-Au-Au4 (10 atoms) with one s function (and ECP) on the 8 layer atoms and LANL2DZ on the two
apex atoms, between the 1-orbital `FCC-s` Bethe-lattice electrodes (`BLDAT/fits/Au-FCC-s`).

`.ant`: `BLPAR1 = BLPAR2 = FCC-s`, `OVERLAP = -1.0`. The route has `gfinput` and
`Scf=(dsymm,conver=6,nodamp)`; `%NProcShared=1`.

This input was in `examples/Au-chains`. It has not been run: the stored `Au4-Au-Au-Au4.log` (July
2022) stops at once with `No executable for file ../bin/l101.exe`. The gapped structure with the same
basis and electrodes, [../../Au4-Au-Au--Au4/s-FCC-s](../../Au4-Au-Au--Au4/s-FCC-s/README.md), has been
run with 3.0.0.
