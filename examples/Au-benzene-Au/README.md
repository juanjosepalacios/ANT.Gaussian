Benzene between Au electrodes, LSDA.

- `input_Electrodo_Au_3_benzene`: the benzene junction. It uses `lsda-6s` Au electrodes and the
  keywords `MOLGAP`, `MOLBLOCKS`, `IP` and `EA`, which ANT 3.0.0 does not have: the run stops at once
  with `ERROR - Undefined keyword: MOLGAP`. The stored outputs (`.log`, `T.*.dat`, `.MOA.dat`, the
  `.jpg`) come from ANT 2.5.2 and were not rerun.
- `Au9-Au4-Au-Au4-Au9`: an Au contact without benzene and without `.ant` (ANT defaults, default
  electrode set `Papacon` of `BLDAT/BL079.dat`). Rerun with ANT.Gaussian 3.0.0 (commit 7709e49
  sources), 8 threads, 13 November 2026: 256 basis functions, E = -1041.523349 a.u. after 37 cycles,
  device E_F = -8.0632 eV, T(E_F) = 2.333 e^2/h, normal termination. Outputs:
  `Au9-Au4-Au-Au4-Au9.log` and `T.Au9-Au4-Au-Au4-Au9.dat`.
