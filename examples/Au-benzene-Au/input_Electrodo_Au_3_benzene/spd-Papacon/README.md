Benzene between Au tips, `Papacon` electrodes (no `BLPAR`, so the default set), LSDA.

`.ant`: `ALPHA = 0.015`, `NPULAY = 4`, `SL = 1.0d-2`, accuracies 1.0d-5/1.0d-7, `EW1 = -5`,
`EW2 = 5`, and `MOLGAP`, `MOLBLOCKS`, `IP = 0.1`, `EA = 0.7`.

ANT 3.0.0 does not have `MOLGAP`, `MOLBLOCKS`, `IP` and `EA`: the run stops at once with
`ERROR - Undefined keyword: MOLGAP`. The stored outputs (`.log`, `T.*.dat`, `Q.*.dat`, `.MOA.dat`,
the `.jpg` and `.ps` plots) come from ANT 2.5.2 and were not rerun.
