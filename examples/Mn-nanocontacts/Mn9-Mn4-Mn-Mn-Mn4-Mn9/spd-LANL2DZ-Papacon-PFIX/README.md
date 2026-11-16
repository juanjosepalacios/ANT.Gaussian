Mn9-Mn4-Mn-Mn-Mn4-Mn9 (28 atoms, multiplicity 141), UBPBE, with LANL2DZ on atoms 10-19 and the spd
basis with ECP on the others, `Papacon` electrodes. The job name is `Mn9-Mn4-Mn-Mn-Mn4-Mn9.lan`.

Route: `Scf=(maxcycle=200,conver=4,dsymm,nodamp)`. `.ant`: as
[../spd-Papacon](../spd-Papacon/README.md) without `SL`, plus `PFIX` with
`'../spd-Papacon/P.Mn9-Mn4-Mn-Mn-Mn4-Mn9.dat'` and `NFIX = 10` (atoms 10-19).

Run [../spd-Papacon](../spd-Papacon/README.md) first: it writes the density matrix that `PFIX`
reads. The path in the `.ant` was `P.Mn9-Mn4-Mn-Mn-Mn4-Mn9.dat` when both jobs shared a directory.

## Outputs

`Mn9-Mn4-Mn-Mn-Mn4-Mn9.lan.log`, `T.Mn9-Mn4-Mn-Mn-Mn4-Mn9.lan.dat` (spin up and down).

## Results

ANT.Gaussian 3.0.0 (commit 7709e49 sources), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads, 13 November 2026, from scratch.

| Quantity | Value |
|---|---|
| Basis functions | 382 |
| SCF energy (a.u.) | -1321.088073 |
| SCF cycles | 28 |
| Device Fermi energy (eV) | -5.5398 |
| T(E_F) up / down (e^2/h) | 0.061 / 0.046 |
| Unsuccessful density integrations | 7 |
| Wall time | 401 s |
