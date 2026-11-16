Mn9-Mn4-Mn-Mn-Mn4-Mn9 (28 atoms, multiplicity 141), UBPBE with `guess=mix`, spd basis and ECP on all
atoms, `Papacon` electrodes.

Route: `Scf=(maxcycle=200,conver=4,dsymm)`. `.ant`: `ALPHA = 0.020`, `NPULAY = 4`, `FERMIACC` and
`CHARGEACC` 1.0d-5, `PACC = 1.0d-8`, `SL = 2.0d-2`, `MAX = 5`. The `.gjf` uses `%NProcShared=12`.

The run writes `P.Mn9-Mn4-Mn-Mn-Mn4-Mn9.dat` (3.5 MB, not stored), which
[../spd-LANL2DZ-Papacon-PFIX](../spd-LANL2DZ-Papacon-PFIX/README.md) reads.

## Outputs

`Mn9-Mn4-Mn-Mn-Mn4-Mn9.log`, `T.Mn9-Mn4-Mn-Mn-Mn4-Mn9.dat` (spin up and down, E - E_F from -3 to
3 eV), `Mn9-Mn4-Mn-Mn-Mn4-Mn9.xyz`.

## Results

ANT.Gaussian 3.0.0 (commit 7709e49 sources), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads, 13 November 2026, from scratch.

| Quantity | Value |
|---|---|
| Basis functions | 252 |
| SCF energy (a.u.) | -438.722578 |
| SCF cycles | 22 |
| Device Fermi energy (eV) | -5.5320 |
| T(E_F) up / down (e^2/h) | 0.069 / 0.001 |
| Unsuccessful density integrations | 9 |
| Wall time | 152 s |
