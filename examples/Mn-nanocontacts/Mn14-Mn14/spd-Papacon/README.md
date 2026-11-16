Mn14-Mn14 (28 atoms, multiplicity 141), UBPBE with `guess=mix`, spd basis and ECP on all atoms,
`Papacon` electrodes.

Route: `Scf=(maxcycle=200,conver=6,dsymm)`. `.ant`: `ALPHA = 0.03`, `FERMIACC` and `CHARGEACC`
1.0d-5, `PACC = 1.0d-8`, `MAX = 4`, `NSPINLOCK = 4`, `LDOS_BEG = 14`, `LDOS_END = 15`.

## Outputs

`Mn14-Mn14.log`, `T.Mn14-Mn14.dat` (spin up and down), `Mn14-Mn14.xyz`.

## Results

ANT.Gaussian 3.0.0 (commit 7709e49 sources), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads, 13 November 2026, from scratch.

| Quantity | Value |
|---|---|
| Basis functions | 252 |
| SCF energy (a.u.) | -438.562168 |
| SCF cycles | 65 |
| Device Fermi energy (eV) | -5.3747 |
| T(E_F) up / down (e^2/h) | 0.070 / 0.002 |
| Unsuccessful density integrations | 157 |
| Wall time | 240 s |

The 157 unsuccessful density integrations (point limit) and the small transmissions call for care
with these numbers.
