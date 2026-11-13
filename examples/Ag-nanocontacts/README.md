Ag nanocontacts between Bethe-lattice electrodes, LSDA, with the book `Papacon` set of
`BLDAT/BL047.dat` (the default electrode set).

## Jobs

| Job | Structure | Settings |
|---|---|---|
| `Ag9-Ag4-Ag-Ag-Ag4-Ag9` | 28 atoms | `scf=(conver=5,nodamp,dsymm)`; `ALPHA = 0.05`, `SELFACC = 1.0d-5`, `PACC = 1.0d-6`, `CHARGEACC` and `FERMIACC` 1.0d-4, `EW1 = -0.1`, `EW2 = 0.1` |
| `Ag16-Ag9-Ag4-Ag-Ag-Ag4-Ag9-Ag16` | 60 atoms, apex atoms 8 A apart (tunnel junction) | `scf=(conver=5,nodamp) symmetry=none`; `ALPHA = 0.02`, `SOC` with `SOC_FAC_P/D/F = 100`, `ETA = 0.030`, `PACC = 1.0d-6`, `CHARGEACC` and `FERMIACC` 1.0d-4, `BIASVOLT = 1.0` |
| `Ag8-chain` | Ag chain | 1D-lead electrodes; not rerun |

`Ag9-Ag4-Ag-Ag-Ag4-Ag9` starts from `P.Ag9-Ag4-Ag-Ag-Ag4-Ag9.dat` if present. The `.gjf` files
substitute links from `../../bin`; change this and `%NProcShared` to your installation. ANT takes the
job name from the title line of the `.gjf`, which matches the file names.

## Results

ANT.Gaussian 3.0.0, Gaussian 09 D.01, NVIDIA HPC SDK build, 13 November 2026, from scratch.

| | Ag9-Ag4-Ag-Ag-Ag4-Ag9 | Ag16-Ag9-Ag4-Ag-Ag-Ag4-Ag9-Ag16 |
|---|---|---|
| Threads | 4 | 8 |
| Basis functions | 252 | 540 |
| SCF energy (a.u.) | -1089.835247 | -2335.954314 |
| SCF cycles | 12 | 17 |
| Bethe-lattice Fermi energy (eV) | 6.527 | 6.560 |
| Device Fermi energy (eV) | -7.2933 | -6.9195 |
| T(E_F) (e^2/h) | 1.169 | 0.00026 |
| Unsuccessful density integrations | 4 | 0 |
| Wall time | 16 s | 501 s |

Outputs: `<job>.log`, `<job>.xyz`, and the transmission in `T.Ag9-Ag4-Ag-Ag-Ag4-Ag9.dat`
(E - E_F = -0.1 to 0.1 eV, from `EW1` and `EW2`) and `T.Ag16-Ag9-Ag4-Ag-Ag-Ag4-Ag9-Ag16.SOC.dat`
(-3 to 3 eV, with SOC).

## Notes

- Ag9-Ag4-Ag-Ag-Ag4-Ag9 changes little from its earlier ANT 2.5.2 output (E = -1089.835708 a.u.,
  T(E_F) = 1.1687), although the earlier `Papacon` set had `edd` and `edt` swapped.
- `Ag16-Ag9-Ag4-Ag-Ag-Ag4-Ag9-Ag16` was called `Ag25-...`, but its geometry and title are
  Ag16-Ag9-Ag4-Ag-Ag-Ag4-Ag9-Ag16. With the old file names ANT looked for an `.ant` named after the
  title, did not find it and ran without SOC and bias; the files are renamed. SOC enters the
  transport calculation (`SCFSOC = F`).
- `Ag9-Ag4-Ag-Ag-Ag4-Ag9.job`, `ERRORS`, `OUTPUT` and the copies of `l101.exe` and `l502.exe` in this
  directory are left from an old cluster run and are not used.
