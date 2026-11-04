Four cases for an Al nanocontact (Al9-Al4-Al-Al-Al4-Al9, 28 atoms, Bethe-lattice electrodes):

1) Minimal sp basis with the ECP10CE pseudopotential on all atoms (`.gbssmall`)
2) LANL2DZ basis and pseudopotential on the central atoms 10-19 (`.gbsdz`)
3) LANL2DZ plus one d shell on atoms 10-19 (`.gbsdzp`) and tighter accuracy settings (`.ant-dzp`)
4) Case 3 at a bias voltage of 1 V (`.ant-dzp-bias`), restarted from the zero-bias density matrix

The outer planes (atoms 1-9 and 20-28) connect to the Bethe lattice and always keep the minimal
basis and ECP10CE, which match the Al Bethe-lattice parameters (BLDAT/BL013.dat). Only the central
atoms change basis, as in the Pt example.

## Running a case

The `.gjf` includes the basis with `@Al9-Al4-Al-Al-Al4-Al9.gbssmall`. For another case, change that
line to `.gbsdz` or `.gbsdzp`. For case 3, also copy `Al9-Al4-Al-Al-Al4-Al9.ant-dzp` to
`Al9-Al4-Al-Al-Al4-Al9.ant`.

| Case | Basis file | `.ant` |
|---|---|---|
| small | `.gbssmall` | `.ant` |
| dz | `.gbsdz` | `.ant` |
| dzp | `.gbsdzp` | `.ant-dzp` (`CHARGEACC`/`FERMIACC` 1.0d-5, `FULLACC`, `DEVBRENT`) |
| dzp-bias | `.gbsdzp` | `.ant-dzp-bias` (`.ant-dzp` plus `BIASVOLT = 1.0`) |

For case 4, first run case 3 and keep its `P.Al9-Al4-Al-Al-Al4-Al9.dat` in the directory: ANT reads
it at start-up and connects the leads from the first cycle. Also replace `dsymm` with `nosymm` in the
`Scf=(...)` options of the `.gjf`, as required for finite bias.

Common settings (`.ant`): `PACC = 1.0d-6`, `NPULAY = 5`, `SL = 0.001`, LDOS on atom 14, energy
window `EW1 = -5` to `EW2 = 5` eV. Default `ETA` and Bethe-lattice parameters.

Outputs of each case carry its suffix: `.log-<case>`, `T.Al9-Al4-Al-Al-Al4-Al9.dat-<case>`
(transmission) and `DOS.Al9-Al4-Al-Al-Al4-Al9.dat-<case>` (energy, total DOS, LDOS of atom 14,
then its orbital-resolved LDOS).

## Results

ANT.Gaussian 2.8.4 (commit 5703257), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads,
3 November 2026.

| | small | dz | dzp | dzp-bias (1 V) |
|---|---|---|---|---|
| Basis functions | 112 | 152 | 202 | 202 |
| SCF energy (a.u.) | -55.302369 | -55.699845 | -55.874762 | -56.420194 |
| SCF cycles (device Fermi searches) | 24 (20) | 20 (16) | 21 (16) | 18 (20) |
| Device Fermi energy (eV) | -6.1245 | -5.9234 | -5.8927 | -5.8768 |
| T(E_F) | 1.919 | 1.760 | 1.782 | 1.855 |
| T(E_F - 1 eV) / T(E_F + 1 eV) | 0.980 / 1.953 | 1.214 / 1.934 | 1.107 / 1.845 | 0.718 / 1.784 |
| Maximum T in the window | 2.46 at +5.00 eV | 3.38 at +3.91 eV | 3.70 at +3.58 eV | 3.42 at +3.70 eV |
| LDOS of atom 14 at E_F (states/eV) | 0.321 | 0.361 | 0.392 | 0.522 |
| Total device DOS at E_F (states/eV) | 7.06 | 7.50 | 7.76 | 8.44 |
| Current at 1 V (uA) | | | 136.6 (from T at 0 V) | 139.0 |
| Wall time | 10 s | 17 s | 99 s | 190 s |

All runs terminate normally and pass the spectral-weight check on the first attempt. The first three
report no integration warnings. The bias run reports one unsuccessful non-equilibrium (real-axis)
integration at 127 points in ANT cycle 7 of 18; the later cycles integrate within the limit.

## Notes

- The conductance is about 1.8-1.9 G0 in all cases, consistent with a few open channels through a
  single-atom Al contact. Moving from the minimal basis to LANL2DZ lowers T(E_F) by about 8%; adding
  d functions changes it by about 1%.
- Above E_F the larger bases bring transmission maxima closer to the Fermi level (3.6-3.9 eV
  instead of the window edge), so the minimal basis is less reliable away from E_F.
- The LDOS on the central atom grows with basis size (0.32 to 0.39 states/eV at E_F), as the
  added functions give the contact atom more flexibility.
- Case 3 changes both the basis and the numerical settings, so its cost increase (99 s) combines the
  extra d functions with the tighter tolerances, FULLACC and Brent root search. Its results differ
  from case 2 mainly through the d functions; the tighter settings are a convergence check.
- SCF energies of different basis sets are not comparable with each other. At finite bias the
  energy comes from a non-equilibrium density and is not comparable with the zero-bias energies.
- The current is I = G0 * integral of T(E) over the bias window, E_F - 0.5 eV to E_F + 0.5 eV, with
  G0 = 77.48 uS; ANT does not print it. Using the self-consistent T(E) at 1 V gives 139.0 uA, an
  average transmission of 1.794 in the window. The zero-bias T(E) of case 3 gives 136.6 uA, so the
  I-V curve is still nearly linear at 1 V (about 1.8 G0).
- At 1 V the transmission at E_F rises from 1.78 to 1.86, while the transmission 1 eV below E_F
  drops from 1.11 to 0.72. The LDOS on the central atom at E_F increases from 0.39 to 0.52 states/eV
  as the bias redistributes charge in the contact.
- The average on-site energies of the atoms, which show the potential drop, are written to
  V.job_name.dat only when the HAMILTON keyword is present; this case does not set it.
