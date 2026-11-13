Al nanocontact Al9-Al4-Al-Al-Al4-Al9 (28 atoms) between Bethe-lattice electrodes, LSDA.

Two 9-atom outer layers at z = -6 and +6 A, two 4-atom layers at z = -4 and +4 A and two apex atoms
on the axis at z = -2 and +2 A (fcc stacking with a = 4.0 A in the layers). The charge is 0 and the
calculation is spin-restricted.

The outer layers (atoms 1-9 and 20-28) connect to the Bethe lattice and always keep the minimal sp
basis and ECP10CE, which match the sp Bethe-lattice sets. Only the central atoms 10-19 change basis.
Each subdirectory is one variation, with its own inputs, outputs and README. The files keep the name
of the system, `Al9-Al4-Al-Al-Al4-Al9`.

- [small-Papacon/](small-Papacon/README.md): minimal sp basis and ECP10CE on all atoms (`.gbssmall`)
- [dz-Papacon/](dz-Papacon/README.md): LANL2DZ basis and ECP on atoms 10-19 (`.gbsdz`)
- [dzp-Papacon/](dzp-Papacon/README.md): LANL2DZ plus one d shell on atoms 10-19 (`.gbsdzp`), tighter
  accuracies, `FULLACC` and `DEVBRENT`. Includes the density matrix that starts the bias case.
- [dzp-Papacon-bias/](dzp-Papacon-bias/README.md): dzp-Papacon at a bias of 1 V, restarted from its
  zero-bias density matrix
- [small-Papacon-nonortho/](small-Papacon-nonortho/README.md): small-Papacon with the non-orthogonal
  set `Papacon-nonortho` at `OVERLAP = 0.5`. Its SCF does not reach a physical solution (see below).

All but the last use the default orthogonal `Papacon` set of `BLDAT/BL013.dat`.

Common `.ant` settings: `PACC = 1.0d-6`, `NPULAY = 5`, `SL = 0.001`, LDOS on atom 14, energy window
`EW1 = -5` to `EW2 = 5` eV; `CHARGEACC` and `FERMIACC` 1.0d-4 (1.0d-5 for dzp and dzp-bias).

## Results

ANT.Gaussian commit 4152255, Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads, 13 November 2026.

| | small | dz | dzp | dzp-bias (1 V) | small-nonortho |
|---|---|---|---|---|---|
| Basis functions | 112 | 152 | 202 | 202 | 112 |
| SCF energy (a.u.) | -55.141462 | -55.552249 | -55.771503 | -56.404968 | -53.153962 |
| SCF cycles (device Fermi searches) | 24 (20) | 20 (16) | 23 (18) | 18 (20) | 36 (32) |
| Bethe-lattice Fermi energy (eV) | 9.937 | 9.937 | 9.937 | 9.937 | 6.904 |
| Device Fermi energy (eV) | -6.5513 | -6.3619 | -6.3123 | -6.3188 | 13.4979 |
| T(E_F) | 1.953 | 1.885 | 1.870 | 1.931 | 1.923 |
| T(E_F - 1 eV) / T(E_F + 1 eV) | 1.031 / 1.973 | 1.057 / 1.788 | 1.114 / 1.677 | 0.948 / 1.626 | 1.947 / 1.628 |
| Maximum T in the window | 2.00 at +0.74 eV | 2.74 at +4.00 eV | 3.06 at +3.78 eV | 2.94 at +3.94 eV | 3.97 at +4.95 eV |
| LDOS of atom 14 at E_F (states/eV) | 0.308 | 0.387 | 0.410 | 0.521 | 0.446 |
| Total device DOS at E_F (states/eV) | 7.07 | 7.53 | 7.71 | 8.02 | 9.82 |
| Current at 1 V (uA) | | | 147.6 (from T at 0 V) | 145.7 | |
| Wall time | 11 s | 18 s | 48 s | 193 s | 9 s |

All runs terminate normally, pass the spectral-weight check without widening the energy window and
report no unsuccessful density integrations.

## Notes

- With the orthogonal book set the conductance is 1.87-1.95 G0 in all cases, a few open channels
  through a single-atom Al contact. LANL2DZ lowers T(E_F) by about 3% from the minimal basis; the d
  shell changes it by 1%.
- Above E_F the larger bases bring transmission maxima to 3.8-4.0 eV, where the minimal basis has a
  plateau near 2 up to the window edge, so the minimal basis is less reliable away from E_F.
- The LDOS on the central atom grows with the basis (0.31 to 0.41 states/eV at E_F).
- SCF energies of different basis sets are not comparable with each other. At finite bias the energy
  comes from a non-equilibrium density and is not comparable with the zero-bias energies.
- The current is I = G0 * integral of T(E) over the bias window, E_F - 0.5 eV to E_F + 0.5 eV, with
  G0 = 77.48 uS; ANT does not print it. The self-consistent T(E) at 1 V gives 145.7 uA and the
  zero-bias T(E) of dzp 147.6 uA, so the I-V curve is still nearly linear at 1 V (about 1.9 G0).
- small-Papacon-nonortho does not converge to a physical state. The device Fermi energy climbs from
  -5.7 eV to about +21 eV and settles at 13.5 eV, 20 eV above the orthogonal case, while the energy
  rises from -56.97 to -53.15 a.u. Gaussian stops on `Density matrix is not changing` (DIIS error
  3e-4). Its T(E_F) agrees with the orthogonal case by coincidence; its T(E) drops to 0 and jumps to
  1.45 at about -1.1 eV. With `OVERLAP = 1.0` ANT stays more than 7 minutes in the Bethe-lattice
  energy-window search and was stopped.
- These results replace those of 3 November 2026, computed with an earlier `Papacon` set of unknown
  origin (es 0.0884, ep 0.4741 Ry, with overlaps unused at the default `OVERLAP`): T(E_F) was 1.919
  (small), 1.760 (dz), 1.782 (dzp) and 1.855 (dzp-bias).
