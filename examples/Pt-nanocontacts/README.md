Pt nanocontact (Pt9-Pt4-Pt-Pt-Pt4-Pt9, 28 atoms, Bethe-lattice electrodes) in five cases.

Each subdirectory is one case, with its own inputs, outputs and README. Every case is run as
`Pt9-Pt4-Pt-Pt-Pt4-Pt9.gjf` with `Pt9-Pt4-Pt-Pt-Pt4-Pt9.ant` in its own directory; the `.gjf` files
substitute links from `../../../bin`.

Spin-orbit coupling:

1) [SOC-small/](SOC-small/README.md): a Pt nanocontact with SOC and a small basis set (`.gbssmall`)
2) [SOC-large/](SOC-large/README.md): a Pt nanocontact with a large basis set and modified values of
   SOC integrals (`.gbslarge`)

Spin rotation test on a Pt contact forced to be magnetic (no SOC, small basis):

3) [mag/](mag/README.md): collinear spin-polarized calculation with multiplicity 15 and spin locking
4) [mag-rotz/](mag-rotz/README.md): case 3 restarted with every atom listed in `SPINROTATOM` at
   theta = 0
5) [mag-rotx/](mag-rotx/README.md): case 3 restarted with every spin rotated to the x axis

The outputs of cases 1 and 2 come from earlier versions of ANT.

## Spin rotation test (cases 3-5)

Pt is not magnetic, but a fixed multiplicity and spin locking keep it spin-polarized. The Pt
Bethe-lattice parameters (`Papacon` in `BLDAT/BL078.dat`) are spin-independent, so the electrodes do
not change when the spins rotate. Rotating every spin of the device by the same angle is then a rigid
rotation of the whole system, and the Fermi energy, charges, moment sizes and transmission must not
change. This tests the implementation of `SPINROTATOM`.

| Case | Settings |
|---|---|
| 3 (`mag`) | Charge 0, multiplicity 15 (147 alpha and 133 beta electrons), `nosymm`, `NSPINLOCK = 1000` |
| 4 (`mag-rotz`) | Case 3 plus `SPINROTATOM` with atoms 1-28 at theta = 0, phi = 0 |
| 5 (`mag-rotx`) | Case 3 plus `SPINROTATOM` with atoms 1-28 at theta = 90, phi = 0 |

The three cases share the same `.gjf` and basis. Run case 3 first. For cases 4 and 5, copy its density
matrix `mag/P.Pt9-Pt4-Pt-Pt-Pt4-Pt9.dat` into `mag-rotz/` and `mag-rotx/`, so both restart from the
same collinear solution. Delete `mag/P.Pt9-Pt4-Pt-Pt-Pt4-Pt9.dat` before repeating case 3 from
scratch; ANT reads that file at start-up as the initial density matrix.

Case 4 does not rotate anything, but listing atoms in `SPINROTATOM` sends the calculation through the
same spinor code as case 5. Comparing cases 4 and 5 therefore isolates the rotation. Cases 4 and 5
write `T.Pt9-Pt4-Pt-Pt-Pt4-Pt9.SOC.dat` with the total transmission in units of e^2/h. Case 3
writes `T.Pt9-Pt4-Pt-Pt-Pt4-Pt9.dat` with a spin-up block and a spin-down block. All transmission
files cover E - E_F = -3 to +3 eV in 0.01 eV steps.

### Results

ANT.Gaussian 2.8.4 (commit 658a5e0), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads, 4 November
2026. Cases 4 and 5 ran at the same time.

| | 3: collinear | 4: all along z | 5: all along x |
|---|---|---|---|
| SCF energy (a.u.) | -770.696875 | -770.695765 | -770.695765 |
| SCF cycles (device Fermi searches) | 46 (68) | 2 (6) | 2 (6) |
| Device Fermi energy in the SCF (eV) | -8.5989 (up), -9.2530 (down) | locked, as case 3 | locked, as case 3 |
| Device Fermi energy for transport (eV) | -9.088631 | -9.090146 | -9.090146 |
| T(E_F) (e^2/h) | 0.571 (up) + 2.483 (down) = 3.053 | 3.0487 | 3.0487 |
| T at E_F - 1 eV / E_F + 1 eV | 5.090 / 1.939 | 5.0986 / 1.9395 | 5.0986 / 1.9395 |
| Electrons in left / right half | 139.977 / 140.023 | 139.98595 / 140.01405 | 139.98595 / 140.01405 |
| Spin of left / right half | 3.054 / 3.101 | 3.04875 / 3.09836 | 3.04875 / 3.09836 |
| Spin of atoms 14 / 15 | 0.2798 / 0.2795 (along z) | 0.27863 / 0.27790 (along z) | 0.27863 / 0.27790 (along x) |
| Wall time | 504 s | 105 s | 128 s |

All runs terminate normally with no unsuccessful integrations. Charges and spins come from ANT's
Mulliken population analysis.

### Notes

- Cases 4 and 5 agree exactly. Their transmission files are byte-identical over the whole window,
  and the Fermi energy, charges and spin of each half agree to 1e-12. Every atomic moment moves from z
  to x with the same size, leaving no y or z component. The rotation is therefore a proper rigid
  rotation of the spin space.
- The SCF energies of cases 4 and 5 are equal because the rotation is applied after
  self-consistency and does not enter the SCF energy. The test of the rotation is the equality of
  the Fermi energy, populations and transmission, which are recomputed after the rotation.
- Spin locking acts only during the SCF. There the two spin channels have separate Fermi energies
  and the 14 unpaired electrons set by the multiplicity (Gaussian's Mulliken spin sum is 13.84). After
  convergence ANT finds a single device Fermi energy for the transport and population analysis, in
  all three cases, so electrons move between the spin channels: ANT's population analysis gives a
  total spin of 6.15, about 0.22-0.28 per atom. Spin-down electrons dominate the transmission at E_F
  (2.48 against 0.57).
- Case 3 differs slightly from cases 4 and 5 (E_F -9.0886 against -9.0901 eV, T(E_F) 3.053 against
  3.049) because the restart of cases 4 and 5 stops after 2 SCF cycles at an energy 1.1 mHa above case
  3, so their collinear density is not exactly that of case 3. This is unrelated to the rotation.
- `SPINROTATOM` is not self-consistent; see `doc/keywords.txt`. A true spin flip is better done with
  `SPINEDIT`.
