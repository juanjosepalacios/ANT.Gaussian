Pt nanocontact Pt9-Pt4-Pt-Pt-Pt4-Pt9 (28 atoms) between Bethe-lattice electrodes, LSDA.

Two 9-atom outer layers at z = 9.46 and -0.84 A, two 4-atom layers at z = 7.50 and 1.12 A and two
apex atoms near the axis at z = 5.54 and 3.08 A (fcc (001) stacking with a = 3.92 A, layers 1.96 A
apart, apex atoms 2.46 A apart). The upper apex atom is displaced 0.1 A off the axis along y. The
charge is 0; the starting multiplicity is 1 in the SOC cases and 15 in the magnetic ones. All cases
use the default Pt Bethe-lattice set `Papacon` in `BLDAT/BL078.dat` and an ECP basis read from a
`.gbssmall` or `.gbslarge` file next to the `.gjf`.

Each subdirectory is one variation of the calculation for this contact, with its own inputs, outputs
and README. Every case is run as `Pt9-Pt4-Pt-Pt-Pt4-Pt9.gjf` with `Pt9-Pt4-Pt-Pt-Pt4-Pt9.ant` in its
own directory.

Spin-orbit coupling:

1) [SOC-small/](SOC-small/README.md): SOC with the small basis set (`.gbssmall`)
2) [SOC-large/](SOC-large/README.md): SOC with a large basis set and modified values of the SOC
   integrals (`.gbslarge`)

Spin rotation test on the contact forced to be magnetic (no SOC, small basis):

3) [mag/](mag/README.md): collinear spin-polarized calculation with multiplicity 15 and spin locking
4) [mag-rotz/](mag-rotz/README.md): case 3 restarted with every atom listed in `SPINROTATOM` at
   theta = 0
5) [mag-rotx/](mag-rotx/README.md): case 3 restarted with every spin rotated to the x axis

All five cases were rerun on 11 November 2026 with the corrected Pt d on-site energies (`edt` on the
t2g orbitals, `edd` on eg; see `CHANGELOG.md`): ANT.Gaussian master (commit 35eefc7), Gaussian 09
D.01, NVIDIA HPC SDK build, 4 threads (`%nproc=4`), all started from scratch without a density
matrix file except cases 4 and 5.

## SOC results (cases 1 and 2)

| | 1: SOC, small basis | 2: SOC, large basis |
|---|---|---|
| Basis functions | 252 | 442 |
| SCF energy (a.u.) | -770.749547 | -1691.297233 |
| SCF cycles (device Fermi searches) | 46 (34) | 56 (42) |
| Device Fermi energy (eV) | -8.9861 | -8.7871 |
| T(E_F) (e^2/h) | 4.231 | 5.898 |
| T at E_F - 1 eV / E_F + 1 eV | 4.798 / 1.939 | 3.857 / 2.186 |
| Electrons in left / right half | 139.988 / 140.012 | 179.982 / 180.018 |
| Wall time | 356 s | 2021 s |

Both runs terminate normally with no unsuccessful integrations and no magnetic moment. With the
swapped d energies (runs of 10 April 2024 with an earlier ANT) T(E_F) was 4.252 and 5.914 and E_F
-8.9753 and -8.7824 eV; the large-basis run then also had 86 unsuccessful integrations.

## Warning: Gaussian reorients the contact

None of the `.gjf` files has `nosymm` on the route line: cases 1 and 2 give `scf=(dsymm)` and cases
3-5 `scf=(nosymm)`, which only switches off symmetry inside the SCF. Gaussian therefore moves the
contact to its standard orientation, and ANT works in that frame, not in the input one. In all five
logs the standard coordinates are

    X = -y_in + 0.0036 A,   Y = z_in - 4.31 A,   Z = -x_in

so the transport axis is Y, not z, and the off-axis displacement of the upper apex atom is along -X.

This is a proper rotation by whole axes, so the cubic axes of the contact still lie along X, Y and Z
and the Bethe-lattice d on-site energies (`edd` on x^2-y^2 and 3z^2-r^2, `edt` on xy, xz and yz,
both set in the fixed X, Y, Z frame) keep their meaning. A geometry that Gaussian rotates by a
general angle would break this. To keep the input frame, give `nosymm` on the route line, as in the
Ni9 examples.

The spin directions of cases 4 and 5 refer to the standard frame too. "Along z" (case 4) is along
-x_in and "along x" (case 5) along -y_in, both perpendicular to the transport axis; neither is along
the contact.

## Spin rotation test (cases 3-5)

Pt is not magnetic, but a fixed multiplicity and spin locking keep case 3 spin-polarized. The Pt
Bethe-lattice parameters (`Papacon` in `BLDAT/BL078.dat`) are spin-independent, so the electrodes do
not change when the spins rotate. Rotating every spin of the device by the same angle is then a rigid
rotation of the whole system, and the Fermi energy, charges, moment sizes and transmission must not
change. Comparing cases 4 and 5 tests the implementation of `SPINROTATOM`.

| Case | Settings |
|---|---|
| 3 (`mag`) | Charge 0, multiplicity 15 (147 alpha and 133 beta electrons), `scf=(nosymm)`, `NSPINLOCK = 1000` |
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

Run on 11 November 2026 with the setup given above. Case 3 ran alongside case 2, and cases 4 and 5
one after the other, each with 4 threads.

| | 3: collinear, locked | 4: all along z | 5: all along x |
|---|---|---|---|
| SCF energy (a.u.) | -770.690678 | -770.706875 | -770.706875 |
| SCF cycles (device Fermi searches) | 46 (68) | 30 (32) | 30 (32) |
| Device Fermi energy in the SCF (eV) | -8.6050 (up), -9.2677 (down) | single, not locked | single, not locked |
| Device Fermi energy for transport (eV) | -9.100902 | -9.153535 | -9.153534 |
| T(E_F) (e^2/h) | 0.563 (up) + 2.470 (down) = 3.033 | 4.050 | 4.050 |
| T at E_F - 1 eV / E_F + 1 eV | 5.091 / 1.946 | 4.800 / 1.943 | 4.800 / 1.943 |
| Electrons in left / right half | 139.976 / 140.024 | 140.00029 / 139.99971 | 140.00029 / 139.99971 |
| Spin of left / right half | 3.024 / 3.071 | 0.67872 / 0.68396 | 0.67872 / 0.68396 |
| Spin of atoms 14 / 15 | 0.2788 / 0.2784 (along z) | 0.13600 / 0.13528 (along z) | 0.13600 / 0.13528 (along x) |
| Wall time | 515 s | 500 s | 525 s |

All runs terminate normally with no unsuccessful integrations. Charges and spins come from ANT's
Mulliken population analysis.

### Notes

- Cases 4 and 5 agree. The charge and spin of each half agree to 1e-9, the transport Fermi energies
  to 1e-6 eV and the transmission to 1e-4 e^2/h (the last printed digit) over the whole window. The
  SCF energies differ by 2e-7 a.u. Every atomic moment moves from z to x with the same size, leaving
  no y or z component. The rotation is therefore a proper rigid rotation of the spin space.
- Since commit f027077 (5 November 2026) ANT ignores `NSPINLOCK` when it starts from a density
  matrix file ("NSPINLOCK = 1000 ignored" in the logs of cases 4 and 5). Cases 4 and 5 therefore do
  not stay on the locked solution of case 3, as they did in the runs of 4 November (2 SCF cycles): they
  relax for 30 cycles with a single Fermi energy to a weakly magnetic state 16 mHa below case 3, with
  a spin of 0.68 per half instead of 3.05. The energy still fell by 1e-4 a.u. per cycle when the SCF
  met `conver=5`, so the moments may shrink further with a tighter threshold. Cases 4 and 5 are
  therefore no longer the rotated case 3, but they remain a valid test of the rotation because both
  follow the same path from the same density matrix.
- Spin locking acts only during the SCF of case 3. There the two spin channels have separate Fermi
  energies and the 14 unpaired electrons set by the multiplicity (Gaussian's Mulliken spin sum is
  13.84). After convergence ANT finds a single device Fermi energy for the transport and population
  analysis, so electrons move between the spin channels: ANT's population analysis gives a total spin
  of 6.10, 0.19-0.28 per atom. Spin-down electrons dominate the transmission at E_F (2.47 against
  0.56).
- Compared with the runs of 4 November with the swapped d energies, case 3 changes little: E_F
  -9.1009 against -9.0886 eV, T(E_F) 3.033 against 3.053, spin per half 3.02-3.07 against 3.05-3.10.
- `SPINROTATOM` is not self-consistent; see `doc/keywords.txt`. A true spin flip is better done with
  `SPINEDIT`.
- The `.gjf` files substitute links from `../../../../bin` and use `%nproc=4` (cases 3-5) or
  `%nproc=8` (cases 1 and 2); change it to the number of cores available. The runs above all used 4.
