Al nanocontact Al16-Al9-Al4-Al-Al-Al4-Al9-Al16 (60 atoms) between Bethe-lattice electrodes, LSDA.

The [Al9-Al4-Al-Al-Al4-Al9](../Al9-Al4-Al-Al-Al4-Al9/README.md) contact with one more layer of 16
atoms on each side, in the same fcc (001) stacking (a = 4.0 A, layers 2 A apart): 16-atom layers at
z = -8 and +8 A (4 atoms at 2 A from the axis, 8 at 4.47 A and 4 at 6 A), 9-atom layers at -6 and
+6 A, 4-atom layers at -4 and +4 A and two apex atoms on the axis at -2 and +2 A (atoms 30 and 31).
The charge is 0 and the calculation is spin-restricted.

The Bethe lattices attach to the outer 16-atom layers (atoms 1-16 and 45-60), which keep the minimal
sp basis and ECP10CE of the sp Bethe-lattice sets. The 28 inner atoms 17-44 have LANL2DZ.

- [dz-Papacon/](dz-Papacon/README.md): LANL2DZ inside, orthogonal Bethe-lattice set `Papacon`,
  `conver=6`

## Results

ANT.Gaussian commit 29d1a76, Gaussian 09 D.01, NVIDIA HPC SDK build, 8 threads, 13 November 2026.

| | dz-Papacon | Al9-Al4-Al-Al-Al4-Al9 dz-Papacon |
|---|---|---|
| Atoms (Bethe-lattice atoms) | 60 (32) | 28 (18) |
| Basis functions | 352 | 152 |
| SCF energy (a.u.) | -120.441810 | -55.552249 |
| SCF cycles (device Fermi searches) | 32 (20) | 20 (16) |
| Device Fermi energy (eV) | -6.3020 | -6.3619 |
| T(E_F) | 2.007 | 1.885 |
| T(E_F - 1 eV) / T(E_F + 1 eV) | 1.258 / 1.250 | 1.057 / 1.788 |
| Maximum T in the window | 3.70 at +3.90 eV | 2.74 at +4.00 eV |
| LDOS of apex atom at E_F (states/eV) | 0.292 | 0.387 |
| Total device DOS at E_F (states/eV) | 15.96 | 7.53 |

The run terminates normally (Gaussian stops on `Density matrix is not changing`, which with
`conver=6` means an RMS density change below 1e-6), passes the spectral-weight check without
widening the energy window and reports no unsuccessful density integrations.

## Timing with the number of cores

The same input on 1, 2, 4 and 8 cores (`%NProcShared` and `OMP_NUM_THREADS` set to the number of
cores), one job at a time on an Intel Core i9-13900KF, of which 8 cores were available.

| Cores | Wall time (s) | Speedup | Efficiency | CPU time (s) |
|---|---|---|---|---|
| 1 | 569.7 | 1.00 | 100% | 569 |
| 2 | 319.8 | 1.78 | 89% | 624 |
| 4 | 189.2 | 3.01 | 75% | 712 |
| 8 | 135.1 | 4.21 | 53% | 976 |

- Nearly all the time is spent in link 502 (the SCF with ANT). A second 8-core run took 135.2 s.
- The times follow Amdahl's law with about 73 s of serial and 497 s of parallel work (a serial
  fraction of about 13%), so the speedup cannot exceed about 7.8 and would be about 6.1 on 24 cores
  (about 94 s). For several jobs at once, 2-4 cores per job use the machine better than 8.
- All core counts give 32 SCF cycles, E_F = -6.302 eV (within 4e-5 eV) and T(E_F) = 2.007. The total
  energy differs by up to 5e-5 a.u. (-120.44178 to -120.44183) because of the order of the parallel
  sums.

## Notes

- With the extra layer T(E_F) rises from 1.885 to 2.007 and the transmission becomes nearly
  symmetric about E_F within 1 eV, while the LDOS on the apex atom drops from 0.39 to 0.29 states/eV.
- The route has no `nosymm`, so Gaussian reorients the contact (X = -y_in, Y = -x_in, Z = -z_in), as
  in the 28-atom cases. The sets are sp and the contact is non-magnetic, so the results do not
  depend on it.
