Fe nanowire (Fe64, bcc, 64 atoms) between Bethe-lattice electrodes with the spd basis, spin-polarized
LSDA.

This is the [BCC-sd](../BCC-sd/README.md) calculation with a 4p shell added to the Fe basis and
`BCC-spd` electrodes. The geometry, charge (0), starting multiplicity (163), two-step procedure and
all other settings are the same.

## Differences from the sd example

- Basis: each Fe has the CRENBS s+d basis and ECP plus the STO-3G 4p shell scaled by 0.586 (the
  `P 3` shell in both `.gjf` files), 576 basis functions instead of 384. It is the basis of
  `BLDAT/fits/Fe-BCC-spd/Fe-crenbs-spd.gbs`.
- `.ant`: `BLPAR1 = BLPAR2 = BCC-spd` (s+p+d, fitted to the CRYSTAL23 LSDA DOS of bulk bcc Fe with
  this basis; see `BLDAT/fits/Fe-BCC-spd/README.md`).
- Step 2 substitutes links from `../../../../bin`.

The running instructions of the BCC-sd README apply unchanged: step 1 must end normally, and
`P.Fe64-bcc-wire.dat` must be deleted before repeating step 2 from the checkpoint.

## Outputs

- `Fe64-bcc-wire.qc.log`: step 1
- `Fe64-bcc-wire.log`: step 2
- `T.Fe64-bcc-wire.dat`: transmission for E - E_F = -3 to +3 eV in 0.01 eV steps, in units of e^2/h,
  spin up and then spin down, separated by a blank line
- `Fe64-bcc-wire.xyz`: the geometry written by ANT; the 18 `Ar` atoms mark the Bethe-lattice
  attachment sites

## Results

ANT.Gaussian 2.8.4 (commit 027b6f8), Gaussian 09 D.01, NVIDIA HPC SDK build, 8 threads,
7-8 November 2026.

| | Step 1 (QC) | Step 2 (ANT) |
|---|---|---|
| Basis functions | 576 | 576 |
| SCF energy (a.u.) | -1426.735171 | -1428.134667 |
| SCF cycles (device Fermi searches) | 23 | 116 (107) |
| Device Fermi energy (eV) | | -5.5834 |
| T(E_F) up / down (e^2/h) | | 2.893 / 2.440 |
| T(E_F) total (e^2/h) | | 5.333 |
| Conductance G (G0) | | 2.67 |
| Spin polarization of T(E_F), (up - down)/(up + down) | | 0.08 |
| T total at E_F - 1 eV / E_F + 1 eV | | 2.564 / 3.497 |
| Electrons in left / right half | | 255.86 / 256.14 |
| Spin of left / right half | | 105.49 / 105.24 |
| Unsuccessful density integrations | | 0 |
| Time | 1 h 13 min 50 s CPU (9 min 47 s wall) | 6 h 38 min 28 s wall |

Both steps terminate normally. Charges and spins are from ANT's Mulliken population analysis.

## Notes

- Step 2 is not tightly converged. It stops at cycle 116 with "Density matrix is not changing but
  DIIS error= 3.65D-02": RMSDP is 9.9e-6, just under `conver=5`, but MaxDP is 3.7e-4 and the last
  cycle moves the device Fermi energy by 0.14 eV (from -5.447 to -5.583 eV). With `ALPHA = 0.01`
  a small RMSDP does not by itself show self-consistency.
- A restart of step 2 from this run's `P.Fe64-bcc-wire.dat` (not included) stopped after 2 cycles
  for the same reason, with E_F = -5.5845 eV, T(E_F) = 2.898 / 2.488 and G = 2.69 G0. Its T(E)
  differs from `T.Fe64-bcc-wire.dat` by at most 0.03 (up) and 0.08 (down) within E_F +- 0.5 eV.
- The device Fermi energy ranges from -4.2 to +6.6 eV in the first 45 device Fermi searches and
  stays between -5.0 and -5.6 eV from the 70th search on. The Bethe-lattice electrode Fermi level is 6.8517 eV.
- Moments: 3.66-3.67 on the outermost 4-atom layers, 3.35-3.44 on the inner layers and 2.67-2.91 on
  the axis atoms. The two halves differ by 0.14 electrons and 0.25 in spin. The run ends with
  361.36 alpha and 150.64 beta electrons, a total spin of 210.7.
- Spin-down T is zero between about E - E_F = +1.27 and +1.40 eV and between +2.81 and +2.88 eV;
  spin-up T is zero between -2.94 and -2.92 eV.
- Gaussian prints "Spurious integrated density or basis function" 58 times (grid integration error
  just above its 1e-3 relative tolerance); the sd example has none.
- The device Fermi searches (`DEVBRENT`) use 1646 charge evaluations, about 22260 s of the 23908 s
  of step 2.
- Comparison with the sd basis for the same wire:

  | Basis and electrodes | E_F (eV) | T(E_F) up + down (e^2/h) | G (G0) |
  |---|---|---|---|
  | sd, Bethe `BCC-sd` ([BCC-sd](../BCC-sd/README.md)) | -2.114 | 2.801 + 1.658 = 4.458 | 2.23 |
  | spd, Bethe `BCC-spd` (this example) | -5.583 | 2.893 + 2.440 = 5.333 | 2.67 |

  Almost all of the difference is in spin down. With spd, the minority T is suppressed from
  E - E_F = -3 to -1.5 eV (mean 0.31 against 1.30) and higher within E_F +- 0.5 eV (mean 1.92 against
  1.25), with a peak of about 3.2 near +0.2 eV. Spin up changes by 3% at E_F.
