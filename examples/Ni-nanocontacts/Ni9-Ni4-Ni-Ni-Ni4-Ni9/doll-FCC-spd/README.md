Ni nanocontact (Ni9-Ni4-Ni-Ni-Ni4-Ni9, 28 atoms) between parallel Bethe-lattice electrodes,
spin-polarized LSDA, with an all-electron basis on the central atoms and `FCC-spd` electrodes.

The case of [doll-Papacon-sd](../doll-Papacon-sd/README.md) (same geometry, all-electron basis on the
central atoms 10-19) with two changes:

- Atoms 1-9 and 20-28 (the outer 9-atom layers) have the CRENBS s and d shells and ECP plus the
  p shell of `BLDAT/fits/Ni-FCC-spd/Ni-crenbs-spd.gbs` (STO-3G Ni 4p with the exponents scaled by
  0.60, optimized for the LSDA energy of bulk fcc Ni), the basis the `FCC-spd` set was fitted to.
- The Bethe-lattice electrodes are `FCC-spd` (`BLDAT/fits/Ni-FCC-spd/README.md`).

The outer layers are at z = -4.83 and +4.83 A, the 4-atom layers at -3.065 and +3.065 A and the apex
atoms at -1.30 and +1.30 A. The charge is 0 and the starting multiplicity is 37.

## Running the example

Run `Ni9-Ni4-Ni-Ni-Ni4-Ni9.gjf` with `Ni9-Ni4-Ni-Ni-Ni4-Ni9.ant`, with
`scf=(conver=5,maxcycle=300,dsymm) pop=min`.

`.ant` settings: Bethe-lattice electrodes `FCC-spd` (`BLPAR1`, `BLPAR2`), `NSPINLOCK = 10`,
`NIP = 7`, `ALPHA = 0.05`, `SL = 1.0d-2`, `PACC = 1.0d-7`, `CHARGEACC` and `FERMIACC` 1.0d-5. The
commented SOC lines of doll-Papacon-sd are kept.

## Outputs

- `Ni9-Ni4-Ni-Ni-Ni4-Ni9.log`: the ANT run. It was run under the job name
  `Ni9-Ni4-Ni-Ni-Ni4-Ni9.doll`, which appears in the log.
- `T.Ni9-Ni4-Ni-Ni-Ni4-Ni9.dat`: transmission for E - E_F = -3 to +3 eV in 0.01 eV steps, in
  units of e^2/h, in two blocks separated by a blank line, spin up and then spin down
- `Ni9-Ni4-Ni-Ni-Ni4-Ni9.xyz`: the geometry written by ANT. The `Ar` atoms mark the sites where
  the Bethe lattices attach to the cluster; they are not part of the calculation.

## Results

ANT.Gaussian `bin/l502.exe` built from the source of commit 9f1ad3f, Gaussian 09 D.01, NVIDIA HPC
SDK build, 4 threads, 12 November 2026.

| Quantity | Value | doll-Papacon-sd |
|---|---|---|
| Basis functions | 472 | 418 |
| SCF energy (a.u.) | -15783.436667 | -15786.572053 |
| SCF cycles | 135 | 41 |
| Device Fermi energy (eV) | -8.06 | -6.4278 |
| Bethe-lattice Fermi energy (eV) | 2.818 | 1.037 |
| T(E_F) up / down (e^2/h) | 1.022 / 2.663 | 0.871 / 2.548 |
| T(E_F) total (e^2/h) | 3.685 | 3.418 |
| Conductance G (G0) | 1.84 | 1.71 |
| Spin polarization of T(E_F), (up - down)/(up + down) | -0.45 | -0.49 |
| T total at E_F - 1 eV / E_F + 1 eV | 4.528 / 2.056 | 4.365 / 3.478 |
| Spin of the two central atoms | 0.83 / 0.83 | 0.37 / 0.37 |
| Spin of atoms 10-13 and 16-19 | 0.80 | 0.58 |
| Electrons in left / right half | 230.00 / 230.00 | 230.00 / 230.00 |
| Spin of left / right half | 11.96 / 11.96 | 11.47 / 11.47 |
| Unsuccessful density integrations (point limit) | 56 (255) | 10 (255) |
| Wall time | 25897 s | 4375 s |

The run terminates normally and passes the spectral-weight check 8 times without widening the energy
window. Charges and spins come from ANT's Mulliken population analysis
(`--- Mulliken population analysis ---` in the log). ANT splits the cluster into atoms 1-14 and 15-28;
the central atoms are 14 and 15.

## Notes

- The SCF is not fully converged. Gaussian stops on `Density matrix is not changing` (RMSDP 8.3e-6)
  with a DIIS error of 0.26, while the energy still falls by about 0.13 a.u. per cycle. A restart of
  an earlier run of this case from its own density (`ALPHA = 0.03`, `NSPINLOCK = 5`, the accuracies of
  doll-Papacon-sd and `conver=5`) stopped at -15785.30 a.u. with the same T(E_F), 1.021 / 2.664, and
  moments within 0.03, so the transport and moments are stable even though the energy is not.
- When the spin lock is released (after the 10 locked cycles), the energy rises from about -15765 to
  -13004 a.u. and the device Fermi energy drops to about -14 eV; most of the 135 cycles go into
  recovering from this. doll-Papacon-sd does not show it. A run with `ALPHA = 0.03` and
  `NSPINLOCK = 5` did the same.
- Compared with doll-Papacon-sd, the central moments are larger (0.83 against 0.37) and the spin-up
  transmission is nearly flat at about 1.0 from E_F to +1.2 eV, where doll-Papacon-sd has a peak of
  2.1 at +0.8 eV. Both the electrodes and the outer basis differ, so the comparison does not separate
  their effects.
- Warning: the route has no `nosymm` (`scf=(dsymm)`), so Gaussian reorients the contact as in
  doll-Papacon-sd: X = y_in, Y = x_in, Z = -z_in. The moments here are collinear and without SOC, so
  the results do not depend on it.
- The `.gjf` uses `%NProcShared=4`; change it to the number of cores available. It substitutes links
  from `../../../../bin`.
