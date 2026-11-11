Ni nanocontact (Ni9-Ni4-Ni-Ni-Ni4-Ni9, 28 atoms) between parallel Bethe-lattice electrodes,
spin-polarized LSDA, with an all-electron basis on the central atoms.

Atoms 1-9 and 20-28 (the outer 9-atom layers) have the CRENBS s and d shells and ECP; atoms 10-19
(the two 4-atom layers and the two apex atoms) have an all-electron basis. Both are given in the
input through `gen` and `pseudo=read`. The geometry is slightly more open than in the other Ni9
cases: the outer layers are at z = -4.83 and +4.83 A, the 4-atom layers at -3.065 and +3.065 A and
the apex atoms at -1.30 and +1.30 A. The charge is 0 and the starting multiplicity is 37.

## Running the example

Run `Ni9-Ni4-Ni-Ni-Ni4-Ni9.gjf` with `Ni9-Ni4-Ni-Ni-Ni4-Ni9.ant`, with
`scf=(conver=4,maxcycle=300,dsymm) pop=min`.

`.ant` settings: Bethe-lattice electrodes `Papacon-sd` (`BLPAR1`, `BLPAR2`), `NSPINLOCK = 5`,
`NIP = 7`, `ALPHA = 0.03`, `SL = 1.0d-2`, `PACC = 1.0d-6`, `CHARGEACC` and `FERMIACC` 1.0d-4.

The case was set up with spin-orbit coupling (`SOC`, `SOCFACATOM`). With SOC and the default
`NIP = 6`, its density integration failed 249 times at the 127-point limit in the first 18 ANT
cycles. The run below therefore uses `NIP = 7` (up to 255 points) without SOC. The SOC lines remain
in the `.ant`, commented out.

## Outputs

- `Ni9-Ni4-Ni-Ni-Ni4-Ni9.log`: the ANT run
- `T.Ni9-Ni4-Ni-Ni-Ni4-Ni9.dat`: transmission for E - E_F = -3 to +3 eV in 0.01 eV steps, in
  units of e^2/h, in two blocks separated by a blank line, spin up and then spin down
- `Ni9-Ni4-Ni-Ni-Ni4-Ni9.xyz`: the geometry written by ANT. The `Ar` atoms mark the sites where
  the Bethe lattices attach to the cluster; they are not part of the calculation.

## Results

ANT.Gaussian 2.8.4 (commit 1e4a005), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads,
3 November 2026.

| Quantity | Value |
|---|---|
| Basis functions | 418 |
| SCF energy (a.u.) | -15786.572053 |
| SCF cycles (device Fermi searches) | 41 (35) |
| Device Fermi energy (eV) | -6.4278 |
| T(E_F) up / down (e^2/h) | 0.871 / 2.548 |
| T(E_F) total (e^2/h) | 3.418 |
| Conductance G (G0) | 1.71 |
| Spin polarization of T(E_F), (up - down)/(up + down) | -0.49 |
| T total at E_F - 1 eV / E_F + 1 eV | 4.365 / 3.478 |
| Spin of the two central atoms | 0.37 / 0.37 |
| Electrons in left / right half | 230.00 / 230.00 |
| Spin of left / right half | 11.47 / 11.47 |
| Unsuccessful density integrations (point limit) | 10 (255) |
| Wall time | 4375 s |

The run terminates normally, passes the spectral-weight check without widening the energy window and
reports no Fermi-search errors. Charges and spins come from ANT's Mulliken population analysis
(`--- Mulliken population analysis ---` in the log). ANT splits the cluster into atoms 1-14 and 15-28;
the central atoms are 14 and 15.

## Notes

- Warning: the route has no `nosymm` (`scf=(dsymm)`), so Gaussian reorients the contact and ANT
  works in its standard frame: X = y_in, Y = x_in, Z = -z_in. The cubic axes stay along X, Y and Z,
  so the Bethe-lattice d on-site energies keep their meaning. The moments here are collinear and
  without SOC, so the results do not depend on it. With SOC, or with the commented `SPINROTATOM`
  block, the spin directions refer to this frame, not to the input one. Give `nosymm` on the route
  line to keep the input frame.
- The moments of the central atoms are small: 0.58 on atoms 10-13 and 16-19 and 0.37 on atoms 14 and
  15. Most outer atoms keep about 1.0 each, so each half carries a spin of 11.47.
- Spin-down transmission dominates at E_F, as with the spd basis on all atoms.
- The `.gjf` uses `%NProcShared=4`; change it to the number of cores available. It substitutes links
  from `../../../../bin`.
