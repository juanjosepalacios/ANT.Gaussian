Ni nanocontact (Ni9-Ni4-Ni-Ni-Ni4-Ni9, 28 atoms) between antiparallel Bethe-lattice electrodes,
spin-polarized LSDA, with a domain wall in the constriction. Basis, ECP and geometry are those of
[../Ni9-spd-F](../Ni9-spd-F/README.md).

**Warning:** as for [../Ni9-spd-F](../Ni9-spd-F/README.md), the results below were computed with an
earlier non-orthogonal `Papacon` set that has since been removed from `BLDAT/BL028.dat`, so a rerun
uses the current orthogonal `Papacon` and gives different numbers.

## Running the example

The run starts from the converged parallel solution:

1. Copy `../Ni9-spd-F/P.Ni9-Ni4-Ni-Ni-Ni4-Ni9.F.dat` to `P.Ni9-Ni4-Ni-Ni-Ni4-Ni9.AF.dat` in this
   folder.
2. Run `Ni9-Ni4-Ni-Ni-Ni4-Ni9.AF.gjf` with `Ni9-Ni4-Ni-Ni-Ni4-Ni9.AF.ant`.

ANT reads that file at start-up as the initial density matrix. `MRSTART = 15` reverses the spin of
atoms 15-28 in it, and `UD` makes the right Bethe-lattice electrode spin-down. Without the copy, the
run does not start from the converged parallel state.

`.ant` settings: `MRSTART = 15`, `UD`, `LEADDOS`, `ALPHA = 0.02`, `NPULAY = 5`, `SL = 1.0d-3`,
`PACC = 1.0d-8`, `CHARGEACC` and `FERMIACC` 1.0d-5.

## Outputs

- `Ni9-Ni4-Ni-Ni-Ni4-Ni9.AF.log`: the ANT run
- `T.Ni9-Ni4-Ni-Ni-Ni4-Ni9.AF.dat`: transmission for E - E_F = -3 to +3 eV in 0.01 eV steps, in units
  of e^2/h, in two blocks separated by a blank line, spin up and then spin down
- `Lead1DOS.dat`, `Lead2DOS.dat`: electrode DOS (`LEADDOS`)

## Results

ANT.Gaussian 2.8.4 (commit 1e4a005), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads,
3 November 2026, with the earlier `Papacon` electrodes.

| Quantity | Value |
|---|---|
| Basis functions | 252 |
| SCF energy (a.u.) | -1126.768283 |
| SCF cycles (device Fermi searches) | 19 (21) |
| Device Fermi energy (eV) | -4.2798 |
| T(E_F) up / down (e^2/h) | 1.202 / 1.227 |
| T(E_F) total (e^2/h) | 2.429 |
| Conductance G (G0) | 1.21 |
| Spin polarization of T(E_F), (up - down)/(up + down) | -0.01 |
| T total at E_F - 1 eV / E_F + 1 eV | 0.900 / 2.938 |
| Spin of the two central atoms | 0.87 / -0.86 |
| Electrons in left / right half | 139.76 / 140.24 |
| Spin of left / right half | 13.16 / -12.78 |
| Unsuccessful density integrations (point limit) | 88 (127) |
| Wall time | 376 s |

The run terminates normally, passes the spectral-weight check without widening the energy window and
reports no Fermi-search errors. Charges and spins come from ANT's Mulliken population analysis
(`--- Mulliken population analysis ---` in the log). ANT splits the cluster into atoms 1-14 and 15-28;
the central atoms are 14 and 15.

## Notes

- The antiparallel configuration is lower in energy than the parallel one by 0.0100 a.u. (0.27 eV)
  with these electrode parameters. Its moments are nearly antisymmetric across the contact (+0.87
  and -0.86 on atoms 14 and 15), so the domain wall sits between the two central atoms.
- Both spin channels carry the same transmission (1.20 and 1.23), and the conductance drops from
  1.85 G0 in the parallel case to 1.21 G0. The magnetoresistance (G_F - G_AF) / G_AF is 52%.
- The run needs 19 cycles against 111 for the parallel case, because it starts from the converged
  parallel density with the spins of one half reversed.
- The unsuccessful density integrations are at the 127-point limit (`NIP = 6`). The run converges
  normally and passes the spectral-weight check; `NIP = 7` would allow 255 points at a higher cost.
- The `.gjf` uses `%nproc=4`; change it to the number of cores available. It substitutes links from
  `../../../bin`.
