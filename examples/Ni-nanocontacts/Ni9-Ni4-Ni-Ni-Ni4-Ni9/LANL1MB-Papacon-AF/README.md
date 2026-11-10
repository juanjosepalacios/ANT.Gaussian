Ni nanocontact (Ni9-Ni4-Ni-Ni-Ni4-Ni9, 28 atoms) between antiparallel Bethe-lattice electrodes,
spin-polarized LSDA, with a domain wall in the constriction. Basis (LANL1MB), ECP (LANL1), geometry
and electrodes (`Papacon`) are those of [../LANL1MB-Papacon-F](../LANL1MB-Papacon-F/README.md).

## Running the example

The run starts from the converged parallel solution:

1. Copy `../LANL1MB-Papacon-F/P.Ni9-Ni4-Ni-Ni-Ni4-Ni9.F.dat` to
   `P.Ni9-Ni4-Ni-Ni-Ni4-Ni9.AF.dat` in this folder.
2. Run `Ni9-Ni4-Ni-Ni-Ni4-Ni9.AF.gjf` with `Ni9-Ni4-Ni-Ni-Ni4-Ni9.AF.ant`, with
   `scf=(maxcycle=300,conver=5,nodamp) nosymm`.

ANT reads that file at start-up as the initial density matrix. `MRSTART = 15` reverses the spin of
atoms 15-28 in it, and `UD` makes the right Bethe-lattice electrode spin-down. Without the copy, the
run does not start from the converged parallel state. Because the run restarts from a density matrix
file, ANT sets `NSPINLOCK` to 0 and the spin is free from the first cycle.

`.ant` settings: `MRSTART = 15`, `UD`, `LEADDOS`, `ALPHA = 0.02`, `NPULAY = 5`, `SL = 1.0d-3`,
`PACC = 1.0d-8`, `CHARGEACC` and `FERMIACC` 1.0d-5.

## Outputs

- `Ni9-Ni4-Ni-Ni-Ni4-Ni9.AF.log`: the ANT run
- `T.Ni9-Ni4-Ni-Ni-Ni4-Ni9.AF.dat`: transmission for E - E_F = -3 to +3 eV in 0.01 eV steps, in units
  of e^2/h, in two blocks separated by a blank line, spin up and then spin down
- `Lead1DOS.dat`, `Lead2DOS.dat`: electrode DOS (`LEADDOS`) from -18.6 to +31.4 eV in 0.1 eV steps,
  spin up and then spin down (with negative sign), in two blocks
- `Ni9-Ni4-Ni-Ni-Ni4-Ni9.AF.xyz`: the geometry written by ANT. The `Ar` atoms mark the sites where
  the Bethe lattices attach to the cluster; they are not part of the calculation.

## Results

ANT.Gaussian 2.8.4 (commit 46e8a71 sources), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads,
10 November 2026.

| Quantity | Antiparallel (this run) | Parallel ([../LANL1MB-Papacon-F](../LANL1MB-Papacon-F/README.md)) |
|---|---|---|
| SCF energy (a.u.) | -1143.033738 | -1143.037136 |
| SCF cycles (device Fermi searches) | 6 (8) | 76 (72) |
| Device Fermi energy (eV) | -8.7592 | -8.7518 |
| T(E_F) up / down (e^2/h) | 0.622 / 0.622 | 0.391 / 2.528 |
| T(E_F) total (e^2/h) | 1.244 | 2.919 |
| Conductance G (G0) | 0.62 | 1.46 |
| Spin polarization of T(E_F), (up - down)/(up + down) | 0.00 | -0.73 |
| T total at E_F - 1 eV / E_F + 1 eV | 4.539 / 1.163 | 6.063 / 1.158 |
| Spin of the two central atoms | 0.48 / -0.48 | 0.52 / 0.52 |
| Electrons in left / right half | 140.01 / 139.99 | 140.00 / 140.00 |
| Spin of left / right half | 8.09 / -8.10 | 8.15 / 8.16 |
| Unsuccessful density integrations (point limit) | 10 (127) | 67 (127) |
| Time | 22 min 8 s CPU, 349 s wall | 3 h 13 min 37 s CPU, 2985 s wall |

The run terminates normally on the `conver=5` criterion (last RMSDP 5.8e-6) and passes the
spectral-weight check. Charges and spins come from ANT's Mulliken population analysis
(`--- Mulliken population analysis ---` in the log). ANT splits the cluster into atoms 1-14 and 15-28;
the central atoms are 14 and 15.

## Notes

- The moments are antisymmetric across the contact: every atom of the right half has the moment of
  its mirror atom in the left half with the opposite sign, to within 0.004, so the domain wall sits
  between the two central atoms. The charges are those of the parallel case to within 0.02.
- Both spin channels carry the same transmission (0.622) and the conductance drops from 1.46 G0 in
  the parallel case to 0.62 G0. The magnetoresistance (G_F - G_AF) / G_AF is 135%.
- The SCF energy rises by about 0.00025 a.u. per cycle and by 0.0105 a.u. in the last cycle, while
  the density meets `conver=5`; the energy of ANT runs is not variational, because the electrodes
  exchange charge with the cluster. The parallel run also ended with a large DIIS error. The two
  energies differ by 0.0034 a.u. (0.09 eV), which is within these changes, so they do not tell which
  configuration is lower.
- The run needs 6 cycles because it starts from the converged parallel density with the spins of one
  half reversed.
- Gaussian's integration grid reports no electron-count error in this run, unlike the parallel one.
- The `.gjf` uses `%nproc=4`; change it to the number of cores available. It substitutes links from
  `../../../../bin`.
