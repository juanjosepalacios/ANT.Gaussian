Ni nanocontact (Ni9-Ni4-Ni-Ni-Ni4-Ni9, 28 atoms) between parallel (ferromagnetic) Bethe-lattice
electrodes, spin-polarized LSDA, with a minimal spd basis.

Each atom has the CRENBS s and d shells and ECP (argon core, 10 valence electrons) plus one contracted
4p shell (exponents 0.631, 0.233, 0.103), given in the input through `gen` and `pseudo=read`. The
charge is 0 and the starting multiplicity is 31.

**Warning:** the results below were computed with an earlier non-orthogonal `Papacon` set that has
since been removed from `BLDAT/BL028.dat`. The `.ant` gives no `BLPAR`, so a rerun now uses the
current default `Papacon`, the orthogonal first-neighbour fit from Papaconstantopoulos's handbook,
and gives different numbers. Runs of this input with the current `Papacon` and `Papacon-nonortho`
have so far not converged: the device Fermi energy runs away. With the current `Papacon`, use the
LANL1MB basis instead ([../LANL1MB-Papacon-F](../LANL1MB-Papacon-F/README.md)).

## Running the example

Run `Ni9-Ni4-Ni-Ni-Ni4-Ni9.gjf` with `Ni9-Ni4-Ni-Ni-Ni4-Ni9.ant`, with
`scf=(maxcycle=300,conver=5,nodamp) nosymm`. It writes `P.Ni9-Ni4-Ni-Ni-Ni4-Ni9.dat`.

`P.Ni9-Ni4-Ni-Ni-Ni4-Ni9.dat` from the run reported below is included. ANT reads it at start-up as
the initial density matrix, so running the input in this folder restarts from the converged density;
delete or move it to repeat the run from scratch.

`.ant` settings: `ALPHA = 0.02`, `NPULAY = 4`, `SL = 0.02`, `MAX = 10`, `PACC = 1.0d-7`, `CHARGEACC`
and `FERMIACC` 1.0d-5.

## Outputs

- `Ni9-Ni4-Ni-Ni-Ni4-Ni9.log`: the ANT run
- `T.Ni9-Ni4-Ni-Ni-Ni4-Ni9.dat`: transmission for E - E_F = -3 to +3 eV in 0.01 eV steps, in units
  of e^2/h, in two blocks separated by a blank line, spin up and then spin down
- `Ni9-Ni4-Ni-Ni-Ni4-Ni9.xyz`: the geometry written by ANT. The `Ar` atoms mark the sites where the
  Bethe lattices attach to the cluster; they are not part of the calculation.
- `P.Ni9-Ni4-Ni-Ni-Ni4-Ni9.dat`: the converged density matrix

## Results

ANT.Gaussian 2.8.4 (commit 1e4a005), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads,
3 November 2026, with the earlier `Papacon` electrodes.

| Quantity | Value |
|---|---|
| Basis functions | 252 |
| SCF energy (a.u.) | -1126.758239 |
| SCF cycles (device Fermi searches) | 111 (99) |
| Device Fermi energy (eV) | -4.3635 |
| T(E_F) up / down (e^2/h) | 0.837 / 2.864 |
| T(E_F) total (e^2/h) | 3.700 |
| Conductance G (G0) | 1.85 |
| Spin polarization of T(E_F), (up - down)/(up + down) | -0.55 |
| T total at E_F - 1 eV / E_F + 1 eV | 2.419 / 3.880 |
| Spin of the two central atoms | 0.87 / 0.85 |
| Electrons in left / right half | 139.73 / 140.27 |
| Spin of left / right half | 12.96 / 12.50 |
| Unsuccessful density integrations | 0 |
| Wall time | 1935 s |

The run terminates normally, passes the spectral-weight check without widening the energy window and
reports no Fermi-search errors. Charges and spins come from ANT's Mulliken population analysis
(`--- Mulliken population analysis ---` in the log). ANT splits the cluster into atoms 1-14 and 15-28;
the central atoms are 14 and 15.

## Notes

- The conductance is carried mostly by the spin-down channel, which has partly filled d states at
  E_F (T down / T up = 3.4).
- The `.gjf` asks for `%nproc=16`; lower it to the number of cores available. It substitutes links
  from `../../../../bin`.
