Ni nanocontact (Ni9-Ni4-Ni-Ni-Ni4-Ni9, 28 atoms) between parallel Bethe-lattice electrodes,
spin-polarized LSDA, with a minimal spd basis.

Each atom has the LANL1MB basis, one contracted shell each of 4s, 4p and 3d (9 functions), with the
Hay-Wadt large-core ECP `LANL1` (argon core, 10 valence electrons). Both are taken from Gaussian's
library through `gen` and `pseudo=read`. The charge is 0 and the starting multiplicity is 31 (155
alpha and 125 beta electrons).

## Running the example

Run `Ni9-Ni4-Ni-Ni-Ni4-Ni9.F.gjf` with `Ni9-Ni4-Ni-Ni-Ni4-Ni9.F.ant`. It starts from Gaussian's
default guess with `scf=(maxcycle=300,conver=5,nodamp) nosymm`; no QC seed is needed.

The run writes `P.Ni9-Ni4-Ni-Ni-Ni4-Ni9.F.dat`. ANT reads that file at start-up as the initial
density matrix, so delete it before repeating the run from scratch. It is not included here.

`.ant` settings: no `BLPAR`, so both electrodes use the default Ni set `Papacon` (s, p and d, the
orthogonal first-neighbour Slater-Koster fit of ferromagnetic fcc Ni from Papaconstantopoulos,
Handbook of the Band Structure of Elemental Solids, 2nd ed., 2015). `ALPHA = 0.02`, `NPULAY = 4`,
`SL = 0.02`, `MAX = 10`, `PACC = 1.0d-7`, `CHARGEACC` and `FERMIACC` 1.0d-5. `NSPINLOCK` is not given,
so spin locking is switched off once RMSDP < SL/10.

## Outputs

- `Ni9-Ni4-Ni-Ni-Ni4-Ni9.F.log`: the ANT run
- `T.Ni9-Ni4-Ni-Ni-Ni4-Ni9.F.dat`: transmission for E - E_F = -3 to +3 eV in 0.01 eV steps, in units
  of e^2/h, in two blocks separated by a blank line, spin up and then spin down
- `Ni9-Ni4-Ni-Ni-Ni4-Ni9.F.xyz`: the geometry written by ANT. The `Ar` atoms mark the sites where the
  Bethe lattices attach to the cluster; they are not part of the calculation.

## Results

ANT.Gaussian 2.8.4 (commit 46e8a71 sources), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads,
10 November 2026.

| Quantity | Value |
|---|---|
| Basis functions | 252 |
| SCF energy (a.u.) | -1143.037136 |
| SCF cycles (device Fermi searches) | 76 (72) |
| Bethe-lattice electrode Fermi energy (eV) | 8.609 |
| Device Fermi energy (eV) | -8.7518 |
| T(E_F) up / down (e^2/h) | 0.391 / 2.528 |
| T(E_F) total (e^2/h) | 2.919 |
| Conductance G (G0) | 1.46 |
| Spin polarization of T(E_F), (up - down)/(up + down) | -0.73 |
| T total at E_F - 1 eV / E_F + 1 eV | 6.063 / 1.158 |
| Spin of the two central atoms | 0.52 / 0.52 |
| Electrons in left / right half | 140.00 / 140.00 |
| Spin of left / right half | 8.15 / 8.16 |
| Unsuccessful density integrations (point limit) | 67 (127) |
| Time | 3 h 13 min 37 s CPU, 2985 s wall |

The run terminates normally. Charges and spins come from ANT's Mulliken population analysis
(`--- Mulliken population analysis ---` in the log), not from Gaussian's population analysis at the
end of the log. ANT splits the cluster into atoms 1-14 and 15-28; the central atoms are 14 and 15.

## Notes

- Gaussian stops the SCF with "Density matrix is not changing but DIIS error= 3.00D-01" after 76
  cycles, not on the `conver=5` criterion, so treat the last digits of the results with caution.
- The device Fermi energy stays between -8.4 and -13.0 eV during the run and settles at -8.72 to
  -8.75 eV over the last 8 searches, which use Muller with 6 charge evaluations each. The
  spectral-weight check passes each time it is made; the window was widened once, from
  [-111.7, 108.7] to [-121.7, 118.7] eV.
- The apex atoms 14 and 15 lose about 1.2 electrons each (8.78 electrons) and carry 0.52 muB. The
  atoms of the 4-atom layers have 9.87 electrons and 0.44 muB, and the outer atoms 10.17-10.32
  electrons and 0.50-0.69 muB. The same loss of charge on the apex atoms appears in the Ni4-Ni-Ni-Ni4
  contact with this basis ([../Ni4-LANL1MB-Papacon](../Ni4-LANL1MB-Papacon/README.md)).
- Spin-down transmission dominates at E_F. Spin-up transmission reaches 3.3 near E - E_F = -1 eV and
  falls to about 0.3 between -0.9 and -0.6 eV; spin-down transmission drops from about 2.5 to 0.55
  between +0.1 and +0.4 eV. Above +1 eV both spins carry about 1.1.
- Gaussian's integration grid misses 1.0-1.4 of the 280 electrons in 25 cycles ("El error",
  relative error about 1.3-1.8e-3, tolerance 1e-3). The most likely cause is the very diffuse 4p
  primitive (exponent 0.024) of LANL1MB. `int=ultrafine` should reduce it.
- No Bethe-lattice parameter set in `BL028.dat` was fitted with this basis.
- The wall time was measured with two other calculations running on the same 8 cores.
- The `.gjf` uses `%nproc=4`; change it to the number of cores available. It substitutes links
  from `../../../../bin`.
