Small Ni nanocontact (Ni4-Ni-Ni-Ni4, 10 atoms) between Bethe-lattice electrodes, spin-polarized LSDA,
with a minimal spd basis.

Each atom has the LANL1MB basis, one contracted shell each of 4s, 4p and 3d (9 functions), with the
Hay-Wadt large-core ECP `LANL1` (argon core, 10 valence electrons). Both are taken from Gaussian's
library through `gen` and `pseudo=read`. The charge is 0 and the starting multiplicity is 15 (57
alpha and 43 beta electrons). A Zeeman field of 20 T along z is applied.

## Running the example

Run `Ni4-Ni-Ni-Ni4.gjf` with `Ni4-Ni-Ni-Ni4.ant`. It starts from Gaussian's default guess with
`scf=(conver=6,nodamp,nosymm)`; no QC seed is needed.

The run writes `P.Ni4-Ni-Ni-Ni4.dat`. ANT reads that file at start-up as the initial density matrix,
so delete it before repeating the run from scratch. It is not included here.

`.ant` settings: Bethe-lattice electrodes with the Ni parameters `Papacon` (`BLPAR1`, `BLPAR2`;
s, p and d, the orthogonal first-neighbour Slater-Koster fit of ferromagnetic fcc Ni from
Papaconstantopoulos, Handbook of the Band Structure of Elemental Solids, 2nd ed., 2015), `ZM` with
`Bz = 20.0`, `NSPINLOCK = 0`, `SL = 1.0d-3`, `ALPHA = 0.01`, `NPULAY = 5`, `PACC = 1.0d-8`,
`CHARGEACC` and `FERMIACC` 1.0d-5. The commented `ZMATOM`, `SPINROTATOM` and `SPINEDIT` blocks
show per-atom fields, spin rotations and spin edits.

## Outputs

- `Ni4-Ni-Ni-Ni4.log`: the ANT run
- `T.Ni4-Ni-Ni-Ni4.dat`: total transmission for E - E_F = -3 to +3 eV in 0.01 eV steps, in units of
  e^2/h, in one block (spinor basis, because of the Zeeman field)
- `Ni4-Ni-Ni-Ni4.xyz`: the geometry written by ANT. The `Ar` atoms mark the sites where the Bethe
  lattices attach to the cluster; they are not part of the calculation.

## Results

ANT.Gaussian 2.8.4 (commit 46e8a71 sources), Gaussian 09 D.01, NVIDIA HPC SDK build, 4 threads,
10 November 2026.

| Quantity | Value |
|---|---|
| Basis functions | 90 |
| SCF energy (a.u.) | -407.755067 |
| SCF cycles (device Fermi searches) | 63 (50) |
| Bethe-lattice electrode Fermi energy (eV) | 8.609 |
| Device Fermi energy (eV) | -9.0523 |
| T(E_F) total (e^2/h) | 2.212 |
| Conductance G (G0) | 1.11 |
| T total at E_F - 1 eV / E_F + 1 eV | 3.911 / 1.815 |
| <Sz> of the two central atoms | 0.44 / 0.44 |
| Electrons in left / right half | 50.00 / 50.00 |
| Spin of left / right half | 3.18 / 3.18 |
| Unsuccessful density integrations (point limit) | 19 (127) |
| Time | 6 min 40 s CPU, 111 s wall |

The run terminates normally and passes the spectral-weight check each time it is made. Charges and
spins come from ANT's Mulliken population analysis (`--- Mulliken population analysis ---` in the
log). ANT splits the cluster into atoms 1-5 and 6-10.

## Notes

- The central atoms lose charge to the outer layers: 9.31 electrons on atoms 5 and 6, 10.17 on the
  outer atoms. Their moments (0.44) are smaller than those of the outer atoms (0.68).
- The device Fermi energy moves from -6.5 to -6.8 eV in the first 5 searches, jumps to -14.0 eV at
  the first Pulay step and settles between -9.0 and -9.2 eV from search 24. The final searches use
  Muller with 6 charge evaluations each.
- Gaussian's integration grid misses 0.53-0.95 of the 100 electrons in 12 cycles ("El error",
  relative error about 2-3e-3, tolerance 1e-3). The most likely cause is the very diffuse 4p
  primitive (exponent 0.024) of LANL1MB. `int=ultrafine` should reduce it.
- No Bethe-lattice parameter set in `BL028.dat` was fitted with this basis.
- The wall time was measured with two other calculations running on the same 8 cores.
- The `.gjf` uses `%NProcShared=4`; change it to the number of cores available. It substitutes links
  from `../../../../bin`.
