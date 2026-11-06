Fe nanowire (Fe64, bcc, 64 atoms) between Bethe-lattice electrodes, spin-polarized LSDA.

The wire is a bcc Fe column along z with a = 2.87 A: 13 square layers of 4 atoms at
(+-1.435, +-1.435) alternate with 12 single atoms on the axis, from z = -17.22 to +17.22 A. Each atom
has a minimal sd basis and an ECP (`gen`, `pseudo=read`). The charge is 0 and the starting
multiplicity is 163 (337 alpha and 175 beta electrons).

## Running the example

The calculation runs in two steps that share the checkpoint `Fe64-bcc-wire.chk`:

1. Run `Fe64-bcc-wire.qc.gjf` with plain Gaussian (no ANT links). It runs quadratically convergent
   SCF to a loose threshold, `scf=(qc,intrep,maxcycle=512,conver=3)` (26 cycles), and writes
   `Fe64-bcc-wire.chk`.
2. Run `Fe64-bcc-wire.gjf` with `Fe64-bcc-wire.ant`. It reads the starting orbitals from the
   checkpoint (`guess=read`) and does the ANT calculation with
   `scf=(maxcycle=300,conver=5,dsymm,fermi,nodiis,vshift=500,nodamp)`.

Step 1 must end normally: Gaussian writes the basis set to the checkpoint only on normal termination,
and step 2 stops with "Basis set data is not on the checkpoint file" otherwise. Step 2 overwrites the
checkpoint, so run step 1 again before repeating step 2.

Step 2 also writes `P.Fe64-bcc-wire.dat`. ANT reads that file at start-up as the initial density
matrix, so delete it before repeating step 2 from the checkpoint. It is not included here.

`.ant` settings: Bethe-lattice electrodes with the Fe parameters `BCC-sd` (`BLPAR1`, `BLPAR2`,
s+d, matching the sd basis, fitted to the CRYSTAL23 LSDA DOS of bulk bcc Fe with this basis), `DEVBRENT`, `SL = 0.03`, `NSPINLOCK = 10`, `ETA = 0.01`,
`ALPHA = 0.01`, `PACC = 1.0d-7`, `CHARGEACC` and `FERMIACC` 1.0d-4.

### Why two steps

Started from Gaussian's default guess with DIIS, the SCF of this wire diverges at cycle 3 (the energy
rises by about 700 a.u.). ANT connects the electrodes when Gaussian's RMSDP falls below `SL`, so they
are attached to that broken density, and ANT then stops with "ERROR: device Fermi bracket failed".
Runs from scratch failed this way with 1D electrodes (`1DLEAD`, `SL = 0.05` with `NATPC1/2 = 10`;
`SL = 0.02` with `NATPC1/2 = 5`, `NPC = 3`) and with Bethe electrodes (`SL = 0.05`; `SL = 0.03` with
`DEVBRENT`). Runs with `SL = 0.01` or `0.001`, or with the Fermi energy fixed at -5 eV
(`FIXEFERMI`), did not converge either, and `SCF=XQC` stalls at RMSDP 0.02. With `BCC-sd`
electrodes, `SL = 0.01` failed the same way from scratch and from the `conver=3` seed, with and
without `nodamp`.

The two steps avoid this. The QC SCF gives a reasonable starting density without the DIIS blow-up,
and in step 2 Fermi broadening with level shifting and no DIIS (`fermi,nodiis,vshift=500`) keeps
the early cycles stable until the electrodes connect at `SL = 0.03`. `nodamp` switches off Gaussian's
dynamic damping of the density updates. `SCF=QC` cannot be used in step 2 itself,
because it runs in link 508 and ANT replaces only link 502. `intrep` is needed in step 1; without it
QC stops with "Density matrix breaks symmetry".

## Outputs

- `Fe64-bcc-wire.qc.log`: step 1
- `Fe64-bcc-wire.log`: step 2
- `T.Fe64-bcc-wire.dat`: transmission for E - E_F = -3 to +3 eV in 0.01 eV steps, in units of e^2/h,
  in two blocks separated by a blank line, spin up and then spin down
- `Fe64-bcc-wire.xyz`: the geometry written by ANT, in Gaussian's standard orientation. The 18 `Ar`
  atoms mark the sites where the Bethe lattices attach to the cluster; they are not part of the
  calculation.

## Results

ANT.Gaussian 2.8.4 (commit 92cd8ae), Gaussian 09 D.01, NVIDIA HPC SDK build, 8 threads,
5 November 2026.

| | Step 1 (QC) | Step 2 (ANT) |
|---|---|---|
| Basis functions | 384 | 384 |
| SCF energy (a.u.) | -1421.003146 | -1423.363367 |
| SCF cycles (device Fermi searches) | 26 | 151 (162) |
| Device Fermi energy (eV) | | -1.1835 |
| T(E_F) up / down (e^2/h) | | 2.796 / 1.815 |
| T(E_F) total (e^2/h) | | 4.611 |
| Conductance G (G0) | | 2.31 |
| Spin polarization of T(E_F), (up - down)/(up + down) | | 0.21 |
| T total at E_F - 1 eV / E_F + 1 eV | | 2.678 / 4.629 |
| Electrons in left / right half | | 256.00 / 256.00 |
| Spin of left / right half | | 102.59 / 102.59 |
| Unsuccessful density integrations | | 3 |
| Time | 26 min 27 s CPU | 2 h 13 min 50 s wall |

Both steps terminate normally. Charges and spins come from ANT's Mulliken population analysis
(`--- Mulliken population analysis ---` in the log), not from Gaussian's population analysis at the
end of the log. ANT splits the cluster into atoms 1-32 and 33-64, which cuts the middle layer
(atoms 31-34, z = 0) in two.

## Notes

- The moments are 3.30-3.33 on the atoms of the inner 4-atom layers, 2.78-2.89 on the axis atoms and
  3.19 on the outermost layer. The solution is symmetric: atoms i and 65 - i have the same charge
  and spin. Each half carries 256 electrons, so the wire is neutral.
- After `NSPINLOCK = 10` cycles the spin is free; the run ends with 358.59 alpha and 153.41 beta
  electrons, a total spin of 205.2 against 162 at the start.
- The wire has 3 spin-up and 2 spin-down open channels at E_F (from the band structure of the
  periodic wire), so T(E_F) is 93% and 91% of the maximum in the two spins.
- Spin-down T is zero between about E - E_F = -0.55 and -0.48 eV. This is a gap in the minority
  DOS of the `BCC-sd` Bethe lattice, not a feature of the wire.
- The device Fermi energy wanders between about -1.4 and +13 eV for most of the run and settles
  near -1.18 eV only in the last 20 or so cycles. The 3 unsuccessful integrations are in cycles
  43, 79 and 113; the run continues past them.
- The SCF energy of step 2 is not variational, because the electrodes exchange charge with the
  cluster. The last cycle lowers it by 0.024 a.u. after a series of steps of about 0.0007 a.u.;
  the density still meets `conver=5`.
- Comparison with other electrodes for the same wire:

  | Electrodes | E_F (eV) | T(E_F) up + down (e^2/h) | G (G0) |
  |---|---|---|---|
  | Bethe `BCC-sd` (this example) | -1.184 | 2.796 + 1.815 = 4.611 | 2.31 |
  | 1D (`1DLEAD`, `NATPC1/2 = 5`, `NPC = 3`), restarted from a converged Bethe-lattice density | -2.181 | 2.679 + 1.572 = 4.251 | 2.13 |

- The `.gjf` files use `%NProcShared=8`; change it to the number of cores available. Step 2
  substitutes links from `../../bin`.
