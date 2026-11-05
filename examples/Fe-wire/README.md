Fe nanowire (Fe64, bcc, 64 atoms) between Bethe-lattice electrodes, spin-polarized LSDA.

The wire is a bcc Fe column along z with a = 2.87 A: 13 square layers of 4 atoms at
(+-1.435, +-1.435) alternate with 12 single atoms on the axis, from z = -17.22 to +17.22 A. Each atom
has a minimal sd basis and an ECP (`gen`, `pseudo=read`). The charge is 0 and the starting
multiplicity is 163 (337 alpha and 175 beta electrons).

## Running the example

The calculation runs in two steps that share the checkpoint `Fe64-bcc-wire.chk`:

1. Run `Fe64-bcc-wire.qc.gjf` with plain Gaussian (no ANT links). It does two cycles of quadratically
   convergent SCF, `scf=(qc,intrep,maxcycle=512,conver=2)`, and writes `Fe64-bcc-wire.chk`.
2. Run `Fe64-bcc-wire.gjf` with `Fe64-bcc-wire.ant`. It reads the starting orbitals from the
   checkpoint (`guess=read`) and does the ANT calculation with
   `scf=(maxcycle=300,conver=5,dsymm,fermi,nodiis,vshift=500)`.

Step 1 must end normally: Gaussian writes the basis set to the checkpoint only on normal termination,
and step 2 stops with "Basis set data is not on the checkpoint file" otherwise. Step 2 overwrites the
checkpoint, so run step 1 again before repeating step 2.

Step 2 also writes `P.Fe64-bcc-wire.dat`. ANT reads that file at start-up as the initial density
matrix, so delete it before repeating step 2 from the checkpoint. It is not included here.

`.ant` settings: Bethe-lattice electrodes with the Fe parameters `lsda-3d4s` (`BLPAR1`, `BLPAR2`,
s+d, matching the sd basis), `DEVBRENT`, `SL = 0.03`, `NSPINLOCK = 10`, `ETA = 0.01`,
`ALPHA = 0.01`, `PACC = 1.0d-7`, `CHARGEACC` and `FERMIACC` 1.0d-4.

### Why two steps

Started from Gaussian's default guess with DIIS, the SCF of this wire diverges at cycle 3 (the energy
rises by about 700 a.u.). ANT connects the electrodes when Gaussian's RMSDP falls below `SL`, so they
are attached to that broken density, and ANT then stops with "ERROR: device Fermi bracket failed".
Runs from scratch failed this way with 1D electrodes (`1DLEAD`, `SL = 0.05` with `NATPC1/2 = 10`;
`SL = 0.02` with `NATPC1/2 = 5`, `NPC = 3`) and with Bethe electrodes (`SL = 0.05`; `SL = 0.03` with
`DEVBRENT`). Runs with `SL = 0.01` or `0.001`, or with the Fermi energy fixed at -5 eV
(`FIXEFERMI`), did not converge either, and `SCF=XQC` stalls at RMSDP 0.02.

The two steps avoid this. The QC SCF gives a reasonable starting density without the DIIS blow-up,
and in step 2 Fermi broadening with level shifting and no DIIS (`fermi,nodiis,vshift=500`) keeps
the early cycles stable until the electrodes connect. `SCF=QC` cannot be used in step 2 itself,
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

ANT.Gaussian 2.8.4 (commit e059d98), Gaussian 09 D.01, NVIDIA HPC SDK build, 8 threads,
4 November 2026.

| | Step 1 (QC) | Step 2 (ANT) |
|---|---|---|
| Basis functions | 384 | 384 |
| SCF energy (a.u.) | -1407.573973 | -1423.369196 |
| SCF cycles (device Fermi searches) | 2 | 91 (102) |
| Device Fermi energy (eV) | | -2.0219 |
| T(E_F) up / down (e^2/h) | | 2.114 / 1.568 |
| T(E_F) total (e^2/h) | | 3.682 |
| Conductance G (G0) | | 1.84 |
| Spin polarization of T(E_F), (up - down)/(up + down) | | 0.15 |
| T total at E_F - 1 eV / E_F + 1 eV | | 1.689 / 3.682 |
| Electrons in left / right half | | 256.00 / 256.00 |
| Spin of left / right half | | 102.22 / 102.22 |
| Unsuccessful density integrations | | 0 |
| Wall time | 13 s | 58 min 28 s |

Both steps terminate normally. Charges and spins come from ANT's Mulliken population analysis
(`--- Mulliken population analysis ---` in the log), not from Gaussian's population analysis at the
end of the log. ANT splits the cluster into atoms 1-32 and 33-64, which cuts the middle layer
(atoms 31-34, z = 0) in two.

## Notes

- The moments are about 3.3 on the atoms of the 4-atom layers and 2.8-2.9 on the axis atoms, and
  3.10 on the outermost layer. The solution is symmetric: atoms i and 65 - i have the same charge and
  spin. Each half carries 256 electrons, so the wire is neutral.
- After `NSPINLOCK = 10` cycles the spin is free; the run ends with 358.22 alpha and 153.78 beta
  electrons, a total spin of 204.4 against 162 at the start.
- Spin up carries slightly more transmission than spin down at E_F (2.11 against 1.57).
- The SCF energy of step 2 is not variational, because the electrodes exchange charge with the
  cluster. The last cycle raises it by 0.033 a.u. after a series of steps of about 0.0004 a.u.; the
  density still meets `conver=5`.
- An earlier run of this input with 1D electrodes (`1DLEAD`, `NATPC1/2 = 10`), started from a
  density matrix left over from a previous calculation, gave E_F = -2.094 eV, T(E_F) = 2.70 + 2.24 =
  4.94 e^2/h and 102.67 spin per half. The electrodes differ, so the numbers are not expected to agree.
- The `.gjf` files use `%NProcShared=8`; change it to the number of cores available. Step 2
  substitutes links from `../../bin`.
