Fe nanowire (Fe64, bcc, 64 atoms) between 1D electrodes with the spd basis, spin-polarized LSDA.

This is the [BCC-spd](../BCC-spd/README.md) calculation with 1D electrodes in place of the Bethe
lattices. The basis, geometry, charge (0), starting multiplicity (163) and `.gjf` are the same.

**Warning: the density-matrix integration is not converged in this run.** ANT prints "Unsuccessful
integration of the density matrix using up to 127 points. Continue at your own risk ..." 455 times,
in 49 of the 67 SCF cycles, including 14 times in the last one. 127 points is the maximum allowed by
the default `NIP = 6` (at most 2^(NIP+1) - 1 points). The integrations probably need a larger `NIP`
(`NIP = 7` allows 255 points), at a higher cost per cycle. Treat the results below as provisional.

## Running the example

The run starts from the converged `BCC-spd` solution, not from a QC seed:

1. Run the [BCC-spd](../BCC-spd/README.md) example (both steps).
2. Copy its `Fe64-bcc-wire.chk` and `P.Fe64-bcc-wire.dat` into this directory.
3. Run `Fe64-bcc-wire.gjf` with `Fe64-bcc-wire.ant`. ANT reads `P.Fe64-bcc-wire.dat` at start-up as
   the initial density matrix and Gaussian reads the orbitals from the checkpoint (`guess=read`).

As in the sd case, runs of this wire from scratch with 1D electrodes fail (see the
[BCC-sd README](../BCC-sd/README.md#why-two-steps)). Step 3 overwrites both copied files.

`.ant` settings: `TYPE1 = TYPE2 = 1DLEAD`, `NATPC1 = NATPC2 = 5`, `NPC = 3`, and otherwise the same
as `BCC-spd`: `DEVBRENT`, `SL = 0.03`, `NSPINLOCK = 10`, `ETA = 0.01`, `ALPHA = 0.01`,
`PACC = 1.0d-7`, `CHARGEACC` and `FERMIACC` 1.0d-4. Each 1D electrode uses unit cells of 5 atoms
(one 4-atom layer and one axis atom) and is rebuilt from the device density in every cycle. Because the restart sets `NSPINLOCK` to 0, the spin is free
from the first cycle.

## Outputs

- `Fe64-bcc-wire.log`: the ANT run
- `T.Fe64-bcc-wire.dat`: transmission for E - E_F = -3 to +3 eV in 0.01 eV steps, in units of e^2/h,
  spin up and then spin down, separated by a blank line

## Results

ANT.Gaussian 2.8.4, Gaussian 09 D.01, NVIDIA HPC SDK build, 8 threads, 8 November 2026.

| | 1D electrodes (this example) | Bethe `BCC-spd` (restart) |
|---|---|---|
| SCF energy (a.u.) | -1428.282269 | |
| SCF cycles | 67 | |
| Device Fermi energy (eV) | -5.8295 | -5.5845 |
| T(E_F) up / down (e^2/h) | 2.579 / 1.090 | 2.898 / 2.488 |
| Conductance G (G0) | 1.83 | 2.69 |
| Electrons in left / right half | 255.63 / 255.55 | 255.86 / 256.14 |
| Spin of left / right half | 101.72 / 101.79 | 105.49 / 105.24 |
| Unsuccessful density integrations | 455 | 0 |
| Time | 2 h 51 min wall (21 h 52 min CPU) | |

The run terminates normally. Charges and spins are from ANT's Mulliken population analysis. The
Bethe-lattice column is the restart described in the BCC-spd README notes.

## Notes

- Convergence is loose as well: the SCF stops at RMSDP 1.00e-5 (`conver=5`) with MaxDP 4.3e-4.
  Every 5 cycles ANT's Pulay mixing (default `NPULAY = 4`) takes a large step, which moves the
  device Fermi energy and the density periodically; a larger `NPULAY` switches these steps off.
- The two halves hold 511.17 electrons in ANT's Mulliken analysis instead of 512, while the device
  charge is 272.00003. The missing 0.83 electrons are not explained and may come from the failed
  integrations.
- Almost all the change in G is in spin down. The minority peak just above E_F (about 3.3) sits at
  +0.26 eV against +0.18 eV with Bethe electrodes, so E_F falls on its foot (T down 1.09 against
  2.49). Below E_F - 1.5 eV the 1D minority T stays at about 1.2-1.5, where the Bethe lattice
  suppresses it below 0.5. The spin-down near-gap at +0.9 to +1.3 eV drops to about 0.001 but is
  not an exact zero, unlike the gaps with Bethe electrodes.
- With the sd basis the change from Bethe to 1D electrodes is 0.10 G0; here it is 0.86 G0. See
  [../README.md](../README.md).
