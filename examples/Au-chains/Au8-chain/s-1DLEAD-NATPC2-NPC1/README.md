Au8-chain (8 atoms, 3 A spacing) between one-dimensional leads, LSDA, with one s orbital per atom.

Each atom has Gaussian's `SHF` ECP, which leaves one valence electron (6s), and one s shell of three
primitives (exponents 0.4409, 0.2626, 0.0617; the s shell of the
[../spd-1DLEAD](../spd-1DLEAD/README.md) basis). This is the basis of the `FCC-s` Bethe-lattice fit
(`BLDAT/fits/Au-FCC-s`).

`.ant`: `TYPE1 = TYPE2 = 1DLEAD`, `NATPC1 = NATPC2 = 2`, `NPC = 1`, `CHARGEACC` and `FERMIACC`
1.0d-5, `PACC = 1.0d-7`, `EW1 = -5`, `EW2 = 5`, `ETA = 0.001`, `ALPHA = 0.01`, `NPULAY = 5`,
`SL = 0.01`. Route: `lsda/gen pseudo=read pop=reg Scf=(conver=5,maxcycle=100,dsymm,nodamp)`,
`%NProcShared=8`.

## Outputs

- `Au8-chain.log`
- `T.Au8-chain.dat`: transmission for E - E_F = -5 to 5 eV in 0.01 eV steps, in units of e^2/h
- `T.Au8-chain.png`: plot of the transmission (log scale), shown below

![Transmission](T.Au8-chain.png)

## Results

ANT.Gaussian 3.0.0 with the bounded device Fermi search (see `BRENT.md`), Gaussian 09 D.01, NVIDIA
HPC SDK build, 8 threads, 16 November 2026.

| Quantity | Value |
|---|---|
| Basis functions | 8 |
| SCF energy (a.u.) | -1.791412 |
| SCF cycles | 27 |
| Lead Fermi energy (eV) | -2.5664 |
| Device Fermi energy (eV) | -2.9269 |
| T(E_F) (e^2/h) | 0.143 |
| T at E_F - 1 eV / E_F + 1 eV (e^2/h) | 1.890 / 1.990 |
| Unsuccessful density integrations | 33 |
| Charge of the end atoms (e) | +0.054 |

The run terminates normally in about 1 s.

## Notes

- T is below 0.01 e^2/h from E - E_F = -0.10 to -0.02 eV (minimum 6e-4 near -0.06 eV), and E_F
  lies on the rising edge of this gap. The 2-atom unit taken from the finite cluster is not exactly
  periodic, so the half-filled s band of the extracted lead opens a gap at E_F. Outside it the
  transmission matches [../s-1DLEAD-NATPC1-NPC2](../s-1DLEAD-NATPC1-NPC2/README.md).
- It also converges worse: 27 SCF cycles and 33 unsuccessful density integrations, against 17 and
  none with the 1-atom unit.
- In the first cycle Muller and the secant step out of the device energy window (shifts of
  +38 and -928 eV); the device Fermi search catches this and falls back to bisection. Without that
  guard the run stops with a failed 1D Dyson solve.
