Al-[4]triangulene-Al (see [../README.md](../README.md)), spin-unrestricted BLYP from a quartet guess,
no gap correction.

`.ant`: `ALPHA = 0.02`, `NIP = 7`, `EW1 = -5`, `EW2 = 5`, `MOLBLOCKS` with atoms 21-68. Route:
`ublyp/gen pseudo=read scf=(maxcycle=300,conver=5,dsymm)`, charge 0, multiplicity 4,
`%NProcShared=12`, `%mem=4GB`. Without `UPLUS`, `MOLBLOCKS` has no effect.

## Outputs

- `Al-triangulene-Al.log`
- `T.Al-triangulene-Al.dat`: transmission for E - E_F = -5 to 5 eV in 0.01 eV steps, in units of
  e^2/h, spin up, a blank line, then spin down
- `T.Al-triangulene-Al.sum.dat`: E - E_F, T up, T down and their sum
- `Al-triangulene-Al.xyz`: the geometry, with the Bethe-lattice sites

## Results

78 SCF cycles, energy -1325.301449 a.u., device Fermi energy -5.3685 eV; 526 s. Net spin 1.27 on
the molecule and 0.15 on the Al tips. T(E_F) = 0.29 (spin up, on a narrow resonance at E_F) and
6.3e-6 (spin down). 362 unsuccessful density integrations with up to 255 points.
