Au-triangulene-Au (see [../README.md](../README.md)), spin-unrestricted BLYP from a quartet, no
gap correction.

`.ant`: `ALPHA = 0.02`, `NIP = 7`, `EW1 = -5`, `EW2 = 5`, `MOLBLOCKS` with atoms 21-68. Route:
`ublyp/gen pseudo=read scf=(maxcycle=300,conver=5,dsymm)`, charge 0, multiplicity 4,
`%NProcShared=12`, `%mem=8GB`. Without `UPLUS`, `MOLBLOCKS` has no effect.

## Outputs

- `Au-triangulene-Au.log`
- `T.Au-triangulene-Au.dat`: transmission for E - E_F = -5 to 5 eV in 0.01 eV steps, in units of e^2/h,
  spin up, a blank line, then spin down
- `T.Au-triangulene-Au.sum.dat`: E - E_F, T up, T down and their sum
- `Au-triangulene-Au.xyz`: the geometry, with the Bethe-lattice sites

## Results

73 SCF cycles, energy -4667.072342 a.u., device Fermi energy -6.0798 eV; 16589 s. Charge +1.18 and net
spin -0.04 on the flake, 0.06 on the Au tips (Mulliken). T(E_F) = 8.4e-4 (up) and 2.2e-4 (down). 244 unsuccessful density integrations with up to 255 points.
