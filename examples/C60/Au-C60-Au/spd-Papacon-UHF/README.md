Au-C60-Au (see [../README.md](../README.md)), spin-unrestricted BLYP from a singlet with `guess=mix`, no
gap correction.

`.ant`: `ALPHA = 0.02`, `NIP = 7`, `EW1 = -5`, `EW2 = 5`, `MOLBLOCKS` with atoms 21-80. Route:
`ublyp/gen pseudo=read guess=mix scf=(maxcycle=300,conver=5,dsymm)`, charge 0, multiplicity 1,
`%NProcShared=12`, `%mem=8GB`. Without `UPLUS`, `MOLBLOCKS` has no effect.

## Outputs

- `Au-C60-Au.log`
- `T.Au-C60-Au.dat`: transmission for E - E_F = -5 to 5 eV in 0.01 eV steps, in units of e^2/h,
  spin up, a blank line, then spin down
- `T.Au-C60-Au.sum.dat`: E - E_F, T up, T down and their sum
- `Au-C60-Au.xyz`: the geometry, with the Bethe-lattice sites

## Results

49 SCF cycles, energy -5686.028654 a.u., device Fermi energy -6.3725 eV; 11961 s. Charge +0.78 and net
spin 0.00 on the C60, 0.00 on the Au tips (Mulliken). T(E_F) = 2.5e-2 (up) and 2.5e-2 (down). 9 unsuccessful density integrations with up to 255 points.
