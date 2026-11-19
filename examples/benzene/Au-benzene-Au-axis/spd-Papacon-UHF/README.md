Au-benzene-Au (see [../README.md](../README.md)), spin-unrestricted BLYP from a singlet with `guess=mix`, no
gap correction.

`.ant`: `ALPHA = 0.02`, `NIP = 7`, `EW1 = -5`, `EW2 = 5`, `MOLBLOCKS` with atoms 21-32. Route:
`ublyp/gen pseudo=read guess=mix scf=(maxcycle=300,conver=5,dsymm)`, charge 0, multiplicity 1,
`%NProcShared=12`, `%mem=8GB`. Without `UPLUS`, `MOLBLOCKS` has no effect.

## Outputs

- `Au-benzene-Au.log`
- `T.Au-benzene-Au.dat`: transmission for E - E_F = -5 to 5 eV in 0.01 eV steps, in units of e^2/h,
  spin up, a blank line, then spin down
- `T.Au-benzene-Au.sum.dat`: E - E_F, T up, T down and their sum
- `Au-benzene-Au.xyz`: the geometry, with the Bethe-lattice sites

## Results

37 SCF cycles, energy -3633.324407 a.u., device Fermi energy -6.7229 eV; 1428 s. Charge +0.38 and net
spin 0.00 on the benzene, 0.00 on the Au tips (Mulliken). T(E_F) = 3.8e-3 (up) and 3.8e-3 (down). 4 unsuccessful density integrations with up to 255 points.
