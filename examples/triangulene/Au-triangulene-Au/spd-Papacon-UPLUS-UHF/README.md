Au-triangulene-Au (see [../README.md](../README.md)), spin-unrestricted BLYP from a quartet, with
the molecular gap correction `UPLUS = 2.0` on atoms 21-68.

`.ant`: `ALPHA = 0.02`, `NIP = 7`, `EW1 = -5`, `EW2 = 5`, `MOLBLOCKS` with atoms 21-68,
`UPLUS = 2.0`. Route:
`ublyp/gen pseudo=read scf=(maxcycle=300,conver=5,dsymm)`, charge 0, multiplicity 4,
`%NProcShared=12`, `%mem=8GB`.

## Outputs

- `Au-triangulene-Au.log`: the `Molecular block` lines give the HOMO and LUMO of the molecule for each spin after the shift, with the DFT values in brackets
- `T.Au-triangulene-Au.dat`: transmission for E - E_F = -5 to 5 eV in 0.01 eV steps, in units of e^2/h,
  spin up, a blank line, then spin down
- `T.Au-triangulene-Au.sum.dat`: E - E_F, T up, T down and their sum
- `Au-triangulene-Au.xyz`: the geometry, with the Bethe-lattice sites

## Results

66 SCF cycles, energy -4667.135131 a.u., device Fermi energy -6.1426 eV; 16383 s. Charge +0.91 and net
spin 2.01 on the flake, -0.05 on the Au tips (Mulliken). 107.02 + 105.02 electrons on the flake (Loewdin). Spin up: HOMO (MO 108) at -0.03 eV (DFT +0.37) and LUMO (MO 106) at +0.00 eV (DFT +0.32), a partly filled level at E_F; spin down: HOMO (MO 105) at -2.44 eV (DFT -1.45), LUMO (MO 106) at +1.64 eV (DFT +0.70). T(E_F) = 2.8e-2 (up) and 4.9e-7 (down). 200 unsuccessful density integrations with up to 255 points.
