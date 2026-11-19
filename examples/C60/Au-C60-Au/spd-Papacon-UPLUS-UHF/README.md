Au-C60-Au (see [../README.md](../README.md)), spin-unrestricted BLYP from a singlet with `guess=mix`, with
the molecular gap correction `UPLUS = 2.0` on atoms 21-80.

`.ant`: `ALPHA = 0.02`, `NIP = 7`, `EW1 = -5`, `EW2 = 5`, `MOLBLOCKS` with atoms 21-80,
`UPLUS = 2.0`. Route:
`ublyp/gen pseudo=read guess=mix scf=(maxcycle=300,conver=5,dsymm)`, charge 0, multiplicity 1,
`%NProcShared=12`, `%mem=8GB`.

## Outputs

- `Au-C60-Au.log`: the `Molecular block` lines give the HOMO and LUMO of the molecule for each spin after the shift, with the DFT values in brackets
- `T.Au-C60-Au.dat`: transmission for E - E_F = -5 to 5 eV in 0.01 eV steps, in units of e^2/h,
  spin up, a blank line, then spin down
- `T.Au-C60-Au.sum.dat`: E - E_F, T up, T down and their sum
- `Au-C60-Au.xyz`: the geometry, with the Bethe-lattice sites

## Results

44 SCF cycles, energy -5685.995164 a.u., device Fermi energy -6.4320 eV; 10861 s. Charge +0.41 and net
spin 0.00 on the C60, 0.00 on the Au tips (Mulliken). 179.49 + 179.49 electrons on the C60 (Loewdin). HOMO (MO 179) at -1.72 eV (DFT -0.86), LUMO (MO 182) at +1.28 eV (DFT +0.34), the same for both spins. T(E_F) = 3.8e-3 (up) and 3.8e-3 (down). 0 unsuccessful density integrations with up to 255 points.
