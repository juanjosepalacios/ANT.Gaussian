Au-benzene-Au (see [../README.md](../README.md)), spin-unrestricted BLYP from a singlet with `guess=mix`, with
the molecular gap correction `UPLUS = 2.0` on atoms 21-32.

`.ant`: `ALPHA = 0.02`, `NIP = 7`, `EW1 = -5`, `EW2 = 5`, `MOLBLOCKS` with atoms 21-32,
`UPLUS = 2.0`. Route:
`ublyp/gen pseudo=read guess=mix scf=(maxcycle=300,conver=5,dsymm)`, charge 0, multiplicity 1,
`%NProcShared=12`, `%mem=8GB`.

## Outputs

- `Au-benzene-Au.log`: the `Molecular block` lines give the HOMO and LUMO of the molecule for each spin after the shift, with the DFT values in brackets
- `T.Au-benzene-Au.dat`: transmission for E - E_F = -5 to 5 eV in 0.01 eV steps, in units of e^2/h,
  spin up, a blank line, then spin down
- `T.Au-benzene-Au.sum.dat`: E - E_F, T up, T down and their sum
- `Au-benzene-Au.xyz`: the geometry, with the Bethe-lattice sites

## Results

37 SCF cycles, energy -3633.319560 a.u., device Fermi energy -6.7349 eV; 1385 s. Charge +0.32 and net
spin 0.00 on the benzene, 0.00 on the Au tips (Mulliken). 20.85 + 20.85 electrons on the benzene (Loewdin). HOMO (MO 20) at -1.79 eV (DFT -0.82), LUMO (MO 22) at +5.07 eV (DFT +4.07), the same for both spins. T(E_F) = 3.0e-3 (up) and 3.0e-3 (down). 0 unsuccessful density integrations with up to 255 points.
