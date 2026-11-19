Al-[4]triangulene-Al (see [../README.md](../README.md)), spin-unrestricted BLYP from a quartet guess,
with the self-consistent molecular-gap correction `UPLUS = 2.0` on the molecule.

`.ant`: `ALPHA = 0.02`, `NIP = 7`, `EW1 = -5`, `EW2 = 5`, `MOLBLOCKS` with atoms 21-68,
`UPLUS = 2.0`. Route: `ublyp/gen pseudo=read scf=(maxcycle=300,conver=5,dsymm)`, charge 0,
multiplicity 4, `%NProcShared=12`, `%mem=4GB`.

## Outputs

- `Al-triangulene-Al.log`: each cycle logs the HOMO and LUMO of each spin after and before (`DFT`)
  the shift; the final analysis lists the molecular orbitals between EW1 and EW2
- `T.Al-triangulene-Al.dat`: transmission for E - E_F = -5 to 5 eV in 0.01 eV steps, in units of
  e^2/h, spin up, a blank line, then spin down
- `T.Al-triangulene-Al.sum.dat`: E - E_F, T up, T down and their sum
- `Al-triangulene-Al.xyz`: the geometry, with the Bethe-lattice sites

## Results

43 SCF cycles, energy -1325.296751 a.u., device Fermi energy -5.1103 eV; 234 s. 106.02 spin-up and
104.99 spin-down electrons on the molecule (Loewdin), net spin 1.03 on the molecule and 0.18 on the
Al tips. Spin up: HOMO (MO 108) at -1.52 eV (DFT -0.54 eV), lowest empty level (MO 106, partly
filled) at +0.001 eV. Spin down: HOMO (MO 105) at -3.48 eV (DFT -2.48 eV), LUMO (MO 106) at
+0.14 eV. T(E_F) = 3.6e-5 (up) and 3.7e-6 (down). 62 unsuccessful density integrations with up to
255 points.
