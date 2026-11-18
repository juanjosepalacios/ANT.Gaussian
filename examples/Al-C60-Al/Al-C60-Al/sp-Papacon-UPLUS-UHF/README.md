Al-C60-Al (see [../README.md](../README.md)), spin-unrestricted BLYP from a `guess=mix` singlet, with
the molecular gap correction `UPLUS = 2.0` on atoms 21-80 (the C60).

`.ant`: `ALPHA = 0.02`, `NIP = 7`, `EW1 = -5`, `EW2 = 5`, `MOLBLOCKS` with atoms 21-80,
`UPLUS = 2.0`. Route: `ublyp/gen pseudo=read guess=mix scf=(maxcycle=300,conver=5,dsymm)`, charge
0, multiplicity 1, `%NProcShared=12`, `%mem=8GB`.

## Outputs

- `Al-C60-Al.log`: the `Molecular block` lines give the HOMO and LUMO of the C60 for each spin
  after the shift, with the DFT values in brackets
- `T.Al-C60-Al.dat`: transmission for E - E_F = -5 to 5 eV in 0.01 eV steps, in units of e^2/h,
  spin up, a blank line, then spin down
- `T.Al-C60-Al.sum.dat`: E - E_F, T up, T down and their sum
- `Al-C60-Al.xyz`: the geometry, with the Bethe-lattice sites

## Results

67 SCF cycles, energy -2360.975460 a.u., device Fermi energy -5.3955 eV; 3827 s. Charge -1.15 and
net spin 0.91 on the C60, -0.03 on the Al tips; 180.68 + 179.76 electrons on the C60 (Loewdin).
Spin up: HOMO (MO 182) at -1.26 eV (DFT -0.31), LUMO (MO 181) at +0.63 eV. Spin down: HOMO (MO 180)
at -2.55 eV, LUMO (MO 182) at +0.65 eV. T(E_F) = 0.038 (up) and 0.023 (down). 5 unsuccessful density
integrations with up to 255 points.
