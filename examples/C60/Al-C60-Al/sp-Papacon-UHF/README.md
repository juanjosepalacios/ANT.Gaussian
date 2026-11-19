Al-C60-Al (see [../README.md](../README.md)), spin-unrestricted BLYP from a `guess=mix` singlet, no
gap correction.

`.ant`: `ALPHA = 0.02`, `NIP = 7`, `EW1 = -5`, `EW2 = 5`, `MOLBLOCKS` with atoms 21-80. Route:
`ublyp/gen pseudo=read guess=mix scf=(maxcycle=300,conver=5,dsymm)`, charge 0, multiplicity 1,
`%NProcShared=12`, `%mem=8GB`. Without `UPLUS`, `MOLBLOCKS` has no effect.

## Outputs

- `Al-C60-Al.log`
- `T.Al-C60-Al.dat`: transmission for E - E_F = -5 to 5 eV in 0.01 eV steps, in units of e^2/h,
  spin up, a blank line, then spin down
- `T.Al-C60-Al.sum.dat`: E - E_F, T up, T down and their sum
- `Al-C60-Al.xyz`: the geometry, with the Bethe-lattice sites

## Results

52 SCF cycles, energy -2360.981015 a.u., device Fermi energy -5.4918 eV; 3131 s. Charge -1.13 on
the C60, no net spin. T(E_F) = 1.38 per spin, on the t1u LUMO resonance. No unsuccessful density
integrations.
