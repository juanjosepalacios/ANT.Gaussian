Al-benzene-Al with the benzene plane containing the transport axis (see [../README.md](../README.md)),
spin-unrestricted BLYP with a symmetry-breaking initial guess and the self-consistent
molecular-gap correction `UPLUS = 2.0` on the benzene.

`.ant`: `ALPHA = 0.02`, `EW1 = -5`, `EW2 = 5`, `MOLBLOCKS` with atoms 21-32, `UPLUS = 2.0`. Route:
`ublyp/gen pseudo=read guess=mix scf=(maxcycle=300,conver=5,dsymm)`, `%NProcShared=12`, `%mem=4GB`.

## Outputs

- `Al-benzene-Al.log`: each cycle logs the benzene HOMO and LUMO of each spin after and before
  (`DFT`) the shift
- `T.Al-benzene-Al.dat`: transmission for E - E_F = -5 to 5 eV in 0.01 eV steps, in units of e^2/h,
  spin up, a blank line, then spin down
- `T.Al-benzene-Al.sum.dat`: E - E_F, T up, T down and their sum
- `Al-benzene-Al.xyz`: the geometry, with the Bethe-lattice sites

## Results

139 SCF cycles, energy -303.987185 a.u., device Fermi energy -6.2786 eV, T(E_F) = 0.048 + 0.048 =
0.095; 70 s. Both spins: HOMO (MO 21) at -0.53 eV (DFT +0.47 eV), LUMO (MO 22) at +6.97 eV (DFT
+5.97 eV), 20.86 electrons each on the benzene. The spin channels differ by less than 0.03 in T.
