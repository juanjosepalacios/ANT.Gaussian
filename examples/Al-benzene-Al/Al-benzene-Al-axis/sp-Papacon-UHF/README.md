Al-benzene-Al with the benzene plane containing the transport axis (see [../README.md](../README.md)),
spin-unrestricted BLYP with a symmetry-breaking initial guess, no gap correction.

`.ant`: `ALPHA = 0.02`, `EW1 = -5`, `EW2 = 5`, `MOLBLOCKS` with atoms 21-32. Route:
`ublyp/gen pseudo=read guess=mix scf=(maxcycle=300,conver=5,dsymm)`, `%NProcShared=12`, `%mem=4GB`.

## Outputs

- `Al-benzene-Al.log`
- `T.Al-benzene-Al.dat`: transmission for E - E_F = -5 to 5 eV in 0.01 eV steps, in units of e^2/h,
  spin up, a blank line, then spin down
- `T.Al-benzene-Al.sum.dat`: E - E_F, T up, T down and their sum
- `Al-benzene-Al.xyz`: the geometry, with the Bethe-lattice sites

## Results

135 SCF cycles, energy -303.994777 a.u., device Fermi energy -6.2380 eV, T(E_F) = 0.135 + 0.131 =
0.266; 72 s. No spin density on the benzene; small opposite moments on the two Al tips (net -0.09)
from `guess=mix`, which split the narrow resonances near -2.9 eV by spin (up to 0.44 in T).
