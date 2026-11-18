Al-benzene-Al with the benzene plane containing the transport axis (see [../README.md](../README.md)),
spin-restricted BLYP, no gap correction.

`.ant`: `ALPHA = 0.02`, `EW1 = -5`, `EW2 = 5`, `MOLBLOCKS` with atoms 21-32. Route:
`blyp/gen pseudo=read scf=(maxcycle=300,conver=5,dsymm)`, `%NProcShared=12`, `%mem=4GB`. Without
`UPLUS`, `MOLBLOCKS` has no effect.

## Outputs

- `Al-benzene-Al.log`
- `T.Al-benzene-Al.dat`: transmission for E - E_F = -5 to 5 eV in 0.01 eV steps, in units of e^2/h
- `Al-benzene-Al.xyz`: the geometry, with the Bethe-lattice sites

## Results

140 SCF cycles, energy -303.995173 a.u., device Fermi energy -6.2499 eV, T(E_F) = 0.240; 53 s. The
occupied benzene pi pair gives the resonance at -0.30 eV.
