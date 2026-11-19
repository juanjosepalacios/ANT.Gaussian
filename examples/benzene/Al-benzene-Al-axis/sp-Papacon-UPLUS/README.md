Al-benzene-Al with the benzene plane containing the transport axis (see [../README.md](../README.md)),
spin-restricted BLYP, with the self-consistent molecular-gap correction `UPLUS = 2.0` on the benzene.

`.ant`: `ALPHA = 0.02`, `EW1 = -5`, `EW2 = 5`, `MOLBLOCKS` with atoms 21-32, `UPLUS = 2.0`. Route:
`blyp/gen pseudo=read scf=(maxcycle=300,conver=5,dsymm)`, `%NProcShared=12`, `%mem=4GB`.

## Outputs

- `Al-benzene-Al.log`: each cycle logs the benzene HOMO and LUMO after and before (`DFT`) the shift;
  the final analysis lists the molecular orbitals between EW1 and EW2
- `T.Al-benzene-Al.dat`: transmission for E - E_F = -5 to 5 eV in 0.01 eV steps, in units of e^2/h
- `Al-benzene-Al.xyz`: the geometry, with the Bethe-lattice sites

## Results

114 SCF cycles, energy -304.003482 a.u., device Fermi energy -6.2779 eV, T(E_F) = 0.094; 35 s.
HOMO (MO 21) at -0.53 eV (DFT +0.47 eV) and LUMO (MO 22)
at +6.97 eV (DFT +5.97 eV); 41.72 electrons on the benzene. The HOMO resonance moves from -0.30 to
-0.68 eV.
