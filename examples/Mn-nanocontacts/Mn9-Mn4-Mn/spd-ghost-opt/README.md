Mn9-Mn4-Mn (multiplicity 71), UBPBE with `guess=mix`, spd basis and ECP (LANL2DZ on part of the
atoms), `Papacon` electrode on one side and a ghost on the other, with a geometry optimization
(`opt(z-matrix)`, `nosymm`).

`.ant`: `TYPE2 = GHOST`, `ALPHA = 0.030`, `NPULAY = 3`.

With ANT.Gaussian 3.0.0 (13 November 2026) it does not finish: after 101 SCF cycles of the first
geometry ANT stops with `Bisec/Error: Root is not bracketed`. No outputs are stored.
