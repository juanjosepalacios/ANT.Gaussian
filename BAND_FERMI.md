# Periodic 1D lead Fermi levels

Periodic 1D leads always find their Fermi levels by band filling; no switch
is required. LEADBANDS is accepted for compatibility and has no effect.
For each lead independently, solve H(k)c = epsilon S(k)c with all stored
couplings through NPC, then fill to NPCEl electrons per chosen repeated cell.
ETA gives occupations 1/2 + atan((EF-epsilon)/ETA)/pi, including the full
tail from minus infinity. No lower-bound search or Dyson charge integration
is performed. The uniform midpoint k mesh doubles from 256 to at most 65536
points. Root charge residual and energy bracket width must meet CHARGEACC
and FERMIACC; consecutive meshes must also agree in Fermi energy and in the
charge evaluated at the previous mesh's Fermi energy. Failure to converge,
or a nonpositive Bloch overlap, stops the calculation. These mesh checks
are empirical error estimates. Device density, transmission and optional
LEADDOS still use Green-function solvers; ADAPTQUAD affects device contour
integration, not the lead band search.

## Validation

Initial tests (2026-11-01, eight threads, ETA=0.01 eV):

- Au10, NATPC1/2=1, NPC=3, saved starting density, FULLACC and input
  FERMIACC=CHARGEACC=1e-4: both leads EF=-6.77104161 eV, 1024 k points;
  lead searches took 0.0803 and 0.0731 seconds. The earlier contour search
  started at -6.76930 eV, with a finite lower cutoff and looser initial
  effective tolerances, so this is not an identical-convergence comparison.
  The test deliberately stopped after the two lead searches; it is not a
  completed device SCF or transmission run.
- Au30 triangular wire, NATPC1/2=3, NPC=3, fresh start, FULLACC and input
  FERMIACC=CHARGEACC=0.01: rejected the lead-1 Bloch overlap. On the initial
  256-point mesh its minimum eigenvalue was -0.00549569012 at
  k=5.436427912 radians. This is larger than numerical roundoff; no physical
  band filling was attempted. Finite-range overlap truncation is a possible
  cause and requires separate investigation. Do not clip these eigenvalues
  or interpret this failure as a converged Fermi level.

The extracted production routines pass analytic tight-binding tests with
NVIDIA and Intel, one and eight threads: complex hopping, nonorthogonal
second-neighbour overlap, spin filling, independent leads and explicit
failure on nonpositive overlap or invalid ETA. The full NVIDIA executable
build and parameter parser checks also pass.

With NPC=4, the Au30 overlap check passes on the tested meshes. Both leads
give EF=-6.41861128 eV at 512 k points, CHARGEACC=FERMIACC=0.01
and ETA=0.01 eV. Their searches take 0.2134 and 0.2162 seconds.
This lead-only test stops after both searches; it is not a full device SCF.
