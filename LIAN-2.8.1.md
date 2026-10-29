# Lian 2.8.1 source snapshot

Imported from /home/palacios/ANT/ANT.Gaussian09/v.2.8.1/nvidia
on lian.fmc.uam.es on 2026-10-29, onto upstream commit 7c2a625.

All canonical installed src/*.F90, src/*.F, src/*.c and
src/lapack_blas.f, plus BLDAT/*.dat, were compared byte for byte.
Only differing files were copied. Upstream-only files were retained.
The installed Makefile.nvidia.g09 and make.nvidia.in were also copied;
the latter contains Lian-specific compiler flags and Gaussian paths.

The source differences include molecular-modification placement,
charge quadrature and experimental QTot2 changes in device.F90,
molecular-analysis changes in MolMod.F90, and diagnostic-output
changes in SpinOrbit.F90 and cluster.F90. Electrode data differ for
Co (027), Pt (078), and Au (079). This records the installed source
state, without claiming scientific validation or binary reproducibility.

No installed executables, objects, module files, generated lowercase
.f90 files, or calculation outputs were imported. Any such files
already tracked by upstream remain part of upstream history.

Builds and calculations were not run. The existing installation,
older clone, shell configuration, and GitHub remote were not modified.
The new NVIDIA Makefile is intended to be selected explicitly from src;
supporting Gaussian build products are still required.

## Static review before development

This branch is an archival baseline, not a validated release.

- In src/device.F90, TotCharge starts with n=31 but qq is not
  initialized before the first abs(q-qq) convergence comparison.
- TotCharge and QTot2 check n > nmax only after passing x(1:n)
  and w(1:n) to gauleg. If convergence fails at n=2047, the next
  n=4095 iteration exceeds the arrays declared with nmax=2047.
  These findings require a separate fix and numerical validation.
- TotCharge's absolute convergence threshold changes from 0.01
  upstream to 0.1 in this snapshot. Mol_Sub moves from before
  density computation into the Evaluation block. Both affect
  scientific behavior and merit explicit validation.
- QTot2 is described in the installed source as experimental.
- git diff --check reports inherited trailing whitespace in
  imported files. It is retained to preserve byte identity.
- Upstream already tracks bin/l101.exe, bin/l502.exe and two
  example executables. The added ignore rules do not untrack them.
  They were not replaced and must not be assumed to match this snapshot.

## Build preparation still required

No build was performed. The clone needs obj/ and lib/ plus matching
Gaussian support products in bin/: ml502.o, l502.a, bdam1.o, caldsu.o,
and the appropriate l101.exe (and l103.exe when needed).
The installed nvidia build contains these products, but none were copied.

The imported NVIDIA Makefile reads make.nvidia.in and selects pgf90.
Its paths are specific to Lian. Fortran module dependencies are not
declared, so a future initial build should be serial.
Changing the clone's directory name changes the generated archive name.

Continue development on a separate branch from this snapshot, fixing
and validating numerical changes before treating it as a production build.
