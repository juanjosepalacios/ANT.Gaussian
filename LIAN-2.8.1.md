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
