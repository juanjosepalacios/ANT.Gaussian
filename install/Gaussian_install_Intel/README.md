# ANT.Gaussian
When installing Gaussian 09 with the Intel compilers, replace i386.make and mdutil.c under g09/bsd by
the files in this folder and remove the suffix.

- `i386.make.ifx`, `mdutil.c.ifx`: ifx/icx (oneAPI 2024) with MKL, used for the Intel build of
  Gaussian 09 D.01 that ANT is tested with. Set `MKLPATH` to your MKL `lib` directory.
- `i386.make.intel`, `mdutil.c.intel`: the earlier ifort/icc (Parallel Studio 2016) setup. With ifx
  it needs the flag below.

Keep `-fp-speculation=safe` in `RUNF77`. With ifx's default `-fp-speculation=fast`, the LAPACK
routines that Gaussian bundles in `util.a` give wrong results: `dlaed4` breaks DSYEVD (`DSYEVD-2
returned Info=...`, NaN for symmetry blocks) and `dlasd4` breaks DGESDD (`Error in orthogonal basis`
in l302), and the CNDO, INDO and `int=am1/pm3/pm3mm` energies of Gaussian test 0189 are wrong. With
the flag, Gaussian tests 0000-0199 match the NVIDIA build. `-axCORE-AVX2` instead of `-axAVX` was
slower on an i9-13900KF.

See http://lichingtao.blogspot.com/2014/03/gaussian-09-for-intel-compiler-with.html for more details.
