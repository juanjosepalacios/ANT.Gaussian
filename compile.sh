# compile source code
gfortran -c -fPIC gaudata.f90 -lblas -llapack

# Make static library
ar crs libgd.a *.o

f2py3 --overwrite-signature -m gaudata -h sgnFile.pyf gaudata.f90

f2py3 -c --fcompiler=gnu95 sgnFile.pyf gaudata.f90 -L. -lgd


