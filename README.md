# ANT.Gaussian
New version with implementations of SOC and spin rotations by Wynand Dednam
It also works both with intel and PGI compilers. No need to maintain two versions. Just set your enviroments variables.

********************************************************************************
                           For Intel installations:
********************************************************************************

Set the following environment variables in our .bashrc, .profile, .bash_profile:

export MKL_NUM_THREADS=1
ulimit -s unlimited
export OMP_STACKSIZE=2G **

** May a require larger value as your input system and basis sets become larger **

ANT.G-2.5.1 is intended to contain the latest changes of Wynand and Palacios already in 2.5.0, now including the evaluation of the polarization. 

I think the issue of the integration of the advanced density matrix is solved. Some care needed when evaluating the Dyson equation

********************************************************************************
Just added test version for the manual opening of the gap of molecules
********************************************************************************
********************************************************************************
                              Release 2.7.0                  
Main new Features:

- Implementation of 1D electrodes (at least preliminary tests are satisfactory). 
  The 1D Hamiltonians H0, V1, S0, S1 are read from HD and SD. 
  Still to check with spin and SOC
- I had to recover previous implementations of MolMod, SpinOrbit and SpinRot (please check)
- I cleaned up modules from previous 1D stuff not used (geom.F90, atomdata.F90)
********************************************************************************
