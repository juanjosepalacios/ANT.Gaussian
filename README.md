# ANT.Gaussian

********************************************************************************
Versions notation: The main module ANT.F90 must contains in its header the label
with the latest version (now 2.7.2). This should coincide with the latest version 
of the  modules recently changed even though the module ANT itself may have not 
changed.  This way, a quick look tells you that some module has been changed. 
Most modules should retain old headers if they have not been changed.
********************************************************************************

It works both with intel and NVIDIA (PG) compilers. No need to maintain two versions. 
Just set your enviroment variables.

********************************************************************************
                           For Intel installations:
********************************************************************************

Set the following environment variables in our .bashrc, .profile, .bash_profile:

export MKL_NUM_THREADS=1
ulimit -s unlimited
export OMP_STACKSIZE=2G (*)

* May a require larger value as your input system and basis sets become larger **
********************************************************************************
                              Release 2.7.1                  
Main new Features:

- Implementation of 1D electrodes (at least preliminary tests are satisfactory). 
  The 1D Hamiltonians H0, V1, S0, S1 are read from HD and SD. 
  Still to check with spin and SOC
- I had to recover previous implementations of MolMod, SpinOrbit and SpinRot (please check)
- I cleaned up modules from previous 1D stuff not used (geom.F90, atomdata.F90)
- Added the possibility of fixing the Fermi energy and let the charge be
- whatever it may be
********************************************************************************
********************************************************************************
                              Release 2.7.2                  
Main new Features:

- SCFSOC implementation only for the real part of the spin diagonal blocks
- Added feature to fully reverse the magnetization based on the Hamiltonian,
  not on the density matrix
********************************************************************************
