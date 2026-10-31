# ANT.Gaussian (2.8.4)

********************************************************************************
Release versions are recorded in VERSION. See CHANGELOG.md for release notes.
For exact source history use Git commits; release tags identify snapshots.
The master header in src/ANT.F90 records the current release. Changed source
files carry the release in which they last changed; unchanged modules retain
their existing headers. Record the commit and any uncommitted changes with
calculation results. Build-time Git metadata is not yet embedded in the binary.
********************************************************************************

## Periodic 1D leads in 2.8.4

Fermi searches use safeguarded Brent iteration, independent warm brackets
for each lead, and a lower-bound seed from repeated-supercell generalized
eigenvalues. Parallel charge integration rejects unconverged results.

The default Dyson mixing is L1DALPHA=0.75; override it in the .ant input.
SELFACC=1.0d-6 controls the mean absolute unmixed Dyson residual.
L1DALPHA is separate from density mixing ALPHA.
FERMIACC (eV) and CHARGEACC (electrons per primitive cell for 1D leads)
have input defaults of 1.0d-5. SCF accuracy adaptation can loosen these
targets; FULLACC retains input accuracy from the first cycle.

See [Brent search](BRENT.md), [warm brackets](WARM_START.md),
[Dyson mixing benchmarks](DYSON_MIXING.md) and
[the keyword reference](doc/keywords.txt). Historical release notes follow.

It works both with intel (ifx) and NVIDIA (PG) compilers. No need to maintain 
two versions.  Just set your enviroment variables.

********************************************************************************
                           For Intel installations:
********************************************************************************

Set the following environment variables in our .bashrc, .profile, .bash_profile:

export MKL_NUM_THREADS=1
ulimit -s unlimited
export OMP_STACKSIZE=2G (*)

* May a require larger value as your input system and basis sets become larger **
********************************************************************************
                              Release 2.8.1 
Main new Features:

- Molecular DFT+U selfconsistency
- Various attempts at making 1DLeads work, but are unsatisfying at this point
********************************************************************************
                              Release 2.8.0 
Main new Features:

- SOC implementation has beed revised and fixed. Now Yukawa potential actually works
- SOC keywords have been changed for clarity
- SCFSOC seems to work, at least if the magentization points in the z direction since
  the magnetic density matrix is real in this case.
********************************************************************************
                             Release 2.7.2                  
Main new Features:

- SCFSOC implementation only for the real part of the spin diagonal blocks (added by Wynand)
- Added feature to fully reverse the magnetization based on the Hamiltonian,
  not on the density matrix (not fully tested)
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
