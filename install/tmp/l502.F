*Deck UHFOpn
      Subroutine UHFOpn(V,NGot)
      Implicit Real*8(A-H,O-Z)
C
C     Copyright (c) 1988,1990,1992,1993,1995,1998,2003,2009,2013
C                Gaussian, Inc.  All Rights Reserved.
C
C     This is part of the Gaussian(R) 09 program.  It is based on
C     the Gaussian(R) 03 system (copyright 2003, Gaussian, Inc.),
C     the Gaussian(R) 98 system (copyright 1998, Gaussian, Inc.),
C     the Gaussian(R) 94 system (copyright 1995, Gaussian, Inc.),
C     the Gaussian 92(TM) system (copyright 1992, Gaussian, Inc.),
C     the Gaussian 90(TM) system (copyright 1990, Gaussian, Inc.),
C     the Gaussian 88(TM) system (copyright 1988, Gaussian, Inc.),
C     the Gaussian 86(TM) system (copyright 1986, Carnegie Mellon
C     University), and the Gaussian 82(TM) system (copyright 1983,
C     Carnegie Mellon University). Gaussian is a federally registered
C     trademark of Gaussian, Inc.
C
C     This software contains proprietary and confidential information,
C     including trade secrets, belonging to Gaussian, Inc.
C
C     This software is provided under written license and may be
C     used, copied, transmitted, or stored only in accord with that
C     written license.
C
C     The following legend is applicable only to US Government contracts
C     under DFARS:
C
C                        RESTRICTED RIGHTS LEGEND
C
C     Use, duplication or disclosure by the US Government is subject to
C     restrictions as set forth in subparagraph (c)(1)(ii) of the Rights
C     in Technical Data and Computer Software clause at DFARS
C     252.227-7013.
C
C     Gaussian, Inc., 340 Quinnipiac St., Bldg. 40, Wallingford CT 06492
C
C     The following legend is applicable only to US Government contracts
C     under FAR:
C
C                        RESTRICTED RIGHTS LEGEND
C
C     Use, reproduction and disclosure by the US Government is subject
C     to restrictions as set forth in subparagraph (c) of the Commercial
C     Computer Software - Restricted Rights clause at FAR 52.227-19.
C
C     Gaussian, Inc., 340 Quinnipiac St., Bldg. 40, Wallingford CT 06492
C
C1L502
C     LINK 502 ... PROGRAM UHFOPN
C
C     Solution of the UHF (Pople-Nesbet) and ROHF equations using
C     repeated diagonalizations, with DIIS extrapolation.
C     J.CHEM.PHYS. 22, 571 (1954)
C
C1Ov5
C
C     Overlay 5 ... Gaussian system SCF routines.
C
C     In the Gaussian system, Overlay 5 consists of a series of
C  programs that perform various types of SCF calculations.
C  Currently, the following links are available:
C
C      L501 ... Closed shell SCF.  This program performs the
C               "regular" repeated diagonalization procedure
C               of Roothaan.
C
C      L502 ... Closed and open shell SCF.  This program implements the
C               Roothaan repeated diagonalization scheme to solve the
C               Roothaan-Hall (RHF), Pople-Nesbit (UHF), or
C               ROHF equations.  This link can use integrals stored in
C               memory, on disk, or recomputed (direct SCF).
C
C      L503 ... Steepest descent SCF.  This link is capable of solving
C               the Pople-Nesbit equations (either closed shell or
C               open shell) using either the method of repeated
C               diagonalizations or the method of steepest descents.
C               When convergence difficulties are encountered with
C               the regular SCF, this procedure will almost always
C               solve the problem.  This link is also capable of
C               handling real or complex orbitals.  Note however,
C               that the steepest descent method can take considerably
C               longer than the regular SCF and that this link is only
C               dimensioned for 256 basis functions.
C
C      L504 ... Handles wavefunctions involving a closed-shell and two
C               special orbitals.  These include triplets, orthogonal
C               open-shell singlets, GVB(1/2) singlets, and 2 electron,
C               2 orbital (2 X 2) CAS wavefunctions.  The latter
C               include the (20, 11, 02) representation of GVB and
C               excited singlet (hole) states without the requirement
C               of orthogonality between the open-shell orbitals.
C
C      L505 ... Obsolete restricted Hartree-Fock method for open shells.
C               This link implements the restricted Hartree-Fock
C               method of Binkley, Pople and Dobosh.
C
C      L506 ... ROHF and generalized valence-bond program.
C
C      L508 ... Quadratically convergent SCF program for RHF and UHF.
C
C      L509 ... VAX specific program which runs the Caltech MCSCF
C               system via subprocesses.
C
C      L510 ... Mike Robb's MCSCF code.
C
C2IOp(5)
C   IOp(5) ... Direct SCF control (L502, L508).
C        0 ... Default (same as 1).
C        1 ... Read the integrals off disk.
C        2 ... Compute 2e integrals.
C        3 ... Compute 2e integrals and store in-core.
C        4 ... Compute 2e integrals and forbid in-core.
C    NNNNx ... use option NNNNN in control of 2e integral calculation.
C   000000 ... Default -- incremental Fock matrix formation only for
C              direct SCF.
C   100000 ... Form full Fock matrix every time.
C   200000 ... Form delta-F each iteration
C  1000000 ... Clear in-core integrals for testing.
C
C2IOp(6)
C IOp(6)  ... Convergence (RMS density except in L506 (SQCDF),
C             L508(rms rotation gradient), and L510 (Energy)).
C
C       0 ... 10**-8, except 10**-7 for PBC.
C       N ... 10**-N.
C
C2IOp(7)
C IOp(7) ... Maximum number of iterations.
C      0 ... 128, except 512 in L503 and L508.
C
C2IOp(8)
C
C IOp(8) ... SELECTION OF THE PROCEDURE OF DIRECT MINIMIZATION (L503).
C     0 ... STEEPEST DESCENT WITH SEARCH PARAMETERS DEFAULT
C     1 ... STEEPEST DESCENT WITH SEARCH PARAMETERS READ (SEE BELOW)
C     2 ... CLASSICAL SCF (ROOTHAAN'S METHOD OF REPEATED DIAGONALIZATION
C     4 ... CONJUGATE GRADIENTS WITH SEARCH PARAMETERS DEFAULT
C     5 ... CONJUGATE GRADIENTS WITH SEARCH PARAMETERS READ (SEE BELOW)
C     THE SEARCH PARAMETERS ARE  MAX. NUMBER OF SEARCH POINTS (I1)
C                                MIN. NUMBER OF SEARCH POINTS (I1)
C                                INITIAL STEPSIZE, TAU (G18.5)
C                                SCALING FACTOR FOR SUBSEQ. TAU (G20.5)
C                                Q (G20.5)
C
C     IOp(8) ... Search method (L508).
C          0 ... Default (1123).
C          1 ... Steepest descent.
C          2 ... Scaled steepest descent.
C          3 ... Quadratic convergence (after rotation gradient is
C                sufficiently small).
C          4 ... Exit when NR point is reached, so L502 can take over.
C         00 ... Default linear search (full search).
C         10 ... Do a full linear search to locate a minimum.
C         20 ... Do a linear search only if the energy goes up after the
C                initial step.
C        000 ... Default handling of wrong curvature (switch direction).
C        100 ... Reverse direction if curvature in NR step direction is
C                wrong.
C        200 ... Take pure NR steps, even if curvature is wrong.
C       1000 ... Turn on linear search and variable step logic.
C       2000 ... Turn off linear search and variable step logic.
C
C     IOp(8) ... Flags for L510:
C          1 ... IRdF2, read damping coefficients.
C         10 ... IFrzCI, freeze CI coefficients after 1st iteration.
C        100 ... Read unformatted symbolic matrix elements from NDATA
C                instead of rwf.
C       1000 ... Read in damping factors from cards.
C      10000 ... Use Levy damping.
C     100000 ... Read Fock matrix restriction matrix.
C
C2IOp(9)
C IOp(9) ... SWITCH TO CLASSICAL SCF AFTER DENSITY MATRIX HAS
C            ACHIEVED A CERTAIN CONVERGENCY (L503 only).
C      0 ... NO
C      1 ... YES, CRITERION DEFAULT  10(**-3)
C      2 ... YES, CRITERION READ IN (FORMAT G16.10)
C
C IOp(9) ... Number of pair iterations (L504, L506).
C     -1 ... None; coefficients are frozen at initial values (L504 only,
C            causes coefficients to be read in in order 11 12 22.
C      0 ... 5.
C
C IOp(9) ... IDoV in Harris functional in L511.
C
C2IOp(10)
C
C  IVShft  Level shifting:
C          -N ... Dynamic level shifting to achieve a gap of -0.001*N
C          -2 ... Dynamic level shift to a default goal (same as -200)
C          -1 ... No level shifting.
C           0 ... Default:  -200 for diagonalization calculations,
C                           -1 for sparse diagonalization replacements,
C                           and if energy DIIS is turned on.
C           N ... Shift by 0.001*N
C
C2IOp(11)
C  IOp(11) ... 3 and 4 Point Density extrapolation control (L502,
C              L503 has only 4 point, L505).
C        0 ... BOTH THREE-POINT AND FOUR-POINT EXTRAPOLATION ARE
C              PERFORMED WHEN APPLICABLE.
C        1 ... THREE-POINT EXTRAPOLATION IS INHIBITED, BUT THE
C              PROGRAM WILL STILL PERFORM FOUR-POINT EXTRAPOLATION
C              WHEN POSSIBLE.
C        2 ... BOTH THREE-POINT AND FOUR-POINT EXTRAPOLATION SCHEMES
C              ARE 'LOCKED OUT' (IE. DISABLED).
C       00 ... Default (20).
C       10 ... Do Camp-King always, taking one step using if
C              |Lpred-1|>=Thresh.
C       20 ... Do not do Camp-King.
C       30 ... Do Camp-King only if the energy rises by at least Thresh.
C       40 ... Same as 1, but CK also done if the energy goes up.
C    NNN00 ... Threshold for CK is NNN/1000 (step for 10, Hartrees for
C              30).  Default is 0.3 for 10, 0.001 for 30; with 30, 999
C              implies 10^-10.
C
C2IOp(12)
C
C     IOp(12) ... Whether to allocate only two N**2 arrays for RHF:
C           0 ... Default (No).
C           1 ... Yes.
C           2 ... No.
C
C IOp(12) ... Number of GVB pairs (L506).  If non zero, the number of
C orbitals in each pair is read in format (30I2).  Each pair consists of
C the highest available occupied from the guess (after high spin orbs
C are accounted for) and the lowest available virtuals.  If <0, pair
C coefficients are read; otherwise standard initial values are used.
C
C2IOp(13)
C IOp(13) ... Action on convergence failure (L502):
C       0 ... Default (2).
C       1 ... Continue the run even on non-convergence.  The ILSW flag
C             for convergence failure is set.
C       2 ... Terminate on non-convergence.
C
C2IOp(14)
C     IOp(14) ... Special functions in L502:
C           0 ... None.
C           1 ... Turn the current RHF run into a uhf run at the end of
C                 this link.
C          10 ... Terminate after computing the 2e terms at the first iteration.
C          20 ... Just recompute band structure from stored real-space Fock matrix.
C         100 ... ADMP/FOSimult, later cycles:  transform the density from L103/L121
C                 before calculating the energy and Fock matrices.
C         200 ... ADMP/FOSimult, first cycle:  use initial AO densities.
C        1000 ... Use Generalized energy-weighted density routines regardless.
C        2000 ... Do not use GEW routines even for CP.
C       10000 ... Fit the converged density even if fitting is not in use during
C                 the SCF.  Also redoes the fit at the end even if using fits
C                 during SCF.
C      100000 ... Compute the 2e and XC quadrature contributions to the gradient
C                 along with the Fock matrix (NYI).
C     1000000 ... Do not call CycOpn at all; instead W, etc. are computed from the
C                 input orbitals and density.
C     2000000 ... Call CycOpn.
C     3000000 ... Call StuCyc instead of CycOpn for debugging.
C     4000000 ... Test in-core spinorbital integrals by doing 2e
C                 contribution both ways.
C    10000000 ... Do simple-minded CPHF for GHF polarizability here.
C   100000000 ... Do numerical differentiation of the Fock matrix and
C                 compare with analytic one. Requires extra input ... see
C                 routine NFmnIn. Strictly for debugging purposes.
C   200000000 ... Do numerical differentiation of the density functional
C                 compare with analytic derivatives. Requires extra input
C                 ... see routine NDXCIn. Strictly for debugging purposes.
C
C  IOp(14) ... CONTROL OF ANNIHILATION OF SPIN CONTAMINANTS (L502).
C     0 ... CALCULATION IS PERFORMED  (PROVIDED OF COURSE THAT
C           ENOUGH SPACE EXISTS IN THE RW-FILES).
C     1 ... CALCULATION IS BYPASSED.
C     2 ... CALCULATION IS PERFORMED, CONTINGENT ON SPACE, AND
C           THE SYSTEM RW-FILES FOR THE APPROPRIATE DENSITY MATRICES
C           ARE UPDATED (USEFUL IF ONE WANTS A POPULATION ANALYSIS).
C
C IOp(14) ... REORDERING OF THE ORBITALS (MAINTAINING CONTINUITY
C             OF THE WAVEFUNCTION ALONG THE SEARCH PATH, L503).
C     0 ... ON  BESSEL CRITERION
C     1 ... ON  STRONGER INDIVIDUAL-OVERLAP CRITERION
C     2 ... OFF
C
C IOp(14) ... Flags for MCSCF:
C       1 ... Skip valence-valence Fock matrix elements.
C      10 ... Skip core-valence Fock matrix elements.
C     100 ... Skip valence-virtual Fock matrix elements.
C    1000 ... Skip core-valence Fock matrix elements.
C   10000 ... Use full diagonalization method rather than Lanczos.
C             (Obsolete; use IOp(17)).
C  100000 ... State average density matrices.
C
C2IOp(15)
C IOp(15) ... Apply Abelian symmetry constraints on orbitals.
C       0 ... Default (1 for L502, 2 for L501 and L506).
C       1 ... No.
C       2 ... Yes, keep occupation of each irrep the same as the
C             initial guess.
C       3 ... Yes, keep overall wavefunction the same as the initial
C             guess, but doing the minimal amount of orbital switching
C             to accomplish this.
C      00 ... Default (use Abelian symmetry in diagonalization).
C      10 ... Use Abelian symmetry in diagonalization.
C      20 ... Do not use Abelian symmetry in diagonalization.
C      30 ... Use the blocking in the diagonalization without checking
C             off-diagonal elements (for energy decomposition).
C
C IOp(15) ... CONTROLS THE AUTOADJUSTMENT OF TAU (L503).
C     0 ... DONE
C     1 ... TAU IS KEPT FIXED
C
C2IOp(16)
C IOp(16) ... Diagonalization method (L502):
C       0 ... Default (1 for full matrices HF/DFT, -30 for semi-empirical, 4 for sparse).
C      -N ... Pseudo-diagonalization with real diagonalization every Nth cycle.
C      -1 ... Same as 3.
C       1 ... DiagD.
C       2 ... KyDiag.
C       3 ... Pseudo-diagonalization whenever allowed by internal tests.
C       4 ... CGDMS.
C       5 ... PDM.
C       6 ... CEM.
C       7 ... Sign Matrix Method.
C       8 ... SNRDMS.
C       9 ... Unused.
C      10 ... Jacobi diagonalization.
C     1xx ... Force formation of the Fock matrix using full storage.
C     2xx ... Force formation of the Fock matrix using sparse storage.
C     3xx ... Used by sparse ADMP?
C     4xx ... Force formation of the Fock matrix using sparse storage if ZDO.
C
C IOp(16) ... INHIBIT PERFORMANCE OF MINIMIZATION OF ALTERNATE
C             WAVEFUNCTION PROVIDED BY SECOND ORDER PROCEDURES (L503).
C     0 ... NO
C     1 ... YES
C
C IOp(16) ... Selection of OCBSE vectors (L506).
C       0 ... by eigenvalue.
C       1 ... by energy least change.
C       2 ... by orbital least change.
C
C IOp(16) ... Lanczos starting vector in L510:
C      -1 ... Read in eigenvector.
C       0 ... C(1) = 1.0
C       N ... C(N) = 1.0
C
C2IOp(17)
C IOp(17) ... CONDITION OFF-DIAGONAL TERMS OF THE Fock MATRIX (L503).
C               -SET TO ZERO IF Abs(F(I,J)).LE.FUZZY
C               -DELETE COUPLING TERMS BETWEEN ALMOST DEGENERATE
C                (DELTA E .LE. DEGEN) M.O. VECTORS
C     0 ... FUZZY=1.D-10, DEGEN=2.D-5
C     1 ... FUZZY AND DEGEN READ IN (2D20.14)
C
C IOp(17) ... Selection of virtual orbitals (L506).
C       0 ... Virtuals obtained by diagonalization of hamiltonians.
C       1 ... Virtuals obtained by Schmidt orthogonalization to
C             occupieds.
C
C IOp(17) ... Use of symmetry (in L502 and L508)
C             and linear equation convergence (in L508):
C       0 ... Default (1032 for 502, 1012 for 508).
C       1 ... Choose LinEq convergence based on orbital gradient.
C       2 ... Always use tight convergence.
C       3 ... Tighten convergence by an extra factor of 10.
C      10 ... If 2E symmetry is on, symmetrize Fock matrices and require
C             proper density matrix symmetry.
C      20 ... If 2E symmetry is on, replicate integrals so that density
C             matrices and wavefunctions need not be symmetric.
C      30 ... If 2E symmetry is on, choose between replicating integrals
C             and symmetrizing the Fock matrix based on whether the current
C             density matrix is symmetric.
C      40 ... Same as 30 in 502 but 20 in 508.
C     100 ... Force the density matrix to have full symmetry at the first
C             iteration.
C     200 ... Force the density matrix to have full symmetry at every
C             iteration.
C    0000 ... Default (1000)
C    1000 ... If the density matrices pass the symmetry test, symmetrize
C             them to ensure that they are exactly symmetric.
C    2000 ... Do not symmetrize the density matrices.
C   00000 ... Default (20000, except if IOp(8) requests old algorithm).
C   10000 ... Always pseudocanonicalize in L508.
C   20000 ... Only pseudocanonicalize in L508 when doing a
C             Newton-Raphson step.
C
C IOp(17)  ... MCSCF flags:
C        0 ... Orthogonalize C,O,V by separate Lowdin, then schmidt.
C        1 ... Lowdin orthogonalize C+O and V, then schmidt.
C        2 ... Just schmidt.
C       10 ... Don't use natural orbitals each iteration.  Bad for 1st
C              order method.
C      100 ... Use full 2nd order convergence.
C      200 ... 2nd order iteration at end, in preparation for CPMCSCF.
C     1000 ... Generate data for multi-reference MP2?
C    10000 ... Attempt to control root flipping in CI.
C   100000 ... Read CI vector and use it every iteration.
C  1000000 ... Use full diagonalization method rather than Lanczos.
C 10000000 ... Use State Average density matrices.(the weights 8F10.8)
C 20000000 ... Do SA and prepare for SA-CPMCSCF.
C 30000000 ... Do SA and prepare for Gradient of Energy difference.
C 40000000 ... Do SA and prepare for SA Second Derivative Computation
C                                    (terms involving 2nd order orbital
C                                     rotation derivatives not included)
C
C
C2IOp(18)
C
C  IOp(18) ... L502:  Mixing when doing damping:
C       -3 ... MO damping at all iterations.
C       -2 ... Turn off damping.
C       -1 ... Dynamic selection of density damping based on band gap and DIIS error.
C        0 ... Default (-1 unless reoptimizing during Stable=Opt).
C        N ... N/100 new density, (100-N)/100 old density.
C
C  IOp(18) ... L503: CUTOFF CRITERIA IN SYMMETRY DETERMINATION OF M.O.S (L503)
C             -SYMMETRY IS DETERMINED IF LARGEST OFF-DIAGONAL M.O.
C              FOCK-MATRIX ELEMENT Abs(F(I,J)).GE.STHRS
C             -ELEMENTS Abs(F(I,J)).LE.SPAN ARE CONSIDERED TO BE ZERO
C     0 ... STHRS=1.D-4, SPAN=5.D-7
C     1 ... STHRS AND SPAN READ IN (2D20.14)
C
C IOp(18) ... Damping (L506)
C
C2IOp(19)
C
C IOp(19) ... Over-ride integral storage control (L501, L502, L506, L508):
C      -1 ... Choose the best given amount of memory available.
C       0 ... 2 if possible, otherwise 1.
C       1 ... Forbid in-core:  force re-reading of integrals even if they
C             fit in 2 buffers if conventional, do not convert to in-core
C             if direct and enough memory for in-core is available.
C       2 ... Force allocation for 1 or 2 buffer case conventional case
C             (VV.ne.IBuf2E).
C       3 ... Force Lower-triangular in-memory storage.
C       4 ... Obsolute.
C      1x ... Save generated integrals on disk (file 610).
C      2x ... Force computation of raff 1 and 2 integrals even for RHF.
C      3x ... Do not save integrals (same as 0x).
C
C IOp(19) ... PRINT F(1),T. (READ ONE CARD WITH START,END 2I2) (L503).
C       0 ... NO
C       1 ... YES
C
C2IOp(20)
C IOp(20) ... Final non-DIIS iteration (L501, L502, L504).
C       0 ... Default (only if doing pseudodiagonalization or QNDMS).
C       1 ... Yes, do a final unextrapolated diagonalization after
C             convergence is reached.
C       2 ... No, just quit when extrapolated convergence is reached.
C       3 ... Do a full diagonalization at the end without recomputing
C             a new Fock matrix.
C
C IOp(20) ... Orbital rotation control (L506).
C       N ... Rotations are turned on when SQCDF is below 10**(-N).
C
C
C2IOp(21)
C
C IOp(21) ... DIIS error for density damping, maximum virtual mixing
C             for MO damping:
C
C       For density damping:
C       0 ... Default (Damp if error > 0.001)
C       N ... Damp if error > 10**-N
C
C       For MO damping:
C       0 ... Default, no more than 1/3 virtual component for any occupied
C             at each iteration.
C       N ... Maximum N/1000 virtual component.
C
C IOp(21) ... ACTION IF OTEST DETECTS PROBLEMS (L503).
C     0 ... ABORT RUN VIA LNK1E.
C     1 ... CONTINUE RUN.
C     2 ... Check orthonormality at first iteration.
C
C IOp(21) ... Extrapolation control in L506.
C
C IOp(21) ... MCSCF flags:
C       2 ... Generate MOs using UHF natural orbitals.
C      10 ... IRdNLp.
C
C2IOp(22)
C     IOp(22) ... Use of DIIS extrapolation (L501, L502, L504).
C           0 ... Default (1042) for calculations using diagonalization
C                   (2) for calculations using sparse diagonalization
C                   replacements.
C           1 ... No.
C           2 ... Yes.
C           3 ... Yes, with Fermi broadening as well, deciding on the fly between
C                 the two forms.
C           4 ... Yes, with "pFON" version of Fermi broadening.
C           5 ... Yes, with "FON" version of Fermi broadening.
C          10 ... Regular DIIS
C          20 ... Energy-based mixing
C          30 ... Energy DIIS when DIIS error has increased significantly or
C                 is above threshold
C          40 ... Energy DIIS when DIIS error has increased significantly,
C                 otherwise, mixture of energy and commutator.
C         1xx ... Use energy DIIS when commutator gives huge coefficients.
C        Nxxx ... Switch from energy to commutator when error is 10^(-N)
C                 in method 3; use (DIIS error/10^-N) for weight of energy
C                 DIIS in method 4.
C       Mxxxx ... Use print level M in DIIS.
C      Lxxxxx ... Whether to include PCM charges in DIIS (default 2).
C      100000 ... Yes.
C      200000 ... No
C     Kxxxxxx ... Number of PCM-only iterations to do per Fock formation
C                 when including PCM in DIIS.  Default is 0.
C
C     IOp(22) ... Orbital mixing control in L506.
C
C2IOp(23)
C IOp(23) ... Flag for later points of an optimization, so that pair
C             and hamiltonian information can be reused (L506, L509).
C       0 ... Read from input stream.
C       1 ... Read from rwf.
C       2 ... Read from chk.
C
C2IOp(24)
C IOp(24) ... Orbital freezing (L506).
C       0 ... Optimize all orbitals.
C       1 ... Freeze all closed, high spin and first natural orbitals.
C             Optimize only 2nd and higher naturals.
C
C2IOp(25)
C
C IOp(25) ... Rotation application (L506).
C       0 ... Default (exponentiate rotation angles).
C       1 ... Apply rotations sequentially.
C
C2IOp(26)
C IOp(26) ... Type of calculation (L504).
C       3 ... 3rd root of CAS(2,2)
C       2 ... Excited singlet as 2nd root of CAS(2,2).
C       1 ... GVB as CAS(2,2)
C       0 ... GVB(1/2)
C      -1 ... Orthogonal open-shell singlet.
C      -2 ... ROHF Triplet (a debugging option).
C
C IOp(26) ... Number of hamiltonians to read in (L506).  If zero, the
C             unpaired orbitals are assumed to be high spin.  If
C             -1, an open-shell singlet is assumed.
C
C IOp(26) ... Closed/Open control for L511:
C       0 ... Default, closed for Multip=1.
C       1 ... Force closed shell, error if Multip>1.
C       2 ... Force UHF.
C
C2IOp(27)
C     IOp(27) ... Spinor testing:
C           0 ... Use default IOpCl.
C           N ... Set generalized IOpCl value N.
C
C2IOp(28)
C IOp(28) ... Root of CI to use in MCSCF.
C       0 ... defaults to 1.
C
C2IOp(29)
C IOp(29) ... Use of rafinetti integrals during direct SCF.
C      -1 ... All integrals done as Raffenetti.
C       0 ... Default: let FoFJK decide.  It will never use
C             Raffenetti for SCF.
C       1 ... All integrals are done as regular integrals.
C       N ... Integrals with degree of contraction greater than
C             or equal to N are done are regular integrals.
C
C2IOp(30)
C IOp(30) ... whether to symmetrize final orbitals using abelian
C             symmetry operations (L502, L505, not needed in L506).
C     0 ... Default (Yes).
C     1 ... Yes.
C     2 ... No.
C     3 ... Symmetrize even if symmetry blocking was done, and print
C           symmetries.
C
C2IOp(31)
C IOp(31) ... How many vectors to form at a time during microiterations
C             in L508 (NYI) and L509:
C       0 ... Default (3 in L509).
C       N ... N.
C
C2IOp(32)
C IOp(32) ... Sleazy SCF (L502, L510):
C       0 ... Default (No).
C       1 ... Yes, use loose integral cutoffs, convergence on either
C             energy or density and always do incremental Fock formation.
C       2 ... No.
C       3 ... Thresholds similar to DGauss for convergence and integrals.
C       4 ... Yes, doing an inexpensive pass 0 and then full accuracy in pass 1.
C       5 ... Decide between 1 and 4 based on details of the calculation.
C       6 ... Do iterations with sleazy XC grid, then one iteration with
C             next grid up.  The default is CoarseGrid for iterations and
C             SG1 for final energy.
C      00 ... No longer used.
C     N00 ... No longer used.
C    I000 ... Use approximation I, 0=normal 1=Linear approximation to Xc.
C   00000 ... Use general DBF logic only if the DBF rwf is present.
C   10000 ... Force use of 1c instead of general DBF logic.
C   20000 ... Force use of general DBF logic.
C
C2IOp(33)
C  PRINT    IOp(33) ... PRINT OPTION.
C           0 ... ONLY SUMMARY RESULTS ARE PRINTED (WITH POSSIBLE
C                 CONTROL FROM THE 'NO-PRINT' OPTION).
C           1 ... THE EIGENVALUES AND THE M. O. COEFFICIENTS ARE
C                 PRINTED AT THE END OF THE SCF.
C           2 ... SAME AS IOp(33)=1, BUT ADDITIONALLY THE DENSITY MATRIX
C                 IS PRINTED.
C           3 ... SAME AS IOp(33)=2, BUT AT THE END OF EACH ITERATION.
C           4 ... SAME AS IOp(33)=3, BUT ALL MATRIX TRANSACTIONS ARE
C                 PRINTED (BEWARE!!!   MUCH OUTPUT EVEN ON SMALL
C                 MOLECULES.)
C
C2IOp(34)
C  DUMP     IOp(34) ... DUMP OPTION.
C           REGULAR SYSTEM DEFAULTS APPLY HERE.
C
C2IOp(35)
C     IOp(35) ... Whether basis is orthonormal (L501, L502).
C           0 ... Default (No).
C           1 ... Yes.
C           2 ... No.
C
C2IOp(36)
C     IOp(36) ... Whether to checkpoint after every SCF cycle.
C           0 ... Default (checkpoint only if direct).
C           1 ... Checkpoint.
C           2 ... Don't checkpoint.
C
C2IOp(37)
C     IOp(37) ... Frequency at which to do full Fock formation
C                 instead of incremental (L502, L508).
C          -1 ... Do not do incremental Fock formation.
C           0 ... Default (every 20 for direct, except 40 if CK is on).
C           N ... Every Nth cycle.
C
C2IOp(38)
C     IOp(38) ... Whether to vary integral cutoffs during direct SCF:
C           0 ... Default (5).
C           1 ... No.
C           2 ... Yes, do integrals 3 digits more accurately than
C                 current convergence.
C           3 ... Yes, do integrals at same accuracy as convergence
C                 until final iteration, then 2 digits more accurately.
C           4 ... Converge to 10**-5 with integrals good to 10**-6 first,
C                 then full convergence.
C           5 ... VarAcc allowed, decide based on details of problems.
C           6 ... VarAcc forbidden because of guess=read; allows different
C                 default actions for PBC.
C           7 ... Full accuracy for 2e part, but do pass 0 with cheaper XC
C                 grid.
C           8 ... Full grid throughout, but do pass 0 with cheaper integrals.
C
C2IOp(39)
C     IOp(39) ... On the fly symbolic matrix element generator.
C
C2IOp(40)
C     IOp(40) ... Use of reaction field; only used now for Onsager
C                 and control of details of SCIPCM
C          -N ... Multipoles of order N, increment field in Gen(2-4)
C           0 ... No.
C           N ... Multipoles of order N, store field in Gen(2-4)
C       00000 ... Default (same as 10000).
C       10000 ... Update surface every iteration.
C       20000 ... Update surface every iteration in pass 1 only.
C       30000 ... Update surface on pass 2 iterations only.
C       40000 ... Same as 3, but re-use 1e matrix instead of surface terms.
C       50000 ... Update surface and restart DIIS when within 10**-2 of
C                 convergence.
C
C2IOp(41)
C     IOp(41) ... Whether to converge on maximum density change as well
C                 or instead of RMS:
C           0 ... 22.
C           N ... Maximum allowed change is 10**N larger than RMS.
C          -1 ... Maximum allowed changed is same as RMS (i.e.,
C                 convergence only on maximum).
C          -2 ... Converge only on RMS density change.
C          N0 ... Converge on energy to 10**(N)*RMS-density-accuracy
C
C                 Also control of iterative diagonalization in L510.
C
C2IOp(42)
C     IOp(42) ... Number of orbitals to localize in L510.
C
C2IOp(43)
C
C     IOp(43): L509:  whether 5th order terms are treated explicitly
C           0 ... default: set to 1
C           1 ... all 5th order terms are treated implicitly
C           2 ... (debug option) 5th order GG and WG terms are
C                 explicitly computed in L715
C
C     IOp(43): L510:  DFT corrections to MCSCF.
C           0 ... No.
C           1 ... Compute energy using final MCSCF density.
C           2 ... Correct diagonal CI matrix elements using DFT.
C
C2IOp(45)
C     IOp(45) ... Some debugging option in L510.
C
C2IOp(47)
C
C     In L510, 1 to set up for CAS-MP2 or 2 to do spin-orbit calculation.
C
C2IOp(48)
C     IOp(48) ... Options to be passed to CalDFT:
C           N ... Control flag for CalDFT is N.
C
C                 Also contradictory use in L510.
C2IOp(49)
C     IOp(49) ... Use of sparse storage and Conjugate Gradient optimization
C                 instead of N**2 memory and diagonalization.
C           0 ... Default (11, or 22 if sparse is set in ILSW).
C           1 ... Diagonalization
C           2 ... Conjugate gradient.
C          10 ... Square storage (only in Fock formation if CG).
C          20 ... Linear storage (only in Fock formation if diagonalization).
C
C                 Also some option for Lanczos in L510.
C
C2IOp(52)
C     IOp(52) ... Amount of memory to allocate to stashing integrals (NYI).
C
C2IOp(53)
C     IOp(53) ... PCM input and solvent type.
C         N>0 ... Solvent type N, default parameters
C         N<0 ... Dielectric constant |N|/1000
C
C2IOp(55)
C     IOp(55) ... How many HOMOs and LUMOs to solve for after CG:
C           0 ... None.
C           N ... N of each.
C
C2IOp(56)
C     IOp(56) ... A0 for Onsager SCRF.
C           N ... N/1000 Bohr.
C
C2IOp(57)
C     IOp(57) ... First iteration at which to level shift and do FON.
C           0 ... Default: 1 unless doing stable=opt, then start after
C                 instability searches.
C
C2IOp(60-62)
C     IOp(60-62) ... Over-ride standard values of IRadAn, IRanWt, and IRanGd.
C
C2IOp(63)
C     IOp(63) ... Whether to do FMM.
C           0 ... Use global default.
C           1 ... Turn off FMM here regardless.
C         100 ... Turn off both FMM and FoFCou here.
C
C2IOp(64)
C     IOp(64) ... Over-ride default value of FMFlags
C           0 ... No.
C           N ... Yes, use N.
C
C2IOp(65)
C     IOp(65) ... Over-ride NFx parameter.
C           0 ... No.
C           N ... Yes, use N.
C
C2IOp(66)
C     IOp(66) ... Over-ride the choice of XC functional.
C           0 ... Use global values.
C           N ... Use functional N, with the same values as for IOp(74)
C                 in overlay 3.
C
C2L510_Junk
C     IOp(65-72) are used for various undocumented features in L510.
C
C2IOp(70)
C     IOp(70) ... Maximum Initial temperature for FON (non-PBC), or temperature
C                 for broadening (PBC and IOp(74)=[1-4]xx).
C          -2 ... None
C          -1 ... Start at a high temperature (limited only by DIIS error).
C           0 ... Default (3000K = 10 milliHartrees for non-PBC, 600K for PBC)
C           N ... N degrees
C
C2IOp(71)
C     IOp(71) ... Number of steps to apply simulated annealing (L502):
C           0 ... Default -- 10 steps FON / 20 steps PFON
C           N ... N steps
C
C2IOp(73)
C     IOp(73) ... Options for ADMP/QNDMS:
C           0 ... Default (2 for ADMP, 1 for QNDMS)
C           1 ... Use Lowdin basis for CP orthonormal transform.
C           2 ... Use Cholesky basis for CP orthonormal transform.
C
C2IOp(74)
C     IOp(74) ... Type of k-point integration:
C           0 ... Default (911, should be 193 for metals)
C           1 ... use LT method (interpolation)
C           2 ... occupy entire points (used together with broadening)
C           3 ... Full points for insulators, temperature broadening
C                 for metals.
C           9 ... occupy lowest NE at each k point regardless of the energies
C          10 ... improved LT with quadratic corrections
C          20 ... original LT method
C          90 ... no concern for corrections
C         100 ... smearing Marzari method I
C         200 ... smearing Marzari method II
C         300 ... first order Hermite-Gaussian of Paxton and Methfessel
C         400 ... Gaussian smearing
C         500 ... classical Fermi-Dirac broadening
C         900 ... no broadening (this will be Gaussian broadening with small T)
C
C2IOp(75-78)
C     IOp(75-78) ... Number of alpha electrons, alpha orbitals, beta electrons,
C                    and beta orbitals for fractional occupation.
C
C2IOp(79)
C     IOp(79) ... Range around Fermi level around which temperature distribution
C                 will be applied if broadening is turned on for PBC.
C           0 ... Default, a value will be chosen in ZInLT1.
C
C2IOp(80)
C     IOp(80) ... The maximum conjugate gradient step size
C          -1 ... No maximum step size
C           0 ... Default maximum (.8)
C        MMNN ... Step size of MM.NN
C
C2IOp(81)
C     IOp(81) ... Conjugate-Gradient Parameters
C          MM ... Maximum Number of CG cycles per SCF iteration.
C                 (defaults to 4 CG cycles).
C        NN00 ... Maximum Number of purification cycles per CG iteration.
C                 (defaults to 3 cycles).
C       00000 ... Don't use CG DIIS
C       10000 ... Use CG DIIS
C      000000 ... Polak-Ribiere CG minimization
C      100000 ... Fletcher-Reeves CG minimization
C     0000000 ... Use diagonal preconditioning in Conjugate-Gradient.
C     1000000 ... No preconditioning.
C
C2IOp(82)
C     IOp(82) ... C.G. Convergence criterion
C           0 ... Defaults to 10**(-7)
C           N ... 10**(-N)
C
C2IOp(83)
C     IOp(83) ... Maximum SCF DIIS vectors
C           0 ... Default (20, except 40 if CK is on).
C           N ... Use SCF DIIS with N vectors
C
C2IOp(84)
C     IOp(84) ... Restart in L509, Restart using Fock matrix in L502.
C           0 ... Default (Yes for PBC and sparse with guess Fock).
C           1 ... Yes.
C           2 ... No.
C
C2IOp(85)
C     IOp(85) ... Over-riding of maximum cycles for XQC.
C          -1 ... Default for first step (128).
C           0 ... No.
C           N ... Limit is N cycles.
C
C2IOp(86)
C     IOp(86) ... Strategy options (L502):
C      000000 ... Default (101100).
C           0 ... Default (1 except during Stable=Opt, then 4).
C           1 ... Just continue as usual if energy goes up.
C           2 ... Reduce DIIS space when energy rises from previous cycle.
C           3 ... Reduce DIIS space when energy goes above the lowest energy.
C           4 ... Reduce DIIS space whenever energy is above the lowest energy.
C          10 ... Turn on dynamic level shift from the beginning
C          20 ... Turn on dynamic level shift only after FON is over.
C         100 ... Keep level shift after energy rises.
C         200 ... Turn off level shift after energy rises.
C        1000 ... Level shift to a maximum of the Goal.
C        2000 ... Level shift to a maximum of 2*Goal.
C        3000 ... Level shift as much as necessary for HOMO>LUMO.
C        4000 ... Level shift only if the HOMO-LUMO gap is zero.
C        5000 ... Level shift only if the HOMO-LUMO gap is zero or
C                 insignificant (>-0.1)
C        6000 ... Level shift only if the HOMO-LUMO gap is zero or
C                 insignificant (>-0.1), up to twice the goal
C       N0000 ... No longer used.
C      100000 ... Turn off 3 and 4 point extrapolation if DIIS is on.
C      200000 ... Retain 3 and 4 point extrapolation if DIIS is on.
C
C     The energy is only checked after FON has been turned off.
C
C2IOp(87)
C     IOp(87) ... Accuracy criterion in Fock matrix formation:
C           0 ... Default, set in FoFCou/CalDSu based on accuracy
C                 part of IOp(5).  Typically 10^-10 for molecules and
C                 10^-12 for periodic systems.
C           N ... 10**-N.
C
C2IOp(88)
C     IOp(88) ... No longer used.
C
C2IOp(89)
C     IOp(89) ... Linearly dependent basis control for PBC; this and
C                 ZFormV should be moved to L302.
C
C2IOp(90)
C     IOp(90) ... Whether to generate sparse guess here.
C           1 ... Yes, do preliminary AM1 calculation.
C           2 ... Yes, do preliminary AM1 calculation and compare
C                 with guess from previous step in geometry optimization.
C
C2IOp(91)
C     IOp(91) ... Control option for chebyshev sparse control.
C
C2IOp(92)
C     IOp(92) ... How to do exact exchange:
C           0 ... Default (Normal processing based on FMM for non-PBC,
C                 separate Coulomb and NFx exchange for PBC).
C           1 ... FoFCou for Coulomb, separate FoFCou/NFx for exchange.
C
C2IOp(93)
C     IOp(93) ... Number of initial iterations for which damping is allowed:
C           0 ... Default (10).
C           N ... N iterations.
C
C2IOp(94)
C     IOp(94) ... Another undocumented option in L510 (step size in
C                 trajectories?)
C
C2IOp(95)
C     IOp(95) ... Whether to do Davidson during RFO in 510:
C           0 ... No.
C           1 ... Yes.
C
C2IOp(96)
C     IOp(96) ... Over-ride IRadAn for CPHF-like step in L509, and
C                 for pass 0 grid in L502.
C           0 ... Use default,
C           N ... Use grid N.
C
C                 RAS control in L510.
C
C2IOp(97)
C     IOp(97) ... Whether to update precomputed grid data with timing
C                 information.
C           0 ... Default (Yes, if available).
C           1 ... Yes.
C           2 ... No.
C
C                 Hopping threshold during trajectories with L510.
C
C2IOp(98)
C     IOp(98) ... Whether to save eigenvalues and orbitals at all k-points.
C           0 ... Default (Yes).
C           1 ... Yes.
C           2 ... No.
C
C                 Some undocumented surface hopping option in L510.
C
C2IOp(99)
C     IOp(99) ... Grid for numerical k-integration in FT-LT method.
C           0 ... Default:  32,12,8 for 1,2,3d
C
C2IOp(100)
C     IOp(100) ... Tight convergence during CGDMS:
C            0 ... Default (No).
C            1 ... Yes.
C            2 ... No.
C
C2IOp(101)
C     IOp(101) ... SDif test on numerical accuracy of PBC diagonalization.
C            0 ... Default (10)
C           -1 ... No test
C          N>0 ... Abort of SDif is larger than N.
C
C2IOp(102)
C     IOp(102) ... Maximum number of configurations for CAS-MP2:
C            0 ... Default (1000).
C            N ... N.
C
C2IOp(103)
C     IOp(103) ... Number of occupied and virtual orbitals to print for
C                  each k-point.
C           -1 ... Default of 5 occupieds and 5 virtuals.
C            0 ... Default is 5 if printing turned up; otherwise 0.
C            N ... N occupieds and N virtuals.
C
C2IOp(104)
C     IOp(104) ... GF testing
C            0 ... No.
C            1 ... Yes.
C            2 ... Yes, and save biorthogonalized orbitals over canonical MOs.
C
C2IOp(105)
C     IOp(105) ... Convergence for Davidson iterations in L510.
C            0 ... Default (6 digits on coefficients).
C            N ... 10^-N on coefficients.
C
C2IOp(106)
C     IOp(106) ... Saving and restart for iterative CI in L510:
C            0 ... No.
C            1 ... Save CI vectors.
C            2 ... Restart CI, possibly adding states.
C
C2IOp(107)
C     IOp(107) ... Maximum number of Davidson iterations in L510:
C            0 ... Default (huge, number of CI configurations).
C           -1 ... No limit
C            N ... N iterations.
C
C2IOp(108)
C     IOp(108) ... Minimum number of iterations at which to damp density.
C            0 ... Default (1 if transition metals present, otherwise 0).
C            N ... N.
C
C2IOp(109)
C     IOp(109) ... Whether to rotate final spin of GHF wavefunctions.
C            0 ... Default (No).
C            1 ... Yes, rotate to align <S> with z.
C            2 ... No.
C
C2IOp(110)
C     IOp(110) ... Debugging options for StuCyc:
C          Bit     Action
C          0-1     Algorithm in StuFo3 0=DB ints 1=SB spinorb ints
C                  2=square reg ints (nyi) 3=lt reg ints (nyi)
C            2  Conjugate computed G(P)* before use.
C            3  Conjugate computed F* before use.
C            4  Compute variations on the energy traces.
C            5  Whether to do CPHF using simple iterations or
C               an in-memory version of lineq1.
C            6  Damping in CPHF.
C            7  3/4 point extrapolation in CPHF.
C            8  permute permutations in CPHF.
C
C2IOp(111)
C     IOp(111) ... QNDMS/SimOpt control:
C            1 ... First-order SimOpt
C           11 ... Extrapolate Fock matrices using sim-DIIS coefficients.
C                  (default if SimOpt).
C           21 ... Do not extrapolate Fock matrices.
C          111 ... (Full matrix) QNDMS iterations here.
C
C2IOp(120)
C     IOp(120) ... Store nuclear repulsion energy as total energy?
C            0 ... Default (No)
C            1 ... Yes.
C
C2IOp(121)
C     IOp(121) ... IDoVI for HarFok in test calculations.
C          xx0 ... Ones digit always set to 4 here.
C
C2IOp(122)
C     IOp(122) ... Whether to do Hirshfeld analysis of spin orientations.
C            0 ... Default (yes for GHF/GKS, no otherwise).
C            1 ... Yes.
C            2 ... No.
C
C2IOp(123)
C     IOp(123) ... Number of iterations for Pseudodiagonalization.
C            0 ... Default (1 for semi-empirical, -1 for HF/DFT).
C           -1 ... Sweep until variable convergence is reached.
C            N ... Do a maximum of N sweeps.
C
C2IOp(124)
C     IOp(124) ... Variable convergence in pseudodiagonalization:
C            0 ... Default (1 for ab initio, 2 for semi-empirical).
C            N ... Off-diagonals larger than the initial max/10^N
C                  are swept.
C           -N ... Off-diagonals larger than OVMax*10^-N are swept,
C                  with OVMax updated each iteration.
C
C2IOp(125)
C     IOp(125) ... Scale factor for Diag/PseudoDiag tradeoff.  Roughly
C                  related to the ratio of Diag to Fock formation for
C                  large systems.
C           -1 ... No test.
C            0 ... Default (30 for ab intio, 15 for semi-empirical).
C            N ... Pseudodiag must be estimated to be at least N/10 times
C                  faster than full diag to be used.
C
C2IOp(126)
C     IOp(126) ... Initial step size for linear searches in L508:
C            0 ... Default (300).
C            N ... N/10000.
C
C2IOp(127)
C     IOp(127) ... Maximum rotation gradient for Newton-Raphson in L508
C                  (above this value, scaled steepest descent is used):
C            0 ... Default (1.d-2).
C            N ... 10**-N.
C
C2IOp(128)
C     IOp(128) ... Diagonalization algorithm (L502).
C            0 ... Default (GSYEVD if memory permits, otherwise GSPEV).
C            N ... Algorithm N in DiagDS.
C
C2IOp(129)
C     IOp(129) ... Threshold for trying alternate initial guess:
C            0 ... Default (2).
C           -1 ... Same as 0.
C           -2 ... Ignore the alternate guess.
C            N ... Try the alternate if the DIIS error >=10^-N.
C?
#include "commonb2.inc"
      Parameter (NSkew=7,LRwGen=1000,NSCF=5,ISEGes=1,ISEEn=2,IGesEn=3,
     $  IVac=4,NProgC=5,NFFDef=20,NDIDef=20,MxCDef=128,MxCXQC=64,
     $  MaxFR=4,MaxLvW=2,MaxDVS=4,IRwGen=501,IRwIBF=508,IRwS=514,
     $  IRwH=515,IRwT=516,IOFInF=521,IRwEig=522,IRwSym=523,IRwCA=524,
     $  IRwCB=526,IRwPA=528,IRwPB=530,IRwPT=532,IRwPS=534,IRwFA=536,
     $  IRwFB=538,IRwUR=540,IRwNCA=554,IRwNCB=556,NEqShl=565,IRwW=571,
     $  IOSBlk=577,NEq=580,IRwFx=584,IOSymA=606,IOEDrv=619,IRwSM=672,
     $  IRwRAO=685,IRwV=686,IRwPOA=688,IRwPOB=699,IOFDen=702,IODFit=723,
     $  NEqSh2=726,IRwAGs=729,IOXCGd=742,IRwFRA=743,IRwFRB=744,
     $  IRwECK=746,IRwSO=751,IRwSMO=769,IRAOFr=785,MxGWIn=32,IPOA1=2001,
     $  IPOB1=2002,IPDDef=-30,MxXCNm=DEFMAXXCNAME,MOmega=DEFMOMEGA,
     $  NRestI=200,NOmega=DEFNOMEGA,IHFX=1,IHFJ=6,ILFrag=67,NRestR=200)
      Dimension V(*), ScaDFX(4), DGen(LRwGen), IProgC(NProgC), NTot(3),
     $  DGenC(LRwGen), Omega(MOmega,NOmega), SVec(3), EField(0:34),
     $  XXA(1), JJA(1), GFac(0:24), MulTab(8,8), CelVec(3,3),
     $  NS(9), NS2(9), NS4(9), NBF(1), NBFU(1), IRData(NRestI),
     $  RRData(NRestR)
      Integer IBfSrt(10), FilNum, VarAcc, NOccSy(2), FMFlag,
     $  FMFlg1, VarAc1, DirSc0, FixRis, IGWInf(MxGWIn), ChrTbl(8,0:8)
      Logical ISort,DIIS,IRSort,ROHF,ZDO,DirSCF,RHF,MemInt,FrcJEn,
     $  IncFoc,SaveMO,SvSort,Sleazy,DFT,AllowP(50),Increm,PChek4,AM1,
     $  FMM,RComp,JustJ,AbDiag,SpFock,SpDiag,RdPCM,SCIPCM,UsSpar,PCM,
     $  Mixed,DoDFJ,DoPrnt,DoPrn1,JustEn,DoSCF(NSCF),ZDOUse,OK,SymmP,
     $  Onsagr,UsKPnt,FkGues,OKSym,HaveTM,UseGEW,UseLow,LBit,FrcNuc,
     $  UsePBC,BndOnl,UpdTim,SavPBC,FixT,Halt,If0K,TigPCV,CBfn,Spinor,
     $  No2E,DoCycO,DoUpd,CFPCM,SEHam,DFTB,DoRotS,DoStCy,TstSpI,DoCyc,
     $  DoStCP,DoNFmn,DoNDXC,DoSim,DoSimP,DoDMS,UseNuc,AM1Ham,NoInCo,
     $  JONEKO,FrcSDg,DISPCM,DFTBHm,HarTHm,HarTst,DoSPI,DoHirS,DefPsu,
     $  Error
      Integer MaxIOp, IOp
      Parameter (MaxIOp=DEFMAXIOP)
      Common /IOp/ IOp(MaxIOp)
#include "commonmol.inc"
#include "repall.inc"
      Common /Sym1/ NBlksS,MaxCmp,NBFS(8),NBFUS(8),NamRep(34)
      Common/IO/In,IOut,IPunch
      Common/IBF/IBF(30)
      Common/MUnit/IUnit(20)
      Common/PhyCon/PhyCon(30)
      Character MethNm*80, NameFr(MaxFR)*58, OnOff(2)*3, SCFNam*80,
     $  ExName*(MxXCNm), CoName*(MxXCNm)
      Equivalence (IFil,IBF(24)), (NWIIB,IBF(29)), (ISym2e,IBF(30)),
     $  (ISMode,IBF(1)),(Hartre,PhyCon(8)),(Boltz,PhyCon(10))
      Save Zero,One,Two,SlzAcc,JJA,XXA,IRwIBS,IOSort,Eight,IDefT,DefVSG,
     $  NameFR,OnOff,Pt1,Big,DGAccP,DGAccE,Pt5,IDefTP,IDefSD,DGen,
     $  IBfSrt,DGenC,NFrag,SmalCK,ErPMax
      Data SlzAcc/1.0d-4/,JJA/0/,XXA/0.0d0/,IOSort/610/,Zero/0.0d0/,
     $  One/1.0d0/,Two/2.0d0/,IRwIBS/608/,Eight/8.0d0/,IDefT/3000/,
     $  DefVSG/0.2d0/,OnOff/'on','off'/,Pt1/0.1d0/,Big/1.d10/,
     $  DGAccP/5.0d-7/,DGAccE/5.0d-7/,Pt5/0.5d0/,IDefTP/600/,IDefSD/10/,
     $  DGen/LRwGen*0.0d0/,IBfSrt/10*0/,DGenC/LRwGen*0.0d0/,NFrag/0/,
     $  SmalCK/1.0d-10/,ErPMax/1.d-3/
      Data NameFr/'No special actions if energy rises.',
     $  'Reduce DIIS space if energy rises from previous iteration.',
     $  'Reduce DIIS space if new energy goes above lowest so far.',
     $  'Reduce DIIS space if new energy is above lowest so far.'/
 1010 Format(' Error in input density:  TNEA=',F12.4,' NAE=',I6,
     $  ' TNEB=',F12.4,' NBE=',I6,'.')
 1070 Format(1X,A,' SCF:')
 1110 Format(' SCF Done:  E(',A,') = ',G19.12,12H A.U. after ,I4,
     $  ' cycles',/,'            NFock=',I3,'  Conv=',D8.2,'     -V/T=',
     $  F7.4)
 1120 Format(' <Sx>=',F7.4,' <Sy>=',F7.4,' <Sz>=',F7.4,' <S**2>=',F7.4,
     $  ' S=',F7.4,/,' <L.S>=',G19.12)
 1160 Format(' Requested convergence on RMS density matrix=',1PD8.2,
     $' within',i4,' cycles.')
 1170 Format(' Requested convergence on MAX density matrix=',1PD8.2,
     $' within',i4,' cycles.')
 1180 Format(' Requested convergence on MAX density matrix=',1PD8.2,'.')
 1190 Format(' Requested convergence on             energy=',1PD8.2,'.')
 1200 Format(' Virtual orbitals will be shifted by ',F7.3,' hartree.')
 1210 Format(' NBasis=',I9,' NAllTT=',I9,' NAllSq=',I9,' NWIIB= ',I9,/,
     $       ' INEq=  ',I9,' LenNEq=',I9,' IIJ=   ',I9,' ID=    ',I9,/,
     $       ' IDD=   ',I9,' IDA=   ',I9,' IDB=   ',I9,' IF=    ',I9,/,
     $       ' IFF=   ',I9,' IFA=   ',I9,' IFB=   ',I9,' IScr2= ',I9,/,
     $       ' IV=    ',I9,' IBuf2E=',I9,' IMin1= ',I9,' IMin=  ',I9,/,
     $       ' IV=    ',I9,' IOpCl= ',I9)
 1220 Format(' Unsorted integral processing.')
 1230 Format(' Sorted integral processing.')
 1250 Format(' IVT=',I12,' IEndB=',I12,' NGot=',I12,' MDV=',I12,/,
     $       ' LenX=',I12,' LenY=',I12)
 1260 Format(' Using DIIS extrapolation, IDIIS=',I6,'.')
 1270 Format(' Overlap will be assumed to be unity.')
 1280 Format(' A Direct SCF calculation will be performed.')
 1290 Format(' Two-electron integral symmetry not used.')
 1300 Format(' Two-electron integral symmetry used by symmetrizing',
     $       ' Fock matrices.')
 1310 Format(' Two-electron integrals replicated using symmetry.')
 1320 Format(' Two-electron integrals will be kept in memory, NGetB=',
     $       I12,'.')
 1330 Format(' Fock matrices will be formed incrementally for',I4,
     $  ' cycles.')
 1340 Format(' KE=',1PD19.12,' PE=',1PD19.12,' EE=',1PD19.12)
 1350 Format(' Failed to select integral storage method, IntTyp=',I10,
     $       '.')
 1390 Format(' Density fitting will be used, IDenFit=',I2,' NucTyp=',I3,
     $  '.')
 1430 Format(' Warning!  Cutoffs for single-point calculations used.')
 1440 Format(' Variable integral accuracy in use.')
 1450 Format(' Integral accuracy reduced to',1PD8.1,
     $  ' until final iterations.')
 1460 Format(' Illegal IOp(38)=',I12,'.')
 1480 Format(' Unrecognized value for IOp(41)=',I12,'.')
 1490 Format(' Integral symmetry usage will be decided dynamically.')
 1500 Format(' Density matrix is not symmetric, so final MOs not',
     $  ' symmetrized.')
 1540 Format(' Level shift goal',F7.3,' maximum shift',F7.3,'.')
 1541 Format(' Level shift goal',F7.3,' unlimited maximum shift.')
 1542 Format(' Level shift goal',F7.3,
     $  ', shift only if positive HOMO-LUMO gap.')
 1543 Format(' Level shift goal',F7.3,' maximum shift',F7.3,
     $  ', only if HOMO-LUMO gap >',F8.3,'.')
 1550 Format(' The wavefunction is already stable.')
 1560 Format(1X,A)
 1570 Format(' Dynamic level shift is ',A,' during FON iterations.')
 1580 Format(' Dynamic level shift is ',A,' after energy rises.')
 1590 Format(' Integral accuracy and convergence thresholds reduced,',
     $  ' IAccFD=',I2,'.')
 1600 Format(' Converting RHF to UHF.')
 1610 Format(' NCelOv=',I6,' NClRep=',I6,' NMtPBC=',I6,' NRecip=',I6,
     $  ' NKPnt=',I8,'.')
 1620 Format(1X,I12,' words used for storage of precomputed grid.')
 1630 Format(' Selected VarAcc=',I2,'.')
 1640 Format(' Forcing IOpCl=',I2,'.')
 1650 Format(' Clearing in-memory integrals.')
C
      Call Drum(V,NGot)
      IIRDat = 1
      IRRDat = IIRDat + InToWP(NRestI)
      NRest = IRRDat + NRestR - 1
      Call TstCor(NRest,NGot,'UHFOpn-Rest')
      Call AClear(NRestR,RRData)
      Call IClear(NRestI,IRData)
      Call SetILR(IOut,.False.,ILRest,NRestI,IRData,NRestR,RRData,NRest,
     $  V(IIRDat))
      IFMM = IOp(63)
      IPrint = IOp(33)
      ISEAll = NRest + 1
      Call SetPFl(IOut,IPrint,IFMM,IPFlag,AllowP,FMM,FMFlag,FMFlg1,
     $  NFxFlg,IHMeth,Omega,ISEAll,LSEAll,jStrt,V,NGot)
      Call GetPBC(NPDir,CelVec)
      Call ILSW(2,ILFrag,IFrgTp)
      If(IFrgTp.ne.0) Call GauErr('IFrgTp in UHFOpn.')
      jIPBC = jStrt
      Call RdIPBC(1,NPDir,0,jIPBC,jStrt,V,NGot)
      LIPBC = (jStrt-jIPBC)*IntPWP(0)
      LenV = NGot - jStrt + 1
      NPDir = 0
      NMtPBC = 0
      NCelOv = 0
      NRecip = 0
      NKPnt = 0
      Call ObtCel(2,V(jIPBC),1,NPDir)
      Call ObtCel(2,V(jIPBC),2,NMtPBC)
      Call ObtCel(2,V(jIPBC),3,NCelOv)
      Call ObtCel(2,V(jIPBC),10,NRecip)
      Call ObtCel(2,V(jIPBC),11,NKPnt)
      UsePBC = NPDir.gt.0
      UsKpnt = UsePBC.and.(NRecip.gt.0.or.NRecip.lt.-1)
      Call GetNB6(NBas6D)
      Call LClear(NSCF,DoSCF)
      DoSCF(NSCF) = .True.
      IRWIBS = FilNum(IRWIBS,IUnit(13))
      Call ILSW(2,1,IOpCl)
      Call ILSW(2,22,IROHF)
      Call ILSW(2,6,IStab)
      If(IStab.eq.1) then
        Write(IOut,1550)
        Goto 999
        endIf
      Call RetDFT(IOp(60),IOp(61),IOp(62),DFT,IExCor,IXCFnc,IExCoW,
     $  IRadAn,ScaHFX,ScaDFX,IRanWt,IRanGd,ICorTp,IEmpDi)
      INExCo = IOp(66)
      If(INExCo.ne.0) then
        IOp79 = IRanWt
        IOp75 = IRadAn
        Call GetDFT(IOut,.False.,IHMeth,IOp79,INExCo,IOp75,0,0,0,0,0,0,
     $    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,IExCor,IXCFnc,IExCoW,
     $    IRadAn,DFT,ExName,CoName,ScaHFX,ScaDFX,IRanWt,IRanGd,ICorTp,
     $    IEmpDi)
        endIf
      ISleaz = Mod(IOp(32),10)
      IAprXC = Mod(IOp(32),10000) / 1000
      Call CkDFit(IHMeth,.True.,IDenFit,NucTyp,NDBF6D)
      HarTst = HarTHm(IHMeth)
      If(HarTst) IAprXC = 1
      If(ISleaz.eq.0) ISleaz = 2
      Sleazy = ISleaz.eq.1.or.ISleaz.eq.6
      MaxCyc = IOp(7)
      FrcJEn = MaxCyc.eq.-1
      Halt = IOp(13).ne.1
      IDamp = IOp(18)
      IDampC = IOp(21)
      MxCAlt = IOp(85)
      If(FrcJEn) then
        MaxCyc = 1
      else if(MxCAlt.eq.-1) then
        MaxCyc = MxCDef
      else if(MxCAlt.gt.0) then
        MaxCyc = MxCAlt
      else if(MaxCyc.eq.0.and.Halt) then
        MaxCyc = MxCDef
      else if(MaxCyc.eq.0) then
        MaxCyc = MxCXQC
      else if(.not.Halt) then
        MaxCyc = Min(MaxCyc,MxCDef)
        endIf
      IStrat = IOp(86)
      FixRis = Mod(IStrat,10)
      If(FixRis.le.0) then
        If(IStab.eq.2) then
          FixRis = 4
        else
          FixRis = 1
          endIf
        endIf
      FixRis = Min(FixRis,MaxFR)
      LvWFON = Mod(IStrat,100) / 10
      If(LvWFON.le.0) LvWFON = 1
      LvWFON = Min(LvWFON,MaxLvW)
      LvRise = Mod(IStrat,1000) / 100
      IDyVSh = Mod(IStrat,10000) / 1000
      If(IDyVSh.le.0) IDyVSh = 1
      IDyVSh = Min(IDyVSh,MaxDVS)
      MaxDamp = IOp(93)
      MinDamp = IOp(108)
      IRotS = IOp(109)
      DoRotS = IRotS.eq.1.and.IOpCl.eq.6
      IStuFl = IOp(110)
      I34DIS = Mod(IStrat,1000000) / 100000
      IAccFD = IOp(87)
      IOpClN = IOp(27)
      If(IOpClN.ne.0) then
        IOpCl = IOpClN
        Write(IOut,1640) IOpCl
        Call ILSW(1,1,IOpCl)
        endIf
      If(IROHF.eq.1) IOpCl = 1
      Call DecoSC(IOpCl,NSpBlk,NRI,NDimBl,CBfn,NSpBlX,Spinor)
      IOpClX = IOpCl - Mod(IOpCl,2)
      NBD = NBasis*NDimBl
      NBasUD = NBsUse*NDimBl
C     RHF is really whether to diagonalize one Fock matrix rather than 2
      RHF = NSpBlk.eq.1.and.IROHF.le.0
      ROHF = (IOpCl.le.1.and.IROHF.eq.1).or.RHF
      If(IDamp.eq.0) then
        If(IStab.eq.2.or.IHMeth.eq.14) then
          IDamp = -2
        else
          IDamp = -1
          endIf
        endIf
      If(IDamp.eq.-3) then
        If(IDampC.eq.0) then
          DampCut = GFloat(1) / GFloat(3)
        else
          DampCut = GFloat(IDampC) / GFloat(10000)
          endIf
      else
        If(IDampC.eq.0) IDampC = 3
        DampCut = GFloat(10)**(-IDampC)
        endIf
      DirSc0 = Mod(IOp(5),10)
C     No conventional Fock formation for complex here.
      If((IOpCl.eq.2.or.IOpCl.eq.3).and.DirSc0.le.1) DirSc0 = 2
      DirSCF = DirSc0.ge.2
      DoPrnt = IPSave(0).eq.0.or.IPrint.gt.0
      DoPrn1 = IPSave(0).le.1.or.IPrint.gt.0
      If(DirSCF.and.IPrint.gt.0) Write(IOut,1280)
      ICntrl = Mod(IOp(5),100000) / 10
      IIFock = Mod(IOp(5),1000000) / 100000
      If(HarTst) IIFock = 1
      No2E = (Mod(IOp(5),10000000)/1000000).eq.1
      IVShft = IOp(10)
      IExtp = Mod(IOp(11),10)
      ICampK = Mod(IOp(11),100)/10
      If(ICampK.eq.2) ICampK = 0
      ICKThr = Mod(IOp(11),100000)/100
      If(ICampK.eq.1.or.ICampK.eq.4) then
        If(ICKThr.eq.0) ICKThr = 300
        CKThr = GFloat(ICKThr) / GFloat(1000)
      else if(ICampK.eq.3) then
        If(ICKThr.eq.0) then
          CKThr = Pt1
        else if(ICKThr.eq.999) then
          CKThr = SmalCK
        else
          CKThr = GFloat(ICKThr) / GFloat(1000)
          endIf
      else
        CKThr = Zero
        endIf
      IUHFT = Mod(IOp(14),10)
      JustJ = (Mod(IOp(14),100)/10).eq.1
      BndOnl = (Mod(IOp(14),100)/10).eq.2
      If(JustJ.or.BndOnl) IPrint = Max(IPrint,1)
      IADMP = Mod(IOp(14),1000) / 100
      IGEWD = Mod(IOp(14),10000) / 1000
      UseGEW = (IADMP.gt.0.and..not.UsKPnt.and.IGEWD.eq.0).or.IGEWD.eq.1
      IFitEn = Mod(IOp(14),100000) / 10000
      IDoGrd = Mod(IOp(14),1000000) / 100000
      IUpd = Mod(IOp(14),10000000)/1000000
      If(IUpd.eq.0) IUpd = 2
      IStuCP = Mod(IOp(14),100000000)/10000000
      DoUpd = IUpd.eq.2.or.IUpd.eq.4
      DoCycO = DoUpd.or.UsePBC
      DoStCy = IUpd.eq.3
      DoCyc = DoCycO.or.DoStCy
      TstSpI = IUpd.eq.4
      DoStCP = IStuCP.eq.1.and.(DoStCy.or.TstSpI)
      INumD = Mod(IOp(14),1000000000)/100000000
      DoNFmn = INumD.eq.1
      If(DoNFmn) Call GauErr('DoNFmn NYI.')
      DoNDXC = INumD.eq.2
      IADBas = Mod(IOp(73),10)
      ISyOcc = Mod(IOp(15),10)
      If(IFrgTp.gt.0) then
        ISyOcc = 2
      else if(ISyOcc.eq.0) then
        ISyOcc = 1
        endIf
      Call GetNOp(NOp1,NOp2)
      Call PcckI(0,Junk,NEqAll,MxAtSO,NOpAll,NAtoms,NOpAll)
      Call ILSW(2,26,NoSymm)
      IAbDag = Mod(IOp(15),100)/10
C     There is some bug with use of symmetry in diagonalization for GHF.
      FrcSDg = IAbDag.eq.3.or.IFrgTp.gt.0
      AbDiag = (NoSymm.eq.0.and.IStab.ne.2.and.(IAbDag.eq.1.or.
     $  (IAbDag.eq.0.and.NDimBl.eq.1.and.(NOp2.gt.1.or.ISyOcc.ge.2))))
     $  .or.FrcSDg
      ISolve = IOp(16)
      IDiag = Mod(ISolve,100)
      Call ILSW(2,40,ISpCut)
      ISpFok = Mod(ISolve,1000)/100
      ZDO = IOp(35).eq.1
      DFTB = DFTBHm(IHMeth)
      DefPsu = DFTB.or.(IHMeth.ne.5.and.SEHam(IHMeth))
      If(IDiag.eq.0) then
        If(UsSpar(0)) then
          IDiag = 4
        else if(DefPsu) then
          IDiag = IPDDef
        else
          IDiag = 1
          endIf
      else if(IDiag.eq.-1.or.IDiag.eq.3) then
        IDiag = -(MaxCyc+10)
        endIf
      NPassI = IOp(123)
      If(NPassI.eq.0) then
        If(DefPsu) then
          NPassI = 1
        else
          NPassI = -1
          endIf
        endIf
      ISwCut = IOp(124)
      If(ISwCut.eq.0) then
        If(DefPsu) then
          ISwCut = 2
        else
          ISwCut = 1
          endIf
        endIf
      IPsCut = IOp(125)
      If(IPsCut.eq.0) then
        If(DefPsu) then
          IPsCut = 15
        else
          IPsCut = 30
          endIf
        endIf
      IAlDNI = IOp(128)
      FixT = (IDiag.eq.4.or.IDiag.eq.8).and.ISpCut.eq.0
      SpDiag = (IDiag.ge.4.and.IDiag.le.8).or.UsSpar(0)
      CnvDef = SCFAcc(0)
      If(UsePBC) CnvDef = CnvDef*GFloat(10)
      If(IOp(6).ne.0) then
        Acurcy = SCFAcc(IOp(6))
      else if(Sleazy) then
        Acurcy = SlzAcc
      else if(ISleaz.eq.3) then
        Acurcy = DGAccP
      else
        Acurcy = CnvDef
        endIf
      IAccMx = IOp(41)
      Mixed = DFT.and.(ScaHFX.ne.Zero.or..not.DoDFJ(0))
      If(IHMeth.eq.14.and.Mixed)
     $  Call GauErr('Int=HCXJ and Mixed in UHFOpn.')
      If(ISpFok.eq.4) then
        If(ZDO) then
          ISpFok = 2
        else
          ISpFok = 1
          endIf
        endIf
      SpFock = ISpFok.eq.2.or.(ISpFok.eq.0.and.SpDiag.and.
     $  (DFT.and..not.Mixed))
      If(IVShft.eq.0.and.(SpDiag.or.UsKpnt)) IVShft = -1
      JSym2E = Mod(IOp(17),100)/10
      If(JSym2E.eq.4) JSym2E = 2
      IFrcSy = Mod(IOp(17),1000)/100
      If(IStab.eq.2.or.FrcSDg) IFrcSy = 0
      SymmP = ((Mod(IOp(17),10000)/1000).le.1).and..not.FrcSDg
      IDIIS = IOp(22)
      IDIIS0 = Mod(IDIIS,10)
      IDIS14 = Mod(IDIIS,100000)/10
      IDIIS5 = Mod(IDIIS,1000000)/100000
      NOnlyS = Mod(IDIIS,10000000)/1000000
      If(IDIS14.eq.0) then
        IDIIS = 1040 + IDIIS0
      else
        IDIIS = IDIS14*10 + IDIIS0
        endIf
      IDIIS0 = Mod(IDIIS,10)
      IDIIS1 = Mod(IDIIS,100)/10
      If(IStab.eq.2) then
        IDIIS0 = 2
      else if(IDIIS0.eq.0) then
        If(SpDiag.or.IDIIS1.gt.1) then
          IDIIS0 = 2
        else
          IDIIS0 = 3
          endIf
        endIf
      DIIS = IDIIS0.ne.1
      IFBrod = IDIIS0 - 3
      If(DIIS.and.I34DIS.le.1) IExtp = 2
      DISPCM = IDIIS5.eq.1
      IDIIS = (IDIIS/10)*10
      If(.not.DISPCM) NOnlyS = 0
      If(NOnlyS.gt.0) MaxCyc = MaxCyc*NOnlyS
      IRaf = IOp(29)
      IAbSym = IOp(30)
      If(IAbSym.eq.0) then
        If(NoSymm.eq.0.and.NDimBl.eq.1) then
          IAbSym = 1
        else
          IAbSym = 2
          endIf
        endIf
      SaveMO = IOp(36).eq.1.or.(IOp(36).eq.0.and.DirSCF)
      NFull = IOp(37)
      NFullF = NFull
      If(NFullF.eq.-1) NFullF = 2
      If(NFullF.eq.0) then
        If(ICampK.ne.0) then
          NFullF = Min(2*NFFDef,MaxCyc)
        else
          NFullF = Min(NFFDef,MaxCyc)
          endIf
        endIf
      Call ILSW(2,4,ISmStp)
      NTT = (NBasis*(NBasis+1))/2
      NTTSp = (NBD*(NBD+1))/2
      NTTX = NRI*NTTSp*NMtPBC
      IFkGes = IOp(84)
      If(IFkGes.eq.0) then
        If(NMtPBC.gt.1.and.IADMP.ne.1) then
          IFkGes = 1
        else
          IFkGes = 2
          endIf
        endIf
      FkGues = ITqry(IRwFA).eq.NTTX.and.IFkGes.eq.1
      If(.not.ROHF) FkGues = FkGues.and.ITqry(IRwFB).eq.NTTX
      VarAcc = IOp(38)
      If(IStab.eq.2.or.HarTst) then
        VarAcc = 1
      else if(VarAcc.eq.0) then
        If(FkGues) then
          VarAcc = 4
        else
          VarAcc = 1
          endIf
        endIf
      IThAlG = IOp(129)
      If(IThAlG.eq.0.or.IThAlG.eq.-1) IThAlG = 2
      If(IThAlG.eq.-2) then
        AlGsMn = -One
      else
        AlGsMn = GFloat(10)**(-IThAlG)
        endIf
      Call ILSW(1,4,0)
      Call GetPCM(IOut,IPrint,.False.,0,1,ISCRF,RJunk)
      ISCRF0 = Mod(ISCRF,10000)/1000
      ISCRFT = Mod(ISCRF,1000)/100
      If(PCM(ISCRF)) then
        Call PCMIOp('load',' ',IJunk,RJunk)
        Call PCMPut('ICalTy',1,RJunk)
        Call PCMGet('NTs',NTs,RJunk)
      else
        NTs = 0
        endIf
      ItSCRF = IAbs(IOp(40))/10000
      Increm = IOp(40).lt.0
      IFOSim = IOp(111)
      DoSim = Mod(IFOSim,10).ge.1
      ISimP = Mod(IFOSim,100)/10
      If(DoSim.and.ISimP.eq.0) ISimP = 1
      DoSimP = DoSim.and.ISimP.eq.1
      DoDMS = (Mod(IFOSim,1000)/100).ge.1
      Call TstVAC(IOut,IPrint,UsePBC,NAtoms,IAn,IAtTyp,AtmChg,ItSCRF,
     $  ISmStp,DoSim,NShell,ShellN,ShellA,EXX,DirSCF,ISCRF0,DFT,Mixed,
     $  IDenFit,IRadAn,IDIIS,VarAcc,Sleazy,IVShft,FkGues,NExtDmp,HaveTM,
     $  MaxDamp,DampLim,LvRise,SpDiag,IHMeth,MinDamp,DoDMS)
      If(IVShft.eq.-2) then
        VSGoal = DefVSG
      else
        VSGoal = Abs(GFloat(IVShft)/GFloat(1000))
        endIf
      If(IDyVSh.eq.2.or.IDyVSh.eq.6) then
        VSMax = VSGoal + VSGoal
      else if(IDyVSh.eq.3) then
        VSMax = Big
      else
        VSMax = VSGoal
        endIf
      If(IDyVSh.eq.4) then
        VSGpMn = Zero
      else if(IDyVSh.eq.5) then
        VSGpMn = -Pt1
      else
        VSGpMn = -Big
        endIf
      ISvInt = IAbs(IOp(19)) / 10
      IntTyp = Mod(IOp(19),10)
      ICtDFT = IOp(48)
      ICtD01 = Mod(ICtDFT,100)
      If(ICtD01.eq.0) ICtDFT = ICtDFT + 11
      ICtDF5 = Mod(ICtDFT,1000000) / 100000
      If(ICtDF5.eq.0.and.IAprXC.eq.0) then
        If(IStab.eq.2.or.IADMP.eq.1.or.DoDMS.or.DoSim) then
          ICtDFT = ICtDFT + 100000
        else
          ICtDFT = ICtDFT + 400000
          endIf
        endIf
      ICtDF6 = Mod(ICtDFT,10000000) / 1000000
      If(ICtDF6.eq.0.and.IAprXC.ne.0) ICtDFT = ICtDFT + 3000000
      If(SpFock.or.SpDiag) then
        AbDiag = .False.
        IAbSym = 2
        endIf
      NAtomP = NAtoms + NPDir
      If(IDoGrd.eq.1) then
        LFXYZ = 3*(NAtomP+1)
        ICtDFT = ICtDFT + 100
        IntTyp = 1
      else
        LFXYZ = 3
        endIf
      IPCM = IOp(53)
      IA0 = IOp(56)
      LvFirst = IOp(57)
      If(IStab.eq.2) then
        If(LvFirst.eq.0) LvFirst = 4
      else
        If(LvFirst.eq.0) LvFirst = 1
        endIf
      If(IRadAn.eq.0) then
        If(Sleazy) then
          If(ISleaz.eq.1) then
            IRadAn = 3
          else
            IRadAn = 2
            endIf
        else if(ISleaz.eq.6) then
          IRadAn = 3
          endIf
        endIf
      If(JSym2E.eq.0) JSym2E = 3
      NFock = 0
      Call NameSC(3,IOpCl,IROHF,IHMeth,IExCor,IXCFnc,IEmpDi,ScaHFX,
     $  MethNm)
C
C     SCRF parameters:
C     a0 = cavity radius of the solute molecule
C     epsi = dielectric constant of the solvent
C     gFac = scale factor of the image dipole (reaction field)
C
      Junk = IPSave(0)
      Call GetRF(IOut,IPrint,Junk,ISCRF,ICharg,NE,PhyCon(1),IPCM,IA0,
     $  gFac,a0,Epsi,fuFac)
      Call FileIO(2,-IRwGen,LRwGen,DGen,0)
      ENR = DGen(41)
      RdPCM = SCIPCM(ISCRF).and.Epsi.eq.Zero
      If(CFPCM(ISCRF)) Call AMove(3,DGen(2),EField(1))
      If(DoPrnt) then
        Call NameSC(2,IOpCl,IROHF,IHMeth,IExCor,IXCFnc,IEmpDi,ScaHFX,
     $    SCFNam)
        LSCFNm = LinEnd(SCFNam)
        Write(IOut,1070) SCFNam(1:LSCFNm)
        endIf
      SvSort = ISvInt.eq.1.or.ISvInt.eq.2
C     Turn off saving in-core integrals until consumers can handle pure
      SvSort = .False.
      IFinIt = IOp(20)
      If(IFinIt.eq.0) then
        If(IDiag.eq.3.or.IDiag.lt.0.or.(DoDMS.and..not.DoSim)) then
          IFinIt = 3
        else if(ISleaz.eq.6.and.Sleazy) then
          IFinIt = 1
        else
          IFinIt = 2
          endIf
        endIf
      If(SpDiag) then
        If(ISCRF.ne.0) Call GauErr('No SCRF with CG yet.')
        If(ROHF.and..not.RHF) Call GauErr('No ROHF so far.')
        endIf
      If(SpFock) IntTyp = 1
      If((SpFock.or.SpDiag).and..not.DirSCF)
     $  Call GauErr('No O(N) storage with conventional SCF.')
      If(IHMeth.le.1.and.(IFitEn.eq.1.or.Mixed.or..not.DFT)) then
        If(IFitEn.eq.1.and.IDenFit.gt.1) then
          IDenFit = -IDenFit
        else
          IDenFit = 1
          endIf
        endIf
      If(DoPrnt.and.IDenFit.ne.1) Write(IOut,1390) IDenFit, NucTyp
      If(IntTyp.eq.0) then
        If(IDenFit.gt.1.or.(PCM(ISCRF).and..not.DirSCF).or.
     $    DirSc0.eq.4) then
          IntTyp = 1
        else if(DirSc0.eq.3.or.DoStCy.or.TstSpI) then
          IntTyp = 3
          endIf
        endIf
      AM1 = AM1Ham(IHMeth)
      RComp = RHF.and.IOp(12).eq.1.and..not.PCM(ISCRF).and.IDenFit.eq.1
      KAlg = IOp(92)
      If(IntTyp.eq.0.and.NoInCo(DFT,Mixed,IPFlag,FMFlg1,IHMeth,IDenFit,
     $  IRaf,IHFX,Omega,LSEAll)) IntTyp = 1
      If(IOp(64).ne.0) FMFlag = IOp(64)
      If(IOp(65).ne.0) NFxFlg = IOp(65)
      NSAStep = IOp(71)
      NAEBr = IOp(75)
      NOABr = IOp(76)
      NBEBr = IOp(77)
      NOBBr = IOp(78)
      If(NAEBr.gt.0.or.NBEBr.gt.0) IFBrod = -1
      NClRep = NCelOv
      ITemp = IOp(70)
      IZMthd = IOp(74)
      If(IZMthd.eq.0) IZMthd = 911
C     If(IZMthd.eq.0) IZMthd = 193
      TempSt = Zero
      FKTemp = Zero
      Deg2AU = Boltz / Hartre
      kTLim = IOp(79)
      If(UsePBC) then
        If(ITemp.eq.0.and.(Mod(IZMthd,10).eq.2.or.Mod(IZMthd,10).eq.3))
     $    ITemp = IDefTP
        If(ITemp.eq.-2) ITemp = 0
        FKTemp = GFloat(ITemp)*Deg2AU
      else if(IFBrod.ge.0) then
        If(ITemp.eq.-1) ITemp = 1000000
        If(ITemp.eq.0) ITemp = IDefT
        If(ITemp.eq.-2) ITemp = 0
        TempSt = GFloat(ITemp)
        endIf
      If0K = ITemp.le.0
      NSPnt = IOp(99)
      If(NSPnt.eq.0) then
        If(NPDir.eq.1) then
          NSPnt = 32
        else if(NPDir.eq.2) then
          NSPnt = 12
        else if(NPDir.eq.3) then
          NSPnt = 8
          endIf
        endIf
      TigPCV = IOp(100).eq.1
      ITstSD = IOp(101)
      If(ITstSD.eq.0) ITstSD = IDefSD
      ThrSDi = GFloat(ITstSD)
C
C     Sparse matrix diagonalization replacement parameters
C
      DoSCF(1) = IOp(90).ge.1
      DoSCF(2) = IOp(90).eq.2.and.ISpTst(IRwPOA,NBasis).eq.2
      DoSCF(3) = IOp(90).eq.2.and.ISpTst(IRwPOA,NBasis).eq.2
      IDLaMx = IOp(80)
      IFxFM = 0
      IMemOp = 0
      If(IDLaMx.eq.0) IDLaMx = 80
      If(IDLaMx.lt.0) IDLaMx = 1000
      SECut = GFloat(10)
      DCut = -One
      IOpDMe = IOp(81)
      IOpDM2 = IOp(91)
      IConCG = IOp(82)
      If(IConCG.eq.0) IConCG = 7
      DThr = GFloat(10)**(-IConCG)
      MxSCFC = IOp(83)
      If(DIIS.and.MxSCFC.eq.0) then
        If(NOnlyS.gt.0.or.ICampK.ne.0) then
          MxSCFC = Min(2*NDIDef,MaxCyc)
        else
          MxSCFC = Min(NDIDef,MaxCyc)
          endIf
        endIf
      DIIS = DIIS.and.MxSCFC.gt.1
      If(DIIS.and.DoPrnt) Write(IOut,1260) IDIIS
      MxSCF1 = MxSCFC + 1
      IOrPBC = IOp(89)
      IRANP0 = IOp(96)
      UpdTim = IOp(97).ne.2
      SavPBC = IOp(98).ne.2.and.UsKPnt
      NPrnKE = IOp(103)
      ITstGF = IOp(104)
      If(ITstGF.ne.0.and.RHF.and.IUHFT.eq.0) IUHFT = 1
      FrcNuc = UsePBC.or.LBit(IPFlag,19)
      If(FrcNuc) ENR = Zero
      UseGEW = UseGEW.or.DoSim
      If(IADBas.eq.0.and.(DoSim.or.DoDMS)) then
        IADBas = 1
      else if(IADBas.eq.0) then
        IADBas = 2
        endIf
      UseNuc = IOp(120).eq.1
      IDoVI = IOp(121)
      IDoHS = IOp(122)
      DoHirS = IDoHS.eq.1.or.
     $  (IDoHS.eq.0.and.IPSave(0).eq.0.and.NDimBl.eq.2)
      IncFoc = (IIFock.eq.2.or.(IIFock.eq.0.and.DirSCF)).and.NFullF.ge.0
     $  .and..not.SpFock.and..not.SpDiag.and.IDenFit.le.2
     $  .and.VarAcc.ne.2
      If(NSAStep.eq.-1) then
        NSAStep = 1000000
      else if(NSAStep.eq.0) then
        NSAStep = 10
        endIf
      Call IClear(NProgC,IProgC)
      IOSort = FilNum(IOSort,IUnit(13))
      NBSq = NBasis*NBasis
      NBNK = NBasis*NKPnt
      NBSqK = 2*NBSq*NKPnt
      NBNU = NRI*NBD*NBasUD
      If(ZDO) then
        Write(IOut,1270)
        IOSUse = 0
      else if(DFTB) then
C       CGDMS clobbers IOSUse and the full S is needed for FoFSE
C       for EHTSC and SC-DFTB.
        IOSUse = IGetRW(NTT,0)
        Call Trsfr(NTT,IRwS,IOSUse,V(JStrt),LenV)
      else
        IOSUse = IRwS
        endIf
      ICntr1 = Mod(ICntrl,100)/10
      If(ICntr1.eq.1) then
        IChuck = 2
      else
        IChuck = 0
        endIf
      If(ICntr1.eq.1.or.ICntr1.eq.3) then
        IAcc0 = 1
      else if(Sleazy) then
        IAcc0 = 2
      else
        IAcc0 = 0
        endIf
      If(IAccFD.eq.0) then
        AccDes = DecAcc(IAcc0,Zero,1)
        If(Acurcy.lt.CnvDef) then
          AccDes = AccDes * Acurcy / CnvDef
        else
          AccDes = Zero
          endIf
      else
        AccDes = GFloat(10)**(-IAccFD)
        endIf
      NumChg = 0
      If(FrcNuc) then
        IRdChg = 31
      else
        IRdChg = 3
        endIf
      Call RdChg(IRdChg,0,NumChg,NumDBF,LChg,LDBF,jStrt-1,jChg,jDBF,IV,
     $  V,NGot)
      jStrt = IV
      IPrtF = IPrint - 2
      Call PrmDim(IOut,IPrtF,IHMeth,IDenFit,0,1,0,0,0,AccDes,IChuck,0,
     $  IPFlag,FMM,FMFlag,FMFlg1,NFxFlg,0,JJA,JJA,NPrDim,LenDup,.True.,
     $  V(jIPBC),CelVec,NAtoms,C,ScaHFX,Zero,0,0,XXA,XXA,NBT,NFixJE,
     $  LFDJE,NClLst,IRaf,V(jStrt),V(jStrt),LenV)
      If(FrcNuc.or.NucTyp.eq.3) NBT = NBT + NAtoms + NumChg
      If(IDenFit.gt.1) NBT = NBT + NDBF6D
      ISort = ITqry(IRwIBS).gt.0
      If(ISort) then
        Write(IOut,1230)
        Call IClear(30,IBF)
        Len = InToWP(10)
        Call FileIO(2,-IRwIBS,Len,IBFSrt,0)
        IRSort = (IBFSrt(9).ge.1.and.IBFSrt(9).le.3)
        If(.not.IRSort)
     $    Call GauErr('UHFOpn cannot handle sorted regular integrals.')
        ISort = .not.IRSort
        nwiib = ibfsrt(3)
        isym2e = ibfsrt(10)
        ibf(1) = ibfsrt(9)
        ibf(2) = ibfsrt(9) + 1
      else
        Call IClear(10,IBFSrt)
        IRSort = .False.
        Call FileIO(2,-IRwIBF,InToWP(30),IBF,0)
        If(DoPrnt.and..not.DirSCF) Write(IOut,1220)
        endIf
      Call GetMT(MulTab)
      NBlks = 1
      NOATot = NAE + (NOABr-NAEBr)
      NOBTot = NBE + (NOBBr-NBEBr)
      NOccSy(1) = NOATot + (NDimBl-1)*NOBTot
      NOccSy(2) = NOBTot
      NBF(1) = NBasis
      NBFU(1) = NBsUse*NDimBl
      MaxCmp = 0
      Call FileIO(2,-IOSBlk,InToWP(52),NBlksS,0)
      Call GetFrg(IOut,IPrint,0,IFrgTp,.True.,NAtoms,NBasis,IAtFrg,
     $  jStrt,NFrag,IHavMO,NACore,NAVal2,NAVal1,NBCore,NBVal2,NBVal1,
     $  IFrStp,NDNBsC,NumNBC,IIChOn,IMulOn,INBCom,IMpF2M,IMpM2F,IIndFr,
     $  IIFrSt,IIFrAt,Next,V,NGot)
      jStrt = Next
      If(IFrgTp.gt.0) then
        NBlks = NFrag
      else if(AbDiag) then
        NBlks = NBlksS
        endIf
      INBF = jStrt
      INBFU = jStrt + InToWP(NBlks)
      INocSy = INBFU + InToWP(NBlks)
      jStrt = INocSy + InToWP(2*NBlks)
      Call TstCor(jStrt,NGot,'UHFOpn-NBF')
      If(AbDiag) then
        Call IMove(NBlks,NBFS,V(INBF))
        Call IMove(NBlks,NBFUS,V(INBFU))
      else
        Call IMove(NBlks,NBF,V(INBF))
        Call IMove(NBlks,NBFU,V(INBFU))
        endIf
C
C     Allocate space.  First set up k-space arrays.
      jWKPnt = 1
      iNBFK = 1
      jCmlxF = 1
      jICoor = 1
      jKCoor = 1
      jIOToK = 1
      jEnd = jStrt
      If(UsKpnt) then
        iNBFK = jEnd
        jEnd = iNBFK + intowp(3*NKPnt)
        Call GnKArr(IOut,IPrint,V(jIPBC),.False.,NKPnt,NDimAK,.True.,
     $    jICoor,jKPntF,jWKPnt,NTot,jIOToK,jCmlxF,jKCoor,jEnd,V,NGot)
        endIf
      If(UsePBC.and.(IPSave(0).eq.0.or.IPrint.gt.0))
     $  Write(IOut,1610) NCelOv, NClRep, NMtPBC, NRecip, NKPnt
C
C     NAtom to basis function mapping arrays.
C
      If(NumChg.gt.0) ISym2E = 0
      IFXYZ = jEnd
      ILLim = IFXYZ + 2*LFXYZ
      IULim = ILLim + InToWP(NAtoms+1)
C     NEqBas(NOp1,NBasis) (only if int symmetry).
      INEq = IULim + InToWP(NAtoms)
      LenNEq = 0
      If(ISym2E.eq.1) LenNEq = InToWP(NOp1*NBasis)
      IIJMap = INEq + LenNEq
      LenIJM = 0
      If(ISym2E.eq.1.and.(JSym2E.eq.2.or.JSym2E.eq.3).and.NOp1.gt.1.and.
     $  .not.DirSCF) LenIJM = InToWP((NOp1-1)*NTT)
      ICP = IIJMap + LenIJM
C
C     NEqShl(MaxShl,NOp1) (only if int symmetry and direct SCF).
      INEqSh = ICP + (MaxCyc+1)*3
      If(ISym2E.eq.1) then
        LenNES = Max(ITqry(NEqShl),0)
        If(HarTst) then
          LenNE2 = InToWP(NOpAll*NAtoms)
        else
          LenNE2 = Max(ITqry(NEqSh2),0)
          endIf
      else
        LenNES = 0
        LenNE2 = 0
        endIf
      If(ISym2E.eq.0) then
        If(IPSave(0).eq.0) Write(IOut,1290)
        JSym2E = 0
      else if(JSym2E.eq.1) then
        If(IPSave(0).eq.0) Write(IOut,1300)
      else if(JSym2E.eq.2) then
        If(IPSave(0).eq.0) Write(IOut,1310)
      else if(JSym2E.eq.3) then
        If(IPSave(0).eq.0) Write(IOut,1490)
      else
        Call GauErr('Illegal value of JSym2E in UHFOpn.')
        endIf
      INEqS2 = INEqSh + LenNES
C
C     IJ(NBasis+1).
      IIJ = INEqS2 + LenNE2
C
C     D(N,N) and DD(N) -- DA(NTT) and DB(NTT)
      NAllBT = IGFix(Pt5*(One+Sqrt(One+Eight*GFloat(NBT))))
      NAlloc = Max(NFullF+1,MxSCF1,8)
      If(DoStCy.or.TstSpI) then
        NBsSp6 = NBas6D*2
        NRIX = 2
      else
        NBsSp6 = NBas6D*NDimBl
        NRIX = NRI
        endIf
      If(.not.SpFock.or..not.SpDiag) NAlloc = Max(NAlloc,NBsSp6)
      If(SpFock.or.SpDiag.or.FrcNuc.or.NucTyp.eq.3)
     $  NAlloc = Max(NAlloc,NAllBT)
      NTT6DS = (NBsSp6*(NBsSp6+1))/2
      NAllSq = NAlloc**2
      NAllTT = Max((NAlloc*(NAlloc+1))/2,NDBF6D)
      If(NMtPBC.gt.1) NAllTT = Max(NAllTT,NTT6DS*NMtPBC)
      If(DISPCM) NAllTT = NAllTT + 2*NTs
      NAlloc = NRIX*NAlloc
      NAllTT = NRIX*NAllTT
      NAllSq = NRIX*NAllSq
      NBT = Max(NBT,NAllTT)
      IEnCyc = IIJ + InToWP(NAlloc+1)
      IErCyc = IEnCyc + MxSCF1
      IIndFk = IErCyc + MxSCF1
      NSkewF = NSkew + MDPgSz(0)
      ID = IIndFk + InToWP(MxSCF1)
      IDA = ID
      IDD = ID + NAllSq
C
C     F(N,N) and FF(N) -- FA(NAllTT) and FB(NAllTT) and Scr(N,5)
      If(RComp) then
        IDB = IDA
        IF = ID
        IFA = IDA + NAllTT
        IFB = IFA
        IFF = IDD
      else
        IDB = IDA + NAllTT
        IF = Max(IDD+NAlloc+NSkewF,IDB+NAllTT)
        IFF = IF + NAllSq
        IFA = IF
        IFB = IF + NAllTT
        endIf
      IScr2 = Max(IFF+NAlloc+NSkewF,IFB+NAllTT)
C
C     V(N,N) and VV(N) -- integral buffer.
      LScr2 = Max(8,NProc(0)+2)
      LenScr = LScr2*NAlloc + NSkew
      If(Onsagr(ISCRF)) LenScr = Max(LenScr,2*NTT)
      IBuf2E = IScr2 + LenScr
      If(DirSCF) NWIIB = 0
      Call ReadGW(IOut,IPrint,DFT,1,0,IBuf2E,MxGWIn,IGWInf,IAtBtD,
     $  IIRBtD,IRRBtD,IGWBtD,Next,V,LenV)
      If(IPSave(0).eq.0.and.IGWInf(1).gt.0.and.Next.gt.IBuf2E)
     $  Write(IOut,1620) (Next-IBuf2E)
      IBuf2E = Next
      IV = IBuf2E
      IMin1 = IBuf2E + 2*NAllTT - 1
      If(DIIS.and..not.SpDiag) IMin1 = IMin1 + (MaxCyc+1)**2
      If(IDiag.eq.2) IMin1 = Max(IMin1,NBasis*(NAE+8),NAE*(NBasis+NAE))
      IMin = Max(IMin1,IBuf2E+2*InToWP(NWIIB)-1)
      IEnd = IMin
      NGetB = IEnd
      If(.not.DirSCF.and..not.ISort.and..not.IRSort.and.IFil.le.2)
     $  NGetB = IMin1 + IFil*InToWP(NWIIB)
      IVT = Max(IEnd,NGetB)
      MDVT = NGot - IVT + 1
      MnInMem = IntMem(-1,0) + NPrDim
      If(PCM(ISCRF)) then
        Call MemPCM(IOut,IPrint,'ExactFock',1,LenPCM,IJunk,IJunk,IJunk,
     $    IPFlag,FMM,FMFlag,FMFlg1,NFxFlg,IHMeth,NAtoms,C,1,1,V(IVT),
     $    MDVT)
        MnInMem = Max(MnInMem,LenPCM)
        endIf
      If(.not.AM1) MnInMem = Max(MnInMem,MemRaf(NBas6D))
      If(DFT) then
        Call DecPrn(IOut,IPrint,IRadAn,0,IPrune,NTheta,NPhi,MaxRad,
     $    MaxAng,ECTolR)
        Call GetExC(IExCor,IExChn,ICorr)
        Call DFTMem(IOut,IPrint,1,10,IExChn,ICorr,MaxRad,MaxAng,NShell,
     $    NAtoms,NClRep,1,NBas6D,NBas6D,0,0,Mem0,Mem1,Mem2,MemM0,MemM1,
     $    MemM2)
        MnInMem = Max(MnInMem,Mem0+MemM0)
        endIf
      Call SCIOCM(IOut,IPrint,.True.,IHMeth,IPFlag,IRaf,DFT,Mixed,IOpCl,
     $  IROHF,0,0,0,SvSort,ISym2E,NOp2,NBasis,NBas6D,NShell,ShellT,
     $  ShellC,JAn,MOmega,Omega(1,IHFJ),Omega(1,IHFX),IVT,0,MnInMem,
     $  JONEKO,IOpClM,IOpClF,LR1,LR2,LR3,NRn,NX,LenR,LRAOJ,LSAOP,NIR2E,
     $  NSAOBk,MulTab,ChrTbl,NS,NS2,NS4,jNSeqB,jISABF,jRSABF,jMap,jNext,
     $  NGetL,V,NGot)
      If((ISMode.gt.0.or.DirSCF).and..not.RComp.and.IOpClM.ge.-1) then
        If(DoStCy.or.TstSpI) then
          If(SvSort) Call GauErr('SvSort and DoStCy in UHFOpn.')
          LR1 = 8
          LR2 = 16*LR1
          LR3 = LR2
          NRn = 3
          IOpClM = 3
          NGetL = IMin + Max((LR1+LR2+LR3)*LenR,LRAOJ) + MnInMem + LSAOP
          endIf
      else
        LR1 = 0
        LR2 = 0
        LR3 = 0
        NGetL = -1
        endIf
      If(IPrint.gt.0) Write(IOut,1210) NBasis, NAllTT, NAllSq, NWIIB,
     $  INEq, LenNEq, IIJ, ID, IDD, IDA, IDB, IF, IFF, IFA, IFB, IScr2,
     $  IV, IBuf2E, IMin1, IMin, IVT, IOpCl
      IR1 = IBuf2E + NSkew
      IR2 = IR1
      IR3 = IR2
      NSAOBk = 0
      jNEqB = IR1
      jISABF = jNEqB
      jRSABF = jISABF
      jMap = jRSABF
      If(.not.ISort.and..not.IRSort.and.NGetL.gt.0.and.DoCyc.and.
     $  NGot.ge.NGetL.and.(IntTyp.lt.0.or.(IntTyp.eq.0.and.DirSCF).or.
     $  IntTyp.eq.3)) then
        Call SetICI(IOut,IPrint,IHMeth,IOpClM,LSAOP,NOp2,NBasis,NBas6D,
     $    LR1,LR2,LR3,LRAOJ,NGetL,NIR2E,IVT,InStor,NSAOBk,MulTab,ChrTbl,
     $    NS,NS2,NS4,jNEqB,jISABF,jRSABF,jMap,IR1,IR2,IR3,jVInt,jVSAO,
     $    jVPerm,V,NGot)
        IBuf2E = jVPerm
        IV = IBuf2E
      else if(IntTyp.eq.3) then
        Call ICErro(IOut,NGetL,NGot)
      else if((IntTyp.eq.0.or.IntTyp.eq.2).and..not.ISort.and.
     $  .not.IRSort.and.IFil.le.2.and.NGot.ge.NGetB.and..not.DirSCF)then
        InStor = 2
        IBuf2E = IMin1 + 1
        If(ROHF.and..not.RHF) IV = IBuf2E + 2*InToWP(NWIIB)
        MemInt = .True.
        If(IPSave(0).le.1) Write(IOut,1320) NGetB
      else if(IntTyp.eq.0.or.IntTyp.eq.1) then
        InStor = 1
        Call TstCor(IMin,NGot,'UHFOPN')
        MemInt = .False.
      else
        Write(IOut,1350) IntTyp
        Call Lnk1E(0)
        endIf
      IVB = IV + NAllTT
      IVV = IV + Max(NAllSq,NAllTT)
      LenVV = NAlloc
      If(DIIS.and..not.SpDiag) LenVV = Max(LenVV,MxSCF1**2)
      IBMat = IVV + LenVV
      Call TstCor(IBMat+MxSCF1**2,NGot,'UHFOpn-Test1')
      MDV = NGot - IBuf2E + 1
      IndX = IBuf2E + 2*InToWP(NWIIB)
      LenX = MDV - 2*InToWP(NWIIB)
      IndY = IBMat
      LenY = NGot - IndY + 1
      If(DoPrnt) Write(IOut,1250) IVT, NGetB, NGot, MDV, LenX, LenY
      If(LenNEq.ne.0) Call FileIO(2,-NEq,LenNEq,V(INEq),0)
      If(LenNES.ne.0) Call FileIO(2,-NEqShl,LenNES,V(INEqSh),0)
      If(LenNE2.ne.0) then
        If(HarTst) then
          Call IMove(NAtoms*NOpAll,NEqAll,V(INEqS2))
        else
          Call FileIO(2,-NEqSh2,LenNE2,V(INEqS2),0)
          endIf
        endIf
      If(LenIJM.gt.0)
     $  Call LIJMap(NBasis,NTT,NOp1,V(INEq),V(IIJMap))
C
C     Make sure input density is consistent with fraction occ
C
      If(NAEBr.gt.0.or.NBEBr.gt.0) then
        Call AClear(NBsUse,V(IVV))
        Call FONDen(IOut,IPrint,'alpha',.True.,NRI,NBasis,NAE-NAEBr,
     $    NAEBr,NOABr,0,0,IRwCA,IRwPA,V(IVV),V(IF),V(IV))
        Call FONDen(IOut,IPrint,'beta',.True.,NRI,NBasis,NBE-NBEBr,
     $    NBEBr,NOBBr,0,0,IRwCB,IRwPB,V(IVV),V(IF),V(IV))
        endIf
C
C     Do initialization for ADMP, full-matrix QNDMS, and SimOpt.
C
      IRwPAN = 0
      IRwPBN = 0
      If(UsKPnt) then
        JD  = IndX
        JDD = JD  + Max(NBSq,2*NTT*NKPnt)
        JVX = JDD + Max(NBSq,NTTX)
        JVV = JVX + NBSqK
        JVY = JVV + 4*NBSq
      else if(ISpFok.eq.3) then
        JD  = ID
        JDD = JD  + NBSq
        JVX = JDD + NBSq
        JVV = JVX + NBSq
        JVY = IBuf2E
      else
        JD  = ID
        JDD = IDD
        JVX = IV
        JVV = IVV
        JVY = IBuf2E
        endIf
      MDVY = NGot - JVY + 1
      If(IADMP.gt.0) then
        Call TstCor(JVY,NGot,'Before ADMIni')
        Call ADMIni(IOut,SpDiag,DoPrnt,UsKPnt,IADMP,NMtPBC,NKPnt,NPDir,
     $    NTTX,IADBas,ZDO,RHF,NBasis,NBsUse,IRwS,IRwPOA,IRwPOB,IRwPA,
     $    IRwPB,IRwRAO,IRwV,IRwPAN,IRwPBN,V(jIPBC),V(jICoor),V(iNBFK),
     $    V(jCmlxF),V(jWKPnt),V(JD),V(JDD),V(JVX),V(JVV),V(JVY),MDVY)
        endIf
      If(DoSim.or.DoDMS) Call QNDIni(IOut,DoPrnt,DoDMS,DoSim,DoSimP,
     $  UseGEW,RHF,ZDO,IADBas,NBasis,NBsUse,IRwS,IRwRAO,MaxCmp,IRwV,
     $  IRwSM,IRwPA,IRwPB,IRwFRA,IRwFRB,NAE,NBE,Acurcy,V(ID),V(IDD),
     $  V(IF),V(IFF),V(IndY),LenY)
      If(IAccMx.eq.0) then
        If(DoSim) then
          IAccEn = 0
        else
          IAccEn = 2
          endIf
        IAccMx = 2
      else if(IAccMx.gt.0) then
        IAccEn = Mod(IAccMx,100)/10
        IAccMx = Mod(IAccMx,10)
      else
        IAccEn = IAccMx
        endIf
      If(IAccMx.gt.0) then
        AccMax = (GFloat(10)**IAccMx)*Acurcy
      else if(IAccMx.eq.-1) then
        AccMax = Acurcy
      else if(IAccMx.eq.-2) then
        AccMax = GFloat(100)
      else
        Write(IOut,1480) IAccMx
        Call Lnk1E(0)
        endIf
      If(Sleazy.or.ISleaz.eq.6) then
        Write(IOut,1430)
        EngCon = Acurcy / Two
      else if(ISleaz.eq.3) then
        EngCon = DGAccE
        If(IAccFD.eq.0) IAccFD = 8
        Write(IOut,1590) IAccFD
      else if(IAccEn.lt.0) then
        EngCon = Zero
      else if(IAccEn.eq.0) then
        EngCon = Acurcy
      else
        EngCon = Acurcy * GFloat(10)**(IAccEn)
        endIf
      If(DoPrn1) then
        If(IAccMx.eq.-1) then
          Write(IOut,1170) AccMax, MaxCyc
        else
          Write(IOut,1160) Acurcy, MaxCyc
          If(IAccMx.ne.-2) Write(IOut,1180) AccMax
          endIf
        If(EngCon.ne.Zero) Write(IOut,1190) EngCon
        If(IVShft.gt.0) then
          Write(IOut,1200) VSGoal
        else if(IVShft.le.-2) then
          If(IDyVSh.le.2) Write(IOut,1540) VSGoal, VSMax
          If(IDyVSh.eq.3) Write(IOut,1541) VSGoal
          If(IDyVSh.eq.4) Write(IOut,1542) VSGoal
          If(IDyVSh.ge.5) Write(IOut,1543) VSGoal, VSMax, VSGpMn
          endIf
        Write(IOut,1560) NameFr(FixRis)(1:LinEnd(NameFr(FixRis)))
        If(IFBrod.ge.0)
     $    Write(IOut,1570) OnOff(LvWFON)(1:LinEnd(OnOff(LvWFON)))
        If(FixRis.ge.2)
     $    Write(IOut,1580) OnOff(LvRise)(1:LinEnd(OnOff(LvRise)))
        endIf
C
C     Initialize files.
C
      IRwFOA = 0
      IRwFOB = 0
      SEEn = Big
      GesEn = Big
      If(SpDiag) then
        IRwFOA =  IGetRW(1,0)
        If(.not.RHF) IRwFOB = IGetRW(1,0)
        endIf
      ZThr = SPCutO(0)
      IOffS = NRI*NTTSp*Min(NMtPBC-1,1)
      If(SpDiag) then
        LnC123 = 10
      else
        Call FilInd(NAlloc+1,V(IIJ))
        Call AClear(4*NAllTT,V(ID))
        If(.not.RComp) Call AClear(2*NAllTT,V(IV))
        If(NMtPBC.gt.1) then
          Call DupMat(NMtPBC,NTT,IRwPA,NMtPBC,V(ID))
          If(.not.RHF) Call DupMat(NMtPBC,NTT,IRwPB,NMtPBC,V(ID))
          endIf
        Call FileIO(2,-IRwPA,NTTX,V(ID),0)
        Call FileIO(2,-IRwS,NTTX,V(IV),IOffS)
        TNEA = SCFTrX(.False.,.False.,IOpClX,NMtPBC,NBasis,V(ID),XXA,
     $    V(IV),XXA)
        If(RHF.and.NDimBl.eq.1) then
          If(Abs(TNEA-GFloat(NE)).lt.Pt5) then
            TNEA = TNEA*Pt5
            Call AScale(NTTX,Pt5,V(ID),V(ID))
            endIf
          Call FileIO(1,-IRwPA,NTTX,V(ID),0)
          Call ConDDF(IRwPB,NTTX)
          Call FileIO(1,-IRwPB,NTTX,V(ID),0)
          Call FileIO(2,-IRwEig,NBsUse,V(ID),0)
          Call ConDDF(IRwEig,2*NBsUse*NDimBl)
          Call FileIO(1,-IRwEig,NBsUse,V(ID),0)
          Call FileIO(1,IRwEig,NBsUse,V(ID),0)
          Call ConDDF(IRwCB,NBNU)
          TNEB = TNEA
        else if(NDimBl.eq.1) then
          Call FileIO(2,-IRwPB,NTTX,V(ID),0)
          TNEB = SCFTrX(.False.,.False.,IOpClX,NMtPBC,NBasis,V(ID),XXA,
     $      V(IV),XXA)
        else
          TNEB = Zero
          endIf
        If(NDimBl.eq.1) then
          Error = Abs(TNEA-NAE).ge.ErPMax.or.Abs(TNEB-NBE).ge.ErPMax
        else
          Error = Abs(TNEA-NAE-NBE).ge.ErPMax
          endIf
        If(Error) then
          Write(IOut,1010) TNEA, NAE, TNEB, NBE
          Call GauErr('Input densites are not normalized.')
          endIf
        Call ConDDF(IRwPT,NTTX)
        Call ConDDF(IRwPS,NTTX)
        Call ConDDF(IRwFA,NTTX)
        Call ConDDF(IRwFB,NTTX)
        If(.not.RHF) Call ConDDF(IRwUR,NBsUse**2)
        LnC123 = 2*NTTX
        endIf
      IRwC1 = IGetRZ(LnC123,0,V(ID),LnC123)
      IRwC2 = IGetRZ(LnC123,0,V(ID),LnC123)
      IRwC3 = IGetRZ(LnC123,0,V(ID),LnC123)
C
C     Allocate RWF's for the K space matrices, and form the transformation matrices
C     to orthonormal coordinates if these have not been provided.  Also set up the
C     matrices for density fitting if necessary.
C
      IRwCAC = 0
      IRwCBC = 0
      IRwEiC = 0
      If(UsKpnt) then
        LenEig = NBasis*NKPnt*NSpBlk
        If(SavPBC) then
          Call ConDD2(IRwECK,NBSqK*NSpBlk+2*LenEig)
          IRwEiC = IGetSF(IRwECK,LenEig,0)
          IRwWgt = IGetSF(IRwECK,LenEig,LenEig)
          IRwCAC = IGetSF(IRwECK,NBSqK,2*LenEig)
          If(NSpBlk.eq.2) IRwCBC = IGetSF(IRwECK,NBSqK,2*LenEig+NBSqK)
        else
          IRwCAC = IGetRW(NBSqK,0)
          If(NSpBlk.eq.2) IRwCBC = IGetRW(NBSqK,0)
          IRwEiC = IGetRW(LenEig,0)
          IRwWgt = IGetRW(LenEig,0)
          endIf
        If(IADMP.gt.0) then
          IRwSk = IGetRW(2*NTT*NKPnt,0)
          Call ConDDF(IRwV,2*NBSqK+NBNK)
          IRwS12 = IGetSF(IRwV,NBSqK+NBNK,NBSqK)
        else
          IRwSk = 0
          IRwS12 = 0
          endIf
        Call ZFormV(IOut,IPrint,V(jIPBC),IOrPBC,NBasis,NKPnt,V(iNBFK),
     $    V(jICoor),V(jCmlxF),IRwS,IRwCAC,IRwSk,IRwS12,V(IndX),LenX)
        If(IADMP.gt.0) Call Trsfr(NBSqK,IRwCAC,IRwV,V(IndX),LenX)
        Call SavKPn(0,NKPnt,NDimAK,V(iNBFK),V(jICoor),V(jKPntF),
     $    V(jCmlxF),V(jWKPnt),V(IndX),LenX)
        If(IRwCBC.ne.0) Call CpyFil(IRwCAC,IRwCBC,V(IndX),LenX)
        endIf
      If(DIIS) IDIIS = IDIIS + 1
      If((DIIS.or.IDiag.eq.2.or.SpDiag).and..not.ZDO) then
        If(IFrgTp.eq.1) then
          IRwTM = IRAOFr
        else if(IFrgTp.gt.0) then
          IRwTM = IGetRw(NSpBlk*NRI*NBD**2,0)
        else if(ITqry(IRwRAO).gt.0) then
          IRwTM = IRwRAO
        else if(SpDiag) then
          IRwTM = 0
          Call GauErr('U matrix not found for sparse.')
        else
          IDiag1 = 1
          If(IDiag.eq.2) IDiag1 = 2
          IRwTM = IGetRW(IDiag1*NBSq,0)
          Call TStamp(1,'Bef FormV')
          Call FormV(IDiag1,NBasis,NBasis,IRwS,IRwTM,V(ID),V(IV),V(IDD),
     $      V(IVV))
          Call TStamp(1,'Aft FormV')
          endIf
      else
        IRwTM = 0
        endIf
C
C     If using the first-order approximation to XC, update a copy of the
C     core Hamiltonian.
C
      ENJ = Zero
      If(IAprXC.ne.0) then
        OKSym = ISym2E.ne.0
        Call RdWrB(-18,0)
        Call HarFok(IOut,IPrint,IDoVI,.True.,.True.,.False.,OKSym,
     $    .False.,.True.,IExCor,IRadAn,IRanWt,IRanGd,AccDes,IOpCl,1,
     $    NAtoms,NBasis,NBas6D,NBT,ICharg,Multip,NE,IAn,IAtTyp,AtmChg,C,
     $    ScaDFX,Junk,ENRXC0,ENJ,V(IDA),V(IFA),V(jIPBC),CelVec,IGWInf,
     $    V(IAtBtD),V(IIRBtD),V(IRRBtD),V(IGWBtD),XXA,V(IBuf2E),MDV)
C       Restore regular fitting functions if necessary.
        If(IDenFit.gt.1.and.IHMeth.ne.14) Call RdWrB(-10,0)
        IRwHXC = IGetRw(2*NTT,0)
        Call FileIO(1,-IRwHXC,NTT,V(IDA),0)
        If(IOpCl.eq.1) then
          Call FileIO(1,IRwHXC,NTT,V(IFA),0)
        else
          Call FileIO(1,IRwHXC,NTT,V(IDA),0)
          endIf
        IExCor = -2
      else if(FrcNuc.and..not.UsePBC) then
        IRwHXC = IRwT
        ENRXC0 = ENR
      else
        IRwHXC = IRwH
        ENRXC0 = ENR
        endIf
C
C     Load the integrals if they are to be kept in memory:
C
      If(InStor.ge.3) then
        IType = -2
      else
        IType = 2
        endIf
      If(RHF) IType = IType/2
      JSym2M = JSym2E
      If(.not.DirSCF) then
        NOpUse = NOp1
      else
        NOpUse = NOpAll
        endIf
      If(InStor.eq.3.or.InStor.eq.4) then
        If(DirSCF) then
          If(Sleazy.and.VarAcc.eq.2) AccDes = Acurcy / GFloat(100)
          If(.not.DFT.and..not.SCIPCM(ISCRF).and..not.FkGues) VarAcc = 1
          NOpUse = NOp1
          IncFoc = IncFoc.and.(NFull.gt.0.or.IIFock.eq.2)
          If(IRaf.eq.-1.or.(JSym2E.ge.2.and.InStor.ne.4)) JSym2M = 0
          JSym2M = Min(JSym2M,1)
          If(JSym2M.eq.0) NOpUse = 1
          If(No2E) then
            Write(IOut,1650)
            Call AClear((LR1+LR2+LR3)*(LenR+NSkew),V(IR1))
          else
            DoSPI = DoStCy.or.TstSpI
            Call StIC2E(IOut,IPrint,IHMeth,IOpClF,IOpClM,ICntrl,InStor,
     $        JONEKO,NX,DoSPI,NAtoms,NBasis,NBas6D,IAn,IAtTyp,AtmChg,C,
     $        AccDes,ScaHFX,JSym2M,NOpUse,NOp1,V(INEqSh),IHFJ,IHFX,
     $        Omega,IPFlag,AllowP,LSEAll,V(ISEAll),MulTab,NSAOBk,NIR2E,
     $        NS,NS2,NS4,V(jISABF),V(jRSABF),V(jMap),LRAOJ,LR1,LR2,LR3,
     $        V(IR1),V(IR2),V(IR3),V(jVInt),NGot-JVInt+1,V(jVSAO),
     $        NGot-jVSAO+1)
            endIf
        else
          If(JSym2E.eq.1) then
            JSym2M = 1
          else
            JSym2M = 0
            endIf
          Call LodRaf(IOut,IPrint,IType,JSym2M,NOp1,NBasis,NTT,V(IIJ),
     $      V(IIJMap),V(IBuf2E),V(IBuf2E),V(IR1),V(IR1),V(IR2),V(IR2))
          endIf
        If(SvSort.and.InStor.eq.3) then
          Call ConDDF(IOSort,(LR1+LR2)*LRAOJ)
          Call FileIO(1,-IOSort,LRAOJ,V(IR1),0)
          If(IOpClM.eq.1) Call FileIO(1,IOSort,LRAOJ,V(IR2),0)
          endIf
        endIf
C
      If(IncFoc.and.IPSave(0).eq.0) Write(IOut,1330) NFullF
      If(IPrint.gt.0) Write(IOut,1630) VarAcc
      If(VarAcc.eq.1.or.(VarAcc.eq.4.and.FkGues)) then
        AccFac = Acurcy / GFloat(100)
      else if(VarAcc.eq.2) then
        AccFac = Acurcy
        Write(IOut,1440)
      else if(VarAcc.eq.3) then
        AccFac = Acurcy
        Write(IOut,1450) AccFac
      else if(VarAcc.eq.4.or.VarAcc.eq.7.or.VarAcc.eq.8) then
        AccFac = GFloat(10)**(-5)
        If(VarAcc.eq.4.and..not.SpFock) Write(IOut,1450) AccFac
      else
        Write(IOut,1460) VarAcc
        Call Lnk1E(0)
        endIf
      If(AbDiag.or.IFrgTp.gt.0) then
        IOCMOS = IGetRW(NSpBlk*NBNU,0)
        Call TstCor(2*NBsUse,LenY,'UHFOpn-GetSOc')
        Call GetSOc(IOut,ROHF,NRI,NDimBl,ISyOcc,IOSymA,NOATot,NOBTot,
     $    NBlks,V(INBFU),NBasis,NBsUse,IRwCA,IRwCB,IOCMOS,MulTab,IWfnSy,
     $    V(INOcSy),V(IScr2),V(ID),V(IndY),V(IV))
      else
        If(NBlks.ne.1) Call GauErr('NBlks.ne.1 w/o AbDiag.')
        Call IMove(2*NBlks,NOccSy,V(INOcSy))
        IOCMOS = 0
        endIf
      If(IStab.eq.2) then
        IRwCS = IGetRW(NSpBlk*NBNU,0)
        Call FileIO(2,-IRwCA,NBNU,V(ID),0)
        Call FileIO(2,-IRwCB,NBNU*(NSpBlk-1),V(ID+NBNU),0)
        Call FileIO(1,-IRwCS,NBNU*NSpBlk,V(ID),0)
      else
        IRwCS = 0
        endIf
      If(SpFock.or..not.DirSCF) then
        NOpUse = NOp1
        NOp2U = NOp1
      else
        NOpUse = NOpAll
        NOp2U = NOp2
        endIf
      IPrPFS = 0
      If(IPrint.gt.0) IPrPFS = Max(IPrint,2)
      LCA = ITqry(IRwCA)
      If(ROHF.and.LCA.gt.0.and.ITqry(IRwCB).le.0) then
        Call ConDDF(IRwCB,LCA)
        Call Trsfr(LCA,IRwCA,IRwCB,V(IndX),LenX)
        endIf
      Call SetILR(IOut,.True.,1,NRestI,IRData,NRestR,RRData,NRest,
     $  V(IIRDat))
C
C     Enter iteration routine CycOpn, possibly for semi-empirical guess
C     (ISCF=1), test energy of SE guess (ISCF=2), test energy of old guess
C     (ISCF=3), Vacuum SCF (ISCF=4), real calculation (ISCF=5).
C
      EVac = Zero
      If(PCM(ISCRF)) GNEl = PCMGNE(IOut,IPrint)
      Do 100 ISCF = 1, NSCF
        If(.not.DoSCF(ISCF)) goto 100
        If(ISCF.eq.NSCF) then
          ISCRF1 = ISCRF
        else
          ISCRF1 = 0
          endIf
        If(ISCF.eq.NSCF.and.DoSCF(IVac)) then
          VarAc1 = 1
        else
          VarAc1 = VarAcc
          endIf
        JustEn = ISCF.eq.ISEEn.or.ISCF.eq.IGesEn.or.IADMP.eq.1.or.FrcJEn
        If(ISCF.eq.ISEGes) then
          IRwHUs = IGetRw(1,0)
          DLaMax = GFloat(100)
          OCut = GFloat(15)
          Call WNDOIn(IOut,IPrint,ICharg,Multip,7,0,0,0,IRwHUs,0,NAtoms,
     $      NBasis,NOrbs,NAEUse,NBEUse,IAn,IAtTyp,AtmChg,C,XXA,XXA,XXA,
     $      XXA,.True.,OCut,DCut,V(ILLim),V(IULim),ENRUse,V(IndX),LenX)
          ZDOUse = .True.
          IExCoU = 0
          IHMeUs = 8
          NBasUs = NOrbs
          NBs6Us = NOrbs
          IDFit = 0
          IOFitX = 0
          IOGesA = 0
        else
          DLaMax = GFloat(IDLaMx) / GFloat(100)
          Call MkLULm(-1,0,NAtoms,V(ILLim),V(IULim),OK)
          If(SpDiag.and..not.OK)
     $      Call GauErr('Sparse requires sorted shells.')
          IHMeUs = IHMeth
          NAEUse = NAE
          NBEUse = NBE
          ZDOUse = ZDO
          IExCoU = IExCor
          NBasUs = NBasis
          NBs6Us = NBas6D
          ENRUse = ENRXC0
          IDFit = IDenFit
          IOFitX = IODFit
          IOGesA = IRwAGs
          IRwHUs = IRwHXC
          If(ISCF.eq.IVac.and.PCM(ISCRF)) then
            ENRUse = ENRUse - GNEl
            If(CFPCM(ISCRF)) then
              If(SProd(3,EField(1),EField(1)).ne.Zero) Call FixHCF(IOut,
     $          IPrint,IRwHXC,IRwHUs,IOpCl,ISCRF,EField(1),IPFlag,
     $          AllowP,FMM,FMFlag,FMFlg1,NFxFlg,IHMeth,NAtoms,IAn,
     $          IAtTyp,AtmChg,C,NBasis,NBas6D,V(IndX),LenX)
              endIf
            endIf
          endIf
C       For sparse guesses here, figure out which density files to use in this pass.
        Call SelDFi(IOut,DoPrnt,RHF,ISCF,NSCF,DoSCF,IPOA1,IPOB1,IRwPOA,
     $    IRwPOB,SEEn,GesEn,IPOAUs,IPOBUs,V(IndX),LenX)
        Call AClear(LFXYZ,V(IFXYZ))
        VShftA = Zero
        VShftB = Zero
        If(DoCycO) then
          Call CycOpn(IOut,Energy,ECoul,Acurcy,MaxCyc,ROHF,JCycle,IExtp,
     $      NAtoms,NAEUse,NBEUse,IAn,IAtTyp,AtmChg,C,IAtFrg,NOpUse,NOp1,
     $      NOp2U,V(INEq),NBasUs,NBD,NBsUse,NBs6Us,NTT,V(IIJ),V(ID),
     $      V(IF),V(IFF),V(IV),V(IDA),V(IDB),V(IFA),V(IFB),V(IV),V(IVB),
     $      V(IVV),V(IBuf2E),ISort,IRSort,IBfSrt,MDV,LScr2,V(IScr2),
     $      IDIIS,IFBrod,V(IBMat),IPrint,IPrPFS,ZDOUse,DirSCF,IncFoc,
     $      V(INEqSh),V(INEqS2),RotAll,NEqAll,JSym2E,JSym2M,V(IIJMap),
     $      MemInt,IFinIt,ICntrl,EOne,V(IndX),LenX,RHF,IRaf,IPFlag,
     $      AllowP,InStor,V(IR1),V(IR2),V(IR3),SaveMO,NFullF,VarAc1,
     $      AccFac,Sleazy,EngCon,AccMax,V(ICP),IVShft,VSGoal,VSMax,
     $      VSGpMn,ISCRF1,ItSCRF,gFac,fuFac,EPol,EField,IExCoU,IRadAn,
     $      IRanWt,IRanGd,ICorTp,ScaHFX,ScaDFX,IDiag,NPassI,IRwEig,
     $      IRwCA,IRwCB,IRwPA,IRwPB,IRwFA,IRwFB,IOSUse,IRwHUs,IRwTM,
     $      IRwC1,IRwC2,IRwC3,ISyOcc,IWfnSy,NBlks,V(INBF),V(INBFU),
     $      V(INOcSy),ICtDFT,FMM,RComp,FMFlag,FMFlg1,NFxFlg,JustJ,
     $      IOCMOS,AbDiag,IFrcSy,ICharg,Multip,JustEn,EPert,SpFock,
     $      SpDiag,NBT,ZThr,MxSCFC,TempSt,NSAStep,Epsi,RdPCM,IHMeUs,
     $      V(jIPBC),LIPBC,CelVec,ENRUse,VShftA,VShftB,LvFirst,IStab,
     $      IRwNCA,IRwNCB,V(IIndFk),V(IEnCyc),V(IErCyc),FixRis,LvWFON,
     $      LvRise,AccDes,MaxTyp,DLaMax,V(ILLim),V(IULim),IOpDME,SECut,
     $      IPOAUs,IRwV,DThr,IProgC,DCut,IFxFM,IMemOp,IOpDM2,IPOBUs,
     $      SymmP,IRwFOA,IRwFOB,NKPnt,V(iNBFK),V(jICoor),V(jWKPnt),
     $      V(jCmlxF),IRwCAC,IRwCBC,IRwEiC,NAEBr,NOABr,NBEBr,NOBBr,
     $      MaxDamp,IDamp,DampCut,DampLim,IDFit,IOFitX,IOFDen,KAlg,
     $      FkGues,NExtDmp,IOGesA,IOSymA,IRwWgt,V(jIOToK),NTot,IADMP,
     $      BndOnl,IGWInf,V(IAtBtD),V(IIRBtD),V(IRRBtD),V(IGWBtD),
     $      FrcNuc,Halt,IZMthd,FKTemp,kTLim,NSPnt,If0K,TigPCV,IRwPAN,
     $      IRwPBN,IRANP0,ThrSDi,Omega,LChg,NumChg,V(jChg),IDoGrd,LFXYZ,
     $      V(IFXYZ),IOpCl,IRwSO,DoUpd,IRWUR,TstSpI,IStuFl,DoNDXC,ErfP,
     $      DoDMS,IADBas,LSEAll,V(ISEAll),V(IndY),LenY,NSAOBk,NS,NS2,
     $      NS4,V(jISABF),V(jRSABF),MulTab,V(jMap),FrcSDg,DISPCM,IAprXC,
     $      ENJ,NDBF6D,ICampK,CKThr,NOnlyS,NFock,DoHirS,ISwCut,IPsCut,
     $      IFrgTp,V(IMpF2M),V(IIndFr),EIter1,IAlDNI,AlGsMn)
        else if(DoStCy) then
          Call GauErr('StuCyc NYI')
C         Call StuCyc(IOut,IPrint,IOpCl,IStuFl,MaxCyc,NBasis,NBsUse,NAE,
C    $      NBE,IRwS,IRwT,IRwH,IRwEig,IRwCA,IRwCB,IRwPA,IRwPB,IRwPT,
C    $      IRwFA,IRwDip,JCycle,Acurcy,AccMax,ENR,EOne,Energy,V(IDA),
C    $      V(IDB),V(IFA),V(IR1),V(IR2),V(IR3),V(IndX),LenX)
C         StuCyc resets IOpCl.
          Call DecoSC(IOpCl,NSpBlk,NRI,NDimBl,CBfn,NSpBlX,Spinor)
          Call ILSW(1,1,IOpCl)
          RHF = .True.
          NBD = NBasis*NDimBl
          NBasUD = NBsUse*NDimBl
          endIf
        If(ISCF.eq.ISEGes) then
          Call DoProj(IOut,IPrint,V(IndX),MDV,V(ILLim),V(IULim),IAn,
     $      AtmChg,C,IAtTyp,IRWPa,IPOAUs,IRWSYM,IRWEIG,0,NAE,NBE,
     $      NAtoms,NBasis,NOrbs,IMemOp)
        else if(ISCF.eq.ISEEn) then
          SEEn = Energy
        else if(ISCF.eq.IGesEn) then
          GesEn = Energy
        else if(ISCF.eq.IVac) then
          EVac = Energy
          endIf
        If(IDoGrd.eq.1) Call IncFX(IOut,IPrPFS,'Forces','L502',IOEDrv,
     $    IRwFX,NAtomP,V(IFXYZ+LFXYZ),V(IFXYZ))
        If(ISCF.eq.IVac.and.CFPCM(ISCRF).and.IRwHUs.ne.IRwHXC)
     $    Call FreeRW(IRwHUs)
        Call TStamp(1,'After CycOpn')
  100   Continue
C
C     Test CPHF.
C
C     MaxCCP = Max(MaxCyc,128)
      If(FrcJEn) goto 999
      If(DoStCP) Call GauErr('StuCP NYI.')
C     If(DoStCP) Call StuCP(IOut,IPrint,IOpCl,IStuFl,MaxCCP,NBasis,
C    $  NBsUse,NAE,NBE,IRwS,IRwDip,IRwEig,IRwCA,IRwCB,IRwPA,IRwPB,IRwFA,
C    $  IRwFB,V(IDA),V(IDB),V(IFA),V(IFB),V(IR1),V(IR2),V(IR3),V(IndX),
C    $  LenX)
C
C     Save back timing information for DFT parallelization.
C
      If(IGWInf(1).gt.0.and.IGWInf(21).ne.0.and.UpdTim) then
        Call FileIO(1,-IOXCGd,InToWP(IGWInf(1)),IGWInf,0)
        Call FileIO(1,IOXCGd,IGWInf(12)*IGWInf(7),V(IRRBtD),
     $    IRRBtD-IAtBtD)
        endIf
C
C     Make the orbitals conform the the largest Abelian subgroup.
C
      If(IStab.eq.2) then
        If(.not.Halt) then
          Call ILSW(1,6,0)
        else
          Call FileIO(2,-IRwCS,NBNU*NSpBlk,V,0)
          Call FileIO(1,-IRwCA,NBNU,V,0)
          Call FileIO(1,-IRwCB,NBNU*(NSpBlk-1),V,0)
          endIf
        endIf
C     ISymMO = 0/1/2 leave as is/default/compute orb syms
      If(SpDiag.or.IAbSym.eq.2.or.NBlks.eq.1) then
        ISymMO = 1
      else if((IDiag.eq.1.or.IDiag.eq.3).and.IAbSym.ne.3) then
        ISymMO = 0
      else if(PChek4(IOut,IPrint,NBasis,NBas6D,V,NGot)) then
        ISymMO = 2
      else
        Write(IOut,1500)
        ISymMO = 1
        endIf
      If(ISymMO.eq.1) then
        Call DefSym(NBasUD,1,IOSymA,V(IF))
      else if(ISymMO.eq.2) then
        IISym = 1
        IX = IISym + InToWP(NBasis)
        LenX = NGot - IX + 1
        Call ConDDF(IOSymA,2*NBsUse)
        If(IAbSym.eq.3) then
          IPrAbS = Max(IPrint,1)
        else
          IPrAbS = IPrint
          endIf
        Call AbSyMO(IOut,IPrAbS,.False.,NRI,NDimBl,NBasis,NBsUse,NBasUD,
     $    1,IRwCA,IOSymA,V(IISym),V(IX),LenX)
        If(NSpBlk.eq.2) Call AbSyMO(IOut,IPrint,.False.,NRI,NDimBl,
     $    NBasis,NBsUse,NBasUD,2,IRwCB,IOSymA,V(IISym),V(IX),LenX)
        endIf
      Call TStamp(1,'After AbSyMO')
C
C     Convergence has been met.  Save the densities.
C
      If(SpDiag) then
        ISpX = IV
        INext = ISpX + NBasis
        Call TstCor(INext,MDV,'UHFOpn-SpPT')
        If(RHF) then
          Call SpGet(V,V,NGot,IRwPA,INext,JIP,IPA,IStart,NBasis,NZP,1,
     $      .True.,0)
          IPB = IPA
        else
          Call MkSmFm(V,NGot,IRwPA,0,0,IRwPB,0,0,INext,INext,0,JIP,IPA,
     $      IPB,NZ,IStart,NBasis,2)
          endIf
        If(ZDO) Call SpStor(IRwPOA,V(JIP),V(ipa),NBasis,1,.true.,0)
        Call SpStor(IRwPA,V(JIP),V(IPA),NBasis,1,.True.,0)
        Call SpStor(IRwPB,V(JIP),V(IPB),NBasis,1,.True.,0)
        If(IPrint.ge.3) then
          Call SPDisp(IOut,v(JIP),v(ipa),NBasis,'PA:',2)
          Call SPDisp(IOut,v(JIP),v(ipb),NBasis,'PB:',2)
          endIf
        Call SpAddM(V,MDV,V(JIP),V(ipa),V(JIP),V(ipb),V(ISpX),IStart,
     $    jipt,ipt,INext,NBasis,NZC,One,One,-One,2)
        If(IPrint.ge.3) Call SPDisp(IOut,v(jipt),v(ipt),NBasis,'PT:',2)
        Call SpStor(IRwPT,V(jipt),V(ipt),NBasis,1,.true.,0)
        Call SpAddM(V,MDV,V(JIP),V(ipa),V(JIP),V(ipb),V(ISpX),IStart,
     $    jips,ips,INext,NBasis,NZC,One,-One,-One,2)
        If(IPrint.ge.3) Call SPDisp(IOut,v(jips),v(ips),NBasis,'PS:',2)
        Call SpStor(IRwPS,V(jips),V(ips),NBasis,1,.true.,0)
        S2 = Zero
        If(FixT) Call MkSprL(V,NGot,IRwT,0,IRwT,1,jia,ian1,iend1,NBasis,
     $    NZA,'T Matrix',zthr,IPrint)
        Call Viral2(V,NGot,IRwPT,IRwT,Vir,T,Energy,NBasis)
      else
        If(DoRotS) then
          Call Spin(IOpCl,IOSUse,IRwCA,IRwCB,IRwUR,NBasis,NBsUse,NAE,
     $      NBE,SVec(1),SVec(2),SVec(3),S2,S,S2Cor,V(ID),NGot-ID+1)
          Call FileIO(2,-IRwCA,NBNU,V(ID),0)
          Call RotSpn(IOut,IPrint,.False.,IOSUse,1,SVec,NBasis,NBasUD,
     $      NAE+NBE,0,V(ID),V(IF),NGot-IF+1)
          Call FileIO(1,-IRwCA,NBNU,V(ID),0)
          endIf
        Call FixPhR(NRI,NBD,NBasUD,IRwCA,IRwSMO,0,V(IFA),V(IndX),LenX)
        If(IOpCl.le.3) Call FixPhR(NRI,NBD,NBasUD,IRwCB,IRwSMO,
     $    NRI*NBD*NBasUD,V(IFA),V(IndX),LenX)
        If(UsePBC.and.UsKPnt) then
          Call FixPhC(NBasis,NKPnt,V(iNBFK),IRwCAC,V(IFA),NGot-IFA+1)
          If(IRwCBC.ne.0) Call FixPhC(NBasis,NKPnt,V(iNBFK),IRwCBC,
     $      V(IFA),NGot-IFA+1)
          endIf
        Call FileIO(2,-IRwPA,NTTX,V(IDA),0)
        If(IOpCl.le.3) then
          Call FileIO(2,-IRwPB,NTTX,V(IDB),0)
          Call AAdd(NTTX,V(IDA),V(IDB),V(IDA))
          Call FileIO(1,-IRwPT,NTTX,V(IDA),0)
          Call ACasB(NTTX,V(IDA),V(IDB),V(IDA),-Two)
          Call FileIO(1,-IRwPS,NTTX,V(IDA),0)
        else
          Call FileIO(1,-IRwPT,NTTX,V(IDA),0)
          Call AClear(NTTX,V(IDA))
          Call FileIO(1,-IRwPS,NTTX,V(IDA),0)
          endIf
        If(RHF.and.IOpCl.eq.0) then
          S2 = Zero
        else
          Call Spin(IOpCl,IOSUse,IRwCA,IRwCB,IRwUR,NBasis,NBsUse,NAE,
     $      NBE,SVec(1),SVec(2),SVec(3),S2,S,S2Cor,V(ID),NGot-ID+1)
          If(UsKPnt) Call ZSpin(IOut,IPrint,V(jIPBC),NBasis,NKPnt,
     $      V(iNBFK),V(jICoor),V(jCmlxF),V(jWKPnt),IRwS,IRwSk,IRwCAC,
     $      IRwCBC,S2,V(IndX),LenX)
          endIf
        If(DoCyc) Call Virial(IOpCl,NBasis,Energy,V(IDA),V(IDB),IRwPT,
     $    IRwT,Vir,T)
        endIf
      Call TStamp(1,'After Virial')
C
      If(NDimBl.eq.2.and.ITqry(IRwSO).gt.0) then
        Call FileIO(2,-IRwPT,NTTX,V(IDA),0)
        Call FileIO(2,-IRwSO,NTTX,V(IFA),0)
        ELS = SCFTrX(.False.,.False.,IOpCl,NMtPBC,NBasis,V(IDA),V(IDA),
     $    V(IFA),V(IFA))
      else
        ELS = Zero
        endIf
      If(DoCyc) then
        Write(IOut,1110) MethNm(1:LinEnd(MethNm)), Energy, JCycle,
     $    NFock, Acurcy, Vir
        If(IOpCl.ne.0.and.IOpCl.ne.2.and.IOpCl.ne.10.and..not.SpDiag)
     $    Write(IOut,1120) SVec, S2, S, ELS
        If(Onsagr(ISCRF)) then
          EPol = -EPol / GFloat(2)
          Call RF2Gen(IOut,Increm,ICharg,gFac,Energy,EPol,EField,DGen)
          If(ISCRFT.gt.1) Call FileIO(1,-IOFInF,35,EField,0)
          endIf
        PE = EOne - T
        EE = Energy - T - PE - ENR
        If(IPSave(0).eq.0) Write(IOut,1340) T, PE, EE
        If(PCM(ISCRF)) then
          Call PCMGet('NSph',NSph,RJunk)
          Call PCMGet('NSphG',NSphG,RJunk)
          Call PCMGet('NSphNE',NSphNE,RJunk)
          jNord  = 1
          jNEord = jNord  + InToWP(NSphG)
          jContS = jNEord + InToWP(NSphNE)
          jContA = jContS + 7*(NAtoms+1)
          jGElSp = jContA + 4*NAtoms
          jEnd   = jGElSp + Min(NSph,(NSphG+1)) - 1
          Call TstCor(jEnd,NGot,'UHFOPN-PCMPr1')
          ENRvac = ENR - GNEl
          Evac1  = Evac - ENRvac
          EPert1 = EPert - ENR
          Energ1 = Energy - ENR
          Energ2 = Energy - GNEl
          Call PCMPr1(IOut,IPrint,2,Evac1,EPert1,Energ1,ENRvac,Energ2,
     $      NAtoms,IAn,IAtTyp,V(jNord),V(jNEord),V(jContS),V(jContA),
     $      V(jGElSp))
          Call PCMGet('ICosTh',ICosTh,RJunk)
          If(ICosTh.eq.1) Call WrtCos(IOut,IPrint,IPFlag,AllowP,Energ1,
     $      ENRvac,NAtoms,IAn,C,V,NGot)
          Call PCMIOp('save',' ',IJunk,RJunk)
          endIf
C
C       Attempt to perform annihilation of the largest spin contaminant.
C
        SS4 = Zero
        S2ARet = Zero
        If(.not.RHF.and..not.SpDiag) Call Annil(IOut,IPrint,IOpCl,IUHFT,
     $    NAE,NBE,NBasis,NBsUse,S2ARet,SS4,V(IDA),V(IFA),V(IndX),LenX)
        Call TStamp(1,'After Annil')
C
C       Store calculated energies.
C
        DGen(1) = Vir
        DGen(5) = ECoul
        DGen(23) = Acurcy
        DGen(24) = S2ARet
        DGen(29) = SS4
        DGen(32) = Energy
        If(UseNuc) then
          DGen(43) = DGen(41)
        else
          DGen(43) = Energy
          endIf
        DGen(44) = S2
        DGen(56) = ErfP
        DGen(60) = EIter1
      else
        IOGenC = FilNum(IRwGen,IUnit(9))
        If(ITqry(IOGenC).gt.0) then
          Call FileIO(2,-IOGenC,LRwGen,DGenC,0)
          DGen(1) = DGenC(1)
          DGen(5) = DGenC(5)
          DGen(23) = DGenC(23)
          DGen(24) = DGenC(24)
          DGen(29) = DGenC(29)
          DGen(32) = DGenC(32)
          DGen(43) = DGenC(43)
          DGen(44) = DGenC(44)
          DGen(55) = DGenC(55)
          DGen(56) = DGenC(56)
          DGen(60) = DGenC(60)
          endIf
        endIf
      Call FileIO(1,-IRwGen,LRwGen,DGen,0)
C
C ANT BEGIN modification
C
C         Perform the conductance calculation at end of scf (1000)
C
      call ANT (.NOT.RHF,1000,IRwH,IRwPA,IRwPB,IRwFA,IRwFB,
     $                IRwS,IRwEig,Acurcy,Crit,.true.,NBasis)
C
C ANT END modification
C
C
C     Print requested output.
C
      Call UHFPrt(IOut,IPrint,SpDiag,IOpCl,IRwEig,IRwCA,IRwCB,IRwPA,
     $  IRwPB,IRwPT,IRwPS,IRwUR,NOATot,NOBTot,NBasis,NBsUse,VShftA,
     $  VShftB,V(IScr2),V(ID))
C
C     Dump out the PBC eigenvalues if requested
      If(DoPrnt.and.UsKPnt) then
        LenEi1 = NBasis*NKPnt
        Call TstCor(ID+NSpBlk*(LenEi1+NBasis)+2*NBsUse,NGot,
     $    'UhfOpn-DmpEig')
        Call FileIO(2,-IRwEig,NBsUse,V(ID),0)
        Call FileIO(2,-IRwEiC,LenEi1,V(ID+NBasis),0)
        If(NSpBlk.eq.2) then
          Call FileIO(2,IRwEig,NBsUse,V(ID+LenEi1+NBasis),0)
          Call FileIO(2,IRwEiC,LenEi1,V(ID+LenEi1+2*NBasis),0)
          endIf
        Call DmpEig(IOut,IPrint,NPrnKE,NBasis,NBsUse,NSpBlk-1,NAE,NBE,
     $    NKPnt,V(ID),V(iNBFK),V(jKCoor))
        endIf
C
C     Calculate the energy-weighted density matrix for possible
C     gradient links later in the route.
C
      If(RHF.and.IOpCl.le.3) then
        If(SpDiag) then
          Call SpGet(V,V,NGot,IRwPT,1,jipt,ipt,iend,NBasis,NZPT,1,
     $      .true.,0)
          Call SpStor(IRwPA,V(jipt),V(ipt),NBasis,1,.true.,0)
        else
          Call FileIO(2,-IRwPT,NTTX,V(IDA),0)
          Call FileIO(1,-IRwPA,NTTX,V(IDA),0)
          endIf
        FactW = Pt5
      else
        FactW = One
        endIf
      Call ILSW(2,5,Kill)
      Call TStamp(1,'Before W')
      If(UsKPnt) then
        Call ConDDF(IRwW,NTTX)
        If(UseGEW) then
          UseLow = IADBas.eq.1
          IFAO  = IndX
          IDAO  = IFAO  + NBSqK
          IVO   = IDAO  + NBSqK
          IScr2 = IVO   + NBSqK
          IScr3 = IScr2 + NBNK
          IScr4 = IScr3 + NBSqK
          Call TstCor(IScr4,LenX,'UHFOpn-ZGEWDn')
          Call ConDDF(IRwFRA,NBSqK)
          Call ConDDF(IRwFRB,NBSqK)
          Call ZGEWDn(NBasis,NTTX,NAE,NBE,.False.,RHF,One,UseLow,UsePBC,
     $      UsKPnt,NKPnt,NMtPBC,NPDir,IRwCA,IRwCB,IRwPA,IRwPB,IRwPAN,
     $      IRwPBN,IRwFA,IRwFB,IRwFRA,IRwFRB,IRwS,IRwV,IRwW,V(IFAO),
     $      V(IDAO),V(IVO),V(IScr2),V(iNBFK),V(jIPBC),V(jICoor),
     $      V(jCmlxF),V(IScr3),V(IScr4),LenX-IScr4)
        else
          Call ZFormD(IOut,IPrint,V(jIPBC),.True.,NSpBlk-1,NBasis,11,
     $      IRwCAC,IRwCBC,IRwW,IRwW,IRwEiC,IRwWgt,NKPnt,V(iNBFK),
     $      0,0,XXA,.False.,V(jICoor),V(jCmlxF),V(IndX),LenX)
          endIf
      else if(SpDiag.and..not.ZDO) then
        If(IADMP.gt.0) then
          If(RHF) FactW = One
          Call SpGEWD(V,NGot,IRwS,IRwV,IRwFA,IRwPA,IRwW,FactW,NBasis,
     $      ZThr,IPrint)
        else
          Call FormWSp(IOut,IPrint,NSpBlk-1,V,NGot,IRWPA,IRwPB,IRWFA,
     $      IRwFB,IRwW,NBasis,ZThr)
          endIf
      else if(.not.ZDO) then
        If(UseGEW) then
          If(RComp) Call GauErr('GEWDen and RComp not allowed')
          If(NRI.eq.2) Call GauErr('No GEW for complex yet.')
          Call GEWDen(IOut,NBasis,NBsUse,NAE,NBE,IRwCA,IRwCB,IRwPA,
     $      IRwPB,IRwFA,IRwFB,IRwS,IRwV,IRwRAO,IRwW,V(IFA),V(IDA),V(IV),
     $      V(IScr2),.False.,RHF,One,IADBas,V(IV+NBSq),NGot-IV-NBSq+1)
        else
          Call FormW(.True.,IRwPA,IRwFA,IRwW,NRI,NBD,FactW,V(IFA),
     $      V(IDA),V(IV))
          If(.not.RHF) Call FormW(.False.,IRwPB,IRwFB,IRwW,NRI,NBD,
     $      One,V(IFA),V(IDA),V(IV))
          endIf
        If(IPrint.ge.2) Call LTOutX(IOut,0,IOpClX,1,0,'-W:',NBasis,
     $    NBasis,V(IF),V(IF),0)
        endIf
      Call TStamp(1,'After W')
C
C     For ADMP and SimOpt, transform Fock matrices to orthogonal basis.
C
      If(IADMP.gt.0) Call ADMFin(IOut,SpDiag,DoPrnt,UsKPnt,IADMP,IADBas,
     $  ZDO,RHF,NAE,NBE,NBasis,NBsUse,NTTX,NMtPBC,NKpnt,NPDir,IRwV,
     $  IRwPOA,IRwPOB,IRwPA,IRwPB,IRwFOA,IRwFOB,IRwFA,IRwFB,IRwFRA,
     $  IRwFRB,V(jIPBC),V(jICoor),V(iNBFK),V(jCmlxF),V(JD),V(JVX),
     $  V(JVV),V(JVY),MDVY)
C
C     For SimOpt, prepare wavefunction DIIS and E-DIIS matrices
C
      If(DoSim) Call FOSimF(IOut,IPrint,RHF,IADBas,NBasis,IRwSM,IRwV,
     $  IRwPA,IRwPB,IRwPOA,IRwPOB,IRwFA,IRwFB,IRwFRA,IRwFRB,
     $  V(JD),V,NGot)
C
      If(Kill.eq.1.and.Halt)
     $  Call GauErr('Convergence failure -- run terminated.')
      Call FreeRW(IRwC1)
      Call FreeRW(IRwC2)
      Call FreeRW(IRwC3)
      Call FreeRW(IRwFOA)
      Call FreeRW(IRwFOB)
      If(.not.SavPBC) then
        Call FreeRW(IRwCAC)
        Call FreeRW(IRwCBC)
        Call FreeRW(IRwEiC)
        endIf
C
C     Duplicate density matrices for PBC testing.
C
      If(NMtPBC.gt.1.and..not.UsKPnt) then
        Call TstCor(NMtPBC*NTT,NGot,'UHFOpn-Dup')
        Call DupMat(NMtPBC,NTT,IRwPA,NMtPBC,V)
        Call DupMat(NMtPBC,NTT,IRwPB,NMtPBC,V)
        Call DupMat(NMtPBC,NTT,IRwPT,NMtPBC,V)
        Call DupMat(NMtPBC,NTT,IRwPS,NMtPBC,V)
        Call DupMat(NMtPBC,NTT,IRwW,NMtPBC,V)
        endIf
C
C     Possibly flag the run as UHF for testing.
C
      If(IUHFT.eq.1) then
        Write(IOut,1600)
        Call ILSW(1,1,1)
        L1 = InToWP(NBsUse)
        Call FileIO(2,-IOSymA,L1,V,0)
        Call ConDDF(IOSymA,2*L1)
        Call FileIO(1,-IOSymA,L1,V,0)
        Call FileIO(1,IOSymA,L1,V,0)
        Call FileIO(2,-IRwPA,NTT,V,0)
        Call AScale(NTT,Pt5,V,V)
        Call FileIO(1,-IRwPA,NTT,V,0)
        Call FileIO(1,-IRwPB,NTT,V,0)
        endIf
C
C     GF testing.
C
      If(ITstGF.ne.0) Call GFTest(IOut,IPrint,ITstGF,NBas6D,NBasis,
     $  NBsUse,NAE,NBE,IPFlag,AllowP,V,MDV)
  999 Call ChainX(0)
      Return
      End
*Deck ADMFin
      Subroutine ADMFin(IOut,SpDiag,DoPrnt,UsKPnt,IADMP,IADBas,ZDO,RHF,
     $  NAE,NBE,NBasis,NBsUse,NTTX,NMtPBC,NKPnt,NPDir,IRwV,IRwPOA,
     $  IRwPOB,IRwPA,IRwPB,IRwFOA,IRwFOB,IRwFA,IRwFB,IRwFRA,IRwFRB,IPBC,
     $  ICoor,NBFK,CmlxF,D,VX,VV,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Transform the Fock matrices to orthogonal basis for ADMP.  At the
C     first point we must also transform the density matrices to the
C     orthonormal basis.
C
      Logical SpDiag, DoPrnt, ZDO, RHF, UsKPnt, Error
      Complex*16 CmlxF(3,*)
      Dimension D(NBasis,NBasis),VX(NBasis,NBasis),VV(NBasis),V(MDV),
     $  ICoor(3,*),XXA(1),IPBC(*),IBckDim(3),NBFK(NKPnt,*)
      Save One, Two, XXA, Small
      Data One/1.0d0/, Two/2.0d0/, XXA/0.0d0/, Small/1.d-8/
 1000 Format(' Tr(PA) in orthogonal basis=',F12.6)
 1010 Format(' Tr(PB) in orthogonal basis=',F12.6)
C
      NTT = (NBasis*(NBasis+1))/2
      NBSq = NBasis**2
      NBSqK = 2*NBSq*NKPnt
      Call ConDDF(IRwFRA,NBSqK)
      Call ConDDF(IRwFRB,NBSqK)
      If(UsKPnt) then
        jStrt = 1
        Call ExtrCl(2,NPDir,NMtPBC,XXA,IPBC,IBckDim,iInCl,jInCl,jCntCl,
     $    iBckCl,jStrt,V,MDV)
        NReal = max(IBckDim(1),IBckDim(2),IBckDim(3))
        iAngAr = jStrt
        iKpnt = iAngAr + 6*(2*NReal+1)
        Call TstCor(iKpnt+NMtPBC,MDV,'ADMFin-ExtrCl')
        Call FileIO(2,-IRwPA,NTTX,D,0)
C       Read in current real Fock matrix
        Call FileIO(2,-IRwFA,NTTX,D,0)
        Call FileIO(2,-IRwV,NBSqK,VX,0)
        Call ARe2Or(NBasis,NKPnt,NMtPBC,NReal,IRwFRA,NBFK,V(iInCl),
     $    ICoor,VV,VX,CmlxF,V(iAngAr),V(iKpnt),D,VV(2*NBSq+1))
        If(.not.RHF) then
          Call FileIO(2,-IRwFB,NTTX,D,0)
          Call ARe2Or(NBasis,NKPnt,NMtPBC,NReal,IRwFRB,NBFK,V(iInCl),
     $      ICoor,VV,VX,CmlxF,V(iAngAr),V(iKpnt),D,VV(2*NBSq+1))
          endIf
      else
        If(SpDiag) then
          If(RHF) Call SpFCpy(RHF,NBasis,IRwPOA,0,IRwPOA,0,Two,One,V,
     $      MDV)
          If(ZDO) then
            Call SpFCpy(RHF,NBasis,IRwFA,IRwFB,IRwFRA,IRwFRB,One,One,V,
     $        MDV)
          else
            Call SpFCpy(RHF,NBasis,IRwFOA,IRwFOB,IRwFA,IRwFB,One,One,V,
     $        MDV)
            Call SpFCpy(RHF,NBasis,IRwFOA,IRwFOB,IRwFRA,IRwFRB,One,One,
     $        V,MDV)
            endIf
        else
          IX = 1
          IT = IX + NBSq
          IScr = IT + NBSq
          MDV1 = MDV - IScr + 1
          Call TstCor(IScr+NBasis,MDV,'ADMFin-OTran')
          If(ZDO) then
            Call Trsfr(NTT,IRwFA,IRwFRA,V,MDV)
            If(.not.RHF) Call Trsfr(NTT,IRwFB,IRwFRB,V,MDV)
          else
            Call FileIO(2,-IRwV,NBSq,V(IX),0)
            Call FileIO(2,-IRwFA,NTT,D,0)
            Call OTran(D,V(IX),V(IT),NBasis,NBsUse,0,IADBas,V(IScr),
     $        MDV1)
            Call FileIO(1,-IRwFRA,NTT,D,0)
            Call FileIO(1,-IRwFA,NTT,D,0)
            If(.not.RHF) then
              Call FileIO(2,-IRwFB,NTT,D,0)
              Call OTran(D,V(IX),V(IT),NBasis,NBsUse,0,IADBas,V(IScr),
     $          MDV1)
              Call FileIO(1,-IRwFRB,NTT,D,0)
              Call FileIO(1,-IRwFB,NTT,D,0)
              endIf
            endIf
          If(IADMP.eq.2) then
            Call FileIO(2,-IRwPA,NTT,D,0)
            If(.not.ZDO) then
              Call FileIO(2,-IRwV,NBSq,V(IX),NBSq)
              Call OTran(D,V(IX),V(IT),NBasis,NBsUse,1,IADBas,V(IScr),
     $          MDV1)
              endIf
            TrPA = TraLT(1,NBasis,D)
            If(RHF) then
              RNE = GFloat(2*NAE)
            else
              RNE = GFloat(NAE)
              endIf
            Error = Abs(TrPA-RNE).ge.Small
            If(DoPrnt.or.Error) Write(IOut,1000) TrPA
            If(Error) Call GauErr('Bad alpha density matrix.')
            Call ConWrt(IRwPOA,NTT,D)
            If(.not.RHF) then
              Call FileIO(2,-IRwPB,NTT,D,0)
              If(.not.ZDO) Call OTran(D,V(IX),V(IT),NBasis,NBsUse,1,
     $          IADBas,V(IScr),MDV1)
              TrPB = TraLT(1,NBasis,D)
              Error = Abs(TrPB-GFloat(NBE)).ge.Small
              If(DoPrnt.or.Error) Write(IOut,1010) TrPB
              If(Error) Call GauErr('Bad beta density matrix.')
              Call ConWrt(IRwPOB,NTT,D)
              endIf
            endIf
          endIf
        endIf
      Return
      End
*Deck ADMIni
      Subroutine ADMIni(IOut,SpDiag,DoPrnt,UsKPnt,IADMP,NMtPBC,NKPnt,
     $  NPDir,NTTX,IADBas,ZDO,RHF,NBasis,NBsUse,IRwS,IRwPOA,IRwPOB,
     $  IRwPA,IRwPB,IRwRAO,IRwV,IRwPAN,IRwPBN,IPBC,ICoor,NBFK,CmlxF,
     $  WKPnt,D,DD,VX,VV,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Initialization of the orthonormal transformation and possibly
C     initial density for ADMP.
C
      Logical SpDiag, DoPrnt, ZDO, RHF, UsKPnt
      Complex*16 CmlxF(3,*)
      Dimension D(NBasis,NBasis), DD(NBasis), VX(NBasis,NBasis),
     $  VV(NBasis),V(MDV),NBFk(NKPnt,*),WKPnt(*),ICoor(3,*),XXA(1),
     $  IPBC(*),IBckDim(3)
      Save Half, One, XXA
      Data Half/0.5d0/, One/1.0d0/, XXA/0.0d0/
 1000 Format(' Using Lowdin orthonormal basis.')
 1010 Format(' Using Cholesky orthonormal basis.')
 1020 Format(' Using SVD-based orthonormal basis.')
C
      If(SpDiag) Call GauErr('Sparse and ADMP NYI.')
      NTT = (NBasis*(NBasis+1))/2
      NBSq = NBasis**2
      NBSqK = 2*NBSq*NKPnt
      NTTK = 2*NTT*NKPnt
      If(UsKPnt) then
        IRwPAN = IGetRW(NTTK,0)
        IRwPBN = IGetRW(NTTK,0)
        If(IADMP.eq.1) then
          jStrt = 1
          Call ExtrCl(2,NPDir,NMtPBC,XXA,IPBC,IBckDim,iInCl,jInCl,
     $      jCntCl,iBckCl,jStrt,V,MDV)
          NReal = Max(IBckDim(1),IBckDim(2),IBckDim(3))
          iAngAr = jStrt
          iKpnt = iAngAr + 6*(2*NReal+1)
          Call TstCor(iKpnt+NMtPBC,MDV,'CPFIni-ExtrCl')
          Call FileIO(2,-IRwPOA,NTTK,D,0)
          Call FileIO(2,-IRwV,NBSqK,VX,0)
          Call AOr2Re(NBasis,NKPnt,NMtPBC,NReal,IRwPAN,NBFK,WKPnt,
     $      V(iInCl),ICoor,D,VX,CmlxF,V(iAngAr),V(iKpnt),DD,VV,
     $      VV(2*NBSq+1))
          Call FileIO(1,-IRwPA,NTTX,DD,0)
          If(RHF) then
            Call FileIO(1,-IRwPB,NTTX,DD,0)
          else
            Call FileIO(2,-IRwPOB,NTTK,D,0)
            Call AOr2Re(NBasis,NKPnt,NMtPBC,NReal,IRwPBN,NBFK,WKPnt,
     $        V(iInCl),ICoor,D,VX,CmlxF,V(iAngAr),V(iKpnt),DD,VV,
     $        VV(2*NBSq+1))
            Call FileIO(1,-IRwPB,NTTX,DD,0)
            endIf
        else
          Call ConDDF(IRwPOA,NTTK)
          If(.not.RHF) Call ConDDF(IRwPOB,NTTK)
          endIf
      else
C
C       Form the transformation matrices to/from orthonormal
C       (IRwV contains U**-1 or S**(-1/2) first, followed by U or S**(1/2))
C
        If(.not.SpDiag) then
          If(.not.ZDO) then
            Call ConDDF(IRwV,2*NBSq)
            If(IADBas.eq.1) then
              If(DoPrnt) Write(IOut,1000)
              Call FormV(2,NBasis,NBasis,IRwS,IRwV,D,VX,DD,VV)
            else if(IADBas.eq.2) then
              If(DoPrnt) Write(IOut,1010)
              Call FileIO(2,-IRwS,NTT,D,0)
              Call CholIn(NBasis,D,VX)
              Call FileIO(1,-IRwV,NBSq,VX,0)
              Call FileIO(1,IRwV,NBSq,D,0)
            else if(IADBas.eq.3) then
              If(DoPrnt) Write(IOut,1020)
              Call FileIO(2,-IRwRAO,NBasis*NBsUse,VX,0)
              Call AClear(NBasis*(NBasis-NBsUse),VX(1,NBsUse+1))
              Call FileIO(1,-IRwV,NBSq,VX,0)
            else
              Call GauErr('Illegal IADBas in ADMIni.')
              endIf
            endIf
          endIf
C
C       After the first point, we must transform the densities produced by l121
C       to the AO basis.
C
        If(IADMP.eq.1) then
          If(SpDiag) then
            If(ZDO) then
              If(RHF) then
                ScaleA = Half
              else
                ScaleA = One
                endIf
              Call SpFCpy(RHF,NBasis,IRwPOA,IRwPOB,IRwPA,IRwPB,ScaleA,
     $          One,V,MDV)
            endIf
          else
            IX = 1
            IT = IX + NBSq
            IScr = IT + NBSq
            MDV1 = MDV - IScr + 1
            Call TstCor(IScr+NBasis,MDV,'ADMIni-OTran')
            If(.not.ZDO) Call FileIO(2,-IRwV,NBSq,V(IX),0)
            Call FileIO(2,-IRwPOA,NTT,D,0)
            If(.not.ZDO) Call OTran(D,V(IX),V(IT),NBasis,NBsUse,1,
     $        IADBas,V(IScr),MDV1)
            Call FileIO(1,-IRwPA,NTT,D,0)
            If(.not.RHF) then
              Call FileIO(2,-IRwPOB,NTT,D,0)
              If(.not.ZDO) Call OTran(D,V(IX),V(IT),NBasis,NBsUse,1,
     $          IADBas,V(IScr),MDV1)
              Call FileIO(1,-IRwPB,NTT,D,0)
              endIf
            endIf
          endIf
        endIf
      Return
      End
*Deck AExpnd
      Subroutine AExpnd(N,IDim1,NDim1,NDim2,A,B)
      Implicit Real*8(A-H,O-Z)
C
C     Copy N elements of each of NDim2 columns of A into B.
C
      Dimension A(IDim1,NDim2), B(NDim1,NDim2)
C
      Do 10 J = 1, NDim2
        Do 10 I = 1, N
   10     B(I,J) = A(I,J)
      Return
      End
*Deck CEMCyc
      Subroutine CEMCyc(IH,H,CC,C,TT,PRODHT,RT,XX,IX,IMap,JMap,NOcc,
     $  xmax,xmin,N,NZ,NP,NumC,Lfail,numbig,numsml,B,IOut,IPrint,
     $  StoreT,LenRT,Ademp,BnFac,IBin,IBT,BT,IScr,ISizCV,NBin,Method,
     $  IFunc)
      Implicit Real*8(A-H,O-Z)
C
C     This code determines the upper and lower limits for eigenvalues
C     of the Hamiltonian
C
      Logical StoreT,fail
      Parameter(MaxIt=500)
      Dimension XX(N*NP,2),IX(N*NP,3),JMap(*),IMap(*),IH(N+1,2),CC(*),
     $  PRODHT(N,2),IBin(NBin+1,N,*),IBT(*),BT(*),ISizCV(*),
     $  IScr(*),BnFac(*),C(NUMC),TT(N+NZ,3),H(N,2),RT(LenRT,numc)
      Save One,Two,xacc
      Data One/1.0d0/,Two/2.0d0/,xacc/1.d-9/
C
C     Shift the Hamiltonian so that its centered at zero
      Diff = (Abs(XMin-XMax)/Two)
      SMiddle = (XMin+XMax)/Two
C     Store the original Hamiltonian in 'H' for later
C     Shift the Hamiltonian so that the middle is at zero
C$OMP Parallel Do Default(Shared) Private(I)
      Do 40 I = 1, N
   40   H(I,1) = H(I,1) - SMiddle
C
C     Scale the Hamiltonian now such that its eigenvalues span
C     from [-1,1]
      call ascale((N+NZ),One/DIFF,H,H)
C     Set initial characteristics for the Heavyside operator
c     this works       B = 1.65d0
C
C     Get the column vectors in RT
C
      If(StoreT) then
        LRT = (N+NZ)*(NumC+1)
        IRWRT = IGetRW(LRT,0)
      else
        IRwRT = 0
        endIf
      Call ChbExp(IH,H,TT,PRODHT,XX,IX,RT,IMap,JMap,N,NZ,NP,NumC,LenRT,
     $  StoreT,IRWRT,BnFac,IBin,IBT,BT,IScr,ISizCV,NBin,Method)
C
C     Newton-Raphson with bisection
C
      Ademp = EnerNR(C,ProdHT,RT,CC,B,N,NOcc,NZ,NumC,Maxit,xacc,IPrint,
     $  iout,NP,XX,LenRT,StoreT,IRWRT,H,fail,IFunc)
      If(Fail) Call GauErr('EnerNR Failed in CEMCyc.')
      Call FreeRW(IRWRT)
      Return
      End
*Deck ChbCof
      Subroutine ChbCof(XPnts,FncPts,CbCofs,aMu,Beta,n,ifunc)
      Implicit Real*8(A-H,O-Z)
C
C     Form the coeffcients for the Chebyshev polynomial expansion (CbCofs)
C       of the complementary error function (ifunc=0) or the
C       Fermi-Dirac distribution (ifunc=1)
C
C     aMu is the fermi level
C     Beta is the ficticious inverse temperature for ifunc=1, for ifunc=0 the
C        spread in the Fermi level.
C
C     N: Length of the Chebyshev expansion
C     XPnts, FncPts:  scratch arrays of length N
C
C     Chebyshev polynomial is evaluated using a uniform (linear) transformation
C     of y in [-1,1] to y in [a,b] via y = .5*(b-a)*y + .5*(b+a)
C
      dimension XPnts(n),FncPts(n),CbCofs(n)
      Save Zero,Pt5,One,Two,Four
      Data Zero/0.d0/,Pt5/0.5d0/,One/1.d0/,Two/2.d0/,Four/4.d0/
C
      Pi = Four*ATan(One)
C     Determine the necessary points to evaluate the function at.
      PIn = Pi/GFloat(N)
      Do 5 I = 1, N
    5   XPnts(I) = Cos(PIn*(GFloat(N-I)+Pt5))
C     Evaluate the Function:
      If(IFunc.eq.0) then
C     Complementary Error Function
        Do 10 I = 1, N
   10     FncPts(I) = Pt5*GErfC((XPnts(I)-aMu)*Beta)
      else
C     Fermi Dirac Distribution
        Do 20 I = 1, N
   20     FncPts(I) = One/(One+Exp(Beta*(XPnts(I)-aMu)))
        endIf
C     Chebyshev fitting to a function
      Scale = Two/GFloat(Max(N,1))
      PidN = Pi/GFloat(N)
      Do 40 I = 0, N-1
        Coef = Zero
        Do 50 J = 1, N
   50     Coef = Coef + FncPts(J)*Cos(PidN*GFloat(I)*(GFloat(N-J)+Pt5))
   40   CbCofs(I+1)=Coef*Scale
      Return
      End
*Deck ChbExp
      Subroutine ChbExp(IH,H,TT,PRODHT,XX,IX,RT,IMap,JMap,N,NZ,NP,NumC,
     $  LenRT,StoreT,IRWRT,BnFac,IBin,IBT,BT,IScr,ISizCV,NBin,Method)
      Implicit Real*8 (A-H,O-Z)
      Logical InitA, StoreT
      Real*8 JJA(1)
      Dimension IH(N+1,2),H(N,2),TT(N+NZ,3),IMAP(N+1,2),JMAP(NZ),
     $  PRODHT(*),RT(LenRT,numc),BnFac(*),XX(N*NP,*),IX(N*NP,*),
     $  IBin(NBin+1,N,*),IBT(*),BT(*),IScr(*),ISizCV(*)
      Save One,Zero,Two
      Data One,Zero,Two/1.d0,0.d0,2.d0/
C
      Call AClear(3*(N+NZ),TT)
      Call AClear(numc*LenRT,RT)
C
C     Make an identity matrix with the same form as the Hamiltonian
C
      If(Method.eq.2) call imove(N+1+NZ,IH,IScr)
      NpNZ=N+NZ
      JJA(1) = zero
      InitA = .true.
      Call ASet(N,One,TT)
C
C     Construct the second recursion matrix which is H
C
      Call AMove(NpNZ,H,TT(1,2))
C
C     Fill first 2 columns of RT with TT(*,1-2)
C
      Call AMove(LenRT,TT(1,1),RT)
      Call AMove(LenRT,TT(1,2),RT(1,2))
      If(StoreT) Call FileIO(1,-IRWRT,2*NpNZ,TT,0)
      ICnter = 2
      ICnt1 = 2
      ICnt2 = 3
      ICnt3 = 1
      Do 20 ICnter = 2, NumC
C       reg: XX(1,2), IX(1,3)
C       bin sort: XX(1,4), IX(1,6)
        Call FSMpy(IH,TT(1,ICnt1),IMap,jMap,JJA(1),JJA(1),IScr,H,JJA(1),
     $    JJA(1),IBT,BT,IH,ProdHt,IX(1,1),XX(1,1),JJA(1),JJA(1),IX(1,4),
     $    JJA(1),JJA(1),IBin(1,1,1),IBin(1,1,2),JJA(1),IX(1,2),XX(1,2),
     $    JJA(1),IX(1,5),IX(1,3),XX(1,3),BnFac,ISizCV,XX(1,4),IX(1,6),N,
     $    NP,NBin,NZC,ZThr,.true.,.true.,.true.,.false.,InitA,InitA,
     $    InitA,InitA,InitA,InitA,InitA,InitA,1,3,Method)
        InitA = .false.
C$OMP Parallel Do Default(Shared) Schedule(Static,1)
C$OMP+ Private(I)
        Do 10 I = 1, NpNZ
   10     TT(I,ICnt2) = Two*PRODHT(I) - TT(I,ICnt3)
        Call AMove(LenRT,TT(1,ICnt2),RT(1,ICnter+1))
        If(StoreT) Call FileIO(1,IRWRT,NpNZ,TT(1,ICnt2),0)
        ICnt1 = Mod(ICnt1,3) + 1
        ICnt2 = Mod(ICnt2,3) + 1
   20   ICnt3 = Mod(ICnt3,3) + 1
      Return
      End
*Deck ChbFnc
      Function ChbFnc(C,FM,RT,X,XX,A,B,N,NOcc,NumC,np,LenRT,ifunc)
      Implicit Real*8(A-H,O-Z)
C
C     Evaluate the Density Matrix given A, B and return the
C     Tr(DM)-NOcc in ChbFnc
C
      Dimension C(NumC),RT(LenRT,NumC),FM(N,2),X(NumC,*),XX(N,NP)
      Save Two
      Data Two/2.d0/
C
C     Get the coefficients for the Chebyshev approximation of the Heavyside
C
      Call chbcof(X(1,1),X(1,2),C,A,B,NumC,ifunc)
C
C     Zero out some matrices
C
      Call AClear(N*NP,XX)
      Call AClear(N,FM)
C
C     Calculate the Fermi matrix
C
      C(1) = C(1)/Two
      Call MatMP3(1,1,1,N,LenRt,NumC,NumC,1,RT,C,FM)
      ChbFnc = ArrSum(N,FM) - GFloat(NOcc)
      Return
      End
*Deck ChbFnf
      Function ChbFnf(C,FM,RT,X,A,B,N,NOcc,NZ,NumC,LenRT,StoreT,IRWRT,
     $  TT,ifunc)
      Implicit Real*8(A-H,O-Z)
C
C     Evaluate the Density Matrix given A, B and return the
C     Tr(DM)-NOcc in ChbFnc
C
      Logical StoreT
      Dimension C(NUMC),RT(LenRT,NumC),FM(N,2),X(NumC,4),TT(N,2)
      Save Two
      Data Two/2.d0/
C
C     Get the coefficients for the Chebyshev approximation of the Heavyside
C
      Call chbcof(X(1,1),X(1,2),C,A,B,NumC,ifunc)
C
C     Zero out some matrices
C
      NpNZ = N + NZ
      Call AClear(NpNZ,FM)
C
C     Calculate the Fermi matrix
C
      If(StoreT) Call FileIO(2,-IRWRT,NpNZ,TT,0)
      Do 10 I = 2, NumC
        If(StoreT) then
          Call FileIO(2,IRWRT,NpNZ,TT,0)
          Call Acasb(NpNZ,FM,TT,FM,C(I))
        else
          Call Acasb(NpNZ,FM,RT(1,I),FM,C(I))
          endIf
   10   Continue
      CC = C(1)/Two
      Do 20 II = 1, N
   20   FM(II,1)=FM(II,1)+CC
C
C     Return the trace of FM
C
      ChbFnf = ArrSum(N,FM) - GFloat(NOcc)
      Return
      End
*Deck CmpXYZ
      Subroutine CmpXYZ(NStart,NDimX,Keep,CurCnt,AngWgt,NAct,DistMx,
     $  DistMn,XYZS,Wgt)
      Implicit Real*8(A-H,O-Z)
C
C     Compress an array of coordinates depending on a logical array.
C
      Logical Keep(NStart)
      Dimension CurCnt(3), AngWgt(*), XYZS(NDimX,3), Wgt(*)
C
      NAct = 0
      DistMx = GFloat(0)
      DistMn = GFloat(10000)
      Do 10 K = 1, NStart
        If(Keep(K)) then
          NAct = NAct + 1
          Dist2 = (XYZS(K,1)-CurCnt(1))**2+(XYZS(K,2)-CurCnt(2))**2+
     $      (XYZS(K,3)-CurCnt(3))**2
          Wgt(NAct) = AngWgt(K)*Dist2
          XYZS(NAct,1) = XYZS(K,1)
          XYZS(NAct,2) = XYZS(K,2)
          XYZS(NAct,3) = XYZS(K,3)
          DistMx = Max(DistMx,Dist2)
          DistMn = Min(DistMn,Dist2)
          endIf
   10   Continue
      DistMx = Sqrt(DistMx)
      DistMn = Sqrt(DistMn)
      Return
      End
*Deck CombRM
      Function CombRM(IR,M,J1,J2,RLSq,RLTmp,L)
      Implicit Real*8(A-H,O-Z)
C
C     CombRM of products of each combinations of IR elements of RLSq,
C     excluding J1 and J2.
C
      Dimension RLSq(*), RLTmp(*), L(*)
      Save Zero, One
      Data Zero/0.0d0/, One/1.0d0/
C
      If(IR.eq.0) then
        CombRM = One
        Return
        endIf
      CombRM = Zero
      I1 = Min(J1,J2)
      I2 = Max(J1,J2)
      If(I1.gt.M) I1 = 0
      If(I2.gt.M) I2 = 0
      If(I2.eq.0) then
        MX = M
        Call AMove(M,RLSq,RLTmp)
      else if(I1.eq.0) then
        MX = M - 1
        Do 10 I = 1, (I2-1)
   10     RLTmp(I) = RLSq(I)
        Do 20 I = (I2+1), M
   20     RLTmp(I-1) = RLSq(I)
      else
        MX = M - 2
        If(IR.gt.(M-2)) Return
        Do 30 I = 1, (I1-1)
   30     RLTmp(I) = RLSq(I)
        Do 40 I = (I1+1), (I2-1)
   40     RLTmp(I-1) = RLSq(I)
        Do 50 I = (I2+1), M
   50     RLTmp(I-2) = RLSq(I)
        endIf
      If(IR.gt.MX) Return
      P = One
      Do 60 I = 1, IR
        L(I) = IR - I + 1
   60   P = P*RLTmp(L(I))
      L(IR+1) = 0
      CombRM = P
      NCom = 1
C     Write(6,1000) NCom, P, CombRM, (l(kk),kk=1,ir)
   70 Do 80 I = IR, 2, -1
        L(I) = L(I) + 1
        If(L(I).lt.L(I-1)) goto 90
   80   L(I) = L(I+1) + 1
      L(1) = L(1) + 1
      If(L(1).gt.MX) Return
   90 NCom = NCom + 1
      P = One
      Do 100 I = 1, IR
  100   P = P*RLTmp(L(I))
      CombRM = CombRM + P
C     Write(6,1000) NCom, P, combrm, (l(kk),kk=1,ir)
      Goto 70
      End
*Deck ConOpn
      Subroutine ConOpn(IOut,IOpCl,N,NU,Acurcy,IFlag,IExtp,A1,A2,A3,
     $  Crit,CritM,CritE,DelEn,DPMax,ICount,Done,IRwPA,IRwPB,IRwC1,
     $  IRwC2,IRwC3,SmalEn,IPass,NMtPBC,Damp,DampFc)
      Implicit Real*8(A-H,O-Z)
C
C     SCF convergence testing.
C
      Dimension A1(*), A2(*), A3(*), XX(1)
      Character*30 ExtMsg(4)
      Logical Done, Damp, Do34, CBfn, Spinor
      Real*8 MDCutO
      Save Zero, Pt99, Pt995, One, OnePt9, Two, Four, SP12, SP22,
     $  ExtMsg, XX
      Data Zero/0.0d0/, Pt99/0.99d0/, Pt995/0.995d0/, One/1.0d0/,
     $  OnePt9/1.9d0/, Two/2.0d0/, Four/4.0d0/, SP12/0.0d0/,
     $  SP22/0.0d0/, ExtMsg/'3-Point extrapolation',
     $  '4-Point extrapolation','3-Point extrapolation skipped',
     $  '4-Point extrapolation skipped'/, XX/0.0d0/
 2000 Format(' Damping current iteration by',1PD9.2)
 2010 Format(1X,A,'.')
C
      Small = MDCutO(1)
      Call DecoSC(IOpCl,NSpBlk,NRI,NDimBl,CBfn,NSpBlX,Spinor)
      IOpClX = IOpCl - Mod(IOpCl,2)
      NMat = NSpBlk*NMtPBC
      NB = N*NDimBl
      NTT = (NB*(NB+1))/2
      IType = 0
      IFlag = 0
      Done = .False.
      ROp = GFloat(3-NSpBlk)
      NTTX = NMtPBC*NTT
      NTT2 = NSpBlk*NTTX
      NTTXR = NTTX*NRI
      NTT2R = NTT2*NRI
      ICount = ICount + 1
      Loc1 = IRwC3
      If(Mod(ICount,2).eq.1) Loc1 = IRwC2
      Loc2 = IRwC2 + IRwC3 - Loc1
      Call FileIO(2,-IRwPA,NTTXR,A3,0)
      If(NSpBlk.eq.2) Call FileIO(2,-IRwPB,NTTXR,A3(1+NTTXR),0)
      If(ICount.gt.1) Call FileIO(2,-IRwC1,NTT2R,A1,0)
      Call FileIO(1,-IRwC1,NTT2R,A3,0)
C     Nothing else to do in first cycle or immediately after extrapolation
      If(ICount.eq.1) Return
C     P(N) now in A3 and P(N-1) in A1; form P(N)-P(N-1) in A2.
      Call ASub(NTT2R,A3,A1,A2)
      DPMax = ArMaxC(NRI,A2,NTT2)
C     Find length DP1
      SP11 = ROp*SCFTrX(.False.,.False.,IOpClX,NMat,N,A2,XX,A2,XX)
      DP1 = Sqrt(SP11/Two)
C     Test for convergence by finding if root mean square dp is less
C     than Acurcy
      RMSDP = DP1 / (GFloat(NU)*Sqrt(GFloat(NDimBl*NMtPBC)))
      Acurcy = RMSDP
C     Convergence criterion met ... exit
C
C ANT BEGIN modification
C
C         Only converge on density matrix criterions
C
C     If(RMSDP.lt.Crit.and.DPMax.lt.CritM.and.
C    $  (DelEn.le.SmalEn.or.IPass.eq.0).and.
C    $  (Abs(DelEn).lt.CritE.or.CritE.eq.Zero)) then
      If(RMSDP.lt.Crit) then
C
C ANT END modification
C
        Done = .True.
        Return
        endIf
      If(Damp) then
        If(IPSave(0).eq.0) Write(IOut,2000) DampFc
        Call ACASB(NTT2R,A3,A2,A3,DampFc-One)
        Call FileIO(1,-IRwC1,NTT2R,A3,0)
        Call FileIO(1,-IRwPA,NTTXR,A3,0)
        If(NSpBlk.eq.2) Call FileIO(1,-IRwPB,NTTXR,A3(1+NTTXR),0)
        Return
        endIf
      Call FileIO(1,-Loc1,NTT2R,A2,0)
      If(ICount.le.3) Return
      Call FileIO(2,-Loc1,NTT2R,A1,0)
      SP23 = SP12
      SP33 = SP22
      SP13 = ROp*SCFTrX(.False.,.False.,IOpClX,NMat,N,A1,XX,A2,XX)
C     Find length DP3
      DP3 = Sqrt(SP33/Two)
C     Read P(N-1)-P(N-2) into A1
      Call FileIO(2,-Loc2,NTT2R,A1,0)
      SP12 = ROp*SCFTrX(.False.,.False.,IOpClX,NMat,N,A1,XX,A2,XX)
      SP22 = ROp*SCFTrX(.False.,.False.,IOpClX,NMat,N,A1,XX,A1,XX)
C     Find length DP2
      DP2 = Sqrt(SP22/Two)
C     Find cosine of angle between successive displacements
      Denom = Two*DP1*DP2
      CosPhi = Zero
      Do34 = Denom.gt.Small
      If(Do34) CosPhi = SP12 / Denom
C     Find cosine of angle between DP3 and plane of DP1 and DP2
      Denom = SP11*SP22 - SP12*SP12
      Do34 = Do34.and.DP3.gt.Small.and.Denom.gt.Small
C     Do not extrapolate unless 4 consecutive points are nearly coplanar
      If(Do34) then
        X = (SP13*SP22-SP12*SP23) / Denom
        Y = (SP23*SP11-SP12*SP13) / Denom
        CosPsi = Sqrt((X*X*SP11+Y*Y*SP22+Two*X*Y*SP12)/Two)/DP3
        ArSqrt = (X*X*SP11+Y*Y*SP22+Two*X*Y*SP12)/Two
        Do34 = ArSqrt.gt.Zero
        If (Do34) then
          CosPsi = Sqrt(ArSqrt)/DP3
          Do34 = Abs(X).gt.Small.and.CosPsi.gt.Pt99
        endif
      else
        CosPsi = Zero
        endIf
      If(Do34) then
C       Express vector DP1 as X*DP3(Projected)+Y*DP2
        Y = -Y/X
        X = One/X
C       Test if 2*2 matrix has real eigenvalues between -.95 and +.95
        XY = Y*Y + Four*X
        XY1 = Abs(Y) + Sqrt(Abs(XY))
        If(XY.ge.Zero.and.XY1.le.OnePt9) then
          If(IExtp.eq.2) then
            IType = 4
          else
            XXX = X / (One-X-Y)
            YYY = (X+Y) / (One-X-Y)
            Call ACASB2(NTT2R,A3,A1,A2,A3,XXX,YYY)
            IType = 2
            IFlag = 1
            ICount = 0
            endIf
        else if(Abs(CosPhi).gt.Pt995) then
          If(IExtp.ge.1) then
            IType = 3
          else
            X = DP1 / (DP2*CosPhi-DP1)
            Call ACASB(NTT2R,A3,A2,A3,X)
            IType = 1
            IFlag = 1
            ICount = 0
            endIf
          endIf
        endIf
      If(IFlag.eq.1) then
        Call FileIO(1,-IRwC1,NTT2R,A3,0)
        Call FileIO(1,-IRwPA,NTTXR,A3,0)
        If(NSpBlk.eq.2) then
          Call FileIO(1,-IRwPB,NTTXR,A3(1+NTTXR),0)
        else
          Call FileIO(1,-IRwPB,NTTXR,A3,0)
          endIf
        endIf
      If(IPSave(0).eq.0.and.IType.ge.1.and.IType.le.2)
     $  Write(IOut,2010) ExtMsg(IType)(1:LinEnd(ExtMsg(IType)))
      Return
      End
*Deck CycOpn
      Subroutine CycOpn(IOut,Energy,ETwoO,Acurcy,MaxCyc,ROHF,JCycle,
     $  IExtp,NAtoms,NAE,NBE,IAn,IAtTyp,AtmChg,C,IAtFrg,NSymOp,NOp1,
     $  NOp2,NEqBas,NBasis,NBD,NBsUse,NBas6D,NTT,IJ,D,F,FF,V,DA,DB,FA,
     $  FB,VA,VB,VV,Buf2E,ISort,IRSort,IBfSrt,MDV,LScr2,Scr2,IDIIS,
     $  IFBrod,BMat,IPrint,IPrPFS,ZDO,DirSCF,IncFoc,NEqShl,NEqSh2,RotOp,
     $  NEqAtm,JSym2E,JSym2M,IJMap,MemInt,IFinIt,ICntrl,EOneO,X,LenX,
     $  RHF,IRaf,IPFlag,AllowP,InStor,R1,R2,R3,SaveMO,NFullF,VarAcc,
     $  AccFac,Sleazy,EngCon,AccMax,CP,IVShft,VSGoal,VSMax,VSGpMn,ISCRF,
     $  ItSCRF,gFac,fuFac,EPol,EField,IExCor,IRadAn,IRanWt,IRanGd,
     $  ICorTp,ScaHFX,ScaDFX,IDiag,NPassI,IRwEig,IRwCA,IRwCB,IRwPA,
     $  IRwPB,IRwFA,IRwFB,IRwS,IRwH,IRwTM,IRwC1,IRwC2,IRwC3,ISyOcc,
     $  IWfnSy,NBlks,NBF,NBFU,NOccSy,ICtDFT,FMM,RComp,FMFlag,FMFlg1,
     $  NFxFlg,JustJ,IOCMOS,AbDiag,IFrcSy,ICharg,Multip,JustEn,EPert,
     $  SpFock,SpDiag,NBTI,ZThr,MxSCFC,TempSt,NSAStep,Epsi,RdPCM,IHMeth,
     $  IPBC,LIPBC,CelVec,ENR,VShftA,VShftB,LvFirst,IStab,IRwNCA,IRwNCB,
     $  IndFok,EnCyc,ErrCyc,FixRis,LvWFON,LvRise,AccDsI,MaxTyp,DLaMax,
     $  LLim,IULim,IOpDME,SECut,IRwPOA,IRwV,DThr,IProgC,DCut,IFxFM,
     $  IMemOp,IOpDM2,IRwPOB,SymmP,IRwFOA,IRwFOB,NKPnt,NBFK,ICoor,WKPnt,
     $  CmlxF,IRwCAC,IRwCBC,IRwEiC,NAEBr,NOABr,NBEBr,NOBBr,MxDmpI,IDamp,
     $  DampCut,DampLim,IDenFit,IODFit,IOFDen,KAlg,FkGues,NExDmI,IOGesA,
     $  IOSymA,IRwWgt,IOToK,NTot,IADMP,BndOnl,IGWInf,IAtBtD,IRdBtD,
     $  RRdBtD,RGWBtD,FrcNuc,Halt,IZMthd,FKTemp,kTLim,NSPnt,If0K,TigPCV,
     $  IRwPAN,IRwPBN,IRANP0,ThrSDi,Omega,LChg,NumChg,Chg,IDoGrd,LFXYZ,
     $  FXYZ,IOpCl,IRwSO,DoUpd,IRwUR,TstSpI,IStuFl,DoNDXC,ErfPO,DoDMS,
     $  IADBas,LSEAll,RSEAll,Y,LenY,NSAOBk,NS,NS2,NS4,ISABF,RSABF,
     $  MulTab,MapSAO,FrcSDg,DISPCM,IAprXC,ENJ,NDBF6D,ICampK,CKThr,
     $  NOnlyS,NFock,DoHirS,ISwCut,IPsCut,IFrgTp,MapF2M,IndFrg,EIter1,
     $  IAlDNI,AlGsMn)
      Implicit Real*8(A-H,O-Z)
C
C     Iteration routine:
C     Energy ... returned with total energy.
C     Acurcy ... Convergence required on density matrix.
C                Returned with actual convergence achieved.
C     MaxCyc ... Maximum number of iterations allowed
C
C     Note possible equivalences in calling routine:
C     F=D=(DA,FA=DB) (RComp)
C     D=(DA,DB) (F,FF)=(FA,FB) (.not.RComp)
C     V=(VA,VB) (and =Buf2E unless MemInt)
C     Y=BMat
C
C     For SpDiag:
C     IPStrt = pointer to beginning of density matrices
C     IFStrt = pointer to beginning of Fock/Hamiltonian matrices
C     IAftrF = pointer to free memory after Fock matrices
C
      Real*8 MDCutO
C
C ANT BEGIN modification
C
      Parameter (IRwGen=501,LenGen=55)
      Dimension DGen(LenGen)
C
C     Status variable for switching on ANT
C
      Logical ANTOn
      SAVE ANTOn
C
C ANT END modification
C
      Parameter (NExtra=5,MaxNU=2,NErSav=2,IHFJ=6)
      Character*1 JGC
      Integer VarAcc, FMFlag, FMFlg1, FixRis
      Logical ISort, DIIS, IRSort, ROHF, ZDO, DirSCF, MemInt, IncFoc,
     $  DoInc, RHF, SaveMO, Sleazy, P1MaxC, DFT, P1Clos, OKAltG, LJunk,
     $  P1Fail, P1Chek, Done, DoDIIS, IfDFT, MixedK, AllowP(*), DoDMS,
     $  Mixed, Exit, DFTSym, FMM, AbOnly, RComp, HaveB, DoPrnt, Confus,
     $  JustJ, AbDiag, ForcSy, DoSurf, StPass, DoSurO, DoRest, HavSrf,
     $  PCM, SCIPCM, TestB, SpFock, SpDiag, Sp1Chk, DoDiag, RdPCM, Damp,
     $  Rises, OldRis, JustEn, DoProg, SlezOK, SymmP, UHF, If0K,
     $  UsKPnt, DoPurP, FkGues, DoFkG, DidAlG, BndOnl, NewExF, DoNDXC,
     $  AboveM, OKExit, Halt, FrcNuc, FirstC, LastC, Var2P, Onsagr,
     $  SolEPt, CBfn, GHF, HavBRI, Spinor, DoUpd, TstSpI, FrcSDg,
     $  DoROHF, DoneQ, DISPCM, DidUpN, DidPCM, SkipEn, IniPCM, DoHirS,
     $  DgAllA, DgAllB, IOTMEx
      Dimension D(NBasis,NBasis), F(NBasis,NBasis), FF(*), DA(*),
     $  BMat(*), V(NBasis,NBasis), FA(*), FB(*), Scr2(NBD,LScr2),
     $  Buf2E(*), IJ(*), NEqBas(NBasis,*), IBfSrt(*), CP(MaxCyc+1,3),
     $  DB(*), NEqShl(*), VA(*), VB(*), IJMap(NTT,*), X(*), R1(*),
     $  R2(*), R3(*), AtmChg(*), C(3,*), EField(0:34), IAn(*), XX(1),
     $  VV(*), gFac(0:24), ScaDFX(*), NBFU(*), NOccSy(*), IAtFrg(*),
     $  RotOp(*), EnCyc(*), IAtTyp(*), IndFok(*), ErrMax(0:NErSav),
     $  ErrCyc(*), CelVec(3,3), LLim(*), IULim(*), IProgC(*), NBFK(*),
     $  CmlxF(*), WKpnt(*), NEqSh2(*), ECDFit(4), IPBC(*), ICoor(*),
     $  IOToK(*), NTot(*), IGWInf(*), IAtBtD(*), IRdBtD(*), RRdBtD(*),
     $  RGWBtD(*), NEqAtm(*), OldCof(2), Omega(*), Chg(LChg,*), Y(*),
     $  RNEFit(4), FXYZ(LFXYZ,2), RSEAll(*), NS(*), NS2(*), NS4(*),
     $  ISABF(*), RSABF(*), MulTab(*), MapSAO(*), EnJunk(1), EOneA(1),
     $  ETwoA(1), ExA(1), EcA(1), RNEA(1), ErfPA(1), EnerCK(0:1),
     $  MapF2M(*), IndFrg(*), NBF(*)
      Equivalence (EOne,EOneA(1)), (ETwo,ETwoA(1)), (Ex,ExA(1)),
     $  (Ec,EcA(1)), (RNE,RNEA(1)), (ErfP,ErfPA(1))
      Save Zero,PT5,One,Two,Pt001,Pt01,Three,Ten,F100,SymCut,Big,XX,
     $  GapMin,GapMax,SmalEn,Pt25,GapDamp,DmpCf1,DmpCf2,RMedEn,OneP1,
     $  PDEMax,PDRMax
      Data Zero,PT5,Two/0.D0,0.5D0,2.D0/,Pt001/0.001d0/,Pt01/0.01d0/,
     $  Three/3.0d0/,Ten/10.0d0/,F100/100.0d0/,SymCut/1.0d-7/,
     $  SyCDFT/1.d-4/,Big/1.0d10/,XX/0.0d0/,One/1.0d0/,PDEMax/2.d-2/,
     $  GapMin/0.03d0/,GapMax/0.15d0/,SmalEn/1.d-9/,Pt25/0.25d0/,
     $  GapDamp/0.2d0/,DmpCf1/0.1d0/,DmpCf2/0.3d0/,RMedEn/1.d-6/,
     $  OneP1/1.1d0/,PDRMax/3.0d-2/
 1000 Format(' Possible problem detected with PBC diagonalization:',/,
     $  ' FLenRc=',1PD9.2,' SDif=',1PD9.2,' STest=',1PD9.2,' ThrSDi=',
     $  1PD9.2,'.')
 1060 Format(' One-electron energy=',D19.12,' <L.S>=',D19.12)
 1070 Format(' One-electron energy=',D19.12)
 1110 Format(' Cycle',I4,'  Pass',I2,'  IDiag',I3,':')
 1120 Format(' E=',G22.15,' Delta-E=',F22.12,' Rises=',L1,' Damp=',L1)
 1130 Format(' E=',G22.15)
 1140 Format(' RMSDP=',1PD8.2,' MaxDP=',1PD8.2,'    ',9X,' OVMax=',
     $  1PD9.2)
 1145 Format(' RMSDP=',1PD8.2,' MaxDP=',1PD8.2,' DE=',1PD9.2,' OVMax=',
     $  1PD9.2)
 1150 Format(1X,10(1H>),' Convergence criterion not met.')
 1160 Format(1X,10(1H>),' Final iteration diverged -- try SCF=NoFinal.')
 1170 Format(' E=','   (Non-Variational)  ')
 1200 Format(' Convergence on energy, delta-E=',1PD9.2)
 1210 Format(' Initial convergence to',1PD8.1,
     $  ' achieved.  Increase integral accuracy.')
 1220 Format(' Failed to achieve initial convergence to',1PD8.1,
     $' in',I3,' cycles.',/,' Switch to full accuracy anyway.')
 1230 Format(' Problem detected with inexpensive integrals.',/,
     $  ' Switching to full accuracy and repeating last cycle.')
 1260 Format(' 1/2 <P',A,'(P)>=',F25.12)
 1270 Format(' Ex=',F25.12,' Ec=',F25.12,/,
     $  ' 1/2 <P',A,'(P)> + E(ex) + E(corr)=',F25.12)
 1275 Format(' Erf(P)=',F25.12)
 1280 Format(' Diagonalized old Fock matrix for initial guess.')
 1290 Format(' GapD=',F9.3,' DampG=',F5.3,' DampE=',F5.3,' DampFc=',
     $  F6.4,' IDamp=',I2,'.')
 1300 Format(' Recover alternate guess density for next cycle.')
 1310 Format(' Switch densities from cycles 1 and 2 for lowest energy.')
 1320 Format(' TrPSA=',F14.4,' TrPSB=',F14.4)
 1330 Format(' Exit after recomputing band structure.')
 1340 Format(' <P*S>=',F20.12,' NE=',F20.12,' Err=',1PD9.2,' Lim=',
     $  1PD9.2,/,' Possibly not enough cells included in PBC.')
 1370 Format(' Density matrix is not changing but DIIS error=',1PD9.2,
     $  ' CofLast=',1PD9.2,'.')
 1380 Format(' <Sx>=',F7.4,' <Sy>=',F7.4,' <Sz>=',F7.4,' <S**2>=',F7.4,
     $  ' S=',F7.4)
C
      If(SpFock.and.NBasis.ne.NBas6D)
     $  Call GauErr('No pure functions with sparsity yet.')
      Call DecoSC(IOpCl,NSpBlk,NRI,NDimBl,CBfn,NSpBlX,Spinor)
      DoROHF = ROHF.and.NDimBl.eq.1.and..not.Spinor.and.NAE.ne.NBE
      RNE = GFloat(NAE+NBE)
      DgAllA = .False.
      DgAllB = .False.
      RNEFit(1) = RNE
      NExDmp = NExDmI
      MaxDamp = MxDmpI
      UHF = .not.RHF
      NTT6D = (NBas6D*NDimBl*(NBas6D*NDimBl+1))/2
      IniPCM = .True.
      Call AClear(3*(MaxCyc+1),CP)
      NMat = 1
      NMatS = 1
      NMatT = 0
      NPDir = 0
      NMtPBC = 0
      IFStrt = 1
      IAftrF = 1
      IConvr = 1
      NOpOld = 0
      Call ObtCel(2,IPBC,1,NPDir)
      Call ObtCel(2,IPBC,2,NMtPBC)
      FLenRc = GnMxRc(NPDir,CelVec)
      NTT2 = (NBD*(NBD+1))/2
      NTTX = NRI*NTT2*NMtPBC
      NTTXB = NTTX*NSpBlk
      NOrbs = NBsUse*NDimBl
      NBNO = NRI*NBD*NOrbs
      UsKPnt = NMtPBC.gt.1.and.NKPnt.gt.0
      IPrt1 = IPrint
      If(IPSave(0).eq.0) IPrt1 = Max(IPrt1,1)
      DIIS = Mod(IDIIS,10).eq.1
      IDITyp = Mod(IDIIS,10000) / 10
      IPtDIS = Max(Mod(IDIIS,100000)/10000,IPrt1)
      AccDes = Zero
      AccXC = AccDes
      VShftA = Zero
      VShftB = Zero
      NSpChk = NSpBlk
      NAEX = NAE + NBE*(NDimBl-1)
      IXScr = 1
      IPStrt = IXScr + NBasis
      DoSurO = .False.
      HavSrf = .False.
      Rises = .False.
      AboveM = .False.
      OldCof(1) = Zero
      OldCof(2) = One
      Damp = .False.
      DidAlG = .False.
      GapD = Ten*GapDamp
      GapDB = GapD
      NSvDIS = 0
      CofLast = One
      IRwSrf = 0
      IRwSCH = 0
      DoPrnt = IPrt1.gt.0
      GHF = Spinor.or.NDimBl.eq.2
      HavBRI = (.not.RComp.and..not.SpDiag).or.(SpDiag.and..not.RHF)
      HaveB = HavBRI.and..not.GHF
      HavBRI = HavBRI.or.(NDimBl.eq.2.and..not.CBfn)
      TestB = .not.RHF
      ISCRF0 = Mod(ISCRF,10000)/1000
      ISCRFT = Mod(ISCRF,1000)/100
      If(ISCRF0.eq.0) ISCRFT = ISCRFT*100
      SCIPCM = ISCRF0.eq.1
      PCM = ISCRF0.eq.2
      DFT = IfDFT(IExCor)
      If(FrcSDg) then
        CutOff = GFloat(1000000)
      else if(DFT.or.SCIPCM.or.PCM) then
        CutOff = SyCDFT
      else
        CutOff = SymCut
        endIf
      IFSelI = 0
      Temp = TempSt
      NIter = -1
      LastDg = 1
      IOBDIS = 0
      IOPADI = 0
      IOPBDI = 0
      IOFADI = 0
      IOFBDI = 0
      IOErDI = 0
      Beta = Zero
      If(AbDiag) then
        IOCAS = IOCMOS
        IOCBS = IOCMOS
        If(ROHF) then
          IndCBS = 0
        else
          IndCBS = NBNO
          endIf
      else
        IOCAS = IRwCA
        IOCBS = IRwCB
        IndCBS = 0
        endIf
      If(AbDiag) then
        IOTMA = IOCAS
        IOTMB = IOCBS
      else
        IOTMA = IRwCA
        IOTMB = IRwCB
        endIf
      IndTMB = IndCBS
      If(ROHF) then
        IOCB = IRwCA
      else
        IOCB = IRwCB
        endIf
      IOSBlk = 0
      IOSInv = 0
      IndSBB = 0
      IOTMEx = IFrgTp.le.1
      IRwS1 = IRwS
      If(ZDO) IRwS1 = 0
      IOpClX = IOpCl - Mod(IOpCl,2)
      IOpClY = IOpClX
      If(HaveB) IOpClY = IOpClY + 1
      IOpClD = IOpCl
      If(ROHF.and..not.RHF) IOpClD = 0
      NDimCP = NFullF + 1
      If(IncFoc) then
        IOPSav = IGetRW(NSpBlk*NFullF*NTTX,0)
        IOFSav = IGetRW(NSpBlk*NFullF*NTTX,0)
      else
        IOPSav = 0
        IOFSav = 0
        endIf
      If(IPrint.ne.0) then
        IPrtF = IPrint - 1
      else if(DoPrnt) then
        IPrtF = 0
      else
        IPrtF = -1
        endIf
      Var2P = VarAcc.eq.3.or.VarAcc.eq.4.or.VarAcc.eq.7.or.VarAcc.eq.8
      If(Var2P.and..not.SpDiag) then
        Call PCOld(-1,NTTX,NDimBl,NSpBlk,NBNO,IRwPA,IRwPB,IRwCA,IRwCB,
     $    IOCMOS,IRwOld,X,LenX)
      else if(SpDiag) then
        IRwOld = IGetRW(10,0)
      else
        IRwOld = 0
        endIf
      NBSq = NBasis*NBasis
      Crit = Acurcy
      CritM = AccMax
      EnOld = GFloat(10000)
      EnMin = EnOld
      Energy = EnOld
      EIter1 = EnOld
      EPert = Zero
      ECDFit(1) = Zero
      Call AClear(NErSav+1,ErrMax)
      NOpUse = 0
      Call AClear(2*LFXYZ,FXYZ)
      ICount = 0
      IFlag = 0
      IFlagO = 0
      Acurcy = Two
      JCycle = 0
      JCyc1 = 0
      NSaved = 0
      BestDP = Big
      DPMax = Big
      OVMax = Zero
      OVMaxB = Zero
      ELS = Zero
      If(IDamp.eq.-3) then
        OVLim = DampCut
      else
        OVLim = -One
        endIf
      NumUp = 0
      OkAltG = ITqry(IOGesA).gt.0
C
C     IPass = 0 ... preliminary low-accuracy pass for direct, or
C                   diagonalization of read-in Fock matrix.
C     IPass = 1 ... main iterations.
C     IPass = 2 ... Final clean-up iterations.
C
      If(FkGues.or.(Var2P.and..not.SpDiag)) then
        IPass = 0
      else
        IPass = 1
        endIf
      If(Var2P.and..not.FkGues.and..not.SpDiag) then
        Crit1 = AccFac
        Crit1M = F100*Crit1
        Crit1E = Zero
      else
        Crit1 = Crit
        Crit1M = CritM
        Crit1E = EngCon
        endIf
      LimP0 = MaxCyc/3
      JLast = 0
      NCycP0 = 0
      DoDiag = .False.
      DoFkG = .False.
      P1Fail = .False.
C
C     Figure out which files with Fock matrices (normal, RO merged
C     Fock, or Gamma-point, etc.) to use in extrapolation.
C
      IRwF1 = IRwFA
      IRwP1 = IRwPA
      If(ROHF) then
        IRwF2 = 0
        IRwP2 = 0
      else
        IRwF2 = IRwFB
        IRwP2 = IRwPB
        endIf
      IRwF1X = IRwF1
      IRwF2X = IRwF2
      IRwP1X = IRwP1
      IRwP2X = IRwP2
      If(UsKPnt) then
        IRwF1X = IGetRW(NTT,0)
        IRwP1 = IGetRW(NTT,0)
        Call Trsfr(NTT,IRwPA,IRwP1,X,LenX)
        If(.not.RHF) then
          IRwF2X = IGetRw(NTT,0)
          IRwP2 = IGetRw(NTT,0)
          Call Trsfr(NTT,IRwPB,IRwP2,X,LenX)
          endIf
      else if(ROHF.and..not.RHF) then
        IRwF1 = IGetRW(NTTX,0)
        IRwF1X = IRwF1
        endIf
      IndCK = -1
      Len1CK = NSpBlk*(NBNO+NTTX)
      EnerCK(0) = EnOld
      EnerCK(1) = EnOld
      If(ICampK.eq.0) then
        IRwCK = 0
      else
        IRwCK = IGetRW(2*Len1CK,0)
        endIf
      DidUpN = .False.
      IOnlyS = NOnlyS - 1
      If(NOnlyS.gt.0) then
        ISavGF = IGetRw(NSpBlk*NTTX,0)
      else
        ISavGF = 0
        endIf
C
C     Initialize progressive convergence and some pointers for full matrix Fock
C     formation but sparse diagonalization.
C
      ZThrP = Zero
      DoProg = SpDiag.and.VarAcc.ge.3.and..not.FkGues
      If(SpDiag) Call ProgCv(IOut,IPrint,IProgC,Acurcy,DoProg,ZThr,
     $  AccDes,ZThrP,AccDsC,AccDsD,TKos3,TKos4,pfthr,.True.,jcycle,
     $  Crit1,SpDiag,DPMax,TigPCV)
      If(.not.SpFock) ZThrP = Zero
      IRwFM1 = 0
      IRwFM2 = 0
      If(SpDiag) Call CGInit(IOut,IPrint,IRwPOA,IRWPA,NOp1,NEqBas,
     $  NBasis,IRwS,IRwTM,IRwV,IRwH,ZThr,ZThrP,MxSCFC,NAEX,ZDO,IMemOp,
     $  UHF,LLim,IULim,NAtoms,X,LenX)
      If(SpDiag.and..not.RHF) Call CGInit(IOut,IPrint,IRwPOB,IRWPB,NOp1,
     $  NEqBas,NBasis,IRwS,IRwTM,IRwV,IRwH,ZThr,ZThrP,MxSCFC,NBE,ZDO,
     $  IMemOp,.True.,LLim,IULim,NAtoms,X,LenX)
C
C     Top of SCF iteration loop.  Test for convergence (Key=0) and 3 or 4
C     point extrapolation (IFlag=1).
C
      DoneQ = .not.PCM
   10 DelEn = Energy - EnOld
      If(JCycle.eq.1) EIter1 = Energy
      Confus = .False.
      If(JCycle.le.NExDmp) then
        FactI1 = Pt5
      else
        FactI1 = One
        endIf
      If(ErrMax(0).lt.(Ten*DampCut)) then
        DampE = One
      else if(ErrMax(0).le.(Ten*Ten*DampCut)) then
        DampE = Pt5*FactI1
      else
        DampE = Pt25*FactI1
        endIf
      If(GapD.le.(Pt25*GapDamp)) then
        DampG = Pt25
      else if(GapD.le.(Pt5*GapDamp)) then
        DampG = Pt5
      else if(GapD.le.(Two*GapDamp)) then
        DampG = One
      else
        DampG = Two
        endIf
C     Decide whether DIIS is giving us a new fock matrix or mostly stepping back.
      NewExF = ((OldCof(1).lt.DmpCf1.or.OldCof(2).gt.DmpCf2)
     $  .and..not.AboveM)
     $  .or.Abs(OldCof(1)).gt.One.or.Abs(OldCof(2)).gt.One.or.P1Fail
      Damp = JCycle.le.MaxDamp.and.IDamp.ge.-1.and.ErrMax(0).ge.DampCut
     $  .and.NewExF.and..not.DidAlG
      If(.not.Damp) then
        DampFc = One
      else if(IDamp.gt.0) then
        DampFc = GFloat(IDamp) / (Ten*Ten)
      else
        DampFc = Max(DampE*DampG,FactI1*DampLim)
        endIf
      If(DoPrnt.and.JCycle.gt.0.and.Damp)
     $  Write(IOut,1290) GapD, DampG, DampE, DampFc, IDamp
      DampFc = Min(DampFc,One)
      Damp = Damp.and.DampFC.ne.One
      If(JCycle.eq.1.or.IHMeth.eq.14) then
        Crt1EX = Zero
      else
        Crt1EX = Crit1E
        endIf
      SkipEn = FkGues.and.IPass.eq.0.and..not.DoFkG
      If(SkipEn) then
C       Nothing to do
      else if(SpDiag) then
        If(RHF) then
          Call SpGet(X,X,LenX,IRWPA,IPStrt,jip,ipa,IFStrt,NBasis,NZP,1,
     $      .true.,0)
          ipb = ipa
        else
          Call MkSmFm(X,LenX,IRwPA,0,0,IRwPB,0,0,IPStrt,IPStrt,
     $      0,jip,ipa,ipb,NZP,IFStrt,NBasis,2)
          endIf
        Call SpConv(X,LenX,X(IXScr),jip,ipa,ipb,IFStrt,NBasis,Acurcy,
     $    IOut,IOpCl,IFlag,IExtp,Crit1,Crit1M,Sleazy,DelEn,DPMax,ICount,
     $    Done,IRwC1,IRwC2,IRwC3)
      else if(IOnlyS.eq.(NOnlyS-1)) then
        Call ConOpn(IOut,IOpCl,NBasis,NBsUse,Acurcy,IFlag,IExtp,V,DA,FA,
     $    Crit1,Crit1M,Crt1EX,DelEn,DPMax,ICount,Done,IRwPA,IRwPB,IRwC1,
     $    IRwC2,IRwC3,SmalEn,IPass,NMtPBC,Damp,DampFc)
        Done = Done.or.IPass.eq.2
        endIf
      Done = Done.and..not.DidAlg
      If(DoFkG) then
        Energy = Zero
      else
        JCycle = JCycle + 1
        If(NOnlyS.gt.0) IOnlyS = Mod(IOnlyS+1,NOnlyS)
        endIf
      Call SetPN(JCycle)
      DoDiag = DoDiag.and.(NOnlyS.eq.0.or.IOnlyS.eq.0)
      If(DoPrnt) then
        If(JCycle.eq.2.and.DoDiag) then
          Write(IOut,1140) Acurcy, DPMax, OVMax
        else if(JCycle.gt.2.and.DoDiag) then
          Write(IOut,1145) Acurcy, DPMax, DelEn, OVMax
          endIf
        Write(IOut,'(1X)')
        endIf
      Call TStamp(1,'Start of iteration.')
      If(IFlag.eq.0) JCyc1 = JCyc1 + 1
      If(Done.and.DoDiag.and..not.Damp.and..not.DidUpN.and.NOnlyS.gt.0)
     $  then
        NOnlyS = 2*NOnlyS
        DidUpN = .True.
        endIf
      Done = (Done.and.DoDiag.and..not.Damp.and.DoneQ).or.DoFkG.or.
     $  (Sleazy.and.IPass.eq.2)
      OKExit = IPass.eq.2.or..not.DoUpd.or.
     $  (IPass.eq.1.and.((IFinIt.eq.2.and.ItSCRF.le.2).or.
     $  (IFinIt.eq.3.and.DgAllA.and.DgAllB)))
      If(Done.and.OKExit) then
        Exit = .True.
      else if(Sleazy.and.JCycle.gt.1.and.IFlag.eq.0.and.OKExit.and.
     $  Energy.lt.EnOld.and.(-DelEn).lt.EngCon) then
        Write(IOut,1200) DelEn
        Exit = .True.
      else if(BndOnl.and.FkGues.and.DoFkG) then
        Write(IOut,1330)
        Exit = .True.
      else
        Exit = .False.
        endIf
  100 If(Exit) then
        If(ErrMax(0).gt.(Crit*GFloat(100)).and..not.JustEn.and.
     $    .not.Sleazy.and.NSvDIS.gt.0.and.CofLast.lt.Pt5.and.
     $    (.not.PCM.or.(PCM.and..not.DoneQ))) then
          Write(IOut,1370) ErrMax(0), CofLast
          Done = .False.
          Exit = .False.
C
C ANT BEGIN modification
C
C         Done = .False.
C         Exit = .False.
C         IDITyp = 1
C         IFSelI = 1
C         IPass = 1
C         Goto 110
C
C ANT END modification
C
C         Call GauErr('The SCF is confused.')
          endIf
        IConvr = 2
        If(.not.JustEn) JCycle = JCycle - 1
        If(SpDiag) then
          Call SpCk2n(IOut,JSym2E,NOp1,JSym2X,NOpUse)
        else
          Call FileIO(2,-IRwPA,NTTX,DA,0)
          Call FileIO(2,-IRwPB,NTTX,DB,0)
          If(NSymOp.eq.1) then
            NOpUse = 1
            JSym2X = 0
          else
            Call RIComp(.True.,TestB,.False.,.True.,IOpCl,NBasis,1,1,1,
     $        0,NMat,NMatS,NMatT,DA,DB,X,LenX)
            Call PChk2N(IOut,IPrPCh,.True.,TestB,SymmP,JSym2E,NOp2,
     $        NSymOp,NBasis,NBas6D,NMat,NMatS+NMatT,CutOff,NEqShl,RotOp,
     $        DA,DB,JSym2X,NOpUse,X,LenX)
            Call RIComp(.False.,TestB,.False.,.True.,IOpCl,NBasis,1,1,1,
     $        0,NMat,NMatS,NMatT,DA,DB,X,LenX)
            endIf
          If(HaveB.and..not.TestB) Call AMove(NTTX,DA,DB)
          Call FileIO(1,-IRwPA,NTTX,DA,0)
          Call FileIO(1,-IRwPB,NTTX,DB,0)
          If(DirSCF.and.InStor.eq.3) then
            If(JSym2X.lt.JSym2M)
     $        Call GauErr('No symm in density but symm in integrals.')
            JSym2X = Min(JSym2X,JSym2M)
            endIf
          endIf
        Goto 999
        endIf
  110 If(Var2P.and.IPass.eq.0.and.DoDiag.and..not.Confus) then
        If(Acurcy.lt.BestDP.or.
     $    (Acurcy.le.(OneP1*BestDP).and.DelEn.lt.Zero)) then
          BestDP = Min(BestDP,Acurcy)
          NumUp = 0
        else if(JCycle.gt.1) then
          If(DelEn.ge.Zero) then
            NumUp = NumUp + 1
          else
            NumUp = 0
            endIf
          endIf
        NumUp = Max(NumUp,0)
        P1MaxC = JCycle.ge.LimP0
        P1Clos = IncFoc.and.(Acurcy.le.(Crit1*Three))
     $    .and.Mod(NSaved,NFullF).eq.0
        If(SpDiag) then
          P1Fail = Sp1Chk(X,LenX,IPStrt,jip,ipa,IFStrt,NBasis,NZP,
     $      JCycle,Energy,DelEn,Acurcy,IRwOld,IOpCl)
        else
          P1Fail = P1Chek(JCycle,Energy,DelEn,Acurcy,NTTX,NDimBl,NSpBlk,
     $      NBNO,IRwPA,IRwPB,IRwCA,IRwCB,IOCMOS,IRwOld,NumUp,MaxNU,X,
     $      LenX)
          endIf
C       If restarting with good integrals, reset damping.
        If(P1Fail) then
          MaxDamp = JCycle + MxDmpI
          ICount = 0
          Damp = .False.
          NumUp = 0
          Call ConOpn(IOut,IOpCl,NBasis,NBsUse,Acurcy,IFlag,IExtp,V,DA,
     $      FA,Crit1,Crit1M,Crt1EX,DelEn,DPMax,ICount,LJunk,IRwPA,IRwPB,
     $      IRwC1,IRwC2,IRwC3,SmalEn,IPass,NMtPBC,Damp,DampFc)
          endIf
      else
        P1MaxC = .False.
        P1Clos = .False.
        P1Fail = Confus.and.IPass.eq.0
        endIf
      If(IPass.eq.0.and.(Done.or.P1MaxC.or.P1Clos.or.P1Fail)) then
        NCycP0 = JCyc1 - 1
        Crit1 = Crit
        Crit1M = CritM
        Crit1E = EngCon
        NIter = -1
        If(DoFkG) then
          Write(IOut,1280)
        else if(Done.and.Var2P) then
          Write(IOut,1210) AccFac
        else if(P1Clos) then
          Write(IOut,1210) Acurcy
        else if(P1MaxC) then
          Write(IOut,1220) Crit1, LimP0
        else if(P1Fail) then
          Write(IOut,1230)
          endIf
      else if(IPass.eq.2.and..not.Done.and.JCycle.ge.(JLast+NExtra))then
        Write(IOut,1160)
        Goto 999
      else if(Confus) then
        NIter = -1
      else
        NIter = JCyc1 - NCycP0
        endIf
      NCycPD = IAbs(IDiag)
      IDiagX = IDiag
      If(IDiagX.lt.0) IDiagX = 3
      If(IDiagX.eq.3) then
        If(Mod(NIter-LastDg,NCycPD).eq.0) IDiagX = 1
        If(Damp.or.ErrMax(0).gt.PDEMax.or.Acurcy.gt.PDRMax) IDiagX = 1
        endIf
      If(Done.or.P1MaxC.or.P1Clos.or.P1Fail) then
        JLast = JCycle - 1
        IPass = IPass + 1
        If(IDiagX.eq.3) IDiagX = 1
        If(.not.P1Fail.and..not.DoFkG) ICount = 0
        If(IPass.ne.2) EnOld = GFloat(10000)
        endIf
      If(IPass.eq.2.and.IDiagX.eq.3) IDiagX = 1
      StPass = JCycle.eq.(JLast+1)
      DoDIIS = DIIS.and.((IPass.eq.0.and..not.FkGues).or.
     $  (IPass.eq.1.and.(VarAcc.eq.1.or.VarAcc.ge.3.or..not.DirSCF)).or.
     $  (IPass.eq.2.and.IFinIt.eq.2))
      DoSurf = ItSCRF.le.1.or.StPass.or.(ItSCRF.eq.2.and.IPass.ge.1)
     $  .or.(ItSCRF.ge.3.and.IPass.eq.2).or.
     $  (ItSCRF.ge.4.and.IPass.eq.1.and.Acurcy.lt.(GFloat(1000)*Crit))
      If(JCycle.gt.MaxCyc) then
        Write(IOut,1150)
        Goto 999
        endIf
      If(DoPrnt) Write(IOut,1110) JCycle, IPass, IDiagX
      If(IPass.eq.2.and.IFinIt.eq.3) goto 500
C
C     Recover core hamiltonian and density matrices.  If two-electron
C     integral symmetry is being handled by symmetrizing the Fock
C     matrices, check that the density matrices have the full (abelian)
C     molecular symmetry.
C
C
C ANT BEGIN modification
C
C         Call ANT before SCF begins for initialization (Cycle "0")
C
      IF(jcycle.eq.1) call ANT (UHF,0,IRwH,IRwPA,IRwPB,IRwFA,
     $                IRwFB,IRwS1,IRwEig,Acurcy,Crit,ANTOn,NBasis)
C
C ANT END modification
C
      If(SpDiag) then
        Call SpGet(X,X,LenX,IRWH,IFStrt,jih,iha,IAftrF,NBasis,NZH,
     $    1,.True.,0)
        ihb = iha
      else
        If(IStab.eq.2.and.JCycle.le.LvFirst) Call StbUpd(IOut,DoPrnt,
     $    IOpCl,JCycle,NBasis,NBsUse,NTT,NAEX,NBE,IRwCA,IRwCB,IRwNCA,
     $    IRwNCB,IRwPA,IRwPB,DA,FA,V)
        Call FileIO(2,-IRwPA,NTTX,DA,0)
        If(HaveB) Call FileIO(2,-IRwPB,NTTX,DB,0)
        Call FileIO(2,-IRwH,NTTX,FA,0)
        If(HaveB) Call FileIO(2,-IRwH,NTTX,FB,NTTX)
        endIf
      If(IPrint.gt.2) then
        If(SpDiag) then
          Call SpDisp(IOut,X(jip),X(ipa),NBasis,'Alpha Density Matrix:',
     $      2)
          Call SpDisp(IOut,X(jih),X(iha),NBasis,
     $      'Core Hamiltonian Matrix:',2)
        else
          Call PrtFok(IOut,IOpClY,'Core hamiltonian',0,1,NMtPBC,NBasis,
     $      FA,FB,X,LenX)
          Call LTOutX(IOut,0,IOpClY,NMtPBC,0,'Density matrix',NBasis,
     $      NBasis,DA,DB,0)
          endIf
        endIf
      IPrPCh = IPrint
      If(IPSave(0).ge.2.and.IPrPCh.eq.0) IPrPCh = -1
      ForcSy = NSymOp.gt.1.and.
     $  (IFrcSy.eq.2.or.(JCycle.eq.1.and.IFrcSy.eq.1))
      If(ForcSy) then
        CutPC = GFloat(1000000)
      else if(JCycle.eq.1) then
        CutPC = SyCDFT
      else
        CutPC = CutOff
        endIf
      If(SpDiag) then
        Call SpCk2n(IOut,JSym2E,NOp1,JSym2X,NOpUse)
        JSym2D = JSym2X
      else
        If(NSymOp.eq.1) then
          NOpUse = 1
          JSym2X = 0
          JSym2D = 0
        else
          Call RIComp(.True.,TestB,.False.,.True.,IOpCl,NBasis,1,1,1,0,
     $      NMat,NMatS,NMatT,DA,DB,X,LenX)
          Call PChk2N(IOut,IPrPCh,.True.,TestB,SymmP,JSym2E,NOp2,NSymOp,
     $      NBasis,NBas6D,NMat,NMatS+NMatT,CutPC,NEqShl,RotOp,DA,DB,
     $      JSym2X,NOpUse,X,LenX)
          Call RIComp(.False.,TestB,.False.,.True.,IOpCl,NBasis,1,1,1,0,
     $      NMat,NMatS,NMatT,DA,DB,X,LenX)
          endIf
        If(HaveB.and..not.TestB) Call AMove(NTTX,DA,DB)
        JSym2D = JSym2X
        If(DirSCF.and.InStor.eq.3) then
          If(JSym2X.lt.JSym2M)
     $      Call GauErr('No symm in density but symm in integrals.')
          JSym2X = Min(JSym2X,JSym2M)
          endIf
        If(ForcSy) then
          If(IRwS1.eq.0) then
            TrA = TraLT(NMtPBC,NBasis,DA)
            TrB = TraLT(NMtPBC,NBasis,DB)
          else
            Call FileIO(2,-IRwS,NTTX,X,0)
            TrA = SCFTrX(.False.,.False.,IOpClX,NMtPBC,NBasis,DA,XX,X,
     $        XX)
            TrB = SCFTrX(.False.,.False.,IOpClX,NMtPBC,NBasis,DB,XX,X,
     $        XX)
            endIf
C
C ANT BEGIN modification
C
C         Allows for personal control of NAE and NBE (ignores multiplicity)
C
      if (ANTOn.or.jcycle.eq.1) goto 666
       Call AScale(NTTX,GFloat(NAEX)/TrA,DA,DA)
         If(NBE.gt.0.and..not.RComp) 
     $       Call AScale(NTTX,GFloat(NBE)/TrB,DB,DB)
 666   continue
      endif
C
C ANT END modification
C
        Call FileIO(1,-IRwPA,NTTX,DA,0)
        Call FileIO(1,-IRwPB,NTTX,DB,0)
        endIf
      If(ABOnly(0).or.SpFock.or.SpDiag) NOpUse = Min(NOpUse,NOp1)
      DFTSym = JSym2D.eq.1.and.NOpUse.gt.1
      Call TStamp(1,'After PChk2N.')
C
C     Compute one-electron contribution to the energy.
      If(SkipEn) then
        EOne = Zero
      else if(SpDiag) then
        EOne = SpTrc2(X(jih),X(iha),X(jip),X(ipa),X(ixscr),NBasis)
        If(HaveB) EOne = EOne +
     $    SpTrc2(X(jih),X(ihb),X(jip),X(ipb),X(ixscr),NBasis)
      else
        EOne = SCFTrX(.False.,.False.,IOpClX,NMtPBC,NBasis,DA,XX,FA,XX)
        If(HaveB) then
          EOne = EOne + SCFTrX(.False.,.False.,IOpClX,NMtPBC,NBasis,DB,
     $      XX,FB,XX)
        else if(NDimBl.eq.2.and.ITqry(IRwSO).gt.0) then
          Call FileIO(2,-IRwSO,NTTX,FA,0)
          ELS = SCFTrX(.False.,.False.,IOpClX,NMtPBC,NBasis,DA,XX,FA,XX)
          endIf
        endIf
      If(.not.HaveB.and..not.GHF) EOne = EOne + EOne
      If(IPrint.gt.0.or.JustJ) then
        If(ELS.ne.Zero) then
          Write(IOut,1060) EOne, ELS
        else
          Write(IOut,1070) EOne
          endIf
        endIf
C
C     Form two-electron part of the Fock matrices.
      DoRest = StPass.or.(NIter.gt.1.and.NOpUse.gt.NOpOld)
      NOpOld = NOpUse
      If(DoRest) NSaved = 0
      DoInc = IncFoc.and.IPass.ne.2
      If(NIter.gt.1.and.IPass.ne.2.and.(((IncFoc.or.DIIS).and.DoRest)
     $  .or.(ItSCRF.ge.4.and.(DoSurf.and..not.DoSurO)))) then
        If(IPass.eq.0) then
          NIter = -1
          NCycP0 = JCyc1 - 1
          endIf
        endIf
      ICutO = Mod(ICntrl,100)/10
      If(Sleazy.and.ICutO.eq.0) ICntrl = ICntrl + 40
      If((VarAcc.eq.2.and.JCycle.eq.1).or.
     $   (VarAcc.ge.3.and.IPass.eq.0.and.IDenFit.le.1)) then
        If(VarAcc.eq.7) then
          AccDes = AccDsI
        else
          AccDes = AccFac/Ten
          endIf
        AccXC = Zero
      else if(VarAcc.eq.2) then
        AccDes = Max(Pt001*Acurcy,Pt01*Crit)
        AccXC = AccDes
      else if(VarAcc.eq.3) then
        AccDes = Pt01*Crit
        AccXC = AccDes
      else
        AccDes = AccDsI
        AccXC = AccDes
        endIf
      If(Sleazy.and.AccDes.eq.Zero) then
        AccDag = AccFac/Ten
      else
        AccDag = AccDes
        endIf
      DoProg = SpDiag.and.VarAcc.ge.3.and..not.DoFkG
      If(SpDiag) Call ProgCv(IOut,IPrint,IProgC,Acurcy,DoProg,ZThr,
     $  AccDes,ZThrP,AccDsC,AccDsD,TKos3,TKos4,pfthr,.False.,JCycle,
     $  Crit1,SpDiag,DPMax,TigPCV)
      Call TStamp(1,'Before 2e terms.')
      Mixed = MixedK(IHMeth,DFT,ScaHFX)
      If(DFT.and..not.Mixed) then
        JGC = 'J'
      else
        JGC = 'G'
        endIf
      DoFkG = FkGues.and.IPass.eq.0
      SlezOK = IPass.eq.0.and.VarAcc.ne.8
      IRAP = IDFTGd(SlezOK,IRadAn)
      If(IPass.eq.2.and.Sleazy.and.IRadAn.eq.2) IRAP = 3
      If(IPass.eq.0.and.IRANP0.ne.0) IRAP = IRANP0
      ICntr5 = Mod(ICtDFT,1000000)/100000
      If(ForcSy.or.Temp.ne.Zero.or.Damp.or.IADMP.eq.1.or.NRI.eq.2) then
        ICtDF1 = ICtDFT - ICntr5*100000 + 100000
      else
        ICtDF1 = ICtDFT
        endIf
      If(DoNDXC) ICtDF1 = ICtDF1 + 100000000
      DidPCM = .False.
      Ex = Zero
      Ec = Zero
      If(DoFkG) then
        ETwo = Zero
      else if(SpFock.or.SpDiag) then
        NFock = NFock + 1
        Call MkFkSp(IOut,IPrint,IPrtF,IOpCl,SpFock,SpDiag,ZDO,JustJ,
     $    IRwH,IRwCA,IOCB,IRwPA,IRwPB,NAtoms,NBasis,NBsUse,NBas6D,NAE,
     $    NBE,ICharg,Multip,IAn,IAtTyp,AtmChg,C,IAtFrg,NMtPBC,IPBC,
     $    LIPBC,CelVec,RHF,ROHF,HaveB,HavBRI,DirSCF,FrcNuc,IHMeth,
     $    IDenFit,IODFit,ICntrl,IRaf,KAlg,IHFJ,Omega,FMM,FMFlag,FMFlg1,
     $    NFxFlg,IPFlag,AllowP,LSEAll,RSEAll,AccDes,AccXC,DFT,ICtDF1,
     $    IExCor,IRAP,IRanWt,IRanGd,ICorTp,ScaHFX,ScaDFX,IGWInf,IAtBtD,
     $    IRdBtD,RRdBtD,RGWBtD,DFTSym,JSym2X,NOpUse,NOp1,NOp2,NEqAtm,
     $    NEqShl,NEqBas,RotOp,NEqSh2,NSAOBk,NS,NS2,NS4,ISABF,RSABF,
     $    MulTab,MapSAO,DCut,EOne,JGC,ETwo,Ex,Ec,ECDFit,RNEFit,IXScr,
     $    IPStrt,IFStrt,JIP,IPA,IPB,JIFA,JIFB,IFA,IFB,IAftrF,ZThrP,NZP,
     $    NZF,NBTI,DA,DB,FA,FB,InStor,R1,R2,R3,Buf2E,IJ,IJMap,X,LenX)
      else
        NMatI = 1
        NMatSI = 1
        NMatTI = 0
        If(IOnlyS.gt.0) then
          Call FileIO(2,-ISavGF,NTTX,FA,0)
          If(NSpBlk.eq.2) then
            Call FileIO(2,ISavGF,NTTX,FB,0)
          else
            Call AMove(NTTX,FA,FB)
            endIf
          Call FileIO(1,-IRwFA,NTTX,FA,0)
          Call FileIO(1,-IRwFB,NTTX,FB,0)
        else
          NFock = NFock + 1
          DidPCM = .not.SCIPCM
          Call MkFok1(IOut,IPrint,IOpCl,DoROHF,RComp,JustJ,IniPCM,
     $      JSym2X,NMatI,NMatSI,NMatTI,NMtPBC,NAE,NBE,ICharg,Multip,
     $      NBasis,NBas6D,NBsUse,DA,DB,FA,FB,NOpAll,NOpUse,NOp1,NOp2,
     $      DirSCF,ICntrl,IPFlag,AllowP,FMM,FMFlag,FMFlg1,NFxFlg,IRaf,
     $      KAlg,IHMeth,LSEAll,RSEAll,Omega,NBTI,AccDes,AccXC,ICtDF1,
     $      IExCor,IRAP,ScaHFX,ScaDFX,IRanWt,IRanGd,ICorTp,ISCRF,IAn,
     $      IAtTyp,AtmChg,C,IAtFrg,NAtoms,IRwH,IRwCA,IOCB,IRwPA,IRwPB,
     $      .True.,InStor,IDenFit,IODFit,IOFDen,DoInc,NSaved,NDimCP,CP,
     $      IOPSav,IOFSav,FrcNuc,NSAOBk,NS,NS2,NS4,ISABF,RSABF,MulTab,
     $      MapSAO,R1,R2,R3,IPBC,LIPBC,CelVec,JGC,EnJunk,ECDFit,EOneA,
     $      ETwoA,ExA,EcA,ErfPA,DoneQ,gFac,fuFac,EField,NEqShl,NEqSh2,
     $      NEqBas,NAtoms,NEqAtm,RotOp,Buf2E,RNEFit,IJ,IJMap,IGWInf,
     $      IAtBtD,IRdBtD,RRdBtD,RGWBtD,X,LenX)
          EPol = -Two*ErfP
          IniPCM = .False.
          If(NOnlyS.gt.0) then
            Call FileIO(1,-ISavGF,NTTX,FA,0)
            If(NSpBlk.eq.2) Call FileIO(1,ISavGF,NTTX,FB,0)
            endIf
          endIf
        Call FileIO(2,-IRwPA,NTTX,DA,0)
        If(HaveB) Call FileIO(2,-IRwPB,NTTX,DB,0)
        endIf
      If(IHMeth.eq.14) ETwo = Two*ETwo + ENJ
      If(IPrint.gt.0) then
        EG2 = ETwo - Ec - Ex
        Write(IOut,1260) JGC, EG2
        If(DFT.and..not.DoFkG) Write(IOut,1270) Ex, Ec, JGC, ETwo
        endIf
C
C     Self-Consistent Reaction Field (SCRF) calculations, either using
C     SCI-PCM or Onsager.
C
      If(ISCRF.ne.0.and..not.DoFkG.and..not.DidPCM) then
        If(RComp) Call GauErr('RComp and SCRF.')
        If(SCIPCM) then
          If(DoSurf.or.ItSCRF.le.3) then
            Call SCISCF(IOut,IPrint,IRwPA,IRwPB,IRwCA,IRwCB,JCycle,
     $        NBasis,NBas6D,DA,DB,FA,FB,NAtoms,AtmChg,C,(NSpBlk-1),NAEX,
     $        NBE,ERFP,IPass,DoSurf,IRwSrf,RdPCM,Epsi,AccDes,SymCut,
     $        IPFlag,AllowP,MaxTyp,Buf2E,MDV)
            If(ItSCRF.ge.4) then
              If(ItSCRF.eq.5.and.HavSrf) then
                If(DoPrnt)
     $            Write(IOut,'('' Average SCRF with previous value.'')')
                Call FileIO(2,-IRwSCH,NTT,DA,0)
                Call ACASB(NTT,FA,DB,FA,-Pt5)
                Call ACASB(NTT,FB,DB,FB,-Pt5)
                Call ACASB(NTT,FA,DA,FA,Pt5)
                Call ACASB(NTT,FB,DA,FB,Pt5)
                Call AAdd(NTT,DA,DB,DB)
                Call AScale(NTT,Pt5,DB,DB)
                endIf
              If(IRwSCH.eq.0) IRwSCH = IGetRW(-NTT,0)
              Call FileIO(1,-IRwSCH,NTT,DB,0)
              HavSrf = .True.
              endIf
          else
            If(DoPrnt) Write(IOut,'('' Recover SCRF 1e matrix.'')')
            Call FileIO(2,-IRwSCH,NTT,DB,0)
            Call AAdd(NTT,FA,DB,FA)
            Call AAdd(NTT,FB,DB,FB)
            endIf
        else if(PCM) then
          If(NRI.eq.2) Call GauErr('No PCM with complex.')
          Call FileIO(2,-IRwPA,NTT,DA,0)
          If(HaveB) Call FileIO(2,-IRwPB,NTT,DB,0)
          If(RComp) Call GauErr('PCMQM not ready for RComp.')
          FirstC = JCycle.eq.1
          LastC  = JCycle.eq.MaxCyc
          Call PCMQM(IOut,IPrint,(NSpBlk-1),FirstC,.True.,.True.,LastC,
     $      FirstC,NOpUse,IPFlag,AllowP,FMM,FMFlag,FMFlg1,NFxFlg,IHMeth,
     $      NAEX,NBE,NAtoms,NAtoms,NBasis,NBas6D,ICharg,.False.,RJunk,
     $      RJunk,ErfP,DoneQ,IAn,IAtTyp,AtmChg,C,DA,DB,FA,FB,X,LenX)
        else if(Onsagr(ISCRF).or.SolEPt(ISCRF)) then
          If(NRI.eq.2) Call GauErr('No Onsager/SCS with complex.')
          Call FileIO(2,-IRwPA,NTT,DA,0)
          If(HaveB) Call FileIO(2,-IRwPB,NTT,DB,0)
          If(Onsagr(ISCRF)) then
            Call Onsag2(IOut,IPrint,1,ISCRFT,1,NAtoms,AtmChg,C,NBasis,1,
     $        NTT,gFac,fuFac,EPol,DA,DB,FA,FB,VA,EField)
          else
            Call SolPMt(IOut,IPrint,1,1,IRwH,NBasis,1,NTT,EPol,DA,DB,FA,
     $        FB,VA)
            endIf
          ErfP = -Pt5*EPol
          endIf
        Call TStamp(1,'After SCRF')
      else if(.not.DidPCM) then
        ErfP = Zero
        endIf
      If(IPrint.gt.0.and.ErfP.ne.Zero) Write(IOut,1275) ErfP
C
C     Finish calculating the energy.
C
      If(IFlagO.eq.0) EnOld = Energy
      Energy = EOne + ETwo + ErfP + ENR
      If(JCycle.eq.1) EPert = Energy
      DelEn = Energy - EnOld
      OldRis = Rises
      AboveM = Energy.gt.(EnMin+RMedEn)
      Rises = JCycle.gt.1.and.Temp.eq.Zero.and.
     $  ((FixRis.eq.2.and.Energy.gt.(EnOld+SmalEn)).or.
     $   (FixRis.ge.3.and.Energy.gt.(EnMin+SmalEn)))
      EnMin = Min(Energy,EnMin)
      If(SpDiag) then
        If(DoFkG) then
          IFockA = IFStrt
          IFockB = IFockA + IOpCl*NTT6D
          IEndFG = IFockB + NTT6D
          Call TstCor(IEndFG,LenX,'SpDiag-FkGues')
          Call FileIO(2,-IRwFA,NTTX,X(IFockA),0)
          If(IOpCl.eq.1) then
            Call FileIO(2,-IRwFB,NTTX,X(IFockB),0)
            Call MkSprL(X,LenX,0,0,0,IFockB,jifb,ifb,IFATmp,NBasis,
     $        NZFB,'FoFkG',ZThrP,1)
            Call AMove(NTT,X(IFockA),X(IFATmp))
            Call SpAMvM(X,LenX,IFStrt,IFockA,jifb,ifb,NBasis,2)
            Call AMove(NTT,X(IFATmp),X(IFockA))
            endIf
          Call MkSprL(X,LenX,0,0,0,IFockA,jifa,ifa,IAftrF,NBasis,NZFA,
     $      'FoFkG',ZThrP,1)
          If(IOpCl.eq.0) then
            jifb = jifa
            ifb = ifa
          else
            Call MkSmFm(X,LenX,0,jifa,ifa,0,jifb,ifb,IFStrt,IFStrt,0,
     $        jifa,ifa,ifb,NZF,IAftrF,NBasis,2)
            endIf
        else
          Call Cmpct2(IOpCl,.True.,NBasis,ZThrP,X(jifa),IFA,IFB,NZFA,X)
          endIf
        Call SpStor(IRWFA,X(jifa),X(ifa),NBasis,1,.True.,0)
        If(.not.RHF) Call SpStor(IRWFB,X(jifa),X(ifb),NBasis,1,.True.,0)
      else if(DoFkG) then
        Call FileIO(2,-IRwFA,NTTX,FA,0)
        If(HaveB) Call FileIO(2,-IRwFB,NTTX,FB,0)
      else
        Call FileIO(1,-IRwFA,NTTX,FA,0)
        If(HaveB) Call FileIO(1,-IRwFB,NTTX,FB,0)
        endIf
      If(IPrint.gt.2) then
        If(SpDiag) then
          Call SpDisp(IOut,X(jifa),X(ifa),NBasis,'Alpha Fock Matrix:',2)
          If(.not.RHF) Call SpDisp(IOut,X(jifa),X(ifb),NBasis,
     $      'Beta Fock Matrix:',2)
        else
          Call LTOutX(IOut,0,IOpClY,NMtPBC,0,'Fock matrix',NBasis,
     $      NBasis,FA,FB,0)
          endIf
        endIf
      If(DoPrnt) then
        If(FkGues.and.IPass.eq.0) then
          Write(IOut,1170)
        else if(IFlag.eq.0.and.JCycle.gt.1) then
          Write(IOut,1120) Energy, DelEn, Rises, Damp
        else if(IFlag.eq.0) then
          Write(IOut,1130) Energy
        else
          Write(IOut,1170)
          endIf
        endIf
      If(ICampK.ne.0) Call SaveCK(NSpBlk,IndCK,Len1CK,NBNO,NTTX,IRwCA,
     $  IRwCB,IRwFA,IRwFB,IRwCK,Energy,EnerCK,X,LenX)
      IFlagO = IFlag
      If(DoHirS.or.(IPrt1.ge.1.and.IOpCl.gt.3.or.
     $  (IPrint.ge.1.and.NSpBlk.eq.2))) then
        Call Spin(IOpCl,IRwS,IRwCA,IRwCB,IRwUR,NBasis,NBsUse,NAE,
     $    NBE,Sx,Sy,Sz,S2,S,S2Cor,X,LenX)
        Write(IOut,1380) Sx, Sy, Sz, S2, S
        If(DoHirS) then
          Call FileIO(2,-IRwPA,NTTX,DA,0)
          If(HaveB) then
            Call FileIO(2,-IRwPB,NTTX,DB,0)
          else
            Call AMove(NTTX,DA,DB)
            endIf
          Call AtMMoH(IOut,IPrint,IOpCl,IRadAn,IRanWt,IRanGd,NAtoms,
     $      IAtTyp,IAn,AtmChg,C,NBasis,DA,DB,X,LenX)
          endIf
        endIf
      If(JCycle.gt.1.and.IHMeth.eq.14) then
        Exit = .True.
        Goto 100
        endIf
C
C ANT BEGIN modification
C
C         Call ANT for evaluating density matrix
C
      call ANT (UHF,JCycle,IRwH,IRwPA,IRwPB,IRwFA,IRwFB,
     $                IRwS1,IRwEig,Acurcy,Crit,ANTOn,NBasis)

      if(ANTOn) DIIS=.false.
      if(ANTOn) IncFoc=.false.
      if(ANTOn) IDIIS=1
      if(ANTOn) Damp=.false.
C
C ANT END modification
C
C
C     Generate gamma point Fock matrices for extrapolation.
C
      If(UsKPnt) Call MkGamF(IOut,IPrint,'raw',NMtPBC,NBasis,NTT,IRwFA,
     $  IRwFB,IRwF1X,IRwF2X,FA,FB)
C
C     Turn alpha and beta Fock matrices into new coefficients.  This
C     is the only place where we care whether we're doing ROHF vs UHF.
C
      If(ROHF.and..not.RHF) Call ROFock(IOut,IPrint,NBasis,NOrbs,NAEX,
     $  NBE,IRwF1,IRwCA,IRwS1,FA,FB,D,V,Scr2)
      If(JCycle.eq.LvFirst.and.IStab.ne.2) then
        IV1 = IVShft
        Temp = TempSt
        NIter = -1
        NCycP0 = JCyc1 - 1
      else if(JCycle.ge.LvFirst) then
        IV1 = IVShft
      else
        IV1 = 0
        endIf
      If((Temp.ne.Zero.and.LvWFON.ne.1).or.(AboveM.and.LvRise.ne.1))
     $  IV1 = 0
      DoDiag = .not.(IStab.eq.2.and.JCycle.lt.LvFirst)
      If(SpDiag) Call FOrthx(IOut,IPrint,X,LenX,1,IRwF1,IRwF2,IRwFOA,
     $  IRwFOB,IRwPA,IRwP2,IRwPOA,IRwPOB,IRwV,IRwH,IRwFM1,IRwFM2,NBasis,
     $  ZThrP,IMemOp,IFxFm,NAEX,NBE,ZDO,JCycle)
C     A hack to just compute the energy and exit.
      If(JustEn) then
        Exit = .True.
        Goto 100
        endIf
C
C     If we picked up the alternate initial guess and it was worse than
C     the regular (projected) guess, then switch the data for ConOpn so
C     that damping is with respect to changes from the better density.
C     Also, turn off EDIIS, which seems to be confused by the alternate
C     guess.
C
      If(DidAlG.and.Energy.gt.EnOld) then
        If(IPSave(0).eq.0) Write(IOut,1310)
        IDITyp = Min(IDITyp,1)
        NumUp = -1
        Call FileIO(2,-IRwC3,NTTXB,VA,0)
        Call ANeg(NTTXB,VA,VA)
        Call FileIO(1,-IRwC3,NTTXB,VA,0)
        Call FileIO(2,-IRwC1,NTTXB,VB,0)
        Call AAdd(NTTXB,VA,VB,VA)
        Call FileIO(1,-IRwC1,NTTXB,VA,0)
        endIf
      If(IFlag.eq.0.and.DoDIIS) Call SCFDIX(IOut,IPtDIS,IDITyp,IFSelI,
     $  RComp,SpDiag,JCycle,LvFirst,Rises,FixRis,OldRis,IStab,MxSCFC,
     $  NIter,IOpClD,NBasis,NBsUse,ISCRF,IRwP1,IRwP1X,IRwF1,IRwF1X,
     $  IRwP2,IRwP2X,IRwF2,IRwF2X,IRwS1,IRwTM,IRwFOA,IRwFOB,IRwPOA,
     $  IRwPOB,IRWFM1,IRWFM2,ZThrP,NSvDIS,IOBDIS,IOPADI,IOPBDI,IOFADI,
     $  IOFBDI,IOErDI,D,F,V,VV,BMat,Scr2,ErrMax,NErSav,IndFok,Energy,
     $  EnCyc,ErrCyc,IMemOp,NMtPBC,OldCof,DISPCM,IFrgTp,NBlks,NBF,NBFU,
     $  MapF2M,NOccSy,IOTMA,IOTMB,IndTMB,IOCAS,IOCBS,IndCBS,IOSBlk,
     $  IndSBB,X,LenX)
      If(NSvDIS.gt.0.and..not.SpDiag) CofLast = Scr2(NSvDIS+1,1)
      IFSelI = 0
      If(.not.RHF) Call FreeRW(IRWFM2)
      If(SpDiag) then
        ITestX = 0
      else if(ISyOcc.ge.2) then
        ITestX = 2
      else
        ITestX = 1
        endIf
      Call TstCor(ITestX*NBSq,LenX,'CycOpn-DgFock')
      If(IFBrod.eq.0) then
        If(JCycle.lt.LvFirst) then
          IPFON = 0
        else
          IPFON = 2
          If(JCycle.gt.1.and.NAEX.lt.NOrbs) then
            Call FileIO(2,-IRwEig,NSpBlk*NOrbs,Scr2,0)
            GapA = Scr2(NAEX+1,1) - Scr2(NAEX,1)
            If(NDimBl.eq.1) then
              GapB = Scr2(NAE+1+(NSpBlk-1)*NOrbs,1)
     $          - Scr2(NAE+(NSpBlk-1)*NOrbs,1)
            else
              GapB = GapA
              endIf
            If(GapA.ge.GapMin.and.GapB.ge.GapMin.and.
     $        (GapA.le.GapMax.or.GapB.le.GapMax)) IPFON = 1
            endIf
          endIf
      else
        IPFON = IFBrod
        endIf
      NS1 = NSAStep
      If(IPFON.eq.1) NS1 = 2*NS1
      If(NS1.gt.0)
     $  Temp = Min(Temp,TempSt*GFloat(NS1-JCycle+LvFirst)/GFloat(NS1))
      Temp = Max(Temp,Zero)
C
C     Form gamma-F for diagonalization.
C
  500 If(UsKPnt) Call MkGamF(IOut,IPrint,'for diagonalization',NMtPBC,
     $  NBasis,NTT,IRwFA,IRwFB,IRwF1X,IRwF2X,FA,FB)
C
C     Now apply the level shift and diagonalize either one (ROHF) or
C     two (UHF) Fock matrices.
C
C
C ANT BEGIN modification
C
C         Only diagonalize Fock matrix when leads not yet connected
C
      if (.NOT.ANTOn) then
C
C ANT END modification
C
      DidAlG = .False.
      If(DoneQ.and.NOnlyS.gt.0) IOnlyS = NOnlyS - 1
      If(ICampK.ne.0.and.JCycle.ge.2) Call CamKin(IOut,IPrint,IOpCl,
     $  ICampK,NAE,NBE,NBsUse,NBasis,IndCK,IRwS,IRwCK,Len1CK,IRwCA,
     $  IRwCB,CKThr,EnerCK,CKLam,X,LenX)
      If(JCycle.eq.1.and.OKAltG.and.ErrMax(0).ge.AlGsMn) then
        If(DoPrnt) Write(IOut,1300)
        DidAlG = .True.
        If(NExDmp.ge.JCycle) NExDmp = NExDmp + 1
        Call FileIO(2,-IOGesA,NBNO,DA,0)
        Call FileIO(1,-IRwCA,NBNO,DA,0)
        If(NDimBl.eq.1) then
          Call FileIO(2,IOGesA,NBNO,DA,0)
          Call FileIO(1,-IRwCB,NBNO,DA,0)
          endIf
        Call FileIO(2,IOGesA,NTTX,DA,0)
        Call FileIO(1,-IRwPA,NTTX,DA,0)
        If(NDimBl.eq.1) then
          Call FileIO(2,IOGesA,NTTX,DA,0)
          Call FileIO(1,-IRwPB,NTTX,DA,0)
          endIf
      else if(UsKpnt) then
        If(IPrint.gt.0) then
C         Check trace of overlap with Fock
          Call FileIO(2,-IRwS,NTTX,FA,NTT)
          Call FileIO(2,-IRwPA,NTTX,DA,0)
          TrPFA =
     $      SCFTrX(.False.,.False.,IOpClX,NMtPBC,NBasis,DA,XX,FA,XX)
          If(HaveB) Call FileIO(2,-IRwPB,NTTX,DA,0)
          TrPFB =
     $      SCFTrX(.False.,.False.,IOpClX,NMtPBC,NBasis,DA,XX,FA,XX)
          Write(IOut,1320) TrPFA, TrPFB
          endIf
        If(DoFkG.or.IVShft.eq.-1.or..not.FkGues) IV1 = 0
        Call ZDFock(IOut,IPrint,IPBC,NSpBlk-1,NBasis,NAEX,NBE,IRwFA,
     $    IRwFB,IRwPA,IRwPB,IRwPOA,IRwPOB,IRwCAC,IRwCBC,IRwEiC,IRwV,
     $    IRwPAN,IRwPBN,NKPnt,NBFK,IADMP,ICoor,CmlxF,WKPnt,IV1,VSGoal,
     $    VSMax,SDif,X,LenX)
        SDifT = FLenRc*SDif
        If(SDifT.gt.ThrSDi.and..not.DoFkG) then
          Write(IOut,1000) FLenRc, SDif, SDifT, ThrSDi
C         Call GauErr('Diagonalization produced ghost orbitals.')
          endIf
        Call ZIntLT(IOut,IPrint,NSpBlk-1,NBasis,NAEX,NBE,IRwEiC,NKPnt,
     $    ICoor,IOToK,NPDir,NTot,WKPnt,NBFK,IRwWgt,IZMthd,FKTemp,kTLim,
     $    NSPnt,If0K,IZWrk,X,LenX)
        If(IZWrk.gt.0) Call ZFormD(IOut,IPrint,IPBC,.False.,NSpBlk-1,
     $    NBasis,IZWrk,IRwCAC,IRwCBC,IRwPA,IRwPB,IRwEiC,IRwWgt,NKPnt,
     $    NBFK,NAEX,NBE,WKPnt,.True.,ICoor,CmlxF,X,LenX)
        If(NSpBlk.eq.1.and.HaveB) then
          Call FileIO(2,-IRwPA,NTTX,DA,0)
          Call FileIO(1,-IRwPB,NTTX,DA,0)
        else if(NSpBlk.eq.1) then
C         Save beta Fock matrix for possible later use
          Call FileIO(2,-IRwFA,NTTX,DA,0)
          Call FileIO(1,-IRwFB,NTTX,DA,0)
          endIf
C       Check density the total density for linear dependencies
        Call FileIO(2,-IRwS,NTTX,FA,NTT)
        Call FileIO(2,-IRwPA,NTTX,DA,0)
        TrPSA =
     $    SCFTrX(.False.,.False.,IOpClX,NMtPBC,NBasis,DA,XX,FA,XX)
        If(HaveB) Call FileIO(2,-IRwPB,NTTX,DA,0)
        TrPSB =
     $    SCFTrX(.False.,.False.,IOpClX,NMtPBC,NBasis,DA,XX,FA,XX)
C       Crash the job if the density is screwed up somehow
        TrTot = TrPSA + TrPSB
        TrErr = TrTot - RNE
        TrCut = GFloat(100)*MDCutO(0)*RNE
        If(Abs(TrErr).gt.TrCut) then
          Write(IOut,1340) TrTot, RNE, TrErr, TrCut
          Call GauErr('Numerical problems in PBC.')
          endIf
        Call DagGam(NSpBlk-1,NBasis,NAEX,NBE,IRwF1X,IRwF2X,IRwP1,IRwP2,
     $    IRwCA,IRwCB,IRwEig,NOrbs,IV1,VSGoal,VSMax,X,LenX)
      else if(DoDMS.and.IPass.ne.2) then
        Call QNDMS(IOut,2,NBasis,IRwPA,IRwFA,0,0,IRwV,NAe,X,LenX,DMSErA,
     $    IADBas)
        If(.not.RHF) Call QNDMS(IOut,2,NBasis,IRwPB,IRwFB,0,0,IRwV,NBE,
     $    X,LenX,DMSErB,IADBas)
      else if(ROHF.and.IOnlyS.eq.(NOnlyS-1)) then
        Call DgFock(IOut,IPrint,0,NRI,NDimBl,ROHF,ZDO,RComp,IV1,VSGoal,
     $    VSMax,VSGpMn,IDiagX,IAlDNI,NPassI,ISyOcc,IWfnSy,DoDiag,NBlks,
     $    NBasis,NBD,NOrbs,NAEX,IRwF1,IRwCA,IOTMEx,IOTMA,0,IOCAS,0,
     $    IRwEig,0,IRwPA,AccDag,NBF,NBFU,NOccSy,CutOff,'Alpha',JCycle,
     $    IPFON,.True.,TempSt,ErrMax,DLaMax,NAtoms,C,LLim,IULim,IOpDME,
     $    SECut,DThr,ZThrP,PFThr,NEqBas,NOp1,IRwPOA,IRwFOA,IRwS,IRwV,
     $    Temp,Beta,EF,VShftA,IRWFM2,IMemOp,IFxFm,IOpDM2,F,D,V,Scr2,
     $    LScr2-1,Scr2(1,2),X,LenX,SpDiag,NAEBr,NOABr,GapD,OVLim,OVMax,
     $    IOSymA,0,ICampK,DgAllA,ISwCut,IPsCut,IFrgTp,MapF2M,Y,LenY)
        If(IDiagX.eq.1) LastDg = IAbs(NIter)
        VShftB = VShftA
        DgAllB = DgAllA
        If(.not.RHF) then
          Call FrmPSA(IOut,IPrint,'Beta',.True.,NRI,IPFON,NBD,NOrbs,
     $      NBE-NBEBr,NBEBr,NOBBr,Temp,Beta,ErrMax,EF,Scr2,V,D,
     $      Scr2(1,2))
          If(IPrint.ge.3) Call LTOutS(IOut,'Beta density at cycle',
     $      JCycle,NBasis,D,0)
          endIf
        If(.not.SpDiag) then
          Call FileIO(1,-IRwEig,NOrbs,Scr2,NOrbs)
          Call FileIO(1,-IRwCB,NBNO,V,0)
          Call FileIO(1,-IRwPB,NTTX,D,0)
          endIf
      else if(IOnlyS.eq.(NOnlyS-1)) then
        Call DgFock(IOut,IPrint,1,NRI,NDimBl,ROHF,ZDO,.False.,IV1,
     $    VSGoal,VSMax,VSGpMn,IDiagX,IAlDNI,NPassI,ISyOcc,IWfnSy,DoDiag,
     $    NBlks,NBasis,NBD,NOrbs,NAE,IRwFA,IRwCA,IOTMEx,IOTMA,0,IOCAS,
     $    0,IRwEig,0,IRwPA,AccDag,NBF,NBFU,NOccSy,CutOff,'Alpha',JCycle,
     $    IPFON,.True.,TempSt,ErrMax,DLaMax,NAtoms,C,LLim,IULim,IOpDME,
     $    SECut,DThr,ZThrP,PFThr,NEqBas,NOp1,IRwPOA,IRwFOA,IRwS,IRwV,
     $    Temp,Beta,EF,VShftA,IRWFM2,IMemOp,IFxFm,IOpDM2,F,D,V,FF,LScr2,
     $    Scr2,X,LenX,SpDiag,NAEBr,NOABr,GapD,OVLim,OVMax,IOSymA,0,
     $    ICampK,DgAllA,ISwCut,IPsCut,IFrgTp,MapF2M,Y,LenY)
        Call DgFock(IOut,IPrint,2,NRI,NDimBl,ROHF,ZDO,.False.,IV1,
     $    VSGoal,VSMax,VSGpMn,IDiagX,IAlDNI,NPassI,ISyOcc,IWfnSy,DoDiag,
     $    NBlks,NBasis,NBD,NOrbs,NBE,IRwFB,IRwCB,IOTMEx,IOTMB,IndTMB,
     $    IOCBS,IndCBS,IRwEig,NOrbs,IRwPB,AccDag,NBF,NBFU,NOccSy,
     $    CutOff,'Beta',JCycle,IPFON,.False.,TempSt,ErrMax,DLaMax,
     $    NAtoms,C,LLim,IULim,IOpDME,SECut,DThr,ZThrP,PFThr,NEqBas,NOp1,
     $    IRwPOB,IRwFOB,IRwS,IRwV,Temp,Beta,EF,VShftB,IRWFM2,IMemOp,
     $    IFxFm,IOpDM2,F,D,V,FF,LScr2,Scr2,X,LenX,SpDiag,NBEBr,NOBBr,
     $    GapDB,OVLim,OVMaxB,IOSymA,InToWP(NBsUse),ICampK,DgAllB,ISwCut,
     $    IPsCut,IFrgTp,MapF2M,Y,LenY)
        If(IDiagX.eq.1) LastDg = IAbs(NIter)
        GapD = Min(GapD,GapDB)
        OVMax = Max(OVMax,OVMaxB)
        endIf
      If(SaveMO.and..not.SpDiag) Call ChkSCF(1,NSpChk,NRI*NBD,NOrbs,
     $  IRwEig,IRwCA,IRwCB,0,0,X,LenX)
      DoSurO = DoSurf
      Call FreeRW(IRWFM2)
      Call TStamp(1,'End of iteration.')
C
C ANT BEGIN modification
C
      Call FileIO(2,-IRwGen,LenGen,DGen,0)
C
      end if
C
C ANT END modification
C
      If(IPass.eq.2.and.IFinIt.eq.3) then
        IConvr = 2
      else
        Goto 10
        endIf
C
  999 If(PCM) then
        Call FileIO(2,-IRwPA,NTTX,DA,0)
        Call FileIO(2,-IRwPB,NTTX,DB,0)
        Call PCMQM(IOut,IPrint,(NSpBlk-1),.False.,.True.,.True.,.True.,
     $    .False.,NOpUse,IPFlag,AllowP,FMM,FMFlag,FMFlg1,NFxFlg,IHMeth,
     $    NAEX,NBE,NAtoms,NAtoms,NBasis,NBas6D,ICharg,.False.,RJunk,
     $    RJunk,ErfP,DoneQ,IAn,IAtTyp,AtmChg,C,DA,DB,FA,FB,X,LenX)
        endIf
      Call FreeRW(IOPSav)
      Call FreeRW(IOFSav)
      Call FreeRw(IRWFM1)
      Call FreeRw(IRWFM2)
      Call FreeRw(IRwCK)
      Call FreeRw(ISavGF)
      Call FreeRw(IOSBlk)
      Call FreeRw(IOSInv)
      Call FreeRw(IRwOld)
      Call ILSW(1,5,IConvr)
      If(UsKPnt) then
        Call FreeRw(IRwF1X)
        Call FreeRw(IRwF2X)
        Call FreeRw(IRwP1)
        Call FreeRw(IRwP2)
        endIf
C
C     Make sure the fitting density left behind is for the full
C     and not an incremental density.
C
      If(IHMeth.eq.14) then
        Call ASet(NDBF6D,Pt5,DA)
        Call ConWrt(IOFDen,NDBF6D,DA)
      else if(.not.SpDiag.and.(IncFoc.or.IDenFit.lt.0).and.
     $  (IConvr.eq.2.or..not.Halt)) then
        Call FileIO(2,-IRwPA,NTTX,DA,0)
        Call FileIO(2,-IRwPB,NTTX,DB,0)
        Call ACSASB(NTTX,DA,Pt5,DB,Pt5,DA)
        IDFitA = Abs(IDenFit)
        TolFit = Crit1 / (Ten*Ten)
        RNEFit(1) = RNE
        Call DenFit(IOut,IPrint,IDFitA,.True.,FrcNuc,0,1,NAtoms,NBasis,
     $    NBas6D,FMM,FMFlag,FMFlg1,NFxFlg,IODFit,IOFDen,IPFlag,AllowP,
     $    LSEAll,RSEAll,IAn,IAtTyp,AtmChg,C,DoPurP,I1Cent,ECDFit,AccDes,
     $    JSym2X,NOpUse,NOp1,NEqShl,NEqSh2,RotOp,NEqAtm,NEqBas,TolFit,
     $    NBTI,RNEFit,DA,FB,IPBC,LIPBC,CelVec,Buf2E,MDV)
        endIf
      If(IDenFit.lt.0.or.IPrint.ge.2) then
        NDBF = ITqry(IOFDen)
        If(NDBF.gt.0) then
          Call FileIO(2,-IOFDen,NDBF,DA,0)
          Call OutMtS(IOut,'Final P fit:',0,0,DA,NDBF,1,NDBF,1)
          endIf
        endIf
      EOneO = EOne
      ETwoO = ETwo
      ErfPO = ErfP
      Return
      End
*Deck DavGts
      Subroutine DavGts(IOut,IPrint,V,MDV,istart,jif,iF,N,xmax,xmin)
      Implicit Real*8(A-H,O-Z)
C
C     Calculate the lowest and highest eigenvectors in the
C     Spectrum of the sparse, symmetric matrix F (jif,iF)
C
      Dimension V(*),IPars(3),JJ(1)
      External SpAXMt, GenQ
      Save Zero,Conv,TenD4,JJ
      Data Zero/0.d0/,Conv/0.01d0/,TenD4/10000.d0/,JJ/0/
C
      NZ = IGetVl(V(jif),N+1)-1
      IPars(1) = jif
      IPars(2) = iF
      MaxMI = 20
      MSek1 = 2
      iWhenS1 = 0
      ifact1 = 1
      IHDiag = IGetRW(N,0)
      IEigen = IGetRW(MSek1,0)
      IGuess = IGetRW(N*MSek1,0)
      xmin  = Zero
      xmax  = Zero
      IRow  = IStart
      I1    = IRow + InToWP(N+1)
      NGott = MDV - I1 + 1
      call tstcor((i1+MSek1*N),MDV,'DavGts-DavDg2-1')
      Call ISeq(N,0,1,V(IRow))
      IPars(3) = IRow
      Do 10 i=0, (MSek1*N-1)
   10   V(I1+i)=gfloat(IRNumb(10000))/TenD4
      Call Fileio(1,-IGuess,MSek1*N,V(I1),0)
      Call Fileio(1,-IHDiag,N,V(iF),0)
      Call DavDag(IOut,IPrint,0,-1,N,MaxMI,MSek1,MSek1,MSek1,0,iWhenS1,
     $  IFact1,.False.,IHDiag,IEigen,IGuess,0,Conv,Conv,SpAXMt,GenQ,
     $  IPars,V,0,0,JJ,V(I1),NGott)
      Call Fileio(2,-IEigen,1,XMin,0)
      Do 50 i=0, (MSek1*N-1)
   50   V(I1+i) = gfloat(IRNumb(10000))/TenD4
      Call Fileio(1,-IGuess,MSek1*N,V(I1),0)
      Call ANeg((N+NZ),V(iF),V(iF))
      Call AInc(N,XMin,V(iF),V(iF))
      Call DavDag(IOut,IPrint,0,-1,N,MaxMI,MSek1,MSek1,MSek1,0,iWhenS1,
     $  IFact1,.False.,IHDiag,IEigen,IGuess,0,Conv,Conv,SpAXMt,GenQ,
     $  IPars,V,0,0,JJ,V(I1),NGott)
      Call AInc(N,-XMin,V(iF),V(iF))
      Call ANeg((N+NZ),V(iF),V(iF))
      Call Fileio(2,-IEigen,1,XMax,0)
      xmax = -(-xmin+xmax)
      Return
      End
*Deck dBkdRk
      Subroutine dbkdrk(NTot,IDid,NDo,GNorm,EFG,dndr,dbdr,dphedr,
     $  rhatx,rhaty,rhatz,dbndr,dphndr,ISlMap,IScal)
      Implicit Real*8 (A-H,O-Z)
      Dimension GNorm(NTot,3),EFG(NTot,10,2),dndr(NDo,3),dbdr(*),
     $  dphedr(*),rhatx(*),rhaty(*),rhatz(*),dbndr(*),dphndr(*),
     $  ISlMap(NTot)
C
c     db(k)/dr(k) and dV(k)/dr(k). V = Phi = potential.
c     b = -E.n
c
c     If IScal=0 then EFGe contains fields etc from nuc and elec.
c     If IScal=2 then do each part seperately.
c
      Do 10 I = 1, NDo
        IIDid=I+IDid
        ISDid = ISlMap(IIDid)
        If(ISDid.ne.0) then
          part1 = rhatx(IIDid)*EFG(ISDid,10,1)*GNorm(ISDid,1)+
     $            rhatx(IIDid)*EFG(ISDid,9,1)*GNorm(ISDid,2)+
     $            rhatx(IIDid)*EFG(ISDid,8,1)*GNorm(ISDid,3)+
     $            rhaty(IIDid)*EFG(ISDid,9,1)*GNorm(ISDid,1)+
     $            rhaty(IIDid)*EFG(ISDid,7,1)*GNorm(ISDid,2)+
     $            rhaty(IIDid)*EFG(ISDid,6,1)*GNorm(ISDid,3)+
     $            rhatz(IIDid)*EFG(ISDid,8,1)*GNorm(ISDid,1)+
     $            rhatz(IIDid)*EFG(ISDid,6,1)*GNorm(ISDid,2)+
     $            rhatz(IIDid)*EFG(ISDid,5,1)*GNorm(ISDid,3)
          part2 = EFG(ISDid,4,1)*dndr(i,1)+EFG(ISDid,3,1)*dndr(i,2)+
     $            EFG(ISDid,2,1)*dndr(i,3)
          dbdr(I)=-(part1+part2)
          dphedr(I)=-(EFG(ISDid,4,1)*rhatx(IIDid)+
     $      EFG(ISDid,3,1)*rhaty(IIDid)+EFG(ISDid,2,1)*rhatz(IIDid))
c
c
          If(IScal.eq.2) then
            parta = rhatx(IIDid)*EFG(ISDid,10,2)*GNorm(ISDid,1)+
     $              rhatx(IIDid)*EFG(ISDid,9,2)*GNorm(ISDid,2)+
     $              rhatx(IIDid)*EFG(ISDid,8,2)*GNorm(ISDid,3)+
     $              rhaty(IIDid)*EFG(ISDid,9,2)*GNorm(ISDid,1)+
     $              rhaty(IIDid)*EFG(ISDid,7,2)*GNorm(ISDid,2)+
     $              rhaty(IIDid)*EFG(ISDid,6,2)*GNorm(ISDid,3)+
     $              rhatz(IIDid)*EFG(ISDid,8,2)*GNorm(ISDid,1)+
     $              rhatz(IIDid)*EFG(ISDid,6,2)*GNorm(ISDid,2)+
     $              rhatz(IIDid)*EFG(ISDid,5,2)*GNorm(ISDid,3)
            partb = EFG(ISDid,4,2)*dndr(i,1)+EFG(ISDid,3,2)*dndr(i,2)+
     $              EFG(ISDid,2,2)*dndr(i,3)
            dbndr(I)=-(parta+partb)
            dphndr(I)=-(EFG(ISDid,4,2)*rhatx(IIDid)+
     $        EFG(ISDid,3,2)*rhaty(IIDid)+EFG(ISDid,2,2)*rhatz(IIDid))
            endif
          endif
   10   Continue
      Return
      End
*Deck DetFit
       Subroutine DetFit(ItMax,IOut,Iprint,CC,Ncore,Norbs,NElect,jif,if,
     $  NZf,NZpA,jipA,ipA,IEnd,Niter,Trace,Thresh,ConvCr,xIdemp,ImemOp,
     $  OuterF,Fixed,FixTyp,CNSI)
      Implicit Real*8(A-H,O-Z)
C
C     This routine calculates the Fermi level of the input Fockian.
C
      Parameter(maxf=20)
      Logical DoHalf,Fixed,FixTyp,CNSI,OuterF
      Real*8 CC(Ncore),fmax,EFermi,TUpper,TLower,XUpper,Xlower
      Save fmax,EFermi,TUpper,TLower,XUpper,Xlower,Half,Zero,Two
      Data Half/0.5d0/,Zero/0.d0/,Two/2.d0/
   34 Format(' Canonical NSI')
   35 Format(' Grand canonical NSI')
   37 Format(' Looking for new Fermi level')
   38 Format(' XUpper = ',f20.10,' Xlower = ',f20.10)
   39 Format(' TUpper = ',f20.10,' TLower = ',f20.10)
   40 Format(' Half distance')
   41 Format(' Fermi level found in',i4,' steps')
 1234 Format(' Fermi level at ',f10.5)
 1111 Format(' istep = ',i3,' EFermi = ',f12.4,' Trace p = ',f12.4,
     $ ' NElect = ',i6,' Dens = ',f8.4)
C
C     Bounds from Gershgorin - algorythm
C
      IX = IEnd
      IEnd = IX + Norbs
      Call GBnds(CC(jif),CC(if),Norbs,NZf,fmax,ZU,ZL,IOut,CC(IX))
      IEnd = IX
      fTrace = gfloat(NElect)
      fNorbs = gfloat(Norbs)
      XUpper = ZU
      TUpper = fNorbs
      Xlower = ZL
      TLower = Zero
C
C     Iterate for Fermi level
C
      If(Niter.eq.1) EFermi = Zero
      DoHalf = .True.
      If(IPrint.ge.1) then
        If(CNSI) then
          Write(IOut,34)
        else
          Write(IOut,35)
          endIf
        endIf
      III = 0
      IEndii = IEnd
      Dens = Zero
      Do 100 II = 1, MaxF
        If(.not.CNSI) then
          If(IPrint.ge.1) Write(IOut,1234) EFermi
C         Make a copy of fock mtr.
          jifc = jif
          ifc = if
          NZfc = NZf
          Istart = IEndii
          Call SpAMvM(CC,Ncore,Istart,IEnd,jifc,ifc,Norbs,2)
C
C         Shift f by Fermi level
C
          Do 723 i = 0, Norbs - 1
  723       CC(ifc+i) = CC(ifc+i) - EFermi
        else
          jifc = jif
          ifc = if
          NZfc = NZf
          endIf
        Traceo = Trace
        Call NSIter(IOut,Iprint,CC,Ncore,Norbs,jifc,ifc,NZfc,IEnd,
     $    Thresh,jipA,ipA,NZpA,ConvCr,ItMax,Niter,NElect,xIdemp,OuterF,
     $    Fixed,FixTyp,CNSI,ImemOp)
        Trace = ArrSum(Norbs,CC(ipa))
C
C       Check convergence
C
        If(Abs(Trace-fTrace).le.Half) then
          If(IPrint.ge.1.and.II.gt.1) then
            Write(IOut,1111) (II-1),EFermi,Trace,NElect,Dens
            Write(IOut,41) (II-1)
            endIf
          Return
          endIf
        If(II.eq.1.and.IPrint.ge.1) Write(IOut,37)
        Dens = (TUpper - TLower) / (XUpper - Xlower)
        Write(IOut,1111) II - 1,EFermi,Trace,NElect,Dens
C
C       Modify bounds of gap and interpolate Fermi level
C
        If(Trace.gt.fTrace) then
          XUpper = EFermi
          TUpper = Trace
        else
          Xlower = EFermi
          TLower = Trace
          endIf
        If(IPrint.ge.1) then
          Write(IOut,38) XUpper,Xlower
          Write(IOut,39) TUpper,TLower
          endIf
        x1 = Abs(Trace - fTrace)
        x2 = Abs(Traceo - fTrace)
        If(II.gt.1.and.x2/x1.lt.Two.and.DoHalf) then
          EFermi = (XUpper + Xlower)/Two
          If(IPrint.ge.1) Write(IOut,40)
          III = III + 1
          If(III.ge.2) DoHalf = .False.
        else
          DoHalf = .True.
          III = 0
          If(Trace.gt.fTrace) then
            EFermi = XUpper +
     $        (XUpper-Xlower)/(TUpper-TLower)*(fTrace+Half-Trace)
          else
            EFermi = Xlower +
     $        (XUpper-Xlower)/(TUpper-TLower)*(fTrace+Half-Trace)
            endIf
          endIf
  100   Continue
      Call Gauerr('Fermi level not found in DetFit.')
      Return
      End
*Deck DgFock
      Subroutine DgFock(IOut,IPrint,ISpin,NRI,NDimBl,ROHF,ZDO,RComp,
     $  IVShft,VSGoal,VSMax,VSGpMn,IDiagI,IAlDNI,NPassI,ISyOcc,IWfnSy,
     $  DoDiag,NBlks,NBasis,NBD,NOrbs,NE,IOF,IRwC,IOTMEx,IOTM,IndTM,IOC,
     $  IndIOC,IOEig,IndEig,IOP,AccDes,NBF,NBFU,NOccSy,CutOff,NameSp,
     $  JCycle,IPFON,UpdTmp,TempSt,ErrMax,DLaMax,NAtoms,C,LLim,IULim,
     $  IOpDME,SECut,DThr,ZThr,PFThr,NEqBas,NOp1,IRWPO,IRWFO,IRwS,IRwV,
     $  Temp,Beta,EF,VShift,IRWFM,IMemOp,IFxFm,IOpCEM,F,D,V,Eig,LScr,
     $  Scr,X,LenX,SpDiag,NEBr,NOBr,Gap,OVLim,OVMax,IOSymA,IndSyA,
     $  ICampK,DgAll,ISwCut,IPsCut,IFrgTp,MapF2M,Y,LenY)
      Implicit Real*8(A-H,O-Z)
C
C     Diagonalize the Fock matrix is IOF, after level-shifting if IVShft
C     is non-zero.  Store the resulting MO coefficients in IOC and the
C     density matrix in IOP.
C
C     IDiag ... 1 -- DiagD.
C               2 -- Krylov iterations.
C               3 -- Pseudo-diagonalization.
C               4 -- CGDMS.
C               5 -- PDM.
C               6 -- CEM.
C               7 -- Sign Matrix Method.
C               8 -- SNRDMS.
C               9 -- Converted to 3 before call to this routine.
C              10 -- Jacobi diagonalizaton.
C               Currrently 1-3,9-10 assume full matrix storage and 4-8 sparse.
C     ISyOcc ... 1 for no orbital symmetry, 2 to maintain occupations,
C                3 to maintain overall wavefunction symmetry (NYI).
C     IWfnSy ... Overall representation of wavefunction to maintain
C                if ISyOcc=3 (NYI).
C     NBlks  ... Number of irreps.
C     NE     ... Number of electrons for the current spin-set.
C     AccDes ... Convergence accuracy for Krylov if IDiag=2.
C     NOccSy ... Number of occupied orbitals of each irrep.
C     F      ... length 2*NTT.
C     D      ... length 2*NTT.  Density is returned here.
C     V      ... length 2*NTT, expect if Krylov is requested must be
C                Max(2*NTT,NBasis*(NE+8)) and if orbital symmetry is
C                constrained must be 2*NBasis**2 as well.
C                MOs are returned here.
C     Eig    ... Receives eigenvalues (NE returned if Krylov, otherwise
C                NBasis).
C     Scr    ... Scratch of length (NBasis,LScr), LScr >= 5.
C     ICampK ... <0 if Camp-King extrapolation was done.  In this case,
C                new orbitals are already in IRwC and we only need to
C                compute new densities here.
C
C     IOTM,IndTM File and offset with transformation to orthonormal basis.
C                IOTMEx flags whether this is a full matrix or blocked.
C
C     IOC,IndIOC File and offset which receive output orbitals.  Normally
C                the same as IOTM, except for fragment methods which set
C                a different space for each fragment each iteration.
C
      Logical ROHF, ZDO, RComp, DoDiag, DoFBrd, UpdTmp, OuterF, PrntDM,
     $  SpDiag, RegDag, DgAll, NonDag, UpdOcc, BlockV, IOTMEx
      Dimension NBF(*), NBFU(*), NOccSy(*), F(NRI,*), D(*), V(NBasis,*),
     $  X(*), Eig(*), Scr(NBD,LScr), ErrMax(0:*), NEqBas(NBasis,*),
     $  LLim(*), IULim(*), C(*), Y(*), MapF2M(*)
      Character*(*) NameSP, NamSP1*100
      Save Big, Zero, One, SmallG
      Data Big/1.0d3/, Zero/0.0d0/, One/1.0d0/, SmallG/0.001d0/
 1010 Format(I4)
 1020 Format(' NOrbs=',I6,' LenY=',I12,' NE2=',I12,' LScr=',I2)
 1030 Format(' Shifting Virtual Orbitals by:',f9.2)
 1040 Format(' Skip diagonalization as ',A,
     $  ' Fock matrix is already diagonal.')
C
      Call TStamp(1,'Top of DgFock.')
      UpdOcc = ISyOcc.le.1
      BlockV = IFrgTp.gt.0
      DgAll = .False.
      If(ISpin.le.1) then
        IndSy = 0
      else
        IndSy = NBlks
        endIf
      IOpClX = (NDimBl-1)*4 + (NRI-1)*2
      IDiag = IDiagI
      NBD = NBasis*NDimBl
      NTT = (NBD*(NBD+1))/2
      NBSq = NBasis*NBasis
      PrntDM = IPrint.ge.2.or.(JCycle.eq.1.and.IPSave(0).eq.0)
      RegDag = IDiag.eq.1.or.IDiag.eq.3.or.IDiag.eq.10
      IAlgDN = IAlDNI
      If(IAlgDN.eq.0.and.IDiag.eq.10) IAlgDN = 6
      NE2 = LenY/Max(NOrbs,1)
      If(ICampK.lt.0) then
C       Nothing to do if pickup up CK extrapolated orbitals.
      else if(SpDiag.and.IDiag.eq.1) then
        IStart = 1
        Call SpGet(X,X,LenX,IRWFO,IStart,jif,if,IAftrF,NBasis,NZF,1,
     $    .true.,0)
        Call Ucmprs(X(IAftrF),X(jif),X(if),NBasis)
        Call AMove(NTT,X(IAftrF),X(IStart))
        IVScr = IStart + NTT
        IScr1 = IVScr + NBSq
        IScr2 = IScr1 + NBasis
        IEig  = IScr2 + NBasis*10
        IEnd  = IEig + NBasis
        Call TstCor(IEnd,LenX,'DiagDN')
        Call DiagDN(IOut,IPrint,IAlgDN,X(IStart),X(IVScr),X(IEig),
     $    NBasis,X(IScr1),10,X(IScr2),NBasis,.False.)
        Call FormP2(.True.,1,NBasis,NE,0,One,X(IVScr),X(IStart))
        Call MkSprL(X,LenX,0,0,0,IStart,jip,ip,IAftrP,NBasis,NZP,
     $    'DiagDN',zthr,1)
      else if(IVShft.eq.0.and.ZDO.and.IDiag.eq.1) then
        Call FileIO(2,-IOF,NTT,F,0)
        If(IPrint.ge.2.or.NE2.lt.LScr)
     $    Write(IOut,1020) NOrbs, LenY, NE2, LScr
        If(NE2.gt.LScr) then
          Call DiagDN(IOut,IPrint,IAlgDN,F,V,Eig,NBasis,Scr,NE2,Y,
     $      NBasis,.False.)
        else
          Call DiagDN(IOut,IPrint,IAlgDN,F,V,Eig,NBasis,Scr,LScr-1,
     $      Scr(1,2),NBasis,.False.)
          endIf
      else if(RegDag.and.RComp) then
        If(ISpin.ne.0) Call GauErr('ISpin error in DgFock.')
        If(BlockV) Call GauErr('RComp and fragment in DgFock.')
        NonDag = JCycle.gt.1
        Call SftOpn(IOut,IPrint,ISpin,NRI,NDimBl,RComp,UpdOcc,CutOff,
     $    NBasis,NOrbs,IFrgTp,NBlks,NOccSy,IndSy,NBF,NBFU,MapF2M,IVShft,
     $    VSGoal,VSMax,VSGpMn,.True.,IOF,IOTMEx,IOTM,IndTM,D(NTT+1),D,V,
     $    Scr,VShift,Gap,NonDag)
        If(Gap.le.SmallG.and.IDiag.eq.3) IDiag = 1
        Call GtSEig(NRI,NBlks,NBFU,D(NTT+1),Eig)
        If(DoDiag.and.NonDag) then
          If(IDiag.eq.1.or.IDiag.eq.10) then
            Call DiagDS(IOut,IPrint,IAlgDN,BlockV,NRI,NBlks,NBFU,NOrbs,
     $        NOccSy,OVLim,OVMax,D(NTT+1),V,Eig,LScr,Scr,D)
            If(IPrint.ge.2) Call PrtSEig(IOut,NBlks,NBFU,NOccSy,Eig)
            Call FileIO(2,-IOTM,NRI*NBasis*NOrbs,D,IndTM)
            If(NRI.eq.2.or.NDimBl.eq.2)
     $        Call GauErr('Missing NRI=2 code for RComp in DgFock.')
            If(NBasis.gt.NOrbs)
     $        Call UnPcck(0,Junk,V,NBasis,NBasis,NOrbs,NOrbs)
            Call MatRec(D,V,Scr(1,2),NBasis,NBasis,NOrbs,NOrbs,4)
            DgAll = .True.
          else if(IDiag.eq.3) then
            Call GauErr('No Pseudo with RComp.')
c           Call PsuDgS(IOut,IPrint,NRI,NPassI,.False.,Zero,NBlks,NBFU,
c    $        NOccSy,NBasis,D(NTT+1),V,Eig,ANorm,OVMax)
          else
            Call GauErr('Inconsistent IDiag #1 in DgFock.')
            endIf
        else if(DoDiag) then
          Write(IOut,1040) NameSp
          endIf
      else if(RegDag) then
        NonDag = JCycle.gt.1
        Call SftOpn(IOut,IPrint,ISpin,NRI,NDimBl,RComp,UpdOcc,CutOff,
     $    NBasis,NOrbs,IFrgTp,NBlks,NOccSy,IndSy,NBF,NBFU,MapF2M,IVShft,
     $    VSGoal,VSMax,VSGpMn,.True.,IOF,IOTMEx,IOTM,IndTM,F,D,V,Scr,
     $    VShift,Gap,NonDag)
        If(Gap.le.SmallG.and.IDiag.eq.3) IDiag = 1
        Call GtSEig(NRI,NBlks,NBFU,F,Eig)
        If(DoDiag.and.NonDag) then
          If(IDiag.eq.1.or.IDiag.eq.10) then
            If(IPrint.ge.2.or.NE2.lt.LScr)
     $        Write(IOut,1020) NOrbs, LenY, NE2, LScr
            If(NE2.gt.LScr) then
              Call DiagDS(IOut,IPrint,IAlgDN,BlockV,NRI,NBlks,NBFU,
     $          NOrbs,NOccSy(IndSy+1),OVLim,OVMax,F,D,Eig,NE2,Y,F)
            else
              Call DiagDS(IOut,IPrint,IAlgDN,BlockV,NRI,NBlks,NBFU,
     $          NOrbs,NOccSy(IndSy+1),OVLim,OVMax,F,D,Eig,LScr,Scr,F)
              endIf
            If(IPrint.ge.2)
     $        Call PrtSEig(IOut,NBlks,NBFU,NOccSy(IndSy+1),Eig)
            If(BlockV) then
              Call BlkMpy(NRI,'N','N',NBlks,NBF,NBFU,NBFU,One,V,NBF,
     $          NBFU,D,NBFU,NBFU,Zero,F,NBF,NBFU)
              Call RecF2M(NRI,NBlks,NBF,NBFU,MapF2M,F,NBD,V)
            else
              Call XGEMM(NRI,'N','N',NBD,NOrbs,NOrbs,One,V,NBD,D,NOrbs,
     $          Zero,F,NBD)
              Call AMove(NRI*NBD*NOrbs,F,V)
              endIf
            DgAll = .True.
          else if(IDiag.eq.3) then
            Call PsuDgS(IOut,IPrint,IAlgDN,NRI,NPassI,ISwCut,IPsCut,
     $        .False.,Zero,NBlks,NBFU,NOccSy(IndSy+1),NBasis,F,V,Eig,D,
     $        LScr,Scr,OVMax,DgAll,Y,LenY)
          else
            Call GauErr('Inconsistent IDiag #2 in DgFock.')
            endIf
        else if(DoDiag) then
          Write(IOut,1040) NameSp
          endIf
      else if(IDiag.eq.2) then
        If(IVShft.ne.0) Call GauErr('VShift with KyDiag.')
        If(RComp) Call GauErr('RComp with KyDiag.')
        Call GauErr('KyDiag broken.')
c       Call FileIO(2,-IOV,NRI*NE*NBasis,V,IndTM)
        Call FileIO(2,-IOTM,NBSq,F,NBSq)
        Call MatMpy(NBasis,NBasis,NE,F,V,D)
        Call FileIO(2,-IOF,NRI*NTT,F(1,NBSq-NTT+1),0)
c       Call FileIO(2,-IOV,NBSq,V,0)
        Call VFVMul(.False.,.False.,.False.,NE,NBasis,NBasis,NBasis,
     $    F(1,NBSq-NTT+1),V,F,F,Scr,NBasis)
        Call KyDiag(IOut,IPrint,NBasis,NE,AccDes,F,D,Eig,V,V(1,NE+1),
     $    F(1,NTT+1),Scr,Scr(1,2))
c       Call FileIO(2,-IOV,NBSq,F,0)
        Call MatMpy(NBasis,NBasis,NE,F,D,V)
        Do 10 I = (NE+1), NBasis
   10     Eig(I) = Big
        Call AClear(NBasis*(NBasis-NE),V(1,NE+1))
      else if(IDiag.le.8) then
        IVShf1 = Abs(IVShft)
        OuterF = IDiag.ne.4.or.Mod((IOpDMe/10000),10).eq.0
        Temp = Zero
        ix = 1
        istart = ix + NBasis
        If(OuterF) then
          Cut = -Big
          If(JCycle.eq.1.and.SECut.gt.Zero.and.ZDO) Cut = SECut
          If(IPrint.ge.1) Write(IOut,1030) Gfloat(IVShf1)
          Call FxInit(IOut,IPrint,X,LenX,IRWPO,IRWFO,IRWFM,NBasis,
     $      ZThr,istart,jip,ip,if,NZ,IAftrP,iend,Cut,C,LLim,IULim,
     $      NAtoms,IMemOp,IVShf1,IFxFm,NE)
          istart = iend
          endIf
        If(IDiag.eq.4) then
C         CGDMS
          Call SCGDMS(IOut,IPrint,PrntDM,X,LenX,NBasis,ZThr,DThr,NE,
     $      DLaMax,IOpDMe,pfthr,JCycle,jip,ip,if,NZ,istart,iend,IMemOp)
        else if(IDiag.eq.5) then
C         Purification of the Density matrix
          Call SCnPDM(IOut,IPrint,PrntDM,X,LenX,X(ix),istart,jip,ip,if,
     $      NZ,iend,nbasis,JCycle,NE,IOpDMe,Tcg2a,npf,xmax,xmin,demp,
     $      pfthr,IMemOp,zthr)
        else if(IDiag.eq.6) then
C         Chebyshev Expansion Method
          Call SChbEM(IOut,IPrint,PrntDM,X,LenX,X(ix),istart,jip,ip,if,
     $      iend,nz,nbasis,JCycle,npfcy,NE,IOpCEM,IOpDMe,xmax,xmin,demp,
     $      pfthr,IMemOp,zthr)
        else if(IDiag.eq.7) then
C         Sign Matrix Method
          If(.not.OuterF) then
            Call SpGet(X,X,LENX,IRWPO,istart,jip,ip,iend,NBasis,NZP,1,
     $        .True.,0)
            istart = iend
            Call SpGet(X,X,LENX,IRWFO,istart,jif,if,iend,NBasis,NZF,1,
     $        .True.,0)
            Call Compct(X(jip),X(ip),zthr,NBasis,NZp,1,.True.)
            Call Compct(X(jif),X(if),zthr,NBasis,NZF,1,.True.)
          else
            NZP = NZ
            NZF = NZ
            endIf
          Call SignMt(IOut,IPrint,PrntDM,X,LenX,X(ix),nbasis,NE,jip,ip,
     $      nzp,jif,if,nzf,istart,iend,jcycle,IOpDMe,zthr,demp,IMemOp,
     $      OuterF)
          If(.not.OuterF) then
            ipa = ip
            nzpa = nzp
            IAftrP = ipa + NOrbs + NZPA
            endIf
        else if(IDiag.eq.8) then
          Call SNRDMS(IOut,IPrint,PrntDM,X,LenX,NBasis,ZThr,DThr,NE,
     $      IOpDMe,JCycle,jip,ip,if,NZ,istart,iend,IMemOp,ifail)
          If(IFail.ge.1) Call SCGDMS(IOut,IPrint,PrntDM,X,LenX,NBasis,
     $      ZThr,DThr,NE,DLaMax,IOpDMe,pfthr,JCycle,jip,ip,if,NZ,istart,
     $      iend,IMemOp)
        else
          Call GauErr('Unknown sparse diagonalization replacement.')
          endIf
        Call Compct(X(jip),X(ip),ZThr,NBasis,NZP,1,.True.)
        If(IPrint.ge.1) Call SpRprt(NZP,(NBasis*(NBasis-1))/2,zthr,
     $    'Ortho Density Matrix')
        Call SpStor(IRWPO,X(jip),X(ip),NBasis,1,.True.,0)
      else
        Call GauErr('Unrecognized diagonalization strategy.')
        endIf
      If(SpDiag) then
        Call POrtho(IOut,IPrint,X,LenX,IAftrP,jip,ip,IRwPO,IOP,IRWS,
     $    IRWV,NBasis,NAtoms,ZThr,NOp1,NEqBas,IMemOp,NE,ZDO,LLim,IULim)
      else if(RegDag.or.IDiag.eq.2) then
        If(ICampK.lt.0) then
          Call FileIO(2,-IRwC,NRI*NBD*NOrbs,V,0)
          DoFBrd = .False.
        else
          Call FixPha(NRI,NBD,NOrbs,0,V)
          Call FileIO(1,-IOC,NRI*NBD*NOrbs,V,IndIOC)
          If(LScr.lt.3) Call GauErr('LScr too small in DgFock.')
          Call UnPerS(IOut,IPrint,ROHF,ISyOcc,IWfnSy,NBlks,NBFU,
     $      NOccSy(IndSy+1),NRI*NBD,NOrbs,V,Eig,D,Scr,Scr(1,2),Scr,
     $      Scr(1,3))
          If(IOSymA.ne.0)
     $      Call FileIO(1,-IOSymA,InToWP(NOrbs),Scr(1,2),IndSyA)
          If(IOC.ne.IRwC) Call FileIO(1,-IRwC,NRI*NBD*NOrbs,V,0)
          Call FileIO(1,-IOEig,NOrbs,Eig,IndEig)
          DoFBrd = IPFON.gt.0
          endIf
        Call EFermi(DoFBrd,UpdTmp,NOrbs,NE,Eig,TempSt,ErrMax,EF,Temp,
     $    Beta)
        Call FrmPSA(IOut,IPrint,NameSp,.True.,NRI,IPFON,NBD,NOrbs,
     $    NE-NEBr,NEBr,NOBr,Temp,Beta,ErrMax,EF,Eig,V,D,Scr)
        Call FileIO(1,-IOP,NRI*NTT,D,0)
        If(IPrint.ge.3) then
          NamSp1 = ' '
          If(NDimBl.eq.1) then
            NamSp1 = NameSp
            Ind1 = LinEnd(NamSp1) + 2
          else
            Ind1 = 1
            endIf
          NamSp1(Ind1:) = 'MO coefficients at cycle'
          Ind = LinEnd(NamSp1) + 2
          Write(NamSp1(Ind:),1010) JCycle
          L = LinEnd(NamSp1)
          Call LTOutX(IOut,3,IOpClX,1,0,NamSp1(1:L),NBasis,NOrbs,V,V,0)
          NamSp1(Ind1:Ind-2) = 'density matrix  at cycle'
          Call LTOutX(IOut,0,IOpClX,1,0,NamSp1(1:L),NBasis,NOrbs,D,D,0)
          endIf
        endIf
      Call TStamp(1,'End of DgFock.')
      Return
      End
*Deck DMPEFG
      Subroutine DmpEfg(Ntot,xyz,Natoms,C,AtmChg,EFG)
      Implicit Real*8 (A-H,O-Z)
C
C     Nuclear potential, field, and field gradient contributions.
C
      Dimension xyz(Ntot,3),C(3,Natoms),AtmChg(Natoms),EFG(Ntot,10)
      Save Three, One
      Data Three/3.0d0/, One/1.0d0/
C
      Do 10 J=1,Natoms
        Do 10 I=1,Ntot
          distx = xyz(i,1)-C(1,J)
          disty = xyz(i,2)-C(2,J)
          distz = xyz(i,3)-C(3,J)
          disti = One/(sqrt(distx**2+disty**2+distz**2))
          d3 = Atmchg(j)*disti**3
          d5 = d3*disti**2
          efg(i,1) = efg(i,1) + atmchg(j)*disti
          efg(i,2) = efg(i,2) + distz*d3
          efg(i,3) = efg(i,3) + disty*d3
          efg(i,4) = efg(i,4) + distx*d3
          efg(i,5) = efg(i,5) + d3 - three*d5*distz**2
          efg(i,7) = efg(i,7) + d3 - three*d5*disty**2
          efg(i,10) = efg(i,10) + d3 - three*d5*distx**2
          efg(i,9) = efg(i,9) - three*d5*distx*disty
          efg(i,8) = efg(i,8) - three*d5*distx*distz
   10     efg(i,6) = efg(i,6) - three*d5*disty*distz
      Return
      End
*Deck dNkdRk
      Subroutine dnkdrk(IDo,IDid,D1Rho,D2Rho,dNdR,rhatx,rhaty,rhatz)
      Implicit Real*8 (A-H,O-Z)
      Dimension D1Rho(IDo,3),D2Rho(IDo,6),dNdR(IDo,3),rhatx(*),
     $  rhaty(*), rhatz(*)
      Save One
      Data One/1.d0/
C
c     Computes dnkdrk.
c
      Do 10 I = 1, IDo
        IIDid = I + IDid
        KDid  = IIDid
        GMag = Sqrt(D1Rho(I,1)**2+D1Rho(I,2)**2+D1Rho(I,3)**2)
        GMagi = One/GMag
        GMagi3 = GMagi**3
        part1x = GMagi*(D2Rho(I,1)*rhatx(KDid)+
     $    D2Rho(I,4)*rhaty(KDid)+D2Rho(I,5)*rhatz(KDid))
        part1y = GMagi*(D2Rho(I,2)*rhaty(KDid)+
     $    D2Rho(I,4)*rhatx(KDid)+D2Rho(I,6)*rhatz(KDid))
        part1z = GMagi*(D2Rho(I,3)*rhatz(KDid)+
     $    D2Rho(I,5)*rhatx(KDid)+D2Rho(I,6)*rhaty(KDid))
        scalar=D2Rho(I,1)*rhatx(KDid)*D1Rho(I,1)+
     $         D2Rho(I,2)*rhaty(KDid)*D1Rho(I,2)+
     $         D2Rho(I,3)*rhatz(KDid)*D1Rho(I,3)+
     $         D2Rho(I,4)*rhatx(KDid)*D1Rho(I,2)+
     $         D2Rho(I,4)*rhaty(KDid)*D1Rho(I,1)+
     $         D2Rho(I,5)*rhatx(KDid)*D1Rho(I,3)+
     $         D2Rho(I,5)*rhatz(KDid)*D1Rho(I,1)+
     $         D2Rho(I,6)*rhaty(KDid)*D1Rho(I,3)+
     $         D2Rho(I,6)*rhatz(KDid)*D1Rho(I,2)
        part2x = GMagi3*scalar*D1Rho(I,1)
        part2y = GMagi3*scalar*D1Rho(I,2)
        part2z = GMagi3*scalar*D1Rho(I,3)
        dNdR(I,1) = -part1x+part2x
        dNdR(I,2) = -part1y+part2y
        dNdR(I,3) = -part1z+part2z
   10   Continue
      Return
      End
*Deck DoRhGr
      Subroutine DoRhGr(IOut,IPrint,IOpCl,IOpt,MOBas,NTot,NDimX,NBasis,
     $  NMO,CutOff,XYZ,P,CMO,Rho,NDimD1,D1Rho,NDimD2,D2Rho,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Evaluate the density and its gradient, possible working in chunks
C     to save space.
C
      Parameter (ISafe=100000)
      Logical MOBas
      Dimension XYZ(NDimX,3), P(*), CMO(NBasis,NMO), Rho(NTot),
     $  D1Rho(NDimD1,3), D2Rho(NTot,6), V(MDV)
 1000 Format(' MDV=',I12,' ISafe=',I8,' Len1=',I8,' NBasis=',I5,' L12=',
     $  2I2,' NTot=',I6)
C
      If(NTot.eq.0) Return
      NTT = (NBasis*(NBasis+1))/2
      L1 = 0
      If(IOpt.ge.2) L1 = 1
      L2 = 0
      If(IOpt.ge.3) L2 = 1
      Len1 = 1 + 3*L1 + 6*L2
      MaxPt = Min((MDV-ISafe)/(Len1*NBasis+3*L1+6*L2),NTot)
      If(MaxPt.lt.1) then
        Write(IOut,1000) MDV, ISafe, Len1, NBasis, L1, L2, NTot
        Call GauErr('Out of space in DoRhGr.')
        endIf
      ID1Rho = 1
      ID2Rho = ID1Rho + L1*3*MaxPt
      IChi = ID2Rho + L2*6*MaxPt
      ID1Chi = IChi + MaxPt*NBasis
      ID2Chi = ID1Chi + L1*3*MaxPt*NBasis
      IV = ID2Chi + L2*6*MaxPt*NBasis
      MDV1 = MDV - IV + 1
      Do 10 IPt = 1, NTot, MaxPt
        NDo = Min(NTot-IPt+1,MaxPt)
        Call RhGrRh(IPrint,NDo,NDimX,NTT,IOpt,NBasis,XYZ(IPt,1),Rho,
     $    V(ID1Rho),V(ID2Rho),V(IChi),V(ID1Chi),V(ID2Chi),P,CutOff,NMO,
     $    IOpCl,CMO,MOBas,V(IV),MDV1)
        If(IOpt.ge.2)
     $    Call AExpnd(NDo,NDo,NDimD1,3,V(ID1Rho),D1Rho(IPt,1))
        If(IOpt.ge.3)
     $    Call AExpnd(NDo,NDo,NDimD2,6,V(ID2Rho),D2Rho(IPt,1))
   10   Continue
      Return
      End
*Deck dWidR1
      Subroutine dwidr1(NTot,XYZ,Gnorm,weight,dwidrn,ICen,NOrg,COrg,
     $  ISurf,IBWFnc,Inbck,PBeck,DPBeck,rhatx,rhaty,rhatz,DMag)
      Implicit Real*8 (A-H,O-Z)
      Parameter (Expent=0.8d0)
      Dimension XYZ(NTot,3),Gnorm(NTot,3),Weight(NTot),dwidrn(NTot,2),
     $  Icen(NTot),Corg(3,NOrg),PBeck(Inbck,NOrg),dPbeck(Inbck,NOrg),
     $  rhatx(*),rhaty(*),rhatz(*),DMag(*)
      Real*8 MDCutO
      Save Zero,Pt5,One,OnePt5,Two
      Data Zero/0.d0/,Pt5/0.5d0/,One/1.d0/,OnePt5/1.5d0/,Two/2.d0/
C
C     Rewrite of dwidri to compute quantities for derivs of
C     weights in large space.
C
C     dwi/dri = dwidrn(I,1) + dwidrn(I,2)*(rhat(I) . dndr(I))   (scalar)
C
C     we don't have dndr here yet.
C
C     dwi/dni = dwidrn(I,2) * rhat(I)                           (vector)
C
c     Deriv of weights wrt r, and n.
c
      VSmall = MDCutO(0)
      Do 10 I = 1, NTot
        rdotgr = One/(rhatx(i)*gnorm(I,1)+
     $    rhaty(i)*gnorm(I,2)+rhatz(i)*gnorm(I,3))
        dwidrn(I,1) = Two*Weight(I)*DMag(I)
        dwidrn(I,2) = -rdotgr*Weight(I)
   10   Continue
c
c     That's all if single origin surface. Otherwise include
c     contributions from Becke terms.
c
      If(ISurf.gt.2) then
        If(IBWFnc.eq.1) Then
          Do 70 I = 1, NOrg
            Do 20 K = 1, NTot
              Pbeck(K,I) = One
   20         dPbeck(K,I) = Zero
            Do 40 J = 1, Norg
              If(J.ne.I) Then
                Distij = Sqrt((corg(1,i)-corg(1,j))**2+
     $                        (corg(2,i)-corg(2,j))**2+
     $                        (corg(3,i)-corg(3,j))**2)
                Do 30 K = 1, NTot
                  Disti = Sqrt((XYZ(K,1)-COrg(1,I))**2+
     $                         (XYZ(K,2)-COrg(2,I))**2+
     $                         (XYZ(K,3)-COrg(3,I))**2)
                  Distj = Sqrt((XYZ(K,1)-COrg(1,J))**2+
     $                         (XYZ(K,2)-COrg(2,J))**2+
     $                         (XYZ(K,3)-COrg(3,J))**2)
                  dmu=(disti-distj)/distij
                  f1=onept5*dmu-pt5*(dmu**3)
                  f2=onept5*f1-pt5*(f1**3)
                  f3=onept5*f2-pt5*(f2**3)
                  ess=Pt5*(One-f3)
   30             Pbeck(K,I)=Pbeck(K,I)*ess
                EndIf
   40         Continue
            Do 60 J = 1, NOrg
              If(J.ne.I) Then
                Distij = Sqrt((corg(1,i)-corg(1,j))**2+
     $                        (corg(2,i)-corg(2,j))**2+
     $                        (corg(3,i)-corg(3,j))**2)
                Do 50 K = 1, NTot
                  If(Abs(PBeck(K,I)).gt.VSmall) then
                    Centx = XYZ(K,1) - COrg(1,ICen(K))
                    Centy = XYZ(K,2) - COrg(2,ICen(K))
                    Centz = XYZ(K,3) - COrg(3,ICen(K))
                    distix = XYZ(K,1)-COrg(1,I)
                    distiy = XYZ(K,2)-COrg(2,I)
                    distiz = XYZ(K,3)-COrg(3,I)
                    distjx = XYZ(K,1)-COrg(1,J)
                    distjy = XYZ(K,2)-COrg(2,J)
                    distjz = XYZ(K,3)-COrg(3,J)
                    Disti = Sqrt(distix**2+distiy**2+distiz**2)
                    Distj = Sqrt(distjx**2+distjy**2+distjz**2)
                    dmu = (disti-distj)/distij
                    f1 = onept5*dmu-pt5*(dmu**3)
                    f2 = onept5*f1-pt5*(f1**3)
                    f3 = onept5*f2-pt5*(f2**3)
                    ess = Pt5*(One-f3)
                    discen = One/(sqrt(centx**2+centy**2+centz**2))
                    doti = (distix*centx+distiy*centy+
     $                      distiz*centz)*discen
                    dotj = (distjx*centx+distjy*centy+
     $                      distjz*centz)*discen
                    dmudrk = (doti/Disti-dotj/Distj)/Distij
                    f1der = (onept5*dmudrk-onept5*dmu**2*dmudrk)
                    f2der = onept5*f1der-onept5*f1**2*f1der
                    f3der = onept5*f2der-onept5*f2**2*f2der
                    essder = -pt5*f3der
                    dPbeck(K,I) = dpbeck(K,I)+essder*Pbeck(K,I)/ess
                    endIf
   50             Continue
                EndIf
   60         Continue
   70       Continue
        Else If(IBWFnc.eq.2) then
          Do 90 I=1,NOrg
            Do 80 K=1,NTot
              DistX=XYZ(K,1)-COrg(1,I)
              DistY=XYZ(K,2)-COrg(2,I)
              DistZ=XYZ(K,3)-COrg(3,I)
              Dist=Sqrt(DistX**2+DistY**2+DistZ**2)
              Dist2=Dist**2
              PBeck(K,I)=Exp(-Expent*Dist2)
              Dot=DistX*RHatX(K)+DistY*RHatY(K)+
     $            DistZ*RHatZ(K)
   80         dPBeck(K,I)=-Two*Expent*Dot*PBeck(K,I)
   90       Continue
        Else
          Call GauErr('dwidr1: Invalid function for Becke weights.')
          EndIf
C
        Do 110 K = 1, NTot
          Denom1 = Zero
          Denum = Zero
          Do 100 I = 1, NOrg
            Denum = Denum + Dpbeck(K,I)
  100       Denom1 = Denom1 + Pbeck(K,I)
          Denom1 = One/Denom1
  110     dwidrn(k,1) = dwidrn(k,1) + Weight(k)*
     $      (dPbeck(K,ICen(k))/PBeck(k,ICen(k))-Denum*Denom1)
      EndIf
      Return
      End
*Deck dwidr2
      Subroutine dwidr2(NDo,IDid,NTot,dndr,dwidrn,rhatx,rhaty,
     $  rhatz,dwdr,dwdn)
      Implicit Real*8(A-H,O-Z)
      Dimension dndr(NDo,3), dwidrn(NTot,2), rhatx(*), rhaty(*),
     $  rhatz(*), dwdr(NDo), dwdn(NDo,3)
C
C     This routine used to do a lot more. Now it just computes
C     the remaining terms in the weight derivatives.
C
      Do 20 I = 1, NDo
        KDid = I+IDid
        dwdr(I) = dwidrn(KDid,1) + dwidrn(KDid,2)*
     $    (rhatx(KDid)*dndr(I,1) + rhaty(KDid)*dndr(I,2) +
     $       rhatz(KDid)*dndr(I,3))
        dwdn(I,1) = dwidrn(KDid,2) * rhatx(KDid)
        dwdn(I,2) = dwidrn(KDid,2) * rhaty(KDid)
        dwdn(I,3) = dwidrn(KDid,2) * rhatz(KDid)
   20   Continue
      Return
      End
*Deck dWidRi
      Subroutine dwidri(NTot,IDid,IDo,XYZ,Gnorm,weight,dndr,dwdr,Cent,
     $  ICen,NOrg,COrg,ISurf,Inbck,PBeck,DPBeck,rhatx,rhaty,rhatz,DMag)
      Implicit Real*8 (A-H,O-Z)
      Logical OldFnc
      Parameter (OldFnc=.True.,Expent=0.8d0)
      Dimension XYZ(NTot,3),Gnorm(NTot,3),Weight(NTot),dndr(IDo,3),
     $  dwdr(IDo),Cent(NTot,3),Icen(NTot),Corg(3,NOrg),
     $  PBeck(Inbck,NOrg),dPbeck(Inbck,NOrg),rhatx(*),rhaty(*),rhatz(*),
     $  DMag(*)
      Real*8 MDCutO
      Save Zero,Pt5,One,OnePt5,Two
      Data Zero/0.d0/,Pt5/0.5d0/,One/1.d0/,OnePt5/1.5d0/,Two/2.d0/
C
c     Deriv of weights wrt r
c
      VSmall = MDCutO(0)
      Do 10 I = 1, IDo
        rdotgr = One/(rhatx(i)*gnorm(I+IDid,1)+
     $    rhaty(i)*gnorm(I+IDid,2)+rhatz(i)*gnorm(I+IDid,3))
        dwdrs = Two*DMag(i) - rdotgr*(
     $    rhatx(i)*dndr(I,1)+rhaty(i)*dndr(I,2)+rhatz(i)*dndr(I,3))
   10   dwdr(I) = weight(i+IDid)*dwdrs
c
c     That's all if single origin surface. Otherwise include
c     contributions from Becke terms.
c
      If(ISurf.gt.2) then
        If(OldFnc) Then
          Do 20 I = 1, NOrg
            Do 30 K = 1, IDo
   30         Pbeck(K,I) = One
            Do 40 J = 1, Norg
              If(J.ne.I) Then
                Distij = Sqrt((corg(1,i)-corg(1,j))**2+
     $                        (corg(2,i)-corg(2,j))**2+
     $                        (corg(3,i)-corg(3,j))**2)
                Do 50 K = 1, IDo
                  Disti = Sqrt((XYZ(K+IDid,1)-COrg(1,I))**2+
     $                         (XYZ(K+IDid,2)-COrg(2,I))**2+
     $                         (XYZ(K+IDid,3)-COrg(3,I))**2)
                  Distj = Sqrt((XYZ(K+IDid,1)-COrg(1,J))**2+
     $                         (XYZ(K+IDid,2)-COrg(2,J))**2+
     $                         (XYZ(K+IDid,3)-COrg(3,J))**2)
                  dmu=(disti-distj)/distij
                    f1=onept5*dmu-pt5*(dmu**3)
                    f2=onept5*f1-pt5*(f1**3)
                    f3=onept5*f2-pt5*(f2**3)
                    ess=Pt5*(One-f3)
   50             Pbeck(K,I)=Pbeck(K,I)*ess
                EndIf
   40         Continue
   20       Continue
          Do 60 I = 1, NOrg
            Do 70 K = 1, IDo
   70         dPbeck(K,I) = Zero
            Do 80 J = 1, NOrg
              If(J.ne.I) Then
                Distij = Sqrt((COrg(1,i)-COrg(1,j))**2+
     $                        (COrg(2,i)-COrg(2,j))**2+
     $                        (COrg(3,i)-COrg(3,j))**2)
                Do 90 K = 1, IDo
                  If(Abs(PBeck(K,I)).gt.vsmall) then
                    Disti = Sqrt((XYZ(K+IDid,1)-COrg(1,I))**2+
     $                           (XYZ(K+IDid,2)-COrg(2,I))**2+
     $                           (XYZ(K+IDid,3)-COrg(3,I))**2)
                    Distj = Sqrt((XYZ(K+IDid,1)-COrg(1,J))**2+
     $                           (XYZ(K+IDid,2)-COrg(2,J))**2+
     $                           (XYZ(K+IDid,3)-COrg(3,J))**2)
                    dmu = (disti-distj)/distij
                      f1 = onept5*dmu-pt5*(dmu**3)
                      f2 = onept5*f1-pt5*(f1**3)
                      f3 = onept5*f2-pt5*(f2**3)
                      ess = Pt5*(One-f3)
                    distix = XYZ(k+IDid,1)-COrg(1,i)
                    distiy = XYZ(k+IDid,2)-COrg(2,i)
                    distiz = XYZ(k+IDid,3)-COrg(3,i)
                    distjx = XYZ(k+IDid,1)-COrg(1,j)
                    distjy = XYZ(k+IDid,2)-COrg(2,j)
                    distjz = XYZ(k+IDid,3)-COrg(3,j)
                    discen = One/(sqrt(cent(k+IDid,1)**2+
     $                                 cent(k+IDid,2)**2+
     $                                 cent(k+IDid,3)**2))
                    doti = (distix*cent(k+IDid,1)+distiy*cent(k+IDid,2)+
     $                      distiz*cent(k+IDid,3))*discen
                    dotj = (distjx*cent(k+IDid,1)+distjy*cent(k+IDid,2)+
     $                      distjz*cent(k+IDid,3))*discen
                    dmudrk = (doti/Disti-dotj/Distj)/Distij
                    f1der = (onept5*dmudrk-onept5*dmu**2*dmudrk)
                    f2der = onept5*f1der-onept5*f1**2*f1der
                    f3der = onept5*f2der-onept5*f2**2*f2der
                    essder = -pt5*f3der
                    dPbeck(K,I) = dpbeck(K,I)+essder*Pbeck(K,I)/ess
                  endIf
   90           Continue
              EndIf
   80       Continue
   60     Continue
        Else
          Do 100 I=1,NOrg
            Do 110 K=1,IDo
              DistX=XYZ(K+IDid,1)-COrg(1,I)
              DistY=XYZ(K+IDid,2)-COrg(2,I)
              DistZ=XYZ(K+IDid,3)-COrg(3,I)
              Dist=Sqrt(DistX**2+DistY**2+DistZ**2)
              Dist2=Dist**2
              PBeck(K,I)=Exp(-Expent*Dist2)
              Dot=DistX*RHatX(K+IDid)+DistY*RHatY(K+IDid)+
     $            DistZ*RHatZ(K+IDid)
  110         dPBeck(K,I)=-Two*Expent*Dot*PBeck(K,I)
  100     Continue
        EndIf
C
        Do 120 K = 1, IDo
          Denom1 = Zero
          Denum = Zero
          Do 130 I = 1, NOrg
            Denum = Denum + Dpbeck(K,I)
  130       Denom1 = Denom1 + Pbeck(K,I)
          Denom1 = One/Denom1
  120     dwdr(k) = dwdr(k) + Weight(k+IDid)*
     $      (dPbeck(K,ICen(k+IDid))/PBeck(k,ICen(k+IDid))-Denum*Denom1)
      EndIf
C
      Return
      End
*Deck EnerNR
      Function EnerNR(C,FM,RT,X,B,N,NOcc,NZ,NumC,NIts,aMuAcc,iprint,
     $  iout,NP,XX,LenRT,StoreT,IRWRT,TT,fail,IFunc)
      Implicit Real*8(A-H,O-Z)
C
C     Find the Mu which minimizes the Trace(P)-Nelecs using
C           Newton-Raphson with bisection
C
      Dimension C(NUMC),RT(LenRT,NumC),FM(N,2),X(NumC,4),XX(N,NP),
     $  TT(*)
C      IH(N+1,*)
      Character*13 Chart
      Logical StoreT,fail,bisect
      Save Zero, pt5, pres, pt99, Two, High
      Data Zero,pt5,pres,pt99,Two,High/0.d0,0.5d0,1.d-16,0.99d0,2.d0,
     $ 10000.d0/
 1315 Format(' CEM Density Deviation from Idempotency:  ',
     $  1pd20.6)
 1000 Format(' ',a,' Tr(P)-N= ',1PD9.2,'  da= ',1PD9.2,'  a= ',1PD9.2)
 2000 Format(' Final Trace:',1PD12.2,' Root:',1PD12.2)
C
C     Values of Mu should always be between -0.99 and 0.99 since F is scaled
C        and shifted
C
      fail=.false.
      aNmLow = ChbFnc(C,FM,RT,X,XX,-pt99,B,N,NOcc,NumC,NP,LenRT,IFunc)
      aNmHi = ChbFnc(C,FM,RT,X,XX,pt99,B,N,NOcc,NumC,NP,LenRT,IFunc)
      If(IPrint.ge.2) Write(IOut,*)'bracketing values',aNmLow,aNmHi
      If((aNmLow.gt.Zero.and.aNmHi.gt.Zero).or.(aNmLow.lt.Zero.and.
     $  aNmHi.lt.Zero)) then
        Write(IOut,*) 'root not bracketed',aNmLow,aNmHi
        fail = .true.
        EnerNR = High
        Return
        endIf
      If(aNmLow.eq.Zero) then
        aMuNxt = -pt99
        goto 999
      else if(aNmHi.eq.Zero) then
        aMuNxt = pt99
        goto 999
      else if(aNmLow.lt.Zero) then
        aMuLw = -pt99
        aMuHi = pt99
      else
        aMuHi = -pt99
        aMuLw = pt99
        endIf
      aMuNxt = Zero
      DMuOld = Two*pt99
      aMuOld = aMuLw
      DMu = DMuOld
      aNmOld = aNmLow
      aNm = ChbFnc(C,FM,RT,X,XX,aMuNxt,B,N,NOcc,NumC,NP,LenRT,IFunc)
      DNm = (aNm-aNmOld)/(aMuNxt-aMuOld)
      aNmOld = aNm
      If(IPrint.ge.2) Write(IOut,1000) '1st bisection',aNm,DMu,aMuNxt
      if(aNm.lt.zero) then
        aMuLw = aMuNxt
      else
        aMuHi = aMuNxt
        endIf
      bisect = .true.
      ibisct = 0
      do 300 i = 1, NIts
        if(Abs(DMuOld).lt.Abs(DMu).or.Abs(DMu).gt.4.d0.or.
     $    Abs(aNm).gt.1.d0.or.Abs(aNm/DNm).gt.3.d0.or.
     $    (bisect.and.ibisct.le.4)) then
          bisect = .true.
          ibisct = ibisct + 1
          DMuOld = DMu
          DMu = pt5*(aMuHi-aMuLw)
          aMuOld = aMuNxt
          aMuNxt = aMuLw + DMu
          write(chart,'(a)')'bi -  section'
          if(aMuLw.eq.aMuNxt) goto 999
        else
          ibisct = 0
          bisect = .false.
          DMuOld = DMu
          DMu = aNm/DNm
          aMuOld = aMuNxt
          aMuNxt = aMuNxt - DMu
          write(chart,'(a)')'Newton-Rapson'
          if(aMuOld.eq.aMuNxt) goto 999
          endIf
        aNm = ChbFnc(C,FM,RT,X,XX,aMuNxt,B,N,NOcc,NumC,NP,LenRT,IFunc)
        If(IPrint.ge.2) write(IOut,1000) chart,aNm,DMu,aMuNxt
        If((Abs(aNm).lt.aMuAcc.and.Abs(DMu).lt.aMuAcc).or.Abs(DNm-
     $    aNmOld).lt.pres.or.Abs(aMuNxt-aMuOld).lt.pres) goto 999
        DNm = (aNm-aNmOld)/(aMuNxt-aMuOld)
        If(Abs(DNm).lt.pres) goto 999
        aNmOld = aNm
        if(aNm.lt.zero) then
          aMuLw = aMuNxt
        else
          aMuHi = aMuNxt
          endIf
  300   Continue
      Write(IOut,*) 'EnerNR exceeded maximum iterations'
      fail = .true.
      EnerNR = High
      Return
C
C     Calculate idempotency
C
  999 aNm = ChbFnf(C,FM,RT,X,aMuNxt,B,N,NOcc,NZ,NumC,LenRT,StoreT,IRWRT,
     $  TT,IFunc)
      If(IPrint.ge.1) write(IOut,2000) aNm,aMuNxt
      dt2=(SpTrac(FM,FM,N,NZ)-arrsum(N,FM))/gfloat(N)
      If(IPrint.ge.1) write(IOut,1315) dt2
      EnerNR = Abs(dt2)
      Return
      End
*Deck FixHCF
      Subroutine FixHCF(IOut,IPrint,IRwH,IRwHU,IOpCl,ISCRF,EField,
     $  IPFlag,AllowP,FMM,FMFlag,FMFlg1,NFxFlg,IHMeth,NAtoms,IAn,IAtTyp,
     $  AtmChg,C,NBasis,NBas6D,V,MDV)
      Implicit Real*8 (A-H,O-Z)
C
C     Fix the core Hamiltonian to remove the PCM cavity-field interaction
C     in the presence of an external electric field.
C
      Logical FMM,AllowP(*)
      Integer FMFlag,FMFlg1,IAn(*),IAtTyp(*),JJ(1)
      Real*8 EField(*),AtmChg(*),C(3,*),V(*),XX(1)
      Save Zero,One,XX,JJ
      Data Zero/0.0d0/,One/1.0d0/,XX/0.0d0/,JJ/0/
C
C     Make a copy the RWF.
C
      LRwH = ItQry(IRwH)
      IRwHU = IGetRW(LRwH,0)
      Call FCopyP(IRwH,0,IRwHU,0,1,LRwH,V,MDV)
C
C     Allocate a little memory.
C
      NTT = (NBasis*(NBasis+1))/2
      NTT6D = (NBas6D*(NBas6D+1))/2
      jDME = 1
      jHA  = jDME + 3*NTT6D
      jHB  = jHA  + NTT
      jEnd = jHB  + NTT - 1
      Call TstCor(jEnd,MDV,'FixHCF')
      jV = jEnd + 1
      MDV1  = MDV - jV + 1
C
C     Compute the cavity-field iteraction.
C
      Call DipInt(IOut,IPrint,11014,0,1,1,1,3,EField,ISCRF,0,XX,IPFlag,
     $  AllowP,FMM,FMFlag,FMFlg1,NFxFlg,IHMeth,NAtoms,IAn,.False.,
     $  IAtTyp,0,AtmChg,C,3,NBasis,JJ,XX,Zero,0,1,JJ,XX,.False.,JJ,XX,
     $  -One,V(jDME),V(jV),MDV1)
C
C     Subtract from the core Hamiltonian.
C
      Call FileIO(2,-IRwHU,NTT,V(jHA),0)
      Call FileIO(2, IRwHU,NTT,V(jHB),0)
      Call AAdd(NTT,V(jHA),V(jDME),V(jHA))
      Call AAdd(NTT,V(jHB),V(jDME),V(jHB))
      Call FileIO(1,-IRwHU,NTT,V(jHA),0)
      Call FileIO(1, IRwHU,NTT,V(jHB),0)
      IOpClX = IOpCl - Mod(IOpCl,2)
      If(IPrint.gt.2) Call PrtFok(IOut,IOpClX,
     $  'Core hamiltonian after FixHCF',0,1,1,NBasis,V(jHA),V(jHB),
     $  V(jV),MDV1)
      Return
      End
*Deck FormBM
      Subroutine FormBM(NSaved,NBasis,BMat,MaxItr,IOErr,A,NAMax)
      Implicit Real*8(A-H,O-Z)
      Dimension BMat(MaxItr,MaxItr),A(NAMax)
      Save Zero,One
      Data Zero/0.D0/,One/1.0D0/
C
C     Do some initilization
C
      NTT = NBasis*(NBasis+1)/2
      Call AClear(MaxItr*MaxItr,BMat)
      IScr1 = 1
      IScr2 = IScr1 + NTT
      Do 10 I = 0, NSaved
        Do 10 J = 0, I
          If(I.eq.0) then
            BMat(1,1) = Zero
          else if(J.eq.0) then
            BMat(I+1, J+1) = -One
            BMat(J+1, I+1) = -One
          else
            Call FileIO(2,-IOErr,NTT,A(IScr1),(I-1)*NTT)
            Call FileIO(2,-IOErr,NTT,A(IScr2),(J-1)*NTT)
            BMat(I+1, J+1) = SCFTrc(A(IScr1),A(IScr2),NBasis,1)
            BMat(J+1, I+1) = BMat(I+1, J+1)
            endIf
 10     Continue
      Return
      End
*Deck FOrthx
      Subroutine FOrthx(IOut,IPrint,X,LenX,istart,IRWFA,IRwFB,IRWFOA,
     $  IRwFOB,IRwPA,IRwPB,IRWPOA,IRwPOB,IRwV,IRWH,IRWFM1,IRWFM2,NBasis,
     $  ZThr,IMemOp,IFxFm,NAE,NBE,ZDO,JCycle)
      Implicit real*8(a-h,o-z)
C
C     Retrieve and transform the Fock matrix
C
      Dimension X(LenX)
      Logical ZDO
C
      If(ZDO) then
        Call SpGet(X,X,LenX,IRWH,istart,jih,ih,iend,NBasis,NZH,0,
     $    .True.,0)
        If(IRWFM1.eq.0) IRWFM1 = IGetRW(1,0)
        Call SpStor(IRWFM1,X(jih),X(ih),NBasis,0,.True.,0)
        Call SpGet(X,X,LenX,IRWPA,istart,jip,ip,iend,NBasis,NZP,1,
     $    .True.,0)
        Call SpStor(IRWPOA,X(jip),X(ip),NBasis,1,.True.,0)
        Call SpGet(X,X,LenX,IRWFA,istart,jif,if,iend,NBasis,NZF,1,
     $    .True.,0)
        Call SpStor(IRWFOA,X(jif),X(if),NBasis,1,.True.,0)
        If(IRwPB.ne.0) Then
          Call SpGet(X,X,LenX,IRWPB,istart,jip,ip,iend,NBasis,NZP,1,
     $      .True.,0)
          Call SpStor(IRWPOB,X(jip),X(ip),NBasis,1,.True.,0)
          Call SpGet(X,X,LenX,IRWFB,istart,jif,if,iend,NBasis,NZF,1,
     $      .True.,0)
          Call SpStor(IRWFOB,X(jif),X(if),NBasis,1,.True.,0)
          endif
      else
        IRWFMA = 0
        Call FOrtho(IOut,IPrint,X,LenX,istart,IRWFA,IRWFOA,IRWPOA,IRwV,
     $    IRWFMA,NBasis,ZThr,IMemOp,IFxFm,NAE)
        If(IRwPB.ne.0) Then
          IRwFMB = 0
          Call FOrtho(IOut,IPrint,X,LenX,istart,IRWFB,IRWFOB,IRWPOB,
     $      IRwV,IRWFMB,NBasis,ZThr,IMemOp,IFxFm,NBE)
          Call MkSmFm(X,LenX,IRwFMA,0,0,IRwFMB,0,0,istart,istart,
     $      IRwFMA,jifm,ifm,ifm2,NZFM,iend,NBasis,0)
          endif
        If(JCycle.le.1) then
          IRWFM1 = IRWFMA
        else
          IRWFM2 = IRWFMA
          endIf
        endIf
      Return
      End
*Deck FOSimF
      Subroutine FOSimF(IOut,IPrint,RHF,IBas,NBasis,IRwSim,IRwV,IRwPA,
     $  IRwPB,IRwPOA,IRwPOB,IRwFA,IRwFB,IRwFRA,IRwFRB,Scr,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Routines to prepare SimOpt wave function error vectors
C
      Logical RHF
      Dimension V(MDV),Scr(*)
      Save Zero, Half, Two
      Data Zero/0.D0/,Half/0.5D0/, Two/2.0D0/
C
      Call FileIO(2,-IRwSim,1,Temp,0)
      Call BNUPak(Temp,NMax,MaxGDI)
      Call FileIO(2,-IRwSim,1,Temp,3)
      Call BNUPak(Temp,NPDIIS,NGDIIS)
      Call RwFGED(1,IRwSim,NMax,MaxGDI,IEn,INNEn,ICoeff,IQNDIS,IEnDIS,
     $  IMask,ICnvQM)
      NPDIIS = NPDIIS + 1
      If(NPDIIS.gt.MaxGDI) NPDIIS = MaxGDI
      Call BNPack(NPDIIS,NGDIIS,Temp)
      Call FileIO(1,-IRwSim,1,Temp,3)
      NTT = NBasis*(NBasis+1)/2
      If(RHF) then
        Call FileIO(2,-IRwPA,NTT,Scr,0)
        Call AScale(NTT,Half,Scr,Scr)
        Call FileIO(1,-IRwPA,NTT,Scr,0)
        endIf
      IDIS   = 1
      IEDIS  = IDIS  + NPDIIS**2
      IStart = IEDIS + NPDIIS**2
      Call AClear(NPDIIS**2,V(IDIS))
      Call AClear(NPDIIS**2,V(IEDIS))
      RMSdP = Zero
      Call SimVec(IOut,IPrint,IBas,NBasis,IRwV,IRwPA,IRwPOA,IRwFA,
     $  IRwFRA,MaxGDI,NPDIIS,RMSdP,V(IDIS),V(IEDIS),V(IStart),
     $  MDV-IStart+1)
      If(.not.RHF) Call SimVec(IOut,IPrint,IBas,NBasis,IRwV,IRwPB,
     $  IRwPOB,IRwFB,IRwFRB,MaxGDI,NPDIIS,RMSdP,V(IDIS),V(IEDIS),
     $  V(IStart),MDV-IStart+1)
      If(RHF) then
        RMSdP = Two*RMSdP
        Call AScale(NPDIIS**2,Two,V(IDIS),V(IDIS))
        Call AScale(NPDIIS**2,Two,V(IEDIS),V(IEDIS))
        endIf
      Call FileIO(1,-IRwSim,1,RMSdP,2)
      Call FileIO(1,-IRwSim,NPDIIS**2,V(IDIS),IQNDIS)
      Call FileIO(1,-IRwSim,NPDIIS**2,V(IEDIS),IEnDIS)
      If(RHF) then
        Call FileIO(2,-IRwPA,NTT,Scr,0)
        Call AScale(NTT,Two,Scr,Scr)
        Call FileIO(1,-IRwPA,NTT,Scr,0)
        endIf
      Return
      End
*Deck FOSmF1
      Subroutine FOSmF1(IOut,IBas,NBasis,NBsUse,IRwV,IRwP,IRwFP,IRwF,
     $  IRwPO,NMax,NP,RMSdP,PDIS,EDIS,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Transform the Fock matrices to orthogonal basis for FoSimult.
C
      Dimension V(MDV),VDummy(1),IDummy(1),PDIS(NP,NP),EDIS(NP,NP)
      Save Half, One
      Data Half/0.5D0/, One/1.0D0/
C
C     Initilization
C
      NTT = (NBasis*(NBasis+1))/2
      NBSq = NBasis**2
      If(NP.eq.1) then
        Call ConDDF(IRwPO,3*NMax*NTT)
        Call ConDDF(IRwFP,NMax*NTT)
        endIf
      ID    = 1
      IVX   = ID  + NBSq
      IVV   = IVX + NBSq
      IScr1 = IVV + NBSq
      IScr2 = IScr1 + NBSq
      IEnd  = IScr2 + NBSq
      Call TstCor(IEnd,MDV,'FOSmF1')
C
C     Transform the Fock matrix into the orthonormal basis
C     Store Fock matrices in orthonormal basis
C
      Call FileIO(2,-IRwF,NTT,V(ID),0)
      Call FileIO(2,-IRwV,NBSq,V(IVX),0)
      Call OTran(V(ID),V(IVX),V(IVV),NBasis,NBsUse,0,IBas,V(IEnd),
     $  MDV-IEnd)
      Call AMove(NTT,V(ID),V(IVV))
      Call TstCor(IScr1+NP*NTT,MDV,'FOSmF1')
      If(NP.ge.2) then
        Call FileIO(2,-IRwPO,NTT*(NP-1),V(IScr1),0)
        Call AMove1((NP-1)*NTT,0,NTT,V(IScr1))
        endIf
      Call AMove(NTT,V(ID),V(IScr1))
      Call FileIO(1,-IRwPO,NTT*NP,V(IScr1),0)
      IRwTmp = IGetRW(NTT,0)
      Call FileIO(1,-IRwTmp,NTT,V(IVV),0)
C
C     Transform the Hamiltonian matrix into the orthonormal basis
C     Store Hamiltonian matrices in orthonormal basis
C
      Call FileIO(2,-IRwF,NTT,V(ID),0)
      Call FileIO(2,-515,NTT,V(IVX),0)
      Call AAdd(NTT,V(ID),V(IVX),V(ID))
      Call AScale(NTT,Half,V(ID),V(ID))
      Call FileIO(2,-IRwV,NBSq,V(IVX),0)
      Call OTran(V(ID),V(IVX),V(IVV),NBasis,NBsUse,0,IBas,V(IEnd),
     $  MDV-IEnd)
      Call TstCor(IScr1+NP*NTT,MDV,'FOSmF1')
      If(NP.ge.2) then
        Call FileIO(2,-IRwPO,NTT*(NP-1),V(IScr1),2*NMax*NTT)
        Call AMove1((NP-1)*NTT,0,NTT,V(IScr1))
        endIf
      Call AMove(NTT,V(ID),V(IScr1))
      Call FileIO(1,-IRwPO,NTT*NP,V(IScr1),2*NMax*NTT)
C
C     Transform the density matrix into the orthonormal basis
C     Store density matrices in orthonormal basis
C
      Call FileIO(2,-IRwP,NTT,V(ID),0)
      Call FileIO(2,-IRwV,NBSq,V(IVX),NBSq)
      Call OTran(V(ID),V(IVX),V(IScr1),NBasis,NBsUse,0,IBas,V(IEnd),
     $  MDV-IEnd)
      If(NP.ge.2) then
        Call FileIO(2,-IRwPO,NTT*(NP-1),V(IScr1),NMax*NTT)
        Call AMove1((NP-1)*NTT,0,NTT,V(IScr1))
        endIf
      Call AMove(NTT,V(ID),V(IScr1))
      Call FileIO(1,-IRwPO,NTT*NP,V(IScr1),NMax*NTT)
C
C     Calculate and store the force with respect to the
C     density matrix.
C
      Call FileIO(2,-IRwTmp,NTT,V(IVV),0)
      Call dEdP(VDummy,1,Junk,IDummy,V(ID),V(IVV),V(IVX),V(IScr1),
     $  V(IScr2),NBasis,NTT,.False.,.True.)
      Call AScale(NTT,-One,V(IVX),V(IVX))
      If(NP.ge.2) then
        Call FileIO(2,-IRwFP,NTT*(NP-1),V(IVV),0)
        Call AMove1((NP-1)*NTT,0,NTT,V(IVV))
        endIf
      Call AMove(NTT,V(IVX),V(IVV))
      Call FileIO(1,-IRwFP,NTT*NP,V(IVV),0)
C
C     Begin to Construct the DIIS matrix for the density step
C     Compute approximate Hessian using only the diagonal elements
C     of P and F, dP is calculated using ForceP/SQRT(H)
C
      Call FileIO(2,-IRwTmp,NTT,V(IVV),0)
      Call FreeRW(IRwTmp)
      Call PHess(IOut,NBasis,V(ID),V(IVV),V(IVX),V(IScr1))
      Call ASqrt(0,NTT,V(IVX),V(IVX))
      Call FileIO(2,-IRwFP,NTT*NP,V(IVV),0)
      Do 15 I = 1, NP
        IShift = IVV + (I-1)*NTT
        Call ADiv(NTT,V(IShift),V(IVX),V(IShift))
   15   Continue
      Do 25 I = 1, NP
        Do 25 J = 1, I
          IShftI = IVV + (I-1)*NTT
          IShftJ = IVV + (J-1)*NTT
          PDIS(I,J) = SCFTrc(V(IShftI),V(IShftJ),NBasis,1) + PDIS(I,J)
   25     PDIS(J,I) = PDIS(I,J)
C
C     Prepare EDIIS e matrix
C
      IScr = 1
      IPN  = IScr1 + NTT
      IFN  = IPN   + NP*NTT
      Call FileIO(2,-IRwFP,NP*NTT,V(IFN),0)
      Call AScale(NP*NTT,-One,V(IFN),V(IFN))
      Call FileIO(2,-IRwPO,NP*NTT,V(IPN),NMax*NTT)
      Call ASub(NTT,V(IPN),V(IPN+NTT),V(IScr))
      RMSdP = RMSdP+Sqrt(SCFTrc(V(IScr),V(IScr),NBasis,1)/
     $  GFloat(NBasis))
      Do 30 I = 1, NP
        IFShft = IFN + (I-1)*NTT
        Do 30 J = 1, NP
          JPShft = IPN + (J-1)*NTT
          EDIS(I,J) = EDIS(I,J) + SCFTrc(V(IFShft),V(JPShft),NBasis,1)
   30   Continue
      Return
      End
*Deck FrmBMS
      Subroutine FrmBMS(NSaved,NBasis,Len,NZP,BMat,NDim,IOErr,A,B)
      Implicit Real*8(A-H,O-Z)
C
      Dimension BMat(NDim,NDim),A(*),B(*)
      Save Zero,One
      Data Zero/0.D0/,One/1.0D0/
C
C     Construct DIIS B-Matrix
C
      Call TStamp(2,'FormBM Start')
      Do 10 I = 0, NSaved
        If(I.ge.1) Call FileIO(2,-IOErr,Len,A,(I-1)*Len)
        Do 10 J = 0, I
          If(I.eq.0) then
            BMat(1,1) = Zero
          else if(J.eq.0) then
            BMat(I+1, J+1) = -One
            BMat(J+1, I+1) = -One
          else
            Call FileIO(2,-IOErr,Len,B,(J-1)*Len)
            BMat(I+1, J+1) = FxSprd(A,B,NBasis,NZP,.true.)
            BMat(J+1, I+1) = BMat(I+1, J+1)
          endif
 10     Continue
      Call TStamp(2,'FormBM Stop')
      Return
      End
*Deck FrmTSB
      Subroutine FrmTSB(IOut,IPrint,NRI,NBD,NBlks,NBFU,NOccSy,IRwS,IOC,
     $  IndC,IOSInv,InSInv,IRwP,CMO,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Form the blocked T and S matrices over fragments.
C
      Parameter (MinPrt=5)
      Logical DoPrnt
      Dimension NBFU(*), NOccSy(NBlks,*), CMO(*), V(*)
C
      If(NRI.ne.1) Call GauErr('Complex in FrmTSB.')
      DoPrnt = IPrint.ge.MinPrt
      NOrb = IArSum(1,NBlks,NOccSy)
      NTTR = NRI*(NBD*(NBD+1))/2
      IBlS = 1
      IBlR = IBlS + NTTR
      IBlT = IBlR + NTTR
      IBlInv = IBlT + NBD*NOrb
      IV = IBlInv + NOrb*NOrb
      Call TstCor(IV,MDV,'FrmTSB')
      MDV1 = MDV - IV + 1
      Call FileIO(2,-IOC,NBD*NBD,CMO,IndC)
      If(DoPrnt) Call LTOutX(IOut,1,0,1,0,'New C',NBD,NBD,CMO,CMO,0)
      Call PcckB(NBD,NBlks,NBFU,NOccSy,CMO,V(IBlT))
      Call AClear(NTTR,V(IBlS))
      Call FileIO(2,-IRwS,NTTR,V(IBlR),0)
      If(DoPrnt) Call LTOutX(IOut,0,0,1,0,'Overlap',NBD,NBD,V(IBlR),
     $  V(IBlR),0)
      Call FrmBlS(IOut,DoPrnt,IRwP,NBD,NOrb,V(IBlR),V(IBlT),V(IBlS),
     $  V(IBlInv),V(IV),MDV1)
      Call FileIO(1,-IOSInv,NOrb**2,V(IBlInv),InSInv)
      Return
      End
*Deck FsPurf
      Subroutine FsPurf(IOut,IPrint,IP,IScr,IBT,Pa,P2,P3,Scr,BT,IMap,
     $  JMap,BnFac,ISizCV,IRowNo,IRowA,JNew,IX,iRABin,iScr2,XRowA,ANew,
     $  Scr1,RowSr,iBinB,iBinBt,NBin,NOcc,N,NZ,NP,ZThr,Method,IniMA)
      Implicit Real*8(A-H,O-Z)
      Logical IniMA
      Dimension IP(N+1,*),IMap(*),JMap(*),IRowNo(N)
      Dimension Pa(N,*),P2(N,*),P3(N,*),Scr(N,*),JJA(1)
      Dimension IRowA(N,NP),JNew(N,NP),XRowA(N,NP),ANew(N,NP)
      Dimension IScr(N+1,*)
      Dimension IBT(N+1,*),BT(N,*),ISizCV(NBin),BnFac(NBin)
      Dimension IX(N,NP),iRABin(N,NP),iScr2(N,NP)
      Dimension Scr1(N,NP),RowSr(N,NP),iBinB(NBin+1,N),iBinBt(NBin+1,N)
      Save JJA
      Data JJA/0.d0/
C
C     McWeeny Purification
C
      If(Method.eq.2) Call IMove(N+1+NZ,IP,IScr)
      Call AMove(N+NZ,Pa,P3)
      Call FSMpy(IP,P3,IMap,jMap,JJA(1),JJA(1),IScr,PA,JJA(1),JJA(1),
     $  IBT,BT,IP,P2,IRowA,XRowA,JJA(1),JJA(1),iRABin,JJA(1),JJA(1),
     $  iBinB,iBinBt,JJA(1),JNew,ANew,JJA(1),IX,iRowNo,Scr1,BnFac,
     $  ISizCV,RowSr,iScr2,N,NP,NBin,NZC,ZThr,.true.,.true.,.true.,
     $  .false.,.true.,.true.,.true.,.true.,.true.,.true.,IniMA,
     $  .true.,1,3,Method)
      IniMA = .false.
      Call FSMpy(IP,P2,IMap,jMap,JJA(1),JJA(1),IScr,PA,JJA(1),JJA(1),
     $  IBT,BT,IP,P3,IRowA,XRowA,JJA(1),JJA(1),iRABin,JJA(1),JJA(1),
     $  iBinB,iBinBt,JJA(1),JNew,ANew,JJA(1),IX,iRowNo,Scr1,BnFac,
     $  ISizCV,RowSr,iScr2,N,NP,NBin,NZC,ZThr,.true.,.true.,.true.,
     $  .false.,.false.,.false.,.false.,.false.,.false.,.false.,IniMA,
     $  .true.,1,3,Method)
      If(Method.eq.2) Then
        Call SpTrp(IScr,Pa(1,2),IBT,BT,N)
        Call SpTrp(IBT,BT,IScr,Pa(1,2),N)
      endIf
      Call AScale((N+NZ),(3.d0),P2,P2)
      Call AScale((N+NZ),(-2.d0),P3,P3)
      Call AAdd((N+NZ),P3,P2,Pa)
      Return
      End
*Deck FuToSp
      Subroutine FuToSp(IPrint,SpDiag,SpFock,ZDO,NBasis,ZThrP,IRWPA2,
     $  IRWFA2,IRWS2,IRWV2,IRWVi2,IRWFO2,IRWPO2,X,LenX,IRWFM1,IRWFM2)
      Implicit Real*8(A-H,O-Z)
C
C     Set up for use of full matrix Fock formation followed by sparse
C     diagonalization.
C
      Logical SpDiag, SpFock, ZDO
      Dimension X(LenX)
C
      If(SpDiag) then
        IRwFM1 = IGetRW(1,0)
        IRwFM2 = IGetRW(1,0)
      else
        IRwFM1 = 0
        IRwFM2 = 0
        endif
      If(SpDiag.and..not.SpFock) then
        IRWPA2 = IGetRW(1,0)
        IRWFA2 = IGetRW(1,0)
        IRWS2 = IGetRW(1,0)
        IRWV2 = IGetRW(1,0)
        IRWVi2 = IGetRW(1,0)
        IRWFO2 = IGetRW(1,0)
        IRWPO2 = IGetRW(1,0)
        IndB = 1 + NBasis**2
        IndAA = IndB + NBasis**2
        IndBB = IndAA + NBasis
        IEnd  =  IndBB + NBasis
        Call TstCor(IEnd,LenX,'FuToSp')
        If(.not.ZDO) then
          Call XRead(IRWS,X,NBasis,NBasis,NBasis,NBasis,1)
          Call RootMt(X,X(IndB),X(IndAA),X(IndBB),NBasis,NBasis,0)
          Call ConDDF(IRWV2,NBasis*NBasis)
          Call ConDDF(IRWVi2,NBasis*NBasis)
          Call XWrite(IRWV2,X,NBasis,NBasis,NBasis,NBasis,0)
          Call RootMt(X,X(IndB),X(IndAA),X(IndBB),NBasis,NBasis,-1)
          Call XWrite(IRWVi2,X,NBasis,NBasis,NBasis,NBasis,0)
          Call MkSprL(X,LenX,IRWS,0,IRWS2,1,jia,ian1,iend,NBasis,NZA,
     $      'A.O Overlap Matrix',zthrP,IPrint)
          Call MkSprQ(X,LenX,IRWVi2,0,IRWVi2,1,jia,ian1,iend,NBasis,
     $      NZA,'A.O Inv 1/2 Overlap',zthrP,IPrint)
          Call MkSprQ(X,LenX,IRWV2,0,IRWV2,1,jia,ian1,iend,NBasis,NZA,
     $      'A.O Inv 1/2 Overlap',zthrP,IPrint)
          endIf
      else
        IRWPA2 = 0
        IRWFA2 = 0
        IRWS2 = 0
        IRWV2 = 0
        IRWVi2 = 0
        IRWFO2 = 0
        IRWPO2 = 0
        endIf
      Return
      End
*Deck FxInit
      Subroutine FxInit(IOut,IPrint,cc,maxcor,IRWPO,IRWFO,IRWFM,N,ZThr,
     $  ibegin,jip,ip,if,NZ,IAftrP,iend,cut,Geo,NFirst,NLast,NAtoms,
     $  IMemOp,IVShft,IFxFm,NOcc)
      Implicit Real*8(A-H,O-Z)
      Dimension CC(maxcor),Geo(*),NFirst(*),NLast(*)
      Logical NewF
      Save Zero,One
      Data Zero,One/0.d0,1.d0/
C
C     Set up Fock and Density matrix for Fixed
C
      istart = ibegin
      call SpGet(cc,cc,maxcor,IRWPO,istart,jip,ip,iend,N,NZP,1,.true.,0)
      istart = iend
      call SpGet(cc,cc,maxcor,IRWFO,istart,jif,if,iend,N,NZF,1,.true.,0)
      ix = iend
      iend = ix + N
      call TstCor(iend,maxcor,'FxInit 1')
      If(IPrint.ge.1) Then
        call SpDiff(cc(jip),'P',cc(jif),'F',cc(ix),N)
        call SpDiff(cc(jif),'F',cc(jip),'P',cc(ix),N)
        endif
      If(Cut.gt.Zero) then
        ISupp = 0
        Call ThrCut(Geo,NFirst,NLast,cc(ix),cc(jif),cc(if+n),cut,N,
     $    NAtoms,nzf)
      else
        ISupp = IFxFm + 1
        endif
C
C     Make Fock and Density share the same form
C
      istart = iend
      NewF = IRWFM.eq.0
      If(NewF) then
        Call FPenh(cc,maxcor,cc(jif),cc(if),cc(jip),cc(ip),cc(ix),zthr,
     $    istart,jipf,ipp,iff,iend,N,NZ,ISupp,.true.,.true.,iprint,
     $    IMemOp,NOcc)
      else
        Call SpGet(cc,cc,maxcor,IRWFM,istart,jipf,ipp,iend,N,NZ,0,
     $    .true.,0)
        ipp = iend
        iff = ipp + NZ + N
        iend = iff + NZ + N
        call AMove(N,cc(if),cc(iff))
        call AMove(N,cc(ip),cc(ipp))
        call NwForm(cc(jif),cc(if+N),cc(jipf),cc(iff+N),cc(ix),N)
        call NwForm(cc(jip),cc(ip+N),cc(jipf),cc(ipp+N),cc(ix),N)
        call TstCor(iend,maxcor,'FxInit 2')
        endIf
      If(iprint.ge.1) Call SpRprt(NZ,(N*(N-1))/2,zthr,'All matrices')
C
C     Allocate memory for all matrices used by CG-DMS
C
      jip   = jipf
      ip    = ipp
      Call SpAMvM(cc,MaxCor,ibegin,if,jip,ip,N,2)
      iend  = if   + NZ+N
      ix    = iend
      itt   = ix    + N
      last  = itt + NZ+N
      IAftrP = if
      Call TstCor(last,maxcor,'FxInit 3')
      call amove1(NZ+N,iff-1,if-1,cc)
      If(IVShft.ne.0) then
        VShift = GFloat(IVShft)
        Call ASet(N,-One,cc(iX))
        Call AAdd(N,cc(ip),cc(iX),cc(itt))
        Call AMove(NZ,cc(ip+N),cc(itt+N))
        Call AScale((N+NZ),-VShift,cc(itt),cc(itt))
        Call AAdd((N+NZ),cc(itt),cc(if),cc(if))
        endif
      call CheckF(cc(jip),cc(ip),cc(ix),N,0,NZ,'Ortho Density',5,IPrint)
      call CheckF(cc(jip),cc(if),cc(ix),N,0,NZ,'Ortho Fock',5,IPrint)
      Return
      End
*Deck GBnds
      Subroutine GBnds(IA,A,Norbs,NZA,AAbs,AMax,AMin,IOut,RDisc)
      Implicit Real*8(A-H,O-Z)
C
C     Find maximum (AMax) , minimum (AMin) and
C     maximum abs. value eigenvalue (AAbs) of a symmetric matrix
C     stored in lower triangle with diagonals
C     using gershgorin's eigenvalue estim.
C
      Dimension A(Norbs,2),RDisc(Norbs),IA(Norbs+1,2)
      Save Zero
      Data Zero/0.d0/
C
C     Sumc : contains maximum of the estimated eigenvalues
C     Sumcc : contains minimum of the estimated eigenvalues
C
      Call AClear(Norbs,RDisc)
      Do 20 i = 1,Norbs
        Sumr = Zero
        Do 10 j = IA(i,1), (IA(I+1,1)-1)
          jj = IA(j,2)
          x = Abs(A(j,2))
C         Row's contr. to maximum of i - th eig
          Sumr = Sumr + x
C         Colmn's contr. to the radius of gershgorin's disc
   10     RDisc(jj) = RDisc(jj) + x
   20   RDisc(i) = RDisc(i) + Sumr
      AMax = A(1,1) + RDisc(1)
      AMin = A(1,1) - RDisc(1)
      Do 30 i = 2,NOrbs
        XMax = A(i,1) + RDisc(i)
        XMin = A(i,1) - RDisc(i)
        If(XMax.gt.AMax) AMax = XMax
        If(XMin.lt.AMin) AMin = XMin
   30   Continue
      AAbs = Max(Abs(AMax),Abs(AMin))
      Return
      End
*Deck GetCGX
      Subroutine GetCGX(IOut,IPrint,NTot,NTotR,A,b,q,Thresh,ffact,
     $  weight,EFG,Phi,xyz,gnorm,SFact,GFact,Beta,bmf,Alpha,qtot,ICen,
     $  Scr,ISubSt,ISlMap,Neighb,DoPrun,PtThr,PtInc,UseOld,MaxSFa,
     $  IRwSl2,LScSl2,SymCut,V,MDV,Const1,Constc,IScal,IMem,ISc1,WgtOld,
     $  XYZOld)
      Implicit Real*8 (A-H,O-Z)
      Parameter (MaxInv=0)
      Logical UseOld, DoPrun
      Real*8 MDCutO
      Dimension A(*),b(NTot,2),q(NTot,ISc1*ISc1),ffact(NTot),
     $  weight(NTot),EFG(Ntot,10,ISc1),Phi(NTot),xyz(NTot,3),
     $  Gnorm(NTot,3),SFact(*),GFact(*),Beta(*),bmf(NTot,ISc1),
     $  Alpha(NTot,2),QTot(NTot),ICen(*),Scr(NTot,10),ISubSt(NTot),
     $  ISlMap(NTot),Neighb(NTot),V(MDV),WgtOld(NTot),
     $  XYZOld(NTot,3)
      Dimension IBBuc(2),IXBuc(2),Local1(2),Local2(2),Minp(2),Nin(2),
     $  Ninp(2),NDone(2),IInBuc(2),IZero(2)
      Save Zero, One, Two, Four, PtMax
      Data Zero/0.d0/, One/1.d0/, Two/2.0d0/, Four/4.0d0/,
     $  PtMax/20.0d0/
 9000 Format(' Using previous values as initial guesses in getcgx.')
 9010 Format(1x,'Could not find pruning parameters that gave all +ve',
     $          ' SFacts.',/,1x,'PtThr: ',f5.3)
 9020 Format(3(1x,'Subspace ',i2,1x,i4,' passes.',:))
 9030 format(1x,i4,4(d17.10,1x))
 9050 Format(1x,'Using a subset of ',i8,' of the ',i8,
     $  ' points (',f5.1,'%).')
C
C     If IMem=0 then A is empty. All values will be computed on the fly.
C        IMem=1 then A has space for the distance matrix for the solvent points.
C        IMem=2 then A has space for the whole A or C matrix.
C        IMem=3 then A has space for A and C, in that order.
C          (Still have to create A and C separately because C depends
C           on q which requires A).
C
C     All solvent linear equations are solved here.
C
C
      TwoPi = Two*Four*ATan(One)
      Cnvc1 = Thresh
      Cnvc2 = Cnvc1*GFloat(10)
      Cnvc3 = MDCutO(0)
      Cnvc4 = Cnvc3*GFloat(10)
      Cnvc5 = Cnvc1
      Cnvc6 = Cnvc5*GFloat(10)
      MaxIt = NTot
      Call GetRFA(IMem,11,NTot,NTot,TwoPi,XYZ,GNorm,SFact,FFact,Weight,
     $  q,A)
      Do 20 ISub = 1, 2
        Ninp(ISub) = 1
        Minp(ISub) = NTot
        NDone(ISub) = 0
        IBBuc(ISub) = IGetRW(ISc1*NTot,0)
        IXBuc(ISub) = IGetRW(ISc1*NTot,0)
        Local1(ISub) = IGetRW(ISc1*5*NTot,0)
        Local2(ISub) = IGetRW(ISc1*5*NTot,0)
   20   IInBuc(ISub) = IGetRW(2*ISc1*NTot,0)
      If(UseOld.and.IPrint.ge.1) Write(IOut,9000)
      Call FileIO(1,-IRwSl2,0,Junk,LScSl2)
      Call FileIO(2,-IRwSl2,0,Junk,LScSl2)
C
C     No scaling case.
C
C     Compute Q, bmf
C
      If(IScal.eq.0) then
C
C       Pruning code in unscaled SCIPCM purely for testing. Set PtThr
C       to some non-zero value to prune points. NTotR = NTot if no pruning.
C
        Call AMove(NTot,Weight,WgtOld)
        Call AMove(NTot*3,XYZ,XYZOld)
        MDVI = IntPWP(1) * MDV
        If(.not.DoPrun) PtThr = Zero
        Call SolPru(IOut,IPrint,SymCut,NTot,NTotR,PtThr,XYZ,Weight,
     $    ISubSt,ISlMap,Neighb,Scr(1,1),Scr(1,2),Scr(1,2),Scr(1,3),V,
     $    MDVI)
C
C     Create rest of smallspace.
C
        If(NTotR.ne.NTot) then
          Do 40 I = 1, NTot
            M = ISlMap(I)
            If(M.ne.0) then
              Scr(M,1) = XYZ(I,1)
              Scr(M,2) = XYZ(I,2)
              Scr(M,3) = XYZ(I,3)
              Endif
   40       Continue
            Do 50 II = 1, NTotR
              XYZ(II,1) = Scr(II,1)
              XYZ(II,2) = Scr(II,2)
   50         XYZ(II,3) = Scr(II,3)
C
            Do 70 I = 1, NTot
              M = ISlMap(I)
              If(M.ne.0) then
                Do 60 Lp = 1, 10
   60             Scr(M,Lp) = EFG(I,Lp,1)
                Endif
   70         Continue
            Do 80 Lp = 1, 10
              Do 80 II = 1, NTotR
                EFG(II,Lp,1) = Scr(II,Lp)
   80           Continue
C
            Do 90 I = 1, NTot
              M = ISlMap(I)
              If(M.ne.0) then
                Scr(M,1) = b(I,1)
                Scr(M,2) = b(I,2)
                Scr(M,3) = GFloat(ICen(I))
                endif
   90         Continue
            Do 100 II = 1, NTotR
              b(II,1) = Scr(II,1)
              b(II,2) = Scr(II,2)
              ICen(II) = Int(Scr(II,3))
  100         Continue
C
            Call AClear(5*NTot,Scr)
            Do 110 I = 1, NTot
              M = ISlMap(I)
              If(M.ne.0) then
                Scr(M,6) = GNorm(I,1)
                Scr(M,7) = GNorm(I,2)
                Scr(M,8) = GNorm(I,3)
                Scr(M,9) = Phi(I)
                Else
                  M = ISlMap(ISubst(I))
                Endif
              Scr(M,5) = Scr(M,5) + Weight(I)
  110         Continue
            Do 120 II = 1, NTotR
              GNorm(II,1) = Scr(II,6)
              GNorm(II,2) = Scr(II,7)
              GNorm(II,3) = Scr(II,8)
              Phi(II) = Scr(II,9)
              Weight(II) = Scr(II,5)
  120         Continue
          endif
C         endif matches NTotR.ne.NTot
        Call GetRFA(IMem,11,NTot,NTotR,TwoPi,XYZ,GNorm,SFact,FFact,
     $    Weight,q,A)
C
        Do 130 I = 1, NTotR
          Scr(I,9) = Weight(I)/(Const1+Ffact(I))
          Scr(I,10) = Weight(I)/(Const1+Ffact(I))
          q(I,1) = -b(I,1)*Scr(I,9)
  130     bmf(I,1) = Phi(I)*Scr(I,10)
        Call Fileio(1,-Ibbuc(1),NTotR,q,0)
        Call Fileio(1,-Ibbuc(2),NTotR,bmf,0)
        If(UseOld) then
          Call Fileio(2,IRwSl2,NTotR,q,0)
          Call Fileio(2,IRwSl2,NTotR,bmf,0)
          endif
        Do 140 ISub = 1, 2
          Minp(ISub) = NTotR
  140     Nin(ISub) = 1
        Call Fileio(1,-IInbuc(1),NTotR,q,0)
        Call Fileio(1,-IInbuc(2),NTotR,bmf,0)
        Call LinEqS(0,IMem,A,NTot,XYZ,GNorm,q,Scr(1,9),IOut,IPrint,
     $    IBbuc,IXbuc,IInBuc,Local1,Local2,2,Minp,Nin,Ninp,NDone,MaxInv,
     $    IZero,Cnvc1,Cnvc2,Cnvc3,Cnvc4,MaxIt,V,MDV)
        If(IPrint.ge.1) Write(IOut,9020)(J,NDone(J),J = 1, 2)
        If(IPrint.ge.1) Write(IOut,*)'Solved q, omega'
        Call Fileio(2,-Ixbuc(1),NTotR,q(1,1),0)
        Call Fileio(2,-Ixbuc(2),NTotR,bmf,0)
        Call Fileio(1,IRwSl2,NTotR,q,0)
        Call Fileio(1,IRwSl2,NTotR,bmf,0)
        Call AMove(NTotR,Q(1,1),Qtot)
C
C     Conductor scaling. Only scaling currently working.
C
C     Compute Q1, Q2, Q3, Q4
C
      else if(IScal.eq.2) then
        Call AMove(NTot,Weight,WgtOld)
        Call AMove(NTot*3,XYZ,XYZOld)
        Do 150 I = 1, NTot
          Scr(I,9) = Weight(I) / (Const1+Ffact(I))
          Scr(I,10) = Weight(I) / (Constc+Ffact(I))
C         Charges due to:
C           electronic field in dielectric
          Q(I,1) = -b(I,1)*Scr(I,9)
C           nuclear field in dielectric
          Q(I,2) = -b(I,2)*Scr(I,9)
C           electronic field in conductor
          Q(I,3) = -b(I,1)*Scr(I,10)
C           nuclear field in conductor
  150     Q(I,4) = -b(I,2)*Scr(I,10)
        Do 160 ISub = 1, 2
          Ninp(ISub) = 2
          NDone(ISub) = 0
          Call Fileio(1,-Ibbuc(ISub),2*NTot,Q(1,1+2*(ISub-1)),0)
          If(UseOld) then
            Call Fileio(2,IRwSl2,2*NTot,Q(1,1+2*(ISub-1)),0)
            endif
          Call Fileio(1,-IInbuc(ISub),2*NTot,Q(1,1+2*(ISub-1)),0)
  160     Nin(ISub) = 2
        Call LinEqS(1,IMem,A,NTot,XYZ,GNorm,q,Scr(1,9),IOut,IPrint,
     $    IBbuc,IXbuc,IInBuc,Local1,Local2,2,Minp,Nin,Ninp,Ndone,MaxInv,
     $    IZero,Cnvc1,Cnvc2,Cnvc3,Cnvc4,MaxIt,V,MDV)
        If(IPrint.ge.1) Write(IOut,9020)(J,NDone(J),J = 1, 2)
        If(IPrint.ge.1) Write(IOut,*)'Solved q(1-4)'
        Do 170 ISub = 1, 2
          Call Fileio(2,-Ixbuc(ISub),2*NTot,q(1,1+(ISub-1)*2),0)
  170     Continue
        Call Fileio(1,-IRwSl2,4*NTot,Q,LScSl2)
        If(IPrint.ge.3) then
          write(IOut,*)'New way: Q(i,1-4)'
          do 180 i = 1, ntot
          write(IOut,9030)i,q(i,1),q(i,2),q(i,3),q(i,4)
  180     continue
          endif
C
C       Now generate GFact (actual scale factors for electronic charges)
C       and Beta (pre-mult factor).
C
C       Create map for coincident points.
C
        MDVI = IntPWP(1) * MDV
        If(.not.DoPrun) PtThr = Zero
C       Do Until(there are less than MaxSFa -ve SFact)
  190   Call SolPru(IOut,IPrint,SymCut,NTot,NTotR,PtThr,XYZ,Weight,
     $    ISubSt,ISlMap,Neighb,Scr(1,1),Scr(1,2),Scr(1,2),Scr(1,3),V,
     $    MDVI)
C
C       Create smallspace of not-close points.
C       Compress q's into small space - then set NTot to NTotR and
C         continue in small space for rest of calculation.
C
C       Compute SFacts in this pruned grid.
        Call AClear(5*NTot,Scr)
        Do 200 I = 1, NTot
          M = ISlMap(I)
          If(M.ne.0) then
            Scr(M,1) = XYZ(I,1)
            Scr(M,2) = XYZ(I,2)
            Scr(M,3) = XYZ(I,3)
            Scr(M,4) = EFG(I,1,2)
          Else
            M = ISlMap(ISubst(I))
            Endif
          Scr(M,5) = Scr(M,5) + Q(I,4)
  200     Continue
        Call MkSFac(NTot,NTotR,Scr,Scr(1,5),Scr(1,4),SFact,Negs)
        If(IPrint.ge.1) Write(IOut,*)'Negative SFacts: ', Negs
        If(Negs.gt.MaxSFa.and.DoPrun) then
          If(IPrint.ge.1) Write(IOut,*)
     $      'Found negative SFact. Re-pruning solvent surface grid.'
          PtThr = PtThr + PtInc
          If(PtThr.gt.PtMax) then
            Write(IOut,9010) PtThr
            Call GauErr('Something wrong with solvent grid.')
            Endif
          Goto 190
          endIf
C
C     Create rest of smallspace.
C
        If(NTotR.ne.NTot.and..false.) then
          Do 210 II = 1, NTotR
            XYZ(II,1) = Scr(II,1)
            XYZ(II,2) = Scr(II,2)
  210       XYZ(II,3) = Scr(II,3)
C
          Do 230 I = 1, NTot
            M = ISlMap(I)
            If(M.ne.0) then
              Do 220 Lp = 1, 10
  220           Scr(M,Lp) = EFG(I,Lp,1)
              Endif
  230       Continue
          Do 240 Lp = 1, 10
            Do 240 II = 1, NTotR
              EFG(II,Lp,1) = Scr(II,Lp)
  240         Continue
C
          Do 260 I = 1, NTot
            M = ISlMap(I)
            If(M.ne.0) then
              Do 250 Lp = 1, 10
  250           Scr(M,Lp) = EFG(I,Lp,2)
              Endif
  260       Continue
          Do 270 Lp = 1, 10
            Do 270 II = 1, NTotR
              EFG(II,Lp,2) = Scr(II,Lp)
  270         Continue
C
          Do 280 I = 1, NTot
            M = ISlMap(I)
            If(M.ne.0) then
              Scr(M,1) = b(I,1)
              Scr(M,2) = b(I,2)
              Scr(M,3) = GFloat(ICen(I))
              endif
  280       Continue
          Do 290 II = 1, NTotR
            b(II,1) = Scr(II,1)
            b(II,2) = Scr(II,2)
            ICen(II) = Int(Scr(II,3))
  290       Continue
C
          Call AClear(5*NTot,Scr)
          Do 310 I = 1, NTot
            M = ISlMap(I)
            If(M.ne.0) then
              Scr(M,6) = GNorm(I,1)
              Scr(M,7) = GNorm(I,2)
              Scr(M,8) = GNorm(I,3)
              Scr(M,9) = Phi(I)
              Else
                M = ISlMap(ISubst(I))
              Endif
            Do 300 Lp = 1, 4
  300         Scr(M,Lp) = Scr(M,Lp) + Q(I,Lp)
            Scr(M,5) = Scr(M,5) + Weight(I)
  310       Continue
          Do 320 II = 1, NTotR
            GNorm(II,1) = Scr(II,6)
            GNorm(II,2) = Scr(II,7)
            GNorm(II,3) = Scr(II,8)
            Phi(II) = Scr(II,9)
            Weight(II) = Scr(II,5)
  320       Continue
C
C         Re-compute q, SFact, and A or dist matrix, in small space
C
          If(IPrint.ge.1)
     $      Write(IOut,9050)NTotR,NTot,GFloat(NTotR*100)/GFloat(NTot)
          Call GetRFA(IMem,11,NTot,NTotR,TwoPi,XYZ,GNorm,SFact,FFact,
     $      Weight,q,A)
          Call Fileio(1,-IInbuc(1),NTotR,Q(1,1),0)
          Call Fileio(1,IInbuc(1),NTotR,Q(1,2),0)
          Call Fileio(1,-IInbuc(2),NTotR,Q(1,3),0)
          Call Fileio(1,IInbuc(2),NTotR,Q(1,4),0)
          Do 330 II = 1, NTotR
            Scr(II,9) = Weight(II) / (Const1+Ffact(II))
            Scr(II,10) = Weight(II) / (Constc+Ffact(II))
            Q(II,1) = -b(II,1)*Scr(II,9)
            Q(II,2) = -b(II,2)*Scr(II,9)
            Q(II,3) = -b(II,1)*Scr(II,10)
  330       Q(II,4) = -b(II,2)*Scr(II,10)
          Do 340 ISub = 1, 2
            Ninp(ISub) = 2
            NDone(ISub) = 0
            Minp(ISub) = NTotR
  340       Nin(ISub) = 2
          Call Fileio(1,-Ibbuc(1),NTotR,Q(1,1),0)
          Call Fileio(1,Ibbuc(1),NTotR,Q(1,2),0)
          Call Fileio(1,-Ibbuc(2),NTotR,Q(1,3),0)
          Call Fileio(1,Ibbuc(2),NTotR,Q(1,4),0)
          Call LinEqS(1,IMem,A,NTot,XYZ,GNorm,q,Scr(1,9),IOut,IPrint,
     $      IBbuc,IXbuc,IInBuc,Local1,Local2,2,Minp,Nin,Ninp,Ndone,
     $      MaxInv,IZero,Cnvc1,Cnvc2,Cnvc3,Cnvc4,MaxIt,V,MDV)
          If(IPrint.ge.1) Write(IOut,9020)(J,NDone(J),J = 1, 2)
          If(IPrint.ge.1) Write(IOut,*)'Solved q(1-4) again'
          Call Fileio(2,-Ixbuc(1),NTotR,q(1,1),0)
          Call Fileio(2,Ixbuc(1),NTotR,q(1,2),0)
          Call Fileio(2,-Ixbuc(2),NTotR,q(1,3),0)
          Call Fileio(2,Ixbuc(2),NTotR,q(1,4),0)
          If(IPrint.ge.3) then
            write(IOut,*)'New way: Q(i,1-4) again'
            do 350 II = 1, NTotR
            write(IOut,9030)II,q(II,1),q(II,2),q(II,3),q(II,4)
  350       continue
            endif
          Call MkSFac(NTot,NTotR,XYZ,Q(1,4),EFG(1,1,2),SFact,Negs)
          If(IPrint.ge.2)
     $      Write(IOut,*)'Negative SFacts in second q: ', Negs
C         this Endif matches If(NTotR.ne.NTot)
          Endif
C
C       Compute Gfact, Beta
C
        ICMat = 1
        If(IMem.eq.3) ICMat = 1 + NTotR**2
        If(IMem.eq.2.or.IMem.eq.3)
     $    Call GetRFA(IMem,2,NTot,NTotR,TwoPi,XYZ,GNorm,SFact,FFact,
     $      Weight,Q(1,3),A(ICMat))
        Do 360 II = 1, NTotR
          Scr(II,9) = One / SFact(II)
          Scr(II,10) = One / (SFact(II) * Q(II,3))
          Scr(II,2) = Phi(II) * Q(II,1) * Scr(II,10)
  360     Scr(II,1) = -EFG(II,1,1) * Scr(II,9)
C
        Call Fileio(1,-IBBuc(1),NTotR,Scr(1,1),0)
        Call Fileio(1,-IBBuc(2),NTotR,Scr(1,2),0)
        If(UseOld) then
          Call Fileio(2,IRwSl2,NTot,GFact,0)
          Call Fileio(2,IRwSl2,NTot,Beta,0)
          Call AClear(NTotR,Scr(1,2))
          Do 370 I = 1, NTot
            M = ISlMap(I)
            If(M.ne.0) then
              Scr(M,1) = GFact(I)
            Else
              M = ISlMap(ISubst(I))
              Endif
            Scr(M,2) = Scr(M,2) + Beta(I)
  370       Continue
        else
          Do 380 II = 1, NTotR
  380       Scr(II,1) = One
          endif
        Do 390 II = 1, NTotR
  390     Scr(II,1) = Scr(II,1) * Q(II,3)
        Call Fileio(1,-IInBuc(1),NTotR,Scr(1,1),0)
        Call Fileio(1,-IInBuc(2),NTotR,Scr(1,2),0)
        Do 400 ISub = 1, 2
          Minp(ISub) = NTotR
          Ninp(ISub) = 1
          NDone(ISub) = 0
  400     Nin(ISub) = 1
        Call LinEqS(2,IMem,A(ICMat),NTot,XYZ,GNorm,Q(1,3),Scr(1,9),
     $    IOut,IPrint,IBbuc,IXbuc,IInBuc,Local1,Local2,2,Minp,Nin,Ninp,
     $    NDone,MaxInv,IZero,Cnvc5,Cnvc6,Cnvc3,Cnvc4,MaxIt,V,MDV)
        If(IPrint.ge.1) Write(IOut,9020)(J,NDone(J),J = 1, 2)
        If(IPrint.ge.1) Write(IOut,*)'Solved gfact, beta'
        Call Fileio(2,-IXBuc(1),NTotR,GFact,0)
        Call Fileio(2,-IXBuc(2),NTotR,Beta,0)
        Do 410 II = 1, NTotR
  410     GFact(II) = GFact(II) / Q(II,3)
C
C       Compute Alpha1, Alpha2, bmf1, bmf2
C
        Do 420 II = 1, NTotR
          Scr(II,3) = Beta(II)*Q(II,3)*GFact(II)*SFact(II)
     $                 / Q(II,4)
          Do 420 JJ = 1, (II-1)
            distx = XYZ(II,1) - XYZ(JJ,1)
            disty = XYZ(II,2) - XYZ(JJ,2)
            distz = XYZ(II,3) - XYZ(JJ,3)
            d1 = One/Sqrt(distx**2+disty**2+distz**2)
            Scr(II,3) = Scr(II,3) + Beta(JJ) * Q(JJ,3) * GFact(JJ)
     $        * d1 / Q(JJ,4)
  420       Scr(JJ,3) = Scr(JJ,3) + Beta(II) * Q(II,3) * GFact(II)
     $        * d1 / Q(II,4)
C
        Do 430 II = 1, NTotR
          Scr(II,9) = Weight(II) / (Const1 + FFact(II))
          Scr(II,10) = Weight(II) / (Constc + FFact(II))
          bmf(II,1) = Phi(II) * GFact(II) * Scr(II,9)
          bmf(II,2) = Phi(II) * Scr(II,9)
          Alpha(II,1) = GFact(II) * Q(II,1) * Phi(II)
     $      * Scr(II,10) / Q(II,3)
  430     Alpha(II,2) = Scr(II,3) * Scr(II,10)
        Call FileIo(1,-IBBuc(1),NTotR,bmf,0)
        Call FileIo(1,IBBuc(1),NTotR,bmf(1,2),0)
        Call FileIo(1,-IBBuc(2),NTotR,Alpha,0)
        Call FileIo(1,IBBuc(2),NTotR,Alpha(1,2),0)
        If(UseOld) then
          Call Fileio(2,IRwSl2,2*NTot,Scr(1,1),0)
          Call Fileio(2,IRwSl2,2*NTot,Scr(1,3),0)
          Call AClear(NtotR,bmf(1,1))
          Call AClear(NtotR,bmf(1,2))
          Call AClear(NtotR,Alpha(1,1))
          Call AClear(NtotR,Alpha(1,2))
          Do 440 I = 1, NTot
            M = ISlMap(I)
            If(M.eq.0) M = ISlMap(ISubst(I))
            bmf(M,1) = bmf(M,1) + Scr(I,1)
            bmf(M,2) = bmf(M,2) + Scr(I,2)
            Alpha(M,1) = Alpha(M,1) + Scr(I,3)
            Alpha(M,2) = Alpha(M,2) + Scr(I,4)
  440       Continue
          endif
        Call FileIo(1,-IInBuc(1),NTotR,bmf,0)
        Call FileIo(1,IInBuc(1),NTotR,bmf(1,2),0)
        Call FileIo(1,-IInBuc(2),NTotR,Alpha,0)
        Call FileIo(1,IInBuc(2),NTotR,Alpha(1,2),0)
        If(IMem.eq.2) Call GetRFA(IMem,1,NTot,NTotR,TwoPi,XYZ,GNorm,
     $    SFact,FFact,Weight,q,A)
        Do 450 ISub = 1, 2
          Minp(ISub) = NTotR
          Ninp(ISub) = 2
          NDone(ISub) = 0
  450     Nin(ISub) = 2
        Call LinEqS(3,IMem,A,NTot,XYZ,GNorm,q,Scr(1,9),IOut,IPrint,
     $    IBbuc,IXbuc,IInBuc,Local1,Local2,2,Minp,Nin,Ninp,Ndone,MaxInv,
     $    IZero,Cnvc1,Cnvc2,Cnvc3,Cnvc4,MaxIt,V,MDV)
        If(IPrint.ge.1) Write(IOut,9020)(J,NDone(J),J = 1, 2)
        If(IPrint.ge.1) Write(IOut,*)'Solved omega, alpha'
        Call FileIo(2,-IXBuc(1),NTotR,bmf,0)
        Call FileIo(2,IXBuc(1),NTotR,bmf(1,2),0)
        Call FileIo(2,-IXBuc(2),NTotR,Alpha,0)
        Call FileIo(2,IXBuc(2),NTotR,Alpha(1,2),0)
        Call AMulA(NTotR,Q(1,2),Q(1,1),GFact,qtot)
C
C       Expand the small space vectors back into large space when saving
C       for UseOld.
C
        Do 460 I = 1, NTot
          JJ = ISlMap(ISubst(I))
          Scr(I,1) = GFact(JJ)
          Ratio = WgtOld(I) / Weight(JJ)
          Scr(I,2) = Beta(JJ) * Ratio
          Scr(I,3) = bmf(JJ,1) * Ratio
          Scr(I,4) = bmf(JJ,2) * Ratio
          Scr(I,5) = Alpha(JJ,1) * Ratio
  460     Scr(I,6) = Alpha(JJ,2) * Ratio
        Call Fileio(1,IRwSl2,6*NTot,Scr(1,1),0)
      else
        Call GauErr('Illegal Iscal in GetCGX')
        endif
      Do 470 ISub = 1, 2
        Call FreeRW(IBBuc(ISub))
        Call FreeRW(IXBuc(ISub))
        Call FreeRW(IInBuc(ISub))
        Call FreeRw(Local1(ISub))
        Call FreeRw(Local2(ISub))
  470   Continue
      If(IPrint.ge.3) then
        write(IOut,*)'New way: Sfact, GFact, bmf1, bmf2'
        GMax = Zero
        IGMax = 0
        do 480 i = 1, NTotR
          GF = Abs(GFact(i))
          If(GF.gt.GMax.or.One/GF.gt.GMax) then
            GMax = Max(GF,One/GF)
            IGMax = i
            Endif
  480     write(IOut,9030)i,SFact(i),GFact(i),bmf(i,1),bmf(i,2)
        write(IOut,'('' GFact('',i5,'') = '',d17.10)')IGmax,GFact(IGMax)
        write(IOut,*)'New way: Beta, Alpha1, Alpha2, Phi'
        do 490 i = 1, NTotR
  490   write(IOut,9030)i,Beta(i),Alpha(i,1),Alpha(i,2),Phi(i)
        write(IOut,*)'New way: EFG1, EFG2, Qtot'
        do 500 i = 1, NTotR
  500   write(IOut,9030)i, EFG(I,1,1), EFG(I,1,2), QTot(I)
        endif
      Return
      End
*Deck GetDat
      Subroutine GetDat(IOut,RdInpt,IRwSl2,IDump,NAtoms,IScal,Epsi,Cont,
     $  ITheta,IPhi,IBWFnc,RadT,NUMR,ISurf,MoBas,IPass,CutPt,NTot0,
     $  RdPCM,AtmChg,IOrg,NOrg,DoPrun,PtThr,PtInc)
      Implicit Real*8 (A-H,O-Z)
C
C     IRwSl2 stores information such as Epsilon, NTot, etc, and
C     the q's etc from previous iterations. Its existence means
C     that there is no need to read from the input file.
C
      Logical NoPrDf, RdPCM, RdInpt, MoBas, DoPrun, NoPrun, LFFGet
      Dimension AtmChg(NAtoms), IOrg(NAtoms)
      Parameter (IPhiDf=302,ITheDf=1,ISurDf=3,NoPrDf=.True.,IScaDf=0,
     $  IBWFDf=1)
      Save Zero, One, DfCont, DfCutP, DfRadT, DfPtTh, DfPtI
      Data Zero/0.0d0/, One/1.0d0/, DfCont/0.0004d0/, DfCutP/1.d-6/,
     $  DfRadT/0.1d0/, DfPtTh/0.10d0/, DfPtI/0.02d0/
C
 9000 Format(1x,'No such atomic center: ',I6,'!')
      MoBas = .True.
      CutPt = DfCutP
      PtThr = DfPtTh
C     Disable pruning: non-convergent SCF. PtThr = Zero, NoPrDf = .true.
      PtThr = Zero
      PtInc = DfPtI
      IBWFnc = IBWFDf
      RadT = DfRadT
      NTot0 = 0
      If(IPass.eq.0) then
        Numr = 19
      else
        Numr = 26
        endIf
      If(RdInpt) then
        Cont = DfCont
        IPhi = IPhiDf
        ITheta = ITheDf
        ISurf = ISurDf
        NoPrun = NoPrDf
        IScal = IScaDf
        If(RdPCM) then
          Call FFSet(IDump)
          Call FFRead(IEof)
          If(IEOF.ne.0) Call GauErr(
     $      'End of File Reading Reaction Field Parameters.')
          Epsil = RFFGet(.True.,.False.,Zero,.False.)
          If(Epsil.eq.Zero) then
            If(Epsi.eq.Zero) Call GauErr(
     $ 'Epsi is zero in input, but no built-in Epsi for this solvent.')
          else
C           Read in Epsi silently overrides built-in Epsi.
            Epsi = Epsil
            endif
          Cont = RFFGet(.True.,.True.,DfCont,.False.)
          IPhi = IFFGet(.True.,IPhiDf,.False.)
          ITheta = IFFGet(.True.,ITheDf,.False.)
          ISurf = IFFGet(.True.,ISurDf,.False.)
          NoPrun = LFFGet(.True.,NoPrDf,.False.)
          IScal = IFFGet(.True.,IScaDf,.False.)
          endIf
        Call ISet(1,NAtoms,IOrg)
        If(ISurf.eq.2) then
C         Single origin
          NOrg = 1
        else if(ISurf.eq.3.or.ISurf.eq.4) then
          If(ISurf.eq.4) then
C           All non-H atoms as origins
            NOrg = 0
            Do 10 I = 1, NAtoms
              If(AtmChg(I).gt.One) then
                NOrg = NOrg + 1
                IOrg(NOrg) = I
                endIf
   10         Continue
            endif
          If(ISurf.eq.3.or.NOrg.eq.0) then
C         All atoms as origins
            ISurf = 3
            NOrg = NAtoms
            Do 20 I = 1, NAtoms
   20         IOrg(I) = I
            endif
        else if(ISurf.eq.5) then
C         Read in atom numbers for origins
          Call FFSet(IDump)
          Call FFRead(IEof)
          Call IClear(NAtoms,IOrg)
          NOrg = 0
          Do 30 I = 1, NAtoms
            IDo = IFFGet(.True.,0,.False.)
            If(IDo.ne.0) then
              If(IDo.lt.1.or.IDo.gt.NAtoms) then
                Write(IOut,9000) IDo
                Call GauErr('Invalid nuclear center for solvent grid.')
                endif
              NOrg = NOrg + 1
              IOrg(NOrg) = IDo
              endIf
   30       Continue
          If(NOrg.eq.0) Call GauErr('There are no solvent origins.')
        else
          Call GauErr('Invalid ISurf in GetDat (SCIPCM)')
          endIf
      else
        Call Fileio(2,-IRwSl2,1,XMask,0)
        Call Bnupak(XMask,NoPrun,Iscal)
        Call Fileio(2,IRwSl2,1,epsi,0)
        Call Fileio(2,IRwSl2,1,Cont,0)
        Call Fileio(2,IRwSl2,1,XMask,0)
        Call Bnupak(XMask,IPhi,ITheta)
        Call Fileio(2,IRwSl2,1,XMask,0)
        Call Bnupak(XMask,Isurf,NTot0)
        Call Fileio(2,IRwSl2,1,PtThr,0)
        PtThr = Max(PtThr - GFloat(2) * PtInc, Zero)
        Call Fileio(2,IRwSl2,1,XMask,0)
        Call Bnupak(XMask,NOrg,Junk)
        Call FileIO(2,IRwSl2,InToWp(NOrg),IOrg,0)
        endIf
      DoPrun = .not.NoPrun
      Return
      End
*Deck GetFlx
      Subroutine GetFlx(IOut,IPrint,IPFlag,NBas6D,NAtoms,CutOff,NTot,
     $  ISc1,IScal,C,AtmChg,PTot,AllowP,XYZ,Coord,EFG,EFlux,GNorm,Phi,V,
     $  MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Compute electrostatic properties for SCIPCM solvation.
C
      Logical AllowP(*)
      Dimension C(3,NAtoms),AtmChg(NAtoms),PTot(*),JJ(1),XX(1),V(MDV)
      Dimension XYZ(NTot,3),Coord(3,NTot),EFG(NTot,10,Isc1),
     $  EFlux(NTot,2),GNorm(NTot,3),Phi(NTot)
      Save JJ, XX
      Data JJ/0/, XX/0/
C
      Do 10 J = 1, 3
        Do 10 I = 1, NTot
   10     Coord(J,I) = xyz(I,J)
C
C     IScal = 0 for no scaling, 2 for yes.
C     ISc1 = 2 if scaling, 1 if no.
C
C     If IScal=0
C         EFG(1,n,2) = Field stuff only from nuclear.
C         EFG(1,n,1) = Total field stuff.
C         EFlux(NTot,2) = -bnuc
C         EFlux(NTot,1) = Total Eflux
C         Phi = total potential
C
C     If IScal=2 as above but:
C         EFG(1,n,1) = Field stuff only from electronic part.
C         EFlux(Ntot,1) = Eflux for electrons
C
      Call AClear(ISc1*10*NTot,EFG)
      Call DmpEFG(NTot,xyz,NAtoms,C,AtmChg,EFG(1,1,ISc1))
      Do 20 I = 1, NTot
   20   EFlux(i,2) = EFG(i,4,ISc1)*GNorm(i,1) +
     $    EFG(i,3,ISc1)*GNorm(i,2) + EFG(i,2,ISc1)*GNorm(i,3)
      Call TStamp(1,'Before OneElI No. 1')
      Call OneElI(IOut,IPrint,1,6,202,0,0,0,2,1,NAtoms,JJ,PTot,PTot,
     $  NBas6D,1,JJ(1),JJ,C,.False.,JJ,XX,NAtoms,JJ,Coord,NTot,
     $  EFG(1,1,1),XX,CutOff,0,1,JJ,XX,IPFlag,AllowP,.False.,JJ,XX,V,V,
     $  MDV)
      Call TStamp(1,'After OneElI No. 1')
      If(IScal.eq.2) then
        Call AAdd(NTot,EFG(1,1,1),EFG(1,1,2),Phi)
      Else
        Call AMove(NTot,EFG(1,1,1),Phi)
        endIf
      Do 30 I = 1, NTot
   30   EFlux(i,1) = EFG(i,4,1)*GNorm(i,1) + EFG(i,3,1)*GNorm(i,2) +
     $    EFG(i,2,1)*GNorm(i,3)
      Return
      End
*Deck GFDen
      Subroutine GFDen(IOut,IPrint,N,M,NMO,ENR,H,R,D1,D2,RL,RLSq,RLTmp,
     $  LL,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
      Dimension H(NMO,NMO), R(NMO,NMO,NMO,NMO), D1(NMO,NMO),
     $  D2(NMO,NMO,NMO,NMO), RL(*), RLSq(*), RLTmp(*), LL(*), V(MDV)
 1000 Format(' GFDen:  D=',F12.8)
 1010 Format(' GFDen:  ',A,'=',G19.12)
 1020 Format(' Non-zero IJKL=',4I3,' D2IJKL=',1PD13.6,' D2JIKL=',
     $  1PD13.6,' D2KLIJ=',1PD13.6,' D2KLJI=',1PD13.6)
C
      Call AClear(NMO**2,D1)
      Call AMul(N,RL,RL,RLSq)
      D = TIJ(0,0,M,N,0,0,RLSq,RLTmp,LL)
      If(IPrint.ge.0) then
        Write(IOut,1000) D
        Call OutMtS(IOut,'RLSq:',0,0,RLSq,N,1,N,1)
        endIf
      Do 10 J = 1, N
   10   D1(J,J) = TI(0,M,N-1,M,N,J,0,RLSq,RLTmp,LL)/D
      Do 20 J = 1, M
        D1(J+N,J+N) = TIJ(0,1,M,N,J,0,RLSq,RLTmp,LL)/D
        D1(J,J+N) = RL(J)*TIJ(1,1,M,N,J,0,RLSq,RLTmp,LL)/D
   20   D1(J+N,J) = D1(J,J+N)
      If(IPrint.ge.0) Call OutMtS(IOut,'D1:',0,0,D1,NMO,NMO,NMO,NMO)
      E1 = SProd(NMO**2,D1,H)
      Write(IOut,1010) 'E1', E1
C
      Call AClear(NMO**4,D2)
      Do 40 I = 1, N
C Diaja,iaja
        Do 30 J = 1, (I-1)
          X = TI(0,M,N-2,M,N,I,J,RLSq,RLTmp,LL)/D
          D2(I,J,I,J) = X
          D2(J,I,I,J) = -X
          D2(J,I,J,I) = X
   30     D2(I,J,J,I) = -X
        Do 40 J = 1, M
          JB = J + N
C Diaja,iajb
          If(I.ne.J) then
            X = RL(J)*TI(1,M-1,N-2,M,N,I,J,RLSq,RLTmp,LL)/D
            D2(I,J,I,JB) = X
            D2(J,I,I,JB) = -X
            D2(I,JB,I,J) = X
            D2(I,JB,J,I) = -X
            endIf
C Diajb,iajb
          X = TIJ(0,1,M,N,I,J,RLSq,RLTmp,LL)/D
          D2(I,JB,I,JB) = X
C Djbia,iajb
          If(I.ne.J) then
            X = TIJ(1,1,M,N,I,J,RLSq,RLTmp,LL)/D
            D2(JB,I,I,JB) = X
            D2(I,JB,JB,I) = X
            endIf
   40     Continue
      Do 50 I = 1, M
        IB = I + N
        Do 50 J = 1, M
          JB = J + N
C Diaja,ibjb
          If(I.ne.J) then
            X = RL(I)*RL(J)*TIJ(2,2,M,N,I,J,RLSq,RLTmp,LL)/D
            D2(I,J,IB,JB) = X
            D2(I,J,JB,IB) = -X
            D2(J,I,IB,JB) = -X
            D2(IB,JB,I,J) = X
            D2(JB,IB,I,J) = -X
            D2(IB,JB,J,I) = -X
            D2(J,I,JB,IB) = X
            D2(JB,IB,J,I) = X
            endIf
C Djaib,iajb
          X = RL(I)*RL(J)*TIJ(1,2,M,N,I,J,RLSq,RLTmp,LL)/D
          D2(J,IB,I,JB) = X
          D2(I,JB,J,IB) = X
C Dibja,iajb
          If(I.ne.J) then
            X = RL(I)*RL(J)*TIJ(2,2,M,N,I,J,RLSq,RLTmp,LL)/D
            D2(IB,J,I,JB) = X
            D2(I,JB,IB,J) = X
C Dibja,ibjb
            X = RL(J)*TIJ(1,2,M,N,I,J,RLSq,RLTmp,LL)/D
            D2(IB,J,IB,JB) = X
            D2(IB,J,JB,IB) = -X
            D2(IB,JB,IB,J) = X
            D2(JB,IB,IB,J) = -X
C Dibjb,ibjb
            X = TIJ(0,2,M,N,I,J,RLSq,RLTmp,LL)/D
            D2(IB,JB,IB,JB) = X
            D2(IB,JB,JB,IB) = -X
            D2(JB,IB,IB,JB) = -X
            D2(JB,IB,JB,IB) = X
            endIf
   50     Continue
      E2 = SProd(NMO**4,D2,R)
      Write(IOut,1010) 'E2', E2
      Do 100 I = 1, N
        Do 100 J = 1, I
          Do 100 K = 1, NMO
            Do 100 L = 1, NMO
              If(D2(I,J,K,L).ne.(-D2(J,I,K,L)).or.
     $          D2(K,L,I,J).ne.(-D2(K,L,J,I))) Write(IOut,1020) I,J,K,L,
     $          D2(I,J,K,L), D2(J,I,K,L), D2(K,L,I,J), D2(K,L,J,I)
  100         Continue
      Do 110 I = (N+1), (N+M)
        Do 110 J = (N+1), I
          Do 110 K = 1, NMO
            Do 110 L = 1, NMO
              If(D2(I,J,K,L).ne.(-D2(J,I,K,L)).or.
     $          D2(K,L,I,J).ne.(-D2(K,L,J,I))) Write(IOut,1020) I,J,K,L,
     $          D2(I,J,K,L), D2(J,I,K,L), D2(K,L,I,J), D2(K,L,J,I)
  110         Continue
      Do 120 I = 1, NMO
        Do 120 J = 1, NMO
          Do 120 K = 1, NMO
            Do 120 L = 1, NMO
              If(D2(I,J,K,L).ne.D2(K,L,I,J)) Write(IOut,1020) I,J,K,L,
     $          D2(I,J,K,L), D2(J,I,K,L), D2(K,L,I,J), D2(K,L,J,I)
  120         Continue
      Return
      End
*Deck GFInts
      Subroutine GFInts(IOut,IPrint,NBas6D,NBasis,NAE,NBE,NMO,IPFlag,
     $  AllowP,CMOA,CMOB,CMOO,NTT6D,RAO,ENR,HMO,RMO,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Generate AO and occupied MO integrals for GF.
C
      Parameter (IRwGen=501,IRwH=515)
      Logical AllowP(*)
      Dimension CMOA(NBasis,*), CMOB(NBasis,*), CMOO(NBas6D,*),
     $  RAO(NTT6D,*), HMO(NMO,NMO), RMO(NMO,NMO,NMO,NMO), V(MDV),
     $  XX(1), JJA(1)
      Save Zero, One, XX, JJA
      Data Zero/0.0d0/, One/1.0d0/, XX/0.0d0/, JJA/0/
 1000 Format(' GFInts:  E1E=',G19.12,' E2E=',G19.12,' ETot=',G19.12)
C
      NTT = (NBasis*(NBasis+1))/2
      Call FileIO(2,-IRwGen,1,ENR,40)
      Call FileIO(2,-IRwH,NTT,HMO,0)
      IRaf = 0
      IOpClF = 0
      ICntrl = 600
      IPrtF = IPrint
      IPrtM1 = Max(IPrint-1,0)
      NBTX = 0
      NOpUse = 1
      Call FoFCou(IOut,IPrtF,IHMeth,IOpClF,ICntrl,IRaf,0,XX,.False.,0,
     $  0,0,IPFlag,AllowP,0,XX,0,0,.False.,.False.,.False.,.False.,
     $  .False.,0,Zero,One,0,1,1,0,1,0,0,NBasis,0,NOpUse,1,JJA,JJA,JJA,
     $  XX,JJA,NBTX,JJA,XX,XX,XX,XX,XX,XX,XX,XX,XX,JJA,XX,XX,JJA,XX,XX,
     $  0,JJA,XX,XX,JJA,JJA,XX,XX,JJA,XX,EJunk,One,RAO,0,0,XX,XX,JJA,
     $  JJA,XX,V,V,MDV)
      Call Square(RAO,RAO,NTT6D,NTT6D,0)
      Do 10 I = 1, NAE
        Call AMove(NBasis,CMOA(1,I),CMOO(1,I))
   10   Continue
      Do 20 I = 1, NBE
        Call AMove(NBasis,CMOB(1,I),CMOO(1,I+NAE))
   20   Continue
      IT = 1
      IV = IT + NBas6D**2
      MDV1 = MDV - IV + 1
      Call VFVMul(.False.,.False.,.False.,0,NBasis,NMO,NBasis,HMO,CMOO,
     $  V(IT),HMO,V(IV),MDV1)
      Call Square(HMO,HMO,NMO,NMO,0)
      If(IPrint.ge.0) Call OutMtS(IOut,'HMO:',0,0,HMO,NMO,NMO,NMO,NMO)
      Call RedoB(NBasis,V,IPureD,IPureF,IPrtM1)
      Call RedoC(IPrtM1,NBasis,NBas6D,NMO,CMOO,V)
      Do 30 IJ = 1, NTT6D
        Call VFVMul(.False.,.False.,.False.,0,NBas6D,NMO,NBas6D,
     $    RAO(1,IJ),CMOO,V(IT),RAO(1,IJ),V(IV),MDV1)
   30   Continue
      Call Trspn2(NTT6D,RAO)
      NMOTT = (NMO*(NMO+1))/2
      Do 40 IJ = 1, NMOTT
        Call VFVMul(.False.,.False.,.False.,0,NBas6D,NMO,NBas6D,
     $    RAO(1,IJ),CMOO,V(IT),RAO(1,IJ),V(IV),MDV1)
   40   Continue
      Do 50 I = 1, NMO
        Do 50 K = 1, I
          IK = (I*(I-1))/2 + K
          Do 50 J = 1, NMO
            JJ = (J*(J-1))/2
            Do 50 L = 1, J
              RMO(L,K,J,I) = RAO(JJ+L,IK)
              RMO(L,I,J,K) = RAO(JJ+L,IK)
              RMO(J,K,L,I) = RAO(JJ+L,IK)
   50         RMO(J,I,L,K) = RAO(JJ+L,IK)
      X1A = Zero
      X1B = Zero
      X2AA = Zero
      X2AB = Zero
      X2BB = Zero
      Do 70 I = 1, NAE
        X1A = X1A + HMO(I,I)
        Do 60 J = 1, (I-1)
   60     X2AA = X2AA + RMO(J,I,J,I) - RMO(J,I,I,J)
        Do 70 J = 1, NBE
   70     X2AB = X2AB + RMO(J+NAE,I,J+NAE,I)
      Do 80 I = 1, NBE
        X1B = X1B + HMO(I+NAE,I+NAE)
        Do 80 J = 1, (I-1)
   80     X2BB = X2BB + RMO(J+NAE,I+NAE,J+NAE,I+NAE)
     $      - RMO(J+NAE,I+NAE,I+NAE,J+NAE)
      X1 = X1A + X1B
      X2 = X2AA + X2AB + X2BB
      ETot = ENR + X1 + X2
      Write(IOut,1000) X1, X2, ETot
      Return
      End
*Deck GFTest
      Subroutine GFTest(IOut,IPrint,ITstGF,NBas6D,NBasis,NBsUse,NAE,NBE,
     $  IPFlag,AllowP,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Convert to biorthogonal orbitals and compute the projected HF
C     energy using the GF energy expression.
C
      Parameter (IRwMOA=524,IRwMOB=526,IRwUR=540)
      Dimension V(MDV)
      Logical AllowP(*)
C
      NBNU = NBasis*NBsUse
      NBNU6 = NBas6D*NBsUse
      NBSq = NBas6D**2
      LS = NBasis*Max(NBasis,5*NAE)
      NTT = (NBasis*(NBasis+1))/2
      NTT6D = (NBas6D*(NBas6D+1))/2
      NMO = NAE + NBE
      IRAO   = 1
      IHMO   = IRAO + NTT6D**2
      IRMO   = IHMO + Max(NMO**2,NTT)
      IRL    = IRMO + NMO**4
      ICMOA  = IRL + NAE
      ICMOB  = ICMOA + NBNU6
      IUA    = ICMOB + NBNU6
      IUB    = IUA + NAE**2
      ICMOAC = IUB + NBE**2
      ICMOBC = ICMOAC + NBSq
      ISABC  = ICMOBC + NBSq
      IT     = ISABC + NBSq
      IS     = IT + NBSq
      ISST   = IS + LS
      IV     = Max(ISST+NAE**2,IT+NMO*NBas6D)
      MDV1   = MDV - IV + 1
      Call TstCor(IV,MDV,'GFTest')
      Call FileIO(2,-IRwMOA,NBNU,V(ICMOAC),0)
      Call FileIO(2,-IRwMOB,NBNU,V(ICMOBC),0)
      Call FileIO(2,-IRwUR,NBsUse**2,V(ISABC),0)
      If(IPrint.ge.0) Call OutMtS(IOut,'SAB-canonical:',0,0,V(ISABC),
     $  NBsUse,NBsUse,NBsUse,NBsUse)
      Call BiOrMX(IOut,IPrint,1,NBasis,NBsUse,NAE,NBE,V(ICMOAC),
     $  V(ISABC),V(ICMOA),V(ISST),V(IUA),LS/NAE,V(IS),V(IT))
      Call BiOrMX(IOut,IPrint,2,NBasis,NBsUse,NBE,NAE,V(ICMOBC),
     $  V(ISABC),V(ICMOB),V(ISST),V(IUB),LS/NBE,V(IS),V(IT))
C
      Call CkBiMO(IOut,IPrint,NBasis,NBsUse,NAE,NBE,V(ICMOA),V(ICMOB),
     $  V(IRL),V(IUB),V(IS),V(IT))
      If(IPrint.ge.0) Call OutMtS(IOut,'Lambda:',0,0,V(IRL),NAE,1,NAE,1)
      If(ITstGF.eq.2) then
        Call FileIO(1,-IRwMOA,NBNU,V(ICMOA),0)
        Call FileIO(1,-IRwMOB,NBNU,V(ICMOB),0)
        endIf
C
      Call GFInts(IOut,IPrint,NBas6D,NBasis,NAE,NBE,NMO,IPFlag,AllowP,
     $  V(ICMOA),V(ICMOB),V(IT),NTT6D,V(IRAO),ENR,V(IHMO),V(IRMO),
     $  V(IV),MDV1)
C
      ID1 = IT
      ID2 = IT + NMO**2
      IRLSq = ID2 + NMO**4
      IRLT = IRLSq + NAE
      ILL = IRLT + NAE
      IV = ILL + InToWP(NMO+1)
      MDV1 = MDV - IV + 1
      Call TstCor(IV,MDV,'GFTest-GFDen')
      Call GFDen(IOut,IPrint,NAE,NBE,NMO,ENR,V(IHMO),V(IRMO),V(ID1),
     $  V(ID2),V(IRL),V(IRLSq),V(IRLT),V(ILL),V(IV),MDV1)
      Return
      End
*Deck MkFkSp
      Subroutine MkFkSp(IOut,IPrint,IPrtF,IOpCl,SpFock,SpDiag,ZDO,JustJ,
     $  IRwH,IRwCA,IOCB,IRwPA,IRwPB,NAtoms,NBasis,NBsUse,NBas6D,NAE,NBE,
     $  ICharg,Multip,IAn,IAtTyp,AtmChg,C,IAtFrg,NMtPBC,IPBC,LIPBC,
     $  CelVec,RHF,ROHF,HaveB,HavBRI,DirSCF,FrcNuc,IHMeth,IDenFit,
     $  IODFit,ICntrl,IRaf,KAlg,IHFJ,Omega,FMM,FMFlag,FMFlg1,NFxFlg,
     $  IPFlag,AllowP,LSEAll,RSEAll,AccDes,AccXC,DFT,ICtDFT,IExCor,
     $  IRadAn,IRanWt,IRanGd,ICorTp,ScaHFX,ScaDFX,IGWInf,IAtBtD,IRdBtD,
     $  RRdBtD,RGWBtD,DFTSym,JSym2X,NOpUse,NOp1,NOp2,NEqAtm,NEqShl,
     $  NEqBas,RotOp,NEqSh2,NSAOBk,NS,NS2,NS4,ISABF,RSABF,MulTab,MapSAO,
     $  DCut,EOne,JGC,ETwo,Ex,Ec,ECDFit,RNEFit,IXScr,IPStrt,IFStrt,JIP,
     $  IPA,IPB,JIFA,JIFB,IFA,IFB,IAftrF,ZThrP,NZP,NZF,NBTI,DA,DB,FA,FB,
     $  InStor,R1,R2,R3,Buf2E,IJ,IJMap,X,LenX)
      Implicit Real*8(A-H,O-Z)
C
C     Form a sparse Fock matrix, or form a full Fock matrix and convert
C     to sparse.
C
      Logical JustJ, SpFock, SpDiag, DFT, Mixed, HavBRI, CBfn, Spinor,
     $  FMM, AllowP(*), MixedK, HaveB, RHF, ZDO, ROHF, DFTSym, DirSCF,
     $  FrcNuc
      Integer FMFlag, FMFlg1
      Dimension IAn(*), IAtTyp(*), C(*), IPBC(*), CelVec(*), Omega(*),
     $  RSEAll(*), ScaDFX(*), IGWInf(*), IAtBtD(*), IRdBtD(*),
     $  RRdBtD(*), RGWBtD(*), NS(*), NS2(*), NS4(*), ISABF(*), RSABF(*),
     $  MulTab(*), MapSAO(*), NEqShl(*), NEqBas(*), RotOp(*), DA(*),
     $  DB(*), FA(*), FB(*), R1(*), R2(*), R3(*), X(*), AtmChg(*),
     $  IJ(*), Buf2E(*), IJMap(*), NEqAtm(*), NEqSh2(*), ECDFit(*),
     $  RNEFit(*), IAtFrg(*), JJ(1), XX(1)
      Character*1 JGC
      Save Zero, Pt5, JJ, XX, SclInc
      Data Zero/0.0d0/, Pt5/0.5d0/, JJ/0/, XX/0.0d0/, SclInc/1.0d0/
 1260 Format(' 1/2 <P',A,'(P)>=',F25.12)
C
      Call DecoSC(IOpCl,NSpBlk,NRI,NDimBl,CBfn,NSpBlX,Spinor)
      Mixed = MixedK(IHMeth,DFT,ScaHFX)
      If(DFT.and..not.Mixed) then
        JGC = 'J'
        IOpCl1 = 0
        ICntr1 = 1000*(ICntrl/1000) + 500 + Mod(ICntrl,100)
        If(HavBRI.and.SpDiag) then
          Call ACSASB(NZP+NBasis,X(ipa),Pt5,X(ipb),Pt5,X(ipa))
          Call AMove(NZP+NBasis,X(ipa),X(ipb))
          endIf
      else
        JGC = 'G'
        ICntr1 = ICntrl
        IOpCl1 = NSpBlX - 1
        endIf
      If(NBasis.ne.NBsUse)
     $  Call GauErr('MkFkSp is stupid about lin. deps.')
      NTT = (NBasis*(NBasis+1))/2
      NTT6D = (NBas6D*(NBas6D+1))/2
      If(SpFock) then
        NMat1 = 1
        NMatS1 = 1
        NMatT1 = 0
        If(.not.SpDiag) then
          jip = 1
          ipa = jip + InToWP(NBasis+1+NTT)
          ipb = ipa
          IFStrt = ipb + NTT + NBasis
          Call CmprLt(DA,X(jip),X(ipa),Zero,NBasis,NZP)
          endIf
        If(ZDO) then
          Call SpGet(X,X,LenX,IRWH,IFStrt,jifa,ifa,IAftrF,NBasis,NZF,0,
     $      .true.,0)
          ifb = ifa + (NZF+NBasis)*IOpCl1
          IAftrF = ifb + NZF + NBasis
          IV = IAftrF
          MDV1 = LenX - IV + 1
          Call AClear((NZF+NBasis),X(ifa))
          Call SpFoSE(IOut,IPrtF,IHMeth,IOpCl1,ICntr1,NAtoms,NBasis,IAn,
     $      IAtTyp,C,DCut,X(jip),X(ipa),X(ipb),X(jifa),X(ifa),X(ifb),
     $      NZF,X(IV),X(IV),MDV1)
        else
          jifa = IFStrt
          ifa = jifa + InToWP(NBasis+1+NBTI)
          ifb = ifa + (NBTI+NBasis)*IOpCl1
          IAftrF = ifb + NBTI + NBasis
          IV = IAftrF
          MDV1 = LenX - IV + 1
          Call FoFCou(IOut,IPrtF,IHMeth,IOpCl1,ICntr1,IRaf,IHFJ,Omega,
     $      FMM,FMFlag,FMFlg1,NFxFlg,IPFlag,AllowP,LSEAll,RSEAll,1,
     $      NAtoms,.True.,.True.,.True.,.True.,.True.,0,AccDes,ScaHFX,
     $      0,NMat1,NMatS1,NMatT1,NMat1,0,0,NBasis,JSym2X,NOpUse,NOp1,
     $      JJ,NEqShl,JJ,RotOp,NEqBas,NBTI,X(jip),X(ipa),X(ipb),XX,
     $      XX,XX,XX,XX,XX,XX,X(jifa),X(ifa),X(ifb),JJ,XX,XX,NAtoms,
     $      IAn,AtmChg,C,IAtTyp,IAtFrg,XX,XX,IPBC,CelVec,CoulEn,SclInc,
     $      XX,0,0,XX,XX,JJ,JJ,XX,X(IV),X(IV),MDV1)
          endIf
        If(.not.SpDiag) then
          Call UCmprs(FA,X(jifa),X(ifa),NBasis)
          If(IOpCl1.ne.0) Call UCmprs(FB,X(jifa),X(ifb),NBasis)
          endIf
      else
        iFockA = IFStrt
        iFockB = iFockA + NTT6D
        iDenA  = iFockB + NTT6D*IOpCl1
        iDenB  = iDenA + NTT6D
        iscr1  = iDenB + NTT6D*IOpCl1
        iscr2  = iscr1 + NTT6D
        IV = iscr2 + NTT6D
        MDV1 = LenX - IV + 1
        Call TstCor(IV,LenX,'FoF2e')
        Call Ucmprs(X(iDenA),X(jip),X(ipa),NBasis)
        If(IOpCl1.ne.0) Call Ucmprs(X(iDenB),X(jip),X(ipb),NBasis)
        Call FoF2E(IOut,IPrint,IOpCl1,ROHF,DFT,.True.,.True.,1,1,0,
     $    NMtPBC,NBasis,NBas6D,NAtoms,IAn,IAtTyp,AtmChg,C,IAtFrg,DirSCF,
     $    InStor,SclInc,FrcNuc,IDenFit,.True.,ICntrl,IHMeth,IRaf,ScaHFX,
     $    AccDes,FMM,FMFlag,FMFlg1,NFxFlg,IPFlag,AllowP,LSEAll,RSEAll,
     $    KAlg,Omega,NBTI,IODFit,IOFDen,JSym2X,NOp1,NOp2,NOpUse,NEqAtm,
     $    NEqBas,NEqShl,NEqSh2,RotOp,IJ,IJMap,NSAOBk,NS,NS2,NS4,ISABF,
     $    RSABF,MulTab,MapSAO,JGC,CoulEn,ECDFit,RNEFit,X(IDenA),
     $    X(IDenB),X(IFockA),X(IFockB),R1,R2,R3,IPBC,LIPBC,CelVec,Buf2E,
     $    X(IV),MDV1)
        If(IOpCl1.ne.0) then
          Call MkSprL(X,LenX,0,0,0,IFockB,jifb,ifb,IFATmp,NBasis,NZFB,
     $      'FoF2e',ZThrP,1)
          Call AMove(NTT,X(IFockA),X(IFATmp))
          Call SpAMvM(X,LenX,IFStrt,IFockA,jifb,ifb,NBasis,2)
          Call AMove(NTT,X(IFATmp),X(IFockA))
          endIf
        Call MkSprL(X,LenX,0,0,0,IFockA,jifa,ifa,IAftrF,NBasis,NZFA,
     $    'FoF2e',ZThrP,1)
        If(IOpCl1.eq.0) then
          jifb = jifa
          ifb = ifa
        else
          Call MkSmFm(X,LenX,0,jifa,ifa,0,jifb,ifb,IFStrt,IV,0,jifa,
     $      ifa,ifb,NZF,IAftrF,NBasis,2)
          endIf
        endIf
      If(SpDiag) then
        If(RHF) then
          Call SpRetr(IRWPA,X(jip),X(ipa),NBasis,NZP,1,.true.,0)
        else
          Call MkSmFm(X,LenX,IRwPA,0,0,IRwPB,0,0,IPStrt,IAftrF,
     $      0,jip,ipa,ipb,NZP,IFStrt,NBasis,2)
          endIf
        If(HaveB.and.IOpCl1.eq.0) then
          NZFB = NZFA
          jifb = jifa
          ifb  = ifa
          IStart = IAftrF
          Call SpAMvM(X,LenX,IStart,IAftrF,jifb,ifb,NBasis,2)
          endIf
        Call SpGet(X,X,LenX,IRWH,IAftrF,jiha,iha,iend,NBasis,NZH,
     $    1,.true.,0)
        istart = iend
        Call Compct(X(jiha),X(iha),ZThrP,NBasis,NZHA,1,.true.)
        Call SpAdaM(X,LenX,X(jifa),X(ifa),X(jiha),X(iha),X(ixscr),
     $    istart,jifha,ifha,iend,NBasis,NZFHA,Zero,2)
        If(RHF) then
          jifa = jifha
          ifa = ifha
          Call SpAMvM(X,LenX,IFStrt,IAftrF,jifa,ifa,NBasis,2)
          ifb = ifa
        else
          istart = iend
          Call SpAdaM(X,LenX,X(jifa),X(ifb),X(jiha),X(iha),X(ixscr),
     $      istart,jifhb,ifhb,iend,NBasis,NZFHA,Zero,2)
          Call MkSmFm(X,LenX,0,jifha,ifha,0,jifhb,ifhb,IFStrt,IEnd,0,
     $      jifa,ifa,ifb,NZF,IAftrF,NBasis,2)
          endIf
        Call Cmpct2(IOpCl,.True.,NBasis,ZThrP,X(jifa),IFA,IFB,NZFA,X)
        If(IPrint.ge.2) Call SpRprt(NZFA,NBasis*(NBasis-1)/2,ZThrP,
     $    'AOA Fock matrix')
        If(RHF) then
          ETwo = SpTrc2(X(jifa),X(ifa),X(jip),X(ipa),X(ixscr),NBasis)-
     $      Pt5*EOne
        else
          ETwo = Pt5*(SpTrc2(X(jifa),X(ifa),X(jip),X(ipa),X(ixscr),
     $      NBasis) + SpTrc2(X(jifa),X(ifb),X(jip),X(ipb),X(ixscr),
     $      NBasis) - EOne)
          endIf
        endIf
      If(DFT) then
        If(IPrint.gt.0) Write(IOut,1260) JGC, ETwo
        Call RIComp(.True.,.True.,.False.,.True.,IOpCl,NBasis,1,1,1,0,
     $    NMat,NMatS,NMatT,DA,DB,X,LenX)
        Call RIComp(.True.,.True.,.True.,.True.,IOpCl,NBasis,1,1,1,0,
     $    NMat,NMatS,NMatT,FA,FB,X,LenX)
        If(IPrint.ge.5) Call PrtFok(IOut,NSpBlX-1,'J/G matrix:',0,
     $    NMat,NMtPBC,NBasis,FA,FB,X,LenX)
        If(JustJ) Call GauErr('Stop after computing J.')
        If(SpFock) then
          IV = IAftrF
          MDV1 = LenX - IV + 1
          Call FoFDFT(IOut,IPrtF,IPFlag,IOpCl,HavBRI,ICtDFT,IExCor,
     $      IRadAn,IRanWt,IRanGd,ICorTp,ScaDFX,DFTSym,NOpUse,NOp1,1,
     $      NAtoms,NAE,NBE,NBasis,NBas6D,ICharg,Multip,IRwCA,0,IOCB,
     $      AccXC,IAn,IAtTyp,AtmChg,C,Ex,Ec,ETwo,.True.,IPBC,CelVec,
     $      X(jip),X(ipa),XX,XX,X(ipb),XX,XX,X(jifa),X(ifa),X(ifb),XX,
     $      XX,IGWInf,IAtBtD,IRdBtD,RRdBtD,RGWBtD,Omega,X(IV),MDV1)
        else
          ifocka = IAftrF
          ifockb = ifocka + NTT6D
          idena = ifockb + IOpCl*NTT6D
          idenb  = idena + NTT6D
          ibuf2e = idenb + IOpCl*NTT6D
          lenbuf = LenX - ibuf2e + 1
          Call TstCor(ibuf2e,LenX,'FoFDFT-SpDiag-NoSpFock')
          Call Ucmprs(X(ifocka),X(jifa),X(ifa),NBasis)
          Call Ucmprs(X(idena),X(jip),X(ipa),NBasis)
          If(.not.RHF) then
            Call Ucmprs(X(ifockb),X(jifa),X(ifb),NBasis)
            Call Ucmprs(X(idenb),X(jip),X(ipb),NBasis)
            endIf
          Call FoFDFT(IOut,IPrtF,IPFlag,IOpCl,HavBRI,ICtDFT,IExCor,
     $      IRadAn,IRanWt,IRanGd,ICorTp,ScaDFX,DFTSym,NOpUse,NOp1,1,
     $      NAtoms,NAE,NBE,NBasis,NBas6D,ICharg,Multip,0,0,0,AccXC,IAn,
     $      IAtTyp,AtmChg,C,Ex,Ec,ETwo,.False.,IPBC,CelVec,JJ,X(idena),
     $      XX,XX,X(idenb),XX,XX,JJ,X(ifocka),X(ifockb),XX,XX,IGWInf,
     $      IAtBtD,IRdBtD,RRdBtD,RGWBtD,Omega,X(ibuf2e),lenbuf)
          If(.not.RHF) then
            Call MkSprL(X,LenX,0,0,0,IFockB,jifb,ifb,IFATmp,NBasis,
     $        NZFB,'FoFDFT',ZThrP,1)
            Call AMove(NTT,X(IFockA),X(IFATmp))
            Call SpAMvM(X,LenX,IFStrt,IFockA,jifb,ifb,NBasis,2)
            Call AMove(NTT,X(IFATmp),X(IFockA))
            endIf
          Call MkSprL(X,LenX,0,0,0,IFockA,jifa,ifa,IAftrF,NBasis,NZFA,
     $      'FoFDFT',ZThrP,1)
          If(IOpCl.eq.0) then
            jifb = jifa
            ifb = ifa
          else
            Call MkSmFm(X,LenX,0,jifa,ifa,0,jifb,ifb,IFStrt,IBuf2E,0,
     $        jifa,ifa,ifb,NZF,IAftrF,NBasis,2)
            endIf
          endIf
        endIf
      Call TStamp(1,'After DFT XC')
      Call RIComp(.False.,.True.,.True.,.True.,IOpCl,NBasis,1,1,1,0,
     $  NMat,NMatS,NMatT,FA,FB,X,LenX)
      Return
      End
*Deck MkGamF
      Subroutine MkGamF(IOut,IPrint,Label,NMtPBC,NBasis,NTT,IOF0KA,
     $  IOF0KB,IRwFGA,IRwFGB,FA,FB)
      Implicit Real*8(A-H,O-Z)
C
C     Form the gamma-point Fock matrix.  FA and FB must be big enough
C     to hold all the F0k's.
C
      Character*(*) Label
      Dimension FA(NTT,NMtPBC), FB(NTT,NMtPBC)
 1000 Format(1X,A,' Fock matrix at Gamma point ',A,':')
C
      Call FileIO(2,-IOF0KA,NTT*NMtPBC,FA,0)
      Call SumCol(.False.,NTT,NTT,NMtPBC-1,FA(1,2),FA)
      Call FileIO(1,-IRwFGA,NTT,FA,0)
      If(IPrint.ge.3) then
        Write(IOut,1000) 'Alpha', Label
        Call LTOut(IOut,NBasis,FA,0)
        endIf
      If(IRwFGB.ne.0) then
        Call FileIO(2,-IOF0KB,NTT*NMtPBC,FB,0)
        Call SumCol(.False.,NTT,NTT,NMtPBC-1,FB(1,2),FB)
        Call FileIO(1,-IRwFGB,NTT,FB,0)
        If(IPrint.ge.3) then
          Write(IOut,1000) 'Beta', Label
          Call LTOut(IOut,NBasis,FB,0)
          endIf
        endIf
      Return
      End
*Deck Mtrx2
      Subroutine Mtrx2(CC,Ncore,Norbs,Iout,IEnd,jiz,iz,NZz,jiz2,iz2,
     $  NZz2,Thresh,Fixed,ImemOp)
      Implicit Real*8(A-H,O-Z)
C
C     Calculate the square (z2) of the input matrix z (lower triangle).
C
      Real*8 CC(Ncore)
      Logical DoLttA,DoLttb,DoLttc,Screen,Fixed
C
C     z**2 matrix
C
      Istart = IEnd
      jizc = jiz
      izc = iz
      Call SpAMvM(CC,Ncore,Istart,IEnd,jizc,izc,Norbs,2)
      DoLtta = .True.
      DoLttb = .True.
      DoLttc = .True.
      Screen = .False.
      imul = 1
      If(Fixed) then
        imode = 3
        jiz2 = IEnd
        Call IMove1(Norbs+1+NZz,0,(jiz2-jiz)*IntPWP(0),CC(jiz))
        iz2 = jiz2 + IntoWp(Norbs + 1 + NZz)
        IEnd = iz2 + Norbs + NZz
      else
        imode = 2
        jiz2 = IEnd
        iz2 = IEnd
        endIf
      Istart = IEnd
      Call FsMpyM(CC,Ncore,CC(jiz),CC(iz),CC(jizc),CC(izc),Istart,
     $  Istart,jiz2,iz2,IEnd,Thresh,Norbs,NZz2,CC(jiz2),CC(iz2),DoLttA,
     $  DoLttb,DoLttc,Screen,imul,imode,ImemOp)
      Istart = jizc
      Call SpAMvM(CC,Ncore,Istart,IEnd,jiz2,iz2,Norbs,2)
      Return
      End
*Deck NewP
      Subroutine NewP(IOut,IPrint,NP,NBasis,Ne,C,IOP,IOErr,IRwPO,A,
     $  NAMax)
      Implicit Real*8(A-H,O-Z)
      Dimension A(*),C(*),VDummy(1),JJ(1)
      Save Half,One,Conv,MaxItr,Reject,JJ
      Data Half/5.D-1/,One/1.D0/,Conv/1.D-12/,MaxItr/3/,Reject/1.D-1/,
     $  VDummy/0.0d0/,JJ/0/
C
C     Do some initialization
C
      NTT   = NBasis*(NBasis+1)/2
      NBSq  = NBasis*NBasis
      IScr1 = 1
      IScr2 = IScr1 + NBSq
      IV    = IScr2 + NBSq
      Call TstCor(IV+NBSq,NAMax,'NewP')
      dPScl = One
C
C     P' = Sum Ci*Pi
C     Pnew = P' + Sum Ci*dPi
C
    5 Call AClear(NTT,A(IScr1))
      Call AClear(NTT,A(IScr2))
      Do 10 I = 1, NP
        Call FileIO(2,-IOP,NTT,A(IScr2),(I-1)*NTT)
        Call ACASB(NTT,A(IScr1),A(IScr2),A(IScr1),C(I))
        Call FileIO(2,-IOErr,NTT,A(IScr2),(I-1)*NTT)
        Call ACASB(NTT,A(IScr1),A(IScr2),A(IScr1),dPScl*C(I))
  10    Continue
C
C     Do McWeeny purification to keep idempotency
C
      Do 20 I = 1, MaxItr
        Scal = GFloat(Ne)/TraLt(1,NBasis,A(IScr1))
        Call AScale(NTT,Scal,A(IScr1),A(IScr1))
        Call McWeny(IOut,IPrint,.True.,Junk,JJ,A(IScr1),A(IScr2),
     $    A(IV),NBasis,NTT,ADemp,1,.False.,VDummy,1)
        Call AMove(NTT,A(IScr2),A(IScr1))
        If(ADemp.le.Conv) goto 30
   20   Continue
   30 If(ADemp.ge.Reject.and.dPScl.ge.1.D-2) then
        dPScl = dPScl*Half
        Goto 5
        endIf
      Call FileIO(1,-IRwPO,NTT,A(IScr1),0)
      Return
      End
*Deck NRDIS
      Subroutine NRDIS(IOut,IPrint,P,G,H,X,BMat,MaxD,NSaved,N,NZ,IOP,
     $  IOG,IS,IAD1,IAD2,D,Q,R,QR)
      Implicit Real*8(A-H,O-Z)
C
C     Fixed Form DIIS for sparse matrices
C     DIIS extrapolation for SCF.  See
C     Pulay, J. Comp. Chem. Vol. 3 No. 4 556-560 (1982).
C
C     P      ... Density matrix (On exit containts the extrapolated matrix)
C     G      ... Gradient
C     H      ... Approximate Hessian
C     X      ... Scratch array with the size of N+NZ
C     BMat   ... Scratch array to store the DIIS error matrix
C     MaxD   ... Maximum number of DIIS vectors to be used
C     NSaved ... Number of current iteration, -1 to restart.
C     N      ... Dimension of density and error matrices
C     NZ     ... Length of density and error matrices
C     IOP    ... Storage for density matrices
C     IOG    ... Storage for gradient matrices
C     IOErr  ... Storage for error matrices
C     IPrint ... Print flag:
C                0 ... print nothing.
C                1 ... print a few messages.
C                2 ... print B-matrix.
C                3 ... print everything
C
      Dimension P(*),G(*),H(*),X(*),BMat(*),IS(*),IAD1(*),IAD2(*),D(*),
     $  Q(*),R(*),QR(*)
      Logical Inv2
      Save One
      Data One/1.0D0/
 1000 Format(' DIIS Inversion failure: switch to NR Step')
 1010 Format(' Coeff:     ',6D12.5)
C
C     If this is the first iteration, initialize.  niter=-1 means
C     restart.
C
      Len = N + NZ
      Leng = Len * MaxD
      If(NSaved.le.1) then
        If(IOP.eq.0) IOP=IGetRW(Leng,0)
        Call ConDDF(IOP,Leng)
        If(IOG.eq.0) IOG=IGetRW(Leng,0)
        Call ConDDF(IOG,Leng)
      EndIf
C
C     Store the matrices.
C
      Call FileIO(1,-IOP,Len,P,(NSaved-1)*Len)
      Call FileIO(1,-IOG,Len,G,(NSaved-1)*Len)
C
C     If this is the first iteration, do simple Newton-Raphson step
C
      If(NSaved.le.1) then
        Call ADiv(Len,G,H,G)
        Call AAdd(Len,P,G,P)
        Return
      endif
C
C     Build B-matrix
C
      IOErr = IGetRW(Leng,0)
      Do 10 I = 1, NSaved
        Call FileIO(2,-IOG,Len,X,(I-1)*Len)
        Call ADiv(Len,X,H,X)
        Call FileIO(1,-IOErr,Len,X,(I-1)*Len)
 10   Continue
      NBDim = NSaved + 1
      LBMat = MaxD*MaxD
      Call AClear(LBMat,BMat)
      Call FrmBMS(NSaved,N,Len,NZ,BMat,MaxD,IOErr,X,H)
C
C     Form and invert Pulay's B matrix, using a routine which ensures
C     numerical stability even if the matrix is quite ill-conditioned.
C
      Call Pcck(0,KK,BMat,MaxD,MaxD,NBDim,NBDim)
      If(IPrint.ge.2) Call OutMat(IOut,1,BMat,NBDim,NBDim,NBDim,NBDim)
      If(.not.Inv2(IOut,1,BMat,NBDim,IS,IAD1,IAD2,D,Det,Q,R,QR,0))then
        If(IPrint.ge.2) Write(IOut,1000)
        Return
      Endif
      Call AClear(NBDim,D)
      D(1) = -One
      Call MatMpy(1,NBDim,NBDim,D,BMat,Q)
      If(IPrint.ge.2) Write(IOut,1010) (Q(I+1),I=1,NSaved)
C
C     Compute the extrapolated matrix.
C
      Call AClear(Len,P)
      Do 20 I = 1, NSaved
        Call FileIO(2,-IOP,Len,H,(I-1)*Len)
        Call FileIO(2,-IOErr,Len,X,(I-1)*Len)
        Call ACASB(Len,P,H,P,Q(I+1))
        Call ACASB(Len,P,X,P,Q(I+1))
 20   Continue
      Call FreeRW(IOErr)
      Return
      End
*Deck NSIter
      Subroutine NSIter(IOut,Iprint,CC,Ncore,Norbs,jif,if,NZf,IEnd,
     $  Thresh,jipa,ipa,NZpa,ConvCr,ItMax,Niter,NElect,Errz,OuterF,
     $  Fixedo,FixTyp,CNSI,ImemOp)
      Implicit Real*8(A-H,O-Z)
C
C     Newton-Schultz iteration algorithm
C     CG-NSI or C-NSI
C     for the calculation of sign ( -> density) matrix
C     related to the input fock matrix
C
      Dimension CC(Ncore)
      Logical DoLttA,DoLttb,DoLttc,Screen,Fixed,pCrit,CNSI
      Logical FiniIt,Quadr,Fixedo,FixTyp,OuterF
      Save Zero,One,Four,Half,Hundred
      Data Zero/0.d0/,One/1.d0/,Four/4.d0/,Half/0.5d0/,Hundred/100.d0/
  788 Format(' pn = ',f12.5)
  789 Format(' pn extremal, finish iteration')
  899 Format(' Convergence achieved in',i5,' steps')
  911 Format(' Form of matrices not accurate enough')
  912 Format(' Finishing iteration, step = ',i5)
 4000 Format(' step= ',i2,' idemp= ',E9.1,' size= ',i6,' SPz= ',
     $  f4.0,' off= ',E9.1,' NZz= ',i10)
 5000 Format(' step= ',i2,' idemp= ',E9.1,' size= ',i6,' SPz= ',
     $  f4.0,' SPz2= ',f4.0,' off= ',E9.1,' NZz= ',i10,' NZz2= ',i10)
c
      IEndOl = jif
C
C     Initial guess for sign (z0)
C
      XlgO = Zero
      Errzo = One
      jiz = jif
      iz  = if
      NZz = NZf
C
C     First renormalization
C
      If(CNSI) then
        IX = IEnd
        IEnd = IX + Norbs
        Call RenZ(CC,Ncore,Norbs,IOut,IEnd,Thresh,CC(IX),jiz,iz,NZz,
     $    NElect)
        endIf
C
C     Iterate for sign(f)
C
      FiniIt = .False.
      Quadr = .False.
      OffMaO = One
      IFin = 0
      Do 100 II = 1, ItMax
        Fixed = Fixedo.and.(II.ne.1.or.FixTyp.or.OuterF)
        IX = IEnd
        IEnd = IX + Norbs
C
C       z**2 matrix
C
        Call Mtrx2(CC,Ncore,Norbs,IOut,IEnd,jiz,iz,NZz,jiz2,iz2,NZz2,
     $    Thresh,Fixed,ImemOp)
C
C       Bounds of z2
C
        Call GBnds(CC(jiz2),CC(iz2),Norbs,NZz2,z2abs,z2u,z2l,IOut,
     $    CC(IX))
C
C       Normalize z and z2
C
        If(.not.CNSI) Call RenZ4(CC,Ncore,Norbs,IOut,IEnd,Thresh,II,
     $    CC(IX),jiz,iz,NZz,jiz2,iz2,NZz2,NElect,z2u,z2l)
        SPz = GFloat(2*NZz+Norbs)/gfloat(Norbs*Norbs)*Hundred
        SPz2 = GFloat(2*NZz2+Norbs)/gfloat(Norbs*Norbs)*Hundred
C
C       Error of z2 matrix (z^2 - i)
C
        Errz = Zero
        Do 566 i = 0, (Norbs-1)
          x = One - CC(iz2+i)
  566     Errz = Errz + x*x
        kl = iz2 + Norbs
        OffMax = Zero
        Do 567 i = 0, (NZz2-1)
          x = Abs(CC(kl+i))
          Errz = Errz + x*x
          If(x.gt.OffMax) OffMax = x
  567     Continue
C       Measure idemp.
        Errz = Sqrt(Errz)/gfloat(Norbs)/Four
C
C       Build (1-z**2)
C
        Do 435 i = 0, (Norbs-1)
  435     CC(iz2+i) = CC(iz2+i) - One
        Call aneg(Norbs + NZz2,CC(iz2),CC(iz2))
C
C       Determine Fixed Format and change Format of z and 1 - z2
C
        If(.not.OuterF.and..not.FixTyp.and.Fixedo.and.II.eq.1) then
          IX = IEnd
          IEnd = IX + Norbs
          Call sform3(CC,Ncore,Norbs,IEnd,IOut,Thresh,jipa,ipa,NZpa,jiz,
     $      iz,NZz,jiz2,iz2,NZz2,ix,ImemOp)
          Fixed = .True.
          endIf
C
C       Dz matrix
C
C       (1 - z**2)*z
C
        DoLtta = .True.
        DoLttb = .True.
        DoLttc = .True.
        Screen = .False.
        IMul = 1
        If(Fixed) then
          IMode = 3
          jidz = IEnd
          Call IMove1(Norbs+1+NZz,0,(jidz-jiz)*IntPWP(0),CC(jiz))
          idz = jidz + IntoWp(Norbs + 1 + NZz)
          IEnd = idz + Norbs + NZz
        else
          IMode = 2
          jidz = IEnd
          idz = IEnd
          endIf
        Istart = IEnd
        Call FsMpyM(CC,Ncore,CC(jiz2),CC(iz2),CC(jiz),CC(iz),Istart,
     $    Istart,jidz,idz,IEnd,Thresh,Norbs,NZdz,CC(jidz),CC(idz),
     $    DoLttA,DoLttb,DoLttc,Screen,IMul,IMode,ImemOp)
C
C       Update z matrix
C
        If(CNSI) then
          Call ZUpdt2(CC,Ncore,Norbs,Thresh,IEnd,IOut,pn,pCrit,jiz,iz,
     $      NZz,jiz2,iz2,NZz2,jidz,idz,NZdz,CC(IX),Fixed)
          If(pCrit) then
            If(Iprint.ge.1) then
              Write(IOut,788) pn
              Write(IOut,789)
              endIf
            Goto 200
            endIf
        else if(Fixed) then
          Call ACASB(Norbs+NZdz,CC(iz),CC(idz),CC(iz),Half)
        else
          Istart = IEnd
          itype = 2
          Call SpAddM(CC,Ncore,CC(jiz),CC(iz),CC(jidz),CC(idz),CC(IX),
     $      Istart,jizn,izn,IEnd,Norbs,NZzn,One,Half,Thresh,itype)
          endIf
C
C       Move back IEnd (update z matrix)
C
        If(CNSI) then
          If(Fixed) then
            IEnd = IX
          else
            Istart = IEndOl
            Call SpAMvM(CC,Ncore,Istart,IEnd,jiz,iz,Norbs,2)
            endIf
        else if(Fixed) then
          IEnd = IX
        else
          Istart = IEndOl
          jiz = jizn
          iz = izn
          NZz = NZzn
          Call SpAMvM(CC,Ncore,Istart,IEnd,jiz,iz,Norbs,2)
          endIf
C
C       Report state of convergence
C
        If(IPrint.ge.1) then
          If(Fixed) then
            Write(IOut,4000) II,Errz,Norbs,SPz,OffMax,NZz
          else
            Write(IOut,5000) II,Errz,Norbs,SPz,SPz2,OffMax,NZz,NZz2
            endIf
          endIf
C
C       Check normal termination of NSI
C
        If(OffMax.le.ConvCr.and.OffMax/OffMaO.gt.Half
     $    .or.OffMax.lt.Thresh) then
          If(IPrint.ge.1) Write(IOut,899) II
          Goto 200
          endIf
C
C       Fixed form satisfactory?
C
        Xlg = GLog10(Errz)
        Quadr = Quadr.or.(II.gt.1.and.(XlgO-Xlg).gt.One)
        FiniIt = FiniIt.or.(Quadr.and.Errz/Errzo.ge.Half)
        If(II.ge.5.and.Errz.gt.Errzo) IFin = IFin + 1
        If(II.ge.5.and.Errz.lt.Errzo) IFin = 0
        FiniIt = FiniIt.or.IFin.eq.3.or.(II.ge.20.and.Errz.gt.Errzo)
        If(Fixed.and.FiniIt) then
          If(IPrint.ge.1) then
            Write(IOut,911)
            Write(IOut,912) II
            endIf
          Goto 200
          endIf
        XlgO = Xlg
        Errzo = Errz
        OffMaO = OffMax
  100   Continue
      Call Gauerr('Newton-Schultz iteration not converged in NSIter.')
C
C     Density mat
C
  200 jipa = jiz
      ipa = iz
      NZpa = NZz
      Call ascale(Norbs+NZz,-Half,CC(ipa),CC(ipa))
      Do 436 i = 0, (Norbs-1)
  436   CC(ipa+i) = CC(ipa+i) + Half
      Return
      End
*Deck NwDISF
      Subroutine NwDISF(IOA,IOB,IA,A,IB,B,AD,X,N,NZA,NZB,IType,NMat,
     $  IUHF)
      Implicit Real*8 (A-H,O-Z)
      Dimension IA(N+1,2),A(NZA,IUHF),IB(N+1,2),B(NZB,IUHF),X(N),AD(*)
C
C  Change the form a NMat matrices in IOA from A to B
C
C  iType = 1 No Diagonal
C          2 Diagonal
C
      LenA = (NZA + (itype-1)*N)*IUHF
      LenB = (NZB + (itype-1)*N)*IUHF
      Do 100 imat=1, NMat
        call fileio(2,-IOA,NZA*IUHF,A,(IMat-1)*(LenA))
        If(iType.eq.2) call fileio(2,IOA,N*IUHF,AD,0)
        Call NwForm(IA,A,IB,B,X,N)
        If(IUHF.eq.2) Call NwForm(IA,A(1,2),IB,B(1,2),X,N)
        call fileio(1,-IOB,NZB*IUHF,B,(IMat-1)*(LenB))
  100   If(iType.eq.2) call fileio(1,IOB,N*IUHF,AD,0)
      Return
      End
*Deck P1Chek
      Logical Function P1Chek(JCycle,Energy,DelEn,Acurcy,LenP,NDimBl,
     $  NSpBlk,NBNO,IRwPA,IRwPB,IRwCA,IRwCB,IOCMOS,IRwOld,NumUp,MaxNU,V,
     $  MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Check if the current iteration of the first (low-accuracy) pass of
C     the SCF has gone awry.  This can happen for nearly linearly dependant
C     cases, for which the cheap integrals are just not accurate enough.
C     This routines .true. if there ARE problems and corrective action
C     must be taken.  This routine also takes care of stashing away the
C     density (if the step is found to be OK) and restoring the previous
C     density (if the step has failed).
C
      Dimension V(*)
      Save EUpMax, RMSLim
      Data EUpMax/0.1d0/, RMSLim/2.01d0/
C
      If(NumUp.ge.MaxNU) goto 900
      If(JCycle.gt.1) then
C
C       Check if the energy has gone up by too much, or down by a
C       suspicious amount.
C
        If(JCycle.gt.2) then
          If(DelEn.ge.EUpMax*Abs(Energy)) goto 900
          If(DelEn.le.(-Abs(Energy))) goto 900
          endIf
C
C       Check if the RMS density change is suspcious.
C
        If(Acurcy.ge.RMSLim) goto 900
        endIf
C
C     The step looks OK -- save the densities.
C
      Call PCOld(1,LenP,NDimBl,NSpBlk,NBNO,IRwPA,IRwPB,IRwCA,IRwCB,
     $  IOCMOS,IRwOld,V,MDV)
      P1Chek = .False.
      Return
C
C     The step failed -- use the old densities.
C
  900 P1Chek = .True.
      Call PCOld(2,LenP,NDimBl,NSpBlk,NBNO,IRwPA,IRwPB,IRwCA,IRwCB,
     $  IOCMOS,IRwOld,V,MDV)
      Return
      End
*Deck PCMDAl
      Subroutine PCMDAl(IOpCl,NBasis,NBsUse,NumFS,MxSCFC,NVQ,LenY,jD,jV)
      Implicit Real*8 (A-H,O-Z)
C
C     Little routine to re-allocate D and V in Y before calling SCFDIS
C     when we're extrapolating also the PCM charges. This will no longer
C     be needed when the proper allocation for NVQ will be accounted for
C     at the top of UHFOpn.
C
      Logical CBfn,Spinor
C
      Call DecoSC(IOpCl,NSpBlk,NRI,NDimBl,CBfn,NSpBlX,Spinor)
C
      NBD = NBasis*NDimBl
      NBUD = NBsUse*NDimBl
      NTT = (NBD*(NBD+1))/2
      NTTU = (NBUD*(NBUD+1))/2
      NBDSq = NBD*NBD
      NTTR = NRI*NTT
      NTTUR = NRI*NTTU
C
C     Set the length of D.
      LenD = Max((NumFS*NTTR+2*NVQ),(NSpBlk*NTTUR+NVQ))
      LenD = Max(LenD,Max(NBDSq,(((MxSCFC+1)*(MxSCFC+2))/2)))
C
C     Set the length of V.
      LenV = NSpBlk*NTTUR + 2*NVQ
      LenV = Max(LenV,Max(NBDSq,((MxSCFC+1)**2)))
      LenV = Max(LenV,LenD)
C
C     Set the offset in Y because Y starts where BMat starts.
      LenB = (MxSCFC+1)**2
C
C     Set the pointers in Y.
      jD    = 1  + LenB
      jV    = jD + LenD
      jEnd  = jV + LenV - 1
      Call TstCor(jEnd,LenY,'PCMDAl')
      Return
      End
*Deck PCMDIS
      Subroutine PCMDIS(IOut,IPrint,SetUp,NumFS,NBasis,NRI,NDimBl,NVQ,
     $  IRwErQ,IRwF1,IRwF2,IRwP1X,IRwP2X,IRwF1U,IRwF2U,IRwP1U,IRwP2U,
     $  WFe,WQe,EWFe,EWQe,V,MDV)
      Implicit Real*8 (A-H,O-Z)
C
C     Setup or restore files before or after SCFDIS when extrapolating
C     also the PCM charges.
C
      Logical SetUp,SPCM,YPCM
      Real*8 V(*)
      Save Zero,RJunk,IJunk
      Data Zero/0.0d0/,RJunk/0.0d0/,IJunk/0/
 1000 Format(' PCMDIS: SetUp=',l1,' S,Y=',2l2,' NVQ=',i8)
 1010 Format(' PCMDIS: WFe=',f8.3,' WQe=',f8.3,' EWFe=',f8.3,' EWQe=',
     $  f8.3)
C
      Call PCMGet('IModel',IModel,RJunk)
      Call PCMGet('IEFTyp',IEFTyp,RJunk)
      SPCM = IModel.eq.2
      YPCM = IModel.eq.3.and.(IEFTyp.ge.4.and.IEFTyp.le.6)
C
      If(SetUp) then
        NVQ = 0
        Call PCMGet('VarPCM',IVarG,RJunk)
        Call PCMGet('DIISQ',IDIISQ,RJunk)
        If(IVarG.eq.1.and.IDIISQ.eq.1) Call PCMGet('NVQ',NVQ,RJunk)
        Call PCMGet('WFe',IJunk,WFe)
        Call PCMGet('WQe',IJunk,WQe)
        Call PCMGet('EWFe',IJunk,EWFe)
        Call PCMGet('EWQe',IJunk,EWQe)
        endIf
      If(IPrint.ge.2) then
        Write(IOut,1000) SetUp,SPCM,YPCM,NVQ
        If(SetUp) Write(IOut,1010) WFe,WQe,EWFe,EWQe
        endIf
C     Quick return.
      If(NVQ.eq.0) Return
      If(.not.(SPCM.or.YPCM))
     $  Call GauErr('Only SPCM or YPCM in PCMDIS.')
C
      NBD = NBasis*NDimBl
      NTT = (NBD*(NBD+1))/2
      NTTR = NRI*NTT
      LenPU = NumFS*NTTR + NVQ
      LenFU = NumFS*NTTR + NVQ
      If(SetUp) then
C       Add dE/dQ after F in IRwF1, and Q after P in IRwP1X.
        LenF = ItQry(IRwF1)
        LenP = ItQry(IRwP1X)
        If((LenFU-LenF).ne.NVQ) Call GauErr('Wrong LenF in PCMDIS.')
        If((LenPU-LenP).ne.NVQ) Call GauErr('Wrong LenP in PCMDIS.')
        IRwF1U = IGetRW(LenFU,0)
        IRwP1U = IGetRW(LenPU,0)
        Call CpyFi2(IRwF1,0,IRwF1U,0,LenF,V,MDV)
        Call CpyFi2(IRwP1X,0,IRwP1U,0,LenP,V,MDV)
        Call TstCor(NVQ,MDV,'PCMDIS-0')
        Call PCMFIO('read','EdQbar',V,NVQ,0,0)
        Call FileIO(1,-IRwF1U,NVQ,V,LenF)
        If(SPCM) then
          Call PCMFIO('read','Qtot',V,NVQ,0,0)
        else
          Call PCMFIO('read','Qttot',V,NVQ,0,0)
          endIf
        Call FileIO(1,-IRwP1U,NVQ,V,LenP)
        If(IRwF2.ne.0) then
          IRwF2U = IGetRW(LenFU,0)
          IRwP2U = IGetRW(LenPU,0)
          Call IniBuc(IRwF2U,-1,Zero,V,MDV)
          Call IniBuc(IRwP2U,-1,Zero,V,MDV)
          Call CpyFi2(IRwF2,0,IRwF2U,0,LenF,V,MDV)
          Call CpyFi2(IRwP2X,0,IRwP2U,0,LenP,V,MDV)
        else
          IRwF2U = 0
          IRwP2U = 0
          endIf
C       Write EdQbar as error vector for the Q in IRwErQ.
        IRwErQ = IGetRW(NVQ,0)
        Call PCMFIO('read','EdQbar',V,NVQ,0,0)
        Call FileIO(1,-IRwErQ,NVQ,V,0)
      else
C       Restore the extrapolated F in IRwF1 and extract the extrapolated
C       dE/dQ into 'EdQbar'.
        LenF = ITqry(IRwF1)
        Call CpyFi2(IRwF1U,0,IRwF1,0,LenF,V,MDV)
        If(IRwF2.ne.0) Call CpyFi2(IRwF2U,0,IRwF2,0,LenF,V,MDV)
        Call TstCor(NVQ,MDV,'PCMDIS-1')
        Call FileIO(2,-IRwF1U,NVQ,V,LenF)
        Call PCMFIO('write','EdQbar',V,NVQ,0,0)
        Call FileIO(2,-IRwP1U,NVQ,V,LenF)
        If(SPCM) then
          Call PCMFIO('write','Qele',V,NVQ,0,0)
        else
          Call PCMFIO('write','Qtele',V,NVQ,0,0)
          endIf
        Call FreeRW(IRwF1U)
        Call FreeRW(IRwF2U)
        Call FreeRW(IRwP1U)
        Call FreeRW(IRwP2U)
        Call FreeRW(IRwErQ)
        endIf
      Return
      End
*Deck PCOld
      Subroutine PCOld(IOp,LenP,NDimBl,NSpBlk,NBNO,IRwPA,IRwPB,IRwCA,
     $  IRwCB,IOCMOS,IRwOld,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Save or restore densities and MOs for internal restart of the SCF.
C     IOp = -1 ... Set IRwOld and save.
C            1 ... Save.
C            2 ... Restore.
C
      Dimension V(*)
C
      If(IOCMOS.eq.0) then
        LenCS = 0
      else
        LenCS = NSpBlk*NBNO
        endIf
      LenOld = 2*LenP + (3-NDimBl)*NBNO + LenCS
      If(IOp.eq.-1) then
        IRwOld = IGetRW(LenOld,0)
      else if(LenOld.ne.ITqry(IRwOld)) then
        Call GauErr('Inconsistent length for IRwOld in PCOld.')
        endIf
      If(IAbs(IOp).eq.1) then
        Call FCopyP(IRwPA,0,IRwOld,0,1,LenP,V,MDV)
        If(NDimBl.eq.1.and.ITqry(IRwPB).gt.0)
     $    Call FCopyP(IRwPB,0,IRwOld,-1,1,LenP,V,MDV)
        Call FCopyP(IRwCA,0,IRwOld,-1,1,NBNO,V,MDV)
        If(NDimBl.eq.1.and.ITqry(IRwCB).gt.0)
     $    Call FCopyP(IRwCB,0,IRwOld,-1,1,NBNO,V,MDV)
        If(IOCMOS.ne.0) Call FCopyP(IOCMOS,0,IRwOld,-1,1,LenCS,V,MDV)
      else if(IOp.eq.2) then
        Call FCopyP(IRwOld,0,IRwPA,0,1,LenP,V,MDV)
        If(NDimBl.eq.1.and.ITqry(IRwPB).gt.0)
     $    Call FCopyP(IRwOld,-1,IRwPB,0,1,LenP,V,MDV)
        Call FCopyP(IRwOld,-1,IRwCA,0,1,NBNO,V,MDV)
        If(NDimBl.eq.1.and.ITqry(IRwCB).gt.0)
     $    Call FCopyP(IRwOld,-1,IRwCB,0,1,NBNO,V,MDV)
        If(IOCMOS.ne.0) Call FCopyP(IRwOld,-1,IOCMOS,0,1,LenCS,V,MDV)
      else
        Call GauErr('Illegal IOp in PCOld.')
        endIf
      Return
      End
*Deck PhAPrd
      Function PhAPrd(JSt,JEnd,I,NTot,XYZ,GNorm,DQ0)
      Implicit Real*8(A-H,O-Z)
C
C     Product of a vector with A for PhInvI.
C
      Dimension XYZ(NTot,3), GNorm(NTot,3), DQ0(NTot)
      Save Zero
      Data Zero/0.0d0/
C
      Sum = Zero
      Do 10 J = JSt, JEnd
        dx = xyz(j,1) - xyz(i,1)
        dy = xyz(j,2) - xyz(i,2)
        dz = xyz(j,3) - xyz(i,3)
        d3 = dq0(j)/(Sqrt((dx**2+dy**2+dz**2)**3))
   10   Sum = Sum + d3*(dx*gnorm(J,1)+dy*gnorm(J,2)+
     $    dz*gnorm(J,3))
      PhAPrd = Sum
      Return
      End
*Deck PHess
      Subroutine PHess(IOut,NBasis,P,F,H,Scr)
      Implicit Real*8 (A-H,O-Z)
      Dimension P(*),F(*),H(NBasis,NBasis),Scr(NBasis,NBasis)
      Save One,Two,Three,Four,ThHess,DifMax
      Data One/1.D0/,Two/2.D0/,Three/3.D0/,Four/4.D0/,ThHess/1.D-1/,
     $  DifMax/1.D1/
 1000 Format(' The smallest Hessian is ',1D12.5,' Adjusted to ',1D12.5)
C
C     This subroutine forms an approximate density Hessian using only the
C     diagonal elements of P and F.
C
      Call Square(P,P,NBasis,NBasis,0)
      Call Square(F,F,NBasis,NBasis,0)
      Call MatMpy(NBasis,NBasis,NBasis,P,F,Scr)
      Call Linear(P,P,NBasis,NBasis)
      Call Linear(F,F,NBasis,NBasis)
      K = 0
      Do 10 I = 1, NBasis
        K = K + I
        M = 0
        Do 10 J = 1, I
          M = M + J
          H(I,J) = (Three - Two*P(K))*F(M) + (Three - Two*P(M))*F(K)
     $      - Four*Scr(I,I) - Four*Scr(J,J)
   10     H(J,I) = H(I,J)
C
C     Take care of negative or small eigenvalues
C
      HMin = ArrMin(H,NBasis*NBasis,.False.)
      If(HMin.lt.ThHess) then
        Write(IOut,1000) HMin,ThHess
        Diff = ThHess - HMin
        If(Diff.ge.DifMax) then
          Do 20 I = 1, NBasis
            Do 20 J = 1, NBasis
   20         H(I,J) = One
        else
          Do 30 I = 1, NBasis
            Do 30 J = 1, NBasis
   30         H(I,J) = H(I,J) + Diff
          endIf
        endIf
      Call Linear(H,H,NBasis,NBasis)
      Return
      End
*Deck POrtho
      Subroutine POrtho(IOut,IPrint,X,LenX,ibegn1,jip,ip,IRwPO,IRwPA,
     $  IRWS,IRWV,NBasis,NAtoms,ZThr,NOp1,NEqBas,IMemOp,NE,ZDO,LLim,
     $  IULim)
      Implicit real*8(a-h,o-z)
C
C     Retrieve and transform the Fock matrix
C
      Dimension X(LenX),NEqBas(nbasis,*),LLim(NAtoms),IULim(NAtoms)
      Logical ZDO
C
      Call TStamp(1,'Top of POrtho.')
      ix = ibegn1
      ibegin = ix + NBasis
      istart = ibegin
      If(ZDO) Then
C     Add in block diagonal elements
        jipa = jip
        ipa = ip
        Call DForm(X,LenX,X(jipa),X(ipa),LLim,IULim,X(ix),ibegin,
     $    istart,jip,ip,iend,NAtoms,NBasis,NZP,iprint)
        Call SpStor(IRWPO,X(jip),X(ip),NBasis,1,.True.,0)
        Call SpStor(IRwPA,X(jip),X(ip),NBasis,1,.True.,0)
      else
        NOp = NOp1
        If(NOp.le.1) NOp = 0
C           Transform Density back to AO basis
        Call SpStor(IRWPO,X(jip),X(ip),NBasis,1,.True.,0)
        Call BackT(X,LenX,NEqBas,NOp,NBasis,IBegin,IStart,jip,ip,
     $    iend,0,IRwPA,IOut,ZThr,NE,IRwS,IRwV,IPrint,IMemOp)
        EndIf
      If(IPrint.ge.3)
     $  Call SpDisp(IOut,X(jip),X(ip),NBasis,'Alpha Dens.',2)
      Call TStamp(1,'After U(-1)PU(-T)')
      Return
      End
*Deck ProgCv
      Subroutine ProgCv(IOut,IPrint,IVarAc,Acurcy,VarAcc,ZThr,AccDes,
     $  ZThrP,AccDsC,AccDsD,TKos3,TKos4,pfthr,Init,JCycle,AFinal,SpDiag,
     $  DPMax,Tight)
      implicit Real*8 (A-H,O-Z)
C
C     Determine some progressive convergence parameters
C
      Dimension IVarAc(*)
      Logical VarAcc, Init, SpDiag, DoPrnt, Tight
      Save ZThrPO,AccDsO,Two,Zero,One,Ten,TenM03,TenM02,FvTM04,TenM07
      Data ZThrPO/1.d0/,AccDsO/1.d0/,Zero/0.d0/,One/1.d0/,Ten/10.d0/,
     $  TenM03/1.0d-3/,TenM02/1.d-2/,FvTM04/5.0d-4/,TenM07/1.0d-7/,
     $  Two/2.d0/
 1100 Format(' Cur acc:       ',1PD12.6,' Min acc:       ',1PD12.6)
 1200 Format(' ZThr:          ',1PD12.6,' Final acc:     ',1PD12.6)
 1300 Format(' Coul acc:      ',1PD12.6,' Min acc:       ',1PD12.6)
 1400 Format(' DFT acc:       ',1PD12.6,' Min acc:       ',1PD12.6)
 1500 Format(' Switch for 35: ',1PD12.6,' Switch for 50: ',1PD12.6)
C
      DoPrnt = IPrint.ge.1.or.IPSave(0).eq.0
      AccDsC = AccDes
      AccDsD = AccDes
      ZThrP = Zero
      If(.not.SpDiag) Return
      ZThrP = ZThr
      PfScal = TenM02
      If(IVarAc(5).ne.0) PfScal = Ten**(gfloat(-IVarAc(5))/Ten)
      PfThr = ZThrP*PfScal
      If(.not.VarAcc) Return
      If(Init) then
C       Initialize zero threshold
        ZThrPO = One
        AccDsO = One
        IKos1 = Mod(IVarAc(2),100)
        IKos2 = Mod(IVarAc(2)/100,100)
        TKos1 = Ten**(gfloat(-40)/Ten)
        TKos2 = Ten**(gfloat(-25)/Ten)
        If(IKos1.ne.0) TKos1 = Ten**(gfloat(-IKos1)/Ten)
        If(IKos2.ne.0) TKos2 = Ten**(gfloat(-IKos2)/Ten)
        ZThrP = TKos1
        If(ZThrP.lt.ZThr.or.Tight) ZThrP = ZThr
        ZThrPO = ZThrP
        If(DoPrnt) Write(IOut,1200) ZThrP, ZThr
      else
C       Coulomb accuracy
        TKos1 = TenM03*Two
        TKos2 = TenM03
        IKos1 = mod(IVarAc(1),100)
        IKos2 = mod(IVarAc(1)/100,100)
        If(IKos1.gt.0) TKos1 = Ten**(gfloat(-IKos1)/Ten)
        If(IKos2.gt.0) TKos2 = Ten**(gfloat(-IKos2)/Ten)
        If(.not.Tight) then
          If(JCycle.lt.2) then
            AccDes = TKos2
          else
            AccDes = Min(DPMax*TKos1,TKos2)
            endIf
          If(AccDes.gt.AccDsO) AccDes = AccDsO
          endIf
        If(IPrint.gt.0) Write(IOut,1100) AccDes, TKos2
        AccDsO = AccDes
C       zero threshold
        IKos1 = Mod(IVarAc(2),100)
        IKos2 = Mod(IVarAc(2)/100,100)
        TKos1 = Ten**(gfloat(-40)/Ten)
        TKos2 = Ten**(gfloat(-25)/Ten)
        If(IKos1.ne.0) TKos1 = Ten**(gfloat(-IKos1)/Ten)
        If(IKos2.ne.0) TKos2 = Ten**(gfloat(-IKos2)/Ten)
        If(JCycle.lt.2) then
          ZThrP = TKos1
        else
          ZThrP = Min(DPMax*TKos2,TKos1)
          Endif
        If(ZThrP.gt.ZThrPO) ZThrP = ZThrPO
        If(JCycle.gt.10) ZThrP = Min(ZThrP,Ten*Ten*ZThr)
        If(JCycle.gt.20.or.Tight) ZThrP = ZThr
        ZThrP = Max(ZThrP,ZThr)
        ZThrPO = ZThrP
C       purification convergenc
        PfThr = ZThrP*PfScal
        If(DoPrnt) Write(IOut,1200) ZThrP, ZThr
C       coulomb accuracy
        TKos1 = One
        TKos2 = One
        IKos1 = mod(IVarAc(3),100)
        IKos2 = mod(IVarAc(3)/100,100)
        If(IKos1.gt.0) TKos1 = Ten**(gfloat(-IKos1)/Ten)
        If(IKos2.gt.0) TKos2 = Ten**(gfloat(-IKos2)/Ten)
        If(.not.Tight) AccDsC = Min(AccDes*TKos1,TKos2)
        If(IPrint.gt.0) Write(IOut,1300) AccDsC,TKos2
C       dft accuracy
        AccDsD = AccDes
        TKos1 = One
        TKos2 = One
        IKos1 = Mod(IVarAc(4),100)
        IKos2 = Mod(IVarAc(4)/100,100)
        If(IKos1.gt.0) TKos1 = Ten**(gfloat(-IKos1)/Ten)
        If(IKos2.gt.0) TKos2 = Ten**(gfloat(-IKos2)/Ten)
        If(.not.Tight) AccDsD = Min(AccDes*TKos1,TKos2)
        If(IPrint.gt.0) write(IOut,1400) AccDsD, TKos2
C       grid accuracy
        TKos3 = FvTM04
        TKos4 = TenM07
        IKos3 = mod(IVarAc(4)/10000,100)
        IKos4 = mod(IVarAc(4)/1000000,100)
        If(IKos3.gt.0) TKos3 = Ten**(gfloat(-IKos3)/Ten)
        If(IKos4.gt.0) TKos4 = Ten**(gfloat(-IKos4)/Ten)
        If(IPrint.gt.0) Write(IOut,1500) TKos3, TKos4
        endIf
      Return
      End
*Deck QNDIni
      Subroutine QNDIni(IOut,DoPrnt,DoDMS,DoSim,DoSimP,UseGEW,RHF,ZDO,
     $  IBas,NBasis,NBsUse,IRwS,IRwRAO,MaxCmp,IRwV,IRwSim,IRwPA,IRwPB,
     $  IRwFRA,IRwFRB,NAE,NBE,Acurcy,D,DD,VX,VV,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Initialization of the orthonormal transformation for FoSimult
C
      Logical DoPrnt, DoSim, DoSimP, UseGEW, RHF, DoDMS, ZDO
      Dimension D(NBasis,NBasis),DD(NBasis),VX(NBasis,NBasis),VV(*),
     $  V(*),XX(1)
      Save Zero, One, Two, XX
      Data Zero/0.0d0/, One/1.0d0/, Two/2.0d0/, XX/0.0d0/
 1000 Format(' Using Lowdin orthonormal basis with evecs from IRwRAO.')
 1010 Format(' Using Cholesky orthonormal basis.')
 1020 Format(' Using Lowdin orthonormal basis.')
 1030 Format(' ZDO:  unit orthonormal basis.')
C
      NTT = (NBasis*(NBasis+1))/2
      NBSq = NBasis**2
C
C     Form the transformation matrices to/from orthonormal
C     (IRwV contains U**-1 or S**(-1/2) first, followed by U or S**(1/2))
C
      Call ConDDF(IRwV,2*NBSq)
      If(ZDO) then
        If(DoPrnt) Write(IOut,1030)
        Call AUnitM(.False.,1,NBasis,NBasis,VX)
        Call FileIO(1,-IRwV,NBSq,VX,0)
        Call FileIO(1,IRwV,NBSq,VX,0)
      else if(IBas.eq.1) then
        Call ConDDF(IRwV,3*NBSq+NBasis)
        If(MaxCmp.eq.1) then
          If(DoPrnt) Write(IOut,1000)
          Call FileIO(2,-IRwRAO,NBasis*NBsUse,D,0)
          Do 10 I = 1, NBsUse
   10       DD(I) = One / SProd(NBasis,D(1,I),D(1,I))
          Call ScalD1(NBasis,NBasis,NBsUse,D,DD,.False.,2,IErr)
          Call GenIn1(.False.,.False.,.True.,NBasis,NBsUse,0,NRank,Zero,
     $      DD,D,XX,XX,VX,V,MDV)
          Call TstCor(NBSq,MDV,'QNDIni-Sq')
          Call Square(VX,V,NBasis,NBasis,3)
          Call FileIO(1,-IRwV,NBSq,V,0)
          Call GenIn1(.False.,.True.,.False.,NBasis,NBsUse,0,NRank,Zero,
     $      DD,D,XX,VX,XX,V,MDV)
          Call Square(VX,V,NBasis,NBasis,3)
          Call FileIO(1,IRwV,NBSq,V,0)
          Call FileIO(1,IRwV,NBasis,DD,0)
          Call FileIO(1,IRwV,NBSq,D,0)
        else
          If(DoPrnt) Write(IOut,1020)
          Call FormV(3,NBasis,NBasis,IRwS,IRwV,D,VX,DD,VV)
          endIf
      else if(IBas.eq.2) then
        If(DoPrnt) Write(IOut,1010)
        Call FileIO(2,-IRwS,NTT,D,0)
        Call CholIn(NBasis,D,VX)
        Call FileIO(1,-IRwV,NBSq,VX,0)
        Call FileIO(1,IRwV,NBSq,D,0)
      else
        Call GauErr('Illegal IBas in QNDIni.')
        endIf
      If(DoSim) then
        UseGEW = .True.
        Call FileIO(2,-IRwSim,1,Temp,0)
        Call BNUPak(Temp,NMax,MaxGDI)
        Call FileIO(2,-IRwSim,1,Temp,3)
        Call BNUPak(Temp,NPDIIS,NGDIIS)
        If(NPDIIS.lt.NGDIIS) Call GauErr('Consistency error in QNDINI.')
        Call FileIO(2,-IRwSim,1,SimThr,1)
        Acurcy = Max(SimThr,Acurcy)
        DoSimP = DoSimP.and.(NGDIIS.ge.2)
        If(DoSimP) then
          Call SimFP(IOut,DoPrnt,NAE,NBasis,IRwV,IRwSim,IRwPA,IRwFRA,
     $      NGDIIS,MaxGDI,NMax,DoDMS,IBas,V,MDV)
          If(.not.RHF) Call SimFP(IOut,DoPrnt,NBE,NBasis,IRwV,IRwSim,
     $      IRwPB,IRwFRB,NGDIIS,MaxGDI,NMax,DoDMS,IBas,V,MDV)
          If(RHF) then
            Call FileIO(2,-IRwPA,NTT,D,0)
            Call AScale(NTT,Two,D,D)
            Call FileIO(1,-IRwPA,NTT,D,0)
            endIf
          endIf
        endIf
      Return
      End
*Deck QNDMS
      Subroutine QNDMS(IOut,IPrint,NBasis,IRwP,IRwF,IRwPO,IRwFO,IRwV,Ne,
     $  A,NAMax,DMSErr,IBas)
      Implicit Real*8(A-H,O-Z)
C
C     Local density matrix search Algorithm using QN-DMS
C     with full matrices
C
C     IOut   ... Output unit.
C     IPrint ... Print flag:
C                0 ... print nothing.
C                1 ... print a few messages.
C                2 ... print B-matrix.
C                3 ... print everything
C     IRwPO  ... RWF on which current alpha density matrix in orthonormal
C                basis is stored (eq.0 requests transformation from the
C                AO basis)
C     IRwFO  ... RWF on which current alpha Fock matrix in orthonormal
C                basis is stored (eq.0 requests transformation from the
C                AO basis)
C     IRwP   ... RWF on which current alpha density matrix in AO basis
C                is stored (new density is returned here)
C     IRwF   ... RWF on which current alpha Fock matrix in AO basis
C                is stored
C     A      ... NAMax scratch matrix
C
      Logical Inv2, OK
      Dimension A(NAMax), VDummy(1)
      Save One, QNConv
      Data One/1.0d0/, QNConv/1.D-10/
 1000 Format(' Local density matrix search using QN-DMS')
 1010 Format(' Coeff:     ',6D12.5)
 1020 Format(' Delta-E ',1PD20.6,1PD12.3)
C
C     Initilization
C
      Write(IOut,1000)
      MaxItr = 4
      NTT    = (NBasis*(NBasis+1))/2
      NBSq   = NBasis**2
      MaxIt1 = MaxItr + 1
      IBMat  = 1
      IScr   = IBMat + MaxIt1**2
      IEnd   = IScr  + MaxIt1*5
      Do 100 Iter = 1, MaxItr
        If(Iter.eq.1) then
C
C         At the first iteration, take a Newton-Raphson step
C
          Call QNErr(IOut,IPrint,NBasis,Iter,MaxItr,IRwP,IRwF,IRwV,
     $      IRwFO,IRwPO,IRwFT,IRwPT,IOP,IOG,IOErr,A,NAMax,IBas,DMSErr)
          If(IPrint.ge.1) Write(IOut,1020) DMSErr,QNConv
          If(DMSErr.lt.QNConv) Goto 300
          A(IScr) = One
          Write(IOut,1010) A(IScr)
          Call NewP(IOut,IPrint,Iter,NBasis,Ne,A(IScr),IOP,IOErr,IRwPT,
     $      A(IEnd),NAMax-IEnd)
        else
C
C         Starting from the second loop, do a DIIS extrapolation using
C         the quasi-Newton-Raphson step -dEdP/H as the error vector
C
          Call QNErr(IOut,IPrint,NBasis,Iter,MaxItr,IRwP,IRwF,IRwV,
     $      IRwFO,IRwPO,IRwFT,IRwPT,IOP,IOG,IOErr,A,NAMax,IBas,ErrTmp)
          If(IPrint.ge.1) Write(IOut,1020) ErrTmp,QNConv
          If(ErrTmp.lt.QNConv) Goto 200
          Call FormBM(Iter,NBasis,A(IBMat),MaxIt1,IOErr,A(IEnd),
     $      NAMax-IEnd)
C
C         Invert Pulay's B matrix, using a routine which ensures numerical
C         stability even if the matrix is quite ill-conditioned.
C
          NDim = Iter + 1
          Call Pcck(0,KK,A(IBMat),MaxIt1,MaxIt1,NDim,NDim)
          OK = Inv2(IOut,0,A(IBMat),NDim,A(IScr),A(IScr+MaxIt1*2),
     $         A(IScr+MaxIt1*3),A(IScr+MaxIt1*4),Det,VDummy,VDummy,
     $         VDummy,0)
          Call AClear(NDim,A(IScr))
          A(IScr) = -One
          If(OK) then
            Call MatMpy(1,NDim,NDim,A(IScr),A(IBMat),A(IScr+MaxIt1))
          else
            Call AClear(NDim,A(IScr+MaxIt1))
            A(IScr+MaxIt1+Iter) = One
            endIf
C
C         Extrapolated / Interpolate the density matrices.
C
          Write(IOut,1010) (A(IScr+MaxIt1+I),I=1,Iter)
          Call NewP(IOut,IPrint,Iter,NBasis,Ne,A(IScr+MaxIt1+1),IOP,
     $      IOErr,IRwPT,A(IEnd),NAMax-IEnd)
          endIf
  100   Continue
  200 Call FileIO(2,-IRwV,NBSq,A,0)
      Call FileIO(2,-IRwPT,NTT,A(1+NBSq),0)
      Call OTran(A(1+NBSq),A,A(1+2*NBSq),NBasis,NBasis,1,IBas,
     $  A(1+3*NBSq),NAMax-1-3*NBSq)
      Call FileIO(1,-IRwP,NTT,A(1+NBSq),0)
  300 Call FreeRW(IOErr)
      Call FreeRW(IOP)
      Call FreeRW(IOG)
      Call FreeRW(IRwPT)
      Call FreeRW(IRwFT)
      Return
      End
*Deck QNErr
      Subroutine QNErr(IOut,IPrint,NBasis,NSaved,MaxItr,IRwP,IRwF,IRwV,
     $  IRwFO,IRwPO,IRwFT,IRwPT,IOP,IOG,IOErr,A,NAMax,IBas,DMSErr)
      Implicit Real*8(A-H,O-Z)
C
C     This subroutine forms the error vectors for QN-DIIS
C     e = -g/H in full matrix form
C
      Dimension A(*)
      Save One
      Data One/1.D0/
C
C     Do some memory allocation
C
      NTT   = (NBasis*(NBasis+1))/2
      NBSq  = NBasis**2
      IP    = 1
      IFk   = IP  + NBSq
      IV    = IFk + NBSq
      IH    = IV  + NBSq
      IScr1 = IH  + NBSq
      IEnd  = IScr1 + NBSq
      MDV   = NAMax - IEnd
      If(NSaved.eq.1) then
        Leng  = MaxItr*NTT
        IOP   = IGetRW(Leng,0)
        IOG   = IGetRW(Leng,0)
        IOErr = IGetRW(Leng,0)
        IRwPT = IGetRW(NTT,0)
        IRwFT = IGetRW(NTT,0)
        endIf
C
C     Transform the density matrix and Fock matrix into orthonormal basis
C
      If(NSaved.eq.1) then
C
C       Transform Fock and dnesity matrices into the orthonormal basis
C
        If(IRwFO.ne.0) then
          Call FileIO(2,-IRwFO,NTT,A(IFk),0)
          Call FileIO(1,-IRwFT,NTT,A(IFk),0)
        else
          Call FileIO(2,-IRwV,NBSq,A(IV),0)
          Call FileIO(2,-IRwF,NTT,A(IFk),0)
          Call OTran(A(IFk),A(IV),A(IP),NBasis,NBasis,1,IBas,A(IEnd),
     $      MDV)
          Call FileIO(1,-IRwFT,NTT,A(IFk),0)
          endIf
        If(IRwPO.ne.0) then
          Call FileIO(2,-IRwPO,NTT,A(IP),0)
          Call FileIO(1,-IOP,NTT,A(IP),0)
        else
          Call FileIO(2,-IRwV,NBSq,A(IV),NBSq)
          Call FileIO(2,-IRwP,NTT,A(IP),0)
          Call OTran(A(IP),A(IV),A(IH),NBasis,NBasis,0,IBas,A(IEnd),
     $      MDV)
          Call FileIO(1,-IOP,NTT,A(IP),0)
          endIf
      else
        Call FileIO(2,-IRwFT,NTT,A(IFk),0)
        Call FileIO(2,-IRwPT,NTT,A(IP),0)
        Call FileIO(1,-IOP,NTT,A(IP),(NSaved-1)*NTT)
        endIf
C
C     Calculate the energy gradient with respect to the density
C
      Call dEdP(A(IEnd),MDV,Junk,IDummy,A(IP),A(IFk),A(IV),A(IH),
     $  A(IScr1),NBasis,NTT,.False.,.True.)
      Call AScale(NTT,-One,A(IV),A(IV))
      Call FileIO(1,-IOG,NTT,A(IV),(NSaved-1)*NTT)
C
C     Calculate the Newton-Raphson step
C
      Call PHess(IOut,NBasis,A(IP),A(IFk),A(IH),A(IScr1))
      Do 10 I = 1, NSaved
        Call FileIO(2,-IOG,NTT,A(IP),(I-1)*NTT)
        Call ADiv(NTT,A(IP),A(IH),A(IP))
        Call FileIO(1,-IOErr,NTT,A(IP),(I-1)*NTT)
   10   Continue
      Call FileIO(2,-IOG,NTT,A(IV),(NSaved-1)*NTT)
      Call FileIO(2,-IOErr,NTT,A(IP),(NSaved-1)*NTT)
      DMSErr = SCFTrc(A(IP),A(IV),NBasis,1)
      Return
      End
*Deck RenZ4
      Subroutine RenZ4(CC,Ncore,Norbs,Iout,IEnd,Thresh,IStep,X,jiz,iz,
     $  NZz,jiz2,iz2,NZz2,NElect,z2u,z2l)
      Implicit Real*8(A-H,O-Z)
C
C     Renormalize z and adjust its Trace to the correct number of electrons
C
      Dimension CC(Ncore),X(Norbs)
      Save Zero,Three
      Data Zero/0.d0/,Three/3.d0/
c
      ZU = Sqrt(z2u)
      If(z2l.lt.Zero) z2l = Zero
      ZL = Sqrt(z2l)
C
C     Check spectral width
C
      f2 = Sqrt(Three*(ZU - ZL)/(ZU*ZU*ZU - ZL*ZL*ZL))
C
C     Normalize z
C
      Call AScale(Norbs+NZz,f2,CC(iz),CC(iz))
C
C     Normalize z2
      f3 = f2*f2
      Call AScale(Norbs+NZz2,f3,CC(iz2),CC(iz2))
      Return
      End
*Deck RenZ
      Subroutine RenZ(CC,Ncore,Norbs,Iout,IEnd,Thresh,X,jiz,iz,NZz,
     $  NElect)
      Implicit Real*8(A-H,O-Z)
C
C     Renormalize z and adjust their Trace to the correct number of electrons
C
      Dimension CC(Ncore),X(Norbs)
      Save Zero,XInfty
      Data Zero/0.d0/,Xinfty/1.d+99/
c
      fnne = gfloat(Norbs-2*NElect)
      fn = gfloat(Norbs)
      trz = ArrSum(Norbs,CC(iz))/fn
      do 10 i = 0, (Norbs-1)
   10   CC(iz+i) = CC(iz+i) - trz
C
C     Bounds of z-trz
C
      Call GBnds(CC(jiz),CC(iz),Norbs,NZz,zabs,ZU,ZL,Iout,X)
      If(Abs(ZU).gt.Thresh) then
        xn1 = (fn-fnne)/ZU
      else
        xn1 = Xinfty
        endIf
      If(Abs(ZL).gt.Thresh) then
        xn2 = (-fn-fnne)/ZL
      else
        xn2 = Xinfty
        endIf
      xn = xn1
      If(xn2.gt.Zero.and.xn2.lt.xn) xn = xn2
C
C     Normalize z
C
      f1 = xn/fn
      Call AScale(Norbs+NZz,f1,CC(iz),CC(iz))
      f2 = fnne/fn
      do 20 i = 0, (Norbs-1)
   20   CC(iz+i) = CC(iz+i) + f2
      Return
      End
*Deck RHat
      Subroutine RHat(NTot,NOrg,XYZ,COrg,ICen,Rhatx,Rhaty,Rhatz,DMag)
      Implicit Real*8(A-H,O-Z)
      Dimension XYZ(NTot,3), COrg(3,NOrg), ICen(NTot), Rhatx(NTot),
     $  RHaty(NTot), Rhatz(NTot), DMag(NTot)
      Save One
      Data One/1.0d0/
C
      Do 10 I = 1, NTot
        Xc = XYZ(I,1) - COrg(1,ICen(I))
        Yc = XYZ(I,2) - COrg(2,ICen(I))
        Zc = XYZ(I,3) - COrg(3,ICen(I))
        Dmag(i) = One/(sqrt(Xc*Xc+Yc*Yc+Zc*Zc))
        rhatx(i) = Xc*dmag(i)
        rhaty(i) = Yc*dmag(i)
        rhatz(i) = Zc*dmag(i)
   10   Continue
      Return
      End
*Deck ROFock
      Subroutine ROFock(IOut,IPrint,NBasis,NOrbs,NAE,NBE,IOFRO,IOC,IOS,
     $  FA,FB,D,V,Scr)
      Implicit Real*8(A-H,O-Z)
C
C     Take the AO alpha and beta Fock matrices in FA and FB, and
C     store the AO combined ROHF Fock matrix in IOFRO.  The combined
C     matrix is (McWheeny and Diercksen, JCP 49(11) 4852-4856 (1968)):
C
C     F' = a*(Pc+Pv)Fc(Pc+Pv) + b*(Po+Pv)Fo*(Po+Pv) +
C          c*(Pc+Po)(2Fc-Fo)(Pc+Po)
C     where Fc = (Falpha+Fbeta)/2, Fo = Falpha.
C     a, b, and c are arbitrary non-zero constants.  The choice
C     a=2/3, b=1/3, c=1/3 makes F'=Fc for closed-shell wavefunctions.
C
      Logical JC, JO, JV, IC, IO, IV
      Dimension FA(*), FB(*), D(*), V(NBasis,NBasis), Scr(*)
      Save Zero, One, Two, Three, Pt5
      Data Zero/0.0d0/, One/1.0d0/, Two/2.0d0/, Three/3.0d0/, Pt5/0.5d0/
C
C     Transform to the current MO basis.
C
      A = Two / Three
      B = One / Three
      C = One / Three
      NTT = (NBasis*(NBasis+1))/2
      NBSq = NBasis*NBasis
      Call FileIO(2,-IOC,NBasis*NOrbs,V,0)
      Call VFVMul(.False.,.False.,.False.,0,NBasis,NOrbs,NBasis,FA,V,D,
     $  FA,Scr,NBasis)
      Call VFVMul(.False.,.False.,.False.,0,NBasis,NOrbs,NBasis,FB,V,D,
     $  FB,Scr,NBasis)
C
C     Now form the merged Fock matrix in D in the MO basis.
C
      IJ = 0
      A2 = A / Two
      Do 100 J = 1, NOrbs
        JC = J.le.NBE
        JV = J.gt.NAE
        JO = .not.JC.and..not.JV
        Do 100 I = 1, J
          IJ = IJ + 1
          IC = I.le.NBE
          IV = I.gt.NAE
          IO = .not.IC.and..not.IV
          X = Zero
          If(IC.and.JC) then
            X = A2*FA(IJ) + (A2+C)*FB(IJ)
          else if(IC.and.JO) then
            X = FB(IJ)
          else if(IO.and.JO) then
            X = B*FA(IJ) + C*FB(IJ)
          else if(IC.and.JV) then
            X = Pt5*(FA(IJ)+FB(IJ))
          else if(IO.and.JV) then
            X = FA(IJ)
          else if(IV.and.JV) then
            X = (A2+B)*FA(IJ) + A2*FB(IJ)
            endIf
  100     D(IJ) = X
C
C     Finally, transform back to the AO basis:
C     F(AO) = S*C*F(MO)*C(Dagger)*S.
C
      If(IPrint.ge.4) Call LTOutS(IOut,
     $  'Combined Fock matrix (MO basis):',0,NOrbs,D,1)
      NOTT = (NOrbs*(NOrbs+1))/2
      Call AMove(NOTT,D,FA)
      Call FileIO(2,-IOC,NBasis*NOrbs,D,0)
      If(IOS.eq.0) then
        Call AMove(NBSq,D,V)
      else
        Call FileIO(2,-IOS,NTT,FB,0)
        Call MultLV(.True.,NOrbs,0,NBasis,NBasis,FB,D,V)
        endIf
      Call AClear(NBasis*(NBasis-NOrbs),V(1,NOrbs+1))
      Call Trspn2(NBasis,V)
      Call Pcck(0,Junk,V,NBasis,NBasis,NOrbs,NBasis)
      Call VFVMul(.False.,.False.,.False.,0,NOrbs,NBasis,NBasis,FA,V,D,
     $  FA,Scr,NBasis)
      Call FileIO(1,-IOFRO,NTT,FA,0)
      If(IPrint.ge.3) Call LTOutS(IOut,
     $  'Combined Fock matrix (AO basis):',0,NBasis,FA,1)
      Return
      End
*Deck SCFDIX
      Subroutine SCFDIX(IOut,IPrint,IDITyp,IFSelI,RComp,SpDiag,JCycle,
     $  LvFirst,Rises,FixRis,OldRis,IStab,MxSCFC,NIter,IOpCl,NBasis,
     $  NBsUse,ISCRF,IRwP1,IRwP1X,IRwF1,IRwF1X,IRwP2,IRwP2X,IRwF2,
     $  IRwF2X,IRwS1,IRwTM,IRWFOA,IRwFOB,IRwPOA,IRwPOB,IRWFM1,IRWFM2,
     $  ZThrP,NSvDIS,IOBDIS,IOPADI,IOPBDI,IOFADI,IOFBDI,IOErDI,D,F,V,VV,
     $  BMat,Scr2,ErrMax,NErSav,IndFok,Energy,EnCyc,ErrCyc,IMemOp,NumFS,
     $  OldCof,DISPCM,IFrgTp,NBlks,NBF,NBFU,MapF2M,NOccSy,IOTMA,IOTMB,
     $  IndTMB,IOCAS,IOCBS,IndCBS,IOSBlk,IndSBB,X,LenX)
      Implicit Real*8(A-H,O-Z)
C
C     Call the appropriate routines for DIIS extrapolation.
C     DIIS extrapolation for SCF.  See
C     Pulay, J. Comp. Chem. Vol. 3 No. 4 556-560 (1982).
C
      Dimension D(*), F(*), V(NBasis,*), VV(NBasis,*), BMat(MxSCFC+1,*),
     $  Scr2(MxSCFC+1,5), IndFok(MxSCFC), EnCyc(MxSCFC), ErrCyc(MxSCFC),
     $  ErrMax(0:NErSav), X(*), OldCof(2), NBF(*), NBFU(*), MapF2M(*),
     $  NOccSy(*)
      Integer FixRis
      Logical RComp,Rises,OldRis,SpDiag,CBfn,Spinor,PCM,DISPCM,DbgFrg
      Save Zero, One
      Data Zero/0.0d0/, One/1.0d0/
 1000 Format(' Inconsistency:  DISPCM=',L1,' but NVQ=',I12,'.')
C
      Call DecoSC(IOpCl,NSpBlk,NRI,NDimBl,CBfn,NSpBlX,Spinor)
      NBD = NBasis*NDimBl
      NBUD = NBsUse*NDimBl
      NTTUR = NRI*(NBUD*(NBUD+1))/2
      If(SpDiag) then
        Call SpDIIS(IOut,IPrint,X,LenX,1,IRwFOA,IRwFOB,IRwPOA,IRwPOB,
     $    IRwFM1,IRwFM2,NBasis,JCycle,MxSCFC,ZThrP,IMemOp)
        NSvDIS = 1
        OldCof(1) = Zero
        OldCof(2) = One
      else
        Call AMove1(NErSav,0,1,ErrMax)
        If(IFSelI.ne.0) then
          IFSel = IFSelI
        else if(IStab.eq.2.and.JCycle.lt.LvFirst) then
          IFSel = 3
        else if(IStab.eq.2.and.JCycle.eq.LvFirst) then
          IFSel = 4
        else if(Rises.and.
     $    (FixRis.eq.4.or.(FixRis.ge.2.and..not.OldRis))) then
          IFSel = 5
        else
          IFSel = 0
          endIf
        If(RComp) then
          IndF = NBasis**2
          LenF = NBasis
        else
          IndF = 1
          LenF = NRI*(NBasis*NDimBl)**2
          endIf
        If(PCM(ISCRF)) then
          Call PCMDIS(IOut,IPrint,.True.,NumFS,NBasis,NRI,NDimBl,NVQ,
     $      IRwErQ,IRwF1,IRwF2,IRwP1X,IRwP2X,IRwF1U,IRwF2U,IRwP1U,
     $      IRwP2U,WFe,WQe,EWFe,EWQe,X,LenX)
        else
          NVQ = 0
          endIf
        If(NVQ.gt.0) then
          IRwF1Y = 0
          IRwF2Y = 0
        else
          IRwP1U = IRwP1X
          IRwF1U = IRwF1
          IRwF1Y = IRwF1X
          IRwP2U = IRwP2X
          IRwF2U = IRwF2
          IRwF2Y = IRwF2X
          IRwErQ = 0
          WFe = One
          WQe = Zero
          EWFe = One
          EWQe = Zero
          endIf
        If(DISPCM.neqv.NVQ.ne.0) then
          Write(IOut,1000) DISPCM, NVQ
          Call Lnk1E(0)
          endIf
        DbgFrg = .False.
        IRwErr = 0
        Call SCFDIS(IOut,IPrint,IDITyp,IFSel,MxSCFC,NIter,NumFS,IOpCl,
     $    NBasis,NBsUse,NVQ,IRwP1,IRwP1U,IRwF1U,IRwF1Y,IRwP2,IRwP2U,
     $    IRwF2U,IRwF2Y,IRwS1,IRwTM,IRwErr,IRwErQ,NSvDIS,IOBDIS,IOPADI,
     $    IOPBDI,IOFADI,IOFBDI,IOErDI,Zero,D,V,VV,BMat,Scr2,ErrMax(0),
     $    IndFok,Energy,EnCyc,ErrCyc,WFe,WQe,EWFe,EWQe,F(IndF),LenF)
        If(NVQ.gt.0) Call PCMDIS(IOut,IPrint,.False.,NumFS,NBasis,NRI,
     $    NDimBl,NVQ,IRwErQ,IRwF1,IRwF2,IRwP1X,IRwP2X,IRwF1U,IRwF2U,
     $    IRwP1U,IRwP2U,WFe,WQe,EWFe,EWQe,X,LenX)
        Do 200 I = Abs(NIter), NErSav
  200     ErrMax(I) = ErrMax(I-1)
        If(NSvDIS.ge.2) then
          OldCof(1) = Scr2(NSvDIS,1)
        else
          OldCof(1) = Zero
          endIf
        OldCof(2) = Scr2(NSvDIS+1,1)
        Call TStamp(1,'After DIIS')
        endIf
      Return
      End
*Deck SCGDMS
      Subroutine SCGDMS(IOut,IPrint,PrntDM,cc,maxcor,N,ZThr,Dens,NAE,
     $  DLaMax,IOpDMe,pfthr,JCycle,jip,ip,if,NZ,istrt1,iend,IMemO2)
      Implicit Real*8(A-H,O-Z)
      Dimension CC(maxcor),Z(1)
      Logical FastM, PrntDM
      Save dthrold,LX,LIX,pt1,Z,Zero,AcInit,DLamDF,DempMx,LBin
      Data dthrold/1.d-4/,LX/11/,LIX/0/,pt1/.1d0/,Z/0.0d0/,Zero/0.0d0/,
     $  AcInit/.05d0/,DLamDF/0.1d0/,LBin/0/,DempMx/5.d-3/
 1030 Format(' Zero threshold:                   ',1pd12.2)
 1060 Format(' C.G. Convergence criteria:        ',1pd12.2)
 1070 Format(' Maximum number of C.G. iterations:',i5)
 1075 Format(' Maximum purification cycles:      ',i5)
 1090 Format(' Maximum CGDIIS vectors:           ',i5)
 1100 Format(' Polak-Ribiere minimization')
 1110 Format(' Fletcher-Reeves minimization')
 1120 Format(' Using Jacobi preconditioning')
 1130 Format(' No preconditioning')
 1135 Format(' Maximum C.G. Step Size         ',f11.2)
 1200 Format(/' ********* CGDMS cycle ',i4/)
 1300 Format(' Using Original Sparse Matrix Multiplication')
 1400 Format(' Using Bin Sorting Sparse Matrix Multiplication')
 1500 Format(' Idempotency:',1pd12.2,' Criterion:',1pd12.2,/,
     $  ' Halting: Density matrix is no longer idempotent.'/,
     $  ' Please try running with tighter neglect thresholds',/,
     $  ' and smaller CG maximum step size.')
 1600 Format(' CGDMS It:',I3,' NCG=',i6,' Delta Func:',1PD12.2,' Grad:',
     $  1PD12.2,' Idemp:',1PD12.2)
C
C     Set convergence threshold and other useful stuff
C
      istart = istrt1
      IMemOp = IMemO2
      Method = Mod(IMemOp,10)
      iBnFac = istart
      iSizCV = istart
      If(NAE.eq.0) then
        Call AClear(NZ,CC(IP))
        Return
        endIf
      nbin = 0
      If(Method.ne.1) then
        Call FSIni1(cc,maxcor,istart,iBnFac,iSizCV,iend,cc(jip),cc(ip),
     $      cc(jip),cc(ip),N,NBin,IMemOp,method,ZThr,.true.,.true.)
        If(Method.eq.2) istart = iend
        endIf
      If(Method.ne.2) then
        If(PrntDM) Write(IOut,1300)
        nbin = 0
      else if(PrntDM) then
        Write(IOut,1400)
        endif
      FastM = (Method.eq.2)
      ItMax =  Mod(IOpDMe,100)
      If(ItMax.eq.0) ItMax = 4
      NPMax = Mod(IOpDMe,10000)/100
      If(NPMax.eq.0) NPMax = 3
      MaxCom = 0
      if(Mod(IOpDMe,100000)/10000.ne.0) MaxCom=ItMax
      IPreC = Mod(IOpDMe,1000000)/100000
      IPr = IOpDMe/1000000
      If(maxcom.le.2) maxcom = 0
      If(PrntDM) then
        Write(IOut,1030) ZThr
        Write(IOut,1060) dens
        Write(IOut,1070) itmax
        Write(IOut,1075) npmax
        If(maxcom.gt.0) Write(IOut,1090) maxcom
        If(ipr.eq.0) then
          Write(IOut,1100)
        else
          Write(IOut,1110)
        endIf
        If(iprec.eq.0) then
          Write(IOut,1120)
        else
          Write(IOut,1130)
        endIf
        If(DLAMax.gt.zero) Write(IOut,1135) DLAMax
        endIf
      If(IPrint.ge.1) write(IOut,1200) jcycle
      Accuc = AcInit
      dloose=pt1*accuc
      dtight=min(dloose,dthrold)
      dthr=max(dens,dtight)
c     dthr=dens
      dthrold=dthr
c
c  Allocate memory for all matrices used by CG-DMS
c
      LBin  = 0
      If(Method.eq.2) LBin = 2
      LIX   = 2
      If(Method.eq.2) LIX = 5
      LenQ = NZ+N
      if(FastM) LenQ = 2*NZ + N
      LSL = 0
      LISL = 0
      LISQ = 0
      if(FastM) LSL  = NZ + N
      If(FastM) LISL = InToWP(NZ+N+1)
      If(FastM) LISQ = InToWP(2*NZ+2*N+1)
      NP = NProc(0)
      ig    = istart
      ih    = ig   + NZ+N
      ixi   = ih   + NZ+N
      isr1  = ixi  + NZ+N
      isr2  = isr1  + LenQ
      jimap = isr2  + LenQ
      jmap  = jimap+ INTOWP(N+1+NZ)
      jiscr = jmap + INTOWP(NZ)
      jiscrq= jiscr + LISQ
      jisrq2= jiscrq + LISQ
      jibt  = jisrq2 + LISQ
      ibt    = jibt + LISL
      iRowNo = ibt  + LSL
      ix     = iRowNo + InToWP(N)
      jixx  = ix   + LX*N*NP
      jibin  = jixx + InToWP(LIX*N*NP)
      last  = jibin + InToWP(LBin*N*(NBin+1))
      if(FastM) then
        jiscr1 = last
        last   = jiscr1 + InToWP(N+NZ+1)
      else
        jiscr1 = jip
        endif
      Call TstCor(last,maxcor,'CGPS 3')
c
c     Fixed version - all matrix forms held constant
c     Make Fock and Density share common element
c     [if supp, include form of PF and FP in the common form]
c
      Call MkMap(cc(jip),cc(jimap),cc(jmap),N)
      Call CheckF(cc(jimap),cc(jmap),cc(ix),N,0,NZ,'Map',-6,IPrint)
      Call FxCgCy(IOut,IPrint,cc(jip),cc(ip),cc(if),cc(ig),cc(ih),
     $  cc(ixi),cc(isr1),cc(isr2),cc(jimap),cc(jmap),cc(ix),LX,ifixed,
     $  N,iexit,dt2,npfcy,NP,ItMax,MaxCom,IPr,Z(1),NAE,IPreC,DThr,
     $  DLaMax,DLamDF,NPMax,tcg2,pfthr,cc(jiscrq),LenQ,cc(jiscr),
     $  zthr,IMemOp,cc(jiscr1),cc(iRowNo),cc(jibt),cc(ibt),
     $  cc(iBnFac),cc(iSizCV),NBin,Method,cc(jixx),lix,cc(jibin),LBin,
     $  cc(jisrq2),demp,Ener,RMSGrd)
C
C     Clean up - store new density matrix
C
      If(IPrint.eq.0) Write(IOut,1600) jcycle,ifixed,Ener,RMSGrd,demp
      If(JCycle.gt.4.and.(Abs(Demp).gt.DempMx)) then
        Write(IOut,1500) Demp,DempMx
        Call GauErr('Stopping in SCGDMS.')
        endIf
      Call CheckF(cc(jip),cc(ip),cc(ix),N,0,NZ,'Conv. density',5,IPrint)
      Return
      End
*Deck SChbEM
      Subroutine SChbEM(IOut,IPrint,PrntCG,cc,ncore,X,istrts,jip,ip,if,
     $  iend,nz,N,niter,npfcy,NOcc,iB,IOpMe,xmax,xmin,demp,pfthr,IMemOp,
     $  zthr)
      Implicit Real*8(A-H,O-Z)
      Dimension cc(ncore),X(N)
      Logical StoreT, canpf, PrntCG
      Save DempMx
      Data DempMx/5.d-3/
 1000 Format(/' Chebyshev Expansion Method')
 1010 Format(' ChbDim Real*8 Arrays:  15*n + 6.5*nz ',/,
     $       ' Memory needed for CEM',i14,' Max memory ',i14)
 1100 Format(' Chebyshev Polynomial Order',i8)
 1150 Format(' Ficticious inverse temperature',f12.2)
 1200 Format(' Maximum Number of Purification Cycles',i8)
 1300 Format(' Storing TT on Disk')
 1400 Format(' Storing TT in Memory')
 1500 Format(' Canonical Purification')
 1600 Format(' Grand Canonical Purification')
 1800 Format(' Purification Threshold',1pd12.2)
 1700 Format(' Maximum Eigenvalue:',f9.4,'Minumum:',f9.4)
 2300 Format(' Using Original Sparse Matrix Multiplication')
 2400 Format(' Using Bin Sorting Sparse Matrix Multiplication')
 2500 Format(' Idempotency:',1pd12.2,' Criterion:',1pd12.2,/,
     $  ' Halting: Density matrix is no longer idempotent.'/,
     $  ' Please try running with tighter neglect thresholds,',/,
     $  ' A Larger CPO, and a Smaller Temperature')
 2600 Format(' Chebyshev Expansion Method Idempotency:',1PD12.2)
 2700 Format(' Using expansion of Fermi-Dirac distribution')
 2800 Format(' Using expansion of complementary error function')
C
C     Initialize for conjugate gradient
C
      NP = NProc(0)
      istart = istrts
      ItMax = Mod(IOpMe,1000)
      If(ItMax.eq.0) ItMax = 100
      NPMax = Mod(IOpMe,100000)/1000
      IFunc = Mod(IOpMe,1000000)/100000
      If(IFunc.ge.1) then
        B = GFloat(iB)/GFloat(10)
        If(iB.eq.0) B = 46.0d0
        Write(IOut,2700)
      else
        B = GFloat(iB)/GFloat(10)
        If(iB.eq.0) B = 16.0d0
        Write(IOut,2800)
        endIf
      StoreT = .True.
      CanPf = .True.
      If(IPrint.ge.1) Write(IOut,1000)
      If(PrntCG) then
        Write(IOut,1000)
        Write(IOut,1100) ItMax
        Write(IOut,1150) B
        If(StoreT) then
          Write(IOut,1300)
        else
          Write(IOut,1400)
          endIf
        If(CanPf.and.NPMax.gt.0) then
          Write(IOut,1500)
          Write(IOut,1200) NPMax
          Write(IOut,1800) pfthr
        else if(NPMax.gt.0) then
          Write(IOut,1600)
          Write(IOut,1200) NPMax
          Write(IOut,1800) pfthr
          endIf
        endIf
      numc = itmax
      numbig=5
      numsml=3
C
C     Find the lowest and highest eigenvalues of F
C
      jiX  = istart
      ipp  = jix + intowp(N)
      iend = ipp + N+NZ
      istart = iend
      call tstcor(iend,ncore,'SChbEM')
      call amove(N+NZ,cc(if),cc(ipp))
      Call DavGts(IOut,IPrint-4,cc,ncore,istart,jIP,ipp,N,xmax,xmin)
      If(IPrint.ge.2) Write(IOut,1700) xmax, xmin
C
C     Handle memory for CEM
C
      LenRT = NZ+N
      If(StoreT) LenRT = N
      LIX = 3
      LXX = 2
      LBin = 0
      NpNZ = 0
      Method = Mod(IMemOp,10)
      iend = istart
      NBin = 0
      iBnFac = istart
      iSizCV = istart
      If(Method.ne.1) Call FSIni1(cc,ncore,istart,iBnFac,iSizCV,iend,
     $  cc(jIP),cc(iP),cc(jIP),cc(iF),N,NBin,IMemOp,method,ZThr,.True.,
     $  .True.)
      jiscr = jip
      if(Method.eq.2) then
        If(IPrint.ge.1) Write(IOut,2400)
        LIX = 6
        LXX = 4
        LBin = 2
        NpNZ = N+NZ
        jiscr = iend
        iend  = jiscr + InToWP(N+1+NZ)
      else
        If(IPrint.ge.1) Write(IOut,2300)
        endIf
      irt    = iend
      ixx    = irt    + LenRT*Numc
      iend1  = ixx    + N*NP*LXX
      ic     = ixx    + N*NP
      icc    = ic     + Numc
      iend2  = icc    + Numc*4
      jIBin  = iend1
      jibt   = jIBin + InToWP(LBin*N*(NBin+1))
      ibt    = jibt + InToWP(NpNZ+1)
      itt    = ibt + NpNZ
      jix    = itt    + (nz+N)*3
      jiMap  = jix    + InToWP(N*NP*LIX)
      iMap   = jiMap  + InToWP(N+1+NZ)
      iend3  = iMap   + InToWP(NZ)
      iend = Max(iend3,iend2)
      call aclear((iend-irt+1),cc(irt))
      if(IPrint.ge.2) write(Iout,1010) iend, ncore
      call tstcor(iend,ncore,'ChbDim 2')
C     compress memory
      call aclear(N+NZ,cc(ip))
      call aclear(3*(N+NZ),cc(itt))
      Lfail = 0
      BB=B
      call CEMCyc(cc(jip),cc(if),cc(icc),cc(ic),cc(iTT),cc(ip),cc(irt),
     $  cc(iXX),cc(jiX),cc(jiMap),cc(iMap),NOcc,xmax,xmin,N,NZ,NP,NumC,
     $  Lfail,numbig,numsml,BB,IOut,IPrint,StoreT,LenRT,demp,cc(iBnFac),
     $  cc(jIBin),cc(jibt),cc(ibt),cc(jiscr),cc(ISizCV),NBin,Method,
     $  ifunc)
      call CheckF(cc(jip),cc(ip),X,N,0,NZ,'Pa After frmchb',5,IPrint)
C
C     Purify P if requested
C
      iend = irt
      If(Method.eq.1) then
        LNNP  = 0
        LN1NZ = 0
        LNNZ  = 0
        LNBin = 0
        jiscr1 = jip
      else if(Method.eq.2) then
        LNNP  = N*NP
        LN1NZ = N+1+NZ
        LNNZ  = N+NZ
        LNBin = N*(NBin+1)
        jiscr1 = iend
        iend   = jiscr1 + InToWP(N+NZ+1)
      else
        Call GauErr('Unknown Method in SChbEM')
        endIf
      ip2    = iend
      ip3    = ip2 + N+NZ
      iScr   = ip3 + N+NZ
      jIMap  = iScr + N+NZ
      iMap   = jiMap  + InToWP(N+1+NZ)
      jIRowA = iMap   + InToWP(NZ)
      iXRowA = jIRowA + InToWP(N*NP)
      jJNEW  = iXRowA + N*NP
      iANEW  = jJNEW + InToWP(N*NP)
      jiX    = iANEW + N*NP
      iRowNo = jiX + InToWP(LNNP)
      iscr1  = iRowNo + InToWP(N)
      jibt   = iscr1 + LNNP
      ibt    = jibt + intowp(LN1NZ)
      iRABin = ibt  + LNNZ
      iBinB  = iRABin + intowp(LNNP)
      iBinBt = iBinB + intowp(LNBin)
      iRowSr = iBinBt + intowp(LNBin)
      iscr2  = iRowSr + LNNP
      iend   = iscr2 + intowp(LNNP)
      Call CnPurf(IOut,IPrint,cc(jIP),cc(jiscr1),cc(jibt),cc(iP),
     $  cc(iP2),cc(iP3),cc(iScr),cc(ibt),cc(jIMap),cc(iMap),cc(iBnFac),
     $  cc(iSizCV),cc(iRowNo),cc(jIRowA),cc(jJNEW),cc(jiX),cc(iRABin),
     $  cc(iscr2),cc(iXRowA),cc(iANEW),cc(iscr1),cc(iRowSr),cc(iBinB),
     $  cc(iBinBt),NBin,NOcc,N,NZ,npfcy,NP,NPMax,demp,pfthr,zthr,Method)
      If(IPrint.eq.0) Write(IOut,2600) Demp
      If(NIter.gt.4.and.(Abs(Demp).gt.DempMx)) then
        Write(IOut,2500) Demp, DempMx
        Call GauErr('Stopping in SChbEM.')
        endIf
      Return
      End
*Deck SciCDR
      Subroutine SciCDR(IOut,NTot,NAtoms,ToAng,Avog,AJpCal,Hartre,Boltz,
     $  Rad,Cav1,Cav2,Disp,Rep,C,XYZ,GNorm,Wgt,NATYP,RHO,AKA,RWA,RDIFF,
     $  NT,AKT,RWT,RSolv,VMol)
      Implicit Real*8(A-H,O-Z)
      Logical GiveUp
      Dimension C(3,NAtoms), XYZ(NTot,3), GNorm(NTot,3), Wgt(NTot),
     $  AKA(NAtoms), RWA(NAtoms), AKT(NATYP), RWT(NATYP), RDIFF(NATYP),
     $  NT(NATYP)
      Save Zero, One, Two, Three, Four, Nine, Ten3, Ten8, Ten10
      Save Caill1, Caill2, Caill3, StdPre
      Data Zero/0.0d0/, One/1.0d0/, Two/2.0d0/, Three/3.0d0/,
     $  Four/4.0d0/, Nine/9.0d0/, Ten3/1.0d3/, Ten8/1.0d8/,
     $  Ten10/1.0d10/
      Data Caill1/0.214d0/, Caill2/47000.0d0/, Caill3/12.35d0/,
     $  StdPre/101325.0d0/
 9000 Format(1x,'No dispersion-repulsion data for solute atom ',i5)
C
C     StdPre = one atm in Pascals
C
C     Compute the cavitation and dispersion-repulsion energies.
C
C     Formulae from Cossi, Tomasi, Cammi, Int. J. Quant. Chem. Quant Chem
C       Symposium 29, 695-701 (1995).
C     and Cossi, Mennucci, Cammi, J. Comp. Chem. V17, 57-73, 1996.
C
C     Cavitation energy. [There is a misprint in the Cavitation paper,
C     y is wrongly defined as 4pi rho Rs / 3.]
C
C     y = 4pi rho Rs**3 / 3   rho = number density of solvent,
C                             Rs = radius of solvent as hard sphere.
C
C     Gcav      -ln(1-y) + (3y/(1-y))(Rad/Rs)
C     ---- =    +[ (3y/(1-y)) + (9/2)(y/(1-y))**2 ] * (Rad/Rs)**2
C      kT       + (yP/(rho * kT)) * (Rad/Rs)**3
C
C     Rad = radius of solute as hard sphere
C     T = absolute temp; P = pressure(atm); k = Boltzmann const
C
      Pi = Four * Atan(One)
      ToBohr = One / ToAng
C
C     Cavitation.
C     RSolv in angstrom.
C     VMol in cm**3 / mol
C     Rad (from volume of solute cavity) already in bohr.
C
      RSolv = RSolv * ToBohr
      RhoC = Avog / (VMol * (ToBohr * Ten8)**3)
      Y = Four * Pi * RhoC * RSolv**3 / Three
      YOneY = Y / (One - Y)
      Cav1 = (Y / RhoC) * (Rad / RSolv)**3
      Cav1 = Cav1 * (StdPre / Hartre) * (ToAng / Ten10)**3
      Cav2 = -log(One - Y)
     $         + (Three * YOneY * Rad / RSolv)
     $         + (Three * YOneY
     $             + (Nine / Two) * YOneY**2 ) * (Rad / RSolv)**2
      Cav2 = Cav2 * Boltz  / Hartre
C
C     Dispersion - Repulsion.
C
C     Convert point positions and weights to angstroms,
C     in place since we will just throw values away after.
C
      Do 10 I = 1, NTot
        XYZ(I,1) = XYZ(I,1) * ToAng
        XYZ(I,2) = XYZ(I,2) * ToAng
        XYZ(I,3) = XYZ(I,3) * ToAng
   10   Wgt(I) = Wgt(I) * ToAng * ToAng
      Disp = Zero
      Rep  = Zero
      GiveUp = .false.
      Do 20 IAt = 1, NAtoms
        If(AKA(IAt).eq.Zero.or.RWA(IAt).eq.Zero) then
          Write(IOut,9000)IAt
          GiveUp = .true.
          endif
   20   Continue
      Do 30 ISol = 1, NATyp
        If(AKT(ISol).eq.Zero.or.RWT(ISol).eq.Zero) then
          Write(IOut,*)
     $      'No dispersion-repulsion data for a solvent atom'
          GiveUp = .true.
          endif
   30   Continue
      If(GiveUp) then
        Write(Iout,*)'Missing dispersion-repulsion data, skipping'
        Return
        endif
      Do 50 ISol = 1, NATYP
        TmpDis = Zero
        TmpRep = Zero
        Do 40 IAt = 1, NAtoms
C         Magic incantations from Caillet.
          RWAT = Two * Sqrt(RWA(IAt) * RWT(ISol))
          DAT = Caill1 * AKA(IAt) * AKT(ISol) * RWAT**6
          BAT = Caill2 * AKA(IAt) * AKT(ISol)
          ALPHAT = Caill3 / RWAT
          Cx = C(1,IAt) * ToAng
          Cy = C(2,IAt) * ToAng
          Cz = C(3,IAt) * ToAng
          Do 40 I = 1, NTot
C           Cavity is expanded to 'solvent accessible surface'
C           we just elongate surface along normal.
            Xc = XYZ(I,1) + GNorm(I,1) * RDIFF(ISol)
            Yc = XYZ(I,2) + GNorm(I,2) * RDIFF(ISol)
            Zc = XYZ(I,3) + GNorm(I,3) * RDIFF(ISol)
            Xd = Xc - Cx
            Yd = Yc - Cy
            Zd = Zc - Cz
            Rai2 = Xd**2 + Yd**2 + Zd**2
            Rai  = Sqrt(Rai2)
            Rai6 = Rai2 * Rai2 * Rai2
            RdotN = Xd*GNorm(I,1) + Yd*GNorm(I,2) + Zd*GNorm(I,3)
C
            ADisp = -DAT / (Three * Rai6)
            TmpDis = TmpDis + Wgt(I) * ADisp * RdotN
C
            S1 = ALPHAT * Abs(Rai)
            S  = One/S1 + Two/(S1*S1) + Two/(S1*S1*S1)
            ARep  = BAT * Exp(-S1) * S
   40       TmpRep = TmpRep + Wgt(I) * ARep * RdotN
        Disp = Disp + TmpDis * Rho * NT(ISol)
   50   Rep  = Rep  + TmpRep * Rho * NT(ISol)
C
C     Convert from kcal / mol to au.
C
      Disp = Disp * ((AJpCal * Ten3) / (Hartre * Avog))
      Rep = Rep * ((AJpCal * Ten3) / (Hartre * Avog))
      Return
      End
*Deck SciDat
      Subroutine SciDat(ISlTyp,MxATYp,RSolv,VMol,Rho,NATyp,NT,RDiff,AKT,
     $  RWT)
      Implicit Real*8(A-H,O-Z)
      Logical DonH, AccH
      Dimension NT(*), RDiff(*), AKT(*), RWT(*)
C
C     Database of optical and physical data for various solvent
C     Modified from Datasol routine of Cossi, et. al. by
C     S. Clifford to work with SCIPCM.
C
C     RSolv = radius of solvent as hard sphere, A
C     VMol = molar volume, cm**3 / mol
C     Rho  = number density, 1 / Ang**3
C     NATyp = Different atom types in solvent. (Seems too low for most cases).
C     NT(I) = Number of atoms of type I
C     RDiff = Distance to extend cavity by for dispersion-repulsion, A
C     AKT, RWT = Caillet-Claverie atomic parameters.
C
C     EPS, EPSINF, DEREPS, TCS, STEN, DSTEN, CMF, DonH, AccH, are not
C     currently used from this routine.
C
      If(ISlTyp.eq.1) then
C
C       WATER solvent
C
        EPS = 78.39d0
        EPSINF = 1.776d0
        DEREPS = -0.3562d0
        RSolv = 1.385d0
        VMol = 18.07d0
        TCE = 2.57d-4
        STEN = 71.81d0
        DSTEN = 0.650d0
        CMF = 1.277d0
        DonH = .True.
        AccH = .True.
C       Atomic parameters for dispersion and repulsion
        Rho = 3.348d-2
        NATyp = 2
        If(NATyp.gt.MxATyp) goto 999
C       Atom O:
        NT(1) = 1
        RDiff(1) = 1.5d0
        AKT(1) = DKList(8)
        RWT(1) = RList(8)
C       Atom H:
        NT(2) = 2
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C
      else If(ISlTyp.eq.19) then
C
C       ACETONITRILE solvent
C
        EPS = 36.64d0
        EPSINF = 1.806d0
C Warning: RSolv and VMol for this solvent have to be checked
        RSolv = 2.155d0
        VMol = 53.68d0
        TCE = 1.192d-3
        STEN = 36.47d0
        DSTEN = 1.373d0
        CMF = 0.808d0
C       Atomic parameters for dispersion and repulsion
        Rho = 1.153d-2
        NATyp = 3
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 2
        RDiff(1) = 1.76d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 3
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C       Atom N:
        NT(3) = 1
        RDiff(3) = 1.6d0
        AKT(3) = DKList(7)
        RWT(3) = RList(7)
C
      else If(ISlTyp.eq.4) then
C
C       METHANOL solvent
C
        EPS = 32.63d0
        EPSINF = 1.758d0
        DEREPS = -0.1984d0
        RSolv = 1.855d0
        VMol = 40.7d0
        TCE = 1.182d-3
        STEN = 22.12d0
        DSTEN = 1.154d0
        CMF = 1.776d0
        DonH = .True.
        AccH = .True.
C       Atomic parameters for dispersion and repulsion
        Rho = 1.495d-2
        NATyp = 3
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 1
        RDiff(1) = 1.76d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 4
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C       Atom O:
        NT(3) = 1
        RDiff(3) = 1.5d0
        AKT(3) = DKList(8)
        RWT(3) = RList(8)
C
      else If(ISlTyp.eq.5) then
C
C       ETHANOL solvent
C
        EPS = 24.55d0
        EPSINF = 1.847d0
        DEREPS = -0.1510d0
        RSolv = 2.180d0
        VMol = 58.7d0
        TCE = 1.103d-3
        STEN = 21.89d0
        DSTEN = 1.146d0
        CMF = 1.543d0
        DonH = .True.
        AccH = .True.
C       Atomic parameters for dispersion and repulsion
        Rho = 1.032d-2
        NATyp = 3
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 2
        RDiff(1) = 1.76d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 6
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C       Atom O:
        RDiff(3) = 1.5d0
        AKT(3) = DKList(8)
        RWT(3) = RList(8)
C
      else If(ISlTyp.eq.12) then
C
C       CHLOROFORM Solvent
C
        EPS = 4.90d0
        EPSINF = 2.085d0
        RSolv = 2.48d0
        VMol = 80.7d0
        TCE = 1.255d-3
        STEN = 26.53d0
        DSTEN = 0.d0
        CMF = 0.d0
C       Atomic parameters for dispersion and repulsion
        Rho = 7.482d-3
        NATyp = 3
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 1
        RDiff(1) = 2.82d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 1
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C       Atom Cl:
        NT(3) = 3
        RDiff(3) = 1.79d0
        AKT(3) = DKList(17)
        RWT(3) = RList(17)
C
      else If(ISlTyp.eq.13) then
C
C       ETHYL ETHER solvent
C
        EPS = 4.335d0
        EPSINF = 0.000d0
        RSolv = 2.785d0
        VMol = 103.84d0
        TCE = 1.617d-3
        STEN = 0.00d0
        DSTEN = 0.d0
        CMF = 0.d0
        DonH = .False.
        AccH = .True.
C       Atomic parameters for dispersion and repulsion
        Rho = 5.799d-3
        NATyp = 3
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 4
        RDiff(1) = 1.76d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 10
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C       Atom O:
        NT(3) = 1
        RDiff(3) = 1.50d0
        AKT(3) = DKList(8)
        RWT(3) = RList(8)
C
      else If(ISlTyp.eq.8) then
C
C       CH2Cl2 solvent
C
        EPS = 8.93d0
        EPSINF = 2.020d0
        RSolv = 2.27d0
        VMol = 64.5d0
        TCE = 1.367d-3
        STEN = 27.33d0
        DSTEN = 0.d0
        CMF = 0.d0
C       Atomic parameters for dispersion and repulsion
        Rho = 9.74d-3
        NATyp = 3
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 2
        RDiff(1) = 1.76d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 2
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C       Atom Cl:
        NT(3) = 2
        RDiff(3) = 1.79d0
        AKT(3) = DKList(17)
        RWT(3) = RList(17)
C
      else If(ISlTyp.eq.7) then
C
C       1,2-DICHLOROETHANE solvent
C
        EPS = 10.36d0
        EPSINF = 2.085d0
        RSolv = 2.505d0
        VMol = 79.4d0
        TCE = 1.156d-3
        STEN = 31.54d0
        DSTEN = 0.d0
        CMF = 0.d0
C       Atomic parameters for dispersion and repulsion
        Rho = 7.576d-3
        NATyp = 3
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 2
        RDiff(1) = 1.76d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 4
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C       Atom Cl:
        NT(3) = 2
        RDiff(3) = 1.79d0
        AKT(3) = DKList(17)
        RWT(3) = RList(17)
C
      else If(ISlTyp.eq.16) then
C
C       CCl4 solvent
C
        EPS = 2.228d0
        EPSINF = 2.129d0
        RSolv = 2.685d0
        VMol = 96.5d0
        TCE = 1.270d-3
        STEN = 26.15d0
        DSTEN = 1.436d0
        CMF = 0.629d0
C       Atomic parameters for dispersion and repulsion
        Rho = 6.241d-3
        NATyp = 2
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 1
        RDiff(1) = 2.82d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom Cl:
        NT(2) = 4
        RDiff(2) = 1.79d0
        AKT(2) = DKList(17)
        RWT(2) = RList(17)
C
      else If(ISlTyp.eq.15) then
C
C       BENZENE solvent
C
        EPS = 2.247d0
        EPSINF = 2.244d0
        RSolv = 2.63d0
        VMol = 88.91d0
        TCE = 1.380d-3
        STEN = 28.18d0
        DSTEN = 1.469d0
        CMF = 0.629d0
C       Atomic parameters for dispersion and repulsion
        Rho = 6.773d-3
        NATyp = 2
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 6
        RDiff(1) = 1.5d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 6
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C
      else If(ISlTyp.eq.14) then
C
C       TOLUENE solvent
C
        EPS = 2.379d0
        EPSINF = 2.232d0
        RSolv = 2.82d0
        VMol = 106.3d0
        TCE = 1.08d-3
        STEN = 27.92d0
        DSTEN = 1.391d0
        CMF = 0.679d0
C       Atomic parameters for dispersion and repulsion
        Rho = 5.665d-3
        NATyp = 2
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 7
        RDiff(1) = 1.5d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 8
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C
      else If(ISlTyp.eq.11) then
C
C       CHLOROBENZENE solvent
C
        EPS = 5.621d0
        EPSINF = 2.320d0
        RSolv = 2.805d0
        VMol = 101.79d0
        TCE = 0.981d-3
        STEN = 32.69d0
        DSTEN = 0.d0
        CMF = 0.d0
C       Atomic parameters for dispersion and repulsion
        Rho = 5.916d-3
        NATyp = 3
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 6
        RDiff(1) = 1.5d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 5
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C       Atom Cl:
        NT(3) = 1
        RDiff(3) = 1.79d0
        AKT(3) = DKList(17)
        RWT(3) = RList(17)
C
      else If(ISlTyp.eq.3) then
C
C       NITROMETHANE solvent
C
        EPS = 38.20d0
        EPSINF = 1.904d0
        RSolv = 2.155d0
        VMol = 53.68d0
        TCE = 1.192d-3
        STEN = 36.47d0
        DSTEN = 1.373d0
        CMF = 0.808d0
C       Atomic parameters for dispersion and repulsion
        Rho = 1.122d-2
        NATyp = 4
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 1
        RDiff(1) = 1.76d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 3
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C       Atom N:
        NT(3) = 1
        RDiff(3) = 1.6d0
        AKT(3) = DKList(7)
        RWT(3) = RList(7)
C       Atom O:
        NT(4) = 2
        RDiff(4) = 1.5d0
        AKT(4) = DKList(8)
        RWT(4) = RList(8)
C
      else If(IslTyp.eq.18) then
C
C       N-HEPTANE solvent
C
        EPS = 1.92d0
        EPSINF = 1.918d0
        RSolv = 3.125d0
        VMol = 146.56d0
        TCE = 1.25d-3
        STEN = 19.80d0
        DSTEN = 1.505d0
        CMF = 0.687d0
C       Atomic parameters for dispersion and repulsion
        Rho = 4.109d-3
        NATyp = 2
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 7
        RDiff(1) = 1.76d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 16
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C
      else If(ISlTyp.eq.17) then
C
C       CYCLOHEXANE solvent
C
        EPS = 2.023d0
        EPSINF = 2.028d0
        RSolv = 2.815d0
        VMol = 108.10d0
        TCE = 1.20d-3
        STEN = 24.38d0
        DSTEN = 1.467d0
        CMF = 0.621d0
C       Atomic parameters for dispersion and repulsion
        Rho = 5.571d-3
        NATyp = 2
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 6
        RDiff(1) = 1.76d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 12
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C
      else If(ISlTyp.eq.10) then
C
C       ANILINE solvent
C
        EPS = 6.89d0
        EPSINF = 2.506d0
        RSolv = 2.80d0
        VMol = 91.15d0
        TCE = 0.85d-3
        STEN = 42.79d0
        DSTEN = 0.731d0
        CMF = 0.972d0
C       Atomic parameters for dispersion and repulsion
        Rho = 6.607d-3
        NATyp = 3
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 6
        RDiff(1) = 1.5d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 7
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C       Atom N:
        NT(3) = 1
        RDiff(3) = 1.6d0
        AKT(3) = DKList(7)
        RWT(3) = RList(7)
C
      else If(ISlTyp.eq.6) then
C
C       ACETONE solvent
C
        EPS = 20.7d0
        EPSINF = 1.841d0
        DEREPS = -9.77d-2
        RSolv = 2.38d0
        VMol = 73.52d0
        TCE = 1.42d-3
        STEN = 22.67d0
        DSTEN = 0.d0
        CMF = 0.d0
C       Atomic parameters for dispersion and repulsion
        Rho = 8.190d-3
        NATyp = 3
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 3
        RDiff(1) = 1.76d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 6
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C       Atom O:
        NT(3) = 1
        RDiff(3) = 1.5d0
        AKT(3) = DKList(8)
        RWT(3) = RList(8)
C
      else If(ISlTyp.eq.9) then
C
C       TETRAHYDROFURAN solvent
C
        EPS = 7.58d0
        EPSINF = 1.971d0
        RSolv = 2.56d0
        VMol = 81.11d0
        TCE = 1.142d-3
        STEN = 26.40d0
        DSTEN = 0.d0
        CMF = 0.d0
C       Atomic parameters for dispersion and repulsion
        Rho = 7.425d-3
        NATyp = 3
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 4
        RDiff(1) = 1.76d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 8
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C       Atom O:
        NT(3) = 1
        RDiff(3) = 1.5d0
        AKT(3) = DKList(8)
        RWT(3) = RList(8)
C
      else If(ISlTyp.eq.2) then
C
C       DIMETHYLSULFOXIDE solvent
C
        EPS = 46.7d0
        EPSINF = 2.179d0
        DEREPS = -0.1902d0
        RSolv = 2.455d0
        VMol = 70.94d0
        TCE = 9.82d-2
        STEN = 42.86d0
        DSTEN = 0.d0
        CMF = 0.d0
C       Atomic parameters for dispersion and repulsion
        Rho = 8.49d-3
        NATyp = 4
        If(NATyp.gt.MxATyp) goto 999
C       Atom C:
        NT(1) = 2
        RDiff(1) = 1.76d0
        AKT(1) = DKList(6)
        RWT(1) = RList(6)
C       Atom H:
        NT(2) = 6
        RDiff(2) = 1.2d0
        AKT(2) = DKList(1)
        RWT(2) = RList(1)
C       Atom S:
        NT(3) = 1
        RDiff(3) = 1.80d0
        AKT(3) = DKList(16)
        RWT(3) = RList(16)
C       Atom O:
        NT(4) = 1
        RDiff(4) = 1.5d0
        AKT(4) = DKList(8)
        RWT(4) = RList(8)
      else
        call GauErr('SciDat: built-in solvent is not known.')
        endif
      Return
  999 Call GauErr('Too many atom types in built-in solvent')
      Return
      End
*Deck SCIFoc
      Subroutine SCIFoc(IOut,IPrint,IScal,NTot,NAtoms,NBas6D,C,xyz,
     $  Coord,gnorm,Weight,ffact,EFG,PTot,ISc1,qtot,Alpha,bmf,Beta,q,
     $  GFact,Ell,Em,Scr1,SolF,Cutoff,MoBas,IOpCl,NMO,CMO,IPFlag,AllowP,
     $  MaxTyp,NTotR,dwidrn,XYZOld,Rhatx,Rhaty,Rhatz,ISlMap,ISubst,
     $  Const1,Cons1c,V,MDV)
      Implicit Real*8 (A-H,O-Z)
      Logical MoBas,AllowP(*)
      Dimension C(3,NAtoms),xyz(NTot,3),Coord(3,NTot),
     $  Gnorm(NTot,3),Weight(NTot),FFact(NTot),EFG(NTot,10,ISc1),
     $  PTot(*),bmf(NTot,ISc1),
     $  V(*),qtot(NTot),SolF(*),CMO(NBas6D,NMO),JJ(1),XX(1)
      Dimension q(NTot,ISc1*ISc1),Alpha(NTot,2),Beta(*),
     $  GFact(*),Em(Ntot,3),Ell(NTot),dwidrn(NTot,2),
     $  XYZOld(NTot,3),Rhatx(*),Rhaty(*),Rhatz(*),
     $  Scr1(NTotR),ISlMap(NTot),ISubst(NTot)
      Save JJ,XX
      Data JJ/0/,XX/0/
C
      NTT6D = (NBas6D*(NBas6D+1))/2
      Call AClear(NTT6D,SolF)
C
C     Cavity terms. (dE/dr . dr/dP) & (dE/dn . dn/dP)
C
      Call TStamp(1,'Before SurDer')
      If(IScal.eq.0.or.IScal.eq.2) then
        Call SurDer(IPrint,NTot,NTotR,xyz,coord,Gnorm,qtot,bmf,Nbas6D,
     $    Ntt6D,const1,Rhatx,Rhaty,Rhatz,PTot,CutOff,SolF,Ffact,Weight,
     $    MoBas,IOpCl,NMO,CMO,Q,ISc1,IScal,Alpha,Beta,GFact,Cons1c,EFG,
     $    Ell,Em,dwidrn,XYZOld,ISlMap,ISubst,MaxTyp,V,MDV)
      else
        Call GauErr('Terms missing for scaled charges.')
        endIf
      Call TStamp(1,'After SurDer')
c
c     Fixed cavity terms (dE/dPmunu)c.
c
      Do 40 J = 1,3
        Do 40 I = 1, NTotR
   40     Coord(J,I) = xyz(I,J)
      IWt = 1
      IV = IWt + 4*NTotR
      MDV1 = MDV - IV + 1
      Call AMove(NTotR,BMF,scr1)
      Call AMove(NTotR,qtot,V(IWt))
      If(IScal.eq.2) then
        Call ASub(NTotR,V(IWt),Beta,V(IWt))
        Call ASub(NTotR,scr1,Alpha,Scr1)
        endIf
      Call ANeg(NTotR,Scr1,Scr1)
      Call AMul(NTotR,Scr1,GNorm(1,1),V(IWt+3*NTotR))
      Call AMul(NTotR,Scr1,GNorm(1,2),V(IWt+2*NTotR))
      Call AMul(NTotR,Scr1,GNorm(1,3),V(IWt+1*NTotR))
      Call ANeg(NTotR,Scr1,Scr1)
      Call TStamp(1,'Before OneElI No. 2')
      Call OneElI(IOut,IPrint,1,6,220,0,0,0,1,1,NAtoms,JJ,XX,XX,NBas6D,
     $  1,JJ(1),JJ,C,.False.,JJ,V(IWt),NAtoms,JJ,Coord,NTotR,XX,SolF,
     $  CutOff,0,1,JJ,XX,IPFlag,AllowP,.False.,JJ,XX,V(IV),V(IV),MDV1)
      Call TStamp(1,'After OneElI No. 2')
      Return
      End
*Deck SciPrt
      Subroutine SciPrt(IOut,IDoWot,IPCM,NAtoms,C,IAn,ToAng,Avog,AJpCal,
     $  Hartre,Boltz,V,MDV)
      Implicit Real*8(A-H,O-Z)
      Parameter (IRwSCI=561)
      Parameter (MxATyp=10)
      Dimension C(3,NAtoms),IAn(NAtoms),V(*)
C
 9000 Format(1x,'---- SCIPCM Solvation Results (all in au) ----')
 9010 Format(  ' Total "Solvent Accessible Surface Area" of solute ',
     $  'cavity =    ',D16.8,/,
     $  ' Volume of solute cavity = ',36x,D16.8,/,
     $  ' (1/(4*pi)) * flux of nuclear E-field thru cavity ',
     $  'surface =    ',1D16.8,/,
     $  ' (1/(4*pi)) * flux of solute E-field thru cavity ',
     $  'surface =     ',1D16.8)
 9020 Format(' (1/(1-1/epsi)) * total surface polarization charge = ',
     $  9x,1D16.8,/,' (1/2) * solvent-solute interaction energy = ',18x,
     $  1D16.8)
 9030 Format(' (1/(1-1/epsi)) * unscaled total surface polarization ',
     $  'charge = ',1D16.8,/,
     $  ' (1/(1-1/epsi)) * scaled total surface polarization ',
     $  'charge =   ',1D16.8,/,
     $  ' (1/2) * solvent-solute interaction energy = ',18x,1D16.8)
 9040 Format(1x,'Cavitation term: Cossi, Tomasi, Cammi, IJQC',
     $  ' Q Chem Symp 29, 695-701 1995',/,1x,'G(cav) = ',D16.8,
     $  ' *P(in atm) + ',D16.8,'*T(in K) au',/,1x,'G(cav) = ',D16.8,
     $  ' au at ',F6.3,'atm, ',F7.2,'K.')
 9050 Format(1x,'Dispersion-Repulsion term: Cossi, Mennucci, Cammi,',
     $  ' J Comp Chem 17, 57-73 1996',/,1x,'G(disp) = ',D16.8,
     $  ' au,  G(Rep) = ',D16.8,' au')
 9060 Format(1x,'Total (Cav+Disp+Rep): ',SP,D16.8,' au')
      Save One, Three, Four, Temp, Press
      Data One/1.0d0/, Three/3.0d0/, Four/4.0d0/, Temp/298.15/,
     $  Press/1.00d0/
C
C     IDoWot = 0        Just print energy, flux, etc.
C     IDoWot = 1        Compute and print Dispersion and Cavitation as well.
C
      If(IDoWot.ge.1) Write(IOut,9000)
      Pi = Four * ATan(One)
      ISlTyp = IAbs(IPCM)
      Thrd = One / Three
      Call FileIO(2,-IRwSCI,1,XMask,0)
      Call Bnupak(XMask,NTotR,NTot)
      Call FileIO(2,IRwSCI,1,XMask,0)
      Call Bnupak(XMask,IScal,ISurf)
      Call Fileio(2,IRwSCI,1,ESol,0)
      Call Fileio(2,IRwSCI,1,GaussN,0)
      Call Fileio(2,IRwSCI,1,Gauss,0)
      Call Fileio(2,IRwSCI,1,PolTot,0)
      Call Fileio(2,IRwSCI,1,PolTtU,0)
      Call Fileio(2,IRwSCI,1,SAreaT,0)
      Call Fileio(2,IRwSCI,1,Volume,0)
      Write(IOut,9010) SAreaT, Volume, GaussN, Gauss
      If(IScal.eq.2) Write(IOut,9030) PolTtU, PolTot, ESol
      If(IScal.eq.0) Write(IOut,9020) PolTot, ESol
C
      If(IDoWot.ge.1.and.ISlTyp.ne.0) then
        IXYZ = 1
        IGNorm = IXYZ + 3*NTotR
        IWgt = IGNorm + 3*NTotR
        IiNT = IWgt + NTotR
        IAKT = IiNT + MxATyp
        IRWT = IAKT + MxATyp
        IRDIFF = IRWT + MxATyp
        IV = IRDIFF + MxATyp
        If(MDV-IV+1.lt.0) Call GauErr('Out of memory in SciPrt')
        Call FileIO(2,IRwSCI,3*NTotR,V(IXYZ),0)
        Call FileIO(2,IRwSCI,3*NTotR,V(IGNorm),0)
        Call FileIO(2,IRwSCI,NTotR,V(IWgt),2*NTotR+7*NTot)
C       Parameters for cavitation and dispersion-repulsion.
        Call SciDat(ISlTyp,MxATyp,RSolv,VMol,Rho,NATyp,V(IiNT),
     $    V(IRDiff),V(IAKT),V(IRWT))
C       Parameters for solute
        Rad = ((Three / (Four * Pi)) * Volume) ** Thrd
        IKA = IV
        IRWA = IKA + NAtoms
        IVV = IRWA + NAtoms
        If(MDV-IVV+1.lt.0) Call GauErr('Out of mem in SciPrt')
        Do 10 IA = 1, NAtoms
          V(IKA+IA-1) = DKLIST(Ian(IA))
   10     V(IRWA+IA-1) = RLIST(Ian(IA))
        Call SciCDR(IOut,NTotR,NAtoms,ToAng,Avog,AJpCal,Hartre,Boltz,
     $    Rad,Cav1,Cav2,Disp,Rep,C,V(IXYZ),V(IGNorm),V(IWgt),NATyp,Rho,
     $    V(IKA),V(IRWA),V(IRDiff),V(IiNT),V(IAKT),V(IRWT),RSolv,VMol)
        Cavit = Cav1*Press+Cav2*Temp
        Write(IOut,9040)Cav1,Cav2,Cavit,Press,Temp
        Write(IOut,9050)Disp,Rep
        Write(IOut,9060)Cavit+Disp+Rep
        endif
      Return
      End
*Deck SCISCF
      Subroutine SCISCF(IOut,IPrint,IRwPA,IRwPB,IRwCA,IRwCB,JCycle,
     $  NBasis,NBas6D,DA,DB,FA,FB,NAtoms,AtmChg,C,IOpCl,NAE,NBE,ESol,
     $  IPass,DoSurf,IRwSrf,RdPCM,Epsi,AccDes,SymCut,IPFlag,AllowP,
     $  MaxTyp,V,MDV)
      Implicit Real*8 (A-H,O-Z)
C
C     Main routine for SCF (HF or DFT) polarizable continuum (PC)
C     energy calculations.  Calls SCISrf to get coordinates, normal
C     vectors and weights of integration points on an electron
C     density isosurface - which defines the solute cavity.  Then
C     calls SCIFoc to get PC interaction energy and PC contribution
C     to Fock matrices.
C
C     Input:  Epsi   Cont   IPhi ITheta ISurf NoPrun IScal
C     e.g. :  78.3  0.0004   40    20     3     F      2
C     ---------------------------------------------------
C     Epsi = Dielectric constant of solvent
C     Cont = Value of isodensity surface in atomic units.  0.0004 au
C            is recommended and the default.
C     IPhi = Number of phi values for numerical surface integration,
C            or number of points in a Lebedev grid.  Default is 302.
C     ITheta = Number of theta Values for numerical surface integration,
C              1 to use a Lebedev grid, which is the default.
C     ISurf = Whether to do surface integrals using a single origin (2)
C             or piecewise, using Becke grids (3,4,5).  3 is all nuclei
C             as origins, 4 is non-H atoms as origins, 5 means another
C             card follows this one with a free form list of atom numbers
C             for origins. 4 is the default.
C     NoPrun = Whether to not prune out close points in conductor scaled
C              calculations. Default is False. Internally DoPrun = .not.NoPrun
C     IScal = Whether to scale the polarization charges.
C             0=no.
C             1=Tomasi scaling, not implemented.
C             2=yes, scale the polarization charges so that the
C               TOTAL potential at all points on the cavity surface is
C               zero when the medium is a conductor (epsi=infinity).
C             Default is 2.  Explicit derivatives of the scale factors
C             are available for 0 and 2 only.
C
C     DoSurf is .True. to compute surface.  It is saved on IRwSrf,
C     which may be set or reallocated, unless IRwSrf=-1.
C
C     T.A. Keith and M.J. Frisch 1994
C     S. Clifford 1997-99
C
      Logical MoBas,UseOld,AllowP(*),DoPrnt,DoGrad,RdInpt,DoSurf,RdPCM,
     $  DoPrun,LJJ(1),NDPRun
      Parameter (IRwSCI=561,IRwSC2=573,IRwGen=501,LSc2=7)
      Dimension AtmChg(Natoms),C(3,Natoms),DA(*),DB(*),V(MDV),XX(1),
     $  JJ(1),FA(*),FB(*)
      Equivalence (IDPRun,NDPRun)
      Save Zero,One,Two,Four,XX,JJ,Pt5,Pt01,LJJ
      Data Zero/0.0d0/,One/1.0d0/,Two/2.0d0/,Four/4.0d0/,XX/0.d0/,JJ/0/,
     $  Pt5/0.5d0/,Pt01/0.01d0/,LJJ/.False./
 9000 Format(' Dielectric constant of solvent = ', 1F12.6,/,
     $  ' Solute cavity defined by isodensity surface of value =',
     $  1F12.6,' au')
 9010 Format(' Using a single origin to do surface integrals.')
 9020 Format(' Using Becke method to do surface integrals',/,
     $'  - all nuclei will be used as origins.')
 9030 Format(' Using Becke method to do surface integrals',/,
     $'  - all nuclei except protons will be used as origins.')
 9040 Format(' Using Becke method to do surface integrals',/,
     $'  - user-selected nuclei will be used as origins.')
 9050 Format(' Illegal surface integration method ',I10,' in SCISCF.')
 9055 Format(' Using old-style function for Becke weights.')
 9057 Format(' Using exponential function for Becke weights.')
 9058 Format(' Unknown function type for Becke weights.')
 9060 Format(' Illegal scaling method',I10,' in SCISCF.')
 9070 Format(' Compute SCI-PCM surface.')
 9080 Format(' Recover previous SCI-PCM surface.')
 9090 Format(' Short by ',I10,' words.')
 9100 Format(' Polarization charges will be uniformly scaled.')
 9110 Format(' Polarization charges will be scaled so the total',
     $  ' potential',/,'   on the cavity surface equals zero in a',
     $  ' conducting medium.')
 9120 Format(' Close points will be merged for conductor scaling.')
 9130 Format(' The A matrix elements will be calculated on the fly.')
 9140 Format(' The distance matrix will be stored in memory.')
 9150 Format(' The A matrix will be stored in memory.')
 9160 Format(' The A and C matrices will both be stored in memory.')
 9170 Format(' Need',I6,'MW to store dists',I6,
     $  'MW to store one matrix',I6,'MW for both.')
 9180 Format(/,' WARNING!  Serious error in surface integrals.',/,
     $  ' Nuclear flux = ',f9.2,' Qnuc =',f9.2,
     $  ' Error in int =',f9.2,/,
     $  ' It is probable that some of the solute is outside ',
     $  'the cavity and/or',/,
     $  ' parts of the cavity surface cannot be reached from the ',
     $  'origin.',/,
     $  ' Try more integration points or a different set of ',
     $  'integration origins.',/)
C
      NTT = (NBasis*(NBasis+1))/2
      NTT6D = (NBas6D*(NBas6D+1))/2
      Pi = Four*ATan(One)
      If(AccDes.eq.Zero) then
        CutOff = GFloat(10)**(-10)
      else
        CutOff = AccDes / GFloat(10)
        endIf
      EfCrit = CutOff * GFloat(10)
C     NB - IRwSCI & IRwSC2 used in SciPrt
      RdInpt = ITqry(IRwSC2).le.0
      IIOrg = 1
      Call GetDat(IOut,RdInpt,IRwSC2,0,NAtoms,IScal,Epsi,Cont,ITheta,
     $  IPhi,IBWFnc,RadT,NumR,ISurf,MoBas,IPass,CutPt,NTot0,RdPCM,
     $  AtmChg,V(IIOrg),NOrg,DoPrun,PtThS,PtInc)
      PtThr = PtThS
      LScSl2 = LSc2 + InToWp(NOrg)
      Call FileIO(1,-IRwGen,1,-Epsi,5)
      Call FileIO(1,-IRwGen,1,Cont,6)
      DoPrnt = IPrint.gt.1.or.(JCycle.eq.1.and.IPSave(0).eq.0)
      NAng = IPhi*ITheta
      MTot = NAng
      If(ISurf.gt.2) MTot = MTot * NOrg
      If(DoPrnt) then
        Write(IOut,9000) Epsi, Cont
        If(ISurf.eq.2) then
          Write(IOut,9010)
        else if(ISurf.eq.3) then
          Write(IOut,9020)
        else if(ISurf.eq.4) then
          Write(IOut,9030)
        else if(ISurf.eq.5) then
          Write(IOut,9040)
        else
          Write(IOut,9050) ISurf
          Call Lnk1E(0)
          endIf
        If(ISurf.gt.2) then
          If(IBWFnc.eq.1) then
            Write(IOut,9055)
          else If(IBWFnc.eq.2) then
            Write(IOut,9057)
          else
            Write(IOut,9058) IBWFnc
            endif
          endif
        Endif
      If(IScal.ne.0.and.IScal.ne.2) then
        Write(IOut,9060) IScal
        Call Lnk1E(0)
        endIf
C
C     Prepare Cartesian total density matrix and possibly orbitals.
C
      Call FileIO(2,-IRwPA,NTT,DA,0)
      Call FileIO(2,-IRwPB,NTT,DB,0)
      Call AAdd(NTT,DA,DB,DA)
      MOBas = .False.
      If(MOBas) then
        NMO = NAE + IOpCl*NBE
      else
        NMO = 0
        endIf
      ICMO = IIOrg + InToWp(NOrg)
      IV = ICMO + NMO*NBas6D
      MDV1 = MDV - IV + 1
      Call TstCor(IV,MDV,'SCISCF-1')
      If(MOBas) then
        Call FileIO(2,-IRwCA,NBasis*NAE,V(ICMO),0)
        If(IOpCl.ne.0)
     $    Call FileIO(2,-IRwCB,NBasis*NBE,V(ICMO+NBasis*NAE),0)
        endIf
      If(NMO.gt.0.and.NBasis.ne.NBas6D)
     $  Call UnPcck(0,Junk,V(ICMO),NBas6D,NMO,NBasis,NMO)
      Call P2Cart(IOut,IPrint,.False.,.True.,.False.,.False.,MOBas,
     $  .False.,.False.,.False.,.True.,NBasis,NBas6D,1,1,1,NMO,IPureD,
     $  IPureF,DA,XX,XX,XX,XX,XX,JJ,JJ,V(ICMO),LJJ,LJJ,V(IV),MDV)
C
C     Generate the surface.
C
      DoGrad = .True.
      If(DoSurf) then
        If(IPSave(0).eq.0) Write(IOut,9070)
        IIOrg = 1
        ICOrg = IIOrg + InToWp(NOrg)
        Ixyzk = ICOrg + 3*NOrg
        IGnorm = Ixyzk+3*MTot
        Iwgt = IGnorm+3*MTot
        Icen = Iwgt+MTot
        Irad = Icen+IntoWp(MTot)
        Iunit = irad+NOrg
        Ixyz = Iunit+3*NAng
        Iunitt = Ixyz+3*NAng
        Ixyzs = Iunitt+3*NAng
        Imx = Ixyzs+3*NAng
        Irho = Imx+IntoWp(NAng)
        IaWgt = Irho+NAng
        Iiterm = Iawgt+NAng
        IGRho = Iiterm+IntoWP(NAng)
        IPBeck = IGRho+3*NAng
        ISArea = IPBeck+NOrg*NAng
        IV = ISArea + NOrg
        MDV1 = MDV-IV+1
        Call TstCor(IV,MDV,'SCISCF')
        Call TStamp(1,'Before SCISrf')
        Call SCISrf(IOut,Iprint,DoPrnt,NAtoms,NBas6D,NTT6D,Cont,C,
     $    AtmChg,ITheta,IPhi,NAng,MTot,NTot,RadT,Numr,Cutoff,CutPt,
     $    IBWFnc,MoBas,NMO,IOpCl,ISurf,NOrg,V(IIOrg),V(ICOrg),V(Ixyzk),
     $    V(IGnorm),V(Iwgt),DA,V(Icen),V(Irad),V(Iunit),V(Ixyz),
     $    V(Iunitt),V(Ixyzs),V(Imx),V(Irho),V(Iawgt),V(Iiterm),V(IGRho),
     $    V(IPBeck),V(ICMO),V(ISArea),SAreaT,Volume,V(IV),MDV1)
        Call TStamp(1,'After SCISrf')
        UseOld = (NTot.eq.NTot0)
        MaxSFa = NTot/100
        MaxSFa = 0
        If(IScal.eq.0) KNTot = 2
        If(IScal.eq.2) KNTot = 10
        Call ConDDF(IRwSC2,LScSl2+KNTot*NTot)
C
        IIOrg = 1
        ICOrg = IIOrg + InToWp(NOrg)
        ICMO = ICOrg + 3*NOrg
        Ixyz = ICMO + NMO*NBas6D
        Call AMove1(NTot,Ixyzk-1,Ixyz-1,V)
        Call AMove1(NTot,Ixyzk+MTot-1,NTot+Ixyz-1,V)
        Call AMove1(NTot,Ixyzk+2*MTot-1,2*NTot+Ixyz-1,V)
        Call AMove1(NTot,IGnorm-1,Ixyz+3*NTot-1,V)
        Call AMove1(NTot,IGnorm+MTot-1,Ixyz+3*NTot+NTot-1,V)
        Call AMove1(NTot,IGnorm+2*MTot-1,Ixyz+3*NTot+2*NTot-1,V)
        Ignorm = Ixyz+3*NTot
        Call AMove1(NTot,Iwgt-1,IGnorm+3*NTot-1,V)
        Iwgt = IGnorm+3*NTot
        Call AMove1(InToWp(NTot),Icen-1,Iwgt+NTot-1,V)
        ICen = Iwgt+NTot
        LSave = ICen + InToWp(NTot) - IXYZ
        If(IRwSrf.gt.-1) then
          If(IRwSrf.gt.0.and.ITqry(IRwSrf).ne.LSave) Call FreeRW(IRwSrf)
          If(IRwSrf.eq.0) IRwSrf = IGetRW(-LSave,0)
          Call FileIO(1,-IRwSrf,LSave,V(IXYZ),0)
          endIf
      else
        If(IPSave(0).eq.0) Write(IOut,9080)
        Call Fileio(2,-IRwSC2,1,XMask,4)
        Call BNUPak(XMask,ISurf,NTot)
        IIOrg = 1
        ICOrg = IIOrg + InToWp(NOrg)
        ICMO = ICOrg + 3*NOrg
        Ixyz = ICMO + NMO*NBas6D
        Ignorm = Ixyz+3*NTot
        Iwgt = IGnorm+3*NTot
        Icen = Iwgt+NTot
        LSave = ICen + InToWp(NTot) - IXYZ
        Call FileIO(2,-IRwSrf,LSave,V(IXYZ),0)
        endIf
C
      Const1 = Two*Pi*(Epsi+one)/(Epsi-One)
      Cons1c = Two*Pi
      Const2 = -One/(One-(One/epsi))
      If(IScal.eq.2) then
        ISc1 = 2
      else
        ISc1 = 1
        endif
      LEFG = 10 * NTot * ISc1
      LEFlux = NTot * 2
      LQtot = NTot
      Lbmf = NTot * ISc1
      LGFact = NTot
      Lq = NTot * ISc1 * ISc1
      LAlpha = NTot * 2 * (ISc1-1)
      LBeta = NTot * (ISc1-1)
      LFFact = NTot
      LPhi = NTot
      Ldwidr = NTot * 2
      LRhat  = NTot
      LISlMp = InToWp(NTot)
      LXYZO  = NTot * 3
C
      Icoord = ICen + IntoWp(NTot)
      Iefg = Icoord + 3 * NTot
      IEFlux = Iefg + LEFG
      IQtot = IEFlux + LEFlux
      Ibmf = Iqtot + LQtot
      IGFact = Ibmf + Lbmf
      Iq = IGFact + LGFact
      IAlpha = Iq + Lq
      IBeta = IAlpha + LAlpha
      IFFact = IBeta + LBeta
      IPhix = IFFact + LFFact
      Idwidr = IPhix + LPhi
      IRhatx = Idwidr + Ldwidr
      IRhaty = IRhatx + LRhat
      IRhatz = IRhaty + LRhat
      ISlMap = IRhatz + LRhat
      ISubst  = ISlMap + LISlMp
      IEnd   = ISubst  + IntoWp(NTot)
C     This code commented for testing pruning for unscaled case.
c     If(IScal.eq.2) then
        IXYZO = IEnd
        IEnd  = IXYZO + LXYZO
c     else
c       IXYZO = IXYZ
c       endif
      If(IEnd.gt.MDV) then
        write(IOut,9090) IEnd-MDV
        Call GauErr('Out of memory in SCISCF -1')
        endIf
C
C     Get the electrostatic properties from OneElI, etc.
C
      MDVEnd = MDV - IEnd + 1
      Call GetFlx(IOut,IPrint,IPFlag,NBas6D,NAtoms,CutOff,NTot,ISc1,
     $  IScal,C,AtmChg,DA,AllowP,V(IXYZ),V(ICoord),V(IEFG),V(IEFlux),
     $  V(IGNorm),V(IPhix),V(IEnd),MDVEnd)
c
c     Before calling GetCGX where we may reduce to the small space
c     we want to compute some quantities for the weight derivatives
c     in the large space
c
      IPbck = IEnd
      Inbck = 1
      If(ISurf.gt.2) Inbck = NTot
      Idpbck = IPbck + Inbck*NOrg
      IDMag  = Idpbck + Inbck*NOrg
      IV2    = IDMag + NTot
      MDVXX = MDV - IV2 + 1
      If(MDVXX.le.0) Call GauErr('Out of memory in SCISCF -2')
      Call Rhat(NTot,NOrg,V(IXYZ),V(ICOrg),V(ICen),V(IRhatx),V(IRhaty),
     $  V(IRhatz),V(IDMag))
      Call dwidr1(NTot,V(Ixyz),V(IGnorm),V(Iwgt),V(Idwidr),V(Icen),NOrg,
     $  V(ICOrg),ISurf,IBWFnc,Inbck,V(Ipbck),V(Idpbck),V(IRhatx),
     $  V(IRhaty),V(IRhatz),V(IDMag))
C
C     Set up memory for GetCGX. Decide now whether A matrix can be stored
C
      ISFact = IEnd
      INeigh = ISFact + NTot * (ISc1 - 1)
      IScr   = INeigh + InToWp(NTot)
      IWgtO  = IScr   + 10 * NTot
      IA     = IWgtO  + NTot
      MDVCGX = MDV - IA + 1
      If(MDVCGX.le.0) Call GauErr('Out of memory in SciFoc -3')
      NTotSq = NTot**2
C     Distance matrix has no diagonal.
      NTotTT = (NTot*(NTot-1))/2
      MDVx = MDV - IntMem(-1,0)
      Lim3 = IA + 2*NTotSq
      Lim2 = IA + NTotSq
      Lim1 = IA + NTotTT
      If(IScal.eq.2.and.Lim3.le.MDVx) then
        IMem = 3
        IV = Lim3
      else if(Lim2.le.MDVx) then
        IMem = 2
        IV = Lim2
      else if(Lim1.le.MDVx) then
        IMem = 1
        IV = Lim1
      else
        IMem = 0
        IV = IA
        endif
      MDV1 = MDV - IV + 1
      If(DoPrnt) then
        If(IScal.eq.1) Write(IOut,9100)
        If(IScal.eq.2) then
          Write(IOut,9110)
          If(PtThr.gt.Zero) Write(IOut,9120)
          Endif
        IMb = 1024*1024
        If(IPrint.ge.1)Write(IOut,9170) ((NTotTT+IMb-1)/IMb),
     $    ((NTotSq+IMb-1)/IMb), ((NTotSq*2+IMb-1)/IMb)
        If(IMem.eq.0) Write(IOut,9130)
        If(IMem.eq.1) Write(IOut,9140)
        If(IMem.eq.2) Write(IOut,9150)
        If(IMem.eq.3) Write(IOut,9160)
        endIf
C
C     GetCGX sets up and solves linear equations for SCIPCM.
C
      Call TStamp(1,'Before GetCGX')
      Call GetCGX(IOut,IPrint,NTot,NTotR,V(IA),V(IEFlux),V(Iq),efcrit,
     $  V(Iffact),V(IWgt),V(IEFG),V(IPhix),V(Ixyz),V(Ignorm),V(ISFact),
     $  V(IGFact),V(IBeta),V(Ibmf),V(IAlpha),V(Iqtot),V(ICen),V(IScr),
     $  V(ISubSt),V(ISlMap),V(INeigh),DoPrun,PtThr,PtInc,UseOld,MaxSFa,
     $  IRwSC2,LScSl2,SymCut,V(IV),MDV1,Const1,cons1c,IScal,IMem,ISc1,
     $  V(IWgtO),V(IXYZO))
      Call TStamp(1,'After GetCGX')
C
      ESol = Pt5*(SProd(NTotR,V(Iqtot),V(IPhix)))
      GaussN = SProd(NTotR,V(IEFlux+NTot),V(IWgt)) / (Four*Pi)
      Gauss = SProd(NTotR,V(IEFlux),V(IWgt)) / (Four*Pi)
      PolTot = Const2*ArrSum(NTotR,V(Iqtot))
      If(IScal.eq.2) then
        Gauss = Gauss + GaussN
        PolTtU = Const2*(ArrSum(NTotR,V(Iq))+ArrSum(NTotR,V(Iq+NTot)))
      else
        PolTtU = PolTot
        endIf
      If(Abs(GaussN - ArrSum(NAtoms,AtmChg)).gt.One
     $  .or.Abs(PolTtU - Gauss).gt.Pt01) then
        Write(IOut,9180)GaussN,ArrSum(NAtoms,AtmChg),
     $    Abs(PolTtU - Gauss)
        Call GauErr('Surface Problems in SciFoc')
        Endif
C
C
C     Finally, compute Fock matrix contribution. Put it in DB.
C
      IScr1 = IEnd
      IEll  = IScr1  + NTotR
      IEm   = IEll   + NTot
      IV    = IEm    + 3*NTot
      MDVFoc = MDV - IV + 1
      If(MDVFoc.le.0) Call GauErr('Out of memory in SciSCF -4')
      Call TStamp(1,'Before SCIFoc')
      Call SCIFoc(IOut,IPrint,IScal,NTot,Natoms,NBas6D,C,V(Ixyz),
     $  V(ICoord),V(Ignorm),V(Iwgt),V(Iffact),V(IEfg),DA,ISc1,V(Iqtot),
     $  V(IAlpha),V(Ibmf),V(IBeta),V(Iq),V(IGFact),V(IEll),V(IEm),
     $  V(IScr1),DB,Cutoff,MoBas,IOpCl,NMO,V(ICMO),IPFlag,AllowP,MaxTyp,
     $  NTotR,V(Idwidr),V(IXYZO),V(IRhatx),V(IRhaty),V(IRhatz),
     $  V(ISlMap),V(ISubst),Const1,Cons1c,V(IV),MDVFoc)
      Call TStamp(1,'After SCIFoc')
C
C     Write out terms for gradients and printing (SciPrt).
C
      If(DoGrad) then
        KNTotR = 9
        KNTot  = 7
        If(IScal.eq.2) KNTotR = KNTotR + 7
        ItsLen = 9 + KNTotR*NTotR + KNTot*NTot
        Call ConDDF(IRwSCI,ItsLen)
        Call BNPack(NTotR,NTot,XMask)
        Call FileIO(1,-IRwSCI,1,XMask,0)
        Call BNPack(IScal,ISurf,XMask)
        Call FileIO(1,IRwSCI,1,XMask,0)
        Call Fileio(1,IRwSCI,1,ESol,0)
        Call Fileio(1,IRwSCI,1,GaussN,0)
        Call Fileio(1,IRwSCI,1,Gauss,0)
        Call Fileio(1,IRwSCI,1,PolTot,0)
        Call Fileio(1,IRwSCI,1,PolTtU,0)
        Call Fileio(1,IRwSCI,1,SAreaT,0)
        Call Fileio(1,IRwSCI,1,Volume,0)
        Call Fileio(1,IRwSCI,NTotR,V(Ixyz),0)
        Call Fileio(1,IRwSCI,NTotR,V(Ixyz+NTot),0)
        Call Fileio(1,IRwSCI,NTotR,V(Ixyz+2*NTot),0)
        Call Fileio(1,IRwSCI,NTotR,V(Ignorm),0)
        Call Fileio(1,IRwSCI,NTotR,V(Ignorm+NTot),0)
        Call Fileio(1,IRwSCI,NTotR,V(Ignorm+2*NTot),0)
        Call Fileio(1,IRwSCI,NTotR,V(Iqtot),0)
        Call Fileio(1,IRwSCI,NTotR,V(Ibmf),0)
        Call Fileio(1,IRwSCI,NTot,V(IEll),0)
        Call Fileio(1,IRwSCI,NTot*3,V(IEm),0)
        Call Fileio(1,IRwSCI,NTot*3,V(IXYZO),0)
        Call Fileio(1,IRwSCI,NTotR,V(IWgt),0)
        If(IScal.eq.2) then
          Call FileIO(1,IRwSCI,NTotR,V(IBeta),0)
          Call FileIO(1,IRwSCI,NTotR,V(IAlpha),0)
          Call FileIO(1,IRwSCI,NTotR,V(IAlpha+NTot),0)
          Call FileIO(1,IRwSCI,NTotR,V(IGFact),0)
          Call FileIO(1,IRwSCI,NTotR,V(IQ+3*NTot),0)
          Call FileIO(1,IRwSCI,NTotR,V(IQ+2*NTot),0)
          Call FileIO(1,IRwSCI,NTotR,V(Ibmf+NTot),0)
          endif
        endIf
C
C     Write out info for subsequent iterations
C
      NDPRun = .not.DoPrun
      Call BNPack(IDPRun,Iscal,XMask)
      Call Fileio(1,-IRwSC2,1,XMask,0)
      Call Fileio(1,IRwSC2,1,epsi,0)
      Call Fileio(1,IRwSC2,1,Cont,0)
      Call BNPack(IPhi,ITheta,XMask)
      Call Fileio(1,IRwSC2,1,XMask,0)
      Call BNPack(ISurf,NTot,XMask)
      Call Fileio(1,IRwSC2,1,XMask,0)
      Call Fileio(1,IRwSC2,1,PtThr,0)
      Junk = 0
      Call BNPack(NOrg,Junk,XMask)
      Call FileIO(1,IRwSC2,1,XMask,0)
      Call FileIO(1,IRwSC2,InToWp(NOrg),V(IIOrg),0)
C
C     Purify and add solvation Fock matrix to total Fock matrix.
C
      Call F2Pure(IOut,IPrint,IPureD,IPureF,1,1,0,NBasis,NBas6D,.True.,
     $  .True.,.False.,.False.,.False.,.False.,DB,XX,XX,XX,V,MDV)
      Call ACASB(NTT,FA,DB,FA,Pt5)
      Call ACASB(NTT,FB,DB,FB,Pt5)
      If(DoPrnt) Call SciPrt(IOut,0,0,NAtoms,C,XX,Junk,Junk,Junk,
     $  Junk,Junk,V,MDV)
      Return
      End
*Deck SCISrf
      Subroutine SCISrf(IOut,IPrint,DoPrnt,NAtoms,NBasis,NTT6D,Cont,C,
     $  AtmChg,ITheta,IPhi,NAng,MPts,NTot,RadT,NumR,CutOff,CutPt,IBWFnc,
     $  MoBas,NMO,IOpCl,ISurf,NOrg,IOrg,COrg,XYZk,GNorm,Wgt,P,ICen,Rad,
     $  Unit,XYZ,UnitT,XYZs,Mx,Rho,AngWgt,Keep,GRho,PBeck,cmo,SArea,
     $  SAreaT,Volume,V,MDV)
      Implicit Real*8 (A-H,O-Z)
C
      Logical MOBas, DoPrnt, OneOrg, Keep(NAng)
      Parameter (MaxRIt=1000)
      Dimension Unit(3,NAng),XYZ(NAng,3),UnitT(NAng,3),XYZS(NAng,3),
     $  Mx(NAng),Rho(NAng),AngWgt(NAng),GNorm(MPts,3),GRho(NAng,3),
     $  C(3,NAtoms),AtmChg(NAtoms),V(MDV),P(NTT6D),CMO(NBasis,NMO),
     $  Rad(NOrg),Wgt(MPts),XYZK(MPts,3),
     $  ICen(MPts),PBeck(NAng,NOrg),XX(1),CurCnt(3),IOrg(NOrg),
     $  COrg(3,NOrg),SArea(NOrg)
      Save Tol,Two,Zero,FStp,XX
      Data Tol/1.0d-13/,Two/2.0d0/,Zero/0.0d0/,XX/0.0d0/,
     $  Fstp/4.0d0/
 1010 Format(' All intersections with cavity surface found')
 1020 Format(' Intersections with cavity surface found from nucleus',I5)
 1030 Format(' Solvent accessible surface area for nucleus ',
     $  I4, ' = ',1PE12.6,' au',/,
     $  ' rMax between nucleus ',I4,' and its part of ',
     $  'cavity surface = ',1PE12.6,' au',/,
     $  ' rMin between nucleus ',I4,' and its part of ',
     $  'cavity surface = ',1PE12.6,' au',/,
     $  ' Number of cavity surface points = ',I6,
     $  ' from nucleus ',I4)
 1040 Format(' rMax between center of nuclear charge and ',
     $  'cavity surface = ',1PE12.6,' au',/,
     $  ' rMin between center of nuclear charge and ',
     $  'cavity surface = ',1PE12.6,' au')
 1050 Format(' Total "Solvent Accessible Surface Area" of solute ',
     $  'cavity = ',1PE12.6,' au',/,
     $  ' Volume of solute cavity = ', 1PE12.6,' au',/,
     $  ' Total number of cavity surface points = ',I6)
C
      If(ISurf.gt.2) then
        Do 5 I = 1, NOrg
          COrg(1,I) = C(1,IOrg(I))
          COrg(2,I) = C(2,IOrg(I))
    5     COrg(3,I) = C(3,IOrg(I))
      else
        Call CntMas(NAtoms,AtmChg,C,SumChg,COrg(1,1))
        endIf
      Do 10 IA = 1, NOrg
   10   Rad(IA) = RadT
      OneOrg = ISurf.eq.2
      NTot = 0
      SAreaT = Zero
      Volume = Zero
      If(OneOrg) then
        NAtms = 1
      else
        NAtms = NOrg
        endIf
      If(ITheta.eq.1.or.IPhi.eq.1) then
        IThetz = 0
        IPhz = 0
        NAngZ = ITheta*IPhi
      else
        IThetz = ITheta
        IPhz = IPhi
        NAngZ = IThetz*IPhz
        endIf
      Do 200 IA = 1, NAtms
        Call AngQad(Unit,AngWgt,NAngZ,IThetz,IPhz,Tol)
        If(OneOrg) then
          Call CntMas(NAtoms,AtmChg,C,SumChg,CurCnt)
        else
          Call AMove(3,COrg(1,IA),CurCnt)
          endIf
        XYZ(1,1)=CurCnt(1)
        XYZ(1,2)=CurCnt(2)
        XYZ(1,3)=CurCnt(3)
        Call DoRhGr(IOut,IPrint,IOpCl,1,MOBas,1,NAng,NBasis,NMO,
     $    CutOff,XYZ,P,CMO,Rho,1,XX,1,XX,V,MDV)
        If(Rho(1).lt.Cont)Call GauErr(
     $    ' A solvent integration origin lies outside the cavity!')
        Do 20 IX = 1, 3
          Do 20 I = 1, NAngZ
   20       XYZS(I,IX) = CurCnt(IX) + (Rad(IA)+FStp)*Unit(IX,I)
        Do 100 ICut = 1, NumR
          Stp = FStp / Two**(ICut-1)
          Call IClear(NAngZ,Mx)
          Do 30 IX = 1, 3
            Do 30 I = 1, NAngZ
              UnitT(I,IX) = Unit(IX,I)
   30         XYZ(I,IX) = XYZS(I,IX) - Stp*Unit(IX,I)
          Left = NAngZ
          MaxItR = 1
          If(ICut.eq.1) MaxItR = MaxRIt
          MaxItR = MaxRIt
          Do 90 ItR = 1, MaxItR
            Call DoRhGr(IOut,IPrint,IOpCl,1,MOBas,Left,NAng,NBasis,NMO,
     $        CutOff,XYZ,P,CMO,Rho,1,XX,1,XX,V,MDV)
            L = 0
            Do 50 J = 1, Left
              If(Rho(J).gt.Cont) then
                L = L + 1
                XYZ(L,1) = XYZ(J,1) + Stp*UnitT(J,1)
                XYZ(L,2) = XYZ(J,2) + Stp*UnitT(J,2)
                XYZ(L,3) = XYZ(J,3) + Stp*UnitT(J,3)
                UnitT(L,1) = UnitT(J,1)
                UnitT(L,2) = UnitT(J,2)
                UnitT(L,3) = UnitT(J,3)
                Mx(L) = Mx(J) + J - L
              else
                XYZS(J+Mx(J),1) = XYZ(J,1)
                XYZS(J+Mx(J),2) = XYZ(J,2)
                XYZS(J+Mx(J),3) = XYZ(J,3)
                endIf
   50         Continue
            If(L.le.0) goto 100
   90       Left = L
  100     Continue
        If(OneOrg) then
          If(DoPrnt) Write(IOut,1010)
        else
          If(DoPrnt) Write(IOut,1020) IOrg(IA)
          endIf
        Call LSet(NAngZ,Keep)
        If(.not.OneOrg) Call SurfBW(NOrg,NAngZ,NAng,IA,CutPt,IBWFnc,
     $    COrg,XYZS,PBeck,AngWgt,Keep)
        Call CmpXYZ(NAngZ,NAng,Keep,CurCnt,AngWgt,NAct,DMx,DMn,XYZS,
     $    Wgt(NTot+1))
        If(NAct.gt.0) then
          Call DoRhGr(IOut,IPrint,IOpCl,2,MOBas,NAct,NAng,NBasis,NMO,
     $      CutOff,XYZS,P,CMO,Rho,NAng,GRho,1,XX,V,MDV)
          Call ScSrSt(NTot,NAct,IA,NAng,MPts,CurCnt,XYZS,GRho,
     $      ICen(NTot+1),XYZK,GNorm,Wgt(NTot+1),Volume)
          endIf
        SArea(IA) = ArrSum(NAct,Wgt(NTot+1))
        NTot = NTot + NAct
        SAreaT = SAreaT + SArea(IA)
        If(DoPrnt) then
          If(OneOrg) then
            Write(IOut,1040) DMx, DMn
          else
            Write(IOut,1030) IOrg(IA),SArea(IA),IOrg(IA),DMx,IOrg(IA),
     $        Dmn,NAct,IOrg(IA)
            endIf
          endIf
  200   Continue
      Volume = Volume / GFloat(3)
      If(DoPrnt) Write(IOut,1050) SAreaT,Volume,NTot
      If((DoPrnt.or.IPrint.ge.1).and.NOrg.gt.1) then
        Write(IOut,*)'Percentage of total area from each origin:'
 1060 Format(4(:,1x,'Org: ',I5,': ',f5.1,'% '))
        Write(IOut,1060)(IAA,(GFLoat(100)*SArea(IAA)/SAreaT),
     $    IAA = 1, NOrg)
        endif
      Return
      End
*Deck SCnPDM
      Subroutine SCnPDM(IOut,IPrint,PrntDM,cc,maxcor,X,istart,jip,ip,
     $  ip2,NZ,iend,N,niter,NOcc,IOpMe,Tcg2,npf,xmax,xmin,demp,pfthr,
     $  IMemOp,zthr)
      Implicit Real*8(A-H,O-Z)
C
C     Top level routine performing Canonical Purification of the Density
C     Matrix (PDM)
C
      Real*8 min, cc(maxcor), X(N)
      Logical PrntDM,Update,Gershg
      Save DempMx
      Data DempMx/5.d-3/
 1000 Format(/' Canonical Purification')
 1010 Format(' CnPDmF Real*8 Arrays:  15*n + 6.5*nz ',/,
     $       ' Memory needed for PDM',i14,' Max memory ',i14)
 1200 Format(' Maximum Number of Purification Cycles',i8)
 1500 Format(' Purification Threshold',1pd12.2)
 1300 Format(' Using Gershgorin"s Formulas to Estimate Eigenvalues')
 1400 Format(' Using Lanczos to Estimate Eigenvalues')
 1700 Format(' Maximum Eigenvalue:',f9.4,'Minumum:',f9.4)
 2300 Format(' Using Original Sparse Matrix Multiplication')
 2400 Format(' Using Bin Sorting Sparse Matrix Multiplication')
 2500 Format(' Grand Canonical Purification Idempotency:',1PD12.2)
 2600 Format(' Idempotency:',1pd12.2,' Criterion:',1pd12.2,/,
     $  ' Halting: Density matrix is no longer idempotent.'/,
     $  ' Please try running with tighter neglect thresholds')
c
c  Initialize for conjugate gradient
c
      NP = NProc(0)
      NPMax = Mod(IOpMe,1000)
      If(NPMax.eq.0) NPMax = 100
      Gershg = ((Mod(IOpMe,10000)/1000).eq.1)
      If(IPrint.ge.1) Write(IOut,1000)
      If(PrntDM) then
        If(Gershg) then
          Write(IOut,1300)
        else
          Write(IOut,1400)
          endIf
        Write(IOut,1200) npmax
        Write(IOut,1500) pfthr
        endIf
      numbig=5
      numsml=3
c
c  Allocate memory for all matrices
c
      ip3    = istart
      iscr   = ip3    + NZ + N
      jix    = iscr   + NZ + N
      ixlanc = jix + INTOWP(N)
      jimap  = jix
      jmap   = jimap + INTOWP(N+1+NZ)
      ixx    = jmap  + INTOWP(NZ)
      jixx  = ixx    + 2*N*NP
      iend = jixx   + InToWP(N*NP*3)
      Len = iend - ip2 + 1
      numcpo = Max(numbig,numsml)
      last = ixlanc + N*(5+4*numcpo) + numcpo*(9*numcpo + 24)
      last = Max(last,iend)
      if(iprint.ge.2) write(Iout,1010) last, maxcor
      call tstcor(last,maxcor,'ChbDim 2')
      call amove(N+NZ,cc(ip2),cc(ip))
      call aclear(Len,cc(ip2))
      UpDate = .true.
      istart = last
C
C     Find the lowest and highest eigenvalues of F
C
      If(Gershg) then
        Call Amove(N,cc(iP),cc(iP2))
        Call Amove(N,cc(iP),cc(iP3))
        Call Gersh(N,cc(jip),cc(iP+N),cc(iP2),cc(iP3))
        Call SpTrp(cc(jIP),cc(iP),cc(jIMap),cc(iP2+N),N)
        Call Gersh(N,cc(jIMap),cc(iP2+N),cc(iP2),cc(iP3))
        xmax = ArrMax(cc(iP3),N,.false.)
        xmin = ArrMin(cc(iP2),N,.false.)
      else
        If(Update) Call DavGts(IOut,IPrint-4,cc,MaxCor,ixlanc,jIP,iP,
     $    N,xmax,xmin)
        endIf
      If(IPrint.ge.2) write(IOut,1700) xmax, xmin
C     Shift the Hamiltonian so that its centered at zero
      gmu = ArrSum(N,cc(iP))/GFloat(N)
      ambda = min((GFloat(NOcc)/(xmax-gmu)),
     $  ((GFloat(N)-GFloat(NOcc))/(gmu-xmin)))
      Call ASet(N,gmu,X)
      Call ANeg((N+NZ),cc(iP),cc(iP))
      Call AAdd(N,cc(iP),X,cc(iP))
      Call AScale((N+NZ),(ambda/GFloat(N)),cc(iP),cc(iP))
      Call ASet(N,GFloat(NOcc)/GFloat(N),X)
      Call AAdd(N,cc(iP),X,cc(iP))
      Method = Mod(IMemOp,10)
      iend = istart
      If(Method.ne.1) Call FSIni1(cc,maxcor,istart,iBnFac,iSizCV,iend,
     $  cc(jIP),cc(iP),cc(jIP),cc(iP),N,NBin,IMemOp,method,ZThr,.true.,
     $  .true.)
      If(Method.eq.1) then
        If(IPrint.ge.2) Write(IOut,2300)
        LNNP  = 0
        LN1NZ = 0
        LNNZ  = 0
        LNBin = 0
        jiscr1 = jip
      else if(Method.eq.2) then
        If(IPrint.ge.2) Write(IOut,2400)
        LNNP  = N*NP
        LN1NZ = N+1+NZ
        LNNZ  = N+NZ
        LNBin = N*(NBin+1)
        jiscr1 = iend
        iend   = jiscr1 + InToWP(N+NZ+1)
      else
        Call GauErr('Unknown Method in SCnPDM')
        endIf
      jIRowA = iend
      iXRowA = jIRowA + InToWP(N*NP)
      jJNEW  = iXRowA + N*NP
      iANEW  = jJNEW + InToWP(N*NP)
      jiX    = iANEW + N*NP
      iRowNo = jiX + InToWP(LNNP)
      iscr1  = iRowNo + InToWP(N)
      jibt   = iscr1 + LNNP
      ibt    = jibt + intowp(LN1NZ)
      iRABin = ibt  + LNNZ
      iBinB  = iRABin + intowp(LNNP)
      iBinBt = iBinB + intowp(LNBin)
      iRowSr = iBinBt + intowp(LNBin)
      iscr2  = iRowSr + LNNP
      iend   = iscr2 + intowp(LNNP)
      Call CnPurf(IOut,IPrint+1,cc(jIP),cc(jiscr1),cc(jibt),cc(iP),
     $  cc(iP2),cc(iP3),cc(iScr),cc(ibt),cc(jIMap),cc(JMap),cc(iBnFac),
     $  cc(iSizCV),cc(iRowNo),cc(jIRowA),cc(jJNEW),cc(jiX),cc(iRABin),
     $  cc(iscr2),cc(iXRowA),cc(iANEW),cc(iscr1),cc(iRowSr),cc(iBinB),
     $  cc(iBinBt),NBin,NOcc,N,NZ,npfcy,NP,NPMax,demp,pfthr,zthr,Method)
      If(NIter.gt.4.and.(Abs(Demp).gt.DempMx)) then
        Write(IOut,2600) Demp, DempMx
        Call GauErr('Stopping in ScnPDM')
        EndIf
      If(IPrint.eq.0) Write(IOut,2500) Demp
      Call CheckF(cc(jip),cc(ip),X,N,0,NZ,'P After frmchb',5,IPrint)
      Return
      End
*Deck ScSrSt
      Subroutine ScSrSt(N,NAct,IA,NDimS,NDimK,CurCnt,XYZS,GRho,ICen,
     $  XYZK,GNorm,Wgt,CurVol)
      Implicit Real*8(A-H,O-Z)
C
C     Store the coordinates during generation of SCI surface.
C
      Dimension CurCnt(3), XYZS(NDimS,3), GRho(NDimS,3), XYZK(NDimK,3),
     $  ICen(*), GNorm(NDimK,3), Wgt(*)
C
      Do 10 K = 1, NAct
        KN=K+N
        ICen(K) = IA
        XYZK(KN,1) = XYZS(K,1)
        XYZK(KN,2) = XYZS(K,2)
        XYZK(KN,3) = XYZS(K,3)
        GRhoMg = Sqrt(GRho(K,1)**2+GRho(K,2)**2+GRho(K,3)**2)
        GNorm(KN,1) = -GRho(K,1)/GRhoMg
        GNorm(KN,2) = -GRho(K,2)/GRhoMg
        GNorm(KN,3) = -GRho(K,3)/GRhoMg
        Dist = Sqrt((XYZS(K,1)-CurCnt(1))**2+(XYZS(K,2)-CurCnt(2))**2+
     $    (XYZS(K,3)-CurCnt(3))**2)
        Wgt(K) = Wgt(K) * Dist / ((GNorm(KN,1)*(XYZS(K,1)-CurCnt(1)) +
     $                             GNorm(KN,2)*(XYZS(K,2)-CurCnt(2)) +
     $                             GNorm(KN,3)*(XYZS(K,3)-CurCnt(3))))
   10   CurVol = CurVol + Wgt(K)*(XYZS(K,1)*GNorm(KN,1)+XYZS(K,2)*
     $    GNorm(KN,2)+XYZS(K,3)*GNorm(KN,3))
      Return
      End
*Deck SelDFi
      Subroutine SelDFi(IOut,DoPrnt,RHF,ISCF,NSCF,DoSCF,IPOA1,IPOB1,
     $  IRwPOA,IRwPOB,SEEn,GesEn,IPOAUs,IPOBUs,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Figure out which sparse density files to use in each pass through CycOpn.
C
      Parameter (ISEGes=1,ISEEn=2,IGesEn=3,IVac=4)
      Logical DoPrnt, RHF, DoSCF(NSCF), DoComp, DoReal
      Dimension V(MDV)
      Save One
      Data One/1.0d0/
 1000 Format(' Using semiemprical P for guess.')
 1010 Format(' Using P from previous geometry for guess.')
C
      DoComp = DoSCF(ISEEn).and.DoSCF(IGesEn)
      DoReal = ISCF.eq.IVac.or.(ISCF.eq.NSCF.and..not.DoSCF(IVac))
      IPOAUs = IRwPOA
      IPOBUs = IRwPOB
      If(DoComp) then
        If(ISCF.eq.ISEGes.or.ISCF.eq.ISEEn) then
          IPOAUs = IPOA1
          IPOBUs = IPOB1
        else if(DoReal) then
          If(SEEn.lt.GesEn) then
            If(DoPrnt) Write(IOut,1000)
            Call SpFCpy(RHF,NBasis,IPOA1,IPOB1,IRwPOA,IRwPOB,One,One,V,
     $        MDV)
            Call FileIO(5,IPOA1,0,V,0)
            If(.not.RHF) Call FileIO(5,IPOB1,0,V,0)
          else
            If(DoPrnt) Write(IOut,1010)
            endIf
          endIf
        endIf
      Return
      End
*Deck SForm3
      Subroutine sform3(CC,Ncore,Norbs,IEnd,Iout,Thresh,jipa,ipa,NZpa,
     $  jiz,iz,NZz,jiz2,iz2,NZz2,IX,ImemOp)
      Implicit Real*8(A-H,O-Z)
      Real*8 CC(Ncore)
      Logical DoLttA,DoLttb,DoLttc,Screen
      Save One,Two
      Data One/1.d0/,Two/2.d0/
C
C     Generate old sign matrix : zo = 2p - i
C
      Istart = IEnd
      jizo = jipa
      izo = ipa
      NZzo = NZpa
      Call SpAMvM(CC,Ncore,Istart,IEnd,jizo,izo,Norbs,2)
      Call AScale(Norbs+NZzo,Two,CC(izo),CC(izo))
      do 10 i = 0, Norbs - 1
   10   CC(izo+i) = CC(izo+i) - One
C
C     Generate Format from the Union of z and (1 - z^2)^2 !
C     and change the Format of z and z2 !
C     z2 contains here (1-z2)
C
C     first: lower triangle of (1-z2)*p
C
      DoLtta = .True.
      DoLttb = .True.
      DoLttc = .True.
      Screen = .False.
      Imul = 1
      IMode = 2
      jidz1 = IEnd
      idz1 = IEnd
      Istart = IEnd
      Call FsMpyM(CC,Ncore,CC(jiz2),CC(iz2),CC(jizo),CC(izo),Istart,
     $  Istart,jidz1,idz1,IEnd,Thresh,Norbs,NZdz1,CC(jidz1),CC(idz1),
     $  DoLttA,DoLttb,DoLttc,Screen,Imul,IMode,ImemOp)
C
C     Second : p*(1-z2)
C
      jidz2 = IEnd
      idz2 = IEnd
      Istart = IEnd
      Call FsMpyM(CC,Ncore,CC(jizo),CC(izo),CC(jiz2),CC(iz2),Istart,
     $  Istart,jidz2,idz2,IEnd,Thresh,Norbs,NZdz2,CC(jidz2),CC(idz2),
     $  DoLttA,DoLttb,DoLttc,Screen,Imul,IMode,ImemOp)
C
C     Symbolic addition, generate new form
C
      jisf1 = IEnd
      Call SpAddS(CC(jidz1),CC(jidz2),CC(jisf1),CC(IX),Norbs,NZsf1)
      IEnd = jisf1 + IntoWp(Norbs+1+NZsf1)
      jisf1t = IEnd
      IEnd = jisf1t + IntoWp(Norbs+1+NZsf1)
      Call SpTrpS(CC(jisf1),CC(jisf1t),Norbs)
      Call SpTrpS(CC(jisf1t),CC(jisf1),Norbs)
      jisf2 = IEnd
      Call SpAddS(CC(jisf1),CC(jipa),CC(jisf2),CC(IX),Norbs,NZsf2)
      IEnd = jisf2 + IntoWp(Norbs+1+NZsf2)
      jisf2t = IEnd
      IEnd = jisf2t + IntoWp(Norbs+1+NZsf2)
      Call SpTrpS(CC(jisf2),CC(jisf2t),Norbs)
      Call SpTrpS(CC(jisf2t),CC(jisf2),Norbs)
      jisf = IEnd
      Call SpAddS(CC(jiz),CC(jisf2),CC(jisf),CC(IX),Norbs,NZsf)
      IEnd = jisf + IntoWp(Norbs+1+NZsf)
      jisft = IEnd
      IEnd = jisft + IntoWp(Norbs+1+NZsf)
      Call SpTrpS(CC(jisf),CC(jisft),Norbs)
      Call SpTrpS(CC(jisft),CC(jisf),Norbs)
C
C     Transform z into the new form
C
      NZzn = NZsf
      jizn = IEnd
      IEnd = jizn + IntoWp(Norbs+1+NZzn)
      Call IMove1(Norbs+1+NZzn,0,(jizn-jisf)*IntPWP(0),CC(jisf))
      izn = IEnd
      IEnd = izn + Norbs + NZzn
      Call AMove1(Norbs,(iz-1),(izn-1),CC)
      Call NwForm(CC(jiz),CC(iz+Norbs),CC(jizn),CC(izn+Norbs),CC(IX),
     $  Norbs)
C
C     Transform z2 into the new form
C
      NZz2n = NZsf
      jiz2n = IEnd
      IEnd = jiz2n + IntoWp(Norbs+1+NZz2n)
      Call IMove1(Norbs+1+NZz2n,0,(jiz2n-jisf)*IntPWP(0),CC(jisf))
      iz2n = IEnd
      IEnd = iz2n + Norbs + NZz2n
      Call AMove1(Norbs,(iz2-1),(iz2n-1),CC)
      Call NwForm(CC(jiz2),CC(iz2+Norbs),CC(jiz2n),CC(iz2n+Norbs),
     $  CC(IX),Norbs)
C
C     Move back z, IX and z2
C
      Istart = jiz
      jiz = jizn
      iz = izn
      NZz = NZzn
      Call SpAMvM(CC,Ncore,Istart,IEnd,jiz,iz,Norbs,2)
      IX = IEnd
      IEnd = IX + Norbs
      Istart = IEnd
      jiz2 = jiz2n
      iz2 = iz2n
      NZz2 = NZz2n
      Call SpAMvM(CC,Ncore,Istart,IEnd,jiz2,iz2,Norbs,2)
      Return
      End
*Deck SignMt
      Subroutine SignMt(IOut,IPrint,PrntDM,CC,Ncore,X,N,NElect,jipa,ipa,
     $  NZpa,jif,if,NZf,istrt1,IEnd,Niter,IopDMe,Thresh,xIdemp,ImemOp,
     $  OuterF)
      Implicit Real*8 (a-h,o-z)
      Dimension CC(Ncore), X(N)
      Logical FixTyp,Fixed,CNSI,OuterF,PrntDM
      Save Ten, Trace
      Data Ten/10.d0/, Trace/0.0d0/
  301 Format(/,' Fixed form sparse matrix operations.')
  300 Format(' Default general form used.')
  310 Format(' Form generated in SignMt: Union of P,F,PF+tp')
  320 Format(' Form generated in SignMt:',
     $  ' Union of P,F,(I-F**2)(2P-I)+tp')
  330 Format(/,' Let-It-Grow sparse matrix operations.')
 2222 Format(/,' Trace of density matrix from iterations:',f21.10,
     $ /,'number of electrons = ',i10)
 1000 Format(' Canonical NSI Idempotency:',1PD12.2)
 2000 Format(' Grand Canonical NSI Idempotency:',1PD12.2)
c
c CNSI = .True. : C-NSI iteration is used
c       .False. : GC-NSI iteration
c
c Fixed = .True. : Fixed Form multiplications and additions
c       .False. : Let-It-Grow multiplications and additions
c
c FixTyp = .True. :  U{P,F,PF,FP} type Format is generated here
c        .False. : U{P,F,(I-F**2)*(2P-I)+tp} type Format is used
c
c OuterF = .True. :  Format for Fixed form operations
c                  comes from an outer routine
c        .False. : Format is generated inside the sign-routines
c
c ItMax = IopDMe/1000 or 100 : maximum number of NSI cycles
c
c
      ItMax = mod(IopDMe,10000)
C     default for max number of iterations for signa
      If(ItMax.eq.0) ItMax = 100
      CNSI = Mod((IopDMe/1000),10).lt.1
      IForm = mod((IopDMe/10000),10)
      Fixed = IForm.ne.1
      FixTyp = IForm.eq.3
      If(PrntDM) then
        Write(IOut,*) ' OuterF = ',OuterF,' CNSI = ',CNSI,' FixTyp = ',
     $    FixTyp,' Fixed = ',Fixed
        If(Fixed) then
          Write(IOut,301)
          If(OuterF) then
            Write(IOut,300)
          else if(FixTyp) then
            Write(IOut,310)
          else
            Write(IOut,320)
            endIf
        else
          Write(IOut,330)
          endIf
        endIf
C
C     conv.crit. for Newton-Schultz Iteration
C
      Istart = istrt1
      Ibegin = jipa
      ConvCr = Thresh*Ten
C
C     Generate Format for Fixed Format multiplies
C
      If(Fixed.and.FixTyp.and..not.OuterF) then
        jifo = jif
        IX = IEnd
        IEnd = IX + N
        Istart = IEnd
        Call FPenh(CC,Ncore,CC(jif),CC(if),CC(jipa),CC(ipa),CC(IX),
     $    Thresh,Istart,jipf,ipp,iff,IEnd,N,NZ,1,.True.,.True.,0,ImemOp,
     $    NElect)
        If(IPrint.eq.1) Call SpRprt(NZ,N*(N-1)/2,Thresh,'all matrices')
        NZf = NZ
        jif = jipf
        if  = iff
        istrt1 = iend
        Call SpAMvM(CC,NCore,istrt1,iend,jif,if,N,2)
        Istart = IEnd
        jipa = jipf
        ipa = ipp
        NZpa = NZ
        Call SpAMvM(CC,Ncore,Istart,IEnd,jipa,ipa,N,2)
        Istart = jifo
        Call SpAMvM(CC,Ncore,Istart,IEnd,jif,if,N,2)
        Istart = IEnd
        Call SpAMvM(CC,Ncore,Istart,IEnd,jipa,ipa,N,2)
        endIf
      If(Fixed.and.OuterF) then
        NZ = NZf
        ifa = if
        jif = ifa
        if = jif + IntoWp(N+1+NZ)
        IEnd = if + N + NZ
        Call TstCor(IEnd,Ncore,'SignMt 1')
        Call AMove1(N+NZ,(ifa-1),(if-1),CC)
        Call IMove1(N+1+NZ,0,(jif-jipa)*IntPWP(0),CC(jipa))
        endIf
C
C     Calculate new p matrix
C
      Istart = IEnd
      Call DetFIt(ItMax,IOut,IPrint,CC,Ncore,N,NElect,jif,if,NZf,NZpa,
     $  jipa,ipa,IEnd,Niter,Trace,Thresh,ConvCr,xIdemp,ImemOp,OuterF,
     $  Fixed,FixTyp,CNSI)
      If(IPrint.eq.0) then
        If(CNSI) Then
          Write(IOut,1000) xIDemp
        else
          Write(IOut,2000) xIDemp
          endIf
        endIf
      If(IPrint.ge.1) Write(IOut,2222) Trace,NElect
C
C     Scale Trace to the exact number of electrons
C
      If(Abs(Trace-gfloat(NElect)).gt.Thresh) then
        Fact = GFloat(NElect)/Trace
        Call AScale(N+NZpa,fact,CC(ipa),CC(ipa))
        endIf
      Istart = Ibegin
      Call SpAMvM(CC,Ncore,Istart,IEnd,jipa,ipa,N,2)
      Return
      End
*Deck SimFP
      Subroutine SimFP(IOut,DoPrnt,Ne,NBasis,IRwV,IRwSim,IRwP,IRwPO,NP,
     $  MaxGDI,NMax,DoDMS,IBas,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     This subroutine takes a density step for the simultaneous
C     optimization
C
      Parameter (NUseM=1000)
      Dimension V(*)
      Logical DoDMS, DoPrnt
      Save One, Pt1
      Data One/1.0d0/, Pt1/0.1d0/
 1000 Format(' Forming Sim-DIIS Fock matrix using ',I4,
     $  ' saved points, Tot=',F8.5)
 1010 Format(' DIIS coeff''s:',5F13.5)
C
C     Do some memory allocation
C
      NTT  = NBasis*(NBasis+1)/2
      NBSq = NBasis*NBasis
      Call TStamp(1,'Top of SimFP')
      IC     = 1
      IFNew  = IC   + NP
      IScr1  = IFNew + NBSq
      IScr2  = IScr1 + NBSq
      IScr3  = IScr2 + NBSq
      IScr4  = IScr3 + NBasis
      IEnd   = IScr4 + NBasis
      Call TstCor(IEnd,MDV,'SimFP1')
C
C     Extrapolate F' = Sum Ci*Fi
C
      Call RwFGED(1,IRwSim,NMax,MaxGDI,IEn,INNEn,ICoeff,IQNDIS,IEnDIS,
     $  IMask,ICnvQM)
      Call FileIO(2,-IRwSim,NP,V(IC),ICoeff)
      NUse = Min(NP,NUseM)
      Tot = ArrSum(NUse,V(IC))
      Write(IOut,1000) NUse, Tot
      If(Tot.lt.Pt1) then
        NUse = 1
        V(IC) = One
      else
        Call AScale(NUse,One/Tot,V(IC),V(IC))
        endIf
      Write(IOut,1010) (V(IC+I-1),I=1,NUse)
      Call TstCor(IScr1+NUse*NTT,MDV,'SimFP2')
      Call FileIO(2,-IRwPO,NTT*NUse,V(IScr1),0)
      Call MatMpy(NTT,NUse,1,V(IScr1),V(IC),V(IFNew))
      If(DoDMS) then
C       Do QN-DMS
        IPOTmp = IGetRW(NTT,0)
        IFOTmp = IGetRW(NTT,0)
        Call FileIO(1,-IFOTmp,NTT,V(IFNew),0)
        Call FileIO(2,-IRwPO,NTT,V(IScr1),MaxGDI*NTT)
        Call FileIO(1,-IPOTmp,NTT,V(IScr1),0)
        Call QNDMS(IOut,2,NBasis,IRwP,0,IPOTmp,IFOTmp,IRwV,Ne,V,
     $    MDV,DMSErr,IBas)
        Call FreeRW(IPOTmp)
        Call FreeRW(IFOTmp)
      else
C       Do diagonalization
        IEig = IScr2
        IE = IEig + NBasis
        IE2 = IE + NBasis
        NE2 = (MDV-IE2+1)/NBasis
        Call DiagDN(IOut,0,0,V(IFNew),V(IScr1),V(IEig),NBasis,V(IE),NE2,
     $    V(IE2),NBasis,.False.)
        Call FormP2(.True.,1,NBasis,Ne,0,One,V(IScr1),V(IScr2))
        Call AMove(NTT,V(IScr2),V(IFNew))
        Call FileIO(2,-IRwV,NBSq,V(IScr1),0)
        Call OTran(V(IFNew),V(IScr1),V(IScr2),NBasis,NBasis,0,1,
     $    V(IScr3),MDV-IScr3)
        Call FileIO(1,-IRwP,NTT,V(IFNew),0)
        endIf
      Call TStamp(1,'End of SimFP')
      Return
      End
*Deck SimVec
      Subroutine SimVec(IOut,IPrint,IBas,NBasis,IRwV,IRwP,IRwFP,IRwF,
     $  IRwPO,MaxGDI,NPDIIS,RMSdP,PDIS,EDIS,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Transform the Fock matrices to orthogonal basis for FoSimult.
C     IRwPO  : first MaxGDI*NTT block saves previous Fock matrices in
C              orthonormal basis second MaxGDI*NTT block saves previous
C              density matrices in orthonormal basis
C     IRwFP  : saves the previous gradients
C
      Dimension V(MDV),IDummy(1),PDIS(NPDIIS,NPDIIS),EDIS(NPDIIS,NPDIIS)
      Save One
      Data One/1.0D0/
 1000 Format(' ',A,'DIIS matrix:')
C
C     Initilization
C
      If(IBas.ne.1) Call GauErr('IBas=1 assumed in SimVec.')
      NTT = (NBasis*(NBasis+1))/2
      NBSq = NBasis**2
      If(NPDIIS.eq.1) then
        Call ConDDF(IRwPO,2*MaxGDI*NTT)
        Call ConDDF(IRwFP,MaxGDI*NTT)
        endIf
      IPN   = 1
      IFN   = IPN   + NBSq
      IVX   = IFN   + NBSq
      IScr1 = IVX   + NBSq
      IScr2 = IScr1 + NBSq
      IEnd  = IScr2 + NBSq
      MDV1  = MDV   - IEnd + 1
      Call TstCor(Max(MaxGDI*NTT*2,IEnd),MDV,'SimVec-1')
      IScr  = IScr1
      Call TstCor(IScr+NPDIIS*NTT,MDV-IScr+1,'SimVec-2')
C
C     Transform the Fock matrix into the orthonormal basis
C     Store Fock matrices in orthonormal basis
C
      NMove = NPDIIS - 1
      Call FileIO(2,-IRwF,NTT,V(IFN),0)
      Call FileIO(2,-IRwV,NBSq,V(IVX),0)
      Call OTran(V(IFN),V(IVX),V(IScr1),NBasis,NBasis,0,IBas,V(IEnd),
     $  MDV1)
      If(NPDIIS.ge.2) then
        Call FileIO(2,-IRwPO,NTT*NMove,V(IScr),0)
        Call FileIO(1,-IRwPO,NTT*NMove,V(IScr),NTT)
        endIf
      Call FileIO(1,-IRwPO,NTT,V(IFN),0)
C
C     Transform the density matrix into the orthonormal basis
C     Store density matrices in orthonormal basis
C
      Call FileIO(2,-IRwP,NTT,V(IPN),0)
      Call FileIO(2,-IRwV,NBSq,V(IVX),NBSq)
      Call OTran(V(IPN),V(IVX),V(IScr1),NBasis,NBasis,0,IBas,V(IEnd),
     $  MDV1)
      If(NPDIIS.ge.2) then
        Call FileIO(2,-IRwPO,NTT*NMove,V(IScr),MaxGDI*NTT)
        Call FileIO(1,-IRwPO,NTT*NMove,V(IScr),MaxGDI*NTT+NTT)
        endIf
      Call FileIO(1,-IRwPO,NTT,V(IPN),MaxGDI*NTT)
C
C     Calculate and store the force with respect to the
C     density matrix
C
      Call FileIO(2,-IRwPO,NTT,V(IPN),MaxGDI*NTT)
      Call FileIO(2,-IRwPO,NTT,V(IFN),0)
      Call dEdP(V(IEnd),MDV1,Junk,IDummy,V(IPN),V(IFN),V(IVX),V(IScr1),
     $  V(IScr2),NBasis,NTT,.False.,.True.)
      Call AScale(NTT,-One,V(IVX),V(IVX))
      If(NPDIIS.ge.2) then
        Call FileIO(2,-IRwFP,NTT*NMove,V(IScr),0)
        Call FileIO(1,-IRwFP,NTT*NMove,V(IScr),NTT)
        endIf
      Call FileIO(1,-IRwFP,NTT,V(IVX),0)
C
C     Compute approximate Hessian using only the diagonal elements
C
      Call FileIO(2,-IRwPO,NTT,V(IPN),MaxGDI*NTT)
      Call FileIO(2,-IRwPO,NTT,V(IFN),0)
      Call PHess(IOut,NBasis,V(IPN),V(IFN),V(IVX),V(IScr1))
C
C     Begin to construct the RFO-DIIS matrix for the density step
C
      Call FileIO(2,-IRwFP,NTT*NPDIIS,V(IScr),0)
      Do 15 I = 1, NPDIIS
        IShift = IScr + (I-1)*NTT
   15   Call ADiv(NTT,V(IShift),V(IVX),V(IShift))
      Do 25 I = 1, NPDIIS
        Do 25 J = 1, I
          IShftI = IScr + (I-1)*NTT
          IShftJ = IScr + (J-1)*NTT
          PDIS(I,J) = SCFTrc(V(IShftI),V(IShftJ),NBasis,1) + PDIS(I,J)
   25     PDIS(J,I) = PDIS(I,J)
      If(IPrint.ge.1) then
        Write(IOut,1000) 'Density RFO-'
        Call OutMat(IOut,1,PDIS,NPDIIS,NPDIIS,NPDIIS,NPDIIS)
        endIf
C
C     Prepare En-DIIS vector matrix e(i,j)=Tr[F(i)P(j)]
C
      IFN  = 1
      IPN  = IPN + NPDIIS*NTT
      Call FileIO(2,-IRwPO,NPDIIS*NTT,V(IFN),0)
      Call FileIO(2,-IRwPO,NPDIIS*NTT,V(IPN),MaxGDI*NTT)
      Do 30 I = 1, NPDIIS
        IFShft = IFN + (I-1)*NTT
        Do 30 J = 1, NPDIIS
          JPShft = IPN + (J-1)*NTT
          EDIS(I,J) = EDIS(I,J) + SCFTrc(V(IFShft),V(JPShft),NBasis,1)
   30   Continue
      If(IPrint.ge.1) then
        Write(IOut,1000) 'Density En-'
        Call OutMat(IOut,1,EDIS,NPDIIS,NPDIIS,NPDIIS,NPDIIS)
        endIf
      Return
      End
*Deck SNRCyc
      Subroutine SNRCyc(IOut,IPrint,IP,P,F,G,H,Xi,IMap,JMap,X,LX,N,NP,
     $  ItMax,NOcc,NPMax,BMat,IS,IAD1,IAD2,D,Q,R,QR,FastM,IScr,IScrQ1,
     $  IScrQ2,ScrQ1,ScrQ2,IBT,BT,BnFac,ISizCV,IRowNo,iBin,NBin,LBin,
     $  IXX,LIX,IMemOp,ZThr,Method,IFail)
      Implicit real*8(a-h,o-z)
C
C     Sparse Matrix Newton-Raphson Density Matrix Search Using Direct
C     Inversion in the Interative Space
C     All matrices are lower (or upper) triangular sparse symmetric matrices
C     with diagonal elements seperate linear arrays.  Map Maps elements
C     of the transposed matrices to elements of the untransposed matrices.
C
      Logical FastM,IniMA
      Dimension IP(N+1,2),P(N,2),F(N,2),G(N,2),H(N,2),Xi(N,2),IMap(*),
     $  JMap(*),X(N,NP,LX),BMat(*),IS(*),IAD1(*),IAD2(*),D(*),Q(*),
     $  R(*),QR(*),IScr(*),IScrQ1(*),IScrQ2(*),ScrQ1(*),ScrQ2(*),
     $  IBT(*),BT(*),BnFac(*),ISizCV(*),IRowNo(*),iBin(NBin+1,N,LBin),
     $  IXX(N,NP,LIX)
      Save Stop1, Stop2, Good, Reject
      Data Stop1/1.D-4/,Stop2/1.D-3/,Good/1.D-11/,Reject/1.D-3/
 1010 Format(' Idempotency after purify:  ',1PD20.6)
 1030 Format(/,' Density is no longer idempotent in QN-DMS,',
     $       ' switch to CG-DMS',/)
 1040 Format(' RMS Gradient  ',1PD20.6,1PD12.3)
 1050 Format(' MAX Gradient  ',1PD20.6,1PD12.3)
C
C     Initilization
C
      IOP  = 0
      IOG  = 0
      NZ   = IP(N+1,1) - 1
      Len  = NZ + N
      MaxD = ItMax + 1
C
C     QN-DMS Iterations
C
      Do 100 NSaved = 1, ItMax
C
C       Calculating dEdP
C
        If(FastM) then
          Call IMove(Len+1,IP,IScr)
          Call Fill2(IP,P,IScrQ1,ScrQ1,N,NZSq)
          Call SpTrp(IScrQ1,ScrQ1,IScrQ2,ScrQ2,N)
          Call SpTrp(IScrQ2,ScrQ2,IScrQ1,ScrQ1,N)
          Call CheckF(IScrQ1,ScrQ1,X,N,0,NZSq,'ScrQ1 Matrix',1,1)
          endIf
        If(FastM) then
          RJunk = FxFnc3(IP,IScr,IScrQ1,P,F,Xi,ScrQ1,ScrQ2,X,N,ZThr,
     $      IMemOp,IBT,BT,IScrQ2,IMap,JMap,BnFac,ISizCV,IRowNo,
     $      IXX(1,1,1),IXX(1,1,2),IXX(1,1,3),IXX(1,1,4),IXX(1,1,5),
     $      X(1,1,2),X(1,1,3),X(1,1,4),X(1,1,5),iBin(1,1,1),iBin(1,1,2),
     $      NBin,Method,NP) / GFloat(N)
        else
          RJunk = FxFunc(IP,P,F,Xi,IMap,JMap,X,X(1,1,2),X(1,1,3),
     $      X(1,1,4),X(1,1,5),X(1,1,6),X(1,1,7),X(1,1,8),N,NP)/GFloat(N)
          endIf
        RMSG = Sqrt(SProd(Len,Xi,Xi))/GFloat(Len)
        GMax = ArrMax(Xi,Len,.true.)
        If(IPrint.ge.2) then
          Write(IOut,1040) RMSG,Stop1
          Write(IOut,1050) GMax,Stop2
          endIf
        If(GMax.le.Stop2 .and. RMSG.le.Stop1) Goto 200
        Call ANeg(Len,Xi,G)
C
C       Calculating Approximate Hessian
C
        Call SpHess(IOut,IPrint,IP,P,F,H,X,X(1,1,2),X(1,1,3),N,NP)
C
C       Do Newton-Raphson density matrix search using DIIS
C
        Call NRDIS(IOut,IPrint,P,G,H,Xi,BMat,MaxD,NSaved,N,NZ,IOP,
     $    IOG,IS,IAD1,IAD2,D,Q,R,QR)
C
C       Mcweeny Purification
C
        Demp = Abs(SpTrac(P,P,N,NZ)-GFloat(NOcc))/GFloat(N)
        If(Demp.le.Good) Goto 200
        If(FastM) then
          Do 10 I = 1, NPMax
            IniMA = .True.
            Call FsPurf(IOut,IPrint,IP,IScr,IBT,P,Xi,ScrQ1,ScrQ2,BT,
     $        IMap,JMap,BnFac,ISizCV,IRowNo,IXX(1,1,1),IXX(1,1,2),
     $        IXX(1,1,3),IXX(1,1,4),IXX(1,1,5),X(1,1,1),X(1,1,2),
     $        X(1,1,3),X(1,1,4),iBin(1,1,1),iBin(1,1,2),NBin,NOcc,
     $        N,NZ,NP,ZThr,Method,IniMA)
            Fac = GFloat(NOcc)/ArrSum(N,P)
            Call AScale(Len,Fac,P,P)
            Demp = Abs(SpTrac(P,P,N,NZ)-GFloat(NOcc))/GFloat(N)
            If(IPrint.ge.2) Write(IOut,1010) Demp
            If(Demp.le.Good) Goto 30
   10       Continue
        else
          Do 20 I = 1, NPMax
            Call FxPurf(IP,P,Xi,IMap,JMap,X,X(1,1,2),X(1,1,3),X(1,1,4),
     $        X(1,1,5),N,NP)
            Fac = GFloat(NOcc)/ArrSum(N,Xi)
            Call AScale(Len,Fac,Xi,P)
            Demp = Abs(SpTrac(P,P,N,NZ)-GFloat(NOcc))/GFloat(N)
            If(IPrint.ge.2) Write(IOut,1010) Demp
            If(Demp.le.Good) Goto 30
   20       Continue
          endIf
   30   If(Demp.ge.Reject) then
          IFail = 1
          Write(IOut,1030)
          Call FileIO(2,-IOP,Len,H,(NSaved-1)*Len)
          Call AMove(Len,H,P)
          Goto 200
          endIf
  100   Continue
  200 Call FreeRW(IOP)
      Call FreeRW(IOG)
      Return
      End
*Deck SNRDMS
      Subroutine SNRDMS(IOut,IPrint,PrntDM,cc,maxcor,N,ZThr,Dens,NE,
     $  IOpDMe,JCycle,jip,ip,if,NZ,ibegin,iend,IMemO,IFail)
      Implicit Real*8(A-H,O-Z)
      Dimension cc(maxcor)
      Logical PrntDM,FastM
      Save LX, LIX, LBin
      Data LX/11/, LIX/5/, LBin/2/
 1010 Format(' Zero threshold:                    ',1pd12.2)
 1020 Format(' NRDMS Convergence criteria:        ',1pd12.2)
 1030 Format(' Maximum number of NRDMS iterations:',i5)
 1040 Format(' Maximum purification cycles:       ',i5)
 1050 Format(/' ********* NRDMS cycle ',i4/)
C
C     Set convergence threshold and other useful stuff
C
      IFail = 0
      If(NE.eq.0) then
        Call AClear(NZ,CC(IP))
        Return
        endIf
      istart = ibegin
      ItMax = Mod(IOpDMe,100)
      If(ItMax.eq.0) ItMax = 4
      NPMax = Mod(IOpDMe,10000)/100
      If(NPMax.eq.0) NPMax = 3
      If(PrntDM) then
        Write(IOut,1010) ZThr
        Write(IOut,1020) Dens
        Write(IOut,1030) Itmax
        Write(IOut,1040) NPmax
        endIf
      If(IPrint.ge.1) write(IOut,1050) JCycle
C
C     Decided the multiplication method to use
C
      IMem   = IMemO
      Method = Mod(IMemO,10)
      Call FSIni1(cc,maxcor,istart,iBnFac,iSizCV,iend,cc(jip),cc(ip),
     $  cc(jip),cc(ip),N,NBin,IMem,Method,ZThr,.true.,.true.)
      istart = iend
      If(Method.eq.2) FastM = .true.
C
C     Allocate memory for all matrices used by NR-DMS
C
      NP    = NProc(0)
      Len   = NZ    + N
      LBMat = ItMax + 1
      ig    = istart
      ih    = ig    + Len
      ixi   = ih    + Len
      jimap = ixi   + Len
      jmap  = jimap + INTOWP(Len+1)
      ibmat = jmap  + INTOWP(NZ)
      ix    = ibmat + LBMat*LBMat
      is    = ix    + LX*N*NP
      iad1  = is    + INTOWP(2*LBMat)
      iad2  = iad1  + INTOWP(LBMat)
      id    = iad2  + INTOWP(LBMat)
      iq    = id    + LBMat
      ir    = iq    + LBMat*LBMat
      iqr   = ir    + LBMat*LBMat
      iend  = iqr   + LBMat*LBMat
      If(FastM) then
        LenQ = 2*NZ + N
        LISL = InToWP(NZ+N+1)
        LISQ = InToWP(2*NZ+2*N+1)
        iscrq1 = iend
        iscrq2 = iscrq1 + LenQ
        jiscr  = iscrq2 + LenQ
        jisrq1 = jiscr  + LISQ
        jisrq2 = jisrq1 + LISQ
        jibt   = jisrq2 + LISQ
        ibt    = jibt   + LISL
        iRowNo = ibt    + Len
        jibin  = iRowNo + INTOWP(N)
        jix    = jibin  + INTOWP(LBin*N*(NBin+1))
        iend   = jix    + INTOWP(LIX*N*NP)
        endIf
      Call TstCor(iend,maxcor,'SNRDMS')
C
C     Fixed version - all matrix forms held constant
C     Make Fock and Density share common element
C
      Call MkMap(cc(jip),cc(jimap),cc(jmap),N)
      Call CheckF(cc(jimap),cc(jmap),cc(ix),N,0,NZ,'Map',-6,IPrint)
      Call SNRCyc(IOut,IPrint,cc(jip),cc(ip),cc(if),cc(ig),cc(ih),
     $  cc(ixi),cc(jimap),cc(jmap),cc(ix),LX,N,NP,ItMax,NE,NPMax,
     $  cc(ibmat),cc(is),cc(iad1),cc(iad2),cc(id),cc(iq),cc(ir),
     $  cc(iqr),FastM,cc(jiscr),cc(jisrq1),cc(jisrq2),cc(iscrq1),
     $  cc(iscrq2),cc(jibt),cc(ibt),cc(ibnfac),cc(isizcv),cc(irowno),
     $  cc(jibin),NBin,LBin,cc(jix),LIX,IMem,ZThr,Method,IFail)
      Call CheckF(cc(jip),cc(ip),cc(ix),N,0,NZ,'Conv. density',5,IPrint)
      Return
      End
*Deck SolPru
      Subroutine SolPru(IOut,IPrint,SymCut,NTot,NTotR,PtThr,XYZ,Wgt,
     $  ISubSt,ISlMap,Neighb,IDelIf,Scr,IScr,Scr1,List,MDVI)
      Implicit Real*8(A-H,O-Z)
      Logical TooCls
      Dimension XYZ(NTot,3), Wgt(NTot), ISubSt(NTot), ISlMap(NTot),
     $  IDelIf(NTot), Neighb(NTot), Scr(NTot), IScr(NTot), Scr1(NTot),
     $  List(*)
      Save One, Zero, OneHun
      Data One/1.0d0/, Zero/0.0d0/, OneHun/100.0d0/
C
C     Delete points in input points that are too close together
C     without destroying symmetry, creating list of who represents
C     whom.
C
C     Throughout, r^2 is used instead of r to save time.
C     Scr and IScr may be equivalenced.
C
C     Output:
C       ISubst(I) = point in large space that represents I.
C       ISlMap(I) = point in small space that represents I, or 0.
C       Neighb(II) = point in large space that is represented by II.
C
 9000 Format(1x,'Point ',i8,' (',i8,') has ',i8,' neighbours:',
     $  3(1x,i8),:,/,(7(1x,i8)))
 9010 Format(1x,'Deleted pts ',i8,' and ',i8,': same weight')
 9020 Format(1x,i8,' & ',i8,' will be deleted if ',i8,' isn''t')
 9030 Format(1x,'Looking at pt ',i8,', ne ',i8,' (',i8,')')
 9040 Format(1x,'Confirmed delete for:',i8)
 9050 Format(1x,'Found high weight neighb ',i8,i8)
 9060 Format(1x,'Couldn''t find substitute for point',i8,
     $  ' - this grid is weird.')
 9070 Format(1x,'Found desperation neighb ',i8,i8)
 9080 Format(1x,'Point ',i8,' is represented by point ',i8)
 9090 Format(1x,'Point ',i8,' is subspace point ',i8)
 9100 Format(1x,i8,' points reduced to ',i8,',',' PtThr=',f5.3)
      Call TStamp(1,'Enter SolPru')
C
C     Find all "too close" neighbours for each point. Each point I
C     has Neighb(I) neighbours whose indices are stored at
C     List(ISlMap(I)+1..Neighb(I)) in order of descending Wgt.
C
      PtTh1 = PtThr*PtThr
      LstInd = 0
      Call IClear(NTot,Neighb)
      Do 30 I = 1, NTot
        IDelIf(I) = I
        ISubst(I) = I
        ISlMap(I) = LstInd
        If(PtTh1.ge.SymCut) then
          Do 20 J = 1, NTot
            If(I.ne.J) then
              Distx = XYZ(I,1) - XYZ(J,1)
              Disty = XYZ(I,2) - XYZ(J,2)
              Distz = XYZ(I,3) - XYZ(J,3)
              Dist = Distx**2 + Disty**2 + Distz**2
              TooCls = (Dist.lt.PtTh1)
              If(TooCls) then
                LstInd = LstInd + 1
                If(LstInd.gt.MDVI) Call GauErr('Out of Mem: SolPru')
                Neighb(I) = Neighb(I) + 1
                List(ISlMap(I)+Neighb(I)) = J
                endif
              Endif
   20       Continue
          Endif
   30   Continue
C
C     Create a sort index of canonical order:
C       highest (Neighb*large# + Wgt) = first
C
      WtMax = ArrMax(Wgt,NTot,.false.)
      Mag = 10 ** (Int(Log10(WtMax)) + 1)
      Do 40 I = 1, NTot
   40   Scr1(I) = GFloat(Neighb(I)*Mag) + Wgt(I)
      Call SortVc(NTot,Scr1,IScr)
C     Debug print out neighbour data
      If(IPrint.ge.3) then
        Do 50 I1 = 1, NTot
          I = IScr(NTot+1-I1)
          Write(Iout,9000)
     $            I1,I,Neighb(I),(List(ISlMap(I)+Ne),Ne=1,Neighb(I))
   50     Continue
        Endif
C
C     Preliminary prune. Remove points that must be deleted for symm
C     reasons, and mark those that may be deleted.
C
      Do 80 I = 1, NTot
        Do 70 Ne = 1, Neighb(I)
          JJ = List(ISlMap(I)+Ne)
          If(IDelIf(JJ).ne.JJ) goto 70
C         First delete all points that have a neighbour with the
C         same weight as they do.
          If(Abs(Wgt(JJ)-Wgt(I)).lt.SymCut) then
            IDelIf(I) = 0
            IDelIf(JJ) = 0
            If(IPrint.ge.3)
     $        Write(IOut,9010)I,JJ
            goto 70
            Endif
C         Then mark neighbours with same weight (and distance from I)
C         as each other to be deleted if point I isn't.
          Do 60 Ne1 = 1, Neighb(I)
            If(Ne1.ne.Ne) then
              KK = List(ISlMap(I)+Ne1)
              If(Abs(Wgt(KK)-Wgt(JJ)).lt.SymCut) then
                Dist1x = XYZ(I,1) - XYZ(JJ,1)
                Dist1y = XYZ(I,2) - XYZ(JJ,2)
                Dist1z = XYZ(I,3) - XYZ(JJ,3)
                Dist1 = Dist1x**2 + Dist1y**2 + Dist1z**2
                Dist2x = XYZ(I,1) - XYZ(KK,1)
                Dist2y = XYZ(I,2) - XYZ(KK,2)
                Dist2z = XYZ(I,3) - XYZ(KK,3)
                Dist2 = Dist2x**2 + Dist2y**2 + Dist2z**2
                If(Abs(Dist1-Dist2).lt.SymCut) then
                  IDelIf(JJ) = I
                  IDelIf(KK) = I
                  If(IPrint.ge.3) Write(Iout,9020)JJ,KK,I
                  goto 70
                  Endif
                Endif
              Endif
   60       Continue
   70     Continue
   80 Continue
C
C     Now assign points to where they're finally going.
C
      Do 100 I1 = 1, NTot
        I = IScr(NTot+1-I1)
        Targwt = Zero
        TargDs = Zero
        If(Neighb(I).eq.0) goto 110
C       Find a reason to keep this point (I).
        Do 90 Ne = 1, Neighb(I)
          JJ = List(ISlMap(I)+Ne)
          If(IPrint.ge.3) Write(IOut,9030)I,Ne,JJ
C         Has target neighbour already been deleted or marked for such?
          If(ISubst(JJ).ne.JJ.or.IDelIf(JJ).ne.JJ) then
            If(IPrint.ge.3) Write(Iout,*)'Been deleted'
            goto 90
            Endif
C         Target must have higher Neighb than us, or higher Wgt
C         if equal.
          If(Neighb(JJ).lt.Neighb(I)) then
            If(IPrint.ge.3) Write(Iout,*)'Less Neighbs'
            goto 90
          ElseIf(Neighb(JJ).eq.Neighb(I)) then
            If(Wgt(JJ).lt.Wgt(I)) then
              If(IPrint.ge.3)
     $          Write(Iout,*)'Same Neighbs, less weight'
              goto 90
              Endif
            Endif
C         Must have highest weight of possible candidates.
          If(Targwt-Wgt(JJ).gt.SymCut) then
            If(IPrint.ge.3) Write(Iout,*)'Previous was higher'
            goto 90
            endif
          If(Abs(Wgt(JJ)-Targwt).lt.SymCut) then
C         If Weights are equal (neighbours with equal weights, but
C         not equidistant from I) then choose closest.
            Distx = XYZ(I,1) - XYZ(JJ,1)
            Disty = XYZ(I,2) - XYZ(JJ,2)
            Distz = XYZ(I,3) - XYZ(JJ,3)
            Dist2 = Distx**2 + Disty**2 + Distz**2
            If(Dist2.gt.TargDs.and.TargDs.ne.Zero) then
              If(IPrint.ge.3)
     $          Write(IOut,*)'Previous was same, but closer'
              goto 90
              Endif
            Endif
C         We've found no reason not to choose to this neighbour.
          If(IPrint.ge.3) Write(Iout,*)'Chose this neighbour'
          ISubst(I) = JJ
          Targwt = Wgt(JJ)
          Distx = XYZ(I,1) - XYZ(JJ,1)
          Disty = XYZ(I,2) - XYZ(JJ,2)
          Distz = XYZ(I,3) - XYZ(JJ,3)
          TargDs = Distx**2 + Disty**2 + Distz**2
   90     Continue
  100 Continue
  110 Continue
C
C     Now delete any points that were to be deleted if a neighbour wasn't.
C
      Do 120 I = 1, NTot
        If(IDelIf(I).ne.I.and.IDelIf(I).ne.0
     $      .and.ISubst(IDelIf(I)).ne.IDelIf(I)) then
          If(IPrint.ge.3) Write(IOut,9040) I
          IDelIf(I) = 0
          Endif
  120   Continue
C
C     Delete points marked for deletion that weren't.
C
      Do 170 I = 1, NTot
        If(ISubst(I).eq.I.and.(IDelIf(I).eq.0)) then
          TargWt = Zero
          TargDs = OneHun
          Do 130 Ne = 1, Neighb(I)
            JJ = List(ISlMap(I)+Ne)
            If(ISubst(JJ).ne.JJ.or.IDelIf(JJ).eq.0) goto 130
C           Now only must be nearest highest weight of possibles.
C           (cannot be more than one of these - already removed)
            If(Abs(Wgt(JJ)-TargWt).gt.SymCut.and.Wgt(JJ).lt.TargWt)
     $        goto 130
            Distx = XYZ(I,1) - XYZ(JJ,1)
            Disty = XYZ(I,2) - XYZ(JJ,2)
            Distz = XYZ(I,3) - XYZ(JJ,3)
            Dist = Distx**2 + Disty**2 + Distz**2
            If(Dist.gt.TargDs) goto 130
            If(IPrint.ge.3) Write(IOut,9050)I,JJ
            ISubst(I) = JJ
            TargWt = Wgt(JJ)
            TargDs = Dist
  130       Continue
C         Check to see if we fixed it. If not then final criteria:
C         pick closest, undeleted, unique point or bail.
          If(ISubst(I).eq.I) then
            Call ASet(NTot,OneHun,Scr)
            Do 140 J = 1, NTot
              If(ISubst(J).ne.J.or.IDelIf(J).eq.0) goto 140
              Distx = XYZ(I,1) - XYZ(J,1)
              Disty = XYZ(I,2) - XYZ(J,2)
              Distz = XYZ(I,3) - XYZ(J,3)
              Scr(J) = Distx**2 + Disty**2 + Distz**2
  140         Continue
            LoopCn = 0
            Indi = 0
            Distch = OneHun-One
  150       LoopCn = LoopCn + 1
            DistO = Distch
            Distch = OneHun-One
            If(LoopCn.gt.NTot) then
              Write(IOut,9060) I
              Call GauErr('FUBAR in grid.')
              Endif
            Do 160 J = 1, NTot
              If(Abs(Scr(J)-Distch).lt.SymCut) then
                Scr(J) = OneHun
                Scr(Indi) = OneHun
                goto 160
                Endif
              If(Scr(J).lt.Distch) then
                Indi = J
                Distch = Scr(J)
                Endif
  160         Continue
            If(Distch.ne.DistO) goto 150
            ISubst(I) = Indi
            If(IPrint.ge.3) Write(IOut,9070)I,Indi
            Endif
          Endif
  170   Continue
C
C     Make sure we only have single level indirection. i.e. no i->j->k->l...
C     and create final indices.
C
      Call IClear(NTot,Neighb)
      NTotR = 0
      Do 190 I = 1, NTot
        IPlarp = ISubSt(I)
        JOld = I
  180   J = ISubSt(JOld)
        If(J.ne.JOld) then
          JOld = J
          Goto 180
          Endif
        ISubSt(I) = J
        If(IPrint.ge.3) Write(IOut,9080)I,J
        If(IPlarp.ne.J) Call GauErr
     $    ('Had >1 level indirection. Error in solvent pruning.')
        If(ISubSt(I).eq.I) then
          NTotR = NTotR + 1
          ISlMap(I) = NTotR
          If(IPrint.ge.3) Write(IOut,9090)I,NTotR
        Else
          ISlMap(I) = 0
          Endif
        If(ISlMap(I).ne.0) Neighb(ISlMap(I)) = I
  190   Continue
      If(IPrint.ge.2) Write(IOut,9100)NTot,NTotR,PtThr
      Call TStamp(2,'Exit SolPru')
      Return
      End
*Deck SpAdAM
      Subroutine SpAdaM(cc,maxcor,IA,AN,IB,BN,X,istart,jIC,iCN,iend,N,
     $  NZC,thresh,IType)
      Implicit Real*8(A-H,O-Z)
C
C     Add two sparse matrices:   C = SA*A + SB*B
C     IType = 0    IC,JC only
C     IType = 1    IC,JC,CN only
C     IType = 2    IC,JC,CN, and CD
C
      Dimension cc(maxcor),IA(N+1,2),AN(N,2),IB(N+1,2),BN(N,2),X(N)
      Logical DoDiag
C
      DoDiag = IType.eq.2
      LenD = 0
      if(DoDiag) LenD = N
      jIC  = istart
      call SpAddS(IA,IB,cc(jIC),X,N,NZC)
      jICT = jIC  + INTOWP(N+1+NZC)
      iend = jICT + INTOWP(N+1+NZC)
      call TstCor(iend,maxcor,'SpAddM 1')
      call SpTrpS(cc(jIC),cc(jICT),N)
      call SpTrpS(cc(jICT),cc(jIC),N)
      iend = jICT
      if(IType.le.0) return
      iCN1 = iend
      iend = iCN1 + NZC+LenD
      call TstCor(iend,maxcor,'SpAdaM 2')
      if(IType.eq.2) call AAdd(N,AN,BN,cc(iCN1))
      call SpAdAN(IA,AN(1,2),IB,BN(1,2),cc(jIC),cc(iCN1+LenD),X,N)
      call compct(cc(jIC),cc(iCN1),thresh,N,NZC,1,DoDiag)
      iCN  = jIC  + INTOWP(N+1+NZC)
      IEnd = iCN  + NZC+LenD
      call TstCor(iend,maxcor,'SpAddM 3')
      call Amove1(NZC+LenD,0,(iCN-iCN1),cc(iCN1))
      Return
      End
*Deck SpAdAN
      Subroutine SpAdAN(IA,AN,IB,BN,IC,CN,X,N)
      Implicit Real*8(A-H,O-Z)
C
C     Sparse C = A + B, values only.
C
      Dimension IA(N+1,2), AN(*), IB(N+1,2), BN(*), IC(N+1,2), CN(*),
     $  X(*)
C
      Do 100 I = 1, N
        ISA = IA(I,1)
        NIA = IA(I+1,1) - ISA
        ISB = IB(I,1)
        NIB = IB(I+1,1) - ISB
        ISC = IC(I,1)
        NIC = IC(I+1,1) - ISC
        Call AScatZ(NIC,X,IC(ISC,2))
        Call ASctIn(NIA,X,IA(ISA,2),AN(ISA))
        Call ASctIn(NIB,X,IB(ISB,2),BN(ISB))
        Call AGathr(NIC,CN(ISC),X,IC(ISC,2))
  100   Continue
      Return
      End
*Deck SpAXMt
      Subroutine SpAXMt(N,NNew,IMult,RMult,IRwX,IRwAX,V,MDV)
      Implicit Real*8 (A-H,O-Z)
      Dimension IMult(3),RMult(*),V(MDV)
C
      jIF  = IMult(1)
      iF   = IMult(2)
      IRow = IMult(3)
      IX = 1
      IAX = IX + N*NNew
      IEnd = IAX + N*NNew
      Call TstCor(IEnd,MDV,'SpAXMt')
      Call Fileio(2,-IRwX,N*NNew,V(IX),0)
      Do 20 i = 1,NNew
        Call SpMvsM(RMult(jIF),RMult(iF),V(IX+(i-1)*N),V(IAX+(i-1)*N),
     $    RMult(IRow),N,N,.true.)
   20   Continue
      Call Fileio(1,-IRwAX,N*NNew,V(IAX),0)
      Return
      End
*Deck SpDIIS
      Subroutine SpDIIS(IOut,IPrint,cc,maxcor,ibegin,IRWFOA,IRwFOB,
     $  IRwPOA,IRwPOB,IRwFM1,IRwFM2,N,jcycle,MxSCFC,ZThr,IMemOp)
      Implicit real*8(a-h,o-z)
C
C     Retrieve and transform the Fock matrix
C
      Dimension CC(Maxcor),JJA(1)
      Logical DoDIIS,newf
      Save Ten,IQSiz2
      Save Zero,One,IOD1,IOERR1,IOD2,IOERR2,IOB,IQStrt,IQSize,JJA,iCForm
      Data Ten/10.d0/,Zero/0.d0/,One/1.d0/,JJA/0.d0/,iCForm/1/
C
C     Take a DIIS step if requested
C
      DoDIIS = mxscfc.gt.2
      NewF = IRWFM2.ne.0
      IUHF = 1
      If(IRwFOB.ne.0) IUHF = 2
      If(DoDIIS) then
        Call TStamp(1,'Top of DIIS.')
        If(JCycle.le.1) then
          IQstrt = 0
          IQsize = 0
          IQsiz2 = 0
          IOB = 0
          IOD1 = IGetRW(1,0)
          IOD2 = IGetRW(1,0)
          IOErr1 = IGetRW(1,0)
          IOErr2 = IGetRW(1,0)
          endIf
        If(Mod(iCForm,2).eq.1)Then
          IOD = IOD1
          IODT = IOD2
          IOErr = IOErr1
          IOErt = IOErr2
        else
          IOD = IOD2
          IODT = IOD1
          IOERR = IOERR2
          IOERT = IOERR1
          endif
C
C     If necessary, Change the matrix form of the previous DIIS steps
C
        ix = ibegin
        istrtX = ix + N
        Call SpGet(cc,cc,maxcor,IRWFM1,istrtX,jifm,ifm,iend,N,NZFM,0,
     $    .true.,0)
        istrtN = iend
        istart = iend
        If(NewF) Then
C
C     First form the union of the old and new forms
C
          Call SpGet(cc,cc,maxcor,IRWFM2,istrtN,jinfm,infm,iend,N,NZNFM,
     $          0,.true.,0)
          istart = iend
          Call SpAddM(cc,maxcor,cc(jifm),cc(ifm),cc(jinfm),cc(infm),
     $         cc(ix),istart,jinw,inw,iend,N,NZNW,One,One,Zero,0)
          NZNFM = NZNW
          jinfm = iStrtN
          iend = jinfm + intowp(N+NZNFM+1)
          Call IMove1((N+NZNFM+1),0,(jiNFM-jinw)*InTPWP(0),cc(jinw))
          Call SpStor(IRWFM1,cc(jiNFM),cc(iend),N,0,.true.,0)
C
C     Now Change the old matrices to the new form
C
          iCForm = iCForm + 1
          I = IOD
          IOD = IODT
          IODT = I
          I = IOERR
          IOERR = IOERT
          IOERT = I
          Leng1=mxscfc*(NZNFM+N)*IUHF
          Leng2=mxscfc*NZNFM*IUHF
          call ConDDF(IOD,Leng1)
          call ConDDF(IOERR,Leng2)
          ifm = iend
          infm = ifm  + NZFM*IUHF
          ifmd = infm + NZNFM*IUHF
          iend = ifmd + N*IUHF
          IQSiz2 = Min(IQStrt+IQSize-1,mxscfc)
          Call NwDISF(IODT,IOD,cc(jifm),cc(ifm),cc(jinfm),cc(infm),
     $         cc(ifmd),cc(ix),N,NZFM,NZNFM,2,IQSiz2,IUHF)
          Call NwDISF(IOERT,IOERR,cc(jifm),cc(ifm),cc(jinfm),cc(infm),
     $         cc(ifmd),cc(ix),N,NZFM,NZNFM,1,IQSiz2,IUHF)
          NZFM = NZNFM
          jiFM = iStrtX
          iend = jiFM + intowp(N+NZFM+1)
          istart = iend
          Call IMove1((N+NZFM+1),0,(jiFM-jinfm)*InTPWP(0),cc(jiNFM))
          endif
C
C     Form G = P*F - F*P
C
        iGa = istart
        If(IUHF.eq.2) Then
          iGb = iGa + N + NZFM
        else
          iGb = iGa
          endif
        iaftrG = iGb + N + NZFM
        Call SpGet(cc,cc,maxcor,IRWFOA,iaftrG,jif,if,iend,N,NZF,1,
     $    .true.,0)
        iaftrF = iend
        Call SpGet(cc,cc,maxcor,IRWPOA,iaftrF,jip,ip,iend,N,NZP,1,
     $    .true.,0)
        istart = iend
        Call FSMpyM(cc,maxcor,cc(jIF),cc(iF),cc(jIP),cc(iP),JJA(1),
     $       istart,JJA(1),JJA(1),iend,Zero,N,NZFm,cc(jIFM),cc(iGa),
     $       .true.,.true.,.true.,.false.,6,3,IMemOp)
        If(IUHF.eq.2) Then
          Call SpGet(cc,cc,maxcor,IRWFOB,iaftrF,jifb,ifb,iend,N,NZF,1,
     $      .true.,0)
          iaftrF = iend
          Call SpGet(cc,cc,maxcor,IRWPOB,iaftrF,jip,ip,iend,N,NZP,1,
     $      .true.,0)
          istart = iend
          Call FSMpyM(cc,maxcor,cc(jIFb),cc(iFb),cc(jIP),cc(iP),JJA(1),
     $       istart,JJA(1),JJA(1),iend,Zero,N,NZFm,cc(jIFM),cc(iGb),
     $       .true.,.true.,.true.,.false.,6,3,IMemOp)
          iGa2 = iGa + N
          Call AMove1(NZFm+N,0,(iGa2-iGa),cc(iGa))
          iFNw = iaftrF
          iFNwDA = iFNw
          iFNwDB = iFNwDA + N
          iFNwA = iFNwDB + N
          iFNwB = iFNwA + NZFM
          iend = iFNwB + NZFM
          Call TstCor(iend,maxcor,'SpDIIS 1')
          call AMove(N,cc(iF),cc(iFNwDA))
          call AMove(N,cc(iFb),cc(iFNwDB))
          call NwForm(cc(jIF),cc(iF+N),cc(jiFM),cc(iFNwA),cc(ix),N)
          call NwForm(cc(jIFb),cc(iFb+N),cc(jiFM),cc(iFNwB),cc(ix),N)
        else
          ifnw = iaftrF
          ix = ifnw + N + NZFM
          iend = ix + N
          Call TstCor(iend,maxcor,'SpDIIS 2')
          call AMove(N,cc(iF),cc(iFNw))
          call NwForm(cc(jIF),cc(iF+N),cc(jiFM),cc(iFNw+N),cc(ix),N)
          endif
C
C     Update the Fock matrix using Fixed Form DIIS
C
        Call FxDIIS(IOut,IPrint,cc(iFNw),cc(iGa),Ten,Ten,Errmax,mxscfc,
     $       jcycle,N*IUHF,NZFM*IUHF,IOD,IOErr,IOB,IQstrt,IQsize,'SCF')
        If(IUHF.eq.2) Then
          Call AMove(N,cc(IFNwDA),cc(iGa))
          Call AMove(N,cc(IFNwDB),cc(iGb))
          Call AMove(NZFM,cc(IFNwA),cc(iGa+N))
          Call AMove(NZFM,cc(IFNwB),cc(iGb+N))
          Call IMove(N+1+NZFM,cc(jiFM),cc(iFNw))
          Call Compct(cc(iFNw),cc(iGb),zthr,N,NZFMB,1,.true.)
          Call Compct(cc(jiFM),cc(iGa),zthr,N,NZFMB,1,.true.)
          Call SpStor(IRWFOB,cc(iFNw),cc(iGb),N,1,.true.,0)
          Call SpStor(IRWFOA,cc(jiFM),cc(iGa),N,1,.true.,0)
        else
          Call Compct(cc(jiFM),cc(iFNw),zthr,N,NZFM,1,.true.)
          Call SpStor(IRWFOA,cc(jiFM),cc(iFNw),N,1,.true.,0)
          Call TStamp(1,'After ScDIIS')
          endif
        endIf
      Return
      End
*Deck SpHess
      Subroutine SpHess(IOut,IPrint,IP,P,F,H,A,B,C,N,NP)
      Implicit Real*8(A-H,O-Z)
C
C     This subroutine forms an approximate density Hessian using only the
C     diagnol elements of P and F
C     The resulting Hessian has a same format as the density
C
C     Huvuv = (3-2Puu)Fvv + (3-2Pvv)Fuu - 4(PF)uu - 4(PF)vv
C
      Dimension IP(N+1,2),P(N,2),F(N,2),H(N,2),A(N),B(N),C(N,NP)
      Save Two,Three,Four,ThHess
      Data Two,Three,Four/2.d0,3.d0,4.d0/,ThHess/1.D-1/
 1000 Format(' The lowest eigenvalue is ',1PD20.6)
C
      Call SpDiag(IP,P,F,C,N,NP)
      Call AScale(N,-Four,C,C)
      Do 10 I = 1, N
        A(I) = Three - Two*P(I,1)
        B(I) = F(I,1)
   10   H(I,1) = Two*(A(I)*B(I)+C(I,1))
      NZ = IP(N+1,1) - 1
C$OMP Parallel Do If(N.gt.25) Schedule(Static,1) Default(Shared)
C$OMP+ Private(i,jj,j)
      Do 30 I = 1, N
        Do 30 JJ = IP(I,1), (IP(I+1,1)-1)
          J = IP(JJ,2)
   30     H(JJ,2) = A(I)*B(J) + B(I)*A(J) + C(I,1) + C(J,1)
C
C     Taking care of negative or small eigenvalues
C
      HMin = ArrMin(H,NZ+N,.false.)
      If(IPrint.ge.2) Write(IOut,1000) HMin
      If(HMin.lt.ThHess) then
        Diff = ThHess - HMin
        Do 40 I = 1, N
 40       H(I,1) = H(I,1) + Diff
        Do 50 I = 1, NZ
 50       H(I,2) = H(I,2) + Diff
        endIf
      Return
      End
*Deck SpIcAM
      Subroutine SpIcAM(cc,maxcor,IB,BN,X,istart,jIA,iAN,iend,N,NZA,
     $  Thresh,Iype)
      implicit real*8(a-h,o-z)
c
c     Increment a sparse matrix:  A  <--  A + B
c     Note, the original A and anything in memory after A may be destroyed!
c     Iype = 0    IA,JA only
c     Iype = 1    IA,JA,AN only
c     Iype = 2    IA,JA,AN, and AD
c
      Dimension cc(maxcor),IB(N+1,2),BN(N,2),X(N)
C
      Call SpAdaM(cc,maxcor,cc(jIA),cc(iAN),IB,BN,X,istart,jIC,iCN,iend,
     $  N,NZC,Thresh,Iype)
      Len = 0
      If(Iype.ge.1) Len = NZC
      if(Iype.ge.2) Len = NZC + N
      iAN  = jIA + InToWP(N+1+NZC)
      iend = iAN + Len
      Call TstCor(iend,maxcor,'SpIncM')
      Call SpAmvS(cc(jIC),cc(jIA),N,NZA)
      Call AMove(Len,cc(iCN),cc(iAN))
      Return
      End
*Deck SpMvsM
      Subroutine SpMvsM(IA,A,X,Y,IRow,NR,N,Init)
      Implicit Real*8(A-H,O-Z)
C
C     Sparse Matrix-vector multiply:   Y = A(T)*X, where X,Y are
C     full vectors and A is a square sparse matrix.
C     IROW specifies the nonzero elements of X.
C     Y is initialized to zero if Init=.true.
C
      Dimension IA(N+1,2),A(N,2),X(N),Y(N),IROW(NR)
      Logical Init
C
      Call SpMVC1(IA,A(1,2),X,Y,IRow,NR,N,Init)
      Call SpMVC2(IA,A(1,2),X,Y,IRow,NR,N,.False.)
      Call SpMVC3(A,X,Y,IRow,NR,N,.False.)
      Return
      End
*Deck StbUp1
      Subroutine StbUp1(NRI,NBD,NBasUD,NE,IRwNE,IRwC,IRwP,RotAng,DA,T,V)
      Implicit Real*8(A-H,O-Z)
C
C     Generate a density based on an instability eigenvector for one
C     spin-state.  RotAng is the angle through which the MOs should be
C     rotated.
C
      Dimension DA(*), T(NBasUD,NBasUD), V(NBD,NBD)
C
      NTTU = NRI*(NBasUD*(NBasUD+1))/2
      NTT  = NRI*(NBD*(NBD+1))/2
      NBNU = NRI*NBD*NBasUD
      Call StbDel(IRwNE,NE,NBasUD-NE,NRI,NBasUD,DA,T)
      If(NRI.ne.1) Call GauErr('StbUp1 cannot handle complex.')
      Call AScale(NTTU,RotAng,DA,DA)
      Call MxExp(NBasUD,DA,T,DA(1+NTTU),V)
      Call FileIO(2,-IRwC,NBNU,DA,0)
      Call MatMpy(NBD,NBasUD,NBasUD,DA,T,V)
      Call FormP(.True.,NBD,NBD,NE,V,DA)
      Call FileIO(1,-IRwP,NTT,DA,0)
      Return
      End
*Deck StbUpd
      Subroutine StbUpd(IOut,DoPrnt,IOpCl,JCycle,NBasis,NBsUse,NTT,NAE,
     $  NBE,IRwCA,IRwCB,IRwNCA,IRwNCB,IRwPA,IRwPB,DA,T,V)
      Implicit Real*8(A-H,O-Z)
C
C     Use an instability eigenvector to update the density matrix.
C
      Logical DoPrnt, CBfn, Spinor
      Dimension DA(NTT,2), T(NBsUse,NBsUse), V(NBasis,NBasis)
      Save Pt1, Pt2, Pt4, Pt8, One, OnePt2
      Data Pt1/0.1d0/, Pt2/0.2d0/, Pt4/0.4d0/, Pt8/0.8d0/, One/1.0d0/,
     $  OnePt2/1.2d0/
 1000 Format(' Recover instability eigenvector and rotate through',F8.2,
     $  ' radians.')
C
      Call DecoSC(IOpCl,NSpBlk,NRI,NDimBl,CBfn,NSpBlX,Spinor)
      PiOv4 = ATan(One)
      If(JCycle.le.2) then
        RotAng = Pt8*PiOv4
      else if(JCycle.le.4) then
        RotAng = OnePt2*PiOv4
      else if(JCycle.le.6) then
        RotAng = Pt4*PiOv4
      else if(JCycle.le.8) then
        RotAng = Pt2*PiOv4
      else if(JCycle.le.10) then
        RotAng = Pt1*PiOv4
      else
        Call GauErr('StbUpd is confused.')
        endIf
      If(Mod(JCycle,2).eq.1) RotAng = -RotAng
      If(DoPrnt) Write(IOut,1000) RotAng
      NBD = NBasis*NDimBl
      NBasUD = NBsUse*NDimBl
      NAEX = NAE + NBE*(NDimBl-1)
      Call StbUp1(NRI,NBD,NBasUD,NAEX,IRwNCA,IRwCA,IRwPA,RotAng,DA,T,
     $  V)
      If(NSpBlk.eq.2) Call StbUp1(NRI,NBD,NBasUD,NBE,IRwNCB,IRwCB,IRwPB,
     $  RotAng,DA,T,V)
      Return
      End
*Deck SurDer
      Subroutine SurDer(Iprint,NTot,NTotR,xyz,coord,Gnorm,qtot,bmf,
     $  NBas6D,Ntt6D,const1,rhatx,rhaty,rhatz,P,CutOff,solfa,Ffact,
     $  weight,MoBas,IopCl,NMO,Como,Q,ISc1,IScal,Alpha,Beta,GFact,
     $  Constc,EFG,Ell,Em,dwidrn,XYZOld,ISlMap,ISubst,MaxTyp,V,MDV)
      Implicit Real*8 (A-H,O-Z)
C
C     Surface derivative contributions to SCI-PCM Fock matrix.
C
      Logical MoBas
      Dimension xyz(NTot,3),qtot(NTot),bmf(NTot,ISc1),solfa(Ntt6D),
     $  P(Ntt6D),V(*),Coord(NTot,3),Ffact(NTot),weight(NTot),
     $  rhatx(NTot),rhaty(NTot),rhatz(NTot),Como(Nbas6D,NMO),
     $  GNorm(NTot,3),Q(NTot,ISc1*ISc1),Alpha(NTot,ISc1),Beta(NTot),
     $  GFact(NTot),EFG(NTot,10,2),Ell(NTot),Em(NTot,3),
     $  dwidrn(NTot,2),XYZOld(NTot,3),ISlMap(Ntot),ISubst(NTot)
C
C     Start loop in large space. Use large space values to compute
C     dndr in large space, contract this into weight derivs in dwidrn
C     and then put everything for this batch into the small space.
C     We must use the special order so that we can lump together
C     appropriate values.
C
      ISafe = IntMem(-1,0)
      MaxDer = 2
      iOffW = 1
      Do 5 Lp = MaxTyp, 0, -2
        LTot = LP + MaxDer
        LSize = ((LTot+1)*(LTot+2)*(LTot+3))/6
    5   iOffW = iOffW + LSize
      ITotal = iOffW + MaxTyp + MaxDer
      NVal = ITotal + 35 + 10*NBas6D
      MaxDo = Min((MDV-ISafe)/NVal,NTot)
C
      Do 20 I = 1, NTot, MaxDo
        NDo = Min(NTot-I+1,MaxDo)
        Idndr = 1
        IRho = Idndr+3*NDo
        IGRho = IRho+NDo
        IChi = IGRho+3*NDo
        IDChi = IChi+NDo*NBas6D
        ID2Rho = IDChi+3*NDo*NBas6D
        ID2Chi = ID2Rho+6*NDo
        IV = ID2Chi+6*NDo*NBas6D
        Call TstCor(IV,MDV,'SurDer')
        MDV1 = MDV - IV + 1
C
C       Within this loop, unpruned points are treated fully, while
C       pruned points only give terms for the weight derivatives.
C
        Do 10 InI = 1, NDo
          K = InI+I-1
          Coord(InI,1) = XYZOld(K,1)
          Coord(InI,2) = XYZOld(K,2)
          Coord(InI,3) = XYZOld(K,3)
   10     Continue
        Call TStamp(1,'Before RhGrRh')
        Call RhGrRh(IPrint,NDo,NTot,NTT6D,3,NBas6D,Coord,V(IRho),
     $    V(IGRho),V(ID2Rho),V(IChi),V(IDChi),V(ID2Chi),P,CutOff,NMO,
     $    IOpCl,Como,MoBas,V(IV),MDV1)
C
        Call TStamp(1,'Before dnkdrk')
        Call dnkdrk(NDo,I-1,V(IGRho),V(ID2Rho),V(Idndr),rhatx,rhaty,
     $    rhatz)
C
        Idwdr  = ID2Rho
        IEnd1 = Idwdr+NDo
        Idwdn  = IEnd1
        IEnd2 = Idwdn+3*NDo
        Call TStamp(1,'Before dwidr2')
        Call dwidr2(NDo,I-1,NTot,V(Idndr),dwidrn,rhatx,rhaty,
     $    rhatz,V(Idwdr),V(Idwdn))
C
        IDXYZ = IEnd2
        IGMag3 = IDXYZ+6*NDo
        IGMag = IGMag3+3*NDo
        IDcXYZ = IGMag+NDo
        Call dndpmn(NTot,NTotR,I-1,NDo,bmf,V(IChi),V(IDChi),
     $    V(IGRho),Nbas6D,Solfa,xyz,gnorm,V(IGMag),V(IGMag3),ffact,
     $    const1,weight,V(IDXYZ),Q,ISc1,Alpha,IScal,Constc,EFG,
     $    V(Idwdn),Em,V(IDcXYZ),ISlMap,ISubst)
C
        Idbdr = IEnd1
        Idpedr = Idbdr+NDo
        Idbndr = Idpedr+NDo
        Idpndr = Idbndr+NDo
        Call dbkdrk(NTot,I-1,NDo,gnorm,EFG,V(Idndr),
     $    V(Idbdr),V(Idpedr),rhatX,rhatY,rhatZ,V(Idbndr),
     $    V(Idpndr),ISlMap,IScal)
        IChMax = Idpndr+NDo
        Call SurFoc(NTot,NTotR,I-1,NDo,xyz,gnorm,V(Idndr),
     $    V(Idwdr),V(Idbdr),V(Idpedr),V(IChi),bmf,weight,qtot,ffact,
     $    nbas6D,ntt6D,Const1,V(IGRho),Ell,Solfa,CutOff,
     $    V(IChMax),rhatX,rhatY,rhatZ,GFact,Constc,Alpha,
     $    Beta,Q,ISc1,V(Idbndr),V(Idpndr),ISlMap,ISubst,IScal)
   20   Continue
C
      Return
      End
*Deck SurfBW
      Subroutine SurfBW(NOrg,NAng,NDimX,IA,CutPt,IBWFnc,COrg,XYZS,PBeck,
     $  AngWgt,Keep)
      Implicit Real*8(A-H,O-Z)
C
C     Generate Becke style weights for a surface integration.
C
      Logical Keep(*)
      Parameter (Expent=0.8d0)
      Dimension COrg(3,NOrg), XYZS(NDimX,3), PBeck(NAng,NOrg),
     $  AngWgt(*)
      Save Zero, Pt5, One, OnePt5
      Data Zero/0.0d0/, Pt5/0.5d0/, One/1.0d0/, OnePt5/1.5d0/
C
      If(IBWFnc.eq.1) Then
        Do 40 I = 1, NOrg
          Do 10 K = 1, NAng
   10       PBeck(K,I) = One
          Do 30 J = 1, NOrg
            If(J.ne.I) Then
              Distij=Sqrt((corg(1,i)-corg(1,j))**2+
     $                    (corg(2,i)-corg(2,j))**2+
     $                    (corg(3,i)-corg(3,j))**2)
              Do 20 K = 1, NAng
                Disti=Sqrt((xyzs(K,1)-corg(1,I))**2+
     $                     (xyzs(K,2)-corg(2,I))**2+
     $                     (xyzs(K,3)-corg(3,I))**2)
                Distj=Sqrt((xyzs(K,1)-corg(1,J))**2+
     $                     (xyzs(K,2)-corg(2,J))**2+
     $                     (xyzs(K,3)-corg(3,J))**2)
                DMu = (disti-distj)/distij
                f1 = onept5*dmu-pt5*(dmu**3)
                f2 = onept5*f1-pt5*(f1**3)
                f3 = onept5*f2-pt5*(f2**3)
                ess = Pt5*(One-f3)
   20           PBeck(K,I) = PBeck(K,I)*ess
            EndIf
   30     Continue
   40   Continue
      Else If(IBWFnc.eq.2) then
        Do 50 I=1,NOrg
          Do 60 K=1,NAng
            Dist=Sqrt((XYZS(K,1)-COrg(1,I))**2+
     $                (XYZS(K,2)-COrg(2,I))**2+
     $                (XYZS(K,3)-COrg(3,I))**2)
            Dist2=Dist**2
   60       PBeck(K,I)=Exp(-Expent*Dist2)
   50   Continue
      Else
        Call GauErr('surfbw: Invalid function for Becke weights.')
        EndIf
C
      Do 70 K = 1, NAng
        Denom = Zero
        Do 80 I = 1, NOrg
   80     Denom = Denom + PBeck(K,I)
        Temp = PBeck(K,IA) / Denom
        Keep(K) = Abs(Temp).gt.CutPt
   70   AngWgt(k) = AngWgt(k)*Temp
C
      Return
      End
*Deck SurFoc
      Subroutine SurFoc(NTot,NTotR,IDid,NDo,XYZ,GNorm,dndr,dwdr,
     $  dbdr,dphedr,Chi,bmf,wgt,qtot,fFact,NBasis,NTT,Const1,grho,Ell,
     $  SolFA,Cutoff,ChiMax,rhatx,rhaty,rhatz,GFact,Constc,Alpha,Beta,Q,
     $  ISc1,dbndr,dphndr,ISlMap,ISubst,IScal)
      Implicit Real*8 (A-H,O-Z)
C
c     Contruct all terms to be multiplied by drk/dPMu,Nu.
c     They are placed in EllJ.
c
c     drk/dPMu,Nu = -Chi(K,Mu)*Chi(k,Nu) / (GRho(K) . Rhat(k))
C
C
      Dimension XYZ(NTot,3),GNorm(NTot,3),dndr(NDo,3),dwdr(NDo),
     $  dphedr(*),dbdr(*),Chi(NDo,NBasis),bmf(NTot,ISc1),Ell(NTot),
     $  grho(NDo,3),Wgt(NTot),fFact(NTot),SolFA(NTT),qtot(NTot),
     $  ChiMax(NBasis),rhatx(*),rhaty(*),rhatz(*),GFact(NTot),
     $  Alpha(NTot,ISc1),Beta(NTot),Q(NTot,ISc1*ISc1),dbndr(*),
     $  dphndr(*),ISlMap(Ntot),ISubst(NTot)
      Save Zero,One,Three
      Data Zero/0.d0/,One/1.d0/,Three/3.d0/
C
C     Compute L = dE/drj
C
      Do 20 J = 1, NDo
        JDid=J+IDid
        JaDid = JDid
        JJDid = ISlMap(JDid)
        If(JJDid.eq.0) then
          JaDid = ISubst(JDid)
          JJDid = ISlMap(JaDid)
          endif
        SclFoc=Zero
        SclFcA=Zero
        SclFcC=Zero
        Ell(JDid) = Zero
        rjnj = rhatx(JaDid)*GNorm(JJDid,1)+rhaty(JaDid)*GNorm(JJDid,2)+
     $    rhatz(JaDid)*GNorm(JJDid,3)
        Do 30 K = 1, NTotR
          EllJ = Zero
          If(K.ne.JJDid) then
            Distx = XYZ(JJDid,1)-XYZ(K,1)
            Disty = XYZ(JJDid,2)-XYZ(K,2)
            Distz = XYZ(JJDid,3)-XYZ(K,3)
            Dist = One/(sqrt(distx**2+disty**2+distz**2))
            d3 = dist**3
            d5 = dist**5
            temp = Three*d5*(distx*rhatx(JaDid)+disty*rhaty(JaDid)+
     $             distz*rhatz(JaDid))
            rjnk = rhatx(JaDid)*GNorm(K,1)+rhaty(JaDid)*GNorm(K,2)+
     $             rhatz(JaDid)*GNorm(K,3)
            djknj = distx*GNorm(JJDid,1)+disty*GNorm(JJDid,2)+
     $            distz*GNorm(JJDid,3)
            djknk = distx*GNorm(K,1)+disty*GNorm(K,2)+distz*GNorm(K,3)
            djkdnj = distx*dndr(J,1)+disty*dndr(J,2)+distz*dndr(J,3)
C
            If(ISlMap(JDid).ne.0) then
c
c             2 index dA terms
c
              dAkjrj = -d3*rjnk + temp*djknk
              dAjjrj = -dAkjrj * Wgt(K) / Wgt(JJDid)
              dAjkrj = rjnj*d3 - temp*djknj + d3*djkdnj
              dAkkrj = -(dAjkrj*Wgt(JJDid) + d3*djknj*dwdr(J))/
     $                 Wgt(K)
c
c             Add all dA terms with appropriate factors
c
              EllJ = EllJ + q(JJDid,1) * bmf(K,1) * dAkjrj
              EllJ = EllJ + q(JJDid,1) * bmf(JJDid,1) * dAjjrj
              EllJ = EllJ + q(K,1) * bmf(JJDid,1) * dAjkrj
              EllJ = EllJ + q(K,1) * bmf(K,1) * dAkkrj
c
              If(IScal.eq.2)then
                EllJ = EllJ + Q(JJDid,2) * bmf(K,2) * dAkjrj
                EllJ = EllJ + Q(JJDid,2) * bmf(JJDid,2) * dAjjrj
                EllJ = EllJ + Q(K,2) * bmf(JJDid,2) * dAjkrj
                EllJ = EllJ + Q(K,2) * bmf(K,2) * dAkkrj
                SclFcA = SclFcA +
     $             (Alpha(K,1)*Q(JJDid,3)-Alpha(K,2)*Q(JJDid,4))*dAkjrj
     $           + (Alpha(JJDid,1)*Q(JJDid,3) -
     $                Alpha(JJDid,2)*Q(JJDid,4))*dAjjrj
     $           + (Alpha(JJDid,1)*Q(K,3)-Alpha(JJDid,2)*Q(K,4))*dAjkrj
     $           + (Alpha(K,1)*Q(K,3)-Alpha(K,2)*Q(K,4))*dAkkrj
C               d(|rj-rk|^-1) / drj terms.
                tmp = -d3*(distx*rhatx(JaDid)+disty*rhaty(JaDid)+
     $                 distz*rhatz(JaDid))
                SclFcC = SclFcC + tmp*
     $         ( -Beta(K)*Q(JJDid,3)*GFact(JJDid)
     $           -Beta(JJDid)*Q(K,3)*GFact(K)
     $           +Beta(JJDid)*Q(JJDid,3)*GFact(JJDid)*Q(K,4)/Q(JJDid,4)
     $           +Beta(K)*Q(K,3)*GFact(K)*Q(JJDid,4)/Q(K,4) )
                endif
C           Else ISlMap.eq.0
            Else
              dAjkrj = rjnj*d3 - temp*djknj + d3*djkdnj
              dAkkrj = -(dAjkrj*Wgt(JJDid) + d3*djknj*dwdr(J))/
     $                 Wgt(K)
              EllJ = EllJ + q(K,1) * bmf(K,1) * dAkkrj
              If(IScal.eq.2)then
                EllJ = EllJ + Q(K,2) * bmf(K,2) * dAkkrj
                SclFcA = SclFcA +
     $            (Alpha(K,1)*Q(K,3)-Alpha(K,2)*Q(K,4))*dAkkrj
                endif
              endif
C           endIf to If(K.eq.JDid)
            endIf
c
          EllJ = -EllJ
   30     Ell(JDid)=Ell(JDid)+EllJ
        If(ISlMap(JDid).ne.0) then
          Ell(JDid) = Ell(JDid) + dbdr(J)*bmf(JJDid,1)
          Ell(JDid) = Ell(JDid) + qtot(JJDid)*dphedr(J)
          endif
        dAjj2 = dwdr(J)/(Wgt(JJDid)**2)
        Ell(JDid) = Ell(JDid) + (Const1+Ffact(JJDid))*dAjj2*q(JJDid,1)*
     $    bmf(JJDid,1)
c
        If(IScal.eq.2) then
          If(ISlMap(JDid).ne.0) then
            Ell(JDid) = Ell(JDid) + dbndr(J)*bmf(JJDid,2)
            Ell(JDid) = Ell(JDid) + qtot(JJDid)*dphndr(J)
            Ell(JDid) = Ell(JDid) + (Const1+Ffact(JJDid))*dAjj2*
     $        Q(JJDid,2)*bmf(JJDid,2)
            SclFoc = SclFoc+Alpha(JJDid,2)*dbndr(J)
            SclFoc = SclFoc-Alpha(JJDid,1)*dbdr(J)
            SclFoc = SclFoc-dphedr(J)*Beta(JJDid)
            SclFoc = SclFoc+Q(JJDid,3)*GFact(JJDid)*Beta(JJDid)*
     $        dphndr(J)/Q(JJDid,4)
            endif
          dAjj2c = -(Constc+Ffact(JJDid))*dAjj2
          SclFcA = SclFcA +
     $      (Alpha(JJDid,1)*Q(JJDid,3)
     $         -Alpha(JJDid,2)*Q(JJDid,4))*dAjj2c
C
          ScalFc = Sclfoc + SclFcA + SclFcC
          Ell(JDid) = Ell(JDid)+ScalFc
          endif
C       NB. Yes, use RHat(JDid).
        Ell(JDid) = Ell(JDid)/(grho(J,1)*rhatx(JDid)+
     $    grho(J,2)*rhaty(JDid)+grho(j,3)*rhatz(JDid))
   20   Continue
C
C     Compute L . (the rest of) dr/dPmn, add to Fock matrix
C
      Do 80 I = 1, NBasis
   80   ChiMax(I) = ArrMax(Chi(1,I),NDo,.True.)
      ADerMx = ArrMax(Ell(IDid+1),NDo,.True.)
C$OMP Parallel Do Schedule(Dynamic) Default(Shared)
C$OMP+ Private(I,J,K,Ind1,Indx)
      Do 100 I = 1, NBasis
        Ind1 = (I*(I-1))/2
        Do 100 J = 1, I
          If(Abs(ChiMax(I)*ChiMax(J)*ADermx).gt.Cutoff) then
            Indx = Ind1 + J
            Do 90 K = 1, NDo
   90         SolFA(Indx) = SolFA(Indx) - Chi(K,I)*Chi(K,J)*Ell(K+IDid)
            endIf
  100     Continue
      Return
      End
*Deck TI
      Function TI(I,MMJ,NMK,M,N,J1,J2,RLSq,RLTmp,L)
      Implicit Real*8(A-H,O-Z)
C
C     Compute the TI(I,J,K) function of Goddard ...
C
      Dimension RLSq(*), RLTmp(*), L(*)
      Save Zero
      Data Zero/0.0d0/
C
      Lim = Min(MMJ,NMK,N-I)
      TI = Zero
      Do 10 IR = 0, Lim
        FN = FacLn(N)
        FR = FacLn(IR+I)
        FNR = FacLn(N-IR-I)
        CF = Exp(FR+FNR-FN)
        X = CombRM(IR,M,J1,J2,RLSq,RLTmp,L)
        TI = TI + CF*X
c     write(6,1234) i,mmj,nmk,j1,j2,cf,x,cf*x,ti
c1234 format('ti:  i=',i2,' mmj=',i2,' nmk=',i2,' j12=',2i2,' cf=',
c    $1pd13.6,' x=',1pd13.6,' t=',1pd13.6,' ti=',1pd13.6)
   10   Continue
      Return
      End
*Deck TIJ
      Function TIJ(I,J,M,N,J1,J2,RLSq,RLTmp,L)
      Implicit Real*8(A-H,O-Z)
C
      Dimension RLSq(*), RLTmp(*), L(*)
C
      TIJ = TI(I,M-J,M-J,M,N,J1,J2,RLSq,RLTmp,L)
      Return
      End
*Deck TstVAC
      Subroutine TstVAC(IOut,IPrint,UsePBC,NAtoms,IAn,IAtTyp,AtmChg,
     $  ItSCRF,ISmStp,DoSim,NShell,ShellN,ShellA,EXX,DirSCF,ISCRF0,DFT,
     $  Mixed,IDenFit,IRadAn,IDIIS,VarAcc,Sleazy,IVShft,FkGues,NExtDmp,
     $  HaveTM,MaxDamp,DampLim,LvRise,SpDiag,IHMeth,MinDamp,DoDMS)
      Implicit Real*8(A-H,O-Z)
C
C     Decide whether to turn off Sleazy integrals or variable accuracy.
C     Returns VarAcc = 1 ... No cheap pass.
C                      4 ... Pass with cheap integrals and XC grid.
C                      7 ... Pass with cheap XC grid.
C                      8 ... Pass with cheap integrals.
C
      Logical UsePBC, Sleazy, FkGues, DFT, HaveTM, HavDif, NoVAC2,
     $  NoVACX, SpDiag, SEHam, Mixed, Ck4TM, DoSim, DirSCF, DoDMS
      Integer ShellN(*), ShellA(*), VarAcc
      Dimension IAn(NAtoms), IAtTyp(NAtoms), AtmChg(NAtoms), EXX(*)
      Parameter (IVSPBC=-25,IVSDef=-100,IVSTM=-100,NExDef=0,NExDTM=1)
      Save Two, DifMax, DLEDIIS, DLCDIIS
      Data Two/2.0d0/, DifMax/0.05d0/, DLEDIIS/0.25d0/, DLCDIIS/0.125d0/
C
C     Whether to do VarAcc.  Turned off for transition metals, Sleazy,
C     after a small geometry step and during SimOpt.  Conversely, level
C     shifting is turned on by default for TMs.
C
      NoVACX = ISmStp.eq.1.or.Sleazy.or.ItSCRF.ge.4.or.ISCRF0.eq.2.or.
     $  UsePBC.or.(DFT.and.(IRadAn.eq.1.or.IRadAn.eq.2.or.IRadAn.eq.3))
     $  .or.DoSim.or.(SEHam(IHMeth).and..not.SpDiag).or..not.DirSCF
     $  .or.DoDMS
      NoVAC2 = NoVACX.or.(DFT.and..not.Mixed.and.IDenFit.gt.1)
      HaveTM = Ck4TM(1,NAtoms,IAn,IAtTyp)
      Call DecPrn(IOut,IPrint,IRadAn,0,IPrune,NTheta,NPhi,MaxRad,MaxAng,
     $  ECTolR)
      NoVACX = NoVACX.or.(DFT.and.HaveTM.and.MaxAng.lt.590)
C
C     Adjust some parameters if energy-DIIS is on, then apply
C     TM vs non-TM defaults.
C
      If((Mod(IDIIS,100)/10).gt.1) then
        If(IVShft.eq.0) IVShft = -1
        DampLim = DLEDIIS
        If(MaxDamp.eq.0) MaxDamp = 1
      else
        DampLim = DLCDIIS
        If(MaxDamp.eq.0) MaxDamp = 10
        endIf
      If(LvRise.eq.0) LVRise = 2
      If(HaveTM.and.IVShft.eq.0) IVShft = IVSTM
      If(MinDamp.ne.0) then
        NExtDmp = MinDamp
      else if(HaveTM) then
        NExtDmp = NExDTM
      else
        NExtDmp = NExDef
        endIf
      MaxDamp = Max(MaxDamp,NExtDmp)
      If(IVShft.eq.0) then
        If(UsePBC) then
          IVShft = IVSPBC
        else
          IVShft = IVSDef
          endIf
      else if(IVShft.eq.-1) then
        IVShft = 0
        endIf
C
C     Turn off Sleazy for all-electron heavy atoms.
C
      Do 20 I = 1, NAtoms
   20   Sleazy = Sleazy.and.
     $    .not.(IAn(I).gt.18.and.Abs(AtmChg(I)-IAn(I)).lt.Two)
C
C     Turn off both variable accuracy and sleazy if diffuse functions are
C     present.
C
      HavDif = .False.
      Do 30 ISh = 1, NShell
        Ind = ShellA(ISh) - 1
        Do 30 I = 1, ShellN(ISh)
   30     HavDif = HavDif.or.EXX(Ind+I).lt.DifMax
c
      NoVAC2 = NoVAC2.or.HavDif
      NoVACX = NoVACX.or.(DFT.and.NoVAC2.and.(Mixed.or.IDenFit.le.1))
      Sleazy = Sleazy.and..not.HavDif
      If(FkGues) then
        If(VarAcc.ne.6) NExtDmp = Max(NExtDmp,2)
        VarAcc = 4
      else if(VarAcc.eq.5) then
        If(NoVACX.and.NoVAC2) then
          VarAcc = 1
        else if(NoVAC2) then
          VarAcc = 7
        else if(NoVACX) then
          VarAcc = 8
        else
          VarAcc = 4
          endIf
      else if(VarAcc.eq.6) then
        VarAcc = 1
        endIf
      Return
      End
*Deck UHFPrt
      Subroutine UHFPrt(IOut,IPrint,SpDiag,IOpCl,IRwEig,IRwCA,IRwCB,
     $  IRwPA,IRwPB,IRwPT,IRwPS,IRwUR,NOATot,NOBTot,NBasis,NBsUse,
     $  VShftA,VShftB,Eig,A)
      Implicit Real*8(A-H,O-Z)
C
C     Print final matrices.
C
      Logical SpDiag, CBfn, Spinor
      Dimension Eig(NBsUse,*), A(*)
 1000 Format(' No print of final matrices with Sparse.')
 1040 Format(' <<<<<<<<<<Alpha-Beta overlap>>>>>>>>>>')
C
      Call DecoSC(IOpCl,NSpBlk,NRI,NDimBl,CBfn,NSpBlX,Spinor)
      If(SpDiag) then
        If(IPrint.ge.2) Write(IOut,1000)
      else
        NOrbs = NBsUse*NDimBl
        Call FileIO(2,-IRwEig,NOrbs,Eig,0)
        NEX = NOATot
        If(NDimBl.eq.2) NEX = NOATot + NOBTot
        Do 10 I = (NEX+1), NOrbs
   10     Eig(I,1) = Eig(I,1) - VShftA
        IndB = NDimBl + 1
        If(NSpBlk.eq.1) then
          Call AMove(NOrbs,Eig,Eig(1,IndB))
        else
          Call FileIO(2,IRwEig,NOrbs,Eig(1,IndB),0)
          Do 20 I = (NOBTot+1), NOrbs
   20       Eig(I,IndB) = Eig(I,IndB) - VShftB
          endIf
        Call FileIO(1,-IRwEig,2*NOrbs,Eig,0)
        If(IPrint.ge.2) then
          NBD = NBasis*NDimBl
          NBDR = NRI*NBD
          NBNU = NBDR*NBsUse*NDimBl
          Call FileIO(2,-IRwCA,NBNU,A,0)
          Call MatPtX(IOut,IOpCl,1,NBasis,NOrbs,NEX,A,Eig)
          If(NSpBlk.eq.2) then
            Call FileIO(2,-IRwCB,NBNU,A,0)
            Call MatPtX(IOut,IOpCl,2,NBasis,NOrbs,NOBTot,A,Eig(1,IndB))
            endIf
          If(IPrint.ge.3) then
            NTTD = NRI*(NBD*(NBD+1))/2
            Call FileIO(2,-IRwPA,NTTD,A,0)
            If(NSpBlk.eq.2) Call FileIO(2,-IRwPB,NTTD,A(1+NTTD),0)
            Call LTOutX(IOut,0,IOpCl,1,0,'Final density matrix',NBasis,
     $        NBasis,A,A(1+NTTD),0)
            Call FileIO(2,-IRwPT,NTTD,A,0)
            IOpClX = IOpCl - Mod(IOpCl,2)
            Call LTOutX(IOut,0,IOpClX,1,0,'Total density matrix',NBasis,
     $        NBasis,A,A,0)
            If(NSpBlk.eq.2) then
              Call FileIO(2,-IRwPS,NTTD,A,0)
              Call LTOutX(IOut,0,IOpClX,1,0,'Spin density matrix',
     $          NBasis,NBasis,A,A,0)
              If(NRI.eq.1) then
                Write(IOut,1040)
                Call FileIO(2,-IRwUR,NBsUse**2,A,0)
                Call OutMat(IOut,0,A,NBsUse,NBsUse,NBsUse,NBsUse)
                endIf
              endIf
            endIf
          endIf
        endIf
      Return
      End
*Deck ZUpdt2
      Subroutine ZUpdt2(CC,Ncore,Norbs,Thresh,IEnd,IOut,pn,pCrit,jiz,iz,
     $  NZz,jiz2,iz2,NZz2,jidz,idz,NZdz,X,Fixed)
      Implicit Real*8(A-H,O-Z)
      Dimension CC(Ncore),X(Norbs)
      Logical Fixed,pCrit
      Save Half,One,Zero,Two,Tenth
      Data Half/0.5d0/,One/1.d0/,Zero/0.d0/,Two/2.d0/,Tenth/0.1d0/
C
C     Update z matrix
C
      trup = ArrSum(Norbs,CC(idz))
      trdw = ArrSum(Norbs,CC(iz2))
      pn = trup/trdw
C
C     Check validity of pn
C
      pCrit = pn.le.-One.or.pn.ge.One
      If(pCrit) then
        If(Abs(trup).gt.Tenth) Return
C       Use NSI form, if pn is unphysical
        f1 = Half
        f2 = Zero
      else
        If(pn.le.Zero) then
          f1 = One/Two/(One-pn)
        else
          f1 = One/Two/(One+pn)
        endIf
        f2 = - pn*f1
        endIf
C
C     zn = z + f1*dz + f2*(1-z2)
C
      If(Fixed) then
        Do 10 i = 0, (Norbs+NZz-1)
   10     CC(iz+i) = CC(iz+i) + f1*CC(idz+i) + f2*CC(iz2+i)
      else
        If(Abs(f2).gt.Thresh) then
          Istart = IEnd
          itype = 2
          Call SpAddM(CC,Ncore,CC(jidz),CC(idz),CC(jiz2),CC(iz2),x,
     $      Istart,jidzn,idzn,IEnd,Norbs,NZdzn,f1,f2,Thresh,itype)
          Istart = IEnd
          itype = 2
          Call SpAddM(CC,Ncore,CC(jiz),CC(iz),CC(jidzn),CC(idzn),x,
     $      Istart,jizn,izn,IEnd,Norbs,NZzn,One,One,Thresh,itype)
        else
          Istart = IEnd
          itype = 2
          Call SpAddM(CC,Ncore,CC(jiz),CC(iz),CC(jidz),CC(idz),x,Istart,
     $      jizn,izn,IEnd,Norbs,NZzn,One,f1,Thresh,itype)
          endIf
        jiz = jizn
        iz = izn
        NZz = NZzn
        endIf
      Return
      End
