*Deck RdGeom
      Subroutine RdGeom(V,MDV)
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
C1OV1
C
C     Overlay 1 ... Gaussian system input and optimization control.
C
C      In the gaussian system, overlay one contains those programs
C  that read in geometry and optimization input and those that control
C  optimization and numerical frequency calculations.  Currently, the
C  following links are implemented:
C
C      Link 101 ... Basic input.  This link reads in the title and
C                   atomic specifications, either in the form of a
C                   symbolic Z-Matrix with associated variables
C                   and constants or in the form of model builder
C                   connectivity an orientation input.  The latter
C                   is still in the experimental stage.
C
C      Link 102 ... Fletcher-Powell optimization program.  This link
C                   implements the algorithm of Fletcher and Powell
C                   as modified by Binkley and Pople.  It optimizes
C                   the geometry using only the energy, by numerically
C                   differentiating when necessary.
C
C      Link 103 ... Berny Schlegel gradient optimization program.  This
C                   program is used in conjunction with those links that
C                   produce analytical energy derivatives to perform
C                   geometry optimizations.  This link can search for
C                   local minima or saddle points (transition states)
C
C      Link 105 ... Murtaugh-Sargent optimization program.  This program
C                   also uses analytically determined first derivatives
C                   to optimize the geometry.  Although in general it is
C                   not as efficient as link 103, it has a more stable
C                   algorithm and may perform better when given a poor
C                   initial geometry or unusual structure.  It can only
C                   locate local minima.
C
C      Link 106 ... Controls numerical differentiation of analytic first
C                   derivatives to produce force constants and dipole
C                   derivatives and differentiation of analytic second
C                   derivatives to produce polarizability derivatives,
C                   dipole second derivatives, and nuclear coordinate
C                   third derivatives.
C
C      Link 107 ... Performs Linear Synchronous Transit searches for
C                   approximate transition states.
C
C      Link 108 ... Performs potential surface scans.
C
C      Link 109 ... Optimizes to stationary points, including transition
C                   states, using only energies.  Uses basic Newton-
C                   Raphson with a read-in Hessian which is not updated.
C
C      Link 110 ... Numerically differentiates energies onece or twice to
C                   produce forces and force constants.
C
C      Link 111 ... Numerically differentiaties energies and dipole
C                   moments twice to produce polarizabilities and
C                   hyperpolarizabilities.
C
C      Link 112 ... SCVS.
C
C      Link 113 ... Eigenvalue-following optimization using analytic
C                   derivatives.
C
C      Link 114 ... Eigenvalue-following optimization using numerical
C                   derivatives.
C
C      Link 115 ... Reaction path following using the method of
C                   Gonzalez and Schlegel.
C
C2IOP(5)
C
C   IOP(5)  L103   MODE OF OPTIMIZATION
C             0    FIND LOCAL MINIMUM
C             1    FIND A SADDLE POINT
C             N    FIND A STATIONARY POINT ON THE ENERGY SURFACE WITH
C                  N NEGATIVE EIGENVALUES OF THE 2ND DERIVATIVE MATRIX
C
C           L107   MODE OF SEARCH
C             0    LOCATE THE MAXIMUM IN THE LST PATH.
C             1    SCAN THE LST PATH.
C
C2IOP(6)
C
C   IOP(6)  L102, L103, L105, L107   MAXIMUM NUMBER OF STEPS (OR NUMBER
C           L109, L113, L114           OF STEPS FOR AN LST SCAN).
C             0    NSTEP = Max(20,NVAR+10) (L103)
C                        = Min(20,NVAR+10) (L102, L105, L109)
C                        = Min(40,NVar+20) (L113, L114)
C             N    NSTEP = N
C2IOP(7)
C
C   IOP(7)  L103, L105, L109, L112, L113, L114   CONVERGENCE ON THE
C           FIRST DERIVATIVE AND ESTIMATED DISPLACEMENT FOR THE
C           OPTIMIZATION RMS FIRST DERIVATIVE .LT. CONFV,
C           RMS EST. DISPLACEMENT .LT. CONVX=4*CONVF.  For L123,
C           similar convergence control for the GS integrator.
C            -1    ConvF = 1/600 HARTRE/BOHR OR RADIAN
C             0    CONVF = 0.0003 HARTRE/BOHR OR RADIAN
C             N    CONVF = N*10**-6
C
C     L116, L117:  Convergence on electric field/charges
C            -1    Default value for optimizations  -- 10**-7.
C             0    Default value for single-points:
C                  10**-5 in L116, 10**-7 in L117.
C             N    10**-N.
C
C     L123: Regular force/displacement convergence for GS2.  For
C           EulerPC, N/1000000 displacement convergence.
C
C2IOP(8)
C
C   IOP(8)  L103, L109, L112   MAXIMUM STEP SIZE ALLOWED DURING OPT.
C             0    DXMAXT = 0.1 BOHR OR RADIAN (L103, Estm or UnitFC).
C                         = 0.3 Bohr or Radian (L103, Read or CalcFC).
C                         = 0.2 Bohr or Radian (L105).
C                         = 0.3 Bohr or Radian (L113, L114).
C             N    DXMAXT = 0.01 * N
C
C     IOp(8) ... L117:  General control.
C          0 ... Which type of basin to use to partition the density
C                isosurface.  Default is 4
C          1 ... GradVne
C          2 ... GradRho
C          3 ... Don't Use Basins, Use only the Center of Nuclear
C                Charge
C          4 ... Use Interlocking Spheres
C         N0 ... Order of Adam's-Bashforth-Moulton (ABM) predictor-
C                corrector method to use in solving diff. eqns. for the
C                grad RHO or Vne trajectories.  Default is 4, max is 9.
C        N00 ... Number of small steps per ABM step to be used in
C                starting ABM and when "slow down" is needed in ABM.
C                Default is 5.
C       N000 ... Which approximation to make.  Default is III for Tomasi
C                (interlocking spheres) and IV for general surface.
C       1000 ... Apprx. I - Don't Do Self-Polarization or "Compensation"
C       2000 ... Apprx. II - Do-Self Polarization, But No Compensation.
C       3000 ... Apprx. III - Do Self-Polarization and Compensation.
C       4000 ... Apprx. IV - Do III and Allow Surface To "Relax" in
C                            Solution if no spheres ...
C      N0000 ... Whether to evaluate densities using orbitals or density
C                matrix.  Default is to use density.
C      10000 ... Use MOs.
C      20000 ... Use density.
C
C     IOp(8) ... L121: Time step, N*0.0001 fs, default 0.1
C
C2IOP(9)
C
C   IOp(9)  L103   Use of Trust radius.
C             0    Whether to update trust radius (DXMaxT, default Yes).
C                  Default is Yes for minima, no for TS.
C             1    No.
C             2    Yes.
C            00    Whether to scale or search the sphere when reducing
C                  the step size to the trust radius (Default search for
C                  minima, scale for transition states.).
C            10    Scale.
C            20    Search.
C
C   IOP(9)  L106   Whether to use symmetry to reduce the number of
C           L110   displacements:
C             0    Yes.
C             1    No.
C            10    Do not use symmetry to skip steps back.
C           100    Do not use symmetry to skip equivalent atoms.
C
C   IOP(9)  L107   Whether to maintain symmetry along the search path,
C             0    Yes.
C             1    No.
C
C   IOp(9)  L117   Whether to delete points which are too close together:
C           0 ... No
C           1 ... Yes, using a default criteria (0.05 Angstroms)
C          -N ... Yes, using a (10**-N Angstroms) criteria
C                 How close to get to the isosurface in search.
C           0 ... Approx 1.0D-6 (N=20)
C           N ... 2.0**-N
C
C   IOp(9)  L121   Whether to read in initial velocities:
C               0  Default (same as 1)
C               1  Generate random initial velocity
C               2  Read in initial cartesian velocity (Bohr/sec)
C               3  Read in initial MW cartesian velocity (sqrt(amu)*Bohr/sec)
C
C2IOP(10)
C
C   IOp(10) L103, L105, L109, L112, L113, L114  Input of initial Hessian:
C                  All values must be in atomic units (Hartre, Bohr,
C                  and radians).
C             0    Use defaults (not valid for L109).
C             1    Read ((FC(I,J),J=1,I),I=1,NVAR) (8F10.6) (L103 only).
C             2    Read I,J,FC(I,J)      (5I3,F20.0) (L103 only).
C                  End with a blank card.
C             3    Read from checkpoint file in internal coordinates.
C             4    Second derivative matrix calculated analytically.
C                  (not valid for L109).
C             5    Read cartesian forces and force constants from the
C                  checkpoint file are convert to internal coordinates.
C             6    Read cartesian forces followed by cartesian force
C                  constants (both in format 6F12.8) from input stream.
C                  followed by a blank line.
C             7    Use semiempirical force constants.
C             8    Use unit matrix (default for L105; only recognized
C                  by 103).
C             9    Estimate force constants using valence force field.
C            10    Use unit matrix throughout.
C
C2IOP(11)
C
C   IOP(11)  L103  TEST OF CURVATURE.  BOMB THE JOB IF THE SECOND
C                  SECOND DERIVATIVE MATRIX HAS THE WRONG NUMBER OF
C                  NEGATIVE EIGENVALUES.
C               0  ... DEFAULT (TEST for z-matrix or cartesian TS but not
C                      for LST/QST or for minimum).
C               1  ... DON'T TEST.
C               2  ... TEST.
C
C   IOp(11) L117 ... Scaling Factor for Determining Overlaps of VDW atoms
C             -1 ... Turn off scaling
C              0 ... Default is 1.010
C              N ... 1.000 + N*(0.001)
C
C                    Step size for ABM method in Trudge for isodensity
C                    method.
C              0 ... 0.05 (N=2)
C              N ... 0.1/N
C
C2IOP(12)
C
C   IOP(12)  L103  OPTIMIZATION CONTROL PARAMETERS
C             0    USE DEFAULT VALUES
C             1    READ IN NEW VALUES FOR ALL PARAMETERS (SEE INITBS)
C
C2IOp(13)
C
C   IOp(13) L103,L113,L114,L115,L123 Type of Hessian Update:
C             0 ... Default (9 for L103 minimization, 7 for L103 TS,
C                   D2Corr and L115, Powell for L113 and L114, Bofill
C                   in L123).
C             1 ... Powell (not in L103).
C             2 ... BFGS (not in L103)
C             3 ... BFGS, safeguarding positive definateness (not in
C                   L103 or L115)
C             4 ... D2Corr (New, only in L103 and L115).
C             5 ... D2Corr (Old, only in L103 and L115).
C             6 ... D2Corr (BFGS)
C             7 ... D2Corr (Bofill Powell+MS for transition states).
C             8 ... D2Corr (No update, use initial Hessian).
C             9 ... D2Corr (New if energy rises, otherwise BFGS).
C
C     IOP(13) L121 Multi-time step parameter (NDtrC,NDtrP)
C               0 ... No multi-time stepping
C              NN ... Iterate density constraints NN times per step
C            MM00 ... Do gradient once every MM steps
C
C     IOp(13) L123 Hessian updat:
C             0 ... Default (Bofill).
C             1 ... Murtagh-Sargent (SR1) update
C             2 ... Powell-symmetric-Broyden (PSB) update
C             3 ... Bofill's update
C             4 ... Sqrt(Bofill) update
C             5 ... No update (keep old Hessian)
C
C2IOP(14)
C
C   IOP(14)  L103  Maximum number of bad steps to allow before doing
C                  a linear minimization (i.e., no quadratic step):
C                  will be attempted.
C            -1    0.
C             0    Default (0 for TS, 1 for minima).
C             N    Allow N -- linear only starts with the N+1st.
C
C2IOP(15)
C
C   IOP(15) L103,L109  ABORT IF DERIVATIVES TOO LARGE
C           -1 or 0    No force test at all.
C              N       FMAXT = 0.1 * N
C
C2IOP(16)
C
C IOP(16) L103,L113,L114  MAXIMUM ALLOWABLE MAGNITUDE OF THE EIGENVALUES
C                  OF THE SECOND DERIVATIVE MATRIX. IF THE LIMIT IS
C                  EXCEEDED, THE SIZE OF THE EIGENVALUE IS REDUCED
C                  TO THE MAXIMUM, AND PROCESSING CONTINUES.
C             0    EIGMAX = 25.0 HARTRE / BOHR**2 OR RADIAN**2
C             N    EIGMAX = 0.1 * N
C2IOP(17)
C
C IOP(17) L103,L113,L114  MINIMUM ALLOWABLE MAGNITUDE OF THE EIGENVALUES
C                  OF THE SECOND DERIVATIVE MATRIX. SIMILAR TO IOP(16)
C             0    EIGMIN = 0.0001
C             N    EIGMIN = 1. / N
C
C2IOP(18)
C
C  IOP(18) L103 ... Coordinate system.
C             0 ... Proceed normally
C             1 ... Second derivatives will be computed as directed on
C                   the variable definition cards.  No optimization will
C                   occur.
C            10 ... Do optimization in cartesian coordinates.
C            20 ... Do full optimization in redundant internal coord.
C            30 ... Do full optimization in pruned distance matrix coords.
C            40 ... Do optimization in Z-matrix coordinates.
C            50 ... Do full optimization in redundant internal coords with
C                   large molecular tools.
C           100 ... Read the AddRedundant input section for each structure.
C          1000 ... Do not define H-bonds (default)
C          2000 ... Define H-bonds with no related coordinates
C          3000 ... Define H-bonds and related coordinates
C         10000 ... Reduce the number of redundant internals
C         20000 ... Define all redundant internals
C        100000 ... Old definition of redundant internals.
C       0000000 ... Default (2000000).
C       1000000 ... Skip MM atoms in internal coordinate definitions and
C                   do microiterations the old way, in L103.
C       2000000 ... Include MM atoms in internal coordinate definitions
C                   (no microiterations).
C       3000000 ... Skip MM atoms in internal coordinate definitions and
C                   do microiterations the new way, in L120.
C       4000000 ... Microiterations for pure MM, done in L402.
C
C2IOP(19)
C
C   IOP(19)  L103  SEARCH SELECTION
C             0    Default (same as 6).
C             2    LINEAR AND STEEPEST DESCENT.
C             3    STEEPEST DESCENT AND LINEAR ONLY WHEN ESSENTIAL.
C             4    Quadratic if curvature is correct; RFO if not.
C                  Linear as usual.
C             5    Quadratic if curvature is correct; RFO if not.
C                  No linear search.
C             6    RFO and linear.
C             7    RFO without linear.
C             8    Newton-Raphson and linear.
C             9    Newton-Raphson only.
C            10    GDIIS and linear
C            11    GDIIS only.
C            13    First-order simultaneous optimization.
C            15    GEDIIS.
C            17    FOSimult with fixed geometry, for debugging.
C            19    FOSimult with full SCF convergence, for debugging.
C
C        L113,L114 Search Selection:
C             0    P-RFO OR RFO STEP ONLY               (DEFAULT)
C             1    P-RFO OR RFO STEP FOR "WRONG" HESSIAN
C                  OTHERWISE NEWTON-RAPHSON
C
C2IOP(20)
C
C  IOP(20) L101, L106, L108, L109, L110 ... INPUT UNITS
C      0 ... ANGSTROMS DEGREES
C      1 ... BOHRS     DEGREES
C      2 ... ANGSTROMS RADIANS
C      3 ... BOHRS     RADIANS
C
C2IOP(21)
C
C     IOP(21) L103,L113,L114 ... EXPERT SWITCH.
C     0 ... NORMAL MODE.
C     1 ... EXPERT MODE ... CERTAIN CUTOFFS USED TO CONTROL THE
C           OPTIMIZATION WILL BE RELAXED.  THESE INCLUDE FMAXT,
C           DXMAXT, EIGMAX AND EIGMIN.
C
C2IOP(22)
C
C  IOP(22) L107 ... whether to reorder coordinates for
C                   maximum coincidence.
C             0 ... yes.
C             1 ... assume reactant order equals product order.
C             2 ... read in a re-ordering vector from the input
C
C  IOp(22) L115/L123 KIND OF SEARCH:
C              0 ... BOTH DIRECTIONS AND GENERATE SEARCH VECTOR
C              1 ... FORWARD DIRECTION AND GENERATE S. VECTOR
C              2 ... BACKWARD DIRECTION AND GENERATE S. VECTOR
C              3 ... BOTH DIRECTIONS AND GENERATE S. VECTOR
C              4 ... FORWARD DIRECTION AND READ S. VECTOR 8F10.6
C              5 ... FORWARD DIRECTION AND READ S. VECTOR 8F10.6
C              6 ... BACKWARD DIRECTION AND READ S. VECTOR 8F10.6
C              7 ... BOTH DIRECTIONS AND READ S. VECTOR 8F10.6
C
C2IOp(23)
C      IOp(23) No longer used.
C
C2IOp(24)
C      IOp(24) ... Whether to round tetrahedral angles.
C            0 ... Default (Yes).
C            1 ... Yes, round angles within 0.001 degree.
C            2 ... No.
C
C2IOp(25)
C
C     IOp(25) ... Wether SCRF is used with numerical polarizability:
C           0 ... No.
C           1 ... Yes, the field in /Gen/ must be cleared each time.
C
C2IOp(26)
C
C     IOp(26) ... Accuracy of function being optimized:
C       -NNMM ... Energy 10**-(NN), Gradient 10**-(MM).
C          -1 ... Read in values
C           0 ... Default (same as 1).
C           1 ... Normal accuracy for HF (energy and gradient both 1.d-7).
C           2 ... Accuracy for DFT with SG1 grid (Energy 1.d-5,
C                 gradient 1.d-4).
C           3 ... Fine grid accuracy for DFT (Energy 1.d-7, gradient 1.d-6)
C           4 ... Ultrafine accuracy (E 1.d-7, G 1.d-6).
C           5 ... Superfine accuracy (E 1.d-7, G 1.d-7).
C
C2IOp(27)
C     IOp(27) = IJKL (i.e. 1000*I+100*J+10*K+L)
C                 Transition state searching using QST
C                 and redundant internal coordinates
C     L= 0,1 ...  Input one structure, either initial guess
C                 of the minimizing structure or transition
C                 structure without QST.
C     L= 2 ...    Input 2 structures, the first one is the
C                 reactant, the second one is the product.
C                 The union of the two redundant coordinates
C                 are taken as the redundant coords for the TS.
C                 The values of the TS coord are estimated by
C                 interpolating the sturcture of R and P.
c                 R and P are used to guide the QST optimization
c                 of the TS.
C     L= 3 ...    Input 3 structures.  The first one is reactant
C                 the second one is the product.  The third one
C                 is the initial guess of the transition structure.
c                 R and P are used to guide the QST optimization
c                 of the TS.
C    K = 1-9      Interpolation of initial guess of TS between
C                 R and P (TS=0.1*J*R + 0.1*(10-J)*P, default J=5)
C    J = 1        LST constraint in internals
C    J = 2        QST constraint in internals
C    J = 3        LST constraint in distance matrix space
C    J = 4        QST constraint in distance matrix space
C    I = 0-9      Control parameters for climbing phase of QST
C                 (e.g. QSTRad = 0.01*I, default QSTrad = 0.05)
C2IOp(28)
C     IOp(28) L103 ... Number of translations and rotations to
C                      remove during redundant coordinate transformations:
C               -2 ... 0.
C               -1 ... Normal (6 or 5 for linear molecules).
C                0 ... Default, same as -1.
C                N ... N.
C
C2IOP(29)
C
C     IOP(29) L101 ... SPECIFICATION OF NUCLEAR CENTERS
C        0 ... BY Z-MATRIX
C        1 ... BY DIRECT COORDINATE INPUT (must set iop(29) in L202).
C        2 ... GET Z-MATRIX AND VARIABLES FROM THE CHECKPOINT FILE.
C        3 ... GET CARTESIAN COORDINATES ONLY FROM THE CHECKPOINT FILE.
C        4 ... By model builder, model A.
C        5 ... By model builder, model B.
C        6 ... Get Z-matrix from the checkpoint file, but read new
C              values for some variables from the input stream.
C        7 ... Get all input (title, charge and multiplicity, structure)
C              from the checkpoint file.
C       10 ... Print details of the model building process.
C      000 ... Default (same as 100).
C      100 ... Do not abort job if model builder generates a z-matrix
C              with too many variables.
C      200 ... Abort job if model builder generates a z-matrix with too
C              many variables.
C     1000 ... Read optimization flags in format 50L1 after the z-matrix.
C     2000 ... Set all optimization flags to optimize.
C     3000 ... Purge flags except the frozen variables.
C     4000 ... Rebuild the coordinate system.
C     5000 ... (2+3) Purge all flags but keep the coordinate definition.
C     6000 ... Generate new redundant coordinates, reading an input
C              section selecting frozen and optimized atoms.
C     7000 ... Mark all internal coordinates as frozen before handling
C              modred input.
C    00000 ... Default, same as 10000.
C    10000 ... Mark Z-matrix constants as frozen variables rather than
C              wired-in constants.
C    20000 ... Do not retain symbolic constants.
C   100000 ... Generate a symbolic z-matrix using all Cartesians if none
C              is present on the checkpoint file (a hack to make IRCs
C              work with Cartesian input).
C   200000 ... Same as one, but retain the redundant internal coordinate
C              definitions.
C  1000000 ... Get input type or chk file name to read from input stream;
C              title and charge/multiplicity for each structure read from
C              input.
C  2000000 ... Read input type for each structure from input stream; title
C              and charge/multiplicity are those of last chk file read.
C  9000000 ... Same as 0000000.
C 00000000 ... Default (read one set of charge/multiplicity pairs unless
C              both NFrag and ONIOM are set).
C 10000000 ... Read ONIOM charge/multiplicity pairs if reading any.
C              Fragment values will be defaulted from the supermolecule.
C 20000000 ... Read fragment charge/multiplicity pairs if reading any.
C              ONIOM model system values will be defaulted from the real
C              system.
C 30000000 ... Read two lines, with ONIOM followed by fragment values.
C
C2IOP(30)
C
C     IOP(30)  L103  ARE THE READ-WRITE FILES TO BE UPDATED?
C                  THIS OPTION IS SET FOR THE LAST CALL TO 103 IN
C                  FREQUENCY CALCULATIONS IN ORDER TO PRESERVE THE
C                  VALUES OF THE VARIABLES FOR ARCHIVING.  It also
C                  suppresses error termination on large gradients.
C             0    Yes.
C             1    No, print internal coordinate step but don't set up
C                  for microiterations and don't update the rwf.
C             2    Set up for step but don't update coordinates; for
C                  QM/MM iterative frequencies.
C
C2IOP(32)
C
C     IOP(32)  TITLE CARD PUNCH CONTROL.
C            0     DON'T PUNCH.
C            1     PUNCH.
C
C2IOP(33)
C
C     IOP(33)  L101 L102 L103 L106 L109 L110 L113 L114  ... DEBUG PRINT
C       0 .. OFF
C       1 .. ON
C
C2IOP(34)
C
C     IOP(34)  L101 L102 L103  ... DEBUG + DUMP PRINT
C       0 .. OFF
C       1 .. ON
C
C2IOP(35)
C
C     IOP(35) ... RESTART (L102-L112).
C     0 ... NORMAL OPTIMIZATION.
C     1 ... FIRST POINT OF A RESTART.  GET GEOMETRY, WAVEFUNCTION, ET.
C           FROM THE CHECKPOINT FILE.
C
C2IOP(36)
C
C     IOP(36) ... CHECKPOINT.
C     0 ... NORMAL CHECKPOINT OF OPTIMIZATION.
C     1 ... SUPPRESS CHECKPOINTING.
C
C2IOP(37)
C
C     IOP(37) ... D2E CLEANUP (obsolete)
C     0 ... NO CLEANUP.
C     1 ... THIS IS THE LAST POINT AT WHICH ANALYTIC SECOND DERIVATIVES
C           WILL BE DONE.  DELETE THE D2E FILE AND THE BUCKETS AND
C           TRUNCATE THE READ/WRITE FILES.
C
C2IOP(38)
C
C     IOP(38) ... Entry control option
C           0 ... Continuation of run.
C           1 ... Initial entry.
C           N>1 . In L103: Initial entry of guided optimization using
C                 N levels.
C          N0 ... In L106: differentiate Nth derivatives once.
C                 In L110 and L111: differentiate energy N times.
C         000 ... In L106: differentiate wrt nuclear coordinates.
C         100 ... In L106: differentiate wrt electric field.
C         200 ... In L106: differentiate wrt field and nuclear.
C        1000 ... In L106: Save original forces and force constants.
C       00000 ... In L106: Assume all quantities available at the central
C                 point will also be computed at displaced points.
C       10000 ... In L106: No analytic nuclear coordinate derivatives will
C                 done at displaced points, even though they were
C                 done at the central point.
C      000000 ... L106 control of number of diff steps.
C      100000 ... Do 2-point differentiation (one step each way).
C      200000 ... Do 4-point differentiation (two steps each way).
C     0000000 ... Default (1)
C     1000000 ... Differentiate wrt translation vectors for PBC for elasticity.
C     2000000 ... Do not differentiate wrt translation vectors.
C
C2IOP(39)
C
C     IOP(39) ... Step size control for numerical differentiation.
C                 (L106, L109, L110, L111).  Path step size in L115.
C           0 ... Use internal default (0.001 Angstroms and 0.001/3 au
C                 E-field in L106, 0.005 A in L109, 0.01 Angstrom in
C                 L110, 0.001 au in L111).
C           N ... Use step-size of 0.0001*N (angstroms in L106, L109, L110,
C                 electric field au in L111).  In L106, the electric field
C                 step will be 3x smaller.
C          -1 ... Read stepsize (up to 2 for L106, 1 for others), free-format.
C        -N>1 ... Use step-size of 0.0001*N atomic units everywhere.
C
C                 L123:  Step-size:
C         N<0 ... Supplied step size is in units of sqrt(amu)*bohr.
C         N>0 ... Supplied step size is in units of bohr.
C         0   ... Default (0.1 bohr, except 0.075 bohr for EulerPC,
C                 or original value if restart.  DVV default is 0.25 fs)
C
C2IOp(40)
C
C     IOp(40) ... L113, L114:  Hessian recalculation.
C          -1 ... Pick up analytic second derivatives every time.
C           0 ... Just update.  The default, execpt for CalcAll.
C           N ... Recalculation the Hessian every N steps.
C
C     IOp(40) ... L116:  Whether to read initial E-field:
C           0 ... Start with 0.0.
C           1 ... Read from checkpoint file.
C           2 ... Read from input stream.
C
C     IOP(40) ... L123:  Truncation error threshold for the modified BS
C                 integration routine (EPS)
C         0   ... default (10**-8 bohr)
C         N   ... 10**-N bohr
C
C2IOp(41)
C
C     IOp(41) ... Take previous geometry from checkpoint file:
C                 N > 0 ... Nth point of geometry optimization
C                           (z-matrix only).  Converted to -(N+1)
C                           if no z-matrix.
C                 N < 0 ... Nth geometry on trajectory file.
C
C2IOp(42)
C
C     IOp(42) ... L103, L115, L123:  Number of points along the reaction
C                 path in each direction.  Default is 10.
C
C     IOp(42) ... L117:  Cutoff to be used in evaluating densities.
C           0 ... 1.0D-10
C           N ... 1.0D-N
C
C2IOp(43)
C
C     IOP(43) ... L116: Extent of Reaction Field.
C           0 ... Dipole
C           1 ... Quadrupole
C           2 ... Octapole
C           3 ... Hexadecapole
C
C     IOP(43) ... L117:  How to define Radii
C           0 ... Default is 11
C           1 ... Use internally stored Radii, centers will be on atoms
C           2 ... Read-in centers and radii on cards
C          10 ... Force Merz-Kollman radii (Default)
C          20 ... Force CHELP (Francl) recommended radii.
C          30 ... Force CHELPG (Breneman) recommended radii.
C         100 ... Read in replacement radii for selected atom types as pairs
C                 (IAn,Rad) or (Symbol,Rad), terminated by a blank line.
C         200 ... Read in replacment radii for selected atoms as pairs
C                 (I,Rad), terminated by a blank line.
C               Initial radius of spheres to be placed around attractors
C               to "capture" the gradient trajectories.  The final
C               radius is then automatically optimized separately
C               for each atom.
C         0 ... 0.1
C        NM ... N.M = NM/10
C
C2IOp(44)
C     IOP(44)   IRC Type (L115, L123, L125):
C               0 ... Optimization coordinates (default 3 in L115 and L123),
C                     2 in L125).
C               1 ... Cartesian.
C               2 ... Internal (NYI in L123).
C               3 ... Mass-weighted Cartesian (NYI in L125).
C              M0 ... Initial interpolation in L125 (default 2).
C
C     IOp(44) ... L117: Maximum distance between a nucleus and its portion
C                 of the isosurface - used in Trudge only to
C                 eliminate, from the outset, points which clearly
C                 lie in another basin.  This parameter should be chosen
C                 with the parameter Cont in mind.
C           0 ... 10.0 au
C          NM ... N.M au = NM/10
C
C     IOP(44)   L121  Seed for random number generator (ISeed)
C              -1 ... Use system time initialize iseed
C                     (Note each run will give different results)
C               0 ... Use default seed value (ISeed = 398465)
C               N ... Set random number seed to N
C
C2IOp(45)
C     IOp(45):  Options for Link 123.
C               00  ... Which IRC Integrator to Use. Default=3, except 6
C                       for ONIOM QuadMac/Micro or 1 for GradientOnly.
C               01  ... Euler
C               02  ... LQA
C               03  ... HPC
C               04  ... GS2
C               05  ... DVV
C               06  ... Euler-based PC
C               10  ... Coordinate driving schemes
C              0xx  ... Is the integration being done on an empirical
C                       surface? Default=2.
C              1xx  ... Yes, this is an empirical surface. The energies
C                       and derivatives required for the IRC integration
C                       are NOT evaluated in this link. Instead it is
C                       assumed that an external program communicates the
C                       appropriate values with Link 402, etc. Also, the
C                       force constant matrix, when needed, is simply
C                       diagonalized, i.e. translation and rotation
C                       projections are NOT used. Also, all atomic masses
C                       are set to 1.
C              2xx  ... No.
C             0xxx  ... Order of magnitude for the step size relative to
C                       the integer value given with the StepSize=N option
C                       on the route line -- IOp(1/39). Default=2.
C             Nxxx  ... Integration step size is taken as IOp(1/39)*10^-N.
C            0xxxx  ... Whether or not energies reported in the final
C                       summary table should be given relative to the TS
C                       energy. Default=1.
C            1xxxx  ... Yes.
C            2xxxx  ... No.
C           0xxxxx  ... Whether or not statistics over coordinates should
C                       be converted to Angstoms/Degrees when reported in
C                       the summary table. Default=1.
C           1xxxxx  ... Yes.
C           2xxxxx  ... No.
C          0xxxxxx  ... Should a URVA input file be written? Default=2.
C          1xxxxxx  ... Yes.
C          2xxxxxx  ... No.
C         0xxxxxxx  ... Should IRC data be saved to the PES data structure
C                       file? Default=1.
C         1xxxxxxx  ... Yes.
C         2xxxxxxx  ... No.
C        0xxxxxxxx  ... When second-order methods are employed, should the
C                       Newton-Raphson step test be carried out? Default=1.
C        1xxxxxxxx  ... Yes.
C        2xxxxxxxx  ... No.
C
C2IOp(46)
C     IOp(46):  Order of multipoles in numerical SCRF:
C         0 ... Dipole
C         1 ... Quadrupole
C         2 ... Octapole
C         3 ... Hexadecapole.
C
C2IOp(47)
C     IOp(47):  Number of redundant internal coordinates to allow for.
C         0 ... Default:  Max(50000,MDV/(100*NStruc))
C         N ... N.
C
C2IOp(48)
C     IOp(48):  IRCMax control.
C         1 ... Do IRCMax
C        20 ... Include zero-point energy.
C
CIOp(49)
C     IOp(49):  Options to IRC path relaxation (IJKL)
C         L ... 2/1 dont/do optimize reactant structure.  Default: 1
C         K ... 2/1 dont/do optimize product structure.  Default: 1
C         J ... 3/2/1 dont/QST3/QST2 optimize TS structure (for QST input).  Default: 1
C         I ... 2/1 unimolecular/bimolecular reaction.  Default: unimolecular
C
C2IOp(51)
C     IOp(51):  Control for L125:
C       LMN ... Path second derivatives (default 213).
C         1     No 2nd derivatives, full path gradient only.
C         2     No 2nd derivatives, individual point gradients only.
C         3     Full path second derivatives.
C         4     Individual point second derivatives.
C        10     Identity matrix.
C       200     BFGS update.
C      K000     Path force definition (default 1)
C      1000     Orthogonal gradient
C      2000     Variational forces
C     J0000     Tangent/Curvature definition (default 3)
C     10000     Central difference
C     20000     Uphill difference
C     30000     Arc of Circle
C     40000     Cubic Spline
C
C2IOP(52)
C    IOP(52)  L101 and L120  Type of ONIOM calculation:
C              0/1 one layer, normal calculation
C               2  two layers
C               3  three layers
C              00  default (20)
C              10  include electrostatics in model systems using MM charges.
C                  (in case of three-layer ONIOM, this includes the charges
C                  in both the small model and the intermediate model system)
C              20  no electrostatics included in the model systems
C              30  As 10, but exclude the MM charges in the calculations on the
C                  smallest model system in case of a three-layer calculation.
C              40  As 10, but exclude the MM charges in the calculations on the
C                  intermediate model system in case of a three-layer calculation.
C             100  Do full square for testing.
C            N000  Use atomic charge type N-1 during microiterations.  The
C                  default is MK charges.
C           M0000  Type of link atoms for the MM calculation in QM/MM.
C                  0 Default (2).
C                  1 Conventional (Maseras) style.
C                  2 ONIOM style.
C          L00000  Whether to read additional charges with electronic embedding.
C                  0 Default (1).
C                  1 No.
C                  2 Yes.
C         K000000  Whether to create new entries in Common/Mol for the link atoms.
C                  0 No.
C                  1 Yes.
C
C2IOP(53)
C    IOP(53)  L120  Action of each invocation of L120:
C               0   do nothing
C               1   set up point MM on rwf from initial data
C               2   set up point MM on rwf from initial data and restore
C                   point MM on chk file if ONIOM data is present there.
C               3   restore point M from data on the rwf.
C               4   integrate energy
C               5   integrate energy and gradient
C               6   integrate energy, gradient, and hessian
C               7   Restore point MM from RWF but do not create a new
C                   model system.
C             NN0   Save necessary information (some rwf's,
C                   energy, gradients, hessian) of point NN
C                   of the ONIOM grid.
C                   NN = MaxLev**2 + 1 (currently 17) to restore real system.
C           MM000   Next point to do is MM.
C
C       Calc Level
C          High     4--7--9*
C                   |  |  |
C          Mid      2--5--8
C                   |  |  |
C          Low      1--3--6
C
C                   S  M  L   system size
C
C2IOp(54)
C     IOp(54) ... Whether to recover initial energy during IRCMax from chk file:
C           0 ... No.
C           1 ... Yes.
C
C2IOp(55)
C     IOp(55) ... L103:  Options for GDIIS:
C                 ICos*1000+IChkC*100+IMix*10+Method form.
C                 L115: IRC optimization.
C           0 ... Default, use gradients to find the next geometry.
C           1 ... Use displacements to find the next geometry.
C
C2IOp(56)
C     IOp(56) ... Set of atom type names to parse:
C           0 ... Accept any.
C           1 ... Dreiding/UFF.
C           2 ... Amber.
C           3 ... Amber allowing any symbol, for use with parameters in input stream.
C
C2IOp(57)
C     IOp(57) ... Whether to produce connectivity:
C           0 ... Default (4 if reading geom from chk file and
C                          connectivity is there, otherwise 3).
C           1 ... No.
C           2 ... Yes, read from input stream
C           3 ... Yes, generate connectivity.
C           4 ... Yes, read from checkpoint file.
C           5 ... Yes, read from rwf file.
C          10 ... Read modifications.
C         100 ... Connectivity input is in terms of z-matrix entries,
C                 including dummy atoms.
C
C2IOp(58)
C     IOp(58) ... IRCMax control in L115.
C           1 ... Do IRC rather than IRCMax
C          10 ... Store compound energy; default for IRCMax.
C          20 ... Zero-point energies are available during IRCMax.
C
C2IOp(59)
C     IOp(59) ... Update of coordinates in L103:
C           0 ... Default (1).
C           1 ... New versions (RedCar/RedQ2X); fall back to
C                 ORedCr/RedQX1 if RedCar fails
C           2 ... Old version (ORedCr/RedQX).
C           3 ... Old version (ORedCr/RedQX1).
C           4 ... New version (RedCar/RedQ2X) but fall back to
C                 ORedCr/RedQX if RedCar fails.
C       MMMM0 ... IAprBG in Red2BG.
C          10 ... Re-use eigenvectors of G only if exact.
C          20 ... Re-use eigenvectors of G if they are linearly independent.
C          30 ... Test old eigenvectors of G but don't re-use them.
C          40 ... Don't look at old eigenvectors.
C          50 ... Re-use eigenvectors of G if the RMS of the elements
C                 of the new G in the old null space is less than the
C                 threshold.
C        NN00 ... RMS < 10**(-N), default 4.
C       M0000 ... Default (1).
C       10000 ... Form G-inverse from the B eigen-values/vectors.
C       20000 ... Form G-inverse directly from G.
C       30000 ... Do G- via diagonalization of G (NYI).
C       40000 ... Do G- via SVD on B, returning only the eigenvectors with
C                 nontrivial eigenvalues.
C
C2IOp(60)
C     IOp(60) ... Interpret extra integer and fp values in z-matrix as scan information.
C           0 ... Default (No).
C           1 ... Yes.
C           2 ... No.
C
C2IOp(61)
C     IOp(61) ... How ONIOM should leave the rwf at the end of each geomtry:
C           0 ... Default (1).
C           1 ... Normal:  leave the rwf set up for the low-level
C                 calculation on the real system.
C           2 ... MOMM:  leave the rwf set up for the real system, but with
C                 NBasis and NBsUse for the high-level calculation on the
C                 model system.  Useful for treating the full system as having
C                 electrons only on the QM atoms. This is really a hack for
C                 two layer QM:MM ONIOM ADMP and should probably be generalized
C                 to behave like a ONIOM-PCM-A case.
C           3 ... Lowest level is MO, but normal setup at end.
C          00 ... Default (10, but 20 if doing EE microiterations).
C          10 ... Leave the charges file (605) from the best calculation that
C                 produced one.
C          20 ... Leave 605 in its normal state (present if from real,
C                 low-level).
C
C2IOp(62)
C     IOp(62) ... Counterpoise control.
C          NN ... NN fragments, NN < 50.
C         1NN ... Force use of new ghost atoms (default).
C         2NN ... Force use of old ghost atoms.
C        1xxx ... Counterpoise (default).
C        2xxx ... Fragment guess.
C        3xxx ... Morokuma energy decomposition.
C
C2IOp(63)
C     IOp(63) ... Step in counterpoise calculation:
C         LMM ... L  = order of derivatives (1=Energy, 2=Gradient, ...
C                 MM = 0                 -- supermolecule
C                      1-NFrag           -- fragments with ghost atoms
C                      NFrag+1 - 2*NFrag -- lone fragments
C
C2IOp(64)
C     IOp(64) ... Molecular mechanics force field selection:
C           0 ... None.
C           1 ... Dreiding.
C           2 ... UFF.
C           3 ... AMBER.
C           4 ... MM2 (NYI).
C           5 ... MM3 (NYI).
C           6 ... MMFF (NYI).
C           7 ... Quartic fitting field (NYI).
C         000 ... Use only hard-wired.
C         100 ... Use soft and hard-wired, hard-wired has priority.
C         200 ... Use soft and hard-wired, soft has priority.
C         300 ... Use only soft.  Lowest 2 digits then have no meaning.
C        0000 ... Do not read modifications to parameter set.
C        1000 ... Read modifications to parameter set.
C       00000 ... With soft parameters, abort when different parameters
C                 match to the same degree.
C       10000 ... Use the first when there are equivalent matches.
C       20000 ... Use the last when there are equivalent matches.
C
C     If IOp(67)=3, then the default is to apply soft parameters with
C     higher priority.
C
C2IOp(65)
C     IOp(65) ... Control of which terms are included in MM, corresponding to
C                 the 'classes' in FncInf.
C           0 ... Do all (default)
C           1 ... Non-bonded
C          10 ... Stretching
C         100 ... Bending
C        1000 ... Torsion
C       10000 ... Out-of-plane
C      100000 ... Stretch-bend
C
C2IOp(66)
C     IOp(66) ... Whether to generate QEQ charges, over-written the
C                 values in AtChMM, or to use the values already there.
C           0 ... Default (2, 1==> 221)
C           1 ... Do QEq.
C           2 ... Don't do QEq.
C          00 ... Default (20)
C          10 ... Do for atoms which were not explicitly typed.
C          20 ... Do for all atoms regardless of typing.
C         000 ... Default (200)
C         100 ... Do for atoms which have charge specified or defaulted to 0.
C         200 ... Do for all atoms regardless of initial charge.
C    MMMMM000 ... IType passed to QEq.
C
C2IOp(67)
C     IOp(67) ... Source of MM parameters.
C           0 ... Default:  2 if reading geom from chk file, else 1.
C           1 ... Generate here, reading from input if requested by IOp(64).
C           2 ... Copy from chk file.
C           3 ... Pick up non-standard parameters from chk file.
C
C2IOp(70)
C     IOP(70)    L118 Type of sampling (Nact)
C               0 ... Default (same as 3)
C               1 ... Orthant sampling
C               2 ... Classical microcanonical normal mode sampling
C               3 ... Fixed normal mode energy
C               4 ... Local mode sampling (thermal sampling based on RTemp)
C                 ... [Currently 0, 2, 3 and 4 working]
C              10 ... Read in Hessian from checkpoint for initial sampling
C
C2IOp(71)
C     IOp(71)   L115:  Whether to print out input files for each structure
C                      along an IRC:
C         0 ... No.
C         1 ... Yes.
C
C     IOp(71)   L118:  Hessian update:
C              -1 ... Gradient only
C              -2 ... Gradient+Hessian, but never calculate full H (only updates)
C               0 ... Full Hessian at every step
C              NN ... Try to do NN updates between full Hessians
C             000 ... Default updating (same as 300)
C             100 ... SR1 Hessian updating algorithm
C             200 ... PSB Hessian updating algorithm
C             300 ... Bofill's Hessian updating algorithm
C             400 ... Sqrt(Bofill) Hessian updating algorithm
C             500 ... No update
C            0000 ... Default (same as 1000)
C            1000 ... Reintegrated updated steps
C            2000 ... Suppress reintegration
C
C     IOp(71)   L123:  Hessian calculation:
C              -1 ... Gradient only.
C               0 ... Either only update or CalcAll as determined by IOp(10).
C               N ... Recalculate analytic Hessian every Nth calculation.
C
C2IOp(72)
C     IOp(72)   L121:  Lagrangian constrain method for ADMP (ICType)
C               Half*Gamma*Tr[(P*P-P)**2] + Lambda*[Tr(P)-Ne] + Eta*Tr(P*P-P)
C               0   ... Default
C                       Same as  7 if no Mass-Weighting (IOp(76) < 0)
C                       Same as 10 if Mass-Weighting (IOp(76) > 0)
C               1   ... Use Lambda and Eta only.  (Gamma=0)
C               2   ... Use Lambda, Eta, Gamma.  Gamma = .2
C               3   ... Use Lambda, Eta, Gamma.  Gamma = 1.
C               Constraints for scalar Mass case:
C               4   ... Use exact constraint Sum(ij)[Vij*(P**2-P)ij]
C               5-7 ... Iterative Scheme same as 4. Different initial
C                       guesses.
C                       7 is default for scalar mass case.
C               Constraints for tensorial Mass:
C               8-11... Mass-weighting constraints. Documentation
C                       maybe found in DVelV1.
C                       10 is default.
C
C     IOp(72)   L124:  ISolv in L124??
C
C2IOp(73)
C     IOp(73)   L123:  Whether to compute projected frequencies:
C               0/1/2 Default(No)/No/Yes.
C
C     IOP(73)   L118 and L121:  Initial Kinetic energy of the Nuclei (EStrtC)
C               0 ... Default (.1 Hartree)
C             N>0 ... N*micro-Hartree
C             N<0 ... 0.0 Hartree
C
C2IOp(74)
C     IOp(74) ... Charge scaling for charge embedding in ONIOM.
C      IJKLMN ... 6th through 1st nearest neighbors of current layer
C                 scaled by I*0.2, J*0.2, etc.  0 ==> 5 => IAtTyp=6 (no scaling);
C                 all layers are scaled by at least as much as ones
C                 farther out.
C                 The default is 500.
C           M   Factor for charges one bond away from link atom
C          L0   Factor for charges two bonds away from link atom
C         K00   Factor for charges three bonds away from link atom
C     ..IJ etc.
C               The actual factors used for each value of IAtTyp are:
C               1: 0.0
C               2: 0.2
C               3: 0.4
C               4: 0.6
C               5: 0.8
C               6: 1.0
C
C2IOp(75)
C     IOP(75)   ADMP control flag (ICntrl)
C               0 ... Standard ADMP
C               1 ... Read converged density at every step
C               2 ... Fix the nuclear coordinates
C               3 ... Test time reversability (MaxStp must be even)
C              00 ... Default (20).
C              10 ... Read stopping parameters from input.
C              20 ... Do not read stopping parameters.
C
C2IOp(76)
C     IOP(76)   +/- XXXXZYYYY
C               Ficticous  electron mass (EMass)
C               YYYY ... Default (1000)
C                        IOP(76)>0 ... YYYY*.0001 AMU   MW core functions
C                                      more than valence functions.
C                        IOP(76)<0 ... YYYY*.0001 AMU.  Use uniform scaling
C                                      for all basis functions
C                        (Note YYYY > 9999 makes no sense)
C               Z    ... Mass-weighting option. If IOP(76)<0, Z is meaningless.
C               XXXX ... If PBC: Mass of Box Coordinates
C                        (BoxMas) = XXXX*.0001 AMU
C                        BoxMas=0 ... Box coordinates not propagated (default).
C
C2IOp(77)
C     IOP(77)   Initial Kinetic energy of the density matrix (EStrtP)
C               (For UHF, Alpha and Beta each get half this energy)
C               And Option Number to compute initial kinetic energy.
C               Format of Input: XXYYYY (six digits)
C               IWType = XX
C               N = YYYY
C               (For UHF, Alpha and Beta each get half this energy)
C               0 ... Default (0.0 Hartree)
C             N>0 ... N*micro-Hartree
C               IWType is used to figure out how the initial velocity is
C               is computed (in gnvelp).
C               If XXYYYY < 0 : Initial velocity = 0.0 Hartee
C                               (i.e., currently same as N=0 above)
C
C2IOp(78)
C     IOp(78)   Sparse in L121?
C              -N ... Sparse here with cutoff 10**(-N), full elsewhere
C               0 ... Use full matrices or spase based on standard settings.
C               1 ... Use sparse fixed form
C
C2IOp(79)
C     IOp(79)   IRCMax convergence in L115
C               Stopping criteria in L121.
C               L118:
C     IOP(79)   0 ... Default, do not analyze pure rotation and vibration
C                     for polyatomic molecules
C               1 ... Do pure rotational and harmonic normal mode analysis
C                     for polyatomic molecule; EBK theory for diatomic
C                     vibrational analysis (require equilibrium information
C                     for each of the polyatomic molecules from saved
C                     checkpoint files and morse parameters for diatomic
C                     molecule)
C               2 ... Do pure rotational and harmonic normal mode analysis
C                     for polyatomic molecule; Local harmonic vibrational
C                     analysis for diatomic molecule (require equilibrium
C                     information for each of the fragments from saved
C                     checkpoint files)
C               3 ... Do pure rotational analysis and select the best
C                     normal mode analysis methods (harmonic and
C                     anharmonic) for polyatomic molecule; Local harmonic
C                     vibrational analysis for diatomic molecule (require
C                     equilibrium information for each of the fragments
C                     from saved checkpoint files)
C              00 ... Default, use default stopping criteria
C              10 ... Use user specified stopping criteria
C
C2IOp(80)
C     IOp(80)   L106: 0/1/2 Cartesian/Normal mode/Internal coordinate
C                     differentiation.  2 is NYI.
C
C     IOp(80)   L118: .eq.1 to surpress the 5th order correction after
C                      surface hop has been made in Trajectory Surface
C                      Hopping calculations. Needs also iop(10/80=1)
C
C     IOP(80)   L121: Nuclear Kinetic Energy Thermostat Option.
C               (Currently only Velocity scaling is implemented)
C               0 ... No Thermostat.
C         11XXXXX ... Velocity scaling, but only for the first
C                     XXXXX simulation steps.
C                     (This options is useful, if thermostating in only
C                     required during equilibration.
C         1000000 ... Velocity scaling, all the way through the simulation.
C
C2IOp(81)
C     IOp(81)   Nuclear KE thermostat in ADMP -- temperate is checked and
C               scaled every IOp(81) steps.
C
C2IOp(82)
C     IOp(82)   Temperature for nuclear KE thermostat in L121.
C
C2IOp(83)
C     IOp(83)   Whether to read in frequencies for electric and magnetic
C               perturbations.
C         0 ... Default (No).
C         1 ... Yes.
C         2 ... No.
C        00 ... Default (10).
C        10 ... Ensure that the static case is done, by putting a zero
C               frequency at the beginning.
C        20 ... Do not put a zero frequency at the beginning.
C       000 ... Default (100).
C       100 ... Sort the frequencies into increasing order.
C       200 ... Do not sort the frequencies (for debugging).
C
C2IOp(84)
C     IOp(84) ... Differentiation of frequency-dependent properties.
C           0 ... No.
C           N ... Mask for which properties on file 721 will be
C                 differentiated.
C
C2IOp(85)
C     IOp(85) ... Bad gap calculation in PBC ADMP:
C           0 ... Default (No).
C           1 ... Diagonalizae Fock matrix to get band gap, evolution, etc.
C           2 ... No.
C
C2IOp(86)
C     IOp(86) ... Printing for NMR for ONIOM.
C           0 ... Default (1).
C           1 ... Print tensors and eigenvalues.
C           2 ... Print eigenvectors as well.
C
C2IOp(87)
C     IOp(87) ... ONIOM integration of density.
C           0 ... Do not integrate.
C           1 ... Integrate current densities.
C           2 ... Integrate densities specified by following digits:
C          K0 ... Density to use from gridpoint 1
C         L00 ... Density to use from gridpoint 2
C        M000 ... etc.
C                 K,L,M,etc: 0: SCF
C                            1: MP first order
C                            2: MP2
C                            3: MP3
C                            4: MP4
C                            5: CI one-particle
C                            6: CI
C                            7: QCI/CC
C                            8: Correct to second order
C
C2IOp(88)
C     IOp(88) ... Whether to read in atomic masses (isotopes):
C           0 ... Default (1 if geometry read from input,
C                          4 if geometry read from chk)
C           1 ... Use most abundant isotopes.
C           2 ... Read isotopes from input.  The temperature
C                 and pressure are read first, for backwards
C                 compatibility.
C           3 ... Read isotopes from rwf.
C           4 ... Read isotopes from chk.
C           5 ... (Generated internally) isotopes were read
C                 from chk file during guess=input.
C
C2IOp(89)
C     IOp(89) ... Maximum allowed deviation from average nuclear KE
C                 during ADMP, in Kelvin.
C
C2IOp(90)
C     IOP(90)   To read in the velocity in cartesian coordinates
C
C     IOP(91)   Nuclear Kinetic Energy Thermostat Option.
C               Average energy (in microhartree) to be maintained during
C               Simulation, as required by IOP(80).
C
C     IOP(92)   Thermostat Option.
C               Maximum allowed deviation from average nuclear KE specified
C               in IOP(81). Also in microhartree.
C
C2IOp(93)
C     IOp(93) ... QM/MM TS vector guess.
C         000     Default (112)
C           1     Retrieve from checkpoint file if available, otherwise
C                 diagonalize QM Hessian or read from input
C           2     Do not try to retrieve from checkpoint file.
C          10     Diagonalize QM contribution to Hessian
C          20     Read from input.
C         N00     How to deal with 'suspicious RFO solutions' (default is 1).
C                 1. Just take the step.
C                 2. Check if there is an eigenvector with wrong curvature.
C                    If there is, flip its eigenvalue.
C                 3. Check if there is an eigenvector with wrong curvature.
C                    If there is, take a small step into this direction,
C                    followed by a linear search. This should step out (or
C                    stay in) the wrong region, and fix the eigenvalue.
C
C2IOp(94)
C     IOp(94) ... Davidson control for quadratic micro-iterations.
C          OP     Number of initial guess vectors (4).
C        MN00     Iteration to scale down number of vectors (5).
C      KL0000     Factor to scale down with; 1 for no scaling (2).
C     J000000     Whether to do geometry steps when the CG is done (2).
C                 1 ... Make the CG steps.
C                 2 ... No displacements.
C                 3 ... Only do displacements at first guess.
C    I0000000     Whether to re-use previous RFO solution or to regenerate guess (1).
C                 1 ... After first step, use previous solution as guess.
C                 2 ... Regenerate guess each time.
C                 3 ... Use previous lowest root, and regenerate remainder.
C   H00000000     Whether (1, default) or not (2) to add 0,..,0,1 guess vector.
C
C2IOp(95)
C     IOp(95) ... RFO/Davidson control for quadratic micro-iterations.
C          MM     Convergence (7).
C       LLL00     <0: Regular Davidson.
C                  0: Only check convergence on first vector, and iterate the
C                  lowest root only. Use all the intermediate vectors.
C                 >0: Only check convergence on first vector, and iterate the
C                  new vectors LLL times with the explicit last row/collumn. This
C                  is specifically appropriate for RFO. The last row/collumn
C                  of the Hessian comes after the diagonal elements.
C
C2IOp(96)
C     IOp(96) ... Options for generating intial guess vectors for RFO/
C                 Davidson diagonalization in coupled QM/MM macro steps.
C                 Note that other RFO/Davidson diagonalization conrols for
C                 coupled QM/MM macro steps are available in IOp(97).
C        GHIJKLLMM    MM  Number of initial guess vectors to get from CG
C                         steps. The default is 0.
C                     LL  Number of initial guess vectors from the diagonal
C                         of the QM block (4). The default is 4.
C                      K  Add 0,..,0,1 guess vector?
C                         0 ... Default: K=1.
C                         1 ... Yes.
C                         2 ... No.
C                      J  Add the gradient vector to the guesses?
C                         0 ... Default: J=1.
C                         1 ... Yes.
C                         2 ... No.
C                      I  Pre-diagonalize a Hessian/RFO matrix without non-
C                         bonding contributions? Note that this control is
C                         only valid for IOp(98)>3; otherwise I is ignored.
C                         0 ... Default: I=1.
C                         1 ... Yes.
C                         2 ... No.
C                      H  Scale factor for the size of the Davidson sub-
C                         space in early iterations.
C                         0 ... Default: H=4.
C                         1 ... Same as no scaling.
C                         H ... Use a sub-space in early iterations that is
C                               H times the number of requested vectors.
C                      G  Number of vectors to solve using Davidson
C                         diagonization.
C                         0 ... Default: G=1.
C                         G ... Solve for G vectors.
C
C2IOp(97)
C     IOp(97) ... RFO/Davidson diagonalization control for coupled QM/MM
C                 macro step. Note that other RFO/Davidson diagonalization
C                 conrols for coupled QM/MM macro steps are available in
C                 IOp(96).
C       GHIJKLMM      MM  Convergence in Davidson iterations. Convergence
C                         is set to 10**-MM. The default value is MM=5.
C                      L  What is being diagonalized?  This option is set
C                         explicitly in subroutines before calling the
C                         Davidson diagonalization code. Therefore, the
C                         value set in this IOp is ignored and serves only
C                         as a place holder.
C                         0 ... the Hessian.
C                         1 ... the augmented-Hessian/RFO matrix.
C                      K  Check convergence on which roots?
C                         0 ... Default: For L=0, K=2; For L=1, K=1.
C                         1 ... Check convergence on lowest root only.
C                         2 ... Check convergence on all roots.
C                      J  Appears to be unused.
C                      I  Number of Davidson iterations to store.
C                         0 ... Default: Keep all iterations.
C                         I ... Keep only the last iteration.
C                      H  Number of new vectors to create in each Davidson
C                         iteration.
C                         0 ... Default: For L=0, H=1; For L=1, H=2.
C                         1 ... Iterate all roots/vectors.
C                         2 ... Iterate lowest root/vector only up to the
C                               maximum number of iterations that are
C                               default in the Davidson code (ignores J
C                               above) and keeping vectors from all
C                               iteratations (ignores I above).
C                         3 ... Iterate lowest root/vector only. Note that
C                               this option is essentially the same as H=2,
C                               though J and I option settings are honored.
C                      G  Initial approximation to use for Davidson
C                         diagonalization.
C                         0 ... Default: G=2.
C                         1 ... Use a diagonal Hessian approximation
C                               together with the gradient vector. This is
C                               best used in RFO applications.
C                         2 ... Use the inverted Hessian for the QM block.
C                         3 ... Use a diagonal Hessian approximation.
C
C2IOp(98)
C     IOp(98) ... Control of quadratic micro-iterations and coupled QM/MM
C                 quadratic macro step.
C          <0 ... Do not use dynamic convergence criteria for the micro-iterations.
C           0 ... Default(11).
C           1 ... Regular non-coupled macro step.
C           2 ... Coupled macro step, full diagonalization.
C           3 ... Coupled macro step, direct /w full Hessian incore.
C           4 ... Coupled macro step, direct /w MM Hessian incore.
C           5 ... Coupled macro step, fully direct.
C           6 ... Go through the QuadMac machinery, but use the fully integrated
C                 ONIOM hessian to calculate the Hessian-vector products.  Switch
C                 to regular microiterations at points without analytic second
C                 derivatives.
C           7 ... Fully Quadratic at 2nd derivative points (1st in
C                 CalcFC, all in CalcAll), QuadMac with integrated
C                 Hessian at non-2nd derivative points.
C          10 ... Regular micro-iterations.
C          20 ... Quadratic micro-iterations, full diagonalization.
C          30 ... Quadratic micro-iterations, direct /w prepared Hessian incore.
C          40 ... Quadratic micro-iterations, direct /w raw MM Hessian incore.
C          50 ... Quadratic micro-iterations, fully direct.
C          60 ... No micro-iterations.
C
C2IOp(99)
C     IOp(99) ... Accuracy used for the non-bonded interactions in the Hessian-vector
C                 product calculations in the fully direct algorithms. Settig this IOp
C                 does not affect the MM energy and gradient calculations; only the
C                 direct evaluation in the QuadMac optimization step.
C                 When IOp(99).lt.0, the full accuracy is used.
C          MM ... Maximum multipole level (8)
C       LLL00 ... Boxsize in FMM (12)
C    KKK00000 ... Cutoff in Vanderwaals (30)
C
C2IOp(101)
C2IOp(102)
C2IOp(103)
C2IOp(104)
C
C     IOp(101-104) ... Phase control in L115 and L118
C               N1, N2, N3, N4
C
C2IOp(105)
C     IOp(105)   Reaction direction
C              00 ... Default (Same as 10)
C              10 ... Forward direction
C              20 ... Reverse direction
C
C               Damped-Velocity Verlet (DVV) options for Dynamic Reaction
C               Path Following
C               0 ... Default (Same as 2)
C               1 ... Use DVV
C               2 ... Do not use DVV
C              00 ... Default (Same as 10)
C              10 ... Follow the rxn path in the forward direction
C              20 ... Follow the rxn path in the reverse direction
C             000 ... Default (Same as 200)
C             100 ... Time step correction not used
C             200 ... Time step correction used but not to recalculate
C                     current DVV step
C             300 ... Time step correction used and current DVV step
C                     recalculated
C            0000 ... Default (Same as 1000)
C            1000 ... Use DVV stopping criteria
C            2000 ... Do NOT use DVV stopping criteria
C
C2IOp(106)
C     IOp(106)   Damping constant for DVV Dynamic Rxn Path following (v0)
C               0 ... Default v0=0.04 (N=400)
C               N ... v0 is set to N*0.0001
C
C2IOp(107)
C     IOp(107)   Error tolerance for DVV time step correction (Error)
C               0 ... Default Error=0.003 (N=30)
C               N ... Error=N*0.0001
C
C2IOp(108)
C     IOp(108)   Gradient magnitude for DVV stopping criteria (Crit1)
C               0 ... Default (N=15)
C              <0 ... Turn off this test.
C               N ... N*0.0001
C
C2IOp(109)
C     IOp(109)   Force-Velocity angle for DVV stopping criteria (Crit2)
C               0 ... Default (90 Degrees)
C              <0 ... Turn off this test.
C               N ... Use N Degrees
C
C2IOp(110)
C     IOp(110) ... Scaling of rigid fragment steps during microiterations.
C            0 ... Do not scale
C            1 ... Scale with 1/NRA   (NRA = number of atoms in fragment)
C            2 ... Scale with 1/Sqrt(NRA)
C     -n ... Scale with 1/n
C
C2IOp(111)
C     IOp(111) ... Step-size to use with steepest descent when L103 is
C                  having trouble:
C           -N ... Scale up to RMS step of N/1000 if DXRMS is less.
C           -1 ... Effectively disables the scaling
C            0 ... Default (50)
C            N ... Scale up or down to maximum change in a variable of N/1000
C
C2IOp(112)
C     IOp(112) ... Temperature for thermochemistry.
C            0 ... Default (standard temperature, unless read in).
C            N ... N/1000 degrees.
C           -N ... N/1000000 degrees.
C
C2IOp(113)
C     IOp(113) ... Pressure for thermochemistry.
C            0 ... Default (1 atomosphere, unless read in).
C            N ... N/1000 atmospheres.
C           -N ... N/1000000 atmospheres.
C
C2IOp(114)
C     IOp(114) ... Scale factor for harmonic frequencies for use in
C                  thermochemistry and harmonic vibration-rotation
C                  analysis.
C            0 ... Default (1 unless specified by IOp in overlay 7 or
C                  read in).
C            N ... N/1000000.
C
C2IOp(115)
C     IOp(115) ... Compression for MOMM quadratic steps:
C            4 ... Second derivatives generated over active atoms, with real
C                  system terms done iteratively during microiterations.
C       N.ne.4 ... Full second derivative matrices are used.
C
C2IOp(116)
C     IOp(116) ... Options for ONIOM Conical intersections.
C            0 ... Default (111, whatever that means).
C
C2IOp(117)
C     IOp(117) ... MMVB control:
C            0 ... No MMVB
C            1 ... Do MMVB calculation.
C           11 ... Do MMVB dynamics allowing for surface hop.
C          000 ... Activate Bridging Corrections to Kij and Qij (Default)
C          100 ... Does Not Activate the Bridging Corrections
C          200 ... Activates Marco's Bridging Corrections
C          300 ... Activate GAUCA-MMVB Bridging Corrections (for debugging only)
C          400 ... Activate GAUCA-MMVB Marco's Bridging Corrections (for debugging only)
C         0000 ... No Conical Intersection Search (Default)
C         1000 ... Search For Conical Intersection For (LRoot,LRoot-1)
C        00000 ... Do not Delete any Coulombic Qij.
C        10000 ... Specify Qij to delete.
C       000000 ... E(MMVB)=E(MM)+E(VB), G(MMVB)=G(MM)+G(VB) (Default)
C       100000 ... E(MMVB)=E(VB), G(MMVB)=G(VB)
C       200000 ... Option 2 in the original GAUCA-MMVB code (NYI)
C      0000000 ... Solve VB problem numerically with Lanczos
C      1000000 ... Solve VB problem analytically
C     00000000 ... Hartree-Waller functions for singlets or triplets
C     10000000 ... Slater determinant
C    000000000 ... Use correct conversion factors
C    100000000 ... Use GAUCA-MMVB conversion factors (for debugging only)
C
C2IOp(118)
C     IOp(118) ... Dump structures for each ONIOM system formatted as input.
C            0 ... Default (No).
C            1 ... Yes.
C
C2IOp(119)
C     IOp(119) ... Control Initial Lanczos Vector (ILzVec)
C           -1 ... Read guess by card in input file
C           -2 ... Use the largest elements of H as a guess
C           -3 ... Use the five largest contributions of H as a guess
C           0x ... For Opt, IRC, dynamics read guess from previous cycle
C           1x ... For Opt, IRC, dynamics generate a fresh guess for each cycle
C
C2IOp(120)
C     IOp(120) ... Flags for QM:QM embedding.  NOTE: The standard embedding
C                  flags must also be set in the same way as necessary for
C                  QM:MM embedding calculations.
C            0 ... Default - Same as 1.
C            1 ... Use Mulliken charges.
C           -1 ... Use the nuclear charge stored in array AtmChg.
C           -2 ... Set the charges to zero.
C           00 ... Default (Same as 20).
C           10 ... Just use the charges that already reside in AtChMM.
C           20 ... Update AtChMM using current atomic charges on the RWF.
C                  This option is only employed immediately following
C                  low-level real-system sub-calculations.
C
C2IOp(121)
C     IOp(121) ... Extra items to differentiate (L106):
C            0 ... Default, none.
C            1 ... Differentiate AO density and Fock matrices.
C
C2IOp(122)
C     IOp(122) ... Read secondary structure information:
C            0 ... Default (4 or 3 if reading geom from chk or rwf file,
C                           otherwise 2).
C            1 ... No.
C            2 ... Yes, read from input stream if any residue information
C                  was provided on the atom definition lines.
C            3 ... Yes, read from rwf.
C            4 ... Yes, read from chk.
C
C2IOp(123)
C     IOp(123) ... Version of /Mol/ to save on disk.
C            0 ... Default (current, version 2)
C            1 ... Version 1 (flag -12345).
C            2 ... Version 2 (flag -12346).
C          N<0 ... Flag value N.
C
C2IOp(124)
C     IOp(124) ... Flavor of ONIOM-PCM to use:
C            1 ... "A", reaction field from the ONIOM integrated density;
C            2 ... "B", reaction field from the real system low level;
C            3 ... "C", reaction field in the real system low level only;
C            4 ... "X", reaction field in all subcalculations using the
C                       real system cavity;
C           0x ... Default (same as 1 unless a semiempirical method is
C                  involved);
C           1x ... Integrate the density for ONIOM-PCM-A;
C           2x ... Integrate the potential for ONIOM-PCM-A;
C          1xx ... Flag to indicate ONIOM-PCM-X as first iteration of
C                  ONIOM-PCM-A.
C
C2IOp(125)
C     IOp(125) ... Solvent charge distribution to add to Hamiltonian:
C            0 ... None
C            1 ... Read charges and DBFs from input stream in input
C                  orientation
C            2 ... Read from RWF.
C            3 ... Read from Chk.
C            4 ... Same as 1.
C            5 ... Read charges and DBFs from input stream in standard
C                  orientation
C           10 ... Force units of Angstroms for coordinates.
C           20 ... Force units of Bohr for coordinates.
C
C2IOp(126)
C     IOp(126) ... Whether to read an input section with atom opt/freeze
C                  information.
C            0 ... Default (2).
C            1 ... Yes.
C            2 ... No.
C
C2IOp(127)
C     IOp(127) ... Unused.
C
C2IOp(128)
C     IOp(128) ... Initial microiterations in L120 before first QM step,
C                  and microiterations in L120 during numerical
C                  differentiation in L103:
C            0 ... Default (No).
C            1 ... Yes.
C            2 ... No.
C
C2IOp(129)
C     IOp(129)  L123:  Whether or not the final statistics printed in the job
C               summary should include coordinate values at each IRC step.
C               0   ... Default (none).
C               1   ... Read user-defined stats from the input file.
C              0x   ... Default (do not report all Cartesian coordinates)
C              1x   ... Report all Cartesian coordinates.
C             0xx   ... Unused.
C            0xxx   ... Default (do not report bond coordinates)
C            1xxx   ... If redundant internals are on the RWF, then report
C                       values for bond coordinates along the IRC.
C           0xxxx   ... Default (do not report angle coordinates)
C           1xxxx   ... If redundant internals are on the RWF, then report
C                       values for angle coordinates along the IRC.
C          0xxxxx   ... Default (do not report dihedral coordinates)
C          1xxxxx   ... If redundant internals are on the RWF, then report
C                       values for dihedral coordinates along the IRC.
C
C2IOp(130)
C     IOp(130)  Eigenvector number to be followed during coordinate driving
C               jobs [IOp(1/19=10)].
C               0   ... Default (1)
C               N   ... Follow input structure Hessian eigenvector number N
C
C2IOp(131)
C     IOp(131)  Options for corrector integration in predictor-corrector
C               IRC calulcations (Link 123).
C              00   ... Should the corrector integration scheme be run in
C                       an (macro-cycle) iterative fashion? Default=2.
C             -01   ... After each corrector integration, evaluate the
C                       actual energy and derivatives, but do not actually
C                       use these. The final IRC will be the same as 1.
C              01   ... Never check convergence of the corrector
C                       integration. Always do one corrector integration
C                       for each predictor integration.
C              02   ... Always check for convergence of the corrector
C                       integration end point. Convergence is acheived
C                       when the change in corrector integration end point
C                       geometry is less than the convergence criteria
C                       indicated by IOp(7).
C              03   ... Same as 2, but never accept convergence after the
C                       first corrector integration at a point.
C              10   ... This flag iteratively refines the DWI fitted
C                       surface if multiple corrector integration
C                       macro-cycles are taken by adding the largest
C                       standard deviation point from the previous BS
C                       cycle.
C             -11   ... This flag forces a PES evaluation step after each
C                       corrector integration. This is similar option -1,
C                       except that the actual energy and derivatives are
C                       used for the next predictor step rather than the
C                       values on the DWI fitted surface at the corrector
C                       end-point.
C             0xx   ... Should DWI surfaces employ numerical thrid-order
C                       terms in Taylor series? Default=1.
C             1xx   ... Use DWI surface with Taylor series expansions
C                       truncated at second-order.
C             2xx   ... Use DWI surface with Taylor series expansions
C                       truncated at third-order using numerical
C                       third-derivatives.
C            Nxxx   ... What power should be used for DWI weights that
C                       include delta-x raised to an even power. The value
C                       of this option sets that power to 2*N. Default=N=1.
C           0xxxx   ... How should the first step off of the transition
C                       structure point be handled in corrector integration
C                       cycles?  Note, in all cases, the transition vector
C                       is used to define the tangent at the transition
C                       structure. Default=2.
C           1xxxx   ... Run the requested number of Euler steps in the
C                       standard way. Only the first Euler step taken uses
C                       the transition vector.
C           2xxxx   ... Take a large step off of the transition structure
C                       point along the transition vector. This step is
C                       taken to be half of the total requested step size
C                       given by TotStp.
C           3xxxx   ... This the same as option 2 in concept. The only
C                       difference is that the first step off of the
C                       transition structure is taken as one-third the
C                       total requested step size given by TotStp.
C           4xxxx   ... This the same as option 2 in concept. The only
C                       difference is that the first step off of the
C                       transition structure is taken as one-fourth the
C                       total requested step size given by TotStp.
C          0xxxxx   ... Should update vectors be used in DWI fits if
C                       possible? Default=1.
C          1xxxxx   ... Yes, when possible.
C          2xxxxx   ... Never.
C
C2IOp(132)
C     IOp(132)  Whether to check for divalent link atoms in ONIOM input:
C         0 ... Default(yes).
C         1 ... Yes.
C         2 ... No.
C
C2IOp(133)
C     IOp(133)  Suppress integration/restore of quantities for Polar=Raman
C               and Polar=ROA ONIOM jobs:
C         0 ... Default (1).
C         1 ... Do not restore or integrate forces, force constants,
C               static electric or magnetic field derivatives.
C         2 ... Restore all files.
C
C2IOp(134)
C     IOp(134)  Number of Fourier coefficients for L125.
C
C2IOp(135)
C     IOp(135)  MM Non-bonded switching function.
C
C2IOp(136)
C     IOp(136)  MM vdW outer cutoff in Angstroms.
C
C2IOp(137)
C     IOp(137)  MM Coulomb outer cutoff in Angstroms.
C
C2IOp(138)
C     IOp(138)  MM soft cutoff range (applied to both vdW and Coulomb,
C               in Angstroms.
C
C2IOp(139)
C     IOp(139)  Number of MM microiterations allowed.
C         0 ... Default, based on NAtoms but at least 5000.
C         N ... N.
C
C2IOp(140)
C     IOp(140)  Whether to restore the real system from chk file
C         0 ... Default, yes if ONIOM.
C         1 ... Yes.
C         2 ... No.
C
C2IOp(141)
C     IOp(141)  Control of error choice during GEDIIS:
C         0 ... Default (1).
C         1 ... Use RFO steps as error vectors, using the NR step at
C               a point if the RFO fails or gives a Hessian eigenvector.
C         2 ... If RFO fails for any point, use the gradient for all
C               points.
C         3 ... Always use the gradients as the errors.
C         4 ... Use RFO steps as error vectors unless any RFO fails and
C               unless Hessian is marked as untrustworthy; then use
C               gradients instead.
C         5 ... Drop back to RFO (no DIIS) if Hessian is untrustworthy.
C               NR steps are used if RFO fails during DIIS.
C
C2IOp(142)
C     IOp(142)  Whether to copy MM charges to link atoms:
C         0 ... Default (3 if QEq is done; otherwise 1).
C         1 ... Copy if link atom charge is zero.
C         2 ... Do not copy.
C         3 ... Always copy.
C
C2IOp(143)
C     IOp(143)  Hessian during IRC restarts
C         0 ... No change in when Hessian is done.
C         1 ... Do Hessian at first new point in each direction.
C
C2IOp(144)
C     IOp(144)  Whether to analyze residue geometry.
C         0 ... Default (Yes, if secondary structure present and
C               NAtoms<=10000)
C         1 ... Yes.
C         2 ... No.
C
C2IOp(145)
C     IOp(145)  Controls for the new internal coordinate data structure.
C         0 ... Should the new internal coordinate data structure be set-up
C               for use? (Default=1).
C         1 ... No.
C         2 ... Yes.
C        0x ... Which coordinate derivative terms should be included in the
C               data structure? (Default=3)
C        1x ... Include only values in the definitions list.
C        2x ... Include only values and their first-derivatives for
C               B-matrix formation.
C        3x ... Include values and first- and second-derivatives for
C               B-matrix and B-matrix derivatives.
C       0xx ... How should the new code define the coordinate set used when
C               building the data structure? (Default=1)
C       1xx ... Use the old coordinate generation code to build IZRed, then
C               have CrdDef build the data structure using IZRed for the
C               coordinates and build derivatives wrt Cartesian
C               coordinates.
C       2xx ... Read the coordinate definitions from the input file using
C               "Format 2" defined in routine CrDf01, and then
C               automatically differentiate the active coordinate
C               derivatives as needed.
C       3xx ... Read the coordinate definitions from the input file using
C               "Format 2" defined in routine CrDf01, and then read in the
C               derivative definitions (one set for IOp(145)=x1x; two sets
C               for IOp(145)=x2x) from the input file using "Format 9"
C               defined in Routine CrDf01.
C       4xx ... Read the coordinate definitions from the input file using
C               "Format 2" defined in routine CrDf01, and then read in the
C               first-derivative definitions using "Format 9" as defined in
C               Routine CrDf01. If second-derivatives are requested
C               (IOp(145)=x2x) then read in the second-derivative
C               definitions using "Format 8" as defined in Routine CrDf01.
C       5xx ... Read the coordinate definitions from the input file using
C               "Format 2" defined in routine CrDf01, and then read in the
C               first-derivative definitions using "Format 8" as defined in
C               Routine CrDf01. If second-derivatives are requested
C               (IOp(145)=x2x) then also read in the second-derivative
C               definitions using "Format 8" as defined in Routine CrDf01.
C       6xx ... Read the coordinate definitions from the input file using
C               the symbolic form, and then automatically differentiate the
C               active coordinate derivatives as needed.
C       7xx ... Define (the full set of) distance matrix coordinates. This
C               option can be combined with 2xxx to use inverse distance
C               matrix coordinates.
C      0xxx ... When building the internal coordinate definitions, should
C               any systematic modifications be done?  (Default=1)
C      1xxx ... No. Simply define the coordinates to be the same as they
C               would have been using the old coorindate code.
C      2xxx ... Invert all bond stretch coordinates.
C
C2IOp(146)
C     IOp(146)  Control of SCVS:
C         0 ... Default (03045011).
C         1 ... Include forces in virial ratio.
C         2 ... Do not include forces in virial ratio.
C        1x ... Use Murdoch's extrapolation.
C        2x ... Do not use Murdoch's extrapolation.
C       Nxx ... Apply SCVS when max force on nuclei is below 10**(-N).
C      Mxxx ... The convergence on the virial Eta.  Default is 5.
C     Lxxxx ... The convergence threshhold on |E| * (2*Eta - 2), where
C               |E| is the magnitude of the kinetic energy is 10**(-M).
C               Default is 3.
C    Kxxxxx ... The initial order of extrapolation is 10**L
C JJJxxxxxx ... Maximum of KKK SCVS iterations.
C
C2IOp(147)
C     IOp(147) ... Control of EDA:
C            0 ... Default (2).
C            1 ... Symmetrically orthogonalize before ES(X).
C            2 ... Do not symmetrically orthogonalize before ES(X).
C
C2IOp(148)
C     IOp(148) ... Storage of derivatives by L106:
C            0 ... Default (1).
C            1 ... Normal storage of numerical first derivatives.
C            2 ... Store numerical first and diagonal second derivatives.
C
C2IOp(149)
C     IOp(149) ... Spring constant in L125:
C            0 ... Default (zero).
C            N ... N*0.01 Hartree/Bohr^2
C
C2IOp(150)
C     IOp(150) ... Initial scale factor in L112.
C            0 ... Default (1.0).
C            N ... 1 + N/100000000.
C
C2IOp(151)
C     IOp(151) ... How many vectors with negative curvature to use in
C                  downhill step during minimization:
C            0 ... Default (3).
C           -1 ... None, do regular RFO step.
C            N ... Up to N vectors.
C
C2IOp(152)
C     IOp(152) ... Control of MaxStp (allocated max number of steps in L103):
C            0 ... Default:  compute based on number of variables, NStep, etc.
C          N>0 ... Make MaxStp at least N.
C          N<0 ... Make MaxStp at most -N.
C
C2IOp(153)
C     IOp(153) ... Whether to fill OT file in L103.
C            0 ... Default:  do so for cartesian and red. int. coord. opts.
C            1 ... Yes.
C            2 ... No.
C
C2IOp(154)
C     IOp(154) ... Linear angle test during internal coordinate generation.C
C            0 ... Default (15 degrees, applied to all 3 tests).
C           -N ... Threshold N degrees, applied to the angle only.
C            N ... Threshold N degrees, applied to all 3 tests.
C
C2IOp(155)
C     IOp(155) ... Number of steps to take in guessing a TS during QST2:
C            0 ... Default (10).
C            N ... Divide the overall step into N increments.
C
C2IOp(156)
C     IOp(156) ... Automatic generation of internal coordinates.
C            0 ... Default (1).
C            1 ... Generate bonds, angles, and dihedrals.
C            2 ... Generate bonds and angles but no dihedrals.
C            3 ... Generate bonds only.
C            4 ... Generate no coordinates automatically.
C
C2IOp(157)
C     IOp(157) ... Maximum step when going down an eigenvector:
C            0 ... Default (0.6)
C            N ... N/1000.
C
C2IOp(158)
C     IOp(158) ... PCM defaults for L124.  Details in link 301.
C
C2IOp(159)
C     IOp(159) ... Default maximum step in Cartesians during red. coord
C                  optimizations.
C           -1 ... Unlimited (10^6 Bohr).
C            0 ... Default (3 Bohr).
C            N ... N/10 Bohr.
C
C2IOp(160)
C     IOp(160) ... Type of fragment calculation:
C            0 ... Default (11).
C            1 ... Force use of new ghost atoms in counterpoise.
C            2 ... Force use of old ghost atoms in counterpoise.
C           10 ... Counterpoise.
C           20 ... Fragment guess.
C           30 ... Morokuma energy decomposition.
C           40 ... XPol.
C           50 ... Variational XPol.
C           60 ... For testing methods within normal iterations.
C           70 ... NEDA.
C
C2IOp(161)
C     IOp(161) ... Hack to save results at each IRC point in L123.
C            0 ... Default (No).
C            1 ... Yes.
C
C2IOp(162)
C     IOp(162) ... Frequency of  analytic Hessians during IRC corrector cycles.
C
C2IOp(163)
C     IOp(163) ... Copy of external input section from chk file:
C            0 ... Default, copy if Geom=AllCheck
C            1 ... Copy regardless.
C            2 ... Do not copy.
C
C2IOp(164)
C     IOp(164) ... Read of atomic pair potential.
C            0 ... Default, copy from chk if Geom=AllCheck.
C            1 ... Read from input.
C            2 ... Do not read.
C            3 ... Read from chk.
C
C     Each pair is given on one line of the form:
C     IA,IB,R0,A,B,ROn,ROff,C(n),C(n-1),...,C(0)
C
C     The potential is of the form P(X)*Exp(A*X)*Exp(B*X*X)*Sw(Rab,ROn,ROff)
C
C     If A is zero, the exponential is omitted.  If B is zero, the Gaussian is
C     omitted, and if ROn is zero, the Stratmann-Scuseria switching function Sw
C     is omitted.  P(X) is the polynomial of degree n given by the specified
C     coefficients.  n is determined from the number of C's provided and can be
C     0 for a constant scale factor.
C
C2IOp(165)
C     IOp(165) ... Convergence of MM microiterations.
C            0 ... Default (10x tighter than macro, except 10x for FOSimult
C                  with MM included).
C            N ... Nx tighter.
C
C2IOp(166)
C     IOp(166) ... Maximum number of vectors in DIIS in L103 (MMNN):
C            0 ... Default.
C           NN ... Save N vectors (default 10 for GEDIIS, 10 for SimOpt).
C         MMNN ... Mix up to MM<NN vectors in DIIS when mixing RFO steps.
C                  (Default NN).
C     A negative value requests uses of the Hessian eigenvector basis for
C     the step.  This is the default and only choice for GEDIIS TS
C     optimizations.
C
C2IOp(167)
C     IOp(167) ... Version of GEDIIS.
C            0 ... Default (1).
C            1 ... Old.
C
C2IOp(168)
C     IOp(168) ... GEDIIS switches:
C           NN ... Switch from RFO to DIIS on RMS force
C                  (10^-NN, default 1.d-3).
C         MM00 ... Switch to En-DIIS from RFO-DIIS on RMS step
C                  (10^-MM, default 2.5d-3).
C      LLL0000 ... Maximum coefficient allowed in RFO-DIIS before
C                  space is reduced (LLL/10, default 10.0).
C   KKK0000000 ... Maximum coefficient allowed in RFO-DIIS before
C                  coefficients are adjusted (KKK/10, default 3.0).
C                  (The minimum value is -KKK/10 + 1).
C
CC1L101
C                       Read geometry input
C
C
C
C      Purpose:  Link 101 reads the title section, Z-matrix section,
C                variables section, and constants section of the input
C                deck.  The data are placed on the r/w file for use
C                by subsequent links in the gaussian system.
C                Optionally, L101 may read coordinates directly or read
C                connectivity information used to create a standard
C                model structure and associated z-matrix and variables.
C      Input:    Data on cards in the input deck and options in IOp.
C      Output:   In atomic units (bohrs/radians).
C                /ZMAT/ CONSISTS OF EIGHT ARRAYS DIMENSIONED TO HOLD THE
C                       DATA FROM UP TO 400 Z-MATRIX CARDS AND TWO
C                       VARIABLES CONTAINING COUNTERS.
C                       IANZ-   ATOMIC NUMBERS
C                       IZ-     CONNECTIVITY DATA.  IZ(I,1)=0 or -1
C                               flags center I as being in cartesian
C                               coordinates, with x,y,z coordinates in
C                               Bl(I), Alpha(I), and Beta(I),
C                               respectively.
C                       BL-     BOND LENGTHS
C                       ALPHA-  VALENCE ANGLES
C                       BETA-   DIHEDRAL OR SECOND VALENCE ANGLE
C                       LBL-    AN ARRAY USED TO MAP FROM THE ARRAY OF
C                               VARIABLE VALUES IN /ZSUBST/ TO THE ARRAY
C                               BL IN /ZMAT/.  FOR LBL(N)=0, NO
C                               SUBSTITUTION IS REQUIRED TO GET BL(N)
C                               IT ALREADY CONTAINS THE PROPER VALUE.
C                               FOR LBL(N)=I, BL(N)=VALUES(I); FOR
C                               LBL(N)=-I, BL(N)=-VALUES(I).
C                      LALPHA-  ANALAGOUS TO LBL FOR ALPHA.
C                      LBETA-   ANALAGOUS TO LBL FOR BETA.
C                      NZ-      THE NUMBER OF CARDS IN THE Z-MATRIX.
C                      NVAR-    THE NUMBER OF VARIABLES, EQUAL TO THE
C                               NUMBER OF SYMBOLS DEFINED IN THE
C                               VARIABLES SECTION OF THE INPUT.
C                /ZSUBST/ CONTAINS THE DATA IN THE VARIABLES SECTION OF
C                         THE INPUT.
C                      ANAMES-  ALPHANUMERIC NAMES OF THE VARIABLES.
C                               (LIMITED TO 8 CHARACTERS IN THE CURRENT
C                                VAX VERSION).
C                      VALUES-  THE CORRESPONDING NUMERIC VALUES.  THESE
C                               WILL BE ALTERED IN THE COURSE OF
C                               GEOMETRY OPTIMIZATIONS AND POTENTIAL
C                               SURFACE SCANS.
C                      INTVEC-  AN ARRAY OF THE INTEGER VALUES FOLLOWING
C                               THE SYMBOL AND VALUE ON A LINE IN THE
C                               VARIABLES SECTION.
C                      FPVEC-   A CORRESPONDING ARRAY CONTAINING THE
C                               FLOATING POINT NUMBER.  THE USE OF THESE
C                               TWO VECTORS DEPENDS UPON THE ROUTE.
C
C
C1
C     A DESCRIPTION OF THE MEANS OF SPECIFYING THE NUCLEAR COORDINATES
C     FOLLOWS.
C
C1GEOMETRY
C2Z-MATRIX
C
C      Z-MATRIX SPECIFICATION SECTION.
C
C          THIS SECTION IS  ALWAYS  REQUIRED.   IT  SPECIFIES  THE
C     NUCLEAR  POSITIONS  AND THE NUMBER OF ELECTRONS OF ALPHA AND
C     BETA SPIN.  THE INPUT IS FREE-FIELD;  THE SEVERAL  ITEMS  ON
C     EACH CARD MAY BE SEPARATED BY EITHER BLANKS OR COMMAS.
C
C          THE  FIRST  CARD  OF  THE  SECTION  SPECIFIES  THE  NET
C     ELECTRIC  CHARGE  (SIGNED INTEGER) AND THE SPIN MULTIPLICITY
C     (POSITIVE INTEGER).  THUS,  FOR  A  NEUTRAL  MOLECULE  IN  A
C     SINGLET STATE, THE ENTRY "0 1" IS APPROPRIATE.
C
C          THE REMAINING CARDS ARE USED TO  SPECIFY  THE  RELATIVE
C     POSITIONS OF THE NUCLEI.  MOST OF THESE WILL BE REAL NUCLEI,
C     USED LATER IN THE MOLECULAR ORBITAL  COMPUTATION.   HOWEVER,
C     IT  IS  FREQUENTLY  USEFUL TO INTRODUCE "DUMMY NUCLEI" WHICH
C     HELP SPECIFY THE  GEOMETRY  BUT  ARE  IGNORED  SUBSEQUENTLY.
C     THEIR USE WILL BECOME CLEAR IN EXAMPLES GIVEN BELOW.
C
C          EACH   NUCLEUS   (INCLUDING   DUMMIES)   IS    NUMBERED
C     SEQUENTIALLY  AND SPECIFIED ON A SINGLE CARD.  (THIS DATA IS
C     REFERRED TO AS THE Z-MATRIX).  THUS, THE NATURE AND LOCATION
C     OF  THE N-TH NUCLEUS IS SPECIFIED ON THE (N+1)TH CARD IN THE
C     SECTION IN TERMS OF THE POSITIONS OF THE PREVIOUSLY  DEFINED
C     NUCLEI 1,2,...(N-1).
C
C          THE INFORMATION ABOUT THE N-TH NUCLEUS IS CONTAINED  IN
C     UP TO EIGHT SEPARATED ITEMS ON THE CARD in one of the forms:
C
C                ELEMENT, N1, LENGTH, N2, ANGLE, N3 TWIST, J
C                Element, 0, x-coordinate, y-coordinate, z-coordinate.
C                Element,    x-coordinate, y-coordinate, z-coordinate.
C                Element, -1, x-coordinate, y-coordinate, z-coordinate.
C
C          "ELEMENT" SPECIFIES THE CHEMICAL NATURE OF THE NUCLEUS.
C     IT  MAY  CONSIST  OF JUST THE CHEMICAL SYMBOL SUCH AS "H" OR
C     "C" FOR OF CARBON.  ALTERNATIVELY, IT MAY BE AN ALPHANUMERIC
C     STRING   BEGINNING   WITH   THE  CHEMICAL  SYMBOL,  FOLLOWED
C     IMMEDIATELY BY A SECONDARY IDENTIFYING INTEGER.  THUS,  "C5"
C     CAN  BE  USED TO SPECIFY A CARBON NUCLEUS, IDENTIFIED AS THE
C     FIFTH CARBON IN THE MOLECULE.  THIS IS SOMETIMES  CONVENIENT
C     IN  FOLLOWING CONVENTIONAL CHEMICAL NUMBERING.  DUMMY NUCLEI
C     ARE DENOTED BY THE SYMBOLS "X" OR "-".  THE  ITEM  "ELEMENT"
C     IS  REQUIRED  FOR  EVERY  NUCLEUS.   FOR  THE  FIRST NUCLEUS
C     SPECIFIED (N=1), it may be the only element on the card,
C     causing it to be placed at the origin, or its cartesian coordinates
C     may be specified.  Ghost atoms (centers with no charge and, by
C     default, no basis functions) are specified by "Bq".  Note that
C     either the MASSAGE route option or a general basis is required to
C     place basis functions on ghost atoms, and that no molecular
C     symmetry is used when Bq atoms are present.  Alternatively, an
C     atom with mechanics type Bq (i.e., "O-Bq") is setup as a ghost
C     of the corresponding atom, with its normal basis functions and
C     XC grid points but no nuclear charge or electrons.
C
C     The most general element specification is
C       SS-MM-CC(prop1=value,prop2=value,...)
C     where SS is the symbol, MM is the optional mechanics type string,
C     CC is the mechanics charge, and the optional properties can include
C     Isotope=X, where X is either integer or floating point.
C
C          "N1" SPECIFIES THE  (PREVIOUSLY  DEFINED)  NUCLEUS  FOR
C     WHICH  THE  INTERNUCLEAR LENGTH R(N,N1) WILL BE GIVEN.  THIS
C     ITEM MAY BE EITHER AN INTEGER (THE VALUE OF N1 <  N)  OR  AN
C     ALPHANUMERIC  STRING.   IN  THE LATTER CASE, THE STRING MUST
C     MATCH THE "ELEMENT" FIELD OF A PREVIOUS Z-MATRIX CARD.  If
C     N1 is 0 or omitted, the cartesian coordinates for this center are
C     expected on the same card.  If N1 is -1, cartesian coordinates
C     are also expected, and in addition this center is flagged as
C     "frozen" to the derivative links, causing this center to be
C     ignored in gradient and frequency calculations.
C
C          "LENGTH" IS THE INTERNUCLEAR LENGTH R(N,N1).  THIS  MAY
C     BE EITHER A POSITIVE FLOATING POINT NUMBER GIVING THE LENGTH
C     IN  ANGSTROMS  OR  AN  ALPHANUMERIC   STRING   (MAXIMUM   16
C     CHARACTERS).   IN THE LATTER CASE, THE LENGTH IS REPRESENTED
C     BY A "VARIABLE" FOR WHICH  A  VALUE  WILL  BE  SPECIFIED  IN
C     SECTION 4.  USE OF VARIABLES IN THE Z-MATRIX IS ESSENTIAL IF
C     OPTIMIZATION IS TO BE CARRIED OUT.  HOWEVER, THEY  CAN  ALSO
C     BE  USED  IN SINGLE-POINT RUNS.  THE ITEMS "N1" AND "LENGTH"
C     ARE REQUIRED FOR ALL NUCLEI AFTER THE FIRST.  FOR THE SECOND
C     NUCLEUS, ONLY "ELEMENT", "N1", AND "LENGTH" ARE REQUIRED.
C
C          "N2" SPECIFIES THE NUCLEUS FOR WHICH  THE  INTERNUCLEAR
C     ANGLE  THETA(N,N1,N2)  WILL  BE GIVEN.  AGAIN THIS MAY BE AN
C     INTEGER (THE VALUE OF N2 <  N)  OR  AN  ALPHANUMERIC  STRING
C     WHICH  MATCHES  A  PREVIOUS "ELEMENT" ENTRY.  NOTE THAT "N1"
C     AND "N2" MUST REPRESENT DIFFERENT NUCLEI.
C
C          "ANGLE" IS THE INTERNUCLEAR ANGLE THETA(N,N1,N2).  THIS
C     MAY  BE  A FLOATING POINT NUMBER GIVING THE ANGLE IN DEGREES
C     OR AN ALPHANUMERIC STRING REPRESENTING A VARIABLE.  "N2" AND
C     "ANGLE"  ARE  REQUIRED FOR ALL NUCLEI AFTER THE SECOND.  FOR
C     THE THIRD NUCLEUS, ONLY "ELEMENT", "N1", "LENGTH", "N2", AND
C     "ANGLE" ARE RQUIRED.
C
C          "N3".  THE SIGNIFICANCE OF "N3" AND "TWIST" DEPENDS  ON
C     THE VALUE OF THE LAST ITEM "J".  IF J=0, OR IS OMITTED, "N3"
C     SPECIFIES TH NUCLEUS FOR  WHICH  THE  INTERNUCLEAR  DIHEDRAL
C     ANGLE,  PHI(N,N1,N2,N3)  WILL  BE  GIVEN;   AS WITH "N1" AND
C     "N2", THIS  MAY  BE  EITHER  AN  INTEGER  (N3  <  N)  OR  AN
C     ALPHANUMERIC STRING MATCHING A PREVIOUS "ELEMENT" ENTRY.
C
C          "TWIST" (IF J=0) IS  THE  INTERNUCLEAR  DIHEDRAL  ANGLE
C     PHI(N,N1,N2,N3).  AGAIN, THIS MAY BE A FLOATING POINT NUMBER
C     GIVING THE  ANGLE  IN  DEGREES  OR  AN  ALPHANUMERIC  STRING
C     REPRESENTING  A  VARIABLE  (OR  A  VARIABLE  PRECEEDED  BY A
C     NEGATIVE SIGN).  THE DIHEDRAL ANGLE IS DEFINED AS THE  ANGLE
C     (-180.0  <  PHI  <= +180.0) BETWEEN THE PLANES (N,N1,N2) AND
C     (N1,N2,N3).  THE SIGN IS POSITIVE IF  THE  MOVEMENT  OF  THE
C     DIRECTED   VECTOR   (N1-->N)  TOWARDS  THE  DIRECTED  VECTOR
C     (N2-->N3) INVOLVES A RIGHTHANDED SCREW MOTION.
C
C          "J".  THE ABOVE DESCRIPTIONS OF "N3" AND "TWIST"  APPLY
C     IF  THE  ITEM "J" IS ZERO OR ABSENT.  ALTHOUGH IT IS USUALLY
C     POSSIBLE TO SPECIFY THE NUCLEUS N BY A BOND LENGTH,  A  BOND
C     ANGLE,  AND  A  DIHEDRAL  ANGLE,  IT IS SOMETIMES SIMPLER TO
C     REPLACE TH DIHEDRAL ANGLE BY A SECOND BOND ANGLE.   THIS  IS
C     CALLED FOR BY USING J = +1 OR -1.
C
C          "N3".  IF "J" IS +1 OR -1, "N3" SPECIFIES  THE  NUCLEUS
C     FOR WHICH THE SECOND INTERNUCLEAR ANGLE CHI(N,N1,N3) WILL BE
C     GIVEN.  AS USUAL, THIS MAY BE EITHER AN INTEGER (N3 < N)  OR
C     AN  ALPHANUMERIC  STRING  REPRESENTING  A PREVIOUSLY DEFINED
C     NUCLEUS.
C
C          "TWIST".  IF "J" IS +1 OR -1, THEN THIS ITEM GIVES  THE
C     VALUE  FOR  THE  SECOND INTERNUCLEAR ANGLE CHI(N,N1,N3).  AS
C     BEFORE, THIS MAY BE EITHER A FLOATING POINT NUMBER (VALUE IN
C     DEGREES) OR AN ALPHANUMERIC STRING REPRESENTING A VARIABLE.
C
C          "J".  IN THE EVENT OF SPECIFICATION BY TWO INTERNUCLEAR
C     ANGLES  THETA, CHI, THERE WILL BE TWO POSSIBLE POSITIONS FOR
C     THE NUCLEUS N.  THIS IS FIXED BY THE SIGN OF "J".  THUS J=+1
C     IF THE TRIPLE VECTOR PRODUCT
C
C          (N1-->N) .  ((N1-->N2) X (N1-->N3))  IS  POSITIVE,  AND
C     J=-1 IF THE PRODUCT IS NEGATIVE.
C
C          THE Z-MATRIX IS TERMINATED  BY  THE  BLANK  CARD  WHICH
C     INDICATES THE END OF THE MOLECULE SPECIFICATION SECTION.  IF
C     NO VARIABLES HAVE BEEN INTRODUCED AND IF THE COMMAND "ALTER"
C     HAS  NOT  BEEN INVOKED IN THE JOB-TYPE SECTION, THE INPUT IS
C     COMPLETE  AND  GAUSSIAN  82  WILL  PERFORM   THE   REQUESTED
C     COMPUTATION.
C
C2ModelBuilder
C
C     Link 101 can also accepted connectivity information and use
C     it to construct a standard model geometry and symbolic
C     z-matrix, which can then be used in a geometry optimization.
C     The two available models (A and B) differ in that model A
C     takes into account the type (single, double, triple, etc.) of
C     a bond in assigning bond lengths, while model B bond lengths
C     depend only on the types of the atoms involved.  Model A is
C     only available for Li, H, C, N, O, and F.  Model B is available
C     for all atoms from H to Cl except He and Ne.
C
C     The basic input to the model builder is a called a "Short
C     Formula Matrix."  This is a collection of cards each of which
C     defines an atom (by its atomic symbol) and its connectivity,
C     by up to six more entries.  Each of these can be either
C     an integer, which is the number of the card defining another
C     explicitly specified atom to which the current atom is bonded,
C     or an atomic symbol (e.g. H, F) to which the current atom is
C     connected by a terminal bond, or a symbol for a terminal
C     functional group which is bonded to the current atom.  The
C     functional groups currently available are OH, NH2, Me, Et,
C     NPr, IPr, NBu, IBu, and TBu.
C
C     The short formula matrix also implicitly defines the rotational
C     geometry about each bond in the following manner.  Suppose atoms
C     X and Y are explicity specified.  Then X will appear in row Y
C     and Y will appear in row X.  Let I be the atom to the right of X
C     in row Y and J be the atom to the right of Y in row X.  Then atoms
C     I and J are put in the trans orientation about the X-Y bond.
C
C     The short formula matrix may be followed by optional cards
C     modifying the generated structure.  There are zero or more
C     of each of the following cards, which must be grouped together
C     in the order given here:
C
C     ATOMGEOM,I,Geom --
C     Normally the local geometry about an atom is defined by the
C     number and types of bond about the atom (e.g., carbon in
C     methane is tetrahedral, in ethylene is trigonal).  All bond
C     angles at one center must be are equal.  The AtomGeom card
C     changes the value of the bonds at center I.  Geom may be
C     the angle as a floating point number, or one of the strings
C     TETR, PYRA, TRIG, BENT, or LINE.
C
C     BONDROT,I,J,K,L,Geom --
C     This changes the orientations of the I-J and K-L bonds about
C     the J-K bond.   Geom is either the dihedral angle or one of
C     the strings CIS (=>0), TRANS (=>180), GAUP (=>+60), or GAUM
C     (=>-60).
C
C     BONDLEN,I,J,NewLen --
C     This sets the length of the I-J bond to NewLen (a floating
C     point value).
C
C     The model builder can only build structures with atoms in
C     their normal valencies.  If a radical is desired, its extra
C     valence can be "tied down" using dummy atoms, which are
C     specified by a minus sign before the atomic symbol (e.g., -H).
C     Only terminal atoms can be dummies.
C
C     This section is still experimental, so there are no guarrantees,
C     but feel free to direct inquires and comments to Mike Frisch.
C?
      Common /IO/     In, IOut, IPunch
      Common /MUnit/  IUnit(20)
      Integer MaxIOp
      Parameter (MaxIOp=DEFMAXIOP)
      Common /IOp/    IOp(MaxIOp)
      Integer FilNum
#include "commonmol.inc"
#include "commonz.inc"
      Common /PhyCon/ PhyCon(30)
      Equivalence (ToAng,PhyCon(1))
      Parameter (IOGen=501,IOCon=503,IRwSym=551,IOZRed=665,IOZSav=669,
     $  IOSLEq=671,IRwOni=674,IRwMM=687,IOFDPr=721,IRwCP=730,IRwCrd=775,
     $  IRwInp=786,IRwPPT=787,MaxStr=4,MaxRed=50,MaxRAD=50000,MxLONI=4,
     $  MaxCon=100,MaxFrq=1000,MinEl=-2,MaxEl=200,NMMVB=9,LChg=7,LDBF=8,
     $  MxRdTp=DEFMAXREDTYPE,MxRdIn=DEFMAXREDINDEX,NumNBC=8)
      Dimension V(*), NTRedS(4,MaxStr), IIniNR(4), IElDat(MinEl:MaxEl),
     $  JJ(1), IMMVB(NMMVB)
      Logical ReadFz, ConsFz, NoFrez, FixTet, SaveZ, AllChk, CkRigF,
     $  DoRed, DoAddB, RdZRed, CartZm, AlCart, DoCopy, DoRRS, NewCrd,
     $  DoReac, DoProd, Large, DoQST, DoPath, InitNR, DoDefI, HavRed,
     $  AllCor, DelMod, GetMod, NoFlag, DoDef, RdChk, KeepMM, KeepM1,
     $  ScanDf, RdFreq, HavCon, DoRedI, RdMMP, PosMul, ForcW0, ONIOM,
     $  RdExtC, PrtChg, DoRdFz, ChkDVL, TitInp, Got1Ck, SrtFrq, AllFrz
      Character*80 Mess01, Mess02
      Equivalence (InitNR,IIniNR(1))
      Save IIniNR, One, JJ, IHavMO, NACore, NAVal2, NAVal1, NBCore,
     $  NBVal2, NBVal1
      Data IIniNR/4*0/, One/1.0d0/, JJ/0/, IHavMO/0/, NACore/0/,
     $  NAVal2/0/, NAVal1/0/, NBCore/0/, NBVal2/0/, NBVal1/0/
 1000 Format(' Fast transformation storage of',I4,' MWords has been',
     $       ' initialized.')
 1030 Format(' Limit of',I7,' fragments for counterpoise exceeded.')
 1040 Format(' Recovering MM parameters from chk file.')
 1050 Format(' Generating MM parameters.')
 1060 Format(' Retrieving read-in MM parameters from chk file.')
 1070 Format(' No MM parameters found on chk file.')
 1080 Format(' Background charge distribution read from input stream:')
 1090 Format(' Background charge distribution read from chk.')
 1100 Format(' Unrecognized background charge reading option.')
 2000 Format(1x,'Initializing the coordinate definitions file.')
 2002 Format(4x,'ICrDer=',I5,' ICrDef=',I5,/,
     $  4x,'ICrDf0=',I5,' ICrDf1=',I5,' ICrDf2=',I5)
 2010 Format(1x,'Forming coordinate definitions.')
 2020 Format(1x,'Forming B-matrix definitions.')
 2030 Format(1x,'Forming B-matrix derivative definitions.')
 2040 Format(1x,'Writing the coordinate definitions file.')
 2900 Format(1x,'RdGeom: ICrDef=',I10)
 3000 Format(1x,'NCoord=',I10,' NTRed=',I10)
C
      Call Drum(V,MDV)
      Call FillEl(MinEl,MaxEl,IElDat)
      NWChek = Max(4,IntPWP(0))
C     Check general values for MxRdTp, MxRdIn
      If(Mod(MxRdTp+1,NWChek).ne.0)
     $  Call GauErr('Illegal MxRdTp in RdGeom.')
      If(Mod(MxRdIn,NWChek).ne.0)
     $  Call GauErr('Illegal MxRdIn in RdGeom.')
C     Now check for current values since most code not general.
      If(MxRdTp.ne.3.or.MxRdIn.ne.4)
     $  Call GauErr('Code is not yet general for MxRdTp or MxRdIn.')
      ICType = Mod(IOp(18),100)/10
      CkRigF = ICType.ne.0
      DoRedI = ICType.eq.2.or.ICType.eq.3.or.ICType.eq.5
      DoRed = DoRedI
      DoAddB = (Mod(IOp(18),1000)/100).gt.0
      IHBond = Mod(IOp(18),10000)/1000
      If(IHBond.eq.0) IHBond = 1
      IAlCor = Mod(IOp(18),100000)/10000
      InpUnt = IOp(20)
      Call ILSW(1,14,InpUnt)
      FixTet = IOp(24).eq.0.or.IOp(24).eq.1
      NStruc = Min(Max(Mod(IOp(27),10),1),MaxStr)
      IWeiF = Mod(IOp(27),100) / 10
      NTrRot = IOp(28)
      If(NTrRot.eq.0) NTrRot = -1
      If(NTrRot.eq.-2) NTrRot = 0
      IDoPun = IOp(32)
      IPrint = IOp(33)
      IDump = IOp(34)
      Call FFSet(IDump)
      NPath = IOp(42)
      I49 = IOp(49)
      DoPath = NPath.ne.0.and.DoRed.and.NStruc.gt.1
      Large  = ICType.eq.5.and.NPath.eq.0.and.DoRed
      AllCor = IAlCor.gt.1.or.(IAlcor.eq.0.and..not.Large)
      MicOTy = Mod(IOp(18),10000000) / 1000000
      If(MicOTy.eq.0) MicOTy = 2
      DoRed = DoRed.and.MicOTy.ne.4
      KeepMM = MicOTy.eq.2
      If(ITqry(IOSLEq).gt.0) Call FileIO(5,IOSLEq,0,V,0)
      If(DoPath) then
        DoReac = Mod(I49,10).le.1
        DoProd = (Mod(I49,100)/10).le.1
        NTot = NPath
        If(DoReac) NTot = NTot - 1
        If(DoProd) NTot = NTot - 1
        If(NStruc.gt.2) NTot = (NPath-1)/2
        WeiFac = GFloat(NTot)/(GFloat(NTot+1))
      else
        If(IWeiF.eq.0) IWeiF = 5
        WeiFac = GFloat(IWeiF) / GFloat(10)
        endIf
      IMolFl = IOp(123)
      Call IniMol(1,IMolFl)
      Call IniZMt(.True.)
      IGoI = Mod(IOp(29),10)
      IGo = IGoI
      AllChk = IGoI.eq.7
      RdChk = IGoI.eq.2.or.IGoI.eq.3.or.IGoI.eq.6.or.AllChk
      IRead = Mod(IOp(29),10000000)/1000000
      If(IRead.eq.9) IRead = 0
      If(RdChk.and.NStruc.gt.1) then
        Mess01 = 'Only one structure can be retrieved from the ' //
     $    'checkpoint file.'
        Call GauEr2(0,Mess01,' ',' ','route')
      else if(RdChk.and.IRead.gt.0) then
        Mess01 = 'Geom=Checkpoint and Geom=AllCheck can not be '
     $    // 'combined with Geom=Input'
        Mess02 = 'or Geom=AllInput.'
        Call GauEr2(0,Mess01,Mess02,' ','route')
        endIf
      IPtMod = Mod(IOp(29),100)/10
      IAbort = Mod(IOp(29),1000)/100
      If(IAbort.eq.0) IAbort = 1
      IFlags = Mod(IOp(29),10000)/1000
      ReadFz = IFlags.eq.1
      NoFrez = IFlags.eq.2.or.IFlags.eq.5
      NoFlag = IFlags.eq.3.or.IFlags.eq.5
      DoDef  = IFlags.eq.4.or.IFlags.eq.6
      DoRdFz = IFlags.eq.6
      AllFrz = IFlags.eq.7
      ConsFz = (Mod(IOp(29),100000)/10000).le.1
      IGenCZ = Mod(IOp(29),1000000)/100000
      IStep = IOp(41)
      NAlRed = IOp(47)
      NuOfLa = Mod(IOp(52),10)
      If(NuOfLa.eq.0) NuOfLa = 1
      NFrag = IOp(62)
      PosMul = (Mod(IOp(160),100)/10).ne.2
      InpFF = IOp(56)
      IDoCon = IOp(57)
      IDCon0 = Mod(IDoCon,10)
      IDCon2 = Mod(IDoCon,1000)/100
      IAprBG = Mod(IOp(59),100000)/10
      ScanDf = IOp(60).eq.1
      MOMMFl = Mod(IOp(61),10)
      IParIn = IOp(64)
      IParI2 = Mod(IOp(64),1000)/100
      ICntMM = IOp(65)
      IQEq   = IOp(66)
      IQEq0  = Mod(IQEq,10)
      If(IQEq0.eq.0.and.IQEq.ne.0) IQEq0 = 1
      IGenMM = IOp(67)
      IRwMMC = FilNum(IRwMM,IUnit(9))
      LRwMMC = Max(ITqry(IRwMMC),0)
      If(IGenMM.eq.0) then
        If(RdChk.and.IParI2.ne.0) then
          IGenMM = 3
        else if(RdChk.and.LRwMMC.gt.0) then
          IGenMM = 2
        else
          IGenMM = 1
          endIf
        endIf
      If(IGenMM.eq.3.and.IParI2.eq.0) IParIn = IParIn + 200
      IParI2 = Mod(IParIn,1000)/100
      ICpLCh = Mod(IOp(142),10)
      If(ICpLCh.eq.0) then
        If(IQEq0.eq.1) then
          ICpLCh = 3
        else
          ICpLCh = 1
          endIf
        endIf
      RdFreq = Mod(IOp(83),10).eq.1
      ForcW0 = (Mod(IOp(83),100)/10).ne.2
      SrtFrq = (Mod(IOp(83),1000)/100).ne.2
      Isotop = IOp(88)
      If(Isotop.eq.0) then
        If(RdChk) then
          Isotop = 4
        else
          Isotop = 1
          endIf
        endIf
      IThTem = IOp(112)
      IThPre = IOp(113)
      IThScl = IOp(114)
      Call IOpMVB(IOp(117),IOp(119),NMMVB,IMMVB)
      IDoSec = IOp(122)
      IBkgCh = IOp(125)
      ChkDVL = Mod(IOp(132),10).ne.2
      ISwitc = IOp(135)
      ICutV = IOp(136)
      ICutC = IOp(137)
      IRanCt = IOp(138)
      DoRRS = IOp(140).ne.2
      IRRAn = IOp(144)
      Call IClear(4*MaxStr,NTRedS)
C
C     Check to see if we are supposed to use the new internal coordinate
C     code. The logical flag NewCrd indicates whether or not this code is
C     to be used within this link. The integer flag ICrDef indicates which
C     options are used when building the coordinate definitions, etc.
C
      ICrOpt = IOp(145)
      IRIC = Mod(ICrOpt,10)
      If(IRIC.eq.0) IRIC = 1
      ICrDer = Mod(ICrOpt,100)/10
      If(ICrDer.eq.0) ICrDer = 3
      ICrDef = Mod(ICrOpt,1000)/100
      If(ICrDef.eq.0) ICrDef = 1
      ICrMod = Mod(ICrOpt,10000)/1000
      NewCrd = IRIC.eq.2
      If(NewCrd.and.ICrDef.eq.1) then
        ICrDf0 = 10030
        ICrDf1 = 20010
        ICrDf2 = 30010
      else if(NewCrd.and.ICrDef.eq.2) then
        ICrDf0 = 10202
        ICrDf1 = 20010
        ICrDf2 = 30010
      else if(NewCrd.and.ICrDef.eq.3) then
        ICrDf0 = 10202
        ICrDf1 = 20902
        ICrDf2 = 30902
      else if(NewCrd.and.ICrDef.eq.4) then
        ICrDf0 = 10202
        ICrDf1 = 20902
        ICrDf2 = 30802
      else if(NewCrd.and.ICrDef.eq.5) then
        ICrDf0 = 10202
        ICrDf1 = 20802
        ICrDf2 = 30802
      else if(NewCrd.and.ICrDef.eq.6) then
        ICrDf0 = 10402
        ICrDf1 = 20010
        ICrDf2 = 30010
      else if(NewCrd.and.ICrDef.eq.7) then
        ICrDf0 = 10021
        ICrDf1 = 20010
        ICrDf2 = 30010
      else if(NewCrd) then
        Write(IOut,2900) ICrDef
        Call GauErr('Invalid ICrDef in RdGeom.')
      else
        ICrDf0 = 10000
        ICrDf1 = 20000
        ICrDf2 = 30000
        endIf
      If(ICrMod.gt.0) ICrDf0 = ICrDf0 + ICrMod*1000
      LinTst = IOp(154)
      NStpTS = IOp(155)
      IDefTp = IOp(156)
      If(IDefTp.eq.0) IDefTp = 1
      ICpInp = IOp(163)
      IRdPPT = IOp(164)
C
C     If the checkpoint file is from an ONIOM run, restore the real system
C     there.  This avoids problems with Geom=Check, etc. after an ONIOM job
C     step which has failed in the middle.
C
      If(DoRRS) Call RestRS(IOut,IPrint,2,2,.False.,.True.,V,MDV)
C
C     If MM parameters are to be read from the chk file, pick up the symbol
C     table information for non-standard symbols first, so it can be used
C     in processing the input structure.
C
      LenPar = 0
      RdMMP = IGenMM.eq.3.and.Mod(IParIn,100).gt.0.and.LRwMMC.gt.0
      If(RdMMP) then
        Call CnvMMF(1,IRwMMC,LenPar,V,MDV)
        IPOld = 1
        IV = IPOld + LenPar
        Call TstCor(IV,MDV,'RdGeom-RdSymT')
        Call CnvMMF(0,IRwMMC,LenPar,V(IPOld),MDV)
        Call InOldP(LenPar,V(IPOld),NCOld,IParCt,IndIC,IndRC,NDSymO,
     $    NSOld,IndSTO,IStepD,NC,IIList)
        Call CpMMST(IOut,IPrint,NDSymO,NSOld,V(IPOld+IndSTO),NumTpS,
     $    MaxTpS,IAtTpS)
        endIf
C
C     Allocate space for sets of redundant coordinates.
C
      IOZRC = FilNum(IOZRed,IUnit(9))
      RdZRed = (IGoI.eq.2.or.IGoI.eq.6.or.IGoI.eq.7)
     $  .and.ITqry(IOZRC).gt.0
      DoRed = DoRed.or.DoRdFz.or.(RdZRed.and.(ICType.eq.0.or.DoRedI))
      If(NStruc.gt.1.and..not.DoRed) Call GauErr(
     $  'More than one input structure requires red. int. coordinates.')
      If(NStruc.eq.3.and.DoPath) then
        NStrAl = 4
      else if(NStruc.ge.2) then
        NStrAl = 3
      else
        NStrAl = 2
        endIf
      IPWP = IntPWP(0)
      If(DoRed.or.RdZRed.or.MicOTy.eq.4.or.IRead.gt.0) then
        If(NAlRed.eq.0) then
          MaxRA = Max(MaxRAD,MDV/(100*NStrAl))
          NAlRed = ((Min(MaxRed*MaxAtm,MaxRA)+IPWP-1)/IPWP)*IPWP
          endIf
      else
        NAlRed = IPWP
        endIf
      NAlSav = Max(MaxAtm,NAlRed)
      MxFrag = Max(MaxAtm,MxLONI**2)
      NDChMu = MxFrag + 1
      If(NFrag.gt.MxFrag) then
        Write(IOut,1030) MxFrag
        Call Lnk1E(0)
        endIf
      MNDNBC = 2*MxFrag + 8
      IIChOn = 1
      IMulOn = IIChOn + InToWP(NDChMu)
      IIChFr = IMulOn + InToWP(NDChMu)
      IMulFr = IIChFr + InToWP(NDChMu)
      IIChIn = IMulFr + InToWP(NDChMu)
      IMulIn = IIChIn + InToWP(2*NDChMu)
      IHvLay = IMulIn + InToWP(2*NDChMu)
      INBCom = IHvLay + InToWP(NDChMu)
      IIZRed = INBCom + InToWP(NumNBC*MNDNBC)
      IValue = IIZRed + InToWP(4*NAlRed*NStrAl)
      ICSave = IValue + NAlRed*Max(NStrAl*2,2+NStrAl)
      IIAnSv = ICSave + 3*NAlSav*NStrAl
      IIntVc = IIAnSv + InToWP(NAlSav*NStrAl)
      IFPVec = IIntVc + NStrAl*InToWP(NAlRed)
      ITRead = IFPVec + NAlRed*NStrAl
      IZTemp = ITRead + InToWP(MaxNZ)
      IOptAt = IZTemp + InToWP(MaxNZ*4)
      IKeep = IOptAt + InToWP(MaxNZ)
      IV0 = IKeep + InToWP(MaxZSb)
      Call TstCor(IV0,MDV,'RdGeom-1')
      Call IClear(NAlRed*NStrAl,V(IIntVc))
      Call AClear(NAlRed,V(IFPVec))
      Call AClear(NAlRed*NStrAl*2,V(IValue))
      Call IClear(MaxNZ,V(ITRead))
      NDNBsC = 0
      Call IClear(MNDNBC,V(INBCom))
      Call IClear(NDChMu,V(IIChOn))
      Call IClear(NDChMu,V(IMulOn))
      Call IClear(NDChMu,V(IIChFr))
      Call IClear(NDChMu,V(IMulFr))
      Call IClear(2*NDChMu,V(IIChIn))
      Call IClear(2*NDChMu,V(IMulIn))
      Call AllZAt('RdGeom',.True.,MaxNZ,MaxZSb,IV0,IIAtWZ,IAtmWZ,IAtSpZ,
     $  IAtZEZ,IAtQMZ,IAtGFZ,IAtFrZ,IAtRsZ,ISymbl,IAtPDZ,IAtZNZ,INamCt,
     $  Next,V,MDV)
      Call AllOni('RdGeom',.True.,MxLONI,MaxNZ,Next,ILayer,IIAnMo,
     $  IIAtMo,ILiAto,IACMod,IDista,IV0,V,MDV)
      IV = IV0
      MDV1 = MDV - IV + 1
C
C     Read in the requested number of structures.
C
      IScr2 = IV
      IVX = IScr2 + InToWP(NAlRed)
      MDVX = MDV - IVX + 1
      TitInp = .False.
      Got1Ck = .False.
      Do 100 IStruc = 1, NStruc
        LastOf = (IStruc-1)*NAlRed
        LastIZ = InToWP(4)*LastOf
        IIntVI = IIntVc + InToWP(LastOf)
        Call RdStrC(In,IOut,IPunch,IPrint,IDump,IPtMod,IAbort,IRead,
     $    IGoI,IDoPun,InpUnt,ScanDf,PosMul,IStruc,ConsFz,SaveZ,IStep,
     $    V(IKeep),V(IZTemp),DoRed,HavRed,RdChk,IGo,NAlRed,NTRedS,
     $    V(IIZRed),V(IValue),V(IIntVI),V(IFPVec),NAlSav,V(IIAnSv),
     $    V(ICSave),MxLONI,NuOfLa,NFrag,V(IIChOn),V(IMulOn),V(ILayer),
     $    V(IIAnMo),V(IIAtMo),V(IACMod),V(ILiAto),V(IDista),IRwOni,
     $    IRwCP,InpFF,DoAddB,KeepMM,KeepM1,IElDat,IHavMO,MNDNBC,NDNBsC,
     $    V(INBCom),V(IIAtWz),V(IAtmWZ),V(IAtSpZ),V(IAtZEZ),V(IAtQMZ),
     $    V(IAtGFZ),V(IAtFrZ),V(IAtRsZ),V(IAtPDZ),V(IAtZNZ),FixTet,
     $    V(ISymbl),V(INamCt),ICkUnt,TitInp,NDChMu,V(IIChFr),V(IMulFr),
     $    V(IIChIn),V(IMulIn),V(IVX),MDVX)
        Got1Ck = Got1Ck.or.RdChk
        DoDefI = DoDef.or.(DoRed.and.(DoRdFz.or..not.HavRed)).or.
     $    (IGo.ne.2.and.IGo.ne.6.and.IGo.ne.7)
        ICnstr = Mod(IStruc,3)*NAlRed
        IndCSv = ICSave + 3*NAlSav*(IStruc-1)
        If(IGenCZ.gt.0) then
          CartZm = AlCart(0)
          If((NZ.le.0.and..not.SaveZ).or.CartZm) then
            DoCopy = .not.CartZm
            Call ZCart(DoCopy,NAtoms,IAn,IAtTyp,AtChMM,C)
            If(IGenCZ.eq.1) then
              If(ITqry(IOZRed).gt.0) Call FileIO(5,IOZRed,0,V,0)
              If(ITqry(IOZRC).gt.0) Call FileIO(5,IOZRC,0,V,0)
              endIf
            Call SaveZM(In,IOut,.False.,.False.,.False.,.True.,0,
     $        V(IKeep))
            Call RdWrZ(1,-IUnit(9),ANames)
            endIf
          endIf
        IOConI = IOCon
        If(RdChk.or.IDCon0.eq.4) IOConI = FilNum(IOCon,ICkUnt)
        IDConX = IDoCon
        If(IDCon0.eq.0) then
          If(RdChk.and.ITqry(IOConI).gt.0) then
            IDConX = IDConX + 4
          else
            IDConX = IDConX + 3
            endIf
        else
          IDConX = IDoCon
          endIf
        HavCon = Mod(IDConX,10).gt.1
        If(HavCon) then
          MxConU = MaxCon
        else
          MxConU = 1
          endIf
        If(IStruc.eq.1) then
          L1R = MxConU*NAtoms
          L1I = InToWP(L1R)
          INBon1 = IVX
          IIBon1 = INBon1 + InToWP(NAtoms)
          IIBTy1 = IIBon1 + L1I
          IRBTy1 = IIBTy1 + L1I
          INBon2 = IRBTy1 + L1R
          IIBon2 = INBon2 + InToWP(NAtoms)
          IIBTy2 = IIBon2 + L1I
          IRBTy2 = IIBTy2 + L1I
          IVX = IRBTy2 + L1R
          Call TstCor(IVX,MDV,'RdGeom-DoCon')
          MDVX = MDV - IVX + 1
          INBonG = INBon1
          IIBonG = IIBon1
          IIBTyG = IIBTy1
          IRBTyG = IRBTy1
        else
          INBonG = INBon2
          IIBonG = IIBon2
          IIBTyG = IIBTy2
          IRBTyG = IRBTy2
          endIf
        If(HavCon) then
          IMapZA = IVX
          IV1 = IMapZA + InToWP(Max(NAtoms,NZ))
          Call TstCor(IV1,MDV,'RdGeom-MapZAt')
          If(IDCon2.eq.1) then
            NMapZA = NZ
            Call MkMpZA(NZ,IAnZ,NAtoms,V(IMapZA))
          else
            NMapZA = NAtoms
            Call ISeq(NMapZA,0,1,V(IMapZA))
            endIf
          Call DoCon(IOut,IPrint,IDConX,IOConI,IOCon,MxConU,NAtoms,IAn,
     $      C,NMapZA,V(IMapZA),ToAng,MxBond,V(INBonG),V(IIBonG),
     $      V(IIBTyG),V(IRBTyG),V(IV1),MDV-IV1+1)
        else
          Call IClear(NAtoms,V(INBonG))
          endIf
        If(RdChk) then
          If(ICkUnt.ne.IUnit(9)) Call FClose(ICkUnt,0)
        else
          Call FixCou(IOut,NFrag,NAtoms,IAtFrg,MxFrag,V(IHvLay),NZ,IAnZ,
     $      NuOfLa,V(ILayer))
          If(.not.HavCon) Call GauErr('No HavCon for FixOni in RdGeom.')
          Call TstCor(IVX+NZ,MDV,'RdGeom-FixOni')
          Call FixOni(IOut,ChkDVL,MxLONI,NuOfLa,NZ,V(IHvLay),V(IIAnMo),
     $      IAnZ,V(ILiAto),V(ILayer),V(IIAtMo),V(IACMod),MaxNZ,
     $      V(IDista),MxConU,V(INBonG),V(IIBonG),V(IRBTyG),V(IVX))
          endIf
C
C       Generate internal coordinates if needed.
C
        If(DoRed.and.IStruc.le.2) then
          GetMod = DoAddB.or.IGoI.eq.6
          DelMod = .True.
          NARed  = 0
          NABond = 0
          NAAng  = 0
          NADih  = 0
          NTRed  = NTRedS(1,IStruc)
          NTBond = NTRedS(2,IStruc)
          NTAng  = NTRedS(3,IStruc)
          NTDih  = NTRedS(4,IStruc)
          Call GenITR(In,IOut,IPrint,DoDefI,KeepM1,DoRdFz,CkRigF,
     $      IStruc,MaxNZ,NZ,IZ,IAnZ,NuOfLa,V(ILayer),V(ILiAto),
     $      V(IIAnMo),NAtoms,IAn,IAtTyp,IAtRes,NResid,IResNa,IResNu,
     $      V(ITRead),MicOpt,V(IOptAt),V(IVX),MDVX)
          Call RedCor(IOut,IPrint,IDefTp,InpUnt,GetMod,DoDefI,NoFrez,
     $      NoFlag,AllFrz,IHBond,AllCor,DelMod,LinTst,NTrRot,NAtoms,
     $      NAlRed,LenCon,NARed,NABond,NAAng,NADih,NTRed,NTBond,NTAng,
     $      NTDih,V(IndCSv),V(IIAnSv),V(ITRead),V(IIZRed+LastIZ),
     $      V(IIntVI),V(IValue+LastOf),V(IValue+ICnstr),V(IFPVec),ToAng,
     $      HavCon,MxConU,V(INBonG),V(IIBonG),V(IRBTyG),MicOpt,V(IScr2),
     $      V(IVX),MDVX)
          NTRedS(1,IStruc) = NTRed
          NTRedS(2,IStruc) = NTBond
          NTRedS(3,IStruc) = NTAng
          NTRedS(4,IStruc) = NTDih
          Call FixMcO(IOut,IPrint,NAtoms,NTRedS(2,IStruc),
     $      NTRedS(3,IStruc),NTRedS(4,IStruc),V(IIZRed+LastIZ),MicOpt)
        else
          Call GenITR(In,IOut,IPrint,.True.,KeepM1,DoRdFz,CkRigF,IStruc,
     $      MaxNZ,NZ,IZ,IAnZ,NuOfLa,V(ILayer),V(ILiAto),V(IIAnMo),
     $      NAtoms,IAn,IAtTyp,IAtRes,NResid,IResNa,IResNu,V(ITRead),
     $      MicOpt,V(IOptAt),V(IVX),MDVX)
          endIf
        If(HavCon) then
          WtNew = One / GFloat(IStruc)
          WtOld = One - WtNew
          If(IStruc.gt.1.and.IPrint.ge.4) then
            Call PrtCon(IOut,'Bef Merge Con1',.True.,NAtoms,MxConU,IAn,
     $        V(INBon1),V(IIBon1),JJ,V(IRbTy1),0)
            Call PrtCon(IOut,'Bef Merge Con2',.True.,NAtoms,MxConU,IAn,
     $        V(INBon2),V(IIBon2),JJ,V(IRbTy2),0)
            endIf
          If(IStruc.gt.1) Call MrgCon(IOut,IPrint,IOCon,NAtoms,IAn,
     $      MxConU,WtNew,WtOld,V(INBon2),V(IIBon2),V(IRBTy2),V(INBon1),
     $      V(IIBon1),V(IRBTy1))
          If(IStruc.eq.NStruc) then
            If(IPrint.ge.3) Call PrtCon(IOut,'Saved connectivity',
     $        .True.,NAtoms,MxConU,IAn,V(INBon1),V(IIBon1),JJ,V(IRbTy1),
     $        0)
            Call RWCon(1,IOCon,NAtoms,MxConU,MxBond,V(INBon1),V(IIBon1),
     $        V(IRBTy1))
            endIf
          endIf
  100   Continue
      If(IRead.eq.2.and..not.Got1Ck) then
        Mess01 = 'Geom=AllInput was specified, but no checkpoint files'
     $    // ' were read.'
        Mess02 = 'Therefore no title, charge or multiplicity was' //
     $    ' specified.'
        Call GauEr2(0,Mess01,Mess02,' ','molecule specification')
        endIf
      If(RdChk.and.Isotop.eq.1) Isotop = 5
      If(IPrint.gt.0)
     $  Call PrtITM(IOut,' ',NAtoms,V(ITRead),MicOpt)
C
C     For transition state optimizations, merge the internal coordinates
C     and generate the initial guess transition structure in terms of them.
C
      IV = IV0
      MDV1 = MDV - IV + 1
      If(NStruc.eq.1) then
        IStruc = 1
        ICnstr = NAlRed
      else if(NStruc.eq.3.and.DoPath) then
        IStruc = 4
        ICnstr = 0
      else
        IStruc = 3
        ICnstr = 0
        endIf
      LastOf = (IStruc-1)*NAlRed
      LastIZ = InToWP(4)*LastOf
      NTRed = NTRedS(1,IStruc)
      DoQST = NPath.eq.0
      If(NStruc.gt.1.or.DoPath) then
        Call RedTS(IOut,IPrint,NStpTS,IAprBG,NAtoms,NTrRot,NAlRed,
     $    NStruc,NAlSav,V(IIAnSv),V(ICSave),V(IIZRed),NTRedS,V(IValue),
     $    WeiFac,V(IValue+ICnstr),V(IIntVc),DoReac,DoQST,DoPath,ANames,
     $    V(IV),MDV1)
        If(DoAddB.or.IGoI.eq.6) Call AddRed(IOut,NAtoms,
     $    V(IIZRed+LastIZ),ToAng,NTRedS(1,IStruc),NTRedS(2,IStruc),
     $    NTRedS(3,IStruc),NTRedS(4,IStruc),V(IndCSv),MicOpt,
     $    V(IValue+LastOf),V(IIntVc),V(IFPVec),V(IValue+ICnstr),
     $    V(IVX),MDVX)
        NTRed = NTRedS(1,IStruc)
        endIf
      If(DoRed) then
C
C       Set up the structure as Cartesians in /ZMat/ if doing redun. coords.
C       Initialize the storage for the inverse G2, P2 and H2 matrices
C       for large molecular calculation
C
        If(Large) then
          NTRed  = NTRedS(1,IStruc)
          NNT    = (NTRed*(NTRed+1))/2
          NAt3X = 3*ICouQM(NAtoms,MicOpt)
          NAt3TX = (NAt3X*(NAt3X+1))/2
          LenF   = NNT + NAt3TX*2 + 1 + 1 + 20*NAt3X
          LenFW = (LenF+1024**2-1)/(1024**2)
          If(IPSave(0).eq.0) Write(IOut,1000) LenFW
          Call ConDDF(IOSLEq,LenF)
          IXXX   = IV
          IB     = IXXX + NTRed
          IIB    = IB + 12*NTRed
          ITri   = IIB + InToWP(4*NTRed)
          IEnd   = ITri + Max(NAt3TX,NNT)
          Call TstCor(IEnd,MDV,'BigMol')
          Call RedVal(IOut,IPrint,.True.,C,V(IIZRed+LastIZ),V(IXXX),
     $      NTRed,NTRedS(2,IStruc),NTRedS(3,IStruc),NTRedS(4,IStruc),
     $      V(IEnd+1),MDV-IEnd)
          Call FormB(IOut,IPrint,ICouQM(NAtoms,MicOpt),NTRedS(2,IStruc),
     $      NTRedS(3,IStruc),NTRed,V(IIZRed+LastIZ),C,V(IXXX),V(IIB),
     $      V(IB),V(IEnd+1),MDV-IEnd)
          Call FormG2(IOut,IPrint,NTRed,V(IIB),V(IB),V(ITri))
          Call FileIO(1,-IOSleq,NNT,V(ITri),0)
          Call AClear(NAt3TX,V(ITri))
          Call FileIO(1,IOSLeq,NAt3TX,V(ITri),0)
          InitNR = .True.
          Call FileIO(1,IOSLeq,1,IIniNR,0)
          endIf
      else if(ITqry(IOZRed).eq.0.and.MicOTy.eq.4) then
        DoRed = .True.
        Call ISet(1,NZ,V(ILayer))
        Call IClear(NZ,V(ILiAto))
        Call IClear(NZ,V(IIAnMo))
        Call GenITR(In,IOut,IPrint,.True.,.False.,.False.,CkRigF,0,
     $    MaxNZ,NZ,IZ,IAnZ,NuOfLa,V(ILayer),V(ILiAto),V(IIAnMo),NAtoms,
     $    IAn,IAtTyp,IAtRes,NResid,IResNa,IResNu,V(ITRead),MicOpt,
     $    V(IOptAt),V(IV),MDV1)
        endIf
      Call TStamp(1,'Aft Large/GenITR2')
C
C     If requested, set-up the new internal coordinate definition data
C     structure using IZRed. For now, we do a quick check to ensure that
C     the system is not a periodic one. At this point, the Generalized
C     Internal Coordinate code cannot handle periodic systems.
C
      If(NewCrd) then
        If(NPDirC(NAtoms,IAn).ne.0)
     $    Call GauErr('GICs cannot be used with PBC.')
        ICrdDt = IV
        Write(IOut,2000)
        If(IPrint.ge.2) Write(IOut,2002) ICrDer,ICrDef,ICrDf0,ICrDf1,
     $    ICrDf2
        Call CrdDef(In,IOut,IPrint,IRWCrd,1,0,0,0,0,0,V(IIZRed+LastIZ),
     $    IV,V(ICrdDt),MDV1)
        If(IPrint.ge.1) Write(IOut,2010)
        Call CrdDef(In,IOut,IPrint,IRWCrd,ICrDf0,0,0,0,0,NTRed,
     $    V(IIZRed+LastIZ),IV,V(ICrdDt),MDV1)
        If(ICrDer.ge.2) then
          If(IPrint.ge.1) Write(IOut,2020)
          Call CrdDef(In,IOut,IPrint,IRWCrd,ICrDf1,-1,0,0,0,NTRed,
     $      V(IIZRed),IV,V(ICrdDt),MDV1)
          endIf
        If(ICrDer.ge.3) then
          If(IPrint.ge.1) Write(IOut,2030)
          Call CrdDef(In,IOut,IPrint,IRWCrd,ICrDf2,-1,0,0,0,NTRed,
     $      V(IIZRed),IV,V(ICrdDt),MDV1)
          endIf
        If(IPrint.ge.1) Write(IOut,2040)
        Call CrdDef(In,IOut,IPrint,IRWCrd,99,0,0,0,0,0,V(IIZRed),IV,
     $    V(ICrdDt),MDV1)
        If(IPrint.ge.1) Call CrExst(IOut,IPrint,-1,0,0,0,IExist,Insert,
     $    V,V,V,0,0,0,0,V,V)
        Call CrdDat(IOut,0,IRWCrd,1,NCoord,ICrdDt,IV,V,MDV)
        If(NCoord.ne.NTRed) then
          Write(IOut,3000) NCoord,NTRed
          Call GauErr('RdGeom: NCoord.ne.NTRed')
          endIf
        MDV1 = MDV - IV + 1
        Call CrdVal(IOut,0,IRWCrd,1,-1,C,V(IValue),V(ICrdDt),V(IV),MDV1)
        If(IPrint.ge.2) Call PrtCrd(IOut,ICrDer,NCoord,IRwCrd,C,
     $    V(IValue),V(IV),MDV1)
        endIf
C
C     Save redundant internal information and set up Cartesians in /ZMat/
C     if using redundant internals.
C
      Call ISubsI(NAtoms,MicOpt,1,-3)
      If(DoRed) then
        If(NVarZM.gt.0) Call SaveZM(In,IOut,.False.,.False.,.False.,
     $    .True.,IOZSav,V(IKeep))
        NVarZM = 0
        Call RedZM(IStruc,NAlSav,NAtoms,V(IIAnSv),V(ICSave),MaxAtm,IAn,
     $    C,MicOpt)
        SaveZ = .True.
        Call SaveZR(In,IOut,ReadFz,NoFrez,NStrAl,NAlRed,IStruc,NStruc,
     $    NAtoms,V(IIZRed),NTRedS,V(IValue),V(IIntVc),V(IFPVec),
     $    V(IKeep),NAlSav,V(ICSave),DoPath,MicOpt)
        endIf
      If(IPSave(0).eq.0.or.IPrint.gt.0) Call CntATy(IOut,MOMMFl,NuOfLa,
     $  NAtoms,MicOpt,V(ILayer),V(ILiAto))
C
C     Save the z-matrix on the rwf.
C
      If(NZ.le.0.and.NAtoms.le.0) Call GauErr('There are no atoms!')
      Call SaveZM(In,IOut,ReadFz,NoFrez,FixTet,SaveZ,0,V(IKeep))
      If(IGoI.eq.1) Call RdWrZ(16,0,V)
      NPDMol = NPDirC(NAtoms,IAn)
      NAtoms = NAtoms - NPDMol
      Call CkPBCV(IOut,NPDMol,C(1,NAtoms+1))
      If(NAtoms.gt.(MaxAtm-3)) Call GauErr('There are too many atoms.')
      Call CmplTr(NPDMol,IAn(NAtoms+1),C(1,NAtoms+1))
C
C     Secondary structure.
C
      Call TstCor(IV+NAtoms,MDV,'RdGeom-SecGet')
      If(IDoSec.eq.0) then
        If(RdChk) then
C         Sec. structure info was read with last geom.
          IDoSec = 1
        else
          IDoSec = 2
          endIf
        endIf
      Call SecGet(In,IOut,IPrint,IDoSec,NAtoms,IAn,IAtFrg,IAtRes,IAtPDB,
     $  NResid,IResNa,IResNu,NSecSt,SecTyp,SecID,SecRes,SecDat,SecAtm,
     $  V(IV))
C
C     Set atomic masses.
C
      Call TstCor(IV+2*NAtoms,MDV,'RdGeom-RedIso')
      Call RedIso(In,IOut,IPrint,IDump,Isotop,IThTem,IThPre,IThScl,
     $  IOGen,V(IV))
C
C     Generate molecular mechanics parameters.
      If(IGenMM.eq.2) then
        If(LRwMMC.gt.0) then
          If(IPSave(0).eq.0) Write(IOut,1040)
          Call CnvMMF(0,IRwMMC,LenPar,V(IV),MDV1)
          Call FileIO(1,-IrwMM,LenPar,V(IV),0)
          endIf
      else if(Mod(IParIn,100).gt.0) then
        Call TStamp(1,'Bef MMInit')
        If(IPSave(0).eq.0) Write(IOut,1050)
        Call GetCon(IOCon,NAtoms,IV,MxBond,INBond,IIBond,IRBTyp,IPOld,
     $    V,MDV)
        If(MxBond.eq.0) Call GauErr('MM requires connectivity.')
        If(IGenMM.ne.3) LRwMMC = 0
        IV = IPOld + LenPar
        If(LRwMMC.gt.0) then
          Write(IOut,1060)
          Call TstCor(IV,MDV,'RdGeom-OldPar')
          Call CnvMMF(0,IRwMMC,LenPar,V(IPOld),MDV)
        else if(IGenMM.eq.3) then
          Write(IOut,1070)
          endIf
        Call InOldP(LenPar,V(IPOld),NCOld,IParCt,IndIC,IndRC,NDSymO,
     $    NSOld,IndSTO,IStepD,NC,IIList)
        MDV1 = MDV - IV + 1
        ONIOM = NuOfLa.gt.1
        Call MMInit(IOut,IPrint,IDump,IParIn,ICntMM,IQEq,ISwitc,ICutV,
     $    ICutC,IRanCt,NAtoms,ICharg,Multip,IAn,IAtTyp,AtChMM,C,NumTpS,
     $    IAtTpS,MxBond,V(INBond),V(IIBond),V(IRBTyp),NCOld,
     $    V(IPOld+IndIC),V(IPOld+IndRC),NDSymO,NSOld,V(IPOld+IndSTO),
     $    LenMM,V(IPOld+IStepD),IMMVB,NMMVB,0,JJ,0,0,JJ,IAtRes,NResid,
     $    IResNu,IResNa,NuOfLa,V(ILayer),V(ILiAto),ONIOM,V(IIChOn),
     $    V(IMulOn),V(IV),MDV1)
        Call ConWrt(IRwMM,LenMM,V(IV))
        Call TStamp(1,'Aft MMInit')
        endIf
C
C     Save ONIOM and fragment information.
C
      Call FixACM(ICpLCh,NAtoms,V(IIAnMo),AtChMM,V(IACMod))
      If(NFrag.gt.1) then
        IFrStp = -1
        NDNBsC = 2*NFrag + 8
        Call RdWrFr(1,IRwCP,NAtoms,NFrag,IHavMO,NACore,NAVal2,NAVal1,
     $    NBCore,NBVal2,NBVal1,IFrStp,V(IIChFr),V(IMulFr),NDNBsC,
     $    V(INBCom))
        endIf
      If(IGetVl(V(ILayer),1).eq.0) Call ISet(1,NAtoms,V(ILayer))
      Call RdWrOn(1,IRwOni,MxLONI,NAtoms,V(IIChOn),V(IMulOn),V(ILayer),
     $  V(IIAnMo),V(IIAtMo),V(ILiAto),V(IACMod),MaxNZ,V(IDista))
      Call UpZMCh(IOut,IPrint,IRwMM,IParIn,NAtoms,IAtTyp,V,MDV)
C
C     Read in frequencies for electromagnetic perturbations.
C
      If(RdFreq) then
        IFreq = 1
        IInd = IFreq + MaxFrq
        IV = IInd + InToWP(MaxFrq)
        MDV1 = MDV - IV + 1
        Call InpFrq(IOut,AllChk,ForcW0,SrtFrq,IOFDPr,MaxFrq,NFrqRd,
     $    V(IFreq),V(IInd),V(IV),MDV1)
        endIf
C
C     Read in external point charges.
C
      IExtCh = Mod(IBkgCh,1000)/100
      ISCSCh = Mod(IBkgCh,10000)/1000
      Call ILSW(1,62,(IExtCh+10*ISCSCh))
      IBkgRd = Mod(IBkgCh,10)
      IBkgUn = Mod(IBkgCh,100)/10
      If(IBkgRd.ne.0) then
        jXMask = 1
        jChg = jXMask + 1
        RdExtC = IBkgRd.eq.1.or.IBkgRd.eq.4.or.IBkgRd.eq.5
        If(RdExtC) then
          MaxRd = (MDV-1) / (LChg+LDBF)
          jDBF = jChg + MaxRd*4
          Write(IOut,1080)
          Call RdBkgC(IDump,IBkgRd,IBkgUn,MaxRd,LChg,LDBF,PhyCon(1),
     $      NumChg,NumDBF,V(jChg),V(jDBF))
          Call WrChg(0,1,NumChg,NumDBF,LChg,LDBF,V(jChg),V(jDBF))
          Call WrChg(0,2,NumChg,NumDBF,LChg,LDBF,V(jChg),V(jDBF))
        else if(IBkgRd.eq.3) then
          Write(IOut,1090)
          Call RdChg(21,-IUnit(9),NumChg,NumDBF,LChgX,LDBFX,0,jChg,jDBF,
     $      jV,V,MDV)
          Call WrChg(0,1,NumChg,NumDBF,LChg,LDBF,V(jChg),V(jDBF))
          Call WrChg(0,2,NumChg,NumDBF,LChg,LDBF,V(jChg),V(jDBF))
        else if(IBkgRd.ne.2) then
          Write(IOut,1100) IBkgRd
          Call Lnk1E(0)
          endIf
        If(IPrint.ge.1.or.RdExtC) then
          PrtChg = IPrint.ge.1.or.(RdExtC.and.IPrint.eq.0)
          Call RdChg(21,0,NumChg,NumDBF,LChgX,LDBFX,0,jChg,jDBF,jV,V,
     $      MDV)
          jCh = jV
          jX = jCh + NumChg + NumDBF
          jEnd = jX + 3*(NumChg+NumDBF) - 1
          Call TstCor(jEnd,MDV,'RdGeom-PtExtC')
          Call PtExtC(IOut,IPrint,PrtChg,NumChg,LChg,V(jChg),NumDBF,
     $      LDBF,V(jDBF),V(jCh),V(jX))
          endIf
        endIf
C
C     If the structure came from the chk file, copy the orientation matrix.
C
      IOSymC = FilNum(IRwSym,IUnit(9))
      LSymC = ITqry(IOSymC)
      If(RdChk.and.LSymC.gt.0)
     $  Call RstFil(IOut,IPrint,' ',' ',IOSymC,LSymC,0,IRwSym,V,MDV)
C
C     Possibly copy external input data from the chk file.
C
      IOInp = FilNum(IRwInp,IUnit(9))
      LInp = ITqry(IOInp)
      If(LInp.gt.0.and.(ICpInp.eq.1.or.(ICpInp.eq.0.and.AllChk)))
     $  Call RstFil(IOut,IPrint,'external input','chk file',IOInp,LInp,
     $  0,IRwInp,V,MDV)
C
C     Possibly read an atomic pair potential or copy one.
C
      IOPPTC = FilNum(IRwPPT,IUnit(9))
      LPPTC = ITqry(IOPPTC)
      If(IRdPPT.eq.0.and.AllChk.and.LPPTC.gt.0) then
        IRdPPT = 3
      else if(IRdPPT.eq.0) then
        IRdPPT = 2
        endIf
      If(IRdPPT.eq.3) then
        InpPPT = IOPPTC
      else
        InpPPT = 0
        endIf
      If(IRdPPT.ne.2) Call RdPrPt(IOut,IPrint,InpPPT,IRwPPT,NAtoms,1,
     $  NPot,IHead,IITPot,IATPot,Next,V,MDV)
C
C     Analyze MM charges.
C
      Call AnlMMC(IOut,IPrint,'MM',NAtoms,AtChMM,C)
C
C     Analyze residue structure.
C
      IStart = IIZRed
      Call AllOni('RdGeom',.False.,MxLONI,NAtoms,IStart,ILayer,IIAnMo,
     $  IIAtMo,ILiAto,IACMod,IDista,Next,V,MDV)
      Call RdWrOn(-2,IRwOni,MxLONI,NAtoms,V(IIChOn),V(IMulOn),V(ILayer),
     $  V(IIAnMo),V(IIAtMo),V(ILiAto),V(IACMod),NAtoms,V(IDista))
      Call GetCon(IOCon,NAtoms,Next,MxBond,INBond,IIBond,IRBTyp,IV,V,
     $  MDV)
      MDV1 = MDV - IV + 1
      Call AnlRes(IOut,IPrint,IRRAn,PhyCon(1),NAtoms,IAn,AtmChg,C,
     $  AtChMM,IAtRes,NResid,IResNu,IResNa,MxLONI,NuOfLa,V(ILayer),
     $  V(ILiAto),MxBond,V(INBond),V(IIBond),V(IRBTyp),V(IV),MDV1)
      Call ChainX(0)
      Return
      End
*Deck AddIZ1
      Subroutine AddIZ1(IA,JA,KA,LA,IRed,IType,C,NTBond,NTAng,NTDih,
     $  NTRed,IZRed,Values,IntVec,FPVec)
      Implicit Real*8(A-H,O-Z)
C
C     Insert a coordinate into IZRed and related arrays.
C
      Dimension C(3,*), IZRed(4,*), Values(*), IntVec(*), FPVec(*)
      Save Zero
      Data Zero/0.0d0/
      Dist(I,J) = Sqrt(Abs((C(1,I)-C(1,J))**2+(C(2,I)-C(2,J))**2
     $                 +(C(3,I)-C(3,J))**2))
C
      If(IType.eq.1) then
        If(IRed.ne.(NTBond+1)) Call GauErr('IRed wrong in AddIZ1.')
        NTBond = NTBond + 1
        If(JA.lt.0) then
          Value = C(-JA,IA)
        else
          Value = Dist(IA,JA)
          endIf
      else if(IType.eq.2) then
        If(IRed.ne.(NTBond+NTAng+1))
     $    Call GauErr('IRed wrong in AddIZ1.')
        NTAng = NTAng + 1
        If(LA.eq.0) then
          Call AngTst(IA,JA,KA,C,Value)
        else
          Call LinAng(IA,JA,KA,LA,C,Value)
          endIf
      else if(IType.eq.3) then
        If(IRed.ne.(NTBond+NTAng+NTDih+1))
     $    Call GauErr('IRed wrong in AddIZ1.')
        NTDih = NTDih + 1
        Call RTors(IA,JA,KA,LA,C,Value)
      else
        Value = Zero
        Call GauErr('AddIZ1 is confused about IType.')
        endIf
      Call InsRed(IRed,IA,JA,KA,LA,Value,0,Zero,NTRed,IZRed,Values,
     $  IntVec,FPVec)
      Return
      End
*Deck AddIZR
      Subroutine AddIZR(IOut,ToBohr,ToRad,IA,JA,KA,LA,C,IZRed,Values,
     $  IntVec,FPVec,Cnstr,NTBond,NTAng,NTDih,NTRed,IAct,GotVal,Value,
     $  FC,NStep,Step)
      Implicit Real*8(A-H,O-Z)
C
C     Add a user defined redundant coordinate to IZRed.
C
      Parameter (IActiv=3,IFreez=4,IHDiag=6,IDeriv=7,IScan=8)
      Logical GotVal
      Dimension IZRed(4,*),Values(*),IntVec(*),FPVec(*),C(3,*),Cnstr(*)
 1010 Format(' Adding R(',I5,',',I5,')')
 1020 Format(' Adding A(',I5,',',I5,',',I5,')')
 1030 Format(' Adding L(',I5,',',I5,',',I5,',',I5,')')
 1040 Format(' Adding L(',I5,',',I5,',',I5,',',I5,',',I5,')')
 1050 Format(' Adding D(',I5,',',I5,',',I5,',',I5,')')
C
C     Find the coordinate and add it if necessary.
C
      Call FndCor(IA,JA,KA,LA,NTBond,NTAng,NTDih,IZRed,IRed,IType)
      If(IType.lt.0) then
        If(IType.eq.-1) then
          Write(IOut,1010) IA, JA
        else if(IType.eq.-2) then
          If(LA.eq.0) then
            Write(IOut,1020) IA, JA, KA
C           Regular angle; delete corresponding linear angles.
            Call DelAng(IA,JA,KA,-1,NTBond,NTAng,NTDih,NTRed,IZRed,
     $        Values,IntVec,FPVec)
          else if(LA.ne.0) then
            If(LA.ge.-9) then
              Write(IOut,1030) IA, JA, KA, LA
            else
              IA4A = -LA/10
              IA4T = Mod(LA,10)
              Write(IOut,1040) IA, JA, KA, IA4A, IA4T
              endIf
C           Linear bend; delete regular angles and dihedrals
            Call DelAng(IA,JA,KA,0,NTBond,NTAng,NTDih,NTRed,IZRed,
     $        Values,IntVec,FPVec)
            Call DelAng(IA,JA,KA,1,NTBond,NTAng,NTDih,NTRed,IZRed,
     $        Values,IntVec,FPVec)
            endIf
          Call FndCor(IA,JA,KA,LA,NTBond,NTAng,NTDih,IZRed,IRed,IType)
        else
          Write(IOut,1050) IA, JA, KA, LA
          endIf
        Call AddIZ1(IA,JA,KA,LA,IRed,-IType,C,NTBond,NTAng,NTDih,NTRed,
     $    IZRed,Values,IntVec,FPVec)
        endIf
C
C     Modify values if necessary.
C
      Call FndCor(IA,JA,KA,LA,NTBond,NTAng,NTDih,IZRed,IRed,IType)
      If(IType.lt.0) Call GauErr('Logic error in AddIZR.')
      If(IType.eq.1) then
        Conv = ToBohr
      else
        Conv = ToRad
        endIf
      If(GotVal) Values(IRed) = Value * Conv
      If(IAct.eq.IDeriv) then
        IntVec(IRed) = 3
      else if(IAct.eq.IHDiag) then
        IntVec(IRed) = 1
        FPVec(IRed) = FC
      else if(IAct.eq.IFreez) then
        IntVec(IRed) = -1
        Cnstr(IRed) = Values(IRed)
      else if(IAct.eq.IActiv) then
        IntVec(IRed) = 0
      else if(IAct.eq.IScan) then
        IntVec(IRed) = -1001 - NStep
        Cnstr(IRed) = Values(IRed)
        FPVec(IRed) = Step * Conv
        endIf
      Return
      End
*Deck AddRed
      Subroutine AddRed(IOut,NAtoms,IZRed,ToAng,NTRed,NTBond,NTAng,
     $  NTDih,C,MicOpt,Values,IntVec,FPVec,Cnstr,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Read and add user defined redundant coordinates
C
C     Input (free format)
C      I, J, [ K, [ L ]], [ Val ], [ F|D|R|A [ FPval]]
C                               or [ S isteps stepsize]
C
C     Output
C      IZRed, NTBond, NTAng, NTDih and NTRed are updated
C      IntVec   - code indicating modification of coordinates
C               F = -1  frozen
C               D = 3   differentiate
C               H = 1   diagonal hessian value
C               R or K  delete coordinate
C               blank or A    (default)
C                       activate retaining current value
C               A Val   activate and set current value to Val
C               S = -1000-isteps   scan the variable
C      FPVec    - step size for numerical differentiation or scan
C
      Parameter (IKillB=2,IRem=5,IL3I=5,ILinP=6)
      Logical ModRed, GotVal
      Dimension IZRed(4,*),C(3,*),V(*),Values(*),IntVec(*),FPVec(*),
     $  Cnstr(*),MicOpt(*),IIA(4)
      Equivalence (IA,IIA(1)),(JA,IIA(2)),(KA,IIA(3)),(LA,IIA(4))
      Save One
      Data One/1.D0/
C
      IPrint = 0
      ModRed = .False.
      ToBohr = One/ToAng
      ToRad = ATan(One)/GFloat(45)
   10 Call Rd1Cor(IOut,.False.,NAtoms,C,MicOpt,NTRed,NTBond,IZRed,
     $  IntVec,IAct,ICType,IIA,GotVal,Value,FC,NStep,Step)
C
C     Delete a coordinate
C
      If(IAct.eq.0) then
        Goto 200
      else if(IAct.eq.IKillB.or.IAct.eq.IRem) then
        Call DelIZR(IOut,IA,JA,KA,LA,IZRed,Values,IntVec,FPVec,NTBond,
     $    NTAng,NTDih,NTRed)
C
C     Add a coordinate
C
      else
        If(LA.lt.0.and.LA.gt.-10.and.IA.gt.KA) then
          ITemp = IA
          IA = KA
          KA = ITemp
          endIf
        Call AddIZR(IOut,ToBohr,ToRad,IA,JA,KA,LA,C,IZRed,Values,IntVec,
     $    FPVec,Cnstr,NTBond,NTAng,NTDih,NTRed,IAct,GotVal,Value,FC,
     $    NStep,Step)
        If(ICType.eq.IL3I.or.ICType.eq.ILinP) Call AddIZR(IOut,ToBohr,
     $    ToRad,IA,JA,KA,LA-1,C,IZRed,Values,IntVec,FPVec,Cnstr,NTBond,
     $    NTAng,NTDih,NTRed,IAct,GotVal,Value,FC,NStep,Step)
        ModRed = ModRed.or.GotVal
        endIf
      Goto 10
C
C     Fix the cartesian coordinates if any of the redundant internals
C     were modified.
C
  200 If(ModRed) then
        Call NewRed(IOut,IPrint,NAtoms,NTRed,NTBond,NTAng,NTDih,IZRed,
     $    Values,Cnstr,IntVec,C,V,MDV)
        Call RedVal(IOut,IPrint,.True.,C,IZRed,Values,NTRed,NTBond,
     $    NTAng,NTDih,V,MDV)
        endIf
      Return
      End
*Deck AlCart
      Logical Function AlCart(IDummy)
      Implicit Real*8(A-H,O-Z)
C
C     Return whether an all-cartesian fixed z-matrix is present.
C     An empty z-matrix returns .false.
C
#include "commonz.inc"
C
      AlCart = .False.
      If(NZ.lt.1) Return
      Do 10 I = 1, NZ
        If(IZ(I,1).gt.0.or.LBl(I).ne.0.or.LAlpha(I).ne.0.or.
     $    LBeta(I).ne.0) Return
   10   Continue
      AlCart = .True.
      Return
      End
*Deck AngOK
      Logical Function AngOK(I,J,K,C,Cut)
      Implicit Real*8(A-H,O-Z)
C
C     Check if the I-J-K angle is within Cut of linear, and return
C     .false. if it is.
C
      Dimension C(3,1)
      Save One, Two, Four
      Data One/1.0d0/, Two/2.0d0/, Four/4.0d0/
      Dist(I,J) = Sqrt((C(1,I)-C(1,J))**2+(C(2,I)-C(2,J))**2
     $                +(C(3,I)-C(3,J))**2)
C
      Pi = Four*ATan(One)
      RIJ = Dist(I,J)
      RIK = Dist(I,K)
      RJK = Dist(J,K)
      Arg = (RIJ*RIJ + RJK*RJK - RIK*RIK)  /  (Two * RIJ * RJK)
      Ang = GACos1(Arg)
      AngOK = Ang.ge.Cut .and. Ang.le.(Pi-Cut)
      Return
      End
*Deck AnlRes
      Subroutine AnlRes(IOut,IPrint,IRRAn,ToAng,NAtoms,IAn,AtmChg,C,
     $  AtChMM,IAtRes,NResid,IResNu,IResNa,MaxLay,NuOfLa,Layer,LiAtom,
     $  MxBond,NBond,IBond,RBType,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Analyze residue structure, if present.
C
      Parameter (MxADef=10000)
      Dimension IAn(*), AtmChg(*), C(3,*), AtChMM(*), IAtRes(*),
     $  IResNu(*), IResNa(*), Layer(*), LiAtom(*), NBond(*),
     $  IBond(MxBond,*), RBType(MxBond,*), V(*)
 1000 Format(' Analysis of residues requested, but no residue',
     $  ' information present.')
C
      IRRAn0 = Mod(IRRAn,10)
      If(IRRAn0.eq.0) then
        If(NAtoms.le.MxADef.and.NResid.gt.0) then
          IRRAn0 = 1
        else
          IRRAn0 = 2
          endIf
        endIf
      If(IRRAn0.eq.2) Return
      If(NResid.eq.0) then
        Write(IOut,1000)
        Return
        endIf
      Call TStamp(1,'Top AnlRes')
      IRRAn1 = Mod(IRRAn,100)/10
C     Use residue 0 for atoms not marked as part of a residue.
      NRDim = NResid + 1
      NRDTT = (NRDim*(NRDim+1))/2
      IAtRs1 = 1
      INAtRs = IAtRs1 + InToWP(NAtoms)
      ILRes  = INAtRs + InToWP(NRDim)
      IResCX = ILRes + InToWP(2*NRDim)
      IResCM = IResCX + 3*NRDim
      IResDM = IResCM + NRDim*MaxLay
      IResMD = IResDM + NRDTT
      IResAD = IResMD + NRDim
      IScr = IResAD + NRDim
      IIScr = IResDM + NAtoms
      IEnd = IResDM + InToWP(NAtoms) - 1
      Call TstCor(IEnd,MDV,'AnlRes')
      Call AnlRs1(IOut,IPrint,IRRAn1,ToAng,NAtoms,IAn,AtmChg,C,AtChMM,
     $  IAtRes,NResid,IResNu,IResNa,MaxLay,NuOfLa,Layer,LiAtom,MxBond,
     $  NBond,IBond,RBType,V(IAtRs1),V(IResCM),V(IResDM),V(IResMD),
     $  V(IResAD),V(INAtRs),V(IResCX),V(ILRes),V(IIScr),V(IScr))
      Call TStamp(1,'Bot AnlRes')
      Return
      End
*Deck AnlRs1
      Subroutine AnlRs1(IOut,IPrint,IRRAn1,ToAng,NAtoms,IAn,AtmChg,C,
     $  AtChMM,IAtRes,NResid,IResNu,IResNa,MaxLay,NuOfLa,Layer,LiAtom,
     $  MxBond,NBond,IBond,RBType,IAtRs1,ResACM,ResDis,ResMD,ResAD,
     $  NAtRes,ResCX,LRes,IScr,Scr)
      Implicit Real*8(A-H,O-Z)
C
C     Analyze residue structure for AnlRes.
C
      Logical Have0
      Dimension IAn(*), AtmChg(*), C(3,*), AtChMM(*), IAtRes(*),
     $  IResNu(*), IResNa(*), Layer(*), LiAtom(*), NBond(*), IAtRs1(*),
     $  IBond(MxBond,*), RBType(MxBond,*), ResACM(MaxLay,0:*),
     $  ResDis(*), ResMD(0:*), ResAD(0:*), IScr(*), Scr(*), NAtRes(0:*),
     $  ResCX(3,0:*), LRes(2,0:*), CX(3)
      Character*1 IC, Ch1, JC, JCh1, NamLay, LayNam, RLName*8
      Save Zero, One, Small, Big
      Data Zero/0.0d0/, One/1.0d0/, Small/1.d-4/, Big/1.d20/
 1010 Format(' Atom',I6,' has IResAt=',I20,' but NResid=',I6,'.')
 1020 Format(' Center=',3F16.8)
 1040 Format('        ',I4,'           ',6X,2X,1X,3X,' charge ',F13.8,
     $  ' Dist=',F13.8)
 1050 Format('        ',I4,' PDB Number',I6,2A1,1X,A3,' charge ',F13.8,
     $  ' Dist=',F13.8)
 1060 Format('        ',4X,'           ',6X,2X,1X,3X,' charge ',F13.8,
     $  ' in layer ',A)
 1070 Format(' Residue',I4,'           ',6X,2X,1X,3X,' charge ',
     $  F13.8,' layers ',A,' AvRTP=',3F16.8)
 1080 Format(' Residue',I4,' PDB Number',I6,2A1,1X,A3,' charge ',
     $  F13.8,' layers ',A,' AvRTP=',3F16.8)
 1090 Format('        ',4X,'           ',6X,2X,1X,3X,' atom',I6,' link',
     $  I6,' ch=',f10.6,' corr=',2f10.6)
 1100 Format('        ',4X,'           ',6X,2X,1X,3X,' atom',I6,
     $  ' bonded to',I6,' res',I5,1X,A,' rb=',f10.3,' ch=',f10.6,
     $  ' corr=',f10.6)
C
      ToDeg = GFloat(45) / ATan(One)
      If(NuOfLa.le.0) then
        NLay = 1
      else
        NLay = Max(NuOfLa,IArMax(.False.,NAtoms,Layer))
        endIf
      NRDim = NResid + 1
      NRDTT = (NRDim*(NRDim+1))/2
      Call IClear(NRDim,NAtRes)
      Call AClear(3,CX)
      Call AClear(3*NRDim,ResCX)
      Call AClear(NRDim*MaxLay,ResACM)
      Call ASet(NRDTT,Big,ResDis)
      Do 10 IRes = 1, NRDim
        LRes(1,IRes-1) = MaxLay + 1
        LRes(2,IRes-1) = 0
   10   ResDis((IRes*(IRes+1))/2) = Zero
      Have0 = .False.
      Do 40 I = 1, NAtoms
        IRes = IAtRes(I)
        If(IRes.lt.0.or.IRes.gt.NResid) then
          Write(IOut,1010) I, IRes, NResid
          Call GauErr('Bad data in AnlRs1.')
        else if(IRes.gt.0) then
          If(IResNu(IRes).eq.0) IRes = 0
          endIf
        Have0 = Have0.or.IRes.eq.0
        IAtRs1(I) = IRes
        NAtRes(IRes) = NAtRes(IRes) + 1
        IIR = (IRes*(IRes+1))/2 + 1
        Lay = Min(LayAtm(NLay,I,Layer,LiAtom),MaxLay)
        LRes(1,IRes) = Min(LRes(1,IRes),Lay)
        LRes(2,IRes) = Max(LRes(2,IRes),Lay)
        ResACM(Lay,IRes) = ResACM(Lay,IRes) + AtChMM(I)
        XI = C(1,I)
        YI = C(2,I)
        ZI = C(3,I)
        ResCX(1,IRes) = ResCX(1,IRes) + XI
        ResCX(2,IRes) = ResCX(2,IRes) + YI
        ResCX(3,IRes) = ResCX(3,IRes) + ZI
        CX(1) = CX(1) + XI
        CX(2) = CX(2) + YI
        CX(3) = CX(3) + ZI
        Do 20 J = 1, (I-1)
   20     Scr(J)= ToAng*Sqrt((C(1,J)-XI)**2+(C(2,J)-YI)**2
     $      +(C(3,J)-zI)**2)
        Do 30 J = 1, (I-1)
          JRes = IAtRs1(J)
          If(IRes.ge.JRes) then
            IJR = IIR + JRes
          else
            IJR = (JRes*(JRes+1))/2 + 1 + IRes
            endIf
   30     ResDis(IJR) = Min(ResDis(IJR),Scr(J))
   40   Continue
      If(IRRAn1.eq.1) Call LTOutS(IOut,'Res-Res Distance Matrix',0,
     $  NRDim,ResDis,0)
      CX(1) = CX(1) * ToAng / GFloat(NAtoms)
      CX(2) = CX(2) * ToAng / GFloat(NAtoms)
      CX(3) = CX(3) * ToAng / GFloat(NAtoms)
      If(IPrint.ge.0) Write(IOut,1020) CX
      Do 45 IRes = 0, NResid
        NAR = NAtRes(IRes)
        If(NAR.gt.0) then
          Sc = ToAng / GFloat(NAR)
          ResCX(1,IRes) = ResCX(1,IRes)*Sc - CX(1)
          ResCX(2,IRes) = ResCX(2,IRes)*Sc - CX(2)
          ResCX(3,IRes) = ResCX(3,IRes)*Sc - CX(3)
          endIf
   45   Continue
      Call AClear(NRDim,ResAD)
      Call ASet(NRDim,Big,ResMD)
      Do 60 IRes = 1, NResid
        IIR = (IRes*(IRes+1))/2 + 1
        ResAD(0) = ResAD(0) + ResDis(IIR)
        ResMD(0) = Min(ResMD(0),ResDis(IIR))
        Do 50 JRes = 1, IRes
          RIJ = ResDis(IIR+JRes)
          ResAD(IRes) = ResAD(IRes) + RIJ
          ResAD(JRes) = ResAD(JRes) + RIJ
          ResMD(IRes) = Min(ResMD(IRes),RIJ)
   50     ResMD(JRes) = Min(ResMD(JRes),RIJ)
   60   Continue
      Call AScale(NRDim,One/GFloat(NResid),ResAD,ResAD)
      Call AHpSrt(NResid,IScr,ResAD(1))
      If(Have0) then
        ISt = NResid + 1
      else
        ISt = NResid
        endIf
      Do 80 I = ISt, 1, -1
        If(I.gt.NResid) then
          IRes = 0
          IRN = 0
        else
          IRes = IScr(I)
          IRN = IResNu(IRes)
          endIf
        Call AnlRs2(NLay,ResACM(1,IRes),LMinC,LMaxC,ChTot)
        If(LRes(1,IRes).gt.MaxLay) then
          LMin = 0
          LMax = 0
        else
          LMin = LRes(1,IRes)
          LMax = LRes(2,IRes)
          endIf
        RLName = ' '
        IRLN = 1
        Do 63 L = LMin, LMax
          RLName(IRLN:IRLN) = LayNam(L,NLay)
   63     IRLN = IRLN + 2
        Call SphCor(ResCX(1,IRes),ResCX(2,IRes),ResCX(3,IRes),RRes,
     $    ThRes,PhRes)
        ThRes = ThRes*ToDeg
        PhRes = PhRes*ToDeg
        If(IPrint.ge.0) then
          If(I.eq.0) then
            Write(IOut,1070) IRes, ChTot, RLName, RRes, ThRes, PhRes
          else if(IRN.eq.0) then
            Write(IOut,1070) IRes, ChTot, RLName, RRes, ThRes, PhRes
          else
            Call GResN2(IRN,IR,IC,Ch1)
            Write(IOut,1080) IRes, IR, IC, Ch1, IResNa(IRes),
     $        ChTot, RLName, RRes, ThRes, PhRes
            endIf
          endIf
        If(LMinC.lt.LMaxC.and.IPrint.ge.0) then
          Do 65 L = LMinC, LMaxC
            NamLay = LayNam(L,NLay)
            Write(IOut,1060) ResACM(L,IRes), NamLay
   65       Continue
          ChL1 = ResACM(LMinC,IRes)
          ChL2 = ResACM(LMaxC,IRes)
          RInt1 = GRound(ChL1,0) - ChL1
          RInt2 = GRound(ChL2,0) - ChL2
          N1 = 0
          N2 = 0
          Do 67 IA = 1, NAtoms
            If(IAtRs1(IA).eq.IRes.and.LiAtom(IA).ne.0) then
              LIA = LiAtom(IA)
              If(IAtRs1(LIA).eq.IRes) then
                N1 = N1 + 1
                Do 66 IB = 1, NBond(IA)
                  JA = IBond(IB,IA)
                  JRes = IAtRs1(JA)
                  If(IRes.eq.JRes.and.JA.ne.LIA.and.
     $              Layer(JA).eq.Layer(IA)) N2 = N2 + 1
   66             Continue
                endIf
              endIf
   67       Continue
          If(Abs(RInt1+RInt2).lt.Small.and.N1.gt.0.and.N2.gt.0) then
            RFix1 = GRound(RInt1/GFloat(N1),6)
            RFix2 = GRound(RInt2/GFloat(N2),6)
          else
            RFix1 = Zero
            RFix2 = Zero
            endIf
          ND1 = 0
          ND2 = 0
          Do 69 IA = 1, NAtoms
            If(IAtRs1(IA).eq.IRes.and.LiAtom(IA).ne.0) then
              LIA = LiAtom(IA)
              ND1 = ND1 + 1
              If(ND1.eq.N1) then
                RCorr = RInt1
              else
                RCorr = Min(Abs(RInt1),Abs(RFix2))
                RCorr = GRound(Sign(RCorr,RInt1),6)
                endIf
              ChNew = AtChMM(IA) + RCorr
              RInt1 = RInt1 - RCorr
              Write(IOut,1090) IA, LIA, AtChMM(IA), ChNew
              Do 68 IB = 1, NBond(IA)
                JA = IBond(IB,IA)
                JRes = IAtRs1(JA)
                LayJ = LayAtm(NLay,JA,Layer,LiAtom)
                If(JA.ne.LIA.and.Layer(JA).eq.Layer(IA)) then
                  ND2 = ND2 + 1
                  If(ND2.eq.N2) then
                    RCorr = RInt2
                  else
                    RCorr = Min(Abs(RInt2),Abs(RFix2))
                    RCorr = GRound(Sign(RCorr,RInt2),6)
                    endIf
                  ChNew = AtChMM(JA) + RCorr
                  RInt2 = RInt2 - RCorr
                  Write(IOut,1100) IA, JA, JRes, LayNam(LayJ,NLay),
     $              RBType(IB,IA), AtChMM(JA), ChNew
                  endIf
   68           Continue
              endIf
   69       Continue
          endIf
        If(Abs(ChTot).gt.Small.and..false.) then
          JRes = 0
          IIR = (IRes*(IRes+1))/2 + 1
          RIJMin = Big
          Do 70 KRes = 1, NResid
            Call AnlRs2(NLay,ResACM(1,KRes),LMinK,LMaxK,ChTotK)
            If(KRes.ne.IRes.and.Abs(ChTot+ChTotK).lt.Small)
     $      then
              If(IRes.ge.KRes) then
                IKR = IIR + KRes
              else
                IKR = (KRes*(KRes+1))/2 + 1 + IRes
                endIf
              If(ResDis(IKR).lt.RIJMin) then
                RIJMin = ResDis(IKR)
                JRes = KRes
                endIf
              endIf
   70       Continue
          JRN = IResNu(JRes)
          Call AnlRs2(NLay,ResACM(1,JRes),LMinJ,LMaxJ,ChTotJ)
          If(JRN.eq.0) then
            Write(IOut,1040) JRes, ChTotJ, RIJMin
          else
            Call GResN2(JRN,JR,JC,JCh1)
            Write(IOut,1050) JRes, JR, JC, JCh1, IResNa(JRes),
     $        ChTotJ, RIJMin
            endIf
          endIf
   80   Continue
      Return
      End
*Deck AnlRs2
      Subroutine AnlRs2(NLay,ChLay,LMin,LMax,ChTot)
      Implicit Real*8(A-H,O-Z)
C
C     Return the range of Layers for which a residue has charges and its
C     total charges.
C
      Dimension ChLay(*)
      Save Zero, Small
      Data Zero/0.0d0/, Small/1.0d-6/
C
      LMin = NLay + 1
      LMax = 0
      ChTot = Zero
      Do 10 I = 1, NLay
        If(Abs(ChLay(I)).ge.Small) then
          LMin = Min(LMin,I)
          LMax = Max(LMax,I)
          endIf
   10   ChTot = ChTot + ChLay(I)
      If(LMax.eq.0) then
        LMin = 1
        LMax = 1
        endIf
      Return
      End
*Deck ApFzAt
      Subroutine ApFzAt(NAtoms,OptAtm,ITRead,MicOpt)
      Implicit Integer(A-Z)
C
C     Turn optimization of atoms on or off in ITRead and MicOpt based on
C     OptAtm.  See GenITr for the possibles values in ITRead and MicOpt.
C
      Logical OptAtm(*)
      Dimension ITRead(*), MicOpt(*)
C
      Do 10 I = 1, NAtoms
        MI = MicOpt(I)
        If(OptAtm(I)) then
          MicOpt(I) = MicNoF(MI)
        else if(MI.eq.0) then
          MicOpt(I) = -3
        else if(MI.eq.-1) then
          MicOpt(I) = -2
        else if(MI.ge.2) then
          MicOpt(I) = -(MI+40)
          endIf
   10   Continue
      Return
      End
*Deck BiltIn
      Logical Function BiltIn(J,I,M,NameSh,F,An,Keep)
      Implicit Integer(A-Z)
C
C     Check if NameSh specifies a standard group.  If it does,
C     return .true., add the group to the ends of F and An, and
C     increment M.  F(J,I) is the current position in F, which
C     is loaded appropriately if NameSh is recognized.  Keep
C     is always set .true. for the groups currently built in.
C
      Logical StrEqC, Keep(1)
      Dimension F(6,1), An(1), NameSh(1)
      Dimension NumAt(10), Name(10), FS(5,78), FS1(5,44), FS2(5,34)
      Equivalence (FS(1,1),FS1(1,1)), (FS(1,45),FS2(1,1))
      Save NumAt, NumNam, Name, FS
      Data NumAt/4,7,10,10,13,13,13,2,3,3/, NumNam/10/
      Data Name/3hME ,3hET ,3hNPR,3hIPR,3hNBU,3hIBU,
     $          3hTBU,3hOH ,3hNH2,3hNO2/
C     Methyl
      Data FS1/6, 0, 2, 3, 4,     1, 1, 0, 0, 0,     1, 1, 0, 0, 0,
     $         1, 1, 0, 0, 0,
C     Ethyl
     $         6, 0, 2, 3, 4,     6, 1, 5, 6, 7,     1, 1, 0, 0, 0,
     $         1, 1, 0, 0, 0,     1, 2, 0, 0, 0,     1, 2, 0, 0, 0,
     $         1, 2, 0, 0, 0,
C     Normal Propyl
     $         6, 0, 2, 4, 5,     6, 1, 3, 6, 7,     6, 2, 8, 9,10,
     $         1, 1, 0, 0, 0,     1, 1, 0, 0, 0,     1, 2, 0, 0, 0,
     $         1, 2, 0, 0, 0,     1, 3, 0, 0, 0,     1, 3, 0, 0, 0,
     $         1, 3, 0, 0, 0,
C     Iso Propyl
     $         6, 0, 2, 3, 4,     6, 1, 5, 6, 7,     6, 1, 8, 9,10,
     $         1, 1, 0, 0, 0,     1, 2, 0, 0, 0,     1, 2, 0, 0, 0,
     $         1, 2, 0, 0, 0,     1, 3, 0, 0, 0,     1, 3, 0, 0, 0,
     $         1, 3, 0, 0, 0,
C     Normal Butyl
     $         6, 0, 2, 5, 6,     6, 1, 3, 7, 8,     6, 2, 4, 9,10,
     $         6, 3,11,12,13,     1, 1, 0, 0, 0,     1, 1, 0, 0, 0,
     $         1, 2, 0, 0, 0,     1, 2, 0, 0, 0,     1, 3, 0, 0, 0,
     $         1, 3, 0, 0, 0,     1, 4, 0, 0, 0,     1, 4, 0, 0, 0,
     $         1, 4, 0, 0, 0/
C     Iso Butyl
      Data FS2/6, 0, 2, 5, 6,     6, 1, 3, 4, 7,     6, 2, 8, 9,10,
     $         6, 2,11,12,13,     1, 1, 0, 0, 0,     1, 1, 0, 0, 0,
     $         1, 2, 0, 0, 0,     1, 3, 0, 0, 0,     1, 3, 0, 0, 0,
     $         1, 3, 0, 0, 0,     1, 4, 0, 0, 0,     1, 4, 0, 0, 0,
     $         1, 4, 0, 0, 0,
C     Tertiary Butyl
     $         6, 0, 2, 3, 4,     6, 1, 5, 6, 7,     6, 1, 8, 9,10,
     $         6, 1,11,12,13,     1, 2, 0, 0, 0,     1, 2, 0, 0, 0,
     $         1, 2, 0, 0, 0,     1, 3, 0, 0, 0,     1, 3, 0, 0, 0,
     $         1, 3, 0, 0, 0,     1, 4, 0, 0, 0,     1, 4, 0, 0, 0,
     $         1, 4, 0, 0, 0,
C     Hydroxyl
     $         8, 0, 2, 0, 0,     1, 1, 0, 0, 0,
C     Amino
     $         7, 0, 2, 3, 0,     1, 1, 0, 0, 0,     1, 1, 0, 0, 0,
C     Nitro
     $         7, 0, 2, 3, 0,     8, 1, 0, 0, 0,     8, 1, 0, 0, 0/
C
      IStart = 0
      Biltin = .false.
      Do 20 INam = 1, NumNam
          If(.not.StrEqC(NameSh,Name(INam),3)) goto 20
              Biltin = .true.
              F(J,I) = M + 1
              Lim = NumAt(INam)
              Do 10 IAt = 1, Lim
                  An(M+1+IAt) = FS(1,IAt+IStart)
                  Keep(M+IAt) = .true.
                  Do 10 L = 1, 4
                      If(FS(L+1,IAt+IStart).gt.0)
     $                    F(L,M+IAt) = FS(L+1,IAt+IStart) + M
   10                 Continue
              F(1,M+1) = I
              M = M + NumAt(INam)
              Return
   20     IStart = IStart + NumAt(INam)
      Return
      End
*Deck BlModA
      Function BlModA(IAn,ICn,JAn,JCn,IType)
      Implicit Real*8(A-H,O-Z)
C
C     Return the model A bond length of an atom with atomic number
C     IAn and connectivity ICn bonded to atomic number JAn,
C     connectivity JCn, by a bond of type IType.  The special
C     cases (such as C-N in amides) are not handled here.  If no
C     model A bond length is difined by the specified combination,
C     then 0.0 is returned.
C
C     Bond lengths for the currently defined atoms are stored in the
C     order C2, C3, C4, N2, N3, O1, O2, H, Li, F.
C
      Dimension IStart(9), MinCn(9), MaxCn(9), BLSing(55), BLDoub(55)
      Save MaxAn, Zero, IStart, MinCn, MaxCn, BLSing, BLDoub
      Data MaxAn/9/, Zero/0.0d0/
      Data IStart/8,0,9,0,0,1,4,6,10/
      Data  MinCn/1,0,1,0,0,2,2,1, 1/
      Data  MaxCn/1,0,1,0,0,4,3,2, 1/
      Data BLSing/1.38d0,1.45d0,1.46d0,1.46d0,1.52d0,1.54d0,1.33d0,
     $            1.40d0,1.47d0,1.45d0,1.33d0,1.40d0,1.47d0,1.45d0,
     $            1.45d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,
     $            1.36d0,1.36d0,1.43d0,1.41d0,1.36d0,0.00d0,1.48d0,
     $            1.06d0,1.08d0,1.09d0,0.99d0,1.01d0,0.00d0,0.96d0,
     $            0.74d0,0.00d0,0.00d0,2.10d0,0.00d0,0.00d0,0.00d0,
     $            0.00d0,1.595d0,0.0d0,1.30d0,1.33d0,1.36d0,1.36d0,
     $            1.36d0,0.00d0,1.42d0,0.92d0,0.00d0,1.42d0/
      Data BLDoub/1.28d0,1.31d0,1.34d0,0.00d0,0.00d0,0.00d0,1.32d0,
     $            1.32d0,0.00d0,1.25d0,0.00d0,0.00d0,0.00d0,0.00d0,
     $            0.00d0,1.16d0,1.22d0,0.00d0,1.22d0,1.24d0,1.21d0,
     $            34*0.00d0/
      LInd(I,J) = (((Max(I,J)*(Max(I,J)-1))/2)+Min(I,J))
C
      BlModA = Zero
      If(IAn.lt.1.or.IAn.gt.MaxAn) Goto 10
      If(IStart(IAn).eq.0.or.ICn.lt.MinCn(IAn).or.ICn.gt.MaxCn(IAn))
     $    goto 10
      If(JAn.lt.1.or.JAn.gt.MaxAn) Goto 10
      If(IStart(JAn).eq.0.or.JCn.lt.MinCn(JAn).or.JCn.gt.MaxCn(JAn))
     $    goto 10
      Num1 = IStart(IAn) + ICn - MinCn(IAn)
      Num2 = IStart(JAn) + JCn - MinCn(JAn)
      If(IType.eq.13) BlModA = BLSing(LInd(Num1,Num2))
      If(IType.eq.14) BlModA = BLDoub(LInd(Num1,Num2))
   10 Return
      End
*Deck Bonded
      Function Bonded(MxBond,InvMap,Map,IAtom,IStart,IEnd,NBond,IBond)
      Implicit Integer(A-Z)
C
C     Return the number of the first atom in the range IStart to
C     IEnd in the new numbering scheme bonded to IAtom, or 0 if there
C     is no such.  Map takes new numbers to old ones and InvMap does
C     the reverse.
C
      Dimension Map(1), InvMap(1), NBond(1), IBond(MxBond,1)
C
      IRet = IEnd + 1
      JLim = NBond(IAtom)
      Do 10 J = 1, JLim
        If(InvMap(IBond(J,IAtom)).ge.IStart.and.
     $     InvMap(IBond(J,IAtom)).le.IEnd)
     $    IRet = Min(IRet,InvMap(IBond(J,IAtom)))
   10   Continue
      If(IRet.gt.IEnd) IRet = 0
      If(IRet.eq.0) goto 999
      JLim = NBond(Map(IRet))
      Do 20 J = 1, JLim
        If(IBond(J,Map(IRet)).eq.IAtom) goto 999
   20   Continue
      Call GauErr('Bonded was given an inconsistent IBond array.')
  999 Bonded = IRet
      Return
      End
*Deck CkAxis
      Logical Function CkAxis(NAtoms,I,J,Angle,Cut,IAn,C)
      Implicit Real*8(A-H,O-Z)
C
C     Determine if the molecule is symmetric under a rotation
C     by Angle about the I-J bond.
C
      Dimension IAn(1), C(3,1), Axis(3), VMJ(3), VRot(3), VMM(3)
C
      Call ASUnit(3,C(1,I),C(1,J),Axis)
      Do 20 M = 1, NAtoms
          If(M.eq.I.or.M.eq.J) goto 20
          Call ASub(3,C(1,M),C(1,J),VMJ)
          Call RotVV(Axis,Angle,VMJ,VRot)
          Call AAdd(3,VRot,C(1,J),VRot)
          Do 10 MM = 1, NAtoms
              If(MM.eq.I.or.MM.eq.J.or.IAn(MM).ne.IAn(M)) goto 10
              Call ASub(3,VRot,C(1,MM),VMM)
              If(Sqrt(SProd(3,VMM,VMM)).lt.Cut) goto 20
   10         Continue
          CkAxis = .false.
          Goto 999
   20     Continue
      CkAxis = .true.
  999 Return
      End
*Deck CkConA
      Subroutine CkConA(MaxNZ,NAtoms,IZ,IAn,C,Al,Cut,AUsed,NConst,Const)
      Implicit Real*8(A-H,O-Z)
C
C     Determine which angles are constrained to their present values
C     by symmetry.  The only case which is currently recognized is
C     a bond which lies along a Cn axis, where n=2,3,4,5,6.
C
      Logical Const(1), AUsed(1), CkAxis
      Dimension IZ(MaxNZ,4), C(3,1), Al(1), IAn(1)
      Save One, Eight
      Data One/1.0d0/, Eight/8.0d0/
C
      NConst = 0
      If(NAtoms.lt.1) goto 999
      A360 = Eight * ATan(One)
      Do 20 I = 1, NAtoms
          Const(I) = .false.
          If(.not.AUsed(I)) goto 20
          Do 10 J = 2, 6
              Const(I) = CkAxis(NAtoms,I,IZ(I,1),(A360/GFloat(J)),Cut,
     $                          IAn,C)
              If(Const(I)) goto 20
   10         Continue
   20     Continue
  999 Return
      End
*Deck CkConB
      Subroutine CkConB(MaxNZ,NAtoms,IZ,IAn,C,Be,Cut,BUsed,NConst,Const)
      Implicit Real*8(A-H,O-Z)
C
C     Determine which dihedral angles in Beta are constrained to
C     their present values by symmetry.  Currently this routine
C     only checks values of 0, 90 and 180 degrees, although in high
C     symmetry point groups other values may be constrained.  For
C     each quartet of I,J,K,L which is planar, the plane is tested for
C     whether it is a plane of symmetry for the whole molecule.  If
C     the dihedral angle is constrained, it is rounded off to the
C     correct value.
C
      Logical Const(1), BUsed(1), CkPlan
      Dimension IZ(MaxNZ,4), C(3,1), Be(1), IAn(1)
      Save Zero, One, Two
      Data Zero/0.0d0/, One/1.0d0/, Two/2.0d0/
C
      NConst = 0
      If(NAtoms.lt.1) goto 999
      PiOv2 = Two * ATan(One)
      Pi = PiOv2 * Two
      Do 30 I = 1, NAtoms
        Const(I) = .false.
        If(.not.BUsed(I)) goto 30
        If(Abs(Abs(Be(I))-Pi).gt.Cut) goto 10
          Const(I) = CkPlan(NAtoms,I,IZ(I,1),IZ(I,2),Cut,IAn,C)
          If(Const(I)) Be(I) = Pi
          Goto 29
   10   If(Abs(Be(I)).gt.Cut) goto 20
          Const(I) = CkPlan(NAtoms,I,IZ(I,1),IZ(I,2),Cut,IAn,C)
          If(Const(I)) Be(I) = Zero
          Goto 29
   20   If((Abs(Abs(Be(I))-PiOv2).gt.Cut)) goto 30
          Const(I) = CkPlan(NAtoms,I,IZ(I,1),IZ(I,2),Cut,IAn,C).and.
     $               CkPlan(NAtoms,IZ(I,1),IZ(I,2),IZ(I,3),Cut,IAn,C)
          If(Const(I)) Be(I) = Sign(PiOv2,Be(I))
   29   If(Const(I)) NConst = NConst + 1
   30   Continue
  999 Return
      End
*Deck CkPBCV
      Subroutine CkPBCV(IOut,NPDir,TV)
      Implicit Real*8(A-H,O-Z)
C
C     Check that the translation vectors are acceptable.
C
      Dimension TV(3,*)
      Save TMxNeg
      Data TMxNeg/-1000.3d0/
 1000 Format(' Translation vector component less than',F6.2,':')
C
      TMin = ArrMin(TV,3*NPDir,.False.)
      If(TMin.lt.TMxNeg) then
        Write(IOut,1000) TMxNeg
        Call OutMat(IOut,0,TV,3,NPDir,3,NPDir)
        Call Lnk1E(0)
        endIf
      Return
      End
*Deck CkPlan
      Logical Function CkPlan(NAtoms,I,J,K,Cut,IAn,C)
      Implicit Real*8(A-H,O-Z)
C
C     Determine if the plane defined by atoms I, J, and K is
C     a symmetry plane of the molecule, by testing whether each
C     atom is either in the plane or has a partner on the other
C     side of it.
C
      Dimension C(3,1), VIJ(3), VJK(3), VNorm(3), VMJ(3), VMM(3), IAn(1)
      Save One, Two
      Data One/1.0d0/, Two/2.0d0/
C
      Call ASUnit(3,C(1,I),C(1,J),VIJ)
      Call ASUnit(3,C(1,K),C(1,J),VJK)
      Call VProd(VNorm,VIJ,VJK)
      R = Sqrt(SProd(3,VNorm,VNorm))
      CkPlan = Abs(R).ge.Cut
      If(.not.CkPlan) goto 999
      Call AScale(3,(One/R),VNorm,VNorm)
      If(NAtoms.lt.3) goto 999
      Do 20 M = 1, NAtoms
          If(M.eq.I.or.M.eq.J.or.M.eq.K) goto 20
          Call ASub(3,C(1,M),C(1,J),VMJ)
          Dist = SProd(3,VMJ,VNorm)
          If(Abs(Dist).lt.Cut) goto 20
          Call ACasB(3,C(1,M),VNorm,VMJ,-Two*Dist)
          Do 10 MM = 1, NAtoms
              If(IAn(MM).ne.IAn(M)) goto 10
              Call ASub(3,VMJ,C(1,MM),VMM)
              If(Sqrt(SProd(3,VMM,VMM)).lt.Cut) goto 20
   10         Continue
          CkPlan = .false.
          Goto 999
   20     Continue
  999 Return
      End
*Deck CkSymb
      Subroutine CkSymb(IfVar,NVar,Name,Len,ValIn,Symbls,NCur,Label,
     $  ValOut,NSymb,OK)
      Implicit Real*8(A-H,O-Z)
C
C     Check the next symbol in the list against the specified name.
C
      Logical IfVar, OK, StrEqC
      Dimension Name(*), Symbls(*), Str(10)
C
      Call GetB(2,Str,LenStr,Symbls,NCur)
      If(Len.eq.LenStr.and.StrEqC(Str,Name,Len)) then
        OK = .True.
        NSymb = NSymb - 1
        If(IfVar) then
          Label = ISign(NVar,Label)
        else
          ValOut = ValIn
          If(Label.lt.0) ValOut = -ValOut
          Label = 0
          endIf
        endIf
      Return
      End
*Deck CntATy
      Subroutine CntATy(IOut,MOMMFl,NLay,NAtoms,MicOpt,Layer,LiAtom)
      Implicit Real*8(A-H,O-Z)
C
C     Count and print the number of atoms of various types.
C
      Logical HaveMM
      Dimension MicOpt(*), Layer(*), LiAtom(*)
 1000 Format(' NAtoms=',I7,' NQM=  ',I7,' NQMF= ',I7,' NMMI=',I7,
     $  ' NMMIF=',I7,/,15X,' NMic= ',I7,' NMicF=',I7,'.')
 1010 Format(' NAtoms=',I7,' but NTot=',I7,'.')
C
      NMic = 0
      NMicF = 0
      NQM = 0
      NQMF = 0
      NMMI = 0
      NMMIF = 0
      HaveMM = NLay.gt.1.and.(MOMMFl.eq.2.or.MOMMFl.eq.3)
      Do 100 I = 1, NAtoms
        If(MicOpt(I).lt.-40.or.MicOpt(I).eq.-3) then
          NMicF = NMicF + 1
        else if(MicOpt(I).eq.-2) then
          If(HaveMM.and.Layer(I).eq.NLay.and.LiAtom(I).eq.0) then
            NMMIF = NMMIF + 1
          else
            NQMF = NQMF + 1
            endIf
        else if(MicOpt(I).eq.-1) then
          If(HaveMM.and.Layer(I).eq.NLay.and.LiAtom(I).eq.0) then
            NMMI = NMMI + 1
          else
            NQM = NQM + 1
            endIf
        else if(MicOpt(I).eq.0.or.MicOpt(I).ge.2) then
          NMic = NMic + 1
          endIf
  100   Continue
      NTot = NMic + NMicF + NQM + NQMF + NMMI + NMMIF
      Write(IOut,1000) NAtoms, NQM, NQMF, NMMI, NMMIF, NMic, NMicF
      If(NAtoms.ne.NTot) then
        Write(IOut,1010) NAtoms, NTot
        Call GauErr('CntATy is confused.')
        endIf
      Call ILSW(1,70,NQM+NQMF)
      Return
      End
*Deck CrdAdd
      Subroutine CrdAdd(IOut,IPrint,AddEq,ISkipI,ISkipR,DoFix,Exists,
     $  INwPrm,OK,MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $  MRVars,IType,IAct,ISet,NIVars,NRVars,IAddVr,RAddVr,InfCrd,IVars,
     $  RVars,NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to add one entry to the internal coordinate
C     definitions data structure.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     AddEq       This logical input argument indicates how this routine
C                 should proceed if the requested coordinate already exists
C                 in the coordinate definition data structure. If this
C                 argument is sent as .True., then we add a coordinate that
C                 sets the new coordinate equal to the pre-existing entry.
C                 If this argument is sent as .False., then we fill INwPrm
C                 with the pre-existing entry number and return to the
C                 calling routine.
C     ISkipI      If > 0, this element in IAddVr is ignored during
C                 comparisions with definitions already in the data
C                 structure. This argument is ignored if AddEq is FALSE.
C     ISkipR      If > 0, this element in RAddVr is ignored during
C                 comparisions with definitions already in the data
C                 structure.  This argument is ignored if AddEq is FALSE.
C     DoFix       This logical input argument indicates whether (.True.) of
C                 not (.False.) the input definition information should be
C                 modified if it is not given in canonical form via Routine
C                 CrDfFx.
C     Exists      This logical output argument indicates whether (TRUE) or
C                 not (FALSE) the requested coordinate definition already
C                 exists in the data structure. This argument is not
C                 meaningful if AddEq is FALSE.
C     INwPrm      This output integer argument is filled with the primitive
C                 number of the added definition. If AddEq=.FALSE. and an
C                 equivalent coordinate is found in the data structure,
C                 then this argument is returned with the entry number of
C                 that equivalent coordinate.
C     OK          This logical argument is returned as .True. if this
C                 routine is able to add the new coordinate to the data
C                 structure, or .False. if the new coordinate cannot be
C                 added. If OK is returned as .False., NPrim, LIVars, and
C                 LRVars have NOT been incremented.
C     MaxSet      The maximum number of definition sets allowed.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C     IType       The value to set for the coordinate type flag.
C     IAct        The value to set for the action status flag.
C     ISet        The value to set for the set number flag.
C     NIVars      The number of integer values to add to IVars.  Note that
C                 this argument is modified if the new coordinate already
C                 exists in the data structure.
C     NRVars      The number of real values to add to RVars.  Note that
C                 this argument is modified if the new coordinate already
C                 exists in the data structure.
C     IAddVr      The list of integer values to add to IVars. Its length is
C                 NIVars.  Note that this argument is modified if the new
C                 coordinate already exists in the data structure.
C     RAddVr      The list of real values to add to IVars. Its length is
C                 NRVars.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Parameter (NInfCr=DEFMAXCOORDINFO)
      Integer CrHash
      Logical DEBUG,AddEq,Exists,DoFix,OK,MultM1
      Dimension IAddVr(*),RAddVr(*),InfCrd(NInfCr,*),IVars(*),RVars(*),
     $  NCrdSt(*)
 1000 Format(1x,'Enter CrdAdd.')
 1100 Format(1x,'CrdAdd does not have enough memory to add coordinate:',
     $  /,4x,'IType =',I4,' IAct=',I10,' NIVars=',I10,' NRVars=',I10,
     $  /,4x,'NPrim =',I10,' LIVars=',I10,' LRVars=',I10,
     $  /,4x,'MPrim =',I10,' MCoord=',I10,' MIVars=',I10,
     $  /,4x,'MRVars=',I10)
 2000 Format(5x,'Adding coordinate: IPrim=',I10,' IType=',I10)
 2500 Format(/,5x,'CrdAdd: Setting hash tree values.',/,
     $  13x,'Prim  =',I12,'  IHash =',I12,/,
     $  13x,'LChild=',I12,'  RChild=',I12,'  Next  =',I12)
 2510 Format(13x,'INSERTED at...',/,
     $  13x,'Prim  =',I12,'  IHash =',I12,/,
     $  13x,'LChild=',I12,'  RChild=',I12,'  Next  =',I12,/)
 9000 Format(1x,'CrdAdd: Primitive ',I10,
     $  ' already exists as definition #',I10,'. IType=',I10)
 9010 Format(1x,'CrdAdd: IType=',I5,' NIVars=',I5,' Ni=',I5)
 9020 Format(1x,'CrdAdd: IType=',I5,' NRVars=',I5,' Nr=',I5)
 9030 Format(1x,'CrdAdd: Insert=',I10)
C
      DEBUG = .False.
      If(IPrint.ge.5.or.DEBUG) Write(IOut,1000)
      OK = .True.
      Exists = .False.
      INwPrm = 0
C
C     If requested, begin by ensuring that the requested definition is
C     given in canonical form. Then look to see if the new coordinate
C     already exists in the data structure. If it is already there, either
C     re-set things so that the new coordinate is one that sets itself
C     equal to another or simply return to the calling routine with IExist
C     set to the equivalent coordinate definition already in the data
C     structure.
C
      Call CrNiNr(IOut,IType,Ni,Nr)
      If(NIVars.ne.Ni) then
        Write(IOut,9010) IType,NIVars,Ni
        Call GauErr('CrdAdd at top: NIVars.ne.Ni.')
        endIf
      If(NRVars.ne.Nr) then
        Write(IOut,9020) IType,NRVars,Nr
        Call GauErr('CrdAdd at top: NRVars.ne.Nr.')
        endIf
      MultM1 = .False.
      IPrmEq = 0
      If(DoFix) Call CrDfFx(IOut,IPrint,IType,.True.,.False.,MultM1,
     $  IPrmEq,IAddVr)
      If(MultM1) Call GauErr('CrdAdd: CrDfFx returned MultM1 as TRUE.')
      If(IAct.lt.0) then
        JSet = 0
      else
        JSet = ISet
        endIf
      Call CrExst(IOut,IPrint,0,NPrim,ISkipI,ISkipR,IExist,Insert,
     $  InfCrd,IVars,RVars,IType,JSet,NIVars,NRVars,IAddVr,RAddVr)
      If(IExist.ne.0.and.AddEq) then
        If(IPrint.ge.5.or.DEBUG) Write(IOut,9000) NPrim+1,IExist,IType
        IType = 56
        NIVars = 1
        NRVars = 0
        IAddVr(1) = IExist
      else if(IExist.ne.0) then
        If(IPrint.ge.5.or.DEBUG) Write(IOut,9000) NPrim+1,IExist,IType
        Exists = .True.
        INwPrm = IExist
        Return
        endIf
      If(IPrmEq.ne.0) then
        INwPrm = IPrmEq
        Return
        endIf
C
C     Initialize things and ensure that the new coordinate will fit in the
C     memory space provided. This is also the place where we ensure that
C     the specified definition set is less than the maximum allowed. Also,
C     if everything looks alright, we increment NCrdSt if the new
C     coordinate is active.
C
      JPrim = NPrim+1
      If(ISet.gt.MaxSet) Call GauErr('CrdAdd: ISet > MaxSet.')
      If(ISet.lt.1) Call GauErr('CrdAdd: ISet < 1.')
      If(IAct.ge.0) then
        JCoord = NCoord+1
        NCrdSt(ISet) = NCrdSt(ISet)+1
      else
        JCoord = NCoord
        endIf
      JIVars = LIVars+NIVars
      JRVars = LRVars+NRVars
      If(JPrim.gt.MPrim.or.JCoord.gt.MCoord.or.JIVars.gt.MIVars.or.
     $  JRVars.gt.MRVars) then
        Write(IOut,1100) IType,IAct,NIVars,NRVars,NPrim,LIVars,LRVars,
     $    MPrim,MCoord,MIVars,MRVars
        If(NIVars.gt.0) Call IMtOuS(IOut,'IAddVr=',0,IAddVr,1,NIVars,1,
     $    NIVars)
        If(NRVars.gt.0) Call OutMtS(IOut,'RAddVr=',0,0,RAddVr,1,NRVars,
     $    1,NRVars)
        OK = .False.
        Return
        endIf
C
C     Append the current definition to the data structure. The InfCrd
C     position giving the total number of dependencies is set to zero here,
C     and filled in by CrdDef when it is told to finalize the coordinate
C     data set. The InfCrd position for this primitive coordinate label
C     cursor is initialized to -1.
C
      Call CrNiNr(IOut,IType,Ni,Nr)
      If(NIVars.ne.Ni) then
        Write(IOut,9010) NIVars,Ni
        Call GauErr('CrdAdd at bot: NIVars.ne.Ni.')
        endIf
      If(NRVars.ne.Nr) then
        Write(IOut,9020) NRVars,Nr
        Call GauErr('CrdAdd at bot: NRVars.ne.Nr.')
        endIf
      InfCrd(InfTyp,JPrim) = IType
      InfCrd(InfAct,JPrim) = IAct
      InfCrd(InfIVr,JPrim) = LIVars + 1
      InfCrd(InfNI,JPrim) = NIVars
      InfCrd(InfRVr,JPrim) = LRVars + 1
      InfCrd(InfNR,JPrim) = NRVars
      InfCrd(InfNDp,JPrim) = 0
      InfCrd(InfSet,JPrim) = ISet
      InfCrd(InfHsh,JPrim) = CrHash(NIVars,IType,IAddVr)
      InfCrd(InfLCh,JPrim) = 0
      InfCrd(InfRCh,JPrim) = 0
      InfCrd(InfNxt,JPrim) = 0
      InfCrd(InfLab,JPrim) = -1
      If(NIVars.gt.0) Call IMove(NIVars,IAddVr,IVars(LIVars+1))
      If(NRVars.gt.0) Call AMove(NRVars,RAddVr,RVars(LRVars+1))
      If(IPrint.ge.5.or.DEBUG) then
        Write(IOut,2000) JPrim,IType
        Call IMtOuS(IOut,'IVars=',0,IVars(LIVars+1),1,NIVars,1,NIVars)
        endIf
C
C     Update the hash value binary tree elements based on Insert from the
C     call to CrExst.
C
      If(NPrim.gt.0) then
        If(Insert.le.0) then
          Write(IOut,9030) Insert
          Call GauErr('CrdAdd: Invalid Insert from CrExst.')
          endIf
        If(InfCrd(InfHsh,JPrim).eq.InfCrd(InfHsh,Insert)) then
          InfCrd(InfNxt,JPrim) = InfCrd(InfNxt,Insert)
          InfCrd(InfNxt,Insert) = JPrim
        else if(InfCrd(InfHsh,JPrim).lt.InfCrd(InfHsh,Insert)) then
          InfCrd(InfLCh,Insert) = JPrim
        else if(InfCrd(InfHsh,JPrim).gt.InfCrd(InfHsh,Insert)) then
          InfCrd(InfRCh,Insert) = JPrim
          endIf
        endIf
      If(IPrint.ge.5.or.DEBUG) Write(IOut,2500) JPrim,
     $  InfCrd(InfHsh,JPrim),InfCrd(InfLCh,JPrim),InfCrd(InfRCh,JPrim),
     $  InfCrd(InfNxt,JPrim)
      If(IPrint.ge.5.or.DEBUG) Write(IOut,2510) Insert,
     $  InfCrd(InfHsh,Insert),InfCrd(InfLCh,Insert),
     $  InfCrd(InfRCh,Insert),InfCrd(InfNxt,Insert)
C
C     Update NPrim, LIVars, and LRVars. Then set INwPrm.
C
      NPrim = JPrim
      NCoord = JCoord
      LIVars = JIVars
      LRVars = JRVars
      INwPrm = NPrim
      Return
      End
*Deck CrdDef
      Subroutine CrdDef(In,IOut,IPrint,IRWCrd,IOpCrd,ISet,ISet2,Iq1,Iq2,
     $  NTRed,IZRed,IV,V,MDV)
      Implicit Real*8(A-H,O-Z)
C1CrdDef
C
C DESCRIPTION
C   This routine is the heart of the internal coordinate code in
C   Gaussian. More specifically, this routine is used to prepare the
C   internal coordinate data structure primarily based on data in /Mol/
C   (common /Mol/ should already be filled before calling this routine).
C   The internal coordinate data structure is stored in IRWCrd. To use
C   the default RWF file send 0 for IRWCrd.
C
C VARIABLES:
C   IOut   ... Unit for output
C   IPrint ... Debug print option
C   IRWCrd ... The RWF file number where the coordinate definition data
C	       is stored. If IRWCrd is sent as 0, then the default file
C	       number is used.
C   IOpCrd ... The options flag used to indicating what this routine is
C              supposed to do in this call.
C              --- What are we actually doing? (IOp0) ---
C                 1 ... Initialize the internal coordinate data structure.
C                 2 ... Append new coordinates by reading user
C                       specification from the input stream according to
C                       the option given in IOp1.
C                10 ... Differentiate all active/propagating coordinates in
C                       definition set ISet2. The derivatives are put into
C                       set ISet. NOTE: ISet2 cannot be passed with -1,
C                       though 0 is acceptable to indicate the current set.
C                19 ... Differentiate coordinate Iq1 wrt Iq2. Both, Iq1 and
C                       Iq2 must already be in the coordinate definition
C                       data structure.
C                20 ... Add all Cartesian coordinates.
C                21 ... Add all (inverse) bond distance coordinates. The
C                       type of bond coordinates added is controlled by
C                       IOp2.
C                30 ... Use IZRed to define the internal coordinate system.
C                       Note that IOp2 is still used in this case. So, if
C                       one wishes to start with the standard IZRed-type
C                       definitions and replace all bond stretches with
C                       their inverse, this option can be used for IOp0
C                       together with IOp2=2.
C                99 ... Write the data structure to the file.
C              --- Adding coordinates read from the input file? (IOp1) ---
C               0xx ... Default (same as 1).
C               1xx ... No.
C               2xx ... Yes, read the entire set of primitive coordinate
C                       definitions from the input file. The expected
C                       format is described in Routine CrDf01 and
C                       cooresponds to the case where its input dummy
C                       argument IRdFmt=2.
C               3xx ... Yes, read a list of conventional internal
C                       coordinate definitions from the input file. All of
C                       these coordinates are flagged as active, and any
C                       other necessary primitive definitions are
C                       automatically added. The expected format is
C                       described in Routine CrDf01 and cooresponds to the
C                       case where its input dummy argument IRdFmt=3.
C               4xx ... Yes, read in a user-defined list of coordinate
C                       definitions using the symbolic interpreter driven
C                       by Routine CrDfSy.
C               8xx ... Yes, read a list previously defined internal
C                       coordinates and add the necessary definitions to
C                       yield their derivatives with respect to a given
C                       list of atomic Cartesian coordinates. The expected
C                       format is described in Routine CrDf01 and
C                       cooresponds to the case where its input dummy
C                       argument IRdFmt=8.
C               9xx ... Yes, read a list previously defined internal
C                       coordinates and add the necessary definitions to
C                       yield their derivatives with respect to a given
C                       list of atomic Cartesian coordinates. The expected
C                       format is described in Routine CrDf01 and
C                       cooresponds to the case where its input dummy
C                       argument IRdFmt=9.
C              --- Options to use when adding bond stretch definitions. (IOp2) ---
C              0xxx ... Default (same as 1).
C              1xxx ... Use the standard bond stretch definition in Routine
C                       CrDfBd.
C              2xxx ... Use the standard inverse bond stretch definitions
C                       in Routine CrDfBd.
C
C              --- Assign the current set number to a type of coordinate. (IOp3) ---
C                  NOTE: If this flag is set more than once, the last
C                        call's set number is the one that will be made in
C                        the assignment.
C             0xxxx ... Do not make any assignments in this call.
C             1xxxx ... Assign the current set number to the internal
C                       coordinate definitions.
C             2xxxx ... Assign the current set number to the B-matrix.
C             3xxxx ... Assign the current set number to the B-matrix
C                       derivatives.
C   ISet       This input argument indicates to which coordinate set the
C              current call is adding coordinates. If ISet=0, then the set
C              number used will be the current set number. (For ISet.ne.0,
C              the 'current' set number is updated to ISet.) If ISet=-1,
C              then the current set number is increased by 1.  When IOpCrd
C              initializes the data structure, the current set number is
C              initialized to 1. Note that ISet is only used when adding
C              coordinates. For other mode settings, ISet is essentially
C              ignored.
C   ISet2      This input argument is used when two definition sets are
C              necessary.  Its behaviour is analogous to argument ISet.
C              Except where indicated in the comments for argument IOp,
C              ISet2 is ignored.
C   Iq1    ... This is a coordinate pointer used for certain IOpCrd
C              choices.
C   Iq2    ... This is a coordinate pointer used for certain IOpCrd
C              choices.
C   NTRed  ... This input argument is the number of internal coordinates
C              included in array IZRed. Except where indicated in the
C              comments for argument IOp, NTRed is ignored.
C   IZRed  ... This input array gives the internal coordinates in the
C              conventional form used in Links 101, 103, 716, and
C              elsewhere. Except where indicated in the comments for
C              argument IOp, IZRed is ignored.
C   IV     ... This output argument is returned with the position in V
C              where memory use can continue without destroying the space
C              reserved for the internal coordinate data structure.
C   V      ... Dynamic memory array
C   MDV    ... Length of V in words
C
C2HISTORY
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C?
#include "crdparams.inc"
#include "commonmol.inc"
      Integer CrDep2
      Logical DEBUG,UsSet2,DoInvR,Init,Found0,RotQM,TmpLog
      Parameter (NInfCr=DEFMAXCOORDINFO)
      Parameter (MxRdIn=DEFMAXREDINDEX)
      Parameter (JRWCrd=775)
      Dimension ITmpAr(1),IZRed(MxRdIn,NTRed,*),V(MDV)
      Save NPrim,NCoord,LIVars,LRVars,JSetCr,JSetB,JSetBx,IFlags,InfCrd,
     $  IVars,IRVars,NCrdSt,ICrdSt,Labels,JSet,MaxSet,ITmpAr
      Data Init/.False./,NPrim/0/,NCoord/0/,LIVars/0/,LRVars/0/,
     $  JSetCr/0/,JSetB/0/,JSetBx/0/,IFlags/1/,InfCrd/1/,IVars/1/,
     $  IRVars/1/,NCrdSt/1/,ICrdSt/1/,Labels/0/,JSet/0/,MaxSet/0/,
     $  ITmpAr/0/
 1000 Format(/,1x,'Enter Routine CrdDef',/,10x,'Version number ',I3,'.')
 1010 Format(1x,'CrdDef: IOpCrd=',I10,/,
     $  9x,'IOp0=',I2,' IOp1=',I2,' IOp2=',I2,' IOp3=',I2,/,
     $  9x,'DoInvR=',L1)
 1020 Format(3x,'KRWCrd=',I10)
 1050 Format(1x,'CrdDef: NPrim =',I10,' LIVars=',I10,' LRVars=',I10,/,
     $  9x,'NPrim0=',I10,' LIVar0=',I10,' LRVar0=',I10)
 1100 Format(1x,'CrdDef: NFroz =',I5,' NRF   =',I5)
 1200 Format(1x,'CrdDef: Incrementing JSet.')
 1210 Format(1x,'CrdDef: JSetCr=',I2,' JSetB=',I2,' JSetBx=',I2)
 2001 Format(9x,'Reading coordinates from the input file.')
 3000 Format(1x,'CrdDef: Setting dependency level for primitive ',I10,
     $  ' to ',I10)
 5000 Format(/,1x,'CrdDef: Defining coordinate derivative. Iq1=',I10,
     $  ' Iq2=',I10)
 9000 Format(1x,'CrdDef: Writing coordinates to the RWF file.',/,
     $  9x,'NPrim=',I10,' NCoord=',I10,' MaxSet=',I10,' LIVars=',I10,/,
     $  9x,'LRVars=',I10)
 9050 Format(1x,'CrdDef: ISet=',I4,' JSet=',I4)
 9060 Format(1x,'CrdDef: UsSet2=',L1,'IOp0=',I2,' ISet2=',I5)
 9100 Format(1x,'CrdDef: NPrim=',I10,' NFound=',I10,' NLeft=',I10)
 9102 Format(1x,'CrdDef: Primitive ',I10,' is a problem.')
 9200 Format(1x,'CrdDef: IPrim causing problem=',I10)
C
C     Start by taking care of DEBUG, IVer, IOpCrd, IRWCrd, JSet, etc.
C
      DEBUG = .False.
      IVer = -2
      If(IPrint.ge.1.or.DEBUG) Write(IOut,1000) IVer
      If(IPrint.ge.2.or.DEBUG) Write(IOut,9050) ISet,JSet
      IOp0 = Mod(IOpCrd,100)
      IOp1 = Mod(IOpCrd,1000)/100
      IOp2 = Mod(IOpCrd,10000)/1000
      IOp3 = Mod(IOpCrd,100000)/10000
      If(IOp1.eq.0) IOp1 = 1
      If(IOp2.eq.0) IOp2 = 1
      DoInvR = IOp2.eq.2
      If(IPrint.ge.2.or.DEBUG) Write(IOut,1010) IOpCrd,IOp0,IOp1,IOp2,
     $  IOp3,DoInvR
      If(IOp0.ne.1.and.IOp0.ne.2.and.IOp0.ne.10.and.IOp0.ne.19.and.
     $  IOp0.ne.20.and.IOp0.ne.21.and.IOp0.ne.30.and.IOp0.ne.99)
     $  Call GauErr('CrdDef: Invalid IOp0.')
      If(IRWCrd.eq.0) then
        KRWCrd = JRWCrd
      else if(IRWCrd.lt.0) then
        KRWCrd = -IRWCrd
      else
        KRWCrd = JRWCrd
        endIf
      UsSet2 = IOp0.eq.10
      If(.not.UsSet2) then
        JSet2 = 0
      else if(ISet2.eq.0) then
        JSet2 = JSet
      else if(ISet2.gt.0) then
        JSet2 = ISet2
      else
        Write(IOut,9060) UsSet2,IOp0,ISet2
        Call GauErr('CrdDef: Invalid ISet2.')
        endIf
      If(ISet.eq.-1) then
        If(IPrint.ge.2.or.DEBUG) Write(IOut,1200)
        JSet = JSet+1
      else if(ISet.gt.0) then
        JSet = ISet
      else if(ISet.lt.0) then
        Write(IOut,9050) ISet,JSet
        Call GauErr('CrdDef: Invalid ISet.')
        endIf
      If(IPrint.ge.2.or.DEBUG) then
        Write(IOut,1020) KRWCrd
        Write(IOut,9050) ISet,JSet
        endIf
      MxSet0 = MaxSet
      If(JSet.gt.MaxSet) MaxSet = JSet
      If(IOp3.eq.1) then
        JSetCr = JSet
      else if(IOp3.eq.2) then
        JSetB = JSet
      else if(IOp3.eq.3) then
        JSetBx = JSet
      else if(IOp3.ne.0) then
        Write(IOut,1010) IOpCrd,IOp0,IOp1,IOp2,IOp3,DoInvR
        Call GauErr('CrdDef: Invalid IOp3.')
        endIf
      If(IPrint.ge.2.or.DEBUG) Write(IOut,1210) JSetCr,JSetB,JSetBx
C
C     If requested by IOp0, write out the coordinate definition data
C     structure to the RWF. Before we can write the data structure to the
C     file we need to do two things. First, we re-set the array NCrdSt.
C     Then, we fill in the row of InfCrd that tells the number of primitive
C     coordinates each entry is dependent upon and also fill array ICrdSt.
C
      If(IOp0.eq.99) then
        If(IPrint.ge.4.or.DEBUG) Call CrPrt(IOut,IPrint,1,0,NPrim,
     $    V(InfCrd),V(IVars),V(IRVars))
        If(IPrint.ge.5.or.DEBUG) Call CrPHsh(IOut,NPrim,V(InfCrd))
        NLeft = NPrim
 210    NFound = 0
        Do 220 IPrim = 1,NPrim
          IPos = (IPrim-1)*NInfCr+InfNDp
          NDep0 = IIGet(V(InfCrd),IPos)
          If(NDep0.eq.0) then
            NDep = CrDep2(IOut,IPrim,0,V(InfCrd),V(IVars))
            Found0 = .False.
            NDepMx = 0
            If(NDep.eq.0) then
              NDepMx = 0
            else
              Do 230 JDep = 1,NDep
                INDep1 = CrDep2(IOut,IPrim,JDep,V(InfCrd),V(IVars))
                IPos = (INDep1-1)*NInfCr+InfNDp
                NDep1 = IIGet(V(InfCrd),IPos)
                Found0 = NDep1.eq.0.or.Found0
                NDepMx = Max(NDepMx,NDep1)
 230            Continue
              endIf
            If(.not.Found0) then
              If(IPrint.ge.5.or.DEBUG) Write(IOut,3000) IPrim,NDepMx+1
              IPos = (IPrim-1)*NInfCr+InfNDp
              Call IIPut(V(InfCrd),IPos,NDepMx+1)
              NFound = NFound+1
              NLeft = NLeft-1
              endIf
            endIf
 220      Continue
        If(NFound.eq.0) then
          Call CrPrt(IOut,IPrint,0,0,NPrim,V(InfCrd),V(IVars),V(IRVars))
          Write(IOut,9100) IPrim,NFound,NLeft
          Do 229 IIPrim = 1,NPrim
            IPos = (IIPrim-1)*NInfCr+InfNDp
 229        If(IIGet(V(InfCrd),IPos).eq.0) Write(IOut,9102) IIPrim
          Call GauErr('CrdDef: Failure A')
          endIf
        If(NLeft.gt.0) Goto 210
        If(IPrint.ge.1.or.DEBUG) then
          Write(IOut,9000) NPrim,NCoord,MaxSet,LIVars,LRVars
          Call IMtOuS(IOut,'CrdDef: NCrdSt=',0,V(NCrdSt),1,MaxSet,1,
     $      MaxSet)
          endIf
        If(IPrint.ge.2.or.DEBUG) Call IMtOuS(IOut,'CrdDef: ICrdSt=',0,
     $    V(ICrdSt),NCoord,1,NCoord,1)
        Call CrPtLb(IOut,IPrint,.True.,0,V(Labels),V(InfCrd))
        IIStrt = 1
        IIEnd  = IIStrt+IIGet(V(NCrdSt),1)-1
        Do 150 II = IIStrt,IIEnd
          IIPrim = IIGet(V(ICrdSt),II)
          Call CrPtLb(IOut,IPrint,.False.,IIPrim,V(Labels),V(InfCrd))
 150      Continue
        Call CrdRW(1,KRWCrd,IVer,NPrim,MaxSet,NCoord,LIVars,LRVars,
     $    JSetCr,JSetB,JSetBx,V(IFlags),V(InfCrd),V(IVars),V(IRVars),
     $    V(NCrdSt),V(ICrdSt),V(Labels))
        Return
        endIf
C
C     Estimate sizes for the coordinate definition data structure and get a
C     set of initial pointers for the associated arrays. If this is an
C     initialization call (IOp0=0), we also intialize NPrim, LIVars, and
C     LRVars.
C
      If(IOp0.eq.1.or..not.Init) then
        Init = .True.
        NPrim = 0
        LIVars = 0
        LRVars = 0
        JSet = 1
        MaxSet = 1
        endIf
      If(.not.Init)
     $  Call GauErr('CrdDef called before initialization call.')
      If(JSet.le.0) then
        Write(IOut,9050) ISet,JSet
        Call GauErr('CrdDef: Invalid JSet.')
        endIf
      NPrim0 = NPrim+10000*NAtoms
      NCoor0 = NCoord+10000*NAtoms
      LIVar0 = LIVars+NPrim0*5
      LRVar0 = LRVars+NPrim0*1
      If(IPrint.ge.2.or.DEBUG) Write(IOut,1050) NPrim,LIVars,LRVars,
     $  NPrim0,LIVar0,LRVar0
      If(IOp0.eq.1.or..not.Init) then
        Call CrMem(IOut,IPrint,1,IVer,MaxSet,NPrim0,NCoor0,LIVar0,
     $    LRVar0,MaxLab,LenLab,InfCrd,IVars,IRVars,NCrdSt,ICrdSt,Labels,
     $    IV)
        Call TstCor(IV-1,MDV,'CrdDef-1')
        Call IClear(MaxSet,V(NCrdSt))
        Call IIPut(V(Labels),1,MaxLab)
        Call IIPut(V(Labels),2,LenLab)
      else
        Call CrReSz(IOut,IPrint,TmpLog,IVer,MxSet0,NCoord,NPrim,LIVars,
     $    LRVars,MaxSet,NCoor0,NPrim0,LIVar0,LRVar0,InfCrd,IVars,IRVars,
     $    NCrdSt,ICrdSt,IV,V,MDV)
        If(MaxSet.gt.MxSet0) Call IClr1(MaxSet-MxSet0,MxSet0,V(NCrdSt))
        If(.not.TmpLog) Call GauErr('CrdDef: CrReSz failure 1.')
        endIf
      Call TstCor(IV-1,MDV,'CrdDef-2')
      If(IOp0.eq.1) Return
C
C     Using the IOp flags, add/append coordinates and/or derivatives to the
C     coordinate definition data structure.  Note that if more space is
C     needed than what has already been allocated using estimates, careful
C     reallocation of memory and shifting of data is done.
C
      NFroz = ICnt(1,-2,NAtoms,MicOpt)
      NFroz = NFroz+ICnt(1,1,NAtoms,MicOpt)
      If(NFroz.gt.0) Call GauErr('CrdDef: NFroz>0 NYI.')
      Call GetNRF(NAtoms,NRF,RotQM,MicOpt)
      If(IPrint.ge.1.or.DEBUG) Write(IOut,1100) NFroz,NRF
      If(IOp0.eq.2) then
        If(IPrint.ge.2.or.DEBUG) Write(IOut,2001)
        If(IOp0.ne.2) Call GauErr('IOp1.ne.1, but IOp0.ne.2.')
        If(IOp1.eq.4) then
          Call CrDfSy(In,IOut,IPrint,JSet,MaxSet,NPrim,NCoord,LIVars,
     $      LRVars,NPrim0,NCoor0,LIVar0,LRVar0,V(InfCrd),V(IVars),
     $      V(IRVars),V(NCrdSt))
        else
          Call CrDf01(IOut,IPrint,IOp1,DoInvR,JSet,MaxSet,NPrim,NCoord,
     $      LIVars,LRVars,NPrim0,NCoor0,LIVar0,LRVar0,V(InfCrd),
     $      V(IVars),V(IRVars),V(NCrdSt),V(ICrdSt),V(IV),MDV-IV+1)
          endIf
      else if(IOp0.eq.10) then
        If(IPrint.ge.2.or.DEBUG) Write(IOut,9060) UsSet2,IOp0,JSet2
        Call CrDfDx(IOut,IPrint,JSet,JSet2,NAtoms,MaxSet,NPrim,NCoord,
     $    LIVars,LRVars,NPrim0,NCoor0,LIVar0,LRVar0,V(InfCrd),V(IVars),
     $    V(IRVars),V(NCrdSt),V(ICrdSt),V(IV),MDV-IV+1)
      else if(IOp0.eq.19) then
        If(IPrint.ge.2.or.DEBUG) Write(IOut,5000) Iq1,Iq2
        ITmpAr(1) = Iq2
        Call CrdDq0(IOut,-1,JSet,Iq1,-555,ITmpAr,1,MaxSet,NPrim,NCoord,
     $    LIVars,LRVars,NPrim0,NCoor0,LIVar0,LRVar0,V(InfCrd),V(IVars),
     $    V(IRVars),V(NCrdSt),V(IV),MDV-IV+1)
      else if(IOp0.eq.20) then
        Call CrDfCar(IOut,IPrint,0,JSet,NAtoms,0,MaxSet,NPrim,NCoord,
     $    LIVars,LRVars,NPrim0,NCoor0,LIVar0,LRVar0,V(InfCrd),V(IVars),
     $    V(IRVars),V(NCrdSt))
      else if(IOp0.eq.21) then
        Call CrDf21(IOut,IPrint,DoInvR,NAtoms,JSet,MaxSet,NPrim,NCoord,
     $    LIVars,LRVars,NPrim0,NCoor0,LIVar0,LRVar0,V(InfCrd),V(IVars),
     $    V(IRVars),V(NCrdSt))
      else if(IOp0.eq.30) then
        Call CrDf30(IOut,IPrint,DoInvR,NTRed,IZRed,JSet,MaxSet,NPrim,
     $    NCoord,LIVars,LRVars,NPrim0,NCoor0,LIVar0,LRVar0,V(InfCrd),
     $    V(IVars),V(IRVars),V(NCrdSt))
      else
        Call GauErr('CrdDef: Unknown IOp0.')
        endIf
C
C     Resize the data structure so that we get rid of unused memory, then
C     re-fill array ICrdSt.
C
      Call CrReSz(IOut,IPrint,TmpLog,IVer,MaxSet,NCoor0,NPrim0,LIVar0,
     $  LRVar0,MaxSet,NCoord,NPrim,LIVars,LRVars,InfCrd,IVars,IRVars,
     $  NCrdSt,ICrdSt,IV,V,MDV)
      If(.not.TmpLog) Call GauErr('CrdDef: CrReSz failure 2.')
      ICrTmp = IV
      IV1    = ICrTmp + InToWP(MaxSet)
      Call TstCor(IV1-1,MDV,'CrdDef-1')
      KStCrd = NCoord-IIGet(V(NCrdSt),MaxSet)+1
      Call IIPut(V(ICrTmp),MaxSet,KStCrd)
      Do 900 i = MaxSet-1,1,-1
        KStCrd = KStCrd-IIGet(V(NCrdSt),i)
        Call IIPut(V(ICrTmp),i,KStCrd)
 900    Continue
      If(KStCrd.ne.1)
     $  Call GauErr('CrdDef: After initializing ICrTmp, KStCrd.ne.1.')
      Do 910 IPrim = 1,NPrim
        IPos = (IPrim-1)*NInfCr+InfAct
        If(IIGet(V(InfCrd),IPos).ge.0) then
          IPos = (IPrim-1)*NInfCr+InfSet
          KSet = IIGet(V(InfCrd),IPos)
          IPos = IIGet(V(ICrTmp),KSet)
          If(IPos.gt.NCoord) then
            Call CrPrt(IOut,IPrint,1,0,NPrim,V(InfCrd),V(IVars),
     $        V(IRVars))
            Write(IOut,9200) IPrim
            Call GauErr('CrdDef: Filling ICrdSt, IPos > NCoord.')
            endIf
          Call IIPut(V(ICrdSt),IPos,IPrim)
          Call IIPut(V(ICrTmp),KSet,IPos+1)
          endIf
 910    Continue
      If(IPrint.ge.4.or.DEBUG) then
        Call IMtOuS(IOut,'CrdDef: Bot,NCrdSt=',0,V(NCrdSt),1,MaxSet,1,
     $    MaxSet)
        Call IMtOuS(IOut,'CrdDef: Bot,ICrdSt=',0,V(ICrdSt),NCoord,1,
     $    NCoord,1)
        endIf
      Return
      End
*Deck CrdDq0
      Subroutine CrdDq0(IOut,IPrint,ISet,Iq1,Iq1CNm,Iq2,NumIq2,MaxSet,
     $  NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $  IVars,RVars,NCrdSt,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to add entries necessary to fully differentiate
C     coordinate Iq1 with respect to each coordinate given in the array
C     Iq2.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     ISet        This input argument indicates to which coordinate set
C                 this derivative is added.
C     Iq1         This input argument gives the coordinate number that is
C                 to be differentiated.
C     Iq1CNm      This input argument gives the "active/propagating"
C                 coordinate number that is to be differentiated.
C     Iq2         This input array gives a list of coordinate numbers with
C                 which coordinate Iq1 is to be differentiated with respect
C                 to.
C     NumIq2      This input argument gives the number of elements sent in
C                 array Iq2.
C     MaxSet      The maximum number of definition sets allowed.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C     V           Dynamic memory array
C     MDV         Length of V in words
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Logical AddEq,DoFix,UsSmDq,DidSmt,Exists,OK,Depend
      Parameter (NInfCr=DEFMAXCOORDINFO)
      Dimension Iq2(NumIq2),InfCrd(NInfCr,*),IVars(*),RVars(*),
     $  NCrdSt(*),ITemp(MxI),RTemp(MxR),V(*)
 1000 Format(1x,'Enter CrdDq0.')
 1010 Format(/,1x,'Forming definition list for deriv of coordinate ',
     $  I10,' wrt coordinate ',I10,'.')
 1100 Format(1x,'CrdDq0: UsSmDq is FALSE.')
 2010 Format(1x,'CrdDq0: dq1/dq2 already exists...q1=',I10,' q2=',I10,
     $  ' IdqPrm=',I10)
 2100 Format(1x,'CrdDq0: Forming dq1/dq2 for q1=',I10,', q2=',I10,
     $  ', Primitive #',I10)
 5000 Format(1x,'CrdDq0: CrTree says derivative is 0. IqA=',I10,
     $  '  Iq2=',I10)
 9000 Format(1x,'CrdDq0: Problem detected...',/,
     $  9x,'IPrim=',I10,' IType=',I10,' IAct=',I5)
C
C     Pick-up operation control settings. For now we statically set AddEq
C     and DoFix.
C
      UsSmDq = .False.
      If(IPrint.ge.5.and..not.UsSmDq) Write(IOut,1100)
      Call IClear(MxI,ITemp)
      Call AClear(MxR,RTemp)
      AddEq = .False.
      DoFix = .True.
      Call CrStOp(IOut,IPrint,2,23,ITemp(1),RTemp(1),UsSmDq)
C
C     Initialize arguments and set pointers in V for calls to CrTree.
C
      NEval  = 0
      NEval1 = 0
      NEval2 = 0
      LPos1  = 0
      LPos2  = 0
      IndxAr = 1
      IScr   = IndxAr + InToWP(MPrim)
      LIScr  = MDV-IScr+1
      Call TstCor(IScr-1,MDV,'CrdDq0-1')
C
C     Loop over the elements of Iq2 and add the necessary elements to
C     differentiate Iq1.
C
      If(IPrint.ge.5) Write(IOut,1000)
      Do 100 Iq2Cur = 1,NumIq2
C
C       Start out by defining dq1/dq2. If this coordinate already exists,
C       then we just exit.
C
        If(IPrint.ge.5) Write(IOut,1010) Iq1,Iq2(Iq2Cur)
        IdqPrm = 0
        NPrmSt = NPrim+1
        IType = 500
        IAct = 1000
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = Iq1
        ITemp(3) = Iq2(Iq2Cur)
        ITemp(4) = Iq1CNm
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IdqPrm,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq0: Failure adding derivative.')
        If(NPrim.lt.NPrmSt) then
          If(IPrint.ge.4) Write(IOut,2010) Iq1,Iq2(Iq2Cur),IdqPrm
          Return
          endIf
        IqA = Iq1
        IPrimA = IdqPrm
C
C       Now, look through the most recently added entries in the coordinate
C       definitions data structure and add new derivative definitions until
C       all of the necessary derivatives have been added to fully define
C       dq1/dq2. If 'smart' differentiation has been requested, then check
C       to see if one of the built-in coordinate types is being
C       differentiated, and then use smart differentiation as appropriate.
C
 200    If(IPrint.ge.5) Write(IOut,2100) IqA,Iq2(Iq2Cur),IPrimA
        DidSmt = .False.
        If(UsSmDq) then
          IActA = InfCrd(InfAct,IPrimA)
          IIVarA = InfCrd(InfIVr,IPrimA)
          IdPrim = IVars(IIVarA+1)
          IdType = InfCrd(InfTyp,IdPrim)
          IdIVar = InfCrd(InfIVr,IdPrim)
          If(IdType.eq.1.and.Iq2(Iq2Cur).lt.0) then
            DidSmt = .True.
            Call CrDfBd(IOut,IPrint,1,.False.,IVars(IdIVar+1),
     $        IVars(IdIVar+2),Iq2(Iq2Cur),0,-1,ISet,MaxSet,INwPrm,NPrim,
     $        NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $        IVars,RVars,NCrdSt)
            IVars(IIVarA) = INwPrm
            If(Abs(IActA).ge.1000) then
              INwAct = Abs(IActA)-1000
              INwAct = Sign(INwAct,IActA)
              InfCrd(InfAct,IPrimA) = INwAct
              endIf
          else if(IdType.eq.2.and.Iq2(Iq2Cur).lt.0) then
            DidSmt = .True.
            Call CrDfAn(IOut,IPrint,1,IVars(IdIVar+1),IVars(IdIVar+2),
     $        IVars(IdIVar+3),Iq2(Iq2Cur),0,-1,ISet,MaxSet,INwPrm,NPrim,
     $        NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $        IVars,RVars,NCrdSt)
            IVars(IIVarA) = INwPrm
            If(Abs(IActA).ge.1000) then
              INwAct = Abs(IActA)-1000
              INwAct = Sign(INwAct,IActA)
              InfCrd(InfAct,IPrimA) = INwAct
              endIf
          else if(IdType.eq.3.and.Iq2(Iq2Cur).lt.0) then
            DidSmt = .True.
            Call CrDfDi(IOut,IPrint,1,IVars(IdIVar+1),IVars(IdIVar+2),
     $        IVars(IdIVar+3),IVars(IdIVar+4),Iq2(Iq2Cur),-1,ISet,
     $        MaxSet,INwPrm,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,
     $        MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
            IVars(IIVarA) = INwPrm
            If(Abs(IActA).ge.1000) then
              INwAct = Abs(IActA)-1000
              INwAct = Sign(INwAct,IActA)
              InfCrd(InfAct,IPrimA) = INwAct
              endIf
          else if(IdType.eq.100.and.Iq2(Iq2Cur).lt.0) then
            DidSmt = .True.
            Call CrDfVV(IOut,IPrint,1,.False.,.True.,.True.,.True.,
     $        .False.,IVars(IdIVar),IVars(IdIVar+1),IVars(IdIVar+2),
     $        IVars(IdIVar+3),Iq2(Iq2Cur),0,ISet,MaxSet,IPrmVV,INwPrm,
     $        IVVxy,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $        MRVars,InfCrd,IVars,RVars,NCrdSt)
            IVars(IIVarA) = INwPrm
            If(Abs(IActA).ge.1000) then
              INwAct = Abs(IActA)-1000
              INwAct = Sign(INwAct,IActA)
              InfCrd(InfAct,IPrimA) = INwAct
              endIf
          else if(IdType.eq.500.and.Iq2(Iq2Cur).lt.0) then
            JdPrim = IVars(IdIVar+1)
            JdType = InfCrd(InfTyp,JdPrim)
            JdIVar = InfCrd(InfIVr,JdPrim)
            If(JdType.eq.1) then
              DidSmt = .True.
              Call CrDfBd(IOut,IPrint,2,.False.,IVars(JdIVar+1),
     $          IVars(JdIVar+2),Iq2(Iq2Cur),IVars(IdIVar+2),-1,ISet,
     $          MaxSet,INwPrm,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,
     $          MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
              IVars(IIVarA) = INwPrm
              If(Abs(IActA).ge.1000) then
                INwAct = Abs(IActA)-1000
                INwAct = Sign(INwAct,IActA)
                InfCrd(InfAct,IPrimA) = INwAct
                endIf
              endIf
            If(.False..and.JdType.eq.2) then
              DidSmt = .True.
              Call CrDfAn(IOut,IPrint,2,IVars(JdIVar+1),IVars(JdIVar+2),
     $          IVars(JdIVar+3),Iq2(Iq2Cur),IVars(IdIVar+2),-1,ISet,
     $          MaxSet,INwPrm,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,
     $          MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
              IVars(IIVarA) = INwPrm
              If(Abs(IActA).ge.1000) then
                INwAct = Abs(IActA)-1000
                INwAct = Sign(INwAct,IActA)
                InfCrd(InfAct,IPrimA) = INwAct
                endIf
              endIf
C
            If(JdType.eq.100) then
              DidSmt = .True.
              Call CrDfVV(IOut,IPrint,2,.False.,.True.,.True.,.True.,
     $          .False.,IVars(JdIVar),IVars(JdIVar+1),IVars(JdIVar+2),
     $          IVars(JdIVar+3),Iq2(Iq2Cur),IVars(IdIVar+2),ISet,MaxSet,
     $          IPrmVV,IVVx,INwPrm,NPrim,NCoord,LIVars,LRVars,MPrim,
     $          MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
              IVars(IIVarA) = INwPrm
              If(Abs(IActA).ge.1000) then
                INwAct = Abs(IActA)-1000
                INwAct = Sign(INwAct,IActA)
                InfCrd(InfAct,IPrimA) = INwAct
                endIf
              endIf
            endIf
          endIf
C
C       If the current derivative hasn't been solved using the "smart"
C       codes above, then use the general routine CrdDq1.
C
        If(.not.DidSmt) then
          Depend = .True.
          If(Iq2(Iq2Cur).lt.0) then
            Call CrTree(IOut,IPrint,15,IqA,Iq2(Iq2Cur),ISet,NPrim,
     $        NDepMx,V(IndxAr),V(IScr),V(IScr),NEval,NEval1,NEval2,
     $        LPos1,LPos2,Depend,InfCrd,IVars,V(IScr),LIScr)
            endIf
          If(Depend) then
            Call CrdDq1(IOut,IPrint,ISet,IqA,Iq2(Iq2Cur),IPrimA,MaxSet,
     $        NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,
     $        InfCrd,IVars,RVars,NCrdSt)
          else
            If(IPrint.ge.5.and..not.Depend) Write(IOut,5000) IqA,
     $        Iq2(Iq2Cur)
            Call CrDfCn(IOut,IPrint,0,ISet,MaxSet,INwPrm,NPrim,NCoord,
     $        LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,
     $        RVars,NCrdSt)
            IVars(IIVarA) = INwPrm
            If(Abs(IActA).ge.1000) then
              INwAct = Abs(IActA)-1000
              INwAct = Sign(INwAct,IActA)
              InfCrd(InfAct,IPrimA) = INwAct
              endIf
            endIf
          endIf
C
C       Check to see if there are any remaining derivatives to solve.
C
        IPrimA = 0
        Do 300 IPrim = NPrmSt,NPrim
          IType = InfCrd(InfTyp,IPrim)
          IAct  = InfCrd(InfAct,IPrim)
          IAct1 = ABS(Mod(IAct,10000))/1000
          IIVar = InfCrd(InfIVr,IPrim)
          If(IAct1.eq.1) then
            If(IType.ne.500) then
              Write(IOut,9000) IPrim,IType,IAct
              Call GauErr('CrdDq0: IType should not have IAct.ge.1000')
              endIf
            IPrimA = IPrim
            IqA = IVars(IIVar+1)
            Goto 200
            endIf
 300      Continue
 100    Continue
      Return
      End
*Deck CrdDq1
      Subroutine CrdDq1(IOut,IPrint,ISet,Iq1,Iq2,Idq,MaxSet,NPrim,
     $  NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,
     $  RVars,NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to build a list of coordinate definitions for
C     differentiating one internal coordinate (number Iq1) wrt coordinate
C     Iq2. The list of coordiante definitions added to the coordinate
C     definition data structure do NOT fully differentiate coordinate Iq1.
C     Instead, it applies the chain rule to differentiate Iq1 wrt other
C     internal coordinates. When a derivative of another coordinates is
C     needed but is not yet present in the definitions list, its definition
C     is added with the "incomplete definition" action value of -999. If
C     Iq1=Iq2 (or if coordinate Iq1 is fully defined in terms of Iq2) then
C     the calling program is done differentiating.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     ISet        This input argument indicates to which coordinate set
C                 this derivative is added.
C     Iq1         This input argument gives the coordinate number that is
C                 to be differentiated.
C     Iq2         This input argument gives the coordinate number for which
C                 coordinate Iq1 is to be differentiated with respect to.
C     Idq         The coordinate number of the derivative of coorindate
C                 Iq1. This value needed for this routine to finalize the
C                 actual derivative's definition in the data structure.
C     MaxSet      The maximum number of definition sets allowed.
C     NPrim       The current number of primitive internal coordinates.
C                 NPrim is incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Logical AddEq,DoFix,Exists,OK
      Parameter (NInfCr=DEFMAXCOORDINFO)
      Dimension ITemp(MxI),RTemp(MxR),InfCrd(NInfCr,*),IVars(*),
     $  RVars(*),NCrdSt(*)
 1000 Format(1x,'Enter CrdDq1.')
 1010 Format(3x,'Differentiating coordinate #',I10,' wrt coordinate #',
     $  I10)
 9000 Format(1x,'CrdDq1: Iq1=',I10,' Idq=',I10,' IdqTyp=',I10)
 9010 Format(1x,'CrdDq1: Iq1=',I10,' Idq=',I10,' IdqAct=',I10)
 9020 Format(1x,'CrdDq1: Iq1Typ=',I10)
C
C     Do some initialization. This includes checking on the action flag for
C     coordinate Idq and re-setting it so that the 1xxx flag is removed.
C
      Call IClear(MxI,ITemp)
      Call AClear(MxR,RTemp)
      If(IPrint.ge.5) then
        Write(IOut,1000)
        Write(IOut,1010) Iq1,Iq2
        endIf
      If(Iq1.eq.0) then
        Call CrPrt(IOut,IPrint,0,0,NPrim,InfCrd,IVars,RVars)
        Call GauErr('CrdDq1: Iq1=0.')
        endIf
      If(Iq2.eq.0) then
        Call CrPrt(IOut,IPrint,0,0,NPrim,InfCrd,IVars,RVars)
        Call GauErr('CrdDq1: Iq2=0.')
        endIf
      IdqTyp = InfCrd(InfTyp,Idq)
      IdqAct = InfCrd(InfAct,Idq)
      IdqIVr = InfCrd(InfIVr,Idq)
      Iq1Typ = InfCrd(InfTyp,Iq1)
      Iq1IVr = InfCrd(InfIVr,Iq1)
      Iq1RVr = InfCrd(InfRVr,Iq1)
      If(Iq2.lt.0) then
        IdXYZ = Mod(-Iq2,10)
        IdXYZ = -IdXYZ
        IdAtom = -Iq2/10
      else
        IdAtom = 0
        IdXYZ = 0
        endIf
      If(IdqTyp.ne.500) then
        Call CrPrt(IOut,IPrint,0,0,NPrim,InfCrd,IVars,RVars)
        Write(IOut,9000) Iq1,Idq,IdqTyp
        Call GauErr('CrdDq1: Idq has improper type flag.')
        endIf
      If(Abs(IdqAct).ge.1000) then
        INwAct = Abs(IdqAct)-1000
        INwAct = Sign(INwAct,IdqAct)
        InfCrd(InfAct,Idq) = INwAct
      else
        Call CrPrt(IOut,IPrint,0,0,NPrim,InfCrd,IVars,RVars)
        Write(IOut,9010) Iq1,Idq,IdqAct
        Call GauErr('CrdDq1: Derivative has improper action flag.')
        endIf
C
C     For now we statically set AddEq and DoFix.
C
      AddEq = .False.
      DoFix = .True.
C
C     Test for the special case where Iq1=Iq2.
C
      If(Iq1.eq.Iq2) then
        Call CrDfCn(IOut,IPrint,1,ISet,MaxSet,IPrimA,NPrim,NCoord,
     $    LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $    NCrdSt)
        IVars(IdqIVr) = IPrimA
        Return
        endIf
C
C     Differentiate Iq1 according to its coordinate type, Iq1Typ. Note that
C     if we have made it to this block of code, then Iq1.ne.Iq2.
C
C     q=Cartesian coordinate
      If(Iq1Typ.eq.0) then
        If(IdAtom.eq.IVars(Iq1IVr).and.IdXYZ.eq.IVars(Iq1IVr+1)) then
          Call CrDfCn(IOut,IPrint,1,ISet,MaxSet,IPrimA,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
        else
          Call CrDfCn(IOut,IPrint,0,ISet,MaxSet,IPrimA,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
          endIf
        IVars(IdqIVr) = IPrimA
C
C     q=Bond stretch term
      else if(Iq1Typ.eq.1) then
        IqA = IVars(Iq1IVr)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimA
C
C     q=Angle bend term
      else if(Iq1Typ.eq.2) then
        IqA = IVars(Iq1IVr)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimA
C
C     q=Torsion term
      else if(Iq1Typ.eq.3) then
        IqA = IVars(Iq1IVr)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimA
C
C     q=q(A)+q(B)
      else if(Iq1Typ.eq.50) then
        IqA = IVars(Iq1IVr)
        IqB = IVars(Iq1IVr+1)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqB
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimB,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 50
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimA
        ITemp(2) = IPrimB
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimC,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimC
C
C     q=q(A)-q(B)
      else if(Iq1Typ.eq.51) then
        IqA = IVars(Iq1IVr)
        IqB = IVars(Iq1IVr+1)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqB
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimB,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 51
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimA
        ITemp(2) = IPrimB
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimC,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimC
C
C     q=q(A)*q(B)
      else if(Iq1Typ.eq.ICtPrd) then
        IqA = IVars(Iq1IVr)
        IqB = IVars(Iq1IVr+1)
        ITypA = InfCrd(InfTyp,IqA)
        ITypB = InfCrd(InfTyp,IqB)
        If(ITypA.ne.58) then
          IType = 500
          IAct = -1001
          NIVars = 5
          NRVars = 0
          ITemp(1) = 0
          ITemp(2) = IqA
          ITemp(3) = Iq2
          ITemp(4) = -1
          ITemp(5) = -1
          Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
          IType = 52
          IAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IqB
          ITemp(2) = IPrimA
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimB,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
          endIf
C
        If(ITypB.ne.58) then
          IType = 500
          IAct = -1001
          NIVars = 5
          NRVars = 0
          ITemp(1) = 0
          ITemp(2) = IqB
          ITemp(3) = Iq2
          ITemp(4) = -1
          ITemp(5) = -1
          Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimC,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
          IType = 52
          IAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IqA
          ITemp(2) = IPrimC
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimD,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
          endIf
C
        If(ITypA.ne.58.and.ITypB.ne.58) then
          IType = 50
          IAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPrimB
          ITemp(2) = IPrimD
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimE,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
          IVars(IdqIVr) = IPrimE
        else if(ITypA.eq.58.and.ITypB.eq.58) then
          Call CrDfCn(IOut,IPrint,0,ISet,MaxSet,IPrimA,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
          IVars(IdqIVr) = IPrimA
        else if(ITypA.eq.58) then
          IVars(IdqIVr) = IPrimD
        else
          IVars(IdqIVr) = IPrimB
          endIf
C
C     q=q(A)/q(B)
      else if(Iq1Typ.eq.ICTQuo) then
        IqA = IVars(Iq1IVr)
        IqB = IVars(Iq1IVr+1)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqB
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimB,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 52
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IqA
        ITemp(2) = IPrimB
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimC,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 52
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IqB
        ITemp(2) = IPrimA
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimD,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 51
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimD
        ITemp(2) = IPrimC
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimE,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 63
        IAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IqB
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimF,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = ICTQuo
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimE
        ITemp(2) = IPrimF
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimG,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimG
C
C     q=S*q(A)
      else if(Iq1Typ.eq.54) then
        IqA = IVars(Iq1IVr)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 54
        IAct = -1
        NIVars = 1
        NRVars = 1
        ITemp(1) = IPrimA
        RTemp(1) = RVars(Iq1RVr)
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimB,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimB
C
C     q=q(A)
      else if(Iq1Typ.eq.56) then
        IqA = IVars(Iq1IVr)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimA
C
C     NYI
      else if(Iq1Typ.eq.57) then
        Call GauErr('CrdDq1: Iq1Typ 57 NYI')
C
C     q=C
      else if(Iq1Typ.eq.58) then
        Call CrDfCn(IOut,IPrint,0,ISet,MaxSet,IPrimA,NPrim,NCoord,
     $    LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $    NCrdSt)
        IVars(IdqIVr) = IPrimA
      else if(Iq1Typ.eq.60) then
        Call GauErr('CrdDq1: Iq1Typ 60 NYI')
C
C     q = q(A)**-1
      else if(Iq1Typ.eq.61) then
        IqA = IVars(Iq1IVr)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        Call CrDfCn(IOut,IPrint,-1,ISet,MaxSet,IPrimB,NPrim,NCoord,
     $    LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $    NCrdSt)
C
        IType = 63
        IAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = Iq1
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimC,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 52
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimB
        ITemp(2) = IPrimC
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimD,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 52
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimA
        ITemp(2) = IPrimD
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimE,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimE
C
C     q = SQRT[q(A)]
      else if(Iq1Typ.eq.62) then
        IqA = IVars(Iq1IVr)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 61
        IAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = Iq1
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimB,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        Call CrDfCn(IOut,IPrint,11,ISet,MaxSet,IPrimC,NPrim,NCoord,
     $    LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $    NCrdSt)
C
        IType = 52
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimB
        ITemp(2) = IPrimC
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimB,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 52
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimA
        ITemp(2) = IPrimB
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimC,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimC
C
C     q = q(A)**2
      else if(Iq1Typ.eq.ICTSq) then
        IqA = IVars(Iq1IVr)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        Call CrDfCn(IOut,IPrint,2,ISet,MaxSet,IPrimB,NPrim,NCoord,
     $    LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $    NCrdSt)
C
        IType = 52
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimA
        ITemp(2) = IPrimB
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimC,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 52
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimC
        ITemp(2) = IqA
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimD,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimD
C
C     q=cos(q(A))
      else if(Iq1Typ.eq.70) then
        IqA = IVars(Iq1IVr)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        Call CrDfCn(IOut,IPrint,-1,ISet,MaxSet,IPrimB,NPrim,NCoord,
     $    LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $    NCrdSt)
C
        IType = 71
        IAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IqA
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimC,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 52
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimA
        ITemp(1) = IPrimB
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimD,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 52
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimD
        ITemp(1) = IPrimC
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimE,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimE
C
C     q=sin(q(A))
      else if(Iq1Typ.eq.71) then
        IqA = IVars(Iq1IVr)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 70
        IAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IqA
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimB,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 52
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimA
        ITemp(1) = IPrimB
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimC,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimC
C
C     q=tan(q(A))
      else if(Iq1Typ.eq.72) then
        IqA = IVars(Iq1IVr)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 70
        IAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IqA
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimB,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 63
        IAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IPrimB
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimC,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 61
        IAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IPrimC
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimD,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 52
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimA
        ITemp(1) = IPrimD
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimE,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimE
C
C     q=arccos(q(A))
      else if(Iq1Typ.eq.ICTACos) then
        IqA = IVars(Iq1IVr)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        Call CrDfCn(IOut,IPrint,-1,ISet,MaxSet,IPrimB,NPrim,NCoord,
     $    LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $    NCrdSt)
C
        Call CrDfCn(IOut,IPrint,1,ISet,MaxSet,IPrimC,NPrim,NCoord,
     $    LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $    NCrdSt)
C
        IType = 63
        IAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IqA
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimD,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 51
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimC
        ITemp(2) = IPrimD
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimE,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 62
        IAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IPrimE
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimF,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 61
        IAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IPrimF
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimG,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 52
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimB
        ITemp(2) = IPrimG
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimH,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
C
        IType = 52
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrimA
        ITemp(2) = IPrimH
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimI,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimI
      else if(Iq1Typ.eq.74) then
        Call GauErr('CrdDq1: Iq1Typ 74 NYI')
      else if(Iq1Typ.eq.75) then
        Call GauErr('CrdDq1: Iq1Typ 75 NYI')
      else if(Iq1Typ.eq.ICTSign) then
        Call GauErr('CrdDq1: Iq1Typ ICTSign NYI')
C
C     Conditional (less-than) function. For q=f(q(A),c,z), q=q(A)+z for
C     q(A) < c; otherwise, q=q(A).
      else if(Iq1Typ.eq.ICTCdL) then
        IqA = IVars(Iq1IVr)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimA
C
C     q = [Cart(:,IAtomA)-Cart(:,IAtomB)].[Cart(:,IAtomC)-Cart(:,IAtomD)]
      else if(Iq1Typ.eq.100) then
        IqA = IVars(Iq1IVr)
        IqB = IVars(Iq1IVr+1)
        IqC = IVars(Iq1IVr+2)
        IqD = IVars(Iq1IVr+3)
        IDiff1 = 0
        IDiff2 = 0
        If(IqC.eq.IdAtom) then
          IDiff1 = 1
        else if(IqD.eq.IdAtom) then
          IDiff1 = -1
          endIf
        If(IqA.eq.IdAtom) then
          IDiff2 = 1
        else if(IqB.eq.IdAtom) then
          IDiff2 = -1
          endIf
        IPrimA = 0
        IPrimB = 0
        If(IDiff1.ne.0) then
          IType = 110
          IAct = -1
          NIVars = 3
          NRVars = 0
          ITemp(1) = IdXYZ
          ITemp(2) = IqA
          ITemp(3) = IqB
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimA,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
          If(IDiff1.eq.-1) Call CrDfM1(IOut,IPrint,IPrimA,ISet,MaxSet,
     $      NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,
     $      InfCrd,IVars,RVars,NCrdSt)
          endIf
        If(IDiff2.ne.0) then
          IType = 110
          IAct = -1
          NIVars = 3
          NRVars = 0
          ITemp(1) = IdXYZ
          ITemp(2) = IqC
          ITemp(3) = IqD
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimB,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
          If(IDiff2.eq.-1) Call CrDfM1(IOut,IPrint,IPrimB,ISet,MaxSet,
     $      NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,
     $      InfCrd,IVars,RVars,NCrdSt)
          endIf
        If(IPrimA.ne.0.and.IPrimB.ne.0) then
          IType = 50
          IAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPrimA
          ITemp(2) = IPrimB
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrimC,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        else if(IPrimA.eq.0.and.IPrimB.eq.0) then
          Call CrDfCn(IOut,IPrint,0,ISet,MaxSet,IPrimC,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
        else if(IPrimA.ne.0) then
          IPrimC = IPrimA
        else if(IPrimB.ne.0) then
          IPrimC = IPrimB
          endIf
        IVars(IdqIVr) = IPrimC
C
C     q = Cart(IqXYZ,IqA)-Cart(IqXYZ,IqB)
      else if(Iq1Typ.eq.110) then
        IqXYZ = IVars(Iq1IVr)
        IqA = IVars(Iq1IVr+1)
        IqB = IVars(Iq1IVr+2)
        If(IdXYZ.eq.IqXYZ.and.IdAtom.eq.IqA) then
          Call CrDfCn(IOut,IPrint,1,ISet,MaxSet,IPrimA,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
        else if(IdXYZ.eq.IqXYZ.and.IdAtom.eq.IqB) then
          Call CrDfCn(IOut,IPrint,-1,ISet,MaxSet,IPrimA,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
        else
          Call CrDfCn(IOut,IPrint,0,ISet,MaxSet,IPrimA,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
          endIf
        IVars(IdqIVr) = IPrimA
C
C     q = 2*Cart(IqXYZ,IqA)-Cart(IqXYZ,IqB)-Cart(IqXYZ,IqC)
      else if(Iq1Typ.eq.111) then
        IqXYZ = IVars(Iq1IVr)
        IqA = IVars(Iq1IVr+1)
        IqB = IVars(Iq1IVr+2)
        IqC = IVars(Iq1IVr+3)
        If(IdXYZ.eq.IqXYZ.and.IdAtom.eq.IqA) then
          Call CrDfCn(IOut,IPrint,2,ISet,MaxSet,IPrimA,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
        else if(IdXYZ.eq.IqXYZ.and.
     $    (IdAtom.eq.IqB.or.IdAtom.eq.IqC)) then
          Call CrDfCn(IOut,IPrint,-1,ISet,MaxSet,IPrimA,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
        else
          Call CrDfCn(IOut,IPrint,0,ISet,MaxSet,IPrimA,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
          endIf
        IVars(IdqIVr) = IPrimA
      else if(Iq1Typ.eq.500) then
        IqA = IVars(Iq1IVr)
        IType = 500
        IAct = -1001
        NIVars = 5
        NRVars = 0
        ITemp(1) = 0
        ITemp(2) = IqA
        ITemp(3) = Iq2
        ITemp(4) = -1
        ITemp(5) = -1
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrimA,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrdDq1: Failure adding derivative.')
        IVars(IdqIVr) = IPrimA
      else
        Write(IOut,9020) Iq1Typ
        Call GauErr('CrdDq1: Unknown Iq1Typ.')
        endIf
      Return
      End
*Deck CrDf01
      Subroutine CrDf01(IOut,IPrint,IRdFmt,DoInvR,ISet,MaxSet,NPrim,
     $  NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,
     $  RVars,NCrdSt,ICrdSt,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     This routine adds coordinate definitions read from the input file to
C     the coordinate definition data structure, and is part of the CrDfxx
C     group of routines. Two different formats (indicated by input argument
C     IRdFmt) can be read from the input file. These two formats are meant
C     to make it straight forward for defining coordinate definitions and
C     for definition derivatives of previously defined definitions wrt
C     Cartesian coordinates.
C
C     When IRdFmt=2, the format for defining an internal coordinate in the
C     input file is:
C
C           IType IAction NIVars NRVars [IVar1,IVar2,...]  [RVar1,RVar2,...]
C
C     where IType and IAction have the meanings in array InfCrd defined in
C     Routine CrdDef. NIVars and NRVars give the number of integer and real
C     values that will be read in. Then, the integer and real values are
C     read and loaded into arrays IVars and RVals.
C
C     When IRdFmt=3, the format for defining an internal coordinate in the
C     input file is:
C
C           IAtom1  IAtom2  [IAtom3  [IAtom4]]
C
C     where IAtom1, IAtom2, IAtom3, and IAtom4 are used to define a
C     convention internal coordinate. IAtom3 and IAtom4 are options
C     arguments. Following the standard convention, if only IAtom1 and
C     IAtom2 are given the internal coordinate added is the bond stretch
C     between these two centers; if IAtom3 is also given, then the
C     angle/bend coordinate with IAtom2 at the center is added; if IAtom4
C     is also given then the dihedral/torsion coordinate is added. Note
C     that this mode also adds all of the necessary primitives definitions
C     as well. All of the defined internal coordinates are flagged as
C     active coordinates.
C
C     When IRdFmt=8, the format for defining an internal coordinate, which
C     is the derivative of a previously defined coordinate definition with
C     respect to a set of atom-center Cartesian coordinates, in the input
C     file is:
C
C           IdQ   IQSet  IAtom  Ixyz
C
C     where IdQ is the "active" coordinate definition number of the
C     coordinate being differentiated. IQSet is the definition set number
C     to which IdQ belongs. IAtom and Ixyz give the atomic center and
C     Cartesian coordinate with respect to which IdQ is being
C     differentiated. Ixyz is given as 1, 2, or 3 for x, y, or z.
C
C     When IRdFmt=9, the format for defining an internal coordinate, which
C     is the derivative of a previously defined coordinate definition with
C     respect to a set of atom-center Cartesian coordinates, in the input
C     file is:
C
C           IdQ   IQSet  NAtoms  [IAtom1,IAtom2,...IAtomN]
C
C     where IdQ is the "active" coordinate definition number of the
C     coordinate being differentiated. IQSet is the definition set number
C     to which IdQ belongs. NAtoms is an integer giving the number of
C     atomic centers for which IdQ should be differentiated. The list
C     IAtom1-IAtomN gives these NAtoms atomic center numbers.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     IRdFmt      An input integer flag indicating which reading format is
C                 used. See the above discussion for choices.
C     DoInvR      A logical input argument that indicates whether or not
C                 bond stretch terms should be built as inverse bond
C                 stretch definitions.
C     ISet        The set number to use when adding the set of definitions
C                 we are about to read.
C     MaxSet      The maximum number of definition sets allowed.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C     V           Dynamic memory array
C     MDV         Length of V in words
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Logical AddEq,DoFix,Exists,OK,DoInvR
      Parameter (MxDx=3)
      Dimension ITemp(MxI),RTemp(MxR),IdqLst(MxDx),InfCrd(*),IVars(*),
     $  RVars(*),NCrdSt(*),ICrdSt(*),V(MDV)
      Save Zero
      Data Zero/0.d0/
 1000 Format(1x,'CrDf01 is appending user defined coordinate to the ',
     $  'coordinate definition arrays.')
 9000 Format(1x,'CrDf01: IRdFmt=',I5)
C
C     For now we statically set AddEq and DoFix.
C
      Call IClear(MxI,ITemp)
      Call AClear(MxR,RTemp)
      AddEq = .True.
      DoFix = .True.
C
C     Read the definitions directly from the input file. Coordinates are
C     added one at a time until we hit a blank line delimiter.
C
C     IRdFmt=2 ... Read raw primitive definitions from the input file.
      If(IRdFmt.eq.2) then
C       Do While(IType.ne.-999)
 100    Call FFRead(IEOF)
        If(IEOF.ne.0)
     $    Call GauErr('Coord definitions missing from the input file.')
        IType = IFFGet(.True.,-999,.False.)
        If(IType.eq.-999) Goto 999
        IAct = IFFGet(.False.,0,.True.)
        NIVars = IFFGet(.False.,0,.True.)
        NRVars = IFFGet(.False.,0,.True.)
        Do 110 i = 1,NIVars
 110      ITemp(i) = IFFGet(.False.,0,.True.)
        Do 120 i = 1,NRVars
 120      RTemp(i) = RFFGet(.False.,.True.,Zero,.True.)
        Write(IOut,1000)
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,INwPrm,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDf01: After CrdAdd, OK=FALSE.')
        Goto 100
C
C     IRdFmt=3 ... Read typical valence coordinate definitions from the
C                  input file. All of the read-in coordinates are flagged
C                  as active coordinates.
      else if(IRdFmt.eq.3) then
 300    Call FFRead(IEOF)
        If(IEOF.ne.0)
     $    Call GauErr('Coord definitions missing from the input file.')
        IAtom1 = IFFGet(.True.,-999,.False.)
        If(IAtom1.eq.-999) Goto 999
        IAtom2 = IFFGet(.True.,-999,.False.)
        If(IAtom2.eq.-999)
     $    Call GauErr('CrDf01: Invalid IAtom2 found in the input file.')
        IAtom3 = IFFGet(.True.,0,.False.)
        If(IAtom3.lt.0)
     $    Call GauErr('CrDf01: Invalid IAtom3 found in the input file.')
        IAtom4 = IFFGet(.True.,0,.False.)
        If(IAtom4.lt.0)
     $    Call GauErr('CrDf01: Invalid IAtom4 found in the input file.')
        Call CrDfIC(IOut,IPrint,IAtom1,IAtom2,IAtom3,IAtom4,DoInvR,0,
     $    ISet,MaxSet,INwPrm,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,
     $    MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
        Goto 300
C
C     IRdFmt=8 ... Read a list of active coordinates and atomic coordinates
C                  over which to differentiate.
      else if(IRdFmt.eq.8) then
 800    Call FFRead(IEOF)
        If(IEOF.ne.0)
     $    Call GauErr('Coord definitions missing from the input file.')
        IdQX = IFFGet(.True.,-999,.False.)
        If(IdQX.eq.-999) Goto 999
        IQSet = IFFGet(.False.,0,.True.)
        IPos = IArSum(1,IQSet-1,NCrdSt)+IdQX
        IdQ = ICrdSt(IPos)
        IAtom = IFFGet(.False.,0,.True.)
        Ixyz = IFFGet(.False.,0,.True.)
        IdqLst(1) = -(10*IAtom+Ixyz)
        NumIq2 = 1
        Call CrdDq0(IOut,IPrint,ISet,IdQ,IdQX,IdqLst,NumIq2,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $    IVars,RVars,NCrdSt,V,MDV)
        Goto 800
C
C     IRdFmt=9 ... Read a list of active coordinates and atomic centers
C                  over which to differentiate.
      else if(IRdFmt.eq.9) then
 900    Call FFRead(IEOF)
        If(IEOF.ne.0)
     $    Call GauErr('Coord definitions missing from the input file.')
        IdQX = IFFGet(.True.,-999,.False.)
        If(IdQX.eq.-999) Goto 999
        IQSet = IFFGet(.False.,0,.True.)
        IPos = IArSum(1,IQSet-1,NCrdSt)+IdQX
        IdQ = ICrdSt(IPos)
        NAtoms = IFFGet(.False.,0,.True.)
        Do 910 i = 1,NAtoms
          IAtom = IFFGet(.False.,0,.True.)
          NumIq2 = 3
          IdqLst(1) = -(10*IAtom+1)
          IdqLst(2) = -(10*IAtom+2)
          IdqLst(3) = -(10*IAtom+3)
          Call CrdDq0(IOut,IPrint,ISet,IdQ,IdQX,IdqLst,NumIq2,MaxSet,
     $      NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,
     $      InfCrd,IVars,RVars,NCrdSt,V,MDV)
 910      Continue
        Goto 900
      else
        Write(IOut,9000) IRdFmt
        Call GauErr('CrDf01: Invalid IRdFmt.')
        endIf
 999  Return
      End
*Deck CrDf21
      Subroutine CrDf21(IOut,IPrint,DoInvR,NAtoms,ISet,MaxSet,NPrim,
     $  NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,
     $  RVars,NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine adds bond stretch and associated definitions to the
C     coordinate definition data structure, and is part of the CrDfxx group
C     of routines. If DoInvR is sent as .TRUE., then all inverse bond
C     stretch and associated definitions are added to the coordinate
C     definition data structure.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     DoInvR      A logical input argument indicated whether (TRUE) or not
C                 (FALSE) the inverse bond stretch derinition is desired.
C     NAtoms      The total number of atoms.
C     ISet        This input argument indicates to which coordinate set
C                 this derivative is added.
C     MaxSet      The maximum number of definition sets allowed.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
      Logical DoInvR
      Dimension InfCrd(*),IVars(*),RVars(*),NCrdSt(*)
 2000 Format(1x,'CrDf21 is appending R(',I5,',',I5,')')
 2002 Format(1x,'CrDf21 is appending 1/R(',I5,',',I5,')')
C
C     Wrap successive calls to Routine CrDfIC adding the appropriate bond
C     stretch defintions.
C
      Do 100 IAtom1 = 1,NAtoms-1
        Do 100 IAtom2 = IAtom1+1,NAtoms
          If(DoInvR) then
            Write(IOut,2002) IAtom1,IAtom2
          else
            Write(IOut,2000) IAtom1,IAtom2
            endIf
          Call CrDfIC(IOut,IPrint,IAtom1,IAtom2,0,0,DoInvR,0,ISet,
     $      MaxSet,INwPrm,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,
     $      MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
 100      Continue
      Return
      End
*Deck CrDf30
      Subroutine CrDf30(IOut,IPrint,DoInvR,NTRed,IZRed,ISet,MaxSet,
     $  NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $  IVars,RVars,NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine uses IZRed to add the same coordinate definitions
C     (together with necessary associated definitions) to the coordinate
C     definition data structure, and is part of the CrDfxx group of
C     routines. If DoInvR is sent as .TRUE., then all bond stretch
C     definitions are replaced by inverse bond stretches.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     DoInvR      A logical input argument indicated whether (TRUE) or not
C                 (FALSE) the inverse bond stretch derinition is desired.
C     NTRed       This input argument is the number of internal coordinates
C                 included in array IZRed.
C     IZRed       This input array gives the internal coordinates in the
C                 conventional form used in Links 101, 103, 716, and
C                 elsewhere. It is dimensioned (MxRdIn,NTRed).
C     ISet        This input argument indicates to which coordinate set
C                 this derivative is added.
C     MaxSet      The maximum number of definition sets allowed.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
      Logical DoInvR
      Parameter (MxRdIn=DEFMAXREDINDEX)
      Dimension IZRed(MxRdIn,NTRed),InfCrd(*),IVars(*),RVars(*),
     $  NCrdSt(*)
 2000 Format(1x,'CrDf30 is appending internal coordinate ',I10)
C
C     Wrap successive calls to Routine CrDfIC adding the appropriate terms.
C
      If(IPrint.ge.4) Call IMtOuS(IOut,'CrDf30: Input IZRed',0,IZRed,
     $  MxRdIn,NTRed,MxRdIn,NTRed)
      Do 100 ICoord = 1,NTRed
        If(IPrint.ge.4) Write(IOut,2000) ICoord
        Call CrDfIC(IOut,IPrint,IZRed(1,ICoord),IZRed(2,ICoord),
     $    IZRed(3,ICoord),IZRed(4,ICoord),DoInvR,0,ISet,MaxSet,INwPrm,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $    IVars,RVars,NCrdSt)
 100    Continue
      Return
      End
*Deck CrDfA1
      Subroutine CrDfA1(IOut,IPrint,IDeriv,IAtom1,IAtom2,IAtom3,Idx,
     $  ISet,MaxSet,IA123,IAx,IR12,IR32,NPrim,NCoord,LIVars,LRVars,
     $  MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine adds one standard internal coordinate and associated
C     definitions to the coordinate definition data structure, and is part
C     of the CrDfxx group of routines.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     IDeriv      An integer input argument indicating wheter the
C                 appropriate bond angle derivative should be added instead
C                 of the actual definition itself.
C     IAtom1-3    These are the numbers of the three atomic centers that
C                 define the bond angle being added. IAtom2 gives the
C                 central atom.
C     Idx         This input argument gives the Cartesian coordinate with
C                 respect to which A(IAtom1,IAtom2,IAtom3) is to be
C                 differentiated when IDeriv=1. If IDeriv=0, this argument
C                 is ignored.
C     ISet        This input argument indicates to which coordinate set
C                 this derivative is added.
C     MaxSet      The maximum number of definition sets allowed.
C     INwPrm      This output integer is filled with the primtive
C                 definition number of the final coordinate.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Logical AddEq,DoFix,Exists,OK,MVV,MR12x,MR32x,MR32y,M12y,Hv0
      Dimension ITemp(MxI),RTemp(MxR),InfCrd(*),IVars(*),RVars(*),
     $  NCrdSt(*)
      Save IAtSv1,IAtSv2,IAtSv3,IVVSv,IR12Sv,IR32Sv,IRRISv,IASv
      Data IAtSv1/0/,IAtSv2/0/,IAtSv3/0/IVVSv/0/,IR12Sv/0/,IR32Sv/0/,
     $  IRRISv/0/,IASv/0/
 9000 Format(1x,'CrDfA1: IAtom1=',I6,' IAtom2=',I6,' IAtom3=',I6,
     $  ' IdXAtm=',I6)
C
C     Initialize some flags. For now, we statically set AddEq and DoFix.
C
      Call IClear(MxI,ITemp)
      Call AClear(MxR,RTemp)
      AddEq = .False.
      DoFix = .True.
      Hv0 = .False.
      IA123 = 0
      IAx   = 0
C
C     Set-up IdXAtm, IdX1, and the local value of the derivative flag IDrv
C
      If(IDeriv.gt.0) then
        IdXAtm = -Idx/10
        IdX1 = Mod(-Idx,10)
        IdX1 = -IdX1
      else
        IdXAtm = 0
        IdX1 = 0
        endIf
      IDrv = Min(IDeriv,1)
C
C     Take care of the zero-order terms. The permanent terms coming from
C     this block of code include:
C           IVV    = V12.V32
C           IR12   = R(1,2)
C           IR32   = R(3,2)
C           IR12In = 1/R(1,2)
C           IR32In = 1/R(3,2)
C           IRRInv = IR12In*IR32In
C           IA123  = IRRInv*IVV
C
      If(IAtom1.eq.IAtSv1.and.IAtom2.eq.IAtSv2.and.IAtom3.eq.IAtSv3)
     $  then
        IVV    = IVVSv
        IR12   = IR12Sv
        IR32   = IR32Sv
        IRRInv = IRRISv
        IA123  = IASv
        Hv0 = .True.
        endIf
      Call CrDfVV(IOut,IPrint,IDrv,.False.,.False.,.True.,.True.,Hv0,
     $  IAtom1,IAtom2,IAtom3,IAtom2,Idx,0,ISet,MaxSet,IVV,IVVx,IVVxy,
     $  NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $  IVars,RVars,NCrdSt)
      MVV = IVV.lt.0
      If(MVV) IVV = -IVV
C
      IDrv12 = 0
      IDrv32 = 0
      If(IDrv.gt.0.and.(IdXAtm.eq.IAtom1.or.IdXAtm.eq.IAtom2))
     $  IDrv12 = 1
      If(IDrv.gt.0.and.(IdXAtm.eq.IAtom3.or.IdXAtm.eq.IAtom2))
     $  IDrv32 = 1
      Call CrDfB2(IOut,IPrint,IDrv12,.True.,.True.,IAtom1,IAtom2,Idx,0,
     $  ISet,MaxSet,Ivv12,IR12,Ivv12x,MR12x,IR12x,Ivv12y,M12y,IR12y,
     $  IdXdY,IBRxy,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $  MRVars,InfCrd,IVars,RVars,NCrdSt)
      Call CrDfB2(IOut,IPrint,IDrv32,.True.,.True.,IAtom3,IAtom2,Idx,0,
     $  ISet,MaxSet,Ivv32,IR32,Ivv32x,MR32x,IR32x,Ivv32y,MR32y,IR32y,
     $  IdXdY,IBRxy,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $  MRVars,InfCrd,IVars,RVars,NCrdSt)
C
      If(.not.Hv0) then
        IType = ICTInv
        IAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IR12
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IR12In,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfA1: Failure adding definition.')
        IType = ICTInv
        IAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IR32
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IR32In,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfA1: Failure adding definition.')
        IType = ICTPrd
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IR12In
        ITemp(2) = IR32In
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IRRInv,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfA1: Failure adding definition.')
        IType = ICTPrd
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IRRInv
        ITemp(2) = IVV
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IA123T,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfA1: Failure adding definition.')
        IA123 = IA123T
        If(MVV) Call CrDfM1(IOut,IPrint,IA123,ISet,MaxSet,NPrim,NCoord,
     $    LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $    NCrdSt)
        endIf
C
C     Save the zero-order terms to avoid re-defining them if they can be
C     used in the next call.
C
      If(.not.Hv0) then
        IAtSv1 = IAtom1
        IAtSv2 = IAtom2
        IAtSv3 = IAtom3
        IVVSv  = IVV
        IR12Sv = IR12
        IR32Sv = IR32
        IRRISv = IRRInv
        IASv   = IA123
        endIf
C
C     Take care of the first-order terms. The permanent terms coming from this
C     block of code include:
C           IVVx   = (V12.V32)x   Note: this term is actually built in the
C                                       zero-order block of code above in
C                                       the call to CrDfVV.
C           IRRx   = [R(1,2)*R(3,2)]x
C           IAx    = IRRInv*(IVVx-IA123*IRRx)
C
      If(IDeriv.ge.1) then
        If(MR12x) Call CrDfM1(IOut,IPrint,IR12x,ISet,MaxSet,NPrim,
     $    NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,
     $    RVars,NCrdSt)
        If(MR32x) Call CrDfM1(IOut,IPrint,IR32x,ISet,MaxSet,NPrim,
     $    NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,
     $    RVars,NCrdSt)
        If(IdXAtm.eq.IAtom1.or.IdXAtm.eq.IAtom2) then
          IType = ICTPrd
          IAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IR12x
          ITemp(2) = IR32
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK) Call GauErr('CrDfA1: Failure adding definition.')
          endIf
        If(IdXAtm.eq.IAtom3.or.IdXAtm.eq.IAtom2) then
          IType = ICTPrd
          IAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IR12
          ITemp(2) = IR32x
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp2,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK) Call GauErr('CrDfA1: Failure adding definition.')
          endIf
        If(IdXAtm.eq.IAtom2) then
          IType = ICTSum
          IAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp1
          ITemp(2) = IPTmp2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IRRx,OK,MaxSet,
     $      NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $      IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,
     $      NCrdSt)
          If(.not.OK) Call GauErr('CrDfA1: Failure adding definition.')
        else if(IdXAtm.eq.IAtom1) then
          IRRx = IPTmp1
        else if(IdXAtm.eq.IAtom3) then
          IRRx = IPTmp2
        else
          Write(IOut,9000) IAtom1,IAtom2,IAtom3,IdXAtm
          Call GauErr('CrDfA1: Logic Error A')
          endIf
C
        IType = ICTPrd
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IA123
        ITemp(2) = IRRx
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfA1: Failure adding definition.')
        IType = ICTDif
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IVVx
        ITemp(2) = IPTmp1
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp2,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfA1: Failure adding definition.')
        IType = ICTPrd
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IRRInv
        ITemp(2) = IPTmp2
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IAx,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfA1: Failure adding definition.')
        endIf
      Return
      End
*Deck CrDfA2
      Subroutine CrDfA2(IOut,IPrint,IDeriv,IAtom1,IAtom2,IAtom3,Idx,Idy,
     $  ISet,MaxSet,IA123,IAx,IR12,IR32,NPrim,NCoord,LIVars,LRVars,
     $  MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine adds one standard internal coordinate and associated
C     definitions to the coordinate definition data structure, and is part
C     of the CrDfxx group of routines.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     IDeriv      An integer input argument indicating wheter the
C                 appropriate bond angle derivative should be added instead
C                 of the actual definition itself.
C     IAtom1-3    These are the numbers of the three atomic centers that
C                 define the bond angle being added. IAtom2 gives the
C                 central atom.
C     Idx         This input argument gives the Cartesian coordinate with
C                 respect to which A(IAtom1,IAtom2,IAtom3) is to be
C                 differentiated when IDeriv=1. If IDeriv=0, this argument
C                 is ignored.
C     ISet        This input argument indicates to which coordinate set
C                 this derivative is added.
C     MaxSet      The maximum number of definition sets allowed.
C     INwPrm      This output integer is filled with the primtive
C                 definition number of the final coordinate.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Dimension ITemp(MxI),RTemp(MxR),InfCrd(*),IVars(*),RVars(*),
     $  NCrdSt(*)
C
C     Initialize some flags.
C
      Call IClear(MxI,ITemp)
      Call AClear(MxR,RTemp)
C
C     Get the zero- and first-order terms, as needed.
C
      Call CrDfA1(IOut,IPrint,IDeriv,IAtom1,IAtom2,IAtom3,Idx,ISet,
     $  MaxSet,IA123,IAx,IR12,IR32,NPrim,NCoord,LIVars,LRVars,MPrim,
     $  MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
C
C     Evaluate the second-order terms if requested.
C
      If(IDeriv.ge.2) then
        Call GauErr('CrDfA2: IDeriv=2 NYI.')
        Call CrDfA1(IOut,IPrint,IDeriv,IAtom1,IAtom2,IAtom3,Idy,ISet,
     $    MaxSet,IA123,IAy,IR12,IR32,NPrim,NCoord,LIVars,LRVars,MPrim,
     $    MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
        endIf
      Return
      End
*Deck CrDfAn
      Subroutine CrDfAn(IOut,IPrint,IDeriv,IAtom1,IAtom2,IAtom3,Idx,Idy,
     $  IAct,ISet,MaxSet,INwPrm,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,
     $  MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine adds one standard internal coordinate and associated
C     definitions to the coordinate definition data structure, and is part
C     of the CrDfxx group of routines.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     IDeriv      An integer input argument indicating wheter the
C                 appropriate bond angle derivative should be added instead
C                 of the actual definition itself.
C     IAtom1-3    These are the numbers of the three atomic centers that
C                 define the bond angle being added. IAtom2 gives the
C                 central atom.
C     Idx         This input argument gives the Cartesian coordinate with
C                 respect to which A(IAtom1,IAtom2,IAtom3) is to be
C                 differentiated when IDeriv=1. If IDeriv=0, this argument
C                 is ignored.
C     ISet        This input argument indicates to which coordinate set
C                 this derivative is added.
C     MaxSet      The maximum number of definition sets allowed.
C     INwPrm      This output integer is filled with the primtive
C                 definition number of the final coordinate.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Logical AddEq,DoFix,Exists,OK,ValIs0
      Dimension ITemp(MxI),RTemp(MxR),InfCrd(*),IVars(*),RVars(*),
     $  NCrdSt(*)
 1000 Format(1x,'CrDfAn: Adding bond angle between atoms ',I6,',',I6,
     $  ' and ',I6,'.')
 1010 Format(1x,'CrDfAn is appending dA(',I6,',',I6,',',I6,')/dx')
 1020 Format(1x,'CrDfAn is appending d2A(',I6,',',I6,',',I6,')/dxdy')
 9000 Format(1x,'CrDfAn: IDeriv=',I2)
 9100 Format(1x,'CrDfAn: IAtom1=',I6,' IAtom2=',I6,' IAtom3=',I6)
C
C     Initialize some key flags. For now we statically set AddEq and DoFix.
C
      Call IClear(MxI,ITemp)
      Call AClear(MxR,RTemp)
      AddEq = .False.
      DoFix = .True.
      If(IDeriv.lt.0.or.IDeriv.gt.2) then
        Write(IOut,9000) IDeriv
        Call GauErr('CrDfBd: Invalid IDeriv.')
        endIf
      If(IAtom1.le.0.or.IAtom2.le.0.or.IAtom3.le.0) then
        Write(IOut,9100) IAtom1,IAtom2,IAtom3
        Call GauErr('CrDfAn: Invalid IAtom sent <= 0.')
        endIf
      If(IPrint.ge.5) then
        If(IDeriv.eq.0) then
          Write(IOut,1000) IAtom1,IAtom2,IAtom3
        else if(IDeriv.eq.1) then
          Write(IOut,1010) IAtom1,IAtom2,IAtom3
        else if(IDeriv.eq.2) then
          Write(IOut,1020) IAtom1,IAtom2,IAtom3
          endIf
        endIf
C
C     If we are differentiating an angle coordinate, do a quick check to
C     see if the derivative is zero.
C
      If(IDeriv.ge.1) then
        IdXAtm = -Idx/10
      else
        IdXAtm = 0
        endIf
      If(IDeriv.ge.2) then
        IdYAtm = -IdY/10
      else
        IdYAtm = 0
        endIf
      ValIs0 = (IDeriv.ge.1.and.IdXAtm.ne.IAtom1.and.IdXAtm.ne.IAtom2
     $  .and.IdXAtm.ne.IAtom3).or.
     $  (IDeriv.ge.2.and.IdYAtm.ne.IAtom1.and.IdYAtm.ne.IAtom2.and.
     $  IdYAtm.ne.IAtom3)
C
C     If ValIs0, then simply set the requested primitive to zero and
C     return. Otherwise, call CrDfA2 to set key intermediates and then
C     finalize the request angle value or derivative term.
C
      If(ValIs0) then
        Call CrDfCn(IOut,IPrint,0,ISet,MaxSet,INwPrm,NPrim,NCoord,
     $    LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $    NCrdSt)
      else
        Call CrDfA2(IOut,IPrint,IDeriv,IAtom1,IAtom2,IAtom3,Idx,Idy,
     $    ISet,MaxSet,IA123,IAx,IR12,IR32,NPrim,NCoord,LIVars,LRVars,
     $    MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
        If(IDeriv.ge.0) then
          IType = ICTACos
          JAct = -1
          NIVars = 1
          NRVars = 0
          ITemp(1) = IA123
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmA1,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfAn: Failure adding definition.')
          endIf
        If(IDeriv.eq.0) then
          IType = ICTAng
          JAct = IAct
          NIVars = 4
          NRVars = 0
          ITemp(1) = IPrmA1
          ITemp(2) = IAtom1
          ITemp(3) = IAtom2
          ITemp(4) = IAtom3
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmA,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfAn: Failure adding definition.')
          INwPrm = IPrmA
        else if(IDeriv.eq.1) then
          IType = ICTSq
          JAct = -1
          NIVars = 1
          NRVars = 0
          ITemp(1) = IA123
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfAn: Failure adding definition.')
          Call CrDfCn(IOut,IPrint,1,ISet,MaxSet,IPTmp2,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
          IType = ICTDif
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp2
          ITemp(2) = IPTmp1
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp3,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfAn: Failure adding definition.')
          IType = ICTSqrt
          JAct = -1
          NIVars = 1
          NRVars = 0
          ITemp(1) = IPTmp3
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp4,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfAn: Failure adding definition.')
          IPTmp5 = IPTmp4
          Call CrDfM1(IOut,IPrint,IPTmp5,ISet,MaxSet,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
          IType = ICTQuo
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IAx
          ITemp(2) = IPTmp5
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IAngX,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfAn: Failure adding definition.')
          INwPrm = IAngX
        else
          Write(IOut,9000) IDeriv
          Call GauErr('CrDfAn: Invalid IDeriv!')
          endIf
        endIf
      Return
      End
*Deck CrDfB1
      Subroutine CrDfB1(IOut,IPrint,IDeriv,DoBigR,IAtom1,IAtom2,Idx,
     $  ISet,MaxSet,Ir,IBigR,IRInv,IdrdX,MdrdX,IBgRx,NPrim,NCoord,
     $  LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $  NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to add the basic coordinate definitions needed
C     for bond stretch-related coordinates and coordinate derivatives.
C     Presently, this routine handles up to second-derivatives. The
C     primitive coordinate definitions formed here and returned to the
C     calling program unit, based on the input value of IDeriv, are:
C
C           For IDeriv.ge.0:
C                 Ir      = V12.V12
C             If DoBigR is sent as .True.:
C                 IBigR   = Sqrt[Ir]
C
C           For IDeriv.ge.1:
C                 IdrDx   = 0.5 * d/dx[V12.V12]
C                         = [KD(IAtom1,IdXAtm)-KD(IAtom2,IdXAtm)]*X12
C                 -- The symbol KD above is used to denote the
C                    Kroneker-delta function.
C                 -- Argument Idx gives the Cartesian coordinate x.
C                 -- IdrDx is returned as 0 if the value of this definition
C                    is 0.
C                 -- MdrDx is returned TRUE if definition IdrDx is actually
C                    (-1)*0.5*d/dx[V12.V12].
C
C           For IDeriv.ge.2:
C                 IdrDy   = 0.5 * d/dy[V12.V12]
C                         = [KD(IAtom1,IdYAtm)-KD(IAtom2,IdYAtm)]*Y12
C                 -- Argument Idy gives the Cartesian coordinate y.
C                 -- IdrDy is returned as 0 if the value of this definition
C                    is 0.
C                 -- MdrDy is returned TRUE if definition IdrDy is actually
C                    (-1)*0.5*d/dy[V12.V12].
C                 IdXdY   = 0.5 * d2/dxdy[V12.V12]
C                 -- Note that the value of this derivative is always -1,
C                    0, or -1. Therefore, IdXdY is NOT returned with a
C                    primitive definition number. Instead, it is returned
C                    with -1, 0, or 1.
C
C
C     NOTE: IdrdX returns the primitive number with 0.5 * the absolute
C     value of the derivative of V12.V12. MdrdX is a logical argument that
C     is returned with .TRUE. or .FALSE. to indicate whether (T) or not (F)
C     a factor of (-1) is needed to complete the derivative.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     IDeriv      An integer input argument indicating wheter the
C                 appropriate bond derivative of the bond definition should
C                 be added instead of the actual definition itself.
C     DoBigR      A logical input argument indicating whether (TRUE) or not
C                 (FALSE) the "big-R" terms should be formed.
C     IAtom1      The first atom involved in the bond stretch.
C     IAtom2      The second atom involved in the bond stretch.
C     Idx         This input argument gives the Cartesian coordinate with
C                 respect to which R(IAtom1,IAtom2) is to be differentiated
C                 when IDeriv >= 1. If IDeriv = 0, this argument is
C                 ignored.
C     Idy         This input argument gives the second Cartesian coordinate
C                 with respect to which R(IAtom1,IAtom2) is to be
C                 differentiated when IDeriv >= 2. If IDeriv < 2, this
C                 argument is ignored.
C     ISet        This input argument indicates to which coordinate set
C                 this derivative is added.
C     MaxSet      The maximum number of definition sets allowed.
C     INwPrm      This output integer is filled with the primtive
C                 definition number of the final coordinate.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Logical MdrdX,AddEq,DoFix,Exists,OK,DoBigR
      Dimension ITemp(MxI),RTemp(MxR),InfCrd(*),IVars(*),RVars(*),
     $  NCrdSt(*)
      Save IAtSv1,IAtSv2,IrSv,IBgRSv,IRInSv
      Data IAtSv1/0/,IAtSv2/0/,IrSv/0/,IBgRSv/0/,IRInSv/0/
 9000 Format(1x,'CrDfB1: IDeriv=',I2)
C
C     Print some details and do some basic error checking.
C
      Call IClear(MxI,ITemp)
      Call AClear(MxR,RTemp)
      If(IDeriv.lt.0.or.IDeriv.gt.2) then
        Write(IOut,9000) IDeriv
        Call GauErr('CrDfB1: Invalid IDeriv.')
        endIf
C
C     Initialize key settings and arguments.
C
      AddEq = .False.
      DoFix = .True.
      Ir    = 0
      IBigR = 0
      IRInv = 0
      IdrdX = 0
      IBgRx = 0
C
C     Based on IDeriv, either add the requested bond stretch of inverse
C     bond stretch definitions or add the requested bond stretch or inverse
C     bond stretch derivative definitions.
C
C     Add the bond stretch or inverse bond stretch definition. Three
C     primitives are necessary for a bond strech; inverse bond stretch
C     definitions require four primitives:
C           Ir = V(1,2).V(1,2)
C
      If(IDeriv.ge.0) then
        If(IAtom1.eq.IAtSv1.and.IAtom2.eq.IAtSv2) then
          Ir = IrSv
        else
          Call CrDfVV(IOut,IPrint,0,.False.,.True.,.True.,.True.,
     $      .False.,IAtom1,IAtom2,IAtom1,IAtom2,0,0,ISet,MaxSet,Ir,
     $      Junk1,Junk2,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,InfCrd,IVars,RVars,NCrdSt)
          IAtSv1 = IAtom1
          IAtSv2 = IAtom2
          IrSv = Ir
          IBgRSv = 0
          IRInSv = 0
          endIf
        If(DoBigR) then
          If(IAtom1.eq.IAtSv1.and.IAtom2.eq.IAtSv2.and.IBgRSv.ne.0) then
            IBigR = IBgRSv
          else
            IType = 62
            IAct = -1
            NIVars = 1
            NRVars = 0
            ITemp(1) = Ir
            Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IBigR,OK,
     $        MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $        MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $        IVars,RVars,NCrdSt)
            If(.not.OK)
     $        Call GauErr('CrDfB1: Failure adding definition.')
            IBgRSv = IBigR
            endIf
          endIf
        endIf
C
C     Add definitions that would be necessary for solving for the
C     first-derivative of a stretch-related definition. We return the
C     primitive numbers for these terms as 0 if they are 0 by definition.
C           ABS[0.5 * d(Ir)/dx] = IdrdX
C
      If(IDeriv.ge.1) then
        Call CrDfVV(IOut,IPrint,1,.True.,.True.,.False.,.True.,.False.,
     $    IAtom1,IAtom2,IAtom1,IAtom2,Idx,0,ISet,MaxSet,Junk1,IdrdX,
     $    Junk2,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,
     $    InfCrd,IVars,RVars,NCrdSt)
        MdrdX = IdrdX.lt.0
        If(MdrdX) IdrdX = -IdrdX
        If(IdrdX.ne.0) then
          If(DoBigR) then
            If(IAtom1.eq.IAtSv1.and.IAtom2.eq.IAtSv2.and.
     $        IRInSv.ne.0) then
              IRInv = IRInSv
            else
              IType = 61
              IAct = -1
              NIVars = 1
              NRVars = 0
              ITemp(1) = IBigR
              Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IRInv,OK,
     $          MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $          MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $          IVars,RVars,NCrdSt)
              If(.not.OK)
     $          Call GauErr('CrDfB1: Failure adding definition.')
              IRInSv = IRInv
              endIf
            IType = 52
            IAct = -1
            NIVars = 2
            NRVars = 0
            ITemp(1) = IdrdX
            ITemp(2) = IRInv
            Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IBgRx,OK,
     $        MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $        MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $        IVars,RVars,NCrdSt)
            If(.not.OK) Call
     $        GauErr('CrDfB1: Failure adding definition.')
            endIf
          endIf
        endIf
      Return
      End
*Deck CrDfB2
      Subroutine CrDfB2(IOut,IPrint,IDeriv,DoBigR,FillRx,IAtom1,IAtom2,
     $  Idx,Idy,ISet,MaxSet,Ir,IBigR,IdrdX,MdrdX,IBRx,IdrdY,MdrdY,IBRy,
     $  IdXdY,IBRxy,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $  MRVars,InfCrd,IVars,RVars,NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to add the basic coordinate definitions needed
C     for bond stretch-related coordinates and coordinate derivatives.
C     Presently, this routine handles up to second-derivatives. The
C     primitive coordinate definitions formed here and returned to the
C     calling program unit, based on the input value of IDeriv, are:
C
C           For IDeriv.ge.0:
C                 Ir      = V12.V12
C             If DoBigR is sent as .True.:
C                 IBigR   = Sqrt[Ir]
C
C           For IDeriv.ge.1:
C                 IdrDx   = 0.5 * d/dx[V12.V12]
C                         = [KD(IAtom1,IdXAtm)-KD(IAtom2,IdXAtm)]*X12
C                 -- The symbol KD above is used to denote the
C                    Kroneker-delta function.
C                 -- Argument Idx gives the Cartesian coordinate x.
C                 -- IdrDx is returned as 0 if the value of this definition
C                    is 0.
C                 -- MdrDx is returned TRUE if definition IdrDx is actually
C                    (-1)*0.5*d/dx[V12.V12].
C
C           For IDeriv.ge.2:
C                 IdrDy   = 0.5 * d/dy[V12.V12]
C                         = [KD(IAtom1,IdYAtm)-KD(IAtom2,IdYAtm)]*Y12
C                 -- Argument Idy gives the Cartesian coordinate y.
C                 -- IdrDy is returned as 0 if the value of this definition
C                    is 0.
C                 -- MdrDy is returned TRUE if definition IdrDy is actually
C                    (-1)*0.5*d/dy[V12.V12].
C                 IdXdY   = 0.5 * d2/dxdy[V12.V12]
C                 -- Note that the value of this derivative is always -1,
C                    0, or -1. Therefore, IdXdY is NOT returned with a
C                    primitive definition number. Instead, it is returned
C                    with -1, 0, or 1.
C
C
C     NOTE: IdrdX returns the primitive number with 0.5 * the absolute
C     value of the derivative of V12.V12. MdrdX is a logical argument that
C     is returned with .TRUE. or .FALSE. to indicate whether (T) or not (F)
C     a factor of (-1) is needed to complete the derivative.
C
C     This routine add a bond stretch or inverse bond stretch and necessary
C     associated definitions to the coordinate definition data structure.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     IDeriv      An integer input argument indicating wheter the
C                 appropriate bond derivative of the bond definition should
C                 be added instead of the actual definition itself.
C     DoBigR      A logical input argument indicating whether (TRUE) or not
C                 (FALSE) the "big-R" terms should be formed.
C     IAtom1      The first atom involved in the bond stretch.
C     IAtom2      The second atom involved in the bond stretch.
C     Idx         This input argument gives the Cartesian coordinate with
C                 respect to which R(IAtom1,IAtom2) is to be differentiated
C                 when IDeriv >= 1. If IDeriv = 0, this argument is
C                 ignored.
C     Idy         This input argument gives the second Cartesian coordinate
C                 with respect to which R(IAtom1,IAtom2) is to be
C                 differentiated when IDeriv >= 2. If IDeriv < 2, this
C                 argument is ignored.
C     ISet        This input argument indicates to which coordinate set
C                 this derivative is added.
C     MaxSet      The maximum number of definition sets allowed.
C     INwPrm      This output integer is filled with the primtive
C                 definition number of the final coordinate.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Logical MdrdX,MdrdY,MRxRy,AddEq,DoFix,Exists,OK,DoBigR,FillRx
      Dimension ITemp(MxI),RTemp(MxR),InfCrd(*),IVars(*),RVars(*),
     $  NCrdSt(*)
 9000 Format(1x,'CrDfB2: IDeriv=',I2)
C
C     Print some details and do some basic error checking.
C
      Call IClear(MxI,ITemp)
      Call AClear(MxR,RTemp)
      If(IDeriv.lt.0.or.IDeriv.gt.2) then
        Write(IOut,9000) IDeriv
        Call GauErr('CrDfB2: Invalid IDeriv.')
        endIf
C
C     Initialize key settings and arguments.
C
      AddEq = .False.
      DoFix = .True.
      Ir = 0
      IdrdX = 0
      MdrdX = .False.
      IdrdY = 0
      MdrdY = .False.
      IdXdY = 0
C
C     Begin by calling Routine CrDfB1 to define zero-th and, if requested
C     with IDeriv.ge.1, first-order terms.
C
      Call CrDfB1(IOut,IPrint,IDeriv,DoBigR,IAtom1,IAtom2,Idx,ISet,
     $  MaxSet,Ir,IBigR,IRInv,IdrdX,MdrdX,IBRx,NPrim,NCoord,LIVars,
     $  LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
C
C     Add the second-derivative of either the requested bond stretch or
C     inverse bond stretch.
C
      If(IDeriv.ge.2) then
        Call CrDfB1(IOut,IPrint,IDeriv,DoBigR,IAtom1,IAtom2,Idy,ISet,
     $    MaxSet,Ir,IBigR,IRInv,IdrdY,MdrdY,IBRy,NPrim,NCoord,LIVars,
     $    LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
        Call CrDfVV(IOut,IPrint,2,.True.,.True.,.True.,.False.,.False.,
     $    IAtom1,IAtom2,IAtom1,IAtom2,Idx,Idy,ISet,MaxSet,Junk1,Junk2,
     $    IdXdY,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,
     $    InfCrd,IVars,RVars,NCrdSt)
        IBRxy = 0
        If(DoBigR) then
          If(IBRx.eq.0.or.IBRy.eq.0) then
            IRxRy = 0
          else
            IType = 52
            IAct = -1
            NIVars = 2
            NRVars = 0
            ITemp(1) = IBRx
            ITemp(2) = IBRy
            Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IRxRy,OK,
     $        MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $        MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $        IVars,RVars,NCrdSt)
            If(.not.OK)Call GauErr('CrDfB2: Failure adding definition.')
            MRxRy = MdrdX.neqv.MdrdY
            endIf
          If(IdXdY.eq.0.and.IRxRy.ne.0) then
            IBRxy = IRxRy
            If(.not.MRxRy) Call CrDfM1(IOut,IPrint,IBRxy,ISet,MaxSet,
     $        NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,
     $        InfCrd,IVars,RVars,NCrdSt)
          else if(IdXdY.ne.0.and.IRxRy.ne.0) then
            Call CrDfCn(IOut,IPrint,IdXdY,ISet,MaxSet,IPTmp1,NPrim,
     $        NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $        IVars,RVars,NCrdSt)
            If(MRxRy) then
              IType = 50
            else
              IType = 51
              endIf
            IAct = -1
            NIVars = 2
            NRVars = 0
            ITemp(1) = IPTmp1
            ITemp(2) = IRxRy
            Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IBRxy,OK,
     $        MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $        MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $        IVars,RVars,NCrdSt)
            If(.not.OK)Call GauErr('CrDfB2: Failure adding definition.')
            endIf
          If(IBRxy.ne.0) then
            IType = 52
            IAct = -1
            NIVars = 2
            NRVars = 0
            ITemp(1) = IRInv
            ITemp(2) = IBRxy
            Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IBRxy,OK,
     $        MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $        MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $        IVars,RVars,NCrdSt)
            If(.not.OK)Call GauErr('CrDfB2: Failure adding definition.')
          else
            Call CrDfCn(IOut,IPrint,0,ISet,MaxSet,IBRxy,NPrim,NCoord,
     $        LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,
     $        RVars,NCrdSt)
            endIf
          endIf
        endIf
C
C     If FillRx is TRUE, then set first-derivative terms explicitly to 0.0
C     as appropriate.
C
      If(FillRx.and.IDeriv.ge.1.and.IBRx.eq.0) Call CrDfCn(IOut,IPrint,
     $  0,ISet,MaxSet,IBRx,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,
     $  MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
      If(FillRx.and.IDeriv.ge.2.and.IBRy.eq.0) Call CrDfCn(IOut,IPrint,
     $  0,ISet,MaxSet,IBRy,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,
     $  MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
      Return
      End
*Deck CrDfBd
      Subroutine CrDfBd(IOut,IPrint,IDeriv,DoInvR,IAtom1,IAtom2,Idx,Idy,
     $  IAct,ISet,MaxSet,INwPrm,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,
     $  MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine add a bond stretch or inverse bond stretch and necessary
C     associated definitions to the coordinate definition data structure.
C     The definition added is the (inverse) bond stretch between atoms
C     IAtom1 and IAtom2.  DoInvR is sent as .TRUE. if the inverse bond
C     stretch is desired.  This routine belongs to the CrDfxx group of
C     routines.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     IDeriv      An integer input argument indicating wheter the
C                 appropriate bond derivative of the bond definition should
C                 be added instead of the actual definition itself.
C     DoInvR      A logical input argument indicated whether (TRUE) or not
C                 (FALSE) the inverse bond stretch derinition is desired.
C     IAtom1      The first atom involved in the bond stretch.
C     IAtom2      The second atom involved in the bond stretch.
C     Idx         This input argument gives the Cartesian coordinate with
C                 respect to which R(IAtom1,IAtom2) is to be differentiated
C                 when IDeriv >= 1. If IDeriv = 0, this argument is
C                 ignored.
C     Idy         This input argument gives the second Cartesian coordinate
C                 with respect to which R(IAtom1,IAtom2) is to be
C                 differentiated when IDeriv >= 2. If IDeriv < 2, this
C                 argument is ignored.
C     ISet        This input argument indicates to which coordinate set
C                 this derivative is added.
C     MaxSet      The maximum number of definition sets allowed.
C     INwPrm      This output integer is filled with the primtive
C                 definition number of the final coordinate.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Logical AddEq,DoFix,Exists,MdrdX,MdrdY,OK,DoInvR
      Dimension ITemp(MxI),RTemp(MxR),InfCrd(*),IVars(*),RVars(*),
     $  NCrdSt(*)
      Save Zero
      Data Zero/0.d0/
 1000 Format(1x,'CrDfBd is appending R(',I5,',',I5,')')
 1002 Format(1x,'CrDfBd is appending 1/R(',I5,',',I5,')')
 1010 Format(1x,'CrDfBd is appending dR(',I5,',',I5,')/dx')
 1012 Format(1x,'CrDfBd is appending d[1/R(',I5,',',I5,')]/dx')
 9000 Format(1x,'CrDfBd: IDeriv=',I2)
C
C     Print some details and set AddEq and DoFix.
C
      Call IClear(MxI,ITemp)
      Call AClear(MxR,RTemp)
      AddEq = .False.
      DoFix = .True.
      If(IDeriv.lt.0.or.IDeriv.gt.2) then
        Write(IOut,9000) IDeriv
        Call GauErr('CrDfBd: Invalid IDeriv.')
        endIf
      If(IPrint.ge.5) then
        If(IDeriv.eq.0) then
          If(DoInvR) then
            Write(IOut,1002) IAtom1,IAtom2
          else
            Write(IOut,1000) IAtom1,IAtom2
            endIf
        else if(IDeriv.eq.1) then
          If(DoInvR) then
            Write(IOut,1012) IAtom1,IAtom2
          else
            Write(IOut,1010) IAtom1,IAtom2
            endIf
          endIf
        endIf
C
C     Begin by calling CrDfB2 to set-up some of the key definitions. Then
C     use code blocks over IDeriv to add any other definitions common to
C     derivatives between 0 and IDeriv.
C
      Call CrDfB2(IOut,IPrint,IDeriv,.True.,.True.,IAtom1,IAtom2,Idx,
     $  Idy,ISet,MaxSet,IV1212,IR12,IdrdX,MdrdX,IRx,IdrdY,MdrdY,IRy,
     $  IdXdY,IR12xy,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $  MRVars,InfCrd,IVars,RVars,NCrdSt)
C
C     Finalize the definition list we're building.
C
      If(IDeriv.eq.0) then
        IType = 1
        If(DoInvR) then
          JAct = -1
        else
          JAct = IAct
          endIf
        NIVars = 3
        NRVars = 0
        ITemp(1) = IR12
        ITemp(2) = IAtom1
        ITemp(3) = IAtom2
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrmR,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfBd: Failure adding definition.')
        INwPrm = IPrmR
        If(DoInvR) then
          IType = ICTInv
          JAct = IAct
          NIVars = 1
          NRVars = 0
          ITemp(1) = IPrmR
          Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrmRI,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK) Call GauErr('CrDfBd: Failure adding definition.')
          INwPrm = IPrmRI
          endIf
      else if(IDeriv.eq.1) then
        If(MdrdX.and.IRx.ne.0) Call CrDfM1(IOut,IPrint,IRx,ISet,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $    IVars,RVars,NCrdSt)
        INwPrm = IRx
      else if(IDeriv.eq.2) then
        If(IR12xy.eq.0) then
          IType = 1
          JAct = -1
          NIVars = 0
          NRVars = 1
          RTemp(1) = Zero
          Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IR12xy,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK) Call GauErr('CrDfBd: Failure adding definition.')
        else
          INwPrm = IR12xy
          endIf
      else
        Write(IOut,9000) IDeriv
        Call GauErr('CrDfBd: Invalid IDeriv.')
        endIf
      Return
      End
*Deck CrDfC0
      Subroutine CrDfC0(IOut,IPrint,SetPrm,IConst,INwPrm)
      Implicit Real*8(A-H,O-Z)
C
C     This routine keeps track of the primitive numbers defining standard
C     constants. If SetPrm is .False., then input argument IConst indicates
C     which constant the calling unit seeks, and INwPrim is returned with
C     the primitive corresponding to that constant. If a constant primitive
C     with that value has yet to be defined, INwPrim is returned as 0.
C
C     If SetPrm is .True., then INwPrm is an input argument holding the
C     primitive number where IConst has been added.
C
C     The table given below gives the list of constants known to this
C     routine.
C     ---------------------------------------------
C         IConst       Value
C     ---------------------------------------------
C           0           0.0
C          11           0.5
C           1           1.0
C           2           2.0
C         100           pi
C         101           2*pi
C     ---------------------------------------------
C
C     To request the negative of any Value in the table above, send
C     -IConst.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
      Logical SetPrm
      Save ISM101,ISM100,ISM2,ISM1,ISM11,ISv0,ISv11,ISv1,ISv2,ISv100,
     $  ISv101
      Data ISM101/0/,ISM100/0/,ISM2/0/,ISM1/0/,ISM11/0/,ISv0/0/,
     $  ISv11/0/,ISv1/0/,ISv2/0/,ISv100/0/,ISv101/0/
 1000 Format(1x,'CrDfC0: IConst = ',I3)
C
C     Do the work...
C
      If(IPrint.ge.3) Write(IOut,1000) IConst
      If(SetPrm) then
        If(IConst.eq.-101) then
          ISM101 = INwPrm
        else if(IConst.eq.-100) then
          ISM100 = INwPrm
        else if(IConst.eq.-2) then
          ISM2   = INwPrm
        else if(IConst.eq.-1) then
          ISM1   = INwPrm
        else if(IConst.eq.-11) then
          ISM11  = INwPrm
        else if(IConst.eq.0) then
          ISv0   = INwPrm
        else if(IConst.eq.11) then
          ISv11  = INwPrm
        else if(IConst.eq.1) then
          ISv1   = INwPrm
        else if(IConst.eq.2) then
          ISv2   = INwPrm
        else if(IConst.eq.100) then
          ISv100 = INwPrm
        else if(IConst.eq.101) then
          ISv101 = INwPrm
        else
          Write(IOut,1000) IConst
          Call GauErr('CrDfC0: Unknown IConst.')
          endIf
      else
        If(IConst.eq.-101) then
          INwPrm = ISM101
        else if(IConst.eq.-100) then
          INwPrm = ISM100
        else if(IConst.eq.-2) then
          INwPrm = ISM2
        else if(IConst.eq.-1) then
          INwPrm = ISM1
        else if(IConst.eq.-11) then
          INwPrm = ISM11
        else if(IConst.eq.0) then
          INwPrm = ISv0
        else if(IConst.eq.11) then
          INwPrm = ISv11
        else if(IConst.eq.1) then
          INwPrm = ISv1
        else if(IConst.eq.2) then
          INwPrm = ISv2
        else if(IConst.eq.100) then
          INwPrm = ISv100
        else if(IConst.eq.101) then
          INwPrm = ISv101
        else
          Write(IOut,1000) IConst
          Call GauErr('CrDfC0: Unknown IConst.')
          endIf
        endIf
      Return
      End
*Deck CrDfCar
      Subroutine CrDfCar(IOut,IPrint,IAct,ISet,NAtoms,Ixyz,MaxSet,NPrim,
     $  NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,
     $  RVars,NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine adds Cartesian coordinate definitions to the coordinate
C     definition data structure, and is part of the CrDfxx group of
C     routines.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     IAct        The value to set when filling action status flags for
C                 these Cartesian coordinates in InfCrd.
C     ISet        This input argument indicates to which coordinate set
C                 this derivative is added.
C     NAtoms      The number of atoms. If NAtoms < 0, then only the three
C                 Cartesian coordinates on atom -NAtom are added.
C     Ixyz        An input integer indicating whether all or only one
C                 Cartesian should be added for each atom specified by
C                 NAtoms. Ixyz=0 indicates that all three Cartesians should
C                 be added; Ixyz={1,2,3} indicated that only the {x,y,z}
C                 Cartesian should be added.
C     MaxSet      The maximum number of definition sets allowed.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Logical AddEq,DoFix,Exists,OK
      Dimension ITemp(MxI),RTemp(MxR),InfCrd(*),IVars(*),RVars(*),
     $  NCrdSt(*)
 1000 Format(1x,'CrDf02 is appending Cartesian coordinates to the ',
     $  'coordinate definition arrays.',/,1x,'NAtoms=',I5)
C
C     For now we statically set AddEq and DoFix.
C
      Call IClear(MxI,ITemp)
      Call AClear(MxR,RTemp)
      AddEq = .True.
      DoFix = .False.
C
C     Set up the range of atoms to loop over, then add all of the requested
C     Cartesian coordinates.
C
      If(IPrint.ge.4) Write(IOut,1000) NAtoms
      If(NAtoms.gt.0) then
        IAtSt = 1
        IAtEnd = NAtoms
      else
        IAtSt = -NAtoms
        IAtEnd = -NAtoms
        endIf
      If(Ixyz.lt.0.or.Ixyz.gt.3) Call GauErr('CrDfCar: Invalid Ixyz.')
      If(Ixyz.eq.0) then
        IxyzSt = 1
        IxyzEn = 3
      else
        IxyzSt = Ixyz
        IxyzEn = Ixyz
        endIf
      Do 100 IAtom = IAtSt,IAtEnd
        Do 100 i = IxyzSt,IxyzEn
          IType = 0
          NIVars = 2
          NRVars = 0
          ITemp(1) = IAtom
          ITemp(2) = -i
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,INwPrm,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK) Call GauErr('CrDf02: After CrdAdd, OK=FALSE.')
 100      Continue
      Return
      End
*Deck CrDfCn
      Subroutine CrDfCn(IOut,IPrint,IConst,ISet,MaxSet,INwPrm,NPrim,
     $  NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,
     $  RVars,NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to add a standard constant to the internal
C     coordinate definitions data structure. The argument IConst indicates
C     which constant to add. The table given below gives the list of
C     constants known to this routine.
C
C     ---------------------------------------------
C         IConst       Value
C     ---------------------------------------------
C           0           0.0
C          11           0.5
C           1           1.0
C           2           2.0
C         100           pi
C         101           2*pi
C     ---------------------------------------------
C
C     To request the negative of any Value in the table above, send
C     -IConst.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Logical AddEq,DoFix,Exists,OK
      Dimension ITemp(MxI),RTemp(MxR),InfCrd(*),IVars(*),RVars(*),
     $  NCrdSt(*)
      Save Zero,Half,One,Two
      Data Zero/0.d0/,Half/0.5d0/,One/1.d0/,Two/2.d0/
 1000 Format(1x,'CrDfCn: IConst = ',I3)
 2000 Format(1x,'CrDfCn is adding definition for IConst = ',I3)
C
      If(IPrint.ge.5) Write(IOut,1000) IConst
C
C     Initialize some key flags. For now we statically set AddEq and DoFix.
C
      Call IClear(MxI,ITemp)
      Call AClear(MxR,RTemp)
      AddEq = .False.
      DoFix = .True.
C
C     Add the requested coordinate.
C
      Call CrDfC0(IOut,IPrint,.False.,IConst,INwPrm)
      If(INwPrm.eq.0) then
        If(IPrint.ge.5) Write(IOut,2000) IConst
        IType = ICTCon
        IAct = -1
        NIVars = 0
        NRVars = 1
        If(IConst.eq.0) then
          RTemp(1) = Zero
        else if(IAbs(IConst).eq.11) then
          RTemp(1) = Half
        else if(IAbs(IConst).eq.1) then
          RTemp(1) = One
        else if(IAbs(IConst).eq.2) then
          RTemp(1) = Two
        else if(IAbs(IConst).eq.100) then
          RTemp(1) = GFloat(4)*GAtan(One)
        else if(IAbs(IConst).eq.101) then
          RTemp(1) = GFloat(8)*GAtan(One)
        else
          Write(IOut,1000) IConst
          Call GauErr('CrDfCn: Invalid IConst.')
          endIf
        If(IConst.lt.0) RTemp(1) = -RTemp(1)
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,INwPrm,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfCn: Failure adding definition.')
        Call CrDfC0(IOut,IPrint,.True.,IConst,INwPrm)
        endIf
C
      Return
      End
*Deck CrDfDi
      Subroutine CrDfDi(IOut,IPrint,IDeriv,IAtom1,IAtom2,IAtom3,IAtom4,
     $  Idx,IAct,ISet,MaxSet,INwPrm,NPrim,NCoord,LIVars,LRVars,MPrim,
     $  MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine adds one standard internal coordinate and associated
C     definitions to the coordinate definition data structure, and is part
C     of the CrDfxx group of routines.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     IDeriv      An integer input argument indicating whether the
C                 appropriate dihedral derivative should be added instead
C                 of the actual definition itself.
C     IAtom1-4    These are the numbers of the four atomic centers that
C                 define the dihedral angle being added. IAtom2 and IAtom3
C                 give the central atoms.
C     Idx         This input argument gives the Cartesian coordinate with
C                 respect to which D(IAtom1,IAtom2,IAtom3,IAtom4) is to be
C                 differentiated when IDeriv=1. If IDeriv=0, this argument
C                 is ignored.
C     ISet        This input argument indicates to which coordinate set
C                 this derivative is added.
C     MaxSet      The maximum number of definition sets allowed.
C     INwPrm      This output integer is filled with the primtive
C                 definition number of the final coordinate.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     The dihedral definition used here follows the formulation given by
C     Tamar Schlick in JCC, 10, 951-956 (1989).
C
C     The dihedral angle (1-2-3-4) is evaluated as
C           arccos(Q)
C     where
C           Q  = N/D
C
C           N  = U.V
C              = (V12.V32)(V32.V43) - (V12.V43)(V32.V32)
C
C           D  = SQRT[D0]
C
C           D0 = (U.U)(V.V)
C              = (U.U)(V.V)
C              = [(V12.V12)(V32.V32)-(V12.V32)**2] *
C                  [(V32.V32)(V43.V43)-(V32.V43)**2]
C
C     The first derivative of the dihedral angle (1-2-3-4) is evaluated as
C            d                      1       -dQ
C           ---- arccos(Q) = -------------- ----
C            dx               SQRT[1-Q**2]   dx
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Logical AddEq,DoFix,Exists,OK,M1223,M1234,M2334
      Dimension ITemp(MxI),RTemp(MxR),InfCrd(*),IVars(*),RVars(*),
     $  NCrdSt(*)
 1000 Format(1x,'CrDfDi: Adding dihedral angle between atoms ',I6,',',
     $  I6,',',I6,' and ',I6,'.')
 1010 Format(1x,'CrDfDi is appending dD(',I6,',',I6,',',I6,',',I6,
     $  ')/dx')
 9000 Format(1x,'CrDfAn: IDeriv=',I2)
 9100 Format(1x,'CrDfDi: IAtom1=',I6,' IAtom2=',I6,' IAtom3=',I6,
     $  ' IAtom4=',I6)
 9200 Format(1x,'CrDfDi: IdXAtm=',I6,' IdX1=',I6)
C
C     Initialize some key flags. For now we statically set AddEq and DoFix.
C
      Call IClear(MxI,ITemp)
      Call AClear(MxR,RTemp)
      AddEq = .False.
      DoFix = .False.
      If(IDeriv.lt.0.or.IDeriv.gt.1) then
        Write(IOut,9000) IDeriv
        Call GauErr('CrDfBd: Invalid IDeriv.')
        endIf
      If(IAtom1.le.0.or.IAtom2.le.0.or.IAtom3.le.0.or.IAtom4.le.0) then
        Write(IOut,9100) IAtom1,IAtom2,IAtom3,IAtom4
        Call GauErr('CrDfDi: Invalid IAtom sent <= 0.')
        endIf
      If(IPrint.ge.5) then
        If(IDeriv.eq.0) then
          Write(IOut,1000) IAtom1,IAtom2,IAtom3,IAtom4
        else if(IDeriv.eq.1) then
          Write(IOut,1010) IAtom1,IAtom2,IAtom3,IAtom4
          endIf
        endIf
      INwPrm = 0
C
C     Define terms that are common to IDeriv >= 0 definitions. These
C     include some key dot-products:
C           IV1212 = (V12.V12)
C           IV2323 = (V23.V23)
C           IV3434 = (V34.V34)
C     and
C           IV1223 = (V12.V23)
C           IV1234 = (V12.V34)
C           IV2334 = (V23.V34)
C     Note that when the second batch of dot-products are defined, we also
C     set logical flags -- M1223, M1234, and M2334 -- indicating whether
C     the cannonical form requires multiplication by -1.
C
C     We also define N (IPrmN), D0 (IPrmD0), and D(IPrmD) as defined above.
C     D is built up from lower terms set to primitive numbers IPrmUU and
C     IPrmVV:
C           N   = (V12.V23)(V23.V34) - (V12.V34)(V23.V23)
C               = IV1223*V2334-V1234*V2323
C               = IPrmN
C           U.U = (V12.V12)(V23.V23) - (V12.V23)**2
C               = IV1212*IV2323-V1223**2
C               = IPrmUU
C           V.V = (V23.V23)(V34.V34) - (V23.V34)**2
C               = IV2323*IV3434-IV2334**2
C               = IPrmVV
C     Note that IPrmN, IPrmUU, and IPrmVV are defined such that they
C     include any factors of -1 according to M1223, M1234, and M2334.
C
C     Define terms that are common to all definitions IDeriv >= 0.
C
      If(IDeriv.ge.0) then
C
C       Define IV1212, IV2323, and IV3434.
        Call CrDfVV(IOut,IPrint,0,.False.,.True.,.True.,.True.,.False.,
     $    IAtom1,IAtom2,IAtom1,IAtom2,0,0,ISet,MaxSet,IV1212,IVVx,IVVxy,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $    IVars,RVars,NCrdSt)
        Call CrDfVV(IOut,IPrint,0,.False.,.True.,.True.,.True.,.False.,
     $    IAtom2,IAtom3,IAtom2,IAtom3,0,0,ISet,MaxSet,IV2323,IVVx,IVVxy,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $    IVars,RVars,NCrdSt)
        Call CrDfVV(IOut,IPrint,0,.False.,.True.,.True.,.True.,.False.,
     $    IAtom3,IAtom4,IAtom3,IAtom4,0,0,ISet,MaxSet,IV3434,IVVx,IVVxy,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $    IVars,RVars,NCrdSt)
C
C       Define IV1223, IV1234, and IV2334.
        Call CrDfA2(IOut,IPrint,0,IAtom1,IAtom2,IAtom3,0,0,ISet,MaxSet,
     $    IA123,IAx,IR12,IR32,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,
     $    MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
        Call CrDfVV(IOut,IPrint,0,.False.,.True.,.True.,.True.,.False.,
     $    IAtom1,IAtom2,IAtom2,IAtom3,0,0,ISet,MaxSet,IV1223,IVVx,IVVxy,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $    IVars,RVars,NCrdSt)
        M1223 = IV1223.lt.0
        If(M1223) IV1223 = -IV1223
        Call CrDfVV(IOut,IPrint,0,.False.,.True.,.True.,.True.,.False.,
     $    IAtom1,IAtom2,IAtom3,IAtom4,0,0,ISet,MaxSet,IV1234,IVVx,IVVxy,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $    IVars,RVars,NCrdSt)
        M1234 = IV1234.lt.0
        If(M1234) IV1234 = -IV1234
        Call CrDfVV(IOut,IPrint,0,.False.,.True.,.True.,.True.,.False.,
     $    IAtom2,IAtom3,IAtom3,IAtom4,0,0,ISet,MaxSet,IV2334,IVVx,IVVxy,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $    IVars,RVars,NCrdSt)
        M2334 = IV1234.lt.0
        If(M2334) IV2334 = -IV2334
C
C       Define IA123 and IA234.
        IType = 62
        JAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IV3434
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IR34,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = 52
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IR32
        ITemp(2) = IR34
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = 53
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IV2334
        ITemp(2) = IPTmp1
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IA234,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        If(.not.M2334) Call CrDfM1(IOut,IPrint,IA234,ISet,MaxSet,NPrim,
     $    NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,
     $    RVars,NCrdSt)
C
C       Put the numerator together in IPrmN.
        IType = 52
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IV1223
        ITemp(2) = IV2334
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = 52
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IV2323
        ITemp(2) = IV1234
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp2,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        If((M1223.eqv.M2334).and.(.not.M1234)) then
          IType = 51
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp1
          ITemp(2) = IPTmp2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmN,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        else if(((M1223.eqv.M2334).and.(M1234)).or.
     $    ((M1223.neqv.M2334).and.(.not.M1234))) then
          IType = 50
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp1
          ITemp(2) = IPTmp2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmN,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          If((M1223.neqv.M2334).and.(.not.M1234)) Call CrDfM1(IOut,
     $      IPrint,IPrmN,ISet,MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,
     $      MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
        else if((M1223.neqv.M2334).and.(M1234)) then
          IType = 51
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp2
          ITemp(2) = IPTmp1
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmN,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        else
          Call GauErr('CrDfDi: Logic Error A.')
          endIf
C
C       Put the denominator together in IPrmD.
        IType = 52
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IV1212
        ITemp(2) = IV2323
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = 63
        JAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IV1223
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp2,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = 51
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPTmp1
        ITemp(2) = IPTmp2
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmD,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = 52
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IV2323
        ITemp(2) = IV3434
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = 63
        JAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IV2334
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp2,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = 51
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPTmp1
        ITemp(2) = IPTmp2
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = 52
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrmD
        ITemp(2) = IPTmp1
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmD,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = 62
        JAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IPrmD
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmD,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
C
C       Form IPrmQ = IPrmN/IPrmD.
        IType = 53
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrmN
        ITemp(2) = IPrmD
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmQ,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
C
C       Form IPSign.
        Call CrDfVX(IOut,IPrint,0,IAtom4,IAtom3,IAtom2,IAtom3,ISet,
     $    MaxSet,ICrsXa,ICrsYa,ICrsZa,NPrim,NCoord,LIVars,LRVars,MPrim,
     $    MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
        IType = 110
        JAct = -1
        NIVars = 3
        NRVars = 0
        ITemp(1) = -1
        ITemp(2) = IAtom2
        ITemp(3) = IAtom1
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IX21,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = ICTPrd
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IX21
        ITemp(2) = ICrsXa
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = 110
        JAct = -1
        NIVars = 3
        NRVars = 0
        ITemp(1) = -2
        ITemp(2) = IAtom2
        ITemp(3) = IAtom1
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IY21,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = ICTPrd
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IY21
        ITemp(2) = ICrsYa
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp2,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = ICTSum
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPTmp1
        ITemp(2) = IPTmp2
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp2,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = 110
        JAct = -1
        NIVars = 3
        NRVars = 0
        ITemp(1) = -3
        ITemp(2) = IAtom2
        ITemp(3) = IAtom1
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IZ21,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = ICTPrd
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IZ21
        ITemp(2) = ICrsZa
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = ICTSum
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPTmp1
        ITemp(2) = IPTmp2
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPSign,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        endIf
C
C     Define terms that are common to all definitions IDeriv >= 1.
C
      If(IDeriv.ge.1) then
        If(IdX.ge.0) Call GauErr('CrDfDi: Invalid IdX.')
        IdXAtm = -Idx/10
        IdX1 = Mod(-Idx,10)
        IdX1 = -IdX1
C
C       Fill in the x, y, and/or z definition(s) for A = V12 x V23 and B =
C       V43 x V23. The definitions for A are stored in primitive numbers
C       ICrsXa, ICrdYa, and ICrdZa; the definitions for B are stored in
C       primitive numbers ICrsXb, ICrdYb, and ICrdZb.
        Call CrDfVX(IOut,IPrint,0,IAtom1,IAtom2,IAtom2,IAtom3,ISet,
     $    MaxSet,ICrsXa,ICrsYa,ICrsZa,NPrim,NCoord,LIVars,LRVars,MPrim,
     $    MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
        Call CrDfVX(IOut,IPrint,0,IAtom4,IAtom3,IAtom2,IAtom3,ISet,
     $    MaxSet,ICrsXb,ICrsYb,ICrsZb,NPrim,NCoord,LIVars,LRVars,MPrim,
     $    MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
C
C       Fill the definition for A.A which is stored in IPrmA2.
        IType = 52
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IV1212
        ITemp(2) = IV2323
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmA2,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IPTmp1 = IV1223
        If(M1223) Call CrDfM1(IOut,IPrint,IPTmp1,ISet,MaxSet,NPrim,
     $    NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,
     $    RVars,NCrdSt)
        IType = 63
        JAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IPTmp1
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = 51
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrmA2
        ITemp(2) = IPTmp1
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmA2,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
C
C       Fill the definition for B.B which is stored in IPrmB2.
        IType = 52
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IV2323
        ITemp(2) = IV3434
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmB2,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IPTmp1 = IV2334
        If(.not.M2334) Call CrDfM1(IOut,IPrint,IPTmp1,ISet,MaxSet,NPrim,
     $    NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,
     $    RVars,NCrdSt)
        IType = 63
        JAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IPTmp1
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = 51
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPrmB2
        ITemp(2) = IPTmp1
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmB2,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
C
C       Complete the derivative wrt a Cartesian coordinate on Atom 1.
        If(IdXAtm.eq.IAtom1) then
          IType = 53
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IR32
          ITemp(2) = IPrmA2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmDx,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          If(IdX1.eq.-1) then
            IPTmp1 = ICrsXa
          else if(IdX1.eq.-2) then
            IPTmp1 = ICrsYa
          else if(IdX1.eq.-3) then
            IPTmp1 = ICrsZa
          else
            Write(IOut,9200) IdXAtm,IdX1
            Call GauErr('CrDfDi: Logic Error B')
            endIf
          IType = 52
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPrmDx
          ITemp(2) = IPTmp1
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmDx,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          Call CrDfM1(IOut,IPrint,IPrmDx,ISet,MaxSet,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
C
C       Derivative wrt a Cartesian coordinate on Atom 2.
        else if(IdXAtm.eq.IAtom2) then
          IType = 53
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IV1223
          ITemp(2) = IR32
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          If(M1223) then
            IType = 51
          else
            IType = 50
            endIf
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IR32
          ITemp(2) = IPTmp1
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          IType = 53
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp1
          ITemp(2) = IPrmA2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          If(IdX1.eq.-1) then
            IPTmp2 = ICrsXa
          else if(IdX1.eq.-2) then
            IPTmp2 = ICrsYa
          else if(IdX1.eq.-3) then
            IPTmp2 = ICrsZa
          else
            Call GauErr('CrDfDi: Logic Error C')
            endIf
          IType = 52
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp1
          ITemp(2) = IPTmp2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          IType = 52
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IR32
          ITemp(2) = IPrmB2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp2,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          IType = 53
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IV2334
          ITemp(2) = IPTmp2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp2,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          If(IdX1.eq.-1) then
            IPTmp3 = ICrsXb
          else if(IdX1.eq.-2) then
            IPTmp3 = ICrsYb
          else if(IdX1.eq.-3) then
            IPTmp3 = ICrsZb
          else
            Call GauErr('CrDfDi: Logic Error D')
            endIf
          IType = 52
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp2
          ITemp(2) = IPTmp3
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp2,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          If(M2334) then
            IType = 51
          else
            IType = 50
            endIf
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp1
          ITemp(2) = IPTmp2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmDx,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
C
C       Derivative wrt a Cartesian coordinate on Atom 3.
        else if(IdXAtm.eq.IAtom3) then
          IType = 53
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IV2334
          ITemp(2) = IR32
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          If(M2334) then
            IType = 51
          else
            IType = 50
            endIf
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IR32
          ITemp(2) = IPTmp1
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          IType = 53
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp1
          ITemp(2) = IPrmB2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          If(IdX1.eq.-1) then
            IPTmp2 = ICrsXb
          else if(IdX1.eq.-2) then
            IPTmp2 = ICrsYb
          else if(IdX1.eq.-3) then
            IPTmp2 = ICrsZb
          else
            Call GauErr('CrDfDi: Logic Error E')
            endIf
          IType = 52
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp1
          ITemp(2) = IPTmp2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          IType = 52
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IR32
          ITemp(2) = IPrmA2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp2,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          IType = 53
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IV1223
          ITemp(2) = IPTmp2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp2,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          If(IdX1.eq.-1) then
            IPTmp3 = ICrsXa
          else if(IdX1.eq.-2) then
            IPTmp3 = ICrsYa
          else if(IdX1.eq.-3) then
            IPTmp3 = ICrsZa
          else
            Call GauErr('CrDfDi: Logic Error F')
            endIf
          IType = 52
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp2
          ITemp(2) = IPTmp3
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp2,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          If(M1223) then
            IType = 51
          else
            IType = 50
            endIf
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp1
          ITemp(2) = IPTmp2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmDx,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          Call CrDfM1(IOut,IPrint,IPrmDx,ISet,MaxSet,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
C
C       Complete the derivative wrt a Cartesian coordinate on Atom 4.
        else if(IdXAtm.eq.IAtom4) then
          IType = 53
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IR32
          ITemp(2) = IPrmB2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmDx,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
          If(IdX1.eq.-1) then
            IPTmp1 = ICrsXb
          else if(IdX1.eq.-2) then
            IPTmp1 = ICrsYb
          else if(IdX1.eq.-3) then
            IPTmp1 = ICrsZb
          else
            Call GauErr('CrDfDi: Logic Error G')
            endIf
          IType = 52
          JAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPrmDx
          ITemp(2) = IPTmp1
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmDx,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
C
C       The requested derivative is invalid.
        else
          Write(IOut,9000) IDeriv
          Write(IOut,9100) IAtom1,IAtom2,IAtom3,IAtom4
          Call GauErr('CrDfDi: Unable to form first derivative of D.')
          endIf
        endIf
C
C     Using IDeriv, finish the requested definition.
C
      If(IDeriv.eq.0) then
        IType = 73
        JAct = -1
        NIVars = 1
        NRVars = 0
        ITemp(1) = IPrmQ
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = ICTSign
        JAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPTmp1
        ITemp(2) = IPSign
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        IType = 3
        JAct = IAct
        NIVars = 5
        NRVars = 0
        ITemp(1) = IPTmp1
        ITemp(2) = IAtom1
        ITemp(3) = IAtom2
        ITemp(4) = IAtom3
        ITemp(5) = IAtom4
        Call CrdAdd(IOut,IPrint,AddEq,1,0,DoFix,Exists,IPrmDi,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK)Call GauErr('CrDfDi: Failure adding definition.')
        INwPrm = IPrmDi
      else if(IDeriv.eq.1) then
C
C       Finalize the first-derivative.
        INwPrm = IPrmDx
      else
        Write(IOut,9000) IDeriv
        Call GauErr('CrDfDi: Invalid IDeriv.')
        endIf
      Return
      End
*Deck CrDfDx
      Subroutine CrDfDx(IOut,IPrint,ISet,ISet2,NAtoms,MaxSet,NPrim,
     $  NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,
     $  RVars,NCrdSt,ICrdSt,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to define new derivative coordinates of all
C     active/propagating coordinates in definition set ISet2 with respect
C     to all Cartesians. These new derivatives are added to definition set
C     ISet. For the definition types where enough information is known, we
C     avoid explicitly defining derivatives that we know will be zero.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     NAtoms      The number of atoms.
C     MaxSet      The maximum number of definition sets allowed.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C     IType       The value to set for the coordinate type flag.
C     IAct        The value to set for the action status flag.
C     ISet        The value to set for the set number flag.
C     NIVars      The number of integer values to add to IVars.  Note that
C                 this argument is modified if the new coordinate already
C                 exists in the data structure.
C     NRVars      The number of real values to add to RVars.  Note that
C                 this argument is modified if the new coordinate already
C                 exists in the data structure.
C     IAddVr      The list of integer values to add to IVars. Its length is
C                 NIVars.  Note that this argument is modified if the new
C                 coordinate already exists in the data structure.
C     RAddVr      The list of real values to add to IVars. Its length is
C                 NRVars.
C
C     InfCrd, IVars, RVars, NCrdSt, and ICrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Logical LogTmp
      Parameter (MaxAtm=4,MaxDx=12)
      Parameter (NInfCr=DEFMAXCOORDINFO)
      Dimension IAtom(MaxAtm),IdqLst(MaxDx),InfCrd(NInfCr,*),IVars(*),
     $  RVars(*),NCrdSt(*),ICrdSt(*),V(MDV)
      Save JdAtom,IdqLst
      Data JdAtom/0/,IAtom/MaxAtm*0/,IdqLst/MaxDx*0/
 1000 Format(1x,'Enter CrDfDx.',/,4x,'ISet=',I5,' ISet2=',I5)
C
C     Figure out which coordinates need to be differentiated and then do
C     the work.
C
      If(IPrint.ge.1) Write(IOut,1000) ISet,ISet2
      IPrim1 = IArSum(1,ISet2,NCrdSt)
      IPrim0 = IPrim1-NCrdSt(ISet2)+1
      ITch1  = 1
      ITch2  = ITch1  + InToWP(NPrim)
      IDrive = ITch2  + InToWP(NAtoms)
      IV     = IDrive + InToWP(NAtoms)
      IIScr  = IV
      MDV1 = MDV-IV+1
      Call TstCor(MDV1,MDV,'CrDfDx-1')
      IILen = MDV1*IntPWP(0)
      Do 100 IPrimX = IPrim0,IPrim1
C
C       Get JType, JJPos, and other flags properly set.
        IPrim = ICrdSt(IPrimX)
        IType = InfCrd(InfTyp,IPrim)
        IIPos = InfCrd(InfIVr,IPrim)
        JType = IType
        JJPos = IIPos
        JdAtom = 0
        JdXYZ = 0
        If(IType.eq.500) then
          If(IVars(IIPos+2).ge.0) then
            Call GauErr('CrDfDx: Consistency Failure A.')
          else
            JdAtom = -IVars(IIPos+2)/10
            JdXYZ = -Mod(IVars(IIPos+2),10)
            endIf
          KType = InfCrd(InfTyp,IVars(IIPos+1))
          If(KType.eq.0.or.KType.eq.1.or.KType.eq.2.or.KType.eq.3) then
            JType = KType
            JJPos = InfCrd(InfIVr,IVars(IIPos+1))
            endIf
          endIf
C
C       Set IAtom1-4 according to JType. Then, use this together with other
C       information (depending on the order of the derivative being formed)
C       to form the list of derivatives we will request (IdqLst).
C
        Call IClear(4,IAtom)
        If(JType.eq.0) then
          IAtom(1) = IVars(JJPos)
        else if(JType.eq.1) then
          IAtom(1) = IVars(JJPos+1)
          IAtom(2) = IVars(JJPos+2)
        else if(JType.eq.2) then
          IAtom(1) = IVars(JJPos+1)
          IAtom(2) = IVars(JJPos+2)
          IAtom(3) = IVars(JJPos+3)
        else if(JType.eq.3) then
          IAtom(1) = IVars(JJPos+1)
          IAtom(2) = IVars(JJPos+2)
          IAtom(3) = IVars(JJPos+3)
          IAtom(4) = IVars(JJPos+4)
          endIf
        NumIq2 = 0
        If(JType.ge.0.and.JType.le.4) then
          Do 200 i=1,4
            If(IAtom(i).ne.0.and.
     $        (JdAtom.eq.0.or.IAtom(i).le.JdAtom)) then
              MaxIdx = 3
              If(IAtom(i).eq.JdAtom) MaxIdx = JdXYZ
              Do 210 Idx=1,MaxIdx
                NumIq2 = NumIq2+1
                IdqLst(NumIq2) = -(10*IAtom(i)+Idx)
 210            Continue
              endIf
 200        Continue
          endIf
C
C       If NumIq2 is non-zero, then we have already made the list of
C       derivatives needed and we can call CrdDq0 directly. Otherwise,
C       we're differentiating a coordinate with an unknown type and we call
C       CrTree to generate the list of atomic centers this coordinate
C       depends upon. Using that list, we request appropriate derivatives.
C
        If(NumIq2.gt.0) then
          Call CrdDq0(IOut,IPrint,ISet,IPrim,IPrimX,IdqLst,NumIq2,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,InfCrd,IVars,RVars,NCrdSt,V(IIScr),IILen)
        else
          IPrim = ICrdSt(IPrimX)
          Call IClear(NPrim,V(ITch1))
          Call IClear(NAtoms,V(ITch2))
          Call IClear(NAtoms,V(IDrive))
          NEval  = 0
          NEval1 = 0
          NEval2 = 0
          LPos1  = 0
          LPos2  = 0
          Call CrTree(IOut,IPrint,14,IPrim,0,ISet,NPrim,NDepMx,V(ITch1),
     $      V(ITch2),V(IDrive),NEval,NEval1,NEval2,LPos1,LPos2,LogTmp,
     $      InfCrd,IVars,V(IIScr),IILen)
          Do 300 IEval = 1,NEval
            IAtX = IIGet(V(IDrive),IEval)
            NumIq2 = 3
            IdqLst(1) = -(10*IAtX+1)
            IdqLst(2) = -(10*IAtX+2)
            IdqLst(3) = -(10*IAtX+3)
            Call CrdDq0(IOut,IPrint,ISet,IPrim,IPrimX,IdqLst,NumIq2,
     $        MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $        MRVars,InfCrd,IVars,RVars,NCrdSt,V(IIScr),IILen)
 300        Continue
          endIf
 100    Continue
      Return
      End
*Deck CrDfFx
      Subroutine CrDfFx(IOut,IPrint,IType,Chk0,DoMtM1,MultM1,IPrmEq,
     $  IAddVr)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to fix coordinate definitions so that a
C     canonical definition is used. This helps ensure that the minimum
C     number of definition entries necessary. Given an IType value, IAddVr
C     may be modified to achieve this goal. In the code below, the
C     conventions used are commented.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     IType       The coordinate type flag.
C     DoMtM1      This is an input logical argument. It indicates whether
C                 or not this routine is allowed to put the coordinate
C                 definition in a form that is equal to -1 multiplied by
C                 the input definition.
C     MultM1       This is an output logical argument. It is returned as
C                 .True. if the returned coordinate definition is equal to
C                 -1 multiplied by the input configuration of the
C                 definition. Note that if DoMtM1 is sent as FALSE, this
C                 argument is always returned as FALSE.
C     IAddVr      The list of integer values belonging to the coordinate
C                 being tested.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
      Logical Chk0,DoMtM1,MultM1,DoFix,VcFlip,VcSwap
      Dimension IAddVr(*)
      Save Zero, DoFix
      Data Zero/0.0d0/, DoFix/.False./
 2000 Format(1x,'CrDfFx is putting coordinate def into canonical form.')
 2100 Format(1x,'CrDfFx: Canonical form of coordinate requires ',
     $  'multiplication with -1.')
 2110 Format(1x,'CrDfFx: Cannot put coordinate definition in full ',
     $'canonical form because DoMtM1=.False.')
C
C     Initialize things, then check DoFix via Routine CrStOp. If we are not
C     to do the fix simply exit leaving things unchanged.
C
      ITemp = 0
      RTemp = Zero
      MultM1 = .False.
      IPrmEq = 0
      Call CrStOp(IOut,IPrint,2,13,ITemp,RTemp,DoFix)
      If(.not.DoFix) Return
      I0 = 0
      If(Chk0) Call CrDfC0(IOut,IPrint,.False.,0,I0)
C
C     Using IType, carrying out tests on the ordering of entries in IAddVr
C     and modify them as needed.
C
C     Bond stretch definition: R(IAtom1,IAtom2); IAtom1 < IAtom2.
      If(IType.eq.1) then
        If(IAddVr(2).gt.IAddVr(3)) then
          If(IPrint.ge.5) Write(IOut,2000)
          ITemp1 = IAddVr(2)
          IAddVr(2) = IAddVr(3)
          IAddVr(3) = ITemp1
          endIf
C
C     Angle bend definition: A(IAtom1,IAtom2,IAtom3); IAtom1 < IAtom2.
      else if(IType.eq.2) then
        If(IAddVr(2).gt.IAddVr(4)) then
          If(IPrint.ge.5) Write(IOut,2000)
          ITemp1 = IAddVr(2)
          IAddVr(2) = IAddVr(4)
          IAddVr(4) = ITemp1
          endIf
      else if(IType.eq.3) then
        If(IAddVr(2).gt.IAddVr(5)) then
          If(IPrint.ge.5) Write(IOut,2000)
          ITemp1 = IAddVr(2)
          IAddVr(2) = IAddVr(5)
          IAddVr(5) = ITemp1
          ITemp1 = IAddVr(3)
          IAddVr(3) = IAddVr(4)
          IAddVr(4) = ITemp1
          endIf
      else if(IType.eq.50) then
        If(IAddVr(1).gt.IAddVr(2)) then
          If(IPrint.ge.5) Write(IOut,2000)
          ITemp1 = IAddVr(1)
          IAddVr(1) = IAddVr(2)
          IAddVr(2) = ITemp1
          endIf
      else if(IType.eq.52) then
        If(IAddVr(1).gt.IAddVr(2)) then
          If(IPrint.ge.5) Write(IOut,2000)
          ITemp1 = IAddVr(1)
          IAddVr(1) = IAddVr(2)
          IAddVr(2) = ITemp1
          endIf
        If((IAddVr(1).eq.I0).or.(IAddVr(2).eq.I0)) IPrmEq=I0
      else if(IType.eq.100) then
        If(IAddVr(1).le.IAddVr(2)) then
          IVcLo1 = 1
          IVcHi1 = 2
        else
          IVcLo1 = 2
          IVcHi1 = 1
          endIf
        If(IAddVr(3).le.IAddVr(4)) then
          IVcLo2 = 3
          IVcHi2 = 4
        else
          IVcLo2 = 4
          IVcHi2 = 3
          endIf
        VcFlip = (IAddVr(IVcLo1).eq.IAddVr(IVcLo2).and.
     $    IAddVr(IVcHi1).gt.IAddVr(IVcHi2)).or.
     $    (IAddVr(IVcLo1).gt.IAddVr(IVcLo2))
        If(VcFlip) then
          ITemp1 = IAddVr(1)
          IAddVr(1) = IAddVr(3)
          IAddVr(3) = ITemp1
          ITemp1 = IAddVr(2)
          IAddVr(2) = IAddVr(4)
          IAddVr(4) = ITemp1
          VcSwap = IVcLo2.gt.IVcHi2
        else
          VcSwap = IVcLo1.gt.IVcHi1
          endIf
        If(VcSwap) then
          ITemp1 = IAddVr(1)
          IAddVr(1) = IAddVr(2)
          IAddVr(2) = ITemp1
          ITemp1 = IAddVr(3)
          IAddVr(3) = IAddVr(4)
          IAddVr(4) = ITemp1
          endIf
        If((VcFlip.or.VcSwap).and.IPrint.ge.5) Write(IOut,2000)
        If(DoMtM1.and.IAddVr(3).gt.IAddVr(4)) then
          ITemp1 = IAddVr(3)
          IAddVr(3) = IAddVr(4)
          IAddVr(4) = ITemp1
          MultM1 = .True.
          If(IPrint.ge.5) Write(IOut,2100)
        else if(IAddVr(3).gt.IAddVr(4)) then
          Write(IOut,2110)
          endIf
C
C     Cartesian scalar difference: Cart(Ixyz,IAtom1)-Cart(Ixyz,IAtom2);
C     IAtom1 < IAtom2.
      else if(IType.eq.110) then
        If(DoMtM1.and.IAddVr(2).gt.IAddVr(3)) then
          ITemp1 = IAddVr(2)
          IAddVr(2) = IAddVr(3)
          IAddVr(3) = ITemp1
          MultM1 = .True.
          If(IPrint.ge.5) Write(IOut,2100)
        else if(IAddVr(2).gt.IAddVr(3)) then
          Write(IOut,2110)
          endIf
        endIf
      Return
      End
*Deck CrDfIC
      Subroutine CrDfIC(IOut,IPrint,IAtom1,IAtom2,IAtom3,IAtom4,DoInvR,
     $  IAct,ISet,MaxSet,INwPrm,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,
     $  MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine adds one standard internal coordinate and associated
C     definitions to the coordinate definition data structure, and is part
C     of the CrDfxx group of routines.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     IAtom1-4    These are the numbers of the 2, 3, or 4 atomic centers
C                 that define the internal coordinate to be added.
C     DoInvR      A logical input argument that indicates whether or not
C                 bond stretch terms should be built as inverse bond
C                 stretch definitions.
C     ISet        This input argument indicates to which coordinate set
C                 this derivative is added.
C     MaxSet      The maximum number of definition sets allowed.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Logical DoInvR,AddEq,DoFix,Exists,OK
      Dimension ITemp(MxI),RTemp(MxR),InfCrd(*),IVars(*),RVars(*),
     $  NCrdSt(*)
      Save Zero,One
      Data Zero/0.d0/,One/1.d0/
 2010 Format(1x,'CrDfIC: Adding bond between atoms ',I6,' and ',I6,'.')
 2020 Format(1x,'CrDfIC: Adding bend between atoms ',I6,',',I6,', and ',
     $  I6,'.')
 2030 Format(1x,'CrDfIC: Adding torsion between atoms ',I6,',',I6,',',
     $  I6,' and ',I6,'.')
 2040 Format(1x,'CrDfIC: Adding linear bend for atoms ',I6,',',I6,
     $  ', and ',I6,' using 4th atom ',I6,'.')
 9000 Format(1x,'CrDfIC: IAtom1=',I6,' IAtom2=',I6,' IAtom3=',I6,
     $  ' IAtom4=',I6)
 9010 Format(1x,'CrDfIC: Linear Bend Requested...LAtom4=',I6,' LinTyp=',
     $  I6)
 9100 Format(1x,'CrDfIC: ICTyp=',I3)
C
      AddEq = .False.
      DoFix = .True.
C
C     Do some error checking. Then use IAtom1-IAtom4 to set ICTyp, which encodes
C     the type of internal coordinate being added. ICTyp values include:
C                 1 ... Bond stretch
C                 2 ... Bond bend
C                 3 ... Bond torsion
C                 4 ... Linear bend that uses a fourth atom in the
C                       definition.
      If(IAtom1.le.0.or.IAtom2.le.0.or.IAtom3.lt.0) then
        Write(IOut,9000) IAtom1,IAtom2,IAtom3,IAtom4
        Call GauErr('CrDfIC: Invalid IAtom sent < 0.')
        endIf
      ICTyp = 0
      If(IAtom3.eq.0) then
        If(IAtom1.eq.IAtom2) then
          Write(IOut,9000) IAtom1,IAtom2,IAtom3,IAtom4
          Call GauErr('CrDfIC: IAtom1=IAtom2.')
          endIf
        ICTyp = 1
      else if(IAtom4.eq.0) then
        If(IAtom1.eq.IAtom2.or.IAtom1.eq.IAtom3.or.
     $    IAtom2.eq.IAtom3) then
          Write(IOut,9000) IAtom1,IAtom2,IAtom3,IAtom4
          Call GauErr('CrDfIC: 2 or more of the centers are the same.')
          endIf
        ICTyp = 2
      else if(IAtom4.gt.0) then
        If(IAtom1.eq.IAtom2.or.IAtom1.eq.IAtom3.or.IAtom1.eq.IAtom4.or.
     $    IAtom2.eq.IAtom3.or.IAtom2.eq.IAtom4.or.
     $    IAtom3.eq.IAtom4) then
          Write(IOut,9000) IAtom1,IAtom2,IAtom3,IAtom4
          Call GauErr('CrDfIC: 2 or more of the centers are the same.')
          endIf
        ICTyp = 3
      else if(IAtom4.lt.0) then
        LAtom4 = -IAtom4
        LinTyp = Mod(LAtom4,10)
        LAtom4 = LAtom4/10
        If(IAtom1.eq.IAtom2.or.IAtom1.eq.IAtom3.or.IAtom1.eq.LAtom4.or.
     $    IAtom2.eq.IAtom3.or.IAtom2.eq.LAtom4.or.
     $    IAtom3.eq.LAtom4) then
          Write(IOut,9000) IAtom1,IAtom2,IAtom3,IAtom4
          Write(IOut,9010) LAtom4,LinTyp
          Call GauErr('CrDfIC: 2 or more of the centers are the same.')
          endIf
        ICTyp = 4
      else
        Call GauErr('CrDfIC: IC: Logic Error')
        endIf
C
C     Using ICTyp and the other option flags passed down from the calling
C     routine, add the appropriate coordinate definition(s).
C
C     Add stretch definition...
      If(ICTyp.eq.1) then
        If(IPrint.ge.4) Write(IOut,2010) IAtom1,IAtom2
        Call CrDfBd(IOut,IPrint,0,DoInvR,IAtom1,IAtom2,0,0,IAct,ISet,
     $    MaxSet,INwPrm,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $    MRVars,InfCrd,IVars,RVars,NCrdSt)
C
C     Add bend definition...
      else if(ICTyp.eq.2) then
        If(IPrint.ge.4) Write(IOut,2020) IAtom1,IAtom2,IAtom3
        Call CrDfAn(IOut,IPrint,0,IAtom1,IAtom2,IAtom3,0,0,IAct,ISet,
     $    MaxSet,INwPrm,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $    MRVars,InfCrd,IVars,RVars,NCrdSt)
C
C     Add torsion definition...
      else if(ICTyp.eq.3) then
        If(IPrint.ge.4) Write(IOut,2030) IAtom1,IAtom2,IAtom3,IAtom4
        Call CrDfDi(IOut,IPrint,0,IAtom1,IAtom2,IAtom3,IAtom4,0,IAct,
     $    ISet,MaxSet,INwPrm,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,
     $    MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
C
C     Add linear bend definition...
      else if(ICTyp.eq.4) then
        If(IPrint.ge.4) Write(IOut,2040) IAtom1,IAtom2,IAtom3,LAtom4
        If(LinTyp.eq.1) then
          JAct = -1
          Call CrDfAn(IOut,IPrint,0,IAtom1,IAtom2,LAtom4,0,0,JAct,ISet,
     $      MaxSet,IPTmp1,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,
     $      MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
          Call CrDfAn(IOut,IPrint,0,LAtom4,IAtom2,IAtom3,0,0,JAct,ISet,
     $      MaxSet,IPTmp2,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,
     $      MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
          IType = ICTSum
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp1
          ITemp(2) = IPTmp2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,INwPrm,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfIC: Failure adding definition.')
        else if(LinTyp.eq.2) then
          TwoPi = GFloat(8)*GAtan(One)
          JAct = -1
          Call CrDfDi(IOut,IPrint,0,IAtom1,IAtom2,LAtom4,IAtom3,0,JAct,
     $      ISet,MaxSet,IPTmp1,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,
     $      MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
          IType = ICTCdL
          NIVars = 1
          NRVars = 2
          ITemp(1) = IPTmp1
          RTemp(1) = Zero
          RTemp(2) = TwoPi
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfIC: Failure adding definition.')
          Call CrDfDi(IOut,IPrint,0,IAtom1,LAtom4,IAtom2,IAtom3,0,JAct,
     $      ISet,MaxSet,IPTmp2,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,
     $      MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
          IType = ICTCdL
          NIVars = 1
          NRVars = 2
          ITemp(1) = IPTmp2
          RTemp(1) = Zero
          RTemp(2) = TwoPi
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp2,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfIC: Failure adding definition.')
          IType = ICTDif
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp1
          ITemp(2) = IPTmp2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfIC: Failure adding definition.')
          Call CrDfCn(IOut,IPrint,11,ISet,MaxSet,IPTmp2,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
          IType = ICTPrd
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp1
          ITemp(2) = IPTmp2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfIC: Failure adding definition.')
          Call CrDfCn(IOut,IPrint,100,ISet,MaxSet,IPTmp2,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
          IType = ICTSum
          NIVars = 2
          NRVars = 0
          ITemp(1) = IPTmp1
          ITemp(2) = IPTmp2
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPTmp1,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,JAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfIC: Failure adding definition.')
          IType = ICTCdL
          NIVars = 1
          NRVars = 2
          ITemp(1) = IPTmp1
          RTemp(1) = Zero
          RTemp(2) = TwoPi
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,INwPrm,OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $      IVars,RVars,NCrdSt)
          If(.not.OK)Call GauErr('CrDfIC: Failure adding definition.')
       else
         Call GauErr('CrDfIC: STOP')
         endIf
C
C     Unknown ICTyp...
      else
        Write(IOut,9100) ICTyp
        Call GauErr('CrDfIC: Unknown ICTyp.')
        endIf
      Return
      End
*Deck CrDfIn
      Subroutine CrDfIn(IOut,IPrint,KAct,JType,JDep1,JDep2,JDep3,JDep4,
     $  MaxRec,RecLst,NCol,NRec,RecTab,ISet,MaxSet,NPrim,NCoord,LIVars,
     $  LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
      Implicit Integer(A-Z)
C
C     This routine translates the user-definition table generated in
C     Routine CrDfSy into the necessary calls that build the corresponding
C     internal coordinate data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Real*8 GFloat,RTemp,RVars
      Logical AddEq,DoFix,Exists,OK
      Dimension RecLst(MaxRec),RecTab(NCol,NRec),ITemp(MxI),RTemp(MxR),
     $  InfCrd(*),IVars(*),RVars(*),NCrdSt(*)
 1000 Format(1x,'Enter CrDfIn.  NRec=',I5)
 9000 Format(1x,'CrDfIn: NIVars=',I3,' NRVars=',I3)
C
      If(IPrint.ge.4) Write(IOut,1000) NRec
      AddEq = .False.
      DoFix = .True.
      Call IClear(MxI,ITemp)
      Call AClear(MxR,RTemp)
C
C     Loop over the columns of RecTab and add them one at a time to the
C     data structure.
C
      Call IClear(MaxRec,RecLst)
      Do 100 I = NRec,1,-1
        IAct = -1
        If(I.eq.1) IAct = KAct
C
C       Add a Cartesian coordinate definition.
        If(RecTab(JType,I).eq.ICTxyz) then
          Call CrNiNr(IOut,RecTab(JType,I),NIVars,NRVars)
          If(NIVars.gt.2.or.NRVars.gt.0) then
            Write(IOut,9000) NIVars,NRVars
            Call GauErr('CrDfIn: Problem with NIVars and/or NRVars.')
            endIf
          ITemp(1) = RecTab(JDep1,I)
          ITemp(2) = RecTab(JDep2,I)
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,RecLst(I),OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,RecTab(JType,I),IAct,ISet,NIVars,NRVars,ITemp,RTemp,
     $      InfCrd,IVars,RVars,NCrdSt)
          If(.not.OK) Call GauErr('CrDfIn: Problem adding coordinate.')
C
C       Add a standard internal coordinate definition.
        else if((RecTab(JType,I).eq.ICTBd).or.
     $      (RecTab(JType,I).eq.ICTAng).or.
     $      (RecTab(JType,I).eq.ICTDih)) then
          If(RecTab(JType,I).eq.ICTBd.and.
     $      (RecTab(JDep3,I).ne.0.or.RecTab(JDep4,I).ne.0))
     $      Call GauErr('Bond defined, but more than 2 atoms given.')
          If(RecTab(JType,I).eq.ICTAng.and.RecTab(JDep4,I).ne.0)
     $      Call GauErr('Angle defined, but more than 3 atoms given.')
          Call CrDfIC(IOut,IPrint,RecTab(JDep1,I),RecTab(JDep2,I),
     $      RecTab(JDep3,I),RecTab(JDep4,I),.False.,IAct,ISet,MaxSet,
     $      RecLst(I),NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,InfCrd,IVars,RVars,NCrdSt)
C
C       Add a dot product of atomic centered Cartesian difference vectors.
        else if(RecTab(JType,I).eq.ICTDot) then
          NIVars = 4
          NRVars = 0
          ITemp(1) = RecTab(JDep1,I)
          ITemp(2) = RecTab(JDep2,I)
          ITemp(3) = RecTab(JDep3,I)
          ITemp(4) = RecTab(JDep4,I)
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,RecLst(I),OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,RecTab(JType,I),IAct,ISet,NIVars,NRVars,ITemp,RTemp,
     $      InfCrd,IVars,RVars,NCrdSt)
          If(.not.OK) Call GauErr('CrDfIn: Problem adding coordinate.')
C
C       Add a constant.
        else if(RecTab(JType,I).eq.ICTCon) then
          If((IABS(RecTab(JDep1,I)).eq.0).or.
     $      (IABS(RecTab(JDep1,I)).eq.1).or.
     $      (IABS(RecTab(JDep1,I)).eq.2)) then
            Call CrDfCn(IOut,IPrint,RecTab(JDep1,I),ISet,MaxSet,
     $        RecLst(I),NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $        MRVars,InfCrd,IVars,RVars,NCrdSt)
          else
            NIVars = 0
            NRVars = 1
            RTemp(1) = GFloat(RecTab(JDep1,I))
            Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,RecLst(I),OK,
     $        MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $        MRVars,RecTab(JType,I),IAct,ISet,NIVars,NRVars,ITemp,
     $        RTemp,InfCrd,IVars,RVars,NCrdSt)
            If(.not.OK)
     $        Call GauErr('CrDfIn: Problem adding coordinate.')
            endIf
C
C       Add a generic coordinate definition.
        else
          Call CrNiNr(IOut,RecTab(JType,I),NIVars,NRVars)
          If(NIVars.gt.2.or.NRVars.gt.0) then
            Write(IOut,9000) NIVars,NRVars
            Call GauErr('CrDfIn: Problem with NIVars and/or NRVars.')
            endIf
          If(RecTab(JDep1,I).ge.0) then
            ITemp(1) = RecTab(JDep1,I)
          else
            ITemp(1) = RecLst(-RecTab(JDep1,I))
            If(ITemp(1).le.0) Call GauErr('CrDfIn: Logic Error A1')
            endIf
          If(RecTab(JDep2,I).ge.0) then
            ITemp(2) = RecTab(JDep2,I)
          else
            ITemp(2) = RecLst(-RecTab(JDep2,I))
            If(ITemp(2).le.0) Call GauErr('CrDfIn: Logic Error A2')
            endIf
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,RecLst(I),OK,
     $      MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $      MRVars,RecTab(JType,I),IAct,ISet,NIVars,NRVars,ITemp,RTemp,
     $      InfCrd,IVars,RVars,NCrdSt)
          If(.not.OK) Call GauErr('CrDfIn: Problem adding coordinate.')
        endIf
 100    Continue
C
      Return
      End
*Deck CrDfM1
      Subroutine CrDfM1(IOut,IPrint,IPrim,ISet,MaxSet,NPrim,NCoord,
     $  LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $  NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to replace a definition number (sent as input
C     argument IPrim) by its product with -1.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     ISet        This input argument indicates to which coordinate set
C                 this derivative is added.
C     MaxSet      The maximum number of definition sets allowed.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Logical AddEq,DoFix,Exists,OK
      Dimension ITemp(MxI),RTemp(MxR),InfCrd(*),IVars(*),RVars(*),
     $  NCrdSt(*)
      Save IPrmM1
      Data IPrmM1/0/
 1000 Format(1x,'CrDfM1 is negating primitive number ',I10)
C
C     Print some details and set AddEq and DoFix.
C
      Call IClear(MxI,ITemp)
      Call AClear(MxR,RTemp)
      If(IPrint.ge.5) Write(IOut,1000) IPrim
      AddEq = .False.
      DoFix = .True.
C
C     Do the work.
C
      If(IPrmM1.eq.0) then
        Call CrDfCn(IOut,IPrint,-1,ISet,MaxSet,IPrmM1,NPrim,NCoord,
     $    LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $    NCrdSt)
        endIf
      IType = ICTPrd
      IAct = -1
      NIVars = 2
      NRVars = 0
      ITemp(1) = IPrmM1
      ITemp(2) = IPrim
      Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrim,OK,MaxSet,
     $  NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $  IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
      If(.not.OK)
     $  Call GauErr('CrDfM1: Failure adding product definition.')
      Return
      End
*Deck CrDfSy
      Subroutine CrDfSy(In,IOut,IPrint,ISet,MaxSet,NPrim,NCoord,LIVars,
     $  LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
      Implicit Integer(A-Z)
C
C     This routine reads and interprets coordinate definitions provided by
C     the user in the input file using the symbolic input format.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
      Real*8 RVars
      Parameter (LenLin=1000,MaxTab=3000,MaxLab=32,MaxRec=10,
     $  MaxVar=1000)
      Parameter (IRetrn=1,IAct=IRetrn+1,IType=IAct+1,IDep1=IType+1,
     $  IDep2=IDep1+1,IDep3=IDep2+1,IDep4=IDep3+1,ILabel=IDep4+1,
     $  MxIRes=ILabel+MaxLab)
      Parameter (JType=1,JDep1=JType+1,JDep2=JDep1+1,JDep3=JDep2+1,
     $  JDep4=JDep3+1,MxJRes=JDep4)
      Dimension PrsTab(MaxTab),PrsRes(MxIRes),RecTab(MxJRes,MaxRec),
     $  RecLst(MaxRec),Line(LenLin),
     $  LabNam(MaxVar*(MaxLab+5)),LabTmp(MaxLab),LabNum(MaxVar),
     $  InfCrd(*),IVars(*),RVars(*),NCrdSt(*)
 2000 Format(/,1x,'CrDfSy: After QParse.',/,
     $  3x,'Stat  =',I5,' KRetrn=',I5,' RecCur=',I5,/,
     $  9x,'KType =',I5,' KDep1 =',I5,' KDep2 =',I5,/)
 2200 Format(/,1x,'CrDfSy found definition label:')
 8000 Format(/,1x,'CrDfSy found the following definition labels:')
 8001 Format(/)
 8010 Format(1x,'Label number ',I5)
 8100 Format(1x,'The coordinate defs at the end of CrDfSy:')
 9000 Format(/,1x,'CrDfSy: QParse returned with a failure flag. Stat=',
     $  I5)
 9010 Format(/,1x,'CrDfSy: Maximum number of recursions exceeded.',/,
     $  9x,'RexNxt=',I5,' MaxRec=',I5)
 9100 Format(/,1x,'CrDfSy found coordinate reference:')
 9200 Format(/,1x,'Coordinate definition number ',I5,' has been ',
     $  'reused as coordinate definition number ',I5,'.')
C
C     Set-up the parse table.
C
      Call CrPTab(IOut,MaxTab,PrsTab,IRetrn,IAct,IType,IDep1,IDep2,
     $  IDep3,IDep4,ILabel)
C
C     Parse the user-defined coordinate input taking care to properly
C     handle recurrsions. This code works using a number of return states
C     in the parse table. Before each return to this routine, the parser
C     loads PrsRes(IRetrn) with an integer flag indicating where the return
C     has happened and how the current values in the result array, PrsRes,
C     should be used for building the user-defined coordinate definition
C     table, UsrTab. Known values for PrsRes(IRetrn) include:
C           1     The label for a new definition has been parsed together
C                 with any general options given for this definition.
C           2     A function definition has been found and the user
C                 definition table counters should be set accordingly.
C          10     The first dependency for the current function has been
C                 given symbolically as a label. The symbol is stored in
C                 PrsRes(ILabel+1). The length of the symbol is stored in
C                 PrsRed(ILabel).
C          11     The first dependency for the current function is another
C                 function and the parser is starting a recurrsion
C                 transition.
C          12     The second dependency for the current function is another
C                 function and the parser is starting a recurrsion
C                 transition.
C         999     The end of a function definition has been reached; pop
C                 RecCur.
C
      LCursr = 0
      NLabel = 0
      RecCur = 1
      RecNxt = 1
      DepCur = 0
      KAct = 0
      Call IClear(MxJRes*MaxRec,RecTab)
      Call IClear(MaxRec,RecLst)
      IQPDbg = Max(IPrint-1,0)
      Call StrIn(In,Line,LenLin,IEOF)
      If(IEOF.ne.0) Call GauErr('Symbolic coordinate input is missing.')
 500  Call QPInit(PrsTab,0,0,0,IQPDbg)
 510  Call IClear(MxIRes,PrsRes)
      Stat = QParse(PrsRes,PrsTab,Line,LenLin)
      If(Stat.lt.0.or.Stat.gt.1) then
        Write(IOut,9000) Stat
        Call QPErro(Line,LenLin,IOut)
        Call Lnk1E(0)
        endIf
      KRetrn = PrsRes(IRetrn)
      KType  = PrsRes(IType)
      KDep1  = PrsRes(IDep1)
      KDep2  = PrsRes(IDep2)
      KDep3  = PrsRes(IDep3)
      KDep4  = PrsRes(IDep4)
      If(KDep1.ne.0) RecTab(JDep1,RecCur) = KDep1
      If(KDep2.ne.0) RecTab(JDep2,RecCur) = KDep2
      If(KDep3.ne.0) RecTab(JDep3,RecCur) = KDep3
      If(KDep4.ne.0) RecTab(JDep4,RecCur) = KDep4
      If(KRetrn.eq.1) then
C
C       Call CrDfIn to add the current definitions to the coordinate data
C       structure.
        If(NLabel.gt.0) then
          If(IPrint.ge.2) Call IMtOuS(IOut,'CrDfSy: Final PrsRes A:',0,
     $      RecTab,MxJRes,MaxRec,MxJRes,MaxRec)
          Call CrDfIn(IOut,IPrint,KAct,JType,JDep1,JDep2,JDep3,JDep4,
     $      MaxRec,RecLst,MxJRes,RecNxt-1,RecTab,ISet,MaxSet,NPrim,
     $      NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $      IVars,RVars,NCrdSt)
          LabNum(NLabel) = RecLst(1)
          endIf
        If(IPrint.ge.2)then
          Write(IOut,2200)
          Call StrOut(IOut,PrsRes(ILabel+1),PrsRes(ILabel),1)
          endIf
        Match = LocStN(2,NLabel,PrsRes(ILabel),PrsRes(ILabel+1),LabNam,
     $    LabTmp)
        If(Match.ne.0) then
          Write(IOut,9200) Match,NLabel+1
          Call GauErr('Coordinate labels may only be used once.')
          endIf
        KAct = PrsRes(IAct)
        NLabel = NLabel+1
        Call PutB(PrsRes(ILabel+1),PrsRes(ILabel),LabNam,LCursr)
        Call PutDel(2,LabNam,LCursr)
        RecCur = 1
        RecNxt = 1
        Call IClear(MxJRes*MaxRec,RecTab)
        Call IClear(MaxRec,RecLst)
      else if(KRetrn.eq.2) then
        If(DepCur.ne.0) RecTab(DepCur,RecCur) = -RecNxt
        DepCur = 0
        RecLst(RecNxt) = RecCur
        RecCur = RecNxt
        RecNxt = RecNxt+1
        RecTab(JType,RecCur) = KType
      else if(KRetrn.eq.10) then
        Match = LocStN(2,NLabel,PrsRes(ILabel),PrsRes(ILabel+1),LabNam,
     $    LabTmp)
        If(Match.gt.0) then
          If(DepCur.eq.0)
     $      Call GauErr('Coordinate reference found, but DepCur=0.')
          RecTab(DepCur,RecCur) = LabNum(Match)
          DepCur = 0
        else
          Write(IOut,9100)
          Call StrOut(IOut,PrsRes(ILabel+1),PrsRes(ILabel),1)
          Call GauErr('CrDfSy: Error resolving definition reference.')
          endIf
      else if(KRetrn.eq.11) then
        DepCur = JDep1
      else if(KRetrn.eq.12) then
        DepCur = JDep2
      else if(KRetrn.eq.999) then
        RecCur = RecLst(RecCur)
        endIf
      If(RecNxt.gt.MaxRec+1) then
        Write(IOut,9010) RecNxt,MaxRec
        Call GauErr(
     $    'CrDfSy: Depth of user-defined recurrsion > MaxRec.')
        endIf
      PrsRes(IRetrn) = 0
      If(IPrint.ge.3) Write(IOut,2000) Stat,KRetrn,RecCur,KType,KDep1,
     $  KDep2
C
C     The current definition is done. Check to see if there are still
C     coordinates, so loop back up to read the next user
C     defined coordinate.
C
      If(Stat.eq.1) Goto 510
      Call StrIn(In,Line,LenLin,IEOF)
      If(IEOF.eq.0) Goto 500
C
C     Call CrDfIn to add the current definitions to the coordinate data
C     structure.
C
      If(IPrint.ge.2) Call IMtOuS(IOut,'Final PrsRes B:',0,RecTab,
     $  MxJRes,MaxRec,MxJRes,MaxRec)
      Call CrDfIn(IOut,IPrint,KAct,JType,JDep1,JDep2,JDep3,JDep4,MaxRec,
     $  RecLst,MxJRes,RecNxt-1,RecTab,ISet,MaxSet,NPrim,NCoord,LIVars,
     $  LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
      LabNum(NLabel) = RecLst(1)
C
C     If requested, write out the list of labels and the coordinate
C     definitions data structure.
C
      Write(IOut,8000)
      Do 900 i = 1,NLabel
        Write(IOut,8010) LabNum(i)
        Call GetNSt(.True.,.False.,2,i,LabNam,LabTmp,LabLen)
        Call StrOut(IOut,LabTmp,LabLen,1)
 900    Continue
      Write(IOut,8001)
      If(IPrint.ge.2) then
        Write(IOut,8100)
        Call CrPrt(IOut,IPrint,1,0,NPrim,InfCrd,IVars,RVars)
        endIf
      Return
      End
*Deck CrDfVV
      Subroutine CrDfVV(IOut,IPrint,IDeriv,HlfDrv,Fill0,FillD1,FillD2,
     $  Have0,IAtom1,IAtom2,IAtom3,IAtom4,Idx,Idy,ISet,MaxSet,IPrmVV,
     $  IVVx,IVVxy,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $  MRVars,InfCrd,IVars,RVars,NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to add V12.V34 coordinate definitions and their
C     derivatives wrt to atomic-centered Cartesian coordinates.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     IDeriv      An integer input argument indicating wheter the
C                 appropriate bond derivative of the bond definition should
C                 be added instead of the actual definition itself.
C     IAtom1      The first atom involved in the bond stretch.
C     IAtom2      The second atom involved in the bond stretch.
C     Idy         This input argument gives the second Cartesian coordinate
C                 with respect to which R(IAtom1,IAtom2) is to be
C                 differentiated when IDeriv >= 2. If IDeriv < 2, this
C                 argument is ignored.
C     ISet        This input argument indicates to which coordinate set
C                 this derivative is added.
C     MaxSet      The maximum number of definition sets allowed.
C     IPrmVV      This output integer is filled with the primtive
C                 definition number of the final coordinate. If putting the
C                 definition into canonical form requires multiplication
C                 with -1, then IPrmVV is < 0. Note that this routine does
C                 NOT include the multiplication with -1 in the definition.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Logical HlfDrv,Fill0,FillD1,FillD2,Have0,AddEq,DoFix,Exists,OK,
     $  MultM1,Mx12,Mx34,MVVx
      Dimension ITemp(MxI),RTemp(MxR),InfCrd(*),IVars(*),RVars(*),
     $  NCrdSt(*)
      Save ITemp,RTemp,IAtSv1,IAtSv2,IAtSv3,IAtSv4,IVVSv
      Data ITemp/MxI*0/,RTemp/MxR*0.d0/,IAtSv1/0/,IAtSv2/0/,IAtSv3/0/,
     $  IAtSv4/0/,IVVSv/0/
 1000 Format(1x,'CrDfVV is appending ',
     $  'V(',I5,',',I5,').V(',I5,',',I5,')')
 9000 Format(1x,'CrDfVV: IDeriv=',I2)
C
C     Print some details and do some basic error checking.
C
      If(IDeriv.lt.0.or.IDeriv.gt.2) then
        Write(IOut,9000) IDeriv
        Call GauErr('CrDfVV: Invalid IDeriv.')
        endIf
C
C     Initialize key settings and arguments.
C
      AddEq = .False.
      DoFix = .True.
      MultM1 = .False.
      IVVx = 0
      IVVxy = 0
C
C     Add the requested coordinate. We try to put the indices in canonical
C     order. If doing so changes the sign of the coordinate then the final
C     returned value of the primitive is given as -1 times its real value.
C
      IType = 100
      IAct = -1
      NIVars = 4
      NRVars = 0
      ITemp(1) = IAtom1
      ITemp(2) = IAtom2
      ITemp(3) = IAtom3
      ITemp(4) = IAtom4
      If(DoFix) Call CrDfFx(IOut,IPrint,IType,.False.,.True.,MultM1,
     $  Junk,ITemp)
      If(ITemp(1).eq.IAtSv1.and.ITemp(2).eq.IAtSv2.and.
     $  ITemp(3).eq.IAtSv3.and.ITemp(4).eq.IAtSv4) then
        IPrmVV = IVVSv
        If(IPrmVV.lt.0.and.Fill0) then
          IPrmVV = -IPrmVV
          Call CrDfM1(IOut,IPrint,IPrmVV,ISet,MaxSet,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
          endIf
      else if(.not.Have0) then
        If(IPrint.ge.5) Write(IOut,1000) IAtom1,IAtom2,IAtom3,IAtom4
        Call CrdAdd(IOut,IPrint,AddEq,0,0,.False.,Exists,IPrmVV,OK,
     $    MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,
     $    IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,
     $    NCrdSt)
        If(.not.OK) Call GauErr('CrDfVV: Failure adding def A.')
        If(MultM1.and.Fill0) then
          Call CrDfM1(IOut,IPrint,IPrmVV,ISet,MaxSet,NPrim,NCoord,
     $      LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,
     $      NCrdSt)
        else if(MultM1) then
          IPrmVV = -IPrmVV
          endIf
        IAtSv1 = ITemp(1)
        IAtSv2 = ITemp(2)
        IAtSv3 = ITemp(3)
        IAtSv4 = ITemp(4)
        IVVSv = IPrmVV
        endIf
C
C     Set-up KD terms for first- and second-derivatives, as requested.
C
      If(IDeriv.ge.1) then
        IdXAtm = -Idx/10
        IdX1 = Mod(-Idx,10)
        IdX1 = -IdX1
        IKD1x = 0
        IKD2x = 0
        IKD3x = 0
        IKD4x = 0
        If(IdXAtm.eq.IAtom1) IKD1x = 1
        If(IdXAtm.eq.IAtom2) IKD2x = 1
        If(IdXAtm.eq.IAtom3) IKD3x = 1
        If(IdXAtm.eq.IAtom4) IKD4x = 1
        IDf12x = IKD1x-IKD2x
        IDf34x = IKD3x-IKD4x
        endIf
      If(IDeriv.ge.2) then
        IdYAtm = -Idy/10
        IdY1 = Mod(-Idy,10)
        IdY1 = -IdY1
        IKD1y = 0
        IKD2y = 0
        IKD3y = 0
        IKD4y = 0
        If(IdYAtm.eq.IAtom1) IKD1y = 1
        If(IdYAtm.eq.IAtom2) IKD2y = 1
        If(IdYAtm.eq.IAtom3) IKD3y = 1
        If(IdYAtm.eq.IAtom4) IKD4y = 1
        IDf12y = IKD1y-IKD2y
        IDf34y = IKD3y-IKD4y
        endIf
C
C     Form the first-derivative.
C
      If(IDeriv.eq.1) then
        Ix12 = 0
        Ix34 = 0
        If(IDf34x.ne.0) then
          IType = 110
          IAct = -1
          NIVars = 3
          NRVars = 0
          ITemp(1) = IdX1
          ITemp(2) = IAtom1
          ITemp(3) = IAtom2
          Call CrDfFx(IOut,IPrint,IType,.False.,.True.,MultM1,Junk,
     $      ITemp)
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,Ix12,OK,MaxSet,
     $      NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $      IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,
     $      NCrdSt)
          If(.not.OK) Call GauErr('CrDfVV: Failure adding def B.')
          Mx12 = MultM1.eqv.(IDf34x.eq.1)
          endIf
        If(IDf12x.ne.0) then
          IType = 110
          IAct = -1
          NIVars = 3
          NRVars = 0
          ITemp(1) = IdX1
          ITemp(2) = IAtom3
          ITemp(3) = IAtom4
          Call CrDfFx(IOut,IPrint,IType,.False.,.True.,MultM1,Junk,
     $      ITemp)
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,Ix34,OK,MaxSet,
     $      NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $      IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,
     $      NCrdSt)
          If(.not.OK) Call GauErr('CrDfVV: Failure adding def C.')
          Mx34 = MultM1.eqv.(IDf12x.eq.1)
          endIf
C
        If(Ix12.ne.0.and.Ix34.eq.0) then
          IVVx = Ix12
          MVVx = Mx12
        else if(Ix12.eq.0.and.Ix34.ne.0) then
          IVVx = Ix34
          MVVx = Mx34
        else if(Ix12.eq.0.and.Ix34.eq.0) then
          IVVx = 0
          MVVx = .False.
        else if(Ix12.eq.Ix34) then
          If(.not.HlfDrv) then
            Call CrDfCn(IOut,IPrint,2,ISet,MaxSet,IPrm2,NPrim,NCoord,
     $        LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,
     $        RVars,NCrdSt)
            IType = 52
            IAct = -1
            NIVars = 2
            NRVars = 0
            ITemp(1) = IPrm2
            ITemp(2) = Ix12
            Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IVVx,OK,
     $        MaxSet,NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,
     $        MRVars,IType,IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,
     $        IVars,RVars,NCrdSt)
            If(.not.OK)
     $        Call GauErr('CrDfVV: Failure adding def D.')
          else
            IVVx = Ix12
            endIf
          MVVx = Mx12
        else if(Mx12.eqv.Mx34) then
          IType = 50
          IAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = Ix12
          ITemp(2) = Ix34
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IVVx,OK,MaxSet,
     $      NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $      IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,
     $      NCrdSt)
          If(.not.OK) Call GauErr('CrDfVV: Failure adding def E.')
          MVVx = Mx12
        else if(Mx12) then
          IType = 51
          IAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = Ix34
          ITemp(2) = Ix12
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IVVx,OK,MaxSet,
     $      NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $      IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,
     $      NCrdSt)
          If(.not.OK) Call GauErr('CrDfVV: Failure adding def F.')
          MVVx = .False.
        else if(Mx34) then
          IType = 51
          IAct = -1
          NIVars = 2
          NRVars = 0
          ITemp(1) = Ix12
          ITemp(2) = Ix34
          Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IVVx,OK,MaxSet,
     $      NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $      IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,
     $      NCrdSt)
          If(.not.OK) Call GauErr('CrDfVV: Failure adding def G.')
          MVVx = .False.
        else
          Call GauErr('CrDfVV: Logic Error B')
          endIf
C
        If(FillD1) then
          If(IVVx.eq.0) then
            Call CrDfCn(IOut,IPrint,0,ISet,MaxSet,IVVx,NPrim,NCoord,
     $        LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,
     $        RVars,NCrdSt)
          else if(MVVx) then
            Call CrDfM1(IOut,IPrint,IVVx,ISet,MaxSet,NPrim,NCoord,
     $        LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,
     $        RVars,NCrdSt)
            endIf
        else if(MVVx) then
          IVVx = -IVVx
          endIf
C
C     Form the second-derivative.
C
      else if(IDeriv.eq.2) then
        IVVxy = 0
        If(IdX1.eq.IdY1) IVVxy = IDf12x*IDf34y + IDf34x*IDf12y
        If(HlfDrv) then
          If(IAbs(IVVxy).eq.2) then
            IVVxy = Sign(1,IVVxy)
          else if(IAbs(IVVxy).eq.1) then
            IVVxy = Sign(11,IVVxy)
          else if(IVVxy.ne.0) then
            Call GauErr('CrDfVV: Logic Error A')
            endIf
          endIf
        IVVxyT = IVVxy
        If(FillD2) Call CrDfCn(IOut,IPrint,IVVxyT,ISet,MaxSet,IVVxy,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,
     $    IVars,RVars,NCrdSt)
        endIf
C
      Return
      End
*Deck CrDfVX
      Subroutine CrDfVX(IOut,IPrint,IDeriv,IAtom1,IAtom2,IAtom3,IAtom4,
     $  ISet,MaxSet,IPrmX,IPrmY,IPrmZ,NPrim,NCoord,LIVars,LRVars,MPrim,
     $  MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to add the three elements resulting from the
C     cross-product of V12 and V34. The three primitives corresponding to
C     the x, y, and z elements of the resulting vector are returned in
C     arguments IPrmX, IPrmY, and IPrmZ, respectively.
C
C     At present, IDeriv is ignored.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     IDeriv      An integer input argument indicating wheter the
C                 appropriate bond derivative of the bond definition should
C                 be added instead of the actual definition itself.
C     IAtom1      The first atom involved in the bond stretch.
C     IAtom2      The second atom involved in the bond stretch.
C     Idy         This input argument gives the second Cartesian coordinate
C                 with respect to which R(IAtom1,IAtom2) is to be
C                 differentiated when IDeriv >= 2. If IDeriv < 2, this
C                 argument is ignored.
C     ISet        This input argument indicates to which coordinate set
C                 this derivative is added.
C     MaxSet      The maximum number of definition sets allowed.
C     IPrmX       This output integer is filled with the primtive
C                 definition number of the x element of the final vector.
C     IPrmY       This output integer is filled with the primtive
C                 definition number of the y element of the final vector.
C     IPrmZ       This output integer is filled with the primtive
C                 definition number of the z element of the final vector.
C     NPrim       The current number of primitive coordinates.  NPrim is
C                 incremented by this routine.
C     NCoord      The current number of active/propagating coordinates.
C                 NCoord is incremented by this routine.
C     LIVars      This is the current length of the IVars array. LIVars is
C                 incremented by this routine.
C     LRVars      This is the current length of the RVars array. LRVars is
C                 incremented by this routine.
C     MPrim       This is the maximum allowed number of primitive
C                 coordinates.
C     MCoord      This is the maximum allowed number of active/propagating
C                 coordinates.
C     MIVars      This is the maximum allowed length of the IVars array.
C     MRVars      This is the maximum allowed length of the RVars array.
C
C     InfCrd, IVars, RVars, and NCrdSt are the standard coordinate
C     definition arrays. See Routine CrdDef for a discussion on the
C     coordinate definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Logical AddEq,DoFix,Exists,MX12,MY12,MZ12,MX34,MY34,MZ34,OK
      Dimension ITemp(MxI),RTemp(MxR),InfCrd(*),IVars(*),RVars(*),
     $  NCrdSt(*)
      Save ITemp,RTemp
      Data ITemp/MxI*0/,RTemp/MxR*0.d0/
 1000 Format(1x,'CrDfVX is appending ',
     $  'V(',I5,',',I5,') x V(',I5,',',I5,')')
 9000 Format(1x,'CrDfVX: IDeriv=',I2)
C
C     Print some details and do some basic error checking.
C
      If(IDeriv.lt.0.or.IDeriv.gt.0) then
        Write(IOut,9000) IDeriv
        Call GauErr('CrDfVX: Invalid IDeriv.')
        endIf
      If(IPrint.ge.5) Write(IOut,1000) IAtom1,IAtom2,IAtom3,IAtom4
C
C     Initialize key settings and arguments.
C
      AddEq = .False.
      DoFix = .True.
C
C     Begin by defining the six necessary Cartesian difference vectors.
C
      IType = 110
      IAct = -1
      NIVars = 3
      NRVars = 0
      ITemp(1) = -1
      ITemp(2) = IAtom1
      ITemp(3) = IAtom2
      Call CrDfFx(IOut,IPrint,IType,.False.,.True.,MX12,Junk,ITemp)
      Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPX12,OK,MaxSet,
     $  NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $  IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
      If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      IType = 110
      IAct = -1
      NIVars = 3
      NRVars = 0
      ITemp(1) = -2
      ITemp(2) = IAtom1
      ITemp(3) = IAtom2
      Call CrDfFx(IOut,IPrint,IType,.False.,.True.,MY12,Junk,ITemp)
      Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPY12,OK,MaxSet,
     $  NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $  IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
      If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      IType = 110
      IAct = -1
      NIVars = 3
      NRVars = 0
      ITemp(1) = -3
      ITemp(2) = IAtom1
      ITemp(3) = IAtom2
      Call CrDfFx(IOut,IPrint,IType,.False.,.True.,MZ12,Junk,ITemp)
      Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPZ12,OK,MaxSet,
     $  NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $  IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
      If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      IType = 110
      IAct = -1
      NIVars = 3
      NRVars = 0
      ITemp(1) = -1
      ITemp(2) = IAtom3
      ITemp(3) = IAtom4
      Call CrDfFx(IOut,IPrint,IType,.False.,.True.,MX34,Junk,ITemp)
      Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPX34,OK,MaxSet,
     $  NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $  IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
      If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      IType = 110
      IAct = -1
      NIVars = 3
      NRVars = 0
      ITemp(1) = -2
      ITemp(2) = IAtom3
      ITemp(3) = IAtom4
      Call CrDfFx(IOut,IPrint,IType,.False.,.True.,MY34,Junk,ITemp)
      Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPY34,OK,MaxSet,
     $  NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $  IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
      If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      IType = 110
      IAct = -1
      NIVars = 3
      NRVars = 0
      ITemp(1) = -3
      ITemp(2) = IAtom3
      ITemp(3) = IAtom4
      Call CrDfFx(IOut,IPrint,IType,.False.,.True.,MZ34,Junk,ITemp)
      Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPZ34,OK,MaxSet,
     $  NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $  IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
      If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
C
C     Form the six necessary products.
C
      IType = 52
      IAct = -1
      NIVars = 2
      NRVars = 0
      ITemp(1) = IPX12
      ITemp(2) = IPY34
      Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPXY,OK,MaxSet,
     $  NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $  IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
      If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      IType = 52
      IAct = -1
      NIVars = 2
      NRVars = 0
      ITemp(1) = IPX12
      ITemp(2) = IPZ34
      Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPXZ,OK,MaxSet,
     $  NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $  IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
      If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      IType = 52
      IAct = -1
      NIVars = 2
      NRVars = 0
      ITemp(1) = IPY12
      ITemp(2) = IPX34
      Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPYX,OK,MaxSet,
     $  NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $  IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
      If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      IType = 52
      IAct = -1
      NIVars = 2
      NRVars = 0
      ITemp(1) = IPY12
      ITemp(2) = IPZ34
      Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPYZ,OK,MaxSet,
     $  NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $  IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
      If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      IType = 52
      IAct = -1
      NIVars = 2
      NRVars = 0
      ITemp(1) = IPZ12
      ITemp(2) = IPX34
      Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPZX,OK,MaxSet,
     $  NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $  IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
      If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      IType = 52
      IAct = -1
      NIVars = 2
      NRVars = 0
      ITemp(1) = IPZ12
      ITemp(2) = IPY34
      Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPZY,OK,MaxSet,
     $  NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $  IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
      If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
C
C     Put everything together to form IPrmX, IPrmY, and IPrmZ.
C
      If((MY12.eqv.MZ34).and.(MZ12.eqv.MY34)) then
        IType = 51
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPYZ
        ITemp(2) = IPZY
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmX,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      else if((MY12.neqv.MZ34).and.(MZ12.neqv.MY34)) then
        IType = 51
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPZY
        ITemp(2) = IPYZ
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmX,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      else if((MY12.eqv.MZ34).and.(MZ12.neqv.MY34)) then
        IType = 50
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPYZ
        ITemp(2) = IPZY
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmX,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      else if((MY12.neqv.MZ34).and.(MZ12.eqv.MY34)) then
        IType = 50
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPYZ
        ITemp(2) = IPZY
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmX,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
        Call CrDfM1(IOut,IPrint,IPrmX,ISet,MaxSet,NPrim,NCoord,LIVars,
     $    LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
      else
        Call GauErr('CrDfVX: Logic error A')
        endIf
      If((MZ12.eqv.MX34).and.(MX12.eqv.MZ34)) then
        IType = 51
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPZX
        ITemp(2) = IPXZ
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmY,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      else if((MZ12.neqv.MX34).and.(MX12.neqv.MZ34)) then
        IType = 51
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPXZ
        ITemp(2) = IPZX
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmY,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      else if((MZ12.eqv.MX34).and.(MX12.neqv.MZ34)) then
        IType = 50
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPZX
        ITemp(2) = IPXZ
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmY,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      else if((MZ12.neqv.MX34).and.(MX12.eqv.MZ34)) then
        IType = 50
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPZX
        ITemp(2) = IPXZ
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmY,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
        Call CrDfM1(IOut,IPrint,IPrmY,ISet,MaxSet,NPrim,NCoord,LIVars,
     $    LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
      else
        Call GauErr('CrDfVX: Logic error B')
        endIf
      If((MX12.eqv.MY34).and.(MY12.eqv.MX34)) then
        IType = 51
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPXY
        ITemp(2) = IPYX
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmZ,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      else if((MX12.neqv.MY34).and.(MY12.neqv.MX34)) then
        IType = 51
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPYX
        ITemp(2) = IPXY
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmZ,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      else if((MX12.eqv.MY34).and.(MY12.neqv.MX34)) then
        IType = 50
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPXY
        ITemp(2) = IPYX
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmZ,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
      else if((MX12.neqv.MY34).and.(MY12.eqv.MX34)) then
        IType = 50
        IAct = -1
        NIVars = 2
        NRVars = 0
        ITemp(1) = IPXY
        ITemp(2) = IPYX
        Call CrdAdd(IOut,IPrint,AddEq,0,0,DoFix,Exists,IPrmZ,OK,MaxSet,
     $    NPrim,NCoord,LIVars,LRVars,MPrim,MCoord,MIVars,MRVars,IType,
     $    IAct,ISet,NIVars,NRVars,ITemp,RTemp,InfCrd,IVars,RVars,NCrdSt)
        If(.not.OK) Call GauErr('CrDfVX: Failure adding definition.')
        Call CrDfM1(IOut,IPrint,IPrmZ,ISet,MaxSet,NPrim,NCoord,LIVars,
     $    LRVars,MPrim,MCoord,MIVars,MRVars,InfCrd,IVars,RVars,NCrdSt)
      else
        Call GauErr('CrDfVX: Logic error C')
        endIf
      Return
      End
*Deck CrExst
      Subroutine CrExst(IOut,IPrint,IOp,NPrim,ISkipI,ISkipR,IExist,
     $  Insert,InfCrd,IVars,RVars,JType,JSet,NIVars,NRVars,IAddVr,
     $  RAddVr)
      Implicit Real*8(A-H,O-Z)
C
C     This subroutine is used to check on the existence of a coordinate
C     definition in the current data structure. The argument IExist is
C     returned as 0 if the coordinate does NOT already exist in the
C     coordinate definition data strcture. If, on output, IExist>0 the
C     current coordinate is the same as coordinate definition number
C     IExist.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     IOp         This is an option flag. Two values are permitted:
C                   0 ... Do the work as described in the comments above.
C                  -1 ... Instead of searching the coordinate hash tree for
C                         a match with the current primitive, just print
C                         the current tree search/traversal stats. For this
C                         option, all arguments other than IOut, IPrint,
C                         and IOp are ignored.
C                  -2 ... Instead of searching the coordinate hash tree for
C                         a match with the current primitive, just
C                         re-initialize all of the tree search/traversal
C                         stats.  For this option, all arguments other than
C                         IOut, IPrint, and IOp are ignored.
C     NPrim       The number of primitive internal coordinates.
C     ISkipI      If > 0, this element in IAddVr is ignored during
C                 comparisions with definitions already in the data
C                 structure.
C     ISkipR      If > 0, this element in RAddVr is ignored during
C                 comparisions with definitions already in the data
C                 structure.
C     JType       The coordinate type flag to look for.
C     JSet        The coordinate set number to look for. If JSet<=0, the
C                 set number is ignored.
C     NIVars      The number of integer values in array IAddVr.
C     NRVars      The number of real values in array IAddVr.
C     IAddVr      The list of integer values belonging to the coordinate
C                 being tested.
C     RAddVr      The list of real values belonging to the coordinate being
C                 tested.
C
C     InfCrd, IVars, and RVars are the standard coordinate definition
C     arrays. See Routine CrdDef for a discussion on the coordinate
C     definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Parameter (NInfCr=DEFMAXCOORDINFO)
      Logical DEBUG,Match1,Match2
      Integer CrHash
      Dimension InfCrd(NInfCr,*),IVars(*),RVars(*),IAddVr(*),RAddVr(*)
      Save Tol,NCalls,NAdds,NFinds,CTot,CMax,CMin
      Data Tol/1.d-10/,NCalls/0/,NAdds/0/,NFinds/0/,CTot/0.d0/,
     $  CMax/0.d0/,CMin/9.d5/,NLeft/0/,NRight/0/
 1000 Format(1x,'Enter CrExst. IOp=',I2)
 1500 Format(/,1x,'CrExst: Number of adds   =',I20,/,
     $  9x,'Number of finds  =',I20,/,
     $  9x,'Max cost of adds =',F20.8,/,
     $  9x,'Min cost of adds =',F20.8,/,
     $  9x,'Avg cost of adds =',F20.8,/,
     $  9x,'left    moves    =',I20,/,
     $  9x,'right   moves    =',I20,/,
     $  9x,'% left  moves    =',F20.8,/,
     $  9x,'% right moves    =',F20.8,/)
 2000 Format(/,4x,'CrExst - Loop 100: ICount=',I12,'  IHash =',I12,
     $  '  LeftCh=',I12,/,
     $  23x,'IPrCur=',I12,'  IHshCu=',I12,'  RightCh=',I12)
 2010 Format(23x,'IPrNxt=',I12,'  NPrim =',I12/)
 2100 Format(4x,'CrExst hash binary tree search took ',I10,'steps.')
C
      DEBUG = .False.
      If(IPrint.ge.5.or.DEBUG) Write(IOut,1000) IOp
C
C     If IOp indicates work on the internal tree search stats just do that
C     work and return. Otherwise, do the standard hash value tree search.
C
      If(IOp.lt.-2.or.IOp.gt.0)
     $  Call GauErr('CrExst: IOp out of bounds.')
      If(IOp.eq.-1) then
        Write(IOut,1500) NAdds,NFinds,CMax,CMin,CTot/GFloat(NCalls),
     $    NLeft,NRight,
     $    GFloat(100)*(GFloat(NLeft)/GFloat(Max(1,NLeft+NRight))),
     $    GFloat(100)*(GFloat(NRight)/GFloat(Max(1,NLeft+NRight)))
        Return
      else if(IOp.eq.-2) then
        NCalls = 0
        NAdds = 0
        NFinds = 0
        CTot = GFloat(0)
        CMax = GFloat(0)
        CMin = GFloat(900000)
        NLeft = 0
        NRight = 0
        Return
        endIf
C
C     Look for coordinates in the data structure with a type flag of IType.
C     If we find a match, test the list of integers and reals.
C
C     Note that if we are looking for an existing derivative, we neglect
C     positions 4 and 5 in the integer list.
C
      NCalls = NCalls+1
      IExist = 0
      Insert = 0
      If(NPrim.lt.1) Return
      NIVar0 = NIVars
      If(JType.eq.500) NIVar0 = 3
      IHash = CrHash(NIVars,JType,IAddVr)
C
C     Search the hash binary tree for any existing primitives with the same
C     hash value or, if the current primitive is new to the data structure,
C     the primitive value that is the current primitive's parent in the
C     binary tree.
C
C     Loop 100: Do While(.not.Match1.and..not.IPrNxt.eq.0)
      Match1 = .False.
      IPrNxt = 1
      ICount = 0
  100 IPrCur = IPrNxt
      IHshCu = InfCrd(InfHsh,IPrCur)
      ICount = ICount+1
      If(IPrint.ge.5.or.DEBUG) Write(IOut,2000) ICount,IHash,
     $  InfCrd(InfLCh,IPrCur),IPrCur,IHshCu,InfCrd(InfRCh,IPrCur)
      If(ICount.gt.NPrim) Call GauErr('STOP in 100: ICount.ge.NPrim')
      If(IHash.eq.IHshCu) then
        Match1 = .True.
      else if(IHash.lt.IHshCu) then
        IPrNxt = InfCrd(InfLCh,IPrCur)
        NLeft = NLeft+1
      else if(IHash.gt.IHshCu) then
        IPrNxt = InfCrd(InfRCh,IPrCur)
        NRight = NRight+1
        endIf
      If(IPrint.ge.5.or.DEBUG) Write(IOut,2010) IPrNxt,NPrim
      If(.not.Match1.and.IPrNxt.ne.0) Goto 100
      If(IPrint.ge.5.or.DEBUG) Write(IOut,2100) ICount
C
C     If the search of the hash value binary tree did NOT find a match,
C     return to the calling routine after setting the value of Insert.
C
      If(.not.Match1) then
        Insert = IPrCur
        NAdds = NAdds+1
        Ideal = Int(Log(GFloat(NPrim))/Log(GFloat(2)))+1
        Cost = GFloat(ICount)/GFloat(Ideal)
        CTot = CTot+Cost
        If(Cost.gt.CMax) then
          CMax = Cost
        else if(Cost.lt.CMin) then
          CMin = Cost
          endIf
        Return
        endIf
C
C     This block of code is only reached if the search of the hash value
C     binary tree above sucessfully found a match (Match1=T),
C
C     Loop 200: Do While(Match1.and.IPrNxt.ne.0)
      IPrNxt = IPrCur
  200 IPrCur = IPrNxt
      If(JType.eq.InfCrd(InfTyp,IPrCur).and.
     $  (JSet.le.0.or.InfCrd(InfSet,IPrCur).eq.JSet).and.
     $  InfCrd(InfNI,IPrCur).eq.NIVars.and.
     $  InfCrd(InfNR,IPrCur).eq.NRVars)then
        Match2 = .True.
        IPtRVr = InfCrd(InfRVr,IPrCur)
        IPtIVr = InfCrd(InfIVr,IPrCur)
        Do 210 I = 1, NIVar0
  210     Match2 = Match2.and.(IVars(IPtIVr+I-1).eq.IAddVr(I)
     $      .or.I.eq.ISkipI)
        If(Match2) then
          Do 220 I = 1, NRVars
  220       Match2 = Match2.and.(Abs(RVars(IPtRVr+I-1)-RAddVr(I))
     $        .lt.Tol.or.I.eq.ISkipR)
          endIf
        If(Match2) then
          IExist = IPrCur
          Insert = IPrCur
          NFinds = NFinds+1
          Return
          endIf
        endIf
      IPrNxt = InfCrd(InfNxt,IPrCur)
      If(IPrNxt.ne.0) Goto 200
      Insert = IPrCur
      NAdds = NAdds+1
      Ideal = Int(Log(GFloat(NPrim))/Log(GFloat(2)))+1
      Cost = GFloat(ICount)/GFloat(Ideal)
      CTot = CTot+Cost
      If(Cost.gt.CMax) then
        CMax = Cost
      else if(Cost.lt.CMin) then
        CMin = Cost
        endIf
      Return
      End
*Deck CrHash
      Integer Function CrHash(NIVars,IType,IVars)
      Implicit Integer(A-Z)
C
C     Return a hash function for a coordinate element.
C
      Dimension IVars(*)
C
      If(IType.eq.500) then
        IT0 = 127
      else
        IT0 = Mod(IType,128)
        endIf
      IT1 = 0
      IT2 = 0
      IT3 = 0
      IT4 = 0
      Goto (1,10,20,30,40), (Min(NIVars,4)+1)
   40 IT4 = IVars(4)
   30 IT3 = IVars(3)
   20 IT2 = IVars(2)
   10 IT1 = IVars(1)
C
    1 JT0 = IT0
      JT1 = Mod(IAbs(IT1),2**8)
      JT2 = Mod(IAbs(IT2),2**6)
      JT3 = Mod(IAbs(IT3),2**6)
      JT4 = Mod(IAbs(IT4),2**4)
      CrHash = (((JT1*2**7+JT0)*2**4+JT4)*2**6+JT3)*2**6+JT2
      If(Mod(IT1,2).eq.1) CrHash = -CrHash
      Return
      End
*Deck CrMkLb
      Subroutine CrMkLb(IOut,IPrint,Init,IPrim,ICur,Label,InfCrd)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to make a standard label for an internal
C     coordinate within the new coordinate definition data structure.
C
C     INPUT ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     Init        A logical input argument indicating whether (.TRUE.) or
C                 not (.FALSE.) the numerical flags used in generating
C                 coordinate labels should be initialized.
C     IPrim       The primitive number for which the label is to be
C                 generated.
C     ICur        This is the cursor in Label where the new string should
C                 begin (it's like an offset, so writing to the beginning
C                 of Label requires ICur=0). ICur is incremented in this
C                 routine.
C     Label       This is the argument where the new generated label is
C                 placed.
C     InfCrd      This is the standard information array used in the new
C                 coordinate data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2011.
C
#include "crdparams.inc"
      Logical DEBUG,Init
      Parameter (NInfCr=DEFMAXCOORDINFO)
      Dimension Label(*),InfCrd(NInfCr,*)
      Save NCart,NBond,NAng,NDih,Nq
      Data NCart/0/,NBond/0/,NAng/0/,NDih/0/,Nq/0/
 1000 Format(1x,'Enter CrMkLb.')
C
C     Take care of initialization.
C
      DEBUG = .False.
      If(DEBUG.or.IPrint.ge.4) Write(IOut,1000)
      If(Init) then
        NCart = 0
        NBond = 0
        NAng  = 0
        NDih  = 0
        Nq    = 0
        endIf
      If(IPrim.le.0) Return
      IType = InfCrd(InfTyp,IPrim)
C
C     Build the appropriate label.
C
      If(IType.eq.ICTxyz) then
        NCart = NCart+1
        Call PutCSt('Cart',Label,ICur)
        Call DecChr(NCart,Label,ICur)
      else if(IType.eq.ICTBd) then
        NBond = NBond+1
        Call PutCSt('R',Label,ICur)
        Call DecChr(NBond,Label,ICur)
      else if(IType.eq.ICTAng) then
        NAng = NAng+1
        Call PutCSt('A',Label,ICur)
        Call DecChr(NAng,Label,ICur)
      else if(IType.eq.ICTDih) then
        NDih = NDih+1
        Call PutCSt('D',Label,ICur)
        Call DecChr(NDih,Label,ICur)
      else
        Nq = Nq+1
        Call PutCSt('q',Label,ICur)
        Call DecChr(Nq,Label,ICur)
        endIf
      Return
      End
*Deck CrMSyG
      Subroutine CrMSyG(LayCur,NAtoms,Layer,LiAtom,Distan,C)
      Implicit Real*8(A-H,O-Z)
C
C     Create the model system geometry for point IPoint which is at
C     layer LayCur.
C
      Dimension Layer(*), LiAtom(*), Distan(*), C(3,*)
C
      Do 10 I = 1, NAtoms
        If(Layer(I).gt.LayCur.and.LiAtom(I).ne.0.and.
     $    Layer(LiAtom(I)).le.LayCur) then
          C(1,I) = C(1,LiAtom(I))+Distan(I)*(C(1,I)-C(1,LiAtom(I)))
          C(2,I) = C(2,LiAtom(I))+Distan(I)*(C(2,I)-C(2,LiAtom(I)))
          C(3,I) = C(3,LiAtom(I))+Distan(I)*(C(3,I)-C(3,LiAtom(I)))
          endIf
   10   Continue
      Return
      End
*Deck CrNiNr
      Subroutine CrNiNr(IOut,IType,Ni,Nr)
      Implicit Integer(A-Z)
C
C     Given an primitive internal coordinate type code (IType), this
C     routine fills and returns the number of integer values that go into
C     IVars (Ni) and the number of real values that go into RVars (Nr).
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
 9000 Format(1x,'IType=',I5)
C
      If(IType.eq.ICTxyz) then
        Ni = 2
        Nr = 0
      else if(IType.eq.ICTBd) then
        Ni = 3
        Nr = 0
      else if(IType.eq.ICTAng) then
        Ni = 4
        Nr = 0
      else if(IType.eq.ICTDih) then
        Ni = 5
        Nr = 0
      else if(IType.eq.ICTSum) then
        Ni = 2
        Nr = 0
      else if(IType.eq.ICTDif) then
        Ni = 2
        Nr = 0
      else if(IType.eq.ICTPrd) then
        Ni = 2
        Nr = 0
      else if(IType.eq.ICTQuo) then
        Ni = 2
        Nr = 0
      else if(IType.eq.ICTScP) then
        Ni = 1
        Nr = 1
      else if(IType.eq.ICTEq) then
        Ni = 1
        Nr = 0
      else if(IType.eq.57) then
        Ni = 2
        Nr = 0
      else if(IType.eq.ICTCon) then
        Ni = 0
        Nr = 1
      else if(IType.eq.60) then
        Ni = 2
        Nr = 0
      else if(IType.eq.ICTInv) then
        Ni = 1
        Nr = 0
      else if(IType.eq.ICTSqrt) then
        Ni = 1
        Nr = 0
      else if(IType.eq.ICTSq) then
        Ni = 1
        Nr = 0
      else if(IType.eq.ICTCos) then
        Ni = 1
        Nr = 0
      else if(IType.eq.ICTSin) then
        Ni = 1
        Nr = 0
      else if(IType.eq.ICTTan) then
        Ni = 1
        Nr = 0
      else if(IType.eq.ICTACos) then
        Ni = 1
        Nr = 0
      else if(IType.eq.74) then
        Ni = 1
        Nr = 0
      else if(IType.eq.75) then
        Ni = 1
        Nr = 0
      else if(IType.eq.ICTSign) then
        Ni = 2
        Nr = 0
      else if(IType.eq.ICTCdL) then
        Ni = 1
        Nr = 2
      else if(IType.eq.ICTDot) then
        Ni = 4
        Nr = 0
      else if(IType.eq.ICTXDf) then
        Ni = 3
        Nr = 0
      else if(IType.eq.111) then
        Ni = 4
        Nr = 0
      else if(IType.eq.500) then
        Ni = 5
        Nr = 0
      else
        Write(IOut,9000) IType
        Call GauErr('CrNiNr: Unknown IType.')
        endIf
      Return
      End
*Deck CrPHsh
      Subroutine CrPHsh(IOut,NPrim,InfCrd)
      Implicit Real*8(A-H,O-Z)
C
#include "crdparams.inc"
      Parameter (NInfCr=DEFMAXCOORDINFO)
      Dimension InfCrd(NInfCr,NPrim)
 1000 Format(/)
 2000 Format(1x,'Coordinate Definitions Hash:')
 2010 Format(3x,'IPrim=',I10,'   Hash=',I20,'   LChild=',I10,
     $  '   RChild=',I10,'   Next=',I10)
      Write(IOut,1000)
      Write(IOut,2000)
      Do 100 i = 1,NPrim
        Write(IOut,2010) i,InfCrd(InfHsh,i),InfCrd(InfLCh,i),
     $    InfCrd(InfRCh,i),InfCrd(InfNxt,i)
 100    Continue
      Write(IOut,1000)
      Return
      End
*Deck CrPTab
      Subroutine CrPTab(IOut,MaxTab,Table,IRetrn,IAct,IType,IDep1,IDep2,
     $  IDep3,IDep4,ILabel)
      Implicit Integer(A-Z)
C
C     This routine prepares a parse table for reading symbolic user-defined
C     coorindate definitions from the input file.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
#include "crdparams.inc"
      Dimension Table(MaxTab)
C
C     Set-up the parse table beginning with the recurrsion into State Oper.
C
      Call PTInit(IOut,Table,MaxTab,I)
      Call PTStat(IOut,Table,MaxTab,I,'RdDefs')
      Call PTSpec(IOut,Table,MaxTab,I,'Wrd','GotLab',ILabel,0)
      Call PTSpec(IOut,Table,MaxTab,I,'EOL','Exi',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     GotLab state.
C     This state handles the return to the calling routine with the label
C     for the current definition and the general definition options, which
C     are parsed here. Then we transition to State Opers where the
C     operators and their operands are parsed.
C
      Call PTStat(IOut,Table,MaxTab,I,'GotLab')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@LabOps',' ',0,0)
      Call PTSpec(IOut,Table,MaxTab,I,'Nul',' ',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'=','Ret',IRetrn,1)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@Opers','RdDefs',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     LabOps State.
C     This state takes care of optional settings one can select with a
C     definition label.
      Call PTStat(IOut,Table,MaxTab,I,'LabOps')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'(','LabOps1',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,'LabOps1')
      Call PTTran(IOut,Table,MaxTab,I,.True.,')','Exi',0,0)
      Call PTTran(IOut,Table,MaxTab,I,.True.,',','LabOps1',0,0)
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@LabOps3','LabOps1',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,'LabOps3')
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Active','Exi',0,0)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'InActive','Exi',IAct,-1)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     Opers State.
C     This state runs through the different classes of operators.
C
      Call PTStat(IOut,Table,MaxTab,I,'Opers')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@OperXYZ','GtXYZDep',0,0)
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@OperStdIC','GtStdICDep',
     $  0,0)
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@Oper1','Gt1Dep',0,0)
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@Oper2','Gt2Dep',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     OperXYZ State.
C     This state takes care of operators that have only two operands.
C           NOTE: If "CART" is found we never set the type flag since it is
C                 currently defined to be zero. At some point, we should
C                 change the value in crdparams.inc and then the transition
C                 to State 'RtnOp' should include setting position IType to
C                 the new Cartesian type value.
C
      Call PTStat(IOut,Table,MaxTab,I,'OperXYZ')
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Cart','RtnOp',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     OperStdIC State.
C     This state takes care of operators belonging to the standard internal
C     coordiante types.
C
      Call PTStat(IOut,Table,MaxTab,I,'OperStdIC')
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Bond','RtnOp',IType,
     $  ICTBd)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Stretch','RtnOp',IType,
     $  ICTBd)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Angle','RtnOp',IType,
     $  ICTAng)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Bend','RtnOp',IType,
     $  ICTAng)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Dihedral','RtnOp',IType,
     $  ICTDih)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Torsion','RtnOp',IType,
     $  ICTDih)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'DotDifference','RtnOp',
     $  IType,ICTDot)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Constant','RtnOp',IType,
     $  ICTCon)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     Oper1 State.
C     This state takes care of operators that have only one operand.
C
      Call PTStat(IOut,Table,MaxTab,I,'Oper1')
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Inverse','RtnOp',IType,
     $  ICTInv)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'SQRT','RtnOp',IType,
     $  ICTSqrt)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Squared','RtnOp',IType,
     $  ICTSq)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Cosine','RtnOp',IType,
     $  ICTCos)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Sine','RtnOp',IType,
     $  ICTSin)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Tangent','RtnOp',IType,
     $  ICTTan)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'ACosine','RtnOp',IType,
     $  ICTACos)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'ArcCosine','RtnOp',IType,
     $  ICTACos)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     Oper2 State.
C     This state takes care of operators that have only two operands.
C
      Call PTStat(IOut,Table,MaxTab,I,'Oper2')
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Sum','RtnOp',IType,
     $  ICTSum)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Addition','RtnOp',IType,
     $  ICTSum)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Difference','RtnOp',
     $  IType,ICTDif)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Subtraction','RtnOp',
     $  IType,ICTDif)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Product','RtnOp',IType,
     $  ICTPrd)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Multiply','RtnOp',IType,
     $  ICTPrd)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Quotient','RtnOp',IType,
     $  ICTQuo)
      Call PTTran(IOut,Table,MaxTab,I,.False.,'Divide','RtnOp',IType,
     $  ICTQuo)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     GetXYZDep State.
C     This state takes care of picking up the atom number for Cartesian
C     operations.
C
      Call PTStat(IOut,Table,MaxTab,I,'GtXYZDep ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'(',' ',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'I10',' ',IDep1,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,',',' ',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'X',' ',IDep2,-1)
      Call PTTran(IOut,Table,MaxTab,I,.True.,'Y',' ',IDep2,-2)
      Call PTTran(IOut,Table,MaxTab,I,.True.,'Z',' ',IDep2,-3)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@RtnEnd','Exi',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     GtStdICDep State.
C     This state takes care of picking up the atom numbers for standard
C     internal coordinate definitions. At this point, this set of states
C     loads up to four numbers in the results array. It is up to the
C     calling unit to ensure that the correct number of operands has been
C     specified.
C
      Call PTStat(IOut,Table,MaxTab,I,'GtStdICDep ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'(',' ',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'I10',' ',IDep1,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@RtnEnd','Exi',0,0)
      Call PTTran(IOut,Table,MaxTab,I,.True.,',',' ',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'I10',' ',IDep2,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@RtnEnd','Exi',0,0)
      Call PTTran(IOut,Table,MaxTab,I,.True.,',',' ',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'I10',' ',IDep3,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@RtnEnd','Exi',0,0)
      Call PTTran(IOut,Table,MaxTab,I,.True.,',',' ',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'I10',' ',IDep4,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@RtnEnd','Exi',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     RtnOp State.
C     The return state for operators.
C
      Call PTStat(IOut,Table,MaxTab,I,'RtnOp')
      Call PTSpec(IOut,Table,MaxTab,I,'Nul','Ret',IRetrn,2)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'Nul','Exi',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     Gt1Dep states...
C
      Call PTStat(IOut,Table,MaxTab,I,'Gt1Dep')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'(',' ',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'I10','Gt1Dp1',IDep1,0)
      Call PTSpec(IOut,Table,MaxTab,I,'Nul','Ret',IRetrn,11)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@Opers','Gt1Dp1',0,0)
      Call PTSpec(IOut,Table,MaxTab,I,'Wrd',' ',ILabel,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'Nul','Ret',IRetrn,10)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,'Gt1Dp1')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@RtnEnd','Exi',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     Gt2Dep states...
C
      Call PTStat(IOut,Table,MaxTab,I,'Gt2Dep')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'(',' ',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'I10','Gt2Dp1',IDep1,0)
      Call PTSpec(IOut,Table,MaxTab,I,'Nul','Ret',IRetrn,11)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@Opers','Gt2Dp1',0,0)
      Call PTSpec(IOut,Table,MaxTab,I,'Wrd',' ',ILabel,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'Nul','Ret',IRetrn,10)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,'Gt2Dp1')
      Call PTTran(IOut,Table,MaxTab,I,.True.,',',' ',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'I10','Gt2Dp2',IDep2,0)
      Call PTSpec(IOut,Table,MaxTab,I,'Nul','Ret',IRetrn,12)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@Opers','Gt2Dp2',0,0)
      Call PTSpec(IOut,Table,MaxTab,I,'Wrd',' ',ILabel,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'Nul','Ret',IRetrn,10)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,'Gt2Dp2')
      Call PTTran(IOut,Table,MaxTab,I,.True.,'@RtnEnd','Exi',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     RtnEnd state.
C     Transitions to this state happen when the end of an operator
C     definition is reached and we need to return to the calling program
C     unit so that it can do some bookkeeping.
C
      Call PTStat(IOut,Table,MaxTab,I,'RtnEnd')
      Call PTTran(IOut,Table,MaxTab,I,.True.,')','Ret',IRetrn,999)
      Call PTEOS(IOut,Table,MaxTab,I)
C
      Call PTStat(IOut,Table,MaxTab,I,' ')
      Call PTSpec(IOut,Table,MaxTab,I,'Nul','Exi',0,0)
      Call PTEOS(IOut,Table,MaxTab,I)
C
C     End the parse table.
C
      Call PTEnd(IOut,Table,MaxTab,I)
      Return
      End
*Deck CrPtLb
      Subroutine CrPtLb(IOut,IPrint,Init,IPrim,Labels,InfCrd)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to make a standard label for an internal
C     coordinate and load it into the Labels array within the new
C     coordinate definition data structure for later use by other Cr* and
C     Crd* Routines.
C
C     INPUT ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     Init        A logical input argument indicating whether (.TRUE.) or
C                 not (.FALSE.) the cursor in Labels and the numerical
C                 flags used in Routine CrMkLb when generating standard
C                 coordinate labels should be initialized.
C     IPrim       The primitive number for which the label is to be
C                 generated.
C     Labels      This is the Labels array in the new coordinate data
C                 structure.
C     InfCrd      This is the standard information array used in the new
C                 coordinate data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2011.
C
#include "crdparams.inc"
      Logical DEBUG,Init
      Parameter (NInfCr=DEFMAXCOORDINFO)
      Dimension Labels(*),InfCrd(NInfCr,*)
 1000 Format(1x,'Enter CrPtLb.')
 9000 Format(1x,'CrPtLb: Primitive number ',I10,
     $  ' has label cursor value of ',I10,'.')
C
C     Take care of initialization.
C
      DEBUG = .False.
      If(DEBUG.or.IPrint.ge.4) Write(IOut,1000)
      MaxLab = Labels(1)
      LenLab = Labels(2)
      LabCur = Labels(3)
      If(Init) then
        Call CBlank(LenLab,Labels(4))
        LabCur = 0
        Labels(3) = LabCur
        Call CrMkLb(IOut,IPrint,.True.,0,LabCur,Labels(4),InfCrd)
        endIf
      If(IPrim.le.0) Return
C
C     Generate the appropriate label, store the Hollerith string in Labels,
C     and update InfCrd.
C
      If(InfCrd(InfLab,IPrim).ne.-1) then
        Write(IOut,9000) IPrim,InfCrd(InfLab,IPrim)
        Call GauErr('CrPtLb: This primitive already has a label.')
        endIf
      InfCrd(InfLab,IPrim) = LabCur
      Call CrMkLb(IOut,IPrint,.False.,IPrim,LabCur,Labels(4),InfCrd)
      If(LabCur.gt.MaxLab)
     $  Call GauErr('CrPtLb: Labels is out of room @ 1.')
      Call PutDel(1,Labels(4),LabCur)
      If(LabCur.gt.MaxLab)
     $  Call GauErr('CrPtLb: Labels is out of room @ 2.')
C
C     Update the third element of Labels with LabCur.
C
      Labels(3) = LabCur
      Return
      End
*Deck CrReSz
      Subroutine CrReSz(IOut,IPrint,OK,IVer,MxSet0,NCrd0,NPrim0,LIVar0,
     $  LRVar0,MxSet1,NCrd1,NPrim1,LIVar1,LRVar1,InfCrd,IVars,IRVars,
     $  NCrdSt,ICrdSt,IV,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to re-size the coordinate definition data
C     structure in memory (i.e., in V). Using a the "old" size parameters
C     and a "new" set of desired size paramters, the memory pointers for
C     the key arrays used by the data structure are updated. The data is
C     also shifted in memory, though newly added space is NOT initialized.
C     It is assumed that the coordinate data structure should begin at
C     V(1).  Note that this routine does NOT kill the job if sufficient
C     memory is unavailable; instead, logical flag OK is returned as
C     .FALSE.
C
C     This is one of the few routines that actually knows about the order
C     of internal coordinate definitions arrays in memory. For complete
C     definitions regarding the internal coordinate definitions data
C     structure, see Routine CrdDef.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      Extra print flag.
C     OK          A logical flag that is returned as .TRUE. if re-sizing
C                 went alright. If necessary memory is unavailable, OK is
C                 be returned .FALSE.
C     IVer        The redundant internal coordinate version number.
C     MxSet0      The old maximum number of definition sets.
C     NCrd0       The old number of coordinates.
C     NPrim0      The old number of primitive coordinates.
C     LIVar0      The old length of array IVars.
C     LRVar0      The old length of array RVars.
C     MxSet1      The requested new maximum number of definition sets.
C     NCrd1       The requested new number of coordinates.
C     NPrim1      The requested new number of primitive coordinates.
C     LIVar1      The requested new length of array IVars.
C     LRVar1      The requested new length of array RVars.
C     InfCrd      The position in V for the definitions array. This
C                 argument is changed to reflect the re-sizing.
C     IVars       The position in V for the integer variables array. This
C                 argument is changed to reflect the re-sizing.
C     IRVars      The position in V for the real variables array. This
C                 argument is changed to reflect the re-sizing.
C     NCrdSt      The position in V for the array carrying the number of
C                 active/propagating coordinates in each definition set.
C     ICrdSt      The pointer for the list of active/propagating
C                 coordinates ordered by definition set.
C     IV          This argument is returned with the starting position for
C                 remainging space in memory.
C     V           The standard memory array.
C     MDV         The size of V.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
      Logical DEBUG,Expand,OK
      Dimension V(MDV)
 1000 Format(1x,'CrReSz: IVer  =',I10,' MxSet0=',I10,' NCrd0 =',I10,/,
     $  9x,'NPrim0=',I10,' LIVar0=',I10,' LRVar0=',I10,/,
     $  9x,'MxSet1=',I10,' NCrd1 =',I10,' NPrim1=',I10,/,
     $  9x,'LIVar1=',I10,' LRVar1=',I10)
 1100 Format(1x,'CrReSz: Not enough memory to resize the internal',
     $  ' coordinate data structure.',/,9x,'IV1=',I20,' MDV=',I20)
 2000 Format(1x,'CrReSz: Pointers after re-sizing...',/,
     $  9x,'Old: InfCrd=',I10,' IVars =',I10,' IRVars=',I10,/,
     $  14x,'NCrdSt=',I10,' ICrdSt=',I10,' IV    =',I10,/,
     $  9x,'New: InfCrd=',I10,' IVars =',I10,' IRVars=',I10,/,
     $  14x,'NCrdSt=',I10,' ICrdSt=',I10,' IV    =',I10)
C
      DEBUG = .False.
C
C     Check to see if the size is NOT actually changing, in which case we
C     set OK to .TRUE. and return. Otherwise, initialize OK and set Expand.
C
      If(IPrint.ge.4.or.DEBUG) Write(IOut,1000) IVer,MxSet0,NCrd0,
     $  NPrim0,LIVar0,
     $  LRVar0,MxSet1,NCrd1,NPrim1,LIVar1,LRVar1
      If(NPrim0.eq.NPrim1.and.LIVar0.eq.LIVar1.and.LRVar0.eq.LRVar1)
     $  then
        OK = .True.
        Return
        endIf
      OK = .False.
      Expand = NPrim0.lt.NPrim1
      If(Expand.neqv.(LIVar0.le.LIVar1))
     $  Call GauErr('CrReSz: Logical Error 1')
      If(Expand.neqv.(LRVar0.le.LRVar1))
     $  Call GauErr('CrReSz: Logical Error 2')
C
C     Get pointers for the new positions in memory and check for sufficient
C     memory.
C
      Call CrMem(IOut,IPrint,1,IVer,MxSet1,NPrim1,NCrd1,LIVar1,LRVar1,
     $  MaxLa1,LenLa1,InfCr1,IVar1,IRVar1,NCrdS1,ICrdS1,Label1,IV1)
      If(InfCr1.ne.InfCrd) Call GauErr('CrReSz: InfCr1.ne.InfCrd')
      If(IV1.gt.MDV) then
        If(IPrint.ge.4.or.DEBUG) Write(IOut,1100) IV1-1,MDV
        Return
        endIf
C
C     Move the data structure in memory and then re-set the key memory
C     pointers before leavingi. The order of arrays in memory is assumed to
C     be InfCrd, IVars, IRVars, NCrdSt, and ICrdSt. Since InfCrd is at the
C     start of memory, we never need to touch its data.
C
      If(Expand) then
        Call AMove1(InToWP(NCrd0),ICrdSt-1,ICrdS1-1,V)
        Call AMove1(InToWP(MxSet0),NCrdSt-1,NCrdS1-1,V)
        Call AMove1(LRVar0,IRVars-1,IRVar1-1,V)
        Call AMove1(InToWP(LIVar0),IVars-1,IVar1-1,V)
      else
        Call AMove1(InToWP(LIVar1),IVars-1,IVar1-1,V)
        Call AMove1(LRVar1,IRVars-1,IRVar1-1,V)
        Call AMove1(InToWP(MxSet1),NCrdSt-1,NCrdS1-1,V)
        Call AMove1(InToWP(NCrd1),ICrdSt-1,ICrdS1-1,V)
        endIf
      If(IPrint.ge.4.or.DEBUG) Write(IOut,2000) InfCrd,IVars,IRVars,
     $  NCrdSt,ICrdSt,IV,InfCr1,IVar1,IRVar1,NCrdS1,ICrdS1,IV1
      IVars = IVar1
      IRVars = IRVar1
      NCrdSt = NCrdS1
      ICrdSt = ICrdS1
      IV = IV1
      OK = .True.
      Return
      End
*Deck CrStOp
      Subroutine CrStOp(IOut,IPrint,IOp,IWhat,IVal,RVal,LVal)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to set and retrive options for the coordinate
C     definitions routines (mostly CrdAdd and its family of routines).
C     Dummy argument IOp indicates whether option flags are being sent from
C     the calling program, returned to the calling program, or if some
C     other action is being taken.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     IPrint      The standard print level flag.
C     IOp         This input integer argument indicates what we are doing.
C                 Allowed options include:
C                       1.  Write a value to an options array.
C                       2.  Read a value from an options array and return
C                           it to the calling program unit.
C     IWhat       This input integer argument indicates which options array
C                 and what position in that array to use.
C                       1.  Use the integer array.
C                       2.  Use the real array.
C                       3.  Use the logical array.
C                      Ix.  Use position I in the appropriate array.
C     IVal        The integer value that the calling program sends in or
C                 receives back, depending on IOp. If IWhat indicates that
C                 the array of interest is the real or logical one, this
C                 argument is ignored.
C     RVal        The real value that the calling program sends in or
C                 receives back, depending on IOp. If IWhat indicates that
C                 the array of interest is the integer or logical one, this
C                 argument is ignored.
C     LVal        The logical value that the calling program sends in or
C                 receives back, depending on IOp. If IWhat indicates that
C                 the array of interest is the integer or real one, this
C                 argument is ignored.
C
C     The elements of the Integer Options array are:
C        1-10. Unused.
C
C     The elements of the Real Options array are:
C        1-10. Unused.
C
C     The elements of the Logical Options array are:
C           1. DoFix  - Whether (True) or not (False) Routine CrDfFx should
C                       actually "fix" coordinate definitions.
C           2. UsSmDq - Whether (True) or not (False) available "smart"
C                       codes should be used for derivative definitions.
C        3-10. Unused.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009-2010.
C
      Parameter (LICrOps=10,LRCrOps=10,LLCrOps=10)
      Logical LCrOps,LVal
      Dimension ICrOps(LICrOps),RCrOps(LRCrOps),LCrOps(LLCrOps)
      Save ICrOps,RCrOps,LCrOps
      Data ICrOps/LICrOps*0/,RCrOps/LRCrOps*0.d0/,
     $  LCrOps/LLCrOps*.True./
 1000 Format(1x,'Enter CrStOp.')
 9000 Format(1x,'CrStOp: IWhat=',I10,' IArTyp=',I2,' IPos=',I10)
 9100 Format(1x,'CrStOp: IOp=',I10)
C
      If(IPrint.ge.5) Write(IOut,1000)
C
C     Process IWhat.
C
      IArTyp = Mod(IWhat,10)
      IPos = IWhat/10
      If(IWhat.le.0) then
        Write(IOut,9000) IWhat,IArTyp,IPos
        Call GauErr('CrStOp: IWhat <= 0.')
        endIf
      If(IArTyp.eq.0.or.IArTyp.gt.3) then
        Write(IOut,9000) IWhat,IArTyp,IPos
        Call GauErr('CrStOp: Invalid IArTyp.')
        endIf
      If(IPos.le.0) then
        Write(IOut,9000) IWhat,IArTyp,IPos
        Call GauErr('CrStOp: IPos <= 0.')
        endIf
      If((IArTyp.eq.1.and.IPos.gt.LICrOps).or.
     $  (IArTyp.eq.2.and.IPos.gt.LRCrOps).or.
     $  (IArTyp.eq.3.and.IPos.gt.LLCrOps)) then
        Write(IOut,9000) IWhat,IArTyp,IPos
        Call GauErr('CrStOp: IPos is too large.')
        endIf
C
C     Following IOp, do the requested work.
C
      If(IOp.eq.1) then
        If(IArTyp.eq.1) then
          ICrOps(IPos) = IVal
        else if(IArTyp.eq.2) then
          RCrOps(IPos) = RVal
        else if(IArTyp.eq.3) then
          LCrOps(IPos) = LVal
          endIf
      else if(IOp.eq.2) then
        If(IArTyp.eq.1) then
          IVal = ICrOps(IPos)
        else if(IArTyp.eq.2) then
          RVal = RCrOps(IPos)
        else if(IArTyp.eq.3) then
          LVal = LCrOps(IPos)
          endIf
      else
        Write(IOut,9100) IOp
        Call GauErr('CrStOp: Invalid IOp.')
        endIf
      Return
      End
*Deck CrValu
      Function CrValu(IOut,NPrim,IPrim,Cart,Coord1,Coord2,IDrive,InfCrd,
     $  IVars,RVars)
      Implicit Real*8(A-H,O-Z)
C
C     Given the internal data structure and a driver array, evaluate the
C     value for coordinate definition IPrim.
C
C     ARGUMENTS:
C     IOut        Output file unit number.
C     NPrim       The total number of primitive coordinate definitions.
C     IPrim       The index of the coordinate whose value we want
C                 evaluated.
C     Cart        The Cartesian coordinate array dimensioned (3,*).
C     Coord1      The current primary list of internal coordinate values.
C     Coord2      The current secondary list of internal coordinate values.
C     IDrive      The driver array built by routine CrDriv. It is
C                 dimensioned (NPrim,2).
C
C     InfCrd, IVars, and RVars are the standard coordinate definition
C     arrays. See Routine CrdDef for a discussion on the coordinate
C     definition data structure.
C
C     -J. L. Sonnenberg & H. P. Hratchian, 2009.
C
#include "crdparams.inc"
      Parameter (NInfCr=DEFMAXCOORDINFO)
      Dimension TmpVec(3,2),Cart(3,*),Coord1(*),Coord2(*),
     $  IDrive(NPrim,2),InfCrd(NInfCr,*),IVars(*),RVars(*)
      Save Small,One,Two
      Data Small/1.d-12/,One/1.d0/,Two/2.d0/
 9000 Format(1x,'CrValu: IPrim = ',I10)
 9010 Format(1x,'CrValu: IType  = ',I10)
C
C     Get the coordinate type and then evaluate the value. See Routine
C     CrdDef for a complete list of available types and their corresponding
C     IType values.
C
      IType = InfCrd(InfTyp,IPrim)
      IIVPos = InfCrd(InfIVr,IPrim)
      IRVPos = InfCrd(InfRVr,IPrim)
      If(IType.eq.0) then
        IAtom1 = IVars(IIVPos)
        Ixyz = -IVars(IIVPos+1)
        CrValu = Cart(Ixyz,IAtom1)
      else if(IType.eq.1) then
        IDep1 = IDrive(IVars(IIVPos),2)
        CrValu = CrGtVl(IDep1,Coord1,Coord2)
      else if(IType.eq.2) then
        IAtom1 = IVars(IIVPos+1)
        IAtom2 = IVars(IIVPos+2)
        IAtom3 = IVars(IIVPos+3)
        Call AngTst(IAtom1,IAtom2,IAtom3,Cart,CrValu)
      else if(IType.eq.3) then
        IAtom1 = IVars(IIVPos+1)
        IAtom2 = IVars(IIVPos+2)
        IAtom3 = IVars(IIVPos+3)
        IAtom4 = IVars(IIVPos+4)
        CrValu = Dihed(Cart(1,IAtom1),Cart(1,IAtom2),Cart(1,IAtom3),
     $    Cart(1,IAtom4))
      else if(IType.eq.50) then
        IDep1 = IDrive(IVars(IIVPos),2)
        IDep2 = IDrive(IVars(IIVPos+1),2)
        CrValu = CrGtVl(IDep1,Coord1,Coord2)+CrGtVl(IDep2,Coord1,Coord2)
      else if(IType.eq.51) then
        IDep1 = IDrive(IVars(IIVPos),2)
        IDep2 = IDrive(IVars(IIVPos+1),2)
        CrValu = CrGtVl(IDep1,Coord1,Coord2)-CrGtVl(IDep2,Coord1,Coord2)
      else if(IType.eq.52) then
        IDep1 = IDrive(IVars(IIVPos),2)
        IDep2 = IDrive(IVars(IIVPos+1),2)
        CrValu = CrGtVl(IDep1,Coord1,Coord2)*CrGtVl(IDep2,Coord1,Coord2)
      else if(IType.eq.53) then
        IDep1 = IDrive(IVars(IIVPos),2)
        IDep2 = IDrive(IVars(IIVPos+1),2)
        If(CrGtVl(IDep2,Coord1,Coord2).lt.Small) then
          Write(IOut,9000) IPrim
          Call GauErr('CrValu: Division by zero detected!')
          endIf
        CrValu = CrGtVl(IDep1,Coord1,Coord2)/CrGtVl(IDep2,Coord1,Coord2)
      else if(IType.eq.54) then
        IDep1 = IDrive(IVars(IIVPos),2)
        CrValu = RVars(IRVPos)*CrGtVl(IDep1,Coord1,Coord2)
      else if(IType.eq.56) then
        IDep1 = IDrive(IVars(IIVPos),2)
        CrValu = CrGtVl(IDep1,Coord1,Coord2)
      else if(IType.eq.57) then
        IDep1 = IDrive(IVars(IIVPos),2)
        IDep2 = IDrive(IVars(IIVPos+1),2)
        CrValu = CrGtVl(IDep1,Coord1,Coord2)
      else if(IType.eq.58) then
        CrValu = RVars(IRVPos)
      else if(IType.eq.60) then
        IDep1 = IDrive(IVars(IIVPos),2)
        IExp = IVars(IIVPos+1)
        CrValu = CrGtVl(IDep1,Coord1,Coord2)**GFloat(IExp)
      else if(IType.eq.61) then
        IDep1 = IDrive(IVars(IIVPos),2)
        If(CrGtVl(IDep1,Coord1,Coord2).lt.Small) then
          Write(IOut,9000) IPrim
          Call GauErr('CrValu: Division by zero detected!')
          endIf
        CrValu = One/CrGtVl(IDep1,Coord1,Coord2)
      else if(IType.eq.62) then
        IDep1 = IDrive(IVars(IIVPos),2)
        If(CrGtVl(IDep1,Coord1,Coord2).lt.Small) then
          Write(IOut,9000) IPrim
          Call GauErr('CrValu: SQRT of negative value detected!')
          endIf
        CrValu = SQRT(CrGtVl(IDep1,Coord1,Coord2))
      else if(IType.eq.63) then
        IDep1 = IDrive(IVars(IIVPos),2)
        CrValu = CrGtVl(IDep1,Coord1,Coord2)**2
      else if(IType.eq.70) then
        IDep1 = IDrive(IVars(IIVPos),2)
        CrValu = Cos(CrGtVl(IDep1,Coord1,Coord2))
      else if(IType.eq.71) then
        IDep1 = IDrive(IVars(IIVPos),2)
        CrValu = Sin(CrGtVl(IDep1,Coord1,Coord2))
      else if(IType.eq.72) then
        IDep1 = IDrive(IVars(IIVPos),2)
        CrValu = Tan(CrGtVl(IDep1,Coord1,Coord2))
      else if(IType.eq.73) then
        IDep1 = IDrive(IVars(IIVPos),2)
        CrValu = ACos(CrGtVl(IDep1,Coord1,Coord2))
      else if(IType.eq.74) then
        IDep1 = IDrive(IVars(IIVPos),2)
        CrValu = ASin(CrGtVl(IDep1,Coord1,Coord2))
      else if(IType.eq.75) then
        IDep1 = IDrive(IVars(IIVPos),2)
        CrValu = ATan(CrGtVl(IDep1,Coord1,Coord2))
      else if(IType.eq.100) then
        IDep1 = IVars(IIVPos)
        IDep2 = IVars(IIVPos+1)
        IDep3 = IVars(IIVPos+2)
        IDep4 = IVars(IIVPos+3)
        Call ASub(3,Cart(1,IDep1),Cart(1,IDep2),TmpVec)
        Call ASub(3,Cart(1,IDep3),Cart(1,IDep4),TmpVec(1,2))
        CrValu = SProd(3,TmpVec,TmpVec(1,2))
      else if(IType.eq.101) then
        IDep1 = IVars(IIVPos)
        IDep2 = IVars(IIVPos+1)
        Call ASub(3,Cart(1,IDep1),Cart(1,IDep2),TmpVec)
        CrValu = SProd(3,TmpVec,TmpVec)
      else if(IType.eq.102) then
        IDep1 = IVars(IIVPos)
        IDep2 = IVars(IIVPos+1)
        IDep3 = IVars(IIVPos+2)
        Call ASub(3,Cart(1,IDep1),Cart(1,IDep2),TmpVec)
        Call ASub(3,Cart(1,IDep3),Cart(1,IDep2),TmpVec(1,2))
        CrValu = SProd(3,TmpVec,TmpVec(1,2))
      else if(IType.eq.110) then
        Ixyz = -IVars(IIVPos)
        IAtom1 = IVars(IIVPos+1)
        IAtom2 = IVars(IIVPos+2)
        CrValu = Cart(Ixyz,IAtom1)-Cart(Ixyz,IAtom2)
      else if(IType.eq.111) then
        Ixyz = -IVars(IIVPos)
        IAtom1 = IVars(IIVPos+1)
        IAtom2 = IVars(IIVPos+2)
        IAtom3 = IVars(IIVPos+3)
        CrValu = Two*Cart(Ixyz,IAtom1)-Cart(Ixyz,IAtom2)-Cart(Ixyz,
     $    IAtom3)
      else if(IType.eq.500) then
        IDep1 = IDrive(IVars(IIVPos),2)
        CrValu = CrGtVl(IDep1,Coord1,Coord2)
      else
        Write(IOut,9000) IPrim
        Write(IOut,9010) IType
        Call GauErr('CrValu: Unknown IType.')
        endIf
      Return
      End
*Deck DelAng
      Subroutine DelAng(IA,JA,KA,LX,NTBond,NTAng,NTDih,NTRed,IZRed,
     $  Values,IntVec,FPVec)
      Implicit Real*8(A-H,O-Z)
C
C     Delete regular angles (LX=0), linear bends (LX<0) or dihedrals
C     (LX>0) involving IA-JA-KA.
C
      Dimension IZRed(4,*), Values(*), IntVec(*), FPVec(*)
C
      If(LX.gt.0) then
        ISt = NTBond + NTAng + 1
        IEnd = NTBond + NTAng + NTDih
        Do 10 I = IEnd, ISt, -1
          If((IA.eq.IZRed(1,I).and.JA.eq.IZRed(2,I).and.
     $      KA.eq.IZRed(3,I)).or.(IA.eq.IZRed(3,I).and.JA.eq.IZRed(2,I)
     $      .and.KA.eq.IZRed(1,I)).or.(IA.eq.IZRed(2,I).and.
     $      JA.eq.IZRed(3,I).and.KA.eq.IZRed(4,I)).or.(IA.eq.IZRed(4,I)
     $      .and.JA.eq.IZRed(3,I).and.KA.eq.IZRed(2,I))) then
            NTDih = NTDih - 1
            Call DelIZ1(I,NTRed,IZRed,Values,IntVec,FPVec)
            endIf
   10     Continue
      else if(LX.eq.0) then
        ISt = NTBond + 1
        IEnd = NTBond + NTAng
        Do 20 I = IEnd, ISt, -1
          If(JA.eq.IZRed(2,I).and.IZRed(4,I).eq.0.and.((IA.eq.IZRed(1,I)
     $      .and.KA.eq.IZRed(3,I)).or.(IA.eq.IZRed(3,I).and.
     $      KA.eq.IZRed(1,I)))) then
            NTAng = NTAng - 1
            Call DelIZ1(I,NTRed,IZRed,Values,IntVec,FPVec)
            endIf
   20     Continue
      else
        ISt = NTBond + 1
        IEnd = NTBond + NTAng
        Do 30 I = IEnd, ISt, -1
          If(JA.eq.IZRed(2,I).and.IZRed(4,I).lt.0.and.((IA.eq.IZRed(1,I)
     $      .and.KA.eq.IZRed(3,I)).or.(IA.eq.IZRed(3,I).and.
     $      KA.eq.IZRed(1,I)))) then
            NTAng = NTAng - 1
            Call DelIZ1(I,NTRed,IZRed,Values,IntVec,FPVec)
            endIf
   30     Continue
        endIf
      Return
      End
*Deck DelBnd
      Subroutine DelBnd(I,J,NBond,MaxCon,IBond,RBType)
      Implicit Real*8(A-H,O-Z)
C
C     Remove the bond between atoms I and J from the connectivity arrays.
C
      Dimension NBond(*), IBond(MaxCon,*), RBType(MaxCon,*)
C
      L = 0
      Do 10 K = 1, NBond(I)
        If(IBond(K,I).ne.J) then
          L = L + 1
          IBond(L,I) = IBond(K,I)
          RBType(L,I) = RBType(K,I)
          endIf
   10   Continue
      NBond(I) = L
      Return
      End
*Deck DelIZ1
      Subroutine DelIZ1(IRem,NTRed,IZRed,Values,IntVec,FPVec)
      Implicit Real*8(A-H,O-Z)
C
C     Shift elements of IZRed, etc.
C
      Dimension IZRed(4,*), Values(*), IntVec(*), FPVec(*)
C
      Do 100 I = (IRem+1), NTRed
        IZRed(1,I-1) = IZRed(1,I)
        IZRed(2,I-1) = IZRed(2,I)
        IZRed(3,I-1) = IZRed(3,I)
        IZRed(4,I-1) = IZRed(4,I)
        Values(I-1) = Values(I)
        IntVec(I-1) = IntVec(I)
  100   FPVec(I-1) = FPVec(I)
      NTRed = NTRed - 1
      Return
      End
*Deck DelIZR
      Subroutine DelIZR(IOut,IA,JA,KA,LA,IZRed,Values,IntVec,FPVec,
     $  NTBond,NTAng,NTDih,NTRed)
      Implicit Real*8(A-H,O-Z)
C
C     Delete user defined redundant coordinate from IZRed
C
      Dimension IZRed(4,*),Values(*),IntVec(*),FPVec(*)
 1000 Format(' Removing R(',I5,',',I5,')')
 1010 Format(' Removing A(',I5,',',I5,',',I5,')')
 1020 Format(' Removing D(',I5,',',I5,',',I5,',',I5,')')
 1030 Format(' Removing L(',I5,',',I5,',',I5,',',I5,',',I5,')')
 1040 Format(' Removing L(',I5,',',I5,',',I5,',',I5,')')
C
C     Check bonds
C
      Call FndCor(IA,JA,KA,LA,NTBond,NTAng,NTDih,IZRed,IRed,IType)
      ITypeA = IAbs(IType)
      If(ITypeA.eq.1) then
        Write(IOut,1000) IA, JA
        If(IType.gt.0) NTBond = NTBond - 1
      else if(ITypeA.eq.2) then
        If(LA.lt.-9) then
          IA4A = -LA/10
          IA4T = Mod(LA,10)
          Write(IOut,1030) IA, JA, KA, IA4A, IA4T
        else if(LA.lt.0) then
          Write(IOut,1040) IA, JA, KA, LA
        else
          Write(IOut,1010) IA, JA, KA
          endIf
        If(IType.gt.0) NTAng = NTAng - 1
      else if(ITypeA.eq.3) then
        Write(IOut,1020) IA, JA, KA, LA
        If(IType.gt.0) NTDih = NTDih - 1
      else
        Call GauErr('DelIZR is confused about IType.')
        endIf
C
C     Remove the coordinate from the list
C
      If(IType.gt.0) Call DelIZ1(IRed,NTRed,IZRed,Values,IntVec,FPVec)
      Return
      End
*Deck DelRed
      Subroutine DelRed(IRed,NTRed,IZRed,Values,IntVec,FPVec)
      Implicit Real*8(A-H,O-Z)
C
C     Delete the redundant internal coordinate number IRed
C
      Dimension IZRed(4,*), Values(*), IntVec(*), FPVec(*)
C
      Do 10 I = (IRed+1), NTRed
        IZRed(1,I-1) = IZRed(1,I)
        IZRed(2,I-1) = IZRed(2,I)
        IZRed(3,I-1) = IZRed(3,I)
        IZRed(4,I-1) = IZRed(4,I)
        Values(I-1) = Values(I)
        IntVec(I-1) = IntVec(I)
   10   FPVec(I-1) = FPVec(I)
      NTRed = NTRed - 1
      Return
      End
*Deck DoCon
      Subroutine DoCon(IOut,IPrint,IDoCon,IOConI,IOCon,MaxCon,NAtoms,
     $  IAn,C,NMapZA,MapZAt,ToAng,MxBond,NBond,IBond,IBType,RBType,V,
     $  MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Read or generate connectivity information and store it on IOCon if
C     IOCon is non-zero.
C
C     Input arguments:
C     IDoCon ... 1 - Do nothing.
C                2 - Read connectivity from input stream.
C                3 - Generate connectivity here based on distances.
C                4 - Read from chk file
C                5 - Read from rwf file.
C               10 - Read modifications from input stream.
C     IOConI ... File from which to read connectivity.
C     IOCon  ... If non-zero, file onto which to save connectivity.
C     MaxCon ... Maximum connectivity allocated for in IBond, etc.
C     NMapZA ... Highest allowed value indexing MapZAt.
C     MapZAt ... Maps numbers in the connectivity input to atom numbers.
C                Typically used to handle input in terms of z-matrix
C                numbering, since this routine returns the connectivty
C                over real atoms only.  Should have identity map for
C                normal input over real atoms only.
C
C     Output:
C     MxBond ... Largest number of bonds to an atom.  IBond, and RBType
C                are returned packed (MxBond,NAtoms).
C     NBond  ... Number of bonds to each atom.
C     IBond  ... List of atoms bonded to each atom.
C     RBType ... Bond order of each bond.
C
      Logical OKLab
      Dimension IAn(*), C(3,*), MapZAt(*), NBond(*), IBond(MaxCon,*),
     $  IBType(MaxCon,*), RBType(MaxCon,*), V(*), JJ(1)
      Save OnePt5, TolA, RCC, RCN, RCO, RNN, AmidBO, JJ
      Data OnePt5/1.5d0/, TolA/0.03d0/, RCC/1.4d0/, RCN/1.34d0/,
     $  RCO/1.36d0/, RNN/1.34d0/, AmidBO/1.41d0/, JJ/0/
      Dist(I,J) = Sqrt((C(1,I)-C(1,J))**2 + (C(2,I)-C(2,J))**2+
     $                 (C(3,I)-C(3,J))**2)
 1000 Format(' Recover connectivity data from disk.')
C
      Call IClear(NAtoms,NBond)
      Call IClear(NAtoms*MaxCon,IBond)
      Call AClear(NAtoms*MaxCon,RBType)
      IDCon0 = Mod(IDoCon,10)
      IDCon1 = Mod(IDoCon,100)/10
      If(IDCon0.eq.3) then
        IP2Box = 1
        IndTab = IP2Box + InToWP(4*NAtoms)
        IndBox = IndTab + InToWP(2*NAtoms)
        IBox = IndBox + InToWP(NAtoms)
        IBoxSt = IBox + InToWP(NAtoms+1)
        IEnd = IBoxSt + InToWP(NAtoms+1) - 1
        Call TstCor(IEnd,MDV,'DoCon-GenCBx')
        Call GenCBx(IOut,IPrint,2,ToAng,MaxCon,NAtoms,IAn,C,NBond,IBond,
     $    IBType,V(IP2Box),V(IndTab),V(IndBox),V(IBox),V(IBoxSt))
        Do 10 I = 1, NAtoms
          IA = IAn(I)
          Do 10 IB = 1, NBond(I)
            J = IBond(IB,I)
            RIJ = Dist(J,I)
            RIJA = ToAng * RIJ
            JA = IAn(J)
            K = IsAmid(ToAng,I,J,NAtoms,MaxCon,IAn,NBond,IBond,C) +
     $          IsAmid(ToAng,J,I,NAtoms,MaxCon,IAn,NBond,IBond,C)
            If(K.gt.0) then
              RBType(IB,I) = AmidBO
            else if(IA.eq.6.and.JA.eq.6.and.Abs(RIJA-RCC).lt.TolA) then
C             resonant C-C bond
              RBType(IB,I) = OnePt5
            else if(IA.eq.7.and.JA.eq.7.and.Abs(RIJA-RNN).lt.TolA) then
C             resonant N-N bond
              RBType(IB,I) = OnePt5
            else if(((IA.eq.6.and.JA.eq.7).or.(IA.eq.7.and.JA.eq.6))
     $        .and.Abs(RIJA-RCN).lt.TolA) then
C             resonant C-N bond
              RBType(IB,I) = OnePt5
            else if(((IA.eq.6.and.JA.eq.8).or.(IA.eq.8.and.JA.eq.6))
     $        .and.Abs(RIJA-RCO).lt.TolA) then
C             resonant C-O bond
              RBType(IB,I) = OnePt5
            else
              RBType(IB,I) = IBType(IB,I)
              endIf
   10       Continue
      else if(IDCon0.eq.4.or.IDCon0.eq.5) then
        Write(IOut,1000)
        Call RWCon(2,IOConI,NAtoms,MaxCon,MxBond,NBond,IBond,RBType)
        endIf
      If(IDCon0.eq.2) then
        Do 30 IZ = 1, NMapZA
          Call FFRead(IEOF)
          If(IEOF.ne.0) Call GauErr('End of file reading connectivity.')
          OKLab = .True.
   20     Call Rd1Con(IOut,OKLab,IZ,NAtoms,NMapZA,MapZAt,MaxCon,NBond,
     $      IBond,RBType,JZ)
          OKLab = .False.
          If(JZ.gt.0) goto 20
   30     Continue
        Call FFRead(IEOF)
      else if(IDCon1.eq.1) then
   40   Call FFRead(IEOF)
        If(IEOF.eq.0) then
          IZ = IFFGet(.True.,-1,.False.)
          If(IZ.gt.0) then
   50       Call Rd1Con(IOut,.False.,IZ,NAtoms,NMapZA,MapZAt,MaxCon,
     $        NBond,IBond,RBType,JZ)
            If(JZ.gt.0) goto 50
            Goto 40
            endIf
          endIf
        endIf
      Call SymCon(IOut,IPrint,NAtoms,MaxCon,NBond,IBond,RBType)
      IIBI = 1
      IRBTI = IIBI + InToWP(MaxCon)
      IInd = IRBTI + MaxCon
      IEnd = IInd + InToWP(MaxCon) - 1
      Call TstCor(IEnd,MDV,'DoCon-SrtCon')
      Call SrtCon(.True.,MaxCon,NAtoms,NBond,IBond,JJ,RBType,V(IIBI),JJ,
     $  V(IRBTI),V(IInd))
      If(IOCon.ne.0)
     $  Call RWCon(1,IOCon,NAtoms,MaxCon,MxBond,NBond,IBond,RBType)
      If(IPrint.ge.2) Call PrtCon(IOut,'Connectivity in DoCon:',.True.,
     $  NAtoms,MaxCon,IAn,NBond,IBond,IBond,RBType,0)
      Return
      End
*Deck ElDum
      Subroutine ElDum(MaxNZ,IOut,IPrint,NAtoms,Keep,IAn,C,IZ,Map,
     $    HavCon,MxBond,NBond,IBond)
      Implicit Real*8(A-H,O-Z)
C
C     Eliminate the atoms not requested in Keep and update NAtoms.
C
      Logical HavCon, Keep(1)
      Dimension IAn(1), C(3,1), Map(1), IZ(MaxNZ,1), NBond(1),
     $    IBond(MxBond,1)
 1000 Format(' Non-terminal atom number',I3,' deleted.')
 1010 Format(' Atom number',I3,' deleted.')
C
      M = 0
      Call IClear(NAtoms,Map)
      Do 100 I = 1, NAtoms
        If(.not.Keep(I).and.IPrint.ge.2) Write(IOut,1010) I
        If(Keep(I)) then
          M = M + 1
          Map(I) = M
          IAn(M) = IAn(I)
          Call AMove(3,C(1,I),C(1,M))
          If(M.eq.1) then
            IZ(M,1) = 1
          else if(IZ(I,1).eq.0) then
            IZ(M,1) = 0
          else
            IZ(M,1) = Map(IZ(I,1))
            endIf
          If(IZ(M,1).eq.0.and.IZ(I,1).ne.0) then
            Write(IOut,1000) IZ(I,1)
            Call Lnk1E(0)
            endIf
          endIf
  100   Continue
      If(.not.HavCon) goto 900
      M1 = 0
      Do 200 I = 1, NAtoms
        If(Map(I).eq.0) goto 200
        M1 = M1 + 1
        NB = 0
        LimJ = NBond(I)
        Do 150 J = 1, LimJ
          If(Map(IBond(J,I)).eq.0) goto 150
          NB = NB + 1
          IBond(NB,M1) = Map(IBond(J,I))
  150     Continue
        NBond(M1) = NB
  200   Continue
  900 NAtoms = M
      Return
      End
*Deck FixACM
      Subroutine FixACM(ICpLCh,NAtoms,IAnMod,AtChMM,ACMod)
      Implicit Real*8(A-H,O-Z)
C
C     Fix ACMod entries for link atoms, possibly copying from AtChMM.
C
      Dimension IAnMod(*), AtChMM(*), ACMod(*)
      Save Zero
      Data Zero/0.0d0/
C
      If(ICpLCh.ne.1.and.ICpLCh.ne.3) Return
      Do 10 I = 1, NAtoms
        If(IAnMod(I).ne.0.and.(ICpLCh.eq.3.or.ACMod(I).eq.Zero))
     $    ACMod(I) = AtChMM(I)
   10   Continue
      Return
      End
*Deck FixCou
      Subroutine FixCou(IOut,NFrag,NAtoms,IAtFrg,NDHvLay,HavLay,NZ,
     $  IAnZ,NuOfLa,Layer)
      Implicit Real*8(A-H,O-Z)
C
C     Check that fragment data is consistent.
C
      Logical HavLay(NDHvLay)
      Dimension IAtFrg(*), IAnZ(*), Layer(*)
 1000 Format(' Atoms in',I2,' fragments were given but there should be',
     $  I2,'.')
 1010 Format(' Atom',I7,' is in fragment',I7,' but NFrag=',I7,'.')
C
      If(NFrag.gt.1) then
        If(NuOfLa.eq.1) then
          NAt = 0
          Do 10 I = 1, NZ
            If(IAnZ(I).gt.-1) then
              NAt = NAt + 1
              If(IAtFrg(NAt).eq.0.and.Layer(I).gt.0)
     $          IAtFrg(NAt) = Layer(I)
              endIf
   10       Continue
          Call IClear(NZ,Layer)
          endIf
        Call LClear(NDHvLay,HavLay)
        Do 50 I = 1, NAtoms
          If(IAtFrg(I).lt.1.or.IAtFrg(I).gt.NFrag) then
            Write(IOut,1010) I, IAtFrg(I), NFrag
            Call GauErr('Inconsistent fragment input.')
            endIf
   50     HavLay(IAtFrg(I)) = .True.
        NHave = 0
        Do 60 L = 1, NDHvLay
          If(HavLay(L)) NHave = NHave + 1
   60     Continue
        If(NHave.ne.NFrag) then
          Write(IOut,1000) NHave, NFrag
          Call Lnk1E(0)
          endIf
        endIf
      Return
      End
*Deck FixMc1
      Subroutine FixMc1(IOut,IPrint,NAtoms,Lab,IR,IA,MicOpt)
      Implicit Integer(A-Z)
C
C     Check that one atom IA from internal coordinate IR is flagged as
C     in the "QM" block in MicOpt.
C
      Logical Error
      Character*(*) Lab
      Dimension MicOpt(*)
C
      Call IniStr
      Call BldStC('Coordinate ')
      Call BldStC(Lab)
      Call BldStI(0,IR)
      Call BldStC(' contains atom')
      Call BldStI(1,IA)
      IA1 = Mod(IA-1,NAtoms) + 1
      OMOIA = MicOpt(IA1)
      Error = .False.
      If(MicOpt(IA1).eq.0) then
        Call BldStC(' moved to QM block.')
        MicOpt(IA1) = -1
      else if(MicOpt(IA1).eq.-1) then
        Call BldStC(' which is already in QM block.')
      else if(MicOpt(IA1).eq.-2.or.MicOpt(IA1).eq.-3) then
        Call BldStC(' which is frozen.')
      else
        IBlk = MicOpt(IA1)
        If(IBlk.lt.-40) IBlk = -IBlk - 40
        Call BldStC(' but this is in rigid fragment')
        Call BldStI(1,IBlk)
        Error = .True.
        endIf
      If(IPrint.gt.2.or.(OMOIA.ne.-1.and.OMOIA.ne.-2)) Call DmpStr(IOut)
      If(Error) Call GauErr('Inconsistent input found by FixMc1.')
      Return
      End
*Deck FixMcO
      Subroutine FixMcO(IOut,IPrint,NAtoms,NTBond,NTAng,NTDih,IZRed,
     $  MicOpt)
      Implicit Integer(A-Z)
C
C     Mark atoms included in internal coordinates as being in the "QM"
C     region.  Some might be marked for microiterations but have been
C     mentioned in ModRedund input.
C
      Dimension IZRed(4,*), MicOpt(*)
      Character*(1) Lab
C
      Do 10 I = 1, NTBond
        If(IZRed(2,I).gt.0) then
          Call FixMc1(IOut,IPrint,NAtoms,'R',I,IZRed(1,I),MicOpt)
          Call FixMc1(IOut,IPrint,NAtoms,'R',I,IZRed(2,I),MicOpt)
        else
          Call FixMc1(IOut,IPrint,NAtoms,'C',I,IZRed(1,I),MicOpt)
          endIf
   10   Continue
      Do 20 I = 1, NTAng
        I1 = I + NTBond
        Call RedL45(IZRed(4,I1),ITypeL,LL1,LL2)
        If(ITypeL.eq.0) then
          Lab = 'A'
        else
          Lab = 'L'
          endIf
        Call FixMc1(IOut,IPrint,NAtoms,Lab,I,IZRed(1,I1),MicOpt)
        Call FixMc1(IOut,IPrint,NAtoms,Lab,I,IZRed(2,I1),MicOpt)
        Call FixMc1(IOut,IPrint,NAtoms,Lab,I,IZRed(2,I1),MicOpt)
        If(ITypeL.eq.2) Call FixMc1(IOut,IPrint,NAtoms,Lab,I,LL1,MicOpt)
   20   Continue
      Do 30 I = 1, NTDih
        I1 = I + NTBond + NTAng
        Call FixMc1(IOut,IPrint,NAtoms,'D',I,IZRed(1,I1),MicOpt)
        Call FixMc1(IOut,IPrint,NAtoms,'D',I,IZRed(2,I1),MicOpt)
        Call FixMc1(IOut,IPrint,NAtoms,'D',I,IZRed(3,I1),MicOpt)
        Call FixMc1(IOut,IPrint,NAtoms,'D',I,IZRed(4,I1),MicOpt)
   30   Continue
      Return
      End
*Deck FixOni
      Subroutine FixOni(IOut,ChkDVL,MaxLay,NuOfLa,NZ,HavLay,IAnMod,IAnZ,
     $  LiAtom,Layer,IAtMod,ACMod,NDist,Distan,MxBond,NBond,IBond,
     $  RBType,IMap)
      Implicit Real*8(A-H,O-Z)
C
C     Check that dummy atoms are not substituted by model atoms or linked
C     to model atoms, that the layers are labelled consistently, etc.
C     If NuOfLa is 0 on input, it is set to the correct value; otherwise
C     it is check against the input.
C
      Logical ChkDVL, HavLay(MaxLay), Failed, OKDist
      Character*1 NamLI, NamLJ, LayNam
      Real*8 MDCutO
      Dimension IAnMod(*), IAnZ(*), LiAtom(*), Layer(*), IAtMod(*),
     $  ACMod(*), Distan(NDist,*), NBond(*), IBond(MxBond,*), IMap(*),
     $  RBType(MxBond,*)
      Save Zero, One, BTMin
      Data Zero/0.0d0/, One/1.0d0/, BTMin/0.4d0/
 1000 Format(' Atoms in',I2,' layers were given but there should be',I2,
     $  '.')
 1010 Format(' Atom',I6,' layer ',A,' linked to atom',I6,
     $  ' is also bonded to atom',I6,' in layer ',A)
 1020 Format('            OK because scale factors are 1.')
 1030 Format('            bad because scale factors are defaulted or',
     $  ' not 1.')
 1040 Format(' Atoms are bonded to more than one atom in an inner',
     $  ' layer and have non-unit',/,' scale factors.  ',
     $  ' Correct the layer boundaries, connectivity, or',
     $  ' scale factors.')
C
      Call LClear(MaxLay,HavLay)
      Do 10 I = 1, NZ
        If(IAnMod(I).ne.0) then
          If(IAnZ(I).eq.-1)
     $      Call GauErr('Dummy atom cannot become a model atom!')
          LiI = LiAtom(I)
          If(LiI.ne.0) then
            If(IAnZ(LiI).eq.-1) Call GauErr(
     $        'Model atom cannot be connected to dummy atom!')
            endIf
          endIf
   10   Continue
      NDummy = 0
      Do 20 I = 1, NZ
        If(IAnZ(I).eq.-1) NDummy = NDummy + 1
   20   IMap(I) = NDummy
      Do 30 J = 1, NZ
        LJ = LiAtom(J)
        If(LJ.ne.0) LiAtom(J) = LJ - IMap(LJ)
   30   Continue
      NAt = 0
      Do 40 I = 1, NZ
        If(IAnZ(I).ne.-1) then
          NAt = NAt + 1
          IMap(NAt) = I
          Layer(NAt) = Max(Layer(I),1)
          IAnMod(NAt) = IAnMod(I)
          IAtMod(NAt) = IAtMod(I)
          ACMod(NAt) = ACMod(I)
          LiAtom(NAt) = LiAtom(I)
          Do 35 II = 1, MaxLay
   35       Distan(NAt,II) = Distan(I,II)
          If(Layer(NAt).gt.MaxLay)
     $      Call GauErr('Layer/Fragment number out of range in FixOni.')
          HavLay(Layer(NAt)) = .True.
          endif
   40   Continue
      NHave = 0
      Do 60 L = 1, MaxLay
        If(HavLay(L)) then
          NHave = NHave + 1
          Do 50 I = 1, NAt
            If(Layer(I).eq.L) Layer(I) = NHave
   50       Continue
          endIf
   60   Continue
      If(NuOfLa.eq.0) then
        NuOfLa = NHave
      else if(NHave.ne.NuOfLa.and.NuOfLa.gt.1) then
        Write(IOut,1000) NHave, NuOfLa
        Call Lnk1E(0)
        endIf
      Failed = .False.
      Small = MDCutO(0)
      Do 110 I = 1, NAt
        LiI = LiAtom(I)
        LaI = Layer(I)
        NamLI = LayNam(LaI,NuOfLa)
        If(LiI.eq.0.and.NuOfLa.gt.1) then
          Do 70 J0 = 1, NBond(I)
            J = IBond(J0,I)
            LaJ = Layer(J)
            If(LaJ.lt.LaI.and.RBType(J0,I).ge.BTMin) then
              NamLJ = LayNam(LaJ,NuOfLa)
              Call IniStr
              Call BldStC('Atom ')
              Call BldAtN(IAnZ(IMap(I)),I)
              Call BldStC(' in layer ')
              Call BldStC(NamLI)
              Call BldStC(' is bonded to ')
              Call BldAtN(IAnZ(IMap(J)),J)
              Call BldStC(' in layer ')
              Call BldStC(NamLJ)
              If(IAnMod(I).eq.0) then
                IAnMod(I) = 1
                Call BldStC(' and was made a link atom replaced by H.')
              else
                Call BldStC(' and was linked to it.')
                endIf
              Call DmpStr(IOut)
              LiAtom(I) = J
              LiI = J
              Goto 80
              endIf
   70       Continue
   80     endIf
        If(LiI.gt.0) then
          If(IAnZ(IMap(I)).eq.IAnMod(I)) then
            Call IniStr
            Call BldStC('Warning:  atom ')
            Call BldAtN(IAnZ(IMap(I)),I)
            Call BldStC(' in layer ')
            Call BldStC(NamLI)
            Call BldStC(' is uses the same element as its link atom.')
            Call DmpStr(IOut)
            endIf
          OKDist = Abs(Distan(I,1)-One).lt.Small
          Do 90 L = 2, NuOfLa
   90       OKDist = OKDist.and.(Distan(I,L).eq.Zero.or.
     $        Abs(Distan(I,L)-One).lt.Small)
          Do 100 IB = 1, NBond(I)
            J = IBond(IB,I)
            LaJ = Layer(J)
            If(J.ne.LiI.and.LaJ.lt.LaI.and.RBType(IB,I).ge.BTMin) then
              NamLJ = LayNam(LaJ,NuOfLA)
              Write(IOut,1010) I, NamLI, LiI, J, NamLJ
              If(OKDist) then
                Write(IOut,1020)
              else
                Failed = .True.
                Write(IOut,1030)
                endIf
              endIf
  100       Continue
          endIf
  110   Continue
      If(Failed.and.ChkDVL) then
        Write(IOut,1040)
        Call GauErr('Terminating because of bad ONIOM partitioning.')
        endIf
      Return
      End
*Deck FndCor
      Subroutine FndCor(IA,JA,KA,LA,NTBond,NTAng,NTDih,IZRed,IRed,IType)
      Implicit Integer(A-Z)
C
C     Find a coordinate in IZRed and return its location or where it
C     should be inserted and its type.  The Type is negative if the
C     coordinate was not found.
C
      Dimension IZRed(4,*)
C
      IRed = 0
      If(KA.eq.0.and.LA.eq.0) then
C       Cartesian coordinates and bonds
        Do 10 I = 1, NTBond
          If((IA.eq.IZRed(1,I).and.JA.eq.IZRed(2,I)).or.
     $      (IA.eq.IZRed(2,I).and.JA.eq.IZRed(1,I))) IRed = I
   10     Continue
        If(IRed.eq.0) then
          IType = -1
          IRed = NTBond + 1
        else
          IType = 1
          endIf
      else if(LA.le.0) then
C       Regular angles and linear bends
        IType = 2
        Do 20 I = (NTBond+1), (NTBond+NTAng)
          If(JA.eq.IZRed(2,I).and.LA.eq.IZRed(4,I).and.
     $      ((IA.eq.IZRed(1,I).and.KA.eq.IZRed(3,I)).or.
     $      (LA.ge.-9.and.IA.eq.IZRed(3,I).and.KA.eq.IZRed(1,I))))
     $      IRed = I
   20     Continue
        If(IRed.eq.0) then
          IType = -2
          IRed = NTBond + NTAng + 1
        else
          IType = 2
          endIf
      else
C       Dihedrals
        Do 40 I = (NTBond+NTAng+1), (NTBond+NTAng+NTDih)
          If((IA.eq.IZRed(1,I).and.JA.eq.IZRed(2,I).and.KA.eq.IZRed(3,I)
     $      .and.LA.eq.IZRed(4,I)).or.(IA.eq.IZRed(4,I).and.
     $      JA.eq.IZRed(3,I).and.KA.eq.IZRed(2,I).and.LA.eq.IZRed(1,I)))
     $      IRed = I
   40     Continue
        If(IRed.eq.0) then
          IType = -3
          IRed = NTBond + NTAng + NTDih + 1
        else
          IType = 3
          endIf
        endIf
      Return
      End
*Deck GenITR
      Subroutine GenITR(In,IOut,IPrint,SetMcO,KeepMM,DoRdFz,CkRigF,
     $  IStruc,MaxNZ,NZ,IZ,IAnZ,NLay,Layer,LiAtom,IAnMod,NAtoms,IAn,
     $  IAtTyp,IAtRes,NResid,IResNa,IResNu,ITRead,MicOpt,OptAtm,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Generate the ITRead and MicOpt flag arrays for RedCor.
C
C     ITRead:
C      0 : Regular optimizer atom
C      2 : Frozen atom
C      3 : Micro-iterations atom
C
C     MicOpt:
C  N<-40 : Atom in Fragment -N-41, but currently frozen (corresponding
C          to -N-40 if not frozen).
C    -2N : Flag for added charge 0<N<4 indicates for which system
C    -1N : Flag for added center, N giving the gridpoint.
C     -3 : Frozen atom in lowest level of ONIOM.  Will be changed to 1
C          in a local copy in L103 and L120 if microiterations are done
C          for ONIOM(MO:MM)
C     -2 : Frozen atom in high level (or non-oniom input)
C     -1 : Atom treated with redundant internals.
C      0 : Optimize with micro-iterations.
C    N>1 : Belongs to rigid fragment N-1, to be optimized with micro-iterations.
C
      Logical SetMcO, KeepMM, DoRdFz, CkRigF, OptAtm(*), Exist, MicQM
      Dimension IZ(MaxNZ,*), IAnZ(*), Layer(*), LiAtom(*), IAnMod(*),
     $  IAn(*), IAtTyp(*), IAtRes(*), IResNa(*), IResNu(*), ITRead(*),
     $  MicOpt(*), V(MDV)
C
      If(SetMcO) then
C       Set MicOpt from IZ
        IA = 0
        Call ISet(-1,NZ,MicOpt)
        Do 10 I = 1, NZ
          If(IAnZ(I).ge.0) then
            IA = IA + 1
            If(Layer(IA).ge.NLay.and.LiAtom(IA).eq.0.and.
     $        IAnMod(IA).eq.0.and..not.KeepMM) then
              If(IZ(I,1).eq.-1) then
                MicOpt(IA) = -3
              else if(IZ(I,1).lt.-1) then
                MicOpt(IA) = -IZ(I,1)
              else
                MicOpt(IA) = 0
                endIf
              ITRead(IA) = 3
            else if(IZ(I,1).eq.-1) then
              MicOpt(IA) = -2
              ITRead(IA) = 2
            else if(IZ(I,1).lt.-1) then
              MicOpt(IA) = -IZ(I,1)
              If(CkRigF) Call GauErr(
     $          'Only 1 rigid fragment with regular optimizer')
              endIf
            endIf
   10     Continue
C       Compress the numbers for the rigid fragments.
        NMax = IArMax(.False.,NZ,MicOpt)
        IShift = 0
        Do 30 I = 2, NMax
          Exist = .False.
          Do 20 J = 1, NZ
            If(MicOpt(J).eq.I) then
              Exist = .True.
              MicOpt(J) = MicOpt(J) - IShift
              endif
   20       Continue
          If(.not.Exist) IShift = IShift + 1
   30     Continue
        endIf
C
C       Possibly read modifications for freezing and optimizing.
C
      If(DoRdFz) then
        If(NAtoms.ne.IA) Call GauErr('Input inconsistency to GenITr.')
        If(IStruc.eq.1) then
          IAtLst = 1
          IBlk = IAtLst + InToWP(NAtoms)
          INoBlk = IBlk + InToWP(NAtoms)
          Call TstCor(INoBlk+InToWP(NAtoms)-1,MDV,'GenITr-RdFzAt')
          Call RdFzAt(In,IOut,IPrint,NAtoms,MicOpt,IAn,IAtTyp,IAtRes,
     $      NResid,IResNa,IResNu,Layer,LiAtom,OptAtm,V(IAtLst),
     $      V(IBlk),V(INoBlk))
          endIf
        Call ApFzAt(NAtoms,OptAtm,ITRead,MicOpt)
        endIf
      Do 200 I = 1, NAtoms
        If(.not.MicQM(MicOpt(I))) then
          ITRead(I) = 3
        else if(MicOpt(I).eq.-1) then
          ITRead(I) = 0
        else
          ITRead(I) = 2
          endIf
  200   Continue
C     Set IZ from MicOpt
      IA = 0
      Do 300 I = 1, NZ
        If(IAnZ(I).ge.0) then
          IA = IA + 1
          If(MicOpt(IA).eq.-1.or.MicOpt(IA).eq.0) then
            If(IZ(I,1).eq.-1) IZ(I,1) = 0
          else if(MicOpt(IA).eq.-2.or.MicOpt(IA).eq.1.or.
     $      MicOpt(IA).eq.-3) then
            If(IZ(I,1).eq.0) IZ(I,1) = -1
          else if(IZ(I,1).le.0) then
            IZ(I,1) = MOp2IZ(MicOpt(IA))
            endIf
          endIf
  300   Continue
      If(IPrint.ge.1)
     $  Call PrtITM(IOut,'After GenITR',NAtoms,ITRead,MicOpt)
      Return
      End
*Deck GeoNam
      Subroutine GeoNam(IOut,IGoI,IRead,IGo,ICkUnt,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Given IGoI, possibly read a structure source from the input file,
C     set IGo, and handle reading the title and charge/multiplicity from
C     a chk file.
C
      Dimension V(*)
      Character*(*) ZMtSt1,ZMtSt2, RdStr, ChkSt1, ChkSt2, ModASt,
     $  ModBSt, Label
      Parameter (ZMtSt1='Z-Matrix',ZMtSt2='ZMatrix',RdStr='Read',
     $  ChkSt1='Check',ChkSt2='Checkpoint',ModASt='ModelA',
     $  ModBSt='ModelB',Label='Structure')
      Common /MUnit/ IUnit(20)
      Character*1024 NamFil, NamTst
      Logical StEqCC
C
      ICkUnt = IUnit(9)
      If(IRead.eq.0) then
        IGo = IGoI
        If(IGo.eq.1) IGo = 0
      else
        Call FFGtNm(-1,.True.,NamFil,LNam,NamTst,LTst)
        If(StEqCC(NamFil(1:LNam),ZMtSt1).or.
     $    StEqCC(NamFil(1:LNam),ZMtSt2)) then
          IGo = 0
        else if(StEqCC(NamFil(1:LNam),RdStr).or.
     $    StEqCC(NamFil(1:LNam),ChkSt1).or.
     $    StEqCC(NamFil(1:LNam),ChkSt2)) then
          IGo = 2
        else if(StEqCC(NamFil(1:LNam),ModASt)) then
          IGo = 4
        else if(StEqCC(NamFil(1:LNam),ModBSt)) then
          IGo = 5
        else
          IGo = 2
          ICkUnt = IUnit(5)
          endIf
        endIf
      If(IGo.eq.2.or.IGo.eq.3.or.IGo.eq.6.or.IGo.eq.7)
     $  Call DoChkN(IOut,ICkUnt,Label,NamFil(1:LNam),V,MDV)
      Return
      End
*Deck GTrans
      Subroutine GTrans(IOut,IPrint,NStepI,IAprBG,NAtoms,NTrRot,CP,CR,
     $  CX,XValue,IZRed,NTRed,NTBond,NTAng,NTDih,WeiFac,MStruc,Rstr,
     $  PStr,DQ,DQI,XStr,B,IB,G,Square,Tri,One4D,CnStr,IFlags,XNames,V,
     $  MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Guess a transition structure.
C
      Parameter (NStepD=10)
      Logical Just1
      Dimension B(3,4,NTRed), IB(4,NTRed), G(NTRed,NTRed),
     $  Square(NTRed,NTRed), Tri(*), One4d(NTRed,4), CR(3,NAtoms),
     $  CP(3,NAtoms), CX(3,NAtoms), IZRed(4,NTRed), RStr(NTRed),
     $  PStr(NTRed), XStr(3,NAtoms), XValue(NTRed), Cnstr(NTRed),
     $  IFlags(NTRed), XX(1), XNames(*), DQ(*), DQI(*), V(MDV)
      Save XX, One, Conv
      Data XX/0.d0/, One/1.0d0/, Conv/1.d-4/
C
C     Interpolate in internal coordinates, then find the corresponding
C     Cartesians, or as close as possible.
C
      NStep = NStepI
      If(NStep.eq.0) NStep = NStepD
      IPrt1 = IPrint
      If(IPSave(0).le.1) IPrt1 = Max(IPrint,1)
      NAt3 = 3*NAtoms
C
C     Find the TS stepping in small increments.
C
      IOldRU = 3
      Just1 = .True.
      IGISav = 0
      IGIOld = 0
      IGoal = 1
      IV = IGoal + NTRed
      Call TstCor(IV,MDV,'GTrans')
      MDV1 = MDV - IV + 1
      Call AMove(NAt3,CR,CX)
      WR = WeiFac
      WP = One - WeiFac
      Call ACSASB(NTRed,PStr,WP,RStr,WR,V(IGoal))
      Call FixDih(NTDih,V(NTBond+NTAng+IGoal))
      Call AMove(NTRed,RStr,XValue)
      Do 100 IPass = 1, NStep
        Call ASub(NTRed,V(IGoal),XValue,DQ)
        Call FixDq(IOut,IPrint,1,NTRed,NTBond,NTAng,NTDih,DQ,XX)
        DQMax = ArrMax(DQ,NTRed,.True.)
        If(DQMax.lt.Conv.and.IPass.gt.1) goto 110
        Call AMove(NTRed,DQ,DQI)
        Call RedCar(IOut,IPrt1,0,NTRed,NTBond,NTang,NTdih,IZRed,XValue,
     $    NAtoms,CX,DQI,XStr,B,IB,Square,Tri,One4D,Cnstr,IFlags,.False.,
     $    IFail,0,V(IV),MDV1)
        If(IFail.ne.0) then
          Call AMove(NTRed,DQ,DQI)
          Call ORedCr(IOut,IPrt1,IOldRU,Just1,IAprBG,IGISav,IGIOld,
     $      NTrRot,NTRed,NTBond,NTAng,NTDih,IZRed,XValue,DXMaxT,NAtoms,
     $      CX,DQI,XStr,B,IB,G,Square,Tri,One4D,Cnstr,IFlags,XNames,
     $      IFail,V(IV),MDV1)
          endIf
        If(IFail.ne.0) Call GauErr('RedCar/ORedCr failed for GTrans.')
        Call AMove(NAt3,XStr,CX)
        Call AMove(NTRed,XValue,DQI)
        Call RedVal(IOut,IPrint,.True.,CX,IZRed,XValue,NTRed,NTBond,
     $    NTAng,NTDih,V(IV),MDV1)
        If(IPrint.ge.2) Call PrntDQ(IOut,IPass,NTRed,DQI,XValue,DQ)
  100   Continue
  110 If(IPrint.ge.0) then
        Call ASub(NTRed,PStr,RStr,DQ)
        Call AScale(NTRed,One-WeiFac,DQ,DQ)
        Call PrntDQ(IOut,0,NTRed,RStr,XValue,DQ)
        endIf
      Return
      End
*Deck IIVbDv
      Subroutine IIVbDv(IOp,IOff,II,LenII,IOut,IPrint,IHess)
      Implicit Real*8(A-H,O-Z)
C
C     This routine is used to fill the integer array of flags (II) for the
C     frequency-by-Davidson diagonalization routines.  IOp is passed in as
C     1 when all of the other arguments are to be used to fill II.  IOp is
C     passed in as 2 when II should be decomposed and all of the other
C     arguments returned accordingly.  IOp is passed in as 3 to only get
C     the value of LenII.  IOff is an offset in the pointer values.  For
C     IOp=3, IOff is ignored.
C
C     The elements of II are:
C           II(1)       LenII       The length of the flags array.  This
C                                   element is NOT affected by IOff.
C           II(2)       IOut        The output file unit number.  This
C                                   element is NOT affected by IOff.
C           II(3)       IPrint      The standard print level flag.  This
C                                   element is NOT affected by IOff.
C           II(4)       IHess       The position of the Hessian matrix
C                                   stored in the real vector space.
C
C     - H. P. Hratchian, 2009.
C
      Parameter (LLnII=4)
      Dimension II(LLnII)
C
      If(IOp.eq.1) then
        II(1)  = LLnII
        II(2)  = IOut
        II(3)  = IPrint
        II(4)  = IHess-IOff
      else if(IOp.eq.2) then
        If(II(1).ne.LLnII)
     $    Call GauErr('InVbDv: Problem with LenII!')
        LenII  = II(1)
        IOut   = II(2)
        IPrint = II(3)
        IHess  = II(4)+IOff
      else if(IOp.eq.3) then
        LenII    = LLnII
      else
        Call GauErr('InVbDv: Unknown IOp value!')
        endIf
      Return
      End
*Deck InpFrq
      Subroutine InpFrq(IOut,AllChk,ForcW0,SrtFrq,IOFDPr,MaxFrq,NFrqRd,
     $  FreqRd,Ind,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Read in the frequencies for frequency-dependent perturbations and
C     initialize the file holding FD properties, if IOFDPr is non-zero.
C     The frequencies are left in V(3)..V(NFrqRd+2).  The frequencies are
C     sorted into increasing order (if SrtFrq is true) and if ForcW0 is
C     set then W=0 is added if necessary.
C
      Parameter (IndFr=3)
      Common /MUnit/ IUnit(20)
      Integer FilNum
      Logical AllChk, ForcW0, SrtFrq
      Dimension FreqRd(MaxFrq), Ind(MaxFrq), V(MDV)
      Save Zero
      Data Zero/0.0d0/
 1020 Format(' Using perturbation frequencies:',5F12.6)
C
      If(MDV.lt.3) Call GauErr('Silly MDV in InpFrq.')
      Call TstCor(MaxFrq,MDV,'InpFrq-Read')
      If(AllChk) then
        IOFDPC = FilNum(IOFDPr,IUnit(9))
        If(ITqry(IOFDPC).lt.2)
     $    Call GauErr('No perturbation frequencies on chk file.')
        Call FileIO(2,-IOFDPC,2,V,0)
        Call BNUPak(V,Junk,NFrqRd)
        Call TstCor(NFrqRd,MDV,'InpFrq-Read')
        Call FileIO(2,IOFDPC,NFrqRd,V,0)
      else
        Call FFGetR(.True.,.True.,.True.,0,IEOF,MaxFrq-1,V,NFrqRd)
        If(NFrqRd.lt.1)
     $    Call GauErr('CPHF=RdFreq but no frequencies in input.')
        Call AAbs(NFrqRd,V,V)
        endIf
      If(SrtFrq) then
        Call AHpSrt(NFrqRd,Ind,V)
      else
        Call ISeq(NFrqRd,0,1,Ind)
        endIf
      If(ForcW0.and.V(Ind(1)).ne.Zero) then
        I0 = 1
        FreqRd(I0) = Zero
      else
        I0 = 0
        endIf
      Do 10 I = 1, NFrqRd
   10   FreqRd(I+I0) = V(Ind(I))
      NFrqRd = NFrqRd + I0
      Lst = IndFr + NFrqRd - 1
      Call TstCor(Lst,MDV,'InpFrq-Store')
      Call BNPack(0,NFrqRd,V)
      Call BNPack(1,1,V(2))
      Call AMove(NFrqRd,FreqRd,V(3))
      If(IOFDPr.ne.0) Call ConWrt(IOFDPr,Lst,V)
      Write(IOut,1020) (V(I),I=IndFr,Lst)
      Return
      End
*Deck InsRed
      Subroutine InsRed(IRed,IA,JA,KA,LA,V,IntV,FPV,NTRed,IZRed,Values,
     $  IntVec,FPVec)
      Implicit Real*8(A-H,O-Z)
C
C     Insert a new redundant internal coordinate, to be number IRed.
C
      Dimension IZRed(4,*), Values(*), IntVec(*), FPVec(*)
C
      Do 10 I = NTRed, IRed, -1
        IZRed(1,I+1) = IZRed(1,I)
        IZRed(2,I+1) = IZRed(2,I)
        IZRed(3,I+1) = IZRed(3,I)
        IZRed(4,I+1) = IZRed(4,I)
        Values(I+1) = Values(I)
        IntVec(I+1) = IntVec(I)
   10   FPVec(I+1) = FPVec(I)
      IZRed(1,IRed) = IA
      IZRed(2,IRed) = JA
      IZRed(3,IRed) = KA
      IZRed(4,IRed) = LA
      Values(IRed) = V
      IntVec(IRed) = IntV
      FPVec(IRed) = FPV
      NTRed = NTRed + 1
      Return
      End
*Deck IntStp
      Subroutine IntStp(IOut,IPrint,JStep,IUnit,NVar,Values,MinXX,ICore)
      Implicit Real*8(A-H,O-Z)
C
C     Fetch the variable values from intermediate point JStep.
C     of an optimization.
C
      Integer FilNum
      Dimension Values(*), ICore(*)
      Logical Prnt,Exit,HavDEP,Done,Large
      Common /OptGrd/ ConvF,FMaxT,DXMaxT,RMax,RMin,RLim,EigMax,EigMin,
     $  DEPred,FSwtch,FncErr,GrdErr,FncCnv,XXIRC,IPOCOU,NStep,IStep,
     $  NMax,NP,Neg,Prnt,Exit,IPSav1,IPSav2,HavDEP,ICheck,IRisen,MaxRis,
     $  Large,ILevel,ISStep,IOffs,IMDS,MaxSv,MaxSv2,NVrMMT
      Save IOptGr
      Data IOptGr/575/
 1000 Format(' Optimization data not found on unit',I2,'.')
 1010 Format(' Consistency failure: NVar=',I4,' MaxSv=',I4,' MaxSv2=',
     $  I4,' LOptGr=',I3,' LOptG2=',I8,' LenOG=',I8,'.')
 1020 Format(' Requested step has not been saved.')
 1030 Format(' Available points:',10I5)
 1040 Format(' Step',I4,' is the next step -- use Opt=Restart instead.')
 1050 Format(' Step',I4,' is the last structure,',
     $  ' equivalent to Geom=CheckPoint.')
C
      IOptGr = FilNum(IOptGr,IUnit)
      LenOG = ITqry(IOptGr)
      If(LenOG.lt.1) then
        Write(IOut,1000) IUnit
        Call Lnk1E(0)
        endIf
      LOptGr = 14 + InToWP(22)
      Call FileIO(2,-IOptGr,LOptGr,ConvF,0)
      Call All103(.False.,NVar,NVrMMT,MinXX,MaxSv,MaxSv2,1,IIC,IX,
     $  IXName,IESave,IFS,IFC1,IFC2,IXX,IFF,IEQMS,IXXMM,IFFMM,IFFMMQ,
     $  IITU,LOptG2)
      If(LenOG.lt.(LOptGr+LOptG2-InToWP(MaxSv))) then
        Write(IOut,1010) NVar, MaxSv, MaxSv2, LOptGr, LOptG2, LenOG
        Call Lnk1E(0)
        endIf
C
C     Now recover the variables.
C
      If(JStep.eq.(IStep+1)) then
        Write(IOut,1040) JStep
        Call Lnk1E(0)
      else if(JStep.eq.IStep) then
        Write(IOut,1050) JStep
      else
        Ind = LOptGr + IIC - 1
        Call FileIO(2,-IOptGr,InToWP(MaxSv),ICore,Ind)
        JStep2 = JStep
        If(JStep2.eq.-1) JStep2 = 0
        JStep1 = 0
        ICMax = 0
        Do 10 I = 1, MaxSv
          ICore(I) = ICore(I) - 1
          If(ICore(I).eq.JStep2.and.JStep1.eq.0) JStep1 = I
          Done = ICore(I).lt.0
          If(I.ge.2) Done = Done.or.(ICore(I).eq.ICore(I-1))
          If(Done) goto 20
   10     ICMax = I
   20   If(JStep1.eq.0) Write(IOut,1020)
        If(IPrint.gt.0.or.JStep1.eq.0)
     $    Write(IOut,1030) (ICore(I),I=1,ICMax)
        If(JStep1.eq.0) Call Lnk1E(0)
        Ind = LOptGr + IXX - 1 + (JStep1-1)*NVar
        Call FileIO(2,-IOptGr,NVar,Values,Ind)
        endIf
      Return
      End
*Deck LoadV
      Subroutine LoadV(NAtoms,MaxNZ,BlUsed,AUsed,BUsed,NumBl,NumA,NumB,
     $    Bl,Al,Be,ConstA,ConstB,ANames,Values,NVar,LBl,LAlpha,LBeta,
     $    IScr)
      Implicit Real*8(A-H,O-Z)
C
C     Load the generated variables specified by use flags xUsed,
C     pointers Numx (see Redund for more details), whose values
C     are already in Bl, Al, and Be into the symbolic Z-matrix arrays.
C
      Logical BlUsed(1), AUsed(1), BUsed(1), ConstA(1), ConstB(1)
      Dimension NumBl(1), NumA(1), NumB(1), Bl(1), Al(1), Be(1),
     $          ANames(1), Values(1), LBl(1), LAlpha(1), LBeta(1),
     $          IScr(1)
C
      NVar = 0
      ICur = 0
      Call CBlank(IntPWP(0)*3*MaxNZ,ANames)
      Call IClear(NAtoms,IScr)
      Do 10 I = 1, NAtoms
          If(.not.BlUsed(I)) goto 10
          NVar = NVar + 1
          Call PutCSt('R',ANames,ICur)
          Call DecChr(I,ANames,ICur)
          Call PutDel(2,ANames,ICur)
          Values(NVar) = Bl(I)
          IScr(I) = NVar
   10     LBl(I) = IScr(NumBl(I))
      Call IClear(NAtoms,IScr)
      Do 20 I = 1, NAtoms
          If(.not.AUsed(I)) LAlpha(I) = IScr(NumA(I))
          If(ConstA(NumA(I))) LAlpha(I) = 0
          If(.not.AUsed(I).or.ConstA(NumA(I))) goto 20
          NVar = NVar + 1
          Call PutCSt('A',ANames,ICur)
          Call DecChr(I,ANames,ICur)
          Call PutDel(2,ANames,ICur)
          Values(NVar) = Al(I)
          IScr(I) = NVar
          LAlpha(I) = NVar
   20     Continue
      Call IClear(NAtoms,IScr)
      Do 30 I = 1, NAtoms
          If(.not.BUsed(I).and.NumB(I).gt.0)
     $        LBeta(I) = IScr(NumB(I))
          If(.not.BUsed(I).and.NumB(I).lt.0)
     $        LBeta(I) = -IScr(-NumB(I))
          If(ConstB(IAbs(NumB(I)))) LBeta(I) = 0
          If(.not.BUsed(I).or.ConstB(IAbs(NumB(I)))) goto 30
          NVar = NVar + 1
          Call PutCSt('D',ANames,ICur)
          Call DecChr(I,ANames,ICur)
          Call PutDel(2,ANames,ICur)
          Values(NVar) = Be(I)
          IScr(I) = NVar
          LBeta(I) = NVar
   30     Continue
      Call PutDel(5,ANames,ICur)
      Return
      End
*Deck LocStN
      Integer Function LocStN(IDelim,N,LenFnd,StrFnd,List,StrTmp)
      Implicit Integer(A-Z)
C
      Logical StrEqC,Found
      Dimension List(*),StrTmp(*)
C
C     Do While(LenFnd.ne.LenTmp.and.i.lt.NLabel)
      LocStN = 0
      Found = .False.
      If(N.lt.1) Return
      i = 1
 100  Call GetNSt(.True.,.False.,IDelim,i,List,StrTmp,Junk)
      Found = Junk.eq.LenFnd
      If(Found) Found = StrEqC(StrTmp,StrFnd,LenFnd)
      i = i+1
      If(.not.Found.and.i.lt.N) Goto 100
      If(Found) LocStN = i-1
      Return
      End
*Deck LSubst
      Function LSubst(NamCnt,N,Str,Len)
      Implicit Integer(A-Z)
C
C     A routine to return the sequential number of the center
C     whose name is in "Str" (a hollerith string of length "Len").
C     The list of known names in the delimited string "NamCnt", and
C     there are "N" of these names.
C
      Logical StrEqC
      Common/IO/In,IOut,IPunch
      Dimension TStr(10), NamCnt(*), Str(*)
C
C     Initialize and check each name from the list.
      NCur = 0
      LSubst = 0
      Do 10 I = 1, N
        Call GetB(2,TStr,LenT,NamCnt,NCur)
        If(Len.eq.LenT.and.StrEqC(TStr,Str,Len)) then
          LSubst = I
          Return
          endIf
   10   Continue
C
C     Unknown center.
C
      Call IniStr
      Call BldStC('Unknown center ')
      Call BldStr(Str,Len)
      Call DmpStr(IOut)
      Call Lnk1E(0)
      Return
      End
*Deck MakZmt
      Subroutine MakZMT(IOut,IPrint,UseSym,MaxNZ,MaxAp3,MaxOp,NZ,ICharg,
     $  Multip,NFree,IAn,IZ,BL,AL,BE,NPerm,C,AtmChg,ToAng,Keep,LBl,
     $  LAlpha,LBeta,NVar,ANames,Values,NumBl,NumA,NumB,BLUsed,AUsed,
     $  BUsed,IScr,ConstA,ConstB,Scr,HavCon,KepOrd,MxBond,NBond,IBond,
     $  TolerD,TolerG,IBType,NToO,OToN,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     This routine takes the fixed Z-Matrix whose connectivity is
C     specified in IZ which was may contain many centers with
C     IZ(,4).ne.0 (specified by two bond angles) and builds a symbolic
C     z-matrix using only dihedral angles and incorporating all
C     available symmetry constraints.  The new Z-matrix is written
C     over the old one.  Most of the other arrays are used only for
C     scratch.  HavCon indicates whether the full connectivity
C     matrix is available in NBond and IBond.  Note that Scr
C     is used as both 8*MaxNZ and integer MaxAP3*NOp.  TolerD and
C     TolerG are tolerances used here and in the symmetry package:  Zero
C     values cause the defaults to be used.  IBType contains
C     bond types which are permuted along with NBond and IBond if
C     HavCon is set.  NToO and OToN receive the permutation applied to
C     the atoms and its inverse.  KepOrd forces the order of atoms to
C     be maintained; it effectively turns off most use of symmetry.
C
      Integer OToN(1)
      Logical HavCon, KepOrd, BLUsed(1), AUsed(1), BUsed(1), OK, Error,
     $  ConstB(1), Keep(1), ConstA(1), UseSym
      Integer MaxAtm,NSymOp,NReps,LblRep,IPrmut
      Real*8 ChrTbl,SymOps
      Parameter (MaxAtm=DEFMAXATM)
      Common/RepCom/NSymOp,NReps,LblRep(32),ChrTbl(10,16),SymOps(9,10),
     $  IPrmut(MaxAtm,10)
      Dimension IAn(*), BL(*), AL(*), C(3,*), NBond(*), BE(*), V(*),
     $  TrVec(3), IBond(MxBond,*), AtmChg(*), IScr(*), NumBl(*),
     $  NumA(*), NumB(*), IZ(MaxNZ,4), Scr(MaxNZ,9), ANames(*),
     $  Values(*), LBl(*), LAlpha(*), LBeta(*), RotMat(3,3), NGrp(4),
     $  IBType(MxBond,*), NPerm(MaxAP3,*), NToO(*)
      Save CutD, Zero, One
      Data CutD/1.d-08/, Zero/0.0d0/, One/1.0d0/
 1000 Format(' Z-Matrix conversion failed!')
 1010 Format(' Coordinates after elimination of dummy atoms:')
 1020 Format(' All atoms were deleted!')
 1030 Format(' Connectivity after elimination of dummy atoms:')
 1040 Format(' I=',I3,' NBond=',I2,' IBond=',10I3)
 1050 Format(' New connectivity before FormZ:')
 1060 Format(' I=',I3,' IZ=',4I4)
 1070 Format(' Coordinates after symmetrization:')
C
C     Determine the point group and permutations of the atoms
C     under symmetry operations.
C
      If(NZ.le.0) goto 999
      Cut = TolerD
      If(Cut.eq.Zero) Cut = CutD
      Call IClear(MaxNZ,IScr)
      Call ZToc(MaxNZ,NZ,IAn,IZ,Bl,Al,Be,.False.,NAtoms,Scr(1,1),C,
     $  Scr(1,2),Scr(1,5),Scr(1,6),Scr(1,7),Scr(1,8),Scr(1,9),IOut,
     $  Error,.False.)
      If(Error) Call Lnk1E(0)
C
C     Eliminate dummy atoms.
      Call ElDum(MaxNZ,IOut,IPrint,NAtoms,Keep,IAn,C,IZ,IScr,
     $  HavCon,MxBond,NBond,IBond)
      NZ = NAtoms
      If(NAtoms.le.0) Write(IOut,1020)
      If(NAtoms.le.0) Call Lnk1E(0)
      If(IPrint.gt.4) Write(IOut,1010)
      If(IPrint.gt.4) Call CorPr1(IOut,NAtoms,IAn,C,ToAng)
      If(IPrint.gt.4.and.HavCon) Write(IOut,1030)
      Do 12 I = 1, NAtoms
        LimJ = NBond(I)
        If(HavCon.and.IPrint.gt.4)
     $    Write(IOut,1040) I, NBond(I), (IBond(J,I),J=1,LimJ)
   12   Continue
      Do 15 I = 1, NAtoms
   15   AtmChg(I) = GFloat(IAn(I))
C
C     Determine the resulting molecular symmetry.
      If(UseSym) then
        IMolF = 1
        IFWG = IMolF + InToWP(2*NAtoms)
        ITrans = IFWG + InToWP(10*NAtoms)
        ICNew  = ITrans + Max(9,MaxAP3)*Max(MxBond,MaxOp)
        ICScr1 = ICNew + 3*MaxAP3
        ICScr2 = ICScr1 + 3*MaxAP3
        ICScr3 = ICScr2 + 3*MaxAP3
        IScr1  = ICScr3 + 3*MaxAP3
        IScr2  = IScr1 + (NAtoms*(NAtoms+1))/2
        IScr3  = IScr2 + NAtoms
        IScr4  = IScr3 + NAtoms
        IScr5  = IScr4 + NAtoms
        IScr6  = IScr5 + NAtoms
        IEnd   = IScr6 + NAtoms - 1
        Call TstCor(ITrans+IEnd-1,MDV,'MakZmt')
        IP1 = Max(IPrint-2,0)
        IDump = IP1
        Call Symm(NAtoms,ICharg,Multip,MaxAP3,MaxOp,1,IDump,IP1,ToAng,
     $    IAn,AtmChg,C,TolerD,TolerG,'    ',0,LenFor,V(IMolF),LenFWG,
     $    V(IFWG),NOp,NPerm,V(ITrans),NoSymE,TrVec,RotMat,NGrp,V(ICNew),
     $    V(ICScr1),V(ICScr2),V(ICScr3),V(IScr1),V(IScr2),V(IScr3),
     $    V(IScr4),V(IScr5),V(IScr6))
        If(NoSymE.ne.0) NOp = 0
        NFree = 0
        If(NAtoms.eq.2) NFree = 1
        If(NAtoms.eq.3) NFree = 3*NAtoms - 6
        If(NoSymE.eq.0) NFree = NumDof(V(IFWG),NAtoms,LenFWG)
        N = 0
        Call PutStr(V(IMolF),LenFor,V(IMolF),N)
        N = 0
        Call PutStr(V(IFWG),LenFWG,V(IFWG),N)
        Call IniStr
        Call BldStC('Stoichiometry ')
        Call BldStr(V(IMolF),LenFor)
        If(NoSymE.eq.0) Call BldStC(' Framework group ')
        If(NoSymE.eq.0) Call BldStr(V(IFWG),N)
        If(IPrint.ge.3) Call DmpStr(IOut)
      else
        NoSymE = 1
        NOp = 1
        Do 10 I = 1, NAtoms
   10     NPerm(I,1) = I
        endIf
C
C     Symmetrize the coordinates, in case loose symmetry cutoffs have
C     been used.  Note that NSymOp will be nonsense if there's no
C     symmetry (NOp1=1, C1 point group).
C
      IPrt1 = Min(IPrint-1,0)
      If(NoSymE.eq.0.and.NOp.gt.1) Call SymCor(IOut,IPrt1,.True.,TrVec,
     $  .False.,RotMat,NAtoms,C)
      If(IPrint.gt.4) Write(IOut,1070)
      If(IPrint.gt.4) Call CorPr1(IOut,NAtoms,IAn,C,ToAng)
C
C     Choose a new connectivity scheme.
      If(HavCon.and..not.KepOrd) then
        IScr3 = 1
        IScr4 = IScr3 + InToWP(NAtoms)
        ICScr1 = IScr4 + InToWP(NAtoms)
        ITrans = ICScr1 + 3*NAtoms
        IEnd = ITrans + Max(MaxAP3*NOp,MxBond*NAtoms) - 1
        Call TstCor(IScr3+IEnd-1,MDV,'MaxZMt-NewIZ')
        Call NewIZ(IOut,MaxNZ,MaxAP3,NAtoms,
     $    NOp,MxBond,NPerm,NBond,IBond,IZ,IAn,C,NToO,OToN,V(IScr3),
     $    V(IScr4),V(ICScr1),V(ITrans),IBType,V(ITrans))
        endIf
      Call ReBldZ(IOut,UseSym,MaxNZ,MaxAP3,Cut,NPerm,NAtoms,NOp,
     $  IAn,C,IZ)
C
C     Compute the bond distances, angles and dihedral angles
C     for the new connectivity.
C
      Call IClear(MaxNZ,LBl)
      Call IClear(MaxNZ,LAlpha)
      Call IClear(MaxNZ,LBeta)
      Call AClear(3*MaxNZ,Values)
      If(IPrint.gt.3) Write(IOut,1050)
      Do 100 I = 1, NZ
        If(IPrint.gt.3) Write(IOut,1060) I, (IZ(I,J),J=1,4)
  100   Continue
      Call FormZ(IOut,0,NAtoms,NAtoms,MaxNZ,IAn,C,IZ,Bl,Al,Be,
     $  .false.,.false.,LBl,LAlpha,LBeta,Values,Scr,ToAng,OK)
      If(.not.OK) Write(IOut,1000)
      If(.not.OK) Call Lnk1E(0)
      If(IPrint.ge.3) Call ZPrint(MaxNZ,NAtoms,IAn,IZ,Bl,Al,Be,One)
C
C     Eliminate redundant bond distances and angles, then
C     check for bond and dihedral angles constrained to their
C     present values by symmetry.
C
      Call Redund(MaxNZ,MaxAP3,NAtoms,NOp,Cut,IZ,NPerm,Be,BlUsed,
     $  AUsed,BUsed,NumBl,NumA,NumB)
      If(UseSym.and.NOp.gt.1) then
        Call CkConA(MaxNZ,NAtoms,IZ,IAn,C,Al,Cut,AUsed,NConA,ConstA)
        Call CkConB(MaxNZ,NAtoms,IZ,IAn,C,Be,Cut,BUsed,NConB,ConstB)
      else
        Call LClear(NAtoms,ConstA)
        Call LClear(NAtoms,ConstB)
        endIf
C
C     Load the information on variables.
      Call LoadV(NAtoms,MaxNZ,BlUsed,AUsed,BUsed,NumBl,NumA,NumB,Bl,Al,
     $    Be,ConstA,ConstB,ANames,Values,NVar,LBl,LAlpha,LBeta,IScr)
  999 Return
      End
*Deck MBld2
      Subroutine MBld2(IModel,IPrint,NEl,El,IUConv,M,An,IEof,Keep,FromA,
     $  EV,EVC,Ch,G)
      Implicit Real*8(A-H,O-Z)
C
C     This program builds a qualitative standard geometric model of a
C     molecule.
C
      Parameter (MaxMBA=DEFMAXATM,MaxMB1=(MaxMBA+1))
      Dimension An(*), IStrng(20), IShort(6,MaxMBA), LShort(6,MaxMBA),
     $  NamStr(20), OutLin(80)
      Equivalence (Name,NamStr(1))
      Integer FROMA(*),EV(*),EVC(*),CH(*),AN,CN,BT,G(6,*),EXSVAL,F,BDAT,
     $  BDBD,RGBD,RGAT,El(-2:*),AtBd,ATRG,BDRG,RGARG,RGBRG,VL,RGSIZE,
     $  DIH,X,Y,TESTA,TOTEV,RGTYPE,ATGEOM,BDTYPE,BDROT,TETR,TRIG,PYRA,
     $  BENT,ASYM,CIS,TRAN,GAUP,GAUM,SKEW,CONJ,SING,DOUB,TRIP,AROM,TARO,
     $  DATV,BLANK,IJunk(10),COMMAT(10),TREE(10,10),VLData(18)
      REAL*8 NURTNG, NEWANG
      Logical Debug, BiltIn, StrEqC, NumEl, Keep(1)
      COMMON/MB/MXCN,NRG,NBD,MXRGSZ,CN(MaxMB1),BT(MaxMB1),
     $  EXSVAL(MaxMB1),F(6,MaxMBA),BDAT(6,MaxMBA),BDBD(6,MaxMBA),
     $  RGBD(6,MaxMBA),RGAT(6,MaxMBA),AtBd(2,MaxMBA),ATRG(10,10),
     $  BDRG(10,10),RGARG(10,10),RGBRG(10,10),NRGAT(MaxMBA),
     $  NRGBD(MaxMBA),VL(18),RGSIZE(10),DIH(2,MaxMBA),RGTYPE(10),
     $  ATGEOM(MaxMBA),BDTYPE(MaxMBA),BDROT(MaxMBA),NXCARD(20),IMBPad,
     $  SKRTNG(MaxMBA),NEWANG(MaxMBA)
      COMMON/IO/In,IOUT,IPUNCH
      Save TETR,TRIG,PYRA,BENT,LINE,ASYM,CIS,TRAN,GAUP,GAUM,SKEW,CONJ
      Save SING,DOUB,TRIP,AROM,TARO,DATV,BLANK,IDash,Debug,IFFFP,IFFInt
      Save IFFStr,IFFEnd,Zero,One,F45,VLData
      Data Debug/.false./
      DATA TETR/1/,TRIG/2/,PYRA/3/,BENT/4/,LINE/5/,ASYM/6/,CIS/7/,
     $TRAN/8/,GAUP/9/,GAUM/10/,SKEW/11/,CONJ/12/,SING/13/,DOUB/14/,
     $TRIP/15/,AROM/16/,TARO/17/,DATV/18/,BLANK/1h /,IDash/1h-/
      DATA IFFFP/2HFP/, IFFInt/3HINT/, IFFSTR/3HSTR/, IFFEND/3HEND/
      Data Zero/0.0d0/, One/1.0d0/, F45/45.0d0/
C
C     VLData contains the normal valences of all atoms in the first
C     two rows of the periodic table.
C
      Data VLData/1,0,1,2,3,4,3,2,1,0,1,2,3,4,5,6,1,0/
 1000 Format(' Short formula matrix:')
 1005 Format(' Use Cis bond rotation for ',4I3,'.')
 1010 Format(' Use Trans bond rotation for ',4I3,'.')
 1020 Format(' Use Gauche + bond rotation for ',4I3,'.')
 1030 Format(' Use Gauche - bond rotation for ',4I3,'.')
 1040 Format(' Use Skew bond rotation for ',4I3,' angle ',F12.4,'.')
 1060 Format(3X,I2,7X,A1,A2,'         ',5X,1X,6(I2,2X))
 1070 Format(3X,I2,7X,A1,A2,'   Tetr  ',5X,1X,6(I2,2X))
 1080 Format(3X,I2,7X,A1,A2,'   Trig  ',5X,1X,6(I2,2X))
 1090 Format(3X,I2,7X,A1,A2,'   Pyra  ',5X,1X,6(I2,2X))
 1100 Format(3X,I2,7X,A1,A2,'   Bent  ',5X,1X,6(I2,2X))
 1110 Format(3X,I2,7X,A1,A2,'   Line  ',5X,1X,6(I2,2X))
 1120 Format(3X,I2,7X,A1,A2,'   ASym  ',F5.1,1X,6(I2,2X))
 1130 Format(4X,4HBond,3X,6HAtom,1,3X,6HAtom,2,3X,15HBond definition,3X,
     $9HBond type,3X,13HBond rot geom)
 1140 Format('                           Ring information:',/,
     $       '  Ring   RingType                   AtRg')
 1150 Format('  ',I2,'             ',I2,'-Ring',10I4)
 1160 Format('  ',I2,'       Arom  ',I2,'-Ring',10I4)
 1170 Format('  ',I2,'       Conj  ',I2,'-Ring',10I4)
C
C     Initialize all variables to zero.
C
      MXCN=0
      NRG=0
      NBD=0
      MXRGSZ=0
      Call IClear(MaxMBA,An)
      Call IClear(MaxMB1,CN)
      Call IClear(MaxMB1,BT)
      Call IClear(MaxMB1,EXSVAL)
      Call AClear(MaxMBA,NewAng)
      Call AClear(MaxMBA,SKRTNG)
      Call IClear(MaxMBA,NRGBD)
      Call IClear(MaxMBA,NRGAT)
      Call IClear(6*MaxMBA,F)
      Call IClear(6*MaxMBA,BDAT)
      Call IClear(6*MaxMBA,BDBD)
      Call IClear(6*MaxMBA,RGBD)
      Call IClear(6*MaxMBA,RGAT)
      Call IClear(2*MaxMBA,AtBd)
      Call IClear(2*MaxMBA,DIH)
      Call IClear(10,RgSize)
      Call IClear(10*10,RgBrg)
      Call IClear(10*10,RGARG)
      Call IClear(10*10,BDRG)
      Call IClear(10*10,ATRG)
      Call IMove(18,VLData,VL)
      Call IClear(10,RgType)
      Call IClear(MaxMBA,AtGeom)
      Call IClear(MaxMBA,BdType)
      Call IClear(MaxMBA,BdRot)
      I=0
      J=0
      K=0
      L=0
      N=0
      X=0
      Y=0
      Call CBlank(6*MaxMBA,IShort)
      Call IClear(6*MaxMBA,LShort)
      AngCon = One
      If(IUConv.gt.1) AngCon = F45 / ATan(One)
C
C     Read the short formula matrix from the input stream and convert
C     it to expanded form (with all atoms numbered explicitly) in F.
C
C     Do while (not BlankLine and not command card)
      N = 0
      Call CBlank(20,NxCard)
  180 N = N + 1
      Call CBlank(20,IStrng)
      Call FFRead(IEof)
      If(IEof.eq.1) goto 210
      Call FFGet(IStrng,LenS,Junk,FJunk,IFound)
      If(IFound.eq.IFFEnd) IEof = 1
      If(IFound.eq.IFFEnd) goto 210
      If(IFound.ne.IFFStr) Call FFErr(IFFStr,IFound)
      ICur = 0
      JCur = 0
      If(LenS.gt.4) Call MovChr(LenS,IStrng,ICur,NxCard,JCur,8)
      If(LenS.gt.4) goto 210
      If(NumEl(NEl,El,IStrng,L)) goto 185
      If(An(N+1).ge.0) goto 185
          Call IniStr
          Call BldStC('Unrecognized atomic symbol ')
          Call BldStr(IStrng,LenS)
          Call DmpStr(IOut)
          Call Lnk1E(0)
  185 An(N+1) = IAbs(L)
      Keep(N) = L.ge.0
      Do 190 J = 1, 6
          Call FFGet(IShort(J,N),LShort(J,N),L,FJunk,IFound)
          If(IFound.eq.IFFEnd) goto 180
          If(IFound.ne.IFFStr.and.IFound.ne.IFFInt)
     $        Call FFErr(IFFStr,IFound)
          If(IFound.eq.IFFInt) F(J,N) = L
  190     Continue
C
C     Now that the short F matrix has been read, expand the
C     abbreviations (atomic and functional group symbols) to form
C     the full F matrix.
C
  210 N = N - 1
      M = N
      Do 230 I = 1, N
          Do 230 J = 1, 6
              If(LShort(J,I).eq.0) goto 230
              If(.not.NumEl(NEl,El,IShort(J,I),L)) goto 220
                  M = M + 1
                  F(J,I) = M
                  F(1,M) = I
                  AN(M+1) = IAbs(L)
                  Keep(M) = L.ge.0
                  Goto 230
  220         If(BiltIn(J,I,M,IShort(J,I),F,An,Keep)) goto 230
                  Call IniStr
                  Call BldStC('Unrecognized group ')
                  Call BldStr(IShort(J,I),LShort(J,I))
                  Call DmpStr(IOut)
                  Call Lnk1E(0)
 230          Continue
      N = M
C
C     Read modification information.
C
      If(IPrint.gt.2) Call StrOut(IOut,NxCard,8,1)
C
C     Now define the connectivity CN for each atom I:
C         CN(I) = Number of atoms bonded to I.
C     MXCN = Maximum connectivity present in the molecule.
C
      MXCN=0
      Do 840 I=1,M
         CN(I+1)=0
         Do 800 J=1,6
             If(F(J,I).gt.0) CN(I+1)=CN(I+1)+1
  800        Continue
  840    MxCN=Max(MxCN,CN(I+1))
C
C     Now fill AtBd(J,I) with the Jth atom in bond I.
C
C     First find non-ring-closure bonds.
C
      CH(1)=1
      Do 850 I=1,M
          CH(I+1)=0
          Do 850 J=1,MXCN
  850         G(J,I)=F(J,I)
      AtBd(2,1)=1
      L=1
      Do 940 K=1,M
          I=AtBd(2,K)
          Do 930 J=1,MXCN
              IND1=G(J,I)
              If(CH(IND1+1).ne.0) goto 930
                  L=L+1
                  AtBd(2,L)=G(J,I)
                  AtBd(1,L)=I
                  CH(IND1+1)=1
                  Do 920 X=1,MXCN
                      If(G(X,IND1).eq.I) G(X,IND1)=0
  920                 Continue
                  G(J,I)=0
  930         Continue
  940     Continue
C
C     Next obtain ring closure bonds.
C
      Do 980 K=1,M
          Do 980 J=1,MXCN
              IND1=AtBd(2,K)
              If(G(J,IND1).le.AtBd(2,K)) goto 980
                  L=L+1
                  AtBd(1,L)=AtBd(2,K)
                  AtBd(2,L)=G(J,IND1)
  980         Continue
      Do 990 I = 1, M
          If(Debug) Write(IOut,9900) I, AtBd(1,I), AtBd(2,I)
 9900     Format(' I, AtBd: ',3I3)
  990     Continue
C
C     Now fill BDat(J,I) with the Jth bond containing atom I.
C
      NBD=0
      NBD=L
      Do 1370 K=2,NBD
          Do 1370 J=1,MXCN
              IND1=AtBd(1,K)
              IND2=AtBd(2,K)
              If(F(J,IND1).eq.AtBd(2,K)) BDAT(J,IND1)=K
              If(F(J,IND2).eq.AtBd(1,K)) BDAT(J,IND2)=K
 1370         Continue
      Do 991 I = 1, M
          If(Debug) Write(IOut,9910) I, (BDat(J,I),J=1,6)
 9910     Format(' I, BDat: ',7I3)
  991     Continue
C
C     Now fill BDBD(J,I) with the Jth bond attached to bond I.
C
      Do 1360 K=1,NBD
          X=0
          Do 1360 J=1,MXCN
              IND1=AtBd(1,K)
              IND2=AtBd(2,K)
              If(F(J,IND1).NE.0.AND.F(J,IND1).NE.AtBd(2,K)) X=X+1
              If(F(J,IND1).NE.0.AND.F(J,IND1).NE.AtBd(2,K))
     $            BDBD(X,K)=BDAT(J,IND1)
              If(F(J,IND2).NE.0.AND.F(J,IND2).NE.AtBd(1,K)) X=X+1
              If(F(J,IND2).NE.0.AND.F(J,IND2).NE.AtBd(1,K))
     $            BDBD(X,K)=BDAT(J,IND2)
 1360         Continue
      Do 992 I = 1, M
          If(Debug) Write(IOut,9920) I, (BDBD(J,I),J=1,6)
 9920     Format(' I, BDBD: ',7I3)
  992     Continue
C
C     Now fill AtRg(J,I) with the Jth atom in ring I.
C     Find the ring introduced by each ring-closure bond.
C
      Y = 0
      Do 1180 K = 1,M
          Do 1180 J = 1,MXCN
 1180         G(J,K)=0
      LIM1 = M+1
      Do 1270 K = LIM1, NBD
          Do 1270 J = 1,MXCN
              IND1 = AtBd(1,K)
              If(F(J,IND1).eq.AtBd(2,K)) IND1 = AtBd(1,K)
              If(F(J,IND1).eq.AtBd(2,K)) G(J,IND1) = F(J,IND1)
              IND1 = AtBd(2,K)
              If(F(J,IND1).eq.AtBd(1,K)) IND1 = AtBd(2,K)
              If(F(J,IND1).eq.AtBd(1,K)) G(J,IND1) = F(J,IND1)
 1270         Continue
      CH(1) = 1
      Do 1280 I=1,M
          Do 1280 J=1,MXCN
 1280         G(J,I) = F(J,I) - G(J,I)
      LIM1=M+1
      Do 1610 KK=LIM1,NBD
        TestA = 0
        Do 1290 I = 1,M
            CH(I+1)=0
 1290       FromA(I) = 0
        Do 1300 I = 1,10
          Do 1300 J = 1,10
 1300       Tree(I,J) = 0
        Tree(1,1) = AtBd(1,KK)
        Ind1 = AtBd(1,KK)
        CH(Ind1+1) = 1
        Do 1500 I = 1,10
          LL = 0
          Do 1470 K = 1,10
            If(Tree(I,K).le.0) goto 1470
              Do 1460 J =1,MXCN
                Ind1 = Tree(I,K)
                Ind2 = G(J,Ind1)
                If(CH(Ind2+1).ne.0) goto 1460
                  LL = LL+1
                  Ind1 = Tree(I,K)
                  Tree(I+1,LL) = G(J,Ind1)
                  Ind1 = Tree(I+1,LL)
                  FromA(Ind1) = Tree(I,K)
                  CH(Ind1+1) = 1
                  If(Tree(I+1,LL).ne.AtBd(2,KK)) goto 1450
                    TestA = 1
                    Y = Y+1
                    RGSIZE(Y) = I+1
                    LIM2 = I+2
                    Do 1440 NN = 1,LIM2
                      If(NN.EQ.1.OR.NN.EQ.I+2) AtRg(NN,Y) = AtBd(2,KK)
                      If(NN.ne.1.and.NN.ne.I+2) Ind1 = ATRG(NN-1,Y)
                      If(NN.ne.1.and.NN.ne.I+2) ATRG(NN,Y) = FromA(Ind1)
 1440                 Continue
                    Ind1 = Tree(I+1,LL)
                    CH(Ind1+1) = 0
 1450             Continue
 1460           Continue
 1470       Continue
          If(TestA.eq.1) goto 1510
 1500     Continue
 1510   Do 1600 J = 1,MXCN
          Ind1 = AtBd(1,KK)
          Ind2 = AtBd(2,KK)
          If(F(J,Ind1).eq.Ind2) G(J,Ind1) = Ind2
          If(F(J,Ind2).eq.Ind1) G(J,Ind2) = Ind1
 1600     Continue
 1610   Continue
      Call IClear(10,CommAt)
      LL = Y
      If(LL.le.1) goto 1820
      Do 1810 K = 1,LL
          LIM1 = K+1
          Do 1810 KK = LIM1,LL
              MM = 0
              LIM2 = RGSIZE(K)
              LIM3 = RGSIZE(KK)
              Do 1700 J = 1,LIM2
                  Do 1690 JJ=1,LIM3
                      If(ATRG(J,K).ne.ATRG(JJ,KK)) goto 1690
                          MM = MM + 1
                          COMMAT(MM) = ATRG(J,K)
                          Goto 1700
 1690                 Continue
 1700             Continue
              If(MM.le.2) goto 1800
                  Y = Y + 1
                  NN = 1
                  ATRG(NN,Y) = COMMAT(1)
                  LIM2 = RGSIZE(K)
                  Do 1760 J = 1,LIM2
                      Do 1750 L = 1,MM
                          If(ATRG(J,K).eq.COMMAT(L)) goto 1760
 1750                     Continue
                      NN = NN + 1
                      ATRG(NN,Y) = ATRG(J,K)
 1760                 Continue
                  LIM3 = RGSIZE(KK)
                  Do 1790 JJ = 1,LIM3
                      Do 1780 L = 1,MM
                          If(ATRG(JJ,KK).eq.COMMAT(L)) goto 1790
 1780                     Continue
                      NN = NN+1
                      ATRG(NN,Y) = ATRG(JJ,KK)
 1790                 Continue
                  NN = NN + 1
                  ATRG(NN,Y) = COMMAT(MM)
                  RGSIZE(Y) = NN
                  ATRG(NN+1,Y) = ATRG(1,Y)
 1800             Continue
 1810         Continue
 1820 NRG = Y
      MM = 0
      Do 1870 K = 1,NRG
 1870     MM = Max(RgSize(K),MM)
      MXRGSZ = MM
      Do 1980 K = 1,NRG
          LIM1 = NBD+1
          Do 1880 I = 1,LIM1
 1880         CH(I) = 0
          CH(1) = 1
          CH(2) = 1
          X = 0
          LIM1 = RGSIZE(K)
          Do 1970 I = 1,LIM1
              Do 1960 J = 1,MXCN
                  Ind1 = ATRG(I,K)
                  If(F(J,Ind1).eq.0) goto 1960
                      LIM2 = I+1
                      Do 1950 KK = LIM2,LIM1
                          Ind2 = BDAT(J,Ind1)
                          If(ATRG(KK,K).ne.F(J,Ind1).or.
     $                        CH(Ind2+1).ne.0) Goto 1950
                              X = X + 1
                              Ind1 = ATRG(I,K)
                              BDRG(X,K) = BDAT(J,Ind1)
                              Ind2 = BDAT(J,Ind1)
                              CH(Ind2+1) = 1
 1950                     Continue
 1960             Continue
 1970         Continue
 1980     Continue
      Do 993 I = 1, 10
          If(Debug) Write(IOut,9930) I, (BDRG(J,I),J=1,10)
 9930     Format(' I, BDBD: ',11I3)
  993     Continue
C
C     Now fill RGAT(J,I) with the Jth ring in which atom I participates.
C
      Do 2030 I = 1,M
          X = 0
          Do 2020 K = 1,NRG
              LIM1 = RGSIZE(K)
              Do 2015 J=1,LIM1
                  If(I.ne.ATRG(J,K)) goto 2015
                      X = X+1
                      RGAT(X,I) = K
 2015             Continue
 2020         Continue
 2030     NRGAT(I) = X
      Do 994 I = 1, M
          If(Debug) Write(IOut,9940) I, (RgAt(J,I),J=1,6)
 9940     Format(' I, RgAt: ',7I3)
  994     Continue
C
C     Now fill RGBD(J,I) with the Jth ring in which bond I participates.
C
      Do 2080 I = 2,NBD
          X = 0
          Do 2070 K = 1,NRG
              LIM1 = RGSIZE(K)
              Do 2065 J=1,LIM1
                  If(I.ne.BDRG(J,K)) goto 2065
                      X = X+1
                      RGBD(X,I) = K
 2065             Continue
 2070         Continue
 2080     NRGBD(I) = X
      Do 995 I = 1, M
          If(Debug) Write(IOut,9950) I, (RgBd(J,I),J=1,6)
 9950     Format(' I, RgBd: ',7I3)
  995     Continue
C
C     Now fill RGARG(J,I) with the Jth ring attached to ring I by
C     one atom, and RGBRG(J,I) with the Jth ring attached to ring I
C     by a bond.
C
      Do 2250 K = 1,NRG
          X = 0
          KK = 0
          LIM1 = RGSIZE(K)
          Do 2160 J = 1,LIM1
              Ind1 = BDRG(J,K)
              If(NRGBD(Ind1).le.1) goto 2160
                  LIM2 = NRGBD(Ind1)
                  Do 2150 I = 1,LIM2
                      If(RGBD(I,Ind1).eq.K) goto 2150
                          X = X + 1
                          RGBRG(X,K) = RGBD(I,Ind1)
                          Ind2 = AtBd(1,Ind1)
                          Ind3 = AtBd(2,Ind1)
                          NRGAT(Ind2) = NRGAT(Ind2)-1
                          NRGAT(Ind3) = NRGAT(Ind3)-1
 2150                 Continue
 2160         Continue
          LIM1 = RGSIZE(K)
          Do 2240 J = 1,LIM1
              Ind1 = ATRG(J,K)
              If(NRGAT(Ind1).le.1) goto 2240
                  LIM2 = NRGAT(Ind1)
                  Do 2230 I = 1,LIM2
                      If(RGAT(I,Ind1).ne.K) KK = KK + 1
                      If(RGAT(I,Ind1).ne.K) RGARG(KK,K) = RGAT(I,Ind1)
 2230                 Continue
 2240         Continue
 2250     Continue
C
C     EXSVAL(I) is the excess valence associated with atom I and helps
C     to define the local atom geometry about I as well as the type of
C     each bond containing I.
C
      Do 2260 I=1,M
          Ind1=AN(I+1)
 2260     EXSVAL(I+1)=VL(Ind1)-CN(I+1)
C
C     AtGeom(I) is the local atom geometry about atom I:
C         Line = Linear
C         Bent = Bent (Angle=Tetrahedral)
C         Trig = Trigonal
C         Pyra = Pyramidal (Angles=Tetrahedral)
C         Tetr = Tetrahedral
C
      Do 2440 I=1,M
          L=0
          Do 2270 J=1,MXCN
              Ind1=F(J,I)
 2270         L=L+EXSVAL(Ind1+1)
          J=AN(I+1)
          K=CN(I+1)
          AtGeom(I) = 0
          If(K.eq.4) AtGeom(I) = Tetr
          If(K.eq.3.and.(J.eq.6.or.(L.gt.0.and.(J.eq.7.or.J.eq.8))))
     $            AtGeom(I) = Trig
          If(K.eq.3.and..not.
     $        (J.eq.6.or.(L.gt.0.and.(J.eq.7.or.J.eq.8))))
     $        AtGeom(I) = Pyra
          If(K.eq.2.and.
     $        ((J.EQ.6.AND.L.GT.1).OR.(J.EQ.7.AND.L.GT.2)))
     $        AtGeom(I) = Line
          If(K.eq.2.and..not.
     $        ((J.EQ.6.AND.L.GT.1).OR.(J.EQ.7.AND.L.GT.2)))
     $        AtGeom(I) = Bent
 2440     Continue
C
C     The following is an option to alter the local atom geometry of
C     any atom.
C
 2470 JunkC=0
      Call PutCSt('AtomGeom',IJunk,JunkC)
      If(.not.StrEqC(NxCard,IJunk,8)) Goto 2590
      I = 0
      Call FFGet(IStrng,Junk,I,FJunk,IFound)
      If(IFound.ne.IFFInt) Call FFErr(IFFInt,IFound)
      NewAng(I) = Zero
      Call CBlank(20,Name)
      Call FFGet(Name,LenNam,Junk,NewAng(I),IFound)
      If(IFound.ne.IFFStr.and.IFound.ne.IFFFP) Call FFErr(IFFStr,IFound)
      Call Captlz(Name,Name,LenNam)
      NewAng(I) = NewAng(I) * AngCon
      If(IFound.eq.IFFFP) AtGeom(I) = ASym
      JunkC=0
      Call PutCSt('Tetr',IJunk,JunkC)
      If(StrEqC(Name,IJunk,4)) AtGeom(I) = Tetr
      If(StrEqC(Name,IJunk,4)) Goto 2480
      JunkC=0
      Call PutCSt('Pyra',IJunk,JunkC)
      If(StrEqC(Name,IJunk,4)) AtGeom(I) = Pyra
      If(StrEqC(Name,IJunk,4)) Goto 2480
      JunkC=0
      Call PutCSt('Bent',IJunk,JunkC)
      If(StrEqC(Name,IJunk,4)) AtGeom(I) = Bent
      If(StrEqC(Name,IJunk,4)) Goto 2480
      JunkC=0
      Call PutCSt('Trig',IJunk,JunkC)
      If(StrEqC(Name,IJunk,4)) AtGeom(I) = Trig
      If(StrEqC(Name,IJunk,4)) Goto 2480
      JunkC=0
      Call PutCSt('Line',IJunk,JunkC)
      If(StrEqC(Name,IJunk,4)) AtGeom(I) = Line
      If(StrEqC(Name,IJunk,4)) Goto 2480
      JunkC=0
      Call PutCSt('Asym',IJunk,JunkC)
      If(StrEqC(Name,IJunk,4)) AtGeom(I) = Asym
 2480 Call CBlank(20,NxCard)
      If(IEof.eq.0) Call FFRead(IEof)
      If(IEof.eq.1) goto 2590
      Call FFGet(NxCard,LenNx,Junk,FJunk,IFound)
      If(IFound.ne.IFFStr.and.IFound.ne.IFFEnd)
     $    Call FFErr(IFFStr,IFound)
      If(IFound.eq.IFFStr) goto 2470
      IEof = 1
C
C     Now all bond types are determined for model A:
C         Sing = Single Bond
C         Doub = Double
C         Trip = Triple
C         Arom = Aromatic
C         Taro = Triple aromatic (as in benzyne).
C         Datv = Dative
C
 2590 If(IModel.eq.1) goto 3710
      EV(1) = 0
      EVC(1) = 0
      Do 2620 I=1,M
 2620     EV(I+1) = EXSVAL(I+1)
      Do 2670 I=1,M
          EVC(I+1) = 0
          Do 2660 J=1,MXCN
              Ind1=F(J,I)
              If(EV(Ind1+1).gt.0) EVC(I+1) = EVC(I+1)+1
 2660         Continue
 2670     Continue
      Do 2680 I=2,NBD
 2680     BT(I+1) = 1
C
C     Now identify all dative bonds.
C
      AN(1) = 0
      CN(1) = 0
      Do 2880 I=1,M
          If(AN(I+1).eq.8.AND.CN(I+1).eq.1.AND.EV(I+1).GT.0) Goto 2730
              Do 2870 J=1,MXCN
                  Ind1 = F(J,I)
                  If(AN(Ind1+1).ne.7.or.CN(Ind1+1).ne.3) Goto 2870
                  Do 2860 K=1,MXCN
                      Ind1 = F(J,I)
                      Ind2 = F(K,Ind1)
                      If(AN(Ind2+1).ne.8.or.CN(Ind2+1).ne.1) Goto 2860
                      Ind3 = BDAT(J,I)
                      BT(Ind3+1) = 6
                      Ind4 = BDAT(K,Ind1)
                      BT(Ind4+1) = 6
                      EV(I+1) = EV(I+1)-1
                      EV(Ind2+1) = EV(Ind2+1)-1
                      EVC(Ind1+1) = EVC(Ind1+1) - 2
 2860                 Continue
 2870             Continue
C
C         N=N dative bonds.
C
 2730     If(AN(I+1).eq.7.AND.CN(I+1).eq.2) Goto 2880
              Do 2770 L=1,MXCN
                  Ind1 = F(L,I)
                  If(AN(Ind1+1).ne.7.or.CN(Ind1+1).ne.1) Goto 2880
                  EVC(I+1) = EVC(I+1)-1
                  EV(Ind1+1) = EV(Ind1+1)-2
                  Ind1 = BDAT(L,I)
 2770             BT(Ind1+1) = 6
 2880     Continue
C
C     First find all terminal multiple bods. keep cycling until all
C     of these have been found.
C
 2890 Do 2940 I=1,M
          If(EV(I+1).le.0.or.EVC(I+1).ne.0) Goto 2940
          If( (AN(I+1).eq.5.and.EV(I+1).eq.2).or.
     $        (AN(I+1).eq.7.and.EV(I+1).eq.2).or.
     $        (AN(I+1).eq.6.and.(EV(I+1).eq.2.or.EV(I+1).eq.4)))
     $        EV(I+1)=0
          If(EV(I+1).ne.0) Write(IOut,2930) I, An(I+1), EV(I+1),
     $                                      EVC(I+1)
 2930     Format(' Grave error: Atom',I3,' An=',I2,' has EV=',I4,
     $           ' and EVC=',I4)
          If(EV(I+1).ne.0) Call Lnk1E(0)
 2940     Continue
      NB2 = 0
 2950 NB1 = NB2
      Do 3190 I=1,M
          If(EV(I+1).le.0.or.EVC(I+1).ne.1) Goto 3190
          Do 3180 J=1,MXCN
              Ind1 = F(J,I)
              If(EV(Ind1+1).le.0.or.BDAT(J,I).le.0) Goto 3180
              EV(I+1) = EV(I+1)-1
              Ind1 = BDAT(J,I)
              BT(Ind1+1) = BT(Ind1+1)+1
              NB2 = NB2+1
              Ind1 = F(J,I)
              EV(Ind1+1) = EV(Ind1+1)-1
              If(EV(I+1).ne.0) goto 3090
                  Do 3080 L=1,MXCN
                      K = F(L,I)
                      If(K.GT.0.and.EVC(K+1).GT.0) EVC(K+1) = EVC(K+1)-1
 3080                 Continue
 3090         Ind1 = F(J,I)
              If(EV(Ind1+1).ne.0) goto 3180
                  Do 3160 L=1,MXCN
                      K = F(L,Ind1)
                      If(K.GT.0.and.EVC(K+1).GT.0) EVC(K+1) = EVC(K+1)-1
 3160                 Continue
 3180         Continue
 3190     Continue
      If(NB2.gt.NB1) goto 2950
C
C     Now all bonds have been identified except aromatic bonds.  If
C     TotEV=0 there ar+ no aromatic bonds.
C
      TotEV = 0
      Do 3270 I=1,M
          If(EV(I+1).gt.0) TotEV = TotEV+1
 3270     Continue
      If(TotEV.eq.0) goto 3710
C
C     Now identify aromatic bonds.
C
      Do 3350 K=1,NRG
          LIM1 = RGSIZE(K)
          Do 3310 J=1,LIM1
              Ind1 = ATRG(J,K)
              If(EV(Ind1+1).eq.0) goto 3350
 3310         Continue
          If(Mod((RgSize(K)-2),4).eq.0) RgType(K) = Arom
          If(Mod((RgSize(K)-2),4).ne.0) RgType(K) = Conj
 3350     Continue
      Do 3480 K = 1, NRg
          If(RgType(K).ne.Arom) goto 3480
          Lim1 = rgsize(k)
          Do 3470 J = 1, Lim1
              Ind1 = AtRg(j,k)
              If(EV(Ind1+1).gt.0) EV(Ind1+1) = EV(Ind1+1) - 1
              Ind1 = BdRg(J,K)
              If(BT(Ind1+1).eq.1) BT(Ind1+1) = BT(Ind1+1) + 3
 3470         Continue
 3480     Continue
      Do 3590 K=1,NRG
          If(RgType(K).NE.CONJ) goto 3590
          LIM1=RgSize(K)
          DO 3530 J=1,LIM1
              IND1 = AtRg(J,K)
              If(EV(IND1+1).EQ.0) goto 3590
 3530         CONTINUE
C
C         All aromatic and conjugated double bonds have been identified.
C         Now see if any more terminal multiple bonds have been
C         uncovered.
C
          Do 3580 J=1,LIM1
              Ind1 = BdRg(J,K)
              Ind2 = AtBd(1,Ind1)
              Ind3 = AtBd(2,Ind1)
              IEVCor = Min(EV(Ind2+1),EV(Ind3+1))
              If(IEVCor.le.0) goto 3580
              EV(Ind2+1) = EV(Ind2+1)-IEVCor
              EV(Ind3+1) = EV(Ind3+1)-IEvCor
              BT(Ind1+1) = BT(Ind1+1)+IEVCor
 3580         Continue
 3590     Continue
      Do 3640 I=1,M
          EVC(I+1)=0
          Do 3630 J=1,MXCN
              Ind1 = F(J,I)
              If(EV(Ind1+1).gt.0) EVC(I+1) = EVC(I+1)+1
 3630         Continue
 3640     Continue
      TotEV = 0
      Do 3680 I=1,M
          If(EV(I+1).gt.0) TotEV = TotEV+1
 3680     Continue
      If(TotEV.gt.0) goto 2890
C
C     All bond types have now been defined.  Now determine bond
C     rotational geometries of all pertinent bonds.
C
 3710 If(IPrint.gt.0) Write(IOut,1000)
      Call IniStr
      Call BldStC('  Atom   Symbol  Geom  Angle ')
      Do 3715 I = 1, MxCN
          If(I.gt.2) Call BldStC('  ')
 3715     Continue
      Call BldStC('Array F')
      If(IPrint.gt.0) Call DmpStr(IOut)
      Do 3720 I = 1, M
          Ind1 = An(I+1)
          IKSym = Blank
          If(.not.Keep(I)) IKSym = IDash
          If(IPrint.gt.0.and.AtGeom(I).eq.0)
     $        Write(IOut,1060) I,IKSym,EL(Ind1),(F(J,I),J=1,MXCN)
          If(IPrint.gt.0.and.AtGeom(I).eq.Tetr)
     $        Write(IOut,1070) I,IKSym,EL(Ind1),(F(J,I),J=1,MXCN)
          If(IPrint.gt.0.and.AtGeom(I).eq.Trig)
     $        Write(IOut,1080) I,IKSym,EL(Ind1),(F(J,I),J=1,MXCN)
          If(IPrint.gt.0.and.AtGeom(I).eq.Pyra)
     $        Write(IOut,1090) I,IKSym,EL(Ind1),(F(J,I),J=1,MXCN)
          If(IPrint.gt.0.and.AtGeom(I).eq.Bent)
     $        Write(IOut,1100) I,IKSym,EL(Ind1),(F(J,I),J=1,MXCN)
          If(IPrint.gt.0.and.AtGeom(I).eq.Line)
     $        Write(IOut,1110) I,IKSym,EL(Ind1),(F(J,I),J=1,MXCN)
          If(IPrint.gt.0.and.AtGeom(I).eq.ASym)
     $        Write(IOut,1120) I,IKSym,EL(Ind1),NewAng(I),
     $        (F(J,I),J=1,MXCN)
 3720     Continue
      Call IClear(2*M,Dih)
      Do 3760 I=2,NBD
          Ind1 = AtBd(1,I)
          Ind2 = AtBd(2,I)
          If(CN(Ind1+1).gt.1.and.CN(Ind2+1).gt.1) BdRot(I) = Tran
          If(CN(Ind1+1).le.1.or.CN(Ind2+1).le.1) BdRot(I) = 0
 3760     Continue
C
C     The following is an option to alter the bond rotational geometry
C     of any bond:
C
C     Gaup = Gauche = Dihedral angle of +60.
C     Gaum = Gauche prime = Dihedral angle of -60.
C     Skew = Dihedral angle specified by user.
C
 3790 JunkC=0
      Call PutCSt('BondRot',IJunk,JunkC)
      If(.not.StrEqC(NxCard,IJunk,7)) goto 4030
      Call FFGet(IStrng,Junk,I,FJunk,IFound)
      If(IFound.ne.IFFInt) Call FFErr(IFFInt,IFound)
      Call FFGet(IStrng,Junk,J,FJunk,IFound)
      If(IFound.ne.IFFInt) Call FFErr(IFFInt,IFound)
      Call FFGet(IStrng,Junk,K,FJunk,IFound)
      If(IFound.ne.IFFInt) Call FFErr(IFFInt,IFound)
      Call FFGet(IStrng,Junk,L,FJunk,IFound)
      If(IFound.ne.IFFInt) Call FFErr(IFFInt,IFound)
      Nurtng = Zero
      Call CBlank(20,Name)
      Call FFGet(Name,LenNam,Junk,Nurtng,IFound)
      If(IFound.ne.IFFStr.and.IFound.ne.IFFFP) Call FFErr(IFFFP,IFound)
      Call Captlz(Name,Name,NChrPW(0))
      Nurtng = Nurtng * AngCon
      IRBd = 0
      JunkC = 0
      Call PutCSt('Cis ',IJunk,JunkC)
      If(StrEqC(Name,IJunk,4)) IRBd = Cis
      JunkC = 0
      Call PutCSt('Tran',IJunk,JunkC)
      If(StrEqC(Name,IJunk,4)) IRBd = Tran
      JunkC = 0
      Call PutCSt('GauP',IJunk,JunkC)
      If(StrEqC(Name,IJunk,4)) IRBd = GauP
      JunkC = 0
      Call PutCSt('GauM',IJunk,JunkC)
      If(StrEqC(Name,IJunk,4)) IRBd = GauM
      If(IFound.eq.IFFFP) IRBd = Skew
      If(IPrint.gt.0.and.IRBd.eq.Cis) Write(IOut,1005) I,J,K,L
      If(IPrint.gt.0.and.IRBd.eq.Tran) Write(IOut,1010) I,J,K,L
      If(IPrint.gt.0.and.IRBd.eq.GauP) Write(IOut,1020) I,J,K,L
      If(IPrint.gt.0.and.IRBd.eq.GauM) Write(IOut,1030) I,J,K,L
      If(IPrint.gt.0.and.IRBd.eq.Skew) Write(IOut,1040) I,J,K,L,Nurtng
      Do 4010 X=1,NBD
          If(.not.(AtBd(1,X).eq.J.and.AtBd(2,X).eq.K).and..not.
     $       (AtBd(1,X).eq.K.and.AtBd(2,X).eq.J)) goto 4010
          BdRot(X) = IRBd
          Skrtng(X) = NuRtng
          If(J.eq.AtBd(1,X)) DIH(1,X)=I
          If(J.ne.AtBd(1,X)) DIH(2,X)=I
          If(K.eq.AtBd(1,X)) DIH(1,X)=L
          If(K.ne.AtBd(1,X)) DIH(2,X)=L
 4010     Continue
      Call CBlank(20,NxCard)
      If(IEof.eq.0) Call FFRead(IEof)
      If(IEof.eq.1) goto 4030
      Call FFGet(NxCard,LenNx,Junk,FJunk,IFound)
      If(IFound.ne.IFFStr.and.IFound.ne.IFFEnd)
     $    Call FFErr(IFFStr,IFound)
      If(IFound.eq.IFFStr) goto 3790
      IEof = 1
C
C     For each pertinent bond X define the two atoms DIH(1,X) and
C     DIH(2,X) whose positions are defined by the bond rotational
C     geometry.  The positions of all other atoms attached to X
C     are then fixed by their bondangles.
C
 4030 Do 4350 X=1,M
          Ind1 = AtBd(1,X)
          Ind2 = AtBd(2,X)
          If(DIH(1,X).ne.0.or.CN(Ind1+1).le.1.or.CN(Ind2+1).le.1)
     $        goto 4350
          KK = AtBd(1,X)
          LL = AtBd(2,X)
          Do 4190 J=1,MXCN
              If(F(J,KK).ne.LL) goto 4190
              LIM1 = J+1
              LIM2 = J+MXCN
              Do 4180 Y=LIM1,LIM2
                  If(Y.le.MXCN) JJ = Y
                  If(Y.gt.MXCN) JJ = Y-MXCN
                  If(F(JJ,KK).ne.0) DIH(1,X) = F(JJ,KK)
                  If(F(JJ,KK).ne.0) goto 4200
 4180             Continue
 4190         Continue
 4200     Do 4330 J=1,MXCN
              If(F(J,LL).ne.KK) goto 4330
              LIM1 = J+1
              LIM2 = J+MXCN
              Do 4320 Y=LIM1,LIM2
                  If(Y.le.MXCN) JJ = Y
                  If(Y.gt.MXCN) JJ = Y - MXCN
                  If(F(JJ,LL).ne.0) DIH(2,X) = F(JJ,LL)
                  If(F(JJ,LL).ne.0) goto 4350
 4320             Continue
 4330         Continue
 4350     Continue
      If(IPrint.gt.1) Write(IOut,4410)
 4410 Format(28X,'Bond information:')
      Do 4600 J=2,NBD
          BdType(J) = 0
          If(BT(J+1).eq.1) BdType(J) = Sing
          If(BT(J+1).eq.2) BdType(J) = Doub
          If(BT(J+1).eq.3) BdType(J) = Trip
          If(BT(J+1).eq.4) BdType(J) = Arom
          If(BT(J+1).eq.5) BdType(J) = TAro
          If(BT(J+1).eq.6) BdType(J) = Datv
 4600     Continue
      If(IPrint.gt.1) Write(IOut,1130)
      Do 4630 I=2,NBD
          Ind1 = AtBd(1,I)
          Ind2 = AN(Ind1+1)
          Ind3 = AtBd(2,I)
          Ind4 = AN(Ind3+1)
          ICur = 0
          If(I.le.9) Call PutCSt('    ',OutLin,ICur)
          If(I.ge.10) Call PutCSt('   ',OutLin,ICur)
          Call DecChr(I,OutLin,ICur)
          Call PutCSt('        ',OutLin,ICur)
          Call PutB(El(Ind2),2,OutLin,ICur)
          Call PutCSt('       ',OutLin,ICur)
          Call PutB(El(Ind4),2,OutLin,ICur)
          If(AtBd(1,I).le.9) Call PutCSt('        ',OutLin,ICur)
          If(AtBd(1,I).ge.10) Call PutCSt('        ',OutLin,ICur)
          Call DecChr(AtBd(1,I),OutLin,ICur)
          If(AtBd(2,I).le.9) Call PutCSt('   ',OutLin,ICur)
          If(AtBd(2,I).ge.10) Call PutCSt('  ',OutLin,ICur)
          Call DecChr(AtBd(2,I),OutLin,ICur)
          Call PutCSt('           ',OutLin,ICur)
          JCur = ICur
          If(BdType(I).eq.Conj) Call PutCSt('Conj',OutLin,ICur)
          If(BdType(I).eq.Sing) Call PutCSt('Sing',OutLin,ICur)
          If(BdType(I).eq.Doub) Call PutCSt('Doub',OutLin,ICur)
          If(BdType(I).eq.Trip) Call PutCSt('Trip',OutLin,ICur)
          If(BdType(I).eq.Arom) Call PutCSt('Arom',OutLin,ICur)
          If(BdType(I).eq.TAro) Call PutCSt('TAro',OutLin,ICur)
          If(BdType(I).eq.Datv) Call PutCSt('Datv',OutLin,ICur)
          If(BdType(I).eq.0)    Call PutCSt('    ',OutLin,ICur)
          If(ICur.eq.JCur)      Call PutCSt('----',OutLin,ICur)
          Call PutCSt('          ',OutLin,ICur)
          JCur = ICur
          If(BdRot(I).eq.Cis) Call PutCSt('Cis',OutLin,ICur)
          If(BdRot(I).eq.Tran) Call PutCSt('Tran',OutLin,ICur)
          If(BdRot(I).eq.GauP) Call PutCSt('GauP',OutLin,ICur)
          If(BdRot(I).eq.GauM) Call PutCSt('GauM',OutLin,ICur)
          If(BdRot(I).eq.Skew) Call PutCSt('Skew',OutLin,ICur)
          If(ICur.eq.JCur.and.BdRot(I).ne.0)
     $      Call PutCSt('----',OutLin,ICur)
          If(IPrint.gt.1) Call StrOut(IOut,OutLin,ICur,1)
 4630     Continue
      If(NRG.le.0) goto 9999
      If(IPrint.gt.1) Write(IOut,1140)
      Do 4700 K = 1, NRg
          Lim = RgSize(K)
          If(IPrint.gt.1.and.RgType(K).eq.0)
     $        Write(IOut,1150) K, RgSize(K), (AtRg(J,K),J=1,Lim)
          If(IPrint.gt.1.and.RgType(K).eq.Arom)
     $        Write(IOut,1160) K, RgSize(K), (AtRg(J,K),J=1,Lim)
          If(IPrint.gt.1.and.RgType(K).eq.Conj)
     $        Write(IOut,1170) K, RgSize(K), (AtRg(J,K),J=1,Lim)
 4700     Continue
 9999 Return
      End
*Deck MBld3
      Subroutine MBld3(IModel,IPrint,MaxNZ,El,IUConv,ToAng,M,An,C,Keep,
     $  Z,BL,Alpha,Beta,IAnZ,IEof,OldNo,Test,From,NewNo,BdDef,BndLth,
     $  BdAng,RotAng,A,B,D,TE,CZ,AScr)
      Implicit Real*8(A-H,O-Z)
C
C     This segment of model-builder uses the qualitative model built in
C     MBld2 to pick out standard bond lengths and bond angles and
C     dihedral angles and calculates the cartesian coordinates of all
C     atoms.  Note that Z should be allocated as (MaxNZ,4) even though
C     it is dimensioned here as (MaxAtm,4) -- it is expanded at the end.
C
      Parameter (MaxMBA=DEFMAXATM,MaxMB1=(MaxMBA+1))
      Dimension An(*),C(3,*),IAnZ(*),IJunk(10),BndLth(*),BdAng(*),
     $  RotAng(*),A(*),B(*),D(*),BL(*),Alpha(*),Beta(*),CZ(3,*),
     $  AScr(MaxNZ,5),V1(3),V2(3),V3(3),VP(3),VJ(3)
      Common/IO/In,IOut,IPunch
      Integer OldNo(*),Z(MaxMBA,4),Test(*),From(*),NewNo(*),BdDef(*),
     $  TETR,TRIG,PYRA,BENT,ASYM,TRAN,GAUP,GAUM,SKEW,SING,DOUB,TRIP,
     $  AROM,TARO,DATV,AN,CN,BT,EXSVAL,F,BDAT,BDBD,RGBD,RGAT,ATBD,ATRG,
     $  BDRG,RGARG,RGBRG,VL,RGSIZE,DIH,RGTYPE,ATGEOM,BDTYPE,BDROT,X,
     $  EL(-2:*)
      Common/MB/MXCN,NRG,NBD,MXRGSZ,CN(MaxMB1),BT(MaxMB1),
     $  EXSVAL(MaxMB1),F(6,MaxMBA),BDAT(6,MaxMBA),BDBD(6,MaxMBA),
     $  RGBD(6,MaxMBA),RGAT(6,MaxMBA),ATBD(2,MaxMBA),ATRG(10,10),
     $  BDRG(10,10),RGARG(10,10),RGBRG(10,10),NRGAT(MaxMBA),
     $  NRGBD(MaxMBA),VL(18),RGSIZE(10),DIH(2,MaxMBA),RGTYPE(10),
     $  ATGEOM(MaxMBA),BDTYPE(MaxMBA),BDROT(MaxMBA),NXCARD(20),IMBPad,
     $  SKRTNG(MaxMBA),NEWANG(MaxMBA)
      Logical IfErr,StrEqC,Te(*),Keep(*)
      Real*8 NULNTH,NEWANG,L1(3),L2(3),L3(3),L4(3)
      Save TETR,TRIG,PYRA,BENT,LINE,ASYM,TRAN,GAUP,GAUM,SKEW
      Save SING,DOUB,TRIP,AROM,TARO,DATV,Zero,One,Two,Three
      Save Four,F120,F180,IFFFP,IFFInt,IFFStr,IFFEnd,DatNO,DatNN
      Save OCNLen, C2C2Ln, C2N1Ln,RN1N1L, C3C3Ln, C3N2Ln,RN2N2L, TAroLn
      DATA TETR/1/,TRIG/2/,PYRA/3/,BENT/4/,LINE/5/,ASYM/6/,
     $  TRAN/8/,GAUP/9/,GAUM/10/,SKEW/11/,SING/13/,DOUB/14/,
     $  TRIP/15/,AROM/16/,TARO/17/,DATV/18/
      Data Zero/0.0d0/, One/1.0d0/, Two/2.0d0/, Three/3.0d0/,
     $  Four/4.0d0/, F120/120.0d0/, F180/180.0d0/
      Data IFFFP/2HFP/, IFFInt/3HINT/, IFFSTR/3HSTR/, IFFEND/3HEND/
C
C     Special Model A bonds.  There are special values for C-N in
C     N-C=O, and for a few triple and aromatic bonds.
C
      Data OCNLen/1.32d0/, C2C2Ln/1.20d0/, C2N1Ln/1.16d0/,
     $     RN1N1L/1.10d0/, C3C3Ln/1.40d0/, C3N2Ln/1.34d0/,
     $     RN2N2L/1.35d0/, TAroLn/1.30d0/, DatNO/1.12d0/,
     $     DatNN/1.24d0/
 1000 Format(' Set the ',I2,' -',I2,' bond to ',F12.6)
 1010 Format('               Model Z-Matrix:')
 1020 Format(I3,I4,F7.4,I4,F11.6,I4,F11.6,I4)
 1040 Format(20X,'Coordinates:',/'    No of atom',2X,'Symbol Keep?',
     $8X,1HX,11X,1HY,11X,1HZ)
 1050 Format(8X,I2,8X,A2,4X,L2,4X,3F12.8)
C
      PI = Four * ATan(One)
      TetAng = F180 * ACos(-One/Three) / Pi
      I=0
      J=0
      K=0
      L=0
      X=0
      KK=0
      KK1=0
      KK2=0
      KK3=0
      LL=0
      MM=0
      JJ=0
      JJ2=0
      JJ3=0
      JJ4=0
      Call AClear(MaxMBA,Bl)
      Call AClear(MaxMBA,Alpha)
      Call AClear(MaxMBA,Beta)
      Call AClear(MaxMBA,Bndlth)
      Call AClear(MaxMBA,Bdang)
      Call AClear(MaxMBA,RotAng)
      Call AClear(MaxMBA,A)
      Call AClear(MaxMBA,B)
      Call AClear(MaxMBA,D)
      Call AClear(3*MaxMBA,C)
      Call AClear(3,V1)
      Call AClear(3,L1)
      Call AClear(3,V2)
      Call AClear(3,L2)
      Call AClear(3,V3)
      Call AClear(3,L3)
      Call AClear(3,L4)
      Call AClear(3,VP)
      Call AClear(3,VJ)
      Call IClear(4*MaxMBA,Z)
      BlConv = One
      If(Mod(IUConv,2).eq.1) BlConv = ToAng
C
C     Pick out standard model A or model B bondlengths for each bond
C     unless its length has been explicitly specified.
C
      Call IClear(MaxMB1,Test)
  200 JunkC = 0
      Call PutCSt('BondLen',IJunk,JunkC)
      If(.not.StrEqC(NxCard,IJunk,7)) goto 320
      Call FFGet(IStrng,Junk,I,FJunk,IFound)
      If(IFound.ne.IFFInt) Call FFErr(IFFInt,IFound)
      Call FFGet(IStrng,Junk,J,FJunk,IFound)
      If(IFound.ne.IFFInt) Call FFErr(IFFInt,IFound)
      Call FFGet(IStrng,Junk,Junk,Nulnth,IFound)
      If(IFound.ne.IFFFP) Call FFErr(IFFFP,IFound)
      Nulnth = NuLnth * BlConv
      If(IPrint.gt.0) Write(IOut,1000) I, J, Nulnth
      Do 260 K=1,MXCN
          If(F(K,I).ne.J) goto 260
              IND1=BDAT(K,I)
              BNDLTH(IND1)=NULNTH
              TEST(IND1+1)=1
  260     Continue
      Call CBlank(20,NxCard)
      If(IEof.eq.0) Call FFRead(IEof)
      If(IEof.eq.1) goto 320
      Call FFGet(NxCard,LenNx,Junk,FJunk,IFound)
      If(IFound.ne.IFFEnd.and.IFound.ne.IFFStr)
     $    Call FFErr(IFFStr,IFound)
      If(IFound.eq.IFFStr) goto 200
      IEof = 1
C
  320 Do 3250 I=2,NBD
          If(Test(I+1).ne.0) goto 3250
          Ind1 = AtBd(1,I)
          Ind2 = AtBd(2,I)
          J = Min(CN(Ind1+1),CN(Ind2+1))
          K = Max(CN(Ind1+1),CN(Ind2+1))
          L = Min(An(Ind1+1),An(Ind2+1))
          KK = Max(An(Ind1+1),An(Ind2+1))
          BndLth(I) = Zero
          If(IModel.eq.1) goto 3240
          If(BDType(I).ne.Sing) goto 350
              BndLth(I) = BlModA(An(Ind1+1),Cn(Ind1+1),An(Ind2+1),
     $                           Cn(Ind2+1),BdType(I))
              If(KK.ne.7.or.L.ne.6) goto 3240
              LL = 0
              If(An(AtBd(1,I)+1).eq.6) LL = 1
              If(An(AtBd(2,I)+1).eq.6) LL = 2
              If(LL.eq.0) goto 3240
              Do 340 MM = 1, MXCN
                  Ind3 = BDat(MM,AtBd(LL,I))
                  Ind4 = AtBd(1,Ind3) + 1
                  Ind5 = AtBd(2,Ind3) + 1
                  If(BT(Ind3+1).eq.2.and.
     $                (An(Ind4).eq.8.or.An(Ind5).eq.8))
     $                BndLth(I) = OCNLen
  340             Continue
              Goto 3240
  350     If(BDType(I).ne.Doub) goto 360
              BndLth(I) = BlModA(An(Ind1+1),Cn(Ind1+1),An(Ind2+1),
     $                           Cn(Ind2+1),BdType(I))
              Goto 3240
  360     If(BDType(I).ne.Trip) goto 370
              If(L.eq.6.and.KK.eq.6) Bndlth(I) = C2C2Ln
              If(L.eq.6.and.KK.eq.7) BndLth(I) = C2N1Ln
              If(L.eq.7.and.KK.eq.7) BndLth(I) = RN1N1L
              Goto 3240
  370     If(BDType(I).ne.Arom) goto 380
              If(L.eq.6.and.KK.eq.6) Bndlth(I) = C3C3Ln
              If(L.eq.6.and.KK.eq.7) BndLth(I) = C3N2Ln
              If(L.eq.7.and.KK.eq.7) BndLth(I) = RN2N2L
              Goto 3240
  380     If(BdType(I).ne.TAro) goto 390
              BndLth(I) = TAroLn
              Goto 3240
  390     If(BdType(I).ne.Datv) goto 3240
              If(L.eq.7.and.KK.eq.8) BndLth(I) = DatNO
              If(L.eq.7.and.KK.eq.7) BndLth(I) = DatNN
 3240     If(BndLth(I).eq.Zero) BndLth(I) = BlModB(L,KK)
 3250     Continue
C
C     All bondlengths have now been found. now find values for
C     bond angles.
C
      Do 3450 I=1,M
          BdAng(I) = Zero
          If(AtGeom(I).eq.ASym) BdAng(I) = NewAng(I)
          If(AtGeom(I).eq.Tetr) BdAng(I) = TetAng
          If(AtGeom(I).eq.Pyra) BdAng(I) = TetAng
          If(AtGeom(I).eq.Bent) BdAng(I) = TetAng
          If(AtGeom(I).eq.Trig) BdAng(I) = F120
          If(AtGeom(I).eq.Line) BdAng(I) = F180
 3450     Continue
C
C     Now find dihedral angles.
C
      Do 3620 I=2,NBD
          RotAng(I) = Zero
          If(BdRot(I).eq.Tran) RotAng(I) = F180
          If(BdRot(I).eq.GauP) RotAng(I) = F120 / Two
          If(BdRot(I).eq.GauM) RotAng(I) = -F120 / Two
          If(BdRot(I).eq.Skew) RotAng(I) = Skrtng(I)
 3620     Continue
C
C     Now fill standard Z-matrix.  First put atoms in new numbering
C     system.
C
      Call IClear(MaxMB1,BDDef)
      Call IClear(MaxMB1,NewNo)
      Call IClear(MaxMB1,From)
      ILim = Min(M,3)
      Do 3670 I=1,ILIM
          IND1 = AtBd(2,I)
          NEWNO(IND1+1) = I
          OLDNO(I) = AtBd(2,I)
 3670     FROM(I+1) = AtBd(1,I)
      If(M.le.3) goto 4040
      NRING=0
      Do 4030 I = 4, M
          If(NRing.ne.1.and.AtBd(1,I).eq.AtBd(1,I-1)) goto 3960
          If(NRing.ne.1) LL = I
          If(NRing.eq.1) NRing = 0
          Do 3910 J = 1, MXCN
              L = AtBd(2,I)
              If(BDat(J,L).eq.0) goto 3910
                  IND1 = BDAT(J,L)
                  If(AtBd(1,Ind1).eq.L) K = AtBd(2,Ind1)
                  If(AtBd(1,Ind1).ne.L) K = AtBd(1,Ind1)
                  Do 3890 KK = 1,MXCN
                      If(BDAT(KK,K).eq.0) goto 3890
                      IND1 = BDAT(KK,K)
                      IF(L.ne.DIH(2,IND1)) goto 3890
                      NEWNO(L+1) = LL
                      OLDNO(LL)  = L
                      IND1 = NEWNO(L+1)
                      BDDEF(IND1+1) = BDAT(KK,K)
                      FROM(IND1+1) = AtBd(1,I)
                      Goto 3970
 3890                 Continue
 3910         Continue
          If(AtBd(1,I).eq.AtBd(1,I+1).and.I.ne.M) NRing = 1
          If(AtBd(1,I).eq.AtBd(1,I+1).and.I.ne.M) Goto 4030
          NEWNO(L+1) = LL
          OLDNO(LL)=L
          FROM(LL+1)=AtBd(1,LL)
          Goto 3970
 3960     IND1 = AtBd(2,I)
          NEWNO(IND1+1) = I
          OLDNO(I) = AtBd(2,I)
          FROM(I+1) = AtBd(1,I)
          LL = I
 3970     IF(LL.ge.I) goto 4030
              LIM1 = I-1
              Do 4010 JJ = LL,LIM1
                  IND1 = AtBd(2,JJ)
                  NEWNO(IND1+1) = JJ+1
                  OLDNO(JJ+1) = AtBd(2,JJ)
                  IND2 = NEWNO(IND1+1)
 4010             FROM(IND2+1) = AtBd(1,JJ)
 4030     Continue
C
C     All atoms have been renumbered.  now set up standard z matrix.
C
 4040 DO 4050 I = 2, M
 4050     FROM(I+1) = NEWNO(From(I+1)+1)
      Z(1,1) = 1
      Do 4130 I = 2, M
          If(BDDEF(I+1).eq.0) Z(I,1) = FROM(I+1)
          If(BDDEF(I+1).eq.0) goto 4130
          IND1 = BDDEF(I+1)
          IND2 = AtBd(2,IND1)
          IND3 = AtBd(1,IND1)
          If(NEWNO(IND2+1).LT.I.AND.NEWNO(IND3+1).LT.I)
     $        Z(I,1) = NewNo(Ind2+1)
          If(NEWNO(IND2+1).LT.I.AND.NEWNO(IND3+1).LT.I) goto 4130
          Z(I,1) = From(I+1)
 4130     Continue
      If(M.lt.3) goto 4930
      Do 4250 I = 3,M
          If(BDDEF(I+1).eq.0) goto 4160
              IND1 = BDDEF(I+1)
              IND2 = AtBd(1,IND1)
              IND3 = AtBd(2,IND1)
              If(NEWNO(IND2+1).LT.I.AND.NEWNO(IND3+1).LT.I)
     $            Z(I,2) = NEWNO(IND2+1)
              If(NEWNO(IND2+1).LT.I.AND.NEWNO(IND3+1).LT.I)
     $            Goto 4250
 4160     If(FROM(I+1).eq.1) goto 4190
              IND1 = FROM(I+1)
              Z(I,2) = FROM(IND1+1)
              Goto 4250
 4190     Continue
              LimJ = I - 2
              Do 4200 J = 1, LimJ
                  If(Z(I-J,1).eq.Z(I,1)) Z(I,2) = I - J
 4200             Continue
 4250     Continue
      If(M.lt.4) goto 4930
      Do 4910 I = 4,M
          X = 0
          Do 4300 J = 1,MXCN
              IND1 = Z(I,1)
              IND2 = OLDNO(IND1)
              IND3 = F(J,IND2)
              IF(NEWNO(IND3+1).LT.I.AND.F(J,IND2).GT.0) X = X + 1
 4300         Continue
          Goto 4890
 4305     Do 4400 J = 1,MXCN
              IND1 = Z(I,1)
              IND2 = OLDNO(IND1)
              IND3 = F(J,IND2)
              K = NEWNO(IND3+1)
              IF(K.LT.I.AND.K.NE.Z(I,2).AND.K.GT.0) Z(I,3) = K
              IF(K.LT.I.AND.K.NE.Z(I,2).AND.K.GT.0) Goto 4410
 4400         Continue
 4410     KK = OLDNO(I)
          IND1 = Z(I,1)
          KK1 = OLDNO(IND1)
          IND1 = Z(I,2)
          KK2 = OLDNO(IND1)
          IND1 = Z(I,3)
          KK3 = OLDNO(IND1)
          If(ATGEOM(KK1).eq.TRIG) Z(I,4) = 1
          IF(.not.(ATGEOM(KK1).EQ.TETR.OR.ATGEOM(KK1).EQ.PYRA.OR.
     $       ATGEOM(KK1).EQ.ASYM)) goto 4910
          Do 4590 J = 1,4
              If(F(J,KK1).eq.KK) JJ = J
              If(F(J,KK1).eq.KK2) JJ2 = J
              If(F(J,KK1).eq.KK3) JJ3 = J
              If(F(J,KK1).ne.KK.and.F(J,KK1).ne.KK2.and.F(J,KK1).ne.KK3)
     $            JJ4 = J
 4590         Continue
          LL = 0
          If(JJ.gt.JJ2) LL = LL + 1
          If(JJ.gt.JJ3) LL = LL + 1
          If(JJ.gt.JJ4) LL = LL + 1
          If(JJ2.gt.JJ3) LL = LL + 1
          If(JJ2.gt.JJ4) LL = LL + 1
          If(JJ3.gt.JJ4) LL = LL + 1
          If(Mod(LL,2).eq.1) Z(I,4) = +1
          If(Mod(LL,2).eq.0) Z(I,4) = -1
          Goto 4910
 4890     Z(I,4) = 0
          IND1 = BDDEF(I+1)
          If(Ind1.eq.0) goto 4305
          IND2 = DIH(1,IND1)
          If(Ind2.eq.0) goto 4305
          Z(I,3) = NEWNO(IND2+1)
          If(Z(I,3).eq.0) goto 4305
 4910     Continue
 4930 Continue
C
C     End of z definition now fill Bl.
C
      Do 4990 I = 2, M
          Do 4980 J = 1, MXCN
              Ind1 = OldnO(I)
              Ind2 = Z(I,1)
              Ind3 = BDat(J,Ind1)
              If(F(J,Ind1).eq.OldNo(Ind2)) BL(I)=BNDLTH(IND3)
              If(F(J,Ind1).eq.OldNo(Ind2)) Goto 4990
 4980         Continue
 4990     Continue
C
C     Now fill alpha with bondangles.
C
      If(M.lt.3) Goto 5190
      Do 5010 I = 3, M
          Ind1 = Z(I,1)
          Ind2 = OldNo(IND1)
 5010     Alpha(I) = BdAng(Ind2)
C
C     Fill beta with dihedral angles if Z(I,4)=0 and bond angles
C     otherwise.
C
      If(M.lt.4) goto 5190
      Do 5120 I = 4, M
          Ind1 = Z(I,1)
          If(Z(I,4).eq.0) goto 5050
              Ind2 = OldNo(Ind1)
              Beta(I) = BdAng(Ind2)
              Goto 5120
 5050     Continue
              K = OldNo(Ind1)
              Ind2 = Z(I,2)
              L = OldNo(Ind2)
              Do 5100 J = 1, MXCN
                  Ind3 = BDat(J,K)
                  If((AtBd(2,Ind3).eq.L).or.(AtBd(1,Ind3).eq.L))
     $                Beta(I) = RotAng(Ind3)
                  If((AtBd(2,Ind3).eq.L).or.(AtBd(1,Ind3).eq.L))
     $                goto 5120
 5100             Continue
 5120     Continue
C
C     Print Z matrix.
C
 5190 If(IPrint.gt.1) Write(IOut,1010)
      Do 5012 I = 1, M
          IAnZ(I) = An(OldNo(I)+1)
          If(IPrint.gt.1) Write(IOut,1020) IAnZ(I),Z(I,1),Bl(I),Z(I,2),
     $        Alpha(I),Z(I,3),Beta(I),Z(I,4)
 5012     Continue
      Call AScale(M,(PI/F180),Alpha,Alpha)
      Call AScale(M,(PI/F180),Beta,Beta)
      Call ZToC(MaxMBA,M,IAnZ,Z,BL,Alpha,Beta,.false.,NAtOut,An,C,CZ,
     $  AScr,AScr(1,2),AScr(1,3),AScr(1,4),AScr(1,5),IOut,IfErr,.False.)
      If(IfErr) Call Lnk1E(0)
C
C     Now reorder Keep to the new numbering system.
C
      Do 5700 I = 1,M
 5700     TE(I) = Keep(OldNo(I))
      Do 5720 I = 1, M
 5720     Keep(I) = TE(I)
C
C     Skip to the end of the input section.
C
 5900 If(IEof.eq.1) goto 5910
      Call FFRead(IEof)
      If(IEof.eq.1) goto 5910
      Call FFGet(NxCard,LenNx,Junk,FJunk,IFound)
      If(IFound.eq.IFFEnd) IEof = 1
      Goto 5900
 5910 If(IPrint.gt.1) Write(IOut,1040)
      DO 5920 J=1,M
        If(IPrint.gt.1) Write(IOut,1050)
     $    J,EL(An(J)),Keep(J),(C(I,J),I=1,3)
 5920   Continue
      Call UnPckI(0,Junk,Z,MaxNZ,4,MaxMBA,4)
      Return
      End
*Deck MkMpZA
      Subroutine MkMpZA(NZ,IAnZ,NAtoms,MapZAt)
      Implicit Integer(A-Z)
C
C     Generate a map from Z-matrix numbers to atom numbers.
C
      Dimension IAnZ(*), MapZAt(*)
C
      NAt = 0
      Do 10 IZ = 1, NZ
        If(IAnZ(IZ).ge.0)then
          NAt = NAt + 1
          MapZAt(IZ) = NAt
        else
          MapZAt(IZ) = 0
          endIf
   10   Continue
      If(NAt.ne.NAtoms) Call GauErr('Consistency failure in MkMpZA.')
      Return
      End
*Deck MrgCon
      Subroutine MrgCon(IOut,IPrint,IOCon,NAtoms,IAn,MaxCon,Wt1,Wt2,
     $  NBond1,IBond1,RBTyp1,NBond2,IBond2,RBTyp2)
      Implicit Real*8(A-H,O-Z)
C
C     Merge two connectivity structures.  Note that these are assumed
C     to be sorted (i.e., by SrtCon).
C
      Dimension IAn(*), NBond1(*), IBond1(MaxCon,*), RBTyp1(MaxCon,*),
     $  NBond2(*), IBond2(MaxCon,*), RBTyp2(MaxCon,*)
C
      Do 50 IA = 1, NAtoms
        Do 10 IB2 = 1, NBond2(IA)
   10     RBTyp2(IB2,IA) = Wt2*RBTyp2(IB2,IA)
        IB2N = 1
        Do 40 IB1 = 1, NBond1(IA)
          IB1A = IBond1(IB1,IA)
          ISt2 = IB2N
          IEnd2 = NBond2(IA)
          Do 20 IB2 = ISt2, IEnd2
            If(IB1A.eq.IBond2(IB2,IA)) then
C             Match
              RBTyp2(IB2,IA) = Wt1*RBTyp1(IB1,IA) + RBTyp2(IB2,IA)
              IB2N = IB2 + 1
              Goto 40
            else if(IB1A.lt.IBond2(IB2,IA)) then
C             No match; insert at IB2
              Goto 30
              endIf
   20       Continue
C         Insert at end
          IB2 = IEnd2 + 1
   30     NBond2(IA) = NBond2(IA) + 1
          If(NBond2(IA).gt.MaxCon)
     $      Call GauErr('MaxCon exceeded in MrgCon.')
          NMove = NBond2(IA) - IB2
          Call IMove1(NMove,IB2-1,IB2,IBond2(1,IA))
          IBond2(IB2,IA) = IB1A
          Call AMove1(NMove,IB2-1,IB2,RBTyp2(1,IA))
          RBTyp2(IB2,IA) = Wt1*RBTyp1(IB1,IA)
   40     Continue
   50   Continue
      If(IPrint.ge.2) Call PrtCon(IOut,'Connectivity in MrgCon:',.True.,
     $  NAtoms,MaxCon,IAn,NBond2,IBond2,IBond2,RBTyp2,0)
      If(IOCon.ne.0)
     $  Call RWCon(1,IOCon,NAtoms,MaxCon,MxBond,NBond2,IBond2,RBTyp2)
      Return
      End
*Deck NEquiv
      Function NEquiv(MaxAP3,NOp,IAtom,NPerm)
      Implicit Integer(A-Z)
C
C     Return the number of atoms equivalent to IAtom, counting itself.
C
      Dimension NPerm(MaxAP3,1)
C
      NEq = 1
      If(NOp.le.1) goto 999
      Do 20 I = 2, NOp
        LimJ = I - 1
        Do 10 J = 1, LimJ
          If(NPerm(IAtom,J).eq.NPerm(IAtom,I)) goto 20
   10     Continue
        NEq = NEq + 1
   20   Continue
  999 NEquiv = NEq
      Return
      End
*Deck NewIZ
      Subroutine NewIZ(IOut,MaxNZ,MaxAP3,NAtoms,NOp,MxBond,NPerm,NBond,
     $  IBond,IZ,IAn,C,NToO,OToN,Used,IScr,Scr,NPermS,IBType,IBTS)
      Implicit Real*8(A-H,O-Z)
C
C     Choose new bond pairs for IZ using the symmetry information in
C     NPerm and the connectivity information in NBond and IBond.  The
C     resulting permutation of atoms and its inverse are returned in
C     NToO and OToN.  This routine assumes that NAtoms=NZ.  Used, IScr,
C     Scr, IBTS and NPermS are scratch; Used and IScr can be the same,
C     as can NPermS and IBTS.
C
      Logical Used(NAtoms)
      Integer Bonded, OToN(NAtoms)
      Dimension NPerm(MaxAP3,1), NBond(NAtoms), IBond(MxBond,NAtoms),
     $  IZ(MaxNZ,4), IAn(NAtoms), C(3,NAtoms), NToO(NAtoms),
     $  IScr(NAtoms), Scr(3,NAtoms), NPermS(MaxAP3,1),
     $  IBType(MxBond,1), IBTS(MxBond,1)
C
C     Make atom 1 a unique atom, or at least a member of the smallest
C     possible group.
C
      Call IClear(4*MaxNZ,IZ)
      Call IClear(NAtoms,NToO)
      Call IClear(NAtoms,OToN)
      MinEqv = NAtoms + 1
      MaxEqv = 0
      IMin = 0
      Do 10 I = 1, NAtoms
        Used(I) = .False.
        NEqI = NEquiv(MaxAP3,NOp,I,NPerm)
        MaxEqv = Max(NEqI,MaxEqv)
        If(NEqI.lt.MinEqv) then
          IMin = I
          MinEqv = NEqI
          endIf
   10   Continue
      NDone = 1
      IZ(1,1) = 1
      NToO(NDone) = IMin
      OToN(IMin) = NDone
      Used(IMin) = .True.
C
C     Next add the atoms bonded to 1-NTest, with the least symmetric
C     first.  Repeat this section until all atoms are accounted for.
C
  100 NTest = NDone
      Do 120 IEq = 1, MaxEqv
        Do 110 I = 1, NAtoms
          If(.not.Used(I).and.NEquiv(MaxAP3,NOp,I,NPerm).le.IEq) then
            ICon = Bonded(MxBond,OToN,NToO,I,1,NTest,NBond,IBond)
            If(ICon.ne.0) then
              NDone = NDone + 1
              NToO(NDone) = I
              OToN(I) = NDone
              Used(I) = .True.
              IZ(NDone,1) = ICon
              endIf
            endIf
  110     Continue
  120   Continue
      If(NDone.eq.NTest) then
        Do 130 I = 1, NAtoms
          If(.not.Used(I)) then
            NDone = NDone + 1
            NToO(NDone) = I
            OToN(I) = NDone
            Used(I) = .True.
            IZ(NDone,1) = 1
            Goto 140
            endIf
  130     Continue
  140   Continue
C       Write(IOut,1000) NTest
C       Write(IOut,1010) MxBond, NAtoms, MaxEqv
C       Write(IOut,1020) NToO
C       Write(IOut,1030) OToN
C       Do 125 I = 1, NAtoms
C         L = NBond(I)
C         Write(IOut,1040) I,L,(IBond(J,I),J=1,L)
C 125     Continue
C       Call Lnk1E(0)
        endIf
      If(NDone.lt.NAtoms) goto 100
C
C     Now permute IAn, C, NBond, IBond, and IBType in accord with IZ.
C
      Call IMove(NAtoms,IAn,IScr)
      Call AMove(3*NAtoms,C,Scr)
      Call IMove(MaxAP3*NOp,NPerm,NPermS)
      Do 200 I = 1, NAtoms
        Call AMove(3,Scr(1,NToO(I)),C(1,I))
  200   IAn(I) = IScr(NToO(I))
      Do 220 I = 1, NAtoms
        Do 220 J = 1, NOp
  220     NPerm(I,J) = OToN(NPermS(NToO(I),J))
      Call IMove(NAtoms,NBond,IScr)
      Do 230 I = 1, NAtoms
  230   NBond(I) = IScr(NToO(I))
      N = MxBond*NAtoms
      Call IMove(N,IBond,IBTS)
      Call IClear(N,IBond)
      Do 240 I = 1, NAtoms
        Lim = NBond(I)
        IOld = NToO(I)
        Do 240 J = 1, Lim
  240     IBond(J,I) = OToN(IBTS(J,IOld))
      Call IMove(N,IBType,IBTS)
      Call IClear(N,IBType)
      Do 250 I = 1, NAtoms
        Lim = NBond(I)
        IOld = NToO(I)
        Do 250 J = 1, Lim
  250     IBType(J,I) = IBTS(J,IOld)
      Return
      End
*Deck NewVar
      Subroutine NewVar(IOut,IConv,ToAng,NVar,ANames,Values,IntVec,
     $  FPVec,MaxNZ,IZ,NZ,LBl,LAlpha,LBeta)
      Implicit Real*8(A-H,O-Z)
C
C     Read in new values for some variables.  IConv controls
C     conversion between angstroms/bohr and degrees/radians.
C
      Parameter (MaxNam=12)
      Logical RdZVar, StrEqC, GotVal, GotInt, GotFP, VDist
      Dimension ANames(*), Values(*), IntVec(*), FPVec(*), Name(MaxNam),
     $  IName(MaxNam), IZ(MaxNZ,4), LBl(*), LAlpha(*), LBeta(*), JJ(1)
      Save Zero, One, F45, JJ
      Data Zero/0.0d0/, One/1.0d0/, F45/45.0d0/, JJ/0/
 1000 Format(' Replacement values for variables:')
C
      Call SZPrnt(0,JJ,0,0)
      ConLen = One / ToAng
      If(Mod(IConv,2).eq.1) ConLen = One
      AngCon = ATan(One) / F45
      If(Mod(IConv,4).ge.2) AngCon = One
      Write(IOut,1000)
      Int = 0
      FP = Zero
      LenNam = 0
      Call CBlank(MaxNam,Name)
      Val = Zero
      GotFP = .False.
      GotInt = .False.
      GotVal = .False.
C
C     Loop over variable cards.
C
   10 If(.not.RdZVar(.True.,.False.,JJ,0,Name,MaxNam,LenNam,GotVal,Val,
     $  GotInt,Int,GotFP,FP)) Return
      Call SZPrnt(1,Name,LenNam,4)
      If(GotVal) Call SZPrnt(1,Val,0,2)
      If(GotInt) Call SZPrnt(1,Int,0,1)
      If(GotFP) Call SZPrnt(1,FP,0,2)
      Call SZPrnt(2,JJ,0,0)
      ICur = 0
      Do 100 IVar = 1, NVar
        Call GetB(2,IName,INameL,ANames,ICur)
        If(INameL.eq.LenNam.and.StrEqC(Name,IName,LenNam)) goto 200
  100   Continue
      Call IniStr
      Call BldStC('Attempt to redefine unrecognized symbol "')
      Call BldStr(Name,LenNam)
      Call BldStC('".')
      Call DmpStr(IOut)
      Call Lnk1E(0)
C
C     Check whether to convert units and store results.
C
  200 If(GotVal) then
        If(VDist(IVar,MaxNZ,NZ,LBl,LAlpha,LBeta,IZ)) then
          Values(IVar) = Val * ConLen
        else
          Values(IVar) = Val * AngCon
          endIf
        endIf
      If(GotInt) IntVec(IVar) = Int
      If(GotFP.and.Int.ne.1) then
        If(VDist(IVar,MaxNZ,NZ,LBl,LAlpha,LBeta,IZ)) then
          FPVec(IVar) = FP * ConLen
        else
          FPVec(IVar) = FP * AngCon
          endIf
        endIf
      Goto 10
      End
*DECK NOSIGN
      SUBROUTINE NOSIGN(STR,LEN,ISGN)
      IMPLICIT INTEGER (A-Z)
C
C          A ROUTINE TO LOOK FOR A LEADING MINUS SIGN IN A STRING.
C     THE STRING (OF LENGTH "LEN") IS IN "STR".  IF A MINUS IS FOUND,
C     THEN THIS IS REMOVED FROM THE STRING ("STR" AND "LEN" ARE UPDATED
C     ACCORDINGLY).  IF NO SIGN IS FOUND, THEN "ISGN" IS RETURNED
C     +1, IF A MINUS IS FOUND, THEN IT IS RETURNED -1.
C     R. A. WHITESIDE - FEBRUARY 1980
C     M. J. Frisch - January 1984.
C
      Dimension Str(1)
      Save Minus
      Data Minus/1h-/
C
C     GET THE FIRST CHARACTER IN THE STRING, AND SEE IF IT'S A MINUS.
C
      CURSOR=0
      ISGN=1
      CHR=GETCHR(STR,CURSOR)
      IF(CHR.NE.Minus) RETURN
C
C     COPY THE CHARACTERS OF THE STRING DOWN, AND SET VALUE OF ISGN.
C
      ISGN=-1
      LEN=LEN-1
      OCURSR=0
  100 IF(CURSOR.GT.LEN) RETURN
      CHR=GETCHR(STR,CURSOR)
      CALL PUTCHR(CHR,STR,OCURSR)
      GO TO 100
      END
*Deck OTStp
      Subroutine OTStp(IOut,IPrint,IStep,IJob,IUnit,IType,NStat,MxStat,
     $  Stat,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Fetch the coordinates from intermediate point JStep
C     on the OT file and store them in /Mol/.
C
#include "commonmol.inc"
      Integer FilNum
      Logical DoStat,HvStat,DoCart,HvCart,DoCMem,HvCMem,DoGrad,HvGrad
      Dimension Stat(*), V(*), XX(1)
      Save MaxSt1,MaxJb1,NStat1,NGVar1,IJOff1,IType1,IDone,XX,IStep1,
     $  IJob1
      Data MaxSt1/0/,MaxJb1/0/,NStat1/0/,NGVar1/0/,IJOff1/0/,IType1/0/,
     $  IDone/0/,XX/0.0d0/,IStep1/0/,IJob1/0/
 1000 Format(' Recover coordinates from OT file, step',I6,' job',I6,'.')
C
      Write(IOut,1000) IStep, IJob
      IJob1 = 1
      Call RdWrOT(MaxSt1,MaxJb1,NStat1,NGVar1,IJOff1,IType1,IStep1,
     $  IJob1,IDone,.False.,XX,.False.,XX,.False.,XX,.False.,XX,IOut,
     $  IUnit,1,V,MDV)
      IStat = 1
      IV = IStat + NStat1
      Call TstCor(IV,MDV,'OTStp')
      MDV1 = MDV - IV + 1
      Call RdWrOT(MaxSt1,MaxJb1,NStat1,NGVar1,IJOff1,IType1,IStep,IJob,
     $  IDone,.True.,V(IStat),.True.,C,.False.,XX,.False.,XX,IOut,IUnit,
     $  3,V(IV),MDV1)
      IType = IType1
      NStat = Min(MxStat,NStat1)
      Call AMove(NStat,V(IStat),Stat)
      Return
      End
*Deck PrmOut
      Subroutine PrmOut(OK,IV0,STR,LEN,VALUE,INT,FP)
      Implicit Integer(A-Z)
C
C          A ROUTINE TO PRINT OUT A Z-MATRIX PARAMETER DEFINITION
C     CARD.  THE ARGUMENTS ARE
C          OK ..... A LOGICAL VARIABLE INDICATING WHETHER THIS NAME
C              WAS ACTUALLY LOCATED IN THE Z-MATRIX.
C          IV0 .... A LOGICAL VARIABLE INDICATING WHETHER ANY OF
C              THE OPTIONAL FIELDS ON THE DEFINITION CARD ARE PRESTNT.
C          STR .... THE NAME OF THE PARAMETER (HOLLERITH).
C          LEN .... THE NUMBER OF CHARACTERS IN "STR".
C          VALUE .. THE VALUE OF THE PARAMETER.
C          INT .... THE OPTIONAL INTEGER VALUE IN THE CARD.
C          FP ..... THE OPTIONAL FLOATING POINT VALUE.
C
      REAL*8 FP,VALUE,TEN,HUNDRD,Abs,Zero
      COMMON/IO/IN,IOUT,IPUNCH
      LOGICAL IV0,OK
      Dimension Str(1), Line(80)
      Save Zero, Ten, Hundrd, IBlank
      Data Zero/0.0d0/, Ten/10.0d0/, Hundrd/100.0d0/, IBlank/1h /
C
C            INIT. THE OUTPUT LINE CURSOR AND CLEAR THE LINE.
      NCUR=0
      CALL PAD(LINE,NCUR,80,IBlank)
      NCUR=1
C
C                  PUT THE PARAMETER NAME INTO THE LINE.
      CALL PUTB(STR,LEN,LINE,NCUR)
C
C                  TAB, AND PUT VALUE.
      NCUR=20
      IF(VALUE.GE.ZERO) NCUR = NCUR + 1
      IF(Abs(VALUE).LT.HUNDRD) NCUR=NCUR+1
      IF(Abs(VALUE).LT.TEN) NCUR=NCUR+1
      CALL PUTFP(VALUE,5,LINE,NCUR)
C            ERASE THE COMMA
      NCUR=NCUR-1
      CALL PutCSt(' ',LINE,NCUR)
C
C                  PUT THE OPTIONAL INTEGER.
      IF(.NOT.IV0) GO TO 100
      NCUR=30
      If(Int.ge.-1000) then
        If(Int.lt.10) NCUR=NCUR+1
        Call PutI(Int,Line,NCur)
      else
        Call PutCSt('Scan ',Line,NCUr)
        Call PutI(-1001-Int,Line,NCur)
        endIf
      NCUR=NCUR-1
      CALL PutCSt(' ',LINE,NCUR)
C
C                  PUT THE OPTIONAL FP VALUE.
      NCUR=38
      IF(FP.GE.ZERO) NCUR = NCUR + 1
      IF(Abs(FP).LT.HUNDRD) NCUR=NCUR+1
      IF(Abs(FP).LT.TEN) NCUR=NCUR+1
      CALL PUTFP(FP,5,LINE,NCUR)
      NCUR=NCUR-1
      CALL PutCSt(' ',LINE,NCUR)
C
C                  PERHAPS PUT ERROR MESSAGE.
  100 NCUR=48
      IF(.NOT.OK) CALL PutCSt('NOT FOUND IN Z-MATRIX.',LINE,NCUR)
C
C                  WRITE OUR LINE.
      CALL STROUT(IOUT,LINE,NCUR,1)
      RETURN
      END
*Deck PrtITM
      Subroutine PrtITM(IOut,Label,NAtoms,ITRead,MicOpt)
      Implicit Integer(A-Z)
C
C     Print the optimization control arrays.
C
      Character*(*) Label
      Dimension ITRead(NAtoms), MicOpt(NAtoms)
 1000 Format(' ITRead=',20I3)
 1010 Format(' MicOpt=',20I3)
 1020 Format(1X,A,':')
C
      If(Label(1:1).ne.' ') Write(IOut,1020) Label
      Write(IOut,1000) (ITRead(I),I=1,NAtoms)
      Write(IOut,1010) (MicOpt(I),I=1,NAtoms)
      Return
      End
*Deck Rd1Con
      Subroutine Rd1Con(IOut,OKLab,IZ,NAtoms,NMapZA,MapZAt,MaxCon,NBond,
     $  IBond,RBType,JZ)
      Implicit Real*8(A-H,O-Z)
C
C     Read and process one connectivity specification.
C
      Logical OKLab
      Dimension MapZAt(*), NBond(*), IBond(MaxCon,*), RBType(MaxCon,*)
      Save Zero
      Data Zero/0.0d0/
 1000 Format(1X,A,' Value',I6,' in connective exceeds NMapZA=',I5,'.')
 1010 Format(1X,A,' Atom',I6,' in connective exceeds NAtoms=',I5,'.')
C
      If(IZ.lt.1.or.IZ.gt.NMapZA) then
        Write(IOut,1000) 'I', IZ, NMapZA
        Call GauErr('Entry out of range in Rd1Con.')
        endIf
      I = MapZAt(IZ)
      If(I.lt.0.or.I.gt.NAtoms) then
        Write(IOut,1010) 'I', I, NAtoms
        Call GauErr('Entry out of range in Rd1Con.')
        endIf
      JZ = IFFGet(.True.,-1,.False.)
      If(JZ.eq.IZ.and.OKLab) JZ = IFFGet(.True.,-1,.False.)
      If(JZ.gt.0) then
        If(JZ.lt.1.or.JZ.gt.NMapZA) then
          Write(IOut,1000) 'J', JZ, NMapZA
          Call GauErr('Entry out of range in Rd1Con.')
          endIf
        J = MapZAt(JZ)
        If(J.lt.0.or.J.gt.NAtoms) then
          Write(IOut,1010) 'J', J, NAtoms
          Call GauErr('Atom out of range in Rd1Con.')
          endIf
        BondOr = RFFGet(.False.,.True.,Zero,.False.)
        If(I.gt.0.and.J.gt.0) then
          Do 10 K = 1, NBond(I)
            If(IBond(K,I).eq.J) goto 20
   10       Continue
          K = NBond(I) + 1
          NBond(I) = K
          IBond(K,I) = J
   20     RBType(K,I) = BondOr
          If(RBType(K,I).lt.Zero) then
            Call DelBnd(I,J,NBond,MaxCon,IBond,RBType)
            Call DelBnd(J,I,NBond,MaxCon,IBond,RBType)
            endIf
          endIf
        endIf
      Return
      End
*Deck RdBkgC
      Subroutine RdBkgC(IDump,IBkgCh,IBkgUn,MaxRd,LChg,LDBF,ToAng,
     $  NumChg,NumDBF,Chg,DBF)
      Implicit Real*8(A-H,O-Z)
C
C     Read in the external charge specification.  This is given by
C     a combination of point charges and density basis functions,
C     specified by the format:
C
C     X,Y,Z, Charge, Exponent, AZ, Rho, CP
C
C     The Exponent is 0 for a point charges.  The mechanics term with a
C     QM atom of charge Z and VDW "polarizability" P is
C     E = AZ exp (-r/rho) - CP/r**6
C     If AZ (or CP) is negative, it is scaled by the nuclear charge of
C     each QM atom as it is used.
C
      Dimension Chg(LChg,MaxRd), DBF(LDBF,MaxRd), XYZ(3)
      Save IFFFP, IFFEnd, Zero
      Data IFFFP/3hFP /, IFFEnd/3hEND/, Zero/0.0d0/
C
      Call TStamp(1,'Top RdBkgC')
      Call FFSet(IDump)
      Call ILSW(2,14,InpUnt)
      If(IBkgUn.eq.1.or.(IBkgUn.eq.0.and.Mod(InpUnt,2).eq.0)) then
        Fact = GFloat(1) / ToAng
      else
        Fact = GFloat(1)
        endIf
      NumChg = 0
      NumDBF = 0
C
   10 Call FFRead(IEof)
      If(IEOF.ne.0) Goto 900
      Call FFGet(IStr,LenStr,I,XYZ(1),IFound)
      If(IFound.eq.IFFEnd) Goto 900
      If(IFound.ne.IFFFP) Call FFErr(IFFFP,IFound)
      XYZ(1) = Fact*XYZ(1)
      XYZ(2) = Fact*RFFGet(.False.,.False.,Zero,.False.)
      XYZ(3) = Fact*RFFGet(.False.,.False.,Zero,.False.)
      Q = RFFGet(.False.,.False.,Zero,.False.)
      Exp = RFFGet(.True.,.False.,Zero,.False.)
      A = RFFGet(.True.,.False.,Zero,.False.)
      Rho = RFFGet(.True.,.False.,Zero,.False.)
      C = RFFGet(.True.,.False.,Zero,.False.)
      If(IBkgCh.eq.5)
     $  Call GauErr('No standard orientation for charges.')
      If(Exp.eq.Zero) then
        NumChg = NumChg + 1
        Chg(1,NumChg) = XYZ(1)
        Chg(2,NumChg) = XYZ(2)
        Chg(3,NumChg) = XYZ(3)
        Chg(4,NumChg) = Q
        Chg(5,NumChg) = A
        Chg(6,NumChg) = Rho
        Chg(7,NumChg) = C
      else
        NumDBF = NumDBF + 1
        DBF(1,NumDBF) = XYZ(1)
        DBF(2,NumDBF) = XYZ(2)
        DBF(3,NumDBF) = XYZ(3)
        DBF(4,NumDBF) = Q
        DBF(5,NumDBF) = Exp
        DBF(6,NumDBF) = A
        DBF(7,NumDBF) = Rho
        DBF(8,NumDBF) = C
        endIf
      Goto 10
  900 Call TStamp(1,'Bot RdBkgC')
      Return
      End
*Deck RdChrg
      Subroutine RdChrg(IOut,DoRead,DoPrnt,NuOfLa,NFragI,ICkUnt,MaxLay,
     $  IChOni,MulOni,IChFrg,MulFrg)
      Implicit Integer(A-Z)
C
C     Read the charges and multiplicities when the structure is coming
C     from the checkpoint file and do some cleanup and printing
C     regardless of their origin.
C
#include "commonmol.inc"
      Logical DoRead, DoPrnt
      Dimension IChOni(MaxLay,*), MulOni(MaxLay,*), IChFrg(*), MulFrg(*)
C
      If(DoRead) then
        Call ILSW(-ICkUnt,1,IOpCl)
        Call ILSW(-ICkUnt,22,IROHF)
        Call ILSW(1,1,IOpCl)
        Call ILSW(1,22,IROHF)
        endIf
      If(NFragI.gt.1.and.NuOfLa.le.1) then
        ICharg = IChFrg(NFragI+1)
        Multip = MulFrg(NFragI+1)
      else
        NU = Max(1,NuOfLa)
        ICharg = IChOni(1,NU)
        Multip = MulOni(1,NU)
        endIf
      If(DoPrnt) Call PrtChg(IOut,.False.,MaxLay,NuOfLa,NFragI,IChOni,
     $  MulOni,IChFrg,MulFrg)
      If(NFragI.gt.1.and.(ICharg.ne.IChFrg(NFragI+1).or.
     $  Multip.ne.MulFrg(NFragI+1)))
     $  Call GauErr('Inconsistent ONIOM and fragment charges.')
      Return
      End
*Deck RdFzAt
      Subroutine RdFzAt(In,IOut,IPrint,NAtoms,MicOpt,IAn,IAtTyp,IAtRes,
     $  NResid,IResNa,IResNu,Layer,LiAtom,OptAtm,IAtLst,IBlk,INoBlk)
      Implicit Real*8(A-H,O-Z)
C
C     Read an input section giving atoms to optimize and/or freeze.
C     These are modifications from whatever was specified in the
C     usual input.
C
      Logical OptAtm(*)
      Dimension MicOpt(*), IAn(*), IAtTyp(*), IAtRes(*), IResNa(*),
     $  IResNu(*), Layer(*), LiAtom(*), IAtLst(*), IBlk(*), INoBlk(*)
      Save Zero
      Data Zero/0.0d0/
 1000 Format(' Rigid block',I3,' both frozen and optimized.')
 1010 Format(' Atoms in rigid block',I3,
     $  ' are both frozen and optimized.')
 1020 Format('   IBlk=',20I5)
 1030 Format(' INoBlk=',20I5)
 1040 Format(' MicOpt=',20I5)
C
      NAtLst = 0
      Do 10 I = 1, NAtoms
        If(MicOpt(I).ge.-1) then
          NAtLst = NAtLst + 1
          IAtLst(NAtLst) = I
          endIf
   10   Continue
      NBlk = 0
      NNoBlk = 0
      Thresh = Zero
      Call AtmSel(In,IOut,IPrint,2,2,'Blocks',NAtoms,0,IAn,IAtTyp,
     $  IAtRes,NResid,IResNa,IResNu,Layer,LiAtom,NAtLst,IAtLst,NBlk,
     $  IBlk,NNoBlk,INoBlk,Thresh)
      Call LClear(NAtoms,OptAtm)
      Do 20 I = 1, NAtLst
   20   OptAtm(IAtLst(I)) = .True.
      Do 50 IB = 1, NBlk
        Do 30 JB = 1, NNoBlk
          If(INoBlk(JB).eq.IBlk(IB)) then
            Write(IOut,1000) IBlk(IB)
            Write(IOut,1020) (IBlk(K),K=1,NBlk)
            Write(IOut,1030) (INoBlk(K),K=1,NNoBlk)
            Call GauErr('Inconsistent input to RdFzAt.')
            endIf
   30     Continue
        IB1 = IBlk(IB)
        IB2 = -IBlk(IB) - 40
        Do 40 I = 1, NAtoms
   40     OptAtm(I) = OptAtm(I).or.MicOpt(I).eq.IB1.or.MicOpt(I).eq.IB2
   50   Continue
      Do 70 IB = 1, NNoBlk
        IB1 = INoBlk(IB)
        IB2 = -INoBlk(IB) - 40
        Do 60 I = 1, NAtoms
   60     OptAtm(I) = OptAtm(I).and.MicOpt(I).ne.IB1
     $      .and.MicOpt(I).ne.IB2
   70   Continue
      MaxBlk = IArMax(.False.,NAtoms,MicOpt)
      IB1 = IArMin(.False.,NAtoms,MicOpt)
      If(IB1.lt.-40) MaxBlk = Max(MaxBlk,-IB1-40)
      Call IClear(MaxBlk,IBlk)
      Do 80 I = 1, NAtoms
        If(MicOpt(I).lt.-40) then
          IB = -MicOpt(I) - 40
          IVal = -1
        else if(MicOpt(I).gt.1) then
          IB = MicOpt(I)
          IVal = 1
        else
          IB = 0
          IVal = 0
          endIf
        If(IB.ne.0) then
          If(IBlk(IB).eq.0) then
            IBlk(IB) = IVal
          else if(IBlk(IB).ne.IVal) then
            Write(IOut,1010) IB
            Write(IOut,1040) (MicOpt(K),K=1,NAtoms)
            Call GauErr('Inconsistent input.')
            endIf
          endIf
   80   Continue
      Return
      End
*Deck RdStrC
      Subroutine RdStrC(In,IOut,IPunch,IPrint,IDump,IPtMod,IAbort,IRead,
     $  IGoI,IDoPun,InpUnt,ScanDf,PosMul,IStruc,ConsFz,SaveZ,IStepI,
     $  Keep,IZTemp,DoRed,HavRed,RdChk,IGo,NAlRed,NTRedS,IZRed,Value,
     $  IntVcR,FPVecR,NAlSav,IAnSav,CSave,MaxLay,NuOfLa,NFragI,IChOni,
     $  MulOni,Layer,IAnMod,IAtMod,ACMod,LiAtom,Distan,IRwOni,IRwCP,
     $  InpFF,DoAddB,KeepMM,KeepM1,IElDat,IHavMO,MNDNBC,NDNBsC,NBSCom,
     $  IAtWtZ,AtmWtZ,IAtSpZ,AtZEfZ,AtQMoZ,AtGFaZ,IAtFrZ,IAtRsZ,IAtPDZ,
     $  AtZNuZ,FixTet,Symbls,NamCnt,ICkUnt,TitInp,NDChMu,IChFrg,MulFrg,
     $  IChInp,MulInp,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Read one structure for RdGeom.
C
      Parameter (MaxOp=DEFMAXOP,NEl=104,IOZMat=507,IOZRed=665,MxStat=10)
#include "commonmol.inc"
#include "commonz.inc"
      Logical ConsFz, SaveZ, Keep(*), AllChk, HavZMt, DoRed, RdChk,
     $  HavRed, DoAddB, KeepMM, KeepM1, ScanDf, FixTet, PosMul, DoModl,
     $  ZMtMod, DoZMat, RdTiCM, TitInp, TestIn
      Dimension Value(NAlRed,*), FPVecR(NAlRed), CSave(3,NAlSav,*),
     $  IAnSav(NAlSav,*), IChOni(*), MulOni(*), Layer(*), IAnMod(*),
     $  IAtMod(*), ACMod(*), LiAtom(*), Distan(*),IAtWtZ(*), AtmWtZ(*),
     $  IAtSpZ(*), AtZEfZ(*), AtQMoZ(*),AtGFaZ(*), IAtFrZ(*), IAtRsZ(*),
     $  NTRedS(4,*), IZRed(4,NAlRed,*), IntVcR(NAlRed), NBSCom(*),
     $  Symbls(*), NamCnt(*), IAtPDZ(*), IZTemp(MaxNZ,4), V(MDV), JJ(1),
     $  IElDat(*), CelVec(3,3), IChFrg(*), MulFrg(*), AtZNuZ(*),
     $  Stat(MxStat), IChInp(NDChMu,2), MulInp(NDChMu,2), NChMu(2)
      Integer FilNum
      Common /PhyCon/ PhyCon(30)
      Equivalence (ToAng,PhyCon(1))
      Save Zero, One, JJ, NFrag, NACore, NAVal2, NAVal1, NBCore, NBVal2,
     $  NBVal1, IFrStp
      Data Zero/0.0d0/, One/1.0d0/, JJ/0/, NFrag/0/, NACore/0/,
     $  NAVal2/0/, NAVal1/0/, NBCore/0/, NBVal2/0/, NBVal1/0/, IFrStp/0/
 1000 Format(' Unrecognized structure location ',I5,' in RdStrc.')
 1010 Format(' RdStrC is calling ZToCFu for IStruc=',I1,' NZ=',I6,'.')
 1020 Format(' Z-Matrix found in file.')
 1030 Format(' Model A used to generate z-matrix:')
 1040 Format(' Model B used to generate z-matrix:')
 1050 Format(' Generated Z-Matrix:')
 1060 Format(' Symbolic Z-matrix:')
 1080 Format(' Generated',I3,' variables but there are only',I3,
     $       ' degrees of freedom.')
 1100 Format(' Keep=',50L1)
 1150 Format(' No Z-Matrix found in file; cartesian coordinates used.')
 1170 Format(' Redundant internal coordinates found in file.')
 1180 Format(' NFrag=',I3,' on checkpoint file does not match value of',
     $  I3,' for this calculation.')
 1190 Format(' NDNBsC=',I7,' on checkpoint file exceeds limit of',I7,
     $  ' for this calculation.')
C
C     Set defaults
C
      NZ = 0
      SaveZ = .False.
      KeepM1 = KeepMM
      NBB = 0
      Call PutDel(5,ANames,NBB)
      AllChk = IGoI.eq.7
      RdTiCM = IRead.eq.2.or.AllChk
      NCMOn = MaxLay**2
      NCMFr = NFragI + 1
C
C     Read title and charge/multiplicity if requested.
C
      If(.not.RdTiCM) then
        Call RdTitl(In,IOut,IPunch,0,IDoPun,TitInp)
        TestIn = IGoI.eq.2.or.IGoI.eq.3.or.IGoI.eq.6
        Call TstCor(4*NDChMu,MDV,'RdStrC-RdChMu')
        Call RdChMu(IOut,IDump,TestIn,NDChMu,NChMu,IChInp,MulInp,V)
        Call DefChM(IOut,IPrint,PosMul,ICharg,Multip,MaxLay,NuOfLa,
     $    NFragI,NDChMu,NChMu,IChInp,MulInp,IChOni,MulOni,IChFrg,MulFrg)
        endIf
C
C     Read instructions from input, or use the source specified by
C     the iop.
C
      Call GeoNam(IOut,IGoI,IRead,IGo,ICkUnt,V,MDV)
C
C     Determine what to do
C
      DoZMat = IGo.eq.0
      If(DoZMat) Write(IOut,1060)
      ZMtMod = IGo.eq.6
      RdChk = IGo.eq.2.or.IGo.eq.3.or.ZMtMod.or.IGo.eq.7
      DoModl = IGo.eq.4.or.IGo.eq.5
C
C     Print initial stuff.
C
      If(RdChk) then
        If(RdTiCM) then
          Call RdTitl(In,IOut,IPunch,ICkUnt,IDoPun,TitInp)
        else
          Call IMove(NCMOn,IChOni,IChInp)
          Call IMove(NCMOn,MulOni,MulInp)
          Call IMove(NCMFr,IChFrg,IChInp(1,2))
          Call IMove(NCMFr,MulFrg,MulInp(1,2))
          endIf
        Call RWMol(2,-ICkUnt)
        IRwOCh = FilNum(IRwOni,ICkUnt)
        IOCPCh = FilNum(IRwCP,ICkUnt)
        Call RdWrOn(-2,IRwOCh,MaxLay,NAtoms,IChOni,MulOni,Layer,IAnMod,
     $    IAtMod,LiAtom,ACMod,MaxNZ,Distan)
        If(NFragI.gt.1.and.ITqry(IOCPCh).gt.0) then
          Call RdWrFr(-2,IOCPCh,NAtoms,NFrag,IHavMO,NACore,NAVal2,
     $      NAVal1,NBCore,NBVal2,NBVal1,IFrStp,IChFrg,MulFrg,NDNBsC,
     $      NBsCom)
          If(NFrag.ne.NFragI) then
            Write(IOut,1180) NFrag, NFragI
            Call Lnk1E(0)
            endIf
          If(NDNBsC.gt.MNDNBC) then
            Write(IOut,1190) NDNBsC, MNDNBC
            Call Lnk1E(0)
            endIf
          Call RdWrFr(2,IOCPCh,NAtoms,NFrag,IHavMO,NACore,NAVal2,NAVal1,
     $      NBCore,NBVal2,NBVal1,IFrStp,IChFrg,MulFrg,NDNBsC,NBsCom)
          endIf
        If(RdTiCM) then
C         Default ONIOM charge and multiplicity data if necessary.
          Call DfCMOn(IOut,IPrint,.False.,.False.,ICharg,Multip,MaxLay,
     $      NuOfLa,0,IChInp,MulInp,IChOni,MulOni)
          Call DfCMFr(IOut,IPrint,.False.,PosMul,ICharg,Multip,NFragI,
     $      0,IChInp,MulInp,IChFrg,MulFrg)
        else
C         Restore read-in charge and multiplicity
          Call IMove(NCMOn,IChInp,IChOni)
          Call IMove(NCMOn,MulInp,MulOni)
          Call IMove(NCMFr,IChInp(1,2),IChFrg)
          Call IMove(NCMFr,MulInp(1,2),MulFrg)
          endIf
        endIf
      Call RdChrg(IOut,RdTiCM,.True.,NuOfLa,NFragI,ICkUnt,MaxLay,IChOni,
     $  MulOni,IChFrg,MulFrg)
      Call ILSW(2,1,IOpClO)
      IOpCl = 0
      If(IArMax(.True.,NCMOn,MulOni).gt.1) IOpCl = 1
      If(IArMax(.True.,NCMFr,MulFrg).gt.1) IOpCl = 1
      If(IOpClO.eq.0.and.IOpCl.eq.1) Call ILSW(1,1,IOpCl)
C
C     Setup layer info
C
      If(NuOfLa.gt.1) then
        MxRdLa = NuOfLa
      else
        MxRdLa = MaxLay
        endIf
C
C     Determine ICkUnt contents
C
      IOZRCk = FilNum(IOZRed,ICkUnt)
      HavRed = RdChk.and.ITqry(IOZRCk).gt.0
      IOZChk = FilNum(IOZMat,ICkUnt)
      HavZMt = RdChk.and.IGo.ne.3.and.ITqry(IOZChk).gt.0
      IStep = IStepI
      If(IStep.gt.0.and..not.HavZMt) IStep = -(IStep+1)
      HavZMt = HavZMt.and.IStep.ge.0
C
C     Read a symbolic Z-matrix from the input stream.
C
      If(DoZMat) then
        Call ZGet(IOut,.True.,ScanDf,InpUnt,ToAng,JJ,0,.True.,.False.,
     $    MxRdLa,Junk1,Junk2,ConsFz,MaxNZ,IAnZ,IZ,Bl,Alpha,Beta,LBl,
     $    LAlpha,LBeta,NZ,NVarZM,IAtTyZ,AtmChZ,IAtWtZ,AtmWtZ,IAtSpZ,
     $    AtZEfZ,AtQMoZ,AtGFaZ,IAtFrZ,IAtRsZ,IAtPDZ,AtZNuZ,NumTpS,
     $    MaxTpS,IAtTpS,ANames,Values,IntVec,FPVec,Symbls,NamCnt,Layer,
     $    IAnMod,IAtMod,ACMod,LiAtom,Distan,InpFF,MaxLay,MaxRes,NResid,
     $    IResNa,IResNu)
        SaveZ = .True.
C
C     Read redundant internal coordinates from checkpoint file.
C
      else if(HavRed) then
        Write(IOut,1170)
        Call GetPBC(NPDir,CelVec)
        NAtoms = NAtoms + NPDir
        Call FileIO(2,-IOZRCk,InToWP(4),NTRedS(1,IStruc),0)
        NTRed = NTRedS(1,IStruc)
        Call FileIO(2,IOZRCk,InToWP(4*NTRed),IZRed(1,1,IStruc),0)
        Call FileIO(2,IOZRCk,NTRed,Value(1,IStruc),0)
        Call FileIO(2,IOZRCk,InToWP(NTRed),IntVcR,0)
        Call FileIO(2,IOZRCk,NTRed,FPVecR,0)
        Call FileIO(2,IOZRCk,NTRed,Value(1,IStruc+1),0)
        Call FileIO(2,IOZRCk,NTRed,Value(1,IStruc+2),0)
C       KeepM1 = IArMax(.False.,NAtoms,MicOpt).lt.0
        If(IStep.gt.0) then
C         For QuadMac, XX stores some temporary data when calling Star
          MinXX = 3*NTRed + 18*NAtoms + 1
          Call IntStp(IOut,IPrint,IStep,ICkUnt,NTRed,Value(1,IStruc),
     $      MinXX,V)
          Call NewRed(IOut,IPrint,NAtoms,NTRedS(1,IStruc),
     $      NTRedS(2,IStruc),NTRedS(3,IStruc),NTRedS(4,IStruc),
     $      IZRed(1,1,IStruc),Value(1,IStruc),Value(1,IStruc+1),IntVcR,
     $      C,V,MDV)
        else if(IStep.lt.0) then
          Call OTStp(IOut,IPrint,-IStep,1,ICkUnt,IType,NStat,MxStat,
     $      Stat,V,MDV)
          endIf
        Call RedVal(IOut,IPrint,.True.,C,IZRed(1,1,IStruc),
     $    Value(1,IStruc),NTRedS(1,IStruc),NTRedS(2,IStruc),
     $    NTRedS(3,IStruc),NTRedS(4,IStruc),V,MDV)
        Call RdWrZ(7,-ICkUnt,V)
        If(NVarZM.eq.0) Call ZMCord(NAtoms,MicOpt,IAn,IAtTyp,AtChMM,C)
        Call PtZSym(.True.,.False.,IOut,IElDat,MaxNZ,IAnZ,IZ,Bl,Alpha,
     $    Beta,LBl,LAlpha,LBeta,NZ,NVarZM,ANames,Values,IntVec,FPVec,JJ,
     $    JJ(1))
C
C     Read Z-matrix from the checkpoint file.
C
      else if(RdChk.and.HavZMt) then
        Write(IOut,1020)
        Call IMove(NAtoms,IAtWgt,IAtWtZ)
        Call AMove(NAtoms,AtmWgt,AtmWtZ)
        Call GetPBC(NPDir,CelVec)
        NAtoms = NAtoms + NPDir
        Call RdWrZ(2,-ICkUnt,V)
        If(NVarZM.gt.0) then
C         For QuadMac, XX stores some temporary data when calling Star
          MinXX = 3*NVarZM + 18*NAtoms + 1
          If(IStep.ne.0)
     $      Call IntStp(IOut,IPrint,IStep,ICkUnt,NVarZM,Values,MinXX,V)
          endIf
        Do 310 I = 1, MaxZSb
          If(IntVec(I).gt.0) IntVec(I) = 0
  310     FPVec(I) = Zero
        If(ZMtMod.and..not.DoRed) Call NewVar(IOut,InpUnt,ToAng,
     $    NVarZM,ANames,Values,IntVec,FPVec,MaxNZ,IZ,NZ,LBl,LAlpha,
     $    LBeta)
        Call PtZSym(.True.,.False.,IOut,IElDat,MaxNZ,IAnZ,IZ,Bl,Alpha,
     $    Beta,LBl,LAlpha,LBeta,NZ,NVarZM,ANames,Values,IntVec,FPVec,JJ,
     $    JJ(1))
        SaveZ = .True.
C
C     Read cartesian coordinates from the checkpoint file.
C
      else if(RdChk.or.IGo.eq.3) then
        Write(IOut,1150)
        DoAddB = DoAddB.or.ZMtMod
        Call GetPBC(NPDir,CelVec)
        NAtoms = NAtoms + NPDir
        NVarZM = 3*NAtoms
C       For QuadMac, XX stores some temporary data when calling Star
        MinXX = 3*NVarZM + 18*NAtoms + 1
        If(IStep.gt.0) then
          Call IntStp(IOut,IPrint,IStep,ICkUnt,NVarZM,C,MinXX,V)
        else if(IStep.lt.0) then
          Call OTStp(IOut,IPrint,-IStep,1,ICkUnt,IType,NStat,MxStat,
     $      Stat,V,MDV)
          endIf
        Call ZMCord(NAtoms,MicOpt,IAn,IAtTyp,AtmChg,C)
        Call WXYZMM(IOut,PhyCon,0,NuOfLa,NAtoms,IAn,IAtTyp,AtChMM,C,
     $    Layer,IAnMod,IAtMod,LiAtom,ACMod,MaxNZ,Distan,MicOpt,AtmWgt,
     $    IAtWgt,IAtFrg,IAtRes,IAtSpn,AtZEff,AtQMom,AtGFac,NumTpS,
     $    IAtTpS,IAtPDZ,IResNa,IResNu,AtZNuc)
C
C     Read model builder input from the input stream.
C
      else if(DoModl) then
        IModel = IGo - 4
        If(IModel.eq.0) Write(IOut,1030)
        If(IModel.eq.1) Write(IOut,1040)
        IPrinL = IPtMod + IPrint + IDump + 1
        IFromA = 1
        IEV = IFromA + InToWP(MaxNZ)
        IEVC = IEV + InToWP(MaxNZ+1)
        ICH = IEVC + InToWP(MaxNZ+1)
        IG = ICH + InToWP(MaxNZ+1)
        IEnd = IG + InToWP(6*MaxNZ) - 1
        Call TstCor(IEnd,MDV,'RdStrC-MBld2')
        Call MBld2(IModel,IPrinL,NEl,IElDat,InpUnt,NAtoms,IAn,IEof,Keep,
     $    V(IFromA),V(IEV),V(IEVC),V(ICH),V(IG))
        If(IPrinL.ge.3) Write(IOut,1100) (Keep(I),I=1,NAtoms)
        IOldNo = 1
        ITest = IOldNo + InToWP(MaxNZ)
        IFrom = ITest + InToWP(MaxNZ+1)
        INewNo = IFrom + InToWP(MaxNZ+1)
        IBDDef = INewNo + InToWP(MaxNZ+1)
        IBndLt = IBDDef + InToWP(MaxNZ+1)
        IBDAng = IBndLt + MaxNZ
        IRotAn = IBDAng + MaxNZ
        IA = IRotAn + MaxNZ
        IB = IA + MaxNZ
        ID = IB + MaxNZ
        ITE = ID + MaxNZ
        ICZ = ITE + InToWP(MaxNZ+1)
        IAScr = ICZ + 3*MaxNZ
        IEnd = IAScr + 5*MaxNZ - 1
        Call TstCor(IEnd,MDV,'RdStrC-MBld3')
        Call MBld3(IModel,IPrinL,MaxNZ,IElDat,InpUnt,ToAng,NAtoms,IAn,C,
     $    Keep,IZTemp,Bl,Alpha,Beta,IAnZ,IEof,V(IOldNo),V(ITest),
     $    V(IFrom),V(INewNo),V(IBDDef),V(IBndLt),V(IBDAng),V(IRotAn),
     $    V(IA),V(IB),V(ID),V(ITE),V(ICZ),V(IAScr))
C
C       Allocate arrays for generation of a new z-matrix:  NumBl,
C       NumA,NumB,BLUsed,AUsed,BUsed,ConstA,ConstB,AtmChg(all NAtoms):
        INBL = 1
        INA  = INBL + InToWP(NAtoms)
        INB  = INA  + InToWP(NAtoms)
        IBLU = INB  + InToWP(NAtoms)
        IAU  = IBLU + InToWP(NAtoms)
        IBU  = IAU  + InToWP(NAtoms)
        IConA = IBU + InToWP(NAtoms)
        IConB = IConA + InToWP(NAtoms)
        IAtmC = IConB + InToWP(NAtoms)
        MaxAP3 = NAtoms + 3
C       NPerm(MaxAP3,MaxOp), Scr(MaxNZ,10):
        INPerm = IAtmC + NAtoms
        IScr = INPerm + InToWP(MaxAP3*MaxOp)
        INToO = IScr + 10*Max(MaxNZ,100,NAtoms)
        IOToN = INToO + InToWP(MaxNZ)
        IEnd = IOToN + InToWP(MaxNZ) - 1
        Call TstCor(INBl+IEnd-1,MDV,'RdStrC')
        NZIn = NAtoms
        Call IMove(NZIn,IZTemp(1,1),IZ(1,1))
        Call IMove(NZIn,IZTemp(1,2),IZ(1,2))
        Call IMove(NZIn,IZTemp(1,3),IZ(1,3))
        Call IMove(NZIn,IZTemp(1,4),IZ(1,4))
        Fact = One / ToAng
        Call AScale(NZIn,Fact,Bl,Bl)
        Call MakZmt(IOut,IPrinL,.True.,MaxNZ,MaxAp3,MaxOp,NZIn,ICharg,
     $    Multip,NFree,IAnZ,IZ,BL,Alpha,Beta,V(INPerm),C,V(IAtmC),ToAng,
     $    Keep,LBl,LAlpha,LBeta,NVarZM,ANames,Values,V(INBL),V(INA),
     $    V(INB),V(IBLU),V(IAU),V(IBU),V(IScr),V(IConA),V(IConB),
     $    V(IScr+MaxNZ),.False.,.False.,1,JJ,JJ,Zero,Zero,JJ,V(INToO),
     $    V(IOToN),V(IEnd+1),MDV-IEnd)
        NZ = NZIn
        Call IClear(MaxZSb,IntVec)
        Call AClear(MaxZSb,FPVec)
        Write(IOut,1050)
        Call PtZSym(.True.,.False.,IOut,IElDat,MaxNZ,IAnZ,IZ,Bl,Alpha,
     $    Beta,LBl,LAlpha,LBeta,NZ,NVarZM,ANames,Values,IntVec,FPVec,JJ,
     $    JJ(1))
        If(NVarZM.ne.NFree) Write(IOut,1080) NVarZM, NFree
        If(NVarZM.ne.NFree.and.IAbort.eq.2) Call Lnk1E(0)
        Call VFloat(NZ,IAnZ,AtZNuZ)
        SaveZ = .True.
      else
        Write(IOut,1000) IGo
        Call Lnk1E(0)
        endIf
C
C     If necessary, compute Cartesian coordinates.
C
      If(NZ.gt.0.and.(NAtoms.eq.0.or.IGo.eq.0.or.IGo.eq.1.or.IGo.eq.4
     $  .or.IGo.eq.5.or.IGo.eq.6)) then
        If(IPrint.ge.1) Write(IOut,1010) IStruc, NZ
        Call ZToCFu(IOut,IPrint,FixTet,.True.,IAtWtz,AtmWtZ,IAtSpZ,
     $    AtZEfZ,AtQMoZ,AtGFaZ,IAtFrZ,IAtRsZ,IAtPDZ,AtZNuZ,NAtoms,IAn,
     $    IAtTyp,AtmChg,AtChMM,C,IAtWgt,AtmWgt,IAtSpn,AtZEff,AtQMom,
     $    AtGFac,IAtFrg,IAtRes,IAtPDB,AtZNuc,V,MDV)
        endIf
      If(NAtoms.eq.0) Call GauEr2(0,
     $  'There are no atoms in this input structure !',' ',' ',
     $  'molecule specification')
      Call IMove(NAtoms,IAn,IAnSav(1,IStruc))
      Call AMove(3*NAtoms,C,CSave(1,1,IStruc))
      Return
      End
*Deck RdTitl
      Subroutine RdTitl(In,IOut,IPunch,IUChk,IPch,TitInp)
      Implicit Integer(A-Z)
C
C     This subroutine reads and prints the job title.
C
C     IUChk  ... if non-zero read from chk file on this unit instead of
C                from input
C     IPch   ... Title printing flag
C     TitInp ... Set to true when a title read from input has been saved.
C
      Logical TitInp
      Character*(*) Dash
      Parameter (IOLbl=502,Dash='-',MaxLin=1024)
#include "commonlab.inc"
      Dimension IC(MaxLin), IC1(MaxLin)
 1020 Format(1X,128A1)
C
      LenTit = InToWP(MaxTit)
      OffTit = InToWP(4)
      LTitC = MaxTit*NChrPW(0)
C
C     Read title from IUChk
C
      If(IUChk.ne.0) then
        Call FileIO(2,-FilNum(IOLbl,IUChk),LenTit,ITitle,OffTit)
        Call StrPak(.True.,.True.,ITitle,LTitC,Len)
C
C     Read title from input
C
      else
        Len = 0
        Call CBlank(MaxTit,ITitle)
  100   Call StrIn(In,IC,MaxLin,IEOF)
        If(IEOF.eq.0) then
          LOld = Len
          L1 = 0
          Call PakStr(IC,MaxLin,IC1,L1)
          LenPak = Min(MaxLin,LTitC-LOld)
          Call PakStr(IC,LenPak,ITitle,Len)
          If(Len.lt.LTitC) Call PutCSt(' ',ITitle,Len)
          If(L1.gt.0) goto 100
          endIf
        Len = Len - 2
        endIf
      If(.not.TitInp) Call FileIO(1,-IOLbl,LenTit,ITitle,OffTit)
      TitInp = TitInp.or.IUChk.eq.0
      LEnd = Min(Len,70)
      Write(IOut,1020) (Dash,I=1,LEnd)
      Call SetLnL(70)
*
* ANT BEGIN modification
*
*        READ NAME FROM GAussian
*
       open(unit=10,file='name',status='unknown')
       Call StrOut(10,ITitle,Len,1)
       close(10)
*
* ANT END modification
*
      Call StrOut(IOut,ITitle,Len,1)
      If(IPch.ne.0) Call StrOut(IPunch,ITitle,Len,2)
      Call SetLnl(256)
      Write(IOut,1020) (Dash,I=1,LEnd)
      Return
      End
*Deck RdZVar
      Logical Function RdZVar(RdInpt,ScanDf,BB,NBB,Name,MaxNam,LenNam,
     $  GotVal,Value,GotInt,IntVal,GotFP,FPVal)
      Implicit Integer(A-Z)
C
C     Read one z-matrix variable card, either from the input stream
C     or from BB, depending on RdInpt.  The value of the function is
C     true if a card was read, false for end of file.
C
      Parameter (MaxStr=80)
      Logical RdInpt, GotVal, GotInt, GotFP, StrEqC, ScanDf
      Real*8 Value, FP, FPVal
      Dimension BB(1), Name(MaxNam), IStr(MaxStr), IJunk(10)
      Save IFFEnd,IFFStr,IFFFP,IFFInt,IhF,IhT,IhA,IHStar,IhD,IhS,IHH
      Data IFFEnd/3hEND/, IFFStr/3hSTR/, IFFFP/3hFP /, IFFInt/3hINT/
      Data IHF/1hF/, IHT/1hT/, IHA/1hA/, IHStar/1h*/, IhD/1hD/
      Data IHS/1hS/, IHH/1hH/
C
C     Read the line.
C
      If(RdInpt) then
        Call FFRead(IEof)
      else
        Call FFRMem(BB,NBB,IEof)
        endIf
      If(IEOF.ne.0) then
        RdZVar = .False.
        Return
        endIf
C
C     Fetch the name.  A blank line means EOF.
C
      Call CBlank(MaxNam,Name)
      Call FFGet(Name,LenNam,I,FP,IF1)
      JunkC = 0
      Call PutCSt('Constant',IJunk,JunkC)
      If(IF1.eq.IFFEnd.or.
     $  (IF1.eq.IFFStr.and.StrEqC(Name,IJunk,8))) then
        RdZVar = .False.
        Return
        endIf
      If(IF1.ne.IFFStr) Call FFErr(IFFStr,IF1)
C
C     Fetch the value.  A "*" means don't change the value, and "F" or
C     "A" mean freeze or activate the variable (i.e., apply to the
C     integer field.  A "D" means differentiate (and also activate).
C
      GotInt = .False.
      GotVal = .False.
      Call FFGet(IStr,LenStr,I,FP,IF1)
      If(IF1.eq.IFFFP) then
        GotVal = .True.
        Value = FP
      else if(IF1.eq.IFFStr.and.LenStr.eq.1.and.(StrEqC(IStr,IhF,1).or.
     $  StrEqC(IStr,IhT,1).or.StrEqC(IStr,IhA,1).or.StrEqC(IStr,IhD,1)
     $  .or.StrEqC(IStr,IhS,1).or.StrEqC(IStr,IhH,1))) then
        GotInt = .True.
      else if(IF1.ne.IFFStr.or.LenStr.ne.1.or.
     $  .not.StrEqC(IStr,IHStar,1)) then
        Call FFErr(IFFFP,IF1)
        endIf
C
C     Fetch the optional integer and FP numbers.  Map F,A,D to -1,0,3
C     and handle the case where a letter was read instead of a value.
C
      If(GotInt) then
        IF1 = IFFStr
      else
        Call FFGet(IStr,LenStr,I,FP,IF1)
        GotInt = IF1.ne.IFFEnd
        endIf
      If(GotInt) then
        If(IF1.eq.IFFInt) then
          IntVal = I
          If(ScanDf) IntVal = -1001 - I
        else if(IF1.eq.IFFStr.and.LenStr.eq.1.and.StrEqC(IStr,IhF,1))
     $  then
          IntVal = -1
        else if(IF1.eq.IFFStr.and.LenStr.eq.1.and.
     $  (StrEqC(IStr,IhT,1).or.StrEqC(IStr,IhA,1))) then
          IntVal = 0
        else if(IF1.eq.IFFStr.and.LenStr.eq.1.and.StrEqC(IStr,IhH,1))
     $  then
          IntVal = 1
        else if(IF1.eq.IFFStr.and.LenStr.eq.1.and.StrEqC(IStr,IhD,1))
     $  then
          IntVal = 3
        else if(IF1.eq.IFFStr.and.LenStr.eq.1.and.StrEqC(IStr,IhS,1))
     $  then
          Call FFGet(IStr,LenStr,I,FP,IF1)
          IntVal = -1001
          If(IF1.eq.IFFInt) then
            IntVal = -1001 - I
          else
            Call FFErr(IFFInt,IF1)
            endIf
        else
          Call FFErr(IFFInt,IF1)
          endIf
        Call FFGet(Junk,LJunk,Junk,FP,IF1)
        GotFP = IF1.ne.IFFEnd
        If(GotFP) then
          If(IF1.ne.IFFFP) Call FFErr(IFFFP,IF1)
          FPVal = FP
          endIf
      else
        GotFP = .False.
        endIf
      RdZVar = .True.
      Return
      End
*Deck ReBldZ
      Subroutine ReBldZ(IOut,UseSym,MaxNZ,MaxAP3,Cut,NPerm,NAtoms,NOp,
     $  IAn,C,IZ)
      Implicit Real*8(A-H,O-Z)
C
C     Rebuild the connectivity matrix IZ.  Only the bond pairs in IZ(,1)
C     are retained.  Bond angles are chosen to involve the new bond
C     and the bond used to define atom IZ(,1).  In selecting dihedral
C     angles for atom I bonded to J making an angle with K, the first
C     choice is to use the lowest numbered atom bonded to K to form the
C     dihedral angle.  If there is no previously defined atom bonded
C     to K, the lowest numbered atom bonded to J is used.  Preference
C     is given to an atom in a plane of symmetry.  Map is a scratch
C     array NAtoms long.
C
      Dimension IZ(MaxNZ,4), IAn(1), C(3,1), NPerm(MaxAP3,1), CutAng(3)
      Logical UseSym, CkPlan, AngOK, CkSym
      Save CutAng
      Data CutAng/0.1d0,1.d-6,0.0d0/
 1000 Format(' Unable to find non-linear bond angle for atom',I4,'.')
C
C     Choose atom 2 to be unique among the atoms bonded to atom 1,
C     if possible.
C
      If(NAtoms.lt.3) goto 999
      CkSym = UseSym.and.NOp.gt.1
      If(NAtoms.eq.3) goto 35
      Do 20 I = 2, NAtoms
        If(IZ(I,1).eq.1) then
          Do 15 IP = 1, NOp
            If(NPerm(I,IP).ne.I.and.IZ(NPerm(I,IP),1).eq.1) goto 20
   15       Continue
          Call SwpAtm(MaxNZ,MaxAP3,NAtoms,NOp,I,2,IAn,C,IZ,NPerm)
          Goto 25
          endIf
   20   Continue
C
C     Choose atom 3 to be one in a plane of symmetry, if possible.
C
   25 If(.not.CkSym.or.((IZ(3,1).eq.1.or.IZ(3,1).eq.2).and.
     $  CkPlan(NAtoms,3,2,1,Cut,IAn,C))) goto 35
      Do 30 I = 4, NAtoms
        If(IZ(I,1).eq.1.or.IZ(I,1).eq.2) then
          If(CkPlan(NAtoms,I,2,1,Cut,IAn,C)) then
            Call SwpAtm(MaxNZ,MaxAP3,NAtoms,NOp,I,3,IAn,C,IZ,NPerm)
            Goto 35
            endIf
          endIf
   30   Continue
C
C     Make sure that the remaining atoms are bonded to a previous one.
C
   35 Do 40 I = 3, NAtoms
        If(IZ(I,1).ge.1.and.IZ(I,1).lt.I) goto 40
        I1 = I + 1
        Do 37 J = I1, NAtoms
          If(IZ(J,1).lt.I) then
            Call SwpAtm(MaxNZ,MaxAP3,NAtoms,NOp,J,I,IAn,C,IZ,NPerm)
            Goto 40
            endIf
   37     Continue
   40   Continue
C
C     Now choose new bond angles.  Try to avoid bond angles too close
C     to 0 or 180 degrees.
C
      Do 140 I = 1, NAtoms
        J = IZ(I,1)
        K = 0
        If(I.le.2) goto 46
        Do 44 ICut = 1, 3
          If(I.gt.2.and.J.ne.1.and.J.ne.0) K = IZ(J,1)
          If(I.gt.2.and.J.eq.1) K = 2
          If(ICut.eq.3.or.AngOK(I,J,K,C,CutAng(ICut))) goto 46
          If(I.eq.3) goto 44
          LimKK = I - 1
          Do 42 KK = J, LimKK
            If(IZ(KK,1).ne.J) goto 42
            K = KK
            If(AngOK(I,J,K,C,CutAng(ICut))) goto 46
   42       Continue
   44     Continue
        Write(IOut,1000) I
        Call Lnk1E(0)
C
C       Now do dihedral angles.  First try atoms L bonded to K for
C       which L-K-J is a plane of symmetry (so there will be another
C       atom I' equivalent to I).
C
   46   L = 0
        If(I.le.3) goto 130
        If(K.ne.0.and.K.ne.1) L = IZ(K,1)
        If(CkSym.and.L.ne.J.and.L.ne.0) then
          If(CkPlan(NAtoms,J,K,L,Cut,IAn,C)) goto 130
          endIf
        LSt = K + 1
        LEnd = I - 1
        If(.not.CkSym.or.LSt.gt.LEnd) goto 60
        Do 50 L = LSt, LEnd
          If(L.ne.J.and.IZ(L,1).eq.K) then
            If(CkPlan(NAtoms,J,K,L,Cut,IAn,C)) goto 130
            endIf
   50     Continue
C
C       Next, try atoms L bonded to J such that K-J-L is a plane
C       of symmetry.
C
   60   LSt = J + 1
        L = 0
        If(J.ne.0.and.J.ne.1) L = IZ(J,1)
        If(CkSym.and.L.ne.K.and.L.ne.0) then
          If(CkPlan(NAtoms,J,K,L,Cut,IAn,C)) goto 130
          endIf
        If(.not.CkSym.or.LSt.gt.LEnd) goto 80
        Do 70 L = LSt, LEnd
          If(L.ne.K.and.IZ(L,1).eq.J) then
            If(CkPlan(NAtoms,J,K,L,Cut,IAn,C)) goto 130
            endIf
   70     Continue
C
C       Now give up on symmetry and find any atom L bonded to J and
C       use the angle between K-J-L and K-J-I
C
   80   L = 0
        If(J.ne.0.and.J.ne.1) L = IZ(J,1)
        If(L.ne.K.and.L.ne.0) goto 130
        LSt = J + 1
        LEnd = I - 1
        If(LSt.gt.LEnd) goto 100
        Do 90 L = LSt, LEnd
          If(L.ne.K.and.IZ(L,1).eq.J) goto 130
   90     Continue
C
C       Last try.  Find any atom bonded to K and use the angle between
C       L-K-J and K-J-I.
C
  100   L = 0
        If(K.ne.0.and.K.ne.1) L = IZ(K,1)
        If(L.ne.J.and.L.ne.0) goto 130
        LSt = K + 1
        If(LSt.gt.LEnd) goto 120
        Do 110 L = LSt, LEnd
          If(L.ne.J.and.IZ(L,1).eq.K) goto 130
  110     Continue
C
C       Finally, keep whatever atom was used previously to define the
C       dihedral angle.
C
  120   L = IZ(I,3)
C
  130   IZ(I,2) = K
        IZ(I,3) = L
        IZ(I,4) = 0
  140   Continue
  999 Return
      End
*Deck RedIso
      Subroutine RedIso(In,IOut,IDump,IPrint,Isotop,IThTem,IThPre,
     $  IThScl,IRwGen,BB)
      Implicit Real*8(A-H,O-Z)
C
C     Handle setting up isotope information in /Mol/.
C
#include "commonmol.inc"
#include "commonmol2.inc"
      Parameter (LenTP=3,IndTP=48)
      Integer BB(*)
      Dimension TP(LenTP), XX(1)
      Character*3 NamFil
      Common /MUnit/  IUnit(20)
      Save Zero, XX
      Data Zero/0.0d0/, XX/0.0d0/
 1010 Format(' Cannot read isotopes from ',A,
     $  ' file:  /Mol/ is not present.')
 1020 Format(' Number of atoms on ',A,' file =',I6,
     $  ' does not match current value of',I6,'.')
 1030 Format(' Atomic number for atom',I6,' on ',A3,' file =',I4,
     $  ' does not match current value of',I4,'.')
 1040 Format(' Using temperature',F9.3,' Kelvin and pressure',F9.3,
     $  ' atmospheres.')
 1050 Format(' Using temperature',F9.3,' Kelvin, pressure',F9.3,
     $  ' atmospheres, and frequency scale factor',F9.6,'.')
 1060 Format(' IAtWgt=',10I12)
 1070 Format(' AtmWgt=',10F12.7)
 1090 Format(20X,'Isotopes and Nuclear Properties:',/,
     $  ' (Nuclear quadrupole moments (NQMom) in fm**2,',
     $ ' nuclear magnetic moments (NMagM)',/,'  in nuclear magnetons)')
 1100 Format(' NucSpn=',10I12)
 1110 Format(' AtZEff=',10F12.7)
 1120 Format(' NQMom= ',10F12.7)
 1130 Format(' NMagM= ',10F12.7)
 1140 Format(/,'  Atom  ',I8,9I12)
 1150 Format(' AtZNuc=',10F12.7)
C
      NBB = 0
      Call AClear(LenTP,TP)
      If(IThTem.ge.0) then
        TempDf = GFloat(IThTem) / GFloat(1000)
      else
        TempDf = -GFloat(IThTem) / GFloat(1000000)
        endIf
      If(IThPre.ge.0) then
        PresDf = GFloat(IThPre) / GFloat(1000)
      else
        PresDf = -GFloat(IThPre) / GFloat(1000000)
        endIf
      ScalDf = GFloat(IThScl) / GFloat(1000000)
      TP(1) = TempDf
      TP(2) = PresDf
      TP(3) = ScalDf
      If(Isotop.eq.1.or.Isotop.eq.2) then
        If(Isotop.eq.2)
     $    Call RdTPM(IDump,TempDf,PresDF,ScalDf,TP(1),TP(2),TP(3))
        IOpt = 4 - Isotop
        Call LodMas(IOpt,In,IOut,IAn,NAtoms,.True.,.True.,IAtWgt,AtmWgt,
     $    IAtSpn,AtQMom,AtGFac,BB,NBB)
        If(Isotop.eq.2) Call FFRead(IEOF)
      else if(Isotop.eq.3.or.Isotop.eq.4) then
        If(Isotop.eq.3) then
          MolUnt = 0
          NamFil = 'rwf'
        else
          MolUnt = -IUnit(9)
          NamFil = 'chk'
          endIf
        Call RWMol(-2,MolUnt)
        If(NAtom2.le.0) then
          Write(IOut,1010) NamFil
          Call Lnk1E(0)
          endIf
        If(NAtom2.ne.NAtoms) then
          Write(IOut,1020) NamFil, NAtom2, NAtoms
          Call Lnk1E(0)
          endIf
        Do 10 I = 1, NAtoms
          If(IAn(I).ne.IAn2(I)) then
            Write(IOut,1030) I, NamFil, IAn2(I), IAn(I)
            Call Lnk1E(0)
            endIf
          IAtWgt(I) = IAtWt2(I)
          Call PutIso(.False.,1,I,IAn(I),IAtWgt(I),BB,NBB)
          AtmWgt(I) = AtmWt2(I)
          IAtSpn(I) = IAtSp2(I)
          AtZEff(I) = AtZEf2(I)
          AtZNuc(I) = AtZNu2(I)
          AtQMom(I) = AtQMo2(I)
   10     AtGFac(I) = AtGFa2(I)
        NBB = NBB - 1
      else if(Isotop.eq.5) then
        Do 20 I = 1, NAtoms
          Call PutIso(.False.,1,I,IAn(I),IAtWgt(I),BB,NBB)
   20     Continue
        NBB = NBB - 1
      else
        Call GauErr('Illegal value for Isotop in RedIso.')
        endIf
      Call FileIO(1,-IRwGen,LenTP,TP,IndTP)
      If(IPrint.gt.0.or.IPSave(0).eq.0) then
        Call GetTPM(In,IOut,0,.False.,.False.,.True.,0,NAtoms,IAn,
     $    TP(1),TP(2),TP(3),XX)
        If(Isotop.eq.2) then
C         Store default values if 0.0, etc. was read from input.
          Call FileIO(1,-IRwGen,LenTP,TP,IndTP)
          If(TP(3).eq.Zero) then
            Write(IOut,1040) TP(1), TP(2)
          else
            Write(IOut,1050) TP
            endIf
          endIf
        Write(IOut,1090)
C       Call StrOut(IOut,BB,NBB,1)
        Do 100 ISt = 1, NAtoms, 10
          IEnd = Min(ISt+9,NAtoms)
          Write(IOut,1140) (I,I=ISt,IEnd)
          Write(IOut,1060) (IAtWgt(I),I=ISt,IEnd)
          Write(IOut,1070) (AtmWgt(I),I=ISt,IEnd)
          Write(IOut,1100) (IAtSpn(I),I=ISt,IEnd)
          Write(IOut,1110) (AtZEff(I),I=ISt,IEnd)
          Write(IOut,1120) (AtQMom(I),I=ISt,IEnd)
          Write(IOut,1130) (AtGFac(I),I=ISt,IEnd)
          Write(IOut,1150) (AtZNuc(I),I=ISt,IEnd)
  100     Continue
        endIf
      Return
      End
*Deck RedTS
      Subroutine RedTS(IOut,IPrint,NStpTS,IAprBG,NAtoms,NTrRot,NAlRed,
     $  NStruc,NAlSav,IAnSav,CSave,IZRed,NTRedS,Value,WeiFac,Cnstr,
     $  IntVec,Doreac,DoQST,DoPath,ANames,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Return the structure that will be used by l103 for either QST or
C     Path Relaxation
C
      Logical DoReac,DoQST,DoPath
      Dimension IAnSav(NAlSav,*), CSave(3,NAlSav,*), IZRed(4,NAlRed,*),
     $  NTRedS(4,*), Value(NAlRed,*), Cnstr(*), IntVec(NAlRed,*),
     $  V(MDV), ANames(*)
C
      If(DoQST.and.NStruc.ge.2) II = 3
      If(DoPath) II = NStruc + 1
C
C     First, form the union of the reactant and product coordinate systems,
C     and convert reactants and product coordinate values to the union.
C
      NTot = NTRedS(1,1) + NTRedS(1,2)
      IScr1 = 1
      IScr2 = IScr1 + InToWP(NTot)
      IZTmp = IScr2 + InToWP(NTot)
      IEnd = IZTmp + InToWP(4*NTot) - 1
      Call TstCor(IEnd,MDV,'RedTS-Union2')
      Call Union2(IOut,NAtoms,NAlRed,IZRed(1,1,1),NTRedS(1,1),
     $  NTRedS(2,1),NTRedS(3,1),NTRedS(4,1),IntVec(1,1),IZRed(1,1,2),
     $  NTRedS(1,2),NTRedS(2,2),NTRedS(3,2),NTRedS(4,2),IntVec(1,2),
     $  IZRed(1,1,II),NTRedS(1,II),NTRedS(2,II),NTRedS(3,II),
     $  NTRedS(4,II),IntVec(1,II),V(IScr1),V(IScr2),V(IZTmp))
      If(II.ne.1) Call IMove(NTRedS(1,II),IntVec(1,II),IntVec)
      Call RedVal(IOut,IPrint,.True.,CSave(1,1,1),IZRed(1,1,II),
     $  Value(1,1),NTRedS(1,II),NTRedS(2,II),NTRedS(3,II),NTRedS(4,II),
     $  V,MDV)
      Call RedVal(IOut,IPrint,.True.,CSave(1,1,2),IZRed(1,1,II),
     $  Value(1,2),NTRedS(1,II),NTRedS(2,II),NTRedS(3,II),NTRedS(4,II),
     $  V,MDV)
      If(NStruc.ge.3) Call RedVal(IOut,IPrint,.True.,CSave(1,1,3),
     $  IZRed(1,1,II),Value(1,3),NTRedS(1,II),NTRedS(2,II),NTRedS(3,II),
     $  NTRedS(4,II),V,MDV)
      If(NStruc.eq.3.and.DoQST) then
C       Nothing else to do
      else if(DoPath.and.DoReac) then
        NAt3 = 3*Natoms
        Call IMove(NAtoms,IAnSav(1,II-1),IAnSav(1,II))
        Call AMove(NAt3,CSave(1,1,1),CSave(1,1,II))
        Call Amove(NalRed,Value(1,1),Value(1,II))
      else
C
C       If only two structures were input or if the reactant structure
C       is not going to be optimized, make a guess at the transition
C       structure or at the first point along the path.
C
        If(II.ge.3) then
          LenMin = 20
          NAlloc = Max(3*NAtoms,NTRedS(1,II))
          NAlloc = Max(NAlloc,LenMin)
          IDQ = 1
          IDQI = IDQ + NTRedS(1,II)
          IXStr = IDQI + NTRedS(1,II)
          IXVal = IXStr + 3*NAtoms
          IB = IXVal + NTRedS(1,II)
          IIB = IB + 12*NTRedS(1,II)+9*NAtoms
          IG = IIB + 4*NTRedS(1,II)
          ISq = IG + NAlloc**2
          ITri = ISq + NAlloc**2
          IOne4D = ITri + (NAlloc*(NAlloc+1))/2
          IXStr = IOne4D + 4*NAlloc
          IEnd = IXStr + NAlloc - 1
          Call TstCor(IEnd,MDV,'RdPth3')
          JJ = II - 1
          Call GTrans(IOut,IPrint,NStpTS,IAprBG,NAtoms,NTrRot,
     $      CSave(1,1,JJ),CSave(1,1,1),CSave(1,1,II),V(IXVal),
     $      IZRed(1,1,II),NTRedS(1,II),NTRedS(2,II),NTRedS(3,II),
     $      NTRedS(4,II),WeiFac,2,Value(1,1),Value(1,JJ),V(IDQ),V(IDQI),
     $      V(IXStr),V(IB),V(IIB),V(IG),V(ISq),V(ITri),V(IOne4D),Cnstr,
     $      IntVec,ANames,V(IEnd+1),MDV-IEnd)
          Call IMove(NAtoms,IAnSav(1,JJ),IAnSav(1,II))
          endIf
        Call RedVal(IOut,IPrint,.True.,CSave(1,1,3),IZRed(1,1,II),
     $    Value(1,3),NTRedS(1,II),NTRedS(2,II),NTRedS(3,II),
     $    NTRedS(4,II),V,MDV)
        endIf
      Return
      End
*Deck Redund
      Subroutine Redund(MaxNZ,MaxAP3,NAtoms,NOp,Cut,IZ,NPerm,Be,BlUsed,
     $  AUsed,BUsed,NumBl,NumA,NumB)
      Implicit Real*8(A-H,O-Z)
C
C     Determine which bond distances, bond angles and dihedral angles
C     are unique and which are constrained to be equal by symmetry.
C     The connectivity IZ and the permutation matrix NPerm are input.
C     NPerm(I,J) is the atom into which atom I is mapped by symmetry
C     operation J.  The output is loaded into
C     the remaining arrays.  BlUsed(I) is .true. if bond length I
C     is a variable which is to be kept (and is not equivalent to a
C     previous bond length).  NumBl(I) is the number of the bond length
C     to use for bond length I.  The other arrays have the same function
C     for bond angles and dihedral angles.  Note that NumB(I) is
C     negative if dihedral angle I is constrained to be the negative of
C     some other dihedral angle.
C
      Dimension IZ(MaxNZ,4), NPerm(MaxAP3,1), NumBl(1), NumA(1),
     $  NumB(1), Be(1)
      Logical BLUsed(1), AUsed(1), BUsed(1)
      Save Zero
      Data Zero/0.0d0/
C
      If(NAtoms.lt.1) goto 999
      Do 10 I = 1, NAtoms
          NumBl(I) = I
          NumA(I) = I
          NumB(I) = I
          BlUsed(I) = I.gt.1
          AUsed(I) = I.gt.2
   10     BUsed(I) = I.gt.3
      If(NAtoms.lt.3.or.NOp.lt.1) goto 999
      Do 70 IAtom = 3, NAtoms
        LimJAt = IAtom - 1
        IBond = IZ(IAtom,1)
        IAngle = IZ(IAtom,2)
        IDih = IZ(IAtom,3)
        Do 70 JAtom = 2, LimJAt
          JBond = IZ(JAtom,1)
          JAngle = IZ(JAtom,2)
          JDih = IZ(JAtom,3)
          Do 70 IOper = 1, NOp
            IMap = NPerm(IAtom,IOper)
            IBMap = NPerm(IBond,IOper)
            IAMap = NPerm(IAngle,IOper)
            IDMap = NPerm(IDih,IOper)
            If(.not.(IMap.eq.JAtom.and.IBMap.eq.JBond).and.
     $         .not.(IMap.eq.JBond.and.IBMap.eq.JAtom)) goto 20
              If(BLUsed(IAtom)) NumBl(IAtom) = JAtom
              BLUsed(IAtom) = .false.
   20       If(JAtom.eq.2) goto 70
            If(.not.(IMap.eq.JAtom.and.IBMap.eq.JBond.and.
     $               IAMap.eq.JAngle).and.
     $         .not.(IMap.eq.JAngle.and.IBMap.eq.JBond.and.
     $               IAMap.eq.JAtom)) goto 30
              If(AUsed(IAtom)) NumA(IAtom) = JAtom
              AUsed(IAtom) = .false.
   30       If(JAtom.eq.3.and.IZ(IAtom,4).ne.0) goto 40
            If(IZ(IAtom,4).ne.0) goto 40
            If(JAtom.eq.3) goto 70
            If(.not.(IMap.eq.JAtom.and.IBMap.eq.JBond.and.
     $               IAMap.eq.JAngle.and.IDMap.eq.JDih).and.
     $         .not.(IMap.eq.JDih.and.IBMap.eq.JAngle.and.
     $               IAMap.eq.JBond.and.IDMap.eq.JAtom)) goto 70
              If(BUsed(IAtom)) NumB(IAtom) = JAtom
              If(BUsed(IAtom).and.(Be(IAtom)*Be(JAtom)).lt.Zero)
     $            NumB(IAtom) = -NumB(IAtom)
              BUsed(IAtom) = .false.
              goto 70
   40       If(.not.(IMap.eq.JAtom.and.IAMap.eq.JAngle.and.
     $               IDMap.eq.JDih.and.IZ(JAtom,4).ne.0).and.
     $         .not.(IMap.eq.JDih.and.IAMap.eq.JAngle.and.
     $               IDMap.eq.JAtom.and.IZ(JAtom,4).ne.0)) goto 50
              If(BUsed(IAtom)) NumB(IAtom) = JAtom
              BUsed(IAtom) = .false.
              Goto 70
   50       If(.not.(IMap.eq.JAtom.and.IBMap.eq.JBond.and.
     $               IDMap.eq.JAngle).and.
     $         .not.(IMap.eq.JAngle.and.IBMap.eq.JBond.and.
     $               IDMap.eq.JAtom)) goto 60
              If(BUsed(IAtom)) NumB(IAtom) = -JAtom
              BUsed(IAtom) = .false.
              Goto 70
   60       If(.not.(IMap.eq.IAtom.and.IBMap.eq.IBond.and.
     $               IDMap.eq.IAngle.and.JAtom.eq.LimJAt)) goto 70
              If(BUsed(IAtom)) NumB(IAtom) = -IAtom
              BUsed(IAtom) = .false.
   70       Continue
  999 Return
      End
*Deck RedZM
      Subroutine RedZM(IStruc,NAlSav,NAtoms,IAnSav,CSave,MaxAtm,IAn,C,
     $  MicOpt)
      Implicit Real*8(A-H,O-Z)
C
C     Setup the z-matrix to look like cartesians when using redundant
C     internal coordinates.
C
      Dimension IAnSav(NAlSav,*), CSave(3,NAlSav,*), IAn(*), C(*),
     $  MicOpt(*)
#include "commonz.inc"
C
      Call IClear(MaxAtm,IAn)
      Call IMove(NAtoms,IAnSav(1,IStruc),IAn)
      Call AClear(3*MaxAtm,C)
      Call AMove(3*NAtoms,CSave(1,1,IStruc),C)
C
      NZ = NAtoms
      Call IClear(MaxNZ,IAnZ)
      Call IMove(NAtoms,IAnSav(1,IStruc),IAnZ)
      Call IClear(MaxNZ*4,IZ)
      Call AClear(MaxNZ,Bl)
      Call AClear(MaxNZ,Alpha)
      Call AClear(MaxNZ,Beta)
      Do 10 I = 1, NAtoms
        IZ(I,1) = MOp2IZ(MicOpt(I))
        Bl(I) = CSave(1,I,IStruc)
        Alpha(I) = CSave(2,I,IStruc)
   10   Beta(I) = CSave(3,I,IStruc)
      Call IClear(MaxNZ,LBl)
      Call IClear(MaxNZ,LAlpha)
      Call IClear(MaxNZ,LBeta)
      Call IClear(MaxZSb,IntVec)
      Call AClear(MaxZSb,FPVec)
      Call AClear(MaxZSb,Values)
      Return
      End
*Deck RGaus
      Subroutine RGaus(IFileI,IPrint,PrtInI,Name,IUnits,InpFF,PDBInI,
     $  ConInp,HavCon,MaxPer,NPer,PerCnt,MaxTT,NTitle,Title,IRtCrd,
     $  LenRte,IChOni,MulOni,IChFrg,MulFrg,MaxCon,NBond,IBond,IBType,
     $  RBType,MaxLay,NuOfLa,HavLay,NDist,Layer,IAnMod,IAtMod,ACMod,
     $  LiAtom,Distan,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Read a Gaussian input file and return the structure, connectivity,
C     PDB info, etc.  The file is assumed to be open on unit IFileI.
C     This routine munges /IO/ temporarily because the underlying
C     routines use In for the input unit.  Commons /Mol/ and /ZMat are
C     loaded.
C
#include "commonmol.inc"
#include "commonz.inc"
      Logical PrtInI, PDBInI, PDBInp, ConInp, PrtInp, HavCon, HavLay(*)
      Character*(*) Name, PerCnt(*), Title(*)
      Common /IO/ In, IOut, IPunch
      Common /PhyCon/ PhyCon(30)
      Dimension IRtCrd(*), IChOni(*), MulOni(*), IChFrg(*), MulFrg(*),
     $  NBond(*), IBond(*), IBType(*), RBType(*), Layer(*), IAnMod(*),
     $  IAtMod(*), ACMod(*), LiAtom(*), Distan(*), V(*), NChMu(2), JJ(1)
      Save JJ
      Data JJ/0/
 1150 Format(1X,A,' Z-matrix:')
C
      InSave = In
      In = IFileI
      PrtInp = PrtInI.or.IPrint.gt.0
      If(PrtInp) Write(IOut,1150) Name
      NDChMu = Max(MaxLay**2,MaxAtm)
      IIChIn = 1
      IMulIn = IIChIn + 2*NDChMu
      IT = IMulIn + 2*NDChMu
      Call TstCor(IT+4*NDChMu-1,MDV,'RGaus-RdTop')
      Call RdTop(IFileI,IOut,IPrint,MaxPer,NPer,Percnt,MaxTT,NTitle,
     $  Title,IRtCrd,LenRte,NDChMu,NChMu,V(IIChIn),V(IMulIn),V(IT))
      Call AllZAt('RGauss',.True.,MaxNZ,MaxZSb,IT,IIAtWZ,IAtmWZ,IAtSpZ,
     $  IAtZEZ,IAtQMZ,IAtGFZ,IAtFrZ,IAtRsZ,ISymbl,IAtPDZ,IAtZNZ,INamCt,
     $  IV,V,MDV)
      MDV1 = MDV - IV + 1
      NuOfLa = MaxLay
      Call IniMol(1,0)
      Call IClear(MaxNZ,Layer)
      Call IClear(MaxNZ,IAnMod)
      Call IClear(MaxNZ,IAtMod)
      Call IClear(MaxNZ,LiAtom)
      Call AClear(MaxNZ,ACMod)
      Call AClear(MaxLay*MaxNZ,Distan)
      Call ZGet(IOut,PrtInp,.False.,IUnits,PhyCon(1),JJ,0,.True.,
     $  .False.,NuOfLa,ICharg,Multip,.True.,MaxNZ,IAnZ,IZ,Bl,Alpha,
     $  Beta,LBl,LAlpha,LBeta,NZ,NVarZM,IAtTyZ,AtmChZ,V(IIAtWZ),
     $  V(IAtmWZ),V(IAtSpZ),V(IAtZEZ),V(IAtQMZ),V(IAtGFZ),V(IAtFrZ),
     $  V(IAtRsZ),V(IAtPDZ),V(IAtZNZ),NumTpS,MaxAtm,IAtTpS,ANames,
     $  Values,IntVec,FPVec,V(ISymbl),V(INamCt),Layer,IAnMod,IAtMod,
     $  ACMod,LiAtom,Distan,InpFF,MaxLay,MaxRes,NResid,IResNa,IResNu)
      PDBInp = PDBInI.or.NResid.gt.0
      ConInp = ConInp.or.PDBInp
      NDist = MaxNZ
      Call ZTocFu(IOut,IPrint,.False.,.True.,V(IIAtWZ),V(IAtmWZ),
     $  V(IAtSpZ),V(IAtZEZ),V(IAtQMZ),V(IAtGFZ),V(IAtFrZ),V(IAtRsZ),
     $  V(IAtPDZ),V(IAtZNZ),NAtoms,IAn,IAtTyp,AtmChg,AtChMM,C,IAtWgt,
     $  AtmWgt,IAtSpn,AtZEff,AtQMom,AtGFac,IAtFrg,IAtRes,IAtPDB,AtZNuc,
     $  V(IV),MDV1)
      NFrag = IArMax(.False.,NZ,IAtFrg)
      Call CkLay(NAtoms,Layer,NuOfLa)
      IIScr1 = IV
      IIScr2 = IIScr1 + MaxNZ
      IV = IIScr2 + MaxNZ
      Call TstCor(IV,MDV,'RGauss-DoCon')
      MDV1 = MDV - IV + 1
      If(ConInp) then
        Call MkMpZA(NZ,IAnZ,NAtoms,V(IIScr1))
        Call DoCon(IOut,IPrint,2,0,0,MaxCon,NAtoms,IAn,C,NZ,V(IIScr1),
     $    ToAng,MxBond,NBond,IBond,IBType,RBType,V(IV),MDV1)
        HavCon = .True.
        endIf
      If(PDBInp) then
        Call SecGet(IFileI,IOut,IPrint,2,NAtoms,IAn,IAtFrg,IAtRes,
     $    IAtPDB,NResid,IResNa,IResNu,NSecSt,SecTyp,SecID,SecRes,SecDat,
     $    SecAtm,V(IIScr1))
        If(IPrint.ge.2) Call DmpPDB(IOut,.True.,.True.,NAtoms,IAtPDB,
     $    IAtRes,NResid,IResNa,IResNu)
        endIf
      NuOfLa = 0
      If(NZ.eq.0) then
        NZ = NAtoms
        Call IMove(NAtoms,IAn,IAnZ)
        endIf
      Call FixOni(IOut,.False.,MaxLay,NuOfLa,NZ,HavLay,IAnMod,IAnZ,
     $  LiAtom,Layer,IAtMod,ACMod,NDist,Distan,MaxCon,NBond,IBond,
     $  RBType,V(IIScr1))
      Call DefChM(IOut,IPrint,.False.,ICharg,Multip,MaxLay,NuOfLa,
     $  NFrag,NDChMu,NChMu,V(IIChIn),V(IMulIn),IChOni,MulOni,IChFrg,
     $  MulFrg)
      Call RdChrg(IOut,.False.,.False.,NuOfLa,NFrag,0,MaxLay,IChOni,
     $  MulOni,IChFrg,MulFrg)
      Call IClear(MaxNZ,V(IIScr1))
      Call GenITR(In,IOut,IPrint-1,.True.,.False.,.False.,.False.,0,
     $  MaxNZ,NZ,IZ,IAnZ,NuOfLa,Layer,LiAtom,IAnMod,NAtoms,IAn,IAtTyp,
     $  IAtRes,NResid,IResNa,IResNu,V(IIScr1),MicOpt,V(IIScr2),V(IV),
     $  MDV1)
      In = InSave
      Return
      End
*Deck RotVV
      Subroutine RotVV(Axis,Angle,V,VP)
      Implicit Real*8(A-H,O-Z)
C
C     Rotate vector V by Angle about Axis.  V is decomposed into a
C     component parallel to the Axis, which is invariant, and a
C     component in the plane perpendicular to the axis, which is
C     rotated.  A positive angle corresponds to rotating the vector
C     counterclockwise (or the axes clockwise).  V and VP may be the
C     same.
C
      Dimension Axis(3), V(3), VP(3), VXA(3)
      Save One
      Data One/1.0d0/
C
      R = Sqrt(SProd(3,Axis,Axis))
      DInv = SProd(3,Axis,V) * (One-Cos(Angle)) / (R*R)
      Call VProd(VXA,Axis,V)
      Call AScale(3,Cos(Angle),V,VP)
      Call ACasB(3,VP,Axis,VP,DInv)
      Call ACasB(3,VP,VXA,VP,(Sin(Angle)/R))
      Return
      End
*Deck SaveZM
      Subroutine SaveZM(In,IOut,ReadFz,NoFrez,FixTet,SaveZ,IRwZSv,Keep)
      Implicit Real*8(A-H,O-Z)
C
C     Finish up building /ZMat/ and /ZSubst/ and initializing /Grdnt/
C     and store them on disk.
C
      Parameter (IGrdnt=511)
      Logical ReadFz, NoFrez, FixTet, SaveZ, Keep(*), LJunk, VarRnd,
     $  TetRnd
#include "commonmol.inc"
#include "commonz.inc"
      Common /Grdnt/ Energy, NVarG, IGetFC
      Save Zero
      Data Zero/0.0D0/
 1110 Format(' OptVar=',50L1)
 1120 Format(50L1)
 1130 Format(' Any frozen variables have been unfrozen.')
 1140 Format(1X,I4,' tetrahedral angles replaced.')
C
      If(ReadFz.and.NVarZM.gt.0) then
        Read(In,1120) (Keep(I),I=1,NVarZM)
        Write(IOut,1110) (Keep(I),I=1,NVarZM)
        Call FFRead(IEOF)
        Do 910 I = 1, NVarZM
          If(.not.Keep(I)) IntVec(I) = -1
  910     Continue
        endIf
      If(NoFrez.and.NVarZM.gt.0) then
        Write(IOut,1130)
        Do 920 I = 1, NVarZM
          If(IntVec(I).eq.-1) IntVec(I) = 0
  920     Continue
        endIf
      If(NVarZM.gt.0) then
        Energy = Zero
        NVarG = NVarZM
        IGetFC = 0
        LGrd1 = 1 + InToWP(2)
        LGrd2 = NVarZM + (NVarZM*(NVarZM+1))/2
        Call ConDDF(IGrdnt,LGrd1+LGrd2)
        Call FileIO(1,-IGrdnt,LGrd1,Energy,0)
        If(FixTet) LJunk = VarRnd(MaxNZ,IZ,NZ,LBl,LAlpha,LBeta,NVarZM,
     $    Values)
        Call UpdVar(IOut,NZ,NVarZM,LBl,LAlpha,LBeta,Bl,Alpha,Beta,
     $    Values)
        endIf
      If(FixTet) then
        NumTet = 0
        Do 930 I = 3, NZ
          If(IZ(I,1).gt.0) then
            If(TetRnd(Alpha(I))) NumTet = NumTet + 1
            If(TetRnd(Beta(I))) NumTet = NumTet + 1
            endIf
  930     Continue
        If(NumTet.gt.0) Write(IOut,1140) NumTet
        endIf
      If(SaveZ) Call RdWrZ(1,IRwZSv,ANames)
      Return
      End
*Deck SaveZR
      Subroutine SaveZR(In,IOut,ReadFz,NoFrez,NStrAl,NAlRed,IStruc,
     $  NStruc,NAtoms,IZRed,NTRedS,Value,IntVec,FPVec,Keep,NAlSav,CRP,
     $  DoPath,MicOpt)
      Implicit Real*8(A-H,O-Z)
C
C     Set up /ZRedun/ on disk.
C     CRP: cartesian coordinates for multiple-structure input
C     NKeep : number of structures for which the coordinates are written on
C             Disk (Nkeep=0 for regular input, NKeep=2,3 for path relaxation)
C
      Parameter (IOZRed=665)
      Logical ReadFz,NoFrez,Keep,DoPath
      Dimension IZRed(4,NAlRed,*), NTRedS(4,*), Value(NAlRed,*),
     $  IntVec(*), FPVec(*), Keep(*), CRP(3,NAlSav,*), MicOpt(*)
 1110 Format(' OptVar=',50L1)
 1120 Format(50L1)
 1130 Format(' Any frozen variables have been unfrozen.')
C
      If(IStruc.gt.NStrAl) Call GauErr('Allocation error 1 in SaveZR.')
      NVar = NTRedS(1,IStruc)
      NAt3 = 3*NAtoms
      If(DoPath) then
        NKeep = Max(2,NStruc)
        NNK = NKeep
        If(NKeep.gt.3) Call Gauerr('Allocation error 2 in SaveZR.')
      else
        NKeep = 0
        NNK = 2
        endIf
      LRwZR = InToWP(5*NVar) + InToWP(4) + 2*NVar + NNK*NVar
     $  + NKeep*NAt3
      ICnst1 = Min(Mod(IStruc,3)+1,NStrAl)
      ICnst2 = Min(Mod(IStruc+1,3)+1,NStrAl)
      ICnst1 = 1
      ICnst2 = 2
      If(NKeep.eq.3) ICnst3 = 3
      Call ConDDF(IOZRed,LRwZR)
      If(ReadFz.and.NVar.gt.0) then
        Read(In,1120) (Keep(I),I=1,NVar)
        Write(IOut,1110) (Keep(I),I=1,NVar)
        Call FFRead(IEOF)
        If(DoPath) then
          Write(IOut,1130)
        else
          Do 10 I = 1, NVar
            If(.not.Keep(I)) then
              IntVec(I) = -1
              Value(I,ICnst1) = Value(I,IStruc)
              endif
   10       Continue
          endIf
        endIf
      Call FileIO(1,-IOZRed,InToWP(4),NTRedS(1,IStruc),0)
      Call FileIO(1,IOZRed,InToWP(4*NVar),IZRed(1,1,IStruc),0)
      Call FileIO(1,IOZRed,NVar,Value(1,IStruc),0)
      Call FileIO(1,IOZRed,InToWP(NVar),IntVec,0)
      Call FileIO(1,IOZRed,NVar,FPVec,0)
      Call FileIO(1,IOZRed,NVar,Value(1,ICnst1),0)
      Call FileIO(1,IOZRed,NVar,Value(1,ICnst2),0)
      If(NKeep.eq.3) Call FileIO(1,IOZRed,NVar,Value(1,ICnst3),0)
      Do 40 I = 1, NKeep
        Call FileIO(1,IOZRed,NAt3,CRP(1,1,I),0)
   40   Continue
      Return
      End
*Deck SwpAtm
      Subroutine SwpAtm(MaxNZ,MaxAP3,NAtoms,NOp,I,J,IAn,C,IZ,NPerm)
      Implicit Real*8(A-H,O-Z)
C
C     Interchange atoms I and J.  Note that only IZ(,1) is updated.
C
      Dimension IAn(1), C(3,1), V(3), IZ(MaxNZ,1), NPerm(MaxAP3,1)
C
      Call AMove(3,C(1,I),V)
      Call AMove(3,C(1,J),C(1,I))
      Call AMove(3,V,C(1,J))
      IT = IAn(I)
      IAn(I) = IAn(J)
      IAn(J) = IT
      IT = IZ(I,1)
      IZ(I,1) = IZ(J,1)
      IZ(J,1) = IT
      Do 10 K = 2, NAtoms
          If(IZ(K,1).eq.J) IZ(K,1) = -I
          If(IZ(K,1).eq.I) IZ(K,1) = J
   10     IZ(K,1) = IAbs(IZ(K,1))
      Do 20 K = 1, NOp
          IT = NPerm(I,K)
          NPerm(I,K) = NPerm(J,K)
   20     NPerm(J,K) = IT
      Do 30 K = 1, NAtoms
          Do 30 L = 1, NOp
              If(NPerm(K,L).eq.J) NPerm(K,L) = -I
              If(NPerm(K,L).eq.I) NPerm(K,L) = J
   30         NPerm(K,L) = IAbs(NPerm(K,L))
      Return
      End
*Deck SZPrnt
      Subroutine SZPrnt(Op,Val,Len,Type)
      Implicit Integer(A-Z)
C
C     PURPOSE:  TO BUILD A STRING CONTAINING THE SYMBOLIC Z-MATRIX AND
C               TO PRINT IT.  ITEMS WILL BE ADDED TO THE STRING ONE
C               PER CALL AND FORMATTED FOR OUPUT BY THIS ROUTINE.
C     INPUT:    OP -- OPCODE
C                     0 -- INITIALIZE BUFFER AND COUNTERS
C                     1 -- APPEND VAL TO BUFFER
C                     2 -- PRINT BUFFER AND RE-INITIALIZE
C               VAL -- QUANTITY TO BE PLACED IN THE BUFFER
C               TYPE -- DESCRIBES VAL
C                     0 -- INTEGER VALUE
C                     1 -- INTEGER VALUE
C                     2 -- FLOATING POINT VALUE (5 DECIMAL FIGURES).
C                     3 -- "LEN"  CHARACTER SYMBOL
C                    -3 -- "LEN" CHARACTER SYMBOL PREFIXED WITH '-'
C                   >=4 -- "LEN" CHARACTER SYMBOL in field of width "Type"
C
      Common/IO/In,IOut,IPunch
      Real*8 Val(*)
      Integer String(80)
      Save NB, String, Tab, Zero
      Data Zero/0.0d0/
 1000 Format(' Unrecognized value of Type in SZPrnt:',I7)
 1010 Format(' SZPrnt outputs null string.')
 1020 Format(' Bad op-code in SZPrnt:',I5)
C
      If(Op.eq.0) then
C       Initialize buffer
        NB = 0
        Call CBlank(80,String)
        Tab = 0
      else if(Op.eq.1) then
        NB = Tab
        If(Type.eq.0.or.Type.eq.1) then
C         Append integer Val to string.
          Call DecChr(Val(1),String,NB)
          Tab = Max(Tab+5,NB+1)
        else if(Type.eq.2) then
C         Append floating Val to string
          If(Val(1).ge.Zero) Call PutCSt(' ',String,NB)
          Call PutFP(Val,5,String,NB)
          NB = NB - 1
          Call PutCSt(' ',String,NB)
          Tab = Max(Tab+10,NB+1)
        else if(Abs(Type).eq.3) then
C         Append character string Val to string
          If(Type.eq.-3) Call PutCSt('-',String,NB)
          Call PutB(Val,Len,String,NB)
          Tab = Max(Tab+9,NB+1)
        else if(Type.ge.4) then
          Call PutB(Val,Len,String,NB)
          Tab = Max(Tab+Type,NB) + 1
        else
          Write(IOut,1000) Type
          Call Lnk1E(0)
          endIf
      else if(Op.eq.2) then
C       Flush buffer and reinitialize
        If(NB.eq.0) then
          Write(IOut,1010)
        else
          Call SetLnL(NB+2)
          Call StrOut(IOut,String,NB,1)
          Call SetLnL(0)
          endIf
        Tab = 0
        NB = 0
        Call CBlank(80,String)
      else
        Write(IOut,1020) Type
        Call Lnk1E(0)
        endIf
      Return
      End
*Deck Union2
      Subroutine Union2(IOut,NAtoms,NAlRed,IZRed1,NTRed1,NTBon1,NTAng1,
     $  NTDih1,IntVc1,IZRed2,NTRed2,NTBon2,NTAng2,NTDih2,IntVc2,IZRed,
     $  NTRed,NTBond,NTAng,NTDih,IntVec,IVal,Ind,IZTemp)
      Implicit Real*8(A-H,O-Z)
C
C     Input:  two redundant coordinate systems
C     Output: the union of the two
C     Chun Yang Peng 1993-1994
C
      Logical Error
      Dimension IZRed1(4,*), IntVc1(*), IZRed2(4,*), IntVc2(*),
     $  IZRed(4,*), IntVec(*), IVal(*), Ind(*), IZTemp(4,*)
 1000 Format(' Inconsistency between reactant coordinate',I7,' IZ=',4I7,
     $  ' Flag=',I8,/,'                    and product coordinate',I7,
     $  ' IZ=',4I7,' Flag=',I8,'.')
C
C     Bonds first.
C
      If((NTRed1.ne.(NTBon1+NTAng1+NTDih1)).or.
     $  (NTRed2.ne.(NTBon2+NTAng2+NTDih2)))
     $  Call GauErr('Consistency failure in Union2.')
      Error = .False.
      NTRed = 0
      Do 10 I = 1, NTBon1
        NTRed = NTRed + 1
        IntVec(NTRed) = IntVc1(I)
        IZRed(1,NTRed) = IZRed1(1,I)
        IZRed(2,NTRed) = IZRed1(2,I)
        IZRed(3,NTRed) = IZRed1(3,I)
   10   IZRed(4,NTRed) = IZRed1(4,I)
      Do 30 I = 1, NTBon2
        ia = IZRed2(1,I)
        ja = IZRed2(2,I)
        Do 20 J = 1, NTBon1
          ib = IZRed1(1,J)
          jb = IZRed1(2,J)
          If((ia.eq.ib.and.ja.eq.jb).or.(ia.eq.jb.and.ib.eq.ja)) then
            If(IntVc1(J).ne.IntVc2(I)) then
              Write(IOut,1000) I, (IZRed2(KK,I),KK=1,4), IntVc2(I), J,
     $          (IZRed1(KK,J),KK=1,4), IntVc1(J)
              Error = .True.
              endIf
            Goto 30
            endIf
   20     Continue
        NTRed = NTRed + 1
        If(NTRed.gt.NAlRed) Call GauErr('Allocation failure in Union2.')
        IntVec(NTRed) = IntVc2(I)
        IZRed(1,NTRed) = IZRed2(1,I)
        IZRed(2,NTRed) = IZRed2(2,I)
        IZRed(3,NTRed) = IZRed2(3,I)
        IZRed(4,NTRed) = IZRed2(4,I)
   30   Continue
      NTBond = NTRed
C
C     Angles
C
      Do 40 i = (NTBon1+1), (NTBon1+NTAng1)
        NTRed = NTRed + 1
        If(NTRed.gt.NAlRed) Call GauErr('Allocation failure in Union2.')
        IntVec(NTRed) = IntVc1(I)
        IZRed(1,NTRed) = IZRed1(1,I)
        IZRed(2,NTRed) = IZRed1(2,I)
        IZRed(3,NTRed) = IZRed1(3,I)
   40   IZRed(4,NTRed) = IZRed1(4,I)
      NTAng = NTRed - NTBond
      Do 60 I = (NTBon2+1), (NTBon2+NTAng2)
        IA = IZRed2(1,I)
        JA = IZRed2(2,I)
        KA = IZRed2(3,I)
        LA = IZRed2(4,I)
        Do 50 J = (NTBond+1), (NTBond+NTAng)
          IB = IZRed(1,J)
          JB = IZRed(2,J)
          KB = IZRed(3,J)
          LB = IZRed(4,J)
          If(JA.eq.JB.and.
     $      ((IA.eq.IB.and.KA.eq.KB).or.(IA.eq.KB.and.KA.eq.IB))) then
            J1 = J - NTBond + NTBon1
            If(LA.ne.0.and.LB.eq.0) then
              IZRed(1,J) = IA
              IZRed(2,J) = JA
              IZRed(3,J) = KA
              IZRed(4,J) = LA
              If(IntVc1(J1).ne.IntVc2(I)) then
                Write(IOut,1000) I, (IZRed2(KK,I),KK=1,4), IntVc2(I),
     $            J1, (IZRed1(KK,J1),KK=1,4), IntVc1(J1)
                Error = .True.
                endIf
              Goto 60
            else if(LA.eq.LB.or.(LA.eq.0.and.LB.ne.0)) then
              If(IntVc1(J1).ne.IntVc2(I)) then
                Write(IOut,1000) I, (IZRed2(KK,I),KK=1,4), IntVc2(I),
     $            J1, (IZRed1(KK,J1),KK=1,4), IntVc1(J1)
                Error = .True.
                endIf
              Goto 60
              endIf
            endIf
   50     Continue
        NTRed = NTRed + 1
        If(NTRed.gt.NAlRed) Call GauErr('Allocation failure in Union2.')
        IntVec(NTRed) = IntVc2(I)
        IZRed(1,NTRed) = IA
        IZRed(2,NTRed) = JA
        IZRed(3,NTRed) = KA
        IZRed(4,NTRed) = LA
   60   Continue
      NTAng = NTRed - NTBond
C
C     Dihedrals
C
      Do 80 I = (NTBon1+NTAng1+1), (NTBon1+NTAng1+NTDih1)
C       check that this dihedral does not involve a linear bend.
        IA = IZRed1(1,I)
        JA = IZRed1(2,I)
        KA = IZRed1(3,I)
        LA = IZRed1(4,I)
        Do 70 J = (NTBond+1), (NTBond+NTAng)
          IB = IZRed(1,J)
          JB = IZRed(2,J)
          KB = IZRed(3,J)
          LB = IZRed(4,J)
          If(LB.lt.0.and.((IA.eq.IB.and.JA.eq.JB.and.KA.eq.KB).or.
     $      (JA.eq.IB.and.KA.eq.JB.and.LA.eq.KB))) goto 80
   70     Continue
        NTRed = NTRed + 1
        If(NTRed.gt.NAlRed) Call GauErr('Allocation failure in Union2.')
        IntVec(NTRed) = IntVc1(I)
        IZRed(1,NTRed) = IA
        IZRed(2,NTRed) = JA
        IZRed(3,NTRed) = KA
        IZRed(4,NTRed) = LA
   80   Continue
      NTDih = NTRed - NTAng - NTBond
      Do 110 I = (NTBon2+NTAng2+1), (NTBon2+NTAng2+NTDih2)
        IA = IZRed2(1,I)
        JA = IZRed2(2,I)
        KA = IZRed2(3,I)
        LA = IZRed2(4,I)
        Do 90 J = (NTBond+1), (NTBond+NTAng)
          IB = IZRed(1,J)
          JB = IZRed(2,J)
          KB = IZRed(3,J)
          LB = IZRed(4,J)
          If(LB.LT.0.and.((IA.eq.IB.and.JA.eq.JB.and.KA.eq.KB).or.
     $      (JA.eq.IB.and.KA.eq.JB.and.LA.eq.KB))) goto 110
   90     Continue
        Do 100 J = (NTBond+NTAng+1), (NTBond+NTAng+NTDih)
          IB = IZRed(1,J)
          JB = IZRed(2,J)
          KB = IZRed(3,J)
          LB = IZRed(4,J)
          If((IA.eq.IB.and.JA.eq.JB.and.KA.eq.KB.and.LA.eq.LB).or.
     $       (IA.eq.LB.and.JA.eq.KB.and.KA.eq.JB.and.LA.eq.IB)) then
            J1 = J - NTBond - NTAng + NTBon1 + NTAng1
            If(IntVc1(J1).ne.IntVc2(I)) then
              Write(IOut,1000) I, (IZRed2(KK,I),KK=1,4), IntVc2(I),
     $          J1, (IZRed1(KK,J1),KK=1,4), IntVc1(J1)
              Error = .True.
              endIf
            Goto 110
            endIf
  100     Continue
        NTRed = NTRed + 1
        If(NTRed.gt.NAlRed) Call GauErr('Allocation failure in Union2.')
        IntVec(NTRed) = IntVc2(I)
        IZRed(1,NTRed) = IA
        IZRed(2,NTRed) = JA
        IZRed(3,NTRed) = KA
        IZRed(4,NTRed) = LA
  110   Continue
       NTDih = NTRed - NTAng - NTBond
       Call SrtRed(NAtoms,NTRed,NTBond,NTAng,NTDih,IZRed,IntVec,IVal,
     $   Ind,IZTemp)
       If(Error)
     $   Call GauErr('Inconsistent reactant and product variables.')
       Return
       End
*Deck UpMIAt
      Subroutine UpMIAt(IFFIn,NAtoms,IAType,IAtTyp)
      Implicit Integer(A-Z)
C
C     Update atom types in IAtTyp with those generated during MM
C     parameter setup.
C
      Dimension IAType(NAtoms,2), IAtTyp(NAtoms)
C
      Do 10 I = 1, NAtoms
        If(IAtTyp(I).eq.0) IAtTyp(I) = IFFin*10000000 +
     $    IAType(I,2)*10000 + IAType(I,1) + 1000
   10   Continue
      Return
      End
*Deck UpZMCh
      Subroutine UpZMCh(IOut,IPrint,IRwMM,IParIn,NAtoms,IAtTyp,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Update the MM charges in the Z-matrix from those in the MM parameter
C     file, if both exist.
C
      Dimension IAtTyp(*), V(MDV), JJ(1)
      Character ShrtNm*20, Name*40
#include "commonz.inc"
      Save Zero, JJ
      Data Zero/0.0d0/, JJ/0/
 1000 Format(' No MM parameters found.')
C
      If(ITqry(IRwMM).lt.1) then
        If(IPrint.gt.0) Write(IOut,1000)
      else
        Call FncInf(-1,MaxFnc,Name,ShrtNm,NClass,NCentr,NInts,NDReal,
     $    IStrct,JJ,MaxPt,IntTot,JJ,MaxCnv)
        Call MMAloc(MaxPt,MaxCnv,NCentr,MaxFnc,NAtoms,NClass,IntTot,
     $    NDInt,INCPar,IIPrt,ICnvrt,IICntr,IIList,IAtL,IChL,IQ,IEClas,
     $    IEFunc,IAType,Len0,MDV)
        Call FileIO(2,-IRwMM,Len0,V,0)
        IFFIn = Mod(IParIn,100)
        Call UpMIAt(IFFIn,NAtoms,V(IAType),IAtTyp)
        If(NZ.gt.0) then
          IA = 0
          Do 20 I = 1, NZ
            If(IAnZ(I).ge.0) then
              IA = IA + 1
              IAtTyZ(I) = IAtTyp(IA)
              AtmChZ(I) = V(IQ+IA-1)
            else
              IAtTyZ(I) = 0
              AtmChZ(I) = Zero
              endIf
   20       Continue
          If(IA.ne.NAtoms) Call GauErr('Consistency failure in UpZMCh.')
          If(IPrint.gt.0)
     $      Call OutMtS(IOut,'New AtmChZ:',0,0,AtmChZ,NZ,1,NZ,1)
          Call RdWrZ(1,0,V)
          endIf
        endIf
      Return
      End
*Deck VarRnd
      Logical Function VarRnd(MaxNZ,IZ,NZ,LBl,LAlpha,LBeta,NVar,Values)
      Implicit Real*8(A-H,O-Z)
C
C     Round Z-matrix variables to the exact tetrahedral angles.
C     The value is .true. if
C
      Logical TetRnd, VDist
      Dimension IZ(MaxNZ,4), LBl(*), LAlpha(*), LBeta(*), Values(*)
C
      VarRnd = .False.
      Do 100 IVar = 1, NVar
        If(.not.VDist(IVar,MaxNZ,NZ,LBl,LAlpha,LBeta,IZ)) then
          If(TetRnd(Values(IVar))) VarRnd = .True.
          endIf
  100   Continue
      Return
      End
*Deck ZCart
      Subroutine ZCart(DoCopy,NAtoms,IAn,IAtTyp,AtChMM,C)
      Implicit Real*8(A-H,O-Z)
C
C     Generate a z-matrix with all Cartesian variables.
C
      Logical DoCopy
      Dimension IAn(NAtoms), IAtTyp(NAtoms), AtChMM(NAtoms), C(3,NAtoms)
#include "commonz.inc"
C
      If(DoCopy) then
        Call IClear(MaxNZ,IAnZ)
        Call IMove(NAtoms,IAn,IAnZ)
        Call IClear(MaxNZ,IAtTyZ)
        Call IMove(NAtoms,IAtTyp,IAtTyZ)
        Call AClear(MaxNZ,AtmChZ)
        Call AMove(NAtoms,AtChMM,AtmChZ)
        NZ = NAtoms
        endIf
      NVarZM = 3*NZ
      Call IClear(4*MaxNZ,IZ)
      Call AClear(MaxZSb,Values)
      Call IClear(MaxZSb,IntVec)
      Call AClear(MaxZSb,FPVec)
      Do 10 I = 1, NZ
        If(DoCopy) then
          Bl(I) = C(1,I)
          Alpha(I) = C(2,I)
          Beta(I) = C(3,I)
          endIf
        LBl(I) = 3*I - 2
        LAlpha(I) = 3*I - 1
        LBeta(I) = 3*I
        Values(LBl(I)) = Bl(I)
        Values(LAlpha(I)) = Alpha(I)
   10   Values(LBeta(I)) = Beta(I)
      Call CartVN(NVarZM,ANames)
      Return
      End
*Deck ZCentr
      Subroutine ZCentr(IOut,Prnt,EndOK,FPOK,IZ,NAMCNT,NZ,FP)
      Implicit Integer(A-Z)
C
C     A routine to read (via "FF" routines) a center specification
C     from a Z-matrix card.  This specification may be either an
C     integer (the sequential number of a previous Z-matrix card),
C     or the name of a previously defined center.  "IZ" is returned
C     as the sequential number of the center being referenced.
C     "NamCnt" is a delimited Hollerith string containing the names
C     of the centers.  "NZ" is the sequential number of the current
C     Z-matrix card.  If a floating-point number is found, IZ is set
C     to -2 and the number is returned in FP.
C     M. J. Frisch - October 1984, various dates in 1992, 2001.
C
      Logical Prnt, EndOK, FPOK
      Real*8 FP, RFFGet, Zero
      Dimension NamCnt(*), STR(10)
      Save IFFStr, IFFInt, IFFFP, IFFEnd, Zero
      Data IFFStr/3hSTR/, IFFInt/3hINT/, IFFFP/3hFP /, IFFEnd/3hEND/,
     $  Zero/0.0d0/
 1000 Format(' Center specification must be either integer or string.')
C
C     Get the next symbol or whatever
      Call FFGet(Str,Len,I,FP,IF1)
C
C     If it's a string, check if it's a valid center name.  If it's
C     a positive integer, use it.  If it's a negative integer, get
C     the following FP number (cartesian coordinate).  If it's a
C     real, return it.
C
      If(IF1.eq.IFFStr) then
        If(Len.gt.4) then
          Call IniStr
          Call BldStC('Center name "')
          Call BldStr(Str,Len)
          Call BldStC('" is too long.')
          Call DmpStr(IOut)
          Call Lnk1E(0)
          endIf
        IZ = LSubst(NamCnt,NZ-1,Str,Len)
        If(Prnt) Call SZPrnt(1,Str,Len,4)
      else if(IF1.eq.IFFInt) then
        IZ = I
        If(Prnt) Call SZPrnt(1,IZ,0,1)
        If(IZ.lt.0) then
          FP = RFFGet(.False.,.False.,Zero,.False.)
          If(Prnt) Call SZPrnt(1,FP,0,2)
          endIf
      else if(IF1.eq.IFFFP.and.FPOK) then
        IZ = -999999
        If(Prnt) Call SZPrnt(1,FP,0,2)
      else if(IF1.eq.IFFEnd.and.EndOK) then
        IZ = 1
      else
        Write(IOut,1000)
        Call FFErr(0,IF1)
        endIf
      Return
      End
*Deck ZGet
      Subroutine ZGet(IOut,PrtInp,ScanDf,InpUnt,ToAng,BB,NBB,RdInpt,
     $  RdChar,NuOfLa,ICharg,Multip,ConsFz,MaxNZ,IAnZ,IZ,Bl,Alpha,Beta,
     $  LBl,LAlpha,LBeta,NZ,NVar,IAtTyZ,AtmChZ,IAtWtZ,AtmWtZ,IAtSpZ,
     $  AtZEfZ,AtQMoZ,AtGFaZ,IAtFrZ,IAtRsZ,IAtPDZ,AtZNuZ,NumNam,MaxNam,
     $  NamLst,ANames,Values,IntVec,FPVec,Symbls,NamCnt,Layer,IAnMod,
     $  IAtMod,ACMod,LiAtom,Distan,InpFF,MaxLay,MaxRes,NRes,IResNa,
     $  IResNu)
      Implicit Real*8(A-H,O-Z)
C
C     If(RdChar) read the charge and multiplicity.
C     Reach up to 3 geometry input sections and fill /ZMat/ and
C     /ZSubst/.  RdInpt indicates whether to read from the input
C     stream or from the array BB, starting at cursor NBB.  If
C     BB is used, units are not converted.  If the input stream is
C     read, the z-matrix is printed and units are converted according
C     to InpUnt.  PrtInp indicates whether in input is to be printed.
C     ConsFz indicates whether constants should be frozen variables
C     (.true.) or wired into the Z-Matrix (.False.).
C
      Logical OK, IV0, RdInpt, RdChar, PrtInp, ConsFz, DoingV, GotFP,
     $  RdZVar, GotVal, EOF, BVar, AVar, Varabl, ScanDf
      Dimension IAnZ(MaxNZ),IZ(MaxNZ,4),Bl(MaxNZ),Alpha(MaxNZ),
     $  Beta(MaxNZ),LBl(MaxNZ),LAlpha(MaxNZ),LBeta(MaxNZ),ANames(*),
     $  Values(*),IntVec(*),FPVec(*),Symbls(*),NamCnt(*),TStr(6),
     $  Layer(*),IAnMod(*),IAtMod(*),ACMod(*),LiAtom(*),Distan(*),
     $  IAtTyZ(*),AtmChZ(*),IAtWtZ(*),AtmWtZ(*),IAtSpZ(*),AtZEfZ(*),
     $  AtQMoZ(*),AtGFaZ(*),IAtFrZ(*),IAtRsZ(*),IAtPDZ(*),NamLst(2,*),
     $  IResNa(*),IResNu(*),AtZNuZ(*)
      Integer BB(*)
      Save Zero,One,F45,IV0,GotFP,GotVal,LenT,TStr
      Data Zero/0.0d0/,One/1.D0/,F45/45.D0/,IV0/.False./,GotFP/.False./,
     $  GotVal/.False./,LenT/0/,TStr/6*0/
 1000 Format(' Integer parameters encountered on z-matrix card',I4,'.')
 1020 Format(' The maximum of',I5,' variables has been exceeded.')
 1070 Format(' Undefined symbol, bond length, card',i4,'.')
 1072 Format(' Undefined symbol, angle alpha, card',i4,'.')
 1074 Format(' Undefined symbol, angle beta, card',i4,'.')
 1080 Format('       Variables:')
 1090 Format('       Constants:')
 1100 Format(' Variable',I4,' is used for an angle and a distance.')
 1110 Format(' ')
C
C     Read the Z-matrix section.  Upon return from ZSymb:
C     NSymb = Total number of symbols (non-numbers) in Z-matrix.
C     LBl, LAlpha, LBeta:
C       1 for integer
C       2 for floating point
C       3 for symbol
C      -3 for -symbol
C
      MaxTST = 6*IntPWP(0)
      NSymb = 0
      MaxZSb = 3*MaxNZ
      Call IClear(MaxNZ,IAnZ)
      Call IClear(4*MaxNZ,IZ)
      Call IClear(MaxNZ,LBl)
      Call IClear(MaxNZ,LAlpha)
      Call IClear(MaxNZ,LBeta)
      Call IClear(MaxNZ,IAtTyZ)
      Call IClear(MaxZSb,IntVec)
      Call AClear(MaxZSb,FPVec)
      Call AClear(MaxNZ,BL)
      Call AClear(MaxNZ,ALpha)
      Call AClear(MaxNZ,Beta)
      Call AClear(MaxNZ,AtmChZ)
      Call CBlank(MaxZSb*IntPWP(0),ANames)
      ICur = 0
      Call PutDel(5,ANames,ICur)
      Call ZSymb(IOut,PrtInp,NSymb,Symbls,BB,NBB,RdInpt,RdChar,ICharg,
     $  Multip,MaxNZ,IAnZ,IZ,Bl,Alpha,Beta,LBl,LAlpha,LBeta,NZ,IAtTyZ,
     $  AtmChZ,IAtWtZ,AtmWtZ,IAtSpZ,AtZEfZ,AtQMoZ,AtGFaZ,IAtFrZ,
     $  IAtRsZ,IAtPDZ,AtZNuZ,NumNam,MaxNam,NamLst,NamCnt,Layer,IAnMod,
     $  IAtMod,ACMod,LiAtom,Distan,InpFF,MaxLay,NuOfLa,MaxRes,NRes,
     $  IResNa,IResNu)
C
C     Scan for integers as Z-matrix parameters and convert the code used
C     for constants.
C
      OK = .True.
      Do 40 I = 1, NZ
        If(LBl(I).eq.2) LBl(I) = 0
        If(LAlpha(I).eq.2) LAlpha(I) = 0
        If(LBeta(I).eq.2) LBeta(I) = 0
        If(LBl(I).eq.1.or.LAlpha(I).eq.1.or.LBeta(I).eq.1) then
          OK = .False.
          Write(IOut,1000) I
          endIf
   40   Continue
      If(.not.OK) Call Lnk1E(0)
C
C     Only 0 and +/- 3 remain in LBl, LAlpha, LBeta.  Convert these to
C     +/- 3000 so that they won't be confused as pointers to variable 3.
C
      Do 60 I = 1, NZ
        LBl(I) = LBl(I) * 1000
        LAlpha(I) = LAlpha(I) * 1000
   60   LBeta(I) = LBeta(I) * 1000
C
C     Top of loop over variables, active then frozen.  Pack the
C     variable names into ANames, make sure that they appear in the
C     Z-matrix, and possibly print the variable definitions.
C
      NVar = 0
      NCursr = 0
      Do 100 IPass = 1, 2
        If(NSymb.eq.0) goto 100
        DoingV = IPass.eq.1
        Varabl = DoingV.or.ConsFz
        If(PrtInp.and.IPass.eq.1) then
          Write(IOut,1080)
        else if(PrtInp) then
          Write(IOut,1090)
          endIf
C
C       While (.not.EOF)
   80   Val = Zero
        Int = 0
        FP = Zero
        EOF = .not.RdZVar(RdInpt,ScanDf,BB,NBB,TStr,MaxTST,LenT,GotVal,
     $    Val,IV0,Int,GotFP,FP)
        If(.not.EOF) then
          If(Varabl) then
            NVar = NVar + 1
            If(NVar.gt.MaxZSb) then
              Write(IOut,1020)
              Call Lnk1E(0)
              endIf
            Values(NVar) = Val
            IntVec(NVar) = Int
            FPVec(NVar) = FP
            Call PutB(TStr,LenT,ANames,NCursr)
            Call PutDel(2,ANames,NCursr)
            endIf
          Call ZMatch(IOut,OK,Symbls,TStr,LenT,Val,Varabl,NSymb,MaxNZ,
     $      IAnZ,IZ,Bl,Alpha,Beta,LBl,LAlpha,LBeta,NZ,NVar,ANames,
     $      Values,IntVec,FPVec)
          If(PrtInp) Call PrmOut(OK,IV0,TStr,LenT,Val,Int,FP)
          If(.not.OK) Call Lnk1E(0)
          If(ConsFz.and..not.DoingV.and.intvec(nvar).gt.-1000)
     $       IntVec(NVar) = -1
          Goto 80
          endIf
        If(NSymb.eq.0.or.(IPass.eq.1.and..not.ConsFz))
     $    Call PutDel(5,ANames,NCursr)
  100   Continue
C
C     Variables and constants have now been read.  Ensure that all
C     symbols have been defined.
C
      If(NSymb.lt.0) then
        Call GauErr('Duplicate variable definitions found.')
      else if(NSymb.gt.0) then
        Do 340 I = 1, NZ
          If(IAbs(LBl(I)).eq.3000) Write(IOut,1070) I
          If(IAbs(LAlpha(I)).eq.3000) Write(IOut,1072) I
          If(IAbs(LBeta(I)).eq.3000) Write(IOut,1074) I
  340     Continue
        Call Lnk1E(0)
        endIf
C
C     Convert to atomic units.
C
      If(RdInpt) then
        ToBohr = One
        ToRad = One
        If(InpUnt.eq.0.or.InpUnt.eq.2) ToBohr = One / ToAng
        If(InpUnt.le.1) ToRad = ATan(One) / F45
        Do 420 I = 1, NZ
          Bl(I) = Bl(I) * ToBohr
          If(IZ(I,1).gt.0) then
            Alpha(I) = Alpha(I) * ToRad
          else
            Alpha(I) = Alpha(I) * ToBohr
            endIf
          If(IZ(I,1).gt.0) then
            Beta(I) = Beta(I) * ToRad
          else
            Beta(I) = Beta(I) * ToBohr
            endIf
 420      Continue
        Do 460 I = 1, NVar
          BVar = .False.
          AVar = .False.
          Do 440 J = 1, NZ
            If(IZ(J,1).gt.0) then
              BVar = BVar.or.IAbs(LBl(J)).eq.I
              AVar = AVar.or.IAbs(LAlpha(J)).eq.I.or.IAbs(LBeta(J)).eq.I
            else
              BVar = BVar.or.IAbs(LBl(J)).eq.I
     $                   .or.IAbs(LAlpha(J)).eq.I
     $                   .or.IAbs(LBeta(J)).eq.I
              endif
  440       Continue
          If(BVar.and.AVar) then
            Write(IOut,1100) I
            Call Lnk1E(0)
          else if(BVar) then
            Values(I) = Values(I) * ToBohr
            If(IntVec(i).ne.1) FPVec(i) = FPVec(i) * ToBohr
          else
            Values(I) = Values(I) * ToRad
            If(IntVec(i).ne.1) FPVec(i) = FPVec(i) * ToRad
            endIf
  460     Continue
        endIf
      If(PrtInp) Write(IOut,1110)
      Return
      End
*Deck ZMatch
      Subroutine ZMatch(IOut,OK,Symbls,Name,Len,Value,IfVar,NSymb,MaxNZ,
     $  IAnZ,IZ,Bl,Alpha,Beta,LBl,LAlpha,LBeta,NZ,NVar,ANames,Values,
     $  IntVec,FPVec)
      Implicit Real*8(A-H,O-Z)
C
C     A routine to search for a match in the symbolic Z-matrix
C     for the symbol "symbol" (length "LenSMB").  The calling routine
C     has read the symbol and its value;  This routine does all of
C     the necessary processing of the data.  The arguments are:
C
C     OK ... Whether the symbol Name was found.
C     Symbls ... A delimited hollerith string containing the
C                names of the symbols encountered in the z-matrix.
C     Name ..... The name of the variable to be matched.
C     Len ...... The number of characters in "Name".
C     Value  ... The value of the parameter.
C     IfVar .... A logical variable indicating whether this
C                parameter is a variable (T) or constant (F).
C     NSymb .... Number of symbols remaining in the z-matrix.  This is
C                decremented by this routine as symbols are matched.
C
      Logical OK, IfVar
      Dimension IAnZ(MaxNZ), IZ(MaxNZ,4), Bl(MaxNZ), Alpha(MaxNZ),
     $  Beta(MaxNZ), LBl(MaxNZ), LAlpha(MaxNZ), LBeta(MaxNZ), ANames(*),
     $  Values(*), IntVec(*), FPVec(*), Name(*), Symbls(*)
C
C     Loop over Z-matrix cards searching for a match for "symbol".
C     on each card, look to see whether any of the parameters match.
C
      NCur = 0
      OK = .False.
      Do 180 I = 1, NZ
        Call CkSymb(IfVar,NVar,Name,Len,Value,Symbls,NCur,LBl(I),Bl(I),
     $    NSymb,OK)
        Call CkSymb(IfVar,NVar,Name,Len,Value,Symbls,NCur,LAlpha(I),
     $    Alpha(I),NSymb,OK)
        Call CkSymb(IfVar,NVar,Name,Len,Value,Symbls,NCur,LBeta(I),
     $    Beta(I),NSymb,OK)
  180   Continue
      If(.not.OK) then
        Call IniStr
        Call BldStC('Symbol "')
        Call BldStr(Name,Len)
        Call BldStC('" not found in Z-matrix.')
        Call DmpStr(IOut)
        Call Lnk1E(0)
        endIf
      Return
      End
*Deck ZMCord
      Subroutine ZMCord(NAtoms,MicOpt,IAn,IAtTyp,AtChMM,C)
      Implicit Real*8(A-H,O-Z)
C
C     Set up /ZMat/ for cartesian coordinates.
C
      Dimension MicOpt(*), IAn(*), IAtTyp(*), AtChMM(*), C(3,*)
#include "commonz.inc"
C
      Call IClear(MaxNZ,IAnZ)
      Call IClear(4*MaxNZ,IZ)
      Call AClear(MaxNZ,Bl)
      Call AClear(MaxNZ,Alpha)
      Call AClear(MaxNZ,Beta)
      Call IClear(MaxNZ,LBl)
      Call IClear(MaxNZ,LAlpha)
      Call IClear(MaxNZ,LBeta)
      NZ = NAtoms
      NVarZM = 0
      Call IClear(2,IPdZMt)
      Call IClear(MaxNZ,IAtTyZ)
      Call AClear(MaxNZ,AtmChZ)
      Call IMove(NAtoms,IAn,IAnZ)
      Call IMove(NAtoms,IAtTyp,IAtTyZ)
      Call AMove(NAtoms,AtChMM,AtmChZ)
      Do 10 I = 1, NAtoms
        IZ(I,1) = MOp2IZ(MicOpt(I))
        Bl(I) = C(1,I)
        Alpha(I) = C(2,I)
   10   Beta(I) = C(3,I)
      Return
      End
*DECK ZPARM
      SUBROUTINE ZPARM(Prnt,X,LX,SYMBLS,NSB,ISYMB,Pad,SymOK)
      IMPLICIT REAL*8(A-H,O-Z)
C
C          A ROUTINE TO READ A Z-MATRIX PARAMETER (VARIABLE OR
C     FLOATING POINT NUMBER).  USES THE "FF" ROUTINES.  IF A
C     CONSTANT IS FOUND, THIS IS STORED INTO "X", ANS "LX" IS SET
C     TO 2.  IN THIS CASE, A NAME OF "0" IS STORED INTO "SYMBLS".
C
C          IF A NAME (STRING) IS FOUND, THEN THIS NAME IS APPENDED
C     TO "SYMBLS", AND "ISYMB" IS INCREMENTED.  "X" IS SET TO ZERO,
C     AND "LX" IS SET TO 3 (OR -3 IF "-NAME" IS FOUND).
C
C          THIS ROUTINE ALSO CALLS THE Z-MATRIX OUTPUT ROUTINE FOR
C     EACH PARAMETER READ IN.
C     R. A. WHITESIDE - FEBRUARY 1980, M. J. Frisch -- October 1984.
C
      Logical Prnt, Pad, SymOK
      Integer STR(5), Symbls(*)
      Save Zero,IFFFP,IFFStr
      Data Zero/0.D0/, IFFFP/3hFP /, IFFStr/3hSTR/
C
      Call FFGet(Str,Len,I,FP,IF1)
      If(IF1.eq.IFFFP) then
        X = FP
        LX = 2
        If(Prnt) CALL SZPrnt(1,FP,0,2)
C       Put a phony parameter name into "Symbls".
        If(Pad) then
          Call PutCSt('0',Symbls,NSb)
          Call PutDel(2,Symbls,NSb)
          endIf
      else if(IF1.eq.IFFStr.and.SymOK) then
        Call NoSign(Str,Len,ISgn)
        Call PutB(Str,Len,Symbls,NSB)
        Call PutDel(2,Symbls,NSB)
        X = Zero
        LX = 3*ISgn
        ISymb = ISymb + 1
        If(Prnt) Call SZPrnt(1,Str,Len,LX)
      else
        Call FFErr(IFFFP,IF1)
        endIf
      Return
      End
*Deck ZSymb
      Subroutine ZSymb(IOut,PrtInp,ISymb,Symbls,BB,NBB,RdInpt,RdChar,
     $  ICharg,Multip,MaxNZ,IAnZ,IZ,Bl,Alpha,Beta,LBl,LAlpha,LBeta,NZ,
     $  IAtTyZ,AtmChZ,IAtWtZ,AtmWtZ,IAtSpZ,AtZEfZ,AtQMoZ,AtGFaZ,IAtFrZ,
     $  IAtRsZ,IAtPDZ,AtZNuZ,NumNam,MaxNam,NamLst,NamCnt,Layer,IAnMod,
     $  IAtMod,ACMod,LiAtom,Distan,InpFF,MaxLay,NuOfLa,MaxRes,NRes,
     $  IResNa,IResNu)
      Implicit Real*8(A-H,O-Z)
C
C     Read a symbolic z-matrix, either from the input stream or
C     from array BB.
C
      Logical RdInpt, RdChar, PrtInp, First, SymOK, StEqHC
      Integer BB(*), Symbls(*), Str(800)
      Character*800 CStr
      Dimension NamCnt(*), IAnZ(MaxNZ), IZ(MaxNZ,4), Bl(*), Alpha(*),
     $  Beta(*), LBl(*), LAlpha(*), LBeta(*), Layer(*), IAnMod(*),
     $  IAtMod(*), ACMod(*), LiAtom(*), Distan(*), IAtTyZ(*), AtmChZ(*),
     $  IAtWtZ(*), AtmWtZ(*), IAtSpZ(*), AtZEfZ(*), AtQMoZ(*),
     $  IAtPDZ(*), AtGFaZ(*), IAtFrZ(*), IAtRsZ(*), JJ(1), NamLst(2,*),
     $  IResNa(*), IResNu(*), AtZNuZ(*)
      Save IFFEnd, IFFStr, IFFInt, IFFNul, JJ
      Data IFFEnd/3hEND/, IFFStr/3hSTR/, IFFInt/3hINT/, IFFNul/3hNUL/,
     $  JJ/0/
 1010 Format(/,' Charge, Multiplicity card seems defective:')
 1020 Format(' Multiplicity of',I10,' is bogus.')
 1030 Format(' Charge is bogus.')
 1040 Format(' Multiplicity is bogus.')
 1050 Format(' The maximum of',I5,' Z-matrix cards has been exceeded.')
C
C     Initialize the output routine, and some flags.
C
      If(PrtInp) Call SzPrnt(0,JJ,0,0)
      NZ = 0
      ISymb = 0
      NSB = 0
      NCur = 0
      LastRs = 0
      IRes = 0
C
C     Read the charge and multiplicity, if so requested.
C
      If(RdChar) then
        If(RdInpt) then
          Call FFRead(IEof)
        else
          Call FFRMem(BB,NBB,IEOF)
          endIf
        If(IEOF.ne.0) Call GauErr('End of file in ZSymb.')
        Call FFGet(Str,Len,ICharg,FP,IF1)
        Call FFGet(Str,Len,Multip,FP,IF2)
        If(If1.eq.IFFNul) then
          ICharg = 0
          If1 = IFFInt
          endIf
        If(IF1.eq.IFFInt.and.IF2.eq.IFFInt.and.Multip.gt.0) then
          If(PrtInp) then
            Call SzPrnt(1,ICharg,1,0)
            Call SzPrnt(1,Multip,1,0)
            Call SzPrnt(2,JJ,0,0)
            endIf
        else
          Write(IOut,1010)
          If(IF1.ne.IFFInt) then
            Write(IOut,1030)
            Call FFErr(IFFINT,IF1)
            endIf
          If(IF2.ne.IFFInt) then
            Write(IOut,1040)
            Call FFErr(IFFINT,IF2)
            endIf
          If(Multip.le.0) Write(IOut,1020) Multip
            Call Lnk1E(0)
          endIf
        endIf
C
C     Here is the top of the loop for reading z-matrix cards.
C     -------------------------------------------------------
   40 If(RdInpt) then
        Call FFRead(IEOF)
      else
        Call FFRMem(BB,NBB,IEOF)
        endIf
      If(IEOF.ne.0) Call GauErr('End of file in ZSymb.')
C
C     Get name of this center and make sure it's legal.
      Len = 0
      Call FFGetX(4,Str,Len,I,FP,IF1)
      If(IF1.eq.IFFEnd) then
        Goto 100
      else if(IF1.eq.IFFInt) then
        Call DecChr(I,Str,Len)
      else if(IF1.eq.IFFStr) then
        If(StEqHC(Str,'Variable',8)) goto 100
      else
        Call FFErr(IFFStr,IF1)
        endIf
      If(PrtInp) Call SzPrnt(1,Str,Len,20)
      Call PutHol(Str,Len,CStr)
      If(Len.eq.1.or.CStr(1:1).eq.'-'.or.CStr(2:2).eq.'-') then
        Len = 1
        Call PutCSt(' ',Str,Len)
        endIf
      NZ = NZ + 1
      First = NZ.eq.1
      If(NZ.gt.MaxNZ) then
        Write(IOut,1050) MaxNZ
        Call Lnk1E(0)
        endIf
C
C     Get atomic number, MM types, etc.
C
      If(IF1.eq.IFFInt) then
        IAnZ(NZ) = I
        LenX = Len
        AtZNuZ(NZ) = GFloat(I)
      else
        Call MSubst(IOut,InpFF,NumNam,MaxNam,NamLst,CStr,IAnZ(NZ),
     $    IAtTyZ(NZ),AtmChZ(NZ),IAtWtZ(NZ),AtmWtZ(NZ),IAtSpZ(NZ),
     $    AtZEfZ(NZ),AtQMoZ(NZ),AtGFaZ(NZ),IAtFrZ(NZ),IAtRsZ(NZ),
     $    IAtPDZ(NZ),IResNm,AtZNuZ(NZ),LenX)
C       Update residue table
        If(IAtRsZ(NZ).ne.LastRs) Call AddRes(.False.,MaxRes,NRes,
     $    IAtRsZ(NZ),IResNm,IResNa,IResNu,IRes)
        LastRs = IAtRsZ(NZ)
        IAtRsZ(NZ) = IRes
        endIf
C
C     Append this center to the list (NamCnt).
      Call PutB(Str,LenX,NamCnt,NCur)
      Call PutDel(2,NamCnt,NCur)
C
C     Get the center to which this atom is attached.
C
      Call ZCentr(IOut,PrtInp,First,.True.,IZ(NZ,1),NamCnt,NZ,FP)
      SymOK = (IZ(NZ,1).gt.0.and.NZ.ge.2).or.IZ(NZ,1).eq.0
      If(IZ(NZ,1).lt.0) BL(NZ) = FP
      If(SymOK) then
        Call ZParm(PrtInp,BL(NZ),LBL(NZ),Symbls,NSB,ISymb,.True.,SymOK)
      else
        Call PutCSt('0',Symbls,NSb)
        Call PutDel(2,Symbls,NSb)
        endIf
      If(IZ(NZ,1).eq.-999999) IZ(NZ,1) = 0
      SymOK = IZ(NZ,1).ge.0
C
C     Get the 3rd center and the bond angle/y-coordinate.
      If(NZ.ge.3.and.IZ(NZ,1).gt.0)
     $  Call ZCentr(IOut,PrtInp,.False.,.False.,IZ(NZ,2),NamCnt,NZ,FP)
      If(NZ.ge.3.or.IZ(NZ,1).le.0) then
        Call ZParm(PrtInp,Alpha(NZ),LAlpha(NZ),Symbls,NSB,ISymb,.True.,
     $    SymOK)
      else
        Call PutCSt('0',Symbls,NSb)
        Call PutDel(2,Symbls,NSb)
        endIf
C
C     Get the 4th center and the dihedral angle/z-coordinate.
      If(NZ.ge.4.and.IZ(NZ,1).gt.0)
     $  Call ZCentr(IOut,PrtInp,.False.,.False.,IZ(NZ,3),NamCnt,NZ,FP)
      If(NZ.ge.4.or.IZ(NZ,1).le.0) then
        Call ZParm(PrtInp,Beta(NZ),LBeta(NZ),Symbls,NSB,ISymb,.True.,
     $    SymOK)
      else
        Call PutCSt('0',Symbls,NSb)
        Call PutDel(2,Symbls,NSb)
        endIf
C
C     Get the final integer
      If(IZ(NZ,1).gt.0.and.NZ.ge.4) then
        Call FFGet(Str,Len,I,FP,IF1)
        If(IF1.eq.IFFEnd) then
          IZ(NZ,4) = 0
        else if(IF1.eq.IFFInt) then
          IZ(NZ,4) = I
        else
          Call FFErr(IFFInt,IF1)
          endIf
        If(PrtInp) Call SzPrnt(1,IZ(NZ,4),0,1)
        endIf
      Call ZSyOni(IOut,InpFF,NumNam,MaxNam,NamLst,NZ,.True.,MaxNZ,IZ,
     $  Layer,IAnMod,LiAtom,MaxLay,MaxNZ,NuOfLa,Distan,ACMod,IAtMod,
     $  PrtInp)
      If(PrtInp) Call SzPrnt(2,JJ,0,0)
      Goto 40
C
  100 Call PutDel(5,Symbls,NSb)
      Return
      End
*Deck ZSyOni
      Subroutine ZSyOni(IOut,InpFF,NumNam,MaxNam,NamLst,NZ,UseIZ,MaxNZ,
     $  IZ,Layer,IAnMod,LiAtom,MaxLay,NDist,NuOfLa,Distan,ACMod,IAtMod,
     $  PrtInp)
      Implicit Real*8(A-H,O-Z)
C
C     Fetch ONIOM input from z-matrix line.
C
      Parameter (LStr=10)
      Integer IZ(MaxNZ,4), Layer(*), IAnMod(*), LiAtom(*), IAtMod(*),
     $  NamLst(2,*), Str(LStr), JJ(1)
      Dimension Distan(NDist,*), ACMod(*)
      Character*80 CStr
      Logical StEqHC, PrtInp, UseIZ
      Save Zero, IFFEnd, IFFStr, IFFInt, IFFFP, JJ
      Data Zero/0.0d0/, IFFEnd/3hEND/, IFFStr/3hSTR/, IFFInt/3hINT/,
     $  IFFFP/3hFP /, JJ/0/
 1000 Format(' Link atom not specified for center',I5,'.')
C
      Len = 0
      Call CBlank(LStr,Str)
      Call FFGet(Str,Len,I,FP,IF1)
C
C     Layer
C
      Layer(NZ) = 1
      IAnMod(NZ) = 0
      LiAtom(NZ) = 0
      Do 5 II = 1, MaxLay
    5   Distan(NZ,II) = Zero
      ACMod(NZ) = Zero
      IAtMod(NZ) = 0
      If(IF1.eq.IFFEnd) then
      else
        If(IF1.eq.IFFStr) then
          If(StEqHC(Str,'High',       Len)) then
            Layer(NZ) = 1
          else if(StEqHC(Str,'Medium',Len)) then
            Layer(NZ) = 2
          else if(StEqHC(Str,'Low',   Len)) then
            Layer(NZ) = 3
          else
            Call IniStr
            Call BldStC('Unrecognized layer "')
            Call BldStr(Str,Len)
            Call BldStC('".')
            Call DmpStr(IOut)
            If(PrtInp) Call SzPrnt(2,JJ,0,0)
            Call Lnk1E(0)
            endIf
        else if(IF1.eq.IFFInt) then
          Layer(NZ) = I
        else
          Call FFErr(IFFStr,IF1)
          endIf
        If(PrtInp) Call SzPrnt(1,Str,2,4)
        Call CBlank(LStr,Str)
        Call FFGetX(1,Str,Len,I,FP,IF1)
        LenX = Len
C       Model atom
        If(IF1.eq.IFFEnd) Return
        If(IF1.ne.IFFStr.and.IF1.ne.IFFInt) Call FFErr(IFFStr,IF1)
        If(PrtInp) Call SzPrnt(1,Str,Len,3)
        If(IF1.eq.IFFInt) then
          IAnMod(NZ) = I
          Call FFGet(Str,Len,I,FP,IF2)
          If(IF2.eq.IFFInt) then
            IAtMod(NZ) = I
            Call SzPrnt(1,I,0,1)
          else if(IF2.ne.IFFEnd) then
            Call FFErr(IFFInt,IF2)
            endIf
          If(IF1.ne.IFFEnd) then
            Call FFGet(Str,Len,I,FP,IF2)
            If(IF2.eq.IFFFP) then
              ACMod(NZ) = FP
              Call SzPrnt(1,FP,0,2)
            else if(IF2.ne.IFFEnd) then
              Call FFErr(IFFFP,IF2)
              endIf
            endIf
        else
          Call PutHol(Str,Len,CStr)
          If(Len.ge.3.and.StEqHC(Str,'NoMicro',Len)) then
            IAnMod(NZ) = -1
          else
            Call MSubst(IOut,InpFF,NumNam,MaxNam,NamLst,CStr,IAnMod(NZ),
     $        IAtMod(NZ),ACMod(NZ),IJunk1,RJunk1,IJunk2,RJunk2,RJunk3,
     $        RJunk4,IJunk3,IJunk4,IJunk5,IJunk6,RJunk5,LenX)
            endIf
          endIf
        Call Pad(Str,LenX,2,IBlank)
C       Link atom, defaulted to bonded atom in z-matrix.
        Call FFGet(Str,Len,I,FP,IF1)
        If(IF1.eq.IFFInt) then
          LiAtom(NZ) = I
          Call FFGet(Str,Len,I,FP,IF1)
        else if(IF1.eq.IFFEnd) then
          If(UseIZ) LiAtom(NZ) = Max(IZ(NZ,1),0)
        else
          Write(IOut,1000) NZ
          If(PrtInp) Call SzPrnt(2,JJ,0,0)
          Call Lnk1E(0)
          endIf
        If(PrtInp) Call SzPrnt(1,LiAtom(NZ),0,1)
C       Distance
        Do 10 II = 1, NuOfLa
          If(IF1.eq.IFFFP) then
            Distan(NZ,II) = FP
          else if(IF1.eq.IFFEnd) then
            Distan(NZ,II) = Zero
          else
            Call FFErr(IFFFP,IF1)
            endIf
          If(IF1.ne.IFFEnd) Call FFGet(Str,Len,I,FP,IF1)
   10     Continue
        If(PrtInp) then
          Do 20 II = 1, NuOfLa
   20       Call SzPrnt(1,Distan(NZ,II),0,2)
          endif
        endIf
      Return
      End
