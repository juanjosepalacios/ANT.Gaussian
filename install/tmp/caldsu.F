*Deck CalDSu
      Subroutine CalDSu(IOut,IPrint,IPFlag,NFrqR1,NMat0,ICnBeg,ICnEnd,
     $  DoSpar,Ex,Ec,IIV,VA,VB,D1E,D2E,VxA,VxB,VxAI,VxBI,ICntrl,IExchn,
     $  ICorr,IExCoW,IRadAI,IRanWt,IRanGd,ICorTp,ScaDFX,C,IAn,IAtTyp,
     $  AtmChg,NAtoms,NAtCel,IP,PA,D1PA,GA,PB,D1PB,GB,NMOA,CMOA,NMOB,
     $  CMOB,NBas6D,IOpCl,AccDes,UseSym,DoSymF,NOpSet,IChHar,MulHar,
     $  NMtPBC,NClRep,IntCel,BckCel,BckDim,CntCel,IPrtP1,FreqP1,NE,
     $  IGWInf,IAtBtD,IRdBtD,RRdBtD,RGWBtD,Omega,AtNetC,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     CalDSu is the top-level interface for density functional
C     calculations.  It handles loading common blocks and controlling
C     parallel execution.  However, input matrices should already be
C     cartesian.  /B/ and /RepAll/ should already be loaded.
C
C     Input:
C            IOut   - Output unit
C            IPrint - Print option
C            IPFlag - Usual direct control flag bits:
C                     2 -- Reverse normal choice of whether to precompute
C                          and store the distance matrix
C                     3 -- Skip consistency checks during XC quadrature.
C                     4 -- Do no do extra work to use cutoffs better.
C                     7 -- Force one matrix at a time code in CPKS.
C                    23 -- Turn off use of Gaussian expansions of atomic
C                          densities
C                    25 -- Allocate for parallel operation but run
C                          sequentially
C                    26 -- Make all atoms large
C                    27 -- Make all shells large
C                    28 -- Do not use symmetry on atomic grid points
C                    29 -- Use precomputed XC weights
C            NFrqR1 - Number of frequencies for which derivative densities
C                     are provided if ICntr1=5.
C            NMat0  - Number of input densities or input density
C                     derivatives if doing CPKS, polar derivs, etc.
C            ICnBeg, ICnEnd - Range of atoms for which to generate
C                             F(x) terms.
C            ICntrl - Controls what is computed
C                        1 = XC energy
C                        2 = XC contribution to hyperpolarizability.
C                       10 = XC potential matrix
C                       20 = XC contrib to TD-KS Lagrangian
C                       30 = XC contrib to TD-KS Lagrangian, and also
C                            non-separable contrib to TD-KS gradient
C                       40 = separable XC contrib to post-KS gradient
C                       50 = XC contrib to polarizability derivatives.
C                       60 = XC contrib to GIAO Gmx(P)
C                       70 = XC contrib to GIAO Gm(Px)
C                       80 = XC contrib to GIAO Gm(Pe)
C                       90 = XC contrib to GIAO <Pe*Gmx(P)>
C                      100 = XC 1st derivs wrt nuclear centers
C                      200 = XC 2nd derivs wrt nuclear centers
C                      300 = XC contrib to CPKS
C                      400 = XC contrib to nuclear F(x)
C                      500 = XC contrib to nuclear F(x), and also
C                            forces
C                      600 = XC contrib to GIAO F(x)
C                      700 = XC 2nd derivs wrt magnetic field (GIAOs)
C                      800 = XC Tau contrib to CSGT/IGAIM F(x)
C                      900 = XC contrib to nuclear F(x) and also XC
C                            2nd derivs wrt nuclear centers.
C                     1000 = XC contrib to CPKS(x)
C                     1100 = XC contrib to Gx(P)Tab for polar derivatives.
C                    10000 = Turn off weight derivatives.
C                    20000 = Turn off radial shell pruning.
C                   000000 = Default, same as 400000 if NMOA is
C                            non-zero, same as 100000 otherwise.
C                   100000 = Force use of density matrix.
C                   200000 = Force use of MOs via single matrix multiply
C                   300000 = Force use of MOs built with cutoffs
C                   400000 = Decide on the fly how to compute density
C                            values.  Both PA/PB and CMOA/CMOB must be
C                            provided.
C                   500000 = Use superposition of atomic densities..
C                   600000 = Test superposition of atomic densities.
C                   700000 = Use DBF expansion in /B2/.
C                   800000 = Use superposition of interpolated atomic
C                            densities.  Net charges are in PA.
C                   900000 = Use superposition of old atomic densities.
C                  0000000 = Normal evaluation of functionals.
C                  1000000 = Harris functional, including superposition
C                            of Coulomb terms.
C                  2000000 = Harris, but include full nuclear charge
C                            (initial guess corresponding to DoV=.False.
C                            in 1e integrals).
C                  3000000 = Harris without Coulomb.
C                  4000000 = Do 1st order approximation to Exc.
C                  5000000 = Do 2nd order approximation to Exc.
C                  6000000 = Return Hirshfeld populations in VA.
C                  7000000 = Same as 3000000 but with large atomic
C                            density expansions
C                 10000000 = When doing grid on an atom, include only basis
C                            functions and potential of that atom.
C                100000000 = Do numerical differentiation of the functional
C                            and compare to analytic derivatives. Requires
C                            extra input ... see routine NDXCIn. Strictly
C                            for debugging purposes.
C            IExchn - Exchange method
C                        2 = Hartree-Fock-Slater
C                        3 = X-Alpha, Alpha = 0.7
C                        4 = Becke 88
C                        5 = Exchange 92
C                        6 = PW91
C                        7 = Gill 96
C                        8 = Perdew 86
C                        9 = mPW
C                       10 = PBE
C                       11 = BA
C                       12 = VSXC
C            ICorr  - Correlation method
C                        0 = None
C                        1 = VWN5
C                        2 = Lee, Yang and Parr
C                        3 = Perdew 81
C                        4 = Perdew 81 + Perdew 86
C                        5 = VWN
C                        6 = VWN + Perdew 86
C                        7 = Odom-Scuseria
C                        8 = PW91
C                        9 = PBE
C                       10 = VSXC
C                       11 = Bc96
C            IExCoW - Frequency-dependence of Xc functional; only
C                     relevent if doing CPKS.
C                      0,1 = None, static limit.
C                        2 = Gross-Kohn form.
C            IRadAI - mmmnnn, where mmm = number of radial
C                        quadrature points, nnn = number of
C                        angular quadrature points
C                     If mmm is zero, then a special pruned grid
C                     is requested.  Current pruned grids:
C                        0 ... same as 1.
C                        1 ... original pruned (50,194) (SG-1)
C                        2 ... pruned (30,72) grid for 1 kcal
C                              single-point accuracy
C                        3 ... sleazier pruned (30,72) for SCF pass 0.
C            ScaDFX - Four scale factors for local exchange, non-local
C                     exchange, local correlation, and non-local
C                     correlation.
C            C      - Atomic coordinates
C            IAN    - Atomic numbers
C            NAtoms - Number of atoms
C            PA     - Alpha density matrix
C            D1PA   - First derivs of alpha density matrix
C            GA     - Alpha post-KS density matrix
C            PB     - Beta  density matrix
C            D1PB   - First derivs of beta density matrix
C            GB     - Beta post-KS density matrix
C            NMOA   - Number of alpha MOs provided in CMOA.  See ICntrl.
C            CMOA   - Alpha MOs.
C            NMOB   - Number of beta MOs provided in CMOA.  See ICntrl.
C            CMOB   - Beta MOs.
C            NBas6D - Number of cartesian basis functions
C            IOpCl  - Type of SCF.
C            AccDes - Desired accuracy
C            UseSym - Whether to try to use symmetry.  If .true., this
C                     routine loads all necessary information from the rwf.
C            NOpSet - Number of symmetry operations to use, in case output
C                     quantities are being incremented and contain terms
C                     already computed at a certain level of symmetry.
C            IPrtP1 - Which atom (if any) is displaced for the perturbation
C                     in each D1PA,D1PB. Should be 0 for any perturbations
C                     which are not atomic displacements.
C            FreqP1 - Frequency of applied perturbations when doing CPKS.
C                     Only used for CPKS and IExCoW > 1.
C            NE     - Number of electrons, used for check on integration
C                     accuracy.
C            V      - Working precision mega-array
C            MDV    - Length of V
C                        (if 0, all available memory is used)
C
C     Output:
C            Ex     - Exchange energies
C            Ec     - Correlation energies
C            VA     - Alpha exchange-correlation matrix
C            VB     - Beta  exchange-correlation matrix
C                     Note: XC contribution to F2 (wrt electric field)
C                     which involves third derivatives of the functional
C                     is returned in VA,VB.
C            D1E    - XC 1st derivatives wrt nuclear centers,
C                     not initialized.
C                     Also used for XC direct contribution to
C                     polarizability derivatives.
C            D2E    - XC 2nd derivatives wrt nuclear centers,
C                     or hyperpolarizabilities, not initialized.
C            VxA    - Alpha XC 1st derivative matrices
C            VxB    - Beta  XC 1st derivative matrices
C            VxAI   - Alpha XC derivative matrices, contribution from
C                     imaginary XC response.
C            VxBI   - BetaXC derivative matrices, contribution from
C                     imaginary XC response.
C
C     Refer to comments in routine CalDFT for implementation details.
C
#include "commonb2.inc"
#include "repall.inc"
      Logical DoE, DoVAB, DoD1E, DoD2E, DoCPKS, DoNuFx, UseSym, GotOne,
     $  LStat, AbOnly, DbgLin, TrcLin, LBit, UseB1, DoZCmp, Do1Mat,
     $  DoMgFx, DoFx, DoSymF, DoSpar, UseB2, DoTStm, DoDyPL, DoGxPT,
     $  WrkOMP, AllFkD, DoSuper, DoRho, DoHyp, DoJig, HavOCBf, DoWrt,
     $  DoChi, DubChi, BigAtm, BigShl, ThrOK, SymAtG, DubWgt, DoMicB,
     $  Fatal, Fail, OldAtD, SeqLin, CSGTau, BareAt, ChekEC, AtBlock,
     $  TDKSLg, PKSDen, DoPolD, LinDyn, CBfn, Spinor, DoGmx, DoHir,
     $  DoGmPx, DoGmPe, PeGmx, DoNDXC, DoVarX, DoVarG, DoDPD, DoPrnt,
     $  DoVN2, UseANC, OldDBF, BigDBE, DoHar
      Integer DoDump, AOut
      Character*(*) CountN
      Parameter (IONEqS=565,IONEqB=580,LMax=DEFLMAX,LenInR=(2*LMax+2),
     $  DoDump=0,MapTyp=1,AOut=1,ITpAtD=3,IType=3,MxDnXC=DEFMXDNXC,
     $  MxTyXC=DEFMXTYXC,MxVrXC=(MxDnXC*MxTyXC),MxXCNm=DEFMAXXCNAME,
     $  IDoDPD=1,CountN='Prism/CalDFT -- NxtVal')
      External NProc, LindEv
      Common /SymInf/ NOp1,NOp2,JTrans(3,8),T(3,3),TrVec(3)
      Dimension VA(*), VB(*), D1E(*), D2E(*), VxA(*), VxB(*), C(*),
     $  IAn(*), PA(*), D1PA(*), GA(*), PB(*), D1PB(*), GB(*), CMOA(*),
     $  CMOB(*), V(*), ScaDFX(*), IPrtP1(*), AtmChg(*), JJ(1), IIV(*),
     $  IP(*), IAtTyp(*), FreqP1(*), VxAI(*), VXBI(*), CntCel(3,NClRep),
     $  IGWInf(*), IAtBtD(*), IRdBtD(*), RRdBtD(*), RGWBtD(*),
     $  IndR(LenInR), Omega(*), IVDen(MxVrXC), IVTyp(MxVrXC),
     $  IVLoc(MxVrXC), IVC2T(MxVrXC), IVT2C(MxVrXC), IUDen(MxVrXC),
     $  IUTyp(MxVrXC), IULoc(MxVrXC), IUC2T(MxVrXC), IUT2C(MxVrXC),
     $  ZTT(3), AtNetC(*)
      Integer IntCel(3,*), BckCel(*), BckDim(3), GLinCo
      Character ExName*(MxXCNm), CoName*(MxXCNm), FName*8
      Real*8 MDCutO
      Save Zero, Pt5, One, Two, JJ, MaxV, SizInc, MaxMBI, DoWrt, NCall,
     $  ECTol2
      Data Zero/0.0d0/, Pt5/0.5d0/, One/1.0d0/, Two/2.0d0/, JJ/0/,
     $  MaxV/0/, SizInc/5.0d0/, MaxMBI/0/, DoWrt/.True./, NCall/0/,
     $  ECTol2/0.1d0/
 1000 Format(' CalDSu:  requested number of processors reduced to:',i4,
     $  ' ShMem',i4,' Linda.')
 1010 Format(' Spurious integrated density or basis function:')
 1015 Format(' Requested info about the electron count:')
 1030 Format(' CalDSu exits because no D1Ps are significant.')
 1040 Format(' CalDSu:  NOpAll=',I3,' NOpUse=',I3,' but NOpSet=',I3,'.')
 1050 Format(' NE=',I5,' NElCor=',I5,' El error=',1PD8.2,' rel=',1PD8.2,
     $  ' Tolerance=',1PD8.2,/,' Shell',I6,'     absolute error=',
     $  1PD8.2,'              Tolerance=',1PD8.2,/,' Shell',I6,
     $  '       signed error=',1PD8.2,'              Tolerance=',1PD8.2)
 1055 Format(' NE=',I5,' NElCor=',I5,' El error=',1PD8.2,' rel=',1PD8.2,
     $  ' Tolerance=',1PD8.2,/,' Shell',I6,'     absolute error=',
     $  1PD8.2,'              Tolerance=no limit',/,' Shell',I6,
     $  '       signed error=',1PD8.2,
     $  '              Tolerance=no limit')
 1070 Format(' CalDSu:  MaxMB=',I4,
     $  ' is small; additional memory will improve performance.')
 1080 Format(' CalDSu:  NPrtUS=',I4,' ThrOK=',L1,' IAlg=',I1,' NPAlg=',
     $  I1,' DoDPD=',L1,' LenP=',I12,' LenD1P=',I14,'.')
 1090 Format('c',I1,'p',I1,'.dat')
 1100 Format(' Not enough memory to run CalDSu, short by',I12,' words.')
 1110 Format(' CalDSu would require',I12,' more words of memory to run',
     $  ' efficiently.')
 1120 Format(' CalDSu: Zero Torque Theorem, X=',D13.6,' Y=',D13.6,' Z=',
     $  D13.6)
C
      Call TStamp(1,'Top of CalDSu')
      DoPrnt = IPrint.gt.0.or.DoTStm(0)
      Call DecoSC(IOpCl,NSpBlk,NRI,NDimBl,CBfn,NSpBlX,Spinor)
      If(NDimBl.eq.2) Call GauErr('GKS NYI in CalDSu.')
      IOpClX = NSpBlX - 1
      NCall = NCall + 1
      Fatal = .not.LBit(IPFlag,3)
      DoZCmp = .not.LBit(IPFlag,4)
      Do1Mat = LBit(IPFlag,7)
      DoJig = .False.
      IRType = IRadDf(IPFlag,0)
      ThrOK = .not.LBit(IPFlag,25)
      Call UnPkCD(ICntrl,DoE,DoVAB,DoD1E,DoD2E,DoCPKS,DoNuFx,DoMgFx,
     $  CSGTau,DoHyp,TDKSLg,PKSDen,DoChi,DoPolD,DoGmx,DoGmPx,DoGmPe,
     $  PeGmx,IDerNu,DoSuper,DoRho,DoGxPT,DoHir,DoVN2,DoHar,DoVarG,
     $  DoVarX)
      ICntr4 = Mod(ICntrl,100000)/10000
      ICntr5 = Mod(ICntrl,1000000)/100000
      ICntr6 = Mod(ICntrl,10000000)/1000000
      ICntr7 = Mod(ICntrl,100000000)/10000000
      ICntr8 = Mod(ICntrl,1000000000)/100000000
      DoNDXC = ICntr8.eq.1
      DoMicB = IAbs(IRanGd).lt.1000000000
      AtBlock = ICntr7.eq.1
      NAtomA = NAtAct(MapTyp,NAtoms,IAtTyp,IAn)
      If(DoHir) then
        NRhoS = NAtomA
      else if(DoSuper) then
        NRhoS = 1
      else
        NRhoS = 0
        endIf
      Call LSetup(IType,IPrint,.True.,IPFlag,1,V,SeqLin,DbgLin,TrcLin,
     $  NPSMem,NPLind,NPSt,NPSSav,NPLSav,LStat,MemWrk,WrkOMP)
      If(DoSpar) NPLind = 1
      If(NPLind.gt.1.and..not.LStat)
     $  Call GauErr('NPLind>1 but not LStat in CalDSu.')
      If(DoGmx) then
        NMtFkD = 9*(ICnEnd-ICnBeg+1)
      else if(DoNuFx) then
        NMtFkD = 3*(ICnEnd-ICnBeg+1)
      else if(CSGTau) then
        NMtFkD = 6
      else if(DoMgFx) then
        NMtFkD = 3
      else
        NMtFkD = 0
        endIf
      DoFx = DoNuFx.or.DoMgFx
      NZV = 0
      If(DoSpar) then
        If(DoVAB) NZV = IIV(NBas6D+1) - 1
        NZP = IP(NBas6D+1) - 1
        MemSpr = InToWP(NBas6D+1+NZP) + InToWP(NZP)
      else
        MemSpr = 0
        endIf
      NTT6D = (NBas6D*(NBas6D+1))/2
      Call DecPrn(IOut,IPrint,IRadAI,0,IPrune,NTheta,NPhi,MaxRad,MaxAng,
     $  ECTolR)
      If(IPrune.ne.0) then
        IRadAn = IPrune
      else
        IRadAn = IRadAI
        endIf
      Thresh = ThrDFT(IRadAn,AccDes,NClRep)
      ISavGI = -1
      If(LBit(IPFlag,29)) ISavGI = 0
      Call CkSvGd(IOut,IPrint,ISavGI,NAtomA,IRadAn,IRanWt,IRanGd,IRType,
     $  NClRep,ICntr4,IPFlag,IGWInf,ISavGd)
      NMat1 = NDimBl*NMat0
      NMat = NMat1*NMtPBC
      If(NAtCel.ne.0) then
        NAtomX = NAtCel
      else
        NAtomX = NAtomA
        endIf
      NMtUse = Max(NMtPBC,NClRep)
      Call NDerXC(IOut,IPrint,IOpCl,IExchn,ICorr,IDerNu,DoNuFx,MetExc,
     $  MetCor,ExName,CoName,Alpha,IfTau,IfLap,MaxDrW,MaxDrB,MaxDrR,
     $  IDRhTy,MxDnXC,MxTyXC,IVDen,IVTyp,IVC2T,IVT2C,IVLoc,NVarVT,
     $  NVarVC,IUDen,IUTyp,IUC2T,IUT2C,IULoc,NVarUT,NVarUC)
C     Use precomputed max number of significant basis functions even
C     if it was for a different grid.
      NBsAll = NBas6D*NClRep
      NAtAll = NAtomA
      If(IGWInf(1).gt.0.and..not.DoSpar) then
        If(DoSuper) then
          IFacNB = 18
        else if(MaxDrW.eq.0) then
          IFacNB = 18
        else
          IFacNB = 18
          endIf
        IFacNA = 14
        NBsAll = Min((IFacNB*IGWInf(18+MaxDrB))/10,NBsAll)
        NAtAll = Min((IFacNA*IGWInf(24))/10,NAtAll)
        endIf
      If(DoGmPe) NAtAll = 0
C
C     Load symmetry information.
C
      If(UseSym.and..not.DoSymF) then
        NOpUse = NOpAll
      else if(UseSym.and..not.DoCPKS.and..not.DoMgFx.and.(.not.DoSymF
     $  .or..not.(DoNuFx.and..not.AllFkD(ICnBeg,ICnEnd,NAtoms)))
     $  .and..not.DoGmPe.and..not.DoGmPx) then
        NOpUse = NOpAll
      else
        NOpUse = 1
        endIf
      If(AbOnly(0)) NOpUse = Min(NOpUse,NOp2)
      If(DoSuper.and.NOpUse.gt.0)
     $  Call CkAtEq(NOpAll,MxAtSO,NAtoms,NEqAll,IAn,IAtTyp,NOp2,NOpUse)
      If(UseSym.and.NOpSet.ne.0) then
        If(NOpUse.lt.NOpSet) then
          Write(IOut,1040) NOpAll, NOpUse, NOpSet
          Call GauErr('CalDSu cannot use requested symmetry level.')
        else
          NOpUse = NOpSet
          endIf
        endIf
      MaxTyX = MaxTyp
      UseB2 = ICntr5.eq.7
      If(UseB2.and..not.DoSuper) MaxTyX = Max(MaxTyX,MaxTyB)
      MMax = ((MaxTyX+1)*(MaxTyX+2)*(MaxTyX+3))/6
      ThrshX = Thresh / GFloat(100)
      NDisMt = NQDisM(MaxDrW,IRanWt,NAtCel,IPFlag)
      OldAtD = LBit(IPFlag,23)
      UseANC = ICntr5.eq.8
      OldDBF = ICntr5.eq.9.or.LBit(IPFlag,31)
      BigDBE = ICntr6.eq.4.or.ICntr6.eq.5.or.ICntr6.eq.7
      If(ICorTp.gt.0) then
        IAtCor = 3
      else
        IAtCor = 1
        endIf
      NElR = NE
      If(NAtCel.ne.0) NElR = NElR*(NAtoms/NAtCel)
      BigAtm = LBit(IPFlag,26)
      BigShl = LBit(IPFlag,27)
      SymAtG = .not.LBit(IPFlag,28)
C     Use 3rd derivs of basis functions here for consistency among
C     energies, gradients, and frequencies.
      MaxDBX = 4
      jStart = 1
      Call Quad0(IOut,IPrint,UseB2,DoSuper,DoHir,OldAtD,UseANC,OldDBF,
     $  BigDBE,ITpAtD,IAtCor,BigAtm,BigShl,NDisMt,MaxDBX,MapTyp,NAtoms,
     $  NAtomA,NClRep,NElR,IChHar,IAn,IAtTyp,AtmChg,C,AtNetC,ThrshX,1,
     $  jStart,HavOCBf,BareAt,MaxNG,MaxAnT,NElCor,NShlDB,NDBF,NDBFU,
     $  IWhich,LenANC,ITypAD,IDoInA,IDoAtm,IMpAtm,IMpAt1,IAnMap,IAtTMp,
     $  IAtCMp,ICMap,IMpJAn,IBegSh,jAtmSz,jShlSz,jDisMt,INAtG,jAtEC,
     $  jAtmSD,iBegSB,jShlSB,jMpJAB,Next,V,MDV)
      IIJV = Next
      IIJKV = IIJV  + InToWP(NVarVT**2)
      IIJU = IIJKV + InToWP(NVarVT**3)
      IIJKU = IIJU  + InToWP(NVarUT**2)
      IBNSum = IIJKU + InToWP(NVarUT**3)
      iCnvrt = IBNSum + NBas6D
      iScale = iCnvrt + InToWP(MMax)
      IRL = iScale + MMax
      IV0 = IRL + MaxTyX + 1
      Call TstCor(IV0,MDV,'CalDSu-0')
      Call GenScl(.True.,MaxTyX,V(iRL),V(iScale))
      Call LdCnvr(1,.True.,MaxTyX,V(iCnvrt))
      Call FilInS(NVarVT,V(IIJV))
      Call FilIn3(NVarVT,V(IIJKV))
      Call FilInS(NVarUT,V(IIJU))
      Call FilIn3(NVarUT,V(IIJKU))
      NElTot = IArSum(1,NAtomA,V(IAnMap))
      UseB1 = .True.
      ICntrX = ICntrl
      If(ICntr5.eq.0) then
        If(IDerNu.eq.0.and..not.DoNuFx.and.NMOA.gt.0.and.ICorTp.eq.0
     $    .and..not.DoMgFx.and..not.DoChi.and..not.DoCPKS.and.
     $    .not.DoSpar.and.NClRep.eq.1.and..not.DoGmPe.and.
     $    .not.DoGmPx) then
          ICntr5 = 4
        else
          ICntr5 = 1
          endIf
        ICntrX = ICntrX + ICntr5*100000
        endIf
      DubChi = ICntr5.eq.2.or.ICntr5.eq.3.or.ICntr5.eq.4
      If(DoD2E) then
        LLT2D  = (3*NAtomA+3)**2
      else
        LLT2D  = 0
        endIf
      iFacAt = IV0
      iLT2D  = iFacAt + NAtoms
      iLT2D1 = iLT2D  + InToWP(LLT2D)
      IV0 = iLT2D1 + InToWP(LLT2D)
      Call TstCor(IV0,MDV,'CalDSu-1')
      ChekEC = .not.(AtBlock.and.BareAt)
      If(LLT2D.gt.0) then
        Call MkLT2D(.True.,0,NAtomA,V(iMpAtm),V(iLT2D))
        Call MkLT2D(.False.,1,NAtomA,V(iMpAtm),V(iLT2D1))
        endIf
      If(SymAtG) then
        NOpAtG = Min(NOpUse,NOp2)
      else
        NOpAtG = 1
        endIf
      LNEqAt = NAtoms*NOpUse
      If(NOpUse.gt.1) then
        Call PcckI(0,Junk,NEqAll,MxAtSO,NOpAll,NAtoms,NOpAll)
        Call PetDFT(V(iDoAtm),V(iFacAt),NEqAll,NAtoms,NOpUse,NOpAtG)
      else
        Call ASet(NAtoms,One,V(iFacAt))
        endIf
      If(NAtCel.ne.0) then
        Call LClear(NAtoms,V(iDoAtm))
        Call LSet(NAtCel,V(iDoAtm))
        Call ASet(NAtoms,One,V(iFacAt))
        endIf
C
C     When doing GKS transform the densities from (A,B,X,Y) to (T,Z,X,Y) and
C     scale (X,Y).
C
      If(NDimBl.eq.2) then
        Call SumDif(.True.,NTT6D,PA,PB)
        Call AScale(NTT6D,Two,PA(1+NTT6D),PA(1+NTT6D))
        Call AScale(NTT6D,Two,PB(1+NTT6D),PB(1+NTT6D))
        endIf
      NMatCP = 0
      NMatG = 0
      LenP = NTT6D*NMat
      If(DoVarX) then
        LenP = NTT6D*NDimBl
        NMatCP = NMat0
        Thrsh1 = Thresh / GFloat(100)
        iIPrP1 = IV0
        IEnd = iIPrP1 + NMatCP - 1
        Call TstCor(IEnd,MDV,'CalDSu-1a')
        Call MapInt(NMatCP,V(IMpAt1),IPrtP1,V(iIPrP1))
        If(NDimBl.eq.2)
     $    Call GauErr('GKS and density derivatives in CalDSu.')
        Call CkPSig(NSpBlX,NDimBl*NMtPBC,NMatCP,NTT6D,Thrsh1,D1PA,
     $    D1PB,GotOne,V(iIPrP1))
        If(.not.GotOne) then
          Write(IOut,1030)
          Goto 999
          endIf
        endIf
      If(DoVarG) then
        LenG = NTT6D*NDimBl
        NMatG = 1
        If(NDimBl.eq.2)
     $    Call GauErr('GKS and post-KS density in CalDSu.')
      else
        LenG = 0
        endIf
      If(DoGxPT) then
        NDPol = NMat0
        NMtCP1 = NMat0
      else if(DoPolD) then
        If(NFrqR1.lt.1) Call GauErr('Bad NFrqR1 in CalDSu.')
        NDPol = NFSDPo(NMat0/NFrqR1)*NFrqR1
        NMtCP1 = NDPol
      else
        NDPol = 0
        NMtCP1 = Max(NMatCP,NMatG)*NDimBl
        endIf
      If(DoGmPx) then
        NMtCP2 = 3*NMtCP1
      else
        NMtCP2 = NMtCP1
        endIf
      LenD1P = NTT6D*NMatCP
      NBfAll = NAllBf(MaxDrB,NBsAll,NShell,ShellT)
      Call DFTMem(IOut,IPrint,IOpCl,ICntrl,IExChn,ICorr,MaxRad,MaxAng,
     $  NShell,NAtomX,NMtUse,NMat1,NBas6D,NBsAll,NShlDB,NDBFU,Mem0,Mem1,
     $  Mem2,MemM0,MemM1,MemM2)
      If(DoD2E) then
        MemSet = Mem2
        MxDrW1 = 1
        MinPt  = MaxAng
      else if(DoD1E.or.DoNuFx.or.DoPolD.or.DoGxPT) then
        MemSet = Mem1
        If(NDisMt.eq.0) then
          MxDrW1 = 0
        else
          MxDrW1 = 1
          endIf
        MinPt  = MaxAng
      else
        MemSet = Mem0
        MxDrW1 = 0
        MinPt  = MaxAng
        endIf
      MemSet = MemSet + MemSpr
      MDV1 = MDV - IV0 + 1
      MxDCon = MxCont(IWhich)
      Call NFPSXC(IOpCl,NMtCP2,NBas6D,NBsAll,MaxFPS,NFPS,NFPSI,MaxFTr,
     $  NFPS2)
C     Trying to fit in cache fouls up too many other things, including
C     inner loop lengths and grid raking.
      MDC = 4*MDCach(0)
      MDC = 0
      MinVLn = MDVLen(2)
      DubWgt = DoD2E
C
C     Loop downward over possible number of processors, until we get
C     a number which fit.
C
      Do 400 NPart = NPSt, 1, -1
        If(NPart.gt.1.and.NPSMem.gt.1.and.NMtPBC.eq.1.and.IDoDPD.ge.1)
     $    then
          NPAlg = 2
        else
          NPAlg = 1
          endIf
        Do 390 IAlg = 1, NPAlg
          DoDPD = IAlg.lt.NPAlg
          Call AlCaSu(0,IOpCl,NAtoms,NAtomA,NMat0,NMat,NMtCP1,NMtFkD,
     $      NDPol,NBas6D,LenP,LenD1P,DoE,DoVAB,DoD1E,DoD2E,DoCPKS,DoHyp,
     $      DoChi,DoGmPx,DoGmPe,PeGmx,IExCoW,DoFx,DoDPD,DoHir,LenExc,
     $      LenVAB,LenD1E,LenD2E,LenVx,LenVxI,IElSum,IEx,IEc,IDA,IDB,
     $      IVA,IVB,ID1E,ID2E,IDPA,IDPB,IVxA,IVxB,IVxAI,IVxBI,IBNorm,
     $      IZTT,NZV,LenAll)
          NPartX = Min(NPart,NPSMem)
          If(NPLind.gt.1.and.NPart.gt.1) NPartX = Max(NPartX,2)
          LenVP = (MDV1-(NPartX-1)*LenAll)/NPartX
          MaxPtI = 0
          MaxMBI = 0
          LenAlg = LenVP - MemSet
          If(LenAlg.gt.0) then
            If(MDC.gt.0) Call AlgDDF(NVarVT,IVLoc,NVarUT,DoSpar,DubChi,
     $        IRanGd,MinPt,0,0,1,-Min(LenAlg,MDC),IDerNu,MaxDrB,MaxDrR,
     $        0,IDRhTy,MxDrW1,NBsAll,NFPS2,NAtomA,NAtAll,NMOA,NMat0,
     $        MxDCon,NBfAll,IOpCl,DoNuFx,DoCPKS,DoGmPe,DoGmPx,PeGmx,
     $        DubWgt,DoHyp,TDKSLg,PKSDen,NDPol,0,MaxFTR,NRhoS,0,DoRho,
     $        DoNDXC,DoVN2,MaxPtI,MaxMBI,jGrid,jWeigt,jD1Wgt,jD2Wgt,
     $        jChi,jD1Chi,jD2Chi,jD3Chi,jAlpha,jDCoef,jR,jExpT,jZeroM,
     $        jW,NColW,jVarVS,jVarUS,jVaVSN,jPhiS,jD1PhN,jVarV,jVarU,
     $        jVarVN,jF,jD1FV,jD2FV,jD3FV,jD1FU,jD2FU,jD3FU,jTRA,jTRB,
     $        jTDA,jTDB,jTBA,jTBB,jDXPr1,jDXPr2,jScr,jNIdx,NDBFU,jVarVX,
     $        jVaVXN,jVarVG,jVaVGN,jTTA,jCPA,jCPB,jD4Chi,jCQA,jTLA,LTBB,
     $        jIUV1,jUV1,jIUV2,jUV2,jIUV3,jUV3,jVVN2,jVVGN2,jVVSN2,
     $        LShort,NPerPt)
            If(MaxPtI.eq.0.or.MaxMBI.lt.MinVLn) Call AlgDDF(NVarVT,
     $        IVLoc,NVarUT,DoSpar,DubChi,IRanGd,MinPt,0,0,1,-LenAlg,
     $        IDerNu,MaxDrB,MaxDrR,0,IDRhTy,MxDrW1,NBsAll,NFPS2,NAtomA,
     $        NAtAll,NMOA,NMat0,MxDCon,NBfAll,IOpCl,DoNuFx,DoCPKS,
     $        DoGmPe,DoGmPx,PeGmx,DubWgt,DoHyp,TDKSLg,PKSDen,NDPol,0,
     $        MaxFTR,NRhoS,0,DoRho,DoNDXC,DoVN2,MaxPtI,MaxMBI,jGrid,
     $        jWeigt,jD1Wgt,jD2Wgt,jChi,jD1Chi,jD2Chi,jD3Chi,jAlpha,
     $        jDCoef,jR,jExpT,jZeroM,jW,NColW,jVarVS,jVarUS,jVaVSN,
     $        jPhiS,jD1PhN,jVarV,jVarU,jVarVN,jF,jD1FV,jD2FV,jD3FV,
     $        jD1FU,jD2FU,jD3FU,jTRA,jTRB,jTDA,jTDB,jTBA,jTBB,jDXPr1,
     $        jDXPr2,jScr,jNIdx,NDBFU,jVarVX,jVaVXN,jVarVG,jVaVGN,jTTA,
     $        jCPA,jCPB,jD4Chi,jCQA,jTLA,LTBB,jIUV1,jUV1,jIUV2,jUV2,
     $        jIUV3,jUV3,jVVN2,jVVGN2,jVVSN2,LShort,NPerPt)
          else
            NPerPt = 0
            LShort = 1000000 - LenAlg
            endIf
          If(NPart.eq.1.or.(MaxPtI.ge.MaxAng.and.MaxMBI.ge.MinVLn)) then
            If(LStat.and.NPart.gt.1) then
              NPrtUL = NPLind
            else
              NPrtUL = 1
              endIf
            If(NPart.eq.1) then
              If(MaxMBI.eq.0) then
                Write(IOut,1100) Max(LShort,20*NPerPt)
                Call Lnk1E(0)
              else if(MaxMBI.lt.MinVLn) then
                Write(IOut,1110) (MinVLn-MaxMBI)*NPerPt
                endIf
              endIf
            NPrtUS = Min(NPart,NPSMem)
            NPartT = NPrtUS*NPrtUL
            LenP1 = LenAll + LenVP
            If(NPart.lt.NPSt) Write(IOut,1000) NPrtUS, NPrtUL
            If(DoMicB) MaxMBI = Min(MaxMBI,2*MDVLen(0))
            If(MaxMBI.lt.MinVLn.and.IPrint.ge.0) Write(IOut,1070) MaxMBI
            LenRO = 12*NOpAll
            If(MaxV.gt.0) LenVP = Min(LenVP,MaxV)
            NxtVal = 1
            LinDyn = DoDyPL(IPFlag,NPrtUL)
C
C           Store input variables in tuple space if running with Linda.
C
            If(NPrtUL.gt.1) then
              Call EvLind(IOut,'LindEv',LindEv,1,NPrtUS,NPrtUL,NPSMem,
     $          WrkOMP,DbgLin,TrcLin,IType,MemWrk,LinDyn,NxtVal)
              NBits = IFBSet(.False.,NPrtUL-1)
              Do 170 IB = NBits, 0, -1
                IPartX = 2**IB
                If(IPartX.lt.NPrtUL) then
                  Call CDLin1(AOut,IOut,TrcLin,-IPartX,NPrtUL,UseB1,
     $              UseB2,IPrint,NFrqR1,NMat0,NMatCP,ICnBeg,ICnEnd,
     $              ICntrX,IExchn,ICorr,IRadAn,NAtoms,NAtomA,NMtPBC,
     $              NClRep,NBas6D,IOpCl,NOpAtg,LenRO,LNEqAt,NAtoms,NMOA,
     $              NMOB,ScaDFX,Thresh,DoZCmp,Do1Mat,IRanWt,IRanGd,
     $              MaxPtI,MinPt,NElCor,NDBF,IPFlag,MaxMBI,SizInc,
     $              IExCoW,IChHar,MulHar,DoJig,IRType,BckDim,ISavGd,
     $              NBsAll,JJ(1),IGWInf,NDisMt,NAtAll,ICorTp,LinDyn,
     $              ITypAD,DoDPD,NDPol,LenANC,NE)
                  Call CDLin2(AOut,IOut,TrcLin,-IPartX,SeqLin,IOpClX,
     $              DoRho,NAtoms,NAtomA,NBas6D,NMOA,NMOB,LenP,LenD1P,
     $              LenG,LenRO,LNEqAt,NMatCP,MaxTyX,NShlDB,LenANC,
     $              V(IAnMap),V(IAtTMp),V(IAtCMp),V(ICMap),PA,PB,D1PA,
     $              D1PB,GA,GB,RotAll,NEqAll,V(iDoAtm),V(iIPrP1),FreqP1,
     $              CMOA,CMOB,V(IFacAt),V(IMpAtm),V(IDoInA),NClRep,
     $              IntCel,BckCel,BckDim,CntCel,ISavGd,IGWInf,IAtBtD,
     $              IRdBtD,RRdBtD,RGWBtD,NShell,NDisMt,V(iMpJAn),
     $              V(jAtmSz),V(jShlSz),V(jDisMt),LLT2D,V(ILT2D),
     $              V(ILT2D1),V(iBegSh),V(iCnvrt),V(iScale),V(iBegSB),
     $              V(jShlSB),V(jMpJAB),Omega,AtNetC)
                  endIf
  170           Continue
              endIf
            ElSum = -One
            If(DoE) then
              Ex = Zero
              Ec = Zero
              endIf
            Call AClear(NBas6D,V(IBNSum))
            Call AClear(LenVAB,VA)
            If(.not.DoHir) Call AClear(LenVAB,VB)
            Call AClear(3,ZTT)
            If(ThrOK) then
              ThrOK = ThrOK.and.NPrtUS.gt.1
              IPart = NProc(1)
              If(ThrOK.and.IPart.ne.NPrtUS) IPart = NProc(-NPrtUS)
              NPrtSt = 0
            else
              IPart = NProc(-1)
              NPrtSt = NPrtUS - 1
              endIf
            If(DoPrnt) Write(IOut,1080) NPrtUS, ThrOK, IAlg, NPAlg,
     $        DoDPD, LenP, LenD1P
C$OMP Parallel If(ThrOK) Default(Shared)
C$OMP+ Private(ID1E,ID2E,IEC,IElSum,IEX,IndV,IPart,LenVAB,LenVX,IBNorm)
C$OMP+ Private(IVA,IVB,IVXA,IVXB,IVxAI,IVxBI,IZTT,LenD1E,LenD2E,LenEXC)
C$OMP+ Private(LenVxI,RJunk,FName,IOutD,IPrinD,IPartX,IDA,IDB,IDPA,IDPB)
            Do 200 IPartX = NPrtSt, 0, -1
              If(ThrOK) then
                IPart = MDPart(0)
              else
                IPart = IPartX
                endIf
              If(DoDump.ne.0) then
                FName = ' '
                Write(FName,1090) NCall, IPart
                IOutD = 50 + IPart
                IPrinD = Max(IPrint,10)
                Open(Unit=IOutD,File=FName,Status='Unknown')
              else
                IOutD = IOut
                IPrinD = IPrint
                endIf
              Call TStamp(6,' ')
              If(IPart.eq.0) then
                Call CalDFT(IOutD,IPrinD,DoWrt,NFrqR1,NMat0,NDPol,
     $            ICnBeg,ICnEnd,DoZCmp,Do1Mat,DoSpar,ElSum,Ex,Ec,IIV,VA,
     $            VB,D1E,D2E,VxA,VxB,VxAI,VxBI,V(IBNSum),ICntrX,IExchn,
     $            ICorr,IExCoW,IRadAn,IRanWt,IRanGd,ScaDFX,V(ICMap),
     $            V(IAnMap),V(IAtTMp),V(IAtCMp),NAtoms,NAtomA,NAtAll,
     $            NClRep,NMtPBC,IP,PA,D1PA,GA,PB,D1PB,GB,NMOA,CMOA,NMOB,
     $            CMOB,NBas6D,NBsAll,IOpCl,Thresh,NOpAtG,RotAll,NEqAll,
     $            V(iDoAtm),V(IFacAt),V(iIPrP1),FreqP1,ECTol,MaxPtI,
     $            MinPt,NDBF,V(jMpJAB),V(iBegSB),V(jShlSB),V(IMpAtm),
     $            V(IDoInA),IPFlag,MaxMBI,SizInc,IChHar,MulHar,DoJig,
     $            IRType,IntCel,BckCel,BckDim,CntCel,ISavGd,IGWInf,
     $            IAtBtD,IRdBtD,RRdBtD,RGWBtD,V(iMpJAn),V(iBegSh),
     $            V(jAtmSz),V(jShlSz),V(jAtmSD),NDisMt,V(jDisMt),
     $            V(ILT2D),V(ILT2D1),0,NPrtUL,NPartT,NxtVal,V(iCnvrt),
     $            V(iScale),MaxNG,V(INAtG),V(jAtEC),Omega,LinDyn,
     $            V(IIJV),V(IIJKV),V(IIJU),V(IIJKU),ZTT,V(IV0),
     $            V(IV0),LenVP)
                If(ElSum.eq.(-One))
     $            Call GauErr('Process 0 failed to complete !')
              else
                Call AlCaSu(LenVP+IV0-1+(IPart-1)*LenP1,IOpCl,NAtoms,
     $            NAtomA,NMat0,NMat,NMtCP1,NMtFkD,NDPol,NBas6D,LenP,
     $            LenD1P,DoE,DoVAB,DoD1E,DoD2E,DoCPKS,DoHyp,DoChi,
     $            DoGmPx,DoGmPe,PeGmx,IExCoW,DoFx,DoDPD,DoHir,LenExc,
     $            LenVAB,LenD1E,LenD2E,LenVx,LenVxI,IElSum,IEx,IEc,IDA,
     $            IDB,IVA,IVB,ID1E,ID2E,IDPA,IDPB,IVxA,IVxB,IVxAI,IVxBI,
     $            IBNorm,IZTT,NZV,IndV)
                If(IndV.ne.(IPart*LenP1+IV0-1))
     $            Call GauErr('Error #1 in CalDSu.')
                V(IElSum) = -One
                Call AClear(LenExc,V(IEx))
                Call AClear(LenExc,V(IEc))
                Call AClear(LenVAB,V(IVA))
                Call AClear(LenVAB,V(IVB))
                Call AClear(LenD1E,V(ID1E))
                Call AClear(LenD2E,V(ID2E))
                Call AClear(LenVx,V(IVxA))
                Call AClear(LenVx,V(IVxB))
                Call AClear(LenVxI,V(IVxAI))
                Call AClear(LenVxI,V(IVxBI))
                Call AClear(3,V(IZTT))
                If(DoDPD) then
                  Call AMove(LenP,PA,V(IDA))
                  Call AMove(LenD1P,D1PA,V(IDPA))
                  If(NSpBlX.eq.2) then
                    Call AMove(LenP,PB,V(IDB))
                    Call AMove(LenD1P,D1PB,V(IDPB))
                    endIf
                  Call CalDFT(IOutD,IPrinD,DoWrt,NFrqR1,NMat0,NDPol,
     $              ICnBeg,ICnEnd,DoZCmp,Do1Mat,DoSpar,V(IElSum),V(IEx),
     $              V(IEc),IIV,V(IVA),V(IVB),V(ID1E),V(ID2E),V(IVxA),
     $              V(IVxB),V(IVxAI),V(IVxBI),V(IBNorm),ICntrX,IExchn,
     $              ICorr,IExCoW,IRadAn,IRanWt,IRanGd,ScaDFX,V(ICMap),
     $              V(IAnMap),V(IAtTMp),V(IAtCMp),NAtoms,NAtomA,NAtAll,
     $              NClRep,NMtPBC,IP,V(IDA),V(IDPA),GA,V(IDB),V(IDPB),
     $              GB,NMOA,CMOA,NMOB,CMOB,NBas6D,NBsAll,IOpCl,Thresh,
     $              NOpAtG,RotAll,NEqAll,V(iDoAtm),V(IFacAt),V(iIPrP1),
     $              FreqP1,RJunk,MaxPtI,MinPt,NDBF,V(jMpJAB),V(iBegSB),
     $              V(jShlSB),V(IMpAtm),V(IDoInA),IPFlag,MaxMBI,SizInc,
     $              IChHar,MulHar,DoJig,IRType,IntCel,BckCel,BckDim,
     $              CntCel,ISavGd,IGWInf,IAtBtD,IRdBtD,RRdBtD,RGWBtD,
     $              V(iMpJAn),V(iBegSh),V(jAtmSz),V(jShlSz),V(jAtmSD),
     $              NDisMt,V(jDisMt),V(ILT2D),V(ILT2D1),IPart,NPrtUL,
     $              NPartT,NxtVal,V(iCnvrt),V(iScale),MaxNG,V(INAtG),
     $              V(jAtEC),Omega,LinDyn,V(IIJV),V(IIJKV),V(IIJU),
     $              V(IIJKU),V(IZTT),V(IndV+1),V(IndV+1),LenVP)
                else
                  Call CalDFT(IOutD,IPrinD,DoWrt,NFrqR1,NMat0,NDPol,
     $              ICnBeg,ICnEnd,DoZCmp,Do1Mat,DoSpar,V(IElSum),V(IEx),
     $              V(IEc),IIV,V(IVA),V(IVB),V(ID1E),V(ID2E),V(IVxA),
     $              V(IVxB),V(IVxAI),V(IVxBI),V(IBNorm),ICntrX,IExchn,
     $              ICorr,IExCoW,IRadAn,IRanWt,IRanGd,ScaDFX,V(ICMap),
     $              V(IAnMap),V(IAtTMp),V(IAtCMp),NAtoms,NAtomA,NAtAll,
     $              NClRep,NMtPBC,IP,PA,D1PA,GA,PB,D1PB,GB,NMOA,CMOA,
     $              NMOB,CMOB,NBas6D,NBsAll,IOpCl,Thresh,NOpAtG,RotAll,
     $              NEqAll,V(iDoAtm),V(IFacAt),V(iIPrP1),FreqP1,RJunk,
     $              MaxPtI,MinPt,NDBF,V(jMpJAB),V(iBegSB),V(jShlSB),
     $              V(IMpAtm),V(IDoInA),IPFlag,MaxMBI,SizInc,IChHar,
     $              MulHar,DoJig,IRType,IntCel,BckCel,BckDim,CntCel,
     $              ISavGd,IGWInf,IAtBtD,IRdBtD,RRdBtD,RGWBtD,V(iMpJAn),
     $              V(iBegSh),V(jAtmSz),V(jShlSz),V(jAtmSD),NDisMt,
     $              V(jDisMt),V(ILT2D),V(ILT2D1),IPart,NPrtUL,NPartT,
     $              NxtVal,V(iCnvrt),V(iScale),MaxNG,V(INAtG),V(jAtEC),
     $              Omega,LinDyn,V(IIJV),V(IIJKV),V(IIJU),V(IIJKU),
     $              V(IZTT),V(IndV+1),V(IndV+1),LenVP)
                  endIf
                endIf
              If(DoDump.ne.0) Close(Unit=IOutD)
  200         Continue
C$OMP END PARALLEL
C
C           Sum up contributions.
C
            Call CASumR(IOut,TrcLin,ThrOK,0,NPrtUS,NPrtUL,LenVP,LenP1,
     $        IOpCl,NAtoms,NAtomA,NMat0,NMat,NMtCP1,NMtFkD,NDPol,NBas6D,
     $        LenP,LenD1P,DoE,DoVAB,DoD1E,DoD2E,DoCPKS,DoHyp,DoChi,
     $        DoGmPx,DoGmPe,PeGmx,IExCoW,DoFx,DoDPD,DoHir,LenExc,LenVAB,
     $        LenD1E,LenD2E,LenVx,LenVxI,NZV,ElSum,Ex,Ec,V(IBNSum),VA,
     $        VB,D1E,D2E,VxA,VxB,VxAI,VxBI,ZTT,V(IV0))
            Junk = GLinCo(4,LinDyn,CountN,NxtVal,0,0)
            Goto 500
            endIf
  390     Continue
  400   Continue
      Call GauErr('Logic error in CalDSu.')
C
C     Restore the densities and the XC potential matrices to (A,B,X,Y).
C
  500 If(NDimBl.eq.2) then
C       print zero torque theorem sum.
        Write(IOut,1120) ZTT(1),ZTT(2),ZTT(3)
C       density matrix.
        Call SumDif(.False.,NTT6D,PA,PB)
        Call AScale(NTT6D,Pt5,PA(1+NTT6D),PA(1+NTT6D))
        Call AScale(NTT6D,Pt5,PB(1+NTT6D),PB(1+NTT6D))
C       density matrix derivatives.
        If(DoCPKS.or.DoHyp.or.DoPolD.or.DoGxPT)
     $    Call GauErr('GKS-CPKS NYI in CalDSu.')
        If(TDKSLg.or.PKSDen)
     $    Call GauErr('GKS-Post-KS NYI in CalDSu.')
C       XC matrix.
        If(DoVAB) then
          Call ASub(NTT6D,VA,VB,VB)
          Call ACSASB(NTT6D,VA,Two,VB,-One,VA)
          endIf
C       XC matrix derivatives.
        If(DoNuFx) Call GauErr('NuFx NYI in CalDSu.')
        endIf
C
C     Apply symmetry
C
      If(DoVAB.and.NMtPBC.gt.1) then
        Call SymPBK(2,NBas6D,NTT6D,NMat0,NMat0,NMtPBC,VA)
        If(NSpBlX.eq.2)
     $    Call SymPBK(2,NBas6D,NTT6D,NMat0,NMat0,NMtPBC,VB)
        endIf
      If(NOpUse.gt.1) then
        If(DoE) then
          Ex = Ex * NOpUse
          Ec = Ec * NOpUse
          endIf
        ElSum = ElSum * NOpUse
        LNES = InToWP(NShell*NOpUse)
        Call SetInR(1,0,0,MinL,MaxL,LenRot,IndR)
        IndRI = IV0
        INEqS = IndRI + LenRot
        IScr1 = INEqS + LNES
        IScr2 = IScr1 + NBas6D
        IV0 = IScr2 + InToWP(NShell)
        Call TstCor(IV0,MDV,'CalDSu-BNSymm')
        MDV1 = MDV - IV0 + 1
        Call FileIO(2,-IONEqS,LNES,V(INEqS),0)
        Call BNSymm(NShell,NOpUse,V(INEqS),V(IBNSum),V(IScr1))
        If(DoSymF) then
          If(DoD1E) Call FrcSym(.True.,0,NAtoms,D1E,NOpUse,NEqAll,
     $      RotAll,V(IV0))
          If(DoD2E) Call FFSym(.True.,0,NAtoms,D2E,NOpUse,NEqAll,RotAll,
     $      V(IV0))
          If(DoVAB.or.DoFx) then
            If(DoSpar) then
              LNEB = InToWP(NBas6D*NOp1)
              INEqB = IV0
              IVSym = INEqB + LNEB
              IV0 = IVSym + NZV+NBas6D
              Call TstCor(IV0,MDV,'CalDSu-EnfSym')
              Call FileIO(2,-IONEqB,LNEB,V(INEqB),0)
            else
              IIndSh = IV0
              IV0 = IIndSh + InToWP(NShell)
              Call TstCor(IV0+NTT6D,MDV,'CalDSu-FNAbSy')
              endIf
            MDV1 = MDV - IV0 + 1
            If(DoVAB) then
              If(DoSpar) then
                Call EnfSym(IIV,VA,V(IVSym),V(INEqB),NOp1,NBas6D,
     $            'XC Alpha Matrix',.false.)
                If(NSpBlX.eq.2) Call EnfSym(IIV,VB,V(IVSym),V(INEqB),
     $            NOp1,NBas6D,'XC Beta Matrix',.False.)
              else
                If(DoPolD) then
                  NumMat = NDPol/NFrqR1
                  NumFr = NFrqR1
                else
                  NumMat = NMat
                  NumFr = 1
                  endIf
                MaxMat = (MDV-IV0+1)/(NumMat*NTT6D)
                Call AdDF1N(.False.,.False.,NTT6D,NumMat,VA)
                If(DoPolD) then
                  Call FNAbSX(.True.,.False.,.False.,0,0,NOpUse,NumMat,
     $              NTT6D,NumFr,NShell,NEqAll,V(INEqS),RotAll,ShellT,
     $              ShellC,VA,V(IV0),MaxMat,V(IIndSh),MinL,MaxL,IndR,
     $              LenRot,V(IndRI))
                else
                  Call FNAbSy(.True.,NOpUse,NTT6D,NumMat,NumMat,NShell,
     $              V(INEqS),RotAll,ShellT,ShellC,VA,V(IV0),V(IIndSh),
     $              MinL,MaxL,IndR,LenRot,V(IndRI))
                  endIf
                Call AdDF1N(.False.,.True.,NTT6D,NumMat,VA)
                If(NSpBlX.eq.2) then
                  Call AdDF1N(.False.,.False.,NTT6D,NumMat,VB)
                  If(DoPolD) then
                    Call FNAbSX(.True.,.False.,.False.,0,0,NOpUse,
     $                NumMat,NTT6D,NumFr,NShell,NEqAll,V(INEqS),RotAll,
     $                ShellT,ShellC,VB,V(IV0),MaxMat,V(IIndSh),MinL,
     $                MaxL,IndR,LenRot,V(IndRI))
                  else
                    Call FNAbSy(.True.,NOpUse,NTT6D,NumMat,NumMat,
     $                NShell,V(INEqS),RotAll,ShellT,ShellC,VB,V(IV0),
     $                V(IIndSh),MinL,MaxL,IndR,LenRot,V(IndRI))
                    endIf
                  Call AdDF1N(.False.,.True.,NTT6D,NumMat,VB)
                  endIf
                endIf
              endIf
            If(DoNuFx) then
              MaxMat = (MDV-IV0+1)/(3*NTT6D)
              Call AdDF1N(.False.,.False.,NTT6D,NMtFkD,VxA)
              If(DoGmx) Call GauErr('FNAbSX for DoGmx NYI')
              Call FNAbSX(.True.,.False.,.False.,NAtoms,ICnEnd,NOpUse,3,
     $          NTT6D,1,NShell,NEqAll,V(INEqS),RotAll,ShellT,ShellC,VxA,
     $          V(IV0),MaxMat,V(IIndSh),MinL,MaxL,IndR,LenRot,V(IndRI))
              Call AdDF1N(.False.,.True.,NTT6D,NMtFkD,VxA)
              If(NSpBlX.eq.2) then
                Call AdDF1N(.False.,.False.,NTT6D,NMtFkD,VxB)
                Call FNAbSX(.True.,.False.,.False.,NAtoms,ICnEnd,NOpUse,
     $            3,NTT6D,1,NShell,NEqAll,V(INEqS),RotAll,ShellT,ShellC,
     $            VxB,V(IV0),MaxMat,V(IIndSh),MinL,MaxL,IndR,LenRot,
     $            V(IndRI))
                Call AdDF1N(.False.,.True.,NTT6D,NMtFkD,VxB)
                endIf
              endIf
            endIf
          endIf
        endIf
C
C     Check the electron count...
      If(NE.gt.0) then
        ErrEC = Abs(GFloat(NE+NElCor)-ElSum)
        ECRel = ErrEC / GFloat(NElTot)
        ErrBN = Zero
        ErrBN2 = Zero
        ECTol1 = GFloat(12)*ECTol
        If(NOpUse.ge.NOp2) ECTol1 = ECTol1*GFloat(NOpUse)/GFloat(NOp2)
        IShErr = 0
        IShEr2 = 0
        If(NAtCel.ne.0.and.NClRep.eq.1) then
          NCl1 = NAtoms / NAtCel
          NShCk = NShell / NCl1
          Call SumCol(.False.,NShCk,NShCk,NCl1-1,V(IBNSum+NShCk),
     $      V(IBNSum))
        else
          NShCk = NShell
          endIf
        Do 510 I = 1, NShCk
          NFI = NShFOS(ShellC(I),ShellT(I))
          Err = V(IBNSum+I-1)/GFloat(NFI) - One
          If(Err.gt.ErrBN2) then
            ErrBN2 = Err
            IShEr2 = I
            endIf
          Err = Abs(Err)
          If(Err.gt.ErrBN) then
            ErrBN = Err
            IShErr = I
            endIf
  510     Continue
        Fail = (ECRel.gt.ECTol.and.ChekEC).or.((ErrBN2.gt.ECTol2.or.
     $    ErrBN.gt.ECTol1).and..not.HavOCBf)
        If(Fail.or.IPrint.ge.2) then
          If(Fail) then
            Write(IOut,1010)
          else
            Write(IOut,1015)
            endIf
          If(HavOCBf) then
            Write(IOut,1055) NE, NElCor, ErrEC, ECRel, ECTol, IShErr,
     $        ErrBN, IShEr2, ErrBN2
          else
            Write(IOut,1050) NE, NElCor, ErrEC, ECRel, ECTol, IShErr,
     $        ErrBN, ECTol1, IShEr2, ErrBN2, ECTol2
            endIf
          endIf
*
* ANT BEGIN modification
*
*         Do not abort after inaccurate quadrature of density
*
*         If(Fail.and.Fatal)
*      $    Call GauErr('Inaccurate quadrature in CalDSu.')
*
* ANT END modification
*
        endIf
C
C     Attempt to trap XC matrix elements which are supposed to
C     be identically zero by symmetry, but which came out slightly
C     nonzero due to grid integration error...
C
      DoWrt = .False.
      Tol = MDCutO(0)
      IPart = IPopVc(LenVAB,VA,Tol,VA)
      If(.not.DoHir) IPart = IPopVc((NSpBlX-1)*LenVAB,VB,Tol,VB)
      IPart = IPopVc(LenD1E,D1E,Tol,D1E)
      IPart = IPopVc(LenD2E,D2E,Tol,D2E)
      IPart = IPopVc(LenVx,VxA,Tol,VxA)
      IPart = IPopVc((NSpBlX-1)*LenVx,VxB,Tol,VxB)
      IPart = IPopVc(LenVxI,VxAI,Tol,VxAI)
      IPart = IPopVc((NSpBlX-1)*LenVxI,VxBI,Tol,VxBI)
C
  999 IPart = NProc(-NPSSav)
      IPart = NProcL(-NPLSav)
      If(NOpUse.gt.1)
     $  Call UnPckI(0,Junk,NEqAll,MxAtSO,NOpAll,NAtoms,NOpAll)
      Call TStamp(1,'Bot of CalDSu')
      Return
      End
