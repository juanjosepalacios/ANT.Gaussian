*Deck Optmz
      Subroutine Optmz(V,MDV)
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
C     Driver for gradient optimization program.  This program uses
C     the method of Schlegel (J. Comp. Chem. 3, 214 (1982)) to which a
C     dynamic trust radius (a la Fletcher) and appropriate modification of
C     large steps (as done by Simons and others) have been added.
C     Care is also taken do the linear step between the most recent
C     point and the previous point of lowest energy, rather than
C     between the two most recent points as in the original algorithm.
C
      Dimension V(*)
      Integer MaxIOp
      Parameter (MaxIOp=DEFMAXIOP)
      COMMON /IOP/    IOp(MaxIOp)
      COMMON /MUNIT/  IUNIT(20)
      COMMON /IO/     IN,IOut,IPUNCH
      COMMON /PHYCON/ PhyCon(30)
      Equivalence (PhyCon(1),ToAng)
C
C1GRDNT
C
C     L103 COMMUNICATES WITH Gaussian PRIMARILY THROUGH /GRDNT/ WHICH
C    CONTAINS THE COMPUTED VALUE OF THE FUNCTION, FIRST AND POSSIBLY
C    SECOND DERIVATIVES.
C    COMMON/GRDNT/
C
C     COMMON /GRDNT/ ENERGY,NVar,IGetFC,F(NVar),FRCNST(NVarTT)
C
C     ENERGY -- THE VALUE OF THE FUNCTION
C     NVar   -- THE NUMBER OF VARIABLES
C     IGETFC -- SAVED COPY OF IOp(10) FOR RESTARTING.
C               (0 IF ONLY FIRST DERIVATIVES ARE AVAILABLE,
C               4 or 7 IF FRCNST CONTAINS SECOND DERIVATIVES).
C     F      -- FORCES (NEGATIVE OF FIRST DERIVATIVE -DE/DX)
C     FRCNST -- SECOND DERIVATIVES IF EVALUATED EXTERNALLY
C
      COMMON/GRDNT/ENERGY,NVar,IGETFC
C
C1OPTGRD
C
C    LOCAL OPTIMIZATION DATA FOR L103 IS KEPT IN COMMON/OPTGRD/
C
      Logical Prnt,Exit,HavDEP,DoReac,DoProd,DidQST,DoMid,DoBiRx,
     $  DoIRC,Do601,Do601H,DoClnS,IniFTr,UseOT,OptMic,UseCA,RotQM,
     $  SkpMic,SmlStp,IntHes,DoEle,MMRest
      Common /OptGrd/ ConvF,FMaxT,DXMaxT,RMax,RMin,RLim,EigMax,Eigmin,
     $  DEPred,FSwtch,FncErr,GrdErr,FncCnv,RLam,NRisDI,NStep,IStep,NMax,
     $  NP,Neg,Prnt,Exit,IPSav1,IPSav2,HavDEP,NGDIIS,IRisen,MaxRis,
     $  Large,ILevel,ISStep,IOffs,IMDS,MaxSv,MaxSv2,NVrMMT
C     IC(NVar),X(NVar),XName(NVar),ES,FS(NVar),FC(2*NVarTT),
C     XX(NVar,NVar),FF(NVar,NVar)
C     Equivalence (ES,ESave(1))
C
C     XX     -- Previous coordinates.
C     ConvF  -- Convergence criterion.
C     FMaxT  -- Maximum force allowed.
C     DXMaxT -- Maximum displacement allowed.
C     RMax   -- See D2Corr and Schmdt.
C     RMin   -- See D2Corr and Schmdt.
C     RLim   -- See D2Corr and Schmdt.
C     EigMax -- Maximum allowed eigenvalue of FC.
C     EigMin -- Minimum allowed eigenvalue of FC.
C     DEPred -- Predicted change in energy if HavDEP is set.
C     FSwtch -- Minimum RMS force for linear search.
C     FncErr -- Estimated error in Energy.
C     GrdErr -- Estimated error in F.
C     FncCnv -- Function value unit conversion factor.
C     RLam   -- Saved value of RFO Lambda.
C     NRisDI -- Number of time energy has risen in GEDIIS.
C     NStep  -- Maximum number of steps allowed.
C     IStep  -- Current step number.
C     NMax   -- Maximum number of points to store.
C     NP     -- Number of points stored.
C     Neg    -- Order of stationary point to search for.
C     Prnt   -- Print switch.
C     Exit   -- Usually indicates convergence, see also Star.
C     IPSav1 -- Indicates the point for which analytic second
C               derivatives are in the first part of FC. Zero
C               if there are none.
C     IPSav2 -- Same as IPSav1 but for the second half.
C     HavDEP -- Whether a prediction for Delta-E is available.
C     NGDIIS -- Number of vectors to use in GDIIS/GDPIS.
C     IRisen -- Number of succesive steps at which the energy has risen.
C     MaxRis -- Maximum number of steps at which the energy can rise
C               before only line minimizations are done.
C     Large  -- Indicates the use of fast transformations for large
C               molecules.
C     ILevel -- Actual level in guided optimization; no longer used.
C     ISStep -- Sum of steps used for previous lower level
C               optimizations; no longer used.
C     IOffs  -- Indicates the offset data to use. Zero is none.
C     IMDS   -- Control flag for multi-dimensional search
C               1 = Old GDIIS.
C               2 = First-Order Simultaneous Optimization
C               3 = New GEDIIS/GDIIS.
C               4 = Fixed coordinates?
C               5 = Mixed Opt in new FOSimOpt routines
C     IC     -- Second derivative info on initial entry scratch storage
C               subsequently.
C     X      -- Coordinates.
C     XName  -- Print labels for coordinates.  Delimited string.
C     ES     -- The value of the function at the previous point.
C     FS     -- The value of the function at earlier points.
C     FC     -- Second derivative matrices.
C     FF     -- Previous forces.
C
#include "commonmol.inc"
#include "repall.inc"
#include "commonz.inc"
      Common /SymInf/ NOp1,NOp2,JTrans(3,8),TRot(3,3),TrVec(3)
      Logical UpTrus, Sphere, DoLin, TstEig, SetStp, FCart, DoRed,
     $  LSTQST, StOnly, First, DoChk, DoRest, DoEstm, Large, DidMix,
     $  Failed, AlUnit, UpdMol, DoMicO, LJJ(1), ONLnOK, MMStpd, Convgd,
     $  PCM, JustCk, HaveRB, DoClnG, FMM, AllowP(50), UpdRWF, NewCrd,
     $  FixedC, ConvP, SavLMO, RmTrRo, RawDIS
      Integer FMFlag, FMFlg1
      Parameter (IOGen=501,IOGenA=505,IGrdNT=511,IOptGr=575,IOScan=550,
     $  IRwF1=584,IRwF2=585,IOZRed=665,IOSLEq=671,IRwSM=672,IRwMM=687,
     $  IRwSvC=698,IMMRFO=713,IMMMOD=714,IMMCRS=715,IRwCrd=775,
     $  MStepD=100,LRwGen=1000,IOPath=601,IRFOTS=748,IGIOld=783,
     $  NStatD=2,MaxMDS=4,MOmega=DEFMOMEGA,NOmega=DEFNOMEGA)
      Character*12 NamMDS(MaxMDS)
      Dimension CSave(3,MaxAtm), NTRedS(4), IntFl2(2), JJ(1), XX(1),
     $  CelVec(3,3), Stats(NStatD), EField(0:3), JScan(4), TRotT(3,3),
     $  Omega(MOmega,NOmega), QMVar(30)
      Save Zero, One, Ten, Big, DConEq, JJ, XX, MinMMC, OneM4, NamMDS,
     $  LJJ, Pt1, Pt9, BigE
      Data Zero/0.0d0/, One/1.0d0/, Ten/1.D1/, Big/1000.0d0/,
     $  DConEq/1.d-8/, JJ/0/, XX/0.0d0/, MinMMC/5000/, OneM4/1.0d-4/,
     $  NamMDS/'GDIIS','FOSimult','GEDIIS/GDIIS','FO SCF Only'/,
     $  LJJ/.False./, BigE/1.d20/, Pt1/0.1d0/, Pt9/0.9d0/
 1000 Format(1X,18(4HGrad))
 1005 Format(' Berny optimization.')
 1010 Format(' Initialization pass.')
 1020 Format(' Numerically estimating second derivatives.')
 1030 Format(' Numerical second derivative computation complete.')
 1040 Format(' Using ',A,' optimizer.')
 1050 Format(' Number of steps in this run=',I7,
     $       ' maximum allowed number of steps=',I7,'.')
 1060 Format(' No Z-matrix variables, so optimization will use',
     $  ' Cartesian coordinates.')
 1070 Format(1X)
 1080 Format(' Search for a saddle point of order',I3,'.')
 2999 Format(1X,19(4H++++))
 3000 Format(1x,'+       ',I3,'-point reaction path relaxation',
     $33x,'+')
 3010 Format(1x,'+       ',4x,'reactant structure will be optimized',
     $27x,'+')
 3020 Format(1x,'+       ',4x,'product structure will be optimized',
     $28x,'+')
 3030 Format(1x,'+       ',4x,'point with highest energy will ',
     $'optimize to TS',18x,'+')
 3040 Format(1x,'+       ',4x,'reaction assumed bimolecular',
     $35x,'+')
 3042 Format(1x,'+       ',4x,'reactants assumed bimolecular',
     $34x,'+')
 3044 Format(1x,'+       ',4x,'products assumed bimolecular',
     $35x,'+')
 3050 Format(1x,'+       ',4x,'will use the eigenvector of the first ',
     $'point',20x,'+')
 3060 Format(' Mismatch in second derivative file, Len=',I8,
     $  ' but NVTT=',I8,'.')
 4050 Format(' Incomplete coordinate system.  Try restarting with',/,
     $  ' Geom=Check Guess=Read Opt=(ReadFC,NewRedundant)')
 4070 Format(' Optmz: redundant internal coordinate code flag=',I2)
 4080 Format(' Save coordinates for EBest=',G23.12)
 4090 Format(' Max number of microiteration restarts exceeded.')
C
C     INITIALIZE THIS LINK.  THE ILSW BIT WILL BE SET UNLESS THIS
C     IS THE FIRST ENTRY TO L103, IN WHICH CASE WE'LL JUST DO SOME
C     INITIALIZATION.
C
      Call Drum(V,MDV)
      IPrint = IOp(33)
      Call SetPFl(IOut,IPrint,0,IPFlag,AllowP,FMM,FMFlag,FMFlg1,NFxFlg,
     $  IHMeth,Omega,0,Junk1,Junk2,V,MDV)
      IPrt1 = IPrint
      If(IPSave(0).le.1) IPrt1 = Max(IPrint,1)
      IGISav = 0
      NStat = NStatD
      Call AClear(NStat,Stats)
      Call IClear(4,JScan)
      MaxJb1 = 1
      IJOff1 = 0
      IType1 = 0
      IStep1 = 0
      IJob1 = 1
      MMStpd = .False.
      Do601 = .False.
      Do601H = .False.
      UpdMol = .False.
      LenCB = 0
      LenSB = 0
      IOSLE1 = 0
      IRwSMM = 0
      MaxGDI = 0
      IExit = 1
      ConvEq = DConEq
      LGrdnt = 1 + InToWP(2)
      LOptGr = 14 + InToWP(22)
      Write(IOut,1070)
      Write(IOut,1000)
      If(IPSave(0).le.1) Write(IOut,1005)
      Exit = .False.
      NNeg = IOp(5)
      UpTrus = Mod(IOp(9),10).eq.2.or.
     $  (Mod(IOp(9),10).eq.0.and.NNeg.eq.0)
      Sphere = (Mod(IOp(9),100)/10).eq.2.or.
     $  ((Mod(IOp(9),100)/10).eq.0.and.NNeg.eq.0)
      DoRest = IOp(35).ne.0
      If(DoRest) then
        Call ChkPnt(1,V,MDV)
        Call FileIO(2,-IGrdnt,LGrdnt,Energy,0)
      else
        IGetFC = IOp(10)
        endIf
      IQST = IOp(27)/100
      LSTQST = Mod(IQST,10).ne.0
      TstEig = IOp(11).eq.2.or.
     $  (IOp(11).eq.0.and.NNeg.ne.0.and..not.LSTQST)
      MthUpd = IOp(13) - 3
      StOnly = Mod(IOp(18),10).ne.0
      FCart = ((Mod(IOp(18),100)/10).eq.1)
      MicOTy = Mod(IOp(18),10000000) / 1000000
      SkpMic = .False.
      JustCk = MicOTy.eq.4.or.IOp(30).eq.1.or.IOp(30).eq.2
      UpdRWF = IOp(30).ne.1.and.IOp(30).ne.2
      SetStp = IOp(30).eq.0.or.IOp(30).eq.2
      If(UpdRWF) then
        MStepU = MStepD
        If(IOp(6).lt.0) MStepU = Max(MStepU,-IOp(6))
      else
        MStepU = 2
        endIf
      MethQ = IOp(19)/2
      FixedC = MethQ.eq.8
      ConvP = MethQ.eq.9
      If(FixedC.or.ConvP) MethQ = 6
      IQST = IOp(27)/100
      LSTQST = Mod(IQST,10).ne.0
      NTrRot = IOp(28)
      IFirst = Mod(IOp(38),10)
      NPath= IOp(42)
      DoLin = Mod(IOp(19),2).eq.0.and.(NNeg.eq.0.or.NPath.gt.0)
      First = IFirst.ne.0
      If(First) NGDIIS = 0
      If(MethQ.gt.4) then
        IMDS = MethQ - 4
        If(IMDS.eq.1.and.IOp(55).ne.0) IMDS = -IAbs(IOp(55))
        If(IMDS.gt.MaxMDS) Call GauErr('Illegal IOp(19).')
      else
        IMDS = 0
        endIf
      Large = ((Mod(IOp(18),100)/10).eq.5).and.NPath.eq.0
      DidQST = Mod(IOp(27),10).gt.2
      I49 = IOp(49)
      DoReac = Mod(I49,10).le.1
      DoProd = (Mod(I49,100)/10).le.1
      IMid = Mod(I49,1000)/100
      If(IMid.eq.0) IMid = 3
      Ibirx = Mod(I49,10000)/1000
      DoMid = IMid.eq.1.and.DidQST
      DoBiRx = IBirx.eq.1
      DoIRC = IMid.eq.2.and.DoReac
      EField(0) = Zero
      Call FileIO(2,-IOGen,3,EField(1),1)
      If(SProd(4,EField,EField).gt.Zero) then
        MaxLEF = 1
      else
        MaxLEF = -1
        endIf
C
C     DoReac : Whether to optimize reactant structure
C     DoProd
C     DidQST : a 3 structure input was used
C     DoMid : a 3 structure input was used and the guess TS
C     structure will optimize to TS irrespective of its relative
C     energy
C     DoBiRx : for a bimolecular reaction. The end point will be
C     treated differently (see paper)
C     DoIRC : a 2 or 3 structure input with the "reactant structure
C     being the TS. A frequency calc. is needed in the first pass
C     to get the transition vector.
C
      If(First) then
        If(ITqry(IOSLEq).gt.0.and..not.Large)
     $    Call FileIO(5,IOSLEq,0,V,0)
        ILevel = 0
        IOffs = 0
        ISStep = 0
        NGDIIS = 0
        Failed = .False.
        endIf
      If(NPath.gt.0) then
        Do601 = .false.
        If(First) then
          NCurr = 1
          Write(IOut,2999)
          Write(IOut,3000) NPath
          If(DoReac) Write(IOut,3010)
          If(DoProd) Write(IOut,3020)
          If(.not.DoMid.and..not.DoIRC) Write(IOut,3030)
          If(DoIRC) Write(IOut,3050)
          If(DoBirx.and..not.DoProd.and..not.DoReac) Write(IOut,3040)
          If(DoBirx.and.DoProd) Write(IOut,3042)
          If(DoBirx.and.DoReac) Write(IOut,3044)
          Write(IOut,2999)
          endIf
        LSTQST = .False.
        IQST = 2
C       stores the original step size, just in case ...
        If(IOP(8).gt.0) then
          DXMXT0 = GFloat(IOP(8))/GFloat(100)
        else
          DXMXT0 = GFloat(3)/GFloat(10)
          endIf
        endIf
      If(MethQ.eq.0) then
        If(SetStp) then
          MethQ = 3
        else
          MethQ = 4
          endIf
        endIf
      Call GetPCM(IOut,IPrint,.False.,0,1,ISCRF,RJunk)
      If(PCM(ISCRF)) Call PCMIOp('load',' ',IJunk,RJunk)
      DoChk = IOp(36).eq.0
      Wght = GFloat(IOp(39)) / GFloat(100)
      IOldRC = Mod(IOp(59),10)
      IOldRB = 3
      If(IOldRC.eq.0) then
        IOldRC = 1
      else if(IOldRC.eq.4) then
        IOldRC = 1
        IOldRB = 2
        endIf
      IAprBG = Mod(IOp(59),100000)/10
      ITRScl = IOp(110)
      IPrjct = 0
      IGuess = IOp(96)
      IRFO = IOp(97)
      IQMac = IOp(98)
      IDir = Mod(IAbs(IQMac),10)
      IntHes = IDir.ge.6
      If(IDir.eq.7) then
        IDir = 5
      else if(IDir.eq.6) then
        IDir = 0
        endIf
      IDir0 = IDir
      If(IntHes.and.IGetFC.ne.4) then
C       No analytical H, so switch to standard micro-iterations,
C       and delete the files associated with QuadMac.
        IntHes = .False.
        If(IDir.lt.2) then
          If(ITqry(IMMRFO).gt.0) Call FileIO(5,IMMRFO,0,XX,0)
          If(ITqry(IRFOTS).gt.0) Call FileIO(5,IRFOTS,0,XX,0)
          endIf
      else if(IntHes) then
        IDir = 3
        endIf
      If(IDir.eq.0.and.MicOTy.eq.1) IDir = 1
      If(NNeg.gt.1.and.IDir.gt.1)
     $  Call GauErr('No second or higher order saddle point.')
      IDirMM = Mod(IAbs(IQMac),100)/10
      If(IDirMM.eq.0.and.MicOTy.eq.1) IDirMM = 1
      IHVAcc = IOp(99)
      If(IOp(99).lt.0) IHVAcc = 0
      If(IOp(99).eq.0) IHVAcc = 3001208
      UseCA  = IDir.ge.2.and.IQMac.gt.0
      ISusp = Mod(IOp(93),1000)/100
      If(ISusp.eq.0) ISusp = 1
      IGetTS = Mod(IOp(93),100)/10
      If(IGetTS.eq.0) IGetTS = 1
      IGetTS = IGetTS + 10*ISusp
      ITSChk = Mod(IOp(93),10)
      If(ITSChk.eq.0) ITSChk = 2
      If(ITSChk.eq.2) IGetTS = -IGetTS
      IGssMM = IOp(94)
      IRFOMM = IOp(95)
      IStepM = IOp(111)
      If(IStepM.eq.0) IStepM = 50
      SteepM = GFloat(IStepM) / GFloat(1000)
      If(ITqry(IOScan).ne.0) Call FileIO(2,-IOScan,InToWP(2),JScan,0)
      Call GetPBC(NPDir,CelVec)
      NAt3R = 3*NAtoms
      NA3RTT = (NAt3R*(NAt3R+1))/2
      NAtoms = NAtoms + NPDir
      Call ILSW(1,23,1)
      Call GetNOp(I,J)
      Call AMove(9,TRot,TRotT)
      Call Trspn2(3,TRotT)
      NAt3 = 3*NAtoms
      NAtNew = NAtoms
      Call AMove(NAt3,C,CSave)
      Call ILSW(2,26,NoSym)
      Call ILSW(2,56,IGrSym)
      NZ = 0
      NVrMMT = 0
      If(ITqry(IOZRed).gt.0) then
        Call FileIO(2,-IOZRed,InToWP(4),NTRedS,0)
      else
        Call IClear(4,NTRedS)
        endIf
      NVarZM = NTRedS(1)
      NTBond = NTRedS(2)
      NTAng  = NTRedS(3)
      NTDih  = NTRedS(4)
      DoRed = NVarZM.gt.0
      ICVer = 0
      Call CrdVer(0,ICVer)
      If(IPrint.ge.1) Write(IOut,4070) ICVer
      NewCrd = ICVer.gt.1.and.DoRed
      AlUnit = IGetFC.eq.10.or.(JustCk.and.NPDir.gt.0.and..not.DoRed)
      If(MicOTy.eq.0) then
        If(DoRed) then
          MicOTy = 3
        else
          MicOTy = 2
          endIf
        endIf
      DoMicO = MicOTy.eq.1.and.IDirMM.ne.6
      OptMic = MicOTy.ne.2.and.MicOTy.ne.4
      NFrag = IOp(62)
      SavLMO = NFrag.le.1.and.NPDir.eq.0.and.UpdRWF
      IDoEle = Mod(IOp(52),100)/10
      IChScl = IOp(74)
      DoEle = IDoEle.eq.1.or.IDoEle.eq.3.or.IDoEle.eq.4.or.
     $  (IDoEle.eq.0.and.IChScl.gt.0)
      DoEle = DoEle.and.MicOTy.eq.3
      ONLnOK = IOp(124).eq.0
      IScMMF = 0
      MaxMic = IOp(139)
      IGDRFO = IOp(141)
      If(IGDRFO.eq.0) IGDRFO = 1
      NGoDwn = IOp(151)
      If(NGoDwn.eq.0) NGoDwn = 3
      IMxStp = IOp(152)
      IUseOT = IOp(153)
      UseOT = IUseOT.eq.1.or.(IUseOT.eq.0.and.(FCart.or.DoRed))
      IEStpM = IOp(157)
      IDLimC = IOp(159)
      If(IDLimC.eq.0) IDLimC = 30
      If(IDLimC.eq.-1) IDLimC = 10000000
      DxLimC = GFloat(IDLimC) / GFloat(10)
      IScMMC = IOp(165)
      If(IScMMC.eq.0) IScMMC = 10
      RawDIS = IOp(166).ge.0
      MxGDII = Mod(IAbs(IOp(166)),100)
      MxRFOD = Mod(IAbs(IOp(166)),10000)/100
C     If(MxRFOD.eq.0) MxRFOD = 5
      IGEVer = IOp(167)
      If(IGEVer.eq.0) IGEVer = 1
      If(IGEVer.ne.1) Call GauErr('Bad IOp(167) in L103.')
      ISwDIS = IOp(168)
      ISwGDI = Mod(ISwDIS,100)
      If(ISwGDI.eq.0) then
        SwitMx = GFloat(1) / GFloat(1000)
      else
        SwitMx = Ten**(-ISwGDI)
        endIf
      ISimSw = Mod(ISwDIS,10000)/100
      If(ISimSw.eq.0) then
        If(IMDS.eq.1.and.Large) then
C         Old GDIIS and large
          SimSw = GFloat(75) /GFloat(100000)
        else if(IMDS.eq.1) then
C         Old GDIIS
          SimSw = GFloat(50) / GFloat(10000)
        else
          SimSw = GFloat(25) / GFloat(10000)
          endIf
      else
        SimSw = Ten**(-ISimSw)
        endIf
      ICThr = Mod(ISwDIS,10000000)/10000
      If(ICThr.eq.0) ICThr = 100
      CThrHd = GFloat(ICThr) / GFloat(10)
      ICDAdj = ISwDIS/10000000
      If(ICDAdj.eq.0) ICDAdj = 30
      CDAdj = GFloat(ICDAdj) / GFloat(10)
      NTRed = NVarZM
      NVar = NVarZM
      IMicOp = 1
      ICStOr = IMicOp + InToWP(NAtoms)
      ICA = ICStOr + NAt3
      ICB = ICA + NAt3
      IIAn = ICB + NAt3
      IAtChg = IIAn + InToWP(NAtoms)
      IFroze = IAtChg + NAtoms
      IFrozA = IFroze + InToWP(NAtoms)
      IIntVc = IFrozA + InToWP(NAtoms)
      IFPVec = IIntVc + InToWP(Max(MaxZSb,NVar))
      IEnd = IFPVec + Max(MaxZSb,NVar) - 1
      Call TstCor(IEnd,MDV,'Optmz-0')
      Call AMove(NAt3,C,V(ICA))
      Call IMove(NAtoms,IAn,V(IIAn))
      Call AMove(NAtoms,AtmChg,V(IAtChg))
      Call IMove(NAtoms,MicOpt,V(IMicOp))
      NAtomA = NAtoms
      Call FrzAt1(1,NAtoms,NAtomA,0,IAn,IAtTyp,MicOpt,NFroz,V(IFroze),
     $  NFrozA,V(IFrozA))
      If(OptMic) Call ISubsI(NAtoms,V(IMicOp),-3,1)
      HaveRB = IArMax(.False.,NAtoms,V(IMicOp)).gt.0
      NVarMM = 0
      NRF = 0
      NVrCQM = 0
      RotQM = .False.
      If(IDir.ge.2.or.OptMic) then
        NAtMO = ICount(0,NAtoms,V(IMicOp))
        NAtMic = NAtMO + ICount(1,NAtoms,V(IMicOp))
        OptMic = OptMic.and.NAtMic.gt.0
        If(IDir.ge.2.or.OptMic) then
          NVarMM = 3*NAtMO
          NVrCQM = 3*ICouQM(NAtoms,V(IMicOp))
          Call GetNRF(NAtoms,NRF,RotQM,V(IMicOp))
          endIf
        endIf
      NVarQ = 0
C     If(PCM(ISCRF)) Call PCMGet('NVQ',NVarQ,RJunk)
      NVrMMT = NVarMM + 6*NRF + NVarQ
      NVarT = NVar + NVrMMT
      ICur = 0
      Call PutDel(5,ANames,ICur)
      NAtMic = 0
C
C     Allocate space for the internal coordinates and, if requested, load
C     the data from the RWF file.
C
      If(DoRed) then
        IIZRed = IEnd + 1
        IMicO1 = IIZRed + InToWP(4*NVar)
        IIZRdO = IMicO1 + InToWP(NAtoms)
        IValue = IIZRdO + InToWP(4*NVar)
        ICnst1 = IValue + NVar
        ICnst2 = ICnst1 + NVar
        IMicIn = ICnst2 + NVar
        IEnd = IMicIn + 2*InToWP(NAtoms)-1
C
C       Allocate extra space for copies of Cartesian and internal coordinates
C       of reactant,product and/or TS
C
        If(NPath.gt.0) then
          ICrt1 = ICnst2 + NVar
          If(DidQST) then
            ICnst3 = ICnst2 + NVar
            ICrt1 = ICnst3 + NVar
            endIf
          ICrt2 = ICrt1 + 3*NAtoms
          IEnd = ICrt2 + 3*NAtoms - 1
          If(DidQST) then
            ICrt3 = ICrt2 + 3*NAtoms
            IEnd = ICrt3 + 3*NAtoms - 1
            endIf
          endIf
        Call TstCor(IEnd+3*NAtoms,MDV,'Optmz-Red')
        Call IMove(NAtoms,V(IMicOp),V(IMicO1))
        Call FileIO(2,IOZRed,InToWP(4*NVar),V(IIZRed),0)
        Call FileIO(2,IOZRed,NVar,V(IValue),0)
        Call FileIO(2,IOZRed,InToWP(NVar),V(IIntVc),0)
        Call FileIO(2,IOZRed,NVar,V(IFPVec),0)
        Call FileIO(2,IOZRed,2*NVar,V(ICnst1),0)
        Call IMove(4*NVar,V(IIZRed),V(IIZRdO))
        Call ChkMic(1,NAtoms,NAtNew,NAtMic,V(IMicO1),V(IMicIn),V(ICA),
     $    NTRed,V(IIZRed),V(IIAn),V(IAtChg),V(IEnd+1),V(IEnd+1))
        DoMicO = DoMicO.and.NAtMic.gt.0.and.NVrMMT.gt.0
        Call NDimZR(NVar,V(IIZRed),NDim,Junk)
        If(NPath.gt.0) then
C         load the cartesian coordinates for Reactant and Product
          If(DidQST) Call FileIO(2,IOZRed,NVar,V(ICnst3),0)
          Call FileIO(2,IOZRed,3*NAtNew,V(ICrt1),0)
          Call FileIO(2,IOZRed,3*NAtNew,V(ICrt2),0)
          If(DidQST) Call FileIO(2,IOZRed,3*NAtNew,V(ICrt3),0)
          endIf
      else
        IIZRed = IEnd   + 1
        IIZRdO = IIZRed
        IEnd   = IIZRed + InToWP(4) - 1
        Call ISet(0,4,V(IIZRed))
        Call ISet(-1000,1,V(IIZRed))
        IMicO1 = IEnd   + 1
        IValue = IEnd   + 1
        ICnst1 = IEnd   + 1
        ICnst2 = IEnd   + 1
        IMicIn = IEnd   + 1
        ICrt1  = IEnd   + 1
        ICrt2  = IEnd   + 1
        Call IClear(MaxZSb,V(IIntVc))
        Call AClear(MaxZSb,V(IFPVec))
        endIf
      Call CkTrRo(DoRed,NAtoms,NVar,NTRedS(2),V(IMicOp),V(IIZRdO),
     $  NTrRot,NFrz)
C
C     Get things set-up for the new internal coordinate code.
C
      ICrDat = IEnd+1
      If(NewCrd) then
        If(NPDir.gt.0) Call GauErr('Optmz: New coord code w/PBC NYI.')
        Call CrdDat(IOut,IPrint,IRWCrd,0,NCoord,ICrDat,IEnd,V,MDV)
        If(NCoord.ne.NVar) Call GauErr('Optmz: NCoord.ne.NVar')
        IEnd = IEnd - 1
        MDV1 = MDV - IEnd
        If(IPrint.ge.3) Call CrdPrt(IOut,IPrint,0,IRWCrd,0,V(ICrDat))
        endIf
C
C     Computation of number of independent degrees of freedom
C     isn't correct for PBC
C
      If(NPDir.gt.0) NtrRot = -3
      If(MthUpd.eq.-3) then
        If(Large.and.LSTQST) then
          MthUpd = 6
        else if(NPath.gt.0) then
          MthUpd = 3
        else if(Large) then
          MthUpd = 6
        else if(DoRed.and.NNeg.ne.0) then
          MthUpd = 4
        else if(DoRed) then
          MthUpd = 3
        else
          MthUpd = 1
          endIf
        endIf
C
C     Read in data and initialize optimization.
C
      If(First) then
        Write(IOut,1010)
        DoClnS = .False.
        DoClnG = .False.
        If(NVarZM.eq.0) Call RdWrZ(7,0,V)
        If(NVarZM.eq.0.and.SetStp) Write(IOut,1060)
        FCart = FCart.or.NVarZM.eq.0
        If(FCart) then
          If(DoRed) Call GauErr('FCart and DoRed #1.')
          If(ITqry(IOZRed).gt.0.and.MicOTy.ne.4)
     $      Call FileIO(5,IOZRed,0,V,0)
          If(NAtNew.le.0.and.NZ.ne.0) then
            NAtNew = 0
            Do 15 I = 1, NZ
              If(IAnZ(I).ge.0) NAtNew = NAtNew + 1
   15         Continue
            endIf
          NAtoms = NAtNew
          NVarZM = 3*NAtNew
          NZ = 0
          If(.not.SetStp) MethQ = 4
          endIf
        If(NVarZM.eq.0)
     $    Call GauErr('There are no variables to optimize!')
        If(DoRed.or.FCart) then
          Call CBlank(MaxZSb*IntPWP(0),ANames)
          ICur = 0
          Call PutDel(5,ANames,ICur)
          If(.not.DoRed) then
            Call IClear(MaxZSb,V(IIntVc))
            Call AClear(MaxZSb,V(IFPVec))
            endIf
          Call AClear(MaxZSb,Values)
        else
          Call IMove(NVarZM,IntVec,V(IIntVc))
          Call AMove(NVarZM,FPVec,V(IFPVec))
          endIf
        If(FCart.and.IDir.eq.5) then
          NVarZX = 3*ICouQM(NAtoms,V(IMicOp))
        else
          NVarZX = NVarZM
          endIf
        NVTT = (NVarZX*(NVarZX+1))/2
        If(NPath.eq.0) then
C         For QuadMac, XX stores some temporary data when calling Star
          MinXX = 3*NVarZM + 18*NAtoms + 1
          Call Ind103(IEnd+1,DoRed,AlUnit,JustCk,MicOTy,NVarZM,NVrMMT,
     $      MStepU,IMxStp,NAtNew,MaxStp,MaxSv,MaxSv2,MinXX,IIC,IX,
     $      IXName,IESave,IFS,IFC1,IFC2,IXX,IFF,IEQMS,IXXMM,IFFMM,
     $      IFFMMQ,IITU,LOptG2,IForce,IFrcns,IEnd)
        else
C         Bucket containing information passed from step to step.
C         Needs reorganization and a few arrays can be made volatile instead.
          Call InPath(IEnd+1,DoRed,NVarZM,NPath,MStepU,IMxStp,NAtNew,
     $      MaxStp,IITU,IIC,IX,INCurr,IXP,IESvP,IFSP,IFC1P,IFC2P,IXXP,
     $      IFFP,IICP,IXName,IESave,IFS,IFC1,IFC2,IXX,IFF,LOptG2,IForce,
     $      IFrcns,IFTStar,IVEcStar,INStar,IRisP,ICP,IConv,IRad,ITRot,
     $      IIRis,IINOp,ISym,IIPrA,Itan,IArcP,IPthL,IEnd)
            Call IClear(NPath,V(IRisP))
          endIf
        Call TstCor(IEnd+3*NAtNew,MDV,'Optmz1')
        N = NVarZM
        If(.not.AlUnit) Call AClear(NVarZM*NVarZM,V(IFC1))
        If(DoRed) then
          Call AMove(NVarZM,V(IValue),V(IX))
          Call IMove(NVarZM,V(IIntVc),V(IIC))
          If(.not.AlUnit) Call APDSQ(1,NVarZM,V(IFPVec),V(IFC1))
          Call RedSym(NVarZM*IntPWP(0),NTRedS(2),NTRedS(3),NTRedS(4),
     $      V(IXName),V(IIZRdO))
        else if(FCart) then
          Call CartVN(NVarZM,V(IXName))
        else
          Call IMove(MaxZSb*IntPWP(0),ANames,V(IXName))
          Call AMove(NVarZM,Values,V(IX))
          Call IMove(NVarZM,V(IIntVc),V(IIC))
          Call APDSq(1,NVarZM,V(IFPVec),V(IFC1))
          endIf
        Call InitBS(IOp,ToAng,AlUnit,MaxStp,N,MethQ,IGEVer,FCart,ConvP,
     $    IRwSM,V(IX),V(IFrcns),V(IFC1),V(IXName),V(IIC),V(IITU),DoRed,
     $    IAprBG,LSTQST,NTrRot,NTRedS(2),NTRedS(3),NTRedS(4),V(IIZRed),
     $    V(IValue),V(ICnst1),V(ICnst2),NAtoms,V(IIAn),V(ICA),MxGDII,
     $    MaxGDI,V(IEnd+1),MDV-IEnd)
        If(OptMic) Call ChkMic(2,NAtoms,NAtNew,NAtMic,V(IMicO1),
     $    V(IMicIn),V(ICA),NTRed,V(IIZRed),V(IIAn),V(IAtChg),V(IEnd+1),
     $    V(IEnd+1))
        Call ConDDF(IOptGr,LOptGr+LOptG2)
        If(AlUnit) then
          LFrcns = 0
        else
          LFrcns = NVTT
          endIf
        Call ConDDF(IGrdnt,LGrdnt+NVar+LFrcns)
        MaxSV1 = Max(NVarZM,MaxStp)+1
C       Save reaction path information for later retrieval
        If(NPath.gt.0) then
          If(Do601) then
            LOpath = 2 + NPath*(2+2+3*NAt3+(NAt3*(NAt3+1))/2)
            Call ConDDF(IOPath,LOPath)
            endIf
          Call AClear(NPath,V(IESvp))
          Call AClear(MaxSv1*NPath,V(IFSp))
          Call AClear(2*NVTT*NPath,V(IFC1P))
          Call AClear(2*NVar*MaxSv1*NPath,V(IXXP))
          Call AClear(NVar*NPath,V(IXP))
          Call IClear(1+NPath*MAxSv1,V(INCurr))
          Call AClear((NVTT+1)*NPath,V(IFTStar))
          Call IClear(NPath+1,V(INStar))
          Call AClear(6*NAtNew*NPath,V(ICP))
          Call IClear(NPath,V(IConv))
          Call AClear(10*NPath,V(IRad))
          Call IClear(4,V(IIRis))
          Call AClear(NVar*NPath,V(ITan))
          Call AClear(NPath,V(IArcP))
          Call AClear(2,V(IPthL))
          Call IClear(NPath,V(IINop))
          Call IClear(NAtNew*MaxOp*NPath,V(IIPra))
        else if(Do601) then
          IL2 = InToWP(2)
          LBlock = IL2 + 1+2*NAt3
          If(Do601H) LBlock = LBlock + (Nat3*(Nat3+1))/2
          LOpath = IL2 + NStep*LBlock
          Call ConDDF(IOPath,LOPath)
          IntFL2(1) = 0
          IntFL2(2) = LBlock
          Call FileIO(1,-IOPath,IL2,IntFL2,0)
          endIf
        Energy = Zero
        If(.not.SetStp) Call FileIO(1,-IGrdnt,1+InToWP(2),Energy,0)
        Call IniRSc(IOut,IOScan,IOGenA,LRwGen,NVar,V(IIntVc),JScan)
        Write(IOut,1050) NStep, MaxStp
        If(NNeg.gt.0) Write(IOut,1080) NNeg
        MaxJb1 = Max(JScan(2),1)
        If(UseOT) Call RdWrOT(MaxStp,MaxJb1,NStat,NAt3,IJOff1,IType1,
     $    IStep1,IJob1,IDone,.True.,Stats,.True.,C,.False.,XX,.True.,XX,
     $    IOut,0,0,V(IEnd+1),MDV-IEnd)
C
C       Allocate for GrdOpt.
C
        MinLen = 20
        NAtN3 = 3*NAtNew
        NAtN3T = (NAtN3*(NAtN3+1))/2
        NAlloc = Max(NVarZM,NAtN3,MinLen)
        NAllT2 = Max((NAlloc*(NAlloc+1))/2,NAtN3*(NAtN3+3))
        If(AlUnit) NAllT2 = MinLen*NAlloc
        IXQuad = IEnd + 1
        IXNew = IXQuad + NAlloc + NVrMMT
        IFTemp = IXNew + NAlloc + NVrMMT
        IA = IFTemp + NAlloc + NVrMMT
        If(Large.or.AlUnit) then
          LenA = MinLen*NAlloc
        else
          LenA = NAlloc*Max(NAlloc,MaxStp)
          endIf
        IB = IA + LenA
        If(Large) then
          LenB = MinLen*NAlloc
        else if(DoRed.or.(FCart.and..not.AlUnit)) then
          LenB = NAlloc**2
        else
          LenB = 0
          endIf
        ISvEig = IB + LenB
        IScr = ISvEig + NAlloc
        IDimSc = Max(NAlloc,MaxStp)
        IAPrt = IScr + MinLen*IDimSc
        IFCT = IAPrt + 6*NAlloc
        ICNew = IFCT + NAllT2
        IIBMat = ICNew + NAt3
        If(DoRed) then
          LenB1 = LenRIC(IOut,IPrint,1,NVar,NAtoms,NDim,V(IIZRed))
          LenIB1 = LenRIC(IOut,IPrint,2,NVar,NAtoms,NDim,V(IIZRed))
        else
          LenB1 = 0
          LenIB1 = 0
          endIf
        If(NewCrd) then
          LenB2 = LenRIC(IOut,IPrint,1,NVar,NAtoms,NDim,V(ICrDat))
          LenIB2 = LenRIC(IOut,IPrint,2,NVar,NAtoms,NDim,V(ICrDat))
        else
          LenB2 = 0
          LenIB2 = 0
          endIf
        LenB = Max(LenB1,LenB2)
        LenIB = Max(LenIB1,LenIB2)
        IBMat = IIBMat + IntoWP(LenIB)
        IEnd = IBMat + LenB + 9*NAtNew
        If(NPath.gt.0) then
          NGoal = Min(2*NVarZM,6*NAtNew)
          MaxStp = MxStSt(MStepU,0,NGoal,IMxStp)
          MaxSV1 = Max(NVarZM,MaxStp) + 1
          ICns12 = IEnd
          ICns22 = ICns12 + NVarZM
          ICrt12 = ICns22 + NVarZM
          If(DidQST) then
            ICns32 = ICns22 + NVarZM
            ICrt12 = ICns32 + NVarZM
            endIf
          ICrt22 = ICrt12 + NAt3
          IIC2 = ICrt22 + NAt3
          If(DidQST) then
            ICrt32 = ICrt22 + NAt3
            IIC2 = ICrt32 + NAt3
            endIf
          IQ = IIC2 + 1
          IDQ = IQ + NVarZM*NPath
          IQ0 = IDQ + NVarZM*NPath
          IGrad = IQ0 + NVarZM*NPath
          IAP = IGrad + NVarZM*NPath
          IEigP = IAP + NVarZM*NVarZM*NPath
          IEnd = IEigP + NVarZM*NPath + 1
          endIf
        Call TstCor(IEnd,MDV,'Optmz4')
        NVar = NVarZM
        If(NPath.gt.0) then
          If(NCurr.eq.0) Call IMove(1,V(INCurr),NCurr)
C         Interpolate the first path.
          If(.not.DidQST) then
            Call AMove(NVar,V(ICnst2),V(ICns22))
            Call AMove(NVar,V(ICnst1),V(ICns12))
            Call AMove(NAt3,V(ICrt1),V(ICrt12))
            Call AMove(NAt3,V(ICrt2),V(ICrt22))
            NStop = NPath
C           PthSpd : Spreads equidistant points on a line (equidistant
C           in Redundant internals)
            NStTmp = 1
            Call PthSpd(IOut,IPrint,NAtNew,V(ICrt22),V(ICrt12),
     $        V(IFtemp),V(IXXP),V(IIZRed),NVar,NTRedS(2),
     $        NTRedS(3),NTRedS(4),V(ICns12),V(ICns22),V(IXQuad),
     $        V(IBMat),V(IIBMat),V(IB),V(IFCT),V(IScr),V(ICns12),
     $        V(IIntVc),MaxSv1,NPath,V(IXNew),V(ICP),DoReac,DoProd,
     $        NStTmp,NStop,V(IEnd+1),MDV-IEnd)
          else
            Call AMove(NVar,V(ICnst3),V(ICns32))
            Call AMove(NVar,V(ICnst1),V(ICns12))
            Call AMove(NAt3,V(ICrt1),V(ICrt12))
            Call AMove(NAt3,V(ICrt3),V(ICrt32))
            NPath1 = (NPath+1)/2
C           First : From Reactant to TS
            NStTmp = 1
            Call PthSpd(IOut,IPrint,NAtNew,V(ICrt32),V(ICrt12),
     $        V(IFtemp),V(IXXP),V(IIZRed),NVar,NTRedS(2),
     $        NTRedS(3),NTRedS(4),V(ICns12),V(ICns32),V(IXQuad),
     $        V(IBMat),V(IIBMat),V(IB),V(IFCT),V(IScr),V(ICns12),
     $        V(IIntVc),MaxSv1,NPath,V(IXNew),V(ICP),DoReac,.true.,
     $        NStTmp,NPath1,V(IEnd+1),MDV-IEnd)
            Call AMove(NVar,V(ICnst2),V(ICns22))
            Call AMove(NVar,V(ICnst3),V(ICns32))
            Call AMove(NAt3,V(ICrt3),V(ICrt32))
            Call AMove(NAt3,V(ICrt2),V(ICrt22))
            NPath1 = (NPath+1)/2
            NStop = NPath
C           Second : From product to TS (rather than TS to product to
C           help conserve symmetry).
            Call PthSpd(IOut,IPrint,NAtNew,V(ICrt32),V(ICrt22),
     $        V(IFtemp),V(IXXP),V(IIZRed),NVar,NTRedS(2),
     $        NTRedS(3),NTRedS(4),V(ICns22),V(ICns32),V(IXQuad),
     $        V(IBMat),V(IIBMat),V(IB),V(IFCT),V(IScr),V(ICns22),
     $        V(IIntVc),MaxSv1,NPath,V(IXNew),V(ICP),DoProd,.true.,
     $        NStop,NPath1,V(IEnd+1),MDV-IEnd)
            endIf
          Do 110 Npt=1,NPath
            Call AMove(NVar,V(IXXP+NVar*MaxSv1*(NPt-1)),
     $        V(IXP+NVar*(NPt-1)))
  110       Continue
          Call PtPath(IOut,NVar,ToAng,NPath,V(IXP),V(IEsvP),V(IXName),
     $      NTRedS(2),.False.)
          Call ChkAng(IOut,NPath,NTRedS(1),NTRedS(2),NTReds(3),
     $      V(IIZRed),V(IXP))
          Call ISet(1,1,V(INCurr))
          endIf
C
C     Second and later calls to L103.
C
      else
        DoClnS = NoSym.eq.0.and..not.HaveRB.and.IGrSym.eq.0
     $    .and.NOpAll.gt.1
        DoClnG = DoClnS.and..not.OptMic
        If(.not.FCart.and..not.DoRed) then
          If(NVarZM.eq.0) Call RdWrZ(7,0,V)
          FCart = FCart.or.NVarZM.eq.0
          endIf
        If(FCart.or.DoRed) then
          Call AMove(NAt3,C,V(ICStOr))
          Call RotF1(NAtoms,TRot,V(ICA))
          Call Transl(-1,3,NAtoms,TrVec,V(ICA))
          Call RotF1(NAtoms,TRot,C)
          Call Transl(-1,3,NAtoms,TrVec,C)
          endIf
        If(OptMic) Call ChkMic(2,NAtoms,NAtNew,NAtMic,V(IMicO1),
     $    V(IMicIn),V(ICA),NTRed,V(IIZRed),V(IIAn),V(IAtChg),V(IEnd+1),
     $    V(IEnd+1))
        If(FCart.or.DoRed) Call RdWrZ(16,0,V)
        If(FCart) then
          If(DoRed) Call GauErr('FCart and DoRed #2.')
          NVarZM = 3*NAtNew
          If(.not.SetStp) MethQ = 4
        else if(.not.DoRed) then
          Call IMove(NVarZM,IntVec,V(IIntVc))
          Call AMove(NVarZM,FPVec,V(IFPVec))
          endIf
        If(FCart.and.IDir.eq.5) then
          NVarZX = 3*ICouQM(NAtoms,V(IMicOp))
        else
          NVarZX = NVarZM
          endIf
        NVTT = (NVarZX*(NVarZX+1))/2
        If(NPath.eq.0) then
C         For QuadMac, XX stores some temporary data when calling Star
          MinXX = 3*NVarZM + 18*NAtoms + 1
          Call Ind103(IEnd+1,DoRed,AlUnit,JustCk,MicOTy,NVarZM,NVrMMT,
     $      MStepU,IMxStp,NAtNew,MaxStp,MaxSv,MaxSv2,MinXX,IIC,IX,
     $      IXName,IESave,IFS,IFC1,IFC2,IXX,IFF,IEQMS,IXXMM,IFFMM,
     $      IFFMMQ,IITU,LOptG2,IForce,IFrcns,IEnd)
        else
          Call INPath(IEnd+1,DoRed,NVarZM,NPath,MStepU,IMxStp,NAtNew,
     $      MaxStp,IITU,IIC,IX,INCurr,IXP,IESvP,IFSP,IFC1P,IFC2P,IXXP,
     $      IFFP,IICP,IXName,IESave,IFS,IFC1,IFC2,IXX,IFF,LOptG2,IForce,
     $      IFrcns,IFTStar,IVecStar,INStar,IRisP,ICP,IConv,IRad,ITRot,
     $      IIRis,IINOp,ISym,IIPrA,ITan,IArcP,IPthL,IEnd)
          endIf
        Call FileIO(2,-IGrdnt,LGrdnt,Energy,0)
        Call FileIO(2,-IOptGr,LOptGr,ConvF,0)
        IError = IOp(26)
        Call InFGEr(IOut,IPrint,.False.,DoRest,IError,NAtoms,FncErr,
     $    GrdErr)
        If(UseOT.and.NP.ge.0) then
          IFxOT = IEnd + 1
          IV2 = IFxOT + NAt3
          Call TstCor(IV2,MDV,'Optmz-RdWrOT-Save1')
          Call FileIO(2,-IRwF1,NAt3,V(IFxOT),0)
          Call RotF1(NAtoms,TRotT,V(IFxOT))
          Stats(1) = Energy
          Stats(2) = Zero
          IDone = IStep + 1
          Junk1 = 0
          Junk2 = 0
          MaxJb1 = Max(JScan(2),1)
          IJob1 = Max(JScan(1),1)
          Call RdWrOT(MaxStp,MaxJb1,NStat,NAt3,Junk1,Junk2,IStep+1,
     $      IJob1,IDone,.True.,Stats,.True.,V(ICStOr),.False.,XX,.True.,
     $      V(IFxOT),IOut,0,2,V(IV2),MDV-IV2+1)
          endIf
        If(Large) IOSLE1 = IOSLEq
C       A hack to handle restarting from old format /OptGr/.
        LOptGX = LOptG2 - InToWP(MaxSv)
        If(ITqry(IOptGr).eq.(LOptGr+LOptGX)) then
          Call FileIO(2,IOptGr,LOptGX,V(IIC),0)
          Call IClear(MaxSv,V(IITU))
          Call ConDDF(IOptGr,LOptGr+LOptG2)
          Call FileIO(1,-IOptGr,LOptGr,ConvF,0)
          Call FileIO(1,IOptGr,LOptGr,V(IIC),0)
        else
          Call FileIO(2,IOptGr,LOptG2,V(IIC),0)
          endIf
        If(DoRed) Call AMove(NVarZM,V(IValue),V(IX))
        If(NP.gt.0.and.IMDS.ne.0) then
          IMDSC = Max(IMDS,1)
          Write(IOut,1040) NamMDS(IMDSC)(1:LinEnd(NamMDS(IMDSC)))
          endIf
        IStart = IEnd
        IScGDI = IEnd + 1
        If(IMDS.ne.0) then
          NNP = Max(NP,NGDIIS) + 1
          NNP = Max(NNP,0)
          If(IMDS.eq.1) then
            Call IndGDI(NNP,NVarZM,Large,IStart,IQQS,IDQQS,IIA,IAA,ICC,
     $        ICSave,IHSc,IDQD,IDQR,IW,IEnd)
            Call TstCor(IEnd,MDV,'Optmz-IndGDI')
            IScGDI = IQQS
          else
            Call SimMem(NNP,NVarT,IStart,IDQQS,IIA,IIAInv,IIBB,IICC,
     $        IICM,IIHSc,IIDQR,IIW,IEnd)
            Call TstCor(IEnd,MDV,'Optmz-Mixed/SimOpt')
            IRwSMM = IRwSM
            IScGDI = IDQQS
            endIf
          endIf
        If(DoRed) Call RedSym(NVarZM*IntPWP(0),NTRedS(2),NTRedS(3),
     $    NTRedS(4),V(IXName),V(IIZRdO))
        If(FCart.and..not.DoRest) then
          Call FileIO(2,-IRWF1,NVarZM,V(IForce),0)
          If(.not.AlUnit) then
            Call AClear(NVTT,V(IFrcns))
            If(IGetFC.eq.8) then
              Call APlusI(NVarZM,One,V(IFrcns))
              Call AMove(NVTT,V(IFrcns),V(IFC1))
            else
              LF2 = ITqry(IRwF2)
              If(LF2.gt.0)
     $          Call FileIO(2,-IRwF2,Min(LF2,NA3RTT),V(IFrcns),0)
              If(LF2.eq.NA3RTT) then
                Call AClear(NVTT-NA3RTT,V(NA3RTT+1))
                Do 260 I = (NAt3R+1), NVarZM
                  II = (I*(I+1))/2
  260             V(II) = One
              else if(LF2.ne.NVTT.and.LF2.ne.0) then
                Write(IOut,3060) LF2, NVTT
                Call GauErr('Inconsistency in Optmz.')
                endIf
              endIf
            endIf
          Call ANeg(NVarZM,V(IForce),V(IForce))
        else
          Call FileIO(2,IGrdnt,NVarZM,V(IForce),0)
          If(.not.AlUnit) Call FileIO(2,IGrdnt,NVTT,V(IFrcns),0)
          endIf
        If(DoRed.and.IExit.eq.1) then
          IFX = IEnd + 1
          IScr = IFX + 3*Max(NAtNew,NAtoms)
          Call TstCor(IScr+3*NAtoms,MDV,'Optmz-FXRed')
          Call FileIO(2,-IRWF1,3*NAtoms,V(IFX),0)
          If(OptMic) Call ChkMic(3,NAtoms,NAtNew,NAtMic,V(IMicO1),
     $      V(IMicIn),V(IFX),JJ(1),JJ,JJ,XX,JJ,V(IScr))
          IniFTr = IStep.eq.0
          If(IGISav.eq.0.and.IOSLE1.eq.0) IGISav = -1
          IFlMin = IArMin(.False.,NVar,V(IIntVc))
          If(Large.or.LSTQST.or.IFlMin.lt.0) IAprBG = (IAprBG/10)*10
          If(NewCrd) then
            Call FXRedC(IOut,IPrint,IOSLE1,IniFTr,IAprBG,IGISav,IGIOld,
     $        NAtNew,NTrRot,NTRed,NTBond,NTAng,NTDih,V(ICrDat),
     $        V(IIntVc),V(IX),V(ICA),V(IFX),V(IForce),IFail,LenSB,LenCB,
     $        V(IXName),V(IScr),MDV-IScr+1)
          else
            Call FXRedC(IOut,IPrint,IOSLE1,IniFTr,IAprBG,IGISav,IGIOld,
     $        NAtNew,NTrRot,NTRed,NTBond,NTAng,NTDih,V(IIZRed),
     $        V(IIntVc),V(IX),V(ICA),V(IFX),V(IForce),IFail,LenSB,LenCB,
     $        V(IXName),V(IScr),MDV-IScr+1)
            endIf
          If(IFail.ne.0) then
            Write(IOut,4050)
            Call GauErr('Incomplete coordinate system.')
            endIf
          Call ANeg(NVarZM,V(IForce),V(IForce))
          endIf
        If(FCart) Call AMove(NVarZM,V(ICA),V(IX))
        DoEstm = IGetFC.ne.3.and..not.AlUnit
        If(DoRest) then
          If(IOp(6).ne.0) then
            NStep = IAbs(IOp(6))
            ITemp = MaxStp
            MaxStp = MxStSt(ITemp,0,NStep,IMxStp)
            Write(IOut,1050) NStep, MaxStp
            MaxJb1 = Max(JScan(2),1)
            IJob1 = Max(JScan(1),1)
            IDone = IStep + 1
            If(UseOT) Call RdWrOT(MaxStp,MaxJb1,NStat,NAt3,IJOff1,
     $        IType1,IStep1,IJob1,IDone,.True.,Stats,.True.,C,.False.,
     $        XX,.True.,XX,IOut,0,-1,V(IEnd+1),MDV-IEnd)
            endIf
        else
          IGetFC = IOp(10)
          If(DoChk) then
            Call FileIO(1,-IGrdnt,LGrdnt,Energy,0)
            Call FileIO(1,IGrdnt,NVarZM,V(IForce),0)
            If(.not.AlUnit) Call FileIO(1,IGrdnt,NVTT,V(IFrcns),0)
            Call ChkPnt(0,V(IEnd+1),MDV-IEnd)
            endIf
          endIf
        If(IGetFC.ne.4.and.IGetFC.ne.7.and..not.AlUnit)
     $    Call AMove(NVTT,V(IFC1),V(IFrcns))
        Call ILSW(1,27,0)
C
C       Allocate for GrdOpt.
C
        MinLen = 20
        NAtN3 = 3*NAtNew
        NAtN3T = NAtN3*(NAtN3+3)
        NAlloc = Max(NVarZM,NAtN3,MinLen)
        NAllT2 = Max(NAlloc**2,NAtN3T)
        If(AlUnit) NAllT2 = MinLen*NAlloc
        IXQuad = IEnd + 1
        IXNew = IXQuad + NAlloc + NVrMMT
        IFTemp = IXNew + NAlloc + NVrMMT
        IA = IFTemp + Max(NAlloc+NVrMMT,NAtoms*3)
        If(Large.or.AlUnit) then
          LenA = MinLen*NAlloc
        else
          LenA = NAlloc*Max(NAlloc,MaxStp)
          endIf
        IB = IA + LenA
        If(Large) then
          LenB = MinLen*NAlloc
        else if(DoRed.or.(FCart.and..not.AlUnit)) then
          LenB = NAlloc**2
        else
          LenB = 0
          endIf
        ISvEig = IB + LenB
        IScr = ISvEig + NAlloc
        LBSc = InToWP(2*NTRedS(1))+IntoWP(NTRedS(1)+2+LenCB+NAtNew)+
     $    3*(LenCB+NAtNew)
        IDimSc = Max(NAlloc,MaxStp,3*NAtoms,NVarT)
        IAPrt = IScr + Max(8*IDimSc,LBSc,3*NAtoms,NAlloc+NVrMMT)
        IFCT = IAPrt + 6*(NAlloc+NVrMMT)
        ICNew = IFCT + NAllT2
        IIBMat = ICNew + 3*NAtoms
        If(DoRed) then
          LenB1 = LenRIC(IOut,IPrint,1,NVar,NAtoms,NDim,V(IIZRed))
          LenIB1 = LenRIC(IOut,IPrint,2,NVar,NAtoms,NDim,V(IIZRed))
        else
          LenB1 = 0
          LenIB1 = 0
          endIf
        If(NewCrd) then
          LenB2 = LenRIC(IOut,IPrint,1,NVar,NAtoms,NDim,V(ICrDat))
          LenIB2 = LenRIC(IOut,IPrint,2,NVar,NAtoms,NDim,V(ICrDat))
        else
          LenB2 = 0
          LenIB2 = 0
          endIf
        LenB = Max(LenB1,LenB2)
        LenIB = Max(LenIB1,LenIB2)
        IBMat = IIBMat + InToWP(LenIB) + InToWP(NTRedS(1)+2+LenSB)
        IEnd = IBMat + LenB + 9*NAtNew+LenSB*3
        If(NPath.gt.0) then
          NGoal = 2*NVarZM
          NGoal = Min(NGoal,6*NAtNew)
          MaxStp = MxStSt(MStepU,0,NGoal,IMxStp)
          MaxSV1 = Max(NVarZM,MaxStp) + 1
          ICns12 = IEnd
          ICns22 = ICns12 + NVarZM
          ICrt12 = ICns22 +NVarZM
          If(DidQST) then
            ICns32 = ICns22 + NVarZM
            ICrt12 = ICns32 + NVarZM
            endIf
          ICrt22 = ICrt12 + 3*NAtNew
          IIC2 = ICrt22 + 3*NAtNew
          If(DidQST) then
            ICrt32 = ICrt22 + 3*NAtNew
            IIC2 = ICrt32 + 3*NAtNew
            endIf
          IQ = IIC2 + 1
          IDQ = IQ + NVarZM*NPath
          IQ0 = IDQ + NVarZM*NPath
          IGrad = IQ0 + NVarZM*NPath
          IAP = IGrad + NVarZM*NPath
          IEigP = IAP + NVarZM*NVarZM*NPath
          IScr2 = IEigP + NVarZM*NPath
          IEnd = IScr2 + 9*NVarZM
          Call IMove(1,V(INCurr),NCurr)
          If(SetStp.and.Do601) then
C           Save path information we have so far
            LOpatN = (2+2+3*Nat3+(Nat3*(Nat3+1))/2)
            LOpath = 2+NPath*LOpatN
C           number of points along the path
            WPath = GFloat(NPath)
            Call Fileio(1,-IOPath,1,WPath,0)
C           cuurent point
            WCurr = GFloat(NCurr)
            Call Fileio(1,IOPath,1,WCurr,0)
            IPoint = LOpatN*(NCurr-1)
C           point index. If the point is converged this inetger will be
C           the length of useable data (reaction coordinate+energy+gradient
C           +tangent+(possibly) second derivatives.
            Call Fileio(1,IOPath,1,WCurr,IPoint)
C           To know whether 2nd derivatives are there
            WGetFC = GFloat(IGetFC)
            Call Fileio(1,IOPath,1,WGetFC,0)
C           space reserved for reaction coordinate value
C           the value is set at "Big" until the reaction coordinate and tangent
C           are indeed stored
            Call Fileio(1,IOPath,1,Big,0)
C           store energy
            Call Fileio(1,IOPath,1,Energy,0)
C           cartesian coordinates
            Call Fileio(1,IOPath,Nat3,CSave,0)
C           forces
            Call FileIO(2,-IRWF1,Nat3,V(IFtemp),0)
            Call Fileio(1,IOPath,Nat3,V(IFtemp),0)
            If(IGetFC.eq.4) then
C             Skip space reserved for the unnormalized tangent
              IPoint = NAt3
C             second derivatives
              NAt3TT = (Nat3*(Nat3+1))/2
              Call FileIO(2,-IRWF2,Nat3TT,V(IFCT),0)
              Call Fileio(1,IOPath,Nat3,V(IFCT),IPoint)
              endIf
            endIf
          endIf
        Call TstCor(IEnd,MDV,'Optmz5')
        NVar = NVarZM
C
C       Now we are either doing the "star" (numerically differentiating
C       to get initial force constants), or doing an optimization step,
C       or transitioning from one phase to the other.  Star sets Exit if
C       the differentiation is complete; GrdOpt then resets it.
C
        If(NP.lt.0) then
          SkpMic = .True.
          IndMCR = 2 + 6*NRF + 6*NAtoms
          If(DoRed.and.IStep.eq.0) then
            Call AMove(3*NAtNew,V(ICA),V(IXX+2*NVarZM))
            Call AMove(NVarZM,V(IForce),V(IXX+2*NVarZM+3*NAtNew))
            Call FileIO(2,-IRWF1,3*NAtoms,V(IXX+3*NVarZM+3*NAtNew+
     $        3*NAtoms),0)
            If(ITqry(IMMCRS).gt.0) then
              Call FileIO(2,-IMMCRS,3*NAtoms,V(IXX+3*NVarZM+3*NAtNew),
     $          IndMCR)
              Call FileIO(2,IMMCRS,3*NAtoms,
     $          V(IXX+3*NVarZM+3*NAtNew+6*NAtoms),0)
              Call FileIO(2,-IMMCRS,1+6*NRF,V(IXX+3*NVarZM+3*NAtNew+
     $          9*NAtoms),1)
              endIf
            endIf
          If(IDir.ge.2) then
C           Star must get only the QM contribution to the gradient.
            NAtN3  = 3*NAtNew
            NAlloc = Max(NVarZM,NAtN3)
            IAG    = IEnd + 1
            IVal2  = IAG + NAlloc*NAlloc
            IScra1 = IVal2 + NVarZM
            IScra2 = IScra1 + Max(NVarZM,3*NAtoms)
            IVJQM = IScra2 + Max(NVarZM,3*NAtoms)
            Call TstCor(IVJQM,MDV,'Optmz-JustQM')
            Call JustQM(IOut,IMMCRS,NatNew,NAtoms,V(IMicOp),NVarZM,
     $        V(IIBMat),V(IBMat),V(IAG),V(IForce),NTrRot,V(ICA),
     $        NTRedS(2),NTRedS(3),NTRedS(4),V(IIZRed),V(IVal2),V(IB),
     $        V(IFCT),V(IScr),IPrint,V(IScra1),V(IScra2),V(IXName),
     $        IndMCR,V(IVJQM),MDV-IVJQM+1)
            endIf
          Call Star(DoEstm,AlUnit,NVarZM,NTRedS(2),NTRedS(3),NTRedS(4),
     $      V(IIZRed),V(IIAn),V(ICA),V(IIC),V(IForce),V(IESave),V(IFS),
     $      V(IX),V(IXX),V(IXX+NVarZM),V(IFC1),V(IFrcns),V(IScr))
          If(Exit) then
            SkpMic = .False.
            If(FCart) then
              Call AMove(NVarZM,V(IX),V(ICA))
            else if(DoRed) then
              Call AMove(3*NAtNew,V(IXX+2*NVarZM),V(ICA))
              Call FileIO(1,-IRWF1,3*NAtoms,V(IXX+3*NVarZM+3*NAtNew+
     $          3*NAtoms),0)
              If(ITqry(IMMCRS).gt.0) then
                Call FileIO(1,-IMMCRS,3*NAtoms,V(IXX+3*NVarZM+3*NAtNew),
     $            IndMCR)
                Call FileIO(1,IMMCRS,3*NAtoms,V(IXX+3*NVarZM+3*NAtNew+
     $            6*NAtoms),0)
                Call FileIO(1,-IMMCRS,1+6*NRF,V(IXX+3*NVarZM+3*NAtNew+
     $            9*NAtoms),1)
                J = 0
                Do 11 I = 1, NAtoms
                  If(IGetVl(V(IMicOp),I).lt.0) then
                    Call AMove(3,V(ICA+J),C(1,I))
                    J = J + 3
                    endIf
   11             Continue
                endIf
              Call AMove(NVarZM,V(IXX+2*NVarZM+3*NAtNew),V(IForce))
              Call AMove(NVarZM,V(IXX),V(IX))
            else
              Call UdZSub(V(IX),NVarZM)
              endIf
            If(DoChk) Call ChkPnt(0,V(IEnd+1),MDV-IEnd)
          else
            Write(IOut,1020)
            If(DoRed) then
              Call AMove(3*NAtNew,V(IXX+2*NVarZM),V(ICA))
              Call ASub(NVar,V(IX),V(IXX),V(IXNew))
              If(NewCrd) then
                Call RedCar(IOut,IPrt1,IOSLE1,NVar,NTRedS(2),NTRedS(3),
     $            NTReds(4),V(ICrDat),V(IXX),NAtNew,V(ICA),V(IXNew),
     $            V(ICNew),V(IBMat),V(IIBMat),V(IB),V(IFCT),V(IScr),
     $            V(ICnst1),V(IIntVc),.True.,IFail,0,V(IEnd+1),MDV-IEnd)
              else
                Call RedCar(IOut,IPrt1,IOSLE1,NVar,NTRedS(2),NTRedS(3),
     $            NTReds(4),V(IIZRed),V(IXX),NAtNew,V(ICA),V(IXNew),
     $            V(ICNew),V(IBMat),V(IIBMat),V(IB),V(IFCT),V(IScr),
     $            V(ICnst1),V(IIntVc),.True.,IFail,0,V(IEnd+1),MDV-IEnd)
                endIf
              If(IFail.eq.1) Call GauErr('RedCar failed.')
              Call AMove(3*NAtNew,V(ICNew),V(ICA))
              endIf
            endIf
          endIf
        If(Exit.or.NP.ge.0) then
          If(NPath.gt.0) then
            If(NewCrd)
     $        Call GauErr('Optmz: In Opt=Path w/NewCrd -- NYI.')
            Call PthOpt(IOut,IPrint,UpTrus,Sphere,DoLin,MethQ,FCart,
     $        SetStp,ToAng,IDimSc,V(IForce),V(IFrcns),V(IX),V(IXName),
     $        V(IFC1),V(IFC2),V(IESave),V(IFS),V(IXX),V(IFF),V(IIC),
     $        V(IITU),V(IXQuad),V(IXNew),V(IFTemp),V(IA),V(IB),V(IScr),
     $        V(IFCT),V(IAPrt),V(IIntVc),NAtNew,V(IIAn),V(IAtChg),
     $        V(ICA),DoClnG,TRotT,MthUpd,Wght,IQST,DoRed,NTrRot,
     $        NTRedS(2),NTRedS(3),NTRedS(4),V(IIZRed),V(IValue),
     $        V(ICnst1),V(ICnst2),V(ICNew),V(IIBMat),V(IBMat),NPath,
     $        NCurr,NNext,MaxSV1,NVTT,V(IICP),V(IFFP),V(IXXP),V(IFC1P),
     $        V(IFC2P),V(IESvP),V(IFSP),V(IQ0),V(IXP),V(IQ),V(IDQ),
     $        V(IGrad),V(IFTStar),V(IVecStar),V(IAP),V(IEigP),V(INStar),
     $        V(IRisP),V(ICP),V(IConv),DxMxT0,DoReac,DoProd,DoMid,
     $        DoBirx,DoIRC,V(IRad),V(ITRot),V(IIRis),V(IScr2),V(IINOp),
     $        V(ISym),V(IIPrA),V(ITan),V(IArcP),V(IPThL),IOSLE1,IAprBG,
     $        IEStpM,NFroz,V(IFroze),V(IEnd+1),MDV-IEnd)
            Call IMove(1,NNext,V(INCurr))
          else
            IDIS = IEnd
            If(IMDS.ge.2) IEnd = IEnd + NMax**2
C           Allocation specifically for coupled QM/MM step.
            If(NewCrd) then
              LenDB  = LenRIC(IOut,IPrint,3,NVarZM,NAtoms,1,V(ICrDat))
              LenIDB = LenRIC(IOut,IPrint,4,NVarZM,NAtoms,1,V(ICrDat))
            else if(DoRed) then
              LenDB  = LenRIC(IOut,IPrint,3,NVarZM,NAtoms,1,V(IIZRed))
              LenIDB = LenRIC(IOut,IPrint,4,NVarZM,NAtoms,1,V(IIZRed))
            else
              LenDB  = 0
              LenIDB = 0
              endIf
            IDBM   = IEnd   + 1
            IIDBM  = IDBM   + LenDB
            IFTmp  = IIDBM  + InToWP(LenIDB)
            IXMMT  = IFTmp  + NVarT
            IFMMQ  = IXMMT  + NVrMMT
            IHMod  = IFMMQ  + NVrMMT
            IHModO = IHMod  + ((NVrCQM+1)*NVrCQM)/2 + 1
            IHReal = IHModO + ((NVrCQM+1)*NVrCQM)/2 + 1
            ICNewQ = IHReal + ((NVrCQM+1)*NVrCQM)/2 + 1
            IFTp   = ICNewQ + 3*NAtoms
            IDT    = IFTp   + 3*NAtoms
            ICnt   = IDT    + NRF
            ICOrig = ICnt   + InToWP(NRF)
            IT     = ICOrig + 3*NAtoms
            ICOM   = IT     + 90*NRF
            IFFQ   = ICOM   + 3*NRF
C           Allocation specifically for coupled QM/Q step.
            IXXQ   = IFFQ   + 3 * (NVarQ + NVarZM )
            IFTmpQ = IXXQ   + 3 * NVarQ
            IXQT   = IFTmpQ + NVarQ + NVarZM
            IEndX  = IXQT   + NVarQ
            Call TstCor(IEndX,MDV,'Optmz-GrdOpt')
            MDVX   = MDV - IEndX - 1
            If(Mod(IRFO,100).eq.0.and.ConvF.lt.OneM4) IRFO = IRFO + 8
C
C           Call GrdOpt to take an optimization step.
C
            If(NewCrd) then
              IZCrd = ICrDat
            else
              IZCrd = IIZRed
              endIf
            Call GrdOpt(IPrint,UpTrus,Sphere,DoLin,MethQ,FCart,TstEig,
     $        SetStp,UpdRWF,ToAng,IDimSc,SteepM,V(IForce),V(IFrcns),
     $        V(IX),V(IXName),V(IFC1),V(IFC2),V(IESave),V(IFS),V(IXX),
     $        V(IFF),V(IIC),V(IITU),V(IXQuad),V(IXNew),V(IFTemp),V(IA),
     $        V(IB),V(ISvEig),V(IScr),V(IFCT),V(IAPrt),V(IIntVc),NAtoms,
     $        NAtNew,V(IIAn),V(IAtChg),V(ICA),DoClnG,TRotT,MthUpd,Wght,
     $        IQST,DoRed,NTrRot,NTRedS(2),NTRedS(3),NTRedS(4),V(IZCrd),
     $        V(IValue),V(ICnst1),V(ICnst2),V(ICNew),V(IIBMat),V(IBMat),
     $        V(IScGDI),Failed,IOSLEq,ConvEq,LenSB,LenCB,IOldRC,IOldRB,
     $        AlUnit,FRMS,NPDir,OptMic,IRwF1,V(IMicO1),C,TrVec,V(IDis),
     $        AtChMM,IAtTyp,NumTpS,IAtTpS,MaxLEF,EField,IRwMM,NVarMM,
     $        V(IIDBM),V(IDBM),V(IFTmp),ISCRF,IDir,IGuess,IRFO,V(IEQMS),
     $        V(IXXMM),V(IFFMM),V(IXMMT),V(IHMod),V(IHModO),V(ICNewQ),
     $        IMMRFO,IMMMOD,IMMCRS,V(IFTp),V(IDT),V(ICnt),V(ICOrig),
     $        V(IT),V(ICOM),NRF,RotQM,JScan,V(IIZRdO),IHVAcc,DoMicO,
     $        SkpMic,IGetTS,IAprBG,IGISav,IGIOld,IRwSM,IPFlag,AllowP,
     $        FMM,FMFlag,FMFlg1,NFxFlg,IHMeth,NVarQ,V(IFFQ),V(IXXQ),
     $        V(IFTmpQ),V(IXQT),IntHes,IRwF2,DoEle,V(IFFMMQ),V(IFMMQ),
     $        V(IHReal),NTRedS,IDir0,ONLnOK,IScMMF,DidMix,IGDRFO,FixedC,
     $        ConvP,NGoDwn,IRwCrd,SavLMO,IEStpM,DXLimC,NFroz,V(IFroze),
     $        IGEVer,SwitMx,SimSw,RawDIS,MxRFOD,CThrHd,CDAdj,V(IEndX),
     $        MDVX)
            endIf
          If(StOnly) then
            Call ILSW(1,25,0)
            Write(IOut,1030)
            Exit = .True.
            endIf
          endIf
        endIf
      If(Exit.and.ITqry(IOScan).gt.0) then
        Call TstCor(IScr+LRwGen-1,MDV,'CkScan')
C       Handle special case of convergence at 1st pt.
        If(UseOT.and.NP.eq.1.and.JScan(1).eq.1) then
          IFxOT = IEnd + 1
          IV2 = IFxOT + NAt3
          Call TstCor(IV2,MDV,'Optmz-RdWrOT-Save1a')
          Call FileIO(2,-IRwF1,NAt3,V(IFxOT),0)
          Call RotF1(NAtoms,TRotT,V(IFxOT))
          Stats(1) = Energy
          Stats(2) = Zero
          IDone = IStep
          Junk1 = 0
          Junk2 = 0
          MaxJb1 = Max(JScan(2),1)
          IJob1 = Max(JScan(1),1)
          Call RdWrOT(MaxStp,MaxJb1,NStat,NAt3,Junk1,Junk2,IStep,
     $      IJob1,IDone,.True.,Stats,.True.,V(ICStOr),.False.,XX,.True.,
     $      V(IFxOT),IOut,0,2,V(IV2),MDV-IV2+1)
          endIf
        HavDEP = .False.
        If(DoRed) then
          Call AMove(NVar,V(IX),V(IScr))
          Call ASub(NVar,V(ICnst1),V(IX),V(IXNew))
          Call AMask(NVar,0,1000,V(IIntVc),NOut,V(IXNew),
     $      V(IXNew))
          Call AAdd(NVar,V(IX),V(IXNew),V(IXNew))
          Call CkScan(IOScan,IOGen,IOGenA,LRwGen,NVar,Exit,IStep,NP,
     $      V(IIntVc),V(IFPVec),V(ICnst1),Energy,V(IX),V(IScr+NVar),
     $      JScan)
C
C         Note that Routine CkScan may have changed Exit.
          If(.not.Exit) then
            Call ASub(NVar,V(IX),V(IXNew),V(IXNew))
            Call AMask(NVar,0,1000,V(IIntVc),NOut,V(IXNew),
     $        V(IXNew))
            Call AAdd(NVar,V(ICnst1),V(IXNew),V(ICnst1))
            Call AMove(NVar,V(IScr),V(IXNew))
            Call ASub(NVar,V(IX),V(IXNew),V(IXNew))
            If(NewCrd) then
              Call RedCar(IOut,IPrt1,IOSLE1,NVar,NTRedS(2),NTRedS(3),
     $          NTReds(4),V(ICrDat),V(IX),NAtNew,V(ICA),V(IXNew),
     $          V(ICNew),V(IBMat),V(IIBMat),V(IB),V(IFCT),V(IScr),
     $          V(ICnst1),V(IIntVc),.True.,IFail,0,V(IEnd+1),MDV-IEnd)
            else
              Call RedCar(IOut,IPrt1,IOSLE1,NVar,NTRedS(2),NTRedS(3),
     $          NTReds(4),V(IIZRed),V(IX),NAtNew,V(ICA),V(IXNew),
     $          V(ICNew),V(IBMat),V(IIBMat),V(IB),V(IFCT),V(IScr),
     $          V(ICnst1),V(IIntVc),.True.,IFail,0,V(IEnd+1),MDV-IEnd)
              endIf
            If(IFail.eq.1) Call GauErr('RedCar failed.')
            Call AMove(3*NAtNew,V(ICNew),V(ICA))
            Call FileIO(1,-IOZRed,2*NVar,V(ICnst1),InToWP(4)+
     $        InToWP(4*NVar)+IntoWP(NVar)+2*NVar)
            If(ITqry(IMMCRS).gt.0) then
              Call FileIO(2,-IMMCRS,6*NRF,V(IXMMT+NVarMM),2)
              Call FileIO(2,IMMCRS,3*NAtoms,V(ICOrig),0)
              Call AClear(NVarT,V(IXNew))
              Call UpCOrg(V(ICA),V(ICNew),V(ICOrig),NAtoms,NAtNew,
     $          V(IMicOp),V(IXMMT),V(ICOM),V(IDT),V(ICnt),V(IT),NRF,
     $          NVar,NVarMM,V(IXNew),IMMCRS,RotQM,.True.)
              endIf
            endIf
        else
          Call CkScan(IOScan,IOGen,IOGenA,LRwGen,NVar,Exit,IStep,NP,
     $      V(IIntVc),V(IFPVec),Values,Energy,V(IX),V(IScr),JScan)
          endIf
        endIf
C
C     Update the r/w files and exit the link, unless prohibited
C     by option or the run is done, in which case just leave.
C     use some care in determining jump action as ChainX alters IOp.
C
      Write(IOut,1000)
      Write(IOut,1070)
      If(Exit.and.ITqry(IOScan).gt.0) Call PtScan(IOut,IOScan,NVar,
     $  ToAng,MaxNZ,NZ,LBl,LAlpha,LBeta,IZ,V(IXName),V(IFS),
     $  DoRed,NTReds(2),JScan)
      If(UseOT.and.NP.eq.1.and.JScan(1).le.1) then
        IFxOT = IEnd + 1
        IV2 = IFxOT + NAt3
        Call TstCor(IV2,MDV,'Optmz-RdWrOT-Save2')
        Call FileIO(2,-IRwF1,NAt3,V(IFxOT),0)
        Call RotF1(NAtoms,TRotT,V(IFxOT))
        Stats(1) = Energy
        Stats(2) = Zero
        IDone = IStep
        Junk1 = 0
        Junk2 = 0
        MaxJb1 = Max(JScan(2),1)
        IJob1 = Max(JScan(1),1)
        Call RdWrOT(MaxStp,MaxJb1,NStat,NAt3,Junk1,Junk2,IStep,
     $    IJob1,IDone,.True.,Stats,.True.,V(ICStOr),.False.,XX,.True.,
     $    V(IFxOT),IOut,0,2,V(IV2),MDV-IV2+1)
        endIf
      If(Exit.and.DoRed.and.NPath.gt.0) then
        Call PtPath(IOut,NVar,ToAng,NPath,V(IXP),V(IEsvP),V(IXName),
     $    NTRedS(2),.True.)
        If(SetStp.and.Do601) then
          LOpatN = (2+2+3*Nat3+(Nat3*(Nat3+1))/2)
          LOPath = 2 + NPath*LOpatN
          WoPath = GFloat(LOPath)
          Call Fileio(1,-IOPath,1,WOPath,1)
          endIf
        endIf
      If(SetStp) then
        If(.not.AlUnit) Call AMove(NVTT,V(IFC1),V(IFrcns))
        Call FileIO(1,-IGrdnt,LGrdnt,Energy,0)
        Call FileIO(1,IGrdnt,NVarZM,V(IForce),0)
        If(.not.AlUnit) Call FileIO(1,IGrdnt,NVTT,V(IFrcns),0)
        endIf
      If(.not.Exit) then
        Call FileIO(1,-IOptGr,LOptGr,ConvF,0)
        Call FileIO(1,IOptGr,LOptG2,V(IIC),0)
        endIf
      If(.not.First.and.FCart) Call AMove(NVarZM,V(IX),V(ICA))
      If(OptMic) Call ChkMic(4,NAtoms,NAtNew,NAtMic,V(IMicO1),
     $  V(IMicIn),V(ICA),NTRed,V(IIZRed),V(IIAn),V(IAtChg),V(IScr),
     $  V(IScr))
      If(DoMicO.and..not.First.and..not.Failed.and.SetStp.and..not.Exit.
     $  and..not.SkpMic) then
        Call TStamp(1,'Before Micro')
C       The code in MMOpt seems to have some parallelism but not
C       enough to be efficient, so for the moment we turn off
C       parallelism here.
        NPSave = NProc(1)
        IFX = NProc(-1)
        If(IPrint.gt.0) Call CorPrS(IOut,'C before microiterations:',0,
     $    NAtoms,V(IIAn),.True.,IAtTyp,V(ICA),ToAng)
        ConvMM = Max(Min(ConvF/Ten,FRMS/Ten),ConvF/GFloat(30))
        ConvMM = ConvF/GFloat(IScMMC)
        IParam = IEnd + 1
        LenPar = ITqry(IRwMM)
        IFX    = IParam + LenPar
        MethMM = MethQ
        If(IDirMM.ge.2.and.IDirMM.le.5) then
          MicTyp = 3
          LenFFX = 3*NAtoms
        else
          MicTyp = 1
          LenFFX = 0
          endIf
        IFFX   = IFX + 3*NAtoms
        ICBest = IFFX + LenFFX
        ICOrig = ICBest + 3*NAtoms
        ICRig  = ICOrig + 3*NAtoms
        LCRig  = IArMax(.False.,NAtoms,V(IMicO1))*6
        IV1    = ICRig + LCRig
        MMMaxS = MaxMic
        If(MMMaxS.eq.0) MMMaxS = Max(3*NAtoms+7*Min(NAtoms,150),MinMMC)
        Call TstCor(IV1,MDV,'Optmz-MMMicr')
        MDV1 = MDV - IV1 + 1
        If(LenPar.gt.0) Call FileIO(2,-IRwMM,LenPar,V(IParam),0)
        NRFX = 0
        If(ITqry(IMMCRS).gt.0) then
          Call FileIO(2,-IMMCRS,1,V(ICOrig),0)
          Call BNUPak(V(ICOrig),NRFX,IDum)
          endIf
        MaxMMR = 10
        EBest = BigE
        Do 300 NMMRes = 1, MaxMMR
          Call AMove(3*NAtoms,V(ICA),V(ICB))
          Call AMove(3*NAtoms,V(ICA),V(ICOrig))
          Call AClear(LCRig,V(ICRig))
          If(NRFX.gt.0.and.NMMRes.eq.1) then
            Call FileIO(2,-IMMCRS,6*NRF,V(ICRig),2)
            Call FileIO(2,IMMCRS,3*NAtoms,V(ICOrig),0)
            endIf
          IPr1 = IPrint
C         These "ISubst"s turn the QM region into a regular rigid block,
C         so that the frozen part of the system retains its coordinates.
          If(RotQM) then
            Call ISubsI(NAtoms,V(IMicO1),-1,NRF+1)
            Call ISubsI(NAtoms,V(IMicO1),1,-1)
            endIf
          Call MMMicr(IOut,IPr1,.False.,ISCRF,MicTyp,MMMaxS,ConvMM,
     $      NAtoms,IAtTyp,AtChMM,NumTpS,IAtTpS,V(IMicO1),.False.,V(ICB),
     $      MaxLEF,EField,V(IParam),LenPar,MMStpd,Convgd,MMCyc,EnerMM,
     $      V(ICA),V(IFX),XX,XX,XX,XX,XX,XX,XX,ITRScl,IAn,V(IFFX),IRWMM,
     $      IGssMM,IRFOMM,IPrjct,0,0,JJ,XX,XX,XX,XX,XX,XX,XX,MethMM,
     $      IDirMM,JJ,XX,UseCA,XX,IMMRFO,V(ICOrig),V(ICRig),0,IHVAcc,
     $      .False.,XX,JJ,0,0,0,0,.False.,ICharg,Multip,SmlStp,0,0,0,
     $      IPFlag,AllowP,FMM,FMFlag,FMFlg1,NFxFlg,IHMeth,NAtoms,LJJ,JJ,
     $      XX,JJ,JJ,JJ,NAtoms,.False.,0,QMVar,.False.,XX,.False.,
     $      .True.,MMRest,IRwSMM,V(IV1),MDV1)
          If(RotQM) then
            Call ISubsI(NAtoms,V(IMicO1),-1,1)
            Call ISubsI(NAtoms,V(IMicO1),NRF+1,-1)
            endIf
          If(EnerMM.lt.EBest) then
            EBest = EnerMM
            Write(IOut,4080) EBest
            Call AMove(NAt3,V(ICA),V(ICBest))
            endIf
          If(NRFX.gt.0) Call FileIO(1,-IMMCRS,6*NRF,V(ICRig),2)
          If(IPrint.gt.0) Call CorPrS(IOut,'C after micro:',NMMRes,
     $      NAtoms,V(IIAn),.True.,IAtTyp,V(ICA),ToAng)
          If(MMRest.and.NMMRes.eq.1) then
          else if(MMRest.and.NMMRes.lt.MaxMMR) then
            Call TstCor(NAt3,MDV1,'Optmz-CA')
            Call ACSASB(NAt3,V(ICBest),Pt9,V(ICA),Pt1,V(IV1))
            Call AMove(NAt3,V(IV1),V(ICA))
            Call AClear(LCRig,V(ICRig))
            Call AMove(NAt3,V(ICBest),V(ICOrig))
            Call AMove(NAt3,V(ICBest),V(ICB))
          else
            Call AMove(NAt3,V(ICBest),V(ICA))
            Goto 310
            endIf
  300     Continue
        Convgd = .False.
        Write(IOut,4090)
  310   Exit = Exit.and.Convgd
        UpdMol = MMStpd.and.Convgd
        IFX = NProc(-NPSave)
        endIf
      DoClnS = DoClnS.and..not.DoClnG
      RmTrRo = .False.
      If(RmTrRo.or.DoClnS) then
        Call ClnCor(IOut,IPrint,RmTrRo,DoClnS,NAtoms,V(IIAn),V(IAtChg),
     $    CSave,0,V(IScr),.True.,TRotT,V(ICA),V(IScr))
        If(IPrint.gt.0) Call CorPrS(IOut,'C after MM cleanup:',0,NAtoms,
     $    V(IIAn),.True.,IAtTyp,V(ICA),ToAng)
        endIf
      If(.not.First) Call AMove(NAt3,V(ICA),C)
      If(Exit.or..not.UpdRWF) then
        Call AMove(3*NAtoms,CSave,C)
        If(SetStp) then
          Call RotF1(NAtoms,TRot,C)
          Call Transl(-1,3,NAtoms,TrVec,C)
          endIf
        IChain = 1
      else
        If(.not.First.and..not.FCart.and..not.DoRed) then
          Call UdZSub(V(IX),NVarZM)
          endIf
        If(IExit.ge.0) then
          IChain = 0
        else
          IChain = 1
          endIf
        If(UpdMol) Call UpdMMO(IOut,IPrint,MaxAtm,NAtoms,V(IMicO1),
     $    TrVec,TRot,C,V(ICA))
        endIf
      If(SetStp) then
        If(IPrint.ge.1.or.Exit) Call PrtCDs(IOut,IPrint,IRwSvC,NAtoms,
     $    IAn,IAtTyp,ToAng,TrVec,TRot,C,V)
        endIf
      NAtoms = NAtoms - NPDir
      If(PCM(ISCRF)) Call PCMIOp('save',' ',IJunk,RJunk)
      If(IChain.eq.0) Call ILSW(1,24,0)
      If(IFirst.eq.1) then
        Call ILSW(1,17,-1)
      else if(NP.ge.0) then
        If(UpdRWF.and.(IChain.eq.0.or.Failed)) then
          Call ILSW(1,17,0)
        else
          Call ILSW(1,17,2)
          endIf
      else
        Call ILSW(1,17,1)
        endIf
      L = ITqry(IGISav)
      If(L.gt.0.and.IGISav.ne.IGIOld) then
        Call ConDD2(IGIOld,L)
        Call Trsfr(L,IGISav,IGIOld,V,MDV)
        endIf
      Call ChainX(IChain)
      Return
      End
*Deck ChkAng
      SubRoutine ChkAng(IOut,NPath,NVar,NTBond,NTAng,IZRed,XP)
      Implicit Real*8(a-h,o-z)
C
C     for path relaxation
C     Checks if angles go past 180 degrees along the path.
C
      Logical OK
      Dimension XP(Nvar,NPath),IZred(4,Nvar)
      Save Zero
      Data Zero /0.0d0/
 1000 Format(' Possible problem with angle ',I3)
C
      Do 100 I = (NTBond+1), (NTBond+NTAng)
        IChang = 0
        OK = .True.
        Do 110 N = 2, (NPath-1)
          DA1 = XP(I,N)-XP(I,N-1)
          DA2 = XP(I,N+1)-XP(I,N)
          If(DA1*DA2.lt.Zero) IChang = IChang + 1
  110     Continue
        If(IChang.gt.0) then
          ia = IZred(1,I)
          ja = IZred(2,I)
          ka = IZred(3,I)
          Do 120 II = (NTBond+NTAng+1), NVar
            ib = IZred(1,II)
            jb = IZred(2,II)
            kb = IZred(3,II)
            lb = IZred(4,II)
            If((jb.eq.ia.and.kb.eq.ja.and.lb.eq.ka).or.
     $         (ib.eq.ka.and.jb.eq.ja.and.kb.eq.ia).or.
     $         (ib.eq.ia.and.jb.eq.ja.and.kb.eq.ka).or.
     $         (jb.eq.ka.and.kb.eq.ja.and.lb.eq.ia)) then
              IChang = 0
              Do 130 N = 2, NPath
                If(XP(II,N)*XP(II,N-1).lt.Zero) IChang = IChang + 1
  130           Continue
              If(IChang.gt.0) OK=.False.
              endIf
  120       Continue
          endIf
        If(.not.OK) Write(IOut,1000) (I-NTBond)
  100   Continue
      Return
      End
*Deck ChkDIS
      Subroutine ChkDIS(IOut,Prnt,Prnt1,NP,CEn,CUpper,CLower,CLarge,OK)
      Implicit Real*8(A-H,O-Z)
C
C     Routine to adjust DIIS coefficient using a uniform scaling method
C
C     NP       Number of point orginally in DIIS matrix
C     CLarge   Threshold when DIIS is considered ill-conditioned
C     CUpper   Ajusted upper-bound of DIIS coefficient
C     CLower   Ajusted lower-bound of DIIS coefficient
C
      Logical Prnt,Prnt1,OK,DidPrt
      Dimension CEn(NP)
      Save One, Small
      Data One/1.D0/, Small/1.d-5/
 1000 Format(' Old DIIS coefficients:',5F13.5)
 1010 Format(' New DIIS coefficients:',5F13.5)
C
      OK = .True.
      DidPrt = .False.
      CMax = CLower
      CMin = CUpper
      IMax = 0
      IMin = 0
      Do 100 i = 1, NP
        If(CEn(i).lt.CMin) then
          CMin = CEn(i)
          IMin = i
        else if(CEn(i).gt.CMax) then
          CMax = CEn(i)
          IMax = i
          endIf
  100   Continue
C
C     DIIS coefficients are too large
C
      If(Abs(CMin).gt.CLarge.or.Abs(CMax).gt.CLarge) then
        If(Prnt1) Write(IOut,1000) (CEn(I),I=1,NP)
        Call AClear(NP,CEn)
        CEn(1) = One
        If(Prnt) Write(IOut,1010) (CEn(I),I=1,NP)
        OK = .False.
        Return
        endIf
C
C     Adjust coefficients smaller than the lower bound
C
      If(CMin.lt.CLower.and.Abs(CMin-CLower).gt.Small) then
        If(Prnt1) Write(IOut,1000) (CEn(I),I=1,NP)
        DidPrt = .True.
        flmC = (CLower - CMin)/(One - CMin)
        If((flmC.gt.Small).and.(Abs(One-flmC).gt.Small)) then
          CEn(IMin) = CLower/(One - flmC)
          Do 200 i = 1, NP
  200       CEn(i) = CEn(i)*(One - flmC)
          endIf
        If(Prnt) Write(IOut,1010) (CEn(I),I=1,NP)
        endIf
C
C     Adjust coefficients greater than the upper bound
C
      CMax = CLower
      Do 300 i = 1, NP
        If(CEn(i).gt.CMax) then
          CMax = CEn(i)
          IMax = i
          endIf
  300   Continue
      If((CMax.gt.CUpper).and.(Abs(CMax-CUpper).gt.Small)) then
        If(Prnt1.and..not.DidPrt) Write(IOut,1000) (CEn(I),I=1,NP)
        flmC = (CUpper - CMax)/(One - CMax)
        If((flmC.gt.Small).and.(Abs(One-flmC).gt.Small)) then
          CEn(IMax) = CUpper/(One - flmC)
          Do 400 i = 1, NP
  400       CEn(i) = CEn(i)*(One - flmC)
          endIf
        If(Prnt) Write(IOut,1010) (CEn(I),I=1,NP)
        endIf
      Return
      End
*Deck ChkDSM
      Subroutine ChkDSM(IOut,Print,DISA,DISC,DISCM,Map,NDim,NP,CThrHd,
     $  OK,DidRem)
      Implicit Real*8(A-H,O-Z)
C
C     DISA     DIIS A Matrix (NDim*NDim)
C     DISC     DIIS coefficients (NDim)
C              Note: The first coefficient corresponds to lamda
C     DISCM    Array used to indicate the deleted points
C              (not including the lamda)
C
C     NDim     Dimension of the problem (Number of points + 1)
C              New dimension of the problem on output
C     NP       Number of point orginally in DIIS matrix
C     CThrHd   Threshold of the coefficients.
C
      Logical Print,OK,DidRem
      Dimension DISA(*),DISC(*),DISCM(*),Map(*)
      Save Zero, Half
      Data Zero/0.D0/, Half/0.5D0/
 1000 Format(' Point number',I4,' has GDIIS coefficient',1PD10.2,
     $  ' and has been removed')
C1010 Format(' Point number',I4,' has GDIIS coefficient',1PD10.2,
C    $  ' -- DIIS solution rejected.')
C
      DidRem = .False.
      J = 0
      Do 10 I = 1, NP
        If(DISCM(I).gt.Half) then
          J = J + 1
          Map(J) = I
          endIf
   10   Continue
      If((J+1).ne.NDim) Call GauErr('Consistency failure #1 in ChkDSM.')
C
C     Check the absolute value of the coefficients
C     Start from the second point
C     Delete one point at a time
C
      JSt = NDim - 1
      Do 20 J = JSt, 1, -1
        I = Map(J)
        If(Abs(DISC(I)).gt.CThrHd) then
          If(I.le.1.or.J.le.1) then
            If(I.ne.J) Call GauErr('Consistency failure #2 in ChkDSM.')
C           OK = .False.
C           If(Print) Write(IOut,1010) I, DISC(I)
            Goto 40
          else
            Call DelBas(NDim,NDim,J+1,J+1,DISA,DISA)
            Call Pcck(0,KK,DISA,NDim,NDim,NDim-1,NDim-1)
            DISCM(I) = Zero
            DidRem = .True.
            NDim = NDim - 1
            If(Print) Write(IOut,1000) I, DISC(I)
            endIf
          OK = .False.
          Return
          endIf
   20   Continue
C
C     Pass coefficients check, construct the final set of coefficients with
C     deleted points being masked.
C
   40 Call AClear(NP,DISCM)
      Do 30 J = 1, (NDim-1)
        I = Map(J)
   30   DISCM(I) = DISC(J)
      Return
      End
*Deck ChkMic
      Subroutine ChkMic(Option,NAtoms,NAtNew,NAtMic,MicOpt,MicInd,C,
     $  NTRed,IZRed,IAN,AtmChg,IScr,Scr)
      Implicit Real*8 (A-H,O-Z)
C
C     Chack the atoms for MM microiteration and prepare the index
C     vector for it and reorder the atoms.
C     Option = 1   : initialize index vector only.
C     Option = 2   : reorder the arrays using existing MicInd.
C     Option = 3   : reorder the C array only using existing MicInd.
C     Option = 4,0 : recover the original order.
C
      Integer Option
      Dimension MicOpt(NAtoms), C(3,NAtoms), Scr(3,NAtoms),
     $  MicInd(2,NAtoms), IScr(NAtoms), IAN(NAtoms), IZRed(4,NTRed),
     $  AtmChg(NAtoms)
C
      If(Option.eq.1) then
        NAtNew = 0
        Do 10 I = 1, NAtoms
          If(MicOpt(I).lt.0) then
            NAtNew = NAtNew + 1
            MicInd(1,NAtNew) = I
            MicInd(2,I) = NAtNew
            endIf
   10     Continue
        NA0 = NAtNew
        NAtMic = 0
        Do 20 I = 1, NAtoms
          If(MicOpt(I).ge.0) then
            NAtMic = NAtMic + 1
            MicInd(1,NA0+NAtMic) = I
            MicInd(2,I) = NA0+NAtMic
            endIf
   20     Continue
        endIf
      If((Option.eq.2.or.Option.eq.3).and.NAtMic.ne.0) then
        Do 30 I = 1, NAtoms
          Scr(1,I) = C(1,MicInd(1,I))
          Scr(2,I) = C(2,MicInd(1,I))
          Scr(3,I) = C(3,MicInd(1,I))
   30     Continue
        Call AMove(3*NAtoms,Scr,C)
        endIf
      If(Option.eq.2.and.NAtMic.ne.0) then
        Do 40 I = 1, NTRed
          Do 40 J = 1, 4
            II = IZRed(J,I)
            If(II.gt.0) then
              IZRed(J,I) = MicInd(2,II)
            else if(II.lt.-10) then
              I1 = IAbs(II)/10
              I2 = Mod(IAbs(II),10)
              IZRed(J,I) = -(10*MicInd(2,I1)+I2)
              endIf
   40       Continue
        Do 50 I = 1, NAtoms
   50     IScr(I) = IAN(MicInd(1,I))
        Call IMove(NAtoms,IScr,IAN)
        Do 60 I = 1, NAtoms
   60     Scr(1,I) = AtmChg(MicInd(1,I))
        Do 70 I = 1, NAtoms
   70     AtmChg(I) = Scr(1,I)
        endIf
      If(Option.eq.4.and.NAtMic.ne.0) then
        Do 80 I = 1, NTRed
          Do 80 J = 1, 4
            II = IZRed(J,I)
            If(II.gt.0) then
              IZRed(J,I) = MicInd(1,II)
            else if(II.lt.-10) then
              I1 = IAbs(II)/10
              I2 = Mod(IAbs(II),10)
              IZRed(J,I) = -(10*MicInd(1,I1)+I2)
              endIf
   80       Continue
        Do 90 I = 1, NAtoms
          Scr(1,I) = C(1,MicInd(2,I))
          Scr(2,I) = C(2,MicInd(2,I))
   90     Scr(3,I) = C(3,MicInd(2,I))
        Call AMove(3*NAtoms,Scr,C)
        Do 100 I = 1, NAtoms
  100     IScr(I) = IAN(MicInd(2,I))
        Call IMove(NAtoms,IScr,IAN)
        Do 110 I = 1, NAtoms
  110     Scr(1,I) = AtmChg(MicInd(2,I))
        Do 120 I = 1, NAtoms
  120     AtmChg(I) = Scr(1,I)
        endIf
      Return
      End
*Deck ChkPth
      Subroutine ChkPth(IOut,NPath,IStep,NConv,IConv,NPStar,
     $DoMid,DoIrc,DoReac,DoProd,ESvP,Stot,PathL)
      Implicit Real*8(A-H,O-Z)
C
C For path relaxation.
C Decide which point will optimize to TS
C
      LogiCal DoMId,DoIrc,DoReac,DoProd
      Dimension IConv(NPath),ESvP(NPath),STot(2),PathL(NPath)
      Save Zero
      Data Zero /0.0d0/
 1000 Format(' Use of a three structure input recommended')
 1010 Format(1x,' Monotonic reaction path')
 1020 Format(1x,' Possibility for multiple intermediates')
 3010 Format(1x,I2,' Points out of ',I2,' obtained')
 3020 Format(1x,' Energy maximum found at ',I2)
 3030 Format(1x,I2,'-th Point elected by LST to optimize to TS')
 3040 Format(1x,I2,'-th point chosen by user to optimize to TS')
 4000 Format(1x,'Path length is ',F10.5,' times the starting length')
C
C     Do an inventory of calculation already performed
C
      NConv=0
      Do 119 NN=1,NPath
      NConv=NConv+IConv(NN)
 119  Continue
      Write(IOut,3010)NConv,NPATH
C
C     Get maximum along the path if all calculations have been
C     performed
C
      If(NConv.eq.NPath) Then
      NPTS=0
      Do 121 NPt=1,NPath
      JJ=1
      Do 120 NNPt=1,NPath
 120  If(ESvP(NNPt).gt.ESvP(NPt))JJ=JJ+1
      If(JJ.eq.1)NPTS=NPt
 121  Continue
      IChange=0
      Do 122 Npt=2,NPath-1
      dE1=ESvP(Npt)-ESvP(NPt-1)
      dE2=ESvP(NPt+1)-ESvP(Npt)
      If(de1*de2.lt.Zero)IChange=IChange+1
 122  Continue
      If(IChange.eq.0)then
      Write(IOut,1010)
      Else if(IChange.gt.1)then
      Write(IOut,1020)
      Endif
CPya
C the highest point will be optimizing to the TS unless domid is
C true, in which case the middle point will optimize to TS.
C
      If(IStep.eq.1)then
      If(.not.DoMid)then
      NPStar=NPTS
      If(DoIrc.and.NPTS.ne.1)Then
      Write(IOut,1000)
      NPStar=1
      Endif
      If(NPTS.eq.1.and.DoReac.and.(.not.DoIrc))
     $Call GauErr(' Point 1 is the highest energy point')
      If(NPTS.eq.Npath.and.DoProd)
     $Call GauErr(' Last point is the highest energy point')
      Write(IOut,3020) NPTS
      Write(Iout,3030) NPStar
      Else
      Write(IOut,3020) NPTS
      NPStar=(NPath+1)/2
      Write(Iout,3040) NPStar
      Endif
      Endif
C
C Chek Path Length
      If(IStep.ge.2)then
      Sjunk=Zero
      Do 130 NPt=1,NPath
      Sjunk=Sjunk+Abs(PathL(Npt))
 130  Continue
      STot(2)=Sjunk
      R=Stot(2)/Stot(1)
      Write(IOut,4000)R
      Endif
      Endif
      Return
      End
*Deck ChoNxt
      Subroutine ChoNxt(IOut,NPath,NIStep,NCurPt,IFinish,IConv,NNext)
      Implicit Real*8(a-h,o-z)
C
C  Choose next point for path relaxation
C
      Dimension IConv(NPath)
 1010 Format(' Skip gradient calculation for point ',I3,
     $' -- Point converged')
C
      If(NIStep.eq.0)NNext=NCurPt-1
      If(NIStep.eq.1)NNext=NCurPt+1
      NIStp2=NIStep
      If(NIStep.eq.1.and.NCurPt.eq.NPath)then
      NNext=NPath
      NIStp2=0
      Endif
      If(NIStep.eq.0.and.NCurPt.eq.1)then
      NNext=1
      NIStp2=1
      Endif
      If(IFinish.eq.NPath)Return
 100  If(IConv(NNext).eq.1)then
      Write(IOut,1010)NNEXT
      NNext2=NNext
      If(NIStp2.eq.0)NNext=NNext2-1
      If(NIStp2.eq.1)NNext=NNext2+1
      If(NIStp2.eq.1.and.NNext2.eq.NPath)then
      NNext=NPath
      NIStp2=0
      Endif
      If(NIStp2.eq.0.and.NNext2.eq.1)then
      NNext=1
      NIStp2=1
      Endif
      Endif
      If(IConv(NNext).eq.1)Goto 100
      Return
      End
*Deck CkHesU
      Logical Function CkHesU(IOut,IPrint,NP,ITU)
      Implicit Real*8(A-H,O-Z)
C
C     Return whether the Hessian appears to be trustworthy, based on
C     whether the trust radius test has been good or bad at recent
C     points.
C
      Parameter (MinPrt=0,NTest1=3,NTest2=10)
      Dimension ITU(*)
      Logical OK
 1010 Format(' ITU=',20I3)
C
      If(IPrint.ge.MinPrt) Write(IOut,1010) (ITU(I),I=1,NP)
      OK = .True.
      Do 10 I = 1, Min(NTest1,NP)
   10   OK = OK.and.ITU(I).ge.0
      ISum = 0
      Do 20 I = 1, Min(NTest2,NP)
   20   ISum = ISum + ITU(I)
      OK = OK.and.ISum.ge.0
      CkHesU = OK
      Return
      End
*Deck CkScan
      Subroutine CkScan(IOScan,IOGen,IOGenA,LRwGen,NVar,Exit,IStep,NP,
     $  IntVec,FPVec,Values,Energy,CurVal,Gen,JScan)
      Implicit Real*8(A-H,O-Z)
C
C     Check for relaxed scan at the end of each optimization.
C
      Dimension IntVec(NVar), FPVec(NVar), Values(NVar), CurVal(NVar),
     $  Gen(LRwGen), JScan(4)
      Logical Exit
C
      Call FileIO(2,-IOScan,InToWP(2),JScan,0)
      IPos = InToWP(2) + JScan(1) - 1
      Call FileIO(1,-IOScan,1,Energy,IPos)
      IPos = InToWP(2) + JScan(2) + (JScan(1)-1)*NVar
      Call FileIO(1,-IOScan,NVar,CurVal,IPos)
      KScan1 = JScan(1) - 1
      KScan = JScan(1)
      Call FileIO(2,-IOGen,LRwGen,Gen,0)
      Call FileIO(1,-IOGenA,LRwGen,Gen,LRwGen*KScan1)
      Do 10 I = 1, NVar
        If(IAbs(IntVec(I)).gt.1000) then
          LScan = IAbs(IntVec(i))-1000
          j = (Mod(KScan,LScan)-Mod(KScan1,LScan)+LScan-1)/LScan
          If(Mod(KScan/LScan,2).eq.1) J = -J
          CurVal(I) = Values(i) + GFloat(j)*FPVec(i)
          KScan1 = KScan1/LScan
          KScan = KScan/LScan
          endIf
   10   Continue
      JScan(1) = JScan(1) + 1
      If(JScan(1).le.JScan(2)) then
        IStep = 0
        NP = 0
        Exit = .False.
        endIf
      Call FileIO(1,-IOScan,InToWP(2),JScan,0)
      Return
      End
*Deck CnsST4
      Subroutine CnsST4(IQST0,N,IFlags,Proj,Values,Cnstr1,Cnstr2,F,
     $  GrdTmp,Vec,NTBond,NTAng,NTDih,RP,RR,RPR,arclnR,arclnP,IDist,
     $  Iprint,doeuler)
      Implicit Real*8(A-H,O-Z)
C
C For path relaxation. Modified ConsST.
C Get 3-Point QST tangent and redistribute points if needed
C
C
C     DoEuler ... If true T= X-R, if false T=T(QST)
C                 Uses the pojection of XP onto RX instead of XP
C     IDist .... 1/0 for whether to redistribute the X
C     RP .... If Idist=1 RP=|P-X|
C             If IDist=0 RP=Radius of arc of circle linking P and X
C     ArclnP .... Arc length of the arc of circle or parabola passing
C             through X and P with T the tangent @ X
C     Vec(1,1) and Vec(1,2) contain R-X and P-X on exit if IDist=0
C
C
      Logical Doeuler,DoPar
      Dimension IFlags(N), Proj(*), Values(N), Cnstr1(N), Cnstr2(N),
     $  F(N), Vec(N,5), XX(1)
      Real*8 MDCutO
      Common /IO/     In, IOut, IPunch
      Save Zero,Half,One,Two,XX
      Data Zero /0.0d0/,Half/0.5d0/,One/1.0d0/,Two /2.0d0/,XX/0.d0/
 2000 Format(1x,' Tangent with ',I3,' points ; Overlap T.G ',F10.5)
 2001 Format(1x,' Deviation from mediane:|P-X|-|Q-X| =',F10.5,
     $' Ravg=',F10.5,' |PR|=',F10.5)
 2003 Format(1x,' ARC  R->X  = ',F10.5,' X->P =',F10.5)
 2004 Format(1x,' Use 3-point parabola instead of arc of circle')
C
C     For LST/QST(int) in internal var space, IQST0 = 1/2
C
C     Cnstr1 : P
C     Cnstr2 : R
C
      If(IQST0.eq.0) Call Gauerr(' IQST0 = 0 in CnsST4')
      If(IQST0.eq.1) then
C  Simple T=R-P
        Call ASub(N,Cnstr1,Cnstr2,Vec(1,3))
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,Vec(1,3),XX)
      else if(IQST0.eq.2) then
C  T tangent to the arc of circle
        Call ASub(N,Cnstr2,Values,Vec(1,4))
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,Vec(1,4),XX)
        Call MultLV(.True.,1,0,N,N,Proj,Vec(1,4),Vec(1,2))
        RR = SProd(N,Vec(1,2),Vec(1,2))
        Tmp = One/Rr
        Rr=Sqrt(Rr)
        Call AScale(N,tmp,Vec(1,2),Vec(1,2))
        Call ASub(N,Cnstr1,Values,Vec(1,4))
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,Vec(1,4),XX)
C If Doeuler, the use the projection of XP onto RX
        If(DoEuler)then
          Dot = -Abs(RR*RR*SProd(N,Vec(1,4),Vec(1,2)))
          Call Ascale(N,Dot,Vec(1,2),Vec(1,4))
          endIf
        Call MultLV(.True.,1,0,N,N,Proj,Vec(1,4),Vec(1,1))
        RP = SProd(N,Vec(1,1),Vec(1,1))
        Tmp = One/RP
        RP = Sqrt(RP)
        Call AScale(N,Tmp,Vec(1,1),Vec(1,1))
        Call ASub(N,Vec(1,1),Vec(1,2),Vec(1,3))
        endIf
      If(IQST0.eq.1.or.IQST0.eq.2) then
        Call MultLV(.True.,1,0,N,N,Proj,Vec(1,3),Vec(1,4))
        Call AMask(N,-100000,-1,IFlags,NOut,Vec(1,4),Vec(1,4))
        endIf
C
C     Normalize constraint in space of active variables
C
      Tmp = One/Sqrt(SProd(N,Vec(1,4),Vec(1,4)))
      Call AScale(N,Tmp,Vec(1,4),Vec(1,4))
C
C     Normalized tq(N) in Vec(i,4)
C
      DoPar = .False.
      TestR = -SProd(N,Vec(1,2),Vec(1,4))/
     $  Sqrt(SProd(N,Vec(1,2),Vec(1,2)))
      TestP = SProd(N,Vec(1,1),Vec(1,4))/
     $  Sqrt(SProd(N,Vec(1,1),Vec(1,1)))
C
C     If angle <90 do a 3-point parabola instead of arc of circle
C     useless and should be scraped. This case arise if too few points
C     are used.
C
      If(Abs(TestR).gt.One) TestR = Abs(TestR)/testR
      TestR = ACos(TestR)
      If(Abs(testP).gt.One) testP = Abs(testP)/testP
      testP = ACos(testP)
      Test = TestP + Testr
      XPi2 = GFloat(2)*ATan(GFloat(1))
      If(test.gt.XPi2) DoPar = .True.
        If(IQST0.eq.2.and.DoPar) then
        Call ASub(N,Cnstr2,Values,Vec(1,4))
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,Vec(1,4),XX)
        Call MultLV(.True.,1,0,N,N,Proj,Vec(1,4),Vec(1,2))
        RR = SProd(N,Vec(1,2),Vec(1,2))
        Tmp = One/Rr
        Rr = Sqrt(Rr)
        Call AScale(N,tmp,Vec(1,2),Vec(1,2))
        Call ASub(N,Cnstr1,Values,Vec(1,4))
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,Vec(1,4),XX)
        If(DoEuler) then
          Dot = RR*RR*SProd(N,Vec(1,4),Vec(1,2))
          Call Ascale(N,Dot,Vec(1,2),Vec(1,4))
          endIf
        Call MultLV(.True.,1,0,N,N,Proj,Vec(1,4),Vec(1,1))
        RP = SProd(N,Vec(1,1),Vec(1,1))
        Tmp = One/RP
        RP = Sqrt(RP)
        Call AScale(N,Tmp,Vec(1,1),Vec(1,1))
        If(IPrint.gt.1) Write(IOut,2004)
        Call AScale(N,RR*RR,Vec(1,2),Vec(1,3))
        Call AScale(N,RP*RP,Vec(1,1),Vec(1,4))
        Call Asub(N,Vec(1,4),Vec(1,3),Vec(1,5))
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,Vec(1,5),XX)
        DRP = Sqrt(SProd(N,Vec(1,5),Vec(1,5)))
        DRP2 = DRP*Half
        Tmp = One/DRP
        Call AScale(N,Tmp,Vec(1,5),Vec(1,5))
        Call AAdd(N,Vec(1,3),Vec(1,4),Vec(1,4))
        Tmp = -Half
        Call AScale(N,tmp,Vec(1,4),Vec(1,4))
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,Vec(1,4),XX)
        DotX = SProd(N,Vec(1,5),Vec(1,4))
        Call AScale(N,DotX,Vec(1,5),Vec(1,3))
        Call Asub(N,Vec(1,4),Vec(1,3),Vec(1,3))
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,Vec(1,3),XX)
        Tmp = SProd(N,Vec(1,3),Vec(1,3))
        Small = MDCutO(0)
        If(tmp.gt.Small) then
          Tmp = One/Sqrt(tmp)
        else
          Tmp = Zero
          endIf
        Call Ascale(N,tmp,Vec(1,3),Vec(1,3))
        DotY = SProd(N,Vec(1,4),Vec(1,3))
        Gam = DotY/((DotX-DRP2)*(DotX+DRP2))
        GamP = Two*DotX*Gam
        Call AScale(N,GamP,Vec(1,3),Vec(1,3))
        Call AAdd(N,Vec(1,5),Vec(1,3),Vec(1,4))
        Tmp = One/Sqrt(SProd(N,Vec(1,4),Vec(1,4)))
        Call AScale(N,Tmp,Vec(1,4),Vec(1,3))
        B1 = -Two*gam*DRP2
        B2 = GamP
        B3 = -B1
        A1 = GLog(B1+Sqrt(B1*B1+One))
        A2 = GLog(B2+Sqrt(B2*B2+One))
        A3 = GLog(B3+Sqrt(B3*B3+One))
        BrcLnR = (A2-A1 + Half*Half*(exp(two*A2)-exp(two*A1)
     $    -exp(-two*A2)+exp(-two*A1)))/(two*two*gam)
        BrcLnP = (A3-A2 + Half*Half*(exp(two*A3)-exp(two*A2)
     $    -exp(-two*A3)+exp(-two*A2)))/(two*two*gam)
        Call MultLV(.True.,1,0,N,N,Proj,Vec(1,3),Vec(1,4))
        Call AMask(N,-100000,-1,IFlags,NOut,Vec(1,4),Vec(1,4))
C
C       Normalize constraint in space of active variables
C
        Tmp = One/Sqrt(SProd(N,Vec(1,4),Vec(1,4)))
        Call AScale(N,Tmp,Vec(1,4),Vec(1,4))
C       Normalized tq(N) in Vec(i,4)
        endIf
      GrdTmp = SProd(N,Vec(1,4),F)
      If(IDist.eq.0) then
C       get XR and XP
        Call AScale(N,RP*RP,Vec(1,1),Vec(1,1))
        Call AScale(N,RR*RR,Vec(1,2),Vec(1,2))
C       Get PR
        Call ASub(N,Vec(1,1),Vec(1,2),Vec(1,3))
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,Vec(1,3),XX)
        RPR = Sqrt(SProd(N,Vec(1,3),Vec(1,3)))
        If(Iprint.gt.1) then
          Write(IOut,2000)IQST0+1,GrdTmp
          Write(IOut,2001) (RP-RR),Half*(RR+RP),RPR
          endIf
C
C       Describe the arc of circle passing through R and X
        Dot = -SProd(N,Vec(1,4),Vec(1,2))
        Theta = Dot/RR
        If(Abs(Theta).gt.One) Theta = Abs(theta)/theta
        Theta = ACos(Theta)
        AA = Tan(theta)
        SMall = MDCutO(0)
        If(Abs(theta).lt.small) then
C         Case of straight line
          S2 = RR/two
          RPR = two*S2
          ArcLnR = RR
          RR = One/Small
        else
          S2 = Abs(RR*RR/(Two*Dot))
          RPR = two*S2
          arclnR = two*S2*theta/AA
          RR = Half*arclnR/theta
          endIf
        If(DoPar) ArcLnR = BrclnR
        RR = S2
C       Do the same for X->P
        Dot = SProd(N,Vec(1,4),Vec(1,1))
        Theta = Dot/RP
        If(Abs(theta).gt.One) Theta = Abs(theta)/theta
        Theta = ACos(Theta)
        AA = Tan(theta)
        Small = MDCutO(0)
        If(Abs(Theta).lt.small) then
          S2 = RP/two
          ArcLnP = RP
          RP = One/Small
        else
          S2 = Abs(RP*RP/(two*Dot))
          arclnP = two*S2*theta/AA
          RP = Half*arclnP/Theta
          endIf
        If(DoPar) ArcLnP = BrclnP
        RP = S2
        If(IPrint.gt.1) Write(IOut,2003) ArclnR,ArcLnP
        endIf
C
C     Redistribution
C
      If(IDist.gt.0) then
        Call AScale(N,RP*RP,Vec(1,1),Vec(1,1))
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,Vec(1,1),XX)
        Call AScale(N,RR*RR,Vec(1,2),Vec(1,2))
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,Vec(1,2),XX)
        Call ASub(N,Vec(1,2),Vec(1,1),Vec(1,3))
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,Vec(1,3),XX)
C       Equidistant
        Call AAdd(N,Vec(1,2),Vec(1,1),Vec(1,5))
        Call AScale(N,-Half,Vec(1,5),Vec(1,5))
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,Vec(1,5),XX)
        Dot = SProd(N,Vec(1,3),Vec(1,3))
        Dot = SProd(N,Vec(1,3),Vec(1,5))/Dot
        Call AScale(N,Dot,Vec(1,3),Vec(1,3))
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,Vec(1,3),XX)
        Call MultLV(.True.,1,0,N,N,Proj,Vec(1,3),Vec(1,5))
        Call ASub(N,Values,Vec(1,5),Values)
        Call FixDih(NTDih,Values(NTBond+NTAng+1))
        endIf
      If(Doeuler) then
        Dot = -One/Sqrt(SProd(N,Vec(1,2),Vec(1,2)))
        Call AScale(N,Dot,Vec(1,2),Vec(1,4))
        endIf
      Return
      End
*Deck ConsST
      Subroutine ConsST(IOut,IQST0,IRad,IStep,N,IFlags,Proj,Values,
     $  Cnstr1,Cnstr2,F,FC,Climb,EigTmp,Vec,NTBond,NTAng,NTDih)
      Implicit Real*8(A-H,O-Z)
C
C     For QST and related procedures, obtain constraint vector and
C     project.
C       QST  T = (P-X)/|P-X|**2 - (R-X)/|R-X|**2
C       LST  T = P - R
C
C     Proj and FC are lower triangular.
C
      Logical Climb
      Real*8 MDCutO
      Dimension IFlags(N), Proj(*), Values(N), Cnstr1(N), Cnstr2(N),
     $  F(N), FC(*), Vec(N,5), XX(1)
      Save Half, One, XX
      Data Half/0.5d0/, One/1.0d0/, XX/0.d0/
 1010 Format(' LST/QST climbing along tangent vector')
C
C     For LST/QST(int) in internal var space, IQST0 = 1/2
C
      If(IQST0.eq.1) then
        Call ASub(N,Cnstr1,Cnstr2,Vec(1,3))
        Call FixDq(IOut,0,1,N,NTBond,NTAng,NTDih,Vec(1,3),XX)
      else if(IQST0.eq.2) then
        Call ASub(N,Cnstr1,Values,Vec(1,4))
        Call FixDq(IOut,0,1,N,NTBond,NTAng,NTDih,Vec(1,4),XX)
        Call MultLV(.True.,1,0,N,N,Proj,Vec(1,4),Vec(1,1))
        Tmp = One/SProd(N,Vec(1,1),Vec(1,1))
        Call AScale(N,Tmp,Vec(1,1),Vec(1,1))
        Call ASub(N,Cnstr2,Values,Vec(1,4))
        Call FixDq(IOut,0,1,N,NTBond,NTAng,NTDih,Vec(1,4),XX)
        Call MultLV(.True.,1,0,N,N,Proj,Vec(1,4),Vec(1,2))
        Tmp = One/SProd(N,Vec(1,2),Vec(1,2))
        Call AScale(N,tmp,Vec(1,2),Vec(1,2))
        Call ASub(N,Vec(1,2),Vec(1,1),Vec(1,3))
        endIf
      If(IQST0.eq.1.or.IQST0.eq.2) then
        Call MultLV(.True.,1,0,N,N,Proj,Vec(1,3),Vec(1,4))
        Call AMask(N,-100000,-1,IFlags,NOut,Vec(1,4),Vec(1,4))
        endIf
C
C     Normalize constraint in space of active variables
C
      Tmp = One/Sqrt(SProd(N,Vec(1,4),Vec(1,4)))
      Call AScale(N,Tmp,Vec(1,4),Vec(1,4))
C
C     Normalized tq(N) in Vec(i,4)
C
      If(IRad.gt.0) then
        Call MultLV(.True.,1,0,N,N,FC,Vec(1,4),Vec(1,2))
        EigTmp = SProd(N,Vec(1,4),Vec(1,2))
        GrdTmp = SProd(N,Vec(1,4),F)
        If(Abs(GrdTmp).gt.MDCutO(0)) then
          XTmp = Half*(EigTmp+Sqrt(EigTmp**2+(GrdTmp+GrdTmp)**2))
          XTmp = -GrdTmp/(XTmp-EigTmp)
        else
          XTmp = One
          endIf
        QSTrad = GFloat(IRad)/GFloat(100)
        Climb = IStep.le.2.or.(Abs(XTmp).gt.QSTrad.and.IStep.le.4)
        If(Climb) Write(IOut,1010)
        Call ACASB(N,Vec(1,2),Vec(1,4),Vec(1,5),-EigTmp)
        endIf
      Return
      End
*Deck D2CorL
      Subroutine D2CorL(IOut,Method,NVar,NNP,NP,JSkip,RMin,RMax,GrdErr,
     $  IC,XX,FF,FC,IS,FCOld,IFlags)
      Implicit Real*8(A-H,O-Z)
C
C     phi*Powell + (1-phi)*Murtagh-Sargent  (JCC vol 15 pg 1 1994)
C     Modified by Dr. H. B. Schlegel in 1994.  MthUpd=3 for min, 4 for TS.
C
      Dimension IC(*), XX(NVar,*), FF(NVar,*), FC(*), IS(*), FCOld(*),
     $          IFlags(*)
      Real*8 MDCutO
      Logical ChkPos
      Save Zero, One
      Data Zero/0.0d0/, One/1.d0/
 1000 Format(' Update second derivatives using D2CorL and',
     $       ' points',5I3,(/1X,55X,5I3))
      LInd(I,J) = ((Max(I,J)*(Max(I,J)-1))/2) + Min(I,J)
C
      MthUpd = IAbs(Method)
      If(NNP.eq.0) then
        NNP = NP
        I1  = NP
        I2  = 1
        ISt = -1
      else
        I1  = 1
        I2  = NP
        ISt = 1
        endIf
      ChkPos = Method.eq.-3
      DGMin = GrdErr
      Small = MDCutO(0)
      RMax2 = RMax**2
      RMin2 = RMin**2
      NN = 0
      Do 270 IP = I1, I2, ISt
        If(NN.ge.NNP) goto 270
        If(IP.eq.JSkip) goto 270
        DX = Zero
        DXDF = Zero
        DXFCDX = Zero
        DFF2 = Zero
        DG = Zero
        DO 215 I = 1, NVar
          FCOld(I) = Zero
          If(IFlags(I).lt.0) goto 215
          XXI=(XX(I,IP)-XX(I,JSkip))
          DX = DX + XXI**2
          DXDF = DXDF + XXI*(FF(I,JSkip) - FF(I,IP))
          DG = DG + (FF(I,JSkip) - FF(I,IP))**2
          Do 210 J = 1, NVar
            If(IFlags(J).lt.0) goto 210
            XXJ = (XX(J,IP)-XX(J,JSkip))
            DXFCDX = DXFCDX + XXI*FC(LInd(I,J))*XXJ
            FCOld(I) = FCOld(I) + FC(LInd(I,J))*XXJ
  210       Continue
          DFF2 = DFF2 + (FF(I,JSkip)-FF(I,IP)-FCOld(I))**2
  215     Continue
        DG = Sqrt(DG)
        If((DX.gt.RMax2.and.IP.gt.2).or.
     $    DX.lt.RMin2.or.DG.lt.DGMin) goto 270
        If(ChkPos.and.DXDF.le.Zero) goto 270
        NN = NN + 1
        Phi = Zero
        RNum = DXDF - DXFCDX
        RDenom = DX * DFF2
C
C  Phi, the Bofill weighting factor
C
        If (MthUpd.eq.7.or.MthUpd.eq.6.and.
     $      Abs(RNum).gt.Small) then
          Phi = Sqrt(RNum**2/RDenom)
          end If
        Write(IOut,'(2I4,3G15.7)') IP,JSkip,RNum,Sqrt(RDenom),Phi
C
C  The update loop
C
        IS(NN) = IC(IP)
        Do 220 I = 1, NVar
          If(IFlags(I).lt.0) goto 220
          DI = (XX(I,IP)-XX(I,JSkip))
          DFI = FF(I,JSkip) - FF(I,IP)
          DFFI = DFI - FCOld(I)
          IIX = (I*(I-1))/2
          Do 219 J = 1, I
            If(IFlags(I).lt.0) goto 219
            DJ = (XX(J,IP)-XX(J,JSkip))
            DFJ = FF(J,JSkip) - FF(J,IP)
            DFFJ = DFJ - FCOld(J)
C
C           For minima, iterate using the BFGS update
C           For saddle points, iterate using Bofill's weigted update
C           phi*Powell + (1-phi)*Murtagh-Sargent  (JCC vol 15 pg 1 1994)
C
            If(MthUpd.eq.7.or.MthUpd.eq.6) then
              If(Abs(RNum).gt.Small)
     $          FC(IIX+J) = FC(IIX+J) + Phi*DFFI*DFFJ/RNum
              end If
            If(MthUpd.eq.3.or.MthUpd.eq.6) then
              FC(IIX+J) = FC(IIX+J) + (One-Phi) * (DFI*DFJ/DXDF
     $          - FCOld(I)*FCOld(J)/DXFCDX)
            else if (MthUpd.eq.7) then
              FC(IIX+J) = FC(IIX+J) + (One-Phi) * ((DFFI*DJ+DI*DFFJ)/DX-
     $          (DXDF-DXFCDX)*DI*DJ/(DX*DX))
              endIf
  219       Continue
  220     Continue
  270   Continue
      NN = NN + 1
      IS(NN) = IC(JSkip)
      If(IPSave(0).le.1) Write(IOut,1000) (IS(I),I=1,NN)
      Return
      End
*Deck D2CorP
      Subroutine D2CorP(IOut,MthUpd,NVar,NP,JSkip,RMin,RMax,GrdErr,IC,
     $  XX,FF,FC,IS,FCOld,T)
      Implicit Real*8(A-H,O-Z)
C For Path relaxation. Variant of D2CorX
C
C    MthUpd 3 ... BFGS/Powell Mixture following Powell
C           4 ... Bofill
C           5 ... BFGS
C           6 ... Powell
C           7 ... BFGS/Powell Mixture following Bofill
C
C
      Logical Debug
      Dimension IC(*), XX(NVar,*), FF(NVar,*), FC(*), IS(*), FCOld(*)
      Dimension T(Nvar)
      Real*8 MDCutO
      Save Zero, One
      Data Zero/0.0d0/, One/1.d0/
 1000 Format(' Update second derivatives using information from',
     $       ' points',5I3,(/1X,53X,5I3))
 1010 Format(' Updating Method = ',I3)
      LInd(I,J) = ((Max(I,J)*(Max(I,J)-1))/2) + Min(I,J)
C
      Debug=.False.
      Write(IOut,1010) MthUpd
      DGMin = GFloat(100)*GrdErr
      NVTT=(NVar*(Nvar+1))/2
      Small = MDCutO(0)
      RMax2 = RMax**2
      RMin2 = RMin**2
      NN = 0
      If(MthUpd.eq.3) then
        Call AMove(Nvar,FF(1,1),T)
        Dot = One/Sqrt(SProd(NVar,T,T))
        Call AScale(NVar,Dot,T,T)
        endIf
      Do 270 IP = NP, 1, -1
        If(IP.eq.JSkip) goto 270
        DX = Zero
        DXDF = Zero
        DXFCDX = Zero
        DFF2 = Zero
        DG = Zero
        DXF=Zero
        G=Zero
        DXT=Zero
        TT=Zero
        DO 215 I = 1, NVar
          XXI=(XX(I,IP)-XX(I,JSkip))
          DX = DX + XXI**2
        If(Mthupd.eq.3.or.Mthupd.eq.7)TT=TT+T(I)*T(I)
        If(Mthupd.eq.3.or.MthUpd.eq.7)DXT=DXT+XXI*T(I)
          DXDF = DXDF + XXI*(FF(I,JSkip) - FF(I,IP))
          DXF = DXF + XXI*(FF(I,JSkip) )
          DG = DG + (FF(I,JSkip) - FF(I,IP))**2
          G = G + (FF(I,JSkip) )**2
          FCOld(I) = Zero
          Do 210 J = 1, NVar
            XXJ = (XX(J,IP)-XX(J,JSkip))
            DXFCDX = DXFCDX + XXI*FC(LInd(I,J))*XXJ
  210       FCOld(I) = FCOld(I) + FC(LInd(I,J))*XXJ
  215     DFF2=DFF2+(FF(I,JSkip) - FF(I,IP)-FCOld(I))**2
        DG = Sqrt(DG)
        If(DX.gt.RMax2.or.DX.lt.RMin2.or.DG.lt.DGMin) goto 270
        NN = NN + 1
        IS(NN) = IC(IP)
        If(mthupd.eq.3.or.Mthupd.eq.7)then
        If(Debug)then
        Write(IOut,*)' DX**2 ',DX
        Write(IOut,*)' DXDF ',DXDF
        Write(IOut,*)' DXFCDX ',DXFCDX
        Write(IOut,*)' TT ',TT
        Write(IOut,*)' DXT ',DXT
        Endif
        Phi=DXT*DXT/(DX*TT)
        If(Debug)Write(IOut,*)' Phi=',Phi
       Endif
        Do 220 I = 1, NVar
          DI = (XX(I,IP)-XX(I,JSkip))
          DFI = FF(I,JSkip) - FF(I,IP)
          DFFI = DFI - FCOld(I)
          IIX = (I*(I-1))/2
          Do 220 J = 1, I
            DJ = (XX(J,IP)-XX(J,JSkip))
            DFJ = FF(J,JSkip) - FF(J,IP)
            DFFJ = DFJ - FCOld(J)
C
C           For minima, iterate using the Powell/bfgs update
C           For saddle points, iterate using Bofill's weigted update
C           phi*Powell + (1-phi)*Murtagh-Sargent  (JCC vol 15 pg 1 1994)
C
            If(MthUpd.eq.3.or.MthUpd.eq.5.or.MthUpd.eq.6) then
              If(MthUpd.eq.5)Phi=Zero
              If(MthUpd.eq.6)Phi=One
              If(DXDF.lt.Small.or.DXFCDX.lt.Small
     $        .and.mthupd.ne.6)Phi=One
              If(MthUpd.eq.3.and.IP.gt.2)Phi=One
              FC(IIX+J) = FC(IIX+J) + (Phi)*((DFFI*DJ+DI*DFFJ)/DX-
     $          (DXDF-DXFCDX)*DI*DJ/(DX*DX))
              FC(IIX+J) = FC(IIX+J) + (One-Phi)*(DFI*DFJ/DXDF
     $          - FCOld(I)*FCOld(J)/DXFCDX)
            else if(MthUpd.eq.4.or.(MthUpd.eq.7.and.IP.gt.2)) then
              RNum = DXDF - DXFCDX
              RDenom = DX * DFF2
              If(Abs(RNum).gt.Small.and.Abs(RDenom).gt.Small) then
                Phi = One - RNum**2/RDenom
              else
                Phi = One
                endIf
              FC(IIX+J) = FC(IIX+J) + Phi*((DFFI*DJ+DI*DFFJ)/DX-
     $          (DXDF-DXFCDX)*DI*DJ/(DX*DX))
              If(Abs(RNum).gt.Small)
     $          FC(IIX+J) = FC(IIX+J) + (One-Phi)*DFFI*DFFJ/RNum
            else if(MthUpd.eq.7.and.IP.eq.2) then
             If(DXDF.lt.Small.or.DXFCDX.lt.Small)Phi=One
             Phi=Sqrt(Phi)
             FC(IIX+J) = FC(IIX+J) + (Phi)*((DFFI*DJ+DI*DFFJ)/DX-
     $         (DXDF-DXFCDX)*DI*DJ/(DX*DX))
             FC(IIX+J) = FC(IIX+J) + (One-Phi)*(DFI*DFJ/DXDF
     $         - FCOld(I)*FCOld(J)/DXFCDX)
              endIf
  220       Continue
      XMax=ARRMax(FC,NVTT,.true.)
      If(Debug)Write(Iout,*)' Maximum Hessian matrix element =',XMax
  270     Continue
      NN = NN + 1
      IS(NN) = IC(JSkip)
      If(IPSave(0).le.1) Write(IOut,1000) (IS(I),I=1,NN)
      Return
      End
*Deck DgQST0
      Subroutine DgQST0(IOut,DoRed,FCart,IQST,IStep,NVar,BigEig,Wght,
     $  Climb,EigTmp,IFlags,Values,Cnstr1,Cnstr2,FTemp,Scr,Proj,
     $  FCTemp,NTBond,NTAng,NTDih,FCSq)
      Implicit Real*8(A-H,O-Z)
C
C     Setup the force constant matrix for diagonalization when using
C     redundant coordinates and/or doing QST.
C
      Logical DoRed, FCart, Climb
      Dimension IFlags(NVar), Values(NVar), Cnstr1(NVar), Cnstr2(NVar),
     $  FTemp(NVar), Scr(NVar,*), Proj(*), FCTemp(*), FCSq(NVar,NVar)
      Save One
      Data One/1.0d0/
C
      IQST0 = Mod(IQST,10)
      IRad = Mod(IQST,100)/10
      If(IRad.eq.0) IRad = 5
      NVarTT = (NVar*(NVar+1))/2
      If(DoRed) Call PFPN(.True.,1,NVar,NVarTT,Proj,FCTemp,FCSq,Scr)
      If(IQST.gt.0) then
        If(FCart.or..not.DoRed)
     $    Call GauErr('QST only with redun. coordinates.')
        Call ConsST(IOut,IQST0,IRad,IStep,NVar,IFlags,Proj,Values,
     $    Cnstr1,Cnstr2,FTemp,FCTemp,Climb,EigTmp,Scr(1,2),NTBond,NTAng,
     $    NTDih)
        endIf
      If(DoRed) then
        Call Square(FCTemp,FCSq,NVar,NVar,0)
        Call Square(Proj,Proj,NVar,NVar,0)
        Call Ad1MPS(NVar,BigEig,Proj,FCSq)
        Call Linear(FCSq,FCTemp,NVar,NVar)
        Call Linear(Proj,Proj,NVar,NVar)
        endIf
      If(IQST.gt.0.and.Climb)
     $  Call FCClim(NVar,One-Wght,EigTmp,Scr(1,3),Scr(1,5),FCTemp)
      Return
      End
*Deck DgQST1
      Subroutine DgQST1(IOut,IQST0,NVar,IVSwap,SMax,TSVec,EVec,Eigen,
     $  Scr1)
      Implicit Real*8(A-H,O-Z)
C
C     Do some printing and setup for other printing during QST
C     diagonalization.
C
      Dimension TSVec(NVar), EVec(NVar), Eigen(NVar), Scr1(NVar,*)
 1010 Format(' LST in optimization variable space.')
 1020 Format(' QST in optimization variable space.')
 1030 Format(' LST in distance matrix space.')
 1040 Format(' QST in distance matrix space.')
 1050 Format(' Eigenvectors 1 and',i4,' swapped, overlap=',F8.4)
C
      If(IQST0.eq.1) then
        Write(IOut,1010)
      else if(IQST0.eq.2) then
        Write(IOut,1020)
      else if(IQST0.eq.3) then
        Write(IOut,1030)
      else if(IQST0.eq.4) then
        Write(IOut,1040)
        endIf
      If(IVSwap.ne.1) Write(IOut,1050) IVSwap, SMax
      Call AMove(NVar,TSVec,Scr1)
      Call AMove(NVar,EVec,Scr1(1,2))
      Call AMove(NVar,Eigen,Scr1(1,4))
      Return
      End
*Deck DgQSTk
      Subroutine DgQSTk(DoRed,IQST,NVar,BigEig,GrdTmp,IFlags,Values,
     $  Cnstr1,Cnstr2,FTemp,Scr,Proj,FCTemp,Vec,NTBond,NTAng,NTDih,
     $  doeuler)
      Implicit Real*8(A-H,O-Z)
C
C     For path relaxation. Variant of dgqst0: Cnsst4 called instead of consst.
C
      Logical DoRed,doeuler
      Dimension IFlags(NVar), Values(NVar), Cnstr1(NVar), Cnstr2(NVar),
     $  FTemp(NVar), Scr(NVar,5),  FCTemp(*), Vec(NVar,NVar),
     $  Proj((NVar*(Nvar+1))/2)
      Common /IO/     In, IOut, IPunch
C
      IQST0 = Mod(IQST,10)
      If(DoRed) Call PFP(NVar,Proj,FCTemp,Vec)
      If(IQST.gt.0) Call CnsST4(IQST0,NVar,IFlags,Proj,Values,Cnstr1,
     $  Cnstr2,FTemp,GrdTmp,Scr,NTBond,NTAng,NTDih,RP,RR,RPR,ArcLnR,
     $  ArcLnP,0,1,doeuler)
      If(DoRed) Call Add1MP(NVar,BigEig,Proj,FCTemp)
      Return
      End
*Deck DgQSTL
      Subroutine DgQSTL(IOut,IPrint,IQST,N,NP,X1,X2,X,FF,W,NTBond,NTAng,
     $  NTDih)
      Implicit Real*8(A-H,O-Z)
C
C     For QST and related procedures, obtain constraint vector and
C     project.
C       QST  T = (X1-X)/|X1-X|**2 - (X2-X)/|X2-X|**2
C       LST  T = X1 - X2
C
C     Associated surface approximation by O. Farkas, 1997, WSU.
C     Code is the modification of the ConsST routine
C
      Dimension X1(N),X2(N),X(N),FF(N,NP),W(N,*),XX(1)
      Save One, Two, XX
      Data One/1.0d0/, Two/2.d0/, XX/0.d0/
C
      IQST0 = Mod(IQST,10)
      Call AClear(N,W)
      If(IQST0.eq.1) then
        Call ASub(N,X1,X2,W)
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,W,XX)
      else if(IQST0.eq.2) then
        Call ASub(N,X1,X,W(1,2))
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,W(1,2),XX)
        Tmp = One/SProd(N,W(1,2),W(1,2))
        Call AScale(N,Tmp,W(1,2),W(1,2))
        Call ASub(N,X2,X,W(1,3))
        Call FixDq(IOut,IPrint,1,N,NTBond,NTAng,NTDih,W(1,3),XX)
        Tmp = One/SProd(N,W(1,3),W(1,3))
        Call AScale(N,Tmp,W(1,3),W(1,3))
        Call ASub(N,W(1,2),W(1,3),W)
        endIf
      If(IQST0.eq.1.or.IQST0.eq.2) then
C
C     Normalize tq in W
C
        Tmp = One/Sqrt(SProd(N,W,W))
        Call AScale(N,Tmp,W,W)
C
C     Change the gradients to make it like around a minimum
C
        Do 10 I=1, NP
          S = SProd(N,FF(1,I),W)
          S = -Two * S
          Call ACASB(N,FF(1,I),W,FF(1,I),S)
   10     continue
        endIf
      Return
      End
*Deck DiagF3
      Subroutine DiagF3(IOut,IPrint,FixEig,SclEig,FCart,NVar,eigmin,
     $  EigMax,IFlags,FC,FCTemp,XName,Vec,Eigen,VecInt,NTRrot,Scr,
     $  FTemp,FTempU,Neg,NNegE,NEig,NAtoms,C,AtmChg,IQST,DoRed,Values,
     $  Cnstr1,Cnstr2,NTBond,NTAng,NTDih,SS2,XM,XQuad,doeuler,Proj,
     $  USeSym,ArcL,NFroz,Frozen)
      Implicit Real*8(A-H,O-Z)
C
C     Diagonalize the force constant matrix, count the negative
C     eigenvalues, and transform the forces in FTemp to the eigenvector
C     basis in FTempU, which can be the same as Scr.  FixEig determines
C     whether eigenvalues having the wrong sign are flipped.   When using
C     redundant coordinates, VecInt should contain the lower-triangular
C     projection matrix on entry.
C
C     Routine similar to Diagfc but modified for points along
C     The path. Compute a GS-like step to measure stability of
C     the path
C
      Parameter (ITSV=8)
      Logical DoRed, FixEig, SclEig, FCart, DoEuler, UseSym, Print,
     $  Failed,Frozen(*)
      Dimension IFlags(NVar), FC(*), FCTemp(*), Vec(NVar,NVar),
     $  Eigen(NVar), Scr(NVar,*), FTemp(NVar), FTempU(NVar),
     $  C(3,NAtoms), XName(*), AtmChg(NAtoms), VecInt(NVar,NVar),
     $  Values(NVar), Cnstr1(NVar), Cnstr2(NVar), XM(NVar), XQuad(Nvar),
     $  Proj((NVar*(Nvar+1)/2)),XX(1)
       Save BigEig,Zero,One,Two,FTol,SmlEig
       Data BigEig/1.d3/,Zero/0.0d0/,One/1.0d0/,Two/2.0d0/,
     $   FTol/1.0d-06/,SmlEig/1.d-6/,XX/0.d0/
 1000 Format('  End-Point of Bimolecular Reaction ')
 1010 Format('  Eigenvectors of the second derivative matrix:')
 1020 Format('  Eigenvectors required to have negative eigenvalues:')
 1040 Format('     Eigenvalues --- ',5F10.5)
 1060 Format('  ==Forces against Tangent == ')
 1070 Format('  Maximum :',F10.5,' RMS :',F10.5)
 1160 Format('  == Quad Disp. == ')
 1170 Format('  Maximum :',F10.5,' RMS :',F10.5)
 1180 Format('  GS Constrained optimization Radius = ',F10.5)
 1200 Format('  Point possibly on a ridge')
 1300 Format('  Angle (T,V) = ',F10.5,' Degrees')
C
      Print=.False.
      If(Doeuler) Write(IOut,1000)
      Call FCMask(NVar,-100000,-1,BigEig,IFlags,NOut,FC,FCTemp)
      IQST0 = Mod(IQST,10)
      Call DgQSTk(DoRed,IQST,NVar,BigEig,GrdTmp,IFlags,Values,Cnstr1,
     $  Cnstr2,FTemp,Scr,Proj,FCTemp,Vec,NTBond,NTAng,NTDih,doeuler)
      Call AClear(NVar,Scr(1,3))
      Dot = Sqrt(SProd(Nvar,Ftemp,FTemp))
      If(Dot.gt.FTol) then
        Dot = One / Dot
      else
        Dot = Zero
        endIf
      Call AScale(Nvar,Dot,Ftemp,Scr(1,3))
      Call GetTan(IOut,1,NVar,.True.,Scr(1,2),XM,Scr(1,3),Scr(1,ITSV),S,
     $  Proj,Scr(1,5),ArcL)
C     For cartesian optimizations, freeze translations and rotations
      If(FCart) then
        Call PrjDag(IOut,IPrint,NTrRot,NAtoms,NVar,AtmChg,C,NFroz,
     $    Frozen,NEig,Eigen,Vec,FCTemp,Scr,VecInt,Scr,0)
      else
        NEig = NVar
        Call DiagDN(IOut,IPrint,0,FCTemp,Vec,Eigen,NVar,Scr,6,Scr(1,2),
     $    NVar,.False.)
        endIf
C
C     For QST and related procedures, swap eigenvectors so that
C     the first one has the maximum overlap with the constraint
C     IVSwap marks the final eigenvec swapped to the first position
C
      If(IPrint.gt.0) then
        Write(IOut,1010)
        Call MatPrt(Vec,NVar,NEig,NVar,NEig,2,0,XName,XName,0,Eigen,1)
      else
        Write(IOut,1040) (Eigen(I),I=1,NEig)
        If(Neg.gt.0) then
          Write(IOut,1020)
          Do 10 I = 1, Neg
            Call MatPrt(Vec(1,I),1,NVar,1,NVar,0,1,XName,XName,0,Eigen,
     $        0)
   10       Continue
          endIf
        endIf
C     Test eigenvalues of 2nd derivative matrix.
      Call ChkEig(IOut,.False.,SclEig,FixEig,NEig,NVar,Neg,EigMin,
     $  EigMax,BigEig,SmlEig,Zero,Zero,.True.,XName,NNegE,Eigen,Vec,
     $  Scr(1,5),Failed)
      IVSwap = 0
      If(IQST.ne.0)
     $  Call EigOvl(NVar,NVar,Scr(1,ITSV),Eigen,Vec,IVSwap,SMax,Scr)
      If(IQST.gt.0) Call DgQST1(IOut,IQST0,NVar,IVSwap,SMax,Scr(1,ITSV),
     $  Vec,Eigen,VecInt)
      Call AClear(NVar,Scr(1,5))
      DotT = SProd(NVar,VecINt(1,1),VecInt(1,1))
      DotG = SProd(NVar,VecInt(1,1),FTemp)
      DotF = SProd(Nvar,FTemp,Ftemp)
      Dot = DotG/(Sqrt(DotF*DotT))
      Ang = ACos(Dot)
      XPi4 = ATan(One)
      Ang = 45.0d0*Ang/XPi4
      Write(IOut,1300) Ang
      DotG = DotG / DotT
      Call AScale(NVar,DotG,VecInt(1,1),Scr(1,5))
      Call ASub(NVar,FTemp,Scr(1,5),Scr(1,5))
      DotP = ARRMax(Scr(1,5),NVar,.True.)
      RMSP = SProd(NVar,Scr(1,5),Scr(1,5))
      RMSP = Sqrt(RMSP/GFLoat(NVar))
      If(Print) then
        Write(IOut,1060)
        Write(IOut,1070) DotP, RMSP
        endIf
C     Form tangent and force in eigenvector basis.
      Call MatMP1(1,0,1,1,1,NVar,NEig,VecInt(1,1),Vec,VecInt(1,2))
      Call MatMP1(1,0,1,1,1,NVar,NEig,FTemp,Vec,FTempU)
      Call AClear(NVar,Scr(1,6))
      Call AMove(NVar,VecInt(1,2),Scr(1,6))
      DotT = SProd(NEig,VecINt(1,2),VecInt(1,2))
      DotG = SProd(NEig,VecInt(1,2),FTempU)
      DotG = DotG/DotT
      Call AScale(NEig,DotG,VecInt(1,2),VecInt(1,2))
      Call ASub(NEig,FTempU,VecInt(1,2),VecInt(1,2))
C     Transform back the complentary force
      Call AClear(NVar,Scr(1,5))
      Do 11 II=1,Neig
        Temp =VecInt(II,2)
        Call ACasB(NVar,Scr(1,5),Vec(1,II),SCR(1,5),Temp)
   11   Continue
      Call AMove(NVar,Scr(1,5),VecINt(1,2))
      S = S*Two
      SS2 = S/two
      If(Print) Write(IOut,1180) SS2
      Call AClear(NEig,Scr(1,5))
      DotT = One/Sqrt(SProd(Nvar,Scr(1,6),Scr(1,6)))
      Call Ascale(Nvar,DotT,Scr(1,6),Scr(1,6))
      Call AScale(NEig,-One,FtempU,Scr(1,5))
      Call Ascale(Neig,SS2,Scr(1,6),Scr(1,6))
      Call DxGS(IOut,NVar,Scr(1,5),Scr(1,6),Vec,Eigen,S,
     $  XQuad,FTol,USeSym,.false.)
      Call MultLV(.True.,1,0,NVar,NVar,Proj,XQuad,Scr(1,5))
      Call AMove(Nvar,Scr(1,5),XQuad)
      Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,XQuad,XX)
      DotQ = ARRMax(XQuad,NVar,.True.)
      RMSQ = SProd(NVar,XQuad,XQuad)
      RMSQ = Sqrt(RMSQ/GFloat(NVar))
      If(Print) then
        Write(IOut,1160)
        Write(IOut,1170) DotQ, RMSQ
        endIf
      If(IVSwap.gt.1) Write(IOut,1200)
      If(IQST.gt.0) then
C       Scr(1,4) is the transition vector
        If(ITSV.ne.4) Call AMove(NVar,Scr(1,ITSV),Scr(1,4))
        If(IPrint.ge.2) then
          call Amove(NVar,FTempU,VecInt(1,3))
          Call QSTOut(2,NVar,VecInt,XName)
          endif
        Call AClear(NVar-NEig,FTempU(NEig+1))
        endIf
      Return
      End
*Deck DiagFC
      Subroutine DiagFC(IOut,IPrint,FixEig,SclEig,FCart,AlUnit,NTrRot,
     $  NVar,EigMin,EigMax,IFlags,FC,FCTemp,XName,Vec,IDimSc,Scr,SimPrj,
     $  NRank,VecInt,FTemp,Neg,NNegE,NEig,NAtoms,C,AtmChg,NFroz,Frozen,
     $  Wght,IQST,IStep,DoRed,Values,Cnstr1,Cnstr2,NTBond,NTAng,NTDih,V,
     $  MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Diagonalize the force constant matrix, count the negative
C     eigenvalues, and transform the forces in FTemp to the eigenvector
C     basis in FTempU=Scr(,2), which can be the same as Scr.  FixEig
C     determines whether eigenvalues having the wrong sign are flipped.
C     When using redundant coordinates, VecInt should contain the
C     lower-triangular projection matrix on entry.  Scr(,1) receives
C     the eigenvalues and Scr(,2) the forces in the eigenvector basis.
C
      Logical DoRed, FixEig, SclEig, FCart, Climb, TstUnL, AlUnit,
     $  Failed, SimPrj, Frozen(*)
      Parameter (NPrtDf=10,ITSV=8)
      Dimension IFlags(NVar), FC(*), FCTemp(*), Vec(NVar,NVar),
     $  Scr(NVar,*), FTemp(NVar), C(3,NAtoms), XName(*),
     $  AtmChg(NAtoms), VecInt(NVar,NVar), Values(NVar), Cnstr1(NVar),
     $  Cnstr2(NVar), V(*)
      Save Zero, One, BigEig, SmlEig
      Data Zero/0.0d0/, One/1.0d0/, BigEig/1.d3/, SmlEig/1.d-6/
 1010 Format(' Eigenvectors of the second derivative matrix:')
 1020 Format(' Eigenvectors required to have negative eigenvalues:')
 1040 Format('     Eigenvalues --- ',5F10.5)
C
      Call TStamp(1,'Top DiagFC')
      NOut = NVar
      If(.not.AlUnit)
     $  Call FCMask(NVar,-100000,-1,BigEig,IFlags,NOut,FC,FCTemp)
      IQST0 = Mod(IQST,10)
      IScale = Mod(IQST-IQST0,100000)/10000
      If(IScale.eq.0) IScale = 5
      FScale = GFloat(IScale) / GFloat(10)
      If(SimPrj) then
        If(NOut.ne.NVar)
     $    Call GauErr('NOut.ne.NVar with SimPrj in DiagFC.')
      else
        Call DgQST0(IOut,DoRed,FCart,IQST,IStep,NVar,BigEig,Wght,Climb,
     $    EigTmp,IFlags,Values,Cnstr1,Cnstr2,FTemp,Scr,VecInt,FCTemp,
     $    NTBond,NTAng,NTDih,Vec)
        If(ITSV.ne.5) Call AMove(NVar,Scr(1,5),Scr(1,ITSV))
        endIf
C
C     For cartesian optimizations, freeze translations and rotations
C
      If(FCart) then
        If(AlUnit) then
          NEig = NVar
          Call ASet(NEig,One,Scr)
        else if(TstUnL(NVar,Zero,FCTemp)) then
          NEig = NVar
          Call ASet(NEig,One,Scr)
          Call AUnitM(.False.,1,NEig,NEig,Vec)
        else
          Call PrjDag(IOut,IPrint,NTrRot,NAtoms,NVar,AtmChg,C,NFroz,
     $      Frozen,NEig,Scr,Vec,FCTemp,Scr(1,2),VecInt,V,MDV)
          endIf
      else
        If(SimPrj) then
          NEig = NRank
          Call VFVMul(.False.,.False.,.False.,0,NVar,NEig,NVar,FCTemp,
     $      VecInt,Vec,FCTemp,V,MDV)
        else
          NEig = NVar
          endIf
        NE2 = MDV / NEig
        If(NE2.gt.6) then
          Call DiagDN(IOut,IPrint,0,FCTemp,Vec,Scr,NEig,Scr(1,2),NE2,V,
     $      NEig,.False.)
        else
          Call DiagDN(IOut,IPrint,0,FCTemp,Vec,Scr,NEig,Scr(1,2),6,
     $      Scr(1,3),NEig,.False.)
          endIf
        If(SimPrj) then
          Call XGEMM(1,'N','N',NVar,NEig,NEig,One,VecInt,NVar,Vec,NRank,
     $      Zero,FCTemp,NVar)
          Call AMove(NVar*NEig,FCTemp,Vec)
          NRed = NVar - NEig
          Call AMove(NVar*NRed,VecInt(1,NEig+1),Vec(1,NEig+1))
          Call ASet(NRed,BigEig,Scr(NEig+1,1))
          endIf
        endIf
      If(.not.AlUnit) then
        If(IPrint.gt.0) then
          Write(IOut,1010)
          Call MatPrt(Vec,NVar,NEig,NVar,NEig,2,0,XName,XName,0,Scr,1)
        else
          Write(IOut,1040) (Scr(I,1),I=1,NEig)
          If(Neg.gt.0) Write(IOut,1020)
          NDo = Min(NPrtDf,NVar)
          Do 10 I = 1, Neg
            Call AHpSrA(NVar,Scr(1,2),Vec(1,I))
            Call IRevrs(NVar,Scr(1,2))
            Call MatPt2(IOut,1,Vec(1,I),1,1,1,1,NVar,1,NDo,0,11,XName,
     $        XName,0,Scr,0,I-1,0,Scr(1,2),Scr(1,2))
   10       Continue
          endIf
C       Test eigenvalues of 2nd derivative matrix.
        Call ChkEig(IOut,.False.,SclEig,FixEig,NEig,NVar,Neg,EigMin,
     $    EigMax,BigEig,SmlEig,Zero,Zero,.True.,XName,NNegE,Scr,Vec,
     $    Scr(1,2),Failed)
        endIf
C
C     For QST and related procedures, swap eigenvectors so that
C     the first one has the maximum overlap with the constraint
C     IVSwap marks the final eigenvec swapped to the first position
C
      IVSwap = 0
      If(IQST.ne.0) Call EigOvl(NVar,NVar,Scr(1,ITSV),Scr,Vec,IVSwap,
     $  SMax,Scr(1,2))
      If(IQST.gt.0) Call DgQST1(IOut,IQST0,NVar,IVSwap,SMax,Scr(1,ITSV),
     $  Vec,Scr,VecInt)
      If(AlUnit) then
        Call AMove(NVar,FTemp,Scr(1,2))
      else
C       Form forces in eigenvector basis.
        Call MatMP1(1,0,1,1,1,NVar,NEig,FTemp,Vec,Scr(1,2))
        endIf
      If(IQST.gt.0) then
C
C       Fix forces if in climbing phase of QST and related procedures
C
        If(IQST.gt.0.and.Climb.and.NVar.gt.1) then
          Call AScale(NVar-1,FScale,Scr(2,2),Scr(1,6))
          Call AMove(NVar-1,Scr(1,6),Scr(2,2))
          endIf
        If(IPrint.ge.0) then
          Call AMove(NVar,Scr(1,2),VecInt(1,3))
          Call QSTOut(1,NVar,VecInt,XName)
          endIf
      else
        Call AClear(NVar-NEig,Scr(NEig+1,2))
        If(IPrint.ge.2) Call OutMtS(IOut,'Forces in eigenvector basis:',
     $    0,0,Scr(1,2),NVar,1,NEig,1)
        endIf
      If(IDimSc.gt.NVar) Call UnPcck(0,Junk,Scr,IDimSc,2,NVar,2)
      Call TStamp(1,'Bot DiagFC')
      Return
      End
*Deck DoLMDS
      Subroutine DoLMDS(IOut,IPrint,IMDS,NP,NUsed,NQ,HInv,QQ,FF,
     $  DQDIIS,DQRes,EigV,W,InitNR,OK,RLam,ELow,ConvEq,NG,NA,IB,B,
     $  IFlags)
      Implicit Real*8 (A-H,O-Z)
C
C     Perform multidimensional search using improved GDIIS with large
C     molecular option.
C
C     Parameters:
C     IOut   : Output channel
C     IPrint : Debug print level
C     IMDS   : MDS optimizer option flag
C              Method = Mod(IMDS,10)     - Reserved, not used
C              IMix   = Mod(IMDS,100)/10 - Not yet implemented
C              IChkC  = Mod(IMDS,1000)/100
C              ICos   = Mod(IMDS,10000)/1000
C              See the GDIIS routine for more details about IChkC and ICos
C     NP     : Number of saved variables and error vectors,
C              number of used vectors on output
C     NQ     : Dimension of the problem
C     HInv   : Guess inverse of the Hessian
C     QQ     : Saved variable vectors, last (or prefered) first
C     FF     : Saved forces, same order than in QQ
C     DQDIIS : QQ * C
C     DQRes  : DQQ * C
C     EigV   : Working space for finding the lowest eigenvalue.
C              If InitNR false it contains the guess for the
C              eigenvector related to the lowest eigenvalue.
C     W      : Working space, allocated by a previous IndGDI call
C     InitNR : True if HInv does not contain a guess.
C     OK     : DoMDS could produce a resonable step
C
C     Program code by O. Farkas
      Dimension QQ(NQ,NP), FF(NQ,NP), DQDIIS(NQ), DQRes(NQ), W(*),
     $  HInv(*), IB(4,NQ), B(3,4,NQ), IFlags(*), EigV(*)
      Logical OK, InitNR, UpdLam, UseRFO
      Save Zero, One
      Data Zero/0.0d0/, One/1.0d0/
C
      OK     = .True.
      NAt3   = 3*NA
      NatT   = (NAt3*(NAt3+1))/2
      IHX    = NAtT + 1
      IPrnt  = IPrint
      IFail  = 0
      Scale  = One
      Method = Mod(IMDS,10)
C     This feature is not implemented yet.
C     IMix   = Mod(IMDS,100)/10
      IChkC  = Mod(IMDS,1000)/100
      ICos   = Mod(IMDS,10000)/1000
      IStart = 0
      Call IndGDI(NP,NQ,.True.,IStart,IQQS,IDQQS,IA,IAInv,ICC,ICSave,
     $  IHSc,IDQD,IDQR,IW,IEnd)
      UseRFO = Method.eq.1
      UpdLam = UseRFO
      If(.not.UseRFO) RLam = Zero
      N0 = 0
      N1 = 2
C     Here we go into the main loop
   10 Continue
      If(N1.le.NP.and.IFail.eq.0) then
C       Prepare the error vectors
        Do 210 I = N0+1, N1
          III = IDQQS+(I-1)*NQ
          NOut = NQ
          Call AMask(NQ,-100000,-1,IFlags,NOut,FF(1,I),W(IDQR))
          Call NRStep(IOut,IPrint,InitNR,UpdLam,ConvEq,RLam,ELow,NQ,
     $      HInv(IHX),HInv,W(IDQR),W(III),EigV,NG,NA,IB,B)
          If(NOut.ne.NQ)
     $      Call AMask(NQ,-100000,-1,IFlags,NOut,W(III),W(III))
          UpdLam = .False.
  210     Continue
        If(N0.eq.0) then
C         Calculate the scaling factors
          Do 220 I = 1, NQ
            IIHSc = IHSc + I - 1
  220       W(IIHSc) = One
          endIf
C       Scale the error vectors and the variables
        Do 230 I = N0, (N1-1)
          IShft = I*NQ
          Call AMul(NQ,QQ(1,I+1),W(IHSc),W(IQQS+IShft))
          Call AMul(NQ,W(IDQQS+IShft),W(IHSc),W(IDQQS+IShft))
  230     Continue
        If(N1.eq.NP) IPrnt = Max(1,IPrint)
        Call GDIIS(IOut,IPrnt,NQ,NQ,N1,NUsed,N0,.True.,ICos,IChkC,
     $    W(IQQS),W(IDQQS),W(IA),W(IAInv),Scale,W(ICC),W(ICSave),
     $    W(IDQD),W(IDQR),W(IW),IFail)
        If(IFail.eq.0.and.NUsed.lt.N1) IFail = 1
        N0 = NUsed
        N1 = N0 + 1
C       Go back and check if we have something more to do
        Goto 10
        endIf
C
C     Scale back the DIIS step
C
      Call ADiv(NQ,W(IDQD),W(IHSc),DQDIIS)
C
C     Scale back the residuum
C
      Call ADiv(NQ,W(IDQR),W(IHSc),DQRes)
      Call AMove(NQ,W(IHSc),W)
      OK = .True.
      If(IFail.gt.1) NP = NUsed
      Return
      End
*Deck DoMDS
      Subroutine DoMDS(IOut,IPrint,IMDS,NP,NUsed,NQ,NNRed,NNeg,EH,UH,QQ,
     $  FF,DQDIIS,DQRes,W,OK)
      Implicit Real*8(A-H,O-Z)
C
C     Perform multidimensional search using improved GDIIS.
C
C     Parameters:
C     IOut   : Output channel
C     IPrint : Debug print level
C     IMDS   : MDS optimizer option flag
C              Method = Mod(IMDS,10)     - Reserved, not used
C              IChkC  = Mod(IMDS,1000)/100
C              ICos   = Mod(IMDS,10000)/1000
C              See the GDIIS routine for more details about IChkC and ICos
C     NP     : Number of saved variables and error vectors,
C              number of used vectors on output
C     NQ     : Dimension of the problem
C     NNRed  : Number of non-redundant coordinates
C     NNeg   : Index of critical point to find
C     EH     : Eigenvalues of the Hessian
C     UH     : Eigenvectors of the Hessian
C     QQ     : Saved variable vectors, last (or prefered) first
C     FF     : Saved forces, same order than in QQ
C     DQDIIS : QQ * C
C     DQRes  : DQQ * C
C     W      : Working space, allocated by a previous IndGDI call
C     OK     : DoMDS could produce a resonable step
C
C     Program code by
C     O. Farkas 1991-1997 at
C     Department of Organic Chemistry,
C     Eotvos University of Budapest, Budapest 112, P.O.Box 32, Hungary
C     and
C     Department of Chemistry,
C     Wayne State University, Detroit, Michigan 48202
C
      Dimension QQ(NQ,NP), FF(NQ,NP), DQDIIS(NQ), DQRes(NQ), W(*),
     $  EH(NQ), UH(NQ,NQ)
      Logical OK
      Save One, EigMin
      Data One/1.0d0/, EigMin/1.d-3/
C
      OK     = .True.
      IPrnt  = IPrint
      IFail  = 0
      Scale  = One
      IChkC  = Mod(IMDS,1000)/100
      ICos   = Mod(IMDS,10000)/1000
      IStart = 0
      Call IndGDI(NP,NQ,.False.,IStart,IQQS,IDQQS,IA,IAInv,ICC,ICSave,
     $  IHSc,IDQD,IDQR,IW,IEnd)
C
C   Check and adjust the eigenvalues of the Hessian
C
      Call AAbs(NNRed,EH,W(IHSc))
      Call AChMin(NNRed,1,EigMin,W(IHSc),W(IHSc))
C
C   Prepare the scaling factors
C
      Do 110 I = 0, NNRed-1
        W(IHSc+I) = Sqrt(W(IHSc+I))
  110 Continue
      N0 = 0
      N1 = NP
C  Here we go into the main loop
   10 Continue
      If (N1 .Le. NP .And. IFail .Eq. 0) Then
C  Go to the eigenvector space, scale and prepare the error vectors
        IShft = N0*NNRed
        Call MatMP1(1,0,1,NNRed,NNRed,NQ,N1-N0,UH,QQ(1,N0+1),
     $              W(IQQS+IShft))
        Call MatMP1(1,0,1,NNRed,NNRed,NQ,N1-N0,UH,FF(1,N0+1),
     $              W(IDQQS+IShft))
        Do 210 I = N0, N1-1
          IShft = I*NNRed
          Call AMul(NNRed,W(IQQS+IShft),W(IHSc),W(IQQS+IShft))
  210   Continue
        Call ANeg(NNeg,W(IHSc),W(IHSc))
        Do 220 I = N0, N1-1
          IShft = I*NNRed
          Call ADiv(NNRed,W(IDQQS+IShft),W(IHSc),W(IDQQS+IShft))
  220   Continue
        Call AAbs(NNRed,W(IHSc),W(IHSc))
        If (N1 .Eq. NP) IPrnt = Max(1,IPrint)
        Call GDIIS(IOut,IPrnt,NNRed,NNRed,N1,NUsed,N0,.True.,ICos,
     $    IChkC,W(IQQS),W(IDQQS),W(IA),W(IAInv),Scale,W(ICC),W(ICSave),
     $    W(IDQD),W(IDQR),W(IW),IFail)
        If (IFail .Eq. 0 .And. NUsed .Lt. N1) IFail = 1
        N0 = NUsed
        N1 = N0 + 1
C  Go back and check if we have something more to do
        Goto 10
      End If
C
C  We have the DIIS coefficients, so it is faster to calculate the DIIS
C  step in the original coordinate system than transform it
C
      Call MatMpy(NQ,NUsed,1,QQ,W(ICC),DQDIIS)
      Call ASub(NQ,DQDIIS,QQ,DQDIIS)
C
C  but it is better to scale back and transform the residuum.
C
      Call ADiv(NNRed,W(IDQR),W(IHSc),W(IDQR))
      Call MatMpy(NQ,NNRed,1,UH,W(IDQR),DQRes)
      OK = NUsed .Ge. 2
      If (IFail .Gt. 1) NP = NUsed
      Return
      End
*Deck DoMix
      Subroutine DoMix(IOut,IPrint,IRwPE,IMDS,IScMMF,NP,NMax,MaxGDI,
     $  NVarT,NVar,NVrMMT,NDimMM,NNeg,NEig,EH,UH,XX,FF,XXMM,FFMM,IStep,
     $  XNew,XQuad,W,OK,Switch,FncErr,IGDRFO,OKHess,DidBck,V,MDV)
      Implicit Real*8 (A-H,O-Z)
C
C     Parameters:
C     IOut   : Output channel
C     IPrint : Debug print level
C     IMDS   : MDS optimizer option flag
C              Method = Mod(IMDS,10)     - Reserved, not used
C     IScMMF : Handling of MM coordinates/forces; see FacMMM.
C     NP     : Number of saved variables and error vectors,
C     NVarT  : Dimension of the problem
C     NVar   : Number of internal coordinates
C     NVrMMT : Number of MM coordinates (cartesian and rigid blocks).
C     NDimMM : 1 or NVrMMT, for dimensioning.
C     NNeg   : Index of critical point to find
C     EH     : Eigenvalues of the Hessian
C     UH     : Eigenvectors of the Hessian
C     XX     : Saved internal coordinate vectors, last (or prefered) first
C     FF     : Saved internal coordinate forces, same order than in XX
C     XXMM   : Save MM coordinates
C     FFMM   : Save MM forces
C     XNew   : Output coordinates (XNew)
C     XQuad  : Output quadratic part of step
C     W      : Working space
C     OK     : DoMix could produce a resonable step
C
      Parameter (NDimIP=100)
      Dimension XX(NVar,*),FF(NVar,*),XNew(*),XQuad(*),W(*),EH(NVar),
     $  UH(NVar,NVar),IStep(*),XXMM(NDimMM,*),FFMM(NDimMM,*),V(MDV),
     $  IPoint(NDimIP)
      Logical OK, UpdLam, DidRFO, DidBck, RFOFai, OKHess
      Save Zero, One
      Data Zero/0.D0/, One/1.D0/
 1000 Format(' En-DIIS/RFO-DIIS IScMMF=',I9,' using points: ',20I5)
 1010 Format(' Point # ',I4,' is marked for removal')
 1020 Format(' Using NR instead of RFO step for point',I2,'.')
 1030 Format(' RFO failed for a point during GDIIS;',
     $  ' revert to using forces.')
 1040 Format(' Using forces rather than RFO steps in GEDIIS.')
C
C     Initialization
C     Do some memory allocation
C
      If(NDimIP.lt.MaxGDI) Call GauErr('NDimIP too small in DoMix.')
      Call SimMem(NP,NVarT,0,IDQQS,IA,IAInv,IB,IC,ICM,IHSc,IDQR,IW,IEnd)
      Call RwFGED(1,IRwPE,NMax,MaxGDI,IEn,INNEn,ICoeff,IQNDIS,IEnDIS,
     $  IMask,ICnvQM)
      IFF   = 1
      ICC   = IFF  + NP*NVarT
      Mask  = ICC  + NP*NVarT
      IMap  = Mask + NMax
      NFF   = IMap  + NMax
      NCC   = NFF   + NVarT*NP
      IVEnd = NCC   + NVarT*NP
      Call TstCor(IVEnd,MDV,'Optmz-DoMix')
C
C     Retrieving saved forces and coordinates; bad points are masked.
C
      FacMMF = FacMMM(IOut,IPrint,.True.,IScMMF,NVrMMT,NEig,EH)
      Call FileIO(2,-IRwPE,NMax,V(Mask),IMask)
      NUsed = 0
      Do 10 I = 0, (NMax-1)
        If(V(Mask+I).gt.Zero) then
          IndQM = NUsed*NVarT
          IndMM = IndQM + NVar
          Call AMove(NVar,FF(1,I+1),V(IFF+IndQM))
          Call AMove(NVar,XX(1,I+1),V(ICC+IndQM))
          Call AScale(NVrMMT,FacMMF,FFMM(1,I+1),V(IFF+IndMM))
          Call AMove(NVrMMT,XXMM(1,I+1),V(ICC+IndMM))
          NUsed = NUsed + 1
          IPoint(NUsed) = IStep(I+1)
          endIf
        If(NUsed.eq.NP) Goto 100
   10   Continue
  100 Write(IOut,1000) IScMMF, (IPoint(I),I=1,NUsed)
C
C     Check and adjust the eigenvalues of the Hessian.  Apply level
C     shift when the lowest eigenvalue is too small or negative.
C
      Call ASqrt(1,NVar,EH,W(IHSc))
      If(NNeg.gt.0) Call ANeg(NNeg,W(IHSc),W(IHSc))
C
C     Transform all saved geometric forces into current eigenvector space
C
      Call AClear(NVarT*NUsed,W(IDQQS))
      Call XGEMM(1,'T','N',NVar,NUsed,NVar,One,UH,NVar,V(IFF),NVarT,One,
     $  W(IDQQS),NVarT)
      Call AMove(NVarT*NUsed,W(IDQQS),V(NFF))
      Call XGEMM(1,'T','N',NVar,NUsed,NVar,One,UH,NVar,V(ICC),NVarT,
     $  Zero,V(NCC),NVarT)
C
C     Construct geometric error vectors -- Forces divided by
C     Sqrt(Hessian).
C
      If(IGDRFO.eq.3.or.(IGDRFO.eq.4.and..not.OKHess)) then
        RFOFai = .True.
        Write(IOut,1040)
      else
        RFOFai = .False.
        UpdLam = .True.
        Do 200 I = 0, (NUsed-1)
          IndQM = I*NVarT
          IndMM = IndQM + NVar
          Call AClear(NVar,W(IDQR))
          Call GDIRFO(IOut,IPrint,.True.,NVar,NNeg,NVar,W(IHSc),UH,
     $      W(IDQQS+IndQM),W(IDQR),RLam0,Rlamda,UpdLam,DidRFO)
          RFOFai = RFOFai.or..not.DidRFO
          If(DidRFO) then
            Call AMove(NVar,W(IDQR),W(IDQQS+IndQM))
          else
            Write(IOut,1020) I+1
            Call ADiv(NVar,W(IDQQS+IndQM),W(IHSc),W(IDQQS+IndQM))
            endIf
          Call AMove(NVrMMT,V(IFF+IndMM),W(IDQQS+IndMM))
          UpdLam = .False.
  200     Continue
        endIf
      Call AAbs(NVar,W(IHSc),W(IHSc))
C
C     If RFO failed for any step, possibly revert to using forces.
C
      If(RFOFai.and.IGDRFO.ge.2.and.IGDRFO.le.4) then
        If(IGDRFO.eq.2) then
          Write(IOut,1030)
          Call AClear(NVarT*NUsed,W(IDQQS))
          Call XGEMM(1,'T','N',NVar,NUsed,NVar,One,UH,NVar,V(IFF),NVarT,
     $      One,W(IDQQS),NVarT)
          endIf
        If(NVrMMT.gt.0) then
          Do 210 I = 0, (NUsed-1)
            IndQM = I*NVarT
            IndMM = IndQM + NVar
            Call AMove(NVrMMT,V(IFF+IndMM),W(IDQQS+IndMM))
  210       Continue
          endIf
        endIf
C
C     It is time to solve DIIS equation A*C = B
C     Note: The density part of the DIIS equation has already been
C           calculated and stored in IRwPE in l502
C
      NVRMS = NVarT
      NPIn = NUsed
C     Use transformed coordinates and forces if TS search is requested
      If(NNeg.eq.0) then
        Call MixDIS(IOut,IPrint,IRwPE,NVarT,NVRMS,NPIn,NUsed,NMax,
     $    MaxGDI,NNeg,W(IDQQS),W(IA),W(IAInv),W(IB),W(IC),W(ICM),
     $    W(IDQR),W(IW),V(ICC),V(IFF),V(IMap),Switch,ISkip,FncErr,OK,
     $    DidBck)
      else
        Call MixDIS(IOut,IPrint,IRwPE,NVarT,NVRMS,NPIn,NUsed,NMax,
     $    MaxGDI,NNeg,W(IDQQS),W(IA),W(IAInv),W(IB),W(IC),W(ICM),
     $    W(IDQR),W(IW),V(NCC),V(NFF),V(IMap),Switch,ISkip,FncErr,OK,
     $    DidBck)
        endIf
C
C     Take a geometry step using the DIIS coefficients
C     X' = Sum Ci*Xi
C     r' = Sum Ci*ri (r' is returned from MixDIS)
C     Xnew = X' + r' (This is done outside later)
C
C     Density step will be taken in l502
C
      If(OK) then
C       For the guided DIIS TS search method linear step computed in
C       Hessian eigenvector space and transformed back.
        If(NNeg.gt.0) then
          Call AClear(NVarT,V(NFF))
          Call MatMpy(NVarT,NUsed,1,V(NCC),W(ICM),V(NFF))
          Call ASub(NVarT,V(NFF),V(NCC),V(NFF))
          Call MatMpy(NVar,NVar,1,UH,V(NFF),XNew)
          Call AMove(NVrMMT,V(NFF+NVar),XNew(NVar+1))
        else
          Call MatMpy(NVarT,NUsed,1,V(ICC),W(ICM),XNew)
          Call ASub(NVarT,XNew,V(ICC),XNew)
          endIf
        Call ADiv(NVar,W(IDQR),W(IHSc),W(IDQR))
        Call MatMpy(NVar,NVar,1,UH,W(IDQR),XQuad)
        Call MatMpy(NVrMMT,NP,1,FFMM,W(ICM),XQuad(NVar+1))
        endIf
C
C     Remove a bad point
C
      If(ISkip.gt.1.and.(IPoint(1)-IPoint(NUsed)).le.2*MaxGDI) then
        Write(IOut,1010) ISkip
        NSkip = 0
        Do 300 I = 0, NMax-1
          If(V(Mask+I).gt.Zero) then
            If((NSkip+1).ne.ISkip) then
              NSkip = NSkip + 1
            else
              V(Mask+I) = -One
              Goto 400
              endIf
            endIf
  300     Continue
  400   Call FileIO(1,-IRwPE,NMax,V(Mask),IMask)
        Call FileIO(2,-IRwPE,1,EngTmp,IEn+ISkip-1)
        Call FileIO(2,-IRwPE,MaxGDI-ISkip,W(IW),IEn+ISkip)
        Call FileIO(1,-IRwPE,MaxGDI-ISkip,W(IW),IEn+ISkip-1)
        Call FileIO(2,-IRwPE,1,EngTmp,INNEn+ISkip-1)
        Call FileIO(2,-IRwPE,MaxGDI-ISkip,W(IW),INNEn+ISkip)
        Call FileIO(1,-IRwPE,MaxGDI-ISkip,W(IW),INNEn+ISkip-1)
        endIf
      Return
      End
*Deck DxGS
      Subroutine DxGS(IOut,NVar,P,T,Vec,Eigen,S,DQ,FToler,Freeze,Print)
      Implicit Real*8(A-H,O-Z)
C
C     For path relaxation.  Do a Gonzlez-Schlegel constrianed optimization
C     step.
C
      Logical Freeze,Print
      Dimension P(Nvar),T(Nvar),DQ(Nvar),Eigen(Nvar),Vec(NVar*NVar)
      Save BigEig,Zero,PT05,Facto1,PT5,XTol
      Data BigEig/1.0d3/,Zero/0.0d0/,PT05/0.05d0/,Facto1/0.05d0/,
     $  PT5/0.5d0/,XTol/1.0d-12/
 2000 Format(' Bracket for Lambda obtained in cycles=',I5)
 2010 Format(2X,F10.5)
 2023 Format(' New point DX updated by quad. search.')
 3030 Format(' Number of Inactive Modes detected = ',I4)
C
C     FIRST GUESS OF LAMBDA FOR MINIMUM OPTIMIZATION.
C     XL1=LOWER BOUNDARY OF LAMBDA
C     XL2=UPPER BOUNDARY FOR LAMBDA.
C
      NTRY = 12000
      IMAX = 12000
      Eigmin = BigEig
      IFreez = 0
      Do 10 I = 1, NVar
        If(Abs(P(I)).lt.FToler.and.Abs(T(I)).lt.FToler.and.Freeze) then
          Ifreez = IFreez + 1
          Eigen(I) = BigEig
          endIf
        If(Abs(P(I)).gt.FToler.and.Eigen(I).lt.Eigmin) Eigmin = Eigen(I)
   10   Continue
      If(Print) Write(IOut,3030) Ifreez
      XL1 = EigMin - Pt05*Abs(EigMin)
      XL2 = EigMin - Abs(EigMin)
      XF1 = HLamb(NVar,P,T,Eigen,XL1,S)
      XF2 = HLamb(NVar,P,T,Eigen,XL2,S)
      Do 20 I = 1, NTRY
        If(XF1*XF2.ge.Zero) then
          IF(Abs(XF1).LT.Abs(XF2)) then
            XL1 = XL1 + Facto1*(XL1-XL2)
            XF1 = HLamb(NVar,P,T,Eigen,XL1,S)
          else
            XL2 = XL2 + Facto1*(XL2-XL1)
            XF2 = HLamb(NVar,P,T,Eigen,XL2,S)
            endIf
          endIf
   20   Continue
C
C     Optimization of Lambda.  This optimization is carried out by using
C     a simple binary search  procedure.
C
      If(Print) Write(IOUT,2000) I
      FMid = HLamb(NVar,P,T,Eigen,ZERO,S)
      FMid = HLamb(NVar,P,T,Eigen,XL2,S)
      XF = HLamb(NVar,P,T,Eigen,XL1,S)
      If(XF.lt.zero) then
        XLamb = XL1
        DLamb = XL2 - XL1
      else
        XLamb = XL2
        DLamb = XL1 - XL2
        endIf
      Do 30 I = 1, IMax
        DLamb = Pt5*DLamb
        XLMID = XLamb + DLamb
        FMID = HLamb(NVar,P,T,Eigen,XLMID,S)
        If(FMID.LE.ZERO) XLamb = XLMID
        If(Abs(DLamb).LT.XTOL.OR.FMID.EQ.ZERO) goto 40
   30   Continue
C
C     Calculation of DX and updating of the points
C
   40 If(Print) then
        Write(IOUT,2023)
        Write(IOut,2010) XLamb
        endIf
      Do 60 I = 1, NVar
        DD = Zero
        Do 50 J = 1, NVar
          DDD = -(P(J)-XLamb*T(J))/(Eigen(J)-XLamb)
          IJ = NVar*(J-1)+I
   50     DD = DD + DDD*VEC(IJ)
   60   DQ(I) = DD
      Return
      End
*Deck DXLinr
      Subroutine DXLinr(ICase,INeg,F,XNew,DXRMS,DXMax,NVar,OK,XOut,ES,
     $  FS,FC1,FC2,XX,FF,DeltaX,HDOld,HDNew,FFMM,XXMM,NVarMX,NVarT,
     $  NDimMM,MaxL,IFCase,DoEle,FFMMQ,IDir0)
      Implicit Real*8(A-H,O-Z)
C
C     Fit a polynomial through the last two points and extrapolate
C     to the minimum along the line.  This step is referred to as the
C     linear search.  The quadratic search is then made starting from
C     this extrapolted minimum in subroutine DXQuad.  The output XOld
C     is 0.0 for the latest point (energy ES), -1.0 to step all the
C     way back to the old point (energy FS(1)).  INeg is the requested
C     curvature.  OK is .true. if a step was generated.
C
C     There are three cases in which this routine can perform the
C     linear search:
C
C     1.  INeg=0 and IPSav1 OR IPSav2=0.  In this case the search is for
C         a local minimum and there are only first derivatives available
C         at one or both of the points.  Routine LnStep is used to fit a
C         cubic or constrained quartic to locate the minimum.  Unless
C         there is only a local maximum between the two points, the
C         quartic will have exactly one minimum.  A minimum in the cubic
C         is only considered valid if it is "reasonably" close to the
C         search region.
C
C     2.  INeg=0 and IPSav1,IPSav2>0.  Here both first and second
C         derivatives are available at both points.  A general quintic
C         polynomical is fit.  If this does not have a reasonable
C         minimum, method 1 is used instead.
C
C     3.  INeg>0 and IPSav1,IPSav2>0.  Here the search is for a higher
C         order stationary point.  The linear search attempts to
C         minimize the norm of the gradient, and consequently second
C         derivatives must be available at both points.  A quartic fit
C         to the gradient norm analogous to 1 is performed.  This is
C         inherently inferior to the newton-raphson step and is only
C         performed if the previous step failed to reduce the norm of
C         of the gradient.
C
C     ICase is set to the relevent case.
C
C     When NVarMM > 1 , the MM region is included. In principle the
C     searches that use second order information can be used for this case
C     as well, but this has not yet been implemented.
C
      Logical OK, Prnt, Exit, Debug, HavDep, Large, FuInFF, DoEle
      Parameter (Debug=.False.)
      Dimension XNew(*), F(*), DeltaX(*), HDOld(*), HDNew(*),
     $  FS(*), FC1(*), FC2(*), XX(NVar,*), FF(NVar,*), FFMM(NVarT,3),
     $  XXMM(NDimMM,3), FFMMQ(NDimMM,3)
      Common /OptGrd/ ConvF,FMaxT,DXMaxT,RMax,RMin,RLim,EigMax,EigMin,
     $  DEPred,FSwtch,FncErr,GrdErr,FncCnv,XXIRC,NLevel,NStep,IStep,
     $  NMax,NP,Neg,Prnt,Exit,IPSav1,IPSav2,HavDEP,NGDIIS,IRisen,MaxRis,
     $  Large,ILevel,ISStep,IOffs,IMDS,MaxSv,MaxSv2,NVrMMT
      Common /IO/ In,IOut,IPunch
      Parameter (MaxLin=6,MxLin1=(MaxLin+1))
      Dimension XVal(2), YVAL(MaxLin), XMat(MaxLin,MaxLin),
     $  IS(MaxLin,4), Scr(MxLin1), A(MaxLin), DA(MaxLin), DDA(MaxLin),
     $  XLoc(MxLin1), XStart(MxLin1), Q(MaxLin,MaxLin),
     $  R(MaxLin,MaxLin), QR(MaxLin,MaxLin)
      Save Zero,One,XStart,XOuMax,XOuMx2
      Data Zero/0.0d0/,One/1.D0/,XOuMax/2.0d0/,XOuMx2/1.0d0/
     $  XStart/-2.0D0,-1.0D0,-0.75D0,-0.50D0,-0.25D0,0.0D0,1.0D0/
 1000 Format(' DXLinr:  EOld= ',1PD12.5,' ENew= ',1PD12.5,/,
     $       '          FOld= ',1PD12.5,' FNew= ',1PD12.5,/,
     $       '          GNOld=',1PD12.5,' GNNew=',1PD12.5)
 1010 Format(' DXLinr:  Minimizing gradient -- IPSav =',2I4)
 1020 Format(' DXLinr:  XVal = ',2(1PD12.5),/,
     $       '          YVal = ',3(1PD12.5),/,
     $       '                 ',3(1PD12.5))
 1030 Format(' DXLinr:  Quintic fit failed.')
 1040 Format(' DXLinr:  Quintic coeficients = ',3(1PD12.5),/,
     $       '                                ',3(1PD12.5),/,
     $       '                       values = ',3(1PD12.5),/,
     $       '                                ',3(1PD12.5),/,
     $       '                                ',3(1PD12.5))
 1050 Format(' DXLinr:  No Minimum found in quintic.')
 1060 Format(' DXLinr:  The Minimum in the quintic is at ',1pd12.5)
 1070 Format(' DXLinr:  Using constrained quartic.')
 1081 Format(' DXLinr:  FOld = ',4(1PD12.5),/,
     $      (('                 ',4(1PD12.5))))
 1082 Format(' DXLinr:  FNew = ',4(1PD12.5),/,
     $      (('                 ',4(1PD12.5))))
 1083 Format(' DXLinr:  DelX = ',4(1PD12.5),/,
     $      (('                 ',4(1PD12.5))))
 1084 Format(' DXLinr:  HDO  = ',4(1PD12.5),/,
     $      (('                 ',4(1PD12.5))))
 1085 Format(' DXLinr:  HDN  = ',4(1PD12.5),/,
     $      (('                 ',4(1PD12.5))))
 1090 Format(' DXLinr:  Hessian for point ',I3)
C
C     Evaluate first derivative components old and new along direction
C     of motion.
C
      OK = .False.
      ICase = 0
      DxRMS = Zero
      DXMax = Zero
      NVar1 = NVar+1
      NVarMM = NVarMX
      If(NVarMM.eq.1) NVarMM = 0
      Call ASub(NVar,XX(1,2),XX(1,1),DeltaX)
      Call ASub(NVarMM,XXMM(1,2),XXMM(1,1),DeltaX(NVar1))
      FuInFF=(IFCase.eq.1.and.IDir0.lt.2).or.IFCase.eq.3.or.IFCase.eq.0
      If(FuInFF) then
C       Full forces are in FF
        FOld = SProd(NVar,FF(1,2),DeltaX) +
     $   SProd(NVarMM,FFMM(NVar1,2),DeltaX(NVar1))
        FNew = SProd(NVar,FF(1,1),DeltaX) +
     $   SProd(NVarMM,FFMM(NVar1,1),DeltaX(NVar1))
        GNOld = SProd(NVar,FF(1,2),FF(1,2)) +
     $   SProd(NVarMM,FFMM(NVar1,2),FFMM(NVar1,2))
        GNNew = SProd(NVar,FF(1,1),FF(1,1)) +
     $   SProd(NVarMM,FFMM(NVar1,1),FFMM(NVar1,1))
        If(IPSav2.ne.0) then
          Call MultLV(.True.,1,0,NVar,NVar,FC2,DeltaX,HDOld)
          Call MultLV(.True.,1,0,NVar,NVar,FC1,DeltaX,HDNew)
          endIf
      else
C       Full forces are in FFMM
        FOld = SProd(NVarT,FFMM(1,2),DeltaX)
        FNew = SProd(NVarT,FFMM(1,1),DeltaX)
        GNOld = SProd(NVarT,FFMM(1,2),FFMM(1,2))
        GNNew = SProd(NVarT,FFMM(1,1),FFMM(1,1))
        endif
      If(Debug) then
        Write(IOut,1000) FS(1), ES, FOld, FNew, GNOld, GNNew
        Write(IOut,1081) (FF(I,2),I=1,NVar)
        Write(IOut,1082) (FF(I,1),I=1,NVar)
        Write(IOut,1083) (DeltaX(I),I=1,NVar)
        If(IPSav2.ne.0.and.NVarMM.eq.0) then
          Write(IOut,1090) IPSav2
          Write(IOut,1084) (HDOld(I),I=1,NVar)
          Write(IOut,1085) (HDNew(I),I=1,NVar)
          endIf
        If(IPSav1.ne.0) Write(IOut,1090) IPSav1
        endIf
C
C     Case 3 -- minimize the norm of the gradient.
C
      If(INeg.gt.0) then
        If(Debug) Write(IOut,1010) IPSav1, IPSav2
        If(IPSav1.ne.0.and.IPSav2.ne.0) then
          If(NVarMM.gt.0) Call GauErr('DXLinr#3: Not yet implemented')
          DGNOld = SProd(NVar,HDOld,FF(1,2))
          DGNNew = SProd(NVar,HDNew,FF(1,1))
          Call LnStep(IOut,Debug,.True.,GNOld,GNNew,DGNOld,DGNNew,XOut,
     $      RJunk,OK)
          If(OK) ICase = 3
          endIf
C
C     Case 2 -- quintic fit.  if this fails, do a constrained quartic.
C
      else if(IPSav1.ne.0.and.IPSav2.ne.0.and.NVarMM.eq.0) then
        XVal(1) = -One
        XVal(2) = Zero
        YVal(1) = FS(1) * FncCnv
        YVal(2) = ES * FncCnv
        YVal(3) = FOld
        YVal(4) = FNew
        YVal(5) = SProd(NVar,HDOld,DeltaX)
        YVal(6) = SProd(NVar,HDNew,DeltaX)
        If(Debug) then
          Write(IOut,1020) XVal, YVal
          IPrint = 1
        else
          IPrint = 0
          endIf
        Call PolFit(IOut,IPrint,2,3,XVal,YVal,XMat,IS,Scr,A,OK,Q,R,QR)
        If(Debug.and..not.OK) Write(IOut,1030)
        If(OK) then
          Call PolMin(5,A,7,XStart,DA,DDA,Scr,XLoc,XOut,RJunk,OK)
          If(.not.OK) XOut = Zero
          FXOld  = FX(5,-One,A)
          FXNew  = FX(5,Zero,A)
          FXMin  = FX(5,XOut,A)
          DFXOld = FX(4,-One,DA)
          DFXNew = FX(4,Zero,DA)
          DFXMin = FX(4,XOut,DA)
          DDFXO  = FX(3,-One,DDA)
          DDFXN  = FX(3,Zero,DDA)
          DDFXM  = FX(3,XOut,DDA)
          If(Debug) then
            Write(IOut,1040) A, FXOld, FXNew, FXMin, DFXOld, DFXNew,
     $        DFXMin, DDFXO, DDFXN, DDFXM
            If(OK) then
              Write(IOut,1060) XOut
            else
              Write(IOut,1050)
              endIf
            endIf
          OK = OK.and.XOut.ge.-One.and.(ES.lt.FS(1).or.XOut.le.Zero)
          If(OK) ICase = 2
          endIf
        If(.not.OK) then
          If(Debug) Write(IOut,1070)
          Call LnStep(IOut,Debug,.True.,FS(1)*FncCnv,ES*FncCnv,FOld,
     $      FNew,XOut,RJunk,OK)
          If(OK) ICase = 1
          endIf
C
C     Case 1 -- constrained quartic.
C
      else
        If(Debug) Write(IOut,1070)
        Call LnStep(IOut,Debug,.True.,FS(1)*FncCnv,ES*FncCnv,FOld,FNew,
     $    XOut,RJunk,OK)
        If(OK) ICase = 1
        endIf
C
C     The step XOut has been determined.  Compute XNew.
C
      If(OK.or.MaxL.gt.0) then
        XOuMx = XOuMax
        If(NVarMM.gt.0) XOuMx = XOuMx2
        If(MaxL.gt.0) XOuMx = GFloat(MaxL)
        If(.not.OK.and.MaxL.gt.0) then
          XOut = XOuMx
          OK = .True.
          ICase = 5
          endif
        If(NVar.gt.1) XOut = Min(XOut,XOuMx)
        If(NVarMM.gt.0) then
          NVarX = NVarT
        else
          NVarX = NVar
          endif
        RX = SProd(NVarX,DeltaX,DeltaX)
        Call ACASB(NVarX,XNew,DeltaX,XNew,-XOut)
        DXRMS = XOut*XOut*SProd(NVar,DeltaX,DeltaX)
        DXMax = Abs(XOut)*ArrMax(DeltaX,NVar,.True.)
        If(IFCase.eq.0) then
C         Standard, no QM/MM
          Call ACASB(NVar,F,FF(1,2),F,-XOut)
          Call ACASB(NVar,F,FF(1,1),F,XOut)
        else if(IFCase.ne.2.or.DoEle) then
C         Need full forces in F
          If(FuInFF) then
            Call AScale(NVar,XOut+One,FF,F)
            Call ACASB(NVar,F,FF(1,2),F,-XOut)
            Call AScale(NVarMM,XOut+One,FFMM(NVar1,1),F(NVarT))
            Call ACASB(NVarMM,F(NVarT),FFMM(NVarT,2),F(NVarT),-XOut)
          else
            Call AScale(NVarT,XOut+One,FFMM,F)
            Call ACASB(NVarT,F,FFMM(1,2),F,-XOut)
            endif
        else
C         Need QM contribution in F
          If(FuInFF) then
            Call AScale(NVar,XOut+One,FFMM,F)
            Call ACASB(NVar,F,FFMM(1,2),F,-XOut)
          else
            Call AScale(NVar,XOut+One,FF,F)
            Call ACASB(NVar,F,FF(1,2),F,-XOut)
            endif
          Call AClear(NVarMM,F(NVarT))
          If(DoEle) then
            Call AScale(NVarMM,XOut+One,FFMMQ,F(NVar1))
            Call ACASB(NVarMM,F(NVar1),FFMMQ(1,2),F(NVar1),-XOut)
            endif
          endif
        FNew = SProd(NVarX,F,DeltaX)/RX
        If(IFCase.ne.2.or.DoEle) Call ACASB(NVarX,F,DeltaX,F,-FNew)
        DXRMS = Sqrt(DXRMS/NVar)
        endIf
      Return
      End
*Deck DXLOnm
      Subroutine DXLOnm(IOut,IPrint,Rises,NVar,NVrMMT,NDimMM,NAtoms,
     $  NAtNew,C,XX,FF,XXMM,XNew,NTBond,NTAng,NTDih,IZRed,BMat,IBMAt,
     $  Scr,Cnstr1,IFlags,IWght,NRF,COrig,MicOpt,COM,DT,T,ICnt,RotQM,
     $  IAtTyp,AtChMM,IRwMM,IMMMOD,ENew,EOld,NVarT,FFMM,XOut,NTrR1,Ok,
     $  DXRMS,DXMax,F,IntHes,IFCase,DoEle,FFMMQ,IMMCRS,NAtmNN,IndexQ,
     $  IDir0,XName,IApBGI,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     NVar        Number of redundant internals
C     NVrMMT      Number of MM variables (cartesian and rigid  fragments)
C     NAtoms      Total number of atoms
C     NAtNew      Number of atoms in model system
C     C           Cartesian coordinates of current point
C     XX(1,1)     Redundant internals of current point
C     XX(1,2)     Redundant internals of previous point
C     XXMM(1,1)   MM variables of current point
C     XXMM(1,2)   MM variables of previous point
C     FF(1,1)     Redundant internal forces (QM contrib only) current point
C     FF(1,2)     Redundant internal forces (QM contrib only) previous point
C
      Parameter (MaxLM=4)
      Logical Rises, RotQM, OK, IntHes, Swapped, DoEle, DoPrnt
      Dimension XNew(*), XX(NVar,*), FF(NVar,*), XXMM(NDimMM,3), C(3,*),
     $  Cnstr1(*), IFlags(*), BMat(*), IBMat(*), Scr(*), COrig(3,*),
     $  MicOpt(*), COM(*), DT(*), T(*), ICnt(*), Coors(2), A(4), F(*),
     $  Values(2,2), FFMM(NVarT,*), JJ(1), Frc(2), FFMMQ(NDimMM,3),
     $  AtChMM(*), IAtTyp(*), IZRed(*), XName(*), V(*)
      Save Zero, One, OneP5, Two, Three, IStep, JJ, XBackM
      Data Zero/0.0d0/, One/1.D0/, OneP5/1.5d0/, Two/2.0d0/,
     $  Three/3.0d0/, IStep/10/, JJ/1/, XBackM/-0.9d0/
 1000 Format(' Generate structures for explicit line search')
 1010 Format(7X,'Point    E(MM,Real)   E(MM,Model)   E(QM,Model)',6X,
     $  'E(Total)')
 1020 Format(I4,F8.3,4F14.8)
 1030 Format(' Generation of coordinates for QM/MM linear search',
     $  ' failed.')
 1040 Format(' Explicit MO:MM linear search failed. Do standard')
 1050 Format(' Fit of polynomial for QM energy failed.')
 1060 Format(' Fit of polynomial for integrated energy failed.')
 1070 Format(' XMin=',1PD10.2,' out of range ',A,'MinMax= ',1PD10.2,
     $  ' to ',1PD10.2)
C
C     If the energy went up, we only want to back up; if it went down, we
C     don't want to back up more than to the original point.
C
      DoPrnt = IPrint.gt.0.or.IPSave(0).eq.0
      If(NVarT.ne.(NVar+NVrMMT))
     $  Call GauErr('Thom is an idiot #1 in DXLOnm.')
      If(Rises) then
        IOutMn = -IGFix(OneP5*GFloat(IStep))
        IOutMx = 0
      else
        IOutMn = -IGFix(One*GFloat(IStep))
        IOutMx = IGFix(Two*GFloat(IStep))
        endIf
      NParts = IOutMx - IOutMn
      Step = One / GFloat(IStep)
      OutMin = GFloat(IOutMn) * Step
      OutMax = GFloat(IOutMx) * Step
      OutCur  = Zero
      NPoint  = NParts + 1
      NAt3    = 3*NAtoms
      III     = Max(3*NAtNew,NVar,4)
      NVarT   = NVar + NVrMMT
      NVarMM  = NVrMMT - 6*NRF
      IXNew   = 1
      IStrct  = IXNew + NVarT
      IGrad   = IStrct + NPoint*NAt3
      ICO     = IGrad + NPoint*NAt3
      ICors   = ICO + NPoint*NAt3
      IETot   = ICors + NPoint
      ITemp   = IETot + NPoint
      IScr    = ITemp + NAt3
      ICOTmp  = IScr + Max(NAt3,III**2)
      IVect   = ICOTmp + NAt3
      IFCT    = IVect + Max(NAt3,III**2)
      IG      = IFCT + 2*III**2
      IDelta  = IG  + III**2
      IHave   = IDelta + NVarT
      IGoal   = IHave + NVarT
      IEF     = IGoal + Max(2*NVarT,((III+1)*III)/2)
      Call LoadEF(MaxLM,MaxLEF,V(IEF))
      LEFild = ((MaxLEF+1)*(MaxLEF+2)*(MaxLEF+3))/6
      IPar   = IEF + LEFild
      LenPar = ITqry(IRwMM)
      LenPrM = ITqry(IMMMOD)
      IDat   = IPar + LenPar
      IQ1    = IDat + LenPrM
      IQ2    = IQ1 + Max(0,NAtmNN)
      IEnd1  = IQ2 + Max(0,NAtmNN) - 1
      IV1    = IEnd1 + 1
      MDV1   = MDV - IEnd1
      Call TstCor(IEnd1,MDV,'DXLOnm')
      Call FileIO(2,-IRwMM,LenPar,V(IPar),0)
      Call FileIO(2,-IMMMOD,LenPrM,V(IDat),0)
      If(DoEle) Call FileIO(2,-IMMCRS,NAtmNN,V(IQ2),IndexQ)
C
      Call AMove(NVarT,XNew,V(IXNew))
      Call AMove(NAt3,C,V(ITemp))
      Call ASub(NVar,XX(1,1),XX(1,2),V(IDelta))
      Call ASub(NVrMMT,XXMM(1,1),XXMM(1,2),V(IDelta+NVar))
C
C     Set up polynomial for QM energy
      If((IFCase.eq.1.and.IDir0.le.0).or.IFCase.eq.3) then
        FNew = -SProd(NVar,FFMM(1,1),V(IDelta))
        FOld = -SProd(NVar,FFMM(1,2),V(IDelta))
      else
        FNew = -SProd(NVar,FF(1,1),V(IDelta))
        FOld = -SProd(NVar,FF(1,2),V(IDelta))
        endIf
      If(DoEle) then
        FNew = FNew - SProd(NVrMMT,FFMMQ(1,1),V(IDelta+NVar))
        FOld = FOld - SProd(NVrMMT,FFMMQ(1,2),V(IDelta+NVar))
        endIf
      Coors(1)    = -One
      Coors(2)    = Zero
      Values(1,1) = EOld
      Values(1,2) = FOld
      Values(2,1) = ENew
      Values(2,2) = FNew
      Call PolFit(IOut,IPrint,2,2,Coors,Values,V(IVect),V(IFCT),
     $  V(ICOTmp),A,OK,V(IFCT),V(IStrct),V(IG))
      If(.not.OK) then
        If(DoPrnt) Write(IOut,1050)
        Goto 9999
        endIf
      If(IPrint.ge.0) Write(IOut,1000)
C
C     Generate the series of QM coordinates.
C     Always want to use the cartesians of an adjecent point, so juggle
C     a bit with the order in which things are done.
C     What is the point closest to OutCur (for which we have Cartesians)?
      ICur = (OutCur - OutMin + Step + Step/Two) / Step
      Do 10 I = 1, NPoint
        J = I
        If(I.le.ICur) J = (ICur - I) + 1
        If(J.gt.ICur) then
          Call AMove(NAt3,V(IStrct+(J-2)*NAt3),V(ITemp))
        else if(J.lt.ICur) then
          Call AMove(NAt3,V(IStrct+J*NAt3),V(ITemp))
          endIf
        Call AMove(NAt3,V(ITemp),V(IStrct+(J-1)*NAt3))
        Call RedVal(IOut,IPrint,.True.,V(ITemp),IZRed,V(IHave),NVar,
     $    NTBond,NTAng,NTDih,V(IV1),MDV1)
        OutDo = OutMin + (J-1)*Step
        V(ICors+J-1) = OutDo
        Call ACASB(NVar,XX(1,1),V(IDelta),V(IGoal),OutDo)
        Call ASub(NVar,V(IGoal),V(IHave),XNew)
        Call FixDQ(IOut,IPrint,0,NVar,NTBond,NTAng,NTDih,XNew,IZRed)
        IPrtm1 = Max(IPrint-2,-1)
        Call RedCar(IOut,IPrtm1,0,NVar,NTBond,NTAng,NTDih,IZRed,
     $    V(IHave),NAtNew,V(ITemp),XNew,V(IStrct+(J-1)*NAt3),BMat,IBMat,
     $    V(IVect),V(IFCT),Scr,Cnstr1,IFlags,.True.,IFail,IWght,V(IV1),
     $    MDV1)
        OK = IFail.eq.0
        If(.not.OK) then
          Write(IOut,1030)
          Return
          endIf
   10   Continue
C
C     Get MM coordinates, combined energies, and forces
      Call AMove(NAt3,COrig,V(ICOTmp))
      If(IPrint.ge.0) Write(IOut,1010)
      Do 20 I = 1, NPoint
        X = V(ICors+I-1)
C       Coordinates of MM region (and possible rotate QM region)
        Call AScale(NVrMMT,X,V(IDelta+NVar),XNew(1+NVar))
        Call AMove(3*(NAtoms-NAtNew),C(1,NAtNew+1),
     $    V(IStrct+(I-1)*NAt3+3*NAtNew))
        Call UpCOrg(V(IStrct+(I-1)*NAt3),V(ITemp),V(ICOTmp),NAtoms,
     $    NAtNew,MicOpt,XXMM,COM,DT,ICnt,T,NRF,NVar,NVarMM,XNew,0,RotQM,
     $    .False.)
        Call AMove(NAt3,V(ICOTmp),V(ICO+(I-1)*NAt3))
        Call AMove(NAt3,V(IStrct+(I-1)*NAt3),V(ITemp))
C       Energies and Cartesian forces
        Call EGHMM(IOut,IPrint,0,1,0,Zero,NAtoms,IAtTyp,AtChMM,V(ITemp),
     $    MaxLEF,V(IEF),V(IPar),LenPar,MicOpt,ERL,V(IGrad+(I-1)*NAt3),
     $    V(IHave),V(IGoal),V(IHave),V(IV1),MDV1)
        EML = Zero
        Call EGHMod(IOut,IPrint,0,Zero,NAtoms,MicOpt,V(ITemp),EML,
     $    V(IScr),V(IGoal),V(IDat),V(IDat),V(IV1),MDV1)
        Call ASub(NAt3,V(IGrad+(I-1)*NAt3),V(IScr),V(IGrad+(I-1)*NAt3))
        If(DoEle) then
C         Take out self-interaction of charge cloud.
          ESelf = Zero
          Call AMove(NAtmNN,V(IDat+InToWP(7)),V(IQ1))
          Call AMove(NAtmNN,V(IQ2),V(IDat+InToWP(7)))
          Call EGHMod(IOut,IPrint,3,Zero,NAtoms,MicOpt,V(ITemp),ESelf,
     $      V(IScr),V(IGoal),V(IDat),V(IDat),V(IV1),MDV1)
          Call AMove(NAtmNN,V(IQ1),V(IDat+InToWP(7)))
          Call AAdd(NAt3,V(IGrad+(I-1)*NAt3),V(IScr),
     $      V(IGrad+(I-1)*NAt3))
          EML = EML - ESelf
          endIf
        EQM = A(1)*(X**3) + A(2)*(X**2) + A(3)*X + A(4)
        ETot = EQM + ERL - EML
        V(IETot-1+I) = ETot
        If(I.eq.1) EMin = ETot + One
        If(ETot.lt.EMin) IMin = I
        If(IPrint.ge.0) Write(IOut,1020) I,X,ERL,EML,EQM,ETOT
   20   EMin = V(IETot-1+IMin)
C
C     Determine the minimum between the lowest point and lowest adjacent point.
      If(IMin.eq.1) then
        II = 1
      else if(IMin.eq.NPoint) then
        II = NPoint - 1
      else if(V(IETot-2+IMin).lt.V(IETot+IMin)) then
        II = IMin - 1
      else
        II = IMin
        endIf
      Do 40 J = II, (II+1)
        IP = 0
        IQ = NVar
        Do 30 I = 1, NAtoms
          If(MicOpt(I).lt.0) then
            Call AMove(3,V(IStrct+(J-1)*NAt3+(I-1)*3),V(ITemp+IP))
            Call AMove(3,V(IGrad+(J-1)*NAt3+(I-1)*3),V(IScr+IP))
            IP = IP + 3
          else if(MicOpt(I).eq.0) then
            Call AMove(3,V(IGrad+(J-1)*NAt3+(I-1)*3),V(IGoal+(J-II)*
     $        NVarT+IQ))
            IQ = IQ + 3
            endIf
   30     Continue
        Call RedVal(IOut,IPrint,.True.,V(ITemp),IZRed,V(IHave),NVar,
     $    NTBond,NTAng,NTDih,V(IV1),MDV1)
        IGISav = 0
        IAprBG = (IApBGI/10)*10 + 4
        Call Red2BG(IOut,IPrint,IAprBG,NAtNew,NTrR1,IGISav,0,V(ITemp),
     $    NVar,NTBond,NTAng,NTDih,IZRed,V(IHave),NRank,BMat,IBMat,V(IG),
     $    V(IVect),V(IFCT+(J-II)*III**2),Scr,XName,V(IV1),MDV1)
        Call TranF(.True.,.False.,NVar,NAtNew,JJ,V(IScr),
     $    V(IGoal+(J-II)*NVarT),IBMat,BMat,V(IG),JJ)
        X = V(ICors+J-1)
        Call AScale(NRF*6,X,V(IDelta+NVar+NVarMM),XNew)
        Call AAdd(NRF*6,XXMM(1+NVarMM,1),XNew,XNew)
        Call RigidX(NAtoms,V(ICO+(J-1)*NAt3),XNew,COM,MicOpt,DT,ICnt,T)
        Call RigCnv(.True.,.False.,NAtoms,MicOpt,T,V(ICO+(J-1)*NAt3),
     $    V(IGrad+(J-1)*NAt3),COM,V(IGoal+(J-II)*NVarT+NVar+NVarMM))
        FQM = Three*A(1)*(X**2) + Two*A(2)*X + A(3)
        Frc(1+J-II) = SProd(NVarT,V(IGoal+(J-II)*NVarT),V(IDelta))*Step+
     $    FQM*Step
   40   Continue
C
C     May need to swap because LnStep depends on order of points.
      Swapped = V(IETot-1+II).lt.V(IETot+II)
      If(Swapped) then
        Call Swap(V(IETot-1+II),V(IETot+II))
        Call Swap(Frc(1),Frc(2))
        Call ANeg(2,Frc,Frc)
        endIf
      Call LnStep(IOut,.False.,.True.,V(IETot-1+II),V(IETot+II),Frc(1),
     $  Frc(2),XOut1,V(IScr),OK)
      If(.not.OK) then
        Write(IOut,1060)
        Goto 9999
        endIf
      If(Swapped) XOut1 = -XOut1 - One
      If(II.le.2) then
        XMin1 = -One
      else
        XMin1 = -Two
        endIf
      If(II.eq.NPoint) then
        XMax1 = Zero
      else
        XMax1 = One
        endIf
      OK = OK.and.XOut1.ge.XMin1.and.XOut1.le.XMax1
      If(.not.OK) then
        Write(IOut,1070) XOut1, 'X', XMin1, XMax1
        Goto 9999
        endIf
      XOut = V(ICors+II) + XOut1*Step
      OK = OK.and.XOut.gt.OutMin.and.XOut.lt.OutMax
      If(Rises) OK = OK.and.XOut.ge.XBackM
      If(.not.OK) then
        Write(IOut,1070) XOut, 'O', OutMin, OutMax
        Goto 9999
        endIf
      DXRMS = XOut*XOut*SProd(NVar,V(IDelta),V(IDelta))
      DXMax = Abs(XOut)*ArrMax(V(IDelta),NVar,.True.)
      NVar1 = NVar+1
      If(IFCase.eq.1.or.IFCase.eq.3.or.DoEle) then
C       F must end up with full forces
        If((IFCase.eq.1.and.IDir0.le.0).or.IFCase.eq.3) then
          Call AScale(NVar,XOut+One,FFMM,F)
          Call ACASB(NVar,F,FFMM(1,2),F,-XOut)
        else
          Call AScale(NVar,XOut+One,FF,F)
          Call ACASB(NVar,F,FF(1,2),F,-XOut)
          endif
        Call AClear(NVrMMT,F(NVar1))
        If(DoEle) then
          Call AScale(NVrMMT,XOut+One,FFMMQ,F(NVar1))
          Call ACASB(NVrMMT,F(NVar1),FFMMQ(1,2),F(NVar1),-XOut)
          endif
        Call ACASB(NVarT,F,V(IGoal+NVarT),F,-(One+XOut1))
        Call ACASB(NVarT,F,V(IGoal),F,XOut1)
        Call AMove(NVar,F,V(IScr))
        Call MultLV(.True.,1,0,NVar,NVar,V(IFCT),V(IScr),F)
      else
C       QuadMac: F has QM contribution only.
        Call ACASB(NVar,F,FF(1,2),F,-XOut)
        Call ACASB(NVar,F,FF(1,1),F,XOut)
        endIf
      DXRMS = Sqrt(DXRMS/NVar)
      Call ACASB(NVarT,V(IXNew),V(IDelta),XNew,XOut)
      Return
 9999 Write(IOut,1040)
      OK = .False.
      XOut = Zero
      Call AMove(NVarT,V(IXNew),XNew)
      Return
      End
*Deck DxPath
      Subroutine DxPath(IOut,IPrint,NVar,NPath,NPStar,NStop,N1,NVtt,
     $  NTBond,NTAng,NTDih,PStr,IDiag,Q0,Q,Q1,DQ,Proj,Grad,FC,Vec,Eigen,
     $  Scr,SCr5,IFlags,DXMaxT,T1,T2,ConvF,ConVX,Weight,Fact,PathL,TTS,
     $  TOld,UseSym)
      Implicit Real*8(A-H,O-Z)
C
C Path relaxation routine
C
C     NSeg : number of points in the segment considered
C     NPStar : Uphill point
C     NStop : Downhill point
C     N1=+/- 1 depending on the ordering of points
C     NVTT,NTBond,NTAng,NTDih : Dimensions.
C     PStr : End point coordinates
C     IDiag : Scratch index vector
C     Q0 : On entry Coordinates for which gradient have been computed
C          On exit Relaxed coordinates.
C     Q,Q1,DQ : Scratch
C     Proj : contains the redundancy projector for each point.
C     Grad : Actual Gradient
C     FC : Force constant matrix.
C     Vec,Eigen,Scr,Scr5 : Scratch
C     DXMaxT : Trust Radius
C     SS : Step size for Gonzales-Schlegel Step.
C     T1,T2,ConVF,ConVX : convergence criteria
C     Weight : Weight for micro iteration step size
C     DoEuler: Whether to use the tangent to the arc of circle
C             or not for the last point.
C
      Dimension Q0(Nvar,NPath),Q(NVar,NPath),Scr(NVar)
      Dimension DQ(NVar,NPath),Proj(NVTT,NPath),Q1(NVar,NPath)
      Dimension Vec(NVar*NVar,NPath),Eigen(NVar,NPath)
      Dimension Scr5(NVar,8),PStr(NVar),IDiag(NPath)
      Dimension Grad(Nvar,NPath),FC(NVTT,NPath),IFlags(Nvar)
      DImension PathL(Npath),TTS(Nvar),TOld(nvar),XX(1)
      Logical conv,DoEuler,DidIRC,Print,USeSym,DoLin
      Save Zero,One,two,Toler,FToler,Osc,XX
      Data Zero /0.0d0/,One /1.0d0/,Two /2.0d0/,Toler/1.0d-06/,
     $  FToler/1.0d-06/,Osc/-0.93969d0/,XX/0.d0/
  999 Format(' Micro-iterations ')
 1062 Format(1x,I2,' Max FP=',F8.5,' RMS FP=',F8.5,
     $' Max FA=',F8.5,' RMS FA=',F8.5)
 1063 Format(1x,I2,' Max QP=',F8.5,' RMS QP=',F8.5,
     $' Max QA=',F8.5,' RMS QA=',F8.5)
 2001 Format(2x,'  X Old  ',1x,'   DX    ',1x,'  X New  ')
 3000 Format(1x,I4,' Iterations in Redistribution')
 3010 Format(1x,' GS Constrained Optimization Radius = ',F10.5)
 3040 Format(1x,' Oscillation detected for points ',I3,' and ',I3)
 3050 Format(1x,' Point ',I3,' travelled too far ')
 3060 Format(1x,' Relaxed Path After ',I4,' Iterations')
C
      Print = .false.
      DoLin = .false.
C     Make the convergence criteria a bit more stringent
      FMax = T1*Fact
      FRMS = ConvF*Fact
      DQMax = T2*Fact
      DQRMS = ConvX*Fact
      Write(IOut,999)
      NMin = Min(NPStar,NStop)
      NMax = Max(NPStar,NStop)
      NSeg = NMax-NMin
      NIter = 0
      R = DXMaxT/(Weight*DQMax)
      MaxIter = 1 + IGFix(R)
C
C     IDiag = 1/0 whether the eigenvalues and eigevectors are available
C     The diagonalizations were done in diagf3 and do not need to repeated
C     here.
    1 NIter=NIter+1
C     Redistribute points
      ICount = 0
    2 Dev = Zero
      ICount = ICount+1
      Call AMove(NVar*(NSeg+1),Q(1,NMin),Q1(1,NMin))
      Do 10 NPt = (NPStar+N1), NStop, N1
        Call AClear(NVar,DQ(1,NPt))
C
C       Q contains the redistributed coordinates
C       Q0 contains the coordinates for which the gradient was computed
C       Redistribute point
        DoEuler = PathL(Npt).lt.Zero
        If(NPt.ne.NStop) then
          Call CnssT4(2,NVar,IFlags,Proj(1,NPt),Q(1,NPt),Q1(1,NPt+N1),
     $      Q1(1,NPt-N1),Scr,Grdtmp,Scr5,NTBond,NTAng,NTDih,RP,RR,RPR,
     $      ArcLnR,ArcLnP,1,0,Doeuler)
        else
          Call CnssT4(2,NVar,IFlags,Proj(1,NPt),Q(1,NPt),PStr,
     $      Q1(1,NPt-N1),Scr,grdtmp,Scr5,NTBond,NTAng,NTDih,RP,RR,RPR,
     $      ArcLnR,ArcLnP,1,0,DoEuler)
          endIf
        Call FixDih(NTDih,Q(NTBond+NTAng+1,NPt))
        RP = Abs(RP-RR)
        If(RP.gt.Dev) dev = RP
   10   Continue
C     Redistribute until maximum deviation <Toler
      If(dev.gt.Toler.and.ICount.le.100) Goto 2
      If(Print) Write(Iout,3000) ICount
      If(NIter.eq.MaxIter) Goto 6
      IConv = 0
      Do 30 NPt = (NPStar+N1), NStop, N1
        Call AMove(NVar,Q0(1,NPt),Scr5(1,5))
        Call AMove(NVTT,FC(1,NPt),Vec(1,NPStar))
C       Get Gradient at Q and put it in Scr
        Call GetGrd(IOut,IPrint,NVar,Scr5(1,5),Grad(1,NPt),
     $    Vec(1,NPSTar),Proj(1,NPt),Q(1,NPt),Scr,Scr5,DotG,NTBond,NTAng,
     $    NTDih)
C       Get QST Tangent and put it in DQ(1,NPStar)
        Doeuler = PathL(Npt).lt.Zero
        If(NPt.ne.NStop) then
          Call CnssT4(2,NVar,IFlags,Proj(1,NPt),Q(1,NPt),Q(1,NPt+N1),
     $      Q(1,NPt-N1),Scr,Grdtmp,Scr5,NTBond,NTAng,NTDih,RP,RR,RPR,
     $      ArcLnR,ArcLnP,0,0,Doeuler)
        else
          Call CnssT4(2,NVar,IFlags,Proj(1,NPt),Q(1,NPt),PStr,
     $      Q(1,NPt-N1),Scr,Grdtmp,Scr5,NTBond,NTAng,NTDih,RP,RR,RPR,
     $      ArcLnR,ArcLnP,0,0,DoEuler)
          endIf
        Call AMove(NVar,Scr5(1,4),DQ(1,NPStar))
        If(NPt.eq.(NPStar+N1)) Call AMove(Nvar,TTS,TOld)
        Call GetTan(IOut,0,NVar,.true.,Scr5(1,2),TOld,Scr,Scr5(1,4),S,
     $    Proj(1,NPt),Scr5(1,5),ArcL)
        S = S*two
        S2 = S/two
        Call AScale(NVar,DotG,Scr,Scr)
        Call GetFPA(NVar,Scr,Scr5(1,4),Scr5(1,5),GMax,GRMS,AMax,ARMS)
        If(Print) Write(IOut,1062) NPt,GMax,GRMS,AMax,ARMS
        Dot = SProd(NVar,SCr5,Scr5(1,4))
C       Attempt to use a linear step; Disabled
        If(Dot.lt.Zero.and.DoLin) then
C         get p"/S in Scr5
          Call ACasB(Nvar,Scr5(1,1),Scr5(1,4),Scr5(1,1),S2)
          Dot = One/Sqrt(SProd(NVar,Scr5,Scr5))
          Call AScale(Nvar,Dot,Scr5,Scr5)
C         get g at x"
          Call ASub(Nvar,Scr5,Scr5(1,4),Scr5(1,5))
          Call AScale(NVar,S2,Scr5(1,5),Scr5(1,5))
          Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,DQ(1,NPt),XX)
          Call AAdd(Nvar,Scr5(1,5),Q(1,NPt),Scr5(1,6))
C         x" is in Scr5(1,6)
          Call AMove(NVar,Q0(1,NPt),Scr5(1,5))
          Call AMove(NVTT,FC(1,NPt),Vec(1,NPStar))
C         Get Gradient at x" and put it in Scr5(1,3)
          Call GetGrd(IOut,IPrint,NVar,Scr5(1,5),Grad(1,NPt),
     $      Vec(1,NPSTar),Proj(1,NPt),Scr5(1,6),Scr5(1,3),Scr5(1,7),
     $      DotGP,NTBond,NTAng,NTDih)
          Call AScale(Nvar,DotGP,Scr5(1,3),Scr5(1,3))
          Call AMove(Nvar,Scr5(1,4),TOld)
          Call PTHARC(IOut,S,Scr,NVAR,Scr5(1,6),Scr5(1,7),TOld,
     $      Scr5,Scr5(1,2),Scr5(1,5),Scr5(1,4),DQ(1,NPt),Scr5(1,8),
     $      Scr5(1,3))
          Dot = One/S2
          Call AScale(Nvar,Dot,Scr5(1,4),Scr5(1,4))
          Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,DQ(1,NPt),XX)
          Call AAdd(NVar,Q(1,NPt),DQ(1,NPt),Q(1,NPt))
          Call FixDih(NTDih,Q(NTBond+NTAng+1,NPt))
          endIf
        Call AMove(Nvar,Scr5(1,4),TOld)
C       Do an IRC GS-Step
        Call AMove(NVTT,FC(1,NPt),Vec(1,NPStar))
C       Diagonalize Hessian ,if necessary, and transform into
C        eigenvector basis
        Call GetPEig(NVar,Scr,Vec(1,NPStar),Proj(1,NPt),Scr5,
     $    Scr5(1,4),Vec(1,Npt),Eigen(1,Npt),IDiag(Npt))
        Call AClear(NVar,Scr5)
        Call AClear(NVar,DQ(1,NPt))
        IConv = IConv+1
        If(Print) Write(IOut,3010) S2
        Call AScale(NVar,-One,Scr,Scr)
        Call AScale(NVar,S2,Scr5(1,4),Scr5(1,4))
        Call DxGS(IOut,NVar,Scr,Scr5(1,4),Vec(1,NPt),Eigen(1,NPt),S,
     $    DQ(1,NPt),FToler,USeSym,Print)
C       Project QST Tangent out of step
        Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,DQ(1,NPt),XX)
        Call GetFPA(NVar,DQ(1,NPt),DQ(1,NPStar),Scr5(1,5),XMax,XRMS,
     $    YMax,YRMS)
        If(Print) Write(IOut,1063) NPt,XMax,XRMS,YMax,YRMS
        DIDIRC=.true.
        If(Didirc) then
          IDiag(NPt) = 1
        else
          IDiag(Npt) = 0
          endIf
C       Scale Step
        Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,DQ(1,NPt),XX)
        Call MultLV(.True.,1,0,NVar,NVar,Proj(1,NPt),DQ(1,NPt),Scr)
        Call AMove(NVar,Scr,DQ(1,NPt))
        Dot = Sprod(NVar,DQ(1,NPt),DQ(1,NPt))
        Dot = Sqrt(Dot/GFloat(NVar))
        Scale = One
        DotMax = ArrMax(Dq(1,NPt),NVar,.True.)
        If(DotMax.gt.Weight*DQMax) Scale = Scale*Weight*DQMax/dotMax
        Call AScale(NVar,Scale,DQ(1,NPt),DQ(1,NPt))
        DotRMs = Sqrt(SProd(NVar,DQ(1,NPt),DQ(1,NPt))/GFloat(NVar))
        If(DotRMS.gt.(Weight*DQRMS)) then
          Scale = Weight*DQRMS/DotRMS
          Call AScale(NVar,Scale,DQ(1,NPt),DQ(1,NPt))
          endIf
        If((GMax.lt.FMax.and.GRMS.lt.FRMS).and.
     $     (DotRMS.lt.DQRMS.and.DotMax.lt.DQMax)) then
          IConv = IConv-1
          Call AClear(NVar,DQ(1,NPt))
          endIf
C       Rudimentary check for oscillation
        If(NPt.ne.(NPStar+N1)) then
          Call AClear(NVar,Scr5)
          Dot = SProd(NVar,DQ(1,NPt),DQ(1,NPt-N1))
          DotN = SProd(NVar,DQ(1,NPt),DQ(1,NPt))
          DotN1 = SProd(NVar,DQ(1,NPt-N1),DQ(1,NPt-N1))
          DN = DotN/GFloat(NVar)
          DotN = Sqrt(DotN)
          DN = Sqrt(DN)
          DN1 = DotN1/GFloat(NVar)
          DotN1 = Sqrt(DotN1)
          DN1 = Sqrt(DN1)
          If(DN.gt.Toler.and.DN1.gt.Toler) then
            Dot = Dot/(DotN*DotN1)
            If(Dot.lt.OSc) then
              Write(IOut,3040) NPt,(NPt-N1)
              Call AScale(NVar,Zero,DQ(1,NPt),DQ(1,NPt))
              endIf
            endIf
          endIf
C       Check overall displacement
        Call AClear(NVar,Scr5)
        Call AAdd(NVar,Q(1,NPt),DQ(1,NPt),Scr5)
        Call FixDih(NTDih,Scr5(NTBond+NTAng+1,1))
        Call ASub(NVar,Scr5,Q0(1,NPt),Scr)
        Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,Scr,XX)
        Call MultLV(.True.,1,0,NVar,NVar,Proj(1,NPt),Scr,Scr5)
        Dx = Sqrt(SProd(NVar,Scr5,Scr5))
        DMaxT = DXMaxT
        If(Dx.gt.DMaxT) then
          If(Print) Write(IOut,3050) NPt
          Call AClear(NVar,DQ(1,NPt))
          IConv = IConv-1
          endIf
   30   Continue
C
    6 Conv=.false.
      If(IConv.eq.0.or.NIter.eq.MaxIter) then
        Conv=.true.
        Write(IOut,3060) NIter
        If(Print.and.NVar.lt.20) Write(IOut,2001)
        endIf
      Call AClear(NVar,DQ(1,NPStar))
      Do 50 NPt = (NPStar+N1), NStop, N1
        Call AAdd(NVar,Q(1,NPt),DQ(1,NPt),Q(1,NPt))
        Call FixDih(NTDih,Q(NTBond+NTAng+1,NPt))
        Call ASub(NVar,Q(1,NPt),Q0(1,NPt),DQ(1,NPt))
        Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,DQ(1,NPt),XX)
        Dot = Sqrt(SProd(NVar,DQ(1,NPt),DQ(1,NPt)))
        Call MultLV(.True.,1,0,NVar,NVar,Proj(1,NPt),DQ(1,NPt),Scr)
        Call AMove(NVar,Scr,DQ(1,NPt))
        Call AAdd(NVar,Q0(1,NPt),DQ(1,NPt),Q(1,NPt))
        If(Conv) Call AMove(NVar,Q(1,NPt),Q0(1,NPt))
   50   Continue
      If(.not.Conv) goto 1
      Return
      End
*Deck DXQuad
      Subroutine DXQuad(IPrint,Sphere,AlUnit,FTemp,XNew,Vec,OK,NVar,
     $  NEig,FTempU,Eigen,Scr)
      Implicit Real*8(A-H,O-Z)
C
C     The current second derivative matrix and gradient are used to
C     compute the displacement from the current point needed to find
C     the minimum on a purely quadratic surface.  This is called
C     the quadratic search.  First an attempt is made to take a straight
C     Newton-Raphson step.  If the resulting step exceeds the trust
C     radius DXMaxT and Sphere is set, the minimum value of the quadratic
C     function on the surface of the sphere of radius DXMaxT is located.
C     Note that any eigenvalues of the Hessian which have incorrect sign
C     have their signs flipped before computing the step; this
C     corresponds to stepping away from the nearby stationary point in
C     (hopefully) the direction of one with the opposite curvature in
C     that mode.
C
      Real*8 MDCutO
      Logical OK, Prnt, Exit, Sphere, AlUnit, HavDep, Large, SmlEig
      Dimension FTemp(NVar), XNew(NVar), FTempU(NVar), Scr(NVar),
     $  Vec(NVar,NVar), Eigen(NVar)
      Common /OptGrd/ CONVF,FMAXT,DXMaxT,RMAX,RMIN,RLIM,EIGMAX,EIGMIN,
     $  DEPred,FSWTCH,FNCERR,GRDERR,FNCCNV,XXIRC,NLEVEL,NSTEP,ISTEP,
     $  NMAX,NP,NEG,PRNT,EXIT,IPSAV1,IPSAV2,HavDEP,NGDIIS,IRisen,MaxRis,
     $  Large,ILevel,ISStep,IOffs,IMDS,MaxSv,MaxSv2,NVrMMT
      COMMON /IO/     IN,IOUT,IPUNCH
      Save Zero
      Data Zero/0.0d0/
 1000 Format(' Quadratic step=',1PD9.3,' exceeds max=',1PD9.3,
     $       ' adjusted using Lamda=',1PD10.3,'.')
 1010 Format(' Angle between NR and scaled steps=',F7.2,' degrees.')
 1020 Format(' Angle between quadratic step and forces=',F7.2,
     $       ' degrees.')
C
      OK = .True.
C
C     Form displacements.  If the straight Newton-Raphson step exceeds
C     the trust radius, choose the best step using the quadratic
C     expansion on the surface of the sphere of radius DXMaxT.
C
      RLamda = Zero
      Small = MDCutO(0)
      Do 200 ITry = 1, 2
        DXMax = Zero
        SmlEig = .False.
        Do 120 I = 1, NVar
          XNew(I) = Zero
          If(AlUnit) then
            EI = Eigen(I)
            XNew(I) = XNew(I) + FTempU(I)/(EI-RLamda)
          else
            Do 110 J = 1, NEig
              EJ = Eigen(J)
              If(Abs(EJ).lt.Small) then
                EJ = Sign(Small,EJ)
                SmlEig = .True.
                endIf
  110         XNew(I) = XNew(I) + Vec(I,J)*FTempU(J)/(EJ-RLamda)
            endIf
          DXMax = Max(DXMax,Abs(XNew(I)))
  120     Continue
        DX = Sqrt(SProd(NVar,XNew,XNew))
        If(DX.le.DXMaxT.or..not.Sphere.or.SmlEig) goto 210
        If(ITry.eq.1) then
          Call AMove(NVar,XNew,Scr)
          Call RedStp(IOut,IPrint,NEig,DXMaxT,FTempU,Eigen,RLamda)
          Write(IOut,1000) DX, DXMaxT, RLamda
          endIf
  200   Continue
      Ang = AngDeg(0,NVar,XNew,Scr)
      Write(IOut,1010) Ang
  210 Ang = AngDeg(0,NVar,XNew,FTemp)
      Write(IOut,1020) Ang
      Return
      End
*Deck EigOvl
      Subroutine EigOvl(N,M,X,EigVal,EigVec,ISelec,SMax,Scr)
      Implicit Real*8(A-H,O-Z)
C
C     Put the eigenvector of EigVec which has the maximum overlap with X
C     first in V.  ISelec is set to the selected eigenvector and SMax to
C     its overlap.  An overlap must be greater than OvrMin times the next
C     largest overlap otherwise V is left unchanged.
C
      Dimension X(N), EigVal(M), EigVec(N,M), Scr(N)
      Save OvrMin, Zero, Two
      Data OvrMin/0.8D0/, Zero/0.D0/,Two/2.d0/
C
      ISelec = 1
      spmax = Zero
      spmax1 = Zero
      ispmax = 1
      Do 10 I = 1, M
        SPrd = Abs(SProd(N,X,EigVec(1,i)))
        if(sprd.gt.spmax) then
          spmax1 = spmax
          spmax = sprd
          ispmax = i
        else if(sprd.gt.spmax1) then
          spmax1 = sprd
          endif
   10   Continue
      If(spmax.gt.ovrmin.or.spmax.gt.(Two*spmax1)) then
        i = ispmax
        Call AMove(N,EigVec(1,i),Scr)
        Call AMove(N,EigVec(1,1),EigVec(1,i))
        Call AMove(N,Scr,EigVec(1,1))
        Tmp = EigVal(i)
        EigVal(i) = EigVal(1)
        EigVal(1) = Tmp
        SMax = spmax
        ISelec = i
        endIf
      Return
      End
*Deck FacMMM
      Function FacMMM(IOut,IPrint,DoSqrt,IScMMF,NVrMMT,NEig,EH)
      Implicit Real*8(A-H,O-Z)
C
C     Return a scale factor for MM forces for use in mixed optimizations.
C     EH are the QM Hessian eigenvalues.
C
C     IScMMF : -2 ... Zero MM forces and MM part of step, equivalent to
C                     not passing MM info.
C              -1 ... Zero MM forces but interpolate step in MM coords.
C               0 ... Default (-1).
C               1 ... Scale MM forces based on 1/RMS(QM-H-Eig).
C               N ... Scale MM forces based on 1/(RMS(QM-H-Eig)*N).
C
      Logical DoSqrt
      Dimension EH(*)
      Save Zero
      Data Zero/0.0d0/
 1000 Format(' DoSqrt=',L1,' IScMMF=',I6,' MM Scaling Factor = ',5F13.5)
C
      If(IScMMF.le.0.or.NVrMMT.eq.0) then
        FacMMM = Zero
      else
        If(DoSqrt) then
          FacMMM = Sqrt(GFloat(NEig)/SProd(NEig,EH,EH))
        else
          FacMMM = GFloat(NEig)/Sqrt(SProd(NEig,EH,EH))
          endIf
        FacMMM = FacMMM / GFloat(IScMMF)
        endIf
      If((IScMMF.gt.0.and.NVrMMT.gt.0.and.IPSave(0).eq.0)
     $  .or.IPrint.gt.0) Write(IOut,1000) DoSqrt, IScMMF, FacMMM
      Return
      End
*Deck FCClim
      Subroutine FCClim(N,W1,E,V1,V2,FC)
      Implicit Real*8(A-H,O-Z)
C
C     Adjust the force constant matrix for the climbing phase of LST/QST.
C
      Dimension V1(N), V2(N), FC(*)
C
      IJ = 1
      W1E = W1*(E+E)
      Do 10 I = 1, N
        Do 10 J = 1, I
        FC(IJ) = FC(IJ) + W1E*V2(I)*V2(J) - W1*(V1(I)*V2(J)+V2(I)*V1(J))
   10   IJ = IJ + 1
      Return
      End
*Deck FewEig
      Subroutine FewEig(N,M,NG,NE,E,NB,U,MU,UMU,UU,W,IOpt)
      Implicit Real*8(A-H,O-Z)
C
C     Find a few lower or upper eigenvalues of M using buffered
C     Lanczos method assuming NE << N.
C
C     N        : Dimensions of the problem
C     M        : Symmetric matrix, lower triangular form
C     NG       : Number of guess vectors privided in U
C     NE       : Number of eigenvectors to optimize
C     E        : Eigenvalues (vector of size NB)
C     NB       : Buffer size, number of vectors can be stored in
C                U and MU
C     U        : Eigenvectors NE is optimized on output, size is N*NB
C     MU       : Working space, size is N*NB
C     UMU      : Working space, lower triangular NB dimensional matrix
C     UU       : Working space, size is NB*NB
C     W        : Working space, size is 4*NB
C     IOpt     : Option flag:
C       0        Find the upper NE eigenvalue
C       1        Find the lower NE eigenvalue
C       2        Find the upper positive eigenvalues, use loose criteria
C                if the larger eigenvalue is negative.
C       3        Find the lower negative eigenvalues, use loose criteria
C                if the lowest eigenvalue is positive.
C       4        Find the all positive eigenvalues
C       5        Find the all negative eigenvalues
C       6        Find the all positive eigenvalues, use loose criteria
C                if the larger eigenvalue is negative.
C       7        Find the all negative eigenvalues, use loose criteria
C                if the lowest eigenvalue is positive.
C
C     This is an experimental version, tested for IOpt=7 only
C
C     O. Farkas, WSU, 1998.
C
      Real*8 M,MU,MDCutO
      Dimension M((N*(N+1))/2),E(NB),U(N,NB),MU(N,NB),
     $  UMU((NB*(NB+1))/2), UU(NB*NB), W(N)
      Logical DoLow, Auto, OK, New
      Save Zero, One, ELimT, ELimL, SLim, MaxStp
      Data Zero/0.D0/, One/1.D0/, ELimT/1.D-8/, ElimL/1.D-3/,
     $  SLim/0.71D0/, MaxStp/100/
C
      Eps = MDCutO(1)
      If(N.eq.1) then
        U(1,1) = One
        E(1) = M(1)
        Return
        endIf
      If(NB.lt.2) Call GauErr('No buffer in FewEig.')
C
C     Initialization
C
      DoLow = IAnd(1,IOpt).eq.1
      Auto  = IAnd(4,IOpt).eq.4
      NNB   = Min(NB,N)
      IET   = 3*NNB+1
      NNG   = Min(NG,NNB)
      NAB   = Min(NNB,NNG+1)
      SSign = One
      If (DoLow) SSign = -SSign
      Call AClear(NNB,E)
      Call AClear(NNB,W(IET))
      If(NNG.eq.0) then
        Call GuessE(N,DoLow,NNB,M,U,MU,MU(1,2))
        NAB = NNB
      else
C
C     Schmidt orthogonalization for the given guess set.
C
        RMin = MDCutO(0)
        RMax = One/MDCutO(0)
        Call Schmdt(U,MU,N,N,NNG,RMax,RMin,RMin)
        NAB  = NNG
        endIf
      Iter = 0
      NDid = 0
      NOK  = 0
      If(DoLow) then
        IDo = 1
      else
        IDo = Min(NAB,NNB-1)
        endIf
C
C     Start the main loop
C
  200 Iter = Iter + 1
C
C     Prepare and diagonalize Ut*M*U
C
      Call MultLV(.True.,NAB-NDid,0,N,N,M,U(1,NDid+1),MU(1,NDid+1))
      JI = 0
      Do 220 I = 1, NAB
        Do 220 J = 1, I
          JI = JI + 1
  220     UMU(JI) = SProd(N,U(1,I),MU(1,J))
      J0 = 0
      If(NNB.gt.NAB.and..not.DoLow) J0 = 1
      Call AMove(NAB,E,W(IET))
      Call DiagD(UMU,UU,E,NAB,W,W(NAB+1),NAB,.False.)
C
C     Update the buffered eigenvectors and product vectors
C
      NNAB  = NAB
      NNDid = NDid
      If(NAB.lt.NNB) then
        NAB = NAB+1
        If(.not.DoLow) IDo = IDo + 1
        endIf
      NDid = NAB-1
      NUpd = (NAB-NOK)/2
      If(NUpd.eq.0) Goto 999
      J0 = 1
      If(.not.DoLow) J0 = NNAB - NDid + 1
      If(.not.DoLow) J0 = NNAB - NNDid + 1
      Call ASub(NAB,E,W(IET),W(IET))
      JK0 = (J0-1)*NNAB
      Do 250 I = 1, N
        Call AClear(NDid,W)
        JK = JK0
        Do 240 J = 1, NDid
          Do 240 K = 1, NNAB
            JK = JK + 1
  240       W(J) = W(J) + UU(JK)*U(I,K)
        Do 250 J = 1, NDid
  250     U(I,J) = W(J)
      Do 260 I = 1, N
        Call AClear(NDid,W)
        JK = JK0
        Do 270 J = 1, NDid
          Do 270 K = 1, NNAB
            JK = JK + 1
  270       W(J) = W(J) + UU(JK)*MU(I,K)
        Do 280 J = 1, NDid
  280     MU(I,J) = W(J)
  260   Continue
C
C     Choose the actual eiginvectors to update and check the
C     convergence on the eigenvalue
C
  300 If(DoLow) then
        IDo = NOK + 1
      else
        IDo = NDid - NOK
        endIf
      OK = E(IDo)*SSign.ge.Zero
      Conv = ELimL
      If(OK) Conv = ELimT
      CConv = Max(Abs(E(IDo))*Conv,Conv*Conv)
      If(Abs(W(IET+IDo-1)).lt.CConv) then
        If((OK.and.Auto.and.NOK.lt.NNDid).or.
     $    (.not.OK.and.NOK.lt.Min(NE,NNDid))) then
          NOK = NOK + 1
          Goto 300
        else
          Goto 999
          endIf
        endIf
C
C     Lanczos update for the actual eigenvectors and check its
C     convergence
C
      If(DoLow) then
        I1 = NOK + 1
      else
        I1 = NNAB
        endIf
      NDidN = NDid
      II = NAB + 1
      I  = I1 - 1
      New = .True.
  500 I  = I + 1
      II = II - 1
      If (I.lt.II) then
        Call AMove(N,MU(1,I),MU(1,II))
        EE = SProd(N,U(1,I),MU(1,II))
        If(Abs(EE).lt.Eps) then
          If((OK.and.Auto.and.NOK.lt.NNDid).or.
     $      (.not.OK.and.NOK.lt.Min(NE,NNDid))) then
            II = II + 1
            Goto 500
          else
            Goto 999
            endIf
        else
          C12 = EE/Sqrt(SProd(N,MU(1,II),MU(1,II)))
          If(I.gt.I1.and.C12.gt.SLim.and.New) then
            Call MultLV(.True.,1,0,N,N,M,MU(1,NAB),MU(1,II))
            SS = Sqrt(SProd(N,MU(1,II),MU(1,II)))
            If(SS.Lt.Conv) then
              New = .False.
              Call MultLV(.True.,1,0,N,N,M,U(1,II),MU(1,II))
              II = II + 1
              Goto 500
              endIf
            Call AScale(N,One/SS,MU(1,II),MU(1,II))
            Call AMove(N,MU(1,II),MU(1,NAB))
            Do 520 J = 1, NAB
              EE = SProd(N,U(1,J),MU(1,II))
              If(Abs(EE).gt.Eps) then
                Call ACASB(N,MU(1,II),U(1,J),MU(1,II),-EE)
                SS = Sqrt(SProd(N,MU(1,II),MU(1,II)))
                If(SS.lt.Conv) then
                  New = .False.
                  Call MultLV(.True.,1,0,N,N,M,U(1,II),MU(1,II))
                  II = II + 1
                  Goto 500
                  endIf
                Call AScale(N,One/SS,MU(1,II),MU(1,II))
                endIf
  520         Continue
            Call AMove1(N*(II-I),N*I,N*(I-1),U)
            Call AMove1(N*(II-I),N*I,N*(I-1),MU)
            Call AMove(N,MU(1,II),U(1,II))
            Call MultLV(.True.,1,0,N,N,M,U(1,II-1),MU(1,II-1))
            I = I - 1
            Goto 500
            endIf
          endIf
        Call ACASB(N,MU(1,II),U(1,I),MU(1,II),-EE)
        SS = Sqrt(SProd(N,MU(1,II),MU(1,II)))
        If (SS.lt.Conv) then
          If((OK.and.Auto.and.NOK.lt.NNDid).or.
     $    (.not.OK.and.NOK.lt.Min(NE,NNDid))) then
            II = II + 1
            Goto 500
          else
            Goto 999
            endIf
        else
          Call AScale(N,One/SS,MU(1,II),MU(1,II))
          endIf
        Do 600 J = II + 1, NAB
          EE = SProd(N,U(1,J),MU(1,II))
          If(Abs(EE).lt.Eps) then
            If((OK.and.Auto.and.NOK.lt.NNDid).or.
     $      (.not.OK.and.NOK.lt.Min(NE,NNDid))) then
              II = II + 1
              Goto 500
            else
              Goto 999
              endIf
            endIf
          Call ACASB(N,MU(1,II),U(1,J),MU(1,II),-EE)
          SS = Sqrt(SProd(N,MU(1,II),MU(1,II)))
          If(SS.lt.Conv) then
            If((OK.and.Auto.and.NOK.lt.NNDid).or.
     $      (.not.OK.and.NOK.lt.Min(NE,NNDid))) then
              II = II + 1
              Goto 500
            else
              Goto 999
              endIf
          else
            Call AScale(N,One/SS,MU(1,II),MU(1,II))
            endIf
  600     Continue
          Call AMove(N,MU(1,II),U(1,II))
          Goto 500
        endIf
      NDid = Min(II-1,NDidN)
      If(NDid.eq.NNAB) Goto 999
      If((OK.and.Iter.lt.(10*MaxStp)).or.Iter.lt.MaxStp) Goto 200
C
C     This is the end
C
  999 NE = NOK
      NG = Min(NAB,NNB-1)
      Return
      End
*Deck FixHis
      Subroutine FixHis(IOut,IPrint,NVar,NP,XX,Scr,NTBond,NTAng,NTDih)
C
C     Make sure that point history has proper dihedral angles
C
      Implicit Real*8(A-H,O-Z)
      Dimension XX(NVar,*),Scr(NVar),XXA(1)
      Save XXA
      Data XXA/0.d0/
C
      Do 100 I=NP,2,-1
        Call ASub(NVar,XX(1,I),XX(1,1),Scr)
        Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,Scr,XXA)
        Call AAdd(Nvar,XX(1,1),Scr,XX(1,I))
 100    Continue
      Return
      End
*Deck FixSym
      Subroutine FixSym(RdCom,NAtoms,NOpAlP,SymAlP,IPrAlP)
      Implicit Real*8(A-H,O-Z)
C
C     For path relaxation.
C     If RdCom=.true. the sym. info is saved in NOpAlP,SymAlP and IPrAlP
C     If RdCom=.false. the sym. info is saved in RWF and Common
C
      Logical RdCom
#include "repall.inc"
      Dimension SymAlP(3,3,MaxOp),IPrAlP(NAtoms,MaxOp)
      Save Zero
      Data Zero/0.0d0/
C
      If(RdCom) then
        Call LodSym(.False.,.True.,0)
        NOpAlP = NOpAll
        Do 10 IAt = 1, NAtoms
          Do 10 IOp = 1, NOpAll
   10       IPrAlP(IAt,IOp) = NEqAll(IAt,IOp)
        Do 30 IOp = 1, NOpAll
          Do 30 J = 1, 3
            Do 30 I = 1, 3
   30         SymAlP(I,J,IOp) = RotAll(I,J,IOp)
      else
        NOpAll = NOpAlP
        NOpPtG = NOpAll
        Do 40 IOp = 1, NOpAlP
          Do 40 IAt = 1, NAtoms
   40       NEqAll(IAt,IOp) = IPrAlP(IAt,IOp)
        Do 60 IOp = 1, NOpAlP
          Do 60 I = 1, 3
            Do 50 J = 1, 3
   50         RotAll(I,J,IOp) = SymALP(I,J,IOp)
   60       RotAll(I,4,IOp) = Zero
        Call OutRep(NAtoms)
        endIf
      Return
      End
*Deck GDIRFO
      Subroutine GDIRFO(IOut,IPrint,AlUnit,NVar,Neg,NEig,Eigen,Vec,
     $  FTempU,XNew,RLam0,RLamda,UpdLam,OK)
      Implicit Real*8(A-H,O-Z)
C
C     Determine an RFO step.
C
      Logical AlUnit, UpdLam, OK
      Dimension Eigen(NEig), Vec(NVar,NEig), FTempU(NVar), XNew(NVar)
      Save Zero, Pt5, Four, Step, Big, MaxIt, Conv, SamEig
      Data Zero/0.0d0/, Pt5/0.5d0/, Four/4.0d0/, Step/0.05d0/,
     $  Big/1.d6/, MaxIt/999/, Conv/1.d-8/, SamEig/1.d-4/
 1000 Format(' RLam0=',1PD9.2,'.')
 1010 Format(' Iter=',I4,' RLamda=',1PD15.8,' Sum=',1PD15.8,'.')
 1020 Format(' RFO step:  Lambda=',1PD15.8,'.')
 1025 Format(' RFO step:  Lambda0=',1PD15.9,' Lambda=',1PD15.8,'.')
 1030 Format(' RFO step:  ModMin=',I2,' Lambda=',1PD15.8,' EMin=',
     $  1PD15.8)
 1040 Format(' RFO could not converge Lambda in',I5,' iterations.')
C
C     Setup, depending on whether we are seeking a minimum or saddle
C     point.  ModMin is the first mode to minimize rather than maximize.
C     Care is needed because some modes may be frozen.
C
      OK = .True.
      Call AClear(NVar,XNew)
      ModMax = 1
      If(Neg.eq.0) then
        ModMin = 1
      else
        ModMin = 2
        RLam0 = Pt5*Eigen(ModMax) +
     $          Pt5*Sqrt(Eigen(ModMax)**2+Four*FTempU(ModMax)**2)
        If(Abs(RLam0-Eigen(ModMax)).lt.Conv) RLam0=RLam0+Conv
        If(IPrint.gt.0) Write(IOut,1000) RLam0
        endIf
      If(.not.UpdLam) Goto 200
C
C     Solve for Lamda, unless there is only one mode.  Ensure that
C     RLamda is less than the smallest eigenvalue for mimization.
C
      If(NVar.gt.1) then
        RLamda = Min(Zero,Eigen(ModMin)-Step)
        RLamd1 = Eigen(ModMin)
        RLamd2 = -Big
        Do 100 Iter = 1, MaxIt
          Sum = Zero
          Do 10 I = ModMin, NEig
   10       Sum = Sum + FTempU(I)*FTempU(I)/(RLamda-Eigen(I))
          If(IPrint.gt.0) Write(IOut,1010) Iter, RLamda, Sum
          If(Abs(RLamda-Sum).le.Conv) goto 200
          If(Eigen(ModMin).gt.Zero) then
            RLamda = Sum
          else
            If(Sum.lt.RLamda) then
              RLamd1 = RLamda
            else
              RLamd2 = RLamda
              endIf
            If(RLamd2.eq.(-Big)) then
              RLamda = RLamda - Step
            else
              RLamda = Pt5*(RLamd1+RLamd2)
              If(Abs(RLamd1-RLamd2).le.Conv) goto 200
              endIf
            endIf
  100     Continue
        Write(IOut,1040) MaxIt
        OK = .False.
        Return
      else if(Neg.eq.0) then
        RLamda = Pt5*(Eigen(1)-Sqrt(Eigen(1)**2+Four*FTempU(1)**2))
      else
        RLamda = RLam0
        endIf
C
C     Now apply some final consistency checks and compute the step.
C
  200 If(UpdLam) then
        If(Neg.eq.0) then
          Write(IOut,1020) RLamda
        else
          Write(IOut,1025) RLam0, RLamda
          endIf
        endIf
      If(NVar.gt.1.and.(RLamda.gt.(Eigen(ModMin)-SamEig).or.
     $   (RLamda.gt.Zero.and.Eigen(ModMin).gt.Zero))) then
        Write(IOut,1030) ModMin, RLamda, Eigen(ModMin)
        OK = .False.
        Return
        endIf
      Do 210 I = 1, NEig
        If(I.lt.ModMin) then
          Temp = -FTempU(I) / (RLam0-Eigen(I))
        else
          Temp = -FTempU(I) / (RLamda-Eigen(I))
          endIf
        If(AlUnit) then
          XNew(I) = XNew(I) + Temp
        else
          Call ACasB(NVar,XNew,Vec(1,I),XNew,Temp)
          endIf
  210   Continue
      Return
      End
*Deck GeoDIS
      Subroutine GeoDIS(IOut,IPrint,IRwGDI,NVarT,NVar,NVrMMT,NMax,
     $  MaxGDI,MxRFOD,SqrtH,NGDIIS,NPDIIS,NNegE,NNeg,IScMMF,AMRFO,AMEn,
     $  AMSim,AMTrc,AMTmp,PotE,Coef1,Coef2,Scr,Map,NPoint,AQuad,ACC,AFF,
     $  ANCC,ANFF,PMask,PDIS,SimSw,FncErr,DidBck,FixedC,ConvP,SimOpt,
     $  CThrHd,CDAdj,Rises,EnDIS,DidDIS,NRisDI,DoIntr,V,MDV)
      Implicit Real*8 (A-H,O-Z)
C
C     Parameters:
C     IRwGDI : RWF that save DIIS data, allocated by RwFGED
C     NVarT  : Total number of coordinates (NVar+NVrMMT)
C     NVar   : Number of internal coordinates
C     NVrMMT : Number of MM coordinates (cartesian and rigid blocks).
C     NMax   : Max number of points that can be saved and used by l103 and l502
C     MaxGDI : Max number of points that can be used by DIIS thoughout optimization
C     MxRFOD : Max number of points to mix when mixing RFO steps
C     NGDIIS : Number of DIIS eligible geometries
C     NPDIIS : Number of DIIS eligible electron densities
C     NNegE  : Number of negative eigenvalues of Hessian
C     NNeg   : Index of critical point to find (NNeg=1 for TS search)
C     IScMMF : Handling of MM coordinates/forces; see FacMMM.
C     AMRFO,AMEn,AMSim,AMTrc,AMTmp:
C              Square matrices of the error vectors (MaxGDI+1)*(MaxGDI+1)
C     PotE   : Array of the total energies for the points used
C     Coef1  : DIIS coefficients
C     Coef2  : DIIS coefficients (output)
C     Scr    : Working space of length Max(6*NVarT,NMax)
C     Map    : Used to indicate which working points are being removed (output)
C     NPoint : Save the indices of DIIS points used in OptDIS
C     AQuad  : Array of quadratic steps of NGDIIS number of points
C     ACC    : Array of original coordinates of NGDIIS number of points
C     AFF    : Array of original forces of NGDIIS number of points
C     ANCC   : Array of normal coordinates of NGDIIS number of points
C     ANFF   : Array of forces in normal coordinate of NGDIIS number of points
C              Note: MM forces and coordinates are always in Cartesian
C     DoIntr : True if interpolating because energy rose.
C
      Parameter (MaxCLn=5,MaxRis=1,NInvIt=6)
      Dimension AMRFO(NGDIIS+1,NGDIIS+1), AMEn(NGDIIS+1,NGDIIS+1),
     $  AMSim(NGDIIS+1,NGDIIS+1), AMTrc(NGDIIS+1,NGDIIS+1),
     $  AMTmp(NGDIIS+1,NGDIIS+1), PotE(NGDIIS), Coef1(NGDIIS+1),
     $  Coef2(NGDIIS+1), Scr(*), Map(NGDIIS), NPoint(NGDIIS),
     $  AQuad(NVarT,NGDIIS), ACC(NVarT,NGDIIS), AFF(NVarT,NGDIIS),
     $  ANCC(NVarT,NGDIIS), ANFF(NVarT,NGDIIS), PMask(NMax),
     $  PDIS(NPDIIS,NPDIIS), V(*)
      Logical Prnt, Prnt1, Inv3, Rises, TstBck, DidSkp, DidBck, FixedC,
     $  ConvP, SimOpt, DidDIS, EnDIS, RFODIS, SimDIS, DoIntr, SqrtH,
     $  OKEnD, OK
      Character*4 NmDIIS
      Save Zero, One, ValMin, SmlDif, Half, TSCFTh, Ten
      Data Zero/0.0d0/,One/1.0d0/,ValMin/-1.0d0/,SmlDif/1.d-5/,
     $  Half/0.5d0/,TSCFTh/5.d-3/, Ten/10.0d0/
 1010 Format(' ',A,'DIIS uses ',I4,' points instead of ', I4)
 1020 Format(' Inversion required scale factor of ',1PD7.1)
 1030 Format(' DidBck=',L1,' Rises=',L1,1X,A,'DIIS coefs:',5F11.5)
 1035 Format('        ',1X,'       ',1X,1X,A,'DIIS coefs:',5F11.5)
 1040 Format(' ',A,'DIIS coefs:',5F11.5)
 1050 Format(' Energy Point ',I5,' ',D18.10)
 1060 Format(' RMS Error=',D18.10,' NUsed=',I2,' OKEnD=',L1,' EnDIS=',
     $  L1)
 1070 Format(' Guided ',A,'DIIS for TS search',5X,'SclDIS = ',F8.1)
 1080 Format(' ',A,'DIIS matrix:')
 1090 Format(' NNeg=',I2,' NGDIIS=',I2,' SimSw=',1PD10.2,' Rises=',L1,
     $  ' DC=',1PD10.2,' SmlDif=',1PD10.2,/,' NRisDI=',I2,' DoIntr=',L1,
     $  ' SqrtH=',L1,' CThrHd=',0PF5.1,' CDAdj=',F5.1)
 1110 Format(' SimOpt: fixing nuclei')
 1120 Format(' SimOpt: converging density')
 1130 Format(' SimOpt: GrScal=',F13.5,' PScal=',F13.5)
 1140 Format(' SimOpt: current RMS error in geom=',1PD9.2,' ConvP=',L1,
     $  ' ISCFTh=',I3,' SCF conv=',1PD9.2)
C
C     Start of the first order simultaneous optimization
C
      Prnt  = IPrint.ge.1
      Prnt1 = Prnt.or.IPSave(0).eq.0
      DidDIS= .False.
      NASq  = (NGDIIS+1)*(NGDIIS+1)
      ISkip = 0
      Call RwFGED(1,IRwGDI,NMax,MaxGDI,IEn,INNEn,ICoef,IQNDIS,IEnDIS,
     $  IMask,ICnvQM)
C
C     Set scaling factors
C
      If(SimOpt) Call GauErr('No SimOpt in GeoDIS.')
      GrScal = One
      PScal = Zero
C
C     Evaluate current point and set optimization options
C     Turn off EnDIS when the error is small enough
C
      If(IScMMF.eq.-2) then
        NVRMS = NVar
      else
        NVRMS = NVarT
        endIf
      GError = SProd(NVRMS,AQuad(1,1),AQuad(1,1))
      RMSErr = Sqrt(GError/GFloat(NVRMS))
C
C     EnDIS:  a(i,j) = g(i)R(i) + g(j)R(j) - g(i)R(j) - g(j)R(i)
C     RFODIS: a(i,j) = g(i)g(j)/H^2
C     SimDIS: a(i,j) = 0.5*(g(i)*g(i)/H + g(j)g(j)/H) +
C                      g(i)R(i) + g(j)R(j) - g(i)R(j) - g(j)R(i)
C
C     The following options are used for test purpose
      RFODIS = .True.
      EnDIS  = .False.
      SimDIS = .False.
      Rises  = .False.
      FncE1  = FncErr
      E1F = PotE(1) - FncE1
      Do 20 I = 2, NGDIIS
   20   Rises = Rises.or.E1F.gt.(PotE(I)+FncE1)
      If(Rises) then
        NRisDI = NRisDI + 1
      else
        NRisDI = 0
        endIf
      DoIntr = NRisDI.gt.MaxRis
      If(NNeg.gt.0) then
        OKEnD = NGDIIS.ge.4
        EnDIS = RMSErr.gt.SimSw
        If(RMSErr.le.SimSw) then
          SclDIS = One
        else
          SclDIS = GFloat(NVarT)
          endIf
      else
        OKEnD = NGDIIS.ge.2
        EnDIS = RMSErr.gt.SimSw.or.Rises
        endIf
      If(NGDIIS.gt.1) then
        DC = PotE(1) - PotE(2)
        OKEnD = OKEnD.and.(Abs(DC).ge.SmlDif.or.Rises)
      else
        DC = Zero
        endIf
      EnDIS = EnDIS.and.OKEnD
      If(Prnt1) Write(IOut,1090) NNeg, NGDIIS, SimSw, Rises, DC,
     $  SmlDif, NRisDI, DoIntr, SqrtH, CThrHd, CDAdj
      If(Prnt1) Write(IOut,1060) RMSErr, NGDIIS, OKEnD, EnDIS
      NmDIIS = 'Unk-'
C
C     If the error for the current point is small enough, take this step
C
      If(GError.le.ValMin) then
        DidDIS = .True.
        Call AClear(NGDIIS,Coef2)
        Coef2(1) = One
        Goto 999
        endIf
C
C     EnDIS:  a(i,j) = g(i)R(i) + g(j)R(j) - g(i)R(j) - g(j)R(i)
C
      NDim = NGDIIS + 1
      If(OKEnD.or.SimDIS) then
C
C       Prepare EnDIS elements e(i,j) = g(i)*R(j)
C       For TS serach,  weight error from first reaction coordinate by SclDIS
C
        Call AClear(NASq,AMTrc)
        Do 30 I = 1, NGDIIS
          If(Prnt) Write(IOut,1050) I,PotE(I)
          Do 30 J = 1, NGDIIS
            If(NNeg.gt.0) then
              AMTrc(I+1,J+1) = -SclDIS*ANFF(1,I)*ANCC(1,J)
     $          - SProd(NVarT-1,ANFF(2,I),ANCC(2,J))
            else
              AMTrc(I+1,J+1) = -SProd(NVarT,ANFF(1,I),ANCC(1,J))
              endIf
   30       Continue
      If(IPrint.ge.2) then
        Write(IOut,1080) 'Geometry trace En-'
        Call OutMat(IOut,1,AMTrc,NDim,NDim,NDim,NDim)
        endIf
C
C       Prepare EnDIS A matrix a(i,j)=e(i,i)+e(j,j)-e(i,j)-e(j,i)
C
        Call AClear(NASq,AMEn)
        Do 40 I = 0, NGDIIS
          Do 40 J = I, NGDIIS
            If(I.eq.0 .and. J.eq.0) then
              AMEn(1,1) = Zero
            else if(I.eq.0.or.J.eq.0) then
              AMEn(I+1,J+1) = -One
              AMEn(J+1,I+1) = AMEn(I+1,J+1)
            else
              AMEn(I+1,J+1) = AMTrc(I+1,I+1) + AMTrc(J+1,J+1)
     $          - AMTrc(I+1,J+1) - AMTrc(J+1,I+1)
              AMEn(J+1,I+1) = AMEn(I+1,J+1)
              endIf
   40       Continue
        NmDIIS = 'En-'
        If(IPrint.ge.1) then
          Write(IOut,1080) 'Total En-'
          Call OutMat(IOut,1,AMEn,NDim,NDim,NDim,NDim)
          endIf
C
C       For SimDIS, the energy term is approximated using second order infomation
C       so that each A matrix element has additional 0.5*(g(i)*g(i)/H+g(j)g(j)/H)
C       contribution
C
        If(SimDIS) then
          Call AMove(NASq,AMEn,AMSim)
          Do 50 I = 1, NGDIIS
            Do 50 J = I, NGDIIS
              If(NNeg.gt.0) then
                AMSim(I+1,J+1) = AMSim(I+1,J+1) +
     $            Half*(SclDIS*AQuad(1,I)*AQuad(1,I) +
     $            SclDIS*AQuad(1,J)*AQuad(1,J)+SProd(NVarT-1,AQuad(2,I),
     $            AQuad(2,I)) + SProd(NVarT-1,AQuad(2,J),AQuad(2,J)))
                AMSim(J+1,I+1) = AMSim(I+1,J+1)
              else
                AMSim(I+1,J+1) = AMSim(I+1,J+1) +
     $            Half*(SProd(NVarT,AQuad(1,I),ANFF(1,I)) +
     $            SProd(NVarT,AQuad(1,J),ANFF(1,J)))
                AMSim(J+1,I+1) = AMSim(I+1,J+1)
                endIf
   50         Continue
          NmDIIS = 'Sim-'
          If(IPrint.ge.1) then
            Write(IOut,1080) 'Total Sim-'
            Call OutMat(IOut,1,AMSim,NDim,NDim,NDim,NDim)
            endIf
          endIf
        endIf
C
C     RFO-DIIS
C
      If(RFODIS) then
        Call AClear(NASq,AMRFO)
        Do 70 I = 0, NGDIIS
          Do 70 J = I, NGDIIS
            If(I.eq.0.and.J.eq.0) then
              AMRFO(1,1) = Zero
            else if(I.eq.0.or.J.eq.0) then
              AMRFO(I+1,J+1) = -One
              AMRFO(J+1,I+1) = AMRFO(I+1,J+1)
            else
              If(NNeg.gt.0) then
                AMRFO(I+1,J+1) = SclDIS*AQuad(1,I)*AQuad(1,J) +
     $            SProd(NVarT-1,AQuad(2,I),AQuad(2,J))
                AMRFO(J+1,I+1) = AMRFO(I+1,J+1)
              else
                AMRFO(I+1,J+1) = SProd(NVarT,AQuad(1,I),AQuad(1,J))
                AMRFO(J+1,I+1) = AMRFO(I+1,J+1)
                endIf
              endIf
   70       Continue
        If(IPrint.ge.2) then
          Write(IOut,1080) 'Geometry RFO-'
          Call OutMat(IOut,1,AMRFO,NDim,NDim,NDim,NDim)
          endIf
        NmDIIS = 'RFO-'
        If(IPrint.ge.1) then
          Write(IOut,1080) 'Total RFO-'
          NDim = NGDIIS + 1
          Call OutMat(IOut,1,AMRFO,NDim,NDim,NDim,NDim)
          endIf
        endIf
C
C     Compute the coefficients by inverting the A matrix
C
      If(RFODIS.or.SimDIS) then
        NUsed = NGDIIS
        Call ISeq(NUsed,0,1,Map)
        If(NUsed.gt.MxRFOD.and.MxRFOD.gt.0) then
          NUsed = MxRFOD
          Call Pcck(0,Junk,AMRFO,NGDIIS+1,NGDIIS+1,NUsed+1,NUsed+1)
          Call Pcck(0,Junk,AMSim,NGDIIS+1,NGDIIS+1,NUsed+1,NUsed+1)
          Call IClear(NGDIIS-NUsed,Map(NUsed+1))
          endIf
   80   NDim = NUsed + 1
        FactA = One
        Call AClear(NDim,Coef1)
C       Account for some numerical silliness in LAPack by
C       scaling the matrix if necessary.
        Do 85 InvIt = 0, NInvIt
          If(RFODIS) Call AScale(NDim**2,FactA,AMRFO,AMTmp)
          If(SimDIS) Call AScale(NDim**2,FactA,AMSim,AMTmp)
          IPrtM1 = IPrint - 1
          DidDIS = Inv3(IOut,IPrtM1,130,NDim,AMTmp,V,MDV)
          If(DidDIS) then
            If(Prnt1.and.InvIt.gt.0) Write(IOut,1020) FactA
            Call AScale(NUsed,-FactA,AMTmp(2,1),Coef1)
            Goto 90
            endIf
   85     FactA = FactA*Ten
C
C       DIIS coefficients check, if one of the older points has a very
C       large coefficient, delete related rows and columns and solve
C       the DIIS equation again.
C
   90   If(RFODIS) Call RemPnt(IOut,Prnt1,DidDIS,AMRFO,Coef1,Coef2,Map,
     $    NUsed,NGDIIS,CThrhd,DidSkp)
        If(SimDIS) Call RemPnt(IOut,Prnt1,DidDis,AMSim,Coef1,Coef2,Map,
     $    NUsed,NGDIIS,CThrhd,DidSkp)
        If(DidSkp.and.NUsed.gt.1) Goto 80
        DidDIS = NUsed.gt.1
        CLower = Min(One-CDAdj,-One)
        If(DidDIS) Call ChkDIS(IOut,Prnt,Prnt1,NGDIIS,Coef2,CDAdj,
     $    CLower,CThrhd,DidDIS)
        If(NUsed.lt.NGDIIS.and.DidDIS)
     $    Write(IOut,1010) NmDIIS(1:LinEnd(NmDIIS)), NUsed, NGDIIS
        endIf
C
C     Use EnCoef solver with constraint c(i)>0
C
      EnDIS = EnDIS.or.(OKEnD.and..not.DidDIS)
      If(EnDIS) then
        Call ISeq(NGDIIS,NGDIIS+1,-1,Map)
        NUsed = NGDIIS
        Call AClear(NUsed+1,AMTrc(1,1))
        AMTrc(NUsed+1,1) = One
        If(NUsed.gt.MxRFOD.and.MxRFOD.gt.0.and..not.DoIntr) then
          NUsed = MxRFOD
          Call Pcck(0,Junk,AMEn,NGDIIS+1,NGDIIS+1,NUsed+1,NUsed+1)
          Call Pcck(0,Junk,AMTrc,NGDIIS+1,NGDIIS+1,NUsed+1,NUsed+1)
          Call IClear(NGDIIS-NUsed,Map)
          endIf
        Call RevEnM(NUsed,PotE,AMEn,AMTrc,AMTmp)
  100   Call EnCoef(IOut,IPrint-1,NUsed,NUsed,PotE,AMEn,AMTrc,Scr,
     $    ISkip,IWarn)
        If(NUsed.gt.1) then
          Call RmEnPt(IOut,Prnt1,NUsed,NGDIIS,IWarn,AMEn,AMTrc,PotE,
     $      Map,Coef1,DidSkp)
          If(DidSkp.and.NUsed.gt.1) goto 100
          endIf
        If(NUsed.gt.1.and.(Coef1(NGDIIS).lt.(One-SmlDif)
     $    .or..not.DoIntr)) then
          Do 110 I = 1, NGDIIS
  110       Coef2(I) = Coef1(NGDIIS-I+1)
          DidDIS = .True.
          NmDIIS = 'En-'
        else
          EnDIS = .False.
          endIf
        If(IPrint.ge.1.and.DidDIS) Write(IOut,1040)
     $    NmDIIS(1:LinEnd(NmDIIS)), (Coef2(I),I=1,NGDIIS)
        endIf
C
C     Save coefficients, deleted points have been masked by
C     zero coefficients
C
      Call FileIO(1,-IRwGDI,NGDIIS,Coef2,ICoef)
      If(DidDIS) then
        DidBck = TstBck(NGDIIS,Coef2).or.(EnDIS.and.DoIntr)
        LNm = LinEnd(NmDIIS)
        If(NNeg.gt.0) then
          Write(IOut,1070) NmDIIS(1:LNm), SclDIS
        else
          Write(IOut,1030) DidBck, Rises, NmDIIS(1:LNm),
     $      (Coef2(I),I=1,Min(NGDIIS,MaxCLn))
          Do 200 ISt = (MaxCLn+1), NGDIIS, MaxCLn
            LSt = Min(ISt+MaxCLn-1,NGDIIS)
            Write(IOut,1035) NmDIIS(1:LNm), (Coef2(I),I=ISt,LSt)
  200       Continue
          endIf
      else
        DidBck = .False.
        endIf
  999 Return
      End
*Deck GetFPA
      Subroutine GetFPA(NVar,F,T,Scr,GMax,GRMS,AMax,ARMS)
      Implicit Real*8(A-H,O-Z)
C
C Get the  Max and RMS values for the projection of F on T.
C G : perp. ;A : Along
C
      Dimension  F(Nvar),T(NVar),Scr(Nvar)
C
      DotT=SProd(Nvar,T,T)
      DotTG=SProd(Nvar,T,F)
      Dot=DotTG/DotT
      Call ACasB(Nvar,F,T,Scr,-Dot)
      Call RMSVec(Nvar,Scr,GRMS,GMax)
      Call ASub(Nvar,Scr,F,Scr)
      Call RMSVec(Nvar,Scr,ARMS,AMax)
      Return
      End
*Deck GetGrd
      Subroutine GetGrd(IOut,IPrint,NVar,Q0,F,FC,Proj,Q1,FNew,Scr,FNorm,
     $  NTBond,NTAng,NTDih)
      Implicit Real*8(a-h,o-z)
C
C     Estimate new gradient given DQ=Q1-Q0, F and FC.
C
      Dimension Q0(NVar),Q1(NVar),F(NVar),FNew(Nvar)
      Dimension Proj(*),FC(*),Scr(NVar,2),XX(1)
      Save One,XX
      Data One/1.0d0/,XX/0.d0/
C
C     Q0 is the point at which energy derivatives have been
C     computed
      Call ASub(NVar,Q1,Q0,Scr)
      Call FixDq(IOut,IPrint,1,NVar,NTbond,NTAng,NTDih,Scr,XX)
C     Projects out redundancy out of displacement
      Call MultLV(.True.,1,0,Nvar,NVar,Proj,Scr,Scr(1,2))
      Call AMove(NVar,Scr(1,2),Scr)
C     Get H.Dq
      Call MultLV(.True.,1,0,Nvar,NVar,FC,Scr,Scr(1,2))
      Call AMove(NVar,Scr(1,2),Scr)
C     Get -dE/dq-d2E/dq2.dq
      Call ASub(NVar,F,Scr,Scr)
C     Project out the new gradient
      Call MultLV(.True.,1,0,Nvar,NVar,Proj,Scr,FNew)
C     Normalize
      FNorm = Sqrt(SProd(Nvar,FNew,FNew))
      Dot = One/FNorm
      Call AScale(NVar,Dot,FNew,FNew)
      Return
      End
*Deck GetPEig
      Subroutine GetPEig(NVar,F,FC,Proj,Scr,T,P,Eig,IDiag)
      Implicit Real*8(a-h,o-z)
C
C     For path relaxation.
C     P : eigenvectors, Eig : Eigenvalues.
C     Transform T and F in eigenvector basis
C     Diagonalize FC if IDiag=0
C
      Dimension F(NVar),Proj(*),FC(*),Scr(NVar,3),T(NVar),P(NVar*NVAr),
     $  Eig(NVar)
      Save BigEig,Thres
      Data BigEig/1.0d3/, Thres/1.0d-10/
C
      If(IDiag.eq.0) then
        Call PFP(NVar,Proj,FC,P)
        Call AClear(NVar,Scr)
        Call Add1MP(NVar,BigEig,Proj,FC)
        Call DiagD(FC,P,Eig,NVar,Scr,Scr(1,2),NVar,.False.)
        endIf
C     Form T,Forces in Eigenvector Basis
      Call AClear(Nvar,SCr)
      Call MatMP1(1,0,1,1,1,NVar,NVar,T,P,Scr)
      Call RndVec(NVar,Scr,Thres)
      Call AMove(Nvar,Scr,T)
      Call AClear(Nvar,SCr)
      Call MatMP1(1,0,1,1,1,NVar,NVar,F,P,Scr)
      Call RndVec(NVar,Scr,Thres)
      Call AMove(Nvar,Scr,F)
      Return
      End
*Deck GetPthl
      Subroutine GetPthL(IOut,IPrint,Dist,Ang,ArcL)
      Implicit real*8 (A-H,O-Z)
C
C     For path relaxation.  Get the arc of circle (or parabola) length
C
      Real*8 MDcutO
      Save One,Two,four5
      Data One/1.0d0/,Two/2.0d0/,four5/45.0d0/
 1000 Format(1x,' Angle (RX,Tr) = ',F10.5,' Degrees ')
 1010 Format(1x,' Arc length R->X =',F10.5 )
C
      Pi4 = ATan(GFloat(1))
      Small = MDCutO(0)
      If(Ang.le.Pi4) then
C       Arc of circle
        ArcL = Dist
        If(Abs(Ang).gt.Small) ArcL = Ang*Dist/Sin(Ang)
      else
C       Parabola
        Theta = Ang-Pi4
        u21 = Sqrt(Tan(theta)**2+One)
        ArcL = Tan(Theta)*u21+GLog(Abs(GTan(theta)+u21))
        ArcL = (ArcL/(One+GTan(theta)))*(Dist*Cos(theta))
        ArcL = ArcL+Dist*Pi4*Sqrt(two)
        endIf
      Theta = Ang*four5/Pi4
      If(IPrint.ge.1) then
        Write(IOut,1000) Theta
        Write(IOut,1010) ArcL
        endIf
      Return
      End
*Deck GetTan
      Subroutine GetTan(IOut,IPrint,NVar,DoAvg,RX,TR,V,TX,S,Proj,Scr,
     $  ArcL)
      Implicit Real*8 (A-H,O-Z)
C
C     For path relaxation
C     Construct tangent for Path Relaxation
C
C     RX : R-X
C     TR : Tangent @R
C     TX : Tangent @X
C     V  : normalized force
C     S : Step size for constrained optimization
C     Proj: projector
C     DoAvg: whether to do average of TX and V
C
      Logical DoAvg,DoPar
      Real*8 MDCutO
      Dimension RX(NVar),TR(NVar),TX(NVar),Scr(NVar),Proj(*),V(*)
      Save Zero,One,Two,Four
      Data Zero/0.0d0/,One/1.0d0/,Two/2.0d0/,Four/4.0d0/
C
      Dot1 = SProd(NVar,RX,RX)
      RRX = Sqrt(Dot1)
      Call MultLV(.True.,1,0,NVar,NVar,Proj,TR,Scr)
      Call AMove(Nvar,Scr,TR)
      Dot = One/Sqrt(SProd(Nvar,TR,TR))
      Call AScale(Nvar,dot,TR,TR)
      Dot3 = SProd(NVar,RX,TR)
C     Get arc length
      Test = -Dot3/RRX
      If(Dot3/RRX.gt.One) Test = One
      If(Dot3/RRX.lt.-One) Test = -One
      Ang = ACos(Test)
      XPi4 = ATan(One)
      IPrinL = IPrint
      If((Ang-XPi4).gt.Zero) IPrinL = 1
      Call GetPthL(IOut,IPrinL,RRX,Ang,ArcL)
      Test = -Dot3/Sqrt(Dot1)
      alpha = -Abs(Dot1/(Two*Dot3) )
      Dot = -One
      Call Ascale(Nvar,Dot,RX,RX)
      Call AClear(NVAr,TX)
      Call ACasB(Nvar,RX,TR,TX,Alpha)
      Dot = One/Sqrt(SProd(Nvar,TX,TX))
      Call Ascale(Nvar,Dot,TX,TX)
      DoPar = .False.
      If(test.lt.One/Sqrt(two)) then
C       Do Parabola
        DoPar = .true.
        Theta = ACos(test)
        Theta = Theta-XPi4
        Dot = -SProd(Nvar,RX,TR)
        Call ACasB(Nvar,RX,TR,Scr,Dot)
        Dot = One/Sqrt(SProd(Nvar,Scr,Scr))
        Call AScale(Nvar,Dot,Scr,Scr)
        Call ASub(Nvar,TR,Scr,TX)
        Call AAdd(Nvar,TR,Scr,Scr)
        Dot = One/Sqrt(SProd(Nvar,Scr,Scr))
        Call AScale(Nvar,Dot,Scr,Scr)
        Dot = One/Sqrt(SProd(Nvar,TX,TX))
        Call AScale(Nvar,Dot,TX,TX)
        Dot = -(one+two*GTAn(Theta))
        Call ACasb(Nvar,Scr,TX,TX,Dot)
        Dot = One/Sqrt(SProd(Nvar,TX,TX))
        Call AScale(Nvar,Dot,TX,TX)
        endIf
      If(DoAvg.and.DoPar) Call AAdd(NVar,TX,V,TX)
      Call MultLV(.True.,1,0,NVar,NVar,Proj,TX,Scr)
      Call AMove(Nvar,Scr,TX)
      Dot = One/Sqrt(Sprod(Nvar,TX,TX))
      Call Ascale(Nvar,Dot,TX,TX)
C     Get Step Size
      Dot2 = SProd(Nvar,RX,TX)
      Dot4 = SProd(NVar,TX,TR)
      Small = MDCutO(0)
      AA = (One-Dot4)*Two
      BB = (-Dot3+Dot2)*two
      CC = -Dot1
      If(Abs(AA).gt.Small) then
        DD = Sqrt(BB*BB-Four*AA*CC)
        SS = (-BB+DD)/(Two*AA)
        SSS = (-BB-DD)/(Two*AA)
      else
        SS = -CC/BB
        SSS = -CC/BB
        endIf
      S = SS
      If(SSS.gt.Zero) S = SSS
      Return
      End
*Deck GrdOpt
      Subroutine GrdOpt(IPrint,UpTrus,Sphere,DoLinS,MethQ,FCart,TstEig,
     $  SetStp,UpdRWF,ToAng,IDimSc,SteepM,F,Frcnst,X,XName,FC,FC2,ESave,
     $  FS,XX,FF,IC,ITU,XQuad,XNew,FTemp,A,VecInt,EigSav,Scr,FCTemp,
     $  APrint,IFlags,NAtoms,NAtNew,IAn,AtmChg,C,DoClnS,TRot,MthUpd,
     $  Wght,IQST,DoRed,NTrRot,NTBond,NTAng,NTDih,IZRed,Value,Cnstr1,
     $  Cnstr2,CNew,IBMat,BMat,ScrGDI,Failed,IOSLEq,ConvEq,LenSB,LenCB,
     $  IOldRC,IOldRB,AlUnit,FRMS,NPDir,OptMic,IRwF1,MicOpt,CMM,TrVec,
     $  PDis,AtChMM,IAtTyp,NumTpS,IAtTpS,MaxLEF,EField,IRwMM,NVarMM,
     $  IDBMat,DBMat,FTmp,ISCRF,IDir,IGuess,IRFO,EQMS,XXMM,FFMM,XMMT,
     $  HMod,HModO,CNewQ,IMMRFO,IMMMOD,IMMCRS,FTp,DT,ICnt,COrig,T,COM,
     $  NRF,RotQM,JScan,IZRedO,IHVAcc,DoMicO,SkpMic,IGetTS,IAprBG,
     $  IGISav,IGIOld,IRwSM,IPFlag,AllowP,FMM,FMFlag,FMFlg1,NFxFlg,
     $  IHMeth,NVarQ,FFQ,XXQ,FTmpQ,XQT,IntHes,IRwF2,DoEle,FFMMQ,FMMQ,
     $  HReal,NTRedS,IDir0,ONLnOK,IScMMF,DidMix,IGDRFO,FixedC,ConvP,
     $  NGoDwn,IRwCrd,SavLMO,IEStpM,DXLimC,NFroz,Frozen,IGEVer,SwitMx,
     $  SimSw,RawDIS,MxRFOD,CThrHd,CDAdj,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     This is the main routine for the Berny optimization procedure.
C     It makes decisions concerning what types of search to carry out.
C     It is also responsible for most of the output and for determining
C     when the stationary point has been found.  Note that FS(1) is
C     assumed equivalenced to ESave(2)
C
C     Note that CMM is the regular /Mol/ C here, while C is the temporary
C     copy in Optmz; this should be fixed later.
C
#include "commonz.inc"
      Parameter (MxITry=10)
      Integer GetChr, FMFlag, FMFlg1
      Logical OK, Linr, Steped, Prnt, Exit, Rises, DidLin, Quad, FixEig,
     $  LinMin, LinTS, Despar, LJunk, UpTrus, Sphere, TstEig, Failed,
     $  DoLin, DidSph, SetStp, Finish, SclEig, FCart, HavDep, DoClnS,
     $  DoRed, LSTQST, Swap, DidGI, RegFin, SmlFin, SimMix, Large,
     $  InitNR, DoLinS, DoDIIS, DoMatC, BigLin, AlUnit, Just1, DidSwap,
     $  OptMic, SimOpt, FixedC, OK1, RotQM, LogTmp, DoMicO, SmlStp,
     $  SkpMic, DdTSLn, ConvP, FMM, AllowP(*), Mixed, QCoupl, DidDgF,
     $  IntHes, DoEle, LnOnly, DoSD, Compct, DoDipD, ONLnOK, DidMix,
     $  UpdRWF, DidBck, DidEig, SimPrj, OKHess, CkHesU, SavLMO, SavXLn,
     $  SimMMM, Frozen(*), DidDIS, OldGDI, SqrtH, RawDIS
      Dimension F(*), Frcnst(*), X(NVar), XName(*), FC(*), FS(*),
     $  XX(NVar,*), FF(*), IC(*), XQuad(*), XNew(*), TRot(3,3), IAn(*),
     $  FTemp(*), A(*), ESave(*), Scr(IDimSc,8), FCTemp(*), ITU(*),
     $  APrint(6,*), FC2(*), IFlags(*), C(3,*), AtmChg(*), VecInt(*),
     $  Cnstr1(NVar), Cnstr2(NVar), TrRMat(6), IZRed(*), Value(*),
     $  CNew(3,*), IBMat(4,*), BMat(3,4,*), ScrGDI(*), IIniNR(4),
     $  CMM(3,*), MicOpt(*), TrVec(3), V(MDV), XXA(1), FTmp(*),
     $  PDis(*), AtChMM(*), IAtTyp(*), EField(0:3), IDBMat(*), DBMat(*),
     $  CnvQM(4), FFMM(*), XXMM(*), XMMT(3,*), HMod(*), CNewQ(3,*),
     $  FTp(3,*), DT(*), ICnt(*), COrig(3,*), T(3,3,10,*), COM(3,*),
     $  HModO(*), JScan(4), JJ(1), IZRedO(4,*), What(1), QMVar(30),
     $  IAtTpS(*), FFQ(*) ,XXQ(*) ,FTmpQ(*) ,XQT(*), EigSav(*), EQMS(*),
     $  FFMMQ(*), FMMQ(*), HReal(*), NTRedS(*), RWhat(1)
      Equivalence (InitNR,IIniNR(1))
      Common /IO/     In, IOut, IPunch
      Common /Grdnt/  Energy,NVar,IGetFC
      Common /OptGrd/ ConvF,FMaxT,DXMaxT,RMax,RMin,RLim,EigMax,Eigmin,
     $  DEPred,FSwtch,FncErr,GrdErr,FncCnv,RLam,NRisDI,NStep,IStep,NMax,
     $  NP,Neg,Prnt,Exit,IPSav1,IPSav2,HavDEP,NGDIIS,IRisen,MaxRis,
     $  Large,ILevel,ISStep,IOffs,IMDS,MaxSv,MaxSv2,NVrMMT
      Integer IReslt(3),Tmp(40),IStr(20)
      Save Zero,Two,CnvFMx,CNVFX,NVarPt,F100,SmStMx,ScMin,XXA,JJ,Pt5,
     $  One,RMinLn,Big,SclMMC
      Data Zero,Two,CnvFMx,CNVFX/0.D0,2.D0,1.5D0,4.D0/,XXA/0.0d0/
      Data NVarPt/20/,F100/100.0d0/,SmStMx/0.02d0/,ScMin/0.01d0/,JJ/0/,
     $  Pt5/0.5d0/,One/1.0d0/,RMinLn/1.d-3/,Big/1.d10/,SclMMC/0.1d0/
 1000 Format('         Item               Value     Threshold  ',
     $       'Converged?')
 1010 Format(' Maximum Force            ',F8.6,5X,F8.6,5X,3A1)
 1020 Format(' RMS     Force            ',F8.6,5X,F8.6,5X,3A1)
 1030 Format(' Maximum Displacement     ',F8.6,5X,F8.6,5X,3A1)
 1040 Format(' RMS     Displacement     ',F8.6,5X,F8.6,5X,3A1)
 1050 Format(' Maximum MM Force         ',F8.6,5X,F8.6,5X,3A1)
 1060 Format(' RMS     MM Force         ',F8.6,5X,F8.6,5X,3A1)
 1070 Format(' The second derivative matrix:')
 1080 Format(' Predicted change in Energy=',1PD13.6)
 1090 Format(' Search for a local minimum.')
 1100 Format(' Search for a saddle point.')
 1110 Format(' Search for a stationary point of order',I3)
 1120 Format(' Step number',I4,' out of a maximum of',I5)
 1130 Format(' Step number',I4,' out of a maximum of',I5,
     $  ' on scan point',I6,' out of',I6)
 1140 Format(' All quantities printed in internal units',
     $' (Hartrees-Bohrs-Radians)')
 1150 Format(' Variable       Old X    -DE/DX',3(3X,7HDelta X),
     $ '     New X',/,33X,'(Linear)    (Quad)   (Total)')
 1170 Format(1X,8A1,2X,6F10.5)
 1180 Format(' Optimization completed.',/,
     $       '    -- Stationary point found.')
 1190 Format(' Optimization completed on the basis of negligible',
     $  ' forces.',/,'    -- Stationary point found.')
 1200 Format(' Optimization stopped.',/,
     $       '    -- Number of steps exceeded,  NStep=',I4,/,
     $       '    -- Flag reset to prevent archiving.')
 1210 Format(' Optimization stopped.',/,
     $       '    -- Wrong number of Negative eigenvalues: Desired=',
     $       I3,' Actual=',I3,/,
     $       '    -- Flag reset to prevent archiving.')
 1220 Format(' Linear search step of',F8.3,' exceeds DXMaxT=',F8.3,
     $  ' scaled by ',F7.3)
 1230 Format(' Second derivative matrix not updated -- first step.')
 1250 Format(' Second derivative matrix not updated -- turned off.')
 1260 Format(' Second derivative matrix not updated --',
     $       ' analytic derivatives used.')
 1270 Format(' Linear search not attempted -- option 19 set.')
 1280 Format(' Linear search not attempted -- no second derivatives',
     $       ' available.')
 1290 Format(' Linear search not attempted -- first point.')
 1300 Format(' Linear search not attempted -- second point.')
 1310 Format(' Linear search not attempted -- rms force is less than',
     $       ' FSwtch (',F8.5,').')
 1320 Format(' Linear search not attempted -- Energy or forces',
     $       ' went down.')
 1330 Format(' Energy rises -- skip ',A,' search.')
 1340 Format(' Energy rises -- skip steepest descent search.')
 1350 Format(' Gradient rises -- skip Quadratic search.')
 1360 Format(' Maximum step size (',f8.3,
     $       ') exceeded in Quadratic search.')
 1363 Format('    -- Step size not scaled.')
 1365 Format('    -- Step size scaled by ',F7.3)
 1370 Format(' Optimization aborted.',/,
     $       '    -- Gradient out of range.',/,
     $       '    -- Maximum allowed force = ',G12.3)
 1380 Format(' Optimization aborted.',/,'   -- No acceptable step.')
 1390 Format(' Quartic linear search produced a step of ',F8.5,'.')
 1400 Format(' Quintic linear search produced a step of ',F8.5,'.')
 1410 Format(' Grad norm linear search produced a step of ',F8.5,'.')
 1415 Format(' QM/MM linear search produced a step of ',F8.5,'.')
 1420 Format(' Skip linear search -- no minimum in search direction.')
 1430 Format(' Steepest descent instead of Quadratic search.')
 1440 Format(' Steepest descent step scaled to max of',F8.5,'.')
 1450 Format(' Steepest descent step scaled up to rms of',F8.5,'.')
 1460 Format(' Unrecognized ICase=',I10,'.')
 1470 Format(' MatCor failed.')
 1480 Format(' Swapping is turned off.')
 1490 Format(' Use linear search instead of GDIIS.')
 1500 Format(' Second derivative matrix not updated --',
     $       ' unit matrix used.')
 1510 Format(' Take small step and do not check convergence.')
 1520 Format(' Linear search for TS vector ',F8.5,
     $       ' (skip quadratic step).')
 1530 Format(' Linear search for TS failed, do quadratic step.')
 1540 Format(' Linear search step of',F8.3,' exceeds DXMaxT=',F8.3,
     $  ' but not scaled.')
 1550 Format(' ITry=',I2,' IFail=',I1,' DXMaxC=',1PD9.2,' DCOld=',
     $  1PD9.2,' DXMaxT=',1PD9.2,' DXLimC=',1PD9.2,' Rises=',L1)
C
      XLin = Zero
      SqrtH = IGetFC.ne.4
      DidDIS = .False.
      DidBck = .False.
      DidEig = .False.
      DidDgF = .False.
      SavXLn = .False.
      DoMatC = (IOldRC.ge.2.or.FCart).and.NTrRot.eq.-1.and.NPDir.eq.0
      Call AMove(3*NAtoms,C,CNew)
      If(Large) then
        IOSLE1 = IOSLEq
        IWght  = 0
      else
        IOSLE1 = 0
        IWght  = 0
        endIf
      NAt3 = 3*NAtNew
      NAt3TT = (NAt3*(NAt3+1))/2
      NNT = (NVar*(NVar+1))/2
      LSTQST = Mod(IQST,10).ne.0
      IStep = IStep + 1
      Call SetL0(IStep)
      If(IPSave(0).le.1) then
        If(Neg.eq.0) then
          Write(IOut,1090)
        else if(Neg.eq.1) then
          Write(IOut,1100)
        else
          Write(IOut,1110) Neg
          endIf
        endIf
      If(JScan(2).eq.0) then
        Write(IOut,1120) IStep, NStep
      else
        Write(IOut,1130) IStep, NStep, JScan(1), JScan(2)
        endIf
      Write(IOut,1140)
      NNegE = Neg
      Failed = .False.
C
C     Zap elements of gradient corresponding to frozen variables.
C
      Call AMask(NVar,-100000,-1,IFlags,NVarM,F,FTemp)
C     Determine whether simplified projection onto the linearly
C     independent space can be done.
      SimPrj = DoRed.and..not.Large.and..not.LSTQST.and.NVarM.eq.NVar
      NVrMMY = NVarMM + 6*NRF
      NVarT = NVar + NVrMMY
      NDimMM = Max(NVrMMY,1)
      NVrMMX = 0
      Call AClear(NVrMMY,FMMQ)
      IFCase = 0
      IndGI = -1
      NVSq = NVar**2
      NVTT = (NVar*(NVar+1))/2
      IndEV = 0
      If(DoRed.and..not.Large) then
        NRank = -1
        If(IGISav.gt.0) then
          Call RdWrGI(IOut,IPrint,2,NAtNew,NVar,NRank,IGISav,IBMat,BMat,
     $      Scr,A,VecInt,FCTemp,IndGI,IndEV,IZRed)
        else
          If(AlUnit)
     $      Call GauErr('AlUnit but G-Inverse needed in GrdOpt.')
          Call Red2BG(IOut,IPrint,IAprBG,NAtNew,NTrRot,IGISav,IGIOld,C,
     $      NVar,NTBond,NTAng,NTDih,IZRed,Value,NRank,BMat,IBMat,A,
     $      VecInt,FCTemp,Scr,XName,V,MDV)
          endIf
C
C       Check if we're doing QM/MM. If so, process the IMMCRS file.
        NRFX = 0
        If(ITqry(IMMCRS).gt.0.and.OptMic) then
C         IFCase = 1: IntHes
C                  2: QuadMac
C                  3: Regular macro-micro
          If(IntHes) then
            If(IDir.ge.4.or.IDir.lt.2) Call GauErr('IntHes 1<IDir<4')
            If(IGetFC.ne.4) Call GauErr('IntHes needs analytical H')
            IFCase = 1
          else if(IDir.ge.2) then
            IFCase = 2
          else
            IFCase = 3
            endif
          Call FileIO(2,-IMMCRS,1,COrig,0)
          Call BNUPak(COrig,NRFX,IDum)
          Call FileIO(2,IMMCRS,1,EQM,0)
          If(NRF.ne.IAbs(NRFX)) Call GauErr('NRF ne Abs NRFX')
          Call FileIO(2,IMMCRS,6*NRF,XMMT(1,NVarMM/3+1),0)
          Call FileIO(2,IMMCRS,3*NAtoms,COrig,6*NAtoms)
          Call FileIO(2,IMMCRS,3*NAtoms,CNewQ,0)
C         Skip left-over garbage
          Call FileIO(2,IMMCRS,1,V,3*(1+NVarT+2*NVrMMY))
          Call BNUPak(V,NAtmNN,ISecnd)
          If(NAtmNN.lt.0.and.DoEle.and.UpdRWF)
     $      Call GauErr('GrdOpt needs charges')
C
C         Get analytical QM contribution to Hessian in redundant internals.
          Call FileIO(2,IMMCRS,1,V,IAbs(NAtmNN))
          Call BNUPak(V,NVarH,ISecnd)
          If(IGetFC.eq.4.and.(UpdRWF.or.NVarH.gt.0)) then
            If(NVarH.eq.0) Call GauErr('Expected something in IMMCRS')
            NAtomI = ICouQM(NAtoms,MicOpt)
            NVarHA = IAbs(NVarH)
            LenH = ((NVarHA+1)*NVarHA)/2
            IHC = 1
            ICX = IHC + Max(LenH,NNT)
            IFX = ICX + 3*NAtoms
            IFrc = IFX + Max(3*NAtoms,NVarHA,NVar)
            IFor = IFrc + Max(LenH,NNT)
            IEnd = IFor + Max(3*NAtoms,NVarHA,NVar) - 1
            IVX = IEnd + 1
            MDVX = MDV - IEnd
            Call TstCor(IEnd,MDV,'GrdOpt-FrcRed')
            Call FileIO(2,IMMCRS,LenH,V(IHC),0)
            Compct = NVarH.lt.0
            If(.not.Compct) NAtomI = NAtoms
            Call AMove(3*NAtoms,CMM,V(ICX))
            Call ASub(3*NAtoms,COrig,CNewQ,V(IFX))
            DoDipD = .False.
            If(ICrdVr(IZRed).eq.1) then
              Call FrcRed(IOut,IPrint,Compct,DoDipD,NTrRot,IAprBG,
     $          IGISav,IGIOld,NAtoms,NAtomI,MicOpt,NTRedS,IZRedO,V(ICX),
     $          V(IFX),V(IHC),V(IFor),V(IFrc),XXA,Value,XName,V(IVX),
     $          MDVX)
            else
              Call FrcRed(IOut,IPrint,Compct,DoDipD,NTrRot,IAprBG,
     $          IGISav,IGIOld,NAtoms,NAtomI,MicOpt,NTRedS,IZRed,V(ICX),
     $          V(IFX),V(IHC),V(IFor),V(IFrc),XXA,Value,XName,V(IVX),
     $          MDVX)
              endIf
C
C           Now we have the QM contrib, in which cases do we need to use it?
C           Right now: 1) QuadMac-ME and QuadMac-EE when IGetFc.eq.4 .
C           For IntHes, it does not matter what is in Frcnst, because the H
C           will be picked up from rwf 585. But if this IntHes is followed
C           by standard QuadMac steps, only the QM contrib should be in
C           Frcnst. For standard Mac/Mic, Frcnst should always have the full
C           force constants.
C
            If(IFCase.eq.2.or.(IFCase.eq.1.and.IDir0.ge.2))
     $        Call AMove(NNT,V(IFrc),Frcnst)
            endif
C
          Call AMove(NVar,F,FTmp)
          Call ANeg(3*NAtoms,COrig,Scr)
          If(IntHes) Call AMove(3*NAtoms,CNewQ,COrig)
          If(ITqry(IMMRFO).gt.0) Call FileIO(1,-IMMRFO,3*NAtoms,COrig,0)
          If(DoEle) then
            LenPrM = ITqry(IMMMOD)
            IDat   = 1
            IEnd1  = IDat + LenPrM - 1
            IV1    = IEnd1 + 1
            MDV1   = MDV - IEnd1
            Call TstCor(IEnd1,MDV,'GrdOpt')
            Call FileIO(2,-IMMMOD,LenPrM,V(IDat),0)
            IndexQ = 6 + 6*NRF + 12*NAtoms + 3*(NVarT+2*NVrMMY)
            Call FileIO(2,-IMMCRS,NAtmNN,V(IDat+InToWP(7)),IndexQ)
            Call EGHMod(IOut,IPrint,3,Zero,NAtoms,MicOpt,CMM,ESelf,
     $        COrig,XXA,V(IDat),V(IDat),V(IV1),MDV1)
            Call ASub(3*NAtoms,COrig,Scr,COrig)
            Call ASub(3*NAtoms,COrig,CNewQ,COrig)
            endif
          Call ANeg(3*NAtoms,CNewQ,CNewQ)
          IP = 0
          JP = 0
          Do 100 I = 0, (NAtoms-1)
            If(MicOpt(I+1).lt.0) then
              Call AMove1(3,I*3,IP,Scr)
              IP = IP + 3
            else if(MicOpt(I+1).eq.0) then
              Call AMove(3,CNewQ(1,I+1),FTemp(JP*3+NVar+1))
              If(DoEle) Call AMove(3,COrig(1,I+1),FMMQ(JP*3+1))
              Call AMove(3,CMM(1,I+1),XMMT(1,JP+1))
              JP = JP + 1
              endIf
  100       Continue
          Call TranF(.True.,.False.,NVar,NAtoms,JJ,Scr,FTemp,IBMat,BMat,
     $      A,JJ)
          Call AMove(3*NAtoms,COrig,Scr)
          If(NRFX.gt.0) then
            Call FileIO(2,-IMMCRS,3*NAtoms,COrig,2+6*NRF)
            Call FileIO(1,-IMMCRS,3*NAtoms,COrig,2+6*NRF+3*NAtoms)
            Call RigidX(NAtoms,COrig,XMMT(1,JP+1),COM,MicOpt,DT,ICnt,T)
            Call RigCnv(.True.,.False.,NAtoms,MicOpt,T,COrig,CNewQ,COM,
     $        FTemp(JP*3+NVar+1))
            If(DoEle) Call RigCnv(.True.,.False.,NAtoms,MicOpt,T,COrig,
     $        Scr,COM,FMMQ(JP*3+1))
          else
            Call AMove(3*NAtoms,CMM,COrig)
            endIf
          Call AMove(NVrMMY,FTemp(NVar+1),FTmp(NVar+1))
          Call ASub(NVar,F,FTemp,F)
          Call AMove(NVar,FTmp,FTemp)
C         F now has QM contrib. FTemp, FTmp full
          Do 110 I = 1, NVar
  110       APrint(2,I) = FTemp(I)
          If(IDir.ge.2) NVrMMX = NVrMMY
        else
          If(OptMic.and.SetStp) Call GauErr('MOMM, but no IMMCRS')
          Call AMove(NVar,F,FTemp)
          endif
        If(NVarQ.gt.0) then
          Call GauErr('NVarQ>0 in GrdOpt.')
          NVarT = NVarQ + NVar
          Call PCMGet('IRwEfQ',IRwEfQ,RJunk)
          Call FileIO(2,-IRwEfQ,NVarQ,XXQ,0)
          Call AMove(NVar,F,FTmpQ)
          endIf
        If(SimPrj) then
          Call GGEMV('T',NVar,NRank,One,VecInt,NVar,FTemp,1,Zero,Scr,1)
          Call GGEMV('N',NVar,NRank,One,VecInt,NVar,Scr,1,Zero,FTemp,1)
        else
          Call AMove(NVTT,FCTemp,VecInt)
          Call CnstrP(IOut,IPrint,NVar,IFlags,-100000,-1,VecInt,A,
     $      FCTemp,Scr)
          Call AMove(NVar,FTemp,Scr)
          Call MultLV(.True.,1,0,NVar,NVar,VecInt,Scr,FTemp)
          Call FileIO(1,-IGISav,NVTT,VecInt,IndEV+NVSq)
          endIf
        endIf
      RtNV = GFloat(NVarM)
      RtNV = Sqrt(RtNV)
C
C     Test gradient.
C
      FMax = ARRMax(FTemp,NVar,.True.)
      FRMS = SProd(NVar,FTemp,FTemp)
      If(NVarM.ne.0) then
        FRMS = Sqrt(FRMS/NVarM)
      else
        FRMS = Zero
        endIf
      Finish = FRMS.lt.ConvF.and.FMax.lt.(CnvFMx*ConvF)
      If(FMax.ge.FMaxT.and.SetStp) then
        NCur = 0
        Call PutCSt('-DE/DX',Tmp,NCur)
        Call PutDel(2,Tmp,NCur)
        Call MatPrt(FTemp,1,NVar,1,NVar,1,1,Tmp,XName,0,XXA,0)
        Write(IOut,1370) FMaxT
        Call Lnk1E(0)
        endIf
      Call AClear(NVarT,XQuad)
C
C     Run FixDQ on the earlier points.
C
      If(DoRed) then
        Do 120 I = 1, NP
          Call ASub(NVar,XX(1,I),X,XNew)
          Call FixDQ(IOut,IPrint,0,NVar,NTBond,NTAng,NTDih,XNew,IZRed)
          Call AAdd(NVar,X,XNew,XX(1,I))
  120     Continue
        endIf
C
C     Set up for simultaneous optimization.
C
      Call MixOpt(IOut,IPrint,IMDS,MethQ,Large,NP,Neg,SimOpt,Mixed,
     $  IStep,JScan,SwitMx,DoLinS,DoDIIS,FRMS,NGDIIS)
      If(SimOpt) Call GauErr('SimOpt in GrdOpt.')
C
C     Save the current point.
C
      Swap = DoLinS
      If(.not.Swap.and..not.LSTQST) Write(IOut,1480)
      DidSwap = .False.
      Call AMove(NVar,X,XNew)
C     For regular QuadMac, the QM contrib to the force needs to go in SavePt,
C     but when the integrated Hessian is used, the full forces (just like the
C     force constants).
C
C     For QM/MM, at this point FTemp has the full forces, and F has the QM
C     contribution only.
C     For SavePt, F must have the QM contribution only for QuadMac.and..not.IntHes,
C     and the full gradient for all other cases. FTemp will have the oposite.
C     Leaving this point, F must be zero for IntHes, QM contribution for QuadMac,
C     and full for all other cases. FTemp will always be the full gradient.
      If((IFCase.eq.1.and.IDir0.lt.2).or.IFCase.eq.3)
     $  Call SwapAB(NVar,1,NVar,F,NVar,FTmp)
      INeg = Neg
      Call SavePt(Energy,NVar,IGetFC,AlUnit,INeg,Rises,F,Frcnst,ESave,
     $  IC,ITU,X,FF,XX,FC,FC2,Swap,DidSwap,NVarT,NVrMMY,NDimMM,FTmp,
     $  XMMT,FFMM,XXMM,FMMQ,FFMMQ,EQM,EQMS)
      If((IFCase.eq.1.and.IDir0.lt.2).or.IFCase.eq.3)
     $  Call SwapAB(NVar,1,NVar,F,NVar,FTmp)
      If(IFCase.eq.3) Call AMove(NVar,FTmp,F)
      Call ASub(NVar,X,XNew,XNew)
      If(DoRed.and.Swap.and.SProd(NVar,XNew,XNew).ne.Zero) then
        If(.not.Rises) Call GauErr('Should only happen when Rises')
C
C       The "last" point has been altered in SavePt, so we fix the Cartesians.
        Call FixDQ(IOut,IPrint,0,NVar,NTBond,NTAng,NTDih,XNew,IZRed)
        DidGI = IGISav.ne.0.and..not.Large
        If(DidGI) Call FileIO(2,-IGISav,NVSq,A,IndGI)
        Call RedCar(IOut,IPrint,IOSLE1,NVar,NTBond,NTAng,NTDih,IZRed,
     $    Value,NAtNew,C,XNew,CNew,BMat,IBMat,VecInt,FCTemp,Scr,Cnstr1,
     $    IFlags,.True.,IFail,IWght,V,MDV)
        endIf
C     If Large then perform the associated surface approximation for LSTQST.
      If(LSTQST.and.Large.and.DoRed) Call DgQSTL(IOut,IPrint,IQST,NVar,
     $  NP,Cnstr1,Cnstr2,X,FF,Scr,NTBond,NTAng,NTDih)
C
C     Update the Hessian
C
      Call TStamp(1,'Before D2Corr')
      If(DidSwap) then
        NNP = 0
      else
        NNP = 1
        endIf
      If(MthUpd.eq.9.and.DidSwap) then
        Method = 1
      else if(MthUpd.eq.9) then
        Method = 3
      else
        Method = MthUpd
        endIf
      If(IGetFC.eq.4) then
        Write(IOut,1260)
      else if(IGetFC.eq.10) then
        Write(IOut,1500)
      else if(NP.eq.1) then
        Write(IOut,1230)
      else if(IAbs(Method).eq.5) then
        Write(IOut,1250)
      else if(IAbs(Method).eq.6.or.IAbs(Method).eq.7) then
        Call D2CorL(IOut,Method,NVar,NNP,NP,1,RMin,RMax,GrdErr,IC,XX,
     $    FF,FC,Scr(1,1),Scr(1,2),IFlags)
      else if(.not.Large) then
        Call D2Corr(Method,NVar,IC,ESave,XX,FF,FC,Scr(1,1),Scr(1,2),
     $    Scr(1,3),A,Swap)
      else if(Abs(Method).gt.2) then
        If(DidSwap) then
          Call D2Corr(Method,NVar,IC,ESave,XX,FF,FC,Scr(1,1),Scr(1,2),
     $      Scr(1,3),XXA,Swap)
        else
          Call D2CorL(IOut,Method,NVar,NNP,NP,1,RMin,RMax,GrdErr,IC,XX,
     $      FF,FC,Scr(1,1),Scr(1,2),IFlags)
          endIf
      else
        Call GauErr('Illegal update method in GrdOpt.')
        endIf
      Call TStamp(1,'D2Corr')
C
C     Update trust radius if allowed.
C
      If(IPrint.gt.0.or.IPSave(0).le.1) then
        IPrDXM = 1
      else
        IPrDXM = 0
        endIf
      If(UpTrus.and.HavDEP.and.ISStep.eq.0) then
        ECor = Zero
        If(IDir.ge.2) then
C         For Electronic Embedding, ECor is always zero. No problem,
C         but could be fixed up at some point.
          IInd = ITqry(IMMRFO)
          Call FileIO(2,-IMMRFO,1,ECor,IInd-1)
          ESave(1) = ESave(1) - ECor
          endIf
        Call UpdDXM(IOut,IPrDXM,SimOpt,Neg,DEPred,ESave,DXMaxT,NVar,XX,
     $    XX(1,2),One,IUpDwn)
        ITU(1) = IUpDwn
        ESave(1) = ESave(1) + ECor
        endIf
      If((IPrint.ge.1.or.(NVar.le.NVarPt.and.IPSave(0).le.1)).and.
     $   .not.Large.and..not.AlUnit) then
        Write(IOut,1070)
        Call Square(FC,A,NVar,NVar,0)
        Call MatPrt(A,NVar,NVar,NVar,NVar,1,1,XName,XName,1,XXA,0)
        endIf
      Steped = .False.
C
C     First order simultaneous optimization
C
      OKHess = CkHesU(IOut,IPrint,NP,ITU)
      SimSwU = SimSw
      If(OKHess) then
        SqrtH = .False.
        If(NP.ge.4) SimSwU = F100*SimSw
        endIf
      SqrtH = .False.
      SimMix = SimOpt.or.Mixed
      SimMMM = .False.
      OldGDI = IGEVer.eq.1.or.IMDS.eq.1
      If(SimMix) then
        If(OldGDI) then
          Call IniFSO(IRwSM,NMax,MaxGDI,NP,NGDIIS,ESave,SimOpt,MaxSim,V,
     $      MDV)
        else
          Call IniGDI(IRwSM,NMax,MaxGDI,NP,ESave,SimOpt,V,MDV)
          endIf
      else
        DoDIIS = IMDS.ne.0.and.NP.ge.2
        endIf
      If(SimMMM) then
        IRwSMM = IRwSM
      else
        IRwSMM = 0
        endIf
C
C     Linear step.  This is only done if we are allowed to do it,
C     at least two points have been computed, and either the
C     search is for a local minimum or there are analytic second
C     derivatives at both points and the gradient has increased.
C
C     If GDIIS failed we come back to this point
C
  130 Steped = .False.
      DidMix = .False.
      Call AClear(NVarT,XNew)
      Call AClear(NVarT,XQuad)
      If(.not.DoDIIS.and.IMDS.ne.0.and.NP.ge.2.and.DoLinS)
     $  Write(IOut,1490)
      DoLin  = DoLinS.and..not.DoDIIS
      LinMin = Neg.eq.0.and.(DoLin.or.(MethQ.eq.1.and.Rises))
      LinTS = Neg.ne.0.and.MethQ.eq.4.and.Rises.and.
     $  IPSav1.ne.0.and.IPSav2.ne.0
      Linr = NP.gt.1.and.FRMS.gt.FSwtch.and.(LinMin.or.LinTS)
      Quad = MethQ.ne.1
      BigLin = .False.
C     Now FTemp will contain QM contribution only!
C     IntHes: FTemp will contain the full gradient.
      If(IDir.ge.2.and..not.(IntHes.or.DoEle)) Call AMove(NVar,F,FTemp)
      What(1) = Zero
      LMMRFO = Itqry(IMMRFO)
      If(LMMRFO.gt.0) Call FileIO(2,-IMMRFO,1,What,LMMRFO-6)
      DdTSLn = What(1).gt.Pt5
      DidLin = .False.
      If(Linr.or.DdTSLn) then
        MaxL = 0
        Negx = Neg
        If(DdTSLn) then
          MaxL = 50
          Negx = 0
          endIf
        If(IFCase.ne.0.and..not.DdTSLn.and.ONLnOK) then
          If(ICrdVr(IZRed).ne.1)
     $      Call GauErr('GrdOpt: New IC data structure NYI for DXLOnm.')
          Call DXLOnm(IOut,IPrint,Rises,NVar,NVrMMY,NDimMM,NAtoms,
     $      NAtNew,CNew,XX,FF,XXMM,XNew,NTBond,NTAng,NTDih,IZRed,BMat,
     $      IBMAt,Scr,Cnstr1,IFlags,IWght,NRF,COrig,MicOpt,COM,DT,T,
     $      ICnt,RotQM,IAtTyp,AtChMM,IRwMM,IMMMOD,EQMS(1),EQMS(2),NVarT,
     $      FFMM,XLin,NTRRot,DidLin,DXRMS,DXMax,FTemp,IntHes,IFCase,
     $      DoEle,FFMMQ,IMMCRS,NAtmNN,IndexQ,IDir0,XName,IAprBG,V,MDV)
          If(DidLin) ICase = 6
          endif
        If(.not.DidLin) then
C         If IntHes.eq.True, we need to do a linear search in the full space.
C         Important: the gradient is then modified, which needs to go
C         correctly into QIMOMM
          Call DXLinr(ICase,Negx,FTemp,XNew,DXRMS,DXMax,NVar,DidLin,
     $      XLin,ESave(1),FS,FC,FC2,XX,FF,Scr(1,1),Scr(1,2),Scr(1,3),
     $      FFMM,XXMM,NVrMMY,NVarT,NDimMM,MaxL,IFCase,DoEle,FFMMQ,IDir0)
          endif
        If(DdTSLn.and.DidLin) then
          Write(IOut,1520) XLin
          Quad = .False.
        else if(DdTSLn) then
          Write(IOut,1530)
          endIf
        Finish = Finish.and..not.(DdTSLn.and.DidLin)
        If(DidLin) then
          Steped = Steped.or.Abs(XLin).gt.RMinLn
          DX = OptDX(DoRed,Large,NVar,NTRedS,XNew)
          BigLin = DX.gt.DXMaxT.and.IRisen.eq.0
          If(BigLin) then
            If(DX.gt.(Two*DXMaxT)) then
              Scale = Two*DXMaxT/DX
              Write(IOut,1220) DX, DXMaxT, Scale
              Call AScale(NVarT,Scale,XNew,XNew)
              XLin = XLin*Scale
              DXMax = DXMax * Scale
              DXRMS = DXRMS * Scale
            else
              Write(IOut,1540) DX, DXMaxT
              endIf
            Quad = Quad.and..not.Swap
          else
            DXMax = Zero
            If(Rises) then
              If(LinTS) then
                Write(IOut,1350)
              else if(Quad.and.LinMin) then
                If(IMDS.eq.0) then
                  Write(IOut,1330) 'Quadratic'
                else
                  Write(IOut,1330) 'Quadratic/GDIIS'
                  endIf
              else if(.not.Quad.and.Neg.eq.0) then
                Write(IOut,1340)
                endIf
              endIf
            endIf
        else
          DXMax = Zero
          BigLin = .False.
          endIf
        endIf
      Call TStamp(1,'Aft Linear')
C
C     Quadratic contributions to displacements.  These are done unless
C     things have gone wrong AND a linear search has been taken which
C     should fix the situation; that is unless the Energy has risen
C     when searching for a minimum or the forces have gone up when
C     searching for a saddle point and the conditions for a linear
C     search have been met (vide supra).  Use Scr(,2) for forces in
C     the eigenvector basis and Scr(,1) for the eigenvalues.
C
      DoSD = Neg.eq.0.and.Linr.and..not.DoDIIS.and..not.DidLin.and.Rises
      LnOnly = Rises.and.(DidLin.or.DoSD)
      Quad = Quad.and..not.LnOnly.and..not.DoDIIS
      DidSph = .False.
      SmlStp = .False.
      What(1) = Zero
      LMMRFO = Itqry(IMMRFO)
      If(LMMRFO.gt.0) Call FileIO(1,-IMMRFO,1,What,LMMRFO-6)
      QCoupl = IDir.ge.2.or.NVarQ.gt.0
      If(Quad.and.QCoupl) then
        If(NVarQ.gt.0) then
          If(IDir.ge.2) Call GauErr(
     $      'No QM/MM quadratic coupled step with NVarQ>0.')
C         Call QIMOQ(...)
          Call GauErr('QIMOQ')
          endIf
        If(IDir.ge.2) then
          If(ICrdVr(IZRed).ne.1)
     $      Call GauErr('GrdOpt: New IC data structure NYI for DXLOnm.')
          Call QIMOMM(IOut,IPrint,NAtoms,NatNew,C,CNewQ,NVar,NTBond,
     $      NTAng,NTDih,IZRed,X,XNew,Value,IBMat,Scr,Cnstr1,IFlags,
     $      IWght,RotQM,NRF,T,MicOpt,COM,COrig,DT,ICnt,XMMT,NVarMM,
     $      IMMMOD,HMod,HModO,CMM,IRWMM,IAtTyp,AtChMM,NumTpS,IAtTpS,
     $      MaxLEF,EField,FTp,IMMRFO,NTrRot,BMat,IDBMat,DBMat,FTemp,A,
     $      ISCRF,Steped,IGuess,IRFO,FC,XQuad,MethQ,IDir,Neg,IHVAcc,
     $      TstEig,XName,IZRedO,FCTemp,.False.,SmlStp,IGetTS,IPFlag,
     $      AllowP,FMM,FMFlag,FMFlg1,NFxFlg,IHMeth,IntHes,DoEle,IRwF2,
     $      HReal,IMMCRS,NAtmNN,V,MDV)
          SkpMic = SkpMic.or.SmlStp
          If(SmlStp) then
            What(1) = One
            If(LMMRFO.gt.0) Call FileIO(1,-IMMRFO,1,What,LMMRFO-6)
            Finish = .False.
            Write(IOut,1510)
            Call AClear(NVarT,XNew)
            endIf
          endIf
      else if(Quad.and..not.Large) then
        FixEig = MethQ.ne.2.and.MethQ.ne.3
        SclEig = IGetFC.ne.4.and.FixEig
        If(DidDgF) then
          Call AMove(NVar,EigSav,Scr)
          Call MatMP1(1,0,1,1,1,NVar,NVar,FTemp,A,Scr(1,2))
        else
          If(IGISav.gt.0.and..not.SimPrj)
     $      Call FileIO(2,-IGISav,NVTT,VecInt,IndEV+NVSq)
          Call DiagFC(IOut,IPrint,FixEig,SclEig,FCart,AlUnit,NTrRot,
     $      NVar,EigMin,EigMax,IFlags,FC,FCTemp,XName,A,IDimSc,Scr,
     $      SimPrj,NRank,VecInt,FTemp,Neg,NNegE,NEig,NAtNew,C,AtmChg,
     $      NFroz,Frozen,Wght,IQST,IStep,DoRed,X,Cnstr1,Cnstr2,NTBond,
     $      NTAng,NTDih,V,MDV)
          endIf
        If(MethQ.eq.3.or.(MethQ.eq.2.and.NNegE.ne.Neg)) then
          Call DXRFO(IOut,IPrint,AlUnit,NGoDwn,NVar,Neg,NEig,DXMaxT,
     $      IEStpM,Scr(1,1),A,Scr(1,2),XQuad,OK,DidEig)
        else
          Call DXQuad(IPrint,Sphere,AlUnit,FTemp,XQuad,A,OK,NVar,NEig,
     $      Scr(1,2),Scr(1,1),Scr(1,3))
          DidSph = Sphere
          endIf
        If(FCart) Call ClnCor(IOut,IPrint,.False.,DoClnS,NAtoms,IAn,
     $    AtmChg,C,0,Scr,.True.,TRot,XQuad,Scr)
        Steped = Steped.or.OK
      else if(Quad.and.Large) then
        ISh = NNT + NAt3TT
        RLamR = RLam
        If(DidLin) RLamR = Zero
        IHX = Nat3TT+1
        IISB = NVar+1
        ISB = NVar+1
        If(ICrdVr(IZRed).ne.1)
     $    Call GauErr('GrdOpt: New IC data structure NYI for TrH2Hx.')
        Call TrH2Hx(IOut,IPrint,NTBond,NTAng,NVar,NAtNew,IZRed,C,
     $    Value,IBMat,BMat,IBMat(1,IISB),BMat(1,1,ISB),IFlags,FC,
     $    FCTemp,Scr,LenSB,LenCB,V,MDV)
        Call AMove(NAt3TT,FCTemp,FCTemp(IHX))
        Call FileIO(2,-IOSLEq,1,IIniNR,ISh)
        NG = 0
        If(.not.InitNR) then
          Call FileIO(2,IOSLEq,NAt3TT,FCTemp,0)
          Call FileIO(2,IOSLEq,1,RNG,0)
          Call IMove(1,RNG,NG)
          Call FileIO(2,IOSLEq,NG*NAt3,VecInt,0)
          endIf
C
C       FCTemp has also room for scratch vectors of length 2*NAt3
C
        Call NRStep(IOut,IPrint,InitNR,.True.,ConvEq,RLamR,ELow,NVar,
     $    FCTemp(IHX),FCTemp,FTemp,XQuad,VecInt,NG,NAtNew,IBMat(1,IISB),
     $    BMat(1,1,ISB))
        If(.not.DidLin) RLam = RLamR
        Call FileIO(1,-IOSLEq,1,IIniNR,ISh)
        Call FileIO(1,IOSLEq,NAt3TT,FCTemp,0)
        RNG = Zero
        Call IMove(1,NG,RNG)
        Call FileIO(1,IOSLEq,1,RNG,0)
        Call FileIO(1,IOSLEq,NG*NAt3,VecInt,0)
        Steped = .True.
      else if(DoDIIS) then
        Quad = .True.
        If(OldGDI) then
          If(.not.SimOpt) NGDIIS = Min(NGDIIS+1,NP)
          If(SimMix) NGDIIS = Min(NGDIIS,MaxSim)
          NGDIIS = Max(NGDIIS,2)
          If(SimOpt) then
            Call BNPack(MaxGDI,NGDIIS,Temp)
            Call FileIO(1,-IRwSM,1,Temp,0)
            endIf
          endIf
C       Description for MDSOpt is in DoLMDS and GDIIS
        If(Large) then
          I00 = 400
        else
          I00 = 300
          endIf
        If(IMDS.eq.1) then
          If(MthUpd.eq.5.or.Finish) then
            MDSOpt = 2001 + I00
          else if (Large) then
            MDSOpt = 4001 + I00
          else
            MDSOpt = 3001 + I00
            endIf
        else if(IMDS.lt.0) then
          MDSOpt = -IMDS
        else
          MDsOpt = 1001 + I00
          endIf
        If(Large) then
C         Read the guess inverse for the Hessian
          ISh = NNT + NAt3TT
          IHX = Nat3TT+1
          IISB = NVar+1
          ISB = NVar+1
          If(ICrdVr(IZRed).ne.1)
     $      Call GauErr('GrdOpt: New IC data structure NYI for TrH2Hx.')
          Call TrH2Hx(IOut,IPrint,NTBond,NTAng,NVar,NAtNew,IZRed,C,
     $      Value,IBMat,BMat,IBMat(1,IISB),BMat(1,1,ISB),IFlags,FC,
     $      FCTemp,Scr,LenSB,LenCB,V,MDV)
          Call AMove(NAt3TT,FCTemp,FCTemp(IHX))
          Call FileIO(2,-IOSLEq,1,IIniNR,ISh)
          NG = 0
          RNG = Zero
          Call IMove(1,NG,RNG)
          If(.not.InitNR) then
            Call FileIO(2,IOSLEq,NAt3TT,FCTemp,0)
            Call FileIO(2,IOSLEq,1,RNG,0)
            Call IMove(1,RNG,NG)
            Call FileIO(2,IOSLEq,NG*NAt3,VecInt,0)
            endIf
C
C         Call the improved GDIIS routines for multidimensional search
C         with large molecular option
C
          If(SimOpt) then
            Call GauErr('SimOpt and Large in GrdOpt.')
          else if(Mixed) then
            Call GauErr('Mixed and Large in GrdOpt.')
          else
            Call DoLMDS(IOut,IPrint,MDSOpt,NGDIIS,NUsed,NVar,FCTemp,XX,
     $        FF,XNew,XQuad,VecInt,ScrGDI,InitNR,OK,RLam,ELow,ConvEq,NG,
     $        NAtNew,IBMat(1,IISB),BMat(1,1,ISB),IFlags)
            endIf
          Call FileIO(1,-IOSLEq,1,IIniNR,ISh)
          Call FileIO(1,IOSLEq,NAt3TT,FCTemp,0)
          RNG = Zero
          Call IMove(1,NG,RNG)
          Call FileIO(1,IOSLEq,1,RNG,0)
          Call FileIO(1,IOSLEq,NG*NAt3,VecInt,0)
          DoDIIS = NUsed.gt.2.or..not.DoLinS
          If(.not.DoDIIS) goto 130
          Steped = .True.
        else
          FixEig = .False.
          SclEig = IGetFC.ne.4.and..not.SimOpt
          If(IGISav.gt.0.and..not.SimPrj)
     $      Call FileIO(2,-IGISav,NVTT,VecInt,IndEV+NVSq)
          Call DiagFC(IOut,IPrint,FixEig,SclEig,FCart,AlUnit,NTrRot,
     $      NVar,EigMin,EigMax,IFlags,FC,FCTemp,XName,A,IDimSc,Scr,
     $      SimPrj,NRank,VecInt,FTemp,Neg,NNegE,NEig,NAtNew,C,AtmChg,
     $      NFroz,Frozen,Wght,IQST,IStep,DoRed,X,Cnstr1,Cnstr2,NTBond,
     $      NTAng,NTDih,V,MDV)
C
C         Do SimOpt or call the GDIIS routines for multidimensional search.
C
          NVarU = NVar
          NVMMTU = 0
          If((SimOpt.or.Mixed).and..not.OldGDI) then
            Call GauErr('Logic error in GrdOpt.')
            OK = DidDIS
            If(.not.DidDIS) then
              MethQ = 3
              DoLinS = .True.
              DidDgF = .True.
              Call AMove(NVar,Scr,EigSav)
C
C             If the Hessian has a wrong number of negative eigenvalues,
C             reset DIIS vector counter.
C
              If(NNegE.ne.Neg) NGDIIS = -1
              endIf
          else if(Mixed.and.OldGDI) then
            If(NNegE.eq.Neg.and..not.(IGDRFO.eq.5.and..not.OKHess)) then
              Call DoMix(IOut,IPrint,IRwSM,MDSOpt,IScMMF,NGDIIS,NMax,
     $          MaxGDI,NVarU,NVar,NVMMTU,NDimMM,Neg,NEig,Scr(1,1),A,XX,
     $          FF,XXMM,FFMM,IC,XNew,XQuad,ScrGDI,OK,SimSw,FncErr,
     $          IGDRFO,OKHess,DidBck,V,MDV)
              DidMix = OK
              If(.not.OK) then
                MethQ = 3
                DoLinS = .True.
                DidDgF = .True.
                Call AMove(NVar,Scr,EigSav)
                endIf
            else
              OK = .False.
              MethQ = 3
              NGDIIS = -1
              DoLinS = .True.
              DidDgF = .True.
              Call AMove(NVar,Scr,EigSav)
              endIf
          else
            Call DoMDS(IOut,IPrint,MDSOpt,NGDIIS,NUsed,NVar,NEig,Neg,
     $        Scr(1,1),A,XX,FF,XNew,XQuad,ScrGDI,OK)
            endIf
          DoDIIS = OK.and.(NGDIIS.gt.2.or..not.DoLinS)
          If(.not.DoDIIS) goto 130
          Steped = Steped.or.OK
          If(.not.Steped) then
            Call DXRFO(IOut,IPrint,AlUnit,NGoDwn,NVar,Neg,NEig,DXMaxT,
     $        IEStpM,Scr(1,1),A,Scr(1,2),XQuad,OK,DidEig)
            Steped = OK
            endIf
          endIf
        DXRMS = SProd(NVar,XQuad,XQuad)
        DX = DXRMS * RtNV
        If(FCart) Call ClnCor(IOut,IPrint,.False.,DoClnS,NAtoms,IAn,
     $    AtmChg,C,0,Scr,.True.,TRot,XQuad,Scr)
        endIf
C
C     Print reason for skipping linear search here as DXQuad prints
C     the eigenvalues of second derivative matrix.
C
      If(Linr.and.DidLin) then
        If(ICase.eq.1) then
          Write(IOut,1390) XLin
        else if(ICase.eq.2) then
          Write(IOut,1400) XLin
        else if(ICase.eq.3) then
          Write(IOut,1410) XLin
        else if(ICase.eq.6) then
          Write(IOut,1415) XLin
        else if(ICase.ne.5) then
          Write(IOut,1460) ICase
          Call Lnk1E(0)
          endIf
      else if(Linr) then
        Write(IOut,1420)
      else
        If(.not.DoLin) then
          If(.not.DoDIIS) Write(IOut,1270)
        else if(NP.le.1) then
          Write(IOut,1290)
        else if(IDir.ge.2.and.NP.le.2) then
          Write(IOut,1300)
        else if(Neg.ne.0.and.Rises.and.(IPSav1.eq.0.or.IPSav2.eq.0))
     $  then
          Write(IOut,1280)
        else if(.not.Rises.and.(MethQ.eq.1.or.Neg.ne.0)) then
          Write(IOut,1320)
        else if(FRMS.le.FSwtch) then
          Write(IOut,1310) FSwtch
        else
          Call GauErr('Linear search skipped for unknown reason.')
          endIf
        endIf
C
C     If both the linear and quadratic steps have failed, step in
C     the direction of the gradient.  Test size of displacement
C     and update position vector values stored in a for printing.
C
      Despar = DoSD.and..not.Steped
      If(.not.DdTSLn.and.Neg.eq.0.and.(.not.Steped.or.
     $  (.not.Rises.and..not.Quad.and..not.BigLin))) then
        Write(IOut,1430)
        If(IDir.le.1) then
          Call DXGrad(NVar,FTemp,XQuad,DXRMS,DXMax,OK)
          OK1 = .True.
        else
          Call DXGrad(NVar,FTemp,XQuad,DXRMS,DXMax,OK)
          Call DXGrad(NVrMMX,FTemp(NVar+1),XQuad(NVar+1),DXRMSM,DXMaxM,
     $      OK1)
          endIf
        Steped = Steped.or.(OK.and.OK1)
        Finish = .False.
        endIf
      Call AMask(NVar,-100000,-1,IFlags,NVarM,XQuad,XQuad)
C
C     SimOpt, when the error in density is large, fixing nuclei
C
      NVarX = NVar + NVrMMX
      If(DoDIIS) then
        Call AAdd(NVarT,XNew,XQuad,Scr)
      else
        Call AMove(NVarT,XQuad,Scr)
        endIf
C     DXMaxT is updated based on the internal coordinate Hessian and so
C     should be applied based on just the internal coord part of the step.
      DX = OptDX(DoRed,Large,NVar,NTRedS,Scr)
      DXMaxU = DXMaxT
      If(DidDIS.and.SqrtH) DXMaxU = DXMaxU*Sqrt(Two)
      If(DX.gt.DXMaxU.and..not.DidSph.and..not.DidBck) then
        Scale = DXMaxU / DX
        Write(IOut,1360) DXMaxU
        If(DidEig.and.DoRed) then
          Write(IOut,1363)
        else
          Write(IOut,1365) Scale
          Call AScale(NVarT,Scale,XQuad,XQuad)
          If(DoDIIS) Call AScale(NVar,Scale,XNew,XNew)
          endIf
        endIf
      DXQR = OptDX(DoRed,Large,NVar,NTRedS,XQuad)/Sqrt(GFloat(NVar))
      DXQMax = OptDXM(DoRed,NVar,NTRedS,XQuad)
      DXTest = Cnvfx*ConvF
      If(Despar.and.(.not.Finish.or.DXQR.gt.DXTest.or.
     $  DXQMax.gt.(CnvFMx*Cnvfx*ConvF))) then
        If(SteepM.gt.Zero) then
          Write(IOut,1440) SteepM
          Call AScale(NVar,SteepM/DXQMax,XQuad,XQuad)
        else if(DXQR.lt.(-SteepM)) then
          Write(IOut,1450) -SteepM
          Call AScale(NVar,-SteepM/DXQR,XQuad,XQuad)
          endIf
        endIf
      Do 140 I = 1, NVar
        APrint(1,I) = X(I)
        If(IDir.le.1) then
          APrint(2,I) = F(I)
        else
          APrint(2,I) = FTmp(I)
          endIf
        APrint(3,I) = XNew(I)
        APrint(4,I) = XQuad(I)
        XNew(I) = XNew(I) + XQuad(I)
  140   X(I) = X(I) + XNew(I)
      Do 150 I = (NVar+1), NVarX
        APrint(3,I) = XNew(I)
        APrint(4,I) = XQuad(I)
  150   XNew(I) = XNew(I) + XQuad(I)
      IPrt1 = IPrint
      If(IPSave(0).le.1) IPrt1 = Max(IPrint,1)
      IStDih = NTBond + NTAng + 1
      Call TStamp(1,'Bef Red2X')
      Just1 = .False.
      If(DoRed) then
        DidGI = IGISav.ne.0.and..not.Large
        If(DidGI) Call FileIO(2,-IGISav,NVSq,A,IndGI)
        If(IWght.eq.2) then
          Call AMove(NVar,ScrGDI,VecInt)
        else if(IWght.eq.3) then
          Do 160 I = 1, NVar
  160       VecInt(I) = Max(Sqrt(Abs(FC((I*(I+1))/2))),ScMin)
          endIf
        Call AMove(NVar,XNew,Scr(1,5))
        If(IOldRC.eq.1) then
          ScFact = Two
          DCOld = Big
          Do 200 ITry = 1, MxITry
            Call AMove(3*NAtoms,C,CNew)
            If(DidEig) then
              ScFact = Pt5*ScFact
            else
              ScFact = GFloat(MxITry-ITry+1)/GFloat(MxITry)
              endIf
            Do 170 I = 1, NVar
  170         XNew(I) = APrint(3,I) + ScFact*APrint(4,I)
            Call FixDq(IOut,IPrint,0,NVar,NTBond,NTAng,NTDih,XNew,IZRed)
            Do 180 I = 1, NVar
  180         X(I) = APrint(1,I) + XNew(I)
            Call FixAng(NTAng,X(NTBond+1))
            Call FixDih(NTDih,X(IStDih))
            Do 190 I = 1, NVar
  190         XNew(I) = X(I) - APrint(1,I)
            Call RedCar(IOut,IPrt1,IOSLE1,NVar,NTBond,NTAng,NTDih,IZRed,
     $        Value,NAtNew,C,XNew,CNew,BMat,IBMat,VecInt,FCTemp,Scr,
     $        Cnstr1,IFlags,.True.,IFailR,IWght,V,MDV)
            IFail = IFailR
            Call TstCor(NAt3,MDV,'GrdOpt-DeltaC')
            Call ASub(NAt3,CNew,C,V)
            DXMaxC = ArrMax(V,NAt3,.True.)
            If(IFail.eq.1.or.ITry.gt.1.or.IPSave(0).eq.0)
     $        Write(IOut,1550) ITry, IFail, DXMaxC, DCOld, DXMaxU,
     $        DXLimC, Rises
            If(DXMaxC.gt.Max(DXLimC,DXMaxU).and..not.Rises.and.
     $        ITry.lt.MxITry) IFail = 1
            If(DXMaxC.gt.DCOld) then
              IFail = 1
            else if(IFailR.eq.0) then
              DCOld = DXMaxC
              endIf
            Just1 = .True.
            If(IFail.eq.0) goto 210
            Call TStamp(1,'Aft RedCar')
  200       Continue
          IOldRU = IOldRB
        else
          IFail = 1
          IOldRU = IOldRC
          endIf
  210   If(IFail.eq.1) then
          If(Large)
     $      Call GauErr('Odon got you.  No ORedCr with Opt=Large.')
          Call AMove(3*NAtoms,C,CNew)
          Call AMove(NVar,Scr(1,5),XNew)
          Call ORedCr(IOut,IPrt1,IOldRU,Just1,IAprBG,IGISav,IGIOld,
     $      NTrRot,NVar,NTBond,NTAng,NTDih,IZRed,Value,DXMaxU,NAtNew,C,
     $      XNew,CNew,BMat,IBMat,A,VecInt,FCTemp,Scr,Cnstr1,IFlags,
     $      XName,IFail,V,MDV)
          DoMatC = NPDir.eq.0.and.NTrRot.eq.-1
          endIf
        If(IFCase.ne.0) Call UpCOrg(CNew,CNewQ,COrig,NAtoms,NAtNew,
     $    MicOpt,XMMT,COM,DT,ICnt,T,NRF,NVar,NVarMM,XNew,IMMCRS,RotQM,
     $    .True.)
      else if(FCart) then
        Call AMove(NAt3,X,CNew)
        endIf
      Call TStamp(1,'Bef ClnCor')
      If(DoRed.or.FCart) then
        Call ClnCor(IOut,IPrint,.False.,DoClnS,NAtoms,IAn,AtmChg,C,1,
     $    TrVec,.True.,TRot,CNew,Scr)
        If(DoMatC) then
          Call MatCor(IOut,IPrt1,.True.,NAtoms,C,CNew,TrRMat,Scr,OK)
          If(.not.OK) Write(IOut,1470)
        else
          Call AMove(NAtoms*3,CNew,Scr)
          endIf
        Call ASub(NAtoms*3,C,Scr,C)
        DXMaxC = ArrMax(C,NAt3,.True.)
        DXRMSC = Sqrt(Abs(SProd(NAt3,C,C)/GFloat(NVar)))
        Call AMove(NAtoms*3,Scr,C)
        endIf
      If(DoRed) then
        Call RedVal(IOut,IPrint,.True.,C,IZRed,X,NVar,NTBond,NTAng,
     $    NTDih,V,MDV)
      else if(FCart) then
        Call AMove(NAt3,C,X)
        endIf
      Do 250 I = 1, NVar
  250   XNew(I) = X(I) - APrint(1,I)
      If(DoRed) Call FixDq(IOut,IPrint,0,NVar,NTBond,NTAng,NTDih,XNew,
     $  IZRed)
      DXMax = Zero
      DXRMS = Zero
      Do 260 I = 1, NVar
        X(I) = APrint(1,I) + XNew(I)
        If(DoRed.and.I.ge.IStDih) Call FixDih(1,X(I))
        APrint(5,I) = XNew(I)
        APrint(6,I) = X(I)
        DXRMS = DXRMS + XNew(I)**2
        If(Abs(XNew(I)).gt.DXMax) DXMax = Abs(XNew(I))
  260   Continue
      DXRMS = Sqrt(DXRMS/NVar)
      If(.not.DoRed.and..not.FCart) then
        DXMaxC = DXMax
        DXRMSC = DXRMS
        endIf
C
C     Predict change in energy
C
      If(AlUnit) then
        DE = SProd(NVar,XNew,XNew)/Two
      else if(IDir.ge.2) then
C       Get the product of the Hessian with the displacement in the same
C       way as in the coupled step.
        IParam = 1
        LenPar = ITqry(IRwMM)
        IDip   = IParam + LenPar
        IV     = IDip + 3
        IEnd1  = IV - 1
        MDV1   = MDV - IEnd1
        Call TstCor(IEnd1,MDV,'GrdOpt-EGHMM')
        Call FileIO(2,-IRwMM,LenPar,V(IParam),0)
        If(IntHes) then
          Call FileIO(2,-IMMRFO,3*NAtoms,Scr,0)
        else
          Call EGHMM(IOut,IPrint,0,1,0,Zero,NAtoms,IAtTyp,AtChMM,CMM,
     $      MaxLEF,EField,V(IParam),LenPar,JJ,ETemp,Scr,XXA,V(IDip),XXA,
     $      V(IV),MDV1)
          endif
        IGM = 1
        IHRed2 = IGM + NVSq
        IHDiag = IHRed2 + NVSq
        IndI = IHDiag + NVarT
        LenI = 13+2*NAtoms+NVar*4
        IndR = IndI + InToWP(LenI)
        MDVX = MDV - IndR
        Call FileIO(2,-IGISav,NVSq,V(IGM),IndGI)
        J = 1
        Do 270 I = 1, NAtoms
          If(MicOpt(I).lt.0) then
            Call AMove(3,CMM(1,I),CNewQ(1,J))
            J = J + 1
            endIf
  270     continue
        Call FormDB(IOut,IPrint,IRWCrd,NAtoms,NTBond,NTAng,NVar,IZRed,
     $    CNewQ,Value,IBMat,BMat,IDBMat,DBMat,V(IV),MDV1)
        Call FormB(IOut,IPrint,NAtoms,NTBond,NTAng,NVar,IZRed,CNewQ,
     $    Value,IBMat,BMat,V(IndR),MDVX)
        If(NRF.gt.0) then
          Call FileIO(2,-IMMCRS,3*NAtoms,COrig,2+3*NAtoms+6*NRF)
          Call RigidX(NAtoms,COrig,XMMT(1,(NVarMM/3)+1),COM,MicOpt,DT,
     $      ICnt,T)
          endIf
        Call MMDav0(IOut,IPrint,Zero,NAtoms,IAtTyp,V(IndI),V(IndR),LenR,
     $    LenI,IRWMM,AtChMM,MicOpt,IDir,NVar,IBMat,BMat,V(IGM),FC,
     $    VecInt,IFlags,FTemp,CMM,DBMat,0,0,IndGrd,HModO,IMMRFO,T,COrig,
     $    COM,Scr,0,IPrjTS,V(IHRed2),.False.,XXA,IntHes,IRwF2,V(IHDiag),
     $    HReal,0,MDVX)
        NVRdMM = NVarT
        IAX = IndR + LenR
        ICC = IAX + NVRdMM + 1
        IVX = ICC + NVRdMM + 1
        MDVX = MDV - IndR - LenR - 2*(NVRdMM+1)
        Call AClear(NVRdMM+1,V(IndGrd+IndR-1))
        Call AMove(NVRdMM,XNew,V(ICC))
        V(ICC+NVRdMM) = Zero
        Call MMDav1(NVRdMM+1,1,V(IndI),V(IndR),XNew,V(IAX),V(IVX),MDVX)
        DE = SProd(NVRdMM,V(IAX),XNew)/Two - SProd(NVRdMM,FTmp,XNew)
      else
        DE = Zero
        IJ = 0
        Do 290 I = 1, NVar
          Do 280 J = 1, I
            IJ = IJ + 1
  280       DE = DE + XNew(I)*FC(IJ)*XNew(J)
  290     DE = DE - XNew(I)*FC(IJ)*XNew(I)/Two
        DE = DE - SProd(NVar,F,XNew)
        endIf
      DE = DE / FncCnv
C
C     Compute the predicted change in Energy for trust radius update.
C     Care is needed for two reasons:
C     1)  The change in Energy computed at the next step will use the
C         lowest Energy point and if the Energy went up this time that
C         won't be the now current point.
C     2)  If the curvature is wrong and we switched the sign on an
C         eigenvalue of the Hessian, the Energy change computed with the
C         real Hessian may not even have the sign we expect after the
C         sign change.
C     If we're not going for a minimum there's no easy solution but if
C     we are looking for a minimum, DEPred is made <0 so that steps down
C     in Energy will look good when updating the trust radius.
C
      DEPred = DE
      If(IRisen.ne.0) DEPred = ESave(1) + DE - ESave(2)
      If(Neg.eq.0) DEPred = -Abs(DEPred)
      HavDEP = Quad
      ConvX = CNVFX*ConvF
      Finish = Finish.and.(NNegE.eq.Neg.or..not.TstEig).and..not.DidEig
      RegFin = DXRMSC.lt.ConvX.and.DXMaxC.lt.CnvFMx*ConvX
      SmlFin = (FRMS.lt.(ConvF/F100)).and.(FMax.lt.(CnvFMx*ConvF/F100))
      Finish = Finish.and.(RegFin.or.SmlFin)
C
C     If we do quadratic QM/MM and micro-iterations for a transition
C     state search, we re-evalutate the TS vector for use in the
C     micro-iterations.
C
      If(IDir.ge.2.and.Neg.gt.0.and.OptMic.and..not.SkpMic) then
        Call AClear(NVarT,XNew)
        SmlStp = .False.
        If(ICrdVr(IZRed).ne.1)
     $    Call GauErr('GrdOpt: New IC data structure NYI for QIMOMM.')
        Call QIMOMM(IOut,IPrint,NAtoms,NatNew,CNew,CNewQ,NVar,NTBond,
     $    NTAng,NTDih,IZRed,X,XNew,Value,IBMat,Scr,Cnstr1,IFlags,IWght,
     $    RotQM,NRF,T,MicOpt,COM,COrig,DT,ICnt,XMMT,NVarMM,IMMMOD,HMod,
     $    HModO,CMM,IRWMM,IAtTyp,AtChMM,NumTpS,IAtTpS,MaxLEF,EField,FTp,
     $    IMMRFO,NTrRot,BMat,IDBMat,DBMat,FTemp,A,ISCRF,Steped,IGuess,
     $    IRFO,FC,XQuad,MethQ,IDir,Neg,IHVAcc,TstEig,XName,IZRedO,
     $    FCTemp,.True.,SmlStp,IGetTS,IPFlag,AllowP,FMM,FMFlag,FMFlg1,
     $    NFxFlg,IHMeth,IntHes,DoEle,IRwF2,HReal,IMMCRS,NAtmNN,V,MDV)
        endIf
C
C     Test convergence of MM microiteration part.
C
      If(OptMic) then
        Call GetNRF(NAtoms,NNRF,LogTmp,MicOpt)
        ConMM1 = ConvF
        If(IDir.lt.2) ConMM1 = ConMM1 / GFloat(10)
        Call FileIO(2,-IRWF1,3*NAtoms,FTemp,0)
        Loptx = 0
        IC1 = 1
        IC2 = IC1 + 3*NAtoms
        IVC = IC2 + 6*NNRF
        MDVC = MDV - 3*NAtoms - 6*NNRF
        Call AMove(3*NAtoms,CMM,V(IC1))
        Call AClear(6*NNRF,V(IC2))
        ETemp = Zero
        If(LogTmp) then
          Call ISubsI(NAtoms,MicOpt,-1,NNRF+1)
          Call ISubsI(NAtoms,MicOpt,1,-1)
          endIf
        Call MMCGOp(IOut,IPrint,0,0,NAtoms,MicOpt,CMM,ETemp,FTemp,LOptx,
     $    ConMM1,CMM,.True.,.True.,FMMRMS,FMMMax,0,IStatx,.False.,0,XX,
     $    IMMRFO,V(IC1),V(IC2),.False.,QMVar,RotQM,.False.,.False.,
     $    LJunk,IRwSMM,V(IVC),MDVC)
        If(LogTmp) then
          Call ISubsI(NAtoms,MicOpt,-1,1)
          Call ISubsI(NAtoms,MicOpt,NNRF+1,-1)
          endif
        Finish = Finish.and.FMMMax.lt.CnvFMx*ConvF.and.FMMRMS.lt.ConvF
        endIf
C
C     If LARGE then remove the associated surface approximation for
C     LSTQST.
C
      If(LSTQST.and.Large.and.DoRed) Call DgQSTL(IOut,IPrint,IQST,NVar,
     $  NP,Cnstr1,Cnstr2,X,FF,Scr,NTBond,NTAng,NTDih)
C
C     Output, using variable names in XName.
C
      Write(IOut,1150)
      NCur = 0
      Do 310 I = 1, NVar
        Call GetB(2,Tmp,Len,XName,NCur)
        Call CBlank(8,IStr)
        ICur = 0
        Len = Min(Len,8)
        JPos = (10-Len) / 2 - 1
        Do 300 J = 1, Len
          IStr(JPos+J) = GetChr(Tmp,ICur)
  300     Continue
        Write(IOut,1170) (IStr(K),K=1,8),(APrint(J,I),J=1,6)
  310   Continue
      T1 = ConvF * CnvFMx
      T2 = ConvX * CnvFMx
      Write(IOut,1000)
      Call Convgd(FMax,T1,IReslt)
      Write(IOut,1010) FMax,T1,IReslt
      Call Convgd(FRMS,ConvF,IReslt)
      Write(IOut,1020) FRMS,ConvF,IReslt
      Call Convgd(DXMaxC,T2,IReslt)
      Write(IOut,1030) DXMaxC,T2,IReslt
      Call Convgd(DXRMSC,ConvX,IReslt)
      Write(IOut,1040) DXRMSC,ConvX,IReslt
      If(SimMMM) then
        C1 = FMax / CnvFMx
        C2 = FRMS
        C3 = DXMaxC / CnvFX
        C4 = DXRMSC / (CnvFMx*CnvFX)
        C2 = Min(C1,C2,C3,C4)
        CnvQM(1) = SclMMC*Max(C2*CnvFMx,T1)
        CnvQM(2) = SclMMC*Max(C2,ConvF)
        CnvQM(3) = SclMMC*Max(C2*CnvFX,T2)
        CnvQM(4) = SclMMC*Max(C2*CnvFX*CnvFMx,ConvX)
        Call RwFGED(2,IRwSM,NMax,MaxGDI,IEn,INNEn,ICoeff,IQNDIS,IEnDIS,
     $    IMask,ICnvQM)
        Call FileIO(1,-IRwSM,4,CnvQM,ICnvQM)
        endIf
      Length = ITqry(IMMRFO)
      If(Length.gt.0.and..not.((DdTSLn.and.DidLin).or.SmlStp)) then
        If(SimMMM) Call GauErr('SimMMM and MMRFO in GrdOpt.')
        Call FileIO(2,-IMMRFO,4,CnvQM,Length-5)
        CnvQM(1) = Min(FMax,CnvQM(1))
        CnvQM(2) = Min(FRMS,CnvQM(2))
        CnvQM(3) = Min(DXMaxC,CnvQM(3))
        CnvQM(4) = Min(DXRMSC,CnvQM(4))
        Call FileIO(1,-IMMRFO,4,CnvQM,Length-5)
        endIf
      If(OptMic) then
        Call Convgd(FMMMax,T1,IReslt)
        Write(IOut,1050) FMMMax,T1,IReslt
        Call Convgd(FMMRMS,ConvF,IReslt)
        Write(IOut,1060) FMMRMS,ConvF,IReslt
        endIf
      If(IPSave(0).le.1) Write(IOut,1080) DEPred
C
C     If we've converged then print final parameters.  Other wise check
C     for exceeding maximum number of steps, wrong curvature, etc.
C
      If(IDir.ge.2) Call AMove(NVar,FTmp,F)
      Finish = Finish.and..not.FixedC
      If(Finish) then
        If(RegFin) then
          Write(IOut,1180)
        else
          Write(IOut,1190)
          endIf
        If(.not.FCart.and..not.DoRed) Call RdWrZ(2,0,ANames)
        Do 320 I = 1, NVar
          X(I) = X(I) - XNew(I)
  320     IC(I) = 99
        If(DoRed) then
          Call RedTbl(IOut,1,XName,X,Cnstr1,Cnstr2,IC,F,NVar,IZRedO,
     $      ToAng,NTBond,NTAng,LSTQST,Zero)
        else if(.not.FCart) then
          Call PrmTbl(1,XName,X,IC,F,NVar,MaxNZ,IZ,LBL,LAlpha,LBeta,NZ,
     $      ToAng)
          endIf
      else if(.not.Steped) then
        Write(IOut,1380)
        Failed = .True.
      else if(NNegE.ne.Neg.and.SetStp.and.TstEig) then
        Write(IOut,1210) Neg, NNegE
        Failed = .True.
      else if(IStep.ge.NStep) then
        Write(IOut,1200) NStep
        Failed = .True.
        endIf
      Exit = Finish.or.Failed
      If(Failed) then
        Call ILSW(1,25,0)
        Call ILSW(1,27,1)
        If(.not.FCart.and..not.DoRed) Call RdWrZ(2,0,ANames)
        Do 330 I = 1, NVar
  330     IC(I) = 99
        If(DoRed) then
          Call RedTbl(IOut,2,XName,X,Cnstr1,Cnstr2,IC,F,NVar,IZRedO,
     $      ToAng,NTBond,NTAng,LSTQST,Zero)
        else if(.not.FCart) then
          Call PrmTbl(2,XName,X,IC,F,NVar,MaxNZ,IZ,LBL,LAlpha,LBeta,NZ,
     $      ToAng)
          endIf
      else if(.not.Exit.and.DXMaxC.lt.SmStMx) then
        Call ILSw(1,4,1)
        endIf
C     ISStep tracks whether we last did a step in an eigenvector
C     direction.
      If(DidEig) then
        ISStep = 1
      else if(.not.Rises) then
        ISStep = 0
        endIf
      If(BigLin.and..not.Quad.and.XLin.ge.(-One).and.XLin.le.Zero
     $  .and.SavLMO) then
        SavXLn = .True.
      else if(SavLMO) then
        Call CpySMO(IOut,IPrint,NP,ESave,SavXLn,V,MDV)
        endIf
      If(SavXLn) then
        RWhat(1) = XLin
      else
        RWhat(1) = Zero
        endIf
      Call XLSW(1,6,2,JJ(1),RWhat)
      Return
      End
*Deck GTrPth
      Subroutine GTrPth(IOut,IPrint,NAtoms,CR,CX,XValue,IZRed,NTRed,
     $  NTBond,NTAng,NTDih,WeiFac,MStruc,Rstr,PStr,XStr,B,IB,Square,Tri,
     $  One4D,CnStr,IFlags,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     For path relaxation . variant of GTrans
C
      Dimension B(3,4,NTRed), IB(4,NTRed),
     $  Square(NTRed,NTRed), Tri(*), One4d(NTRed,4), CR(3,NAtoms),
     $  CX(3,NAtoms), IZRed(4,NTRed), RStr(NTRed),
     $  PStr(NTRed), XStr(NTRed), XValue(NTRed), Cnstr(NTRed),
     $  IFlags(NTRed),XX(1),V(MDV)
      Save XX
      Data XX/1.0d0/
 1000 Format(' ***** Convergence failure in GTrPth2 *****',/,
     $  ' Try using 3 structures as input for',/,
     $  ' QST transition state search')
C
C     Translate and rotate CR to maximum coincidence with CP, but then
C     Gaussian version of Coince is unreliable.
c     call Coince(IOut,IPrint,Toang,NAtoms,NOrd,IAnP,CP,IanR,CR,RStr,
c    $  MPtoR,MRtoP)
C
C     Estmate CX by interpolation and then produce internal coordinates.
C
      NAt3 = 3*NAtoms
      Call AMove(Nat3,CR,CX)
C
C     Guess the transition structure using the weighting factor WeiFac,
C     XValue is the trans structure, then compare the two TS coordinate
C     values.
C
      Call ASub(NTRed,RStr,PStr,XValue)
      Call FixDq(IOut,IPrint,1,NTRed,NTBond,NTAng,NTDih,XValue,XX)
      Call ACASB(NTRed,PStr,XValue,XValue,WeiFac)
      Call FixDih(NTDih,XValue(NTBond+NTAng+1))
      Call ASub(NTRed,XValue,RStr,PStr)
      Call FixDq(IOut,IPrint,1,NTRed,NTBond,NTAng,NTDih,PStr,XX)
C
C     Now we have CX, deltaQ, which is PStr, we can find the true CX
C     which fits the redundant internal at XValue.  If RedCor fails,
C     use the old CX as the guess for the TS.
C
      Call RedCar(IOut,IPrint,0,NTRed,NTBond,NTang,NTdih,IZRed,RStr,
     $  NAtoms,CX,PStr,XStr,B,IB,Square,Tri,One4D,Cnstr,IFlags,.False.,
     $  IFail,0,V,MDV)
      If(IFail.eq.0) then
        Call AMove(NAt3,XStr,CX)
      else if(MStruc.eq.2) then
        Write(IOut,1000)
        Call Lnk1E(0)
      else
        Call AMove(NAt3,CR,CX)
        endIf
      Return
      End
*Deck GuessE
      Subroutine GuessE(N,Lower,NE,M,U,E,IE)
      Implicit Real*8(A-H,O-Z)
C
C     The routine creates guess eigenvectors using the diagonal elements
C     of matrix M.
C
C     O. Farkas, WSU, 1998.
C
      Real*8 M,U
      Dimension M((N*(N+1))/2),U(N,NE),E(N),IE(N)
      Logical Lower
C
      Save One
      Data One/1.D0/
C
      NNE = Min(N,NE)
      II  = 0
      Do 100 I = 1, N
        II = II + I
        E(I) = M(II)
  100   Continue
      Call MySort(N,E,IE)
      I0 = 0
      If (.not.Lower) I0 = N - NNE
      Call AClear(NNE*N,U)
      Do 200 I = 1, NNE
        U(IE(I0+I),I) = One
  200   Continue
      Return
      End
*Deck Ind103
      Subroutine Ind103(IStart,DoRed,AlUnit,JustCk,MicOTy,NVarZM,NVrMMT,
     $  MStepD,IMxStp,NAtoms,MaxStp,MaxSv,MaxSv2,MinLXX,IIC,IX,IXName,
     $  IESave,IFS,IFC1,IFC2,IXX,IFF,IEQMS,IXXMM,IFFMM,IFFMMQ,IITU,
     $  LOptG2,IForce,IFrcns,IEnd)
      Implicit Integer(A-Z)
C
C     Allocate the variable-dimension portion of /OptGrd/.  MStepD
C     is the minimum number of steps to allocate for, while MaxStp
C     is returned as the maximum number of steps actually allocated for.
C
      Logical DoRed, AlUnit, JustCk
C
      NGoal = 2*NVarZM
      If(DoRed) NGoal = Min(NGoal,6*NAtoms)
      If(JustCk) then
        ValMin = 2
      else
        ValMin = 0
        endIf
      MaxStp = MxStSt(MStepD,ValMin,NGoal,IMxStp)
      If(MicOTy.eq.4) MaxStp = Min(MaxStp,2)
      MaxSV = Max(NVarZM,MaxStp) + 1
      MaxSV2 = Min(MaxSV,MStepD) + 1
      NVTT = (NVarZM*(NVarZM+1))/2
      If(AlUnit) NVTT = 0
      Call All103(AlUnit,NVarZM,NVrMMT,MinLXX,MaxSv,MaxSv2,IStart,IIC,
     $  IX,IXName,IESave,IFS,IFC1,IFC2,IXX,IFF,IEQMS,IXXMM,IFFMM,IFFMMQ,
     $  IITU,LOptG2)
      IForce = IIC + LOptG2
      IFrcns = IForce + NVarZM
      IEnd = IFrcns + NVTT
      Return
      End
*Deck IndGDI
      Subroutine IndGDI(NP,NQ,Large,IStart,IQQS,IDQQS,IA,IAInv,ICC,
     $  ICSave,IHSc,IDQD,IDQR,IW,IEnd)
      Implicit Integer(A-Z)
C
C     Allocate the scratch space and variables for GDIIS.
C
      Logical Large
C
      IQQS   = IStart + 1
      IDQQS  = IQQS   + NQ*NP
      IA     = IDQQS  + NQ*NP*2
      IAInv  = IA     + (NP*(NP+1))/2
      ICC    = IAInv  + (NP*(NP+1))/2
      ICSave = ICC    + NP
      IHSc   = ICSave + NP
      IDQD   = IHSc   + NQ
      IDQR   = IDQD   + NQ
      IW     = IDQR   + NQ
      If(Large) then
        IEnd = IW + 6*Max(NP,NQ)
      else
        IEnd = IW + 6*NP
        endIf
      Return
      End
*Deck IniFSO
      Subroutine IniFSO(IRwPE,NMax,MaxGDI,NP,NDIIS,ESave,SimOpt,MaxSim,
     $  V,MDV)
      Implicit Real*8(A-H,O-Z)
      Parameter (IRwGen=501)
      Dimension V(*),ESave(*)
      Logical SimOpt
      Save One
      Data One/1.D0/
C
C     Save energies for LFO/SimOpt/GEDIIS
C
      Call RwFGED(2,IRwPE,NMax,MaxGDI,IEn,INNEn,ICoeff,IQNDIS,IEnDIS,
     $  IMask,ICnvQM)
      MaxSim = MaxGDI
      NDIIS = Min(NDIIS,MaxSim)
      NSave = Min(NP,MaxGDI)
      Call TstCor(NSave,MDV,'IRwSM')
      Call FileIO(2,-IRwPE,NSave-1,V,IEn)
      Call FileIO(1,-IRwPE,NSave-1,V,IEn+1)
      Call FileIO(1,-IRwPE,1,ESave(1),IEn)
      Call FileIO(2,-IRwGen,1,RepulN,40)
      Call FileIO(2,-IRwPE,NSave-1,V,INNEn)
      Call FileIO(1,-IRwPE,NSave-1,V,INNEn+1)
      Call FileIO(1,-IRwPE,1,RepulN,INNEn)
      Call FileIO(2,-IRwPE,NP-1,V,IMask)
      Call AMove1(NP-1,0,1,V)
      V(1) = One
      Call FileIO(1,-IRwPE,NP,V,IMask)
      Return
      End
*Deck IniGDI
      Subroutine IniGDI(IRwGDI,NMax,MaxGDI,NP,ESave,SimOpt,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Save energies for LFO/SimOpt/GEDIIS
C     The newest point is saved at the beginning of each array
C
C     NMax   : Max number of points that can be saved and used by l103 and l502
C     MaxGDI : Max number of points that can be used by DIIS thoughout optimization
C     NP     : Number of prvious points saved
C
      Parameter (IRwGen=501)
      Dimension V(*),ESave(*)
      Logical SimOpt
      Save One
      Data One/1.D0/
C
      Call RwFGED(2,IRwGDI,NMax,MaxGDI,IEn,INNEn,ICoeff,IQNDIS,IEnDIS,
     $  IMask,ICnvQM)
      NSave = Min(NP,NMax)
      Call TstCor(NSave,MDV,'IniGDI')
      Call FileIO(2,-IRwGDI,NSave-1,V,IEn)
      Call FileIO(1,-IRwGDI,NSave-1,V,IEn+1)
      Call FileIO(1,-IRwGDI,1,ESave(1),IEn)
      Call FileIO(2,-IRwGen,1,RepulN,40)
      Call FileIO(2,-IRwGDI,NSave-1,V,INNEn)
      Call FileIO(1,-IRwGDI,NSave-1,V,INNEn+1)
      Call FileIO(1,-IRwGDI,1,RepulN,INNEn)
      Call FileIO(2,-IRwGDI,NP-1,V,IMask)
      Call AMove1(NP-1,0,1,V)
      V(1) = One
      Call FileIO(1,-IRwGDI,NP,V,IMask)
      Return
      End
*Deck IniRSc
      Subroutine IniRSc(IOut,IOScan,IOGenA,LRwGen,NVar,IntVec,JScan)
      Implicit Real*8(A-H,O-Z)
C
C     Initialize link 103 for relaxed potential surface scans.
C
      Dimension IntVec(NVar), JScan(4)
      Parameter (MaxInf=DEFMAXINFO)
      Common /Info/ Info(MaxInf)
 1000 Format(' Number of optimizations in scan=',I4)
C
      JScan(1) = 1
      JScan(2) = 1
      Do 20 I = 1, NVar
        If(IAbs(IntVec(I)).gt.1000)
     $    JScan(2) = JScan(2) * (IAbs(IntVec(I))-1000)
   20   Continue
      If(JScan(2).gt.1) then
        Write(IOut,1000) JScan(2)
        Call ConDDF(IOScan,JScan(2)*(NVar+1)+InToWP(2))
        Call FileIO(1,-IOScan,intowp(2),jscan,0)
        Info(7) = 8
        Call ConDDF(IOGenA,JScan(2)*LRwGen)
        endIf
      Return
      End
*Deck InitBS
      Subroutine InitBS(IOp,ToAng,AlUnit,MaxStp,NVarD,MethQ,IGEVer,
     $  FCart,ConvP,IRwSM,X,Frcnst,FC,XName,IC,ITU,DoRed,IAprBG,LSTQST,
     $  NTrRot,NTBond,NTAng,NTDih,IZRed,Value,Cnstr1,Cnstr2,NAtoms,IAn,
     $  C,MxGDII,MaxGDI,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Initialize L103 for geometry optimizations.
C
#include "commonz.inc"
      Parameter (IRdFCD=1,ICalcE=2,ICalcF=3,ICalFF=-2,IFrzFF=-3)
      Logical FCart, DoRed, LSTQST, Large, AlUnit, DoStuf, ConvP
      Common /IO/     In,IOut,IPunch
      Common /Grdnt/  Energy,NVar,NDum
      Logical Prnt, Exit, HavDEP
      Common /OptGrd/ ConvF,FMaxT,DXMaxT,RMax,RMin,RLim,EigMax,EigMin,
     $  DEPred,FSwtch,FncErr,GrdErr,FncCnv,RLam,NRisDI,NStep,IStep,NMax,
     $  NP,Neg,Prnt,Exit,IPSav1,IPSav2,HavDEP,NGDIIS,IRisen,MaxRis,
     $  Large,ILevel,ISStep,IOffs,IMDS,MaxSv,MaxSv2,NVrMMT
      Dimension IOp(*), V(*), X(*), Frcnst(*), FC(NVarD,NVarD),
     $  XName(*), IC(*), IZRed(4,*), Cnstr1(*), Cnstr2(*), Value(*),
     $  IAn(*), C(3,*), ITU(*)
      Save CutOff,Delta,Zero,DefSCT
      Data CutOff/1.0d-6/,Delta/0.005D0/,Zero/0.0d0/,DefSCT/1.d-4/
 1000 Format(' NVar =',I4,' Out of valid range in L103.')
 1010 Format(' Option 10 out of range in L103, IOp(10)=',I6,'.')
 1020 Format(8F10.6)
 1030 Format(2I3,F20.0)
 1040 Format(10F8.6)
 1050 Format(' New optimization control parameters:',/,
     $       ' FMaxT,DXMaxT,RMax,RMin,RLim,EigMax,EigMin:',1X,7F10.6)
 1060 Format(' Expert switch is set: DXMaxT, EigMax, EigMin relaxed.')
 1070 Format(' Maximum allowed force set to',1PD9.2,'.')
 1080 Format(' Trust Radius=',1PD8.2,' FncErr=',1PD8.2,' GrdErr=',
     $  1PD8.2)
C
C     Check for valid NVar and read in /ZSubst/.
C
      IPrint = IOp(33)
      ISStep = 0
      RLam = GFloat(0)
      If(NVarZM.le.0) then
        Write(IOut,1000) NVarZM
        Call Lnk1E(0)
        endIf
      NVTT = (NVarZM*(NVarZM+1)) / 2
      If(IOp(5).ne.0.and.MethQ.eq.1)
     $  Call GauErr('Steepest descent can only locate minima.')
      IPSav1 = 0
      IPSav2 = 0
      HavDEP = .False.
      NGDIIS = 0
      IRisen = 0
      NRisDI = 0
      MaxRis = IOp(14)
      Call IClear(MaxSv,ITU)
C
C     Move relavent data read by link 101 to the optimization commons.
C
      IHOp = IOp(10)
      If(IHOp.eq.0) IHOp = 4
      NVar  = NVarZM
C
C     Read the second derivative matrix.
C
      DoStuf = .not.AlUnit
      If(IHOp.eq.1) then
        Read(In,1020) ((FC(I,J),J=1,I),I=1,NVar)
      else if(IHOp.eq.2) then
  120   Read(In,1030) I,J,FFIN
        If(I.ne.0) then
          FC(I,J) = FFIN
          FC(J,I) = FFIN
          Goto 120
          endIf
      else if(IHOp.eq.3.or.IHOp.eq.5.or.IHOp.eq.6) then
        NVar0 = NVar
        Call ReadFC(In,IOut,IHOp,IPrint,FCart,DoRed,IAprBG,NAtoms,IAn,
     $    C,NTrRot,NVar0,NTBond,NTAng,NTDih,IZRed,Value,XName,Frcnst,V,
     $    MDV)
        Do 138 I = 1, NVar
          IDX = (I*(I+1))/2
          If(IC(I).eq.0) then
            IC(I) = IRdFCD
          else if(IC(I).eq.1) then
            Frcnst(IDX) = FC(I,I)
          else if(IC(I).eq.2.or.IC(I).eq.3) then
            Frcnst(IDX) = FC(I,I)
            IF(FC(I,I).lt.CutOff) Frcnst(IDX) = Delta
            endIf
  138     Continue
        DoStuf = .False.
      else if(IHOp.eq.8) then
        Do 140 I = 1, NVarZM
          If(IC(I).eq.0) then
            FC(I,I) = GFloat(1)
            IC(I) = 1
            endIf
  140     Continue
      else if(IHOp.ne.4.and.IHOp.ne.7.and.IHOp.ne.10) then
        Write(IOut,1010) IHOp
        Call Lnk1E(0)
        endIf
C
C     Stuff read in force constants into FC.
C
      If(DoStuf) then
        IJ = 0
        Do 180 I = 1, NVar
          If(IC(I).eq.ICalcE.or.IC(I).eq.ICalcF) then
            If(Abs(FC(I,I)).le.CutOff) FC(I,I) = Delta
          else if(IAbs(IC(I)).le.1000) then
            If(Abs(FC(I,I)).gt.CutOff) IC(I) = IRdFCD
            endIf
          If(IOp(10).eq.4) then
            If(IC(I).ge.0.and.IC(I).lt.10) then
              IC(I) = ICALFF
            else if(IAbs(IC(I)).lt.1000) then
              IC(I) = IFrzFF
              endIf
            FC(I,I) = GFloat(0)
            endIf
          Do 160 J = 1, I
            IJ = IJ + 1
  160       Frcnst(IJ) = FC(I,J)
  180     Continue
        endIf
C
C     Print info concerning initial second derivatives.
C
      If(.not.AlUnit) Call AMove(NVTT,Frcnst,FC)
      If(DoRed) then
        Call RedTbl(IOut,0,XName,X,Cnstr1,Cnstr2,IC,FC,NVar,IZRed,
     $    ToAng,NTBond,NTAng,LSTQST,Zero)
      else if(.not.FCart) then
        Call PrmTbl(0,XName,X,IC,FC,NVar,MaxNZ,IZ,LBL,LAlpha,LBeta,
     $    NZ,ToAng)
        endIf
C
C     Now that printing has been done, flag the frozen variables so that
C     Estm won't clobber their force constants, if these were read in.
C
      If(IHOp.eq.3.or.IHOp.eq.5.or.IHOp.eq.6) then
        Do 210 I = 1, NVar
          If(IC(I).eq.-1.or.IC(I).eq.1) IC(I) = ICalFF
  210     Continue
        endIf
C
C     Set optimization control parameters.
C
      NP = -1
      If(IOp(10).eq.4) NP=0
      IStep = 0
      Neg = IOp(5)
      If(Neg.eq.0.and.MaxRis.eq.0) MaxRis = 1
      If(MaxRis.lt.0) MaxRis = 0
      NMax = NVar + 1
      MaxGDI = MxGDII
      If(IMDS.ge.2) then
C
C       DIIS vector size of 5 is seen to be optimum at this stage of development
C       Future development involves some dynamic control over the number of saved vectors
C
        NMax = Min(NMax,100)
        If(MaxGDI.eq.0) then
C         MaxGDI controls the maximum number of points that can be used in DIIS
          If(IGEVer.eq.1) then
            MaxGDI = Min(5,NMax)
          else if(IMDS.eq.3) then
            MaxGDI = Min(10,NMax)
          else
            MaxGDI = Min(5,NMax)
            MaxGDI = Min(10,NMax)
            endIf
          endIf
        Call RwFGED(0,IRwSM,NMax,MaxGDI,IEn,INNEn,ICoeff,IQNDIS,IEnDIS,
     $    IMask,ICnvQM)
        Call BNPack(NMax,MaxGDI,Temp)
        Call FileIO(1,-IRwSM,1,Temp,0)
        Call BNPack(0,0,Temp)
        Call FileIO(1,-IRwSM,1,Temp,3)
        If(ConvP) then
          SCFThr = Zero
        else
          SCFThr = DefSCT
          endIf
        Call FileIO(1,-IRwSM,1,SCFThr,1)
        Call TstCor(4,MDV,'InitBS-GED')
        Call AClear(4,V)
        Call FileIO(1,-IRwSM,4,V,ICnvQM)
      else if(DoRed) then
        NMax = Min(NMax,100)
        endIf
      Prnt = IPrint.ne.0.or.IOp(34).ne.0
C
C     Set various control parameters to standard Values.  For the
C     initial trust radius, a small value (corresponding to a maximum
C     step of 0.1 au or radians in one variable) is used.  This is
C     increased if force constants are being read in.
C
      RMax = GFloat(6) / GFloat(10)
      ConvF = GeoCon(IOp(7))
      RMin = GFloat(4) * ConvF
      RLim = GFloat(7) / GFloat(100)
      DEPred = GFloat(5) / GFloat(10000)
      NStep = IOp(6)
      If(NStep.eq.0) then
        NStep = Max(20,NVar+10)
      else
        NStep = IAbs(NStep)
        endIf
      NStep = Min(NStep,MaxStp)
      If(IOp(8).eq.0) then
        DXMaxT = GFloat(1) / GFloat(10)
        If(IOp(10).eq.1.or.IOp(10).eq.3.or.IOp(10).eq.4.or.IOp(10).eq.5
     $    .or.IOp(10).eq.6.or.IOp(10).eq.7.or.DoRed)
     $    DXMaxT = DXMaxT * GFloat(3)
      else
        DXMaxT = GFloat(IOp(8)) / GFloat(100)
        endIf
      FSwtch = GFloat(0)
      IOp15 = IOp(15)
      If(IOp(21).ne.0) IOp15 = -1
      If(IOp15.gt.0) then
        FMaxT = GFloat(IOp15) / GFloat(10)
      else
        FMaxT = GFloat(10)**20
        endIf
      If(IOp15.gt.0.or.Prnt) Write(IOut,1070) FMaxT
      If(IOp(16).eq.0) then
        EigMax = GFloat(250)
      else
        EigMax = GFloat(IOp(16))
        endIf
      If(IOp(17).eq.0) then
        EigMin = GFloat(1) / GFloat(10000)
      else
        EigMin = GFloat(1) / GFloat(IOp(17))
        endIf
      If(IOp(12).ne.0) then
        Read(In,1040) ConvF,FMaxT,DXMaxT,RMax,RMin,RLim,EigMax,
     $    EigMin,DEPred
        Write(IOut,1050) ConvF,FMaxT,DXMaxT,RMax,RMin,RLim,EigMax,
     $    EigMin,DEPred
        endIf
      If(IOp(21).ne.0) then
        Write(IOut,1060)
        DXMaxT = DXMaxT*GFloat(10)
        EigMax = EigMax*GFloat(10)
        EigMin = EigMin/GFloat(10)
        endIf
C
C     FncErr ... Estimate of error in value of function
C     GrdErr ... Estimate of error in value of gradient
C     FncCnv ... Unit conversion factor for function.
C
      IError = IOp(26)
      Call InFGEr(IOut,IPrint,.True.,.False.,IError,NAtoms,FncErr,
     $  GrdErr)
      If(IPSave(0).le.1.or.Prnt.or.IOp(8).ne.0.or.IError.lt.0)
     $  Write(IOut,1080) DXMaxT, FncErr, GrdErr
      FncCnv = GFloat(1)
      Return
      End
*Deck InPath
      Subroutine InPath(IStart,DoRed,NVarZM,Npath,MStepD,IMxStp,NAtoms,
     $  MaxStp,IITU,IIC,IX,INCurr,IXP,IESvP,IFSP,IFC1P,IFC2P,IXXP,IFFP,
     $  IICP,IXName,IESave,IFS,IFC1,IFC2,IXX,IFF,LOptG2,IForce,IFrcns,
     $  IFTStar,IVecStar,INStar,IIRisP,ICP,IConv,IRad,ITRot,IIRis,IINOp,
     $  ISym,IIPrA,ITan,IArcL,IPthL,IEnd)
      Implicit Integer(A-Z)
C
C     Allocate the variable-dimension portion of /OptGrd/.  MStepD
C     is the minimum number of steps to allocate for, while MaxStp
C     is returned as the maximum number of steps actually allocated for.
C
      Logical DoRed
      Parameter (MaxOp=DEFMAXOP)
C
      NGoal = 2*NVarZM
      If(DoRed) NGoal = Min(NGoal,6*NAtoms)
      MaxStp = MxStSt(MStepD,0,NGoal,IMxStp)
      NVTT = (NVarZM*(NVarZM+1))/2
      MaxSV = Max(NVarZM,MaxStp)+1
      IITU = IStart
      IIC = IITU + InToWP(MaxSv)
      IX = IIC + InToWP(MaxSV)
      IXName = IX + MaxSV
      IESave = IXName + MaxSV
      IFS = IESave + 1
      IFC1 = IFS + MaxSV
      IFC2 = IFC1 + NVTT
      IXX = IFC2 + NVTT
      IFF = IXX + NVarZM*MaxSv
      IESvP = IFF + NVarZM*MaxSv
C ESvP store the energy of each point
      IFSP = IESvP + NPath
C FSP store the history of the energy of each point
      IFC1P = IFSP + MaxSV*NPath
C FC1P stores the hessian for each point
      IFC2P = IFC1P + NVTT*NPath
C FC2P stores the previous hessian for each point
      IXXP = IFC2P + NVTT*NPath
C XXP stores the history of coordinates for each point
      IFFP = IXXP + NVarZM*MaxSv*NPath
C FFP stores the history of the gradient for each point
      IXP = IFFP + NVarZM*MaxSv*NPath
C XP stores the coordinates for each point
      INCurr = IXP + NVarZM*NPath
C NCurr stores the index of the current point
      IICP=INCurr+IntoWP(1)
C ICP stores the history of index for each point
      IFTStar= IICP + IntoWP(NPath*MaxSV)
C FTStar stores the redistribution info. for each point
      IVecStar=IFTStar + Npath
C VecStar stores the projector for each point
      INStar= IVecStar+NVTT*NPath
C NStar stores the TS index
      IIRisP= INStar+InToWP(1)
C IRisP stores the number of steps for each point
      ICP=IIRisP+InToWP(NPath)
C CP stores the current and previous cartesians for each point
      IConv=ICP+3*NAtoms*2*NPAth
C Conv stores indexes for the path
      IRad=IConv+IntoWp(Npath)
C Rad stores the trust radius for each point
      ITRot=IRad+NPath
C TRot stores the rotation matrix for each point
      IIRis=ITRot+9*NPAth
C IRis stores the energy rising info for reactant and products
C Symmetry Info
      ITan=IIRis+IntoWP(4)
C ITan stores Tangent
      IArcL=Itan+NVarZM*Npath
C IArcL stores the arc length
      IPthL=IArcL+NPath
C IPthL stores path length
      IINOp= IPthL+2
      ISym=IINOp+IntoWP(NPath)
      IIPrA=ISym+9*MaxOp*NPath
      LOptG2=IIPrA+IntoWP(NAtoms*MaxOp*NPath)-IStart+1
      IForce = IIC + LOptG2
      IFrcns = IForce + NVarZM
      IEnd = IFrcns + NVTT
      Return
      End
*Deck JustQM
      Subroutine JustQM(IOut,IMMCRS,NatNew,NAtoms,MicOpt,NVarZM,IBMat,
     $  BMat,A,F,NTrRot,C,NTBond,NTAng,NTDih,IZRed,Value,SqScr,Tri,
     $  One4D,IPrint,Scr1,Scr2,XName,IndMCR,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Take out the MM contribution from the forces in redundant internals.
C
      Dimension MicOpt(*), Scr1(*), Scr2(*), IBMat(*), BMat(*), A(*),
     $  JJ(1), F(*), C(*), IZRed(*), Value(*), SqScr(*), Tri(*), V(*),
     $  One4D(*), XName(*)
      Save JJ
      Data JJ/0/
C
      Call FileIO(2,-IMMCRS,3*NAtoms,Scr1,IndMCR)
      J = 0
      Do 10 I = 0, (NAtoms-1)
        If(MicOpt(I+1).lt.0) then
          Call AMove1(3,I*3,J,Scr1)
          J = J + 3
          endIf
   10   Continue
      If(NAtNew.eq.NAtoms) then
        NTRR1 = NTrRot
      else
        NTRR1 = -1
        endIf
      IGISav = 0
      Call Red2BG(IOut,IPrint,0,NAtNew,NTRR1,IGISav,0,C,NVarZM,NTBond,
     $  NTAng,NTDih,IZRed,Value,NRank,BMat,IBMat,A,SqScr,Tri,One4D,
     $  XName,V,MDV)
      Call TranF(.True.,.False.,NVarZM,NAtoms,JJ,Scr1,Scr2,IBMat,BMat,A,
     $  JJ)
      Call AAdd(NVarZM,F,Scr2,F)
      Return
      End
*Deck MemDIS
      Subroutine MemDIS(MaxGDI,NMax,NVarT,IStart,IAMRFO,IAMEn,IAMSim,
     $  IAMTm1,IAMTm2,ICoef1,ICoef2,IPotE,IHEig,IScr,IAQuad,IAFF,IACC,
     $  IANFF,IANCC,IQuad,ILinr,Mask,IMap,IPDIS,IEnd,MDV)
      Implicit Integer(A-Z)
C
C     Allocate the scratch space and variables for GDIIS.
C
      NDim   = MaxGDI + 1
      IAMRFO = IStart
      IAMEn  = IAMRFO + NDim*NDim
      IAMSim = IAMEn  + NDim*NDim
      IAMTm1 = IAMSim + NDim*NDim
      IAMTm2 = IAMTm1 + NDim*NDim
      ICoef1 = IAMTm2 + NDim*NDim
      ICoef2 = ICoef1 + NDim
      IPotE  = ICoef2 + NDim
      IHEig  = IPotE  + MaxGDI
      IScr   = IHEig  + NVarT
      IAQuad = IScr   + Max(NVarT*6,NMax,4*(MaxGDI+1))
      IAFF   = IAQuad + NVarT*MaxGDI
      IACC   = IAFF   + NVarT*MaxGDI
      IANFF  = IACC   + NVarT*MaxGDI
      IANCC  = IANFF  + NVarT*MaxGDI
      IQuad  = IANCC  + NVarT*MaxGDI
      ILinr  = IQuad  + NVarT
      Mask   = ILinr  + NVarT
      IMap   = Mask   + NMax
      IPDIS  = IMap   + MaxGDI
      IEnd   = IPDIS  + NDim*NDim
      Call TstCor(IEnd,MDV,'Optmz-GeoDIS')
      Return
      End
*Deck MixDIS
      Subroutine MixDIS(IOut,IPrint,IRwPE,NQ,NQRMS,NP,NUsed,NMax,MaxGDI,
     $  NNeg,DQQ,A,AInv,En,C,CM,DQRes,W,QQ,FF,Map,Switch,ISkip,FncErr,
     $  OK,DidBck)
      Implicit Real*8 (A-H,O-Z)
C
C     Calculate the Simultaneous DIIS coefficients
C
C     Parameters:
C     IOut   : Output channel
C     IPrint : Debug print level
C     NQ     : Dimension of the problem
C     NQRMS  : Dimension to use in RMS calculation.
C     NP     : Number of saved variables and error vectors
C     NUsed  : Number of used vectors in GDIIS on output
C     NBasis : Number of basis functions
C     DQQ    : Geometric Error vectors on input
C     A      : Square matrix of the error vectors (NP+1)*(NP+1)
C     AInv   : Inverse of A matrix
C     En     : Scratch space (NP+1)
C     C      : DIIS coefficients
C     CM     : Array used to indicate deleted points
C              Final DIIS coefficients with ill-conditioned points masked
C              upon output
C     DQRes  : Sum DQQi * Ci (Geometric Residuum on output)
C     W      : Working space of length 5*NP
C     NUsed  : Number of used vectors in GDIIS
C
C     Stored in IRwPE
C     NMax   : Max number of points to be saved
C     NP     : Number of points saved
C     NConvP : SCF convergence threshold (not used here).
C     Energy : Energies saved (NMax)
C     DIISC  : DIIS coefficients (NMax)
C     A1     : SCF-DIIS matrix using quasi-Newton step (NMax*NMax)
C     A2     : SCF-EDIIS matrix (NMax*NMax)
C
      Dimension DQQ(NQ,*),C(NP+1),CM(NP+1),DQRes(NQ),W(5*(NP+1)),
     $  A(NP+1,NP+1),AInv(NP+1,NP+1),En(*),QQ(NQ,NP),FF(NQ,NP),Map(*),
     $  XX(1)
      Logical Print, OK, Inv2, EDIIS, DidBck, Rises, DidRem, Prnt1,
     $  TstBck
      Character*4 NmDIIS
      Save Zero, One, ValMin, CThrhd, XX, SmlDif
      Data Zero/0.0d0/, One/1.0d0/, ValMin/-1.0d0/, CThrhd/10.0d0/
     $  XX/0.0d0/, SmlDif/1.d-5/
      Data CLower, CUpper, CLarge/-0.5D0, 1.5D0, 5.D0/
 1000 Format(' SimOpt point',I4,':',/,' Geom Err:',D18.10)
 1010 Format(' NNeg=',I2,' NP=',I2,' Switch=',1PD10.2,' Rises=',L1,
     $  ' DC=',1PD10.2,' SmlDif=',1PD10.2)
 1020 Format(' DIIS inversion failure, using last point.')
 1030 Format(' DidBck=',L1,' Rises=',L1,1X,A,'DIIS coefs:',20F11.5)
 1050 Format(' Energy Point ',I5,' ',D18.10)
 1060 Format(' RMS Error=',D18.10,' NUsed=',I2,' EDIIS=',L1)
 1070 Format(' Guided ',A,'DIIS for TS search',5X,'SclDIS = ',F8.1)
C
C     Start of the first order simultaneous optimization
C
      Print = IPrint.ge.1
      Prnt1 = Print.or.IPSave(0).eq.0
      OK    = .False.
      NASq  = (NP+1)*(NP+1)
      NUsed = NP
      ISkip = 0
      Call RwFGED(1,IRwPE,NMax,MaxGDI,IEn,INNEn,ICoeff,IQNDIS,IEnDIS,
     $  IMask,ICnvQM)
      Call ASet(NUsed,One,CM)
C
C     Evaluate current point and set optimization options
C
      TError = SProd(NQRMS,DQQ(1,1),DQQ(1,1))
C
C     Turn off EDIIS when the error is small enough
C
      RMSErr = Sqrt(TError/GFloat(NQRMS))
      Call FileIO(2,-IRwPE,NUsed,C,IEn)
      Rises = .False.
      FncE1 = FncErr*GFloat(1000)
      E1F = C(1) - FncE1
      Do 20 I = 2, NUsed
   20   Rises = Rises.or.E1F.gt.(C(I)+FncE1)
      If(NNeg.gt.0) then
        EDIIS = NP.ge.4.and.(RMSErr.gt.Switch)
        If(RMSErr.lt.Switch) then
          SclDIS = One
        else
          SclDIS = GFloat(NQ)
          endIf
      else
        EDIIS = NP.ge.4.and.(RMSErr.gt.Switch.or.Rises)
        endIf
      If(NUsed.gt.1) then
        DC = C(NUsed) - C(NUsed-1)
        EDIIS = EDIIS.and.Abs(DC).ge.SmlDif
      else
        DC = Zero
        endIf
      If(Prnt1) Write(IOut,1010) NNeg, NP, Switch, Rises, DC,
     $  SmlDif
      If(Prnt1) Write(IOut,1060) RMSErr, NUsed, EDIIS
      NmDIIS = ' '
C
C     If the error for the current point is small enough, take this step
C
      If(TError.le.ValMin) then
        NUsed = NP
        OK = .True.
        Call AClear(NP,CM)
        CM(1) = One
      else if(EDIIS) then
C       In EDIIS, the newest point is saved at the end.
        Call AClear(NASq,AInv)
        Call AClear(NUsed,En)
C
C       Prepare EDIIS elements e = -g*R
C
        Do 30 I = 1, NUsed
          En(I) = C(NUsed-I+1)
          If(Print) Write(IOut,1050) NUsed-I+1,En(I)
          Do 30 J = 1, NUsed
C           Weight error from first reaction coordinate by SclDIS
            If(NNeg.gt.0) then
              AInv(I+1,J+1) = -SclDIS*FF(1,NUsed-I+1)*QQ(1,NUsed-J+1)
     $          -SProd(NQ-1,FF(2,NUsed-I+1),QQ(2,NUsed-J+1))
            else
              AInv(I+1,J+1) = -SProd(NQ,FF(1,NUsed-I+1),QQ(1,NUsed-J+1))
              endIf
   30       Continue
        Do 40 I = 1, NUsed
   40     AInv(I+1,1) = CM(NUsed-I+1)
C
C       Prepare EDIIS A matrix a(i,j)=e(i,i)+e(j,j)-e(i,j)-e(j,i)
C
        Call AClear(NASq,A)
        Do 50 I = 0, NUsed
          Do 50 J = I, NUsed
            If(I.eq.0 .and. J.eq.0) then
              A(1,1) = Zero
            else if(I.eq.0.or.J.eq.0) then
              A(I+1,J+1) = -One
              A(J+1,I+1) = A(I+1,J+1)
            else
C             Add in quadratic term and scale
              If(NNeg.gt.0) then
                A(I+1,J+1) = 0.5D0*(SclDIS*DQQ(1,I)*DQQ(1,I) +
     $            SclDIS*DQQ(1,J)*DQQ(1,J) + SProd(NQ-1,DQQ(2,I),
     $            DQQ(2,I)) + SProd(NQ-1,DQQ(2,J),DQQ(2,J))) +
     $            AInv(I+1,I+1) + AInv(J+1,J+1) - AInv(I+1,J+1) -
     $            AInv(J+1,I+1)
                A(J+1,I+1) = A(I+1,J+1)
              else
                A(I+1,J+1) = AInv(I+1,I+1) + AInv(J+1,J+1)
     $            - AInv(I+1,J+1) - AInv(J+1,I+1)
                A(J+1,I+1) = A(I+1,J+1)
                endIf
              endIf
   50       Continue
C
C       Solve EDIIS equation with constraint c(i)>0
C
        Call EnCoef(IOut,IPrint-1,NP,NUsed,En,A,AInv,W,ISkip,IWarn)
        ISkip = 0
        Do 60 I = 1, NUsed
   60     CM(I) = AInv(NUsed-I+2,1)
        If(Abs(CM(I)-One).ge.SmlDif) then
          OK = .True.
          NmDIIS = ' En-'
          endIf
        endIf
      If(.not.OK) then
C       RFO-DIIS: Beginning of the main loop to construct DIIS A matrix
C       (NP+1)*(NP+1)
        Call AClear(NASq,A)
        Do 70 I = 0, NUsed
          Do 70 J = I, NUsed
            If(I.eq.0.and.J.eq.0) then
              A(1,1) = Zero
            else if(I.eq.0.or.J.eq.0) then
              A(I+1,J+1) = -One
              A(J+1,I+1) = A(I+1,J+1)
            else
C             Weight error along TS vector for RFO-DIIS
              If(NNeg.gt.0) then
                A(I+1,J+1) = SclDIS*DQQ(1,I)*DQQ(1,J)+SProd(NQ-1,
     $            DQQ(2,I),DQQ(2,J))
                A(J+1,I+1) = A(I+1,J+1)
              else
                A(I+1,J+1) = SProd(NQ,DQQ(1,I),DQQ(1,J))
                A(J+1,I+1) = A(I+1,J+1)
                If(I.eq.J) En(I) = A(I+1,J+1)
                If(Print.and.I.eq.J) Write(IOut,1000) I,En(I)
                endIf
              endIf
   70       Continue
C
C       The A matrix is ready, then calculate the inverse.
C
        If(NUsed.ne.NP) Call Pcck(0,KK,A,NP+1,NP+1,NUsed+1,NUsed+1)
        NDim = NUsed + 1
   80   If(IPrint.ge.3) Call OutMat(IOut,1,A,NDim,NDim,NDim,NDim)
        Call AMove(NDim*NDim,A,AInv)
        OK = Inv2(IOut,0,AInv,NDim,W(1),W(2*NDim+1),W(3*NDim+1),
     $    W(4*NDim+1),Det,XX,XX,XX,0)
C
C       The inverse matrix is ready. It is time to calculate the DIIS
C       coefficients
C
        Call AClear(NDim,C)
        If(.not.OK) then
          Write(IOut,1020)
          NUsed = 1
          Call AClear(NDim,CM)
          CM(1) = One
          Call AMove(NQ,DQQ,DQRes)
        else
          Call ANeg(NDim,AInv,C)
C
C         DIIS coefficients check, if one of the older points has a very
C         large coefficient, delete related rows and columns and solve
C         the DIIS equation again.
C
          Call ChkDSM(IOut,Prnt1,A,C(2),CM,Map,NDim,NUsed,CThrhd,OK,
     $      DidRem)
          If(DidRem) Goto 80
          If(.not.OK) Return
          NmDIIS = 'RFO-'
          endIf
        endIf
      If(NP.lt.NMax) ISkip = 0
C
C     Calculated residuum, deleted points have been masked by
C     zero coefficients
C
      If(NNeg.gt.0.and.OK) Call ChkDIS(IOut,Prnt1,Prnt1,NUsed,CM,CUpper,
     $  CLower,CLarge,OK)
      Call MatMpy(NQ,NUsed,1,DQQ,CM,DQRes)
      Call FileIO(1,-IRwPE,NUsed,CM,ICoeff)
      DidBck = TstBck(NUsed,CM)
      If(NNeg.gt.0) then
        Write(IOut,1070) NmDIIS,SclDIS
      else
        Write(IOut,1030) DidBck, Rises, NmDIIS, (CM(I),I=1,NUsed)
        endIf
      Return
      End
*Deck MixOpt
      Subroutine MixOpt(IOut,IPrint,IMDS,MethQ,Large,NP,NNeg,SimOpt,
     $  Mixed,IStep,JScan,SwitMx,DoLinS,DoDIIS,FRMS,NDIIS)
      Implicit Real*8 (A-H,O-Z)
      Logical Large,SimOpt,Mixed,DoLins,DoDIIS
      Dimension JScan(2)
 1000 Format(' Mixed Optimization -- RFO/linear search')
 1010 Format(' Mixed Optimization -- En-DIIS/RFO-DIIS')
 1030 Format(' RMS Force = ',D10.5,' SwitMx=',D10.5,' MixMth=',I2)
 1040 Format(' First Order Simultaneous Optimization')
C
      SimOpt = IMDS.eq.2.or.IMDS.eq.5
      Mixed  = IMDS.eq.3.or.IMDS.eq.5
      DoDIIS = .False.
      If(Mixed) then
        MixMth = 1
        If((IStep.gt.1.or.JScan(2).eq.0.or.JScan(1).eq.1)
     $    .and.FRMS.le.SwitMx) MixMth = 2
        If(NNeg.gt.0) MixMth = 2
        If(IPrint.ge.1.or.IPSave(0).eq.0)
     $    Write(IOut,1030) FRMS, SwitMx, MixMth
        If(MixMth.eq.1) then
          MethQ = 3
          NDIIS = -1
          Write(IOut,1000)
        else if(MixMth.eq.2) then
          Write(IOut,1010)
          If(NDIIS.lt.1) then
            NDIIS = 1
            MethQ = 3
            endIf
        else
          Call GauErr('Unrecognized MixMth in MixOpt.')
          endIf
        DoDIIS = MethQ.ne.3
        DoLinS = .not.DoDIIS
        endIf
      Return
      End
*Deck MxStSt
      Function MxStSt(IniVal,MaxVal,MinVal,IMxStp)
      Implicit Integer(A-Z)
C
C     Return MaxStp, doing a min/max with MaxVal/MinVal if they are
C     non-zero and accounting for the IMxStp option.  The min is applied
C     2nd if both are requested.
C
      MaxStp = IniVal
      If(MinVal.gt.0) MaxStp = Max(MaxStp,MinVal)
      If(MaxVal.gt.0) MaxStp = Min(MaxStp,MaxVal)
      If(IMxStp.gt.0) then
        MaxStp = Max(MaxStp,IMxStp)
      else if(IMxStp.lt.0) then
        MaxStp = Min(MaxStp,-IMxStp)
        endIf
      MxStSt = MaxStp
      Return
      End
*Deck NRStep
      Subroutine NRStep(IOut,IPrint,Init,UpdLam,ConvEq,RLam,ELow,N,H,
     $  HInv,F,dX,W,NG,NA,IB,B)
      Implicit Real*8 (A-H,O-Z)
C
C     Calculate the Newton-Raphson or RFO step by solving
C     F = (M+Lambda*I) * dX
C
C     Odon Farkas, WSU, 1997-1998
C
      Parameter (NB=9)
      Dimension H((N*(N+1))/2), HInv(N*(N+1)/2), F(N), dX(N), W(20*N),
     $  IB(*), B(*), UMU((NB*(NB+1))/2), UU(NB*NB), WB(6*NB), E(NB)
      Logical Init, UpdLam
      Real*8 MDCutO
      Save F1, F0
      Data F1/1.0d0/, F0/0.0d0/
 1020 Format(' RFO step:  Lambda=',1PD15.8,'.')
 1030 Format(' lowest eigenvalue of the Hessian is',G12.4,'.')
C
      Call TStamp(1,'Enter NRStep')
      NA3 = 3*NA
      NA3T = (NA3*(NA3+1))/2
C
C     This is a size independent RFO implementation and
C     the number of strongly  coupled coordinates can be adjusted as
C     a constant multiplier.
C     Later this can be a user defined option.
C
      RNSiz = Max(Sqrt(GFloat(NA3)),F1)
      SScal = RNSiz
      Eps = MDCutO(0)
      IFX = NA3T + 1
      IDX = IFX + NA3
      IU  = 1
      IMU = IU  + NB * NA3
      IUW = IMU
C
C     Transform internal forces into cartesian forces
C     in order to sovle the linear equation in cartesian space
C     Note: cartesian space is smaller than internal space (redundency)
C
      If(IB(1).eq.0.or.IB(1).gt.NA) then
        Call MSBtQ(N,NA,IB,B,F,H(IFX))
      else
        Call MulBtQ(N,NA,IB,B,F,H(IFX))
        endIf
C
C     Generate the initial guess for Hessian
C     Initialize it as the inverse of the diagonoal of Hessian
C
      If(Init) then
        Call AClear(NA3T,HInv)
        Do 10 I = 1, NA3
          II = (I*(I+1))/2
          If(Abs(H(II)).gt.Eps) HInv(II) = F1/Abs(H(II))
   10     Continue
        endIf
C
C     Initialize Lamda
C     Note: Logical comparison between real numbers (RLam.eq.F0) should be avoided
C     suggest (RLam-F0).lt.ThrHd
C
      If(UpdLam.and.RLam.eq.F0) RLam = SProd(NA3,H(IFX),H(IFX))/RNSiz
      RLam0 = F0
      If(UpdLam) then
        NE = 1
        NG = 0
        ELow = H(1)
        Call FewEig(NA3,H,NG,NE,E,NB,W(IU),W(IMU),UMU,UU,WB,7)
        ELow = E(1)
        If(ELow.lt.F0) RLam0 = -ELow
        If(IPSave(0).eq.0) Write(IOut,1030) ELow
        endIf
      RLamE = RLam0 + RLam
C
C     Solve linear equation f=(H+RLamE*I)*dX
C
      Call SLEqF1(IOut,IPrint,NA3,SScal,RLam0,RLamE,H,HInv,H(IDX),
     $  H(IFX),W(IUW),ConvEq,0,UpdLam,IFail)
      If(IFail.ne.0) Call SLEqF1(IOut,IPrint,NA3,SScal,RLam0,RLamE,H,
     $  HInv,H(IDX),H(IFX),W(IUW),ConvEq,0,UpdLam,IFail)
C
C     Transform dX from cartesian space into internal space
C
      If(IB(1).eq.0.or.IB(1).gt.NA) then
        Call MSBX(N,IB,B,H(IDX),dX)
      else
        Call MulBX(N,NA,IB,B,H(IDX),dX)
        endIf
      RLam = RLamE
      If(IPrint.ge.0.and.UpdLam) Write(IOut,1020) RLam
      Init = .False.
      Call TStamp(1,'NRStep')
      Return
      End
*Deck OptDat
      Subroutine OptDat(Q0,XXP,GradPt,FFP,FCP,FC2P,ESvP,FSP,VecInP,ICP,
     $  TRotP,CP,Value,XX,F,FF,FC,FC2,ESave,FS,VecInt,IC,TRot,C,NAtoms,
     $  NVar,MaxSv1,NVTT,NPath,NN)
      Implicit Real*8(A-H,O-Z)
C
C     Read or store data for path relaxation
C
      Dimension Q0(Nvar,NPath),XXP(NVar,MaxSv1,NPath),GradPt(NVar,NPath)
      Dimension FFP(Nvar,MaxSv1,NPath),FCP(NVTT,NPath),FC2P(NVTT,NPath)
      Dimension ESvP(NPath),FSP(MaxSv1,NPath),VecInP(NVTT,NPath)
      Dimension ICP(MaxSv1,NPath),TRotP(3,3,NPath),CP(3,NAtoms,2,NPath)
      Dimension Value(*),XX(Nvar,*),F(1),FF(*),FC(1),FC2(1),Esave(*)
      Dimension FS(1),VecInt(*),IC(1),TRot(3,3),C(3,*)
C
      IType=NN
      If(IType.eq.0) Call GauErr('Wrong IType in OptDat')
      If(IType.lt.0) then
      NN=-NN
      Call AMove(NVar,Q0(1,NN),Value)
      Call IMove(MaxSv1,ICP(1,NN),IC)
      Call AMove(NVar*MaxSv1,FFP(1,1,NN),FF)
      Call AMove(NVar*MaxSv1,XXP(1,1,NN),XX)
      Call AMove(NVtt,FCP(1,NN),FC)
      Call AMove(NVtt,FC2P(1,NN),FC2)
      Call AMove(MaxSv1,FSP(1,NN),FS)
      Call AMove(NVar,GradPt(1,NN),F)
      Call AMove(NVTT,VecInP(1,NN),VecInt)
      Call AMove(9,TRotP(1,1,NN),TRot)
      Call AMove(3*NAtoms,CP(1,1,1,NN),C)
      ESave(1)=ESVP(NN)
      NN=-NN
      Endif
      If(IType.gt.0)then
      Call AMove(NVar,Value,Q0(1,NN))
      Call IMove(MaxSv1,IC,ICP(1,NN))
      Call AMove(NVar*MaxSv1,FF,FFP(1,1,NN))
      Call AMove(NVar*MaxSv1,XX,XXP(1,1,NN))
      Call AMove(NVtt,FC,FCP(1,NN))
      Call AMove(NVtt,FC2,FC2P(1,NN))
      Call AMove(MaxSv1,FS,FSP(1,NN))
      Call AMove(NVar,F,GradPt(1,NN))
      Call AMove(NVTT,VecInt,VecInP(1,NN))
      Call AMove(9,TRot,TRotP(1,1,NN))
      Call AMove(3*NAtoms,C,CP(1,1,1,NN))
      ESvP(NN)=ESave(1)
      Endif
      Return
      End
*Deck OptDXM
      Function OptDXM(DoRed,NVar,NTRedS,DeltaX)
      Implicit Real*8(A-H,O-Z)
C
C     Return a maximum step element for trust-radius testing.  For
C     internal coordinates, this may weight dihedrals differently.
C
      Logical DoRed
      Dimension NTRedS(*), DeltaX(*)
      Save Pt1
      Data Pt1/0.1d0/
C
      If(DoRed) then
        NVSB = NVar - NTRedS(4)
        DXSB = ArrMax(DeltaX,NVSB,.True.)
        DXD = ArrMax(DeltaX(1+NVSB),NTRedS(4),.True.)
        DX = Max(DXSB,Pt1*DXD)
      else
        DX = ArrMax(DeltaX,NVar,.True.)
        endIf
      OptDXM = DX
      Return
      End
*Deck OptFlag
      Subroutine OptFlag(Neg,MethQ,DoLin,DoLin2,Rises,NIRis,Steped,IQST,
     $  IQST1,FixEig,SclEig,IType)
      Implicit Real*8(A-H,O-Z)
C  For path relaxation.
C  Set flags for optimizer. IType=1 min; Itype>1 saddle point
C  This a rudimentary routine that can be modified to make the
C  optimizer's options more complexed.
      Logical DoLin,DoLin2,Rises,Steped,FixEig,SclEig
C
      Steped=.false.
      FixEig=.false.
      SclEig=.false.
      Neg=IType-1
      MethQ=3
      Rises=.false.
      If(NIRis.eq.1.and.Neg.lt.1)Rises=.true.
      DoLin2=.false.
      If(Neg.lt.1.and.DoLin)DoLin2=.true.
      IQST1=0
      If(Neg.ge.1)IQST1=IQST
      Return
      End
*Deck OptStp
      Subroutine OptStp(IOut,IPrint,Sphere,DoLin,MethQ,FCart,UpdRWF,
     $  IGetFC,IDimSc,XName,FC,FC2,ESave,FS,XX,FF,XQuad,XNew,FTemp,A,
     $  VecInt,Scr,FCTemp,IFlags,NAtoms,IAn,AtmChg,C,UseSym,TRot,Wght,
     $  IQST,DoRed,NTBond,NTAng,NTDih,NVar,Cnstr1,Cnstr2,NVarM,Rises,
     $  NTrRot,IEStpM,NFroz,Frozen)
      Implicit Real*8(A-H,O-Z)
C
C     For path relaxation.
C     This routine gathers all the stepping routines found in grdopt.
C
      Logical OK, Linr, Steped, Prnt, Exit, Rises, DidLin, Quad, FixEig,
     $  LinMin, LinTS, Despar, Sphere, DoLin, DidSph, UpdRWF, Frozen(*),
     $  Finish, SclEig, FCart, HavDep, UseSym, DoRed, Large, DidEig
      Dimension XName(1), FC(1), FS(1), XX(NVar,*), FF(*), XQuad(1),
     $  XNew(NVar), TRot(3,3), FTemp(1), A(*), ESave(*), Scr(IDimSc,8),
     $  FCTemp(1), FC2(1), IFlags(1), C(3,*), AtmChg(*), VecInt(*),
     $  Cnstr1(NVar), Cnstr2(NVar), XXA(1), IAn(*)
      Common /OptGrd/ ConvF,FMaxT,DXMaxT,RMax,RMin,RLim,EigMax,EigMin,
     $  DEPred,FSwtch,FncErr,GrdErr,FncCnv,XXIRC,NLEVEL,NStep,IStep,
     $  NMax,NP,Neg,Prnt,Exit,IPSav1,IPSav2,HavDEP,NGDIIS,IRisen,MaxRis,
     $  Large,ILevel,ISStep,IOffs,IMDS,MaxSv,MaxSv2,NVrMMT
      Integer Tmp(20)
      Save Zero,CnvFMx,CNVFX,SteepM,XXA
      Data Zero,CnvFMx,CNVFX/0.D0,1.5D0,4.D0/,SteepM/0.05d0/,XXA/0.0d0/
 1130 Format(' Maximum step size (',F8.3,') exceeded in linear search.',
     $     /,'    -- Step size scaled by ',F7.3,/,
     $       '    -- Skip Quadratic or steepest descent search.')
 1220 Format(' Linear search not attempted -- option 19 set.')
 1222 Format(' Linear search not attempted -- no second derivatives',
     $       ' available.')
 1224 Format(' Linear search not attempted -- first point.')
 1226 Format(' Linear search not attempted -- rms force is less than',
     $       ' FSwtch (',F8.5,').')
 1227 Format(' Linear search not attempted -- Energy or forces',
     $       ' went down.')
 1230 Format(' Energy Rises -- skip Quadratic search.')
 1235 Format(' Energy Rises -- skip steepest descent search.')
 1237 Format(' Gradient Rises -- skip Quadratic search.')
 1240 Format(' Maximum step size (',f8.3,
     $       ') exceeded in Quadratic search.',/,
     $       '    -- Step size scaled by ',F7.3)
 2010 Format(' Optimization aborted.',/,
     $       '    -- Gradient out of range.',/,
     $       '    -- Maximum allowed force = ',F8.3)
 2041 Format(' Quartic linear search produced a step of ',F8.5,'.')
 2042 Format(' Quintic linear search produced a step of ',F8.5,'.')
 2043 Format(' Grad norm linear search produced a step of ',F8.5,'.')
 2050 Format(' Skip linear search -- no minimum in search direction.')
 2070 Format(' Steepest descent instead of Quadratic search.')
 2080 Format(' Steepest descent step Scaled up to rms of',F8.5,'.')
 2090 Format(' Unrecognized ICase=',I10,'.')
C
C     Test gradient.
C
      Steped=.false.
      RtNV = GFloat(NVarM)
      RtNV = Sqrt(RtNV)
      FMax = ARRMax(FTemp,NVar,.True.)
      FRMS = SPROD(NVar,FTemp,FTemp)
      FRMS = Sqrt(FRMS/NVarM)
      Finish = (FRMS.lt.ConvF) .and. (FMax.lt.(CnvFMx*ConvF))
c     FSwtch=ConvF
      If(FMax.ge.FMaxT.and.UpdRWF) then
        NCur = 0
        Call PutCSt('-DE/DX',Tmp,NCur)
        Call PutDel(2,Tmp,NCur)
        Call MatPrt(FTemp,1,NVar,1,NVar,1,1,Tmp,XName,0,XXA,0)
        Write(IOut,2010) FMaxT
        Call Lnk1E(0)
        endIf
      Call AClear(NVar,XNew)
      Call AClear(NVar,XQuad)
C
C     Linear step.  This is only done if we are allowed to do it,
C     at least two points have been computed, and either the
C     search is for a local minimum or there are analytic second
C     derivatives at both points and the gradient has increased.
C
      If(DoRed) Call FixHis(IOut,IPrint,NVar,NP,XX,Scr,NTBond,NTAng,
     $  NTDih)
      LinMin = Neg.eq.0.and.(DoLin.or.(MethQ.eq.1.and.Rises))
      LinTS = Neg.ne.0.and.MethQ.eq.4.and.Rises.and.
     $                         IPSav1.ne.0.and.IPSav2.ne.0
      Linr = NP.gt.1.and.FRMS.gt.FSwtch.and.(LinMin.or.LinTS)
      Quad = MethQ.ne.1
      If(Linr) then
        Call DXLinr(ICase,Neg,FTemp,XNew,DXRMS,DXMax,NVar,DidLin,
     $    XLin,ESave(1),FS,FC,FC2,XX,FF,Scr(1,1),Scr(1,2),Scr(1,3),XXA,
     $    XXA,1,1,1,0,0,.False.,XXA,0)
        Steped = Steped .or. DidLin
        DX = DXRMS * RtNV
        If(DX.gt.DXMaxT.and.IRisen.eq.0) then
          Scale=DXMaxT/DX
          XLin=XLin*Scale
          Write(IOut,1130) DXMaxT,Scale
          Call AScale(NVar,Scale,XNew,XNew)
          DXMax = DXMax * Scale
          DXRMS = DXRMS * Scale
          Quad = .False.
        else
          DXMax = Zero
          If(Rises.and.Quad.and.LinMin) Write(IOut,1230)
          If(Rises.and..not.Quad.and.Neg.eq.0) Write(IOut,1235)
          If(Rises.and.LinTS) Write(IOut,1237)
          endIf
        endIf
C
C     Quadratic contributions to displacements.  These are done unless
C     things have gone wrong AND a linear search has been taken which
C     should fix the situation; that is unless the Energy has risen
C     when searching for a minimum or the forces have gone up when
C     searching for a saddle point and the conditions for a linear
C     search have been met (vide supra).  Use Scr(,2) for forces in
C     the eigenvector basis and Scr(,1) for the eigenvalues.
C
      Quad = Quad.and..not.(Rises.and.Linr.and.(DidLin.or.Neg.eq.0))
      DidSph = .False.
      If(Quad) then
        FixEig = MethQ.ne.2.and.MethQ.ne.3
        SclEig = IGetFC.ne.4.and.FixEig
        Call DiagFC(IOut,IPrint,FixEig,SclEig,FCart,.False.,NTrRot,NVar,
     $    EigMin,EigMax,IFlags,FC,FCTemp,XName,A,IDimSc,Scr,.False.,0,
     $    VecInt,FTemp,Neg,NNegE,NEig,NAtoms,C,AtmChg,NFroz,Frozen,Wght,
     $    IQST,IStep,DoRed,XX,Cnstr1,Cnstr2,NTBond,NTAng,NTDih,XXA,0)
        If(MethQ.eq.3.or.(MethQ.eq.2.and.NNegE.ne.Neg)) then
          Call DXRFO(IOut,IPrint,.False.,0,NVar,Neg,NEig,DXMaxT,IEStpM,
     $      Scr(1,1),A,Scr(1,2),XQuad,OK,DidEig)
        else
          Call DXQuad(IPrint,Sphere,.False.,FTemp,XQuad,A,OK,NVar,NEig,
     $      Scr(1,2),Scr(1,1),Scr(1,3))
          DidSph = Sphere
          endIf
        If(FCart) Call ClnCor(IOut,IPrint,.False.,UseSym,NAtoms,IAn,
     $    AtmChg,C,0,Scr,.True.,TRot,XQuad,Scr)
        Steped = Steped .or. OK
        endIf
C
C     Print reason for skipping linear search here as DXQuad prints
C     the eigenvlaues of second derivative matrix.
C
      If(Linr.and.DidLin) then
        If(ICase.eq.1) then
          Write(IOut,2041) XLin
        else if(ICase.eq.2) then
          Write(IOut,2042) XLin
        else if(ICase.eq.3) then
          Write(IOut,2043) XLin
        else
          Write(IOut,2090) ICase
          Call Lnk1E(0)
          endIf
      else if(Linr) then
        Write(IOut,2050)
      else
        If(.not.DoLin) then
          Write(IOut,1220)
        else if(NP.le.1) then
          Write(IOut,1224)
        else if(Neg.ne.0.and.Rises.and.(IPSav1.eq.0.or.IPSav2.eq.0))
     $  then
          Write(IOut,1222)
        else if(.not.Rises.and.(MethQ.eq.1.or.Neg.ne.0)) then
          Write(IOut,1227)
        else if(FRMS.le.FSwtch) then
          Write(IOut,1226) FSwtch
        else
          Call GauErr('Linear search skipped for unknown reason.')
          endIf
        endIf
C
C     If both the linear and quadratic steps have failed, step in
C     the direction of the gradient.  Test size of displacement
C     and update position vector values stored in a for printing.
C
      Despar = Neg.eq.0.and.Linr.and..not.DidLin.and..not.Steped
     $         .and.Rises
      If(Neg.eq.0.and..not.Steped) then
        Write(IOut,2070)
        Call DXGrad(NVar,FTemp,XQuad,DXRMS,DXMax,OK)
        Steped = Steped.or.OK
        endIf
      Call AMask(NVar,-100000,-1,IFlags,NVarM,XQuad,XQuad)
      DXMax = ArrMax(XQuad,NVar,.True.)
      DXSq = SProd(NVar,XQuad,XQuad) / GFloat(NVar)
      DXRMS = Sqrt(DXSq)
      DX = DXRMS * RtNV
      If(DX.gt.DXMaxT.and..not.DidSph) then
        Scale = DXMaxT / DX
        Write(IOut,1240) DXMaxT, Scale
        Call AScale(NVar,Scale,XQuad,XQuad)
        endIf
      DXMax = Zero
      DXRMS = Zero
      DXQR = Sqrt(SProd(NVarM,XQuad,XQuad)/GFloat(NVar))
      DXQMax = ArRMax(XQuad,NVar,.True.)
      DXTest = Cnvfx*ConvF
      If(Despar.and.DXQR.lt.SteepM.and.(.not.Finish.or.
     $   DXQR.gt.DXTest.or.DXQMax.gt.(CnvFMx*Cnvfx*ConvF))) then
        Write(IOut,2080) SteepM
        Call AScale(NVar,(SteepM/DXQR),XQuad,XQuad)
        endIf
      Return
      End
*Deck PathUp
      Subroutine Pathup(IOut,IPrint,NVar,MaxSv1,NVTT,NPath,NPTS,XP,XXP,
     $  FFP,FCP,ICP,ESvP,FSP,IRisP,XX,FF,IC,ESave,Scr,cnstr1,cnstr2,
     $  NTBond,NTAng,NTDih,Proj,DoReac,DoProd,mthUpd,DoIrc)
      Implicit Real*8(A-H,O-Z)
      Dimension XXP(NVar,MaxSv1,NPath),FFP(Nvar,MaxSv1,NPath)
      Dimension FCP(NVTT,NPath),ICP(MaxSv1,NPath),ESvP(NPath)
      Dimension FSP(MaxSv1,NPath),IRisP(NPath),XP(Nvar,NPath)
      Dimension XX(NVar,*),FF(*),IC(*),ESave(*),Scr(NVar,8)
      Dimension Cnstr1(*),Cnstr2(*),Proj(NVTT,NPath)
      Logical Prnt,Exit,HavDEP,Doreac,DoProd,DoIrc,Large
      Common /OptGrd/ ConvF,FMaxT,DXMaxT,RMax,RMin,RLim,EigMax,EigMin,
     $  DEPred,FSwtch,FncErr,GrdErr,FncCnv,XXIRC,NLEVEL,NStep,IStep,
     $  NMax,NP,Neg,Prnt,Exit,IPSav1,IPSav2,HavDEP,NGDIIS,IRisen,MaxRis,
     $  Large,ILevel,ISStep,IOffs,IMDS,MaxSv,MaxSv2,NVrMMT
      Save Zero
      Data Zero/0.0d0/
C
C     Driver for the hessian update for path relaxation.  Update along
C     the path and then update using the history of each point.
C
      MethHi = 0
C     reactant valley
      Do 10 NN = 1, (NPts-1)
C
C       Choose between Powell or mixed update for updating using
C       the point's history. The update along the path is made using
C       Powell's update except for minima for which BFGS is used
C
        If(MthUpd.eq.4) MethHi = 6
        If(MthUpd.eq.3) MethHi = 3
        If(NN.gt.1) then
          Asv0 = EsvP(NN)
          AsvM = EsvP(NN-1)
          AsvP = EsvP(NN+1)
          Call UpDate(IOut,IPrint,NVar,MaxSv1,NVTT,NTBond,NTAng,NTDih,
     $      NPath,NN,NPTS,XXP(1,1,NN),XP(1,NN-1),XP(1,NN+1),XX,ASv0,
     $      AsvM,AsvP,FSP(1,NN),ESave,FFP(1,1,NN),FFP(1,1,NN-1),
     $      FFP(1,1,NN+1),FF,FCP(1,NN),Proj(1,NN),ICP(1,NN),IC,
     $      IRisP(NN),Scr,.False.,.true.,6,MethHi,4)
        else if(NN.eq.1) then
          If(DoReac) MethHi = 5
          RJunk = Zero
          Asv0 = EsvP(NN)
          AsvM = RJunk
          AsvP = EsvP(NN+1)
          If(.not.Doreac) then
            Call UpDate(IOut,IPrint,NVar,MaxSv1,NVTT,NTBond,NTAng,NTDih,
     $        NPath,NN,NPTS,XXP(1,1,NN),Cnstr1,XP(1,NN+1),XX,ASv0,AsvM,
     $        AsvP,FSP(1,NN),ESave,FFP(1,1,NN),Scr(1,8),FFP(1,1,NN+1),
     $        FF,FCP(1,NN),Proj(1,NN),ICP(1,NN),IC,IRisP(NN),Scr,
     $        .false.,.true.,6,MethHi,4)
          else
C           Do BFGS for minimum
            Call UpDate(IOut,IPrint,NVar,MaxSv1,NVTT,NTBond,NTAng,NTDih,
     $        NPath,NN,NPTS,XXP(1,1,NN),Cnstr1,XP(1,NN+1),XX,ASv0,AsvM,
     $        AsvP,FSP(1,NN),ESave,FFP(1,1,NN),Scr(1,8),FFP(1,1,NN+1),
     $        FF,FCP(1,NN),Proj(1,NN),ICP(1,NN),IC,IRisP(NN),Scr,.true.,
     $        .true.,5,MethHi,4)
            endIf
          endIf
   10   Continue
C
C Do TS
      NN = NPTS
      If(MthUpd.eq.4) MethHi = 4
      If(MthUpd.eq.3) MethHi = 7
      If(NN.gt.1.and.NN.lt.Npath) then
        Asv0 = EsvP(NN)
        AsvM = EsvP(NN-1)
        AsvP = EsvP(NN+1)
        MethHi=4
        Call UpDate(IOut,IPrint,NVar,MaxSv1,NVTT,NTBond,NTAng,NTDih,
     $    NPath,NN,NPTS,XXP(1,1,NN),XP(1,NN-1),XP(1,NN+1),XX,ASv0,AsvM,
     $    AsvP,FSP(1,NN),ESave,FFP(1,1,NN),FFP(1,1,NN-1),FFP(1,1,NN+1),
     $    FF,FCP(1,NN),Proj(1,NN),ICP(1,NN),IC,IRisP(NN),Scr,.true.,
     $    .true.,4,MethHi,4)
      else if(NN.eq.Npath) then
        Asv0 = EsvP(NN)
        AsvM = EsvP(NN-1)
        AsvP = RJunk
        MethHi=4
        Call UpDate(IOut,IPrint,NVar,MaxSv1,NVTT,NTBond,NTAng,NTDih,
     $    NPath,NN,NPTS,XXP(1,1,NN),XP(1,NN-1),Cnstr2,XX,ASv0,AsvM,AsvP,
     $    FSP(1,NN),ESave,FFP(1,1,NN),FFP(1,1,NN-1),Scr(1,8),FF,
     $    FCP(1,NN),Proj(1,NN),ICP(1,NN),IC,IRisP(NN),Scr,.true.,
     $    .false.,4,MethHi,4)
      else if(NN.eq.1) then
        Asv0 = EsvP(NN)
        AsvM = RJunk
        AsvP = EsvP(NN+1)
        Call UpDate(IOut,IPrint,NVar,MaxSv1,NVTT,NTBond,NTAng,NTDih,
     $    NPath,NN,NPTS,XXP(1,1,NN),Cnstr1,XP(1,NN+1),XX,ASv0,AsvM,AsvP,
     $    FSP(1,NN),ESave,FFP(1,1,NN),Scr(1,8),FFP(1,1,NN+1),FF,
     $    FCP(1,NN),Proj(1,NN),ICP(1,NN),IC,IRisP(NN),Scr,.true.,
     $    .false.,4,MethHi,4)
        endIf
      NSTop = NPath
      NStart = NPTS + 1
      N1 = 1
C     Do product valley
      Do 20 NN = NStart, NStop, N1
        If(IStep.eq.1.and.DoIrc.and.NN.eq.NStart)
     $    Call AMove(NVTT,FCP(1,NN-N1),FCP(1,NN))
        If(MthUpd.eq.4) MethHi = 6
        If(MthUpd.eq.3) MethHi = 3
        If(NN.lt.NStop) then
          Asv0 = EsvP(NN)
          AsvM = EsvP(NN-1)
          AsvP = EsvP(NN+1)
          Call UpDate(IOut,IPrint,NVar,MaxSv1,NVTT,NTBond,NTAng,NTDih,
     $      NPath,NN,NPTS,XXP(1,1,NN),XP(1,NN-1),XP(1,NN+1),XX,ASv0,
     $      AsvM,AsvP,FSP(1,NN),ESave,FFP(1,1,NN),FFP(1,1,NN-1),
     $      FFP(1,1,NN+1),FF,FCP(1,NN),Proj(1,NN),ICP(1,NN),IC,
     $      IRisP(NN),Scr,.False.,.true.,6,MethHi,4)
        else if(NN.eq.NPATH) then
          If(DoProd) MethHi = 5
          Asv0 = EsvP(NN)
          AsvM = EsvP(NN-1)
          AsvP = RJunk
          If(.not.DoProd) then
            Call UpDate(IOut,IPrint,NVar,MaxSv1,NVTT,NTBond,NTAng,NTDih,
     $        NPath,NN,NPTS,XXP(1,1,NN),XP(1,NN-1),Cnstr2,XX,ASv0,AsvM,
     $        AsvP,FSP(1,NN),ESave,FFP(1,1,NN),FFP(1,1,NN-1),Scr(1,8),
     $        FF,FCP(1,NN),Proj(1,NN),ICP(1,NN),IC,IRisP(NN),Scr,
     $        .false.,.True.,6,MethHi,4)
          else
            Call UpDate(IOut,IPrint,NVar,MaxSv1,NVTT,NTBond,NTAng,NTDih,
     $        NPath,NN,NPTS,XXP(1,1,NN),XP(1,NN-1),Cnstr2,XX,ASv0,AsvM,
     $        AsvP,FSP(1,NN),ESave,FFP(1,1,NN),FFP(1,1,NN-1),Scr(1,8),
     $        FF,FCP(1,NN),Proj(1,NN),ICP(1,NN),IC,IRisP(NN),Scr,.true.,
     $        .True.,5,MethHi,4)
            endIf
          endIf
   20   Continue
      Return
      End
*Deck PrjDag
      Subroutine PrjDag(IOut,IPrint,NTrRot,NAtoms,NVar,AtmChg,C,NFroz,
     $  Frozen,NEig,Eigen,Vec,FCTemp,Scr,VecInt,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Diagonalize the Hessian for cartesian optimization, projecting out
C     translations and rotations.
C
      Logical Frozen(*), UseAtC, SchOr2
      Parameter (UseAtC=.True.)
      Dimension AtmChg(*), C(3,*), Eigen(NVar), Vec(NVar,NVar), JJ(1),
     $  FCTemp(*), Scr(NVar,7), COM(3), PMom(3), PVec(3,3), V(*),
     $  VecInt(NVar,NVar), XX(1)
      Save Zero, One, Tol, JJ, XX
      Data Zero/0.0d0/, One/1.0d0/, Tol/1.d-5/, JJ/0/, XX/0.0d0/
 1000 Format(' DiagFC:  COM=',3D20.10)
 1010 Format(' Eigenvectors of Hessian in projected basis:')
 1020 Format(' Eigenvectors of Hessian in Cartesian basis:')
C
      Call TStamp(1,'Top PrjDag')
      NAt3 = 3*NAtoms
      If(UseAtC) then
        Call AMove(NAtoms,AtmChg,Scr)
      else
        Do 10 I = 1, NAtoms
   10     Scr(I,1) = One
        endIf
      If(NTrRot.eq.-3) then
        NTrRo = 3*NFroz
        NEig = NAt3 - NTrRo
        Call FrezFC(NAtoms,NAt3,NEig,Frozen,VecInt,FCTemp,Scr(1,2))
      else
        Call MofI(NAtoms,C,Scr,COM,PMom,PVec)
        If(IPrint.ge.3) Write(IOut,1000) COM
        Call TRVect(IOut,0,NAtoms,Scr,C,COM,PVec,NTrRo,VecInt)
        Call TStamp(1,'Aft TrVect')
        If(.not.SchOr2(IOut,IPrint,NAt3,NTrRo,NAt3,NRet,.True.,VecInt,V,
     $    MDV)) Call GauErr('SchOr2 failed for PrjDag.')
        Call TStamp(1,'Aft SchOr2')
        If(NRet.ne.NAt3) Call GauErr('SchOr2 failed for PrjDag.')
        If(NTrRot.gt.0) NTrRo = Min(NTrRot,NTrRo)
        If(NVar.ne.NAt3) Call GauErr('PrjDag is confused about NVar.')
        NEig = NVar - NTrRo
        Call AMove(NVar*NTrRo,VecInt,Scr)
        Call AMove1(NVar*NEig,NVar*NTrRo,0,VecInt)
        Call AMove(NVar*NTrRo,Scr,VecInt(1,NEig+1))
        If(IPrint.ge.3) Call OutMtS(IOut,'Proj. vectors in PrjDag:',0,0,
     $    VecInt,NVar,NVar,NVar,NVar)
        Call VFVMul(.False.,.False.,.False.,0,NVar,NVar,NVar,FCTemp,
     $    VecInt,Vec,FCTemp,V,MDV)
        Call TStamp(1,'Aft VFVMul')
        endIf
      If(IPrint.ge.3)
     $  Call LTOutS(IOut,'FC in projected basis:',0,NEig,FCTemp,0)
      NE2 = MDV / NVar
      Call DiagDN(IOut,IPrint,0,FCTemp,Vec,Eigen,NEig,Scr,NE2,V,NVar,
     $  .False.)
      Call TStamp(1,'Aft DiagDN')
      If(IPrint.ge.3) then
        Write(IOut,1010)
        Call MatPt1(IOut,Vec,NVar,NVar,NEig,NEig,0,0,JJ,JJ,0,Eigen,1)
        endIf
      NEigO = NEig
      If(NTrRo.eq.0) then
        N = 0
        Do 20 I = 1, NEig
          If(Abs(Eigen(I)).ge.Tol) then
            N = N + 1
            Eigen(N) = Eigen(I)
            If(N.lt.I) Call AMove(NVar,Vec(1,I),Vec(1,N))
            endIf
   20     Continue
        NEig = N
        endIf
      NVE = NVar*NEig
      Call TstCor(NVE,MDV,'PrjDag-MM')
      Call XGEMM(1,'N','N',NVar,NEig,NEigO,One,VecInt,NVar,Vec,NVar,
     $  Zero,V,NVar)
      Call AMove(NVE,V,Vec)
      Call FixPh1(NVar,NVar,NEig,0,0,1,1,Vec,XX,XX)
      If(IPrint.ge.3) then
        Write(IOut,1020)
        Call MatPt1(IOut,Vec,NVar,NVar,NVar,NEig,0,0,JJ,JJ,0,Eigen,1)
        endIf
      Call TStamp(1,'Bot PrjDag')
      Return
      End
*Deck PrtCDs
      Subroutine PrtCDs(IOut,IPrint,IRwSvC,NAtoms,IAn,IAtTyp,ToAng,
     $  TrVec,TRot,C,V)
      Implicit Real*8(A-H,O-Z)
C
C     Print information about the Cartesian displacements from
C     the original structure.
C
      Parameter (MinPrt=2)
      Dimension IAn(NAtoms), IAtTyp(NAtoms), TrVec(3), TRot(3,3),
     $  C(3,NAtoms), V(3,NAtoms)
      Save Zero
      Data Zero/0.0d0/
 1000 Format(' Largest change from initial coordinates is atom',I5,
     $  F12.3,' Angstoms.')
 1010 Format(' No saved initial coordinates found by PrtCDs.')
C
      If(ITqry(IRwSvC).lt.1) then
        If(IPrint.ge.MinPrt) Write(IOut,1010)
      else
        Call FileIO(2,-IRwSvC,3*NAtoms,V,0)
        If(IPrint.ge.MinPrt) then
          Call CorPrS(IOut,'Initial coordinates:',0,NAtoms,IAn,.True.,
     $      IAtTyp,V,ToAng)
          Call CorPrS(IOut,'Current coordinates before tr-rot:',0,
     $      NAtoms,IAn,.True.,IAtTyp,C,ToAng)
          endIf
        Call Transl(1,3,NAtoms,TrVec,C)
        Call Trspn2(3,TRot)
        Call RotF1(NAtoms,TRot,C)
        Call Trspn2(3,TRot)
        If(IPrint.ge.MinPrt) then
          Call OutMtS(IOut,'Translation vector:',0,1,TrVec,3,1,3,1)
          Call OutMtS(IOut,'Rotation matrix:',0,1,TRot,3,3,3,3)
          Call CorPrS(IOut,'Current coordinates after rotation:',0,
     $      NAtoms,IAn,.True.,IAtTyp,C,ToAng)
          endIf
        Call ASub(3*NAtoms,C,V,V)
        Call RotF1(NAtoms,TRot,C)
        Call Transl(-1,3,NAtoms,TrVec,C)
        If(IPrint.ge.MinPrt) Call CorPrS(IOut,
     $    'Displacements from initial coordinates:',0,NAtoms,IAn,.True.,
     $    IAtTyp,V,ToAng)
        DMax = Zero
        IMax = 0
        Do 10 I = 1, NAtoms
          D = V(1,I)**2 + V(2,I)**2 + V(3,I)**2
          If(D.gt.DMax) then
            DMax = D
            IMax = I
            endIf
   10     Continue
        D = Sqrt(D)*ToAng
        If(IPSave(0).eq.0.or.IPrint.gt.0) Write(IOut,1000) IMax, D
        endIf
      Return
      End
*Deck PTHARC
      Subroutine PTHARC(IOut,S,FTEMP,NVAR,FRE,FINPR,PVEC,PPRIME,PPER,
     $  PIPPR,PTIL,FTIL,PIPTL,FTMP2)
      Implicit Real*8(A-H,O-Z)
C
C     Linear search routine.
C
      Logical Debug
      Dimension FTEMP(*),FRE(*),FINPR(*),PVEC(*),PPRIME(*),PPER(*),
     $  PIPPR(*),PTIL(*),FTIL(*),PIPTL(*),FTMP2(*)
      Save Zero,Pt5,One
      Data Zero,Pt5,One/0.0D0,0.5D0,1.0d0/
 1030 Format(2X,3F10.5)
 1040 Format(' X  ,  DX,  XNEW ')
 1070 Format(' PPRIME, P, PTILDA ')
 1090 Format(' LINEAR SEARCH ALONG THE ARC HAS BEEN TURNED ON ')
 1130 Format(2X,' FPRIME, F, FTILDA ')
 1140 Format(' GFPRIME-(-PPRIME) ANGLE= ',D16.8)
 1150 Format(' G*F-(-P) ANGLE= ',D16.8)
C
      WRITE(IOUT,1090)
      DEBUG=.FALSE.
C
C     CALCULATION OF PPR AND FPR
C
C     PPRIME=SEARCH VECTOR AT THE PREVIOUS POINT.
C     FF(I,2)=THE FORCE AT THE PREVIOUS POINT.
C     PVEC=THE SEARCH VECTOR AT THE CURRENT POINT.
C     FTEMP=THE FORCE AT THE CURRENT POINT.
C
C
      HSTEP=PT5*S
      Call AScale(Nvar,HStep,PVec,PVec)
      Call AScale(Nvar,HStep,PPRime,PPRime)
      Call AScale(Nvar,-One,Ftemp,Ftemp)
      Call AScale(NVar,-One,FTMP2,FTMP2)
C
C     CALCULATION OF THETAPRIME
C
      THETPR=ANGDEG(1,NVAR,PVEC,PPRIME)
      XCOP=Cos(THETPR)
      XSIP=Sin(THETPR)
      FP=ZERO
      FPAR=ZERO
      FPIPR=ZERO
      FPARPR=ZERO
C
C     CALCULATION OF PER. AND PAR. COMPONENTS OF PVEC AND PPRIME.
C
      FIRSTF=ZERO
      SECONF=ZERO
      DO 20 I=1,NVAR
      PPER(I)=(-PVEC(I)*XCOP/XSIP)+PPRIME(I)/XSIP
      PIPPR(I)=PVEC(I)*(-XSIP-XCOP*XCOP/XSIP)+
     $PPRIME(I)*XCOP/XSIP
   20 CONTINUE
      FIRSTF=SPROD(NVAR,PPER,PPER)
      SECONF=SPROD(NVAR,PIPPR,PIPPR)
      PSQ=SPROD(NVAR,PVEC,PVEC)
      PPSQ=SPROD(NVAR,PPRIME,PPRIME)
      FIRSTF=(PT5*S)/(Sqrt(FIRSTF))
      SECONF=(PT5*S)/(Sqrt(SECONF))
      PSQ=(PT5*S)/(Sqrt(PSQ))
      PPSQ=(PT5*S)/(Sqrt(PPSQ))
      DO 22 I=1,NVAR
      FRE(I)=-PVEC(I)
      FINPR(I)=-PPRIME(I)
      FP=FP+FTEMP(I)*PPER(I)*FIRSTF
      FPAR=FPAR+FTEMP(I)*PVEC(I)*PSQ
      FPIPR=FPIPR+FTMP2(I)*PIPPR(I)*SECONF
      FPARPR=FPARPR+FTMP2(I)*PPRIME(I)*PPSQ
   22 CONTINUE
      ANGLE1=ANGDEG(0,NVAR,FTEMP,FRE)
      ANGLE2=ANGDEG(0,NVAR,FTMP2,FINPR)
      THETL=-FP*THETPR/(FPIPR-FP)
      XCOPP=Cos(THETL)
      XSIPP=Sin(THETL)
      U=ZERO
      IF(DEBUG) THEN
      WRITE(IOUT,1150) ANGLE1
      WRITE(IOUT,1140) ANGLE2
      WRITE(IOUT,1070)
      WRITE(IOUT,'('' ****************************** '')')
      WRITE(IOUT,'('' G1*F = '')')
      WRITE(IOUT,23) (FTEMP(I),I=1,NVAR)
   23 Format(2X,3F14.9)
      WRITE(IOUT,'('' G1*FPRIME = '')')
      WRITE(IOUT,23) (FTMP2(I),I=1,NVAR)
      WRITE(IOUT,'(''P PERPENDICULAR = '')')
      WRITE(IOUT,23) (PPER(I),I=1,NVAR)
      WRITE(IOUT,'('' PPRIME PERP. = '')')
      WRITE(IOUT,23) (PIPPR(I),I=1,NVAR)
      WRITE(IOUT,24)FP
   24 Format(2X,'FP =',F14.9)
      WRITE(IOUT,25)FPIPR
   25 Format(2X,'FPIPR=',F14.9)
      WRITE(IOUT,'('' ********************************** '')')
      ENDIF
C
C     CALCULATION OF PTILDA AND PIPTL.
C
      WRITE(IOUT,26)THETPR,THETL
   26 Format(2X,'THETAPR=',F14.9,5X,'THETATIL=',F14.9)
      If(Thetl.gt.ThetPr)then
      Call AClear(Nvar,FTil)
      Call AMove(Nvar,PVec,PTil)
      Write(IOut,*)' Linear step rejected '
      Return
      Endif
      DO 30 I=1,NVAR
      PTIL(I)=PVEC(I)*(XCOPP-XSIPP*XCOP/XSIP)+
     $PPRIME(I)*XSIPP/XSIP
      PIPTL(I)=PVEC(I)*(-XSIPP-XCOPP*XCOP/XSIP)
      PIPTL(I)=PIPTL(I)+PPRIME(I)*XCOPP/XSIP
      U=U+PTIL(I)*PTIL(I)
      IF(DEBUG) WRITE(IOUT,1030)PPRIME(I),PVEC(I),PTIL(I)
   30 CONTINUE
      U=HSTEP/Sqrt(U)
      FIRSTT=SPROD(NVAR,PTIL,PTIL)
      SECONT=SPROD(NVAR,PIPTL,PIPTL)
      FIRSTT=(PT5*S)/(Sqrt(FIRSTT))
      SECONT=(PT5*S)/(Sqrt(SECONT))
C
C     CALCULATION OF FRES., FRESPRIME AND FTILDA.
C
      IF(DEBUG) WRITE(IOUT,1130)
      DO 40 J=1,NVAR
      FRE(J)=FTEMP(J)-FPAR*PVEC(J)*PSQ-FP*PPER(J)*FIRSTF
      FINPR(J)=FTMP2(J)-FPARPR*PPRIME(J)*PPSQ-
     $FPIPR*PIPPR(J)*SECONF
      CONE=FPAR+(FPARPR-FPAR)*THETL/THETPR
      CONE=CONE*FIRSTT
      CTWO=FP+(FPIPR-FP)*THETL/THETPR
      CTWO=CTWO*SECONT
      CTHREE=FRE(J)+(FINPR(J)-FRE(J))*THETL/THETPR
      FTIL(J)=CONE*PTIL(J)+CTWO*PIPTL(J)+CTHREE
      FTEMP(J)=FTIL(J)
   40 CONTINUE
      Call AScale(Nvar,-One,Ftemp,Ftemp)
      IF(DEBUG) WRITE(IOUT,1040)
      DO 50 I=1,NVAR
   50   FTil(I)=PTIL(I)-PVec(I)
      Return
      End
*Deck PthOpt
      Subroutine PthOpt(IOut,IPrint,UpTrus,Sphere,DoLin,MethQ,FCart,
     $  UpdRWF,ToAng,IDimSc,F,Frcnst,X,XName,FC,FC2,ESave,FS,XX,FF,IC,
     $  ITU,XQuad,XNew,FTemp,A,VecInt,Scr,FCTemp,APrint,IFlags,NAtoms,
     $  IAn,AtmChg,C,UseSym,TRot,MthUpd,Wght,IQST,DoRed,NTrRot,NTBond,
     $  NTAng,NTDih,IZRed,Value,Cnstr1,Cnstr2,CNew,IBMat,BMat,NPath,
     $  NCurPt,NNext,MaxSv1,NVTT,ICP,FFP,XXP,FCP,FC2P,ESvP,FSP,Q0,XP,Q,
     $  DQ,GradPt,PathL,VecInP,AP,EigP,NPStar,IRisP,CP,IConv,dxmxt0,
     $  doReac,doprod,doMid,dobirx,DoIrc,SS,TRotP,IRis,Scr2,NOpAlP,
     $  SymAlp,IPrAlP,TanP,ArcP,Stot,IOSLE1,IAprBG,IEStpM,NFroz,Frozen,
     $  V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     NPath .... Number of discrete points
C     NCurPt .... Index of point for which energy and gradient
C                 have just been computed.
C     NNext .... Index of point for which energy and gradient
C                are going to be computed next.
C     MAxSv1 .... Max. # of points saved.
C     NVTT .... Length of Hessian matrix
C     ??P .... Same as ?? but for each discrete point unless
C              specified otherwise.
C     Q,Q0 .... Scratch for path coordinates
C     GradPt .... Scratch for Path forces
C     PathL .... Storage for path length
C     VecInP ... Used to keep a copy of each projector and for
C                 and to store eigenvectors.
C     EigP   .... Used to store eigenvalues.
C     IConv .... keeps track of calculations along the path.
C     DXmxT0 ... Original trust radius.
C     DoReac ... Whether to optimize reactant
C     DoProd ... Whether to optimize products
C     DoMid .... Whether to force the middle point to optimize
C                to TS.
C     DoIrc .... The first point is the TS and the transition vector is
C                given by hessian eigenvector
C     SS .... Storage for Trust Radius
C
      Integer GetChr
      Parameter (MaxOp=DEFMAXOP,IOPath=601)
      Logical OK, Steped, Prnt, Exit, Rises, FixEig, LJunk, UpTrus,
     $  Sphere, Failed, DoLin, UpdRWF, Finish, SclEig, FCart,
     $  HavDep, UseSym, DoRed, Swap, DidGI, Do601, Large, DidSwap,
     $  DoReac,DoProd,DoneTS,DoMid,DoX,DoBiRx,DoIrc,DoLin2,PthSym,
     $  UseSyT, Frozen(*)
      Dimension F(1), Frcnst(1), X(NVar), XName(1), FC(1), FS(1),
     $  XX(NVar,*), FF(*), IC(1), XQuad(1), XNew(NVar), TRot(3,3),
     $  FTemp(1), A(*), ESave(*), Scr(IDimSc,8), FCTemp(1), XXA(1),
     $  APrint(6,*), FC2(1), IFlags(1), C(3,*), AtmChg(*), ITU(*),
     $  VecInt(*),  Cnstr1(NVar), Cnstr2(NVar), TrRMat(6), RR(1),
     $  IZRed(4,*), Value(*), CNew(3,NAtoms), IBMat(4,*), BMat(3,4,*),
     $  IAn(*)
      Dimension ICP(MaxSv1,NPath),FFP(NVAr,MaxSv1,NPAth),
     $  XXP(NVAr,MaxSv1,NPath),FCP(NVTT,NPAth),FC2P(NVTT,NPAth),
     $  ESvP(NPAth),FSP(MAxSv1,NPath),Q0(NVAr,NPath),XP(NVar,NPath),
     $  Q(NVar,NPath),DQ(NVar,NPAth),GradPt(NVar,NPath),PathL(NPath),
     $  VecInP(NVTT,NPath),AP(NVar*NVar,NPath),EigP(NVar,NPath),
     $  IRisP(NPath),CP(3,NAtoms,2,NPath),IConv(NPath),SS(Npath),
     $  TRotP(3,3,NPath),IRis(4),Scr2(Nvar,9),NOpAlP(NPath),
     $  SymAlP(3,3,MaxOp,NPath),IPrAlP(NAtoms,MaxOp,NPath),
     $  TanP(Nvar,NPath),ArcP(Npath),STot(2),V(*)
      Common /Grdnt/  Energy,NVar,IGetFC
      Common /OptGrd/ ConvF,FMaxT,DXMaxT,RMax,RMin,RLim,EigMax,EigMin,
     $  DEPred,FSwtch,FncErr,GrdErr,FncCnv,XXIRC,NLEVEL,NStep,IStep,
     $  NMax,NP,Neg,Prnt,Exit,IPSav1,IPSav2,HavDEP,NGDIIS,IRisen,MaxRis,
     $  Large,ILevel,ISStep,IOffs,IMDS,MaxSv,MaxSv2,NVrMMT
      Dimension IReslt(3),ITmp(20),IStr(20)
      Save Zero,One,CnvFMx,CNVFX,Four,Pt5,XXA
      Data Zero,One,CnvFMx,CNVFX/0.D0,1.0d0,1.5D0,4.D0/,XXA/0.0d0/
     $  Four/4.0d0/,Pt5/0.5d0/
 1000 Format('         Item               Value     Threshold  ',
     $       'Converged?')
 1002 Format(' Maximum Force            ',F8.6,5X,F8.6,5X,3A1)
 1004 Format(' RMS     Force            ',F8.6,5X,F8.6,5X,3A1)
 1006 Format(' Maximum Displacement     ',F8.6,5X,F8.6,5X,3A1)
 1008 Format(' RMS     Displacement     ',F8.6,5X,F8.6,5X,3A1)
 1016 Format(' Search for a local minimum.')
 1117 Format(' Search for a reaction path.')
 1018 Format(' Search for a stationary point of order',I3)
 1120 Format(' Pass number',I4,' out of a maximum of',I4,/,
     $' All quantities printed in internal units',
     $' (Hartrees-Bohrs-Radians)')
 1050 Format(' Variable       Old X    -DE/DX',3(3X,7HDelta X),
     $ '     New X',/,33X,'(Linear)    (Quad)   (Total)')
 1052 Format(' Variable       Old X    -DE/DX    -DE/DX',
     $2(3X,7HDelta X),
     $ '     New X',/,33X,'(Projected) (Quad)   (Total)')
 1060 Format(1X,8A1,2X,6F10.5)
 1200 Format(' Second derivative matrix not updated -- first step.')
 1210 Format(' Second derivative matrix not updated --',
     $       ' analytic derivatives used.')
 2010 Format(' Optimization aborted.',/,
     $       '    -- Gradient out of range.',/,
     $       '    -- Maximum allowed force = ',F8.3)
 2100 Format(' MatCor failed.')
 3000 Format(' No relaxation : TS found and path stable')
 3060 Format(' Energy along the path')
 3064 Format('Min E(',I2,') = ',G19.12)
 3065 Format('    E(',I2,') = ',G19.12)
 3066 Format('TS  E(',I2,') = ',G19.12)
 3070 Format(' Minimization step for reactants')
 3080 Format(' Minimization step for products')
 3090 Format(' QST Step for point ',I2)
 4000 Format(' Relaxation step number ',I3)
 4010 Format(' Point ',I3)
 4012 Format(' Point ',I3,' S = ',F10.5)
 4020 Format('+++++++++++++')
 4022 Format('+++++++++++++++++++++++++++++')
 4030 Format(' Number of gradient calculations for point ',I3,
     $' so far = ',I4)
 4040 Format(' Trust radius = ',F10.5)
 4050 Format(I6,' Points converged ')
 5000 Format(' Prepare gradient calculation for point ',I3)
 5020 Format(1x,'        ',9x,'Path converged')
 5030 Format(1x,'        ',9x,'Path not converged')
 5040 Format(10x,I6,' Steps to optimize ts')
 5050 Format(10x,I6,' Steps to optimize reactants')
 5060 Format(10x,I6,' Steps to optimize products')
 5070 Format(10x,F6.2,' Steps to optimize remaining points')
 5080 Format(10x,I6,' Total steps')
 5999 Format(1X,19(4H++++))
C
      Do601 = .false.
      IGISav = 0
      NAt3 = 3*NAtoms
      IStep = IStep + 1
      Call SetL0(IStep)
      If(IPSave(0).le.1) then
        If(Neg.eq.0) then
          Write(IOut,1016)
          Call Gauerr(' Wrong route for path relaxation')
        else if(Neg.eq.1) then
          Write(IOut,1117)
        else
          Write(IOut,1018) Neg
          Call Gauerr(' Wrong route for path relaxation')
          endIf
        endIf
      Write(IOut,1120) IStep, NStep
      NNegE = Neg
      Failed = .False.
C     Use symmetry for path relaxation
      PthSym=.true.
      Call IClear(NVar,IFlags)
      Call AMask(NVar,-100000,-1,IFlags,NVarM,F,FTemp)
C     Dored is assumed true for the relaxation.
C     intended, for the present version, to stay clear of the cartesian
C     specific routines.
      If(.not.Dored) Call GauErr(' Red. Internal Required ')
C     Make sure Value contains X
      Call AMove(Nvar,X,Value)
      IGISav = 0
      Call Red2BG(IOut,IPrint,0,NAtoms,NTrRot,IGISav,IAprBG,C,NVar,
     $  NTBond,NTAng,NTDih,IZRed,Value,NRank,BMat,IBMat,A,VecInt,FCTemp,
     $  Scr,XName,V,MDV)
      IGISav = IGetRW(NVar**2,0)
      Call FileIO(1,-IGISav,NVar**2,A,0)
      Call AMove(NVTT,FCTemp,VecInt)
      Call CnstrP(IOut,IPrint,NVar,IFlags,-100000,-1,VecInt,A,FCTemp,
     $  Scr)
      Call MultLV(.True.,1,0,NVar,NVar,VecInt,F,FTemp)
C     Save a copy of forces and coordinates
      If(NCurpt.eq.1.and.DoReac) Call AMove(Nvar,X,Cnstr1)
      If(NCurpt.eq.NPath.and.DoProd) Call AMove(Nvar,X,Cnstr2)
C
C     Test gradient.
C
      FMax = ARRMax(FTemp,NVar,.True.)
      FRMS = SPROD(NVar,FTemp,FTemp)
      FRMS = Sqrt(FRMS/NVarM)
      Convx = CNVFX*ConvF
      T1 = ConvF*CnvFMX
      T2 = ConvX*CnvFMX
      Finish = (FRMS.lt.ConvF) .and. (FMax.lt.(CnvFMx*ConvF))
      If(FMax.ge.FMaxT.and.UpdRWF) then
        NCur = 0
        Call PutCSt('-DE/DX',ITmp,NCur)
        Call PutDel(2,ITmp,NCur)
        Call MatPrt(F,1,NVar,1,NVar,1,1,ITmp,XName,0,XXA,0)
        Write(IOut,2010) FMaxT
        Call Lnk1E(0)
        endIf
      Call AClear(NVar,XNew)
      Call AClear(NVar,XQuad)
      If(NPath.gt.0) Swap = .False.
C
C     IRis contains information on whteher energy has risen for
C     reactant and product.
C     IRis is initialized at first .
      If(IStep.eq.1.and.NCurPt.eq.1) Call IClear(4,IRis)
C     set  IRisen info for reactant and/or product.
      If(DoReac.and.NCurpt.eq.1) IRisen = IRIs(2)
      If(DoProd.and.NCurpt.eq.NPath) IRisen=IRIs(4)
C     IRisP contains the number of gradient calc. for this point
      NP = IRisP(NCurPt)
      IStep2 = Istep
      IStep = NP+1
      Call SavePt(Energy,NVar,IGetFC,.False.,0,Rises,F,Frcnst,ESave,IC,
     $  ITU,X,FF,XX,FC,FC2,Swap,DidSwap,1,1,1,RR,RR,RR,RR,RR,RR,RR(1),
     $  RR)
C     Update IRisP for this point
      IRisP(NCurPt) = NP
      Write(IOut,4030) NCurPt,IRisP(NCurPt)
      IStep = IStep2
C     Update IRis for reactant and/or product
      If(DoReac.and.NCurpt.eq.1) then
        If(Rises) then
          IRis(1) = 1
        else
          Iris(1) = 0
          endIf
        IRis(2) = IRisen
      else if(DoProd.and.NCurpt.eq.NPath) then
        If(Rises) then
          IRis(3) = 1
        else
          IRis(3) = 0
          endIf
        IRis(4) = IRisen
        endIf
      If(IGetFC.eq.4) then
        Write(IOut,1210)
      else if(NP.eq.1) then
        Write(IOut,1200)
        endIf
C
C     Update trust radius if allowed.
C
      Steped = .False.
C     Save information of NCurPt
      IConv(NCurPt) = 1
C     Fixsym save or substitute symmetry data pecific to point NCuprPt
C     This is a quick fix that may have a better alternative.
      If(USeSym.and.PthSym) Call FixSym(.True.,NAtoms,NOpAlP(NCurPt),
     $  SymAlP(1,1,1,NCurPt),IPrAlP(1,1,NCurPt))
      Call OptDat(XP,XXP,GradPt,FFP,FCP,FC2P,ESvP,FSP,VecInP,ICP,TRotP,
     $  CP,Value,XX,F,FF,FC,FC2,ESave,FS,VecInt,IC,TRot,C,NAtoms,NVar,
     $  MaxSv1,NVTT,NPath,NCurPt)
      NIStep = Mod(IStep,2)
      Call ChkPth(IOut,NPath,IStep,NConv,IConv,NPStar,DoMid,DoIrc,
     $  DoReac,DoProd,ESvP,STot,PathL)
C Trust radius update and generate empirical hessian
      If(NConv.eq.NPath) then
        If(IStep.gt.1) then
          NStop = 1
          N1 = -1
C         First Valley
          ICount = 0
          Call PthTrs(IOut,IPrint,NVar,NTBond,NTAng,NTDih,NVtt,NPath,
     $      NPStar,NPStar,NStop,N1,MaxSv1,IStep,XNew,XXP,FFP,XQuad,FCP,
     $      ESvP,FSP,DXMxT0,ESave,IRIsP,SS,Doreac,ICount,T2,ConvX)
          NStop = NPath
          NStart = Min(NPStar+1,NPath)
          N1 = 1
C         Second valley
          Call PthTrs(IOut,IPrint,NVar,NTBond,NTAng,NTDih,NVtt,NPath,
     $      NPStar,NStart,NStop,N1,MaxSv1,IStep,XNew,XXP,FFP,XQuad,FCP,
     $      ESvP,FSP,DXMxT0,ESave,IRIsP,SS,DoProd,ICount,T2,ConvX)
          endIf
        If(IStep.lt.2.or..Not.UpTrus) then
          ICount = NPath
          Do 125 II = 1, NPath
  125       SS(II) = DxmxT0
        else
          Call RMSVec(NPath,SS(1),SRMS,SMAx)
C         Get RMS value for the radii of the points that are actually moving.
          If(Icount.gt.0) then
            SRMS = SRMS*Sqrt(GFloat(NPath)/GFloat(ICount))
            Write(IOut,4040) SRMS
          else
            SRMS = DxMxT0
            endIf
          DXMaxT = SRMS
          Do 126 II = 1, NPath
  126       SS(II) = DxMaxT
          endIf
C       Hessian update
      If(IGetFC.ne.4) then
        Call aclear(nvar,DQ(1,NPStar))
C       For the first pass, be generous with a large acceptance criteria
C       for updating along the line.
        If(IStep.eq.1) RMax = RMax*Sqrt(GFLoat(NVar)/GFloat(3*NAtoms-6))
        Call Pathup(IOut,IPrint,NVar,MaxSv1,NVTT,NPath,NPStar,XP,XXP,
     $    FFP,FCP,ICP,ESvP,FSP,IRisP,XX,FF,IC,ESave,Scr2,cnstr1,cnstr2,
     $    NTBond,NTAng,NTDih,VecInP,doreac,doprod,MthUpd,DoIrc)
        If(IStep.eq.1) RMax = RMax*Sqrt(GFloat(3*NAtoms-6)/GFloat(NVar))
        endIf
      Write(Iout,3060)
      Do 130 NPt = 1, NPath
        If((Npt.eq.NPath.and.Doprod).or.(NPt.eq.1.and.DoReac)) then
          Write(Iout,3064)NPt,ESvP(NPT)
        else if(Npt.eq.NPStar) then
          Write(IOut,3066)NPt,ESvP(Npt)
        else
          Write(IOut,3065)NPt,ESvP(NPt)
          endIf
C       Prepare arrays for relaxation
        Call AMove(NVar,XP(1,NPt),Q0(1,NPt))
        Call AMove(NVar,FFP(1,1,NPt),GradPt(1,NPt))
  130   Continue
      Sjunk = Zero
      If(IStep.eq.1) then
        Do 134 Npt = (NPstar-1), 1, -1
          Call ASub(Nvar,XP(1,NPt+1),XP(1,NPt),Scr2(1,1))
          If(DoRed) Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,
     $      Scr2,XXA)
          Call MultLV(.True.,1,0,NVar,NVar,VecInP(1,NPt),Scr2,Scr2(1,2))
          SJunk = SJunk + Sqrt(SProd(Nvar,Scr2,Scr2))
  134     Continue
        Do 136 Npt = (NPstar+1), NPath
          Call ASub(Nvar,XP(1,NPt-1),XP(1,NPt),Scr2(1,1))
          If(DoRed) Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,
     $      Scr2,XXA)
          Call MultLV(.True.,1,0,NVar,NVar,VecInP(1,NPt),Scr2,Scr2(1,2))
          SJunk = SJunk + Sqrt(SProd(Nvar,Scr2,Scr2))
  136     Continue
        STot(1) = SJunk
        endIf
      If(NConv.eq.NPath) then
        If(DoReac.and.NPStar.ne.1) then
C         Optimize reactants; Feed the information for Reactant-like point
          Write(IOut,3070)
          NPtTmp = -1
          Call OptDat(Q0,XXP,GradPt,FFP,FCP,FC2P,ESvP,FSP,VecInP,ICP,
     $      TRotP,CP,Value,XX,F,FF,FC,FC2,ESave,FS,VecInt,IC,TRot,C,
     $      NAtoms,NVar,MaxSv1,NVTT,NPath,NPtTmp)
          Call AMove(NVar,XP(1,1),Cnstr1)
          Call MultLV(.True.,1,0,NVar,NVar,VecInt,F,FTemp)
C         Set Optimization Flags
          Call OptFlag(Neg,MethQ,DoLin,DoLin2,Rises,IRis(1),Steped,IQST,
     $      IQST1,FixEig,SclEig,1)
          IRisen = Iris(2)
          NP = IRisP(1)
          DxMaxT = SS(1)
          Call AClear(Nvar,XQuad)
          Call AClear(Nvar,Xnew)
C         OptStp contains most of all of the stepping routines found in GrdOpt.
          If(IRisP(1).eq.IStep) Call OptStp(IOut,IPrint,Sphere,DoLin2,
     $      MethQ,FCart,UpdRWF,IGetFC,IDimSc,XName,FC,FC2,ESave,FS,XX,
     $      FF,XQuad,XNew,FTemp,A,VecInt,Scr,FCTemp,IFlags,NAtoms,IAn,
     $      AtmChg,C,UseSym,TRot,Wght,IQST1,DoRed,NTBond,NTAng,NTDih,
     $      Nvar,Cnstr1,Cnstr2,NVarM,Rises,NTrRot,IEStpM,NFroz,Frozen)
          Call AClear(NVar,DQ(1,1))
          If(IRisP(1).lt.IStep) then
            Call AClear(Nvar,XNew)
            Call AClear(Nvar,XQuad)
            endIf
          Call AAdd(NVar,XQuad,XNew,DQ(1,1))
          IPrt1 = IPrint
          If(IPSave(0).le.1) IPrt1 = Max(IPrint,1)
          IGISav = 0
C         Update red. internal coordinates
          Call IClear(Nvar,IFlags)
          Call AMove(NAt3,CP(1,1,1,1),C)
          If(USeSym.and.PthSym) Call FixSym(.False.,NAtoms,NOpAlP(1),
     $      SymAlP(1,1,1,1),IPrAlP(1,1,1))
          UseSyT = UseSym.and.PthSym.and.NOpAlP(1).gt.1
          Call UpdRed(IOut,IPrint,IPrt1,Dored,NVar,NTBond,NTAng,NTDih,
     $      IZred,Value,NAtoms,C,DQ(1,1),Cnew,BMat,IBMat,VecInt,FCTemp,
     $      IDimSc,Scr,XNew,IFlags,IAn,AtmChg,TRot,TrRMat,X,Cnstr1,
     $      CnvFX,ConvF,CnvFMX,GradPt(1,1),UseSyT,.False.,CP(1,1,1,1),
     $      FMax,FRMS,DxMaxC,DxRMS,DoneTS,IOSLE1,V,MDV)
          Call AMove(Nvar,Cnstr1,Q(1,1))
          Call AMove(Nvar,Cnstr1,Q0(1,1))
          Call AMove(Nvar,Cnstr1,XP(1,1))
C         Reactant Done
          endIf
      If(DoProd) then
C       Feed the information for Product-like point
C       The following does the same thing as before. This should really be put
C       in one compact routine available for reactant,ts,product.
        Write(IOut,3080)
        NPtTmp = -NPath
        Call OptDat(Q0,XXP,GradPt,FFP,FCP,FC2P,ESvP,FSP,VecInP,ICP,
     $    TRotP,CP,Value,XX,F,FF,FC,FC2,ESave,FS,VecInt,IC,TRot,C,
     $    NAtoms,NVar,MaxSv1,NVTT,NPath,NPtTmp)
        Call AMove(NVar,XP(1,NPath),Cnstr2)
        Call MultLV(.True.,1,0,NVar,NVar,VecInt,F,FTemp)
        Call OptFlag(Neg,MethQ,DoLin,DoLin2,Rises,IRis(3),Steped,IQST,
     $    IQST1,FixEig,SclEig,1)
        IRisen = IRis(4)
        NP = IRisP(NPath)
        DxMaxT = SS(NPath)
        Call AClear(Nvar,XQuad)
        Call AClear(Nvar,Xnew)
        If(IRisP(NPath).eq.IStep) Call OptStp(IOut,IPrint,Sphere,DoLin2,
     $    MethQ,FCart,UpdRWF,IGetFC,IDimSc,XName,FC,FC2,ESave,FS,XX,FF,
     $    XQuad,XNew,FTemp,A,VecInt,Scr,FCTemp,IFlags,NAtoms,IAn,AtmChg,
     $    C,UseSym,TRot,Wght,IQST1,DoRed,NTBond,NTAng,NTDih,Nvar,Cnstr1,
     $    Cnstr2,NVarM,Rises,NTrRot,IEStpM,NFroz,Frozen)
        Call AClear(NVar,DQ(1,NPath))
        If(IRisP(NPath).lt.IStep) then
          Call AClear(Nvar,XNew)
          Call AClear(Nvar,XQuad)
          endIf
        Call AAdd(NVar,XQuad,XNew,DQ(1,NPath))
        IPrt1 = IPrint
        If(IPSave(0).le.1) IPrt1 = Max(IPrint,1)
        IGISav = 0
        Call IClear(Nvar,IFlags)
        Call AMove(NAt3,CP(1,1,1,NPath),C)
        If(USeSym.and.PthSym) Call FixSym(.False.,NAtoms,NOpAlP(NPath),
     $    SymAlP(1,1,1,NPath),IPrAlP(1,1,NPath))
        UseSyT = UseSym.and.PthSym.and.NOpAlP(NPath).gt.1
        Call UpdRed(IOut,IPrint,IPrt1,Dored,NVar,NTBond,NTAng,NTDih,
     $    IZred,Value,NAtoms,C,DQ(1,NPath),Cnew,BMat,IBMat,VecInt,
     $    FCTemp,IDimSc,Scr,XNew,IFlags,IAn,AtmChg,TRot,TrRMat,X,Cnstr2,
     $    CnvFX,ConvF,CnvFMX,GradPt(1,NPath),UseSyT,.False.,
     $    CP(1,1,1,NPAth),FMax,FRMS,DxMaxC,DxRMS,DoneTS,IOSLE1,V,MDV)
        Call AMove(Nvar,Cnstr2,Q(1,NPath))
        Call AMove(Nvar,Cnstr2,Q0(1,NPath))
        Call AMove(Nvar,Cnstr2,XP(1,NPath))
C       Product Done
        endIf
C     Feed the information for QST step on ts-like point
      Write(IOut,3090) NPStar
      NPtTmp = -NPStar
      Call OptDat(Q0,XXP,GradPt,FFP,FCP,FC2P,ESvP,FSP,VecInP,ICP,TRotP,
     $  CP,Value,XX,F,FF,FC,FC2,ESave,FS,VecInt,IC,TRot,C,NAtoms,NVar,
     $  MaxSv1,NVTT,NPath,NPtTmp)
      Call MultLV(.True.,1,0,NVar,NVar,VecInt,F,FTemp)
      Call OptFlag(Neg,MethQ,DoLin,DoLin2,Rises,0,Steped,IQST,IQST1,
     $  FixEig,SclEig,2)
      NP = IRisP(NPStar)
      If(DoIrc) IQST1 = 0
      DxMaxT = SS(NPStar)
      Call AClear(Nvar,XQuad)
      Call AClear(Nvar,Xnew)
C     A few special cases to handle ...
      If(NPStar.gt.1.and.NPStar.lt.NPath) then
        Call AMove(Nvar,Q0(1,NPStar+1),Scr2)
        Call AMove(Nvar,Q0(1,NPStar-1),Scr2(1,2))
      else if(NPStar.eq.1) then
C       the putative TS is the first point of the list.
        Call AMove(Nvar,Q0(1,NPStar+1),Scr2)
        Call AMove(Nvar,Cnstr1,Scr2(1,2))
      else if(NPStar.eq.NPath) then
C       The putative TS is the last point of the list
        Call AMove(Nvar,Cnstr2,Scr2)
        Call AMove(Nvar,Q0(1,NPStar-1),Scr2(1,2))
        endIf
      Call OptStp(IOut,IPrint,Sphere,DoLin2,MethQ,FCart,UpdRWF,IGetFC,
     $  IDimSc,XName,FC,FC2,ESave,FS,XX,FF,XQuad,XNew,FTemp,A,VecInt,
     $  Scr,FCTemp,IFlags,NAtoms,IAn,AtmChg,C,UseSym,TRot,Wght,IQST1,
     $  DoRed,NTBond,NTAng,NTDih,NVar,Scr2,Scr2(1,2),NVarM,Rises,NTrRot,
     $  IEStpM,NFroz,Frozen)
      PathL(NPStar) = Zero
      If(DoIrc) then
C       If we do an IRC type calc. feed the transition vector.
        Call MultLV(.True.,1,0,NVar,NVar,VecInt,FTemp,XNew)
        Call AMove(NVar,XNew,DQ(1,NPStar))
        Dot = One/Sqrt(SProd(Nvar,DQ(1,NPStar),DQ(1,NPStar)))
C       make sure it points in the right direction.
        Call ASub(NVar,XP(1,NPStar+1),XP(1,NPStar),Xnew)
        DotD = SProd(NVar,DQ(1,NPStar),XNew)
        If(DotD.lt.Zero) Dot = -Dot
        Call AScale(NVar,Dot,DQ(1,NPStar),DQ(1,NPStar))
        Call AMove(NVar,DQ(1,NPStar),Scr2(1,4))
        Call AClear(NVar,XNew)
        endIf
C
C     Update TS coordinates
C
      Call AClear(NVar,DQ(1,NPStar))
      Call AClear(NVar,DQ(1,NPStar+1))
      Call AClear(NVar,DQ(1,NPStar-1))
      If(IRisP(NPStar).lt.IStep) then
        Call AClear(Nvar,XNew)
        Call AClear(Nvar,XQuad)
        endIf
      Call AAdd(NVar,XQuad,XNew,DQ(1,NPStar))
      IPrt1 = IPrint
      If(IPSave(0).le.1)IPrt1 = Max(IPrint,1)
      IGISav = 0
      Call IClear(Nvar,IFlags)
      Call AMove(NAt3,CP(1,1,1,NPStar),C)
      If(USeSym.and.PthSym) Call FixSym(.False.,NAtoms,NOpAlP(NPStar),
     $  SymAlP(1,1,1,NPStar),IPrAlP(1,1,NPStar))
      UseSyT = UseSym.and.PthSym.and.NOpAlP(NPStar).gt.1
      Call UpdRed(IOut,IPrint,IPrt1,Dored,NVar,NTBond,NTAng,NTDih,IZred,
     $  Value,NAtoms,C,DQ(1,NPStar),Cnew,BMat,IBMat,VecInt,FCTemp,
     $  IDimSc,Scr,XNew,IFlags,IAn,AtmChg,TRot,TrRMat,X,Q0(1,NPStar),
     $  CnvFX,ConvF,CnvFMX,GradPt(1,NPStar),UseSyT,.False.,
     $  CP(1,1,1,NPstar),FMax,FRMS,DxMaxC,DxRMS,DoneTS,IOSLE1,V,MDV)
      Call AMove(NVar,Q0(1,NPStar),Q(1,NPStar))
      Call AMove(Nvar,Q(1,NPStar),XP(1,NPStar))
      If(.not.DoIrc) then
C       get Tangent with new TS coordinates and old neighbors' coordinates
        If(NPStar.gt.1.and.NPStar.lt.NPath) then
          Call ConsST(IOut,2,0,IStep,NVar,IFlags,VecInP(1,NPstar),
     $      XP(1,NPStar),XP(1,NPstar+1),XP(1,NPStar-1),F,FC,LJunk,
     $      EigTmp,Scr2,NTBond,NTAng,NTDih)
        else if(NPStar.eq.1) then
          Call ConsST(IOut,2,0,IStep,NVar,IFlags,VecInP(1,NPstar),
     $      XP(1,NPStar),XP(1,NPStar+1),Cnstr1,F,FC,LJunk,EigTmp,Scr2,
     $      NTBond,NTAng,NTDih)
        else if(NPStar.eq.NPath) then
          Call ConsST(IOut,2,0,IStep,NVar,IFlags,VecInP(1,NPstar),
     $      XP(1,NPStar),Cnstr2,XP(1,NPstar-1),F,FC,LJunk,EigTmp,Scr2,
     $      NTBond,NTAng,NTDih)
          endIf
        endIf
C     Save a copy of Tangent in XNew
      Call AMove(NVar,Scr2(1,4),XNew)
      Call AMove(NVar,Scr2(1,4),Ftemp)
      Call AMove(NVar,Scr2(1,4),TanP(1,NPStar))
      Call IClear(NPath,IConv)
      Write(IOut,4000) IStep
      FixEig = .false.
      Scleig = .false.
      If(NPStar.gt.1) then
        NStop = 1
        If(DoReac) NStop = 2
C       Check the path approximation using the new TS coordinates
C       and new tangent.  Do reactant valley first.
        Do 150 NPt = (NPStar-1), NStop, -1
          Write(IOut,4010) NPt
          Write(IOut,4020)
C         Prepare IRC Step
          If(NPt.eq.(NPStar-1)) Call AMove(NVar,Xnew,XQuad)
          Call AClear(5*NVar,SCr2)
          Call AMove(NVtt,VecInP(1,NPt),VecInt)
          Call AMove(NVtt,VecInP(1,NPt),A)
C         A contains the eigenvector matrix
C         Scr2(1,1) contains the eigenvalues
          If(NPt.gt.NStop) then
            DoX = .false.
            PathL(Npt) = One
C           DiagF3 checks the force aginst the tangent, get the tangent
C           for the next point downhill, check for stationary point etc...
            Call DiagF3(IOut,IPrint,FixEig,SclEig,FCart,NVar,EigMin,
     $        EigMax,IFlags,FCP(1,NPt),FCTemp,XName,A,Scr2(1,1),VecInt,
     $        NTrRot,Scr2(1,2),GradPt(1,NPt),Scr2(1,2),0,NNegE,NEig,
     $        NAtoms,C,AtmChg,2,DoRed,XP(1,NPt),
     $        XP(1,NPt-1),XP(1,NPt+1),NTBond,NTAng,NTDih,SS2,XQuad,
     $        DQ(1,NPt),DoX,VecInP(1,NPt),USeSym,ArcL,NFroz,Frozen)
            endIf
          If(NPt.eq.NStop) then
            DoX = .false.
            PathL(Npt) = One
            If(DoBirx.and..not.DoReac) then
              DoX = .true.
              PathL(Npt) = -One
              endIf
            Call DiagF3(IOut,IPrint,FixEig,SclEig,FCart,NVar,EigMin,
     $        EigMax,IFlags,FCP(1,NPt),FCTemp,XName,A,Scr2(1,1),VecInt,
     $        NTrRot,Scr2(1,2),GradPt(1,NPt),Scr2(1,2),0,NNegE,NEig,
     $        NAtoms,C,AtmChg,2,DoRed,XP(1,NPt),Cnstr1,
     $        XP(1,NPt+1),NTBond,NTAng,NTDih,SS2,XQuad,DQ(1,NPt),DoX,
     $        VecInP(1,NPt),USeSym,ArcL,NFroz,Frozen)
            endIf
C
C         Save eigenvalues and eigenvectors
C
          Call AMove(NVar*NVar,A,AP(1,NPt))
          Call AMove(NVar,Scr2(1,1),EigP(1,NPt))
          PathL(Npt) = PathL(NPt)*ArcL
          IConv(NPt) = 1
C         Put the Tangent in XQuad for next point IRC Step
          Call AMove(NVar,VecInt,XQuad)
          Call AMove(NVar,VecInt,TanP(1,Npt))
  150     Continue
        endIf
      If(NPStar.lt.NPath) then
C       Now do product valley
        NStop = NPath
        If(DoProd) NStop = NPath-1
        Do 160 NPt = (NPStar+1), NStop
          Write(IOut,4010) NPt
          Write(IOut,4020)
C         Prepare IRC Step
          If(NPt.eq.(NPStar+1)) then
            Call AMove(NVar,Xnew,XQuad)
C           Reverse transition vector
            Call AScale(NVar,-One,XQUad,XQuad)
            If(DoIrc) then
C             Make sure the vector points the rigth way
              Dot = SProd(Nvar,XQuad,XP(1,2))-SProd(Nvar,XQuad,XP(1,1))
              If(Dot.lt.Zero) Call AScale(NVar,-One,XQUad,XQuad)
              endIf
            endIf
          Call AClear(5*NVar,SCr2)
          Call AMove(NVtt,VecInP(1,NPt),VecInt)
          Call AMove(NVtt,VecInP(1,NPt),A)
          If(NPt.lt.NStop) then
            DoX = .false.
            PathL(NPt) = One
            Call DiagF3(IOut,IPrint,FixEig,SclEig,FCart,NVar,EigMin,
     $        EigMax,IFlags,FCP(1,NPt),FCTemp,XName,A,Scr2(1,1),VecInt,
     $        NTrRot,Scr2(1,2),GradPt(1,NPt),Scr2(1,2),0,NNegE,NEig,
     $        NAtoms,C,AtmChg,2,DoRed,XP(1,NPt),
     $        XP(1,NPt+1),XP(1,NPt-1),NTBond,NTAng,NTDih,SS2,XQuad,
     $        DQ(1,NPt),DoX,VecInP(1,NPt),USeSym,ArcL,NFroz,Frozen)
            endIf
          If(NPt.eq.NStop) then
            DoX = .false.
            PathL(Npt) = One
            If(DoBirx.and..not.DoProd) then
              DoX = .true.
              PathL(Npt) = -One
              endIf
            Call DiagF3(IOut,IPrint,FixEig,SclEig,FCart,NVar,EigMin,
     $        EigMax,IFlags,FCP(1,NPt),FCTemp,XName,A,Scr2(1,1),VecInt,
     $        NTrRot,Scr2(1,2),GradPt(1,NPt),Scr2(1,2),0,NNegE,NEig,
     $        NAtoms,C,AtmChg,2,DoRed,XP(1,NPt),Cnstr2,
     $        XP(1,NPt-1),NTBond,NTAng,NTDih,SS2,XQuad,DQ(1,NPt),DoX,
     $        VecInP(1,NPt),USeSym,ArcL,NFroz,Frozen)
            endIf
          Call AMove(NVar*NVar,A,AP(1,NPt))
          Call AMove(NVar,Scr2(1,1),EigP(1,NPt))
          IConv(NPt) = 1
          PathL(Npt) = PathL(NPt)*ArcL
C         Put the Tangent in XQuad for next point IRC Step
          Call AMove(NVar,VecInt,XQuad)
          Call AMove(NVar,VecInt,TanP(1,Npt))
  160     Continue
        endIf
      If(NPStar.gt.1) then
        NStop = 1
        If(DoReac) NStop = 2
        N1 = -1
C
C       Check convergence of path.  If points down-hill from TS are
C       converged there is no need to include them in the relaxation.
C
        NStart = NPStar
        NQuad = NPStar
        Do 170 II = (NPStar+N1), NStop, N1
          Call GetFPA(Nvar,GradPt(1,II),TanP(1,II),XNew,FMax,FRMS,Rjunk,
     $      RJunk2)
          Call RMSVEC(NVar,DQ(1,II),XRMS,XMax)
          If((XRMS.lt.ConvX.and.XMax.lt.T2.and.FMax.lt.T1.and.
     $      FRMS.lt.ConvF).and.(NQuad.eq.II-N1))then
            NQuad = II
          else
            Goto 171
            endIf
  170     Continue
  171   If(NQuad.eq.NStop.and.DoneTS) then
          Write(IOut,3000)
        else
C         Note: There is no need to test for convergence on the
C         reactant. The stability of the path approximation was measured
C         using the new reactant geometry in Diagf3. The path is much
C         less sensitive to reactant than to TS. If the path is stable,
C         then let it be eventhough the reactant geometry changed a bit.
C         This saves a lot of gradient calculations.
C
C         Feed transition vector and relax
C
          Call AMove(NVar,Ftemp,XQuad)
          Rjunk = Zero
          Do 169 II = (NPStar+N1), NStop, N1
            If(SS(II).gt.Rjunk) Rjunk = SS(II)
  169       Continue
C         DxPath relaxes the path approximation in a given valley
          If(DoReac) then
            Call ASub(Nvar,XP(1,2),XP(1,1),Scr2(1,1))
            If(DoRed) Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,
     $        Scr2,XXA)
            Call MultLV(.True.,1,0,NVar,NVar,VecInP(1,1),Scr2,Scr2(1,2))
            Call GetTan(IOut,0,Nvar,.False.,Scr2(1,2),TanP(1,2),Scr2,
     $        Scr2(1,3),SJunk,VecInP(1,1),Scr2(1,4),ArcL)
            PathL(1) = ArcL
            endIf
          Call DxPath(IOut,IPrint,NVar,NPath,NStart,NStop,N1,NVtt,
     $      NTBond,NTAng,NTDih,Cnstr1,IConv,XP,Q0,Q,DQ,VecInP,GradPt,
     $      FCP,AP,EigP,Scr2,SCr2(1,2),IFlags,RJunk,T1,T2,ConvF,ConvX,
     $      Four,Pt5,PathL,XQuad,XNew,USeSym)
          endIf
        endIf
      If(NPStar.lt.NPath) then
        NStop = NPath
        If(DoProd) NStop = NPath - 1
        N1 = 1
        NStart = NPStar
        NQuad = NPStar
        Do 180 II = (NPStar+N1), NStop, N1
          Call GetFPA(Nvar,GradPt(1,II),TanP(1,II),XNew,FMax,FRMS,Rjunk,
     $      RJunk2)
          Call RMSVEC(NVar,DQ(1,II),XRMS,XMax)
          If((XRMS.lt.ConvX.and.XMax.lt.T2.and.FMax.lt.T1.and.
     $      FRMS.lt.ConvF).and.(NQuad.eq.II-N1)) then
            NQuad = II
          else
            Goto 181
            endIf
  180     Continue
  181   If(NQuad.eq.NStop.and.DoneTS) then
          Write(IOut,3000)
        else
          Call AMove(NVar,Ftemp,XQuad)
          Call AScale(NVar,-One,XQUad,XQuad)
          If(DoIrc) then
C           Make sure the vector points the rigth way
            Dot = SProd(Nvar,XQuad,XP(1,2))-SProd(Nvar,XQuad,XP(1,1))
            If(Dot.lt.Zero) Call AScale(NVar,-One,XQUad,XQuad)
            endIf
          Rjunk = Zero
          Do 179 II = (NPStar+N1), NStop, N1
            If(SS(II).gt.Rjunk) Rjunk = SS(II)
  179       Continue
          If(DoProd) then
            Call ASub(Nvar,XP(1,NPath-1),XP(1,NPath),Scr2(1,1))
            If(DoRed) Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,
     $        Scr2,XXA)
            Call MultLV(.True.,1,0,NVar,NVar,VecInP(1,NPath),Scr2,
     $        Scr2(1,2))
            Call GetTan(IOut,0,Nvar,.False.,Scr2(1,2),TanP(1,NPath-1),
     $        Scr2,Scr2(1,3),SJunk,VecInP(1,NPath),Scr2(1,4),ArcL)
            PathL(Npath) = ArcL
            endIf
          Call DxPath(IOut,IPrint,NVar,NPath,NStart,NStop,N1,NVtt,
     $      NTBond,NTAng,NTDih,Cnstr2,IConv,XP,Q0,Q,DQ,VecInP,GradPt,
     $      FCP,AP,EigP,Scr2,SCr2(1,2),IFlags,RJunk,T1,T2,ConvF,ConvX,
     $      Four,Pt5,PathL,XQuad,Xnew,UseSym)
          endIf
        endIf
      Steped=.true.
      endIf
      endIf
      If(NPath.gt.0) then
        IFinish = 0
        If(NConv.eq.NPath) then
C         Get new Path
          Do 480 NPt = 1, NPath
            SJunk = Zero
            If(NPt.lt.NPstar) then
              Do 475 NN = (NPStar-1), NPt, -1
  475           SJunk = SJunk - Abs(PathL(NN))
              endIf
            If(NPt.gt.NPstar) then
              Do 476 NN = (NPStar+1), NPt
  476           SJunk = SJunk + Abs(PathL(NN))
              endIf
            ArcP(Npt) = Sjunk
            Call AClear(NVar,XNew)
            Call AClear(NVar,XQuad)
            Call AMove(NVar,XXP(1,1,NPt),X)
            Call AMove(NVar,X,Value)
            Call AMove(NVar,FFP(1,1,NPt),F)
            Call ASub(NVar,XP(1,NPt),X,XNew)
            If(DoRed) Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,
     $        XNew,XXA)
            Call AMove(NAt3,CP(1,1,1,NPt),C)
            Call AMove(9,TRotP(1,1,NPt),TRot)
            DXMax = Zero
            DXRMS = Zero
            Call RMSVEC(NVAr,F,FRMS,FMax)
            If(Npt.ne.NPstar.and..not.(Npt.eq.1.and.Doreac).and..not.
     $        (Npt.eq.NPath.and.DoProd)) then
              Call GetFPA(Nvar,F,TanP(1,NPt),XQuad,FMax,FRMS,Rjunk,
     $          RJunk2)
              Call AAdd(Nvar,XQuad,F,XQuad)
              Call Asub(Nvar,XNew,Xquad,Xnew)
              endIf
            Write(IOut,4012)NPt,SJunk
            Write(IOut,4022)
            Do 490 I = 1, NVar
              APrint(1,I) = X(I)
              APrint(2,I) = F(I)
              APrint(4,I) = XNew(I)
              APrint(3,I) = XQuad(I)
              XNew(I) = XNew(I) + XQuad(I)
  490         X(I) = X(I) + XNew(I)
            IPrt1 = IPrint
            If(IPSave(0).le.1) IPrt1 = Max(IPrint,1)
            IGISav = 0
            If(DoRed) then
              DidGI = IGISav.ne.0
              If(DidGI) Call FileIO(2,-IGISav,NVar**2,A,0)
              Call IClear(Nvar,IFlags)
              If(UpdRWF.and.Do601.and..not.(Npt.eq.1.and.Doreac).and.
     $          .not.(Npt.eq.NPath.and.DoProd)) then
C               transform tangent before saving to disk
                Call RedCar(IOut,IPrt1,0,NVar,NTBond,NTAng,NTDih,IZRed,
     $            Value,NAtoms,C,TanP(1,NPt),CNew,BMat,IBMat,VecInt,
     $            FCTemp,Scr,XQuad,IFlags,.True.,IFail,0,V,MDV)
                LOpatN = (2+2+3*Nat3+(Nat3*(Nat3+1))/2)
                IPoint = LOpatN*(NCurPt-1)+4
                Call Fileio(1,-IOPath,1,ArcP(Npt),IPoint)
                IPoint = 2*NAtoms*3+1
                Call Trspn2(3,TRot)
                Call RotF1(NAtoms,TRot,CNew)
                Call Fileio(1,IOPath,NAtoms*3,CNew,IPoint)
                Call Trspn2(3,TRot)
                endIf
              Call RedCar(IOut,IPrt1,0,NVar,NTBond,NTAng,NTDih,IZRed,
     $          Value,NAtoms,C,XNew,CNew,BMat,IBMat,VecInt,FCTemp,Scr,
     $          XQuad,IFlags,.True.,IFail,0,V,MDV)
            else if(FCart) then
              Call AMove(NAt3,X,CNew)
              endIf
            If(DoRed.or.FCart) then
              If(USeSym.and.PthSym) Call FixSym(.False.,NAtoms,
     $          NOpAlP(NPt),SymAlP(1,1,1,NPt),IPrAlP(1,1,NPt))
              UseSyT = UseSym.and.PthSym.and.NOpAlP(NPt).gt.1
              Call ClnCor(IOut,IPrint,.False.,UseSyT,NAtoms,IAn,AtmChg,
     $          C,0,Scr,.True.,TRot,CNew,Scr)
              Call MatCor(IOut,IPrt1,.True.,NAtoms,C,CNew,TrRMat,Scr,OK)
              If(.not.OK) Write(IOut,2100)
              Call ASub(NAt3,C,Scr,C)
              DXMaxC = ArrMax(C,NAt3,.True.)
              DXRMSC = Sqrt(Abs(SProd(NAt3,C,C)/GFloat(NVar)))
              Call AMove(NAt3,Scr,C)
              Call AMove(NAt3,CP(1,1,1,Npt),CP(1,1,2,NPt))
              Call AMove(Nat3,C,CP(1,1,1,NPt))
              endIf
            If(DoRed) then
              Call RedVal(IOut,IPrint,.True.,C,IZRed,X,NVar,NTBond,
     $          NTAng,NTDih,V,MDV)
            else if(FCart) then
              Call AMove(NAt3,C,X)
              endIf
            Do 495 I = 1, NVar
  495         XNew(I) = X(I) - APrint(1,I)
            If(DoRed) Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,
     $        XNew,XXA)
            Do 497 I = 1, NVar
              X(I) = APrint(1,I) + XNew(I)
              APrint(5,I) = XNew(I)
              APrint(6,I) = X(I)
              DXRMS = DXRMS + XNew(I)**2
              If(Abs(XNew(I)).gt.DXMax) DXMax = Abs(XNew(I))
  497         Continue
            DXRMS = Sqrt(DXRMS/NVar)
C
C           For points along the path, get the less stringent
C           convergence criteria
C
            If(Npt.ne.NPstar.and..not.(Npt.eq.1.and.DoReac).and..not.
     $        (Npt.eq.NPath.and.DoProd)) then
              Write(IOut,1052)
              If(DXRMS.lt.DxRMSC) DxRMSC = DXRMS
              If(DXmax.lt.DxMaxC) DxMaxC = DxMax
            else
              Write(IOut,1050)
              endIf
            NCur = 0
            Do 515 I = 1, NVar
              Call GetB(2,ITmp,Len,XName,NCur)
              Call CBlank(8,IStr)
              ICur = 0
              Len = Min(Len,8)
              JPos = (10-Len) / 2 - 1
              Do 510 J = 1, Len
  510           IStr(JPos+J) = GetChr(ITmp,ICur)
              Write(IOut,1060) (IStr(K),K=1,8),(APrint(J,I),J=1,6)
  515         Continue
            T1 = ConvF * CnvFMx
            T2 = ConvX * CnvFMx
            Write(IOut,1000)
            Call Convgd(FMax,T1,IReslt)
            Write(IOut,1002) FMax,T1,IReslt
            Call Convgd(FRMS,ConvF,IReslt)
            Write(IOut,1004) FRMS,ConvF,IReslt
            Call Convgd(DXMaxC,T2,IReslt)
            Write(IOut,1006) DXMaxC,T2,IReslt
            Call Convgd(DXRMSC,ConvX,IReslt)
            Write(IOut,1008) DXRMSC,ConvX,IReslt
            Iconv(NPt) = 0
            If(DXRMSC.lt.ConvX.and.DXMaxC.lt.T2) then
              IFInish = IFinish + 1
              IConv(NPt) = 1
              If(NPt.eq.NPStar) then
*
* ANT BEGIN modification
*
*        Avoiding to stop on forces
*
                If(FRMS.gt.ConvF.or.FMax.gt.T1) then
                 !Ifinish = Ifinish - 1
                 !IConv(NPt) = 0
                  endIf
                endIf
              If(NPt.eq.1.and.DoReac) then
                If(FRMS.gt.ConvF.or.FMax.gt.T1) then
                 !Ifinish = Ifinish - 1
                 !IConv(NPt) = 0
                  endIf
                endIf
              If(NPt.eq.NPath.and.DoProd) then
                If(FRMS.gt.ConvF.or.FMax.gt.T1) then
                 !Ifinish = Ifinish - 1
                 !IConv(NPt) = 0
*
* ANT END modification
*
                  IConv(NPt) = 0
                  endIf
                endIf
              If(UpdRWF.and.Do601) then
                If(IConv(Npt).gt.0) then
                  LOpatN = (2+2+3*Nat3+(Nat3*(Nat3+1))/2)
                  IPoint = 2+(Npt-1)*LOPatN+1
                  Call Fileio(2,-IOPath,1,WJunk,IPOint)
                  NJunk = IGFix(WJUnk)
                  LOPatN = LOPAtN - 2
                  If((Npt.eq.1.and.DoReac).or.(NPt.eq.NPAth.and.DoProd))
     $              LOPatN = LOPAtN - NAt3
                  If(NJunk.ne.4) LOPatN = LOPatN-(NAt3*(NAt3+1))/2
                  WOPatN = GFloat(LOPAtN)
                  Call Fileio(1,-IOPath,1,WOPatN,IPOint)
                  If((Npt.eq.1.and.DoReac).or.(NPt.eq.NPAth.and.DoProd)
     $              .and.NJunk.eq.4) then
                    IPoint = IPoint+2+3*NAt3
                    Call Fileio(2,-IOPath,NAtTT,A,IPoint)
                    IPoint = IPoint-NAt3
                    Call Fileio(1,-IOPath,NAtTT,A,IPoint)
                    endIf
                  endIf
                endIf
              endIf
            If(Npt.eq.NPath) Write(IOut,4050) IFinish
C
C           Replace XP by its back-transformed components.
C
            Call FixDih(NTDih,X(NTBOnd+NTAng+1))
            Call AMove(NVar,X,XP(1,NPt))
C
C           If Point converged do not update coordinates.
C
            If(IConv(NPt).eq.1) then
              Do 470 I = 1, NVar
  470           XP(I,NPt) = APrint(1,I)
              Call AMove(Nvar,XXP(1,1,Npt),X)
              NP2 = NP
              NP = IRisP(NPt)
              IStep2 = Istep
              IStep = NP
              Esave(1) = EsvP(Npt)
              Energy = Esave(1)
              Call AMove(MaxSV1,FSP(1,NPt),FS)
C             gradient won't be calculated, so fill the arrays with what
C             we already got.
              Call SavePt(Energy,NVar,IGetFC,.False.,0,Rises,
     $          GradPt(1,Npt),FCP(1,NPt),ESave,ICP(1,Npt),ITU,X,
     $          FFP(1,1,Npt),XXP(1,1,Npt),FCP(1,Npt),FC2P(1,Npt),
     $          .False.,LJunk,1,1,1,RR,RR,RR,RR,RR,RR,RR(1),RR)
              Call AMove(MaxSv1,FS,FSP(1,NPt))
              NP = NP2
              IStep = IStep2
              Call AMove(NAt3,CP(1,1,2,NPt),CP(1,1,1,NPt))
              endIf
  480       Continue
          endIf
C
C       Choose next Point for gradient calculation
C
      If(IFinish.ne.NPath) then
        Call ChoNxt(IOut,NPath,NIStep,NCurPt,IFinish,IConv,NNext)
      else
        NNExt = NPStar
        endIf
      If(UseSym.and.PthSym) Call FixSym(.False.,NAtoms,
     $  NOpAlP(NNext),SymAlP(1,1,1,NNext),IPrAlP(1,1,NNext))
      NPtTmp = -NNext
      Call OptDat(XP,XXP,GradPt,FFP,FCP,FC2P,ESvP,FSP,VecInP,ICP,TRotP,
     $  CP,Value,XX,F,FF,FC,FC2,ESave,FS,VecInt,IC,TRot,C,NAtoms,NVar,
     $  MaxSv1,NVTT,NPath,NPtTmp)
      Write(IOut,5000)NNext
      Call AMove(NVar,XP(1,NNext),X)
      Call AMove(NVar,FFP(1,1,NNext),F)
      Finish=.false.
C     Summarize performance
      If(IFinish.eq.NPath) then
        Write(IOut,5999)
        Write(IOut,5020)
        Write(IOut,5040)ICP(1,NPStar)
        If(DoReac.and.NPStar.ne.1) Write(IOut,5050) ICP(1,1)
        If(DoProd) Write(IOut,5060) ICP(1,NPath)
        NST = 0
        NRem = 0
        Do 701 II = 1, NPath
          If(DoReac.and.II.eq.1) Goto 701
          If(DoProd.and.II.eq.NPath) Goto 701
          If(II.eq.NPStar) Goto 701
          NST = NST + ICP(1,II)
          NRem = NRem + 1
  701     Continue
        StAvg = GFloat(NST)/GFloat(Nrem)
        Write(IOut,5070) StAvg
        If(Doreac.and.NPStar.ne.1) NST = NST + ICP(1,1)
        If(DoProd) NST = NST + ICP(1,NPath)
        NST = NST + ICP(1,NPStar)
        Write(IOut,5080) NST
        Do 702 I = 1, NVar
  702     IC(I) = 99
        Write(IOut,5999)
        Call RedTbl(IOut,1,XName,XP(1,NPStar),Cnstr1,Cnstr2,IC,
     $    GradPt(1,NPStar),NVar,IZRed,ToAng,NTBond,NTAng,.True.,Zero)
        Finish = .true.
      else if(IStep.ge.NStep) then
        Write(IOut,5030)
        Finish = .true.
        Failed = .true.
        endIf
      Exit = Finish.or.Failed
      If(Failed) then
        Call ILSW(1,25,0)
        Call ILSW(1,27,1)
        Call RedTbl(IOut,2,XName,XP(1,NPStar),Cnstr1,Cnstr2,IC,
     $    GradPt(1,NPStar),NVar,IZRed,ToAng,NTBond,NTAng,.True.,Zero)
        endIf
      If(NConv.lt.NPath) then
        IStep = IStep - 1
        NP = NP - 1
        If(IStep.eq.0) NP = -1
        endIf
      Endif
      Call FreeRW(IGISav)
      Return
      End
*Deck PthSpd
      Subroutine pthspd(IOut,IPrint,NAtoms,CP,CR,CX,XValue,IZRed,NTRed,
     $  NTBond,NTAng,NTDih,Rstr,PStr,XStr,B,IB,Square,Tri,One4D,CnStr,
     $  IFlags,MaxSv1,NPath,Scr,CPP,DoReac,DoProd,NStart,NStop,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Get an Npath-point LST path. This routine is basically the same as
C     RedTS but does the interpolation in internals rather than Cartesians
C     to deal better with group migrations and internal rotations.
C
      Logical DoReac,DoProd
      Dimension B(3,4,NTRed),IB(4,NTRed),Square(NTRed,NTRed),Tri(*),
     $  One4d(NTRed,4),CR(3,NAtoms),CP(3,NAtoms),CX(3,NAtoms),
     $  IZRed(4,NTRed),RStr(NTRed),PStr(NTRed),XStr(NTRed),
     $  XValue(NTRed,MaxSv1,NPath),Cnstr(NTRed),IFlags(NTRed),
     $  Scr(NTred),CPP(3,NAtoms,2,NPath),V(*)
C
      NAt3=3*NAtoms
      If(DoProd)then
      Call AMove(NTred,PStr,XValue(1,1,NStop))
      Call Amove(NAt3,CP,CPP(1,1,1,NStop))
      Endif
      If(DoReac)then
      Call AMove(NTred,RStr,XValue(1,1,NStart))
      Call Amove(NAt3,CR,CPP(1,1,1,NStart))
      Endif
      Call Amove(NTred,PStr,Scr)
      N1=1
      If(NStop.lt.NStart)N1=-1
      III=N1*(NStop-NStart)+1
      If(Doreac)Then
      NStart=NStart+N1
      III=III-1
      Endif
      If(DoProd)then
      NStop=NStop-N1
      III=III-1
      Endif
      II=0
      Do 100 Npt=NStart,NStop,N1
      II=II+1
C
C     Go one sgement at a time
      WeiFcN=GFloat(III-II+1)/GFloat(III-II+2)
      Call GTrPth(IOut,IPrint,NAtoms,CR,CX,XValue(1,1,NPt),IZRed,NTRed,
     $  NTBond,NTAng,NTDih,WeiFcN,2,RStr,PStr,XStr,B,IB,Square,Tri,
     $  One4D,CnStr,IFlags,V,MDV)
       Call RedVal(IOut,IPrint,.True.,CX,IZRed,XStr,NTred,NTBond,NTAng,
     $   NTDih,V,MDV)
C
C     Make the new point the old one
      Call AMove(NTred,XStr,XValue(1,1,NPt))
      Call AMove(NTred,XValue(1,1,NPt),RStr)
      Call AMove(NAt3,CX,CR)
C
C     PStr was destroyed , so reload
      Call AMove(NTred,Scr,PStr)
      Call AMove(NAt3,CX,CPP(1,1,1,NPt))
 100  Continue
      Return
      End
*Deck PthTrs
      Subroutine PthTrs(IOut,IPrint,NVar,NTBond,NTAng,NTDih,NVtt,NPath,
     $  NPStar,NStart,NStop,N1,MaxSv1,IStep,XNew,XXP,FFP,XQuad,FCP,ESvP,
     $  FSP,DXMxT0,ESave,IRisP,SS,DoProd,ICount,T2,ConVX)
      Implicit Real*8(A-H,O-Z)
C
C     For Path Relaxation.
C     Update the trust radius for the points along the path.
C
      Logical Doprod
      Real*8 MDCutO
      Dimension XNew(NVar),XXP(Nvar,MaxSv1,NPath),
     $  FFP(Nvar,MaxSv1,NPath),XQuad(Nvar),FCP(NVtt,Npath),
     $  ESvP(NPath),FSP(MaxSv1,NPath),ESave(*),
     $  IRisP(Npath),SS(NPath),XX(1)
      Save Two, Scale, StpMin, RMin, RGood, Zero, StpMax, One, XX
      Data Two, Scale, StpMin, RMin, RGood, Zero, StpMax, One, XX
     $  /2.0d0, 2.0d0, 0.05d0, 0.7d0, 0.85d0, 0.0d0, 0.6d0, 1.0d0, 0.d0/
  900 Format(' Update trust radius for point ',I3)
 1000 Format(' Trust test=',1PD9.2,' RLast=',1PD9.2,' DXMaxT set to',
     $  1PD9.2)
 1010 Format(' No trust radius update -- first pass')
 1020 Format(' Skip trust radius update -- point converged')
C
      Do 100 NPt = NStart, NStop, N1
        Write(IOut,900) NPt
        If(IRisP(Npt).lt.2) then
          Write(IOut,1010)
          Goto 100
          endIf
C       Test if product or ts has already converged or first step
        If(DoProd.and.Npt.eq.NStop.and.IRisP(Npt).lt.IStep) then
          SS(Npt) = Zero
          Write(IOut,1020)
        else if(Npt.eq.NPStar.and.IRisP(Npt).lt.IStep) then
          SS(Npt) = Zero
          Write(IOut,1020)
        else if(IRisP(Npt).ge.2) then
C         Estimate the energy change.  Get Delta X
          Call AClear(NVar,XNew)
          Call Asub(NVar,XXP(1,2,NPt),XXP(1,1,NPt),XNew)
          Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,XNew,XX)
          Call RMSVEC(Nvar,XNew,DXRMS,DXMax)
          DA = SProd(Nvar,FFP(1,2,NPt),XNew)
          Call MultLV(.True.,1,0,NVar,NVar,FCP(1,NPt),XNew,XQuad)
C         Estimate gradient change
          DePred = SProd(NVar,XQuad,XQuad)
          DE = SProd(Nvar,FFP(1,1,NPt),XQuad) -
     $      Sprod(Nvar,FFP(1,2,Npt),XQuad)
          DA = DA + SProd(Nvar,XQuad,Xnew)/Two
          Call Asub(NVar,FFP(1,1,NPt),FFP(1,2,NPt),XQuad)
          Dot1 = SProd(Nvar,XQuad,XNew)
          Call RMSVEC(Nvar,XQuad,DFRMS,DFMax)
          If(DXRMS.lt.ConvX.and.DXMax.lt.T2) then
            SS(Npt) = Zero
            Write(IOut,1020)
          else
            ESave(1) = EsvP(Npt)
            ESave(2) = FSP(1,NPt)
            Call AMove(NVar,XNew,XQuad)
            Call AAdd(NVar,XNew,XXP(1,1,NPt),XNew)
            Neg = 1
            If(Npt.eq.NStop.and.Doprod) Neg = 0
            If(IRisP(NPt).ge.2) then
              DXPT = SS(Npt)
              If(DXPt.lt.StpMin) DXPt = DxMxT0
              If(Npt.eq.NPStar.or.(npt.eq.NStop.and.Doprod)) then
                S = One
C               Test change in curvature; be conservative
                If(Neg.eq.0.and.IRisP(Npt).ge.3) then
                  Call ASub(Nvar,XXP(1,3,NPt),XXP(1,2,Npt),Xquad)
                  Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,
     $              XQuad,XX)
                  Dot2 = SProd(NVar,XQuad,FFP(1,2,NPt))
     $              - SProd(Nvar,Xquad,FFP(1,3,NPt))
                  S = Dot1*Dot2
                  endIf
C               Do update for TS and minima
                Call UpdDXM(IOut,1,.True.,Neg,DA,ESave,DXPT,NVar,
     $            XXP(1,1,NPt),XNew,S,IUpDwn)
              else
                Small = MDCutO(0)
                Rlast = Zero
                Do 10 I = 1, NVar
   10             RLast = RLast + (XNew(I)-XXP(I,1,Npt))**2
                RLast = Sqrt(RLast)
                R = Zero
                If(Abs(DEPred).ge.Small) then
                  R = DE / DEPred
                  DAA = ESave(1) - ESave(2)
C                 take the average of the two diagnostics to stay on
C                 the safe side.
                  R = R + DAA/DA
                  R = R / Two
                  If(R.le.RMin.or.R.ge.(Two-RMin)) DXPT = DXPT / Scale
                  If(R.ge.RGood.and.R.le.(Two-RGood))
     $              DXPT = DXPT*Sqrt(Scale)
                  endIf
C               Ensure that the steps aren't absurdly small or large.
                DXPT = Max(DXPT,StpMin)
                DXPT = Min(DXPT,StpMax)
                Write(IOut,1000) R, RLast, DXPT
                endIf
              ICount = ICount + 1
              SS(Npt) = DXPT
              endIf
            endIf
          endIf
  100   Continue
      Return
      End
*Deck PtPath
      Subroutine PtPath(IOut,NVar,ToAng,NPath,XP,EsvP,XName,NTBond,Exit)
      Implicit Real*8(A-H,O-Z)
C
C     Print out a summary of the results of a relaxed reaction path
C     scan.  Note that Res should be dimensioned (NVar+1)*NScan.
C
      Dimension   XName(*), XP(NVar,NPath),ESvP(NPath)
      Logical   ConvBd,EXit
      Save One
      Data One/1.0d0/
 1000 Format(12x,'  Converged reaction path structures ')
 1010 Format(12x,' Interpolated reaction path structures')
 1020 Format(26X,'  Angstrom and degree ')
C
      ToDeg = GFloat(45) / ATan(GFloat(1))
      Do 10 I = 1, NVar
          ConvBd = I.le.NTBond
        If(ConvBd) then
          Conver = ToAng
        else
          Conver = ToDeg
          endIf
        Do 10 J = 1, NPath
   10     XP(I,J) = Conver*XP(I,J)
      If(Exit)then
      Write(IOut,1000)
      Write(IOut,1020)
        Call MatPt1(IOut,XP,NVar,NPath,NVar,NPath,1,0,XName,XName,
     $  0,ESvP,1)
      Else
      Write(IOut,1010)
      Write(IOut,1020)
        Call MatPt1(IOut,XP,NVar,NPath,NVar,NPath,1,0,XName,XName,
     $  0,ESvP,0)
      Endif
      Do 20 I = 1, NVar
          ConvBd = I.le.NTBond
        If(ConvBd) then
          Conver = One/ToAng
        else
          Conver = One/ToDeg
          endIf
        Do 20 J = 1, NPath
   20     XP(I,J) = Conver*XP(I,J)
      Return
      End
*Deck PtScan
      Subroutine PtScan(IOut,IOScan,NVar,ToAng,MaxNZ,NZ,LBl,LAlpha,
     $  LBeta,IZ,XName,Res,DoRed,NTBond,JScan)
      Implicit Real*8(A-H,O-Z)
C
C     Print out a summary of the results of a relaxed potential surface
C     scan.  Note that Res should be dimensioned (NVar+1)*NScan.
C
      Dimension LBl(MaxNZ), LAlpha(MaxNZ), LBeta(MaxNZ), IZ(MaxNZ,4),
     $  XName(*), Res(NVar,*), JScan(4)
      Logical VDist, DoRed, ConvBd
 1000 Format(' Summary of Optimized Potential Surface Scan')
C
      ToDeg = GFloat(45) / ATan(GFloat(1))
      Call FileIO(2,-IOScan,InToWP(2),JScan,0)
      IPos = InToWP(2)
      Call FileIO(2,-IOScan,JScan(2),Res(1,JScan(2)+1),IPos)
      IPos = InToWP(2) + JScan(2)
      Call FileIO(2,-IOScan,JScan(2)*NVar,Res,IPos)
      Do 10 I = 1, NVar
        If(DoRed) then
          ConvBd = I.le.NTBond
        else
          ConvBd = VDist(I,MaxNZ,NZ,LBl,LAlpha,LBeta,IZ)
          endIf
        If(ConvBd) then
          Conver = ToAng
        else
          Conver = ToDeg
          endIf
        Do 10 J = 1, JScan(2)
   10     Res(I,J) = Conver*Res(I,J)
      Write(IOut,1000)
      Call MatPt1(IOut,Res,NVar,JScan(2),NVar,JScan(2),1,0,XName,XName,
     $  0,Res(1,JScan(2)+1),1)
      Return
      End
*Deck QIMOMM
      Subroutine QIMOMM(IOut,IPrint,NAtoms,NatNew,C,CNewQ,NVar,NTBond,
     $  NTAng,NTDih,IZRed,X,XNew,Value,IBMat,Scr,Cnstr1,IFlags,IWght,
     $  RotQM,NRF,T,MicOpt,COM,COrig,DT,ICnt,XMMT,NVarMM,IMMMOD,HMod,
     $  HModO,CMM,IRWMM,IAtTyp,AtChMM,NumTpS,IAtTpS,MaxLEF,EField,FTp,
     $  IMMRFO,NTrR1,BMat,IDBMat,DBMat,FTemp,A,ISCRF,Steped,IGuess,IRFO,
     $  FC,XQuad,MethQ,IDir,Neg,IHVAcc,TstEig,XName,IZRedO,FCTemp,
     $  TSOnly,SmlStp,IGetTS,IPFlag,AllowP,FMM,FMFlag,FMFlg1,NFxFlg,
     $  IHMeth,IntHes,DoEle,IRwF2,HReal,IMMCRS,NAtmNN,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Do the fully coupled QM/MM macro step.
C     B-matrices and related quantities will be changed!
C
      Logical RotQM, Steped, StepD, Convg, TstEig, OK, TSOnly, SmlStp,
     $  FMM, AllowP(*), LJJ(1), IntHes, DoEle, MMRest
      Integer FMFlag, FMFlg1
      Dimension V(*), C(3,*), CNewQ(3,*), IZRed(4,*), IZRedO(4,*), X(*),
     $  XNew(*), Value(*), IBMat(4,*), Scr(*), Cnstr1(*), IFlags(*),
     $  JJ(1), XX(1), MicOpt(*), DT(*), ICnt(*), A(*), IAtTpS(*),
     $  COrig(3,*), T(3,3,10,*), COM(3,*), XMMT(3,*), CMM(3,*), HMod(*),
     $  HModO(*), IAtTyp(*), AtChMM(*), EField(0:3), FTp(3,*), BMat(*),
     $  IDBMat(*), DBMat(*), FTemp(*), FC(*), XQuad(*), XName(*),
     $  FCTemp(*), QMVar(30), HReal(*)
      Save Zero, JJ, XX, LJJ
      Data Zero/0.D0/, JJ/0/, XX/0.0d0/, LJJ/.False./
C
C     First get the cartesian coordinates from the internal coordinates and
C     the original coordinates. This needs to be done because the linear
C     search and/or swapping in GrdOpt can have modified only the internal
C     coordinates, while we need cartesians to evaluate the various MM
C     contributions.
C
      Call TStamp(1,'Top QIMOMM')
      III   = Max(3*NAtNew,NVar)
      IVect = 1
      IFCT  = IVect + III**2
      ITTmp = IFCT + III**2
      ITTmp2= ITTmp + III
      IEnd  = ITTmp2 + III - 1
      Call TstCor(IEnd,MDV,'GrdOpt-RedCar')
      Call AMove(3*NAtoms,C,CNewQ)
      Call AMove(NVar,XNew,V(ITTmp2))
      Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,XNew,XX)
      Call AAdd(NVar,X,XNew,V(ITTmp))
      Call FixDih(NTDih,V(ITTmp+NTBond+NTAng))
      Call ASub(NVar,V(ITTmp),X,XNew)
      Call RedCar(IOut,IPrint,0,NVar,NTBond,NTAng,NTDih,IZRed,Value,
     $  NAtNew,C,XNew,CNewQ,BMat,IBMat,V(IVect),V(IFCT),Scr,Cnstr1,
     $  IFlags,.True.,IFail,IWght,V(IEnd+1),MDV-IEnd)
      If(IFail.eq.1) Call AMove(3*NAtoms,C,CNewQ)
      Call AMove(NVar,V(ITTmp2),XNew)
      Call UpCOrg(CNewQ,V,COrig,NAtoms,NAtNew,MicOpt,XMMT,COM,DT,ICnt,
     $  T,NRF,NVar,NVarMM,XNew,0,RotQM,.False.)
C
C     Do the MM calculations with the new coordinates.
C
      IDat = 1
      LenDat = ITqry(IMMMOD)
      IVX = IDat + LenDat
      IEndX = IVX - 1
      Call TstCor(IEndX,MDV,'GrdOpt-EGHMod')
      Call FileIO(2,-IMMMOD,LenDat,V,0)
      MDVX = MDV - LenDat
      NAtomI = ICount(-1,NAtoms,MicOpt)
      If(IntHes) then
        Call AClear(3*NAtoms,FTP)
        NClear = ((NAtomI*3+1)*NAtomI*3)/2
        Call AClear(NClear,HMod)
        Call AClear(NClear,HModO)
        Call AClear(NClear,HReal)
      else
C       Now all has been set up, get the one for the original coordinates
C       (in CMM) to be used later on for the prediction of the energy change.
        Call EGHMod(IOut,IPrint,1,Zero,NAtoms,MicOpt,CMM,ETemp,Scr,
     $    HModO,V(IDat),V(IDat),V(IVX),MDVX)
        Call EGHMod(IOut,IPrint,1,Zero,NAtoms,MicOpt,CNewQ,ETemp,Scr,
     $    HMod,V(IDat),V(IDat),V(IVX),MDVX)
        If(DoEle) then
          IndexQ = 6 + 6*NRF + 12*NAtoms + 3*NVar + 9*(NVarMM+6*NRF)
          Call FileIO(2,-IMMCRS,NAtmNN,V(IDat+InToWP(7)),IndexQ)
          Call EGHMod(IOut,IPrint,3,Zero,NAtoms,MicOpt,CNewQ,ETemp,Ftp,
     $     XX,V(IDat),V(IDat),V(IVX),MDVX)
          Call ASub(3*NAtoms,Scr,Ftp,Scr)
          endif
        IParam = 1
        LenPar = ITqry(IRwMM)
        IDDip  = IParam + LenPar
        IDip   = IDDip + 9*NAtoms
        IV     = IDip + 3
        IEnd1  = IV - 1
        MDV1   = MDV - IEnd1
        Call TstCor(IEnd1,MDV,'GrdOpt-EGHMM')
        Call FileIO(2,-IRwMM,LenPar,V(IParam),0)
        Call EGHMM(IOut,IPrint,0,52,0,Zero,NAtoms,IAtTyp,AtChMM,CNewQ,
     $    MaxLEF,EField,V(IParam),LenPar,MicOpt,ETemp,FTp,HReal,V(IDip),
     $    V(IDDip),V(IV),MDV1)
        Call ASub(3*NAtoms,FTp,Scr,FTp)
C       Zap the elements referring to frozen QM centers.
        Call FrzAtm(.True.,NAtoms,MicOpt,NFroz,V(IV))
        If(NFroz.gt.0) Call ZeroD(NAtoms,V(IV),1,1,Ftp)
        Call FileIO(1,-IMMRFO,3*NAtoms,FTp,0)
        I1 = 1
        I2 = I1 + NAtoms
        Call IMove(NAtoms,MicOpt,V(I1))
        Do 35 I = 1, NAtoms
          If(MicOpt(I).ge.0) MicOpt(I) = 1
   35     continue
        Call FrzAtm(.True.,-NAtoms,MicOpt,NFroz,V(I2))
        Call IMove(NAtoms,V(I1),MicOpt)
        If(NFroz.gt.0) then
          Call ZeroD(NAtomI,V(I2),1,2,HModO)
          Call ZeroD(NAtomI,V(I2),1,2,HMod)
          Call ZeroD(NAtomI,V(I2),1,2,HReal)
          endif
        endif
C
      IQM = 1
      Do 40 I = 1, NAtoms
        If(MicOpt(I).lt.0) then
          Call ANeg(3,FTp(1,I),FTp(1,IQM))
          IQM = IQM + 1
          endif
   40   Continue
C
C     Set up the remaining data for the coupled step.
C
      IPrjct = 0
      MaxCyc = 50
      ConvMM = 8.0d80
      ITrScl = 0
      LenPar = ITqry(IRwMM)
      NAt3 = 3*NAtNew
      NAt3TT = (NAt3*(NAt3+1))/2
      IParam = 1
      IFX    = IParam + LenPar
      IFFX   = IFX + 3*NAtoms + NAt3TT
      LenFFX = 3*NAtoms
      IGM    = IFFX + LenFFX
      IFFtmp = IGM + NVar**2
      NVarT  = NVar + NVarMM + 6*NRF
      ICRig  = IFFtmp + NVarT
      IVecIn = ICRig + NRF*6
      IV     = IVecIn + NVar**2
      MDV1   = MDV - IV
      IQM = IFX
      Do 50 I = 1, NAtoms
        If(MicOpt(I).lt.0) then
          Call AMove(3,CNewQ(1,I),V(IQM))
          IQM = IQM + 3
          endif
   50   continue
      Call FormDB(IOut,IPrint,0,NAtoms,NTBond,NTAng,NVar,IZRed,V(IFX),
     $  Value,IBMat,BMat,IDBMat,DBMat,V(IV),MDV1)
      IGISav = 0
      Call Red2BG(IOut,IPrint,0,NAtNew,NTrR1,IGISav,0,V(IFX),NVar,
     $  NTBond,NTAng,NTDih,IZRed,Value,NRank,BMat,IBMat,V(IGM),
     $  V(IVecIn),FCTemp,Scr,XName,V(IV),MDV1)
      NVTT = (NVar*(NVar+1))/2
      Call AMove(NVTT,FCTemp,V(IVecIn))
      Call CnstrP(IOut,IPrint,NVar,IFlags,-100000,-1,V(IVecIn),A,
     $  FCTemp,Scr)
      Call TranF(.True.,.False.,NVar,NAtoms,JJ,FTp,Scr,IBMat,BMat,
     $  V(IGM),JJ)
C     FTemp will end up in RMult(IGAll) and be used for the RFO step.
C     With electronic embedding (IntHes or regular), FTemp already
C     contains the full gradient.
      If(DoEle.or.IntHes) Call AClear(NVar,Scr)
      Call AAdd(NVar,FTemp,Scr,Scr)
      Call AMove(NVarT,FTemp,V(IFFTmp))
      Call MultLV(.True.,1,0,NVar,NVar,V(IVecIn),Scr,V(IFFTmp))
C     CRig is only used for the gradient contribution to the Hessian. For
C     IntHes, this must correspond to the original Hessian.
      Call AAdd(6*NRF,XNew(NVar+NVarMM+1),XMMT(1,1+NVarMM/3),V(ICRig))
      If(IntHes) Call AMove(6*NRF,XMMT(1,1+NVarMM/3),V(ICRig))
      Call FileIO(2,-IRwMM,LenPar,V(IParam),0)
      Call AMove(3*NAtoms,CNewQ,Scr)
      NNT = (NVar*(NVar+1))/2
      Call AMove(NNT,FC,FCTemp)
      If(IntHes) Call AClear(NNT,FC)
      Call MMMicr(IOut,IPrint,.False.,ISCRF,3,MaxCyc,ConvMM,NAtoms,
     $  IAtTyp,AtChMM,NumTpS,IAtTpS,MicOpt,.True.,Scr,MaxLEF,EField,
     $  V(IParam),LenPar,Stepd,Convg,MMCyc,EMM,Scr,V(IFX),XX,XX,XX,XX,
     $  XX,XX,XX,ITrScl,JJ,V(IFFX),IRwMM,IGuess,IRFO,IPrjct,NVarMM,NVar,
     $  IBMat,FC,V(IFFTmp),BMat,V(IGM),XQuad,DBMat,X,MethQ,IDir,IFlags,
     $  V(IVecIn),.False.,HMod,IMMRFO,COrig,V(ICRig),Neg,IHVAcc,TstEig,
     $  XName,IZRedO,NVar,NTBond,NTAng,NTDih,TSOnly,0,0,SmlStp,IGetTS,0,
     $  0,IPFlag,AllowP,FMM,FMFlag,FMFlg1,NFxFlg,IHMeth,NAtoms,LJJ,JJ,
     $  XX,JJ,JJ,JJ,NAtoms,IntHes,IRwF2,QMVar,.False.,HReal,DoEle,
     $  .False.,MMRest,0,V(IV),MDV1)
      Call AMove(NNT,FCTemp,FC)
      OK = .True.
      Steped = Steped .or. OK
      Call TStamp(1,'Bot QIMOMM')
      Return
      End
*Deck QSTOut
      Subroutine QSTOut(IType,NVar,TV,XName)
      Implicit Real*8(A-H,O-Z)
C
C     Print out QST parameters.
C     Chun Yang Peng 1992
C
      Dimension TV(NVar,4), XName(*), ColLab(80)
C
      JC = 0
      Call PutCst('Tangent',ColLab,JC)
      Call PutDel(2,ColLab,JC)
      If(IType.eq.1) then
        Call PutCst('TS vect',ColLab,JC)
      else
        Call PutCst('F(1-TT)',ColLab,JC)
        endIf
      Call PutDel(2,ColLab,JC)
      Call PutCst('// Eig Forces',ColLab,JC)
      Call PutDel(2,ColLab,JC)
      Call PutCst('Eigenvalues',ColLab,JC)
      Call PutDel(2,ColLab,JC)
      Call MatPrt(TV,NVar,4,NVar,4,2,1,XName,ColLab,0,TV,0)
      Return
      End
*Deck RedStp
      Subroutine RedStp(IOut,IPrint,NEig,h,F,Eig,RLamda)
      Implicit Real*8(A-H,O-Z)
C
C     Given the trust radius h, forces F (in basis of Hessian
C     eigenvectors) and Hessian eigenvalues Eig, compute RLamda
C     such that the resulting Newton-Raphson step with eigenvalues
C     level-shifted by RLamda locates a minimum on the sphere
C     of radius h.  This is done by minimizing the error in the
C     function (R-h)**2 with respect to the scaling parameter Lamda,
C     using the fact that the step can be computed explicitly from
C     F and Eig given Lamda.
C
C     Mike Frisch, July 1986.
C
      Parameter (MaxIt=1000)
      Dimension F(*), Eig(*)
      Save Zero, Conv, Two, Three
      Data Zero/0.0d0/, Two/2.0d0/, Three/3.0d0/, Conv/1.d-6/
 1000 Format(' RedStp:  It=',I2,' RLamda=',1PD10.3,' R=',1PD9.3,/,
     $       '          Err=',1PD9.3,' DErr=',1PD10.3,' DDErr=',1PD9.3)
 1010 Format(' Maximum of',I3,' iterations exceeded in RedStp.')
C
      RLamda = Zero
      Do 100 It = 1, MaxIt
        RSq = Zero
        DR = Zero
        DDR = Zero
        Do 10 I = 1, NEig
          RSq = RSq + F(I)**2/(RLamda-Eig(I))**2
          DR = DR - F(I)**2/(RLamda-Eig(I))**3
   10     DDR = DDR + Three*F(I)**2/(RLamda-Eig(I))**4
        R = Sqrt(RSq)
        E = (h-R)**2
        DR = DR / R
        DDR = DDR/R - DR**2/R
        DE = DR * Two * (R-h)
        DDE = Two*(R-h)*DDR + Two*DR**2
        If(IPrint.gt.0) Write(IOut,1000) It, RLamda, R, E, DE, DDE
        If(Abs(E).lt.Conv) Return
        RLamda = RLamda - DE/DDE
  100   Continue
      Write(IOut,1010) MaxIt
      Call Lnk1E(0)
      Return
      End
*Deck RemPnt
      Subroutine RemPnt(IOut,Prnt,DidInv,DISA,CIn,COut,Map,NUsed,NDIIS,
     $  CThrHd,DidRem)
      Implicit Real*8(A-H,O-Z)
C
C     DISA   : DIIS A Matrix. If a point is removed, this matrix will be repacked
C     CIn    : DIIS coefficients (input)
C     COut   : DIIS coefficients (output)
C     Map    : Used to record the original indices of the DIIS points
C     NUsed  : Number of points used. If a point is removed, NUsed will be modfied
C     NDIIS  : Starting number of DIIS points
C     CThrHd : Threshold of the coefficients.
C
      Logical Prnt,DidInv,DidRem
      Dimension DISA(*),CIn(*),COut(*),Map(*)
 1000 Format(' Point number',I4,' has DIIS coefficient',1PD10.2,
     $  ' and has been removed')
 1020 Format(' DIIS inversion failure, remove point',I4,'.')
C
C     Check the absolute value of the coefficients
C
      DidRem = .False.
      NDim = NUsed + 1
      If(DidInv) then
        CMax = ArMax1(CIn,NUsed,.True.,IMax)
      else
        IMax = 0
        CMax = CThrHd
        endIf
      If(CMax.lt.CThrHd) then
        If(NUsed.eq.NDIIS) then
          Call AMove(NUsed,CIn,COut)
        else
          Call AClear(NDIIS,COut)
          ICount = 0
          Do 10 I = 1, NDIIS
            If(Map(I).gt.0) then
              ICount = ICount + 1
              COut(I) = CIn(ICount)
              endIf
   10       Continue
          endIf
      else
        If(IMax.eq.0) then
          IMax = NUsed
          Write(IOut,1020) IMax
        else
          ISt = IMax + 1
          Do 20 I = ISt, NUsed
            If(Abs(CIn(I)).gt.CThrHd) IMax = I
   20       Continue
          If(IMax.eq.1) IMax = NUsed
          If(Prnt) Write(IOut,1000) IMax, CIn(IMax)
          endIf
        If(NUsed.ne.NDIIS) then
          ICount = 0
          Do 30 I = 1, NDIIS
            If(Map(I).gt.0) ICount = ICount + 1
            If(ICount.eq.IMax) then
              Map(I) = 0
              Goto 40
              endIf
   30       Continue
        else
          Map(IMax) = 0
          endIf
   40   Call DelBas(NDim,NDim,IMax+1,IMax+1,DISA,DISA)
        Call Pcck(0,KK,DISA,NDim,NDim,NDim-1,NDim-1)
        NUsed  = NUsed - 1
        DidRem = .True.
        endIf
      Return
      End
*Deck RevEnM
      Subroutine RevEnM(NGDIIS,PotE,AMEn,AMTrc,AMTmp)
      Implicit Real*8 (A-H,O-Z)
      Dimension AMEn(0:NGDIIS,0:NGDIIS), AMTrc(0:NGDIIS,0:NGDIIS),
     $  AMTmp(0:NGDIIS,0:NGDIIS),PotE(NGDIIS)
C
      NASq = (NGDIIS+1)*(NGDIIS+1)
      Call AMove(NASq,AMEn,AMTmp)
      Do 100 I = 1, NGDIIS
        Do 100 J = I, NGDIIS
          AMTmp(I,J) = AMEn(NGDIIS-I+1,NGDIIS-J+1)
          AMTmp(J,I) = AMTmp(I,J)
  100     Continue
      Call AMove(NASq,AMTmp,AMEn)
      Call AMove(NASq,AMTrc,AMTmp)
      Do 200 I = 1, NGDIIS
        Do 200 J = 1, NGDIIS
          AMTmp(I,J) = AMTrc(NGDIIS-I+1,NGDIIS-J+1)
  200     Continue
      Call AMove(NASq,AMTmp,AMTrc)
      Do 300 I = 1, NGDIIS
  300   AMTmp(I,0) = PotE(NGDIIS-I+1)
      Call AMove(NGDIIS,AMTmp(1,0),PotE)
      Return
      End
*Deck RmEnPt
      Subroutine RmEnPt(IOut,Prnt,NUsed,NDIIS,IWarn,AMEn,AMTrc,PotE,Map,
     $  COut,DidSkp)
      Implicit Real*8(A-H,O-Z)
C
C     DISA   : DIIS A Matrix. If a point is removed, this matrix will be repacked
C     PotE    : DIIS coefficients (input)
C     COut   : DIIS coefficients (output)
C     Map    : Used to record the original indices of the DIIS points
C     NUsed  : Number of points used. If a point is removed, NUsed will be modfied
C     NDIIS  : Starting number of DIIS points
C     CThrHd : Threshold of the coefficients.
C
      Logical Prnt,DidSkp
      Dimension AMEn(*), AMTrc(*), Map(*), COut(*), PotE(*)
      Save One
      Data One/1.d0/
 1000 Format(' Point number',I4,' has energy ',D18.10,
     $  ' and has been removed')
 1010 Format(' Point number',I4,' has energy ',D18.10,
     $  ' -- DIIS solution rejected.')
 1020 Format(' IWarn = ', I4)
 1030 Format(' Energies = ',5G20.12)
C
      If(Prnt) Write(IOut,1020) IWarn
      DidSkp = .False.
      If(IWarn.eq.0) then
        Call AClear(NDIIS,COut)
        ICount = 0
        Do 300 I = 1, NDIIS
          If(Map(I).gt.0) then
            ICount = ICount + 1
            COut(I) = AMTrc(ICount+1)
            endIf
  300     Continue
        Return
        endIf
C
C     Remove the point with the highest energy
C
      NDim = NUsed + 1
      IMax = 0
      RJunk = ArMax1(PotE,NUsed,.False.,IMax)
      If(Prnt) Write(IOut,1030) (PotE(I),I=1,NUsed)
      If(IMax.eq.NUsed.and.NUsed.gt.2)
     $  RJunk = ArMax1(PotE,NUsed-1,.False.,IMax)
      If(IMax.eq.NUsed) then
        Write(IOut,1010) IMax, PotE(IMax)
        NUsed = 1
        Return
      else
        If(NUsed.ne.NDIIS) then
          ICount = 0
          Do 100 I = 1, NDIIS
            If(Map(I).gt.0) ICount = ICount + 1
            If(ICount.eq.IMax) then
              Map(I) = 0
              Goto 200
              endIf
  100       Continue
        else
          Map(IMax) = 0
          endIf
  200   Call DelBas(NDim,NDim,IMax+1,IMax+1,AMEn,AMEn)
        Call Pcck(0,KK,AMEn,NDim,NDim,NDim-1,NDim-1)
        Call DelBas(NDim,NDim,IMax+1,IMax+1,AMTrc,AMTrc)
        Call Pcck(0,KK,AMTrc,NDim,NDim,NDim-1,NDim-1)
        If(Prnt) Write(IOut,1000) IMax, PotE(IMax)
        Do 400 I = IMax, (NUsed-1)
  400     PotE(I) = PotE(I+1)
        NUsed  = NUsed - 1
        Call AClear(NUsed+1,AMTrc)
        AMTrc(NUsed+1) = One
        DidSkp = .True.
        endIf
      Return
      End
*Deck SavePt
      Subroutine SavePt(Energy,NVar,IGetFC,AlUnit,INeg,Rises,F,Frcnst,
     $  ESave,IC,ITU,X,FF,XX,FC,FC2,SwapL,DidSwap,NVarT,NVrMMX,NDimMM,
     $  FMM,XMM,FFMM,XXMM,FMMQ,FFMMQ,EQM,EQMS)
      Implicit Real*8(A-H,O-Z)
C
C     Original version by H. B. Schlegel.  Rewritten by M. Frisch.
C     Push the new energy, position, forces, and possibly force
C     constants onto the stacks.  Be careful to keep the best point
C     at the top of the stack.  Rises is returned .True. if the
C     Energy (INeg=0) or gradient (INeg=1) of the new point is not
C     the lowest so far, and MaxRis tries to lower it have been made.
C
C     When climbing, the linear search should be carried out with the
C     current point (IStep) and previous point (IStep-1); when not
C     climbing, the linear search should be done with the current point
C     (IStep) and the best previous points (not including the climbing
C     points.  CYP Mar, 1993.
C
      Logical Prnt, Exit, Rises, SwapL, DidSwap, AlUnit
      Dimension F(*), Frcnst(*), ESave(*), IC(*), X(*), FC2(*),
     $  FF(NVar,*), XX(NVar,*), FC(*), FMM(*), XMM(*), FFMM(NVarT,*),
     $  XXMM(NDimMM,*), EQMS(3), FMMQ(*), FFMMQ(NDimMM,*), ITU(*)
      Logical HavDEP, Large
      Common /OptGrd/ ConvF,FMaxT,DXMaxT,RMax,RMin,RLim,EigMax,Eigmin,
     $  DEPred,FSwtch,FncErr,GrdErr,FncCnv,XXIRC,NLevel,NStep,IStep,
     $  NMax,NP,Neg,Prnt,Exit,IPSav1,IPSav2,HavDEP,NGDIIS,IRisen,MaxRis,
     $  Large,ILevel,ISStep,IOffs,IMDS,MaxSV,MaxSV2,NVrMMT
C
      DidSwap = .False.
      NP = Min(NP+1,NMax)
      NVrMMY = NVrMMX
      If(NVrMMY.eq.1) NVrMMY = 0
      If(NP.ge.2) then
        NN = NP - 1
        Call AMove1(NN,0,1,ESave)
        Call IMove1(NN,0,1,IC)
        Call IMove1(NN,0,1,ITU)
        Call AMove1(NN*NVar,0,NVar,FF)
        Call AMove1(NN*NVar,0,NVar,XX)
        If(NVrMMY.gt.0) then
          Call AMove1(NN*NVarT,0,NVarT,FFMM)
          Call AMove1(NN*NVrMMY,0,NVrMMY,XXMM)
          Call AMove1(NN*NVrMMY,0,NVrMMY,FFMMQ)
          Call AMove1(NN,0,1,EQMS)
          endIf
        endIf
      ESave(1) = Energy
      IC(1) = IStep
      ITU(1) = 0
      Call AMove(NVar,X,XX)
      Call AMove(NVar,F,FF)
      If(NVrMMY.gt.0) then
        Call AMove(NVarT,FMM,FFMM)
        Call AMove(NVrMMY,XMM,XXMM)
        Call AMove(NVrMMY,FMMQ,FFMMQ)
        EQMS(1) = EQM
        endIf
C
C     Make sure that the best point other than the current
C     one is second from the top of the stack.
C
      If(NP.gt.2.and.SwapL) then
        IP = 2
        GNSave = SProd(NVar,FF(1,IP),FF(1,IP))
        Do 10 I = 3, NP
          GNI = SProd(NVar,FF(1,I),FF(1,I))
          If(INeg.eq.0.and.ESave(I).lt.ESave(IP)) IP = I
          If(INeg.ne.0.and.GNI.lt.GNSave) IP = I
          GNSave = Min(GNI,GNSave)
   10     Continue
        DidSwap = IP.ge.3
        Do 20 I = IP, 3, -1
          Call SwapPt(I,I-1,NVar,NVrMMY,ESave,IC,XX,FF,EQMS,XXMM,FFMM,
     $      FFMMQ)
   20     Continue
        endIf
      NVarTT = (NVar*(NVar+1))/2
      If(IC(2).eq.IPSav1.and..not.AlUnit) then
        Call AMove(NVarTT,FC,FC2)
        IPSav2 = IPSav1
        endIf
      If(IC(2).ne.IPSav2) IPSav2 = 0
      If(IGetFC.ne.0.and..not.AlUnit)
     $  Call AMove(NVarTT,Frcnst,FC)
      IPSav1 = 0
      IF(IGetFC.eq.4) IPSav1 = NP
      Rises = .False.
      If(NP.eq.1) Return
      GNNew = SProd(NVar,FF(1,1),FF(1,1))
      GNOld = SProd(NVar,FF(1,2),FF(1,2))
C
C     The Rises flag really indicates which things have gotten so
C     bad that only linear minimization is to be done until a lower
C     point is reached.  First check if the point is really worse.
C
      Rises = (INeg.eq.0.and.(ESave(1)-FncErr).gt.(ESave(2)+FncErr))
     $  .or.(INeg.ne.0.and.(GNNew.gt.GNOld))
      If(Rises.and.SwapL) then
        IRisen = IRisen + 1
      else
        IRisen = 0
        endIf
C
C     Now check whether this has happened the allowed number of times
C     in a row.  If so, make sure that the linear search is done using
C     the best point (#2 on the stack) and the second-best point.  This
C     logic currently only works for MaxRis=1 (linear search after second
C     failure).
C
      Rises = Rises.and.IRisen.gt.MaxRis
      If(Rises) then
        If(MaxRis.gt.0.and.IRisen.eq.(MaxRis+1)) then
          If(INeg.eq.0) then
            ROld = ESave(3)
            RNew = ESave(1)
          else
            ROld = SProd(NVar,FF(1,3),FF(1,3))
            RNew = SProd(NVar,FF(1,1),FF(1,1))
            endIf
          If(RNew.gt.ROld) then
            Call SwapPt(1,3,NVar,NVrMMY,ESave,IC,XX,FF,EQMS,XXMM,FFMM,
     $        FFMMQ)
            Energy = ESave(1)
            Call AMove(NVar,XX,X)
            Call AMove(NVar,FF,F)
            If(NVrMMY.gt.0) then
              Call AMove(NVarT,FFMM,FMM)
              Call AMove(NVrMMY,XXMM,XMM)
              Call AMove(NVrMMY,FFMMQ,FMMQ)
              EQM = EQMS(1)
              endif
            endIf
          endIf
        endIf
      Return
      End
*Deck SimMem
      Subroutine SimMem(NP,NQ,IStart,IDQQS,IA,IAInv,IB,IC,ICM,IHSc,
     $  IDQR,IW,IEnd)
      Implicit Integer(A-Z)
C
C     Allocate the scratch space and variables for GDIIS.
C
      NDim   = NP + 1
      IDQQS  = IStart+ 1
      IA     = IDQQS + NQ*NP
      IAInv  = IA    + NDim*NDim
      IB     = IAInv + NDim*NDim
      IC     = IB    + NDim
      ICM    = IC    + NDim
      IHSc   = ICM   + NDim
      IDQR   = IHSc  + NQ
      IW     = IDQR  + NQ
      IEnd   = IW    + NQ*6
      Return
      End
*Deck TrH2Hx
      Subroutine TrH2Hx(IOut,IPrint,NTBond,NTAng,NTRed,NAtoms,IZRed,C,
     $  Value,IB,B,ISB,SB,IntVec,H,HX,Scr,LenSB,LenCB,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     Transform the Hessian represented in the redundant coordinate
C     space to a smaller, non-redundant space.
C
C     Constraints and rotational part by O. Farkas, 1997, WSU.
C
      Dimension IZRed(4,*), C(3,*), Value(*), IB(4,*), B(*),
     $  IntVec(*), H(*), HX(*), Scr(*), IRB(3), ISB(*), SB(3,*), V(MDV)
      Save Zero, One, Half, Two, EigMax
      Data Zero/0.0d0/,One/1.D0/,Half/0.5D0/,Two/2.D0/,EigMax/1.D4/
      LInd(I,J)  = ((Max(I,J))*(Max(I,J)-1))/2 + Min(I,J)
C
      If(NTRed.le.0) Return
      NAt3 = 3*NAtoms
      NAtT = (NAt3*(NAt3+1))/2
      IRB(1) = 12*NTRed
      IRB(2) = IRB(1) + NAt3
      IRB(3) = IRB(2) + NAt3
      II = 0
      Do 10 I = 1, NTRed
        II = II + I
        If(IntVec(I).lt.0) H(II) = EigMax
   10   Continue
      Call AClear(NAtT,HX)
      Call FormSB(IOut,IPrint,NTBond,NTAng,NTRed,IZRed,C,Value,IB,B,ISB,
     $  SB,IntVec,NAtoms,Scr,LenSB,LenCB,V,MDV)
      I0  = NTRed + 1
      ISh = ISB(1) - 1
      Do 100 I = 1, NTRed
        Call AClear(NAt3,Scr)
        Do 120 J = 1, NTRed
          HJI = H(LInd(J,I))
          IP  = ISB(J+1)
          Len = ISB(J+2) - IP
          IIP = I0 + IP
          IP  = ISh + IP
          Do 120 K = 1, Len
            KK = ISB(IIP+K)
            L0 = (KK-1)*3
            Do 120 L = 1, 3
  120         Scr(L0+L) = Scr(L0+L) + SB(L,IP+K) * HJI
        IP  = ISB(I+1)
        Len = ISB(I+2) - IP
        IIP = I0+IP
        IP  = ISh+IP
        Do 130 K = 1, Len
          KK = ISB(IIP+K)
          L0 = (KK-1)*3
          Do 130 L = 1, 3
            LL = L0 + L
            SBLK = SB(L,IP+K) * Half
            Do 130 LLL = 1, NAt3
              DHX = SBLK * Scr(LLL)
              If(LL.eq.LLL) DHX = DHX*Two
              LX = LInd(LLL,LL)
  130         HX(LX) = HX(LX) + DHX
  100   Continue
      If(IZRed(2,1).lt.0) Return
C
C     Put matrix elements related to translation and rotation
C
      RNA  = GFloat(NAtoms)**2
      If(RNA.ne.Zero) RNA = One / RNA
      Do 150 IX = 1, 3
        Do 150 I = IX, NAt3, 3
          II = (I*(I-1))/2
          Do 150 J = IX, I, 3
  150       HX(II+J) = HX(II+J) + RNA
      IRB1 = IRB(1)
      IRB2 = IRB(2)
      IRB3 = IRB(3)
      Do 170 I = 1, NAt3
        RBI1 = B(IRB1+I)
        RBI2 = B(IRB2+I)
        RBI3 = B(IRB3+I)
        IX = (I*(I-1))/2
        Do 170 J = 1, I
  170     HX(IX+J) = HX(IX+J) + B(IRB1+J)*RBI1 + B(IRB2+J)*RBI2 +
     $      B(IRB3+J)*RBI3
      Return
      End
*Deck TstBck
      Logical Function TstBck(NUsed,CM)
      Implicit Real*8(A-H,O-Z)
C
C     Check whether DIIS stepped back.
C
      Dimension CM(*)
      Save ThrBak, One
      Data ThrBak/0.8d0/, One/1.0d0/
C
      TB1 = One - ThrBak
      TstBck = .False.
      Do 10 I = 2, NUsed
   10   TstBck = TstBck.or.CM(I).ge.TB1
      TstBck = TstBck.and.CM(1).lt.ThrBak
      Return
      End
*Deck UpCOrg
      Subroutine UpCOrg(CNew,Temp,COrig,NAtoms,NAtNew,MicOpt,XMMT,COM,
     $  DT,ICnt,T,NRF,NVar,NVarMM,XNew,IMMCRS,RotQM,ReOrdr)
      Implicit Real*8(A-H,O-Z)
C
C     This routine takes care of several changes in geometry steps for QuadMac.
C     When the internal coordinates have changed, it can (optionally) update COrig
C     accordingly (including the rotation back to the original geometry).
C     It can then take care of a step in the coordinates for the MM rigid fragments.
C
      Logical RotQM, ReOrdr
      Dimension CB(3), CC(3), Temp(*), CNew(3,NAtoms), COrig(3,*),
     $  MicOpt(*), XMMT(3,*), COM(3,*), DT(*), ICnt(*), T(3,3,10,*),
     $  XNew(*)
C
C     Orig is in conventional ordering, while CNew has the QM centers first. Reorder.
      Call AMove(3*NAtoms,CNew,Temp)
      IQM = 1
      IMM = 1 + 3*NAtNew
      Do 10 I = 1, NAtoms
        If(MicOpt(I).ge.0) then
          Call AMove(3,Temp(IMM),CNew(1,I))
          IMM = IMM + 3
        else
          Call AMove(3,Temp(IQM),CNew(1,I))
          IQM = IQM + 3
          endIf
   10   Continue
C
C     Rotate QM region back to COrig orientation
      If(RotQM) then
        Call RigidX(NAtoms,CNew,XMMT(1,1+NVarMM/3),COM,MicOpt,DT,ICnt,T)
        Call Trspn2(3,T(1,1,1,NRF))
        Do 20 I = 1, NAtoms
          If(MicOpt(I).lt.0) then
            Call ASub(3,CNew(1,I),COM(1,NRF),CB)
            Call RotF(1,T(1,1,1,NRF),CB,CC)
            Call ASub(3,CC,XMMT(1,NRF*2-1+NVarMM/3),CB)
            Call AAdd(3,CB,COM(1,NRF),COrig(1,I))
            endIf
   20     Continue
        endIf
C
C     Update MMCRS file
      Call AAdd(6*NRF,XNew(NVar+NVarMM+1),XMMT(1,1+NVarMM/3),Temp)
      If(ITqry(IMMCRS).gt.0) then
        Call FileIO(1,-IMMCRS,6*NRF,Temp,2)
        Call FileIO(1,IMMCRS,3*NAtoms,COrig,0)
        endIf
C
C     Process translation/rotation of rigid fragments.
      Call RigidX(NAtoms,COrig,Temp,COM,MicOpt,DT,ICnt,T)
      IP = NVar + 1
      Do 30 I = 1, NAtoms
        IFrag = MicOpt(I) - 1
        If(RotQM.and.MicOpt(I).lt.0) IFrag = NRF
        If(MicOpt(I).eq.0) then
          Call AAdd(3,XNew(IP),CNew(1,I),CNew(1,I))
          IP = IP + 3
        else if(IFrag.gt.0) then
          Call ASub(3,COrig(1,I),COM(1,IFrag),CB)
          Call RotF(1,T(1,1,1,IFrag),CB,CC)
          Call AAdd(3,CC,Temp(1+(IFrag-1)*6),CB)
          Call AAdd(3,CB,COM(1,IFrag),CNew(1,I))
          endIf
   30   Continue
      If(.not.ReOrdr) Return
C
C     Reorder so that the QM centers come first.
      IQM = 1
      IMM = 1 + 3*NAtNew
      Do 40 I = 1, NAtoms
        If(MicOpt(I).ge.0) then
          Call AMove(3,CNew(1,I),Temp(IMM))
          IMM = IMM + 3
        else
          Call AMove(3,CNew(1,I),Temp(IQM))
          IQM = IQM + 3
          endIf
   40   Continue
      If(IQM.ne.(1+3*NAtNew).or.IMM.ne.(1+3*NAtoms))
     $  Call GauErr('Thom is an idiot in UpCOrg.')
      Call AMove(3*NAtoms,Temp,CNew)
      Return
      End
*Deck UpDate
      Subroutine UpDate(IOut,IPrint,NVar,NDim,NVTT,NTBond,NTAng,NTDih,
     $  NPath,NN,NPStar,XX0,XXM,XXP,XX,ESv0,EsvM,EsvP,FS0,FS,FF0,FFM,
     $  FFP,FF,FC0,Proj,IC0,IC,MStep,Scr,DoStop,DoAlong,MethAl,MethHi,
     $  IMax)
      Implicit real*8(a-h,o-z)
C
C     Path relaxation.  Updating routine for a point on the path.
C
      Dimension XX0(Nvar,NDim),XXM(Nvar),XXP(Nvar),XX(Nvar,NDim),
     $  FF0(Nvar,NDim),FFM(Nvar),FFP(Nvar),FF(Nvar,NDim),FS0(NDim),
     $  FS(NDim+1),FC0(NVTT),Proj(NVTT),IC0(NDim),IC(NDim),Scr(Nvar,7),
     $  XS(3),Values(3,2),XMat(6,6),IS(4,6),YS(6),FScr(6),Q(36),R(36),
     $  QR(36)
      Logical Prnt,Exit,HavDEP,OK,DoStop,DoAlong,DOSwap,Debug,Large
      Common /OptGrd/ ConvF,FMaxT,DXMaxT,RMax,RMin,RLim,EigMax,EigMin,
     $  DEPred,FSwtch,FncErr,GrdErr,FncCnv,XXIRC,NLEVEL,NStep,IStep,
     $  NMax,NP,Neg,Prnt,Exit,IPSav1,IPSav2,HavDEP,NGDIIS,IRisen,MaxRis,
     $  Large,ILevel,ISStep,IOffs,IMDS,MaxSv,MaxSv2,NVrMMT
      Save Zero,One,Two,Thresh
      Data Zero /0.0d0/, One /1.0d0/, Two /2.0d0/,Thresh/1.0d-05/
 1000 Format(' Hessian update for path ')
 1010 Format('++++++++++++++++++++++++++')
 1020 Format(' Point ',I2)
 1030 Format(' Update along the path ')
 1045 Format(' Second derivative along the tangent : Old = ',F10.5,
     $  ' New = ',F10.5)
 1050 Format(' Update history ')
 2000 Format(' No Update for point ',I3,': Point converged')
 5000 Format(' Maximum Hessian matrix element =',I10)
 5010 Format(' DM=',F10.5,'  DP=',F10.5)
 5020 Format(' DotGP=',F10.5)
 5030 Format(' DotGR=',F10.5)
 5050 Format(' Maximum Hessian matrix element = ',F10.5)
C
      Debug=.false.
      If(DoStop.and.MStep.ne.IStep) then
        Write(IOut,2000) NN
        Return
        endIf
      Write(IOut,1000)
      Write(IOut,1020) NN
      Write(IOut,1010)
      If(DoAlong) then
        Call AMove(NVar,XX0(1,1),XX(1,1))
        Call AMove(NVar,FF0(1,1),FF(1,1))
        FS(1)=Esv0
        Write(IOut,1030)
        Write(IOut,1010)
        II=1
        IC(1)=NN
        If(NN-1.ge.1) then
          MM=NN-1
          II=II+1
          Call AMove(NVar,XXM(1),XX(1,II))
          Call AMove(NVar,FFM(1),FF(1,II))
          FS(II)=ESvM
          IC(II)=MM
          endIf
        If(NN+1.le.NPath) then
          MM=NN+1
          II=II+1
          Call AMove(NVar,XXP(1),XX(1,II))
          Call AMove(NVar,FFP(1),FF(1,II))
          FS(II)=ESvP
          IC(II)=MM
          endIf
        MyUpd=MethAL
        NP=II
        If(NP.gt.1) then
          Call FixHis(IOut,IPrint,Nvar,NP,XX,Scr,NTBond,NTAng,NTDih)
          If((ESvP.gt.ESvM).and.NP.eq.3)
     $      Call SwapPt(3,2,NVar,0,FS,IC,XX,FF,XX,XX,XX,XX)
          XMax=ARRMax(FC0(1),NVTT,.true.)
          endIf
        If(Debug) Write(Iout,5000) XMax
        If(NP.gt.1.and.MStep.le.IMax) Call D2CorP(IOut,MyUpd,NVar,NP,1,
     $    RMin,Rmax,GrdErr,IC,XX,FF,FC0(1),scr(1,1),Scr(1,2),Scr(1,4))
        Call AClear(NVar,Scr)
        Call AClear(NVar,Scr(1,2))
        Call AClear(NVar,Scr(1,3))
        Call AClear(NVar,Scr(1,4))
        Call AClear(NVar,Scr(1,5))
        If(NN.gt.1.and.NN.lt.NPath) Call Cnsst4(2,NVar,Scr,Proj(1),
     $    XX0(1,1),XXM(1),XXP(1),FF0(1,1),GrdTmp,Scr,NTBond,NTAng,NTDih,
     $    RP,RR,RPR,ArclnR,ArclnP,0,1,.false.)
        If(NN.gt.1.and.NN.lt.NPath) then
          DM=-ArclnR
          DP=ArclnP
          If(Debug) Write(IOut,5010) DM,DP
          DotTG=-Grdtmp
C         Get Tangent at Next Point T(P)= -sT(X)+XP
          Call ACasB(Nvar,Scr(1,1),Scr(1,4),Scr(1,1),-RP)
          Tmp=One/Sqrt(SProd(Nvar,Scr(1,1),Scr(1,1)))
          Call Ascale(Nvar,tmp,Scr(1,1),Scr(1,1))
          DotGP=-SProd(Nvar,Scr(1,1),FFM(1))
          If(Debug) Write(IOut,5020) DotGP
C         Get Tangent at Previous Point
          Call AScale(NVar,-One,Scr(1,2),Scr(1,2))
          Call ACasB(Nvar,Scr(1,2),Scr(1,4),Scr(1,2),-RR)
          Tmp=One/Sqrt(SProd(Nvar,Scr(1,2),Scr(1,2)))
          Call Ascale(Nvar,tmp,Scr(1,2),Scr(1,2))
          DotGR=-SProd(Nvar,Scr(1,2),FFP(1))
          If(Debug) Write(IOut,5030) DotGR
          XS(1)=-ArcLnR
          XS(2)=Zero
          XS(3)=ArcLnP
          Values(1,1)=EsvP
          Values(2,1)=Esv0
          Values(3,1)=EsvM
          If(Abs(DotGR).lt.Thresh) DotGR=Zero
          If(Abs(DotGP).lt.Thresh) DotGP=Zero
          If(Abs(DotTG).lt.Thresh) DotTG=Zero
          Values(1,2)=DotGR
          Values(2,2)=DotTG
          Values(3,2)=DotGP
          Call PolFit(Iout,1,3,2,XS,values,XMat,IS,FScr,YS,OK,Q,R,QR)
          If(Debug) Call OutMtS(IOut,'YS in Update:',0,0,YS,1,6,1,6)
          B=Two*YS(4)
          Test=DotGR*DotTG*DotGP
          ITest=Max(NN,NPstar)-Min(NN,NPstar)
          If(test.lt.Zero.and.Itest.ne.1.and..not.DoStop) OK=.false.
          Call MultLV(.True.,1,0,NVar,NVar,FC0(1),Scr(1,4),XX(1,2))
          EigTmp=SProd(NVar,Scr(1,4),XX(1,2))
          Test=GRdTmp/Sqrt(GFloat(NVar))
          XTmp =(EigTmp+Sqrt(EigTmp**2+(GrdTmp+GrdTmp)**2))/Two
          XTmp = -GrdTmp/(XTmp-EigTmp)
          If((DM.lt.-RMax/Two.or.DP.gt.RMax/Two).and.OK.and.
     $      MStep.le.IMax) then
            Write(IOut,1045) Eigtmp,B
            IJ = 1
            W1E = (Eigtmp+B)
            Do 100 I = 1, NVar
              Do 100 J = 1, I
                FC0(IJ) = FC0(IJ) + W1E*Scr(I,4)*Scr(J,4) -
     $            (XX(I,2)*Scr(J,4)+Scr(I,4)*XX(J,2))
  100           IJ = IJ + 1
            endIf
          endIf
        endIf
      Write(IOut,1010)
      II=MStep
      NP=II
      If(NP.gt.1) then
C       Feed history of point NN
        Call AMove(NVar*MStep,XX0(1,1),XX(1,1))
        Call AMove(NVar*MStep,FF0(1,1),FF(1,1))
        FS(1)=Esv0
        Call AMove(MStep,FS0(1),FS(2))
        Call IMove(IStep,IC0(1),IC)
        Write(IOut,1050)
        Write(IOut,1010)
        DoSwap=.false.
        If(MethHi.eq.7.and.NP.lt.5) MethHi=4
        If((FS(1).gt.FS(2)).and.DoSwap)
     $    Call SwapPt(1,2,NVar,0,FS,IC,XX,FF,XX,XX,XX,XX)
        Call FixHis(IOut,IPrint,Nvar,NP,XX,Scr,NTBond,NTAng,NTDih)
        XMax=ARRMax(FC0(1),NVTT,.true.)
        If(Debug) Write(IOut,5050) XMax
        MyUpd=MethHi
        If(IC(1).ne.IC(2)) then
          If(MethHi.eq.3)Call AMove(Nvar,XX0(1,NDim),Scr(1,4))
          Call D2CorP(IOut,MyUpd,NVar,NP,1,RMin,Rmax,GrdErr,IC,XX,FF,
     $      FC0(1),scr(1,1),Scr(1,2),Scr(1,4))
          endIf
        endIf
      Return
      End
*Deck UpdMMO
      Subroutine UpdMMO(IOut,IPrint,MaxAtm,NAtoms,OptAt,TrVec,TRot,CNew,
     $  CTemp)
      Implicit Real*8(A-H,O-Z)
C
C     Update /Mol/ with the results of the microiterations.
C
#include "commonmol2.inc"
      Integer OptAt(*)
      Dimension TRot(3,3), CNew(3,NAtoms), CTemp(3,NAtoms), TrVec(*)
      Common /MUnit/ IUnit(20)
 1000 Format(' Updating /Mol/ with results of microiterations.')
 1010 Format(' OptAt=',50I1)
C
      If(IPrint.ge.1) Write(IOut,1000)
      If(IPrint.ge.2) Write(IOut,1010) (OptAt(I),I=1,NAtoms)
      Call RWMol(-2,-IUnit(9))
      Call RotF1(NAtoms,TRot,CB)
      Call Transl(1,3,NAtoms,TrVec,CNew)
      Do 10 I = 1, NAtoms
        If(OptAt(I).ge.0) then
          CB(1,I) = CNew(1,I)
          CB(2,I) = CNew(2,I)
          CB(3,I) = CNew(3,I)
          endIf
   10   Continue
      Call Transl(-1,3,NAtoms,TrVec,CNew)
      Call Trspn2(3,TRot)
      Call RotF1(NAtoms,TRot,CB)
      Call Trspn2(3,TRot)
      Call RWMol(-1,-IUnit(9))
      Return
      End
*Deck UpdRed
      SubRoutine UpdRed(IOut,IPrint,IPrt1,Dored,NVar,NTBond,NTAng,NTDih,
     $  IZred,Value,NAtoms,C,DQ,Cnew,BMat,IBMat,VecInt,FCTemp,IDimSc,
     $  Scr,XNew,IFlags,IAn,AtmChg,TRot,TrRMat,X,Q0,CnvFX,ConvF,CnvFMX,
     $  F,UseSym,UpCart,CP,GMax,GRMS,DxMaxC,DxRMSC,Done,IOSLE1,V,MDV)
      Implicit real*8(A-H,O-Z)
C
C     Transform redundant internal displacement to cartesian and back
C     to red. internal.  Check convergence and update coordinates.
C
      Logical DoRed,UseSym,OK,Done,UpCart,Print
      Dimension IZred(4,*),Value(*),C(3,*),DQ(Nvar),CNew(3,*),IAn(*),
     $  BMat(3,4,*),IBMat(4,*),VecInt(*),FCTemp(*),Scr(IDimSc,8),
     $  XNew(*),IFlags(*),ATmChg(*),TRot(3,3),TrRMat(6),X(*),Q0(NVar),
     $  F(Nvar),CP(3,NAtoms,2),IReslt(3),XX(1),V(MDV)
      Save XX
      Data XX/0.d0/
 1000 Format('         Item               Value     Threshold  ',
     $       'Converged?')
 1002 Format(' Maximum Force            ',F8.6,5X,F8.6,5X,3A1)
 1004 Format(' RMS     Force            ',F8.6,5X,F8.6,5X,3A1)
 1006 Format(' Maximum Displacement     ',F8.6,5X,F8.6,5X,3A1)
 1008 Format(' RMS     Displacement     ',F8.6,5X,F8.6,5X,3A1)
 1100 Format(' MatCor failed.')
C
      Print = .False.
      NAt3 = NAtoms*3
      If(.not.DoRed) Call GauErr('Wrong Type of Coord.')
      Call RedCar(IOut,IPrt1,IOSLE1,NVar,NTBond,NTAng,NTDih,IZRed,Value,
     $  NAtoms,C,DQ,CNew,BMat,IBMat,VecInt,FCTemp,Scr,XNew,IFlags,
     $  .True.,IFail,0,V,MDV)
      Call ClnCor(IOut,IPrint,.False.,UseSym,NAtoms,IAn,AtmChg,C,0,Scr,
     $  .True.,TRot,CNew,Scr)
      Call MatCor(IOut,IPrt1,.True.,NAtoms,C,CNew,TrRMat,Scr,OK)
      If(.not.OK) Write(IOut,1100)
      Call Asub(Nat3,C,Scr,C)
      DXMaxC = ArrMax(C,NAt3,.True.)
      DXRMSC = Sqrt(Abs(SProd(NAt3,C,C)/GFloat(NVar)))
      Call AMove(NAt3,Scr,C)
      If(UpCart) then
        Call AMove(NAt3,CP(1,1,1),CP(1,1,2))
        Call AMove(NAt3,C,CP)
        endIf
      Call RedVal(IOut,IPrint,.True.,C,IZRed,X,NVar,NTBond,NTAng,NTDih,
     $  V,MDV)
      If(NTDih.ge.1) Call FixDih(NTDih,X(NTBond+NTAng+1))
      Call ASub(Nvar,X,Q0,DQ)
      If(NTDih.ge.1) Call FixDq(IOut,IPrint,1,NVar,NTBond,NTAng,NTDih,
     $  DQ,XX)
      Call AClear(NVar,XNew)
      ConvX = CNVFX*ConvF
      T1 = ConvF * CnvFMx
      T2 = ConvX * CnvFMx
      Call RMSVec(Nvar,F,GRMS,GMax)
      If(Print) Write(IOut,1000)
      Call Convgd(GMax,T1,IReslt)
      If(Print) Write(IOut,1002) GMax,T1,IReslt
      Call Convgd(GRMS,ConvF,IReslt)
      If(Print) Write(IOut,1004) GRMS,ConvF,IReslt
      Call Convgd(DXMaxC,T2,IReslt)
      If(Print) Write(IOut,1006) DXMaxC,T2,IReslt
      Call Convgd(DXRMSC,ConvX,IReslt)
      If(Print) Write(IOut,1008) DXRMSC,ConvX,IReslt
      Done = (DXRMSC.lt.ConvX.and.DXMaxC.lt.T2.and.
     $  GMax.lt.T1.and.GRMS.lt.ConvF)
      If(.not.Done) Call AMove(Nvar,X,Q0)
      Return
      End
