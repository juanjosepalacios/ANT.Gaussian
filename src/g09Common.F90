!*********************************************************!
!*********************  ANT.G-2.7.0  *******************!
!*********************************************************!
!                                                         !
!   Copyright (c) by                                      !
!                                                         !
!   Juan Jose Palacios (1)                                !
!   David Jacob (2)                                       !
!   Angel J. Perez-Jimenez (3)                            !
!   Emilio SanFabian (3)                                  !
!                                                         !
!  (1) Departamento de Fisica de la Materia Condensada    !
!      Universidad Autonoma de Madrid                     !      
!      28049 Madrid (SPAIN)                               !
!  (2) Theory Department                                  !
!      Max-Planck-Institute for Microstructure Physics    !
!      Halle, 06120 (GERMANY)                             !
!  (3) Departamento de Quimica Fisica                     !
!      Universidad de Alicante                            !
!      03690 Alicante (SPAIN)                             !
!                                                         !
!*********************************************************!
  MODULE G09Common
!**********************************************************
!  Common blocks for communication with Gaussian09        !
!**********************************************************
  USE parameters, ONLY: Nalpha, Nbeta
  USE preproc, ONLY: DEFMAXSHL => MaxShl, DEFMAXATM => MaxAtm
  IMPLICIT NONE
  
  PRIVATE

  ! ***********************
  ! Gaussian Common block B 
  ! ***********************
   Integer MaxShl,MaxPrm,MaxSh1,MaxS21,JAN,ShellA,ShellN,ShellT,ShellC,ShlADF,AOS,JAnSav,NShell,MaxTyp,I5DB1,I7FB1
      Real*8 EXX,C1,C2,C3,C4,X,Y,Z,RLam,RLamSv
      Parameter (MaxShl=DEFMAXSHL,MaxPrm=(3*MaxShl),MaxSh1=(MaxShl+1),MaxS21=(2*MaxShl+1))
      Common/B/ EXX(MaxPrm),C1(MaxPrm),C2(MaxPrm),C3(MaxPrm),X(MaxShl),&
        Y(MaxShl),Z(MaxShl),JAN(MaxShl),ShellA(MaxShl),ShellN(MaxShl),&
        ShellT(MaxShl),ShellC(MaxShl),AOS(MaxShl),JAnSav(MaxShl),&
        RLam(MaxShl),RLamSv(MaxShl),NShell,MaxTyp,I5DB1,I7FB1
      Dimension C4(MaxShl),ShlADF(MaxShl)
      Equivalence (C4(1),C3(MaxSh1)),(ShlADF(1),C3(MaxS21))

  PUBLIC :: GetAtm4Sh, Get1stAO4Sh, GetNShell, GetShellT, GetShellC, GetShellA, GetShlADF, GetShellN, GetEXX, GetC1, GetC2, GetC3, GetC4

  ! ************************
  ! Gaussian Common block B2 
  ! ************************    
  
   Integer JANB,ShelAB,SHELNB,SHELTB,SHELCB,AOSB,JAnSvB,NShelB,MaxTyB,ShlADB,I5DB2,I7FB2
      Real*8 EXXB,C1B,C2B,C3B,C4B,XB,YB,ZB,RLamB,RLmSvB
      Common/B2/EXXB(MaxPrm),C1B(MaxPrm),C2B(MaxPrm),C3B(MaxPrm),XB(MaxShl),&
      YB(MaxShl),ZB(MaxShl),JANB(MaxShl),ShelAB(MaxShl),ShelNB(MaxShl),&
      ShelTB(MaxShl),ShelCB(MaxShl),AOSB(MaxShl),JAnSvB(MaxShl),RLamB(MaxShl),&
      RLmSvB(MaxShl),NShelB,MaxTyB,I5DB2,I7FB2  
      Dimension C4B(MaxShl),ShlADB(MaxShl)
      Equivalence (C4B(1),C3B(MaxSh1)),(ShlADB(1),C3B(MaxS21))      
   
  PUBLIC :: GetC1B, GetC2B, GetC3B         

  ! *************************
  ! Gaussian Common block Mol
  ! *************************

      Integer MaxAtm,NAtoms,ICharg,Multip,NAE,NBE,NE,NBasis,IAn,NBsUse,&
        IAtWgt,IAtTpR,IAtFrg,IAtRes,NPtMol,NPDMol,NumTpS,MolDum,IAtSpn,&
        IAtTpS,MaxTpS,MicOpt,IAtTyp
      Real*8 AtmChg,C,AtChMM,AtmWgt,AtZEff,AtQMom,AtGFac
      Parameter (MaxAtm=DEFMAXATM,MaxTpS=MaxAtm)
      Common /Mol/ NAtoms,ICharg,Multip,NAE,NBE,NE,NBasis,IAn(MaxAtm),&
        NBsUse,AtmChg(MaxAtm),C(3,MaxAtm),IAtTyp(MaxAtm),AtChMM(MaxAtm),&
        AtmWgt(MaxAtm),IAtWgt(MaxAtm),IAtTpR(MaxAtm),IAtFrg(MaxAtm),&
        IAtRes(MaxAtm),NPtMol,NPDMol,NumTpS,MolDum,IAtSpn(MaxAtm),&
        AtZEff(MaxAtm),AtQMom(MaxAtm),AtGFac(MaxAtm),IAtTpS(2,MaxTpS),&
        MicOpt(MaxAtm)
  
  PUBLIC :: GetNAtoms, GetNAE, GetNBE, GetNE, GetNBasis, GetAN, GetAtmChg, GetAtmCo

  CONTAINS

    ! **********************************
    ! Get Atom number for shell number i
    ! **********************************
    INTEGER FUNCTION GetAtm4Sh( i )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i
      GetAtm4Sh = JAN(i)
    END FUNCTION GetAtm4Sh

    ! *********************************************
    ! Get number of first atomic orbital on shell i
    ! *********************************************
    INTEGER FUNCTION Get1stAO4Sh( i )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i
      Get1stAO4Sh = AOS( i )
    END FUNCTION Get1stAO4Sh

    ! ****************************
    ! Number of shells in molecule
    ! ****************************
    INTEGER FUNCTION GetNShell()
      IMPLICIT NONE
      GetNShell = NShell
    END FUNCTION GetNShell

    ! *******************************
    ! Number of primitives in shell i
    ! *******************************
    INTEGER FUNCTION GetShellN( i )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i
      GetShellN = ShellN(i)
    END FUNCTION GetShellN        

    ! **************************
    ! Get shell type for shell i
    ! **************************
    ! s = 0, p = 1, d = 2 ,f = 3 ...
    ! sp = 1 also but shell constraint is different
    INTEGER FUNCTION GetShellT( i )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i
      GetShellT = ShellT( i )
    END FUNCTION GetShellT

    ! ********************************
    ! Get shell constraint for shell i
    ! ********************************
    ! to distinguish between p and sp-shell
    ! sp = 2 , p = 1  
    INTEGER FUNCTION GetShellC( i )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i
      GetShellC = ShellC( i )
    END FUNCTION GetShellC
    
    ! ********************************
    ! Get starting location within 
    ! (EXX,C1,C2) of the data for shell i
    ! ********************************
    ! Coefficients for shells with L<=1
    ! Exponents for all shells  
    INTEGER FUNCTION GetShellA( i )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i
      GetShellA = ShellA( i )
    END FUNCTION GetShellA
    
    ! ********************************
    ! Get exponent of primitive j  
    ! in shell i
    ! ********************************
    ! Exponents for all shell types  
    REAL*8 FUNCTION GetEXX( i )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i
      GetEXX = EXX( i )
    END FUNCTION GetEXX    
    
    ! ********************************
    ! Get coefficient of primitive j  
    ! in shell i
    ! ********************************
    ! Coefficient for S shells  
    REAL*8 FUNCTION GetC1( i )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i
      GetC1 = C1( i )
    END FUNCTION GetC1 

    ! ********************************
    ! Get coefficient of primitive j  
    ! in shell i
    ! ********************************
    ! Coefficient for S shells  
    REAL*8 FUNCTION GetC1B( i )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i
      GetC1B = C1B( i )
    END FUNCTION GetC1B         

    ! ********************************
    ! Get coefficient of primitive j  
    ! in shell i
    ! ********************************
    ! Coefficient for P shells  
    REAL*8 FUNCTION GetC2( i )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i
      GetC2 = C2( i )
    END FUNCTION GetC2  
    
    ! ********************************
    ! Get coefficient of primitive j  
    ! in shell i
    ! ********************************
    ! Coefficient for P shells  
    REAL*8 FUNCTION GetC2B( i )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i
      GetC2B = C2B( i )
    END FUNCTION GetC2B                               
    
    ! ********************************
    ! Get starting location within 
    ! (C3,C4) of the data for shell i
    ! ********************************
    ! For shells with L>=2  
    INTEGER FUNCTION GetShlADF( i )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i
      GetShlADF = ShlADF( i )
    END FUNCTION GetShlADF

    ! ********************************
    ! Get coefficient of primitive j  
    ! in shell i
    ! ********************************
    ! Coefficient for D shells  
    REAL*8 FUNCTION GetC3( i )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i
      GetC3 = C3( i )
    END FUNCTION GetC3   

    ! ********************************
    ! Get coefficient of primitive j  
    ! in shell i
    ! ********************************
    ! Coefficient for F and higher shells  
    REAL*8 FUNCTION GetC4( i )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i
      GetC4 = C4( i )
    END FUNCTION GetC4   
        
    ! ********************************
    ! Get coefficient of primitive j  
    ! in shell i
    ! ********************************
    ! Coefficient for D shells  
    REAL*8 FUNCTION GetC3B( i )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: i
      GetC3B = C3B( i )
    END FUNCTION GetC3B                                      

    ! ***************************
    ! Number of atoms in molecule
    ! ***************************
    INTEGER FUNCTION GetNAtoms()
      IMPLICIT NONE
      GetNAtoms = NAtoms
    END FUNCTION GetNAtoms

    ! *************************
    ! Number of Alpha electrons
    ! *************************
    INTEGER FUNCTION GetNAE()
      IMPLICIT NONE
      if (Nalpha  < 0) then
      GetNAE = NAE
      else
      GetNAE = Nalpha
      end if
    END FUNCTION GetNAE
    
    ! ************************
    ! Number of Beta electrons
    ! ************************
    INTEGER FUNCTION GetNBE()
      IMPLICIT NONE
      if (Nbeta  < 0) then
      GetNBE = NBE
      else
      GetNBE = Nbeta
      end if
    END FUNCTION GetNBE

    ! *******************
    ! Number of electrons
    ! *******************
    INTEGER FUNCTION GetNE()
      IMPLICIT NONE
      if (Nbeta  >= 0 .and. Nalpha >= 0) then
      GetNE = Nalpha+Nbeta
      else
      GetNE = NE
      end if
    END FUNCTION GetNE

    ! *************************
    ! Number of basis functions
    ! *************************
    INTEGER FUNCTION GetNBasis()
      IMPLICIT NONE
      GetNBasis = NBasis
    END FUNCTION GetNBasis

    ! ************************
    ! Atomic number of atom ia
    ! ************************
    INTEGER FUNCTION GetAN( ia )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ia
      GetAN = IAN(ia)
    END FUNCTION GetAN

    ! ****************************************
    ! Atom core charge = nuke + core electrons
    ! ****************************************
    REAL*8 FUNCTION GetAtmChg( ia )
      IMPLICIT NONE 
      INTEGER, INTENT(in) :: ia
      GetAtmChg = AtmChg(ia)
    END FUNCTION GetAtmChg

    ! **********************
    ! Coordinates of atom ia
    ! **********************
    REAL*8 FUNCTION GetAtmCo( j, ia )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: j, ia
      GetAtmCo = C(j,ia)
    END FUNCTION GetAtmCo

  END MODULE G09Common
