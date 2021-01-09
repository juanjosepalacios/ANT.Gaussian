!*********************************************************!
!*********************  ANT.G03-2.5.0  *******************!
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
  MODULE G03Common
!*********************************************************!
!  Common blocks for communication with Gaussian03        !
!*********************************************************!
  USE parameters, ONLY: Nalpha, Nbeta
  USE preproc, ONLY: DEFMAXSHL => MaxShl, DEFMAXATM => MaxAtm
  IMPLICIT NONE
  
  PRIVATE

  ! ***********************
  ! Gaussian Common block B 
  ! ***********************
  INTEGER MaxPrm,MaxSh1,MaxS21,JAN,ShellA,ShellN,ShellT, ShellC,ShlADF,AOS,JAnSav,NShell,MaxTyp,I5DB1,I7FB1,MaxShl
  REAL*8 EXX,C1,C2,C3,C4,X,Y,Z,RLam,RLamSv
  PARAMETER (MaxShl=DEFMAXSHL,MaxPrm=(3*MaxShl),MaxSh1=(MaxShl+1), MaxS21=(2*MaxShl+1))
  COMMON/B/EXX(MaxPrm),C1(MaxPrm),C2(MaxPrm),C3(MaxPrm),X(MaxShl), &
       &  Y(MaxShl),Z(MaxShl),JAN(MaxShl),ShellA(MaxShl),ShellN(MaxShl), &
       &  ShellT(MaxShl),ShellC(MaxShl),AOS(MaxShl),JAnSav(MaxShl), &
       &  RLam(MaxShl),RLamSv(MaxShl),NShell,MaxTyp,I5DB1,I7FB1
  DIMENSION C4(MaxShl),ShlADF(MaxShl)
  EQUIVALENCE (C4(1),C3(MaxSh1)),(ShlADF(1),C3(MaxS21))

  PUBLIC :: GetAtm4Sh, Get1stAO4Sh, GetNShell, GetShellT, GetShellC

  ! *************************
  ! Gaussian Common block Mol
  ! *************************
  INTEGER NAtoms,ICharg,Multip,NAE,NBE,NE,NBasis,IAn,NBsUse, IAtWgt,IAtTpR,IAtFrg,IAtRes,NPtMol,NPDMol,NumTpS,MolDum,IAtSpn,IAtTpS,MaxTpS,MicOpt,MaxAtm,IAtTyp
  REAL*8 AtmChg,C,AtChMM,AtmWgt,AtZEff,AtQMom,AtGFac
  PARAMETER (MaxAtm=DEFMAXATM,MaxTpS=MaxAtm)
  COMMON /Mol/ NAtoms,ICharg,Multip,NAE,NBE,NE,NBasis,IAn(MaxAtm), &
       &  NBsUse,AtmChg(MaxAtm),C(3,MaxAtm),IAtTyp(MaxAtm),AtChMM(MaxAtm), &
       &  AtmWgt(MaxAtm),IAtWgt(MaxAtm),IAtTpR(MaxAtm),IAtFrg(MaxAtm), &
       &  IAtRes(MaxAtm),NPtMol,NPDMol,NumTpS,MolDum,IAtSpn(MaxAtm), &
       &  AtZEff(MaxAtm),AtQMom(MaxAtm),AtGFac(MaxAtm),IAtTpS(2,MaxTpS), &
       &  MicOpt(MaxAtm)

  
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

    ! *******************
    ! Coordinates of atom
    ! *******************
    REAL*8 FUNCTION GetAtmCo( j, ia )
      IMPLICIT NONE
      INTEGER, INTENT(in) :: j, ia
      GetAtmCo = C(j,ia)
    END FUNCTION GetAtmCo

    END MODULE G03Common
