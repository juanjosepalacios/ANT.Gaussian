!*********************************************************!
!*********************  ANT.G-2.5.2  *********************!
!*********************************************************!
!                                                         !
!   Copyright (c) by                                      !
!                                                         !
!   Juan Jose Palacios (1)                                !
!   David Jacob (2)                                       !
!                                                         !
!  (1) Departamento de Fisica de la Materia Condensada    !
!      Universidad Autonoma de Madrid                     !      
!      28049 Madrid (SPAIN)                               !
!  (2) Theory Department                                  !
!      Max-Planck-Institute for Microstructure Physics    !
!      Halle, 06120 (GERMANY)                             !
!                                                         !
!*********************************************************!
  MODULE constants
!*********************************************************!
!  Module containing mathematical and physical constants  !
!*********************************************************!

  IMPLICIT NONE
  
  ! Precision of numerical types
  integer, parameter :: single = 4
  integer, parameter :: double = 8  

  ! Mathematical constants
  REAL*8, PARAMETER :: d_pi   = 3.14159265358979323846d0
  REAL*8, PARAMETER :: d_zero = 0.0d0
  REAL*8, PARAMETER :: d_one  = 1.0d0
  
  COMPLEX*16, PARAMETER :: c_zero = (0.0d0,0.0d0)
  COMPLEX*16, PARAMETER :: c_one  = (1.0d0,0.0d0)
  COMPLEX*16, PARAMETER :: ui   = (0.0d0,1.0d0)

  ! Physical constants
  REAL*8, PARAMETER :: Hart = 27.2113834d0
  REAL*8, PARAMETER :: Ryd  = 0.5d0*Hart
  REAL*8, PARAMETER :: Bohr  = 0.5291772108d0
  
  REAL*8, PARAMETER :: eleccharge = 1.6d-19
!  REAL*8, PARAMETER :: hbar = 1.05457d-34 ! J*s
  REAL*8, PARAMETER :: hbar = 6.582d-16 ! eV*S ! WITH THIS VALUE THE GibbsY OPERATOR RESULTS HERMITIAN.  
  
  character(LEN=4), dimension(2), parameter :: SpinStr = (/"UP","DOWN"/)

  ! Array limits
  integer, parameter :: MaxAODOS    = 1000
  integer, parameter :: MaxCorrBl   = 100
  integer, parameter :: MaxDAtoms   = 1000
  integer, parameter :: MaxPCAtoms  = 100
  integer, parameter :: MaxLAO      = 1500
  integer, parameter :: MaxAtomData = 100
  integer, parameter :: MaxAOAtom   = 100
  integer, parameter :: MaxShAtom   = 20  

  ! Computational constants
  
  END MODULE constants
