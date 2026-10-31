!*********************************************************!
!*********************  ANT.G-2.8.4  *********************!
!*********************************************************!
!                                                         !
!  Copyright (c) by                                       !
!                                                         !
!  Juan Jose Palacios (1)                                 !
!  David Jacob (2)                                        !
!  Wynand Dednam (5)                                      !
!                                                         !
! (1) Departamento de Fisica de la Materia Condensada     !
!     Universidad Autonoma de Madrid                      !    
!     28049 Madrid (SPAIN)                                !
! (2) Theory Department                                   !
!     Max-Planck-Institute for Microstructure Physics     !
!     Halle, 06120 (GERMANY)                              !
! (5) Departamento de Fisica Aplicada                     !
!     Universidad de Alicante                             !    
!     03690 Alicante (SPAIN)                              !
!                                                         !
!*                                     *
!*  ANT1D - OneDLead.f90               *
!*********************************************************!
!                                                         !
!   Module for description of one-dimensional lead        !
!                                                         !
!*********************************************************!
  module OneDLead
!*********************************************************!
  use constants
  use parameters
  implicit none
  save
  
  private

  !*********** 
  !Public Type
  !***********
  public :: T1DLead

  !****************
  !Public variables
  !****************
  public :: Lead1D

  !************************ 
  !Public module procedures
  !************************
  public :: Init1DLead
  public :: ReadHamiltonian
  public :: CleanUp1DLead
  public :: CompSelfEnergy1D
  public :: L1D_NSpin, L1D_NAOrbs, L1D_NElectrons
  public :: L1D_EMin, L1D_EMax, L1D_EFermi
  public :: L1D_H0, L1D_V1, L1D_V2, L1D_S0, L1D_S1, L1D_S2
 !public :: L1D_NAtoms, L1D_Atom
  public :: PrintMatrices,PrintDOS

  !*******************************
  !One-dimensional lead parameters
  !*******************************
  type T1DLead
!    private

     !Lead number
     integer :: LeadNo
     
     ! Number of atoms per PUC
     integer :: NPCAtoms
     ! Atom data indices
    !type(TAtom), dimension(MaxPCAtoms) :: PCAtoms

     !Number of non-degenerate Spin bands
     integer :: NSpin    
     !Number of atomic orbitals in supercell
     integer :: NAOrbs
     !Spin-Orbit coupling: 0 = off, 1 = on
     integer :: SO
     !Number of electrons in supercell
     integer :: NElectrons

     !Lower and upper bound for non-zero DOS
     real(double) :: EMin = -100.0, EMax = 1000.0
     real(double) :: EFermi = -5.0

     !Number of primitive unit cells
     integer :: NPC
     !Number of atomic orbitals in PRIMITIVE unit cell
     integer :: NPCAO
     !Number of electrons in PRIMITIVE unit cell     
     integer :: NPCEl

     integer, dimension(MaxLAO) :: AOT

     !Coupling between n-th neighbour PUCs
     complex(double), dimension(:,:,:,:,:), pointer :: vn
     !Overlap between n-th neighbour PUCs 
     complex(double), dimension(:,:,:), pointer :: sn     
  end type T1DLead

  type(T1DLead),dimension(2) :: Lead1D

  !Coupling between n-th neighbour PUCs
  complex(double), dimension(:,:,:,:,:),allocatable, target :: vn1, vn2
  !Overlap between n-th neighbour PUCs 
  complex(double), dimension(:,:,:),allocatable, target :: sn1, sn2
  

  !******************
  !Internal variables
  !******************
 
  real(double)    :: ChargeOffset 
  integer :: WhichLead
  integer :: FermiEvaluations=0
  ! Independent guesses in the unshifted Hamiltonian reference, per electrode.
  ! Preserve across CleanUp1DLead during SCF; never retain old charge values.
  real(double) :: WarmLeadFermi(2)=0d0
  logical :: WarmLeadFermiValid(2)=.false.
  ! Initial search half-width (eV), a heuristic rather than a convergence tolerance.
  real(double), parameter :: WarmFermiHalfWidth=0.25d0
 
  !eXtended supercell matrices 
  integer :: NXAO, NXDim, NXSpin
  complex(double), dimension(:,:,:), allocatable :: HX
  complex(double), dimension(:,:),   allocatable :: GX, SX
  complex(double), dimension(:,:),   allocatable :: Sigma1, Sigma2

  ! whether to create a device input file from extended supercell
  logical :: MakeDevice = .false., init=.true.
  integer :: NDevBeg=1,NDevEnd=0

contains

  !***************************
  !*** Initialize 1D Leads ***
  !***************************
  subroutine Init1DLead (ilead,NSpin,HD,SD,NAOrbs)
    use parameters
!   use device, only: LeadsOn
    implicit none
    integer :: inpf,ios
    integer, intent(in) :: ilead,NSpin,NAOrbs
    real(double), dimension(:,:,:), intent(in) :: HD
    real(double), dimension(:,:), intent(in) :: SD
    integer :: ixyzf, spin_index

    print *, "*******************************"
    print *, "*    Initializing Lead No. ", ilead
    print *, "*******************************"
 
    call ReadHamiltonian (ilead,Nspin,HD,SD,NAOrbs)

    Lead1D(ilead)%NAOrbs     = Lead1D(ilead)%NPC * Lead1D(ilead)%NPCAO
    Lead1D(ilead)%NElectrons = Lead1D(ilead)%NPC * Lead1D(ilead)%NPCEl
    Lead1D(ilead)%SO = 0
    Lead1D(ilead)%LeadNo = ilead

   !call PrintMatrices( ilead )
   !
   ! Extend supercell by one or more primitive cells to each side 
   ! in order to compute charge and DOS in the case of 
   ! non-orthogonal basis set
   !
    ! Seed the search independently from this lead's repeated supercell.
    call SetInitialLeadLowerBound(Lead1D(ilead))
    call CreateXSC( Lead1D(ilead) )
   !
   ! If no Energy boundaries have been defined in input find them
   !
    call FindEBounds( Lead1D(ilead) )
   !
   ! Find Fermi level of lead if requested
   !
    call FindFermi( Lead1D(ilead) )
   !
   ! Shift Hamiltonian to adjust Fermi level to zero
   !
    call AdjustFermi( ilead ) 
   !
   ! Print out DOS of bulk Lead if requested
   !
     if( LeadDOS ) then
       ! The self-energies already use the Fermi-shifted lead Hamiltonian.
       ! Shift the temporary bulk Hamiltonian to the same energy reference.
       do spin_index=1,NXSpin
          HX(:,:,spin_index) = HX(:,:,spin_index) - Lead1D(ilead)%EFermi*SX
       end do
      !call CreateXSC( Lead1D(ilead) )
      !call PrintMatrices (ilead)
       call PrintDOS( ilead )
     end if 
    !
    ! Deallocate all extended supercell matrices
    !
    deallocate( HX, SX, GX, Sigma1, Sigma2 )
    !

  end subroutine Init1DLead


  !*****************************************
  !*** Subroutine for Cleaning up module ***
  !*****************************************
  subroutine CleanUp1DLead(ilead)
    integer, intent(in) :: ilead      
    if (ilead == 1) deallocate( vn1, sn1 )
    if (ilead == 2) deallocate( vn2, sn2 )
  end subroutine CleanUp1DLead


  !**********************************************************
  !*** Compute self-energy matrices projected into device ***
  !**********************************************************
  recursive subroutine CompSelfEnergy1D( L1D, spin, z, Sigma, NBasis, idir )
    use parameters, only: eta
    use constants, only: ui, c_zero
    implicit none

    type(T1DLead), intent(in) :: L1D
    integer, intent(in) :: spin
    complex(double), intent(in) :: z
   !complex(double), dimension( NBasis*(L1D%SO+1), NBasis*(L1D%SO+1) ), intent(out) :: Sigma
    complex(double), dimension(:,:), intent(out) :: Sigma
    integer , intent(in) :: idir, NBasis !! AO1, 

    integer :: i,j,AOrb1,AOrb2, ispin, jspin, NAO, ierr,solver_status

   !complex(double), dimension(:,:), allocatable :: Veff1, Heff, Veff2
    complex(double), dimension(L1D%NAOrbs,L1D%NAOrbs) :: Veff1, Heff, Veff2
    complex(double), dimension(L1D%NAOrbs,L1D%NAOrbs) :: sig


    NAO = L1D%NAOrbs
   !allocate( Veff1(NAO,NAO), Heff(NAO,NAO), Veff2(NAO,NAO),  STAT=ierr )
   !if( ierr /= 0 ) then
   !   print *, "OneDLead/CompSelfEnergy1D/ALLOCATION ERROR: Veff1, Heff, Veff2"
   !   call AbortProg
   !end if
    Veff1 = c_zero
    Veff2 = c_zero
    Heff = c_zero
    sig = c_zero
    Sigma = c_zero

    ispin=spin
    if( L1D%NSpin==1 .and. spin == 2 ) ispin = 1

    do i=1,NAO
       do j=1,NAO
          Veff1(i,j) =  L1D_V1(L1D,ispin,i,ispin,j) - (z+ui*eta)*L1D_S1(L1D,i,j)
          Veff2(i,j) =  L1D_V2(L1D,ispin,i,ispin,j) - (z+ui*eta)*L1D_S2(L1D,i,j)
          Heff (i,j) = -L1D_H0(L1D,ispin,i,ispin,j) + (z+ui*eta)*L1D_S0(L1D,i,j) 
       end do
    end do
    if( idir == 1 )then
       call SolveDyson1D(sig,Veff1,Heff,Veff2,NAO,status=solver_status)
       AOrb1 = 1
       AOrb2 = NAO
    else
       call SolveDyson1D(sig,Veff2,Heff,Veff1,NAO,status=solver_status)
       AOrb1=NBasis-NAO+1
       AOrb2=NBasis
    end if

    if (solver_status/=0) then
       print *, 'ERROR: 1D Dyson failed, lead/spin/direction/status=',L1D%LeadNo,spin,idir,solver_status
       print *, 'Energy (eV)=',z,' mixing=',L1DAlpha
       stop 1
    end if
    do i=AOrb1,AOrb2
       do j=AOrb1,AOrb2
          Sigma( i, j ) = sig( i-AOrb1+1, j-AOrb1+1 )
       end do
    end do
   !deallocate( Veff1, Heff, Veff2, STAT=ierr )
   !if( ierr /= 0 ) then
   !   print *, "OneDLead/CompSelfEnergy1D/DEALLOCATION ERROR: Veff1, Heff, Veff2"
   !   call AbortProg
   !end if
  end subroutine CompSelfEnergy1D
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*** Access functions to T1DLead members ***
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  !**********************************************
  !*** Number of atoms in supercell ***
  !**********************************************
  integer function L1D_NAOrbs( L1D )
    type(T1DLead),intent(IN) :: L1D
    L1D_NAOrbs = L1D%NAOrbs
  end function L1D_NAOrbs

  !**********************************************
  !*** Number of atoms in supercell ***
  !**********************************************
 !integer function L1D_NAtoms( L1D )
 !  type(T1DLead),intent(IN) :: L1D
 !  L1D_NAtoms = L1D%NPCAtoms*L1D%NPC
 !end function L1D_NAtoms

  !**********************************************
  !*** Number of atoms in supercell ***
  !**********************************************
 !type(TAtom) function L1D_Atom( L1D, iatom )
 !  type(T1DLead),intent(IN) :: L1D
 !  integer, intent(in) :: iatom
 !  L1D_Atom = L1D%PCAtoms( mod(iatom-1,L1D%NPCAtoms)+1 )
 !end function L1D_Atom

  !**********************************************
  !*** Number of non-degenerate spin-channels ***
  !**********************************************
  integer function L1D_NSpin( L1D )
    type(T1DLead),intent(IN) :: L1D
    L1D_NSpin = L1D%NSpin
  end function L1D_NSpin
  
  !****************************************
  !*** Number of electrons in supercell ***
  !****************************************
  integer function L1D_NElectrons( L1D )
    type(T1DLead),intent(IN) :: L1D
    L1D_NElectrons = L1D%NElectrons
  end function L1D_NElectrons

  !*****************************
  !*** Lower energy boundary ***
  !*****************************
  !for printing DOS and finding EF
  real(double) function L1D_EMin( L1D )
    type(T1DLead),intent(IN) :: L1D
    L1D_EMin = L1D%EMin
  end function L1D_EMin

  !*****************************
  !*** Upper energy boundary ***
  !*****************************
  !for printing DOS and finding EF
  real(double) function L1D_EMax( L1D )
    type(T1DLead),intent(IN) :: L1D
    L1D_EMax = L1d%EMax
  end function L1D_EMax

  !******************************
  !*** Fermi energy of a Lead ***
  !******************************
  real(double) function L1D_EFermi( L1D )
    type(T1DLead),intent(IN) :: L1D
    L1D_EFermi = L1d%EFermi
  end function L1D_EFERMI

  !**********************************
  !*** Lead Supercell Hamiltonian ***
  !**********************************
  recursive complex(double) function L1D_H0( L1D, ispin, i, jspin, j )
    type(T1DLead),intent(IN) :: L1D
    integer, intent(IN) ::  ispin, i, j
    integer :: ipc, jpc, iao, jao, jspin
    ! PC indices
    ipc = (i-1)/L1D%NPCAO + 1
    jpc = (j-1)/L1D%NPCAO + 1
    ! AO indices
    iao = mod( i-1, L1D%NPCAO ) + 1
    jao = mod( j-1, L1D%NPCAO ) + 1
    if( ipc >= jpc )then
       L1D_H0 = conjg( L1D%vn( jspin, ispin, ipc-jpc, jao, iao ) )
    else
       L1D_H0 = L1D%vn( ispin, jspin, jpc-ipc, iao, jao )
    end if
  end function L1D_H0

  !********************************************
  !*** Lead Supercell Hopping (to the left) ***
  !********************************************
  recursive complex(double) function L1D_V1( L1D, ispin, i, jspin, j )
    type(T1DLead),intent(IN) :: L1D
    integer, intent(IN) ::  ispin, i, j, jspin
    L1D_V1 = conjg( L1D_V2( L1D, jspin, j, ispin, i ) )
  end function L1D_V1

  !*********************************************
  !*** Lead Supercell Hopping (to the right) ***
  !*********************************************
  recursive complex(double) function L1D_V2( L1D, ispin, i, jspin, j )
    type(T1DLead),intent(IN) :: L1D
    integer, intent(IN) ::  ispin, i, j, jspin
    integer :: ipc, jpc, iao, jao
    ! PC indices
    ipc = (i-1)/L1D%NPCAO + 1
    jpc = (j-1)/L1D%NPCAO + 1
    ! AO indices
    iao = mod( i-1, L1D%NPCAO ) + 1
    jao = mod( j-1, L1D%NPCAO ) + 1
    L1D_V2 = c_zero
    if( ipc < jpc ) return
    L1D_V2 = L1D%vn( ispin, jspin, L1D%NPC - ipc + jpc, iao, jao )
  end function L1D_V2

  !******************************
  !*** Lead supercell overlap ***
  !******************************
  recursive complex(double) function L1D_S0( L1D, i, j )
    type(T1DLead),intent(IN) :: L1D
    integer, intent(IN) :: i, j
    integer :: ipc, jpc, iao, jao
    ! PC indices
    ipc = (i-1)/L1D%NPCAO + 1
    jpc = (j-1)/L1D%NPCAO + 1
    ! AO indices
    iao = mod( i-1, L1D%NPCAO ) + 1
    jao = mod( j-1, L1D%NPCAO ) + 1
    if( ipc >= jpc )then
       L1D_S0 = conjg( L1D%sn( ipc-jpc, jao, iao ) )
    else
       L1D_S0 = L1D%sn( jpc-ipc, iao, jao )
    end if
  end function L1D_S0

  !********************************
  ! Lead supercell overlap to left 
  !********************************
  recursive complex(double) function L1D_S1( L1D, i, j )
    type(T1DLead),intent(IN) :: L1D
    integer, intent(IN) :: i, j
    L1D_S1 = conjg( L1D_S2( L1D, j, i ) )
  end function L1D_S1

  !*********************************
  ! Lead supercell overlap to right 
  !*********************************
  recursive complex(double) function L1D_S2( L1D, i, j )
    type(T1DLead),intent(IN) :: L1D
    integer, intent(IN) :: i, j
    integer :: ipc, jpc, iao, jao
    ! PC indices
    ipc = (i-1)/L1D%NPCAO + 1
    jpc = (j-1)/L1D%NPCAO + 1
    ! AO indices
    iao = mod( i-1, L1D%NPCAO ) + 1
    jao = mod( j-1, L1D%NPCAO ) + 1
    L1D_S2 = c_zero
    if( ipc < jpc ) return
    L1D_S2 = L1D%sn( L1D%NPC - ipc + jpc, iao, jao )
  end function L1D_S2


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!! Private Subroutines and Functions !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
  !********************************
  ! Read lead parameters from ANT 
  !********************************
  subroutine ReadHamiltonian(ilead,NSpin,HD,SD,NAOrbs )
!   use device, only: DevNSpin, DevFockMat, DevOverlapMat, DevNAOrbs, LeadsOn  
    use parameters, only: NatPC, NPC
    use cluster, only: LoAOrbNo,  HiAOrbNo, LeadAtmNo, LeadNAOrbs, NALead,  NAOAtom
#ifdef G03ROOT
    use g03Common, only: GetNAtoms, GetAtmChg, GetAN, GetAtmCo
#endif
#ifdef G09ROOT
    use g09Common, only: GetNAtoms, GetAtmChg, GetAN, GetAtmCo
#endif  
    use numeric, only: CSetId
    use util
    use messages    
    implicit none
    
    real(double), dimension(:,:,:), intent(in) :: HD
    real(double), dimension(:,:), intent(in) :: SD
    
    ! Number of non-degenerate Spin bands
    integer, intent(in) :: NSpin
    integer, intent(in) :: ilead
    integer, intent(in) :: NAOrbs

    ! Number of atoms per PUC
    integer :: NPCAtoms
    ! Number of atomic orbitals in PRIMITIVE unit cell
    integer :: NPCAO
    ! Number of electrons in PRIMITIVE unit cell     
    integer :: NPCEl
 
    integer :: ierr, ispin, ipc, iatom, jatom, AN, i, j, k, firstao, lastao
    character(len=10) :: ANStr
    logical :: initzero
    
    Lead1D(ilead)%LeadNo   = ilead
    Lead1D(ilead)%NSpin    = NSpin
    Lead1D(ilead)%NPC      = NPC
    Lead1D(ilead)%NPCAtoms  = NatPC(ilead)
    NPCAtoms = Lead1D(ilead)%NPCAtoms
    if (ilead ==1 ) then
      Lead1D(ilead)%NPCEl    = NatPC(ilead)*GetAtmChg(1)
      Lead1D(ilead)%NPCAO    = NatPC(ilead)*NAOAtom(1)
    else 
      Lead1D(ilead)%NPCEl    = NatPC(ilead)*GetAtmChg(GetNAtoms())
      Lead1D(ilead)%NPCAO    = NatPC(ilead)*NAOAtom(GetNatoms())
    end if      
    NPCEl = Lead1D(ilead)%NPCEl
    NPCAO = Lead1D(ilead)%NPCAO

    !
    ! Allocate Hoppings and overlaps
    !
    if( ilead == 1 )then
       allocate( vn1(2,2,0:NPC,NPCAO,NPCAO),sn1(0:NPC,NPCAO,NPCAO), stat=ierr)
       if( ierr /= 0 )call AllocErr("OneDLead/ReadParameter","vn1,sn1")
       Lead1D(1)%vn => vn1
       Lead1D(1)%sn => sn1
       vn1 = 0.0d0
       sn1 = 0.0d0
       do ispin=1,NSpin
        do ipc=0,NPC
            do i=1,NPCAO
              do j=1,NPCAO
               vn1(ispin,ispin,ipc,i,j)=HD(ispin,i ,j+ipc*NPCAO )
               sn1(ipc,i,j)=SD(i ,j+ipc*NPCAO )
              end do  
          end do  
        end do
     end do          
    endif
     
    if( ilead == 2 )then
       allocate( vn2(2,2,0:NPC,NPCAO,NPCAO),sn2(0:NPC,NPCAO,NPCAO), stat=ierr)
       if( ierr /= 0 )call AllocErr("OneDLead/ReadParameter","vn2,sn2")
       Lead1D(2)%vn => vn2
       Lead1D(2)%sn => sn2
       vn2 = 0.0d0
       sn2 = 0.0d0
       do ispin=1,NSpin
        do ipc=0,NPC
            do i=1,NPCAO
              do j=1,NPCAO
               vn2(ispin,ispin,ipc,i,j)=HD(ispin,i+(NAOrbs-NPCAO)-ipc*NPCAO ,j+(NAOrbs-NPCAO) )
               sn2(ipc,i,j)=SD(i+(NAOrbs-NPCAO)-ipc*NPCAO ,j+(NAOrbs-NPCAO))
   !           vn2(ispin,ispin,ipc,i,j)=HD(ispin,i ,j+ipc*NPCAO )
   !           sn2(ipc,i,j)=SD(i ,j+ipc*NPCAO )
              end do  
          end do  
        end do
     end do          
    endif
    
  end subroutine ReadHamiltonian  


  !********************************************
  ! Print out Hamiltonian and Overlap Matrices
  !********************************************
  subroutine PrintMatrices( ilead )
    implicit none

    integer, intent(in) :: ilead
    integer :: NPC, NPCAO, NAO, NSpin, i, j, ispin, proc_id

    NAO   = Lead1D(ilead)%NAOrbs
    NSpin = Lead1D(ilead)%NSpin
   !if(proc_id == 0)then
       print *
       print *, "*** Supercell on-site and coupling matrices ***"
       print *
       do ispin=1,NSpin
          if( NSpin == 1 ) print *, "H0 = "
          if( NSpin == 2 .and. ispin == 1 ) print *, "H0 up = "
          if( NSpin == 2 .and. ispin == 2 ) print *, "H0 down = "       
          do i=1,NAO
            print '(1000(E14.4))', ( DREAL( L1D_H0( Lead1D(ilead), ispin, i, ispin, j)), j=1,NAO )
          end do
          if( NSpin == 1 ) print *, "V1 = "
          if( NSpin == 2 .and. ispin == 1 ) print *, "V1 up = "
          if( NSpin == 2 .and. ispin == 2 ) print *, "V1 down = "       
          do i=1,NAO
             print '(1000(E14.4))', ( DREAL( L1D_V1( Lead1D(ilead), ispin, i, ispin, j)), j=1,NAO )
          end do
       end do
       print *, "S0 = "
       do i=1,NAO
          print '(1000(E14.4))', ( DREAL(L1D_S0( Lead1D(ilead), i,j)), j=1,NAO )
       end do
       print *, "S1 = "
       do i=1,NAO
          print '(1000(E14.4))', (DREAL(L1D_S1( Lead1D(ilead), i,j)), j=1,NAO )
       end do
       call FlUSH(6)
   !end if
  end subroutine PrintMatrices

  !***************************
  ! Create extended supercell
  !***************************
  subroutine CreateXSC( L1D )
    use parameters
    use constants, only: c_zero
    implicit none
    
    type(T1DLead) :: L1D
    
    integer :: ierr, ispin, m, n, i, j
    integer :: NPCAO

!   NPC = L1D%NPC
    NPCAO = L1D%NPCAO
    !
    ! Allocate extended supercell matrices
    !
    ! XSC = SC + PC + SC
    NXAO   = L1D%NAOrbs + L1D%NPCAO + L1D%NAOrbs
    NXSpin = L1D%NSpin
    NXDim  = NXAO
   !print*,'NXDim=', NXDim

    allocate( &
         HX(NXDim,NXDim,NXSpin), SX(NXDim,NXDim), &
         GX(NXDim,NXDim), Sigma1(NXDim,NXDim), &
         Sigma2(NXDim,NXDim), STAT=ierr )
    if( ierr /= 0 )then
       print *, "OneDLead/Warning: HX, SX, GX, Sigma1, Sigma2 could not be allocated."
      !call AbortProg
    end if

    HX = c_zero
    GX = c_zero
    SX = c_zero
    Sigma1 = c_zero
    Sigma2 = c_zero    
    !
    ! Compute Hamiltonian HX and overlap SX
    !
    ! Loop over all primitive cells 
    do m=1,2*NPC+1
       !
       ! Diagonal blocks 
       !
       do ispin=1,NXSpin
          HX( (m-1)*NPCAO+1:m*NPCAO, (m-1)*NPCAO+1:m*NPCAO, ispin ) = L1D%vn( ispin, ispin, 0,:,:)
       end do
       SX( (m-1)*NPCAO+1:m*NPCAO, (m-1)*NPCAO+1:m*NPCAO ) = L1D%sn( 0,:,:)
       !
       !Off diagonal blocks:
       !
       ! 1) hoppings with PCs to the right of PC m: n > m
       do n=m+1,2*NPC+1
          if( n-m <= NPC )then
             do ispin = 1, NXSpin
                HX( (m-1)*NPCAO+1:m*NPCAO, (n-1)*NPCAO+1:n*NPCAO, ispin ) = L1D%vn( ispin, ispin, n-m,:,:)
             end do
             SX((m-1)*NPCAO+1:m*NPCAO, (n-1)*NPCAO+1:n*NPCAO ) = L1D%sn(n-m,:,:)
          end if
       end do
       ! 2) hoppings with PCs to the left of PC m: n < m
       do n=1,m-1
          if( m-n <= NPC )then
             do ispin=1,NXSpin 
                HX( (m-1)*NPCAO+1:m*NPCAO, (n-1)*NPCAO+1:n*NPCAO, ispin ) = conjg(transpose(L1D%vn(  ispin, ispin, m-n,:,:)))
             end do
             SX((m-1)*NPCAO+1:m*NPCAO, (n-1)*NPCAO+1:n*NPCAO ) = conjg(transpose(L1D%sn(m-n,:,:)))
          end if
       end do
    end do
       
   !if( PrintHS == 1 )then       
   !   print *
   !   print *, "*** Extended Supercell Matrices ***"
   !   print *
   !   do ispin=1,NXSpin
   !      if( NXSpin == 2 .and. ispin == 1) print *, "HX spin-up = "
   !      if( NXSpin == 2 .and. ispin == 2) print *, "HX spin-down = "
   !      if( NXSpin == 1) print *, "HX = "
   !      do i=1,NXDim
   !         print '(1000(ES14.4))', ( DREAL(HX( i, j, ispin )), j=1,NXDim )
   !         write(3000,'(1000(ES14.4))'), ( DREAL(HX( i, j, ispin )), j=1,NXDim )
   !      end do
   !   end do
   !   print *, "SX = "
   !   do i=1,NXDim
   !      print '(1000(ES14.4))', ( DREAL(SX(i,j)), j=1,NXDim )
   !      write(3001, '(3001(ES14.4))'), ( DREAL(SX(i,j)), j=1,NXDim )
   !   end do
   !   call FLUSH(6)
   !end if

  end subroutine CreateXSC

  !****************************************
  ! Write Bulk DOS of an electrode to file 
  !****************************************
  subroutine PrintDOS( ilead )
    use parameters, only: gamma, DE => EStep
    use constants, only: ui
    use filemaster
    implicit none

    integer, intent(in) :: ilead
    
    real(double) :: EMin, EMax, energy
    complex(double) :: zenergy
    character(len=12),parameter :: file(2) = (/"Lead1DOS.dat","Lead2DOS.dat"/)
    integer :: ispin, ios, iunit
    
    print '(A,I1,A,A)', "Printing DOS of Lead ", ilead, "to file ", file(ilead)
    
   !EMin = Lead1D(ilead)%EMin
   !EMax = Lead1D(ilead)%EMax
    EMin = -10.0
    EMax = 10.0

    iunit=fopen(file(ilead),'unknown', ios)
    do energy=EMin, EMax, DE 
       zenergy=energy+gamma*ui
      !zenergy=energy+gamma*ui
      !zenergy=energy
       write(UNIT=iunit,FMT='(F10.5,(1000(ES14.6)))'), energy,&
            ( BulkSDOS( Lead1D(ilead), ispin, zenergy ), ispin=1,NXSpin )
       call flush(iunit)
    end do
    call fclose(iunit)
  end subroutine PrintDOS

  !*******************************************************
  !*** Iterative Dyson solver for 1D lead with overlap ***
  !*******************************************************
  !*                          1                          *
  !*       Sigma = Veff1 ------------ Veff2              *
  !*                     Heff - Sigma                    *
  !*******************************************************
  recursive subroutine SolveDyson1D(Sigma,Veff1,Heff,Veff2,NDim,iterations,residual,status)
    use parameters, only: conv=>SelfAcc,alpha=>L1DAlpha,MaxCycle=>L1DMaxCyc
    use constants, only: c_zero,c_one,ui
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifdef PGI
    use lapack_blas, only: zgetrf,zgetri,zgemm
#endif
    implicit none
    integer, intent(in) :: NDim
    complex(double), intent(out) :: Sigma(NDim,NDim)
    complex(double), intent(in) :: Veff1(NDim,NDim),Heff(NDim,NDim),Veff2(NDim,NDim)
    integer, intent(out), optional :: iterations,status
    real(double), intent(out), optional :: residual
    integer :: ipiv(NDim),i,j,info,ncycle,solve_status
    real(double) :: error
    complex(double) :: Sigma_aux(NDim,NDim),temp(NDim,NDim),work(4*NDim)

    Sigma=c_zero
    ncycle=0
    error=huge(0d0)
    solve_status=4
    if (NDim<1 .or. MaxCycle<1) goto 900
    if (.not.ieee_is_finite(alpha) .or. alpha<=0d0 .or. alpha>1d0) goto 900
    if (.not.ieee_is_finite(conv) .or. conv<=0d0) goto 900
    do i=1,NDim
       Sigma(i,i)=-ui
    end do
    solve_status=1
    do ncycle=1,MaxCycle
       Sigma_aux=Heff-Sigma
       call zgetrf(NDim,NDim,Sigma_aux,NDim,ipiv,info)
       if (info/=0) then
          solve_status=2
          exit
       end if
       call zgetri(NDim,Sigma_aux,NDim,ipiv,work,4*NDim,info)
       if (info/=0) then
          solve_status=2
          exit
       end if
       call zgemm('N','N',NDim,NDim,NDim,c_one,Sigma_aux,NDim,Veff2,NDim,c_zero,temp,NDim)
       call zgemm('N','N',NDim,NDim,NDim,c_one,Veff1,NDim,temp,NDim,c_zero,Sigma_aux,NDim)
       ! True fixed-point residual, before mixing and reset every iteration.
       ! It must not be scaled by (1-alpha), which falsely favors large alpha.
       error=0d0
       do i=1,NDim
          do j=1,NDim
             error=error+abs(Sigma_aux(i,j)-Sigma(i,j))
          end do
       end do
       error=error/dble(NDim*NDim)
       if (.not.ieee_is_finite(error)) then
          solve_status=3
          exit
       end if
       if (error<=conv) then
          solve_status=0
          exit
       end if
       ! Leave Sigma at the iterate whose residual is reported on exhaustion.
       if (ncycle==MaxCycle) exit
       Sigma=alpha*Sigma_aux+(1d0-alpha)*Sigma
    end do
900 continue
    if (present(iterations)) iterations=ncycle
    if (present(residual)) residual=error
    if (present(status)) then
       status=solve_status
    else if (solve_status/=0) then
       print *, 'Warning in Dyson selfconsistency, status=',solve_status,' residual=',error
    end if
  end subroutine SolveDyson1D


  !****************************************
  !*** Compute Bulk Green's function G0 *** 
  !**************************************** 
  subroutine CompGreensFunc( L1D, spin, z )
    type(T1DLead), intent(in) :: L1D
    integer, intent(in) :: spin
    complex(double), intent(in) :: z

    ! Serial DOS callers retain the module output; quadrature uses explicit scratch.
    call CompGreensFuncWorkspace(L1D, spin, z, GX)
  end subroutine CompGreensFunc

  recursive subroutine CompGreensFuncWorkspace( L1D, spin, z, green )
    use parameters, only: eta
    use numeric, only: CInv
    use constants, only: ui
    implicit none

    type(T1DLead), intent(in) :: L1D
    integer,intent(in) :: spin
    complex(double), intent(in) :: z
    integer :: i, j, ispin, iinv
    complex(double), intent(out) :: green(NXDim,NXDim)
    complex(double) :: Sigma1(NXDim,NXDim), Sigma2(NXDim,NXDim)

    ispin = spin
    if( spin == 2 .and. L1D%NSpin == 1 ) ispin = spin

    call CompSelfEnergy1D( L1D, ispin, z, Sigma1, NXAO, 1 )
    call CompSelfEnergy1D( L1D, ispin, z, Sigma2, NXAO, 2 )
    !
    ! Compute Greensfunction G Tilde of extended device
    !
    do i=1,NXDim
       do j=1,NXDim
          green( i, j ) = (z+ui*eta)*SX(i,j) - HX(i,j,ispin) - Sigma1(i,j) - Sigma2(i,j)
       end do
    end do

    iinv = CInv( green )
  end subroutine CompGreensFuncWorkspace

  !****************************** 
  !*** Spin-resolved Bulk DOS ***
  !******************************
  real(double) function BulkSDOS( L1D, spin, energy )
   !use parameters, only: SOC
    use numeric, only: CInv
    use constants, only: d_pi

    implicit none

    type(T1DLead), intent(inout) :: L1D
    integer, intent(in) :: spin
    complex(double), intent(in) :: energy

    complex(double), dimension( NXDim, NXDim )  :: GXSX
    integer :: i, NAO, NPCAO

    NPCAO = L1D%NPCAO   
    NAO = L1D%NAOrbs

    call CompGreensFunc( L1D, Spin, energy )

    GXSX = matmul( GX, SX )

    BulkSDOS = d_zero
    do i=NAO+1,NAO+NPCAO
       BulkSDOS = BulkSDOS + DIMAG(GXSX(i,i)) 
      !if( SOC ) BulkSDOS = BulkSDOS + DIMAG(GXSX(i+NXAO,i+NXAO)) 
    end do
    BulkSDOS = -BulkSDOS/d_pi
  end function BulkSDOS


  !*******************************************
  !*** Total charge up to energy E         ***
  !*** up to Energy E, lower bound is EMin ***
  !*******************************************
  recursive real(double) function ChargePoint1D(spin, phi, E0, R) result(qj)
    use constants, only: ui, d_pi
    implicit none
    integer, intent(in) :: spin
    real(double), intent(in) :: phi, E0, R
    complex(double) :: z, green(NXDim,NXDim)
    integer :: k, l, NAO, NPCAO

    NAO = Lead1D(WhichLead)%NAOrbs
    NPCAO = Lead1D(WhichLead)%NPCAO
    z = E0 - R*(cos(phi) - ui*sin(phi))
    call CompGreensFuncWorkspace(Lead1D(WhichLead), spin, z, green)
    qj = d_zero
    do k=NAO+1,NAO+NPCAO
       do l=1,NXAO
          qj = qj -(DIMAG(R*(sin(phi)+ui*cos(phi))*green(k,l))/d_pi)*SX(l,k)
       end do
    end do
  end function ChargePoint1D

  recursive real(double) function TotCharge( E )
    use parameters, only: ChargeAcc, Infty
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
!$  use omp_lib, only: omp_in_parallel
    use numeric, only: gauleg 
    use constants, only: ui, d_pi
    implicit none
    
    real(double), intent(in) :: E
    real(double) :: EMin

    integer, parameter :: nmax = 2047
   !integer, parameter :: nmax = 1023
   !integer, parameter :: nmax = 511
   !integer, parameter :: nmax = 255
   !integer, parameter :: nmax = 127
   !integer, parameter :: nmax = 63
    
    real(double) :: q,qq, E0, R
    integer :: n, j, ispin
    real(double), dimension(nmax) :: x, w, contribution
    logical :: parallel_charge



    if (.not.ieee_is_finite(ChargeAcc) .or. ChargeAcc<=0d0) then
       print *, "ERROR: CHARGEACC must be finite and positive for 1D charge integration"
       stop 1
    end if
    EMin = Lead1D(WhichLead)%EMin
    !
    ! Integration contour parameters:
    !
    If (ChargeOffSet == 0.0d0) then
       E0 = 0.5*(E - Infty)
       R  = 0.5*(E + Infty)
    else   
       E0 = 0.5*(E + EMin)
       R  = 0.5*(E - EMin)
    end if 
    !
    ! Computing integral of Green's function
    ! over complex contour using Gauss-Legendre 
    ! quadrature
    !
    n = 1
    qq = d_zero

    do 
       call gauleg(d_zero,d_pi,x(1:n),w(1:n),n)
       q = d_zero
       parallel_charge = .true.
!$     parallel_charge = .not. omp_in_parallel()
       do ispin=1,NXSpin
          ! Each point owns its Green function and self-energy scratch. Preserve
          ! serial summation order to keep adaptive refinement reproducible.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j) IF(parallel_charge) SCHEDULE(STATIC)
          do j=1,n
             contribution(j) = w(j)*ChargePoint1D(ispin, x(j), E0, R)
          end do
!$OMP END PARALLEL DO
          do j=1,n
             q = q + contribution(j)
          end do
       end do
      !if( NXSpin == 1 .and. .not. SOC ) q = q*2.0d0
       if( NXSpin == 1 ) q = q*2.0d0
      !print *, n, q
       if( n > 1 .and. abs(q-qq) < ChargeAcc ) exit
       n = 2*n+1
       if( n > nmax )then
          print '(A,I5,A)', "WARNING: TotCharge/gaussian quadrature has not converged after", nmax , " steps."
          print *, "E = ", E
          print *, "Error = ", abs(q-qq)
          stop 1
       end if
       qq = q
    end do
    if (n < nmax) print '(A,I5,A)', "GauLeg quadrature converged after", n, " steps."
    call flush(6)
    TotCharge = q-ChargeOffset

  end function TotCharge


  !*****************************************************
  !*** Find upper and lower energy bounds for a lead ***
  !*** such that non-zero DOS lies completely inside ***
  !*****************************************************
  ! Generalized isolated-supercell spectrum H0*c = e*S0*c. This provides
  ! a starting estimate, not a guaranteed lower edge of the periodic spectrum.
  ! FindEBounds must still verify the integrated spectral weight below it.
  subroutine SetInitialLeadLowerBound(L1D)
    use numeric, only: CHDiag
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    type(T1DLead), intent(inout) :: L1D
    complex(double) :: overlap(L1D%NAOrbs,L1D%NAOrbs)
    complex(double) :: transform(L1D%NAOrbs,L1D%NAOrbs)
    complex(double) :: hcell(L1D%NAOrbs,L1D%NAOrbs)
    complex(double) :: horth(L1D%NAOrbs,L1D%NAOrbs)
    real(double) :: seval(L1D%NAOrbs),eval(L1D%NAOrbs),lowest
    integer :: n,i,j,spin,info

    n=L1D%NAOrbs
    if (n<1 .or. L1D%NSpin<1) then
       print *, 'ERROR: invalid lead supercell dimensions',L1D%LeadNo
       stop 1
    end if
    do j=1,n
       do i=1,n
          overlap(i,j)=L1D_S0(L1D,i,j)
       end do
    end do
    call CHDiag(overlap,seval,info)
    if (info/=0) then
       print *, 'ERROR: lead overlap diagonalization failed',L1D%LeadNo,info
       stop 1
    end if
    if (.not.all(ieee_is_finite(seval))) then
       print *, 'ERROR: nonfinite lead overlap eigenvalues',L1D%LeadNo
       stop 1
    end if
    if (any(seval<=0d0)) then
       print *, 'ERROR: lead supercell overlap is not positive definite',L1D%LeadNo,minval(seval)
       stop 1
    end if
    do j=1,n
       transform(:,j)=overlap(:,j)/sqrt(seval(j))
    end do
    lowest=huge(0d0)
    do spin=1,L1D%NSpin
       do j=1,n
          do i=1,n
             hcell(i,j)=L1D_H0(L1D,spin,i,spin,j)
          end do
       end do
       horth=matmul(conjg(transpose(transform)),matmul(hcell,transform))
       call CHDiag(horth,eval,info)
       if (info/=0) then
          print *, 'ERROR: lead supercell diagonalization failed',L1D%LeadNo,spin,info
          stop 1
       end if
       if (.not.all(ieee_is_finite(eval))) then
          print *, 'ERROR: nonfinite lead supercell eigenvalues',L1D%LeadNo,spin
          stop 1
       end if
       if (eval(1)<lowest) lowest=eval(1)
    end do
    L1D%EMin=lowest
    print '(A,I2,A,I6,A,F16.8)', 'Initial lower bound, lead ',L1D%LeadNo, &
          ' supercell orbitals ',n,' lowest eigenvalue (eV): ',lowest
    call flush(6)
  end subroutine SetInitialLeadLowerBound

  subroutine FindEBounds( L1D )
    use parameters, only: ChargeAcc, FermiAcc, Infty
    implicit none

    type(T1DLead), intent(inout) :: L1D

    real(double)    :: Q, qmax,lower_step
    
    qmax = 2.0d0*L1D%NPCAO
   !qmax = 2.0d0*(NXAO-(L1D%NPCAO*2))

    ChargeOffSet = d_zero
    WhichLead = L1D%LeadNo
    
    ! 1. Decrease EMin until charge integration becomes zero
    !
    ! Retain geometric downward expansion, including zero/positive seeds.
    lower_step=abs(L1D%EMin)
    if (lower_step<1d0) lower_step=1d0
    if (L1D%EMin < -Infty) then
       print *, 'ERROR: initial lower bound is below -INFTY',L1D%EMin
       stop 1
    end if
    print '(A)', "Searching optimum Emin such that DOS integrates to zero ..."
    do
       Q = TotCharge( L1D%EMin )
       print *, "EMin = ", L1D%EMin, Q
       call FLUSH(6)
       if( abs(Q) < ChargeAcc ) exit
       L1D%EMin = L1D%EMin - lower_step
       lower_step=2d0*lower_step
       if (L1D%EMin < -Infty) then
           print*,'Warning: Emin < -Infty. Decrease this value'
           stop 1
       end if
    end do
    !
    ! 1. Increase EMax from default until charge integration becomes number of spin orbitals
    !
    print '(A)', "Searching optimum EMax such that the DOS integrates to the total spectral weight in the unit cell ..."
    do 
       Q = TotCharge( L1D%EMax )
       print *, "EMax = ", L1D%EMax, Q
       call FLUSH(6)
       if( abs(Q - qmax ) < ChargeAcc ) exit
       L1D%EMax = L1D%EMax + L1D%EMax
    end do
    !
    ! 2. Try to decrease EMax if possible 
    !
!   if( L1D%EMax == 0.0d0 )then
!      EStep = 1.0d0
!      do 
!         L1D%EMax = L1D%EMax - EStep
!         if( L1D%EMax < -Infty )exit
!         Q = TotCharge( L1D%EMax )
!         print *, "EMax = ", L1D%EMax, Q
!         call FLUSH(6)
!         if( abs(Q - qmax) >  ChargeAcc ) exit
!         EStep = 2.0d0*EStep +1.0d0
!      end do
!      L1D%EMax = L1D%EMax + EStep
!   end if
    !
    print *, "EMin=", L1D%EMin, "  EMax=", L1D%EMax
    call FLUSH(6)
  end subroutine FindEBounds


  !***************************************************
  !***                                             ***
  !*** Subroutine to find Fermi level of Bulk lead ***
  !***                                             ***
  !***************************************************
  real(double) function CountedLeadCharge(E) result(q)
    real(double), intent(in) :: E
    FermiEvaluations=FermiEvaluations+1
    q=TotCharge(E)
  end function CountedLeadCharge

  subroutine FindFermi( L1D )
    use parameters, only: Max, FermiAcc
    use numeric, only: BrentRoot, BracketRoot
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use omp_lib, only: omp_get_wtime
    implicit none
    type(T1DLead), intent(inout) :: L1D
    integer :: iterations,evaluations,info,expansions,bracket_evaluations
    real(double) :: EFermi,Delta,residual,width,started,left,right,fleft,fright,halfwidth

    print *, "Searching Fermi energy..."
    WhichLead=L1D%LeadNo
    if (WhichLead<1 .or. WhichLead>size(WarmLeadFermi)) then
       print *, 'ERROR: invalid lead number in Fermi search',WhichLead
       stop 1
    end if
    ChargeOffset=L1D%NPCEl
    FermiEvaluations=0
    print *, 'ChargeOffset',ChargeOffset
    print *, 'Brent Method'
    call flush(6)
    ! Absolute energy accuracy (eV), configured by FERMIACC.
    Delta=FermiAcc
    if (.not.ieee_is_finite(Delta) .or. Delta<=0d0) then
       print *, 'ERROR: FERMIACC must be finite and positive'
       stop 1
    end if
    halfwidth=WarmFermiHalfWidth
    if (2d0*Delta>halfwidth) halfwidth=2d0*Delta
    started=omp_get_wtime()
    if (WarmLeadFermiValid(WhichLead)) then
       print '(A,I2,A,F14.7)', 'Warm Fermi bracket, lead ',WhichLead, &
                              ' previous energy (eV): ',WarmLeadFermi(WhichLead)
       call BracketRoot(CountedLeadCharge,WarmLeadFermi(WhichLead), &
                        L1D%EMin,L1D%EMax,halfwidth, &
                        left,right,fleft,fright,expansions,bracket_evaluations,info)
       print '(A,I5,A,2F16.7)', 'Lead Fermi bracket expansions: ',expansions, &
                              ' interval (eV): ',left,right
       if (info/=0) then
          print *, 'ERROR: lead Fermi bracket failed, status = ',info
          stop 1
       end if
       ! Endpoint values are current-Hamiltonian evaluations, not stored history.
       EFermi=BrentRoot(CountedLeadCharge,left,right,Delta, &
                       5*Max,iterations,evaluations,info,residual,width,fleft,fright)
    else
       print '(A,I2)', 'Full Fermi bracket, no previous solution for lead ',WhichLead
       EFermi=BrentRoot(CountedLeadCharge,L1D%EMin,L1D%EMax,Delta, &
                       5*Max,iterations,evaluations,info,residual,width)
    end if
    print '(A,I5,A,I5)', 'Lead Fermi iterations: ',iterations, &
                        ' charge evaluations: ',FermiEvaluations
    print '(A,ES16.7,A,ES16.7)', 'Lead Fermi residual: ',residual, &
                                ' bracket width (eV): ',width
    print '(A,F12.6)', 'Lead Fermi wall seconds: ',omp_get_wtime()-started
    if (info/=0) then
       print *, 'ERROR: lead Fermi search failed, Brent status = ',info
       stop 1
    end if
    L1D%EFermi=EFermi
    ! Store only the energy, before AdjustFermi shifts this lead's matrices.
    WarmLeadFermi(WhichLead)=EFermi
    WarmLeadFermiValid(WhichLead)=.true.
    print '(A,I3,A,F10.5)', 'EFermi in bulk electrode ',WhichLead,' = ',L1D%EFermi
    call flush(6)
  end subroutine FindFermi


  !*********************************************
  !*** Routine to adjust Fermi-level to zero ***
  !*********************************************
  !
  ! Shifts Lead Hamiltonian to adjust Fermi level to zero
  !
  subroutine AdjustFermi( ilead )
   !use numeric, only : CsetId
    implicit none
    
    integer, intent(in) :: ilead
    integer :: ispin, i, j, ipc, NSpin, NPC, NPCAO

    NSpin = Lead1D(ilead)%NSpin
    NPC = Lead1D(ilead)%NPC
    NPCAO = Lead1D(ilead)%NPCAO
    !
    ! Shift Lead Hamiltonian so that EFermi becomes = 0
    !
    do ispin = 1, NSpin
       do ipc = 0, NPC
          do i = 1, NPCAO
             do j = 1, NPCAO
                Lead1D(ilead)%vn( ispin, ispin, ipc, i, j ) = &
                     Lead1D(ilead)%vn( ispin, ispin, ipc, i, j ) - Lead1D(ilead)%EFermi * Lead1D(ilead)%sn( ipc, i, j )
             end do
          end do
       end do
    end do
    !
    ! Also shifts EMin, EMax to zero
    !
    Lead1D(ilead)%EMin = Lead1D(ilead)%EMin - Lead1D(ilead)%EFermi
    Lead1D(ilead)%EMax = Lead1D(ilead)%EMax - Lead1D(ilead)%EFermi       
    print *, "Lead Hamiltonian shifted by -EFermi = ", -Lead1D(ilead)%EFermi
    !! Lead1D(ilead)%EFermi = 0.0d0
  end subroutine AdjustFermi


end module OneDLead


!!$    !
!!$    ! *** Hamiltonian of super-cell
!!$    !     
!!$    !       / v0    v1    v2    ... vN-1 \
!!$    ! H0 = |  v'1   v0    v1    ... vN-2  |
!!$    !      |  :     :     :         :     |
!!$    !       \ v'N-1 v'N-2 v`N-3 ... v0   /
!!$    !      
!!$    ! *** Coupling between super-cells (from right to left)
!!$    ! 
!!$    !        / vN       0  \
!!$    !       |  ·  ·         |
!!$    ! V2 =  |  ·    ·       |
!!$    !       |  ·      ·     |
!!$    !        \ v1  ···  vN /
!!$    !
!!$    ! where N =  NPC = Number of primitve cells in super-cell, 
!!$    ! vn hopping from n-th right neighbour primitve cell
!!$    ! v'n = hermitian adjoint of v_n



