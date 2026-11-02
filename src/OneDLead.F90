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
    ! Periodic band filling is the sole lead Fermi-level search.
    call FindBandFermi(Lead1D(ilead))
    call CreateXSC(Lead1D(ilead))
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

  recursive subroutine CompGreensFuncWorkspace( L1D, spin, z, green, inversion_status )
    use parameters, only: eta
    use numeric, only: CInv
    use constants, only: ui
    implicit none

    type(T1DLead), intent(in) :: L1D
    integer,intent(in) :: spin
    integer,intent(out),optional :: inversion_status
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
    if(present(inversion_status)) inversion_status=iinv
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
  recursive real(double) function ChargePoint1D(spin, phi, E0, R, point_status) result(qj)
    use constants, only: ui, d_pi
    implicit none
    integer, intent(in) :: spin
    integer, intent(out),optional :: point_status
    real(double), intent(in) :: phi, E0, R
    complex(double) :: z, green(NXDim,NXDim)
    integer :: k, l, NAO, NPCAO,info

    NAO = Lead1D(WhichLead)%NAOrbs
    NPCAO = Lead1D(WhichLead)%NPCAO
    z = E0 - R*(cos(phi) - ui*sin(phi))
    call CompGreensFuncWorkspace(Lead1D(WhichLead), spin, z, green,info)
    if(present(point_status)) point_status=info
    if(info/=0) then
       if(.not.present(point_status)) then
          print *, 'ERROR: lead Green-function inversion failed',info
          stop 1
       end if
       qj=0d0;return
    end if
    qj = d_zero
    do k=NAO+1,NAO+NPCAO
       do l=1,NXAO
          qj = qj -(DIMAG(R*(sin(phi)+ui*cos(phi))*green(k,l))/d_pi)*SX(l,k)
       end do
    end do
  end function ChargePoint1D

  recursive subroutine AdaptiveLeadCharge(E0,R,accuracy,q,error,status,evaluations)
    use numeric, only: AdaptiveGK15
    use constants, only: d_pi
    use omp_lib
    implicit none
    real*8,intent(in) :: E0,R,accuracy
    real*8,intent(out) :: q,error
    integer,intent(out) :: status,evaluations
    real*8 :: result(1),tolerances(1),ratio
    integer :: panels
    tolerances=accuracy
    call AdaptiveGK15(AdaptiveLeadBatch,1,0d0,d_pi,tolerances,result,ratio,evaluations,panels,status,(/E0,R/))
    q=result(1);error=ratio*accuracy
  end subroutine AdaptiveLeadCharge

  recursive subroutine AdaptiveLeadBatch(ncomp,npoint,x,y,info,context)
      use omp_lib
      integer,intent(in) :: ncomp,npoint
      real*8,intent(in) :: x(npoint),context(:)
      real*8,intent(out) :: y(ncomp,npoint)
      integer,intent(out) :: info
      integer :: j,spin,node_info,point_status(npoint)
      real*8 :: E0,R,value
      E0=context(1);R=context(2)
      info=0
!$omp parallel do default(shared) private(j,spin,node_info,value) if(.not.omp_in_parallel()) schedule(static)
      do j=1,npoint
         y(1,j)=0d0;point_status(j)=0
         do spin=1,NXSpin
            value=ChargePoint1D(spin,x(j),E0,R,node_info)
            if(node_info/=0) point_status(j)=node_info
            y(1,j)=y(1,j)+value
         end do
         if(NXSpin==1) y(1,j)=2d0*y(1,j)
      end do
!$omp end parallel do
      if(any(point_status/=0)) info=1
  end subroutine AdaptiveLeadBatch


  recursive real(double) function TotCharge( E, accuracy, error_estimate, status )
    use parameters, only: ChargeAcc, Infty, AdaptiveQuad
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
!$  use omp_lib, only: omp_in_parallel
    use numeric, only: gauleg 
    use constants, only: ui, d_pi
    implicit none
    
    real(double), intent(in) :: E
    real(double), intent(in), optional :: accuracy
    real(double), intent(out), optional :: error_estimate
    integer, intent(out), optional :: status
    real(double) :: tolerance,quadrature_error
    real(double) :: EMin

    integer, parameter :: nmax = 2047
   !integer, parameter :: nmax = 1023
   !integer, parameter :: nmax = 511
   !integer, parameter :: nmax = 255
   !integer, parameter :: nmax = 127
   !integer, parameter :: nmax = 63
    
    real(double) :: q,qq, E0, R
    integer :: n, j, ispin, adaptive_status, evaluations
    real(double), dimension(nmax) :: x, w, contribution
    logical :: parallel_charge



    tolerance=ChargeAcc
    if(present(accuracy)) tolerance=accuracy
    if(present(status)) status=0
    if(present(error_estimate)) error_estimate=huge(0d0)
    TotCharge=0d0
    if (.not.ieee_is_finite(tolerance) .or. tolerance<=0d0) then
       print *, "ERROR: CHARGEACC must be finite and positive for 1D charge integration"
       if(present(status)) then
          status=1;return
       end if
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
    if(AdaptiveQuad) then
       call AdaptiveLeadCharge(E0,R,tolerance,q,quadrature_error,adaptive_status,evaluations)
       if(present(error_estimate)) error_estimate=quadrature_error
       if(adaptive_status/=0) then
          if(present(status)) then
             status=2;return
          end if
          print *, 'ERROR: adaptive lead quadrature failed, status/evaluations=',adaptive_status,evaluations
          stop 1
       end if
       TotCharge=q-ChargeOffset
       print '(A,I6,A,ES12.4)', 'Adaptive lead quadrature evaluations: ',evaluations, &
            ' estimated charge error: ',quadrature_error
       return
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
       quadrature_error=abs(q-qq)
       if (.not.ieee_is_finite(q).or..not.ieee_is_finite(quadrature_error)) then
          if(present(status)) then
             status=3;return
          end if
          print *, 'ERROR: nonfinite lead charge integral'
          stop 1
       end if
       if(present(error_estimate)) error_estimate=quadrature_error
       if( n > 1 .and. quadrature_error < tolerance ) exit
       n = 2*n+1
       if( n > nmax )then
          print '(A,I5,A)', "WARNING: TotCharge/gaussian quadrature has not converged after", nmax , " steps."
          print *, "E = ", E
          print *, "Error = ", abs(q-qq)
          if(present(status)) then
             status=2;return
          end if
          stop 1
       end if
       qq = q
    end do
    if (n < nmax) print '(A,I5,A)', "GauLeg quadrature converged after", n, " steps."
    call flush(6)
    TotCharge = q-ChargeOffset

  end function TotCharge


  ! Bloch matrices use every stored neighbour block of the chosen unit cell.
  subroutine LeadBandPoint(L1D,k,energies,info,min_overlap)
    use numeric, only: CHDiag
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    type(T1DLead),intent(in) :: L1D
    real(double),intent(in) :: k
    real(double),intent(out) :: energies(L1D%NPCAO,L1D%NSpin)
    integer,intent(out) :: info
    real(double),optional,intent(out) :: min_overlap
    complex(double) :: sk(L1D%NPCAO,L1D%NPCAO),hk(L1D%NPCAO,L1D%NPCAO)
    complex(double) :: t(L1D%NPCAO,L1D%NPCAO),phase
    real(double) :: se(L1D%NPCAO)
    integer :: r,j,spin
    sk=L1D%sn(0,:,:)
    do r=1,L1D%NPC
       phase=cmplx(cos(r*k),sin(r*k),kind=double)
       sk=sk+phase*L1D%sn(r,:,:)+conjg(phase)*conjg(transpose(L1D%sn(r,:,:)))
    end do
    call CHDiag(sk,se,info)
    if(info/=0) return
    if(.not.all(ieee_is_finite(se))) then
       info=-1;return
    end if
    if(present(min_overlap)) min_overlap=minval(se)
    if(any(se<=0d0)) then
       info=-2;return
    end if
    do j=1,L1D%NPCAO
       t(:,j)=sk(:,j)/sqrt(se(j))
    end do
    do spin=1,L1D%NSpin
       hk=L1D%vn(spin,spin,0,:,:)
       do r=1,L1D%NPC
          phase=cmplx(cos(r*k),sin(r*k),kind=double)
          hk=hk+phase*L1D%vn(spin,spin,r,:,:)+ &
               conjg(phase)*conjg(transpose(L1D%vn(spin,spin,r,:,:)))
       end do
       hk=matmul(conjg(transpose(t)),matmul(hk,t))
       call CHDiag(hk,energies(:,spin),info)
       if(info/=0) return
    end do
    if(.not.all(ieee_is_finite(energies))) info=-3
  end subroutine LeadBandPoint

  real(double) function BandPopulation(energies,mu,width,degeneracy) result(q)
    implicit none
    real(double),intent(in) :: energies(:,:,:),mu,width,degeneracy
    real(double) :: pi
    pi=acos(-1d0)
    ! Exact integral of each ETA-broadened pole from minus infinity to mu.
    q=degeneracy*sum(0.5d0+atan((mu-energies)/width)/pi)/size(energies,2)
  end function BandPopulation

  subroutine BandSearchFailure(message)
    character(*),intent(in) :: message
    print *, 'ERROR: ',message
    call flush(6)
    stop 1
  end subroutine BandSearchFailure

  subroutine FindBandFermi(L1D)
    use parameters, only: eta,ChargeAcc,FermiAcc
    use omp_lib, only: omp_get_wtime
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    type(T1DLead),intent(inout) :: L1D
    real(double),allocatable :: bands(:,:,:),overlap_min(:)
    integer,allocatable :: status(:)
    integer :: nk,j,it,level,worst
    real(double) :: pi,k,lo,hi,mu,q,oldmu,oldq,deg,step,started
    logical :: converged
    if(eta<=0d0.or.ChargeAcc<=0d0.or.FermiAcc<=0d0.or. &
       .not.ieee_is_finite(eta+ChargeAcc+FermiAcc)) call BandSearchFailure('Invalid band-search tolerances/ETA')
    if(L1D%NSpin<1.or.L1D%NSpin>2.or.L1D%NPCEl<=0.or. &
       L1D%NPCEl>=2*L1D%NPCAO) call BandSearchFailure('Invalid band-search filling')
    pi=acos(-1d0);deg=1d0
    if(L1D%NSpin==1) deg=2d0
    nk=256;oldmu=huge(0d0);started=omp_get_wtime();converged=.false.
    do level=1,9
       allocate(bands(L1D%NPCAO,nk,L1D%NSpin),status(nk),overlap_min(nk))
       overlap_min=huge(0d0)
!$omp parallel do default(shared) private(j,k) schedule(static)
       do j=1,nk
          k=2d0*pi*(dble(j)-0.5d0)/nk
          call LeadBandPoint(L1D,k,bands(:,j,:),status(j),overlap_min(j))
       end do
!$omp end parallel do
       if(any(status/=0)) then
          print *, 'ERROR: Bloch diagonalization/overlap, lead ',L1D%LeadNo,minval(status),maxval(status)
          worst=minloc(overlap_min,dim=1)
          print '(A,F16.9,A,ES16.8)', 'Worst overlap k (radians): ', &
               2d0*pi*(dble(worst)-0.5d0)/nk,' minimum eigenvalue: ',overlap_min(worst)
          call flush(6)
          stop 1
       end if
       lo=minval(bands);hi=maxval(bands);step=hi-lo+eta
       do it=1,64
          if(BandPopulation(bands,lo,eta,deg)<L1D%NPCEl) exit
          lo=lo-step;step=2d0*step
       end do
       if(it>64) call BandSearchFailure('Band lower bracket failed')
       step=hi-lo+eta
       do it=1,64
          if(BandPopulation(bands,hi,eta,deg)>L1D%NPCEl) exit
          hi=hi+step;step=2d0*step
       end do
       if(it>64) call BandSearchFailure('Band upper bracket failed')
       do it=1,256
          mu=lo+(hi-lo)*0.5d0
          q=BandPopulation(bands,mu,eta,deg)
          if(q<L1D%NPCEl) then
             lo=mu
          else
             hi=mu
          end if
          if(hi-lo<=FermiAcc.and.abs(q-L1D%NPCEl)<=ChargeAcc) exit
       end do
       if(it>256) call BandSearchFailure('Band filling root failed')
       oldq=huge(0d0)
       if(level>1) oldq=abs(BandPopulation(bands,oldmu,eta,deg)-L1D%NPCEl)
       print '(A,I2,A,I8,A,F16.8,A,ES12.4,A,ES12.4)', &
          'Band Fermi lead ',L1D%LeadNo,' k points ',nk,' EF (eV) ',mu, &
          ' charge residual ',q-L1D%NPCEl,' previous-EF charge residual ',oldq
       call flush(6)
       if(level>1) converged=abs(mu-oldmu)<=FermiAcc.and.oldq<=ChargeAcc
       if(converged) then
          L1D%EFermi=mu
          L1D%EMin=minval(bands)-25d0;L1D%EMax=maxval(bands)+25d0
       end if
       deallocate(bands,status,overlap_min)
       if(converged) exit
       oldmu=mu;nk=2*nk
    end do
    if(.not.converged) call BandSearchFailure('Band Fermi k mesh did not converge')
    print '(A,I3,A,F12.6)', 'EFermi in bulk electrode ',L1D%LeadNo,' = ',L1D%EFermi
    print '(A,F12.6)', 'Band Fermi wall seconds: ',omp_get_wtime()-started
    call flush(6)
  end subroutine FindBandFermi




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



