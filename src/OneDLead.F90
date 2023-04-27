!*********************************************************!
!*********************  ANT.G-2.7.0  *********************!
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
     real(double) :: EMin = -10.0, EMax = 10.0
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
    integer :: ixyzf

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
  subroutine CompSelfEnergy1D( L1D, spin, z, Sigma, NBasis, idir )
    use parameters, only: eta
    use constants, only: ui, c_zero
    implicit none

    type(T1DLead), intent(in) :: L1D
    integer, intent(in) :: spin
    complex(double), intent(in) :: z
   !complex(double), dimension( NBasis*(L1D%SO+1), NBasis*(L1D%SO+1) ), intent(out) :: Sigma
    complex(double), dimension(:,:), intent(out) :: Sigma
    integer , intent(in) :: idir, NBasis !! AO1, 

    integer :: i,j,AOrb1,AOrb2, ispin, jspin, NAO, ierr

   !complex(double), dimension(:,:), allocatable :: Veff1, Heff, Veff2
    complex(double), dimension(L1D%NAOrbs,L1D%NAOrbs) :: Veff1, Heff, Veff2
    complex(double), dimension(L1D%NAOrbs,L1D%NAOrbs) :: sig

    save 

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
       call SolveDyson1D( sig, Veff1, Heff, Veff2, NAO )
       AOrb1 = 1
       AOrb2 = NAO
    else
       call SolveDyson1D( sig, Veff2, Heff, Veff1, NAO )
       AOrb1=NBasis-NAO+1
       AOrb2=NBasis
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
  complex(double) function L1D_H0( L1D, ispin, i, jspin, j )
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
  complex(double) function L1D_V1( L1D, ispin, i, jspin, j )
    type(T1DLead),intent(IN) :: L1D
    integer, intent(IN) ::  ispin, i, j, jspin
    L1D_V1 = conjg( L1D_V2( L1D, jspin, j, ispin, i ) )
  end function L1D_V1

  !*********************************************
  !*** Lead Supercell Hopping (to the right) ***
  !*********************************************
  complex(double) function L1D_V2( L1D, ispin, i, jspin, j )
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
  complex(double) function L1D_S0( L1D, i, j )
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
  complex(double) function L1D_S1( L1D, i, j )
    type(T1DLead),intent(IN) :: L1D
    integer, intent(IN) :: i, j
    L1D_S1 = conjg( L1D_S2( L1D, j, i ) )
  end function L1D_S1

  !*********************************
  ! Lead supercell overlap to right 
  !*********************************
  complex(double) function L1D_S2( L1D, i, j )
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
    use parameters, only: NEmbed, NPC
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
    Lead1D(ilead)%NPCAtoms  = NEmbed(ilead)
    NPCAtoms = Lead1D(ilead)%NPCAtoms
    if (ilead ==1 ) then
      Lead1D(ilead)%NPCEl    = NEmbed(ilead)*GetAtmChg(1)
      Lead1D(ilead)%NPCAO    = NEmbed(ilead)*NAOAtom(1)
    else 
      Lead1D(ilead)%NPCEl    = NEmbed(ilead)*GetAtmChg(GetNAtoms())  
      Lead1D(ilead)%NPCAO    = NEmbed(ilead)*NAOAtom(GetNatoms())
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
    if(proc_id == 0)then
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
    end if
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
    
    real(double) :: EMin, EMax, DE, energy
    complex(double) :: zenergy
    character(len=12),parameter :: file(2) = (/"Lead1DOS.dat","Lead2DOS.dat"/)
    integer :: ispin, ios, iunit
    
    print '(A,I1,A,A)', "Printing DOS of Lead ", ilead, "to file ", file(ilead)
    
    EMin = Lead1D(ilead)%EMin
    EMax = Lead1D(ilead)%EMax

    iunit=fopen(file(ilead),'unknown', ios)
    do energy=EMin, EMax, DE 
      !zenergy=energy-Lead1D(ilead)%EFermi+gamma*ui
       zenergy=energy+gamma*ui
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
  subroutine SolveDyson1D( Sigma, Veff1, Heff, Veff2, NDim )
    use parameters, only: conv => selfacc, alpha => L1DALPHA, MaxCycle => L1DMaxCyc
    use constants, only: c_zero, ui
   !use device, only: OnSigma
#ifdef PGI    
   use lapack_blas
#endif
    implicit none

    integer, intent(in) :: NDim

    complex(double), dimension(NDim,NDim),intent(out) :: Sigma 
    complex(double), dimension(NDim,NDim),intent(in) :: Veff1, Heff, Veff2

    integer, dimension(NDim) :: ipiv
    integer :: ispin, i, j, info, ierr
    real(double) :: error
    integer :: ncycle
    
    !auxiliary self energy for self-consistent calculation, temporary matrix
    complex(double), dimension(NDim,NDim) :: Sigma_aux, temp 
    !work array for inversion
    complex(double), dimension( 4*NDim ) :: work

    !
    ! Initializing Sigma
    !
    Sigma = c_zero
    do i=1,NDim
       Sigma(i,i)=-ui
    end do
    !
    ! Iterative solution of Dyson equation
    !
    error = 1.0d0
    ncycle = 0
    do while ( error.gt.conv .AND. ncycle < MaxCycle )
       ncycle = ncycle+1
       do i=1,NDim
          do j=1,NDim
             Sigma_aux(i,j) = Heff(i,j) - Sigma(i,j)
          end do
       end do
       !
       ! Inverting
       !
       call zgetrf(NDIM,NDIM,Sigma_aux,NDIM,ipiv,info)
       call zgetri(NDIM,Sigma_aux,NDIM,ipiv,work,4*NDIM,info)
       !
       ! Matrix multiplication: Veff1 ( Heff - Sigma )^-1 Veff2
       !
       call zgemm('N','N',NDIM,NDIM,NDIM, c_one,Sigma_aux,NDIM,Veff2,NDIM, c_zero,temp,NDIM )
       call zgemm('N','N',NDIM,NDIM,NDIM, c_one,Veff1,NDIM,temp,NDIM,c_zero, Sigma_aux,NDIM )
       !
       ! Mixing with old self-energy matrix
       !
       do i=1,NDIM
          do j=1,NDIM
             Sigma(i,j) = alpha*Sigma_aux(i,j)+(1.0d0-alpha)*Sigma(i,j)
             error = error+2.0d0*abs(Sigma_aux(i,j)-Sigma(i,j))
          end do
       end do     
       error=error/(2.0*NDIM*NDIM)
    enddo 
   !print*,'ncycle Dyson',ncycle
    if (ncycle == MaxCycle) print*, 'Warning in Dyson selfconsistency'

  end subroutine SolveDyson1D


  !****************************************
  !*** Compute Bulk Green's function G0 *** 
  !**************************************** 
  subroutine CompGreensFunc( L1D, spin, z ) 
    use parameters, only: eta
    use numeric, only: CInv
    use constants, only: ui
    implicit none

    type(T1DLead), intent(in) :: L1D
    integer,intent(in) :: spin
    complex(double), intent(in) :: z
    integer :: i, j, ispin, jspin, NPCAO, NAO, leadno, iinv

    ispin = spin
    if( spin == 2 .and. L1D%NSpin == 1 ) ispin = spin

    call CompSelfEnergy1D( L1D, ispin, z, Sigma1, NXAO, 1 )
    call CompSelfEnergy1D( L1D, ispin, z, Sigma2, NXAO, 2 )
    !
    ! Compute Greensfunction G Tilde of extended device
    !
    do i=1,NXDim
       do j=1,NXDim
          GX( i, j ) = (z+ui*eta)*SX(i,j) - HX(i,j,ispin) - Sigma1(i,j) - Sigma2(i,j)
       end do
    end do

    iinv = CInv( GX )       
   !print*,'z=',z,'SigmaL=',Sigma1(1,1),'SigmaR=',Sigma2(NXDim-10,NXDim-10),'Green=', GX(1,1)
  end subroutine CompGreensFunc

  !****************************** 
  !*** Spin-resolved Bulk DOS ***
  !******************************
  real(double) function BulkSDOS( L1D, spin, energy )
    use parameters, only: SOC
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
       if( SOC ) BulkSDOS = BulkSDOS + DIMAG(GXSX(i+NXAO,i+NXAO)) 
    end do
    BulkSDOS = -BulkSDOS/d_pi
  end function BulkSDOS


  !*******************************************
  !*** Total charge up to energy E         ***
  !*** up to Energy E, lower bound is EMin ***
  !*******************************************
  real(double) function TotCharge( E )
    use parameters, only: ChargeAcc, Infty, SOC
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
    
    real(double) :: q,qq, E0, R, phi, qj
    integer :: n, j, ispin, k, l, NAO, NPCAO
    real(double), dimension(nmax) :: x, w

    complex(double) :: z


    NAO = Lead1D(WhichLead)%NAOrbs
    NPCAO = Lead1D(WhichLead)%NPCAO
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

    do 
       call gauleg(d_zero,d_pi,x(1:n),w(1:n),n)
       q = d_zero
       do ispin=1, NXSpin
          do j = 1,n
             phi = x(j)
             z = E0 - R*(cos(phi) - ui*sin(phi))
             
             call CompGreensFunc( Lead1D(WhichLead), ispin, z ) 

             qj= 0.0d0
             do k=NAO+1, NAO+NPCAO
                do l=1,NXAO
                   !
                   ! Charge = Tr[ P S ]
                   !
                  !q = q -(w(j)*DIMAG(R*(sin(phi)+ui*cos(phi))*GX(k,l))/d_pi)*SX(l,k)
                   qj = qj -(DIMAG(R*(sin(phi)+ui*cos(phi))*GX(k,l))/d_pi)*SX(l,k)
                  !if( SOC ) q = q -(w(j)*DIMAG(R*(sin(phi)+ui*cos(phi))*GX(k+NXAO,l+NXAO))/d_pi)*SX(l+NXAO,k+NXAO)
                   if( SOC ) qj = qj -(DIMAG(R*(sin(phi)+ui*cos(phi))*GX(k+NXAO,l+NXAO))/d_pi)*SX(l+NXAO,k+NXAO)
                end do
             end do
            !if (j == 1) print*,'z= ',z,'phi=',x(1),GX(1,1),'DOS=',qj*2.0
             q = q + w(j)*qj
          end do
            !print*,'----------------'
       end do
       if( NXSpin == 1 .and. .not. SOC ) q = q*2.0d0
      !print *, n, q
      !if( n > 1 .and. (q == d_zero .or. abs(q-qq) < ChargeAcc*100 ) ) exit  
       if( n > 1 .and. (q == d_zero .or. abs(q-qq) < 0.1 ) ) exit  
       n = 2*n+1
       if( n > nmax )then
          print '(A,I5,A)', "WARNING: TotCharge/gaussian quadrature has not converged after", nmax , " steps."
         !print *, "E = ", E
         !print *, "deltaq = ", abs(q-qq)
          exit
       end if
       qq = q
    end do
    print '(A,I5,A)', "GauLeg quadrature converged after", n, " steps."
    call flush(6)
    TotCharge = q-ChargeOffset

  end function TotCharge


  !*****************************************************
  !*** Find upper and lower energy bounds for a lead ***
  !*** such that non-zero DOS lies completely inside ***
  !*****************************************************
  subroutine FindEBounds( L1D )
    use parameters, only: ChargeAcc, FermiAcc, Infty
    implicit none

    type(T1DLead), intent(inout) :: L1D

    real(double)    :: EStep, Q, qmax
    
    qmax = 2.0d0*L1D%NPCAO
   !qmax = 2.0d0*(NXAO-(L1D%NPCAO*2))

    ChargeOffSet = d_zero
    WhichLead = L1D%LeadNo
    
    ! 1. Decrease EMin until charge integration becomes zero
    !
    EStep = 10.0d0
    print '(A)', "Searching optimum Emin such that DOS integrates to zero ..."
    do
       Q = TotCharge( L1D%EMin )
       print *, "EMin = ", L1D%EMin, Q
       call FLUSH(6)
       if( abs(Q) <  0.1 ) exit
       L1D%EMin = L1D%EMin - EStep
       if (L1D%EMin < -Infty) then
           print*,'Warning: Emin < -Infty. Decrease this value'
           stop
       end if
    end do
    L1D%EMin = L1D%EMin - 10.0    ! for safety
    L1D%EMax = -L1D%EMin

    !
    ! 1. Increase EMax until charge integration becomes number of spin orbitals
    !
    print '(A)', "Searching optimum EMax such that the DOS integrates to the total spectral weight in the unit cell ..."
    do 
       Q = TotCharge( L1D%EMax )
       print *, "EMax = ", L1D%EMax, Q
       call FLUSH(6)
       if( abs(Q - qmax ) <  0.1 ) exit
       L1D%EMax = L1D%EMax + EStep
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
  subroutine FindFermi( L1D )
    use parameters, only: ChargeAcc, FermiAcc, Infty, Max
    use numeric, only: bisec, Muller
    implicit none

    type(T1DLead), intent(inout) :: L1D

    integer :: cond, nsteps ,k
    real(double) :: EStep, Q, EFermi
    real(double) :: EF0,EF1,EF2,Z, Delta, Epsilon

    print *, "Searching Fermi energy..."
    
    WhichLead = L1D%LeadNo
   !EFermi = L1D%EFermi
    ChargeOffset = L1D%NPCEl
    print*, 'ChargeOffset', ChargeOffset
    
   !if( L1D%EFermi == 0.0d0 )then
       print *, "Bisec Method"
       call FLUSH(6)
       Delta=0.1
       EFermi = bisec(TotCharge,L1D%EMin,L1D%EMax,Delta,5*Max,K) 
       L1D%EFermi = EFermi
   !end if

   !print '(A,F10.5)' , "Initial EFermi = ", L1D%EFermi

   !EF0 = EFermi-1.0
   !EF1 = EFermi
   !EF2 = EFermi+1.0

   !nsteps = 25
   !Delta=FermiAcc
   !Epsilon=ChargeAcc
    
   !print *, "Muller Method"
   !call FLUSH(6)
    
   !call MULLER(TotCharge,EF0,EF1,EF2,Delta,Epsilon,nsteps,EFermi,Z,K,Cond)
   !L1D%EFermi = EFermi

    print '(A,I3,A,F10.5)', 'EFermi in bulk electrode ', WhichLead,' = ', L1D%EFermi 
    call FLUSH(6)
   !init = .true.

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



