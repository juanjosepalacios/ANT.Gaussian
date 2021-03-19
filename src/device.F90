!*********************************************************!
!*********************  ANT.G-2.5.2  *********************!
!*********************************************************!
!                                                         !
!  Copyright (c) by                                       !
!                                                         !
!  Juan Jose Palacios (1)                                 !
!  David Jacob (2)                                        !
!  Maria Soriano (1)                                      !
!  Angel J. Perez-Jimenez (3)                             !
!  Emilio SanFabian (3)                                   !
!  Jose Antonio Antonio Verges (4)                        !
!  Enrique Louis (5)                                      !
!  Wynand Dednam (5)                                      !
!                                                         !
! (1) Departamento de Fisica de la Materia Condensada     !
!     Universidad Autonoma de Madrid                      !    
!     28049 Madrid (SPAIN)                                !
! (2) Theory Department                                   !
!     Max-Planck-Institute for Microstructure Physics     !
!     Halle, 06120 (GERMANY)                              !
! (3) Departamento de Quimica Fisica                      !
!     Universidad de Alicante                             !
!     03690 Alicante (SPAIN)                              !
! (4) Insto. de Ciencias de Materiales de Madrid (ICMM)   !
!     Consejo Superior de Investigacion y Ciencia (CSIC)  !
!     28049 Madrid (SPAIN)                                !
! (5) Departamento de Fisica Aplicada                     !
!     Universidad de Alicante                             !    
!     03690 Alicante (SPAIN)                              !
!                                                         !
!*********************************************************!
  MODULE device
!*********************************************************!
!  Main module for computation of device Green's function !
!  and transport                                          !
!*********************************************************!
  use ANTCommon
  implicit none
  save

  private

  public :: DevNAOrbs, DevNSpin, DevShift, DevFockMat, DevDensMat, SetDevDensMat, SetDevFockMat
  public :: DevHWFockMat, SetDevHWFockMat, DevDGibbsYMat, SetDevDGibbsYMat
  public :: DevDGibbsYKernel1Mat, SetDevDGibbsYKernel1Mat, DevDGibbsYKernel2Mat, SetDevDGibbsYKernel2Mat
  public :: BuildLiouvillian, CompCGibbsY  
  public :: LeadsOn, SwitchOnLeads, SecantOn, SwitchOnSecant, SwitchOffSecant
  public :: EvaluationOn, SwitchOnEvaluation
  public :: SwitchOnChargeCntr, SwitchOffChargeCntr, SwitchOffSpinLock, SwitchOnSpinLock
  public :: InitDevice, ReadDensMat, ReadFockMat, CleanUpDevice, InitElectrodes, Transport
  public :: IntDDOStimesE, DDOStimesE0, WorkFock, WorkEnergy  


  !*****************************
  ! Module's internal variables
  !*****************************

  ! *** Number of atomic orbitals, number of non-degenerate spin-channels ***
  integer :: NAOrbs, NSpin, DNAOrbs

  ! *** Total number of electrons in central device region
  integer :: NCDEl, NCDAO1, NCDAO2

  ! *** Actual electron charge for certain Fermi energy ***
  real*8 :: QAlpha, QBeta, Q_SOC

  ! *** Overlap matrix S of device ***
  real*8, dimension(:,:),allocatable :: SD, InvSD
  real*8, dimension(:,:), allocatable :: S_SOC, InvS_SOC

  ! *** Complex S^+1/2 matrix ***
  complex*16, dimension(:,:),allocatable :: SPH, SNH
  
  ! *** Hamiltonian and density matrix of device
  real*8, dimension(:,:,:),allocatable :: HD
  real*8, dimension(:,:,:),allocatable :: PD
  real*8, dimension(:,:,:),allocatable :: PDGIBBS
  real*8, dimension(:,:,:),allocatable :: HW
  real*8, dimension(:,:,:),allocatable :: DGibbsY, DGibbsYKernel1, DGibbsYKernel2  
  complex*16, dimension(:,:,:),allocatable :: PDOUT
  complex*16, dimension(:,:,:),allocatable :: PDOUTGIBBS
  complex*16, dimension(:,:,:),allocatable :: HWOUT
  complex*16, dimension(:,:,:),allocatable :: LiouvSOp
  complex*16, dimension(:,:,:),allocatable :: CGibbsY, CGibbsYKernel1, CGibbsYKernel2  
  complex*16, dimension(:,:),allocatable :: H_SOC, PD_SOC,PD_SOC_R,PD_SOC_A
  complex*16, dimension(:,:),allocatable :: PDOUT_SOC

  ! *** Orthognalization matrix for device ***
  real*8, dimension(:,:),allocatable :: OD

  ! *** Energy shifts ***
  real*8 :: shift, ShiftUp, ShiftDown

  ! *** spin excess charge ***
  real*8 :: exspin

  ! *** internal spin variable 1=up,2=down ***
  integer  :: ispin
  
  ! *** Lowest and highest eigen value of Fock matrix ***
  real*8 :: LEV,HEV
 
  ! *** lower and upper band edge of electrode ***
  real*8 :: EMinEc, EMaxEc

  ! *** lower and upper energy bound ***
  real*8 :: EMin, EMax

  ! *** Density of states projected on atoms at the Fermi energy
  real*8, dimension(:,:), allocatable :: AtomDOSEF

  ! *** Control Switches ***
  logical :: Leads       = .false.
  logical :: Secant      = .false.
  logical :: Evaluation  = .false.
  logical :: UDTrans     = .false.
  logical :: ChargeCntr  = .false.
  logical :: SpinLock    = .true.

  ! Whether device Hamiltonain has been orthogonalized
  logical :: HDOrtho = .false.
  
  ! *** In which atom to calculate the spinmu ***
  integer :: spinatom
  ! *** Number of electrons projected on atoms at the Fermi energy
  ! *** Actual electron charge for certain Fermi energy ***
  real*8 :: LocalQAlpha, LocalQBeta
  real*8, dimension(:), allocatable :: alphaelec, betaelec, alphalocalshift, betalocalshift
  real*8, dimension(:,:,:),allocatable :: LocPD
  complex*16, dimension(:,:,:),allocatable :: LocPDOUT  
  
!!$OMP THREADPRIVATE(shift)
  contains

  !********************************
  !*** Public module procedures ***
  !********************************

  !**************************************
  ! Access functions to private entities
  !**************************************


  ! *** Total number of atomic orbitals in device Hilbert space ***
  integer function DevNAOrbs()
    implicit none
    DevNAOrbs=NAOrbs
  end function DevNAOrbs

  ! *** Number of non-degenerate spin bands ***
  integer function DevNSpin() 
    implicit none
    DevNSpin=NSpin
  end function DevNSpin

  ! *** Spin excess charge to equilibrate Fermi level up and down ***
  real*8 function DevShift()
    implicit none
    DevShift=-shift
  end function DevShift

  ! *** Get matrix Element of Density matrix ***
  real*8 function DevFockMat( is, i, j )
    implicit none
    integer, intent(in) :: is, i, j
    DevFockMat = HD(is, i, j)
  end function DevFockMat
  
  ! *** Set matrix Element of Fock matrix ***
  subroutine SetDevFockMat( is, i, j, fij )
    implicit none
    integer, intent(in) :: is, i, j
    real, intent(in) :: fij
    HD(is, i, j) = fij
  end subroutine SetDevFockMat  

  ! *** Get matrix Element of Density matrix ***
  real*8 function DevDensMat( is, i, j )
    implicit none
    integer, intent(in) :: is, i, j
    DevDensMat = PD(is, i, j)
  end function DevDensMat

  ! *** Set matrix Element of Density matrix ***
  subroutine SetDevDensMat( is, i, j, pij )
    implicit none
    integer, intent(in) :: is, i, j
    real*8, intent(in) :: pij
    PD(is, i, j) = pij
  end subroutine SetDevDensMat
  
  ! *** Get matrix Element of HWFock matrix ***
  real function DevHWFockMat( is, i, j )
    implicit none
    integer, intent(in) :: is, i, j
    DevHWFockMat = HW(is, i, j)
  end function DevHWFockMat

  ! *** Set matrix Element of HWFock matrix ***
  subroutine SetDevHWFockMat( is, i, j, hwij )
    implicit none
    integer, intent(in) :: is, i, j
    real, intent(in) :: hwij
    HW(is, i, j) = hwij
  end subroutine SetDevHWFockMat

  ! *** Get matrix Element of CGibbsY matrix ***
  real*8 function DevDGibbsYMat( is, i, j )
    implicit none
    integer, intent(in) :: is, i, j
    DevDGibbsYMat = DGibbsY(is, i, j)
  end function DevDGibbsYMat

  ! *** Set matrix Element of CGibbsY matrix ***
  subroutine SetDevDGibbsYMat( is, i, j, DGibbsYij )
    implicit none
    integer, intent(in) :: is, i, j
    real, intent(in) :: DGibbsYij
    DGibbsY(is, i, j) = DGibbsYij
  end subroutine SetDevDGibbsYMat

  ! *** Get matrix Element of CGibbsYKernel1 matrix ***
  real*8 function DevDGibbsYKernel1Mat( is, i, j )
    implicit none
    integer, intent(in) :: is, i, j
    DevDGibbsYKernel1Mat = DGibbsYKernel1(is, i, j)
  end function DevDGibbsYKernel1Mat

  ! *** Set matrix Element of CGibbsYKernel1 matrix ***
  subroutine SetDevDGibbsYKernel1Mat( is, i, j, DGibbsYKernel1ij )
    implicit none
    integer, intent(in) :: is, i, j
    real, intent(in) :: DGibbsYKernel1ij
    DGibbsYKernel1(is, i, j) = DGibbsYKernel1ij
  end subroutine SetDevDGibbsYKernel1Mat

  ! *** Get matrix Element of CGibbsYKernel2 matrix ***
  real*8 function DevDGibbsYKernel2Mat( is, i, j )
    implicit none
    integer, intent(in) :: is, i, j
    DevDGibbsYKernel2Mat = DGibbsYKernel2(is, i, j)
  end function DevDGibbsYKernel2Mat

  ! *** Set matrix Element of CGibbsYKernel2 matrix ***
  subroutine SetDevDGibbsYKernel2Mat( is, i, j, DGibbsYKernel2ij )
    implicit none
    integer, intent(in) :: is, i, j
    real, intent(in) :: DGibbsYKernel2ij
    DGibbsYKernel2(is, i, j) = DGibbsYKernel2ij
  end subroutine SetDevDGibbsYKernel2Mat  

  ! ***
  logical function LeadsOn()
    implicit none
    LeadsOn = Leads
  end function LeadsOn

  ! *** 
  subroutine SwitchOnLeads
    implicit none
    Leads = .true.
  end subroutine SwitchOnLeads

  ! ***
  logical function SecantOn()
    implicit none
    SecantOn = Secant
  end function SecantOn

  ! ***
  subroutine SwitchOnSecant()
    implicit none
    Secant = .true.
  end subroutine SwitchOnSecant

  ! ***
  subroutine SwitchOffSecant()
    implicit none
    Secant = .false.
  end subroutine SwitchOffSecant

  ! ***
  logical function EvaluationOn()
    implicit none
    EvaluationOn = Evaluation
  end function EvaluationOn

  ! ***
  subroutine SwitchOnEvaluation()
    implicit none
    Evaluation = .true.
  end subroutine SwitchOnEvaluation

  subroutine SwitchOnChargeCntr()
    implicit none
    ChargeCntr = .true.
  end subroutine SwitchOnChargeCntr

  ! ***
  subroutine SwitchOffChargeCntr()
    implicit none
    ChargeCntr = .false.
  end subroutine SwitchOffChargeCntr

  ! ***
  subroutine SwitchOffSpinLock()
    implicit none
    SpinLock = .false.
  end subroutine SwitchOffSpinLock

  ! ***
  subroutine SwitchOnSpinLock()
    implicit none
    SpinLock = .true.
  end subroutine SwitchOnSpinLock
  

  !***********************************************
  !* Initialize device for transport calculation *
  !***********************************************
  subroutine InitDevice( NBasis, UHF, S )
    use constants, only: d_zero
    use numeric, only: RMatPow
    use parameters, only: ElType, FermiStart, Overlap, HybFunc, CompFock, SOC, biasvoltage
    use cluster, only: AnalyseCluster, AnalyseClusterElectrodeOne, AnalyseClusterElectrodeTwo, NAOAtom, NEmbedBL
#ifdef G03ROOT
    use g03Common, only: GetNAtoms, GetAtmChg
#endif
#ifdef G09ROOT
    use g09Common, only: GetNAtoms, GetAtmChg
#endif
    use correlation
    use orthogonalization
    use ANTCommon
    
    implicit none

    integer, intent(in) :: NBasis
    logical, intent(in) :: UHF
    real*8, dimension(NBasis,NBasis),intent(in) :: S

    integer :: AllocErr, ios, iatom, NEmbed1, NEmbed2

    real*8, dimension(NBasis,NBasis) :: RSPH 

    write(ifu_log,*) "-------------------"
    write(ifu_log,*) "Initializing device"
    write(ifu_log,*) "-------------------"

    NAOrbs = NBasis
    DNAOrbs = 2*NBasis
    NSpin=1
    if( UHF ) NSpin=2
    exspin=d_zero

    !if (biasvoltage /= 0.0) allocate(PDOUT(NSpin,NAOrbs,NAOrbs), STAT=AllocErr) 
    ! Dynamic arrays 
    allocate( SD(NAOrbs,NAOrbs), InvSD(NAOrbs,NAOrbs),  HD(NSpin,NAOrbs,NAOrbs),  PD(NSpin,NAOrbs,NAOrbs), PDOUT(NSpin,NAOrbs,NAOrbs), STAT=AllocErr )
    if( AllocErr /= 0 ) then
       print *, "DEVICE/Allocation error for SD, InvSD, SMH, SPH, H, P"
       stop
    end if
    if(CompFock)then
      allocate( HW(NSpin,NAOrbs,NAOrbs), &
         HWOUT(NSpin,NAOrbs,NAOrbs), &
         STAT=AllocErr )
      if( AllocErr /= 0 ) then
        print *, "DEVICE/Allocation error for HW, HWOUT"
        stop
      end if
      allocate( CGibbsY(NSpin,NAOrbs,NAOrbs), &
                DGibbsY(NSpin,NAOrbs,NAOrbs), &
                CGibbsYKernel1(NSpin,NAOrbs,NAOrbs), &
                DGibbsYKernel1(NSpin,NAOrbs,NAOrbs), &
                CGibbsYKernel2(NSpin,NAOrbs,NAOrbs), &
                DGibbsYKernel2(NSpin,NAOrbs,NAOrbs), &
                STAT=AllocErr )
      if( AllocErr /= 0 ) then
        print *, "DEVICE/Allocation error for CGibbsY"
        stop
      end if
      allocate( PDGIBBS(NSpin,NAOrbs,NAOrbs), &
                PDOUTGIBBS(NSpin,NAOrbs,NAOrbs), &
                STAT=AllocErr )
      if( AllocErr /= 0 ) then
        print *, "DEVICE/Allocation error for PDGIBBS"
        stop
      end if
    end if    

    SD = S
    call RMatPow( SD, -1.0d0, InvSD )

    if( HybFunc ) call InitCorrelation(NAOrbs,NSpin)

    allocate( SPH(NAorbs,NAOrbs), STAT=AllocErr )
    if( AllocErr /= 0 ) then
       print *, "DEVICE/InitDevice: Allocation error for SPH(:,:)"
       stop
    end if
    ! Computing transformation matrix S^+1/2
    call RMatPow( SD,  0.5d0, RSPH )
    SPH = RSPH

    shift = -FermiStart
    shiftup = shift
    shiftdown = shift

    IF( ElType(1) == "BETHE" .and. ElType(2) == "BETHE" ) THEN 
      call AnalyseCluster
    ELSE IF  (ElType(1) == "BETHE" .and. ElType(2) == "GHOST" ) THEN
      call AnalyseClusterElectrodeOne
    ELSE IF  (ElType(1) == "GHOST" .and. ElType(2) == "BETHE" ) THEN
      call AnalyseClusterElectrodeTwo
    ELSE IF  (ElType(1) == "GHOST" .and. ElType(2) == "GHOST" ) THEN
      continue                           
    ELSE 
      print *, 'These electrodes are not implemented yet !!!'
      stop
    END IF

    call InitElectrodes

    EMin = 0.0d0
    EMax = 0.0d0

    LEV = 0.0d0
    HEV = 0.0d0

    ! Compute number of electrons 
    ! in central device region
    

    IF (Overlap < 0.01) THEN
       NEmbed1=0
       NEmbed2=0
    ELSE
       NEmbed1=NEmbedBL(1)
       NEmbed2=NEmbedBL(2)
    END IF

   
    print *, "---------------------------------------------------"
    print *, " Details on device and contacting atoms -----------"
    print *, "---------------------------------------------------"

    print *, "NEmbed(1) =", NEmbedBL(1)
    print *, "NEmbed(2) =", NEmbedBL(2)

    NCDEl = 0
    do iatom=NEmbed1+1,GetNAtoms()-NEmbed2
      NCDEl = NCDEl + GetAtmChg(iatom)
    end do
    
    print *, "Number of electrons in neutral reduced device"
    print *, "NCDEl = ", NCDEl
    
    print *, "First and last orbital in reduced device"
    IF (Overlap < 0.01) THEN
       NCDAO1 = 1
    ELSE
       NCDAO1 = 0
    END IF
    do iatom=1,NEmbed1
       NCDAO1 = NCDAO1 + NAOAtom(iatom) 
    end do

    print *, "NCDAO1 = ", NCDAO1

    NCDAO2 = NCDAO1-1
    do iatom = NEmbed1+1,GetNAtoms()-NEmbed2
       NCDAO2 = NCDAO2 + NAOAtom(iatom) 
    end do

    IF  (ElType(1) == "GHOST" .and. ElType(2) == "GHOST" ) NCDAO2 = NAOrbs
    print *, "NCDAO2 = ", NCDAO2
    print *, "---------------------------------------------------"

  end subroutine InitDevice
  
  !***********************************************
  !* Read initial density matrix from file P.dat *
  !***********************************************
  subroutine ReadDensMat(densitymatrix)
    use parameters, only: NSpinEdit, SpinEdit, MRStart, SpinDel, PFix, NFix, IFix, densitymatrixx
    use constants, only: d_zero
    use numeric, only: RMatPow
    use cluster, only: LoAOrbNo, HiAOrbNo,NAOAtom
#ifdef G03ROOT
    use g03Common, only: GetNE, GetNAtoms, GetNAE, GetNBE
#endif
#ifdef G09ROOT
    use g09Common, only: GetNE, GetNAtoms, GetNAE, GetNBE
#endif
    use ANTCommon
    implicit none
    
    integer :: norb, ni, nj, isp, n, iatom, jatom, is, i, j, ii, jj, AOStart, AO_BEG, AO_END, iorb, jorb, iato, jato
    real*8 :: density, TrP, xxx
    real*8, dimension(NSpin,NAOrbs,NAOrbs) :: PDMod
    real*8, dimension(NAOrbs,NAOrbs) :: PDAux
    character (len=50) :: densitymatrix
  
    write(ifu_log,*) "-------------------------------------------------------------------------------------------"
    write(ifu_log,*) "Starting calculation with density matrix from file  ", densitymatrix
    write(ifu_log,*) "-------------------------------------------------------------------------------------------"

    PD=0.0d0
    PDMod=0.0d0

    open(ifu_dm,file=densitymatrix,status='old')
    read(ifu_dm,*,end=21) shift
    shift = -shift
    do i=1,2*NAOrbs*NAOrbs
       read(ifu_dm,*,end=21)is,ii,jj,density
       PDMod(is,ii,jj)=density
       PDMod(is,jj,ii)=density
    end do
 21 close (ifu_dm)

    if (PFix) then
       write(ifu_log,*) "-------------------------------------------------------------------------------------------"
       write(ifu_log,*) "  ... and using supplementary density matrix from file  ", densitymatrixx
       write(ifu_log,*) "-------------------------------------------------------------------------------------------"
       if (NFix == 0) print *,'Warning ... NFix = 0'
       open(ifu_dmx,file=densitymatrixx,status='old')
       read(ifu_dmx,*) xxx   
       do norb=1,8*NAOrbs*NAOrbs
          read(ifu_dmx,*,end=22)is,ni,nj,density,iato,iorb,jato,jorb
          do isp=1,NSpin
             i=0
             do iAtom=1,GetNAtoms()
                do n=1,NFix
                   if (iAtom == IFix(n)) then 
                      i=i+NAOAtom(iAtom)
                      goto 11
                   end if
                end do
                do ii=1,NAOAtom(iAtom)
                   i=i+1
                   j=0
                   do jAtom=1,GetNAtoms()
                      do n=1,NFix
                         if (jAtom == IFix(n)) then 
                            j=j+NAOAtom(jAtom)
                            goto 12
                         end if
                      end do
                      do jj=1,NAOAtom(jAtom)
                         j=j+1
                         if (is == isp .and. iato == iAtom .and. jato == jAtom .and. iorb == ii .and. jorb == jj) then
                            PDMod(isp,i,j)=density
                            PDMod(isp,j,i)=density
                         end if
                      end do
12                 end do
                end do
11           end do
          end do
       end do
22     close(ifu_dmx)
       
    end if
  
    !open(111,file='readdensitymatrix',status='unknown')
    !do is=1,NSpin
    !i=0
    !do iAtom=1,GetNAtoms()
    ! do ii=1,NAOAtom(iAtom)
    !    i=i+1
    !    j=0
    !    do jAtom=1,GetNAtoms()
    !       do jj=1,NAOAtom(jAtom)
    !          j=j+1
    !          write(111,'(i3,4i5,e18.10)')is,iAtom,ii,jAtom,jj,PDMod(is,i,j)
    !end do
    !end do
    !end do
    !end do
    !end do
    !close(111)

    !
    ! Manipulate atomic spins of initial guess 
    ! if SpinEdit or MRStart set
    !
    if( NSpinEdit > 0 .or. MRStart > 0 .or. SpinDel )then

       if( MRStart > 0 )then
          do iAtom=MRStart,GetNAtoms()
             SpinEdit( iAtom ) = -1
          end do
       end if

       if( SpinDel )then
          do iAtom=1,GetNAtoms()
             SpinEdit( iAtom ) = 0
          end do
       end if

       PD = PDMod
       PDMod = d_zero

       do iAtom=1,GetNAtoms()
          do jAtom=1,GetNAtoms()
             if( SpinEdit(iAtom) ==  1 .and. SpinEdit(jAtom) ==  1 )then
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)
                      PDMod(1,i,j) = PD(1,i,j)
                      PDMod(2,i,j) = PD(2,i,j)
                   end do
                end do
             else if( SpinEdit(iAtom) == -1 .and. SpinEdit(jAtom) == -1 )then
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)                
                      PDMod(2,i,j) = PD(1,i,j)
                      PDMod(1,i,j) = PD(2,i,j)
                   end do
                end do
             else
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)                
                      PDMod(1,i,j) = 0.5d0*(PD(1,i,j)+PD(2,i,j))            
                      PDMod(2,i,j) = 0.5d0*(PD(1,i,j)+PD(2,i,j))            
                   end do
                end do
             end if
          end do
       end do

    end if

    ! Normalize to correct number of electrons when recommended
     
    TrP = d_zero
    do is=1,NSpin
       PDAux=MATMUL(PDMod(is,:,:),SD)
       do i=1,NAOrbs
          TrP = TrP + PDAux(i,i)
       end do
    end do
    if (NSpin ==1) then
       PRINT *, "Tr[P*S] of initial guess =", TrP*2.0
       if( NSpinEdit > 0 .or. MRStart > 0 .or. PFIX) then
         PD = PDMod*GetNE()/(TrP*2.0)
       else
         PD = PDMod
       end if
    else
       PRINT *, "Tr[P*S] of initial guess =", TrP
       if( NSpinEdit > 0 .or. MRStart > 0 .or. PFIX) then
          PD = PDMod*GetNE()/(TrP)
       else
          PD = PDMod
       end if
    end if
    PRINT *, "--------------------------------"
       
  end subroutine ReadDensMat

  !********************************************
  !* Read initial Fock matrix from file F.dat *
  !********************************************
  subroutine ReadFockMat(fockmatrix)
    use parameters, only: NSpinEdit, SpinEdit, MRStart, SpinDel !!!, PFix, NFix, IFix, densitymatrixx
    use constants, only: d_zero
    use numeric, only: RMatPow
    use cluster, only: LoAOrbNo, HiAOrbNo,NAOAtom
#ifdef G03ROOT
    use g03Common, only: GetNE, GetNAtoms, GetNAE, GetNBE
#endif
#ifdef G09ROOT
    use g09Common, only: GetNE, GetNAtoms, GetNAE, GetNBE
#endif
    use ANTCommon
    implicit none
    
    integer :: norb, ni, nj, isp, n, iatom, jatom, is, i, j, ii, jj, AOStart, AO_BEG, AO_END, iorb, jorb, iato, jato
    real*8 :: fock !, TrP, xxx

    real*8, dimension(NSpin,NAOrbs,NAOrbs) :: HDMod
    character (len=50) :: fockmatrix
  
    write(ifu_log,*) "-------------------------------------------------------------------------------------------"
    write(ifu_log,*) "Starting calculation with fock matrix from file  ", fockmatrix
    write(ifu_log,*) "-------------------------------------------------------------------------------------------"

    HD=0.0d0
    HDMod=0.0d0

    open(ifu_fm,file=fockmatrix,status='old')
    read(ifu_fm,*,end=21) shift
    shift = -shift
    do i=1,2*NAOrbs*NAOrbs
       read(ifu_fm,*,end=21)is,ii,jj,fock
       HD(is,ii,jj)=fock
       HD(is,jj,ii)=fock
    end do
21 close (ifu_fm)

    !
    ! Manipulate atomic spins of initial guess 
    ! if SpinEdit or MRStart set
    !
    if( NSpinEdit > 0 .or. MRStart > 0 .or. SpinDel )then

       if( MRStart > 0 )then
          do iAtom=MRStart,GetNAtoms()
             SpinEdit( iAtom ) = -1
          end do
       end if

       if( SpinDel )then
          do iAtom=1,GetNAtoms()
             SpinEdit( iAtom ) = 0
          end do
       end if

       HDMod = d_zero

       do iAtom=1,GetNAtoms()
          do jAtom=1,GetNAtoms()
             if( SpinEdit(iAtom) ==  1 .and. SpinEdit(jAtom) ==  1 )then
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)
                      HDMod(1,i,j) = HD(1,i,j)
                      HDMod(2,i,j) = HD(2,i,j)
                   end do
                end do
             else if( SpinEdit(iAtom) == -1 .and. SpinEdit(jAtom) == -1 )then
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)                
                      HDMod(2,i,j) = HD(1,i,j)
                      HDMod(1,i,j) = HD(2,i,j)
                   end do
                end do                
             else
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)                
                      HDMod(1,i,j) = 0.5d0*(HD(1,i,j)+HD(2,i,j))            
                      HDMod(2,i,j) = 0.5d0*(HD(1,i,j)+HD(2,i,j))            
                   end do
                end do
             end if
          end do
       end do

       HD = HDMod

    end if

  end subroutine ReadFockMat


  !**************************
  ! Deallocate dynamic arrays
  !**************************
  subroutine CleanUpDevice
    use BetheLattice, only: CleanUpBL, LeadBL
    use parameters, only: ElType, BiasVoltage, CompFock
    implicit none    
    integer :: AllocErr, LeadNo

    deallocate( SD, InvSD, HD, PD, PDOUT, SPH, HW, HWOUT, DGibbsY, CGibbsY, DGibbsYKernel1, CGibbsYKernel1, DGibbsYKernel2, CGibbsYKernel2, PDGIBBS, PDOUTGIBBS, STAT=AllocErr )
    if( AllocErr /= 0 ) then
       print *, "DEVICE/Deallocation error for SD, InvSD, HD, PD, PDOUT, SPH, HW, HWOUT, DGibbsY, CGibbsY, DGibbsYKernel1, CGibbsYKernel1, DGibbsYKernel2, CGibbsYKernel2, PDGIBBS, PDOUTGIBBS"
       stop
    end if
    do LeadNo=1,2
       select case( ElType(LeadNo) )
       case( "BETHE" )
          call CleanUpBL( LeadBL(LeadNo) ) 
       end select
    end do
  end subroutine CleanUpDevice


  !*************************
  !* Initialize electrodes *
  !*************************
  subroutine InitElectrodes
    use BetheLattice, only: InitBetheLattice, LeadBL, BL_EMin, BL_EMax
    use OneDLead, only: Init1DLead, Lead1d, L1D_EMin, L1D_EMax
    use parameters, only: ElType
    implicit none 
    integer :: LeadNo
    real*8, dimension(2) :: EMin, EMax

    do LeadNo=1,2
       select case( ElType(LeadNo) )
       case( "BETHE" )
          call InitBetheLattice( LeadBL(LeadNo), LeadNo )
          EMin(LeadNo) = BL_EMin( LeadBL(LeadNo) )
          EMax(LeadNo) = BL_EMax( LeadBL(LeadNo) )
       case( "1DLEAD" )
          call Init1DLead( Lead1d(LeadNo), LeadNo )
          EMin(LeadNo) = L1D_EMin( Lead1d(LeadNo) )
          EMax(LeadNo) = L1D_EMax( Lead1d(LeadNo) )
       case( "GHOST" )
          EMin(LeadNo) = -100.0                       
          EMax(LeadNo) =  100.0                       
       case DEFAULT
          print *, "Unknown option for ElType:", ElType(LeadNo)
          stop
       end select
    end do
    
    EMinEc = min( Emin(1),EMin(2) )
    EMaxEc = max( EMax(1),EMax(2) )
  end subroutine InitElectrodes

  
  !***************************
  !* Solve transport problem *
  !***************************
  subroutine Transport(F,ADDP) 
    use parameters, only: RedTransmB, RedTransmE, ANT1DInp, ElType, HybFunc, POrtho, DFTU, DiagCorrBl, DMImag, LDOS_Beg, LDOS_End, &
                          NSpinEdit, SpinEdit, SOC, ROT, PrtHatom, IntEnergy, BiasEnergy, DiagFock, SPINMU, BiasVoltage, CompFock, CompGibbsY
    use numeric, only: RMatPow, RSDiag
    use cluster, only: LoAOrbNo, HiAOrbNo
    use correlation
    use orthogonalization
    use constants, only: Hart    
#ifdef G03ROOT
    use g03Common, only: GetNAtoms
#endif
#ifdef G09ROOT
    use g09Common, only: GetNAtoms
#endif    
    implicit none

    logical,intent(out) :: ADDP
    real*8, dimension(NSpin,NAOrbs,NAOrbs),intent(in) :: F

    real*8, dimension(NAOrbs) :: evals
    real*8, dimension(:,:),allocatable :: SPM

    real*8 :: diff !!,TrP,QD
    integer :: i,j,is, info, AllocErr, iatom, jatom, Atom

    HD = F
    !
    ! Estimate upper bound for maximal eigenvalue of HD and use it for upper and lower energy boundaries
    !
    if( NSpin == 1 ) EMax = maxval(sum(abs(HD(1,:,:)),1))
    if( NSpin == 2 ) EMax = max(maxval(sum(abs(HD(1,:,:)),1)),maxval(sum(abs(HD(2,:,:)),1)))
    EMin = -EMax
    print *, "EMin=", EMin
    print *, "EMax=", EMax

    if( .not. DMImag .and. ChargeCntr )then
       print *, "--------------"
       print *, "Charge Control"
       print *, "--------------"
       ! Find upper and lower energy bound 
       call FindEnergyBounds
    end if

    if( DFTU ) call Add_DFT_plus_U_Pot( PD, HD )

    if(.not.DMImag) call CompDensMat(ADDP,.false.,.false.)
    if(DMImag) call CompDensMat2(ADDP)

    if( Evaluation )then
       print *
       print *, "****************************************** "
       print *, "*                                        * "
#ifdef G03ROOT
       print *, "*        ANT.G03 final analysis          * "
#endif
#ifdef G09ROOT
       print *, "*        ANT.G09 final analysis          * "
#endif
       print *, "*                                        * "       
       print *, "****************************************** "
       print *
       
      if( SPINMU )then
        !call Hamiltonian
        !call LDOS
        !call MullPop
        if(.not.DMImag) call CompLocalMu(ADDP)
        if(DMImag) call CompLocalMu(ADDP)
      end if
      !if( IntEnergy ) call WorkEnergy
      if( DiagFock ) call WorkFock
      if( IntEnergy ) call WorkEnergy
      if( BiasEnergy ) call WorkBiasEnergy
      if( CompFock )then
      !if( CompGibbsY )then
        if(BiasVoltage==0.0d0)then
          ! 1st TRUE MEANS COMPUTE FOCK WITH QHWTot. 2nd FALSE MEANS NOT COMPUTE CGibbsY.
          !call CompDensMat(ADDP,.true.,.false.)
          call BuildLiouvillian ! COMMENTED ON 2018-04-24 BECAUSE LIOUVILIAN IS NOT USED ANYMORE.
          call CompDensMat(ADDP,.true.,.true.)
          call DeallocateLiouvillian ! COMMENTED ON 2018-04-24 BECAUSE LIOUVILIAN IS NOT USED ANYMORE.
        elseif(BiasVoltage/=0.0d0)then
          call BuildLiouvillian ! COMMENTED ON 2018-04-24 BECAUSE LIOUVILIAN IS NOT USED ANYMORE.
          ! 1st TRUE MEANS COMPUTE FOCK WITH QHWTot. 2nd TRUE MEANS COMPUTE CGibbsY.
          call CompDensMat(ADDP,.true.,.true.)
          call DeallocateLiouvillian ! COMMENTED ON 2018-04-24 BECAUSE LIOUVILIAN IS NOT USED ANYMORE.
        end if
      end if       

       IF( ANT1DInp ) call WriteANT1DInput

       if( POrtho )then
          allocate( OD(NAorbs,NAOrbs), STAT=AllocErr )
          if( AllocErr /= 0 ) then
             print *, "DEVICE/InitDevice: Allocation error for OD(:,:)"
             stop
          end if
          do ispin=1,NSpin
             PD(ispin,:,:) = matmul( SD, matmul(PD(ispin,:,:), SD ) )
          end do
          HDOrtho = .true.
          call ProjOrtho(cix, SD, OD )
          ! Othogonalize density matrix density matrix and Hamiltonian
          do ispin=1,NSpin
             HD(ispin,:,:) = matmul( transpose(OD), matmul(F(ispin,:,:), OD) )
             PD(ispin,:,:) = matmul( OD, matmul(PD(ispin,:,:), transpose(OD) ) )
          end do
          HDOrtho = .true.
       end if
       if( DiagCorrbl ) call DiagCorrBlocks( HD, SD )
       call Hamiltonian
       IF ( HybFunc ) call CompHybFunc
       IF ((ElType(1) == "GHOST" .or. ElType(2) == "GHOST") .and. LDOS_Beg <= LDOS_End) CALL LDOS
       IF (ElType(1) /= "GHOST" .and. ElType(2) /= "GHOST") THEN            
          IF( RedTransmE >= RedTransmB  ) call EigenChannelAnalysis
          call transmission
       END IF
       
       !
       ! Print hamiltonian of atoms 
       ! with SpinEdit set
       !
       if( NSpinEdit > 0 )then                         
       
       Atom = PrtHatom
 
          do iAtom=1,GetNAtoms()
             do jAtom=1,GetNAtoms()
                !if( SpinEdit(iAtom) == -1 .and. SpinEdit(jAtom) == -1 )then
                if( iAtom == Atom .and. jAtom == Atom )then
                   PRINT *, " Hamiltonian of edited spin (atom) ",iAtom," is: "
                   PRINT *, " Up-Up "  
                   do i=LoAOrbNo(Atom),HiAOrbNo(Atom)                                                               
                      PRINT '(1000(F11.5))', ( (HD(1,i,j)), j=LoAOrbNo(Atom),HiAOrbNo(Atom) )               
                   end do                                                                                                       
                   PRINT *, " Down-Down "                                                                              
                   do i=LoAOrbNo(Atom),HiAOrbNo(Atom)                                                            
                      PRINT '(1000(F11.5))', ( (HD(2,i,j)), j=LoAOrbNo(Atom),HiAOrbNo(Atom) )             
                   end do                                                                              
                end if   
             end do
          end do
       end if   
        
       if (SOC .or. ROT) then 
          call MullPop_SOC
       else 
          call MullPop
       end if   
 
    end if
  end subroutine transport

  subroutine BuildLiouvillian
!    use g09Common, only: GetNAE, GetNBE
    use parameters, only: eta
    use constants, only: d_zero, d_one, c_zero, c_one, ui, eleccharge, hbar
    use numeric, only: CInv
    implicit none

!    complex, dimension(:,:,:),allocatable :: LiouvSOp!, LiouvSOpL, LiouvSOpR
    complex, dimension(NAOrbs,NAOrbs) :: IDD
    complex, dimension(NAOrbs*NAOrbs,NAOrbs*NAOrbs) :: LiouvSOpPiv
!    real ::
!    real ::
    integer :: iSpin, i, j, k, l
    integer :: info, allocerr
    integer :: n, ipiv(NAOrbs*NAOrbs)
    complex*16, DIMENSION( 4*NAOrbs*NAOrbs ) :: work

    if(DebugDev)then
    Write(*,'(A)')"***********************************"
    Write(*,'(A)')"**** ENTER BuildLiouvillian!!! ****"
    Write(*,'(A)')"***********************************"
    end if

    allocate( LiouvSOp(NSpin,NAorbs*NAorbs,NAorbs*NAorbs),&
              STAT=AllocErr )
    if( AllocErr /= 0 ) then
      print *, "DEVICE/BuildLiouvillian: Allocation error for LiouvSOp(:,:,:)"
      stop
    end if
    LiouvSOpPiv = c_zero
    LiouvSOp = c_zero

    IDD = c_zero                           ! Initialize the array.
    forall(j = 1:NAOrbs) IDD(j,j) = c_one     ! Set the diagonal.

    Write(*,'(A)')"DON'T FORGET TO MULTIPLY BY THE ELECTRON CHARGE!!!"
    do iSpin=1,NSpin
      LiouvSOpPiv = c_zero
      do i=1,NAOrbs
        do j=1,NAOrbs
          do k=1,NAOrbs
            do l=1,NAOrbs
              !LiouvSOpL((i-1)*NAOrbs+k,(j-1)*NAOrbs+l)=HD(i,j)*IDD(k,l)
              !LiouvSOpR((i-1)*NAOrbs+k,(j-1)*NAOrbs+l)=IDD(i,j)*HD(k,l)
              ! FOR EACH SUPEROP ELEMENT, 1st SUMANDO BELONGS TO LiouvSOpL & 2nd TO LiouvSOpR, NEGATIVE.
              ! WITH THIS VALUE THE GibbsY OPERATOR RESULTS HERMITIAN.
              !LiouvSOpPiv((i-1)*NAOrbs+k,(j-1)*NAOrbs+l)= -(COMPLEX(HD(iSpin,i,j),0.0D0)*IDD(k,l) - IDD(i,j)*COMPLEX(HD(iSpin,k,l),0.0D0))
              LiouvSOpPiv((i-1)*NAOrbs+k,(j-1)*NAOrbs+l)= - HD(iSpin,i,j)*IDD(k,l) + IDD(i,j)*HD(iSpin,k,l) ! STRICTLY -L.
              ! ADD COMPLEX ETA TO THE DIAGONAL.
              if((i==j) .and. (k==l))then
                LiouvSOpPiv((i-1)*NAOrbs+k,(j-1)*NAOrbs+l) = LiouvSOpPiv((i-1)*NAOrbs+k,(j-1)*NAOrbs+l) + ui*eta ! BECAUSE -L+i*eta
              end if
              !Write(*,'(I4,I4,I4,F12.8,F12.8)')iSpin,i,j,LiouvSOp(iSpin,i,j)
            end do
          end do
        end do
      end do

      if(DebugDev)then
      if(iSpin==1)Write(*,'(A)')"ALPHA LIOUVILLIAN BEFORE INVERSION"
      if(iSpin==2)Write(*,'(A)')"BETA LIOUVILLIAN BEFORE INVERSION"
      call PrintCMatrix(LiouvSOpPiv)
      end if

      !info = Cinv(LiouvSOpPiv)
      n = SIZE( LiouvSOpPiv, 1)
      CALL zgetrf(n,n,LiouvSOpPiv,n,ipiv,info)
      CALL zgetri(n,LiouvSOpPiv,n,ipiv,work,4*n,info)
      if( info /= 0 ) THEN
        WRITE(ifu_log,*)'Device/BuildLiouvillian using CInv (zgetrf,zgetri) in device.f90'
        WRITE(ifu_log,*)'INFO=',info
!       STOP
      end if

      do i=1,NAOrbs*NAOrbs
        do j=1,NAOrbs*NAOrbs
          LiouvSOp(iSpin,i,j)=LiouvSOpPiv(i,j)
        end do
      end do

    end do

    if(DebugDev)then
    do iSpin=1,NSpin
      if(iSpin==1)Write(*,'(A)')"ALPHA INVERSE LIOUVILLIAN"
      if(iSpin==2)Write(*,'(A)')"BETA INVERSE LIOUVILLIAN"
      !do i=1,NAOrbs*NAOrbs
      !  do j=1,NAOrbs*NAOrbs
      !    Write(*,'(I4,I4,I4,A,F12.8,A,F12.8,A)')iSpin,i,j,"(",DREAL(LiouvSOp(iSpin,i,j)),") + i*(",DIMAG(LiouvSOp(iSpin,i,j)),")"
      !  end do
      !end do
      if(iSpin==1)Write(*,'(A)')"ALPHA LIOUVILLIAN AFTER INVERSION"
      if(iSpin==2)Write(*,'(A)')"BETA LIOUVILLIAN AFTER INVERSION"
      call PrintCMatrix(LiouvSOpPiv)
    end do
    end if

    if(DebugDev)then
    Write(*,'(A)')"***********************************"
    Write(*,'(A)')"**** EXIT BuildLiouvillian!!! *****"
    Write(*,'(A)')"***********************************"
    end if

  end subroutine BuildLiouvillian

  subroutine DeallocateLiouvillian
    implicit none
    integer :: allocerr
    !deallocate(LiouvSOp,stat=allocerr)
    if (allocerr /= 0 ) then
      print*,"Problems deallocating LiouvSOp"
      !stop
    end if
  end subroutine DeallocateLiouvillian

  subroutine WorkEnergy
    use g09Common, only: GetNAE, GetNBE
    use constants, only: d_zero, d_pi
    use parameters, only: RedTransmB, RedTransmE, ANT1DInp, ElType, HybFunc, POrtho, DFTU, DiagCorrBl, DMImag, LDOS_Beg, LDOS_End, DiagFock, IntEnergy, QEXCESS
    implicit none
    ! *** Eigenvalues of the Fock-Matrix ***
    real, dimension(:,:,:),allocatable :: OHD
    real :: upIntegerDOSE, downIntegerDOSE, IntegerDOSE, E1
    real :: hartreeupIntegerDOSE, hartreedownIntegerDOSE, hartreeIntegerDOSE,evperhartree
    integer :: i,j,is, info, AllocErr, NROT
    
    evperhartree = 2.721138D1

        !if(DMImag)then
        if(NSpin==2 .and. SPINLOCK )then
        !if(NSpin==2)then
          print'(A,f10.1)',' SPINLOCK: NAlpha-NBeta = ',dble((GetNAE()-GetNBE())*(NCDEl+QExcess))/dble(GetNAE()+GetNBE())
          do ispin=1,NSpin
            !if(shiftup/=shiftdown)then
              if (ispin.eq.1)then
                E1=shiftup
                write(ifu_log,*)'ispin',ispin
                write(ifu_log,*)'- SHIFTUP',-shiftup
              end if
              if (ispin.eq.2)then
                E1=shiftdown
                write(ifu_log,*)'ispin',ispin
                write(ifu_log,*)'- SHIFTDOWN',-shiftdown
              end if
              if (ispin.eq.1)upIntegerDOSE = FullEnergy(d_zero)
              if (ispin.eq.2)downIntegerDOSE = FullEnergy(d_zero)
            !end if
            write(ifu_log,*)'--------------------------------------------------------'
            if (ispin.eq.1) then
              write(ifu_log,'(A,F9.5)') ' Fermi energy for alpha electrons= ', -shiftup
              write(ifu_log,*)
              write(ifu_log,'(A,F16.8)') ' Energy of the alpha electrons: ', upIntegerDOSE
            end if
            if (ispin.eq.2) then
              write(ifu_log,'(A,F9.5)') ' Fermi energy for beta electrons=  ', -shiftdown
              write(ifu_log,*)
              write(ifu_log,'(A,F16.8)') ' Energy of the beta electrons: ', downIntegerDOSE
            end if
            write(ifu_log,*)'--------------------------------------------------------'
          end do
        else if(NSpin==2 .and. .not. SPINLOCK )then
          print'(A,f10.1)',' SPINLOCK: NAlpha-NBeta = ',dble((GetNAE()-GetNBE())*(NCDEl+QExcess))/dble(GetNAE()+GetNBE())
          ispin = 1
            E1=shift
            write(ifu_log,*)'ispin',ispin
            write(ifu_log,*)'- SHIFT',-shift
            IntegerDOSE = FullEnergy(d_zero)
            write(ifu_log,*)'--------------------------------------------------------'
            write(ifu_log,'(A,F9.5)') ' Fermi energy for alpha/beta electrons= ', -shift
            write(ifu_log,*)
            write(ifu_log,'(A,F16.8)') ' Energy of the alpha/beta electrons: ', IntegerDOSE
            write(ifu_log,*)'--------------------------------------------------------'

 !if(.not.DMImag) then
        !if(NSpin==1)then
        else
          print'(A,f10.1)',' SPINLOCK: NAlpha-NBeta = ',dble((GetNAE()-GetNBE())*(NCDEl+QExcess))/dble(GetNAE()+GetNBE())
          do ispin=1,NSpin
            E1=shift
            write(ifu_log,*)'ispin',ispin
            write(ifu_log,*)'- SHIFT',-shift
            IntegerDOSE = FullEnergy(d_zero)
            write(ifu_log,*)'--------------------------------------------------------'
            write(ifu_log,'(A,F9.5)') ' Fermi energy for alpha/beta electrons= ', -shift
            write(ifu_log,*)
            write(ifu_log,'(A,F16.8)') ' Energy of the alpha/beta electrons: ', IntegerDOSE
            write(ifu_log,*)'--------------------------------------------------------'
          end do
        end if
        if(NSpin==2 .and. SPINLOCK )then
          IntegerDOSE = upIntegerDOSE + downIntegerDOSE
          hartreeIntegerDOSE = IntegerDOSE/evperhartree
        else if(NSpin==2 .and. .not. SPINLOCK )then
          IntegerDOSE = 2.0d0*IntegerDOSE
          hartreeIntegerDOSE = IntegerDOSE/evperhartree
        else if(NSpin==1)then
          IntegerDOSE = 2.0d0*IntegerDOSE
          hartreeIntegerDOSE = IntegerDOSE/evperhartree
        end if
        write(ifu_log,*)
        write(ifu_log,*)'--------------------------------------------------------'
        !write(ifu_log,'(A,F10.5)') ' Total energy of electrons:  ', IntegerDOSE
        write(ifu_log,'(A,F16.8)') ' Total energy of electrons (eV):  ', IntegerDOSE
        write(ifu_log,'(A,F10.5)') ' Total energy of electrons (Hartree):  ', hartreeIntegerDOSE
        write(ifu_log,*)'--------------------------------------------------------'
  end subroutine WorkEnergy

  subroutine WorkBiasEnergy
    use parameters, only: RedTransmB, RedTransmE, ANT1DInp, ElType, HybFunc, POrtho, DFTU, DiagCorrBl, DMImag, LDOS_Beg, LDOS_End, DiagFock, BiasEnergy
    implicit none
    ! *** Eigenvalues of the Fock-Matrix ***
    real, dimension(:,:,:),allocatable :: OHD
    real :: upIntegerDOSE, downIntegerDOSE, IntegerDOSE
    real :: hartreeupIntegerDOSE, hartreedownIntegerDOSE, hartreeIntegerDOSE,evperhartree
    integer :: i,j,is, info, AllocErr, NROT
    
    evperhartree = 2.721138D1

 !if(.not.DMImag) then
        if(NSpin==1)then
            write(ifu_log,*)''
            write(ifu_log,*)''
            write(ifu_log,*)'*************************************************************'
            write(ifu_log,*)'ENERGY UP TO WHICH WE INTEGRATE DOS(E)*E'
            write(ifu_log,*)'- SHIFT',-shift
            !write(ifu_log,*)'SHIFT',shift
            write(ifu_log,*)'*************************************************************'
            write(ifu_log,*)'' 

            IntegerDOSE = QXTotEnergy(shift) ! THE - SIGN DOESN'T WORK.
            !IntegerDOSE = CompEnergy(shift) ! THE BEST ONE. INTEGRATES TO MU.
            !IntegerDOSE = CompEnergy(-shift) ! THE BEST ONE. INTEGRATES TO MU.

            hartreeIntegerDOSE = IntegerDOSE/evperhartree
            !print *, "Integral of DOS*E =", upIntegerDOSE
     write(ifu_log,*)'Integral of DOS*E (eV) =',IntegerDOSE
     write(ifu_log,*)'Integral of DOS*E (Hartree) =',hartreeIntegerDOSE
            write(ifu_log,*)''
        end if
        !if(DMImag)then
        if(NSpin==2)then
          write(ifu_log,*)'- SHIFTUP',-shiftup
          write(ifu_log,*)'- SHIFTDOWN',-shiftdown

          if(shiftup/=shiftdown)then
            write(ifu_log,*)'*************************************************************'
            write(ifu_log,*)'ENERGY UP TO WHICH WE INTEGRATE ALPHA DOS(E)*E'
            write(ifu_log,*)'- SHIFTUP',-shiftup
            !write(ifu_log,*)'SHIFTUP',shiftup
            write(ifu_log,*)'ENERGY UP TO WHICH WE INTEGRATE BETA DOS(E)*E'
            write(ifu_log,*)'- SHIFTDOWN',-shiftdown
            !write(ifu_log,*)'SHIFTDOWN',shiftdown
            write(ifu_log,*)'*************************************************************'
            write(ifu_log,*)''

            upIntegerDOSE = QXTotEnergy(shiftup) ! THE - SIGN DOESN'T WORK.
            !upIntegerDOSE = CompEnergy(shiftup) ! THE BEST ONE. INTEGRATES TO MU.
            !upIntegerDOSE = CompEnergy(-shiftup) ! THE BEST ONE. INTEGRATES TO MU.

            hartreeupIntegerDOSE = upIntegerDOSE/evperhartree
            !print *, "Integral of spin-up DOS*E =", upIntegerDOSE
     write(ifu_log,*)'Integral of spin-up DOS*E (eV) =',upIntegerDOSE
     write(ifu_log,*)'Integral of spin-up DOS*E (Hartree) =',hartreeupIntegerDOSE
            write(ifu_log,*)''

            downIntegerDOSE = QXTotEnergy(shiftdown) ! THE - SIGN DOESN'T WORK.
            !downIntegerDOSE = CompEnergy(shiftdown) ! THE BEST ONE. INTEGRATES TO MU.
            !downIntegerDOSE = CompEnergy(-shiftdown) ! THE BEST ONE. INTEGRATES TO MU.

     hartreedownIntegerDOSE = downIntegerDOSE/evperhartree
            !print *, "Integral of spin-down DOS*E =", downIntegerDOSE
     write(ifu_log,*)'Integral of spin-down DOS*E (eV) =',downIntegerDOSE
     write(ifu_log,*)'Integral of spin-down DOS*E (Hartree) =',hartreedownIntegerDOSE
            write(ifu_log,*)''

            IntegerDOSE = upIntegerDOSE + downIntegerDOSE
            hartreeIntegerDOSE = IntegerDOSE/evperhartree
          else
            write(ifu_log,*)''
            write(ifu_log,*)''
            write(ifu_log,*)'*************************************************************'
            write(ifu_log,*)'ENERGY UP TO WHICH WE INTEGRATE DOS(E)*E'
            write(ifu_log,*)'- SHIFT',-shift
            !write(ifu_log,*)'SHIFT',shift
            write(ifu_log,*)'*************************************************************'
            write(ifu_log,*)'' 

            IntegerDOSE = QXTotEnergy(shift) ! THE - SIGN DOESN'T WORK.
            !IntegerDOSE = CompEnergy(shift) ! THE BEST ONE. INTEGRATES TO MU.
            !IntegerDOSE = CompEnergy(-shift) ! THE BEST ONE. INTEGRATES TO MU.

            hartreeIntegerDOSE = IntegerDOSE/evperhartree
            !print *, "Integral of DOS*E =", upIntegerDOSE
     write(ifu_log,*)'Integral of DOS*E (eV) =',IntegerDOSE
     write(ifu_log,*)'Integral of DOS*E (Hartree) =',hartreeIntegerDOSE
            write(ifu_log,*)''
          end if
        end if

        write(ifu_log,*)'*************************************************************'
        write(ifu_log,*)''
 write(ifu_log,*)'Integral of total DOS*E (eV) =',IntegerDOSE
 write(ifu_log,*)'Integral of total DOS*E (Hartree) =',hartreeIntegerDOSE 
        write(ifu_log,*)''
        write(ifu_log,*)'*************************************************************'

  end subroutine WorkBiasEnergy

  subroutine WorkFock
    use parameters, only: RedTransmB, RedTransmE, ANT1DInp, ElType, HybFunc, POrtho, DFTU, DiagCorrBl, DMImag, LDOS_Beg, LDOS_End, DiagFock, IntEnergy, glue
    use numeric, only: RMatPow, RSDiag, CDiag, Jacobi, eigsrt, ordks, balanc, RInv, GAUSSJ
    use correlation
    use util
    use orthogonalization
    !use lapack_blas, only: zgetri, zgetrf
    implicit none
    ! *** Eigenvalues of the Fock-Matrix ***
    real, dimension(:,:,:),allocatable :: OHD
    complex*16, dimension(:,:,:),allocatable :: COHD
    real, dimension(:,:),allocatable :: kseigenvalues
    real, dimension(:),allocatable :: upkseigenvalues, DupOHD
    real, dimension(:),allocatable :: downkseigenvalues, DdownOHD
    real, dimension(:),allocatable :: allkseigenvalues, ordkseigenvalues, hartreekseigenvalues
    real, dimension(:,:),allocatable :: upHD, upOHD, invupOHD, pivupOHD, VupOHD
    real, dimension(:,:),allocatable :: downHD, downOHD, invdownOHD, pivdownOHD, VdownOHD
    complex*16, dimension(:),allocatable :: DupCOHD
    complex*16, dimension(:,:),allocatable :: upCOHD, VupCOHD
    complex*16, dimension(:),allocatable :: DdownCOHD
    complex*16, dimension(:,:),allocatable :: downCOHD, VdownCOHD
    real :: kohnshamenergy, upkohnshamenergy, downkohnshamenergy
    real :: hartreeks, hartreeupks, hartreedownks, evperhartree
    integer :: i,j,is, info, AllocErr, NROT
    integer :: upocckscount, downocckscount, norbks
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl,sigr
    complex*16 :: cshiftup, cshiftdown


! Dynamic arrays 
    allocate( OHD(NSpin,NAOrbs,NAOrbs), &
         kseigenvalues(NAOrbs,NAOrbs), &
         allkseigenvalues(2*NAOrbs), &
         hartreekseigenvalues(2*NAOrbs), &
         ordkseigenvalues(2*NAOrbs), &
         upkseigenvalues(NAOrbs), &
         downkseigenvalues(NAOrbs), &
         upHD(NAOrbs,NAOrbs), &
         downHD(NAOrbs,NAOrbs), &
         upOHD(NAOrbs,NAOrbs), &
         downOHD(NAOrbs,NAOrbs), &
         invupOHD(NAOrbs,NAOrbs), &
         invdownOHD(NAOrbs,NAOrbs), &
         pivupOHD(NAOrbs,NAOrbs), &
         pivdownOHD(NAOrbs,NAOrbs), &
         DupOHD(NAOrbs), &
         DdownOHD(NAOrbs), &
         VupOHD(NAOrbs,NAOrbs), &
         VdownOHD(NAOrbs,NAOrbs), &
         COHD(NSpin,NAOrbs,NAOrbs), &
         upCOHD(NAOrbs,NAOrbs), &
         downCOHD(NAOrbs,NAOrbs), &
         DupCOHD(NAOrbs), &
         DdownCOHD(NAOrbs), &
         VupCOHD(NAOrbs,NAOrbs), &
         VdownCOHD(NAOrbs,NAOrbs), &
         STAT=AllocErr )
    
    if( AllocErr /= 0 ) then
       print *, "DEVICE/Allocation error for SD, InvSD, SMH, SPH, H, P"
       stop
    end if
    
    evperhartree = 2.721138D1   
    upkohnshamenergy = 0.d0
    downkohnshamenergy = 0.d0
    !write(ifu_log,*)''
    !write(ifu_log,*)'******************************************************************'
    !write(ifu_log,*)'Before orthogonalization and diagonalization (looks yet symmetric)'
      !NAOrbs = SIZE(HD(1,:,:),2)
      !NBasis = SIZE(HD(1,:,:),2)
10111 format(a6,i4)
10112 format(f12.6)
!10111 format(a6,i4,f12.6)

        !write(ifu_log,*)'SNH'
        !call PrintCMatrix( SNH(1:4,1:4) )
        !write(ifu_log,*)'SPH'
        !call PrintCMatrix( SPH(1:4,1:4) )

        do ispin=1,NSpin
            !write(ifu_log,*)'HD'
            !call PrintRMatrix( HD(ispin,1:4,1:4) )
        end do

        do ispin=1,NSpin
           !OHD(ispin,:,:) = HD(ispin,:,:)
           ! Looks like this is the correct way.
           !OHD(ispin,:,:) = matmul( SNH, matmul( HD(ispin,:,:), SPH ) )
           ! Orthogonalization by SNH*HD*SNH
           OHD(ispin,:,:) = matmul( SNH, matmul( HD(ispin,:,:), SNH ) ) ! This seems to be the correct expression.
           !OHD(ispin,:,:) = HD(ispin,:,:)
           !OHD(ispin,:,:) = matmul( SNH, matmul( F(ispin,:,:), SPH ) )

           !OHD(ispin,:,:) = matmul( transpose(SD), matmul(HD(ispin,:,:), SD) )
           !OHD(ispin,:,:) = matmul( transpose(OD), matmul(F(ispin,:,:), OD) )
           !PD(ispin,:,:) = matmul( OD, matmul(PD(ispin,:,:), transpose(OD) ) )
        end do

        if(.not.DMIMag)then
          cshiftup = dcmplx(shift)
        else
          cshiftup = dcmplx(shiftup)
          cshiftdown = dcmplx(shiftdown)
        end if
        do ispin=1,NSpin
            !if ((ispin == 1)) upOHD=OHD(ispin,:,:)+sigl+sigr
            if ((ispin == 1))then
              call CompSelfEnergies( ispin, cshiftup, sigl, sigr )
              sigr = glue*sigr
              sigl = glue*sigl
              upCOHD=dcmplx(OHD(ispin,:,:))+sigl+sigr
            end if
            !if(DMImag)
            !if ((ispin == 2)) downOHD=OHD(ispin,:,:)+sigl+sigr
            if ((ispin == 2))then
              call CompSelfEnergies( ispin, cshiftdown, sigl, sigr )
              sigr = glue*sigr
              sigl = glue*sigl
              downCOHD=dcmplx(OHD(ispin,:,:))+sigl+sigr
            end if
        end do
    !write(ifu_log,*)'***************************************************************'
    !write(ifu_log,*)''
    !	write(ifu_log,*)'After orthogonalization'
    !write(ifu_log,*)'***************************************************************'
    !write(ifu_log,*)''

    do ispin=1,NSpin
      if ((ispin == 1))then
        !write(ifu_log,*)'Spin-up Hamiltonian'
        !call PrintRMatrix( upOHD(1:4,1:4) )
        !call PrintCMatrix( upCOHD(1:4,1:4) )
        !write(ifu_log,*)'***************************************************************'
      end if
      !if(DMImag)
      if ((ispin == 2))then
        !write(ifu_log,*)'Spin-down Hamiltonian'
        !call PrintRMatrix( downOHD(1:4,1:4) )
        !call PrintCMatrix( downCOHD(1:4,1:4) )
        !write(ifu_log,*)'***************************************************************'
      end if
    end do
    do ispin=1,NSpin
      if ((ispin == 1)) then
        !call balanc(upOHD)
        !call zgetrf(NAOrbs,NAOrbs,upCOHD,NAOrbs,ipiv,info)
        !write(ifu_log,*)'After balancing'
        !write(ifu_log,*)'***************************************************************'
        !write(ifu_log,*)''
        !write(ifu_log,*)'Spin-up Hamiltonian'
        !call PrintRMatrix( upOHD(1:4,1:4) )
        !call PrintCMatrix( upCOHD(1:4,1:4) )
      end if
      !if(DMImag)
      if ((Nspin == 2)) then
        !call balanc(downOHD)
        !call zgetrf(NAOrbs,NAOrbs,downCOHD,NAOrbs,ipiv,info)
        !write(ifu_log,*)'After balancing'
        !write(ifu_log,*)'***************************************************************'
        !write(ifu_log,*)''
        !write(ifu_log,*)'Spin-down Hamiltonian'
        !call PrintRMatrix( downOHD(1:4,1:4) )
        !call PrintCMatrix( downCOHD(1:4,1:4) )
      end if
    end do

    !write(ifu_log,*)'***************************************************************'
    !write(ifu_log,*)''

    !write(ifu_log,*)'*************************************************************'
    !write(ifu_log,*)''
    !write(ifu_log,*)'Before diagonalization'
    do ispin=1,NSpin
      if ((ispin == 1)) then
        call CDiag(upCOHD,DupCOHD,info)
        ! Sorting eigenvalues and eigenvectors in ascending order of eigenvalues.
        upOHD = real(upCOHD,8)
        DupOHD = real(DupCOHD,8)
        DupOHD = -DupOHD
        call eigsrt(DupOHD,VupOHD)
        DupOHD = -DupOHD
        !print *, ' DupOHD(j)'
        !print *, ( DupOHD(j), j=1,NAOrbs )
      end if
      !if(DMImag)
      if ((ispin == 2)) then
        call CDiag(downCOHD,DdownCOHD,info)
        ! Sorting eigenvalues and eigenvectors in ascending order of eigenvalues.
        downOHD = real(downCOHD,8)
        DdownOHD = real(DdownCOHD,8)
        DdownOHD = -DdownOHD
        call eigsrt(DdownOHD,VdownOHD)
        DdownOHD = -DdownOHD
        !print *, ' DdownOHD(j)'
        !print *, ( DdownOHD(j), j=1,NAOrbs )
      end if
    end do
    
    !if(.not.DMImag)
    if ((Nspin == 1)) DdownOHD = DupOHD !because beta and alpha eigenvalues are equal, overwrite the betas with the alphas. YOU DON'T NEED TO MULTIPLY LATER BY TWO.
    if ((Nspin == 1)) norbks = nint(2*QAlpha)
    if ((Nspin == 2)) norbks = nint(QAlpha+QBeta)

    do i=1,NAOrbs
      do j=i,NAOrbs
        do ispin=1,NSpin
          if ((ispin == 1)) upOHD(j,i)=upOHD(i,j)
          if ((ispin == 2)) downOHD(j,i)=downOHD(i,j)
        end do
      end do
    end do

    !write(ifu_log,*)'************************************************************'
    !write(ifu_log,*)''
    !write(ifu_log,*)'After diagonalization'
        
    allkseigenvalues(1:NAOrbs) = DupOHD(:)
    ordkseigenvalues(1:NAOrbs) = 1.0
    allkseigenvalues(NAOrbs+1:2*NAOrbs) = DdownOHD(:)
    ordkseigenvalues(NAOrbs+1:2*NAOrbs) = 2.0

    allkseigenvalues = -allkseigenvalues
    call ordks(allkseigenvalues,ordkseigenvalues)
    allkseigenvalues = -allkseigenvalues

    !do i=1,NAOrbs+1 !I first wrote this, don't remember why.

!    write(ifu_log,*)'************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Total num of K-S occ. orbitals',norbks
!    write(ifu_log,*)''
!    write(ifu_log,*)'K-S occ. eigen. (eV)',(allkseigenvalues(i),i=1,norbks)
!    write(ifu_log,*)''
!    write(ifu_log,*)'K-S virt. eigen. (eV)',(allkseigenvalues(i),i=norbks+1,2*NAOrbs)
!    write(ifu_log,*)''
!    write(ifu_log,*)'************************************************************'
!    write(ifu_log,*)''
!    !do i=1,2*NAOrbs
!    !   shiftallkseigenvalues(i) = allkseigenvalues(i) - shift
!    !end do
!    write(ifu_log,*)'************************************************************'
!    write(ifu_log,*)'SHIFT = ', -SHIFT
!    write(ifu_log,*)'************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Total num of K-S occ. orbitals',norbks
!    write(ifu_log,*)''
!    write(ifu_log,*)'K-S occ. eigen. (eV) (shifted)',(allkseigenvalues(i)-shift,i=1,norbks)
!    write(ifu_log,*)''
!    write(ifu_log,*)'K-S virt. eigen. (eV) (shifted)',(allkseigenvalues(i)-shift,i=norbks+1,2*NAOrbs)
!    write(ifu_log,*)''
!    write(ifu_log,*)'************************************************************'
!    write(ifu_log,*)''
    do i=1,2*NAOrbs
      hartreekseigenvalues(i) = allkseigenvalues(i)/27.21184
      !shifthartreekseigenvalues(i) = shiftallkseigenvalues(i)/27.21184
    end do

    write(ifu_log,*)'************************************************************'
    write(ifu_log,*)''
    write(ifu_log,*)'Total num of K-S occ. orbitals',norbks
    write(ifu_log,*)''
    write(ifu_log,*)'K-S occ. eigen. (Hartree)',(hartreekseigenvalues(i),i=1,norbks)
    write(ifu_log,*)''
    !write(ifu_log,*)'K-S virt. eigen. (Hartree)',(hartreekseigenvalues(i),i=norbks+1,2*NAOrbs)
    !write(ifu_log,*)''
    write(ifu_log,*)'************************************************************'
    write(ifu_log,*)''
    write(ifu_log,*)'************************************************************'
    write(ifu_log,*)'SHIFT = ', -SHIFT
    write(ifu_log,*)'************************************************************'
    write(ifu_log,*)''
    write(ifu_log,*)'Total num of K-S occ. orbitals',norbks
    write(ifu_log,*)''
    write(ifu_log,*)'K-S occ. eigen. (Hartree) (shifted)',(hartreekseigenvalues(i)-shift/evperhartree,i=1,norbks)
    write(ifu_log,*)''
    !write(ifu_log,*)'K-S virt. eigen. (Hartree) (shifted)',(hartreekseigenvalues(i)-shift/evperhartree,i=norbks+1,2*NAOrbs)
    !write(ifu_log,*)''
    write(ifu_log,*)'************************************************************'
    write(ifu_log,*)''

    upocckscount = 0
    downocckscount = 0
    do i=1,norbks
      if(ordkseigenvalues(i).eq.1.0) upocckscount = upocckscount + 1
      if(ordkseigenvalues(i).eq.2.0) downocckscount = downocckscount + 1
      IF (upocckscount + downocckscount .GE. norbks) EXIT
    end do

      !do i=1,NAOrbs
      !  !if(DupOHD(i).le.0.0d0)then
      !  if(i.le.NAOrbs/2)then
      !    upkohnshamenergy=upkohnshamenergy+DupOHD(i)
      !  end if
      !  !if(DdownOHD(i).le.0.0d0)then
      !  if(i.le.NAOrbs/2)then
      !    downkohnshamenergy=downkohnshamenergy+DdownOHD(i)
      !  end if
      !end do

    !do i=1,2*NAOrbs
    !  shiftDupOHD(i) = DupOHD(i)
    !  shiftDdownOHD(i) = DdownOHD(i)
    !end do

    upkohnshamenergy = sum(DupOHD(1:upocckscount))
    downkohnshamenergy = sum(DdownOHD(1:downocckscount))

    !shiftupkohnshamenergy = upkohnshamenergy - upocckscount*shift
    !shiftdownkohnshamenergy = downkohnshamenergy - downocckscount*shift

    !kohnshamenergy = upkohnshamenergy + downkohnshamenergy
    kohnshamenergy = sum(allkseigenvalues(1:norbks))
    hartreeupks = upkohnshamenergy/evperhartree
    hartreedownks = downkohnshamenergy/evperhartree
    hartreeks = kohnshamenergy/evperhartree
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Occupied alpha eigenvalues',upocckscount
!    !print '(1000(ES14.4))', ( DupOHD(j), j=1,upocckscount )
!    print *, ( DupOHD(j), j=1,upocckscount )
!    write(ifu_log,*)'up Kohn-Sham energy sum',upkohnshamenergy
!    write(ifu_log,*)'upKS in Hartrees',hartreeupks
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Occupied beta eigenvalues',downocckscount
!    !print '(1000(ES14.4))', ( DdownOHD(j), j=1,downocckscount )
!    print *, ( DdownOHD(j), j=1,downocckscount )
!    write(ifu_log,*)'down Kohn-Sham energy sum',downkohnshamenergy
!    write(ifu_log,*)'downKS in Hartrees',hartreedownks
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Total Kohn-Sham eigenvalues',norbks
!    write(ifu_log,*)'Total Kohn-Sham energy sum',kohnshamenergy
!    write(ifu_log,*)'KS in Hartrees',hartreeks
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!
!    write(ifu_log,*)'************************************************************'
!    write(ifu_log,*)'SHIFT = ', -SHIFT
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Occupied alpha eigenvalues (shifted)',upocckscount
!    !print '(1000(ES14.4))', ( DupOHD(j)-shift, j=1,upocckscount )
!    print *, ( DupOHD(j)-shift, j=1,upocckscount )
!    write(ifu_log,*)'up Kohn-Sham energy sum (shifted)',upkohnshamenergy-upocckscount*shift
!    write(ifu_log,*)'upKS in Hartrees (shifted)',hartreeupks-upocckscount*shift/evperhartree
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Occupied beta eigenvalues (shifted)',downocckscount
!    !print '(1000(ES14.4))', ( DdownOHD(j)-shift, j=1,downocckscount )
!    print *, ( DdownOHD(j)-shift, j=1,downocckscount )
!    write(ifu_log,*)'down Kohn-Sham energy sum (shifted)',downkohnshamenergy-downocckscount*shift
!    write(ifu_log,*)'downKS in Hartrees (shifted)',hartreedownks-downocckscount*shift/evperhartree
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Total Kohn-Sham eigenvalues',norbks
!    write(ifu_log,*)'Total Kohn-Sham energy sum (shifted)',kohnshamenergy-norbks*shift
!    write(ifu_log,*)'KS in Hartrees (shifted)',hartreeks-norbks*shift/evperhartree
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''

    write(ifu_log,*)'****************************************************************'
    write(ifu_log,*)''
    write(ifu_log,*)'Occupied alpha eigenvalues in Hartrees',upocckscount
    !print '(1000(ES14.4))', ( DupOHD(j), j=1,upocckscount )
    print *, ( DupOHD(j)/evperhartree, j=1,upocckscount )
    write(ifu_log,*)'up Kohn-Sham energy sum',upkohnshamenergy
    write(ifu_log,*)'upKS in Hartrees',hartreeupks
    write(ifu_log,*)'****************************************************************'
    write(ifu_log,*)''
    write(ifu_log,*)'Occupied beta eigenvalues in Hartrees',downocckscount
    !print '(1000(ES14.4))', ( DdownOHD(j), j=1,downocckscount )
    print *, ( DdownOHD(j)/evperhartree, j=1,downocckscount )
    write(ifu_log,*)'down Kohn-Sham energy sum',downkohnshamenergy
    write(ifu_log,*)'downKS in Hartrees',hartreedownks
    write(ifu_log,*)'****************************************************************'
    write(ifu_log,*)''
    write(ifu_log,*)'Total Kohn-Sham eigenvalues',norbks
    write(ifu_log,*)'Total Kohn-Sham energy sum',kohnshamenergy
    write(ifu_log,*)'KS in Hartrees',hartreeks
    write(ifu_log,*)'****************************************************************'
    write(ifu_log,*)''

!    write(ifu_log,*)'************************************************************'
!    write(ifu_log,*)'SHIFT = ', -SHIFT
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Occupied alpha eigenvalues in Hartrees (shifted)',upocckscount
!    !print '(1000(ES14.4))', ( DupOHD(j)-shift, j=1,upocckscount )
!    print *, ( (DupOHD(j)-shift)/evperhartree, j=1,upocckscount )
!    write(ifu_log,*)'up Kohn-Sham energy sum (shifted)',upkohnshamenergy-upocckscount*shift
!    write(ifu_log,*)'upKS in Hartrees (shifted)',hartreeupks-upocckscount*shift/evperhartree
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Occupied beta eigenvalues in Hartrees (shifted)',downocckscount
!    !print '(1000(ES14.4))', ( DdownOHD(j)-shift, j=1,downocckscount )
!    print *, ( (DdownOHD(j)-shift)/evperhartree, j=1,downocckscount )
!    write(ifu_log,*)'down Kohn-Sham energy sum (shifted)',downkohnshamenergy-downocckscount*shift
!    write(ifu_log,*)'downKS in Hartrees (shifted)',hartreedownks-downocckscount*shift/evperhartree
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Total Kohn-Sham eigenvalues',norbks
!    write(ifu_log,*)'Total Kohn-Sham energy sum (shifted)',kohnshamenergy-norbks*shift
!    write(ifu_log,*)'KS in Hartrees (shifted)',hartreeks-norbks*shift/evperhartree
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''



      !do ispin=1,NSpin
      !  do j=1,6
      !    if(NSpin==1) write(ifu_log,10111)'Elem:',j
      !    if(NSpin==1) write(ifu_log,10112)(upOHD(j,i),i=1,3)
      !    if(NSpin==2) write(ifu_log,10111)'Elem:',j
      !    if(NSpin==2) write(ifu_log,10112)(downOHD(j,i),i=1,3)
      !  end do
      !end do

      !call DiagFullHamiltonian( HD, SD)
      !call DiagHamiltonian
      deallocate( OHD, kseigenvalues, allkseigenvalues, upkseigenvalues, downkseigenvalues, upHD, downHD, upOHD, downOHD, DupOHD, DdownOHD, VupOHD, VdownOHD )
      !deallocate( shiftallkseigenvalues, shiftupkseigenvalues, shiftdownkseigenvalues, shiftDupOHD, shiftDdownOHD )
  end subroutine WorkFock

  subroutine PrintCMatrix( A )
    implicit none
    complex*16,dimension(:,:),intent(in) :: A
    integer :: i,j,dim1,dim2
    dim1 = size(A,1)
    dim2 = size(A,2)
    do i=1,dim1
       !print '(1000(ES14.4))', ( A(i,j), j=1,dim2 )
       print '(100(g15.5,g15.5,2x))', ( real(A(i,j)),&
                                        AIMAG(A(i,j)), j=1,dim2 )
    end do
  end subroutine PrintCMatrix
  
  subroutine WriteANT1DInput
    use parameters, only: eta
    use AntCommon
    implicit none

    real*8 :: dsmall
    integer :: i,j

    CHARACTER(len=55) :: fname

    dsmall = eta
    
    fname='dev.'//trim(ant1dname)//'.dat'
    open(ifu_ant,file=fname,status='unknown')

    write(ifu_ant,'(A)'),       "&DevParams"
    write(ifu_ant,'(A,I1)'),    "NDSpin = ", NSpin
    write(ifu_ant,'(A,I4)'),    "NDAO = ", NAOrbs
    write(ifu_ant,'(A,I4)'),    "NDEl = ", NCDEl
    write(ifu_ant,'(A,F12.8)'), "EFermi = ", -shift
    write(ifu_ant,'(A)'),       "sparse = .true."
    write(ifu_ant,'(A)'),       "/"
    write(ifu_ant,*)
    write(ifu_ant,'(A)'),       "! Hamiltonian"
    do ispin=1,NSpin
       if( NSpin == 2 .and. ispin == 1 ) write(ifu_ant,'(A)'),       "! Spin-up"
       if( NSpin == 2 .and. ispin == 2 ) write(ifu_ant,'(A)'),       "! Spin-down"
       do i=1,NAOrbs
          do j=1,NAOrbs
             if(abs(HD(ispin,i,j))>=dsmall) write(ifu_ant,'(I6,I6,ES20.8)'), i, j, HD(ispin,i,j)
          enddo
       enddo
       write(ifu_ant,'(I6,I6,ES20.8)'), 0, 0, 0.0d0
       write(ifu_ant,*)
    end do
    write(ifu_ant,'(A)'),       "! Overlap"
    do i=1,NAOrbs
       do j=1,NAOrbs
          if(abs(SD(i,j))>=dsmall) write(ifu_ant,'(I6,I6,ES20.8)'), i, j, SD(i,j)
       enddo
    enddo
    write(ifu_ant,'(I6,I6,ES20.8)'), 0, 0, 0.0d0
    write(ifu_ant,*)
    close(ifu_ant)
  end subroutine WriteANT1DInput
  
  !**************************************************************
  !* Subroutine for determining Fermi energy and density matrix *
  !* for some total charge                                      *
  !**************************************************************
  !* Pre-condition:                                             *
  !*   HC: Fock-matrix                                          *
  !*   shiftup,shiftdown: starting values for Fermi-energies    *
  !* Results:                                                   *
  !*   PC: Density-matrix                                       *
  !*   shiftup,shiftdown: Fermi-energies                        *
  !**************************************************************
  subroutine CompDensMat(ADDP, boolComputeFock, boolComputeCGibbsY)
    use numeric, only: Secant, Muller, BISEC
#ifdef G03ROOT
    use g03Common, only: GetNAE, GetNBE
#endif
#ifdef G09ROOT
    use g09Common, only: GetNAE, GetNBE
#endif    
    use parameters, only: FermiAcc,ChargeAcc,Max,QExcess
    !use ieee_arithmetic
    implicit none

    logical,intent(out) :: ADDP
    logical,intent(in) :: boolComputeFock, boolComputeCGibbsY
    integer :: i,j, k,cond
    real*8 :: E0,E1,E2,E3,DE,Z, Delta, Epsilon
    real*8 :: TotChargeFromGlesser
    logical :: root_fail
    
    Write(*,'(A)')"ENTERED Device/CompDensMat"

    Z=10.0d0*FermiAcc
    Delta=FermiAcc
    Epsilon=ChargeAcc*(NCDEl+QExcess)

    if( NSpin == 2 .and. SPINLOCK )then
       print'(A,f10.1)',' SPINLOCK: NAlpha-NBeta = ',dble((GetNAE()-GetNBE())*(NCDEl+QExcess))/dble(GetNAE()+GetNBE())
       do ispin=1,NSpin
          root_fail = .true.
          if (ispin.eq.1) E1=shiftup
          if (ispin.eq.2) E1=shiftdown
          E0=E1-Z
          E2=E1+Z
          if( root_fail )then
             print*,'MULLER method'
             call MULLER(F,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
             !call MULLER_OMP(F,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
             if(k.eq.Max .or. E3<EMin .or. E3>EMax) then
                print *, 'Warning: MULLER method failed to find root. Using SECANT.'
                root_fail = .true.
             else
                if (ispin.eq.1)shiftup=E3
                if (ispin.eq.2)shiftdown=E3
                root_fail = .false.
             end if
          end if
          if( root_fail )then
             print*,'SECANT method'
             call SECANT(F,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
             if(k.eq.Max .or. E3<EMin .or. E3>EMax) then
                print *, 'Warning: SECANT method failed to find root. Using BISEC.'
                root_fail = .true.
             else
                if (ispin.eq.1)shiftup=E3
                if (ispin.eq.2)shiftdown=E3
                root_fail = .false.
             end if
          end if
          if (root_fail) then
             print *, 'BISEC method'
             if (ispin.eq.1) shiftup = BISEC(F,EMin,EMax,Delta,5*Max,K)
             if (ispin.eq.2) shiftdown = BISEC(F,EMin,EMax,Delta,5*Max,K)
             DE=Delta
             if(k.lt.5*Max) root_fail = .false.
             if(k.ge.5*Max) print *, 'Warning: BISECT method failed to find root. Skipping this cycle.'
          end if
          write(ifu_log,*)'--------------------------------------------------------'
          if (ispin.eq.1) then
             write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for alpha electrons= ', -shiftup,'  +/-',dabs(DE)
             write(ifu_log,*)
             write(ifu_log,'(A,F10.5)') ' Number of alpha electrons: ', QAlpha
          end if
          if (ispin.eq.2) then
             write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for beta electrons=  ', -shiftdown,'  +/-',dabs(DE)
             write(ifu_log,*)
             write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', QBeta
             write(ifu_log,*)
             write(ifu_log,'(A,F10.5)') ' Total number of electrons:  ', QAlpha+QBeta
          end if
          write(ifu_log,*)'--------------------------------------------------------'
       end do
    else
       root_fail = .true.
       E0=shift-Z 
       E1=shift
       E2=shift+Z
       if (root_fail) then
          print*,'MULLER method'
          if(boolComputeFock)then
            call MULLER(QHWTot,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
          else
            call MULLER(QXTot,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
          end if
          !call MULLER_OMP(QXTot,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
          if(k .eq. Max .or. E2<EMin .or. E2>EMax) then
             print *, 'Warning: MULLER method failed to find root. Using SECANT.'
             root_fail = .true.
          else
             shift = E3
             root_fail = .false.
          end if
       end if
       if (root_fail) then
          print*,'SECANT method'
          if(boolComputeFock)then
            call SECANT(QHWTot,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
          else
            call SECANT(QXTot,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
          end if
          if(k .eq. Max .or. E3<EMin .or. E3>EMax) then
             print *, 'Warning: SECANT method failed to find root. Using BISEC.'
             root_fail = .true.
          else
             shift = E3
             root_fail = .false.
          end if
       end if
       if (root_fail) then
          print *, 'BISEC method'
          if(boolComputeFock)then
            shift = BISEC(QHWTot,EMin,EMax,Delta,5*Max,K)
          else
            shift = BISEC(QXTot,EMin,EMax,Delta,5*Max,K)
          end if
          DE=Delta
          if(k.lt.5*Max) root_fail = .false.
          if(k.ge.5*Max) print *, 'Warning: BISECT method failed to find root. Skipping this cycle.'
       end if

       if(boolComputeCGibbsY)then
         ! MUST NOT CALL MULLER METHOD RELYING ON  QCGibbsYTot BECAUSE IT USES GLESSER IN THE FULL ENERGY RANGE.
         !call MULLER(QCGibbsYTot,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
         write(ifu_log,*)''
         write(ifu_log,*)''
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)'--- BEFORE QCGibbsYTot(shift) ---'
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)'--- boolComputeCGibbsY IS true ---'
         write(ifu_log,*)'--- ENTERING QCGibbsYTot(shift) ---'
         TotChargeFromGlesser = QCGibbsYTot(shift)
         write(ifu_log,'(A,F10.5)') ' TotChargeFromGlesser:  ', TotChargeFromGlesser
         write(ifu_log,*)'--- AFTER QCGibbsYTot(shift) ---'
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)''
         write(ifu_log,*)''
       else
         write(ifu_log,*)''
         write(ifu_log,*)''
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)'--- boolComputeCGibbsY IS false ---'
         write(ifu_log,*)'--- NOT ENTERING QCGibbsYTot(shift) ---'
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)''
         write(ifu_log,*)''
         TotChargeFromGlesser = 0.0d0
       end if

       write(ifu_log,*)'-----------------------------------------------'
       write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy= ',-shift,'  +/-',dabs(DE)
       write(ifu_log,*)
       write(ifu_log,'(A,F10.5)') ' Number of alpha electrons: ', QAlpha
       write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', QBeta
       write(ifu_log,'(A,F10.5)') ' Total number of electrons:  ', QAlpha+QBeta
       write(ifu_log,*)'-----------------------------------------------'
    end if
    ADDP = .not. root_fail
  end subroutine CompDensMat
  
  
  !**************************************************************
  !* Subroutine for determining Fermi energy and density matrix *
  !* for some total charge                                      *
  !**************************************************************
  !* Pre-condition:                                             *
  !*   HC: Fock-matrix                                          *
  !*   shiftup,shiftdown: starting values for Fermi-energies    *
  !* Results:                                                   *
  !*   PC: Density-matrix                                       *
  !*   shiftup,shiftdown: Fermi-energies                        *
  !**************************************************************
!  subroutine CompDensMat(ADDP)
!    use numeric, only: SECANT, MULLER, BISEC
!#ifdef G03ROOT
!    use g03Common, only: GetNAE, GetNBE
!#endif
!#ifdef G09ROOT
!    use g09Common, only: GetNAE, GetNBE
!#endif
!    use parameters, only: FermiAcc,ChargeAcc,Max,QExcess
!    !use ieee_arithmetic
!    implicit none
!
!    logical,intent(out) :: ADDP
!    integer :: i,j, k,cond
!    real*8 :: E0,E1,E2,E3,DE,Z, Delta, Epsilon
!    
!    logical :: root_fail
!    
!    Z=0.1d0
!    Delta=FermiAcc
!    Epsilon=ChargeAcc*(NCDEl+QExcess)
!
!    if( NSpin == 2 .and. SPINLOCK )then
!       print'(A,f10.1)',' SPINLOCK: NAlpha-NBeta = ',dble((GetNAE()-GetNBE())*(NCDEl+QExcess))/dble(GetNAE()+GetNBE())
!       do ispin=1,NSpin
!          if (ispin == 1) print*,'Spin alpha'
!          if (ispin == 2) print*,'Spin beta'
!          root_fail = .true.
!          if (ispin.eq.1) E1=shiftup
!          if (ispin.eq.2) E1=shiftdown
!          E0=E1-Z
!          E2=E1+Z
!          if( root_fail )then
!             print*,'MULLER method'
!             call MULLER(QTot_SPINLOCK,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
!             if(k.eq.Max .or. E3<EMin .or. E3>EMax) then
!                print *, 'Warning: MULLER method failed to find root. Using SECANT.'
!                root_fail = .true.
!             else
!                if (ispin.eq.1)shiftup=E3
!                if (ispin.eq.2)shiftdown=E3
!                root_fail = .false.
!             end if
!          end if
!          if( root_fail )then
!             print*,'SECANT method'
!             call SECANT(QTot_SPINLOCK,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
!             if(k.eq.Max .or. E3<EMin .or. E3>EMax) then
!                print *, 'Warning: SECANT method failed to find root. Using BISEC.'
!                root_fail = .true.
!             else
!                if (ispin.eq.1)shiftup=E3
!                if (ispin.eq.2)shiftdown=E3
!                root_fail = .false.
!             end if
!          end if
!          if (root_fail) then
!             print *, 'BISEC method'
!             if (ispin.eq.1) shiftup = BISEC(QTot_SPINLOCK,EMin,EMax,Delta,5*Max,K)
!             if (ispin.eq.2) shiftdown = BISEC(QTot_SPINLOCK,EMin,EMax,Delta,5*Max,K)
!             DE=Delta
!             if(k.lt.5*Max) root_fail = .false.
!             if(k.ge.5*Max) print *, 'Warning: BISECT method failed to find root. Skipping this cycle.'
!          end if
!          write(ifu_log,*)'--------------------------------------------------------'
!          if (ispin.eq.1) then
!             write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for alpha electrons= ', -shiftup,'  +/-',dabs(DE)
!             write(ifu_log,*)
!             write(ifu_log,'(A,F10.5)') ' Number of alpha electrons: ', QAlpha
!          end if
!          if (ispin.eq.2) then
!             write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for beta electrons=  ', -shiftdown,'  +/-',dabs(DE)
!             write(ifu_log,*)
!             write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', QBeta
!             write(ifu_log,*)
!             write(ifu_log,'(A,F10.5)') ' Total number of electrons:  ', QAlpha+QBeta
!          end if
!          write(ifu_log,*)'--------------------------------------------------------'
!       end do
!    else
!       root_fail = .true.
!       E0=shift-Z 
!       E1=shift
!       E2=shift+Z
!       if (root_fail) then
!          print*,'MULLER method'
!          call MULLER(QTot,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
!          if(k .eq. Max .or. E2<EMin .or. E2>EMax) then
!             print *, 'Warning: MULLER method failed to find root. Using SECANT.'
!             root_fail = .true.
!          else
!             shift = E3
!             root_fail = .false.
!          end if
!       end if
!       if (root_fail) then
!          print*,'SECANT method'
!          call SECANT(QTot,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
!          if(k .eq. Max .or. E3<EMin .or. E3>EMax) then
!             print *, 'Warning: SECANT method failed to find root. Using BISEC.'
!             root_fail = .true.
!          else
!             shift = E3
!             root_fail = .false.
!          end if
!       end if
!       if (root_fail) then
!          print *, 'BISEC method'
!          shift = BISEC(QTot,EMin,EMax,Delta,5*Max,K)
!          DE=Delta
!          if(k.lt.5*Max) root_fail = .false.
!          if(k.ge.5*Max) print *, 'Warning: BISECT method failed to find root. Skipping this cycle.'
!       end if
!
!       write(ifu_log,*)'-----------------------------------------------'
!       write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy= ',-shift,'  +/-',dabs(DE)
!       write(ifu_log,*)
!       write(ifu_log,'(A,F10.5)') ' Number of alpha electrons: ', QAlpha
!       write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', QBeta
!       write(ifu_log,'(A,F10.5)') ' Total number of electrons:  ', QAlpha+QBeta
!       write(ifu_log,*)'-----------------------------------------------'
!    end if
!    ADDP = .not. root_fail
!    return
!  end subroutine CompDensMat

!**************************************************************
!* Subroutine for determining Fermi energy and density matrix *
!* for some total charge                                      *
!**************************************************************
!* Pre-condition:                                             *
!*   HC: Fock-matrix                                          *
!*   shiftup,shiftdown: starting values for Fermi-energies    *
!* Results:                                                   *
!*   PC: Density-matrix                                       *
!*   shiftup,shiftdown: Fermi-energies                        *
!**************************************************************
  subroutine CompDensMat_SOC(ADDP)
    use numeric, only: SECANT, MULLER, BISEC
#ifdef G03ROOT
    use g03Common, only: GetNAE, GetNBE
#endif
#ifdef G09ROOT
    use g09Common, only: GetNAE, GetNBE
#endif
    use parameters, only: FermiAcc,ChargeAcc,Max,QExcess
    !use ieee_arithmetic
    implicit none

    logical,intent(out) :: ADDP
    integer :: i,j, k,cond
    real*8 :: E0,E1,E2,E3,DE,Z, Delta, Epsilon
    
    logical :: root_fail
    
    Z=0.1d0               
    Delta=FermiAcc
    Epsilon=ChargeAcc*(NCDEl+QExcess)

    root_fail = .true.
    E0=shift-Z 
    E1=shift
    E2=shift+Z  
    if (root_fail) then
        print*,'SECANT method'
        call SECANT(QTot_SOC,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
        if(k .eq. Max .or. E3<EMin .or. E3>EMax) then
           print *, 'Warning: SECANT method failed to find root. Using MULLER.'
           root_fail = .true.
        else
           shift = E3
           root_fail = .false.
        end if
    end if      
    if (root_fail) then
        print*,'MULLER method'
        call MULLER(QTot_SOC,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
        if(k .eq. Max .or. E2<EMin .or. E2>EMax) then
           print *, 'Warning: MULLER method failed to find root. Using BISEC.'
           root_fail = .true.
        else
           shift = E3
           root_fail = .false.
       end if
    end if  
    if (root_fail) then
       print *, 'BISEC method'
       shift = BISEC(QTot_SOC,EMin,EMax,Delta,5*Max,K)
       DE=Delta
       if(k.lt.5*Max) root_fail = .false.
       if(k.ge.5*Max) print *, 'Warning: BISECT method failed to find root. Skipping this cycle.'
    end if

    write(ifu_log,*)'-----------------------------------------------'
    write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy= ',-shift,'  +/-',dabs(DE)
    write(ifu_log,*)
    write(ifu_log,'(A,F10.5)') ' Total number of electrons:  ', Q_SOC
    write(ifu_log,*)'-----------------------------------------------'
    ADDP = .not. root_fail

    return
  end subroutine CompDensMat_SOC


  !**************************************************************
  !* Subroutine for determining Fermi energy and density matrix *
  !* for some total charge                                      *
  !**************************************************************
  !* Pre-condition:                                             *
  !*   HC: Fock-matrix                                          *
  !*   shiftup,shiftdown: starting values for Fermi-energies    *
  !* Results:                                                   *
  !*   PC: Density-matrix                                       *
  !*   shiftup,shiftdown: Fermi-energies                        *
  !**************************************************************
  subroutine CompDensMat2(ADDP)
    use numeric, only: SECANT, MULLER, BISEC
#ifdef G03ROOT
    use g03Common, only: GetNAE, GetNBE
#endif
#ifdef G09ROOT
    use g09Common, only: GetNAE, GetNBE
#endif
    use parameters, only: FermiAcc,ChargeAcc,Max,QExcess
    !use ieee_arithmetic
    implicit none

    logical,intent(out) :: ADDP
    integer :: i,j, k,cond
    real*8 :: E0,E1,E2,E3,DE,Z, Delta, Epsilon
    
    logical :: root_fail
    
    Z=0.1d0           
    Delta=FermiAcc
    Epsilon=ChargeAcc*(NCDEl+QExcess)

    root_fail = .true.
    if( NSpin == 2 .and. SPINLOCK )then
       print'(A,f10.1)',' SPINLOCK: NAlpha-NBeta = ',dble((GetNAE()-GetNBE())*(NCDEl+QExcess))/dble(GetNAE()+GetNBE())
       do ispin=1,NSpin
          if (ispin.eq.1) E0=shiftup
          if (ispin.eq.2) E0=shiftdown
          E1=E0-Z
          E2=E0+Z 
          if( root_fail ) then
             print*,'Secant method'
             call SECANT(CompSpinPD,E0,E1,Delta,Epsilon,Max,E2,DE,Cond,K)
             if(k.eq.Max .or. E2<EMin .or. E2>EMax)then
                print *, 'Warning: Secant method failed to find root. Using BISEC.'
                root_fail = .true.
             else
                if (ispin.eq.1)shiftup=E2
                if (ispin.eq.2)shiftdown=E2
                root_fail = .false.
             end if
          end if
          if( root_fail ) then
             print*,'Muller method'
             call MULLER(CompSpinPD,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
             if(k.eq.Max .or. E3<EMin .or. E3>EMax)then
                print *, 'Warning: Muller method failed to find root. Using BISEC.'
                root_fail = .true.
             else
                if (ispin.eq.1)shiftup=E3
                if (ispin.eq.2)shiftdown=E3
                root_fail = .false.
             end if
          end if
          if( root_fail )then
             print *, 'BISEC method'
             if( ispin.eq.1) shiftup = BISEC(CompSpinPD,EMin,EMax,Delta,Max,K)
             if( ispin.eq.2) shiftdown = BISEC(CompSpinPD,EMin,EMax,Delta,Max,K)
             DE=Delta
             if(k.lt.Max) root_fail = .false.
          end if
          write(ifu_log,*)'--------------------------------------------------------'
          if (ispin.eq.1) then
             write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for alpha electrons= ', -shiftup,'  +/-',dabs(DE)
             write(ifu_log,*)
             write(ifu_log,'(A,F10.5)') ' Number of alpha electrons: ', QAlpha
          end if
          if (ispin.eq.2) then
             write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for beta electrons=  ', -shiftdown,'  +/-',dabs(DE)
             write(ifu_log,*)
             write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', QBeta
          end if
          write(ifu_log,*)'--------------------------------------------------------'
       end do
    else
       E0=shift
       E1=E0-Z 
       E2=E0+Z 
       if( root_fail )then
          print*,'Secant method'
          call SECANT(CompPD,E0,E1,Delta,Epsilon,Max,E2,DE,Cond,K)
          if(k.eq.Max .or. E2<EMin .or. E2>EMax)then
             print *, 'Warning: Secant method failed to find root. Using BISEC.'
             root_fail = .true.
          else
             shift = E2
             root_fail = .false.
          end if
       end if
       if( root_fail )then
          print*,'Muller method'
          call MULLER(CompPD,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
          if(k.eq.Max .or. E3<EMin .or. E3>EMax)then
             print *, 'Warning: Muller method failed to find root. Using BISEC.'
             root_fail = .true.
          else
             shift = E3
             root_fail = .false.
          end if
       end if
       if( root_fail )then
          print *, 'BISEC method'
          print *, EMin, EMax
          shift = BISEC(CompPD,EMin,EMax,Delta,Max,K)
          DE=Delta
          if(k.lt.Max) root_fail = .false.
       end if
       write(ifu_log,*)'-----------------------------------------------'
       write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy= ',-shift,'  +/-',dabs(DE)
       write(ifu_log,*)
       write(ifu_log,'(A,F10.5)') ' Number of alpha electrons: ', QAlpha
       write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', QBeta
       write(ifu_log,*)'-----------------------------------------------'
    end if
    ADDP = .not. root_fail
  end subroutine CompDensMat2

!-----------------------------------------------------------------------------
    subroutine CompDensMat2_SOC(ADDP)
!**************************************************************
!* Subroutine for determining Fermi energy and density matrix *
!* for some total charge with SOC along the imaginary axis    *
!**************************************************************
!* Pre-condition:                                             *
!*   HC: Fock-matrix                                          *
!*   shift: starting values for Fermi-energies                *
!* Results:                                                   *
!*   PC: Density-matrix                                       *
!*   shift: Fermi-energies                                    *
!**************************************************************

    use numeric, only: SECANT, MULLER, BISEC
    use parameters, only: FermiAcc,ChargeAcc,Max,QExcess
    !use ieee_arithmetic
    implicit none

    logical,intent(out) :: ADDP
    integer :: i,j, k,cond
    real*8 :: E0,E1,E2,E3,DE,Z, Delta, Epsilon
    
    logical :: root_fail
    
    Z=0.1d0              
    Delta=FermiAcc
    Epsilon=ChargeAcc*(NCDEl+QExcess)

    root_fail = .true.
    E0=shift
    E1=E0-Z 
    E2=E0+Z 
    if( root_fail )then
       print*,'Secant method'
       call SECANT(CompPD_SOC,E0,E1,Delta,Epsilon,Max,E2,DE,Cond,K)
       if(k.eq.Max .or. E2<EMin .or. E2>EMax)then
          print *, 'Warning: Secant method failed to find root. Using BISEC.'
          root_fail = .true.
       else
          shift = E2
          root_fail = .false.
       end if
    end if
    if( root_fail )then
       print*,'Muller method'
       call MULLER(CompPD_SOC,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
       if(k.eq.Max .or. E3<EMin .or. E3>EMax)then
          print *, 'Warning: Muller method failed to find root. Using BISEC.'
          root_fail = .true.
       else
          shift = E3
          root_fail = .false.
       end if
    end if
    if( root_fail )then
       print *, 'BISEC method'
       print *, EMin, EMax
       shift = BISEC(CompPD_SOC,EMin,EMax,Delta,Max,K)
       DE=Delta
       if(k.lt.Max) root_fail = .false.
    end if
    write(ifu_log,*)'-----------------------------------------------'
    write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy= ',-shift,'  +/-',dabs(DE)
    write(ifu_log,*)
    write(ifu_log,'(A,F10.5)') ' Number of electrons:  ', Q_SOC
    write(ifu_log,*)'-----------------------------------------------'
    ADDP = .not. root_fail

  end subroutine CompDensMat2_SOC
  
subroutine CompLocalMu(ADDP)
    !use numeric, only: LocalSecant, LocalMuller, LocalBISEC
    use numeric, only: Secant, Muller, BISEC
    use g09Common, only: GetNAE, GetNBE
    use cluster, only: LoAOrbNo, HiAOrbNo, LoAOrbNo, HiAOrbNo
    use parameters, only: FermiAcc,ChargeAcc,Max,QExcess,BiasVoltage,SPIN_Beg,SPIN_End
    !use ieee_arithmetic
    implicit none

    logical,intent(out) :: ADDP
    integer :: i,j, k,cond
    real :: E0,E1,E2,E3,DE,Z, Delta, Epsilon
    
    logical :: root_fail

!-------------------------------------------------------
    integer :: n,l
    real :: sdeg, ro_a, ro_b, chargeregion, spinregion
    real, dimension(NAOrbs,NAOrbs) :: rho_a, rho_b, tmp
!-------------------------------------------------------
    
    write(ifu_log,*)'-------------------------------------'
    write(ifu_log,*)'------- I am in CompLocalMu ---------'
    write(ifu_log,*)'-------------------------------------'
    write(ifu_log,*)'-------------------------------------'
    write(ifu_log,*)'---  Mulliken population analysis ---'
    write(ifu_log,*)'-------------------------------------'


    if (NSpin.eq.2) sdeg=1.0d0
    if (NSpin.eq.1) sdeg=2.0d0

    rho_a = matmul( PD(1,:,:), SD )
    if( NSpin == 2 ) rho_b = matmul( PD(2,:,:), SD )

    write(ifu_log,*)'---------------------------------------------------------------------'
    write(ifu_log,*)'Charges in selected region to compute local electrochemical potential'
    write(ifu_log,*)'---------------------------------------------------------------------'

    chargeregion=0.0d0
    spinregion=0.0d0
    do j=SPIN_Beg,SPIN_End
       ro_a=0.0d0
       ro_b=0.0d0
       do i=LoAOrbNo(j),HiAOrbNo(j)
          ro_a=ro_a+rho_a(i,i)
          !if(NSpin==1) ro_b=ro_a
          if(NSpin==2) ro_b=ro_b+rho_b(i,i)
       end do
       alphaelec(j-SPIN_Beg+1)=ro_a
       if(NSpin==1) betaelec(j-SPIN_Beg+1)=ro_a
       if(NSpin==2) betaelec(j-SPIN_Beg+1)=ro_b
       if(NSpin ==1 ) write(ifu_log,'(A,I5,A,f9.5)')'Atom:',j,' El.dens:',ro_a*sdeg
       if(NSpin ==2 ) write(ifu_log,'(A,I5,A,f9.5,A,f9.5)')'Atom:',j,' El.dens:', (ro_a+ro_b),' Sp.dens:', (ro_a-ro_b)

       chargeregion=chargeregion+(ro_a+ro_b)
       if (NSpin == 2) spinregion=spinregion+(ro_a-ro_b)
    end do
    
   write(ifu_log,*)'----------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in selected region:',chargeregion*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in selected region:',spinregion
   write(ifu_log,*)'----------------------------------------------------'

!-------------------------------------------------------------------------------------------------
    Z=10.0d0*FermiAcc
    Delta=FermiAcc
    !Epsilon=ChargeAcc*(NCDEl+QExcess)
    Epsilon=ChargeAcc*(chargeregion)

!-------------------------------------------------------------------------------------------------

    !--- ADDED BY CARLOS -----------------------------------------------
    print*,'BiasVoltage = ', BiasVoltage
    !-------------------------------------------------------------------
!if( NSpin == 2 .and. SPINLOCK )then
if( NSpin == 2 )then
!  print*,'SPINLOCK is true'
  print*,'NSpin = 2'
  print'(A,f10.1)',' SPINLOCK: NAlpha-NBeta = ',dble((GetNAE()-GetNBE())*(NCDEl+QExcess))/dble(GetNAE()+GetNBE())
  do spinatom=SPIN_Beg,SPIN_End
    alphalocalshift(spinatom) = shiftup
    betalocalshift(spinatom) = shiftdown
    do ispin=1,NSpin
       root_fail = .true.
       !if (ispin.eq.1) E1=shiftup
       !if (ispin.eq.2) E1=shiftdown
       if (ispin.eq.1) E1=alphalocalshift(spinatom)
       if (ispin.eq.2) E1=betalocalshift(spinatom)
       E0=E1-Z
       E2=E1+Z
       if( root_fail )then
          print*,'MULLER method'
          call MULLER(FPart,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
          if(k.eq.Max .or. E3<EMin .or. E3>EMax) then
             print *, 'Warning: MULLER method failed to find root. Using SECANT.'
             root_fail = .true.
          else
             !if (ispin.eq.1)shiftup=E3
             !if (ispin.eq.2)shiftdown=E3
             if (ispin.eq.1)alphalocalshift(spinatom)=E3
             if (ispin.eq.2)betalocalshift(spinatom)=E3
             root_fail = .false.
          end if
       end if
       if( root_fail )then
          print*,'SECANT method'
          call SECANT(FPart,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
          if(k.eq.Max .or. E3<EMin .or. E3>EMax) then
             print *, 'Warning: SECANT method failed to find root. Using BISEC.'
             root_fail = .true.
          else
             !if (ispin.eq.1)shiftup=E3
             !if (ispin.eq.2)shiftdown=E3
             if (ispin.eq.1)alphalocalshift(spinatom)=E3
             if (ispin.eq.2)betalocalshift(spinatom)=E3
             root_fail = .false.
          end if
       end if
       if (root_fail) then
          print *, 'BISEC method'
          !if (ispin.eq.1) shiftup = LocalBISEC(FPart,EMin,EMax,Delta,5*Max,K)
          !if (ispin.eq.2) shiftdown = LocalBISEC(FPart,EMin,EMax,Delta,5*Max,K)
          if (ispin.eq.1) alphalocalshift(spinatom) = BISEC(FPart,EMin,EMax,Delta,5*Max,K)
          if (ispin.eq.2) betalocalshift(spinatom) = BISEC(FPart,EMin,EMax,Delta,5*Max,K)          
          DE=Delta
          if(k.lt.5*Max) root_fail = .false.
          if(k.ge.5*Max) print *, 'Warning: BISECT method failed to find root. Skipping this cycle.'
       end if
       write(ifu_log,*)'--------------------------------------------------------'
       if (ispin.eq.1) then
          write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for alpha electrons= ', -alphalocalshift(spinatom),'  +/-',dabs(DE)
          write(ifu_log,*)
          write(ifu_log,'(A,F10.5)') ' Number of alpha electrons: ', LocalQAlpha
       end if
       if (ispin.eq.2) then
          write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for beta electrons=  ', -betalocalshift(spinatom),'  +/-',dabs(DE)
          write(ifu_log,*)
          write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', LocalQBeta
          write(ifu_log,*)
          write(ifu_log,'(A,F10.5)') ' Total number of electrons:  ', LocalQAlpha+LocalQBeta
       end if
       write(ifu_log,*)'--------------------------------------------------------'
    end do
  end do
    write(ifu_log,*)'-----------------------------------------------'
    write(ifu_log,*)'-----------------------------------------------'
  do ispin=1,NSpin 
    do spinatom=SPIN_Beg,SPIN_End
      !write(ifu_log,'(A,F9.5,A,f9.5)') ' Local Mu= ',-alphalocalshift(spinatom),'  +/-',dabs(DE)
      if (ispin.eq.1)write(ifu_log,'(A,I5,A,F9.5)') 'Atom',spinatom,'	Alpha Local Mu:  ',-alphalocalshift(spinatom)
      if (ispin.eq.2)write(ifu_log,'(A,I5,A,F9.5)') 'Atom',spinatom,'	Beta Local Mu:  ',-betalocalshift(spinatom)
    end do
    write(ifu_log,*)'-----------------------------------------------'
  end do
else
  do spinatom=SPIN_Beg,SPIN_End
    !print*,'SPINLOCK is false'
    print*,'NSpin = 1'
    !alphalocalshift(spinatom) = shiftup
    !betalocalshift(spinatom) = shiftdown
    alphalocalshift(spinatom) = shift
    root_fail = .true.
    E0=alphalocalshift(spinatom)-Z 
    E1=alphalocalshift(spinatom)
    E2=alphalocalshift(spinatom)+Z
    if (root_fail) then
      print*,'MULLER method'
      !call MULLER(QXPart,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
      call MULLER(QXPart,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
      if(k .eq. Max .or. E2<EMin .or. E2>EMax) then
         print *, 'Warning: MULLER method failed to find root. Using SECANT.'
         root_fail = .true.
      else
         alphalocalshift(spinatom) = E3
         root_fail = .false.
      end if
    end if
    if (root_fail) then
      print*,'SECANT method'
      !call SECANT(QXPart,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
      call SECANT(QXPart,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
      if(k .eq. Max .or. E3<EMin .or. E3>EMax) then
         print *, 'Warning: SECANT method failed to find root. Using BISEC.'
         root_fail = .true.
      else
         alphalocalshift(spinatom) = E3
         root_fail = .false.
      end if
    end if
    if (root_fail) then
      print *, 'BISEC method'
      !Delta = 10*Delta !Added by Carlos Salgado
      !Max = 10*Max !Added by Carlos Salgado
      !alphalocalshift(spinatom) = BISEC(QXPart,EMin,EMax,Delta,5*Max,K)
      alphalocalshift(spinatom) = BISEC(QXPart,EMin,EMax,Delta,5*Max,K)
      !alphalocalshift(spinatom) = BISEC(QXPart,EMin,EMax,1.0D2*Delta,500*Max,K)
      DE=Delta
      if(k.lt.5*Max) root_fail = .false.
      if(k.ge.5*Max) print *, 'Warning: BISECT method failed to find root. Skipping this cycle.'
    end if

    write(ifu_log,*)'-----------------------------------------------'
    write(ifu_log,'(A,F9.5,A,f9.5)') ' Alpha Local Mu:  ',-alphalocalshift(spinatom),'  +/-',dabs(DE)
    write(ifu_log,*)
    write(ifu_log,'(A,F10.5)') ' Number of alpha electrons:  ', LocalQAlpha
    write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', LocalQBeta
    write(ifu_log,'(A,F10.5)') ' Total number of electrons:  ', LocalQAlpha+LocalQBeta
    write(ifu_log,*)'-----------------------------------------------'

  end do
    write(ifu_log,*)'-----------------------------------------------'
    write(ifu_log,*)'-----------------------------------------------'
  do spinatom=SPIN_Beg,SPIN_End
    !write(ifu_log,'(A,F9.5,A,f9.5)') ' Local Mu= ',-alphalocalshift(spinatom),'  +/-',dabs(DE)
    write(ifu_log,'(A,I5,A,F9.5)') 'Atom',spinatom,'	Local Mu:  ',-alphalocalshift(spinatom)
  end do
    write(ifu_log,*)'-----------------------------------------------'
    write(ifu_log,*)'-----------------------------------------------'

end if
  ADDP = .not. root_fail
end subroutine CompLocalMu  
  
  ! 
  ! Computes the density matrix for a fixed chemical potential mu
  ! by integrating Greens function on matsubara axis. No lower energy
  ! bound required anymore!
  ! - Returns number of electrons in device region
  !
  ! - replaces old code in functions F(x) and QXTot(x)
  !
  real*8 function CompPD( mu )
     use constants
     use util
     use numeric, only: CHDiag, gauleg, RTrace
     use parameters, only: eta, PAcc
!    USE IFLPORT
     use omp_lib
#ifdef G03ROOT
    use g03Common, only: GetNAE, GetNBE
#endif
#ifdef G09ROOT
    use g09Common, only: GetNAE, GetNBE
#endif
     implicit none

     ! chemical potential
     real*8, intent(in) :: mu

     complex*16, dimension(NAOrbs,NAOrbs) :: GD

     integer, parameter :: nmin=1, nmax=10, npmax=2**nmax-1
     real*8, dimension(2*npmax) :: x, w
     real*8 :: Ei, dEdx, Q, QQ, DPD, E0 !, Qi
     integer :: n, np, i, j, k, l !, info, ierr
     
     real*8, parameter :: x0 = 0.5d0
     real*8 :: aa, bb, cc, EM

     EM = EMax
     aa = EM/x0
     bb = aa*(1.0d0-x0)**2
     cc = EM - bb/(1-x0)

     shift = mu

     Q = d_zero

     do n=nmin,nmax

        np=2**n-1
        
        QQ = Q
        PD=d_zero
        QAlpha = d_zero; QBeta = d_zero
        
        ! Compute Gauss-Legendre abcsissas and weights
        call gauleg(0.0d0,2.0d0,x(1:2*np),w(1:2*np),2*np)

        do ispin=1,NSpin
!$OMP PARALLEL PRIVATE(Ei,dEdx,GD,DPD)
!$OMP DO
           do i=1,np
              Ei = 2.0d0*EMax*x(i)
              if( x(i) > 0.5d0 ) Ei = 0.5d0*EMax/(1.0d0-x(i))
              dEdx = 2.0d0*EMax
              if( x(i) > 0.5d0 ) dEdx = 0.5d0*EMax/(1.0d0-x(i))**2   
              call gplus0( ui*Ei, GD )
!$OMP CRITICAL
              do k=1,NAOrbs
                 do l=1,NAOrbs
                    DPD = w(i)*(dEdx*real(GD(k,l))/d_pi + 0.5d0*InvSD(k,l))
                    PD(ispin,k,l) = PD(ispin,k,l) + DPD
                    if(ispin.eq.1) QAlpha = QAlpha + DPD*SD(l,k)
                    if(ispin.eq.2) QBeta = QBeta + DPD*SD(l,k)
                 end do
              end do
!$OMP END CRITICAL
           end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
        end do
        if(NSpin.eq.1) QBeta=QAlpha
        Q = QAlpha + QBeta        
        if( n > nmin .and. abs(Q-QQ) < PAcc ) exit
        if (n == nmax) print*, 'Warning!, not enough integration points'
     end do
        
     print '(A,I4,A)', ' Integration of density matrix along imaginary axis has needed ', np, ' points.'
    !print '(A,F10.5,A,F10.5,A,F8.5)', ' mu =', -mu, '  Num. of electrons =', Q, ' +/-', abs(Q-QQ) 

     CompPD = Q - dble(NCDEl)

  end function CompPD

!------------------------------------------------------------------------------------
   real*8 function CompPD_SOC( mu )
! 
! Computes the density matrix for a fixed chemical potential mu
! by integrating the spectral function. 
!
! Returns number of electrons in device region
!
     use constants
     use util
     use numeric, only: CHDiag, gauleg, RTrace
     use parameters, only: eta, PAcc
!    USE IFLPORT
     use omp_lib
     implicit none

     ! chemical potential
     real*8, intent(in) :: mu

     complex*16, dimension(DNAOrbs,DNAOrbs) :: GDR, GDA

     integer, parameter :: nmin=1, nmax=15, npmax=2**nmax-1
     real*8, dimension(npmax) :: x, w
     real*8 :: Q, QQ, R
     integer :: n, np, i, j, k, l, M !, info, ierr
     
     ! Radius of complex contour integration
     ! add 10eV just in case 
     R  = 0.5*abs(EMin)+10.0

     shift = mu
     PD_SOC=c_zero

! Retarded density matrix
  
     Q=0.0d0
     do n=nmin,nmax

        np=2**n-1
        
        QQ = Q
        PD_SOC_R=c_zero       
                
        ! Compute Gauss-Legendre abcsissas and weights
        call gauleg(0.d0,d_pi,x(1:np),w(1:np),np)

!$OMP PARALLEL PRIVATE(GDR)
!$OMP DO
        do i=1,np
           call gplus0_SOC( R*exp(ui*x(i))-R, GDR, 1 )        
    !print*,'x(i)',i,x(i),GDR(1,1)

!$OMP CRITICAL
           do k=1,DNAOrbs
              do l=1,DNAOrbs
                 PD_SOC_R(k,l) = PD_SOC_R(k,l) + w(i)*(ui*R*exp(ui*x(i))*GDR(k,l))
              end do
           end do
!$OMP END CRITICAL
        end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
        Q=d_zero
        do k=1,DNAOrbs
           do l=1,DNAOrbs
              Q = Q + real(PD_SOC_R(k,l)*S_SOC(l,k))
           end do
        end do
       !print*,'Q=',Q
        if( n > nmin .and. abs(Q-QQ) < PAcc ) exit
        if (n == nmax) print*, 'Warning!, not enough integration points'

     end do
        
     print '(A,I4,A)', 'Integration of retarded density matrix along the upper semicircle on the complex plane has needed ', np, ' points.'

! Advanced density matrix

     Q=0.0d0
     do n=nmin,nmax

        np=2**n-1
        
        QQ = Q
        PD_SOC_A=c_zero
        
        call gauleg(0.0d0,-d_pi,x(1:np),w(1:np),np)

!$OMP PARALLEL PRIVATE(GDA)
!$OMP DO
        do i=1,np           
            call gplus0_SOC(R*exp(ui*(x(i)))-R, GDA, -1)
    !print*,'x(i)',i,x(i),GDA(1,1)
!$OMP CRITICAL
             do k=1,DNAOrbs
                do l=1,DNAOrbs
                   PD_SOC_A(k,l) = PD_SOC_A(k,l) + w(i)*(ui*R*exp(ui*(x(i)))*GDA(k,l))
                end do
             end do
!$OMP END CRITICAL
        end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
        Q=d_zero
        do k=1,DNAOrbs
           do l=1,DNAOrbs
              Q = Q + real(PD_SOC_A(k,l)*S_SOC(l,k))
           end do
        end do
       !print*,'Q=',Q
        if( n > nmin .and. abs(Q-QQ) < PAcc ) exit
        if (n == nmax) print*, 'Warning!, not enough integration points'
     end do
     print '(A,I4,A)', 'Integration of advanced density matrix along the lower semicircle on the complex plane has needed ', np, ' points.'
        
     PD_SOC =-ui*(PD_SOC_R  - PD_SOC_A)/(2.0d0*d_pi)

     Q_SOC= d_zero
     do k=1,DNAOrbs
        do l=1,DNAOrbs
           Q_SOC = Q_SOC + PD_SOC(k,l)*S_SOC(l,k)
        end do
     end do
     print '(A,F10.5,A,F10.5,A,F8.5)', ' mu =', -mu, '  Num. of electrons =', Q_SOC
       !print '(2(F20.5))',PD_SOC(1,1)
       !print '(4(F20.5))',PD_SOC(3,1),PD_SOC(1,3)          

     CompPD_SOC = Q_SOC - dble(NCDEl)

  end function CompPD_SOC

!-------------------------------------------------------------------------------------
  
  ! 
  ! Computes the spin-resolved density matrix for a fixed chemical potential mu
  ! by integrating Greens function on matsubara axis. No lower energy
  ! bound required anymore!
  ! - Returns number of electrons in device region
  !
  ! - replaces old code in functions F(x) and QXTot(x)
  !
  real*8 function CompSpinPD( mu )
     use constants
     use util
     use numeric, only: CHDiag, gauleg, RTrace
     use parameters, only: eta, PAcc
!    USE IFLPORT
     use omp_lib
#ifdef G03ROOT
    use g03Common, only: GetNAE, GetNBE
#endif
#ifdef G09ROOT
    use g09Common, only: GetNAE, GetNBE
#endif
     implicit none

     ! chemical potential
     real*8, intent(in) :: mu

     complex*16, dimension(NAOrbs,NAOrbs) :: GD

     integer, parameter :: nmin=1, nmax=8, npmax=2**nmax-1
     real*8, dimension(2*npmax) :: x, w
     real*8 :: Ei, dEdx, Q, QQ, DPD !, E0, Qi
     integer :: n, np, i, j, k, l !, info, ierr
     real*8 :: NCDAB

     ! Number of alpha OR beta electrons in central region C of device
     if(ispin.eq.1) NCDAB = dble(GetNAE()*NCDEl)/dble(GetNAE()+GetNBE())
     if(ispin.eq.2) NCDAB = dble(GetNBE()*NCDEl)/dble(GetNAE()+GetNBE())
     
     shift = mu

     Q=0.0d0
     do n=nmin,nmax
        np=2**n-1
        
        QQ = Q
        Q = 0.0d0
        PD(ispin,:,:) = 0.0d0
        
        ! Compute Gauss-Legendre abcsissas and weights
        call gauleg(0.0d0,2.0d0,x(1:2*np),w(1:2*np),2*np)

!$OMP PARALLEL PRIVATE(Ei,dEdx,GD,DPD)
!$OMP DO
        do i=1,np
           Ei = 2.0d0*EMax*x(i)
           if( x(i) > 0.5d0 ) Ei = 0.5d0*EMax/(1.0d0-x(i))
           dEdx = 2.0d0*EMax
           if( x(i) > 0.5d0 ) dEdx = 0.5d0*EMax/(1.0d0-x(i))**2   
           call GPlus0( ui*Ei, GD )
!$OMP CRITICAL
           do k=1,NAOrbs
              do l=1,NAOrbs
                 DPD = w(i)*(dEdx*real(GD(k,l))/d_pi + 0.5d0*InvSD(k,l))
                 PD(ispin,k,l) = PD(ispin,k,l) + DPD
                 Q = Q + DPD*SD(l,k)
              end do
           end do
!$OMP END CRITICAL
        end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
        if( n > nmin .and. abs(Q-QQ) < PAcc ) exit
     end do
        
     print '(A,I4,A)', ' Integration of density matrix has needed ', np, ' points.'
     print '(A,F10.5,A,F10.5,A,F8.5)', ' mu =', -mu, '  Num. of electrons =', Q, ' +/-', abs(Q-QQ) 

     if(ispin.eq.1) QAlpha = Q
     if(ispin.eq.2) QBeta = Q

     CompSpinPD = Q - NCDAB

   end function CompSpinPD
   
  ! 
  ! Computes the density matrix for a fixed chemical potential mu
  ! by integrating Greens function on matsubara axis. No lower energy
  ! bound required anymore!
  ! - Returns number of electrons in device region
  !
  ! - replaces old code in functions F(x) and QXTot(x)
  !
  real function CompEnergy( mu )
     use constants
     use util
     use numeric, only: CHDiag, gauleg, RTrace
     use parameters, only: eta, PAcc
    use g09Common, only: GetNAE, GetNBE
     implicit none

     ! chemical potential
     real, intent(in) :: mu

     complex*16, dimension(NAOrbs,NAOrbs) :: GD

     integer, parameter :: nmin=1, nmax=8, npmax=2**nmax-1
     !integer, parameter :: nmin=1, nmax=16, npmax=2**nmax-1
     real, dimension(2*npmax) :: x, w
     real :: Ei, dEdx, Q, QQ, DPD, E0 !, Qi
     integer :: n, np, i, j, k, l !, info, ierr
     
     real, parameter :: x0 = 0.5d0
     real :: aa, bb, cc, EM

     real :: IntDOSE, IntDOSEAlpha, IntDOSEBeta, IntDOSEE, CompPD

     EM = EMax
     aa = EM/x0
     bb = aa*(1.0d0-x0)**2
     cc = EM - bb/(1-x0)
     
     print *, " - SHIFT:", -shift
     shift = mu ! ADDED BECAUSE IT IS IN OTHER FUNCTIONS
     print *, " - mu:", -mu
     !shift = mu ! I HAVE TO CLEAR THIS TO AVOID OVERWRITING SHIFT
     print *, " - SHIFT:", -shift
     Q = d_zero
     IntDOSE = d_zero

     do n=nmin,nmax

        np=2**n-1
        
        QQ = Q
        PD=d_zero
        QAlpha = d_zero; QBeta = d_zero
        IntDOSEE = IntDOSE
        IntDOSEAlpha = d_zero; IntDOSEBeta = d_zero
        
        ! Compute Gauss-Legendre abcsissas and weights
        call gauleg(0.0d0,2.0d0,x(1:2*np),w(1:2*np),2*np)

        do ispin=1,NSpin
!$OMBLABLABLAP PARALLEL PRIVATE(Ei,dEdx,GD,DPD)
!$OMBLABLABLAP DO
           do i=1,np
              Ei = 2.0d0*EMax*x(i)
              if( x(i) > 0.5d0 ) Ei = 0.5d0*EMax/(1.0d0-x(i))
              dEdx = 2.0d0*EMax
              if( x(i) > 0.5d0 ) dEdx = 0.5d0*EMax/(1.0d0-x(i))**2   
              call GPlus0( ui*Ei, GD )
!$OMBLABLABLAP CRITICAL
              do k=1,NAOrbs
                 do l=1,NAOrbs
                    DPD = w(i)*(dEdx*real(GD(k,l))/d_pi + 0.5d0*InvSD(k,l))
                    PD(ispin,k,l) = PD(ispin,k,l) + DPD
                    if(ispin.eq.1) QAlpha = QAlpha + DPD*SD(l,k)
                    if(ispin.eq.2) QBeta = QBeta + DPD*SD(l,k)
                    if(ispin.eq.1) IntDOSEAlpha = IntDOSEAlpha + (x(i)-mu)*DPD*SD(l,k)
                    if(ispin.eq.2) IntDOSEBeta = IntDOSEBeta + (x(i)-mu)*DPD*SD(l,k)
                 end do
              end do
!$OMBLABLABLAP END CRITICAL
           end do
           print *, ' Grid points:',np
           !print *, ' Charge alpha ', QAlpha
           !print *, ' Charge beta ', QBeta
           print *, ' Charge alpha + beta ', QAlpha + QBeta
           !print *, ' IntDOSE alpha ', IntDOSEAlpha
           !print *, ' IntDOSE beta ', IntDOSEBeta
           print *, ' IntDOSE alpha + beta ', IntDOSEAlpha + IntDOSEBeta
!$OMBLABLABLAP END DO
!$OMBLABLABLAP BARRIER
!$OMBLABLABLAP END PARALLEL
        end do
        if(NSpin.eq.1) QBeta=QAlpha
        Q = QAlpha + QBeta 
        if(NSpin.eq.1) IntDOSEBeta=IntDOSEAlpha  
        IntDOSE = IntDOSEAlpha + IntDOSEBeta     
        if( n > nmin .and. abs(Q-QQ) < PAcc*NCDEl ) exit
        !if( n > nmin .and. abs(Q-QQ) < (1/16)*PAcc*NCDEl ) exit
     end do
        
     print '(A,I4,A)', ' Integration of P has needed ', np, ' points.'
     print '(A,F10.5,A,F10.5,A,F8.5)', ' mu =', -mu, '  Num. of electrons =', Q, ' +/-', abs(Q-QQ) 

     print *, ' Computed Charge alpha + beta ', Q

     CompPD = Q - dble(NCDEl)
     CompEnergy = IntDOSE

  end function CompEnergy   

  
  !*************************************
  !* Compute retarded Green's function *
  !*************************************
  subroutine gplus0(z,green)
    use PARAMETERS, only: eta,glue
    use constants, only: c_zero, ui
#ifdef PGI
    use lapack_blas, only: zgetri, zgetrf
#endif 
    implicit none
    external zgetrf, zgetri    

    integer :: i, j, info, omp_get_thread_num
    integer, dimension(NAOrbs) :: ipiv
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl,sigr, temp 
    complex*16 :: work(4*NAOrbs) 

    complex*16, intent(in) :: z 

    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: green
    

    ! Initilization 
    green=c_zero
    sigr=-ui*eta*SD 
    sigl=-ui*eta*SD 

    call CompSelfEnergies( ispin, z, sigl, sigr, 1 )
    sigr=glue*sigr
    sigl=glue*sigl
    
   !print*,'selfenergy right'
   !print*,sigr(1,1),sigr(1,2),sigr(2,1)
   !print*,'selfenergy left'
   !print*,sigl(1,1),sigl(1,2),sigl(2,1)
   !print*,'..........'
    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************
    do i=1,NAOrbs
       do j=1,NAOrbs
          green(i,j)=(z-shift)*SD(i,j)-HD(ispin,i,j)-sigl(i,j)-sigr(i,j)
       enddo
    enddo

    call zgetrf(NAOrbs,NAOrbs,green,NAOrbs,ipiv,info)
    call zgetri(NAOrbs,green,NAOrbs,ipiv,work,4*NAOrbs,info)

  end subroutine gplus0


  !********************************************************
  !* Compute retarded Green's function and gamma matrices *
  !********************************************************
  subroutine gplus(z,green,gammar,gammal)
    use PARAMETERS, only: eta, glue
    use constants, only: c_zero, ui
#ifdef PGI
    use lapack_blas, only: zgetri, zgetrf
#endif

    implicit none
    external zgetrf, zgetri 
    
    integer :: i, j, info, omp_get_thread_num
    integer, dimension(NAOrbs) :: ipiv
    complex*16 :: work(4*NAOrbs)
    complex*16, intent(in) :: z 
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl,sigr, temp 
    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: green
    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: gammar,gammal
    
    ! Initilization 
    green=c_zero
    gammar=c_zero
    gammal=c_zero
    sigr=-ui*eta*SD 
    sigl=-ui*eta*SD 

    call CompSelfEnergies( ispin, z, sigl, sigr, 1 )
    sigr=glue*sigr
    sigl=glue*sigl
    
    !************************************************************************
    !c Coupling matrices                                   
    !************************************************************************

    gammar=ui*(sigr-conjg(transpose(sigr)))
    gammal=ui*(sigl-conjg(transpose(sigl)))

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************
    do i=1,NAOrbs
       do j=1,NAOrbs
          green(i,j)=(z-shift)*SD(i,j)-HD(ispin,i,j)-sigl(i,j)-sigr(i,j)
       enddo
    enddo

    call zgetrf(NAOrbs,NAOrbs,green,NAOrbs,ipiv,info)
    call zgetri(NAOrbs,green,NAOrbs,ipiv,work,4*NAOrbs,info)

    !write(ifu_log,*)omp_get_thread_num(),green(1,1)
  end subroutine gplus

  !*************************************
  !* Compute lesser Green's function *
  !*************************************
  subroutine glesser(z,gless)
    use parameters, only: eta, biasvoltage, glue
    use constants, only: c_zero, ui, c_one
#ifdef PGI
    use lapack_blas, only: zgetri, zgetrf
#endif

    implicit none
    external zgetrf, zgetri, zgemm     

    integer :: i, j, info, error
    integer, dimension(NAOrbs) :: ipiv
    complex*16 :: work(4*NAOrbs)

    complex*16, intent(in) :: z

    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: gless
    complex*16, dimension(:,:), allocatable :: green,gammar,gammal,sigl,sigr,glessR,glessL

     allocate(glessR(NAOrbs,NAOrbs),glessL(NAOrbs,NAOrbs),green(NAOrbs,NAOrbs),stat=error)
     if (error /= 0 ) then 
        print*,"Problems allocating" 
        stop
     end if
     allocate(sigl(NAOrbs,NAOrbs),sigr(NAOrbs,NAOrbs),gammar(NAOrbs,NAOrbs),gammal(NAOrbs,NAOrbs),stat=error)
     if (error /= 0 ) then 
        print*,"Problems allocating" 
        stop
     end if

    gammar=c_zero
    gammal=c_zero
    sigr=-ui*eta*SD 
    sigl=-ui*eta*SD 

    call CompSelfEnergies( ispin, z, sigl, sigr, 1 )
    sigr=glue*sigr
    sigl=glue*sigl
    
    !************************************************************************
    !c Coupling matrices                                   
    !************************************************************************
    gammar=ui*(sigr-conjg(transpose(sigr)))
    gammal=ui*(sigl-conjg(transpose(sigl)))

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************
    do i=1,NAOrbs
       do j=1,NAOrbs
          green(i,j)=(z-shift)*SD(i,j)-HD(ispin,i,j)-sigl(i,j)-sigr(i,j)
       enddo
    enddo

    call zgetrf(NAOrbs,NAOrbs,green,NAOrbs,ipiv,info)
    call zgetri(NAOrbs,green,NAOrbs,ipiv,work,4*NAOrbs,info)

    !************************************************************************
    !* G< (Gless)
    !************************************************************************

     call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one,gammal,NAOrbs,green,NAOrbs, &
          &           c_zero,gless,NAOrbs)
     call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,green,NAOrbs,gless,NAOrbs, &
          &           c_zero,glessL,NAOrbs)
     call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one,gammar,NAOrbs,green,NAOrbs, &
          &           c_zero,gless,NAOrbs)
     call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,green,NAOrbs,gless,NAOrbs, &
          &           c_zero,glessR,NAOrbs)

     gless=ui*(glessL*theta(dble(z)-biasvoltage/2)+glessR*theta(dble(z)+biasvoltage/2))

     deallocate(glessR,glessL,sigl,sigr,gammar,gammal,stat=error)
     if (error /= 0 ) then 
        print*,"Problems deallocating" 
        stop
     end if

  end subroutine glesser

  !*********************************************************************
  !* Compute lesser Green's function  with SOC (doubling the matrices) *
  !*********************************************************************
  subroutine glesser_SOC(z,gless)
    use parameters, only: eta, glue, biasvoltage
    use constants, only: c_zero, ui, c_one
#ifdef PGI
    use lapack_blas, only: zgetri, zgetrf, zgemm
#endif

    implicit none
    external zgetrf, zgetri, zgemm     

    integer :: i, j, info, omp_get_thread_num, error
    integer, dimension(DNAOrbs) :: ipiv
    complex*16, dimension(4*DNAOrbs) :: work
    complex*16, intent(in) :: z 
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl1,sigr1
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl2,sigr2
    complex*16, dimension(DNAOrbs,DNAOrbs), intent(out) :: gless
    complex*16, dimension(DNAOrbs,DNAOrbs) :: gammar,gammal
    complex*16, dimension(DNAOrbs,DNAOrbs) :: sigmar,sigmal
    complex*16, dimension(DNAOrbs,DNAOrbs) :: green,glessL,glessR

!   complex*16, dimension(:,:), allocatable :: green,gammar,gammal,sigl,sigr,glessR,glessL

!    allocate(glessR(NAOrbs,NAOrbs),glessL(NAOrbs,NAOrbs),green(NAOrbs,NAOrbs),stat=error)
!    if (error /= 0 ) then 
!       print*,"Problems allocating" 
!       stop
!    end if
!    allocate(sigl(NAOrbs,NAOrbs),sigr(NAOrbs,NAOrbs),gammar(NAOrbs,NAOrbs),gammal(NAOrbs,NAOrbs),stat=error)
!    if (error /= 0 ) then 
!       print*,"Problems allocating" 
!       stop
!    end if

   
    ! Initilization 
    green=c_zero
    gammar=c_zero
    gammal=c_zero
    sigmar=c_zero
    sigmal=c_zero


    sigr1=-ui*eta*SD 
    sigl1=-ui*eta*SD 
    sigr2=-ui*eta*SD 
    sigl2=-ui*eta*SD 

    call CompSelfEnergies( 1, z, sigl1, sigr1, 1 )
   
    if (NSpin == 2) call CompSelfEnergies( 2, z, sigl2, sigr2, 1 )
    sigr1=glue*sigr1
    sigl1=glue*sigl1
    sigr2=glue*sigr2
    sigl2=glue*sigl2
    
    !************************************************************************
    !c Coupling matrices                                   
    !************************************************************************
    if (NSpin == 2) then
       do i=1,NAOrbs
       do j=1,NAOrbs
          sigmar(i,j)=sigr1(i,j)
          sigmar(i+NAOrbs,j+NAOrbs)=sigr2(i,j)
          sigmal(i,j)=sigl1(i,j)
          sigmal(i+NAOrbs,j+NAOrbs)=sigl2(i,j)
       end do
       end do
    else
       do i=1,NAOrbs
       do j=1,NAOrbs
          sigmar(i,j)=sigr1(i,j)
          sigmar(i+NAOrbs,j+NAOrbs)=sigr1(i,j)
          sigmal(i,j)=sigl1(i,j)
          sigmal(i+NAOrbs,j+NAOrbs)=sigl1(i,j)
       end do
       end do
    end if

    gammar=ui*(sigmar-conjg(transpose(sigmar)))
    gammal=ui*(sigmal-conjg(transpose(sigmal)))

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************

    do i=1,DNAOrbs
       do j=1,DNAOrbs
          green(i,j)=(z-shift)*S_SOC(i,j)-H_SOC(i,j)-sigmal(i,j)-sigmar(i,j)
       enddo
    enddo

    call zgetrf(DNAOrbs,DNAOrbs,green,DNAOrbs,ipiv,info)
    call zgetri(DNAOrbs,green,DNAOrbs,ipiv,work,4*DNAOrbs,info)


    !************************************************************************
    !* G< (Gless)
    !************************************************************************

     call zgemm('N','C',DNAOrbs,DNAOrbs,DNAOrbs,c_one,gammal,DNAOrbs,green,DNAOrbs, &
          &           c_zero,gless,DNAOrbs)
     call zgemm('N','N',DNAOrbs,DNAOrbs,DNAOrbs,c_one,green,DNAOrbs,gless,DNAOrbs, &
          &           c_zero,glessL,DNAOrbs)
     call zgemm('N','C',DNAOrbs,DNAOrbs,DNAOrbs,c_one,gammar,DNAOrbs,green,DNAOrbs, &
          &           c_zero,gless,DNAOrbs)
     call zgemm('N','N',DNAOrbs,DNAOrbs,DNAOrbs,c_one,green,DNAOrbs,gless,DNAOrbs, &
          &           c_zero,glessR,DNAOrbs)

     gless=ui*(glessL*theta(dble(z)-biasvoltage/2)+glessR*theta(dble(z)+biasvoltage/2))

!    deallocate(glessR,glessL,sigl,sigr,gammar,gammal,stat=error)
!    if (error /= 0 ) then 
!       print*,"Problems deallocating" 
!       stop
!    end if

  end subroutine glesser_SOC



  ! *************
  ! Step function
  ! *************
  function theta(x)
    implicit none
    real*8,intent(in) :: x
    real*8 :: theta
    theta=0.0d0
    if(x.lt.0.0d0) theta=1.0
    return
  end function theta

  !****************************************
  ! Function gives excess charge of device
  ! in dependence of Fermi energy x
  ! when spin is locked
  !****************************************
  double precision function QTot_SPINLOCK(x)
    use parameters, only: biasvoltage,QExcess
    use constants, only: d_pi, d_zero, c_zero
#ifdef G03ROOT
    use g03Common, only: GetNAE, GetNBE
#endif
#ifdef G09ROOT
    use g09Common, only: GetNAE, GetNBE
#endif
    implicit none

    real*8, intent(in) :: x

    real*8 :: chargeup, chargedown, rrr, a, b, Q
    integer :: i,j,M

    shift=x
    
    PD(ispin,:,:) = d_zero
    PDOUT(ispin,:,:) = c_zero

    ! Radius of complex contour integration
    ! add 10eV just in case 
    rrr = 0.5*abs(EMin)+10.0d0;

    !c c Integral limits ... (a,b)
    a = 0.d0
    b = d_pi
    M=1000
    call IntCompPlane(rrr,a,b,M,d_zero-dabs(biasvoltage/2.0))
      
    ! Density matrix out of equilibirum
    if (biasvoltage /= 0.0) then
       M=1000
       call IntRealAxis(-dabs(biasvoltage/2.0),dabs(biasvoltage/2.0),M)
       PD(ispin,:,:) = PD(ispin,:,:) + REAL(PDOUT(ispin,:,:))
    end if
    ! Density matrix out of equilibirum
      
    Q=d_zero
    !do i=1,NAOrbs
    do i=NCDAO1, NCDAO2
       do j=1,NAOrbs
          Q=Q+PD(ispin,i,j)*SD(j,i)
       end do
    end do
    if (ispin.eq.1) QAlpha = Q
    if (ispin.eq.2) QBeta  = Q
 !if (ispin.eq.1) f = QAlpha - dble(GetNAE())
 !if (ispin.eq.2) f = QBeta  - dble(GetNBE())
    if (ispin.eq.1) QTot_SPINLOCK = QAlpha - dble(GetNAE()*NCDEl)/dble(GetNAE()+GetNBE()) -QExcess/2.0
    if (ispin.eq.2) QTot_SPINLOCK = QBeta  - dble(GetNBE()*NCDEl)/dble(GetNAE()+GetNBE()) -QExcess/2.0
    return
  end function QTot_SPINLOCK


  !****************************************
  ! Function gives excess charge of device
  ! in dependence of Fermi energy x
  ! when spin is not locked
  !****************************************
  double precision function QTot(x)
    use parameters, only: biasvoltage,QExcess
    use constants, only: d_pi, d_zero, c_zero
#ifdef G03ROOT
    use g03Common, only: GetNAE, GetNBE
#endif
#ifdef G09ROOT
    use g09Common, only: GetNAE, GetNBE
#endif
    implicit none

    real*8, intent(in) :: x

    real*8 :: rrr, a, b, Q
    integer :: i,j,M,omp_get_thread_num


    do ispin=1,NSpin

       PD(ispin,:,:) = d_zero
       PDOUT(ispin,:,:) = c_zero

       shift=x

       Q = d_zero

       ! Radius of complex contour integration
       ! add 10eV just in case 
       rrr = 0.5*abs(EMin)+10.0d0;

       !c c Integral limits ... (a,b)
       a = 0.d0
       b = d_pi
       M=1000

       call IntCompPlane(rrr,a,b,M,d_zero-dabs(biasvoltage/2.0))

       ! Density matrix out of equilibirum
       if (biasvoltage /= 0.0) then
          M=1000
          call IntRealAxis(-dabs(biasvoltage/2.0),dabs(biasvoltage/2.0),M)
          PD(ispin,:,:) = PD(ispin,:,:) + REAL(PDOUT(ispin,:,:))
       end if
       ! Density matrix out of equilibirum

       !do i=1,NAOrbs
       do i=NCDAO1, NCDAO2
          do j=1,NAOrbs
             Q=Q+PD(ispin,i,j)*SD(j,i)
          end do
       end do
       if( ispin == 1 ) QAlpha = Q
       if( ispin == 2 ) QBeta  = Q
    end do

    if( NSpin == 1 ) QBeta = QAlpha
    !QTot = QAlpha + QBeta -dble(GetNAE()) -dble(GetNBE())
    QTot = QAlpha + QBeta - dble(NCDEl) - QExcess
    return
  end function QTot

  !****************************************
  ! Function gives excess charge of device 
  ! dependent of Fermi energy x (shift)
  ! integrating in the complex plane
  ! when SOC is present
  !****************************************
  double precision function QTot_SOC(x)
    use parameters, only: biasvoltage,QExcess
    use constants, only: d_pi, d_zero, c_zero
#ifdef G03ROOT
    use g03Common, only: GetNAE, GetNBE
#endif
#ifdef G09ROOT
    use g09Common, only: GetNAE, GetNBE
#endif
    implicit none

    real*8, intent(in) :: x

    real*8 :: rrr, a, b, Q, E0
    integer :: i,j,M,omp_get_thread_num

       shift=x
    !write(ifu_log,*)omp_get_thread_num(),'in qxtot',shift
       Q = 0.d0
       PD_SOC = c_zero       
       PDOUT_SOC = c_zero

       ! Radius of complex contour integration
       ! add 10eV just in case 
       rrr = 0.5*abs(EMin)+10.0;

       !c c Integral limits ... (a,b)
       M=1000
       a = 0.d0
       b = d_pi
       call IntCompPlane_SOC(1,rrr,a,b,M,-dabs(biasvoltage/2.0))!Retarded
       PD_SOC_R=PD_SOC
	   
       M=1000
       a = -d_pi
       b = 0.d0
       call IntCompPlane_SOC(-1,rrr,a,b,M,-dabs(biasvoltage/2.0))!Advanced
       PD_SOC_A=PD_SOC
	   
       PD_SOC = PD_SOC_R - PD_SOC_A

      !if (NSpin == 2)

       ! Density matrix out of equilibirum
       if (biasvoltage /= 0.0) then
          M=1000
          call IntRealAxis_SOC(-dabs(biasvoltage/2.0),dabs(biasvoltage/2.0),M)
          PD_SOC = PD_SOC + PDOUT_SOC
       end if
       ! Density matrix out of equilibirum

      do i=NCDAO1, NCDAO2
          do j=1,NAOrbs
             Q=Q+PD_SOC(i,j)*S_SOC(j,i)
          end do
       end do
       do i=NCDAO1+NAOrbs, NCDAO2+NAOrbs
          do j=NAOrbs+1,DNAOrbs
             Q=Q+PD_SOC(i,j)*S_SOC(j,i)
          end do
       end do

    Q_SOC = Q
    !print*,Q_SOC
    QTot_SOC = Q - dble(NCDEl) - QExcess
    return
  end function QTot_SOC


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !c                                                                              c
  !c     Change of variable no. 3 for the numerical integration:                  c
  !c                                                                              c
  !c     int_{-infty}^{Eq} DOS(E)dE = int_{-1}^{1} DOS(E(x))*(dE/dx)dx            c
  !c                                                                              c
  !c     E = Em*(1-bx)/(1+x)                                                      c
  !c                                                                              c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  real*8 function edex3(Em,b,x)
    implicit none
    real*8 :: Em, b ,x
    edex3 = 0.5d0*((Em-b)*x + (Em+b))
    return
  end function edex3
  
  !******************************!
  ! Writing the Hamiltonian      !
  !******************************!
  subroutine Hamiltonian
    use cluster, only: NALead, NAMol, NAOAtom, NAOMol, NEmbedBL
    USE parameters, only: Hamilton
#ifdef G03ROOT
    use g03Common, only: GetAtmCo
#endif
#ifdef G09ROOT
    use g09Common, only: GetAtmCo
#endif
    use constants, only: Bohr
    implicit none

    integer :: i,j, I1, is ,n
    real*8 :: ro_a, ro_b

    write(ifu_log,*)'-------------------------------------'
    write(ifu_log,*)'---  Self-consistent Hamiltonian  ---'
    write(ifu_log,*)'-------------------------------------'

    
    !open(ifu_ham,file='V.'//trim(xxx)//'.dat',status='unknown')
    write(ifu_log,*)'---------'
    write(ifu_log,*)'Left lead'
    write(ifu_log,*)'---------'
    I1=0
    do j=1,NALead(1)
       ro_a=0.0d0
       ro_b=0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+hd(1,i,i)
          if(NSpin==2) ro_b=ro_b+hd(2,i,i)
       end do
       ro_a=ro_a/NAOAtom(j)
       ro_b=ro_b/NAOAtom(j)
       if(NSpin ==1 ) write(ifu_log,1011)'Atom:',j, ro_a
       if(NSpin ==2 ) write(ifu_log,1011)'Atom:',j, (ro_a+ro_b)/2.0
       IF (Hamilton) THEN
       if(NSpin ==1 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a
       if(NSpin ==2 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b)/2.0
       END IF
       I1 = I1 + NAOAtom(j)
    end do
    
    if (NAOMol().ge.1) then 
       write(ifu_log,*)'--------'
       write(ifu_log,*)'Molecule'
       write(ifu_log,*)'--------'
       do j = NALead(1)+1,NALead(1)+NAMol()
          ro_a = 0.0d0
          ro_b = 0.0d0
          do i=I1+1,I1+NAOAtom(j)
             ro_a=ro_a+hd(1,i,i)
             if(NSpin==2) ro_b=ro_b+hd(2,i,i)
          end do
          ro_a=ro_a/NAOAtom(j)
          ro_b=ro_b/NAOAtom(j)
          if(NSpin==1) write(ifu_log,1011)'Atom:',j, ro_a
          if(NSpin==2) write(ifu_log,1011)'Atom:',j, (ro_a+ro_b)/2.0
          IF (Hamilton) THEN
          if(NSpin ==1 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a
          if(NSpin ==2 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b)/2.0
          END IF
          I1 = I1 + NAOAtom(j)
       end do
    end if
    
    write(ifu_log,*)'----------'
    write(ifu_log,*)'Right lead'
    write(ifu_log,*)'----------'
    do j=NALead(1)+NAMol()+1,NALead(1)+NAMol()+NALead(2)
       ro_a=0.0d0
       ro_b = 0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+hd(1,i,i)
          if(NSpin==2) ro_b=ro_b+hd(2,i,i)
       end do
       ro_a=ro_a/NAOAtom(j)
       ro_b=ro_b/NAOAtom(j)
       I1 = I1 + NAOAtom(j)
       if(NSpin==1) write(ifu_log,1011)'Atom:',j, ro_a
       if(NSpin==2) write(ifu_log,1011)'Atom:',j, (ro_a+ro_b)/2.0
       IF (Hamilton) THEN
       if(NSpin ==1 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a
       if(NSpin ==2 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b)/2.0
       END IF
    end do
1012 format(4f10.5)
1011 format(a6,i4,f12.6)
     !close(ifu_ham)
  end subroutine Hamiltonian 

  !******************************!
  ! Mulliken population analysis !
  !******************************!
  subroutine MullPop
    use cluster, only: NALead, NAMol, NAOAtom, NAOMol, LoAOrbNo, HiAOrbNo
    USE parameters, only: Mulliken, LDOS_Beg, LDOS_End, PrtHatom
#ifdef G03ROOT
    use g03Common, only: GetAtmCo, GetNAtoms
#endif
#ifdef G09ROOT
    use g09Common, only: GetAtmCo, GetNAtoms
#endif
    use constants, only: Bohr
    implicit none

    integer :: i,j, I1, is ,n, l, iAtom, jAtom
    real*8 :: sdeg, ro_a, ro_b, chargemol, chargelead1, chargelead2, spinlead1, spinlead2, spinmol
    real*8, dimension(NAOrbs,NAOrbs) :: rho_a, rho_b
 
    write(ifu_log,*)'-------------------------------------'
    write(ifu_log,*)'---  Mulliken population analysis ---'
    write(ifu_log,*)'-------------------------------------'
    
    if (PrtHatom > 1) then
       do iAtom=1,GetNAtoms()
          do jAtom=1,GetNAtoms()
             !if( SpinEdit(iAtom) == -1 .and. SpinEdit(jAtom) == -1 )then
             if( iAtom == PrtHatom .and. jAtom == PrtHatom )then
                PRINT *, " Density matrix atom ",iAtom," is: "
                PRINT *, " Up-Up "  
                do i=LoAOrbNo(PrtHatom),HiAOrbNo(PrtHatom)                                                               
                   PRINT '(1000(F11.5))', ( (PD(1,i,j)), j=LoAOrbNo(PrtHatom),HiAOrbNo(PrtHatom) )               
                end do                                                                                                       
                PRINT *, " Down-Down "                                                                              
                do i=LoAOrbNo(PrtHatom),HiAOrbNo(PrtHatom)                                                            
                   PRINT '(1000(F11.5))', ( (PD(2,i,j)), j=LoAOrbNo(PrtHatom),HiAOrbNo(PrtHatom) )             
                end do                                                                              
             end if   
          end do
       end do
    end if          
    
    if (NSpin.eq.2) sdeg=1.0d0
    if (NSpin.eq.1) sdeg=2.0d0

    !rho_a = matmul( PD(1,:,:), SD )
    rho_a=0.0d0
    rho_b=0.0d0
    do i=1,NAOrbs
    do j=1,NAOrbs
    do l=1,NAOrbs
       rho_a(i,j)=rho_a(i,j)+PD(1,i,l)*SD(l,j)
    end do
    end do
    end do
    !if( NSpin == 2 ) rho_b = matmul( PD(2,:,:), SD )
    if( NSpin == 2 ) then             
    do i=1,NAOrbs
    do j=1,NAOrbs
    do l=1,NAOrbs
       rho_b(i,j)=rho_b(i,j)+PD(2,i,l)*SD(l,j)
    end do
    end do
    end do
    end if

    write(ifu_log,*)'----------------------'
    write(ifu_log,*)'Charges in electrode 1'
    write(ifu_log,*)'----------------------'
    I1=0
    chargelead1=0.0d0
    spinlead1=0.0d0
    do j=1,NALead(1)
       ro_a=0.0d0
       ro_b=0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+rho_a(i,i)
         !print*,rho_a(i,i)
          if(NSpin==2) ro_b=ro_b+rho_b(i,i)
       end do
       if(NSpin ==1 ) write(ifu_log,1011)'Atom:',j,' El.dens:',ro_a*sdeg
       if(NSpin ==2 ) write(ifu_log,1012)'Atom:',j,' El.dens:', (ro_a+ro_b),' Sp.dens:', (ro_a-ro_b)
       IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg,AtomDOSEF(1,j)
         if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
       END IF
       IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg
         if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b)
       END IF
       chargelead1=chargelead1+ro_a+ro_b
       if (NSpin == 2) spinlead1=spinlead1+(ro_a-ro_b)
       I1 = I1 + NAOAtom(j)
    end do
    
   write(ifu_log,*)'---------------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in electrode 1:',chargelead1*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in electrode 1:',spinlead1
   write(ifu_log,*)'---------------------------------------------------------'

    chargemol=0.0d0
    spinmol=0.0d0
    if (NAOMol().ge.1) then 
       write(ifu_log,*)'-------------------'
       write(ifu_log,*)'Charges in molecule'
       write(ifu_log,*)'-------------------'
       do j = NALead(1)+1,NALead(1)+NAMol()
          ro_a = 0.0d0
          ro_b = 0.0d0
          do i=I1+1,I1+NAOAtom(j)
             ro_a=ro_a+rho_a(i,i)
             if(NSpin==2) ro_b=ro_b+rho_b(i,i)
          end do
          if(NSpin==1) write(ifu_log,1011)'Atom:',j,' El.dens:',sdeg*ro_a
          if(NSpin==2) write(ifu_log,1012)'Atom:',j,' El.dens:', (ro_a+ro_b),' Sp.dens:', (ro_a-ro_b)
          IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
            if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg,AtomDOSEF(1,j)
            if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
          END IF
          IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
            if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg
            if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b)
          END IF
          chargemol=chargemol+ro_a+ro_b
          if (NSpin == 2) spinmol=spinmol+(ro_a-ro_b)
          I1 = I1 + NAOAtom(j)
       end do
    end if
    
   write(ifu_log,*)'----------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in molecule:',chargemol*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in molecule:',spinmol
   write(ifu_log,*)'----------------------------------------------------'

    write(ifu_log,*)'----------------------'
    write(ifu_log,*)'Charges in electrode 2'
    write(ifu_log,*)'----------------------'
    chargelead2=0.0d0
    spinlead2=0.0d0
    do j=NALead(1)+NAMol()+1,NALead(1)+NAMol()+NALead(2)
       ro_a = 0.0d0
       ro_b = 0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+rho_a(i,i)
          if(NSpin==2) ro_b=ro_b+rho_b(i,i)
       end do
       I1 = I1 + NAOAtom(j)
       if(NSpin==1) write(ifu_log,1011)'Atom:',j,' El.dens:',ro_a*sdeg
       if(NSpin==2) write(ifu_log,1012)'Atom:',j,' El.dens:', (ro_a+ro_b),' Sp.dens:', (ro_a-ro_b)
       IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg,AtomDOSEF(1,j)
         if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
       END IF
       IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg
         if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b)
       END IF
       chargelead2=chargelead2+ro_a+ro_b
       if (NSpin == 2) spinlead2=spinlead2+(ro_a-ro_b)
    end do

   write(ifu_log,*)'---------------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in electrode 2:',chargelead2*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in electrode 2:',spinlead2
   write(ifu_log,*)'---------------------------------------------------------'

1011 format(a6,i4,a10,f8.4)
1012 format(a6,i4,a10,f8.4,a10,2f9.4)
1013 format(7f9.4)
  end subroutine MullPop
  
  !*******************************************************************!
  ! Mulliken population analysis with non-collinear or SOC hamiltonian!
  !*******************************************************************!
  subroutine MullPop_SOC
    use cluster, only: NALead, NAMol, NAOAtom, NAOMol, LoAOrbNo, HiAOrbNo
    USE parameters, only: Mulliken, LDOS_Beg, LDOS_End, biasvoltage, PrtHatom
#ifdef G03ROOT
    use g03Common, only: GetAtmCo, GetNAtoms
#endif
#ifdef G09ROOT
    use g09Common, only: GetAtmCo, GetNAtoms
#endif
    use constants, only: Bohr, d_zero, c_zero
    implicit none

    integer :: i,j, I1, is ,n, l, iAtom, jAtom
    real*8 :: sdeg, ro_a, ro_ab, ro_ba, ro_ab_I, ro_ba_I, ro_b, spindens, chargemol, chargelead1, chargelead2, spinlead1, spinlead2, spinmol
    real*8, dimension(NAOrbs,NAOrbs) :: rho_a, rho_ab, rho_ba, rho_ab_I, rho_ba_I, rho_b, tmp
    complex*16, dimension(DNAOrbs,DNAOrbs) :: rho

    write(ifu_log,*)'-----------------------------------------------'
    write(ifu_log,*)'--- Mulliken population analysis (with SOC) ---'
    write(ifu_log,*)'-----------------------------------------------'

    rho = c_zero
    rho_a = d_zero
    rho_ab = d_zero
    rho_ab_I = d_zero
    rho_ba = d_zero
    rho_ba_I = d_zero
    rho_b = d_zero     
    
    rho =matmul(PD_SOC,S_SOC)
       
    do i=1,NAOrbs
    do j=1,NAOrbs
       rho_a(i,j)=rho_a(i,j)+rho(i,j)
       rho_b(i,j)=rho_b(i,j)+rho(i+NAOrbs,j+NAOrbs)
       rho_ab(i,j)=rho_ab(i,j)+REAL(rho(i,j+NAOrbs))
       rho_ba(i,j)=rho_ba(i,j)+REAL(rho(i+NAOrbs,j))
       rho_ab_I(i,j)=rho_ab_I(i,j)+IMAG(rho(i,j+NAOrbs))
       rho_ba_I(i,j)=rho_ba_I(i,j)+IMAG(rho(i+NAOrbs,j))
    end do
    end do             

    write(ifu_log,*)'----------------------'
    write(ifu_log,*)'Charges in electrode 1'
    write(ifu_log,*)'----------------------'
    I1=0
    chargelead1=0.0
    spinlead1=0.0
    do j=1,NALead(1)
       ro_a=0.0d0
       ro_ab=0.0d0
       ro_ab_I=0.0d0
       ro_ba=0.0d0
       ro_ba_I=0.0d0
       ro_b=0.0d0
       spindens=0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+rho_a(i,i)
         !print*,rho_a(i,i)
          ro_b=ro_b+rho_b(i,i)
       !  IF(NSpin==2 .or. (NSpin == 1 .and. biasvoltage /= 0.0)) THEN 
            ro_ab=ro_ab+rho_ab(i,i)
            ro_ab_I=ro_ab_I+rho_ab_I(i,i)
            ro_ba=ro_ba+rho_ba(i,i)  
            ro_ba_I=ro_ba_I+rho_ba_I(i,i)
       !  END IF 
       end do
      !if(NSpin == 1 .and. biasvoltage == 0.0) write(ifu_log,2011)'Atom:',j,' El.dens:',ro_a+ro_b
      !IF(NSpin == 2 .or. (NSpin == 1 .and. biasvoltage /= 0.0)) THEN     
         spindens = sqrt((ro_ab+ro_ba)**2+(ro_ab_I-ro_ba_I)**2+(ro_a-ro_b)**2)
         write(ifu_log,2012)'Atom:',j,' El.dens:',(ro_a+ro_b),' Sp.dens.x:',(ro_ab+ro_ba),' Sp.dens.y:',(ro_ab_I-ro_ba_I),' Sp.dens.z:',(ro_a-ro_b),' Coll. sp.dens:',spindens
      !END IF
       IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
      !  if(NSpin ==1 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a+ro_b,AtomDOSEF(1,j)
      !  if(NSpin ==2 .or. (NSpin == 1 .and. biasvoltage /= 0.0) ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
         write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
       END IF
       IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
       ! if(NSpin ==1 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a+ro_b
       ! if(NSpin ==2 .or. (NSpin == 1 .and. biasvoltage /= 0.0) ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b)
         write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b)
       END IF
       chargelead1=chargelead1+(ro_a+ro_b)
       spinlead1=spinlead1+spindens
       I1 = I1 + NAOAtom(j)
    end do
    
   write(ifu_log,*)'---------------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in electrode 1:',chargelead1
   if (NSpin == 2) write(ifu_log,*)'Total spin in electrode 1:',spinlead1
   write(ifu_log,*)'---------------------------------------------------------'

    chargemol=0.0d0
    spinmol=0.0d0
    if (NAOMol().ge.1) then 
       write(ifu_log,*)'-------------------'
       write(ifu_log,*)'Charges in molecule'
       write(ifu_log,*)'-------------------'
       do j = NALead(1)+1,NALead(1)+NAMol()
       ro_a=0.0d0
       ro_ab=0.0d0
       ro_ab_I=0.0d0
       ro_ba=0.0d0
       ro_ba_I=0.0d0
       ro_b=0.0d0
       spindens=0.0d0
          do i=I1+1,I1+NAOAtom(j)
             ro_a=ro_a+rho_a(i,i)
             ro_b=ro_b+rho_b(i,i)
         !   IF(NSpin==2 .or. (NSpin == 1 .and. biasvoltage /= 0.0)) THEN 
                ro_ab=ro_ab+rho_ab(i,i)
                ro_ab_I=ro_ab_I+rho_ab_I(i,i)
                ro_ba=ro_ba+rho_ba(i,i)  
                ro_ba_I=ro_ba_I+rho_ba_I(i,i)
         !   END IF 
          end do
      !if(NSpin == 1 .and. biasvoltage == 0.0) write(ifu_log,2011)'Atom:',j,' El.dens:',ro_a+ro_b
      !IF(NSpin == 2 .or. (NSpin == 1 .and. biasvoltage /= 0.0)) THEN     
            spindens = sqrt((ro_ab+ro_ba)**2+(ro_ab_I-ro_ba_I)**2+(ro_a-ro_b)**2)
            write(ifu_log,2012)'Atom:',j,' El.dens:',(ro_a+ro_b),' Sp.dens.x:',(ro_ab+ro_ba),' Sp.dens.y:',(ro_ab_I-ro_ba_I),' Sp.dens.z:',(ro_a-ro_b),' Coll. sp.dens:',spindens
      !   END IF
          IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
      !     if(NSpin ==1 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a+ro_b,AtomDOSEF(1,j)
      !     if(NSpin ==2 .or. (NSpin == 1 .and. biasvoltage /= 0.0) ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
            write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
          END IF
          IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
         !  if(NSpin ==1 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a+ro_b
         !  if(NSpin ==2 .or. (NSpin == 1 .and. biasvoltage /= 0.0) ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b)
            write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b)
          END IF
          chargemol=chargemol+ro_a+ro_b
          spinmol=spinmol+spindens
          I1 = I1 + NAOAtom(j)
       end do
    end if
    
   write(ifu_log,*)'----------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in molecule:',chargemol
   if (NSpin == 2) write(ifu_log,*)'Total spin in molecule:',spinmol
   write(ifu_log,*)'----------------------------------------------------'

    write(ifu_log,*)'----------------------'
    write(ifu_log,*)'Charges in electrode 2'
    write(ifu_log,*)'----------------------'
    chargelead2=0.0
    spinlead2=0.0
    do j=NALead(1)+NAMol()+1,NALead(1)+NAMol()+NALead(2)
       ro_a=0.0d0
       ro_ab=0.0d0
       ro_ab_I=0.0d0
       ro_ba=0.0d0
       ro_ba_I=0.0d0
       ro_b=0.0d0
       spindens=0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+rho_a(i,i)
          ro_b=ro_b+rho_b(i,i)
         !IF(NSpin==2 .or. (NSpin == 1 .and. biasvoltage /= 0.0)) THEN 
            ro_ab=ro_ab+rho_ab(i,i)
            ro_ab_I=ro_ab_I+rho_ab_I(i,i)
            ro_ba=ro_ba+rho_ba(i,i)  
            ro_ba_I=ro_ba_I+rho_ba_I(i,i)
         !END IF 
       end do
       I1 = I1 + NAOAtom(j)
    !  if(NSpin == 1 .and. biasvoltage == 0.0) write(ifu_log,2011)'Atom:',j,' El.dens:',ro_a+ro_b
    !  IF(NSpin == 2 .or. (NSpin == 1 .and. biasvoltage /= 0.0)) THEN     
         spindens = sqrt((ro_ab+ro_ba)**2+(ro_ab_I-ro_ba_I)**2+(ro_a-ro_b)**2)
         write(ifu_log,2012)'Atom:',j,' El.dens:',(ro_a+ro_b),' Sp.dens.x:',(ro_ab+ro_ba),' Sp.dens.y:',(ro_ab_I-ro_ba_I),' Sp.dens.z:',(ro_a-ro_b),' Coll. sp.dens:',spindens
    !  END IF
       IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
    !    if(NSpin ==1 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a+ro_b,AtomDOSEF(1,j)
    !    if(NSpin ==2 .or. (NSpin == 1 .and. biasvoltage /= 0.0) ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
         write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
       END IF
       IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
      !  if(NSpin ==1 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a+ro_b
     !   if(NSpin ==2 .or. (NSpin == 1 .and. biasvoltage /= 0.0) ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b)
         write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b)
       END IF
       chargelead2=chargelead2+ro_a+ro_b
       spinlead2=spinlead2+spindens
    end do

   write(ifu_log,*)'---------------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in electrode 2:',chargelead2
   if (NSpin == 2) write(ifu_log,*)'Total spin in electrode 2:',spinlead2
   write(ifu_log,*)'---------------------------------------------------------'    

2011 format(a6,i4,a10,f8.4)
2012 format(a6,i4,a10,f8.4,3(a12,f9.4),a16,f9.4)
2013 format(10f9.4)
  end subroutine MullPop_SOC  

!*******************************************************************************
!* Subroutine to compute LDOS(E) when it is not computed by Transmission        *
!*******************************************************************************
  SUBROUTINE LDOS
    use Cluster, only : hiaorbno, loaorbno
    use constants, only: c_one, c_zero, d_zero, d_pi
    use parameters, only: LDOS_Beg, LDOS_End, EW1, EW2, EStep, DOSEnergy
    use preproc, only: MaxAtm
!   USE IFLPORT
    use omp_lib
#ifdef G03ROOT
    use g03Common, only: GetNAtoms
#endif
#ifdef G09ROOT
    use g09Common, only: GetNAtoms
#endif

    real*8 :: energy, trans, DOS, energ
    real*8, dimension(10001) :: xxx
    real*8, dimension(MaxAtm) :: AtomDOS
    complex*16 :: cenergy,ctrans
    integer :: n, nsteps, i, imin, imax, info ,j
    
    complex*16, dimension(NAOrbs,NAOrbs) :: GammaL, GammaR, Green, T, temp, SG

    print *
    print *, "-------------------------"
    print *, "--- Calculating  LDOS ---"
    print *, "-------------------------"
    print *

    allocate ( AtomDOSEF(2,MaxAtm) )

    nsteps = (EW2-EW1)/EStep + 1
    do ispin=1,NSpin

       open(333,file='tempDOS',status='unknown')
#ifdef PGI
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,cenergy,energy,green,gammar,gammal) 
!$OMP DO SCHEDULE(STATIC,10)
#endif
       do n=1,nsteps
          energy=EW1+EStep*(n-1)
          cenergy=dcmplx(energy)

          !*********************************************************************
          !* Evaluation of the retarded "Green" function and coupling matrices *
          !*********************************************************************
          call gplus(cenergy,Green,GammaR,GammaL)

          ! Mulliken DOS 
#ifdef PGI
!$OMP CRITICAL
#endif
          SG = matmul( SD, green )
          ! computing total DOS
          DOS=d_zero
          AtomDOS=d_zero
          do j=1,GetNAtoms()
          do i=LoAOrbNo(j),HiAOrbNo(j)
             AtomDOS(j)=AtomDOS(j)-dimag(SG(i,i))/d_pi
             DOS=DOS-dimag(SG(i,i))/d_pi
          end do
          end do

          if (dabs(energy-DOSEnergy) < EStep/2) AtomDOSEF(ispin,:)=AtomDOS

          ! print out DOS and atomic orbital resolved DOS ***
          imin = LoAOrbNo(LDOS_Beg)
          if( imin < 1 ) imin = 1
          imax = HiAOrbNo(LDOS_End)
          if( imax > NAOrbs ) imax = NAOrbs
          call flush(333)
          write(333,3333) energy,DOS*(-1)**(ispin+1),(AtomDOS(j)*(-1)**(ispin+1),j=LDOS_Beg,LDOS_End),(-dimag(SG(i,i))*(-1)**(ispin+1)/d_pi,i=imin,imax)

#ifdef PGI
!$OMP END CRITICAL
#endif

       end do ! End of energy loop
      
#ifdef PGI
!$OMP END DO
!$OMP END PARALLEL
#endif

  ! Reordering in energy for nice output
       do n=1,nsteps
          energy=EW1+EStep*(n-1)
          rewind(333)
          do i=1,10000000000
          read(333,*)energ
          if (dabs(energy-energ) < 0.000001) then
             backspace(333)
             read(333,3333) (xxx(j),j=1,2+(LDOS_End-LDOS_Beg+1)+(imax-imin+1))
             write(ifu_dos,3333) (xxx(j),j=1,2+(LDOS_End-LDOS_Beg+1)+(imax-imin+1))
             exit
          end if
          end do
       end do
 
      close(333,status='delete')

        write(ifu_dos,*) '    '                    
      end do ! End of spin loop


3333 format(f10.5,10000E14.5)

  END SUBROUTINE LDOS

!*******************************
!* Subroutine to evaluate T(E) *
!*******************************
  subroutine transmission
    use Cluster, only : hiaorbno, loaorbno
    use constants, only: c_one, c_zero, d_zero, d_pi
    use parameters, only: NChannels,HTransm,EW1,EW2,EStep,LDOS_Beg,LDOS_End, DOSEnergy, SOC, ROT, FermiAcc, QExcess, ChargeAcc, DMIMAG
    use numeric, only: CMatPow, CHDiag, CDiag, sort, MULLER, RMatPow
    use preproc, only: MaxAtm
#ifdef PGI
    use lapack_blas, only: zgemm
#endif
!   USE IFLPORT
    use omp_lib
#ifdef G03ROOT
    use g03Common, only: GetNAtoms
#endif
#ifdef G09ROOT
    use g09Common, only: GetNAtoms
#endif
    implicit none
    external zgemm    

    real*8 :: energy, trans, DOS, energ, E0, E1, E2, E3, Delta, Epsilon, DE
    real*8, dimension(MaxAtm) :: AtomDOS
    real*8, dimension(10001) :: xxx
    complex*16 :: cenergy,ctrans
    integer :: n, nsteps, i, imin, imax, info, j, AllocErr, cond, k, wcount
    integer :: Max = 20
    complex*16, dimension(:,:), allocatable :: GammaL, GammaR, Green, T, temp, SG
    complex*16, dimension(:,:), allocatable :: Green_UU, Green_DD, Green_UD, Green_DU, GammaL_UU, GammaR_UU, GammaL_DD, GammaR_DD
    complex*16, dimension(:,:), allocatable :: DGammaL, DGammaR,  DGreen, DT, Dtemp, DSG
    complex*16, dimension(:,:),allocatable :: dummy
    complex*16, dimension(:), allocatable :: ctn,Dctn
    real*8, dimension(:), allocatable :: tn,Dtn
    real*8, dimension(:),allocatable   :: tchan1,tchan2
    real*8   :: T_uu,T_dd,T_ud,T_du,polar,trans2
    logical :: ADDP

    print *
    print *, "--------------------------------"
    print *, "--- Calculating Transmission ---"
    print *, "---  (and DOS if required)   ---"
    print *, "--------------------------------"
    print *

    if (SOC .or. ROT) then
       write(ifu_log,*)' Finding new Fermi level after adding SOC or rotating spins or both ..........'

       allocate(H_SOC(DNAOrbs,DNAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(PD_SOC(DNAOrbs,DNAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(PD_SOC_R(DNAOrbs,DNAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(PD_SOC_A(DNAOrbs,DNAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(PDOUT_SOC(DNAOrbs,DNAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(S_SOC(DNAOrbs,DNAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(InvS_SOC(DNAOrbs,DNAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(DGreen(DNAOrbs,DNAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(DGammaR(DNAOrbs,DNAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(DGammaL(DNAOrbs,DNAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(DT(DNAOrbs,DNAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(T(NAOrbs,NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(Dtemp(DNAOrbs,DNAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(temp(NAOrbs,NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(DSG(DNAOrbs,DNAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(Dtn(DNAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(Dctn(DNAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop

       allocate(GammaL_UU(NAOrbs,NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(GammaL_DD(NAOrbs,NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(GammaR_UU(NAOrbs,NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(GammaR_DD(NAOrbs,NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(Green_UU(NAOrbs,NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(Green_UD(NAOrbs,NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(Green_DU(NAOrbs,NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(Green_DD(NAOrbs,NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop

       call spin_orbit
       
    else

       allocate(GammaL(NAOrbs,NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(GammaR(NAOrbs,NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(Green(NAOrbs,NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(T(NAOrbs,NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(temp(NAOrbs,NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(SG(NAOrbs,NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(tn(NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop
       allocate(ctn(NAOrbs), STAT=AllocErr);if( AllocErr /= 0 ) stop

    end if

    if( NChannels > 0 ) then
       print *, "Number of eigen-channels to print out: ", NChannels 
       allocate( tchan1(NChannels), tchan2(NChannels), dummy(NChannels,NChannels), STAT = AllocErr )
       if( AllocErr /= 0 ) stop
    end if

    allocate( AtomDOSEF(2,MaxAtm), STAT=AllocErr)
    if( AllocErr /= 0 ) stop

    if (SOC .or. ROT) then
       if (DMIMAG) then
          call RMatPow( S_SOC, -1.0d0, InvS_SOC )
          call CompDensMat2_SOC(ADDP)
       else   
          call CompDensMat_SOC(ADDP)
       end if 
    else
       if (DMIMAG) then
          call CompDensMat2(ADDP)
       else
          call CompDensMat(ADDP)
       end if 
    end if 

    nsteps = (EW2-EW1)/EStep + 1

    if ((.not. SOC) .and. (.not. ROT)) then

    do ispin=1,NSpin

      open(334,file='tempT',status='unknown')
      if (LDOS_Beg <= LDOS_End ) open(333,file='tempDOS',status='unknown')

#ifdef PGI
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,cenergy,energy,green,gammar,gammal,T,temp) FIRSTPRIVATE(nsteps)
!$OMP DO SCHEDULE(STATIC,10)
#endif
       do n=1,nsteps
          energy=EW1+EStep*(n-1)
          cenergy=dcmplx(energy)

          !*********************************************************************
          !* Evaluation of the retarded "Green" function and coupling matrices *
          !*********************************************************************
          call gplus(cenergy,Green,GammaR,GammaL)

          if( .not. HTransm )then
             !*************************************************************
             !* Here we use the following non-Hermitian expression  for T *
             !* [Gamma_L G^a Gamma_R G^r]                                 *
             !* It works better for large clusters                        *
             !*************************************************************
             call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one, GammaL,NAorbs, Green,  NAOrbs, c_zero, T,    NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, T,     NAOrbs, GammaR, NAOrbs, c_zero, temp, NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, temp,  NAOrbs, Green,  NAOrbs, c_zero, T,    NAOrbs)
          else
             !********************************************************
             !* Here we use the following Hermitian expression for T *
             !* [Gamma_L^1/2 G^a Gamma_R G^r Gamma_L^1/2]            *
             !********************************************************
             call CMatPow(GammaL,0.5d0,temp)
             GammaL=temp
             call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one,GammaL,NAOrbs,Green, NAOrbs,c_zero,temp,NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,temp,  NAOrbs,GammaR,NAOrbs,c_zero,T,   NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,T,     NAOrbs,Green, NAOrbs,c_zero,temp,NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,temp,  NAOrbs,GammaL,NAOrbs,c_zero,T,   NAOrbs)
          end if

#ifdef PGI
!$OMP CRITICAL
#endif
          ! Mulliken DOS 
          if (LDOS_Beg <= LDOS_End ) then
            SG = matmul( SD, green )
            ! computing total DOS
            DOS=d_zero
            AtomDOS=d_zero
            do j=1,GetNAtoms()
            do i=LoAOrbNo(j),HiAOrbNo(j)
               AtomDOS(j)=AtomDOS(j)-dimag(SG(i,i))/d_pi
               DOS=DOS-dimag(SG(i,i))/d_pi
            end do
            end do

            if (dabs(energy-DOSEnergy) < EStep/2) AtomDOSEF(ispin,:)=AtomDOS

            ! print out DOS and atomic orbital resolved DOS ***
            imin = LoAOrbNo(LDOS_Beg)
            if( imin < 1 ) imin = 1
            imax = HiAOrbNo(LDOS_End)
            if( imax > NAOrbs ) imax = NAOrbs
            call flush(333)
            write(333,3333) energy,DOS*(-1)**(ispin+1),(AtomDOS(j)*(-1)**(ispin+1),j=LDOS_Beg,LDOS_End),(-dimag(SG(i,i))*(-1)**(ispin+1)/d_pi,i=imin,imax)
          end if

          ! computing transmission T
          ctrans=c_zero
          do i=1,NAOrbs
             ctrans=ctrans + T(i,i)
          end do

          if (dimag(ctrans).gt.1.0d-5) then
             write(ifu_log,*)'Transmission not real !!!'
             stop
          end if

          !Conductance in units of e^2/h

          if (NSpin == 1) trans=ctrans*2
          if (NSpin == 2) trans=ctrans

          ! Diagonalize the T matrix to get eigen channels
          if( NChannels > 0 )then
             if( HTransm ) then 
                call CHDiag( T, tn, info )
             else
                call CDiag( T, ctn, info )
                do i=1,NAOrbs
                  tn(i) = dble( ctn(i) )
                end do
                ! sort eigenvalues smallest to biggest
                call sort(NAOrbs,tn)
             end if
             if( n > 3 ) call SeparateSpaghettis( tchan1, tchan2, tn(NAOrbs-NChannels+1:NAOrbs), dummy, NChannels)
             tchan1=tchan2
             tchan2=tn(NAOrbs-NChannels+1:NAOrbs)
          end if

          call flush(334)
          write(334,1002)energy,trans,(tn(i),i=NAOrbs,NAOrbs-NChannels+1,-1)
          
#ifdef PGI
!$OMP END CRITICAL
#endif
       end do ! End of energy loop
#ifdef PGI
!$OMP END DO
!$OMP END PARALLEL
#endif

  ! Reordering in energy for nice DOS output
       if (LDOS_Beg <= LDOS_End ) then
       do n=1,nsteps
          energy=EW1+EStep*(n-1)
          rewind(333)
          do i=1,10000000000
          read(333,*)energ
          if (dabs(energy-energ) < 0.000001) then
             backspace(333)
             read(333,3333) (xxx(j),j=1,2+(LDOS_End-LDOS_Beg+1)+(imax-imin+1))
             write(ifu_dos,3333) (xxx(j),j=1,2+(LDOS_End-LDOS_Beg+1)+(imax-imin+1))
             exit
          end if
          end do
       end do
      write(ifu_dos,*)'   '
      close(333,status='delete')
      end if

  ! Reordering in energy for nice T output
      do n=1,nsteps
          energy=EW1+EStep*(n-1)
          rewind(334)
          do i=1,10000000000
          read(334,*)energ
          if (dabs(energy-energ) < 0.000001) then
             backspace(334)
             read(334,1002) (xxx(j),j=1,2+NChannels)
             write(ifu_tra,1002) (xxx(j),j=1,2+NChannels)
             exit
          end if
          end do
       end do
      write(ifu_tra,*)'   '
     !close(334,status='delete')
      close(334)

    end do ! End of spin loop

    else !SOC case
      
      open(334,file='tempT',status='unknown')
      if (LDOS_Beg <= LDOS_End ) open(333,file='tempDOS',status='unknown')
      
      wcount = 0  ! Issue warning of failure in transmission just once

#ifdef PGI
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,cenergy,energy,Dgreen,Dgammar,Dgammal,DT,Dtemp) 
!$OMP DO SCHEDULE(STATIC,10)
#endif
       do n=1,nsteps
          energy=EW1+EStep*(n-1)
          cenergy=dcmplx(energy)

          !*********************************************************************
          !* Evaluation of the retarded "Green" function and coupling matrices *
          !*********************************************************************
             call gplus_SOC(cenergy,DGreen,DGammaR,DGammaL,1)
             call zgemm('N','C',DNAOrbs,DNAOrbs,DNAOrbs,c_one, DGammaL,DNAOrbs, DGreen,  DNAOrbs, c_zero, DT,    DNAOrbs)
             call zgemm('N','N',DNAOrbs,DNAOrbs,DNAOrbs,c_one, DT,     DNAOrbs, DGammaR, DNAOrbs, c_zero, Dtemp, DNAOrbs)
             call zgemm('N','N',DNAOrbs,DNAOrbs,DNAOrbs,c_one, Dtemp,  DNAOrbs, DGreen,  DNAOrbs, c_zero, DT,    DNAOrbs)

#ifdef PGI
!$OMP CRITICAL
#endif

      ! Mulliken DOS 
           if (LDOS_Beg <= LDOS_End ) then
             DSG = matmul( S_SOC, DGreen )
             DOS=d_zero
             AtomDOS=d_zero
             do j=1,GetNAtoms()
             do i=LoAOrbNo(j),HiAOrbNo(j)
                AtomDOS(j)=AtomDOS(j)-dimag(DSG(i,i))/(2*d_pi)-dimag(DSG(i+NAOrbs,i+NAOrbs))/(2*d_pi)
                DOS=DOS-dimag(DSG(i,i))/(2*d_pi)-dimag(DSG(i+NAOrbs,i+NAOrbs))/(2*d_pi)
             end do
             end do

             if (dabs(energy-DOSEnergy) < EStep/2) AtomDOSEF(ispin,:)=AtomDOS

     ! print out DOS and atomic orbital resolved DOS ***
             imin = LoAOrbNo(LDOS_Beg)
             if( imin < 1 ) imin = 1
             imax = HiAOrbNo(LDOS_End)
             if( imax > NAOrbs ) imax = NAOrbs
             call flush(333)
             write(333,3333) energy,DOS,(AtomDOS(j),j=LDOS_Beg,LDOS_End),((-dimag(DSG(i,i))/(2*d_pi)-dimag(DSG(i+NAOrbs,i+NAOrbs)))/(2*d_pi),i=imin,imax)
           end if

          ! computing transmission T
          ctrans=c_zero
          do i=1,DNAOrbs
             ctrans=ctrans + DT(i,i)
          end do

          if (dimag(ctrans).gt.1.0d-5) then
             write(ifu_log,*)'Transmission not real !!!'
             stop
          end if
          trans=ctrans   ! in units of e^2/h

          ! Diagonalize the T matrix 
          ! to get eigen channels
          if( NChannels > 0 )then
             if( HTransm ) then 
                call CHDiag( DT, Dtn, info )
             else
                call CDiag( DT, Dctn, info )
                do i=1,DNAOrbs
                  Dtn(i) = dble( Dctn(i) )
                end do
                ! sort eigenvalues smallest to biggest
                call sort(DNAOrbs,Dtn)
             end if
             if( n > 3 ) call SeparateSpaghettis( tchan1, tchan2, Dtn(DNAOrbs-NChannels+1:DNAOrbs), dummy, NChannels)
             tchan1=tchan2
             tchan2=Dtn(DNAOrbs-NChannels+1:DNAOrbs)
          end if

     ! Computing polarization

          T_uu=0.0d0
          T_ud=0.0d0
          T_du=0.0d0
          T_dd=0.0d0

          do i=1,NAOrbs
          do j=1,NAOrbs
             GammaR_UU(i,j) = DGammaR(i,j)
             GammaR_DD(i,j) = DGammaR(i+NAOrbs,j+NAOrbs)
             GammaL_UU(i,j) = DGammaL(i,j)
             GammaL_DD(i,j) = DGammaL(i+NAOrbs,j+NAOrbs)
             Green_UU(i,j) = DGreen(i,j)
             Green_DD(i,j) = DGreen(i+NAOrbs,j+NAOrbs)
             Green_UD(i,j) = DGreen(i,j+NAOrbs)
             Green_DU(i,j) = DGreen(i+NAOrbs,j)
          end do
          end do

! up-up
          T=0.0d0
          temp=0.0d0
             call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one, GammaL_UU,NAOrbs, Green_UU,  NAOrbs, c_zero, T, NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, T,     NAOrbs, GammaR_UU, NAOrbs, c_zero, temp, NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, temp,  NAOrbs, Green_UU,  NAOrbs, c_zero, T,    NAOrbs)

             do i=1,NAOrbs
                T_uu = T_uu+REAL(T(i,i))
             end do
! up-down
          T=0.0d0
          temp=0.0d0
             call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one, GammaL_UU,NAOrbs, Green_UD,  NAOrbs, c_zero, T, NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, T,     NAOrbs, GammaR_DD, NAOrbs, c_zero, temp, NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, temp,  NAOrbs, Green_UD,  NAOrbs, c_zero, T,    NAOrbs)

             do i=1,NAOrbs
                T_ud = T_ud+REAL(T(i,i))
             end do
! down-up
          T=0.0d0
          temp=0.0d0
             call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one, GammaL_DD,NAOrbs, Green_DU,  NAOrbs, c_zero, T, NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, T,     NAOrbs, GammaR_UU, NAOrbs, c_zero, temp, NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, temp,  NAOrbs, Green_DU,  NAOrbs, c_zero, T,    NAOrbs)

             do i=1,NAOrbs
                T_du = T_du+REAL(T(i,i))
             end do
! down-down
          T=0.0d0
          temp=0.0d0
             call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one, GammaL_DD,NAOrbs, Green_DD,  NAOrbs, c_zero, T, NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, T,     NAOrbs, GammaR_DD, NAOrbs, c_zero, temp, NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, temp,  NAOrbs, Green_DD,  NAOrbs, c_zero, T,    NAOrbs)

          do i=1,NAOrbs
             T_dd = T_dd+REAL(T(i,i))
          end do

          polar = T_uu + T_du - T_dd - T_ud
          trans2 = T_uu + T_du + T_dd + T_ud

          
          if (dabs(trans2-trans) >= 1.0d-5 .and. wcount < 1) then
               if (SOC) print*,'Warning in the transmission with SOC'
               if (ROT) print*,'Warning in the transmission with spin rotations'
               wcount = wcount + 1
               !stop
           end if

          call flush(334)
          write(334,1002)energy,trans,polar,(Dtn(i),i=DNAOrbs,DNAOrbs-NChannels+1,-1)

#ifdef PGI
!$OMP END CRITICAL
#endif
       end do ! End of energy loop
#ifdef PGI
!$OMP END DO
!$OMP END PARALLEL
#endif

  ! Reordering in energy for nice DOS output
       if (LDOS_Beg <= LDOS_End ) then
       do n=1,nsteps
          energy=EW1+EStep*(n-1)
          rewind(333)
          do i=1,10000000000
          read(333,*)energ
          if (dabs(energy-energ) < 0.000001) then
             backspace(333)
             read(333,3333) (xxx(j),j=1,2+(LDOS_End-LDOS_Beg+1)+(imax-imin+1))
             write(ifu_dos,3333) (xxx(j),j=1,2+(LDOS_End-LDOS_Beg+1)+(imax-imin+1))
             exit
          end if
          end do
       end do
      close(333,status='delete')
      end if

  ! Reordering in energy for nice T output
      do n=1,nsteps
          energy=EW1+EStep*(n-1)
          rewind(334)
          do i=1,10000000000
          read(334,*)energ
          if (dabs(energy-energ) < 0.000001) then
             backspace(334)
             read(334,1002) (xxx(j),j=1,3+NChannels)
             write(ifu_tra,1002) (xxx(j),j=1,3+NChannels)
             exit
          end if
          end do
       end do
      close(334,status='delete')

  end if !End of SOC if

      if (SOC .or. ROT) then
         deallocate(DGammaL)
         deallocate(DGammaR)
         deallocate(DGreen)
         deallocate(DT)
         deallocate(Dtemp)
         deallocate(DSG)
         deallocate(Dtn)
         deallocate(Dctn)
         !deallocate(S_SOC)   ! DO NOT deallocate when calculating Mulliken population analysis with S_SOC!!!
         deallocate(H_SOC)
      else 
         deallocate(GammaL)
         deallocate(GammaR)
         deallocate(Green)
         deallocate(T)
         deallocate(temp)
         deallocate(SG)
         deallocate(tn)
         deallocate(ctn)
      end if

    if( NChannels > 0 ) then
       deallocate( tchan1, tchan2, dummy )
    end if

1002 format(f10.5,10000E14.5)
3333 format(f10.5,10000E14.5)

  end subroutine transmission
!--------------------------------------------------------------------------------
  subroutine CompHybFunc
  !
  ! Computes Hybridization functions for all correlated subspaces
  ! for all energy points defined in mesh.dat and writes them to 
  ! file
  !
    use constants
    use parameters, only: NCorrBl, CorrBeg, CorrEnd, eta
    use util
    use correlation
    use antcommon
!   USE IFLPORT
    use omp_lib
    implicit none

    integer :: ios, n, nmax, iblock, i, j, iao,jao
    real*8 :: En, Ep
    complex*16 :: zn
    character(len=3) :: istr
    character(len=100) :: fname

    real*8, dimension(:), allocatable :: EMesh
    
    ! Device Green's function
    complex*16, dimension(NAorbs,NAOrbs) :: GD

    complex*16, dimension(:,:,:,:,:), allocatable :: delta

    print *, "-----------------------------------------"
    print *, "--- Computing Hybridization functions ---"
    print *, "-----------------------------------------"

    call SetHamOvl( HD, SD )
    !
    ! Read mesh file mesh.dat
    ! 
    print *, "Reading energy mesh from file mesh.dat" 
    ! open mesh file
    open(unit=ifu_msh,file='mesh.dat',status='old',iostat=ios)
    if( ios /= 0 )then
       print *, "Device/CompHybFunc/Error: Could not open energy mesh file mesh.dat. Abort."
       STOP
    end if
    ios = 0; nmax=0; Ep=-1e+10
    do while( ios == 0 )
       read(unit=ifu_msh,fmt=*,iostat=ios), En
       if( ios /= 0 ) exit
       if( En <= Ep ) exit
       Ep = En
       nmax=nmax+1
       !print *, nmax, En
    end do
    print *, "Mesh file has", nmax, " data points."
    allocate( EMesh( nmax ) )
    rewind(ifu_msh)
    do n=1,nmax
       read(unit=ifu_msh,fmt=*,iostat=ios), EMesh(n)
       if( ios /= 0 ) exit
    end do
    close(ifu_msh)
    !
    ! calculate hybridization function for all mesh points
    ! 
    allocate( delta(nmax,NSpin,NCorrBl,NMaxCorr,NMaxCorr) ) 
    delta = c_zero

    do ispin=1,NSpin
!$OMP PARALLEL PRIVATE(En,zn,GD)
!$OMP DO
       do n=1,nmax
          En = EMesh(n)
          zn = En+ui*eta
          call gplus0(zn,GD)
!$OMP CRITICAL
          call CompDelta( ispin, En, -shift, GD, delta(n,ispin,:,:,:) )
!$OMP END CRITICAL
       end do
!$OMP END DO
!$OMP END PARALLEL
    end do
    !
    ! write hybridization functions to files
    !
    do iblock=1,NCorrBl
       print *, "Correlated block", iblock
       call int2str( iblock, istr )
       fname='delta.out.'//istr
       print *, "Output file for Hybridization:", fname
       open(unit=ifu_hyb,file=fname,status='unknown',iostat=ios)
       fname='Ac.out.'//istr 
       print *, "Output file for Bath Sepctral function:", fname
       open(unit=ifu_ac,file=fname,status='unknown',iostat=ios)
       do n=1,nmax
          En=EMesh(n)
          !
          ! write hybridization function Delta
          !
          write(ifu_hyb,fmt='(E20.10,1000E14.6)'),& 
               En, ( (delta(n,ispin,iblock,i,i),ispin=1,NDSpin),i=1,ncorrao(iblock) )
          call flush(ifu_hyb)
          !
          ! write bath spectral function = -Im(Delta)/pi
          !
          write(ifu_ac,fmt='(E20.10,1000E14.6)'), &
               En, ( (-AIMAG(delta(n,ispin,iblock,i,i))/d_pi,ispin=1,NDSpin),i=1,ncorrao(iblock) )
          call flush(ifu_ac)
       end do
       close(ifu_hyb)
       close(ifu_ac)
    end do
    print *, "done." 

    deallocate( EMesh, delta )
    
  end subroutine CompHybFunc


  !*******************************************!
  ! Routine for orbital eigen-channel         !
  ! analysis with reduced transmission matrix !
  !*******************************************!
  subroutine EigenChannelAnalysis
    use parameters, only:  RedTransmB,RedTransmE, eta, EW1, EW2,EStep
    use Cluster, only : hiaorbno, loaorbno
    use constants, only: ui, d_pi
    use numeric, only: CInv, CMatPow, CHDiag
    implicit none
    
    real*8 :: rho, phi,DE,energy,delta,mindelta,tsave
    real*8, parameter :: dsmall = 1.0d-10
    integer :: NMaxEV = 10 ! At how many points to print out eign vectors: 2*NMaxEV+1

    integer :: NSD, NLD, NRD, NS1, NS2, N,NMax, info

    complex*16,dimension(NAOrbs,NAOrbs) :: SigmaL, SigmaR
    complex*16,dimension(:,:),allocatable :: GSD, SigmaLD, SigmaRD, GammaLD, GammaRD, GammaLDph, TS, TS0
    complex*16,dimension(:,:),allocatable :: gLD, gRD
    complex*16,dimension(:,:),allocatable :: VLS, VRS

    complex*16 :: cenergy, ci
    integer :: is, i, nchan, j, jmin, k

    real*8,dimension(:),allocatable :: tchan,tchan1,tchan2

    ! Dimension of scattering region inside device
    NSD = HiAOrbNo(RedTransmE)-LoAOrbNo(RedTransmB)+1
    ! Dimension of left region of device
    NLD = LoAOrbNo(RedTransmB)-1
    ! Dimension of right region of device
    NRD = NAOrbs-HiAOrbNo(RedTransmE)
    
    NS1 = LoAOrbNo(RedTransmB)
    NS2 = HiAOrbNo(RedTransmE)

    allocate( GSD(NSD,NSD), SigmaLD(NSD,NSD), SigmaRD(NSD,NSD), &
         GammaLD(NSD,NSD), GammaRD(NSD,NSD), GammaLDph(NSD,NSD), TS(NSD,NSD), TS0(NSD,NSD), &
         gLD(NLD,NLD), gRD(NRD,NRD), VLS(NLD,NSD), VRS(NRD,NSD), tchan(NSD), tchan1(NSD), tchan2(NSD) )

    !open(ifu_red,file='t.dat',status='unknown')

    print *
    print *, "-----------------------------------------------------------------------"
    print *, "--- Orbital Eigen-channel Analysis with reduced Transmission matrix ---"
    print *, "-----------------------------------------------------------------------"
    print *
    print *, "Begin of scattering region: AO ", NS1
    print *, "End of scattering region:   AO ", NS2
    print *

    do is=1,NSpin
       if(NSpin ==2 .and. is==1)print*,"Spin UP"
       if(NSpin ==2 .and. is==2)print*,"Spin DOWN"       

       NMax = int(EW2/EStep)

       !DE = d_zero
       !IF(NMaxE>0) DE = EWindow/DBLE(NMaxE)

       do N = -NMax,NMax!-NMaxE,NMaxE  
          energy = N*EStep
          cenergy = energy

          call CompSelfEnergies(is,cenergy,SigmaL,SigmaR,1)
       
          do i=1,NLD
             do j=1,NLD
                gLD(i,j) = (energy-shift+ui*eta)*SD(i,j) - HD(is,i,j) - SigmaL(i,j)
             end do
          end do

          do i=1,NRD
             do j=1,NRD
                gRD(i,j) = (energy-shift+ui*eta)*SD(NS2+i,NS2+j) - HD(is,NS2+i,NS2+j) - SigmaR(NS2+i,NS2+j)
             end do
          end do

          info = CInv( gLD )
          info = CInv( gRD )
          
          VLS = HD(is,1:NLD, NS1:NS2) - (energy-shift)*SD(1:NLD, NS1:NS2)
          VRS = HD(is,NS2+1:NAOrbs,NS1:NS2) - (energy-shift)*SD(NS2+1:NAOrbs,NS1:NS2)
          
          SigmaLD = matmul( conjg(transpose(VLS)), matmul( gLD, VLS ) )
          SigmaRD = matmul( conjg(transpose(VRS)), matmul( gRD, VRS ) )
          
          GSD = (energy-shift+ui*eta)*SD(NS1:NS2,NS1:NS2) & 
               - HD(is,NS1:NS2,NS1:NS2) - SigmaLD - SigmaRD

          info = CInv( GSD )

          GammaLD = ui*(SigmaLD-conjg(transpose(SigmaLD)))
          GammaRD = ui*(SigmaRD-conjg(transpose(SigmaRD)))

          call CMatPow( GammaLD, 0.5d0, GammaLDph )
          
          ![Gamma_L^1/2 G^a Gamma_R G^r Gamma_L^1/2] 
          TS = matmul( GammaLDph, conjg(transpose(GSD)) )
          TS = matmul( TS, GammaRD )
          TS = matmul( TS, GSD )
          TS = matmul( TS, GammaLDph )

          call CHDiag( TS, tchan, info )
          if( info /= 0 )then
             print*, "Error diagonalizing Reduced transmission matrix: info=", info
             return
          end if

          ! *** Ordering of eigenchannels ***
          if( N >= -NMax+2 ) call SeparateSpaghettis( tchan1, tchan2, tchan, TS, NSD )
          tchan1=tchan2
          tchan2=tchan

          write(ifu_red,'(F10.5,100E14.5)'), energy,(tchan(i),i=NSD,1,-1)

          if( N ==0 )then 
             print *, "Eigenchannel composition at Fermi level:"
             print *, "----------------------------------------"
             do nchan=NSD,1,-1
                print '(A,I2,A,F9.3)',"Channel ", NSD-nchan+1, ": Transmission = ", tchan(nchan)
                do i=1,NSD
                   ci = TS(i,nchan)
                   rho = abs(ci)
                   if( abs(real(ci)) < dsmall * abs(aimag(ci)) )then
                      phi = sign(0.5*d_pi,aimag(ci))
                   else
                      phi = atan( aimag(ci)/real(ci) )
                      if( real(ci) < 0 .and. aimag(ci) > 0 ) phi = phi + 0.5*d_pi
                      if( real(ci) < 0 .and. aimag(ci) < 0 ) phi = phi - 0.5*d_pi
                   end if
                   print '(I4,A,F8.4,F8.4)', i, " : ", rho, phi
                end do
             end do
          end if

       end do
       write(ifu_red,*)
       print *, " "
    end do

    deallocate( GSD, SigmaLD, SigmaRD, GammaLD, GammaRD, GammaLDph, TS, TS0, gLD, gRD, VLS, VRS , tchan, tchan1, tchan2 )

    !close(ifu_red)

  end subroutine EIGENCHANNELANALYSIS


  !
  ! *** Subroutine to separate individual eigen channel     ***
  ! *** transmissions (= spaghettis) using first derivative ***
  !
  subroutine SeparateSpaghettis( evals1, evals2, evals, evecs, N )
    implicit none
    
    ! eigen values at before last, last and actual energy
    real*8, dimension(N) :: evals1, evals2, evals
    ! eigenvectors at actual energy
    complex*16, dimension(N,N) :: evecs
    ! number of eigenvectors
    integer :: N
    
    integer :: i, j, jmin
    real*8 :: yex, delta, mindelta, valsave

    complex*16, dimension(N) :: vecsave

    do i=1,N

       ! extrapolate actual eigenvalue from last 2 eigenvalues
       yex = 2.0d0*evals2(i)-evals1(i)

       ! Find actual eigenvalue which deviates minimally 
       ! from extrapolated value 
       mindelta = abs(yex-evals(i))
       jmin = i
       do j=i+1,N
          delta =  abs(yex-evals(j))
          if( delta < mindelta )then 
             jmin = j 
             mindelta = delta
          end if
       end do

       ! change eigenvector i with eigenvector jmin
       if( jmin /= i )then 
          vecsave = evecs(:,jmin)
          evecs(:,jmin)=evecs(:,i)
          evecs(:,i)=vecsave
          valsave = evals(jmin)
          evals(jmin) = evals(i)
          evals(i) = valsave
       end if
    end do
    
  end subroutine SeparateSpaghettis
  

  ! ***********************************
  ! Compute Self energies of electrodes
  ! ***********************************
  subroutine CompSelfEnergies( spin, cenergy, Sigma1, Sigma2, sgn )
    use parameters, only: ElType, DD, UD, DU, Overlap
    use BetheLattice, only: CompSelfEnergyBL, LeadBL 
    use OneDLead, only: CompSelfEnergy1D, Lead1D
    use constants
#ifdef PGI
    use lapack_blas, only: zgemm
#endif
    implicit none
    external zgemm
    
    integer, intent(in) :: spin,sgn
    complex*16, intent(in) :: cenergy
    complex*16, dimension(:,:),intent(inout) :: Sigma1
    complex*16, dimension(:,:),intent(inout) :: Sigma2
    complex*16, dimension(NAOrbs,NAOrbs) :: temp
    integer :: is,omp_get_thread_num
    
    !Sigma1=(0.0,0.0)
    !Sigma2=(0.0,0.0)

    is = spin
    if( DD .and. spin == 1 ) is=2
    if( DD .and. spin == 2 ) is=1
    if( UD .and. spin == 1 ) is=1
    if( UD .and. spin == 2 ) is=2
    if( DU .and. spin == 1 ) is=2
    if( DU .and. spin == 2 ) is=1
    !write(ifu_log,*)omp_get_thread_num(),'in CompSelfEnergies',shift,cenergy,'lead 1'
    select case( ElType(1) )
    case( "BETHE" ) 
       call CompSelfEnergyBL( LeadBL(1), is, cenergy, Sigma1, sgn )
       ! transform to non-orthogonal basis
       if( Overlap < -0.01 .and. .not. HDOrtho )then
          ! Sigma1 -> S^1/2 * Sigma1 * S^1/2
          call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, Sigma1,NAorbs, SPH,  NAOrbs, c_zero, temp,   NAOrbs)
          call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, SPH,   NAorbs, temp, NAOrbs, c_zero, Sigma1, NAOrbs)
       endif
    case( "1DLEAD" )
       call CompSelfEnergy1D( Lead1D(1), is, cenergy, Sigma1 )
    case( "GHOST" )
        continue
    end select

    if( DD .and. spin == 1 ) is=2
    if( DD .and. spin == 2 ) is=1
    if( UD .and. spin == 1 ) is=2
    if( UD .and. spin == 2 ) is=1
    if( DU .and. spin == 1 ) is=1
    if( DU .and. spin == 2 ) is=2
    
    !write(ifu_log,*)omp_get_thread_num(),'in CompSelfEnergies',shift,cenergy,'lead 2'
    select case( ElType(2) )
    case( "BETHE" )
       call CompSelfEnergyBL( LeadBL(2), is, cenergy, Sigma2, sgn )
       ! transform to non-orthogonal basis 
       if( Overlap < -0.01 .and. .not. HDOrtho )then
          ! Sigma2 -> S^1/2 * Sigma2 * S^1/2
          call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, Sigma2,NAorbs, SPH,  NAOrbs, c_zero, temp,   NAOrbs)
          call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, SPH,   NAorbs, temp, NAOrbs, c_zero, Sigma2, NAOrbs)
       endif
    case( "1DLEAD" )
       call CompSelfEnergy1D( Lead1D(2), is, cenergy, Sigma2 )
    case( "GHOST" )
        continue
    end select
  end subroutine CompSelfEnergies


  !
  ! *** Integrand for charge integration on complex contour ***
  !
  real*8 function DDOS( E0, R, phi )
    use constants, only: c_zero, ui, d_pi

    real*8, intent(in) :: phi, E0, R
    integer :: i, j !!, ispin 
    complex*16,dimension(NAOrbs,NAOrbs) :: green,gammar,gammal
    complex*16 :: TrGS, z

    z = E0 - R*(cos(phi) - ui*sin(phi)) 

    TrGS=c_zero
    do ispin=1,NSpin
       call gplus(z,green,gammar,gammal)
       !do i=1,NAOrbs
       do i=NCDAO1,NCDAO2
          do j=1,NAOrbs
             TrGS = TrGS + green(i,j)*SD(j,i)
          end do
       end do
    end do
    ! Account for spin degeneracy
    if(NSpin==1) TrGS = TrGS * 2.0d0
    DDOS = -DIMAG(R*(sin(phi)+ui*cos(phi))*TrGS)/d_pi 
  end function DDOS

  ! 
  ! *** Total charge up to energy E, lower bound is EMin ***
  ! 
  real*8 function TotCharge( E )
    use constants, only: d_zero, d_pi
    use parameters, only: ChargeAcc
    use numeric, only: gauleg 
    implicit none
    
    real*8, intent(in) :: E

    integer, parameter :: nmax = 2047
    
    real*8 :: q,qq, E0, R, w_j, phi_j
    integer :: n, n1, n2, j, i
    
    real*8, dimension(nmax) :: x, w
    
    ! Integration contour parameters:
    E0 = 0.5*(E + EMin)
    R  = 0.5*(E - EMin)
    ! Computing integral of DOS over 
    ! complex contour using Gauss-Legendre 
    ! quadrature
    n=1
    do  
       q = d_zero
       do j=1,n
          call gauleg(d_zero,d_pi,x(1:n),w(1:n),n)
          q = q + w(j)*DDOS( E0, R, x(j) )
       end do
       if( n > 1 .and. (q == d_zero .or. abs(q-qq) < ChargeAcc*NCDEl ) ) exit  
       n=2*n+1
       if( n > nmax )then
          print *, "TotCharge/gaussian quadrature has not converged after", nmax, " steps."
          TotCharge = 2.0d0*(NCDAO2-NCDAO1+1) - 10.0d0*ChargeAcc*NCDEl
          return
       end if
       qq = q
    end do
    !print *, "gaussian quadrature converged after", n, " steps. Error:", abs(q-qq)

    TotCharge = q
  end function TotCharge

  ! *************************************
  !  Estimates upper/lower energy 
  !  boundary EMin/EMax,
  !  above/below which DOS is gauranteed 
  !  to be zero
  ! *************************************
  subroutine FindEnergyBounds
    use parameters, only: ChargeAcc,eta

    integer :: i, cond,k 
    real*8 :: EStep, Q
    
    print *, "Searching energy boundaries [EMin, EMax] such that"
    print '(A,I4)', " Total Integrated Spectral Weight (TISW)=", 2*(NCDAO2-NCDAO1+1)

    EStep = 10.0d0 +10000.0 *eta

    do
       Q = TotCharge( EMax )
       print '(A,F15.3,A,F15.3,A,F12.5)', " EMin=", EMin, "  EMax=", EMax , "  TISW=", Q
       if( abs(Q - 2.0d0*(NCDAO2-NCDAO1+1)) < ChargeAcc*NCDEl*10.0 ) then
          exit
       end if
       EMin = EMin - EStep
       EMax = EMax + EStep
    end do
    print *, "--------------------------------------------------"

  end subroutine FindEnergyBounds

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
!    second kind                                                               c
!        eps: Tolerance                                                        c
!        M:   On input, maximum number of points allowed                       c
!             On output, 0 for an alleged successfull calculation, 1 otherwise c
!        F(): External function to be integrated.                              c
!        CH:  The value of the integral. Interval [-1,1]                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine IntRealAxis(Er,El,M)

    use constants, only: ui,d_pi,d_zero
    use parameters, only: PAcc 
!   USE IFLPORT
    use omp_lib
   
    real*8,intent(in) :: Er,El
    integer,intent(inout) :: M
    real*8, dimension(M) :: xs,xcc
    complex*16, dimension(NAOrbs,NAOrbs) :: green
    complex*16, dimension(NAOrbs,NAOrbs) :: greenp,greenm 
    complex*16, dimension(NAOrbs,NAOrbs) :: p,q
    complex*16, dimension(NAOrbs,NAOrbs) :: PDP
    complex*16 :: E0,Em,Ep
    integer :: n,i,j,l,k,k1
    real*8 :: pi,S0,c0,rchp,rchq,xp,c1,s1,s,cc,x,xx,achp,achq

      pi=d_pi

! Initializing M, n, S0, C0, CH and p

      M = (M-1)*0.5d0
      n = 1
      S0=1
      C0=0
      E0=edex3(El,Er,d_zero)
      call glesser(E0,green)
      do i=1,NAOrbs
       do j=1,NAOrbs
        PDOUT(ispin,i,j) = -ui*green(i,j)/(2*pi)
        p(i,j) = PDOUT(ispin,i,j)
       enddo
      enddo
! Computing the (2n+1) points quadrature formula ...
! ... updating q, p, C1, S1, C0, S0, s and c
1     continue
      do i=1,NAOrbs
       do j=1,NAOrbs
         q(i,j) = 2*p(i,j)
         p(i,j) = 2*PDOUT(ispin,i,j)
       enddo
      enddo
      C1 = C0
      S1 = S0
      C0 = sqrt((1+C1)*0.5d0)
      S0 = S1/(2*C0)
      !s = S0
      !cc = C0
      xs(1) = S0
      xcc(1) = C0
      do l=1,n,2
         xs(l+2)=xs(l)*C1+xcc(l)*S1
         xcc(l+2)=xcc(l)*C1-xs(l)*S1
      end do
! ... computing F() at the new points
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,xx,Em,Ep,greenp,greenm,i,j,pdp)
       PDP=d_zero
!$OMP DO SCHEDULE(STATIC,1)
      do l=1,n,2
        xx = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
        Em=edex3(El,Er,-xx)
        Ep=edex3(El,Er,xx)
!!$OMP  PARALLEL DEFAULT(SHARED)
!!$OMP  SECTIONS
!!$OMP  SECTION
       call glesser(Em,greenm)
!!$OMP  SECTION
       call glesser(Ep,greenp)
!!$OMP  END SECTIONS
!!$OMP  END PARALLEL
          do i=1,NAOrbs
           do j=1,NAOrbs
            pdp(i,j) = pdp(i,j)-ui*(greenm(i,j)+greenp(i,j))*xs(l)**4/(2*pi)
           enddo
          enddo
      enddo
!$OMP END DO
!$OMP CRITICAL
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             PDOUT(ispin,i,j)=PDOUT(ispin,i,j)+PDP(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL

! ... replacing n by 2n+1
         n = n + n + 1
! Stopping?
      do i=1,NAOrbs
       do j=1,NAOrbs
        rCHp=dble(PDOUT(ispin,i,j)-p(i,j))
        aCHp=dimag(PDOUT(ispin,i,j)-p(i,j))
        rCHq=dble(PDOUT(ispin,i,j)-q(i,j))
        aCHq=dimag(PDOUT(ispin,i,j)-q(i,j))
        if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc.and.n.le.M) goto 1
        if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc.and.n.le.M) goto 1
       enddo
      enddo
! Test for successfullness and integral final value
      M = 0
      do i=1,NAOrbs
      do j=1,NAOrbs
        rCHp=dble(PDOUT(ispin,i,j)-p(i,j))
        aCHp=dimag(PDOUT(ispin,i,j)-p(i,j))
        rCHq=dble(PDOUT(ispin,i,j)-q(i,j))
        aCHq=dimag(PDOUT(ispin,i,j)-q(i,j))
        if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc) M = 1
        if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc) M = 1
        PDOUT(ispin,i,j) = 16*PDOUT(ispin,i,j)/(3*(n+1))
        PDOUT(ispin,i,j) = PDOUT(ispin,i,j)*(El-Er)/2
      enddo
      enddo
      if (M == 1) write(ifu_log,'(A)')'Unsuccesful integration'
      write(ifu_log,'(A,i5,A)')' Integration of the non-equilibrium density matrix has needed ',(((n-1)/2)+1)/2, ' points'

      return
    end subroutine IntRealAxis

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
!    second kind                                                               c
!        eps: Tolerance                                                        c
!        M:   On input, maximum number of points allowed                       c
!             On output, 0 for an alleged successfull calculation, 1 otherwise c
!        F(): External function to be integrated.                              c
!        CH:  The value of the integral. Interval [-1,1]                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine IntRealAxis_SOC(Er,El,M)

    use constants, only: ui,d_pi,d_zero
    use parameters, only: PAcc 
!   USE IFLPORT
    use omp_lib
   
    real*8,intent(in) :: Er,El
    integer,intent(inout) :: M
    real*8, dimension(M) :: xs,xcc
    complex*16, dimension(DNAOrbs,DNAOrbs) :: green
    complex*16, dimension(DNAOrbs,DNAOrbs) :: greenp,greenm 
!   complex*16, dimension(DNAOrbs,DNAOrbs) :: p,q
    complex*16, dimension(DNAOrbs,DNAOrbs) :: PDP
    complex*16 :: E0,Em,Ep
    integer :: n,i,j,l,k,k1
    real*8 :: pi,S0,c0,rchp,rchq,xp,c1,s1,s,cc,x,xx,achp,achq,CH,q,CHI

      pi=d_pi

! Initializing M, n, S0, C0, CH and p

      M = (M-1)*0.5d0
      n = 1
      S0=1
      C0=0
      E0=edex3(El,Er,d_zero)
     call glesser_SOC(E0,green)

    CH = 0.d0
      do i=1,DNAOrbs
       do j=1,DNAOrbs
        PDOUT_SOC(i,j) = -ui*green(i,j)/(2*pi)
        CH = CH + REAL(PDOUT_SOC(i,j)*S_SOC(j,i))
       enddo
      enddo
 !print*,'CH',CH
! Computing the (2n+1) points quadrature formula ...
! ... updating q, p, C1, S1, C0, S0, s and c

    xp = CH
1   q = xp + xp
    xp = CH + CH

      C1 = C0
      S1 = S0
      C0 = sqrt((1+C1)*0.5d0)
      S0 = S1/(2*C0)
      !s = S0
      !cc = C0
      xs(1) = S0
      xcc(1) = C0
      do l=1,n,2
         xs(l+2)=xs(l)*C1+xcc(l)*S1
         xcc(l+2)=xcc(l)*C1-xs(l)*S1
      end do
! ... computing F() at the new points
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,xx,Em,Ep,greenp,greenm,i,j,pdp)
       PDP=d_zero
!$OMP DO SCHEDULE(STATIC,1)
      do l=1,n,2
        xx = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
        Em=edex3(El,Er,-xx)
        Ep=edex3(El,Er,xx)
!!$OMP  PARALLEL DEFAULT(SHARED)
!!$OMP  SECTIONS
!!$OMP  SECTION
       call glesser_SOC(Em,greenm)
!!$OMP  SECTION
       call glesser_SOC(Ep,greenp)
!!$OMP  END SECTIONS
!!$OMP  END PARALLEL
          do i=1,DNAOrbs
           do j=1,DNAOrbs
            pdp(i,j) = pdp(i,j)-ui*(greenm(i,j)+greenp(i,j))*xs(l)**4/(2*pi)
           enddo
          enddo
      enddo
!$OMP END DO
!$OMP CRITICAL
       do i = 1,DNAOrbs
          do j = 1,DNAOrbs
             PDOUT_SOC(i,j)=PDOUT_SOC(i,j)+PDP(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL

! ... replacing n by 2n+1
         n = n + n + 1

    CH = 0.d0
    CHI = 0.d0
    do i=1,DNAOrbs
       do j=1,DNAOrbs
          CH = CH + REAL(PDOUT_SOC(i,j)*S_SOC(j,i))
          CHI = CHI + DIMAG(PDOUT_SOC(k,k1)*S_SOC(k1,k))
       end do
    enddo
   !print*,n,ch,chi
      !print*, n,16*CH*(El-Er)/(6*(n+1))
    ! Stopping?
   !print*, n,(CH-xp)*(CH-xp)*16-3*(n+1)*abs(CH-q)*PAcc
     if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc.and.n.le.M) goto 1
    ! Test for successfullness and integral final value
    M = 0
     if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc) M = 1

    do i=1,DNAOrbs
       do j=1,DNAOrbs
          PDOUT_SOC(i,j) = 16*PDOUT_SOC(i,j)/(3*(n+1))
          PDOUT_SOC(i,j) = PDOUT_SOC(i,j)*(El-Er)/2
       enddo
    enddo


      if (M == 1) write(ifu_log,'(A)')'Unsuccesful integration'
      write(ifu_log,'(A,i5,A)')' Integration of the non-equilibrium density matrix has needed ',(((n-1)/2)+1)/2, ' points'

      return
    end subroutine IntRealAxis_SOC     

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !c    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
  !c second kind                                                                  c
  !c        eps: Tolerance                                                        c
  !c        b: parameter of the change of variable                                c
  !c        Em: maximum value of the energy range                                 c
  !c        M:   On input, maximum number of points allowed                       c
  !c             On output, 0 for an alleged successfull calculation, 1 otherwise c
  !c        dn:  On output, the density matrix                                    c
  !c        CH:  On output, the value of the integral (charge density).           c
  !c             Interval [-1,1]                                                  c
  !c        Eq:  On output, the value of the upper bound of the integral.         c
  !c        The rest of arguments are neeed by the subrtn. gplus                  c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine IntCompPlane(rrr,bi,Emi,M,Eq)
    use parameters, only: PAcc 
    use constants, only: d_pi, d_zero, ui
!   USE IFLPORT
    use omp_lib

    implicit none

    real*8,intent(in) :: rrr, bi, Emi, Eq
    integer,intent(inout) :: M
    real*8, dimension(NAOrbs,NAOrbs) :: PDP

    real*8 :: a,b,Em,S0,c0,x0,er0,der0,ch,xp,q,c1,s1,s,cc,x,erp,erm
    integer :: n,i,j,l,k,k1,chunk!,omp_get_thread_num,omp_get_num_threads
    real*8, dimension(M) :: xs,xcc

    complex*16 :: E0,E
    complex*16, dimension(2) :: EE
    complex*16, dimension(NAOrbs,NAOrbs) :: green
    complex*16, dimension(2,NAOrbs,NAOrbs) :: greenn 

   !logical :: omp_get_nested
 
    a = 1.d0/d_pi
    b = bi
    Em = Emi
    PD(ispin,:,:) = d_zero
    M = (M-1)*0.5
    n = 1
    S0 = 1
    C0 = 0
    x0 = 0.d0
    er0 = edex3(Em,b,x0)
    der0 = 0.5d0*(Em-b)
    E0 = rrr*exp(ui*er0)-rrr+Eq
    call gplus0(E0,green)
   !print*,green(1,1)
   !print*,green(1,3),green(3,1)
   
   !print*,'E0,er0,der0,rrr',E0,er0,der0,rrr            
 CH = 0.d0
    do i = 1,NAOrbs
       do j =1,NAOrbs
          PD(ispin,i,j)= a*dimag(ui*rrr*exp(ui*er0)*green(i,j))*der0
          CH = CH + PD(ispin,i,j)*SD(j,i)
       enddo
    enddo
    !print*,PD(ispin,1,1)
    !print*,PD(ispin,1,3)
    !print*,PD(ispin,3,1)

   !print*,'ch', 16*CH/(3*(n+1))
    xp = CH
1   q = xp + xp
    xp = CH + CH
    C1 = C0
    S1 = S0
    C0 = sqrt((1+C1)*0.5d0)
    S0 = S1/(C0+C0)
    xs(1) = S0
    xcc(1) = C0
    do l=1,n,2
       xs(l+2)=xs(l)*C1+xcc(l)*S1
       xcc(l+2)=xcc(l)*C1-xs(l)*S1
    end do
    !call omp_set_nested(.true.)
    !call omp_set_num_threads(2)
    !print *, omp_get_nested()
    !print *,'--------------------------' 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,x,erp,erm,EE,greenn,i,j,pdp)
    PDP=d_zero
     chunk=max(((n+1)/2)/omp_get_num_threads(),1)
    !chunk=1
!$OMP DO SCHEDULE(STATIC,chunk)
    do l=1,n,2
       !write(ifu_log,*)'thread',omp_get_thread_num(),'l=',l
       x = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
       erp = 0.5d0*((Em-b)*x + (Em+b))
       erm = 0.5d0*(-(Em-b)*x + (Em+b))
       EE(1) = rrr*exp(ui*erp)-rrr+Eq
       EE(2) = rrr*exp(ui*erm)-rrr+Eq
       do k=1,2
          !print *,'l',l,'k',k,omp_get_thread_num()
          call gplus0(EE(k),greenn(k,:,:))
       end do
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             PDP(i,j) = PDP(i,j)+ a*(dimag(ui*rrr*exp(ui*erp)*greenn(1,i,j))*der0 &
                  &   +dimag(ui*rrr*exp(ui*erm)*greenn(2,i,j))*der0)*xs(l)**4
          end do
       end do
    end do
!$OMP END DO
!$OMP CRITICAL
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             PD(ispin,i,j)=PD(ispin,i,j)+PDP(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL
    CH = 0.d0
    do k=1,NAOrbs
       ! Non-orthogonal basis: Ch = Tr[P*SD]   
       do k1=1,NAOrbs
          CH = CH + PD(ispin,k,k1)*SD(k1,k)
       end do
    enddo
   !print*,'ch',n,  16*CH/(3*(n+1))
    ! ... replacing n by 2n+1
    n = n + n + 1
    ! Stopping?
     if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc.and.n.le.M) goto 1
    ! Test for successfullness and integral final value
    M = 0
     if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc) M = 1
    CH = 16*CH/(3*(n+1))
    do k=1,NAOrbs
       do l=1,NAOrbs
          PD(ispin,k,l) = 16*PD(ispin,k,l)/(3*(n+1))
       enddo
    enddo
   !print*,'ch',n, CH
    write(ifu_log,'(A,I4,A)')' Integration of the density matrix on the complex plane has needed ',(((n-1)/2)+1)/2,' points'

    return
  end subroutine IntCompPlane

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !c    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
  !c second kind                                                                  c
  !c        eps: Tolerance                                                        c
  !c        b: parameter of the change of variable                                c
  !c        Em: maximum value of the energy range                                 c
  !c        M:   On input, maximum number of points allowed                       c
  !c             On output, 0 for an alleged successfull calculation, 1 otherwise c
  !c        dn:  On output, the density matrix                                    c
  !c        CH:  On output, the value of the integral (charge density).           c
  !c             Interval [-1,1]                                                  c
  !c        Eq:  On output, the value of the upper bound of the integral.         c
  !c        The rest of arguments are neeed by the subrtn. gplus                  c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine IntCompPlane_SOC(sgn,rrr,bi,Emi,M,Eq)
    use parameters, only: PAcc 
    use constants, only: d_pi, d_zero, ui, c_zero
!   USE IFLPORT
    use omp_lib

    implicit none

    integer,intent(inout) :: M
    integer,intent(in) :: sgn
    integer :: n,i,j,l,k,k1,chunk!,omp_get_thread_num,omp_get_num_threads

    real*8 :: a,b,Em,S0,c0,x0,er0,der0,ch,xp,q,c1,s1,s,cc,x,erp,erm
    real*8,intent(in) :: rrr, bi, Emi, Eq
    real*8, dimension(M) :: xs,xcc

    complex*16 :: E0,E
    complex*16, dimension(2) :: EE
    complex*16, dimension(DNAOrbs,DNAOrbs) :: green,PDP
    complex*16, dimension(2,DNAOrbs,DNAOrbs) :: greenn

   !logical :: omp_get_nested
 
    a = 1.0/(2.0d0*d_pi)
    b = bi
    Em = Emi
    PD_SOC = c_zero
    PDP = c_zero
    M = (M-1)*0.5
    n = 1
    S0 = 1
    C0 = 0
    x0 = 0.d0
    er0 = edex3(Em,b,x0)
    der0 = 0.5d0*(Em-b)
    E0 = rrr*exp(ui*er0)-rrr+Eq
       
    call gplus0_SOC(E0,green,sgn)

    CH = 0.d0

    do i = 1,DNAOrbs
       do j =1,DNAOrbs
          PD_SOC(i,j)= PD_SOC(i,j) + a*rrr*exp(ui*er0)*green(i,j)*der0
          CH = CH + real(PD_SOC(i,j)*S_SOC(j,i))
       enddo
    enddo
    
   !print*,'ch', 16*CH/(3*(n+1))    
    
    xp = CH
1   q = xp + xp
    xp = CH + CH
    C1 = C0
    S1 = S0
    C0 = sqrt((1+C1)*0.5d0)
    S0 = S1/(C0+C0)
    xs(1) = S0
    xcc(1) = C0
    do l=1,n,2
       xs(l+2)=xs(l)*C1+xcc(l)*S1
       xcc(l+2)=xcc(l)*C1-xs(l)*S1
    end do
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,x,erp,erm,EE,greenn,i,j,PDP)
     PDP=c_zero
     chunk=max(((n+1)/2)/omp_get_num_threads(),1)
    !chunk=1
!$OMP DO SCHEDULE(STATIC,chunk)
    do l=1,n,2
    !print*,n
       !write(ifu_log,*)'thread',omp_get_thread_num(),'l=',l
       x = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
       erp = 0.5d0*((Em-b)*x + (Em+b))
       erm = 0.5d0*(-(Em-b)*x + (Em+b))
       EE(1) = rrr*exp(ui*erp)-rrr+Eq
       EE(2) = rrr*exp(ui*erm)-rrr+Eq
       do k=1,2
          call gplus0_SOC(EE(k),greenn(k,:,:),sgn)    
       end do
        
       do i = 1,DNAOrbs
          do j = 1,DNAOrbs
             PDP(i,j) = PDP(i,j)+ a*(rrr*exp(ui*erp)*greenn(1,i,j)*der0+rrr*exp(ui*erm)*greenn(2,i,j)*der0)*xs(l)**4
          end do
       end do
    end do
!$OMP END DO
!$OMP CRITICAL
       do i = 1,DNAOrbs
          do j = 1,DNAOrbs
             PD_SOC(i,j)=PD_SOC(i,j)+PDP(i,j)                    
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL
    CH = 0.d0
    do i=1,DNAOrbs
       do j=1,DNAOrbs
          CH = CH + real(PD_SOC(i,j)*S_SOC(j,i))
       end do
    enddo
    !print*,'ch', 16*CH/(3*(n+1))
    ! ... replacing n by 2n+1
    n = n + n + 1
    ! Stopping?
     if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc.and.n.le.M) goto 1
    ! Test for successfullness and integral final value
    M = 0
     if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc) M = 1
    CH = 16*CH/(3*(n+1))
    do i=1,DNAOrbs
       do j=1,DNAOrbs
          PD_SOC(i,j) = sgn*16*PD_SOC(i,j)/(3*(n+1))
       enddo
    enddo
   !print*,'ch', CH
    if (sgn == 1) write(ifu_log,'(A,I4,A)')' Integration of the retarded density matrix on the complex plane has needed ',(((n-1)/2)+1)/2,' points'
    if (sgn == -1) write(ifu_log,'(A,I4,A)')' Integration of the advanced density matrix on the complex plane has needed ',(((n-1)/2)+1)/2,' points'

    return
  end subroutine IntCompPlane_SOC


  !*********************************************************************
  !* Compute retarded Green's function with SOC (doubling the matrices)*
  !*********************************************************************
  subroutine gplus_SOC(z,green,gammar,gammal,sgn)
    use parameters, only: eta, glue
    use constants, only: c_zero, ui
#ifdef PGI
    use lapack_blas, only: zgetri, zgetrf
#endif

    implicit none
    external zgetri, zgetrf     

    integer :: i, j, info, omp_get_thread_num
    integer, intent(in) :: sgn
    integer, dimension(DNAOrbs) :: ipiv
    complex*16, dimension(4*DNAOrbs) :: work
    complex*16, intent(in) :: z 
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl1,sigr1
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl2,sigr2
    complex*16, dimension(DNAOrbs,DNAOrbs), intent(out) :: green
    complex*16, dimension(DNAOrbs,DNAOrbs), intent(out) :: gammar,gammal
    complex*16, dimension(DNAOrbs,DNAOrbs) :: sigmar,sigmal
    
    ! Initilization 
    green=c_zero
    gammar=c_zero
    gammal=c_zero
    sigmar=c_zero
    sigmal=c_zero


    sigr1=-sgn*ui*eta*SD 
    sigl1=-sgn*ui*eta*SD 
    sigr2=-sgn*ui*eta*SD 
    sigl2=-sgn*ui*eta*SD 

    call CompSelfEnergies( 1, z, sigl1, sigr1, sgn )
   
    if (NSpin == 2) call CompSelfEnergies( 2, z, sigl2, sigr2, sgn )
    sigr1=glue*sigr1
    sigl1=glue*sigl1
    sigr2=glue*sigr2
    sigl2=glue*sigl2
    
    !************************************************************************
    !c Coupling matrices                                   
    !************************************************************************
    if (NSpin == 2) then
       do i=1,NAOrbs
       do j=1,NAOrbs
          sigmar(i,j)=sigr1(i,j)
          sigmar(i+NAOrbs,j+NAOrbs)=sigr2(i,j)
          sigmal(i,j)=sigl1(i,j)
          sigmal(i+NAOrbs,j+NAOrbs)=sigl2(i,j)
       end do
       end do
    else
       do i=1,NAOrbs
       do j=1,NAOrbs
          sigmar(i,j)=sigr1(i,j)
          sigmar(i+NAOrbs,j+NAOrbs)=sigr1(i,j)
          sigmal(i,j)=sigl1(i,j)
          sigmal(i+NAOrbs,j+NAOrbs)=sigl1(i,j)
       end do
       end do
    end if

    gammar=ui*(sigmar-conjg(transpose(sigmar)))
    gammal=ui*(sigmal-conjg(transpose(sigmal)))

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************

    do i=1,DNAOrbs
       do j=1,DNAOrbs
          green(i,j)=(z-shift)*S_SOC(i,j)-H_SOC(i,j)-sigmal(i,j)-sigmar(i,j)
       enddo
    enddo

    call zgetrf(DNAOrbs,DNAOrbs,green,DNAOrbs,ipiv,info)
    call zgetri(DNAOrbs,green,DNAOrbs,ipiv,work,4*DNAOrbs,info)

  end subroutine gplus_SOC

  !****************************************************
  !* Compute retarded Green's function with SOC*
  !****************************************************
  subroutine gplus0_SOC(z,green,sgn)
    use parameters, only: eta, glue
    use constants, only: c_zero, ui
#ifdef PGI
    use lapack_blas, only: zgetri, zgetrf
#endif

    implicit none
    external zgetri, zgetrf 
    
    integer :: i, j, info, omp_get_thread_num 
    integer, dimension(DNAOrbs) :: ipiv
    integer, intent(in) :: sgn
    complex*16, dimension(4*DNAOrbs) :: work
    complex*16, intent(in) :: z 
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl1,sigr1
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl2,sigr2
    complex*16, dimension(DNAOrbs,DNAOrbs), intent(out) :: green
    complex*16, dimension(DNAOrbs,DNAOrbs) :: sigmar,sigmal
    
    ! Initilization 
    green=c_zero
    sigmar=c_zero
    sigmal=c_zero
    sigr1=c_zero    
    sigl1=c_zero     
    sigr2=c_zero     
    sigl2=c_zero     
    sigr1=-sgn*ui*eta*SD
    sigl1=-sgn*ui*eta*SD 
    sigr2=-sgn*ui*eta*SD 
    sigl2=-sgn*ui*eta*SD

    call CompSelfEnergies( 1, z, sigl1, sigr1, sgn )
   
    if (NSpin == 2) call CompSelfEnergies( 2, z, sigl2, sigr2, sgn )
    sigr1=glue*sigr1
    sigl1=glue*sigl1
    sigr2=glue*sigr2
    sigl2=glue*sigl2
    
    !************************************************************************
    !c Coupling matrices                                   
    !************************************************************************
    if (NSpin == 2) then
       do i=1,NAOrbs
       do j=1,NAOrbs
          sigmar(i,j)=sigr1(i,j)
          sigmar(i+NAOrbs,j+NAOrbs)=sigr2(i,j)
          sigmal(i,j)=sigl1(i,j)
          sigmal(i+NAOrbs,j+NAOrbs)=sigl2(i,j)
       end do
       end do
    else
       do i=1,NAOrbs
       do j=1,NAOrbs
          sigmar(i,j)=sigr1(i,j)
          sigmar(i+NAOrbs,j+NAOrbs)=sigr1(i,j)
          sigmal(i,j)=sigl1(i,j)
          sigmal(i+NAOrbs,j+NAOrbs)=sigl1(i,j)
       end do
       end do
    end if

   !print*,'selfenergy right'
   !print*,sigmar(1,1),sigmar(1,2),sigmar(2,1),sigmar(1,3),sigmar(3,1)

   !print*,'selfenergy left'
   !print*,sigmal(1,1),sigmal(1,2),sigmal(2,1),sigmal(1,3),sigmar(3,1)

   !print*,'..........'
    !************************************************************************
    !c Retarded or advanced Green's functions
    !************************************************************************

    if (sgn == 1) then
    do i=1,DNAOrbs
       do j=1,DNAOrbs
          green(i,j)=(z-shift)*S_SOC(i,j)-H_SOC(i,j)-sigmal(i,j)-sigmar(i,j)
       enddo
    enddo
    else if (sgn == -1) then 
    do i=1,DNAOrbs
       do j=1,DNAOrbs
         !green(i,j)=(z-shift)*S_SOC(i,j)-H_SOC(i,j)-conjg(sigmal(i,j))-conjg(sigmar(i,j))
          green(i,j)=(z-shift)*S_SOC(i,j)-H_SOC(i,j)-sigmal(i,j)-sigmar(i,j)
       enddo
    enddo
    else
       print*,'Warning in gplus0_SOC'
       stop
    end if

    call zgetrf(DNAOrbs,DNAOrbs,green,DNAOrbs,ipiv,info)
    call zgetri(DNAOrbs,green,DNAOrbs,ipiv,work,4*DNAOrbs,info)

   !print*,'in gplus0'
   !print*,sgn
   !print*,green(1,1)
   !print*,green(1,3),green(3,1)
   !print*,'in gplus0'
    return
    
  end subroutine gplus0_SOC


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! SOC subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine spin_orbit

    use SpinOrbit, only: CompHSO
    use SpinRotate, only: CompHROT
    use cluster, only: LoAOrbNo, HiAOrbNo
    use parameters, only: PrtHatom, SOC, ROT
    use constants, only: c_zero, d_zero
#ifdef G03ROOT
    use g03Common, only: GetNShell
#endif
#ifdef G09ROOT
    use g09Common, only: GetNShell
#endif    
      
    complex*16, dimension(DNAOrbs,DNAOrbs) :: overlaprot, hamilrot, hamil_SO
    integer :: i,j,totdim,nshell,Atom
    real*8 :: uno
 
 Atom = PrtHatom   
 totdim=NAOrbs   
 hamilrot = c_zero
 overlaprot = c_zero
 hamil_SO = c_zero 
 H_SOC = c_zero
 S_SOC = d_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Duplicate the size of the Hamiltonian and Overlap matrix to include up and down
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if (SOC) then 
   if (NSpin == 2) then
      do i=1,NAOrbs
      do j=1,NAOrbs
         H_SOC(i,j)=dcmplx(HD(1,i,j),0.0d0)
         H_SOC(i+NAOrbs,j+NAOrbs)=dcmplx(HD(2,i,j),0.0d0)
         S_SOC(i,j)=SD(i,j)              
         S_SOC(i+NAOrbs,j+NAOrbs)=SD(i,j)         
      end do
      end do
   else 
      do i=1,NAOrbs
      do j=1,NAOrbs
         H_SOC(i,j)=dcmplx(HD(1,i,j),0.0d0)
         H_SOC(i+NAOrbs,j+NAOrbs)=dcmplx(HD(1,i,j),0.0d0)
         S_SOC(i,j)=SD(i,j)              
         S_SOC(i+NAOrbs,j+NAOrbs)=SD(i,j)         
      end do
      end do
   end if
 end if
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Return the input matrix with double size plus the soc interaction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 nshell = GetNShell()
 
 If (ROT) CALL CompHROT(HD,hamilrot,SD,overlaprot,NAOrbs,nshell)
 If (SOC) CALL CompHSO(hamil_SO,HD,NAOrbs,nshell)
 
!PRINT *, "Hamil matrix for atom ",Atom," : "
!PRINT *, "Up-Up" 
!do i=LoAOrbNo(Atom),HiAOrbNo(Atom)
!    PRINT '(1000(F11.5))',  ( REAL(hamil( i, j )), j=LoAOrbNo(Atom),HiAOrbNo(Atom) ) 
!end do  
!PRINT *, "Up-Down" 
!do i=LoAOrbNo(Atom),HiAOrbNo(Atom)
!    PRINT '(1000(F11.5))',  ( REAL(hamil( i, j )), j=totdim+LoAOrbNo(Atom),totdim+HiAOrbNo(Atom) ) 
!end do  
!PRINT *, "Down-Up"
!do i=totdim+LoAOrbNo(Atom),totdim+HiAOrbNo(Atom)                                                     
!    PRINT '(1000(F11.5))',  ( REAL(hamil( i, j )), j=LoAOrbNo(Atom),HiAOrbNo(Atom) ) 
!end do                                   
!PRINT *, "Down-Down"
!do i=totdim+LoAOrbNo(Atom),totdim+HiAOrbNo(Atom)                                                     
!    PRINT '(1000(F11.5))',  ( REAL(hamil( i, j )), j=totdim+LoAOrbNo(Atom),totdim+HiAOrbNo(Atom) ) 
!end do                                                                                             
!
!PRINT *, "Real part of Hamil_SO matrix for atom ",Atom," : "
!PRINT *, "Up-Up" 
!do i=LoAOrbNo(Atom),HiAOrbNo(Atom)
!   PRINT '(1000(F11.5))',  ( REAL(hamil_SO( i, j )), j=LoAOrbNo(Atom),HiAOrbNo(Atom) ) 
!end do  
!PRINT *, "Up-Down" 
!do i=LoAOrbNo(Atom),HiAOrbNo(Atom)
!    PRINT '(1000(F11.5))',  ( REAL(hamil_SO( i, j )), j=totdim+LoAOrbNo(Atom),totdim+HiAOrbNo(Atom) ) 
!end do  
!PRINT *, "Down-Up"
!do i=totdim+LoAOrbNo(Atom),totdim+HiAOrbNo(Atom)                                                   
!    PRINT '(1000(F11.5))',  ( REAL(hamil_SO( i, j )), j=LoAOrbNo(Atom),HiAOrbNo(Atom) ) 
!end do                                   
!PRINT *, "Down-Down"
!do i=totdim+LoAOrbNo(Atom),totdim+HiAOrbNo(Atom)                                                   
!    PRINT '(1000(F11.5))',  ( REAL(hamil_SO( i, j )), j=totdim+LoAOrbNo(Atom),totdim+HiAOrbNo(Atom) ) 
!end do                                                                                             
!

!PRINT *, "Imaginary part of Hamil_SO matrix for atom ",Atom," : "
!PRINT *, "Up-Up" 
!do i=LoAOrbNo(Atom),HiAOrbNo(Atom)
!   PRINT '(1000(F11.5))',  ( IMAG(hamil_SO( i, j )), j=LoAOrbNo(Atom),HiAOrbNo(Atom) ) 
!end do  
!PRINT *, "Up-Down" 
!do i=LoAOrbNo(Atom),HiAOrbNo(Atom)
!    PRINT '(1000(F11.5))',  ( IMAG(hamil_SO( i, j )), j=totdim+LoAOrbNo(Atom),totdim+HiAOrbNo(Atom) ) 
!end do  
!PRINT *, "Down-Up"
!do i=totdim+LoAOrbNo(Atom),totdim+HiAOrbNo(Atom)                                                   
!    PRINT '(1000(F11.5))',  ( IMAG(hamil_SO( i, j )), j=LoAOrbNo(Atom),HiAOrbNo(Atom) ) 
!end do                                   
!PRINT *, "Down-Down"
!do i=totdim+LoAOrbNo(Atom),totdim+HiAOrbNo(Atom)                                                   
!    PRINT '(1000(F11.5))',  ( IMAG(hamil_SO( i, j )), j=totdim+LoAOrbNo(Atom),totdim+HiAOrbNo(Atom) ) 
!end do                                                                                             

 if (ROT) then
    S_SOC=d_zero
    do i=1, totdim*2
       do j=1, totdim*2
          S_SOC(i,j)=REAL(overlaprot(i,j))         
          H_SOC(i,j)=hamilrot(i,j)
       end do
    end do 
 end if    

 if (SOC) then
   do i=1, totdim*2
      do j=1, totdim*2
         H_SOC(i,j)=H_SOC(i,j) + hamil_SO(i,j)
      end do
   end do
 end if

 return
 end subroutine spin_orbit

END MODULE device
