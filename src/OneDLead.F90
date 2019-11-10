!*********************************************************!
!*********************  ANT.G-2.4.0  *********************!
!*********************************************************!
!                                                         !
!   Copyright (c) by                                      !
!                                                         !
!   David Jacob                                           !
!      Theory Department                                  !
!      Max-Planck-Institute for Microstructure Physics    !
!      Halle, 06120 (GERMANY)                             !
!                                                         !
!*********************************************************!
  MODULE OneDLead
!*********************************************************!
!* Module for description of one-dimensional              !
!* lead parameters calculated from 1st principles         !
!* NOT FINISHED YET                                       !
!*********************************************************!
  IMPLICIT NONE
  SAVE
  
  ! *******************************
  ! One-dimensional lead parameters
  ! *******************************
  TYPE T1DLead
     PRIVATE

     ! *** Lead number
     INTEGER :: LeadNo
     ! *** Number of non-degenerate Spin bands
     INTEGER :: NSpin                       
     ! *** Number of atomic orbitals in supercell
     INTEGER :: NAOrbs
     ! *** Number of electrons in supercell
     INTEGER :: NElectrons

     ! *** Supercell Hamiltonian 
     COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE   :: H0    
     ! *** Coupling between neighbouring super-cells
     COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: V     
     ! *** Self-energy
     COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: Sigma1, Sigma2
     ! *** Overlap within supercell
     COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: S0
     ! *** Overlap between neighbouring supercells
     COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: S1

     ! *** lower and upper bound for non-zero DOS
     REAL*8    :: EMin, EMax

     ! *** Number of primitive unit cells
     INTEGER :: NPC
     ! *** Number of atomic orbitals in PRIMITIVE unit cell
     INTEGER :: NPCAO
     ! *** Number of electrons in PRIMITIVE unit cell
     INTEGER :: NPCEl

     ! *** coupling between n-th neighbour PUCs
     COMPLEX*16, DIMENSION(:,:,:,:), ALLOCATABLE :: vn
     ! *** overlap between n-th neighbour PUCs 
     COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: sn
  END TYPE T1DLead

  TYPE(T1DLead),DIMENSION(2) :: Lead1D

  ! Internal constants for calculation of selfenergies
  INTEGER, PARAMETER, PRIVATE :: MaxCycle = 5000
  REAL*8, PARAMETER, PRIVATE :: Mix = 0.2d0

  ! Internal variables
  REAL*8, PRIVATE :: ChargeOffSet
  INTEGER, PRIVATE :: WhichLead

  !**************************** 
  !* Public module procedures *
  !****************************
  PUBLIC :: L1D_NAOrbs
  PUBLIC :: L1D_NSpin
  PUBLIC :: L1D_NElectrons
  PUBLIC :: Init1DLead
  PUBLIC :: CompSelfEnergy1D

  !*****************************
  !* Private module procedures *
  !*****************************
  PRIVATE :: SolveDyson1D
  PRIVATE :: BulkSDOS
  PRIVATE :: CompGreensFunc
  PRIVATE :: ReadParameters
  PRIVATE :: TotCharge
  PRIVATE :: AdjustFermi
  PRIVATE :: DDOS

CONTAINS

  ! ***********************************
  ! Access functions to T1DLead members
  ! ***********************************
  !
  ! *** Number of atomic orbitals in supercell ***
  !
  INTEGER FUNCTION L1D_NAOrbs( L1D )
    TYPE(T1DLead),INTENT(IN) :: L1D
    L1D_NAOrbs = L1D%NAOrbs
  END FUNCTION L1D_NAOrbs
  !
  ! *** Number of non-degenerate spin-channels ***
  !
  INTEGER FUNCTION L1D_NSpin( L1D )
    TYPE(T1DLead),INTENT(IN) :: L1D
    L1D_NSpin = L1D%NSpin
  END FUNCTION L1D_NSpin
  !
  ! *** Number of electrons in supercell ***
  !
  INTEGER FUNCTION L1D_NElectrons( L1D )
    TYPE(T1DLead),INTENT(IN) :: L1D
    L1D_NElectrons = L1D%NElectrons
  END FUNCTION L1D_NElectrons
  !
  ! *** 
  !
  REAL*8 FUNCTION L1D_EMin( L1D )
    TYPE(T1DLead),INTENT(IN) :: L1D
    L1D_EMin = L1d%EMin
  END FUNCTION L1D_EMin
  !
  ! *** 
  !
  REAL*8 FUNCTION L1D_EMax( L1D )
    TYPE(T1DLead),INTENT(IN) :: L1D
    L1D_EMax = L1d%EMax
  END FUNCTION L1D_EMax


  !
  ! *** Initialize 1D Lead ***
  !
  SUBROUTINE Init1DLead( L1D, LeadNo )
    USE constants
    USE numeric
    USE parameters, ONLY: LEADDOS, ESTEP
    IMPLICIT NONE

    TYPE(T1DLead),INTENT(OUT) :: L1D
    INTEGER,INTENT(IN) :: LeadNo
    INTEGER :: NPC, NPCAO, ispin, n, m, AllocErr, NAO , NSpin
    COMPLEX*16 :: zenergy
    REAL*8 :: energy
    
    ! Read lead parameters from file
    CALL ReadParameters( L1D )

    PRINT *, "Initializing 1D Lead No.", LeadNo
 
    L1D%LeadNo     = LeadNo
    L1D%NAOrbs     = L1D%NPC * L1D%NPCAO
    L1D%NElectrons = L1D%NPC * L1D%NPCEl
    NPC   = L1D%NPC
    NPCAO = L1D%NPCAO
    NAO   = L1D%NAOrbs
    NSpin = L1D%NSpin

    ALLOCATE( &
         L1D%    H0( NSpin, NAO, NAO ),&
         L1D%     V( NSpin, NAO, NAO ), &
         L1D% Sigma1( NSpin, NAO, NAO ), &
         L1D% Sigma2( NSpin, NAO, NAO ), &
         L1D% S0( NAO, NAO ), & 
         L1D% S1( NAO, NAO ), &
         STAT = AllocErr )
    IF( AllocErr /= 0 ) THEN
       PRINT *, "Error: Program could not allocate memory for TElectrode%H, TElectrode%V, and TElectrode%Sigma"
       STOP
    END IF
    ! *** Hamiltonian of super-cell
    !     
    !       / v0    v1    v2    ... vN-1 \
    ! H0 = |  v'1   v0    v1    ... vN-2  |
    !      |  :     :     :         :     |
    !       \ v'N-1 v'N-2 v`N-3 ... v0   /
    !      
    ! *** Coupling between super-cells (from right to left)
    !
    !       / vN       0  \
    !      |  �  �         |
    ! V  = |  �    �       |
    !      |  �      �     |
    !       \ v1  렁�  vN /
    !
    ! where N = NPC = Number of primitve cells in super-cell, 
    ! vn hopping from n-th right neighbour primitve cell
    ! v'n = hermitian adjoint of v_n
    !
    L1D%H0 = c_zero
    L1D%V  = c_zero
    L1D%S0 = c_zero
    L1D%S1 = c_zero
    DO m=1, NPC 
       DO n=1, m
          DO ispin=1, NSpin
             ! *** super-cell Hamiltonian H0
             L1D%H0(ispin,((m-1)*NPCAO+1):(m*NPCAO),((n-1)*NPCAO+1):(n*NPCAO)) = CONJG(TRANSPOSE( L1D%vn( ispin, m-n,:,: ) ))
             IF( n /= m ) L1D%H0(ispin,((n-1)*NPCAO+1):(n*NPCAO),((m-1)*NPCAO+1):(m*NPCAO)) = L1D%vn( ispin, m-n,:,: ) 
             ! *** coupling between neighbouring supercells
             L1D%V(ispin, ((m-1)*NPCAO+1):(m*NPCAO),((n-1)*NPCAO+1):(n*NPCAO) ) = L1D%vn( ispin, NPC-m+n,:,:  )
          END DO
          ! *** super-cell overlap
          L1D%S0(((m-1)*NPCAO+1):(m*NPCAO),((n-1)*NPCAO+1):(n*NPCAO)) =  CONJG(TRANSPOSE( L1D%sn( m-n,:,: ) ))
          IF( n /= m ) L1D%S0( ((n-1)*NPCAO+1):(n*NPCAO),((m-1)*NPCAO+1):(m*NPCAO)) = L1D%sn( m-n,:,: )
          ! *** super-cell/super-cell overlap
          L1D%S1(((m-1)*NPCAO+1):(m*NPCAO),((n-1)*NPCAO+1):(n*NPCAO) ) = L1D%sn( NPC-m+n,:,: )
       END DO
    END DO
 
    ! Left lead: V / S1 = coupling/overlap from left to right
    IF( LeadNo == 1 )THEN
       DO ispin=1,NSpin
          L1D%V(ispin,:,:) = CONJG(TRANSPOSE( L1D%V(ispin,:,:) ))
       END DO
       L1D%S1 = CONJG(TRANSPOSE( L1D%S1 ))
    END IF

!!$    print *, "H0 = "
!!$    print *, L1D%H0(1,:,:)
!!$    print *, "V = "
!!$    print *, L1D%V(1,:,:)
!!$   
!!$    print *, "S0 = "
!!$    print *, L1D%S0(:,:)
!!$    print *, "S1 = "
!!$    print *, L1D%S1(:,:)
!!$   
    
    !!ChargeOffSet = 0.0d0
    !!WhichLead = LeadNo
    !!L1D%EMin = -2.0d0
    !!L1D%EMax =  6.0d0
    !print *, "Charge(-20.0d0,-5.0d0) = ",  TotCharge( -5.0d0  )

    ! Adjust parameters such that EFermi=0
    CALL AdjustFermi(L1D)

    ! Write electrode bulk DOS to file
    IF( LeadDOS )THEN
       IF(LeadNo==1)THEN 
          OPEN(UNIT=10,FILE="Lead1DOS.dat",STATUS="UNKNOWN")
       ELSE
          OPEN(UNIT=10,FILE="Lead2DOS.dat",STATUS="UNKNOWN")
       END IF
       DO energy=L1D%EMin,L1D%EMax,0.5*EStep
          zenergy=energy !!+1.0d-2*ui
          !!WRITE(UNIT=10,FMT=*), energy, BulkSDOS( L1D, 1, zenergy ), -DIMAG(Trace(L1D%Sigma1(1,:,:)))!!, -DIMAG(Trace(L1D%Sigma2(1,:,:)))
          !!, &
          !!     BulkSDOS( L1D, 2, zenergy ), -DIMAG(Trace(L1D%Sigma1(2,:,:))), -DIMAG(Trace(L1D%Sigma2(2,:,:)))
          WRITE(UNIT=10,FMT=*), energy, BulkSDOS( L1D, 1, zenergy ), -DIMAG(CTrace(L1D%Sigma1(1,:,:)))/d_pi
       END DO
       CLOSE(10)
    END IF
    STOP
  END SUBROUTINE Init1DLead


  !
  ! *** Compute self-energy matrices projected into device ***
  !
  SUBROUTINE CompSelfEnergy1D( L1D, spin, cenergy, Sigma_n )
#ifdef G03ROOT
    USE g03Common, ONLY: GetNBasis
#endif
#ifdef G09ROOT
    USE g09Common, ONLY: GetNBasis
#endif
    IMPLICIT NONE
    TYPE(T1DLead), intent(inout) :: L1D
    INTEGER, INTENT(in) :: spin
    COMPLEX*16, INTENT(in) :: cenergy
    COMPLEX*16, DIMENSION(:,:) :: Sigma_n

    INTEGER :: i,j,AOrb1,AOrb2, ispin

    ispin=spin
    IF( L1D%NSpin==1 .AND. spin == 2 )ispin=1

    ! Left lead (=1):  connected to orbitals 1 to NAOrbs in cluster
    ! Right lead (=2): connected to orbitals NBasis-NAOrbs+1 to NBasis 
    IF(L1D%LeadNo==1)THEN
       AOrb1=1
       AOrb2=L1D%NAOrbs
    ELSE
       AOrb1=GetNBasis()-L1D%NAOrbs+1
       AOrb2=GetNBasis()
    END IF
    CALL SolveDyson1D( Sigma_n(AOrb1:AOrb2,AOrb1:AOrb2), cenergy, L1D%S0, L1D%H0(ispin,:,:), L1D%V(ispin,:,:), L1D%S1 )
  END SUBROUTINE CompSelfEnergy1D


  !
  ! *** Self-consisten Dyson solver for 1D lead with overlap ***
  !
  !                                  1
  ! Sigma(E) = (V - E*S1) ------------------------ (V' - E*S1')   
  !                       ( E*S0 - H0 - Sigma(E) )
  !
  SUBROUTINE SolveDyson1D( Sigma, E, S0, H0, V, S1 )
    USE constants
    USE parameters
#ifdef PGI
     USE lapack_blas, ONLY: zgetri,zgetrf
#endif
    IMPLICIT NONE
    external zgetri,zgetrf

    COMPLEX*16, DIMENSION(:,:),INTENT(OUT) :: Sigma 
    COMPLEX*16, DIMENSION(:,:),INTENT(IN) :: V, H0, S0, S1
    COMPLEX*16, INTENT(IN) :: E
    
    INTEGER :: ispin, i, j, N, ncycle, ipiv(SIZE(Sigma,1)), info, AllocErr
    REAL*8 :: DOSS, OldDOSS
    
    ! auxiliary self energy for self-consistent calculation, tempory matrix
    COMPLEX*16, DIMENSION(SIZE(Sigma,1),SIZE(Sigma,1)) :: Sigma_aux, temp
    ! work array for inversion
    COMPLEX*16, DIMENSION( 4*SIZE(Sigma,1) ) :: work

    N = SIZE(Sigma,1)

    !! *** Initialize directional self-energies and energy matrix
    !!Sigma = c_zero
    !!DO i=1,N
    !!   Sigma(i,i) = -10.0d0*ui
    !!END DO
    
    OldDOSS = d_zero

    ncycle = 0
    ! *** Selfconsistency ***
    DO !!ncycle=1,MaxCycle
       ncycle = ncycle+1
       DO i=1,N
          DO j=1,N
             Sigma_aux(i,j) = E*S0(i,j) - H0(i,j) - Sigma(i,j)
          END DO
       END DO
       ! Inverting 
       CALL zgetrf(N,N,Sigma_aux,N,ipiv,info)
       CALL zgetri(N,Sigma_aux,N,ipiv,work,4*N,info)

       ! Now Sigma_aux contains the 
       ! surface Green's function 
       ! -> Compute actual Surface DOS
       !    *  DOSS = -1/pi Tr[GS]   *
       DOSS = d_zero
       DO i=1,N
          DOSS = DOSS + DIMAG( Sigma_aux(i,i) )
       END DO
       DOSS = -DOSS/d_pi
       ! Mix with old DOS
       IF( ncycle > 1 ) DOSS = (1.0d0-Mix)*DOSS + Mix*OldDOSS

       ! Matrix multiplication: (V - E*S1) (E*S0 - H0 - Sigma )^-1 (V' - E*S1')
       CALL zgemm('N','C',N,N,N, c_one,Sigma_aux,N,(V-E*S1), N, c_zero,temp,N )
       CALL zgemm('N','N',N,N,N, c_one,(V-E*S1),N,temp,N, c_zero, Sigma_aux,N )

       ! Mixing with old self-energy matrix
       Sigma = (1.0d0-Mix)*Sigma_aux + Mix*Sigma

       ! test convergence criterion: compare with previous DOS
       IF( ncycle > 1 .AND. (DOSS < 1.0d-10 .OR. ABS(DOSS-OldDOSS) < 1.0d-12*OldDOSS) ) EXIT
       OldDOSS = DOSS
    ENDDO ! End of self-consistency loop    
  END SUBROUTINE SolveDyson1D


  ! *************************************
  ! Routine to adjust Fermi-level to zero
  ! 1. Estimates upper/lower energy 
  !    boundary of lead DOS EMin/EMax,
  !    above/below which DOS is gauranteed 
  !    to be zero
  ! 2. Searches Fermi level
  ! 3. Shifts on-site Hamiltonian H0 to
  !    adjust Fermi level to zero
  ! *************************************
  SUBROUTINE AdjustFermi( L1D )
    USE parameters, ONLY: ChargeAcc, FermiAcc
    USE constants
    USE numeric

    TYPE(T1DLead), INTENT(inout) :: L1D

    INTEGER :: ispin, i, cond, max,k 
    REAL*8 :: EStep, EFermi, Q
    REAL*8 :: E0,E1,E2,E3,Z, Delta, Epsilon, EMin

    PRINT *, "Adjusting Fermi level to zero for Lead ", L1D%LeadNo

    ChargeOffSet = d_zero
    WhichLead = L1D%LeadNo
    
    PRINT *, "Searching boundaries [EMin, EMax]"
    PRINT *, "such that Int[EMin, EMax] DOS(E) =", 2*L1D%NAOrbs, " ( Number of spin orbitals )."

    EStep = 10.0d0
    L1D%EMin = -EStep
    L1D%EMax =  EStep
    
    DO
       Q = TotCharge( L1D%EMax )
       PRINT *, "EMin=", L1D%EMin, "  EMax=", L1D%EMax , "  Charge=", Q
       IF( ABS(Q - 2*L1D%NAOrbs) <  ChargeAcc ) THEN
          EXIT
       END IF
       L1D%EMin = L1D%EMin - EStep
       L1D%EMax = L1D%EMax + EStep
    END DO

    ! Fine tuning for EMin ...
    DO
       L1D%EMin = L1D%EMin + 0.2d0*EStep
       Q = TotCharge( L1D%EMax )
       PRINT *, "EMin=", L1D%EMin, "  EMax=", L1D%EMax , "  Charge=", Q
       IF( ABS(Q - 2*L1D%NAOrbs) >  ChargeAcc ) EXIT       
    END DO
    L1D%EMin = L1D%EMin - 0.2d0*EStep

    ! Fine tuning for EMax ...
    DO
       L1D%EMax = L1D%EMax - 0.2d0*EStep
       Q = TotCharge( L1D%EMax )
       PRINT *, "EMin=", L1D%EMin, "  EMax=", L1D%EMax , "  Charge=", Q
       IF( ABS(Q - 2*L1D%NAOrbs) >  ChargeAcc ) EXIT       
    END DO
    L1D%EMax = L1D%EMax + 0.2d0*EStep
    Q = TotCharge( L1D%EMax )
    PRINT *, "EMin=", L1D%EMin, "  EMax=", L1D%EMax , "  Charge=", Q

    PRINT *, "Searching Fermi energy..."

    E0 = 0.5*(L1D%EMin+L1D%EMax)
    E1 = E0 - 1.0d0
    E2 = E0 + 1.0d0
    Max = 25
    Delta=FermiAcc
    Epsilon=FermiAcc

    EMin = L1D%EMin
    L1D%EMin=2.0d0*EMin

    print *, L1D%EMin

    print*, "Q(Emin)=", TotCharge( EMin )
    print*, "Q(Emax)=", TotCharge( L1D%EMax )
    
    STOP

    ChargeOffset = L1D%NElectrons
    CALL MULLER(TotCharge,E0,E1,E2,Delta,Epsilon,Max,EFermi,Z,K,Cond)
    !!!EFermi = rtbis(TotCharge,EMin,L1D%EMax,FermiAcc) !!-0.218*Hart
    ChargeOffset = d_zero

    !!L1D%EMin=EMin

    PRINT *, "EFermi = ", EFermi

    ! Shift on-site Hamiltonian H0 so that EFermi = 0
    DO ispin=1,L1D%NSpin
       DO i=1,L1D%NAOrbs
          L1D%H0(ispin,i,i)=L1D%H0(ispin,i,i)-EFermi
       END DO
    END DO

    ! Also shift EMin, EMax
    L1D%EMin = L1D%EMin - EFermi
    L1D%EMax = L1D%EMax - EFermi       
  END SUBROUTINE AdjustFermi


  ! ***********************************
  ! Total charge up to energy E
  ! up to Energy E, lower bound is EMin
  ! ***********************************  
  REAL*8 FUNCTION TotCharge( E )
    USE constants
    USE parameters, ONLY: ChargeAcc
    USE numeric
    IMPLICIT NONE
    
    REAL*8, INTENT(in) :: E

    INTEGER, PARAMETER :: nmax = 2047
    
    REAL*8 :: q,qq, E0, R
    INTEGER :: n, n1, n2, j, i,ispin
    REAL*8, DIMENSION((nmax*(nmax+1))/2),SAVE :: x, w
    LOGICAL, SAVE :: InitQGauss = .TRUE.
    
    ! Computing all abscissas and weights
    IF( InitQGauss ) THEN
       DO n=1,nmax
          n1=(n*(n-1))/2+1
          n2=(n*(n+1))/2
          CALL gauleg(d_zero,d_pi,x(n1:n2),w(n1:n2),n)
       END DO
       InitQGauss = .FALSE.
    END IF

    !! *** Initialize directional self-energies and energy matrix
    Lead1D(WhichLead)%Sigma1 = c_zero
    Lead1D(WhichLead)%Sigma2 = c_zero
    DO i=1,Lead1D(WhichLead)%NAOrbs
       DO ispin=1,Lead1D(WhichLead)%NSpin
          Lead1D(WhichLead)%Sigma1(ispin,i,i) = -10.0d0*ui
          Lead1D(WhichLead)%Sigma2(ispin,i,i) = -10.0d0*ui
       END DO
    END DO
    !Integration contour parameters:

    E0 = 0.5*(E + Lead1D(WhichLead)%EMin)
    R  = 0.5*(E - Lead1D(WhichLead)%EMin)

    ! Computing integral of DOS over 
    ! complex contour using Gauss-Legendre 
    ! quadrature
    n=1
    DO 
       n1=(n*(n-1))/2+1
       n2=(n*(n+1))/2 
       q = d_zero
       DO j=n1,n2
          !!print*, "phi=", x(j)
          q = q + w(j)*ddos( Lead1D(WhichLead), E0, R, x(j) )
       END DO

       !!print *, "n =", n, "  q =", q
       
       IF( n > 1 .AND. (q == d_zero .OR. ABS(q-qq) < ChargeAcc ) ) EXIT  
       n=2*n+1
       IF( n > nmax )THEN
          PRINT *, "TotCharge/gaussian quadrature has not converged after", nmax, " steps."
          STOP
       END IF
       qq = q
    END DO

    TotCharge = q-ChargeOffSet
  END FUNCTION TotCharge


  !
  ! *** Integrand for charge integration along complex contour ***
  !
  REAL*8 FUNCTION DDOS( L1D, E0, R, phi )
    USE constants
    USE numeric
    !USE parameters

    TYPE(T1DLead), INTENT(INOUT) :: L1D
    REAL*8, INTENT(in) :: phi, E0, R
    
    INTEGER :: i, j, ispin 
    COMPLEX*16,DIMENSION(L1D%NAOrbs,L1D%NAOrbs) :: green,S0G0, gL, gR, S1GL0, S1GR0, VMES1
    COMPLEX*16 :: TrG, z

    z = E0 - R*(COS(phi) - ui*SIN(phi)) 

    TrG=c_zero
    DO ispin=1,L1D%NSpin
       CALL CompGreensFunc( L1D, ispin, z, green )
       S0G0 = MATMUL( L1D%S0, green )

!!$       gL = (z+ui*eta)*L1D%S0 - L1D%H0(ispin,:,:) - L1D%Sigma1(ispin,:,:)
!!$       gR = (z+ui*eta)*L1D%S0 - L1D%H0(ispin,:,:) - L1D%Sigma2(ispin,:,:)
!!$       
!!$       CALL CInv( gL )
!!$       CALL CInv( gR )
!!$       
!!$       VMES1 = L1D%V(ispin,:,:)-z*L1D%S1
!!$       
!!$       S1GL0 = MATMUL( CONJG(TRANSPOSE(L1D%S1)), gL )
!!$       S1GL0 = MATMUL( S1GL0, VMES1 )
!!$       S1GL0 = MATMUL( S1GL0, green )
!!$       
!!$       S1GR0 = MATMUL( L1D%S1, gL )
!!$       S1GR0 = MATMUL( S1GL0, CONJG(TRANSPOSE(VMES1)) )
!!$       S1GR0 = MATMUL( S1GL0, green )

       DO i=1,L1D%NAOrbs
          TrG = TrG + S0G0(i,i) !!+ S1GL0(i,i) +  S1GR0(i,i)
          !!DO j=1,L1D%NAOrbs
          !!   TrG = TrG + L1D%S0(i,j)*green(j,i)
          !!END DO
       END DO
    END DO
    IF(L1D%NSpin==1) TrG = TrG * 2.0d0

    DDOS = -DIMAG(R*(SIN(phi)+ui*COS(phi))*TrG)/d_pi 

  END FUNCTION DDOS


  ! 
  ! *** Spin-resolved Bulk DOS ***
  !
  REAL*8 FUNCTION BulkSDOS( L1D, spin, energy )
    USE constants
    USE numeric
    USE parameters, ONLY: ETA

    TYPE(T1DLead), intent(inout) :: L1D
    INTEGER, INTENT(in) :: spin
    COMPLEX*16, INTENT(in) :: energy

    COMPLEX*16, DIMENSION( L1D%NAOrbs, L1D%NAOrbs ) :: G0, S0G0, gL, gR, S1GL0, S1GR0, VMES1
    INTEGER :: i,ispin,info

    ! Initialize self-energies
    L1D%Sigma1 = c_zero
    L1D%Sigma2 = c_zero
    DO i=1,L1D%NAOrbs
       DO ispin=1,L1D%NSpin
          L1D%Sigma1(ispin,i,i) = -10.0d0*ui
          L1D%Sigma2(ispin,i,i) = -10.0d0*ui
       END DO
    END DO

    CALL CompGreensFunc( L1D, Spin, energy, G0 )

    ispin = Spin

    gL = (energy+ui*eta)*L1D%S0 - L1D%H0(ispin,:,:) - L1D%Sigma1(ispin,:,:)
    gR = (energy+ui*eta)*L1D%S0 - L1D%H0(ispin,:,:) - L1D%Sigma2(ispin,:,:)

    info = CInv( gL )
    info = CInv( gR )

    VMES1 = L1D%V(ispin,:,:)-energy*L1D%S1

    S0G0 = MATMUL( L1D%S0, G0 )

    S1GL0 = MATMUL( CONJG(TRANSPOSE(L1D%S1)), gL )
    S1GL0 = MATMUL( S1GL0, VMES1 )
    S1GL0 = MATMUL( S1GL0, G0 )

    S1GR0 = MATMUL( L1D%S1, gL )
    S1GR0 = MATMUL( S1GL0, CONJG(TRANSPOSE(VMES1)) )
    S1GR0 = MATMUL( S1GL0, G0 )

    !!S1GL0 = MATMUL( MATMUL( CONJG(TRANSPOSE(L1D%S1)), MATMUL( gL, L1D%V(ispin,:,:)-energy*L1D%S1 )), G0 )
    !!S1GR0 = MATMUL( MATMUL( L1D%S1, MATMUL( gR, CONJG(TRANSPOSE(L1D%V(ispin,:,:)-energy*L1D%S1)) )), G0 )

    BulkSDOS = d_zero
    DO i=1,L1D%NAOrbs
       BulkSDOS = BulkSDOS + DIMAG(S0G0(i,i)+S1GL0(i,i)+S1GR0(i,i))
    END DO
    BulkSDOS = -BulkSDOS/d_pi
  END FUNCTION BulkSDOS


  ! 
  ! *** Compute Bulk Green's function G0 *** 
  ! 
  SUBROUTINE CompGreensFunc( L1D, Spin, energy, G0 )
    USE constants
    USE parameters
#ifdef PGI
    USE lapack_blas, ONLY: zgetri,zgetrf
#endif    
    external zgetri,zgetrf
    !USE lapack_blas, ONLY: zgetri,zgetrf

    TYPE(T1DLead), INTENT(inout) :: L1D
    INTEGER, INTENT(in) :: Spin
    COMPLEX*16, INTENT(in) :: energy
    COMPLEX*16, DIMENSION( L1D%NAOrbs, L1D%NAOrbs ),INTENT(out) :: G0

    COMPLEX*16, DIMENSION( L1D%NAOrbs, L1D%NAOrbs ) :: E
    INTEGER :: n, k, ipiv(L1D%NAOrbs), info, i, j, ispin
    COMPLEX*16, DIMENSION( 4*L1D%NAOrbs  ) :: work

    ispin = spin
    IF( spin == 2 .AND. L1D%NSpin == 1 ) ispin = 1

    IF( ispin < 1 .OR. ispin > 2 ) THEN
       PRINT *, "Undefined value for Spin: ", ispin
    END IF

    N=L1D%NAOrbs

    CALL SolveDyson1D( L1D%Sigma1(ispin,:,:), energy, L1D%S0, & 
         L1D%H0(ispin,:,:),  L1D%V(ispin,:,:), L1D%S1(:,:) )
    CALL SolveDyson1D( L1D%Sigma2(ispin,:,:), energy, L1D%S0, & 
        L1D%H0(ispin,:,:),  CONJG(TRANSPOSE(L1D%V(ispin,:,:))), CONJG(TRANSPOSE(L1D%S1(:,:))) )

    E = (energy+ui*eta)*L1D%S0

    G0 = E - L1D%H0(ispin,:,:)-L1D%Sigma1(ispin,:,:)-L1D%Sigma2(ispin,:,:)

    CALL zgetrf(n,n,G0,n,ipiv,info)
    CALL zgetri(n,G0,n,ipiv,work,4*n,info)
  END SUBROUTINE CompGreensFunc


  !
  ! *** Read lead parameters from file ***
  ! 
  SUBROUTINE ReadParameters( L1D )
    USE constants
    IMPLICIT NONE
    TYPE(T1DLead),INTENT(INOUT) :: L1D

!!$      ! Experimental version: Al with ECP Clarkson
!!$      ODL%NSpin      = 1
!!$      ODL%NPUC       = 3
!!$      ODL%NAOrbs     = 4
!!$      ODL%NElectrons = 1

    ! Experimental version: TB chain
    L1D%NSpin = 1
    L1D%NPC   = 1
    L1D%NPCAO = 1
    L1D%NPCEl = 1
    
    ALLOCATE( L1D%vn(L1D%NSpin,0:L1D%NPC,L1D%NPCAO,L1D%NPCAO),  L1D%sn(0:L1D%NPC,L1D%NPCAO,L1D%NPCAO) );

    L1D%vn = c_zero
    L1D%sn = c_zero
    
    ! v0 = 0
    L1D%vn(1,0,1,1) = 0.0d0
    ! v1 = -1
    L1D%vn(1,1,1,1) = -1.0d0
    
    ! s0 = 1
    L1D%sn(0,1,1) = 1.0d0
    ! s1 = 0
    L1D%sn(1,1,1) = 0.4d0
    
!!$    ! *** Initializing coupling matrices
!!$    
!!$         ODL(LeadNo)%V(1,0,1,1) = -3.2187E-01*Hart
!!$         ODL(LeadNo)%V(1,1,1,1) = -2.8557E-01*Hart
!!$         ODL(LeadNo)%V(1,2,1,1) = -6.7472E-02*Hart
    
!!$    ODL%V(1,0,1,1) = 0.0
!!$    ODL%V(1,1,1,1) = -1.0   
!!$         ODL%V(1,2,1,1) = -0.1   
!!$         ODL%V(1,3,1,1) = -0.01
!!$
!!$    ODL%S(0,1,1) = 1.0d0
!!$    ODL%S(1,1,1) = 0.3
!!$         ODL%S(2,1,1) = 0.03
!!$         ODL%S(3,1,1) = 0.006
    
    !!         ODL%S(2,1,1) = 0.05
    !!         ODL%S(3,1,1) = 0.01
    
    !!         ODL(LeadNo)%V(1,1,1,1) = -1.0d0
    !!         ODL(LeadNo)%V(1,2,1,1) = -0.5d0
    
    ! *** Initializing Al parameters
    
!!$         ODL%V(1,0,1,1) = -5.6825E-01*Hart
!!$         ODL%V(1,0,2,1) = -1.4102E-16*Hart
!!$         ODL%V(1,0,2,2) = -3.0444E-01*Hart
!!$         ODL%V(1,0,3,3) = -2.0515E-01*Hart
!!$         ODL%V(1,0,4,3) =  1.0424E-19*Hart
!!$         ODL%V(1,0,4,4) = -2.0515E-01*Hart
!!$
!!$         DO i=1,ODL%NAOrbs
!!$            DO j=i+1,ODL%NAOrbs
!!$               ODL%V(1,0,i,j) = ODL%V(1,0,j,i) 
!!$            END DO
!!$         END DO
!!$
!!$         ODL%V(1,1,1,1) = -3.5940E-01*Hart
!!$         ODL%V(1,1,2,1) = -3.2842E-01*Hart
!!$         ODL%V(1,1,1,2) =  3.2842E-01*Hart
!!$         ODL%V(1,1,2,2) =  5.1138E-03*Hart
!!$         ODL%V(1,1,3,3) = -1.8051E-01*Hart
!!$         ODL%V(1,1,4,4) = -1.8051E-01*Hart
!!$
!!$         ODL%V(1,2,1,1) = -5.5767E-02*Hart
!!$         ODL%V(1,2,2,1) = -1.1254E-01*Hart
!!$         ODL%V(1,2,1,2) =  1.1254E-01*Hart
!!$         ODL%V(1,2,2,2) =  1.1129E-01*Hart
!!$         ODL%V(1,2,3,3) = -7.0537E-02*Hart
!!$         ODL%V(1,2,4,4) = -7.0537E-02*Hart
!!$
!!$         ODL%V(1,3,1,1) = -7.3145E-03*Hart
!!$         ODL%V(1,3,2,1) = -1.7339E-02*Hart
!!$         ODL%V(1,3,1,2) =  1.7339E-02*Hart
!!$         ODL%V(1,3,2,2) =  3.1674E-02*Hart
!!$         ODL%V(1,3,3,3) = -1.3979E-02*Hart
!!$         ODL%V(1,3,4,4) = -1.3979E-02*Hart

!!$         ODL%V(1,4,1,1) =  4.6893E-05*Hart
!!$         ODL%V(1,4,2,1) = -2.4757E-03*Hart
!!$         ODL%V(1,4,1,2) =  2.4757E-03*Hart
!!$         ODL%V(1,4,2,2) =  1.1153E-02*Hart
!!$         ODL%V(1,4,3,3) =  1.5110E-03*Hart
!!$         ODL%V(1,4,4,4) =  1.5110E-03*Hart
!!$
!!$         ODL%V(1,5,1,1) = -8.2112E-04*Hart
!!$         ODL%V(1,5,2,1) =  1.4406E-03*Hart
!!$         ODL%V(1,5,1,2) = -1.4406E-03*Hart
!!$         ODL%V(1,5,2,2) = -1.2502E-03*Hart
!!$         ODL%V(1,5,3,3) = -1.4941E-03*Hart
!!$         ODL%V(1,5,4,4) = -1.4941E-03*Hart
         
         ! *** Initializing overlap matrices
         
!!$         ODL%S(0,1,1) = 1.0000E+00
!!$         ODL%S(0,2,2) = 1.0000E+00
!!$         ODL%S(0,3,3) = 1.0000E+00
!!$         ODL%S(0,4,4) = 1.0000E+00
!!$         
!!$         ODL%S(1,1,1) =  5.0386E-01
!!$         ODL%S(1,2,1) =  5.8185E-01
!!$         ODL%S(1,1,2) = -5.8185E-01
!!$         ODL%S(1,2,2) = -1.1063E-01
!!$         ODL%S(1,3,3) = 5.1692E-01
!!$         ODL%S(1,4,4) = 5.1692E-01
!!$
!!$         ODL%S(2,1,1) = 8.4678E-02
!!$         ODL%S(2,2,1) = 1.9927E-01
!!$         ODL%S(2,1,2) = -1.9927E-01
!!$         ODL%S(2,2,2) = -3.0313E-01
!!$         ODL%S(2,3,3) = 1.0640E-01
!!$         ODL%S(2,4,4) = 1.0640E-01
!!$
!!$         ODL%S(3,1,1) = 7.5509E-03
!!$         ODL%S(3,2,1) = 2.6183E-02
!!$         ODL%S(3,1,2) = -2.6183E-02
!!$         ODL%S(3,2,2) = -7.5186E-02
!!$         ODL%S(3,3,3) = 1.0999E-02
!!$         ODL%S(3,4,4) = 1.0999E-02

!!$         ODL%S(4,1,1) = 3.7899E-04
!!$         ODL%S(4,2,1) = 1.6695E-03
!!$         ODL%S(4,1,2) =-1.6695E-03
!!$         ODL%S(4,2,2) =-6.9958E-03
!!$         ODL%S(4,3,3) = 5.6944E-04
!!$         ODL%S(4,4,4) = 5.6944E-04
!!$
!!$        ODL%S(5,1,1) = 1.0669E-05
!!$        ODL%S(5,2,1) = 5.4107E-05
!!$        ODL%S(5,1,2) =-5.4107E-05
!!$        ODL%S(5,2,2) =-2.7284E-04
!!$        ODL%S(5,3,3) = 1.3990E-05
!!$        ODL%S(5,4,4) = 1.3990E-05

  END SUBROUTINE ReadParameters

!!$    !
!!$    ! Routina para orthogonalizar el puto overlap
!!$    ! COJONES!
!!$    ! 
!!$    SUBROUTINE OrthogonalizeOverlap( LeadNo )
!!$      USE constants
!!$      USE numeric 
!!$      INTEGER*4, INTENT(IN) :: LeadNo
!!$
!!$      ! Overlap, Fock matrix and diagonalizing matrix in reciprocal space
!!$      COMPLEX*16, DIMENSION(:,:,:), allocatable :: Sk, Fk,SPrime, FPrime
!!$      COMPLEX*16, DIMENSION(:,:), allocatable :: Smhk
!!$      COMPLEX*16, DIMENSION(:,:,:,:), allocatable :: Smh
!!$
!!$      INTEGER*4 :: ik, Nk, NAO, m,n,mu,nu,m1,n1,mu1,nu1, NPUC, NSpin, ispin
!!$      REAL*8 :: k
!!$
!!$      print *, "Ortogonalizando el puto overlap de los cojones!"
!!$
!!$      Nk = 100
!!$      NPUC = ODL_NPUC( LeadNo )
!!$      NAO = ODL_NAOrbs( LeadNo )
!!$      NSpin = ODL_NSpin( LeadNo )
!!$
!!$      ALLOCATE( &
!!$           Sk(0:Nk-1,NAO,NAO), &
!!$           Fk(0:Nk-1,NAO,NAO), &
!!$           Smhk(NAO,NAO), &
!!$           Smh(-NPUC:NPUC,NAO,-NPUC:NPUC,NAO), &
!!$           SPrime(0:NPUC,NAO,NAO), &
!!$           FPrime(0:NPUC,NAO,NAO) )
!!$      
!!$      Sk = c_zero
!!$      Fk = c_zero
!!$      
!!$      DO ispin=1,NSpin
!!$         Smh  = c_zero
!!$         SPrime = c_zero
!!$         FPrime = c_zero
!!$
!!$         DO ik=0,Nk-1
!!$            Smhk = c_zero
!!$            ! compute matrices in reciprocal space
!!$            k = ((2*d_pi)/(Nk-1))*ik-d_pi
!!$
!!$            Fk(ik,:,:)=ODL(LeadNo)%V(ispin,0,:,:)
!!$            Sk(ik,:,:)=ODL(LeadNo)%S(0,:,:)
!!$
!!$            DO n=1,NPUC
!!$               Fk(ik,:,:)=Fk(ik,:,:) +EXP(-ui*k*n)*ODL(LeadNo)%V(ispin,n,:,:) &
!!$                    + EXP(ui*k*n)*TRANSPOSE(CONJG(ODL(LeadNo)%V(ispin,n,:,:)))
!!$               Sk(ik,:,:)=Sk(ik,:,:) +EXP(-ui*k*n)*ODL(LeadNo)%S(n,:,:) &
!!$                    + EXP(ui*k*n)*TRANSPOSE(CONJG(ODL(LeadNo)%S(n,:,:)))
!!$            END DO
!!$
!!$            ! Compute S(k)^-1/2
!!$            !!CALL CMatPow( Sk(ik,:,:), -0.5d0, Smhk(:,:) )
!!$            Smhk=Sk(ik,:,:)
!!$            
!!$            print *,  "ik = ", ik, "  k=", k
!!$            print *, " S(k)^-1/2 S(k) S(k)^-1/2 = "
!!$            print *,  MATMUL( Smhk, MATMUL( Sk(ik,:,:), Smhk ) )
!!$
!!$            ! S^-1/2 = \sum_k U_B(k) * S(k)^-1/2 * U_B(k)^*
!!$            DO m=-NPUC,NPUC
!!$               DO n=-NPUC,NPUC
!!$                  DO mu=1,NAO
!!$                     DO nu=1,NAO
!!$                        IF( ik == 0 .OR. ik == Nk-1 )THEN 
!!$                           Smh(m,mu,n,nu) = Smh(m,mu,n,nu) + (0.5d0*EXP(ui*k*m) * Smhk(mu,nu) * EXP(-ui*k*n))/(Nk-1)
!!$                        ELSE
!!$                        !!IF( IABS( m-n ) <= NPUC ) THEN
!!$                        Smh(m,mu,n,nu) = Smh(m,mu,n,nu) + (EXP(ui*k*(m-n)) * Smhk(mu,nu))/(Nk-1)
!!$                        !!END IF
!!$                        END IF
!!$                     END DO
!!$                  END DO
!!$               END DO
!!$            END DO
!!$           
!!$         END DO ! End Loop over k-vectors
!!$         
!!$         print*, "Smh(m,n)="
!!$         DO m=-NPUC,NPUC
!!$            DO n=-NPUC,NPUC
!!$               DO mu=1,NAO
!!$                  DO nu=1,NAO
!!$                     print *, m, n, dreal(Smh(m,mu,n,nu) )
!!$                  END DO
!!$               END DO
!!$            END DO
!!$         END DO
!!$            
!!$
!!$         ! Orthogonalization:
!!$         ! S' = S^-1/2 S S^-1/2
!!$         ! F' = S^-1/2 F S^-1/2
!!$         DO n=0,NPUC           
!!$            print *, "n =", n, "  SPrime(mu,nu)="
!!$            DO mu=1,NAO
!!$               DO nu=1,NAO
!!$                  DO m1=-NPUC,NPUC
!!$                     DO mu1=1,NAO
!!$                        DO n1=-NPUC,NPUC
!!$                           DO nu1=1,NAO
!!$                              IF( IABS(n1-m1) <= NPUC )THEN
!!$                                 IF( n1 >= m1 )THEN
!!$                                    SPrime(n,mu,nu) = SPrime(n,mu,nu) &
!!$                                         + Smh(0,mu,m1,mu1) &
!!$                                         * CONJG(ODL(LeadNo)%S(n1-m1,nu1,mu1)) &
!!$                                         * CONJG(Smh(n,nu,n1,nu1))
!!$                                 ELSE 
!!$                                    SPrime(n,mu,nu) = SPrime(n,mu,nu) &
!!$                                         + Smh(0,mu,m1,mu1) &
!!$                                         * ODL(LeadNo)%S(m1-n1,mu1,nu1) &
!!$                                         * CONJG(Smh(n,nu,n1,nu1))
!!$                                 END IF
!!$                              END IF
!!$                           END DO
!!$                        END DO
!!$                     END DO
!!$                  END DO
!!$                  print *, mu, nu, SPrime(n,mu,nu)
!!$               END DO
!!$            END DO
!!$
!!$         END DO
!!$
!!$      END DO ! End loop over spin
!!$
!!$    END SUBROUTINE OrthogonalizeOverlap

END MODULE OneDLead
