!*********************************************************!
!*********************  ANT.G-2.5.0  *********************!
!*********************************************************!
!                                                         !
!   Copyright (c) by                                      !
!                                                         !
!   David Jacob                                           !
!                                                         !
!      Theory Department                                  !
!      Max-Planck-Institute for Microstructure Physics    !
!      Halle, 06120 (GERMANY)                             !
!                                                         !
!*********************************************************!
  MODULE Correlation
!*********************************************************!
  use parameters, only: NCorrBl, CorrBeg, CorrEnd, UCoul, JHund
  use util
  implicit none

  ! Dimension of device Hilbert space and number of spin channels
  integer :: NDAO, NDSpin
  ! Indexes of correlated orbitals
  logical, dimension(:), allocatable :: cix
  ! Hamiltonian within each correlated subspace
  real*8, dimension(:,:,:,:),allocatable :: h_d
  ! Overlap within each correlated subspace (usually unity)
  real*8, dimension(:,:,:),allocatable :: s_d
  ! Dimensions of individual correlated subspaces
  integer, dimension(:),allocatable :: ncorrao
  ! Dimension of biggest correlated subspace
  integer :: NMaxCorr

  ! Inverse GF of correlated subspace
  complex*16, dimension(:,:), allocatable :: inv_g_d

contains
  
  subroutine InitCorrelation( NBasis, NSpin ) 
    implicit none
    integer, intent(in) :: nbasis, nspin
    ! Device Hamiltonian and overlap matrix

    integer :: iblock, iao, jao, i, j, ispin

    print *, "--------------------------------------------"
    print *, "--- Initialization of Correlation module ---"
    print *, "--------------------------------------------"

    if( NCorrBl < 1 )then
       print *, "No correlated blocks defined." 
       return
    end if

    NDAO=NBasis
    NDSpin=NSpin

    allocate( cix(nbasis), ncorrao(NCorrBl) )

    cix(:) = .false.
    NMaxCorr=0
    do iblock=1,NCorrBl
       print *, "iblock =", iblock
       print *, "  U=", UCoul(iblock), "  J=", JHund(iblock)
       ncorrao(iblock) = CorrEnd(iblock)-CorrBeg(iblock)+1
       if( ncorrao(iblock) > NMaxCorr ) NMaxCorr = ncorrao(iblock)
       do iao=CorrBeg(iblock),CorrEnd(iblock)
          print *, "  AO", iao
          cix(iao) = .true.
       end do
    end do
    allocate( h_d(NDSpin,NCorrBl,NMaxCorr,NMaxCorr), s_d(NCorrBl,NMaxCorr,NMaxCorr), inv_g_d(NMaxCorr,NMaxCorr) )
  end subroutine InitCorrelation


  subroutine SetHamOvl( HD, SD )
    use parameters, only: CorrBeg, CorrEnd
    use util
    implicit none
    real*8, dimension(NDSpin,NDAO,NDAO), intent(in) :: HD
    real*8, dimension(NDAO,NDAO), intent(in) :: SD
    integer :: iblock, i, j, iao, jao, ispin

    do iblock = 1,NCorrBl
       print '(A,I2,A,I4,A,I4)', "CorrBl #", iblock, ": AO #", CorrBeg(iblock),"- AO #", CorrEnd(iblock)
       i=0
       do iao=CorrBeg(iblock),CorrEnd(iblock)
          i=i+1
          j=0
          do jao=CorrBeg(iblock),CorrEnd(iblock)
             j=j+1
             do ispin=1,NDSpin
                h_d(ispin,iblock,i,j) = HD(ispin,iao,jao)
             end do
             s_d(iblock,i,j) = SD(iao,jao)
          end do
       end do
       do ispin=1,NDSpin
          if( NDSpin == 1 ) print *, "h_d ="
          if( NDSpin == 2 .and. ispin == 1 ) print *, "Spin-up h_d ="
          if( NDSpin == 2 .and. ispin == 2 ) print *, "Spin-down h_d ="
          call PrintRMatrix( h_d(ispin,iblock,:,:) )
       end do
       print *, "s_d ="
       call PrintRMatrix( s_d(iblock,:,:) )
       print *, "----------------------------------------------------------------------------"
    end do
  end subroutine SetHamOvl

  !
  ! Compute the hybridization functions Delta 
  ! for energy <omega> and spin <ispin>
  !
  subroutine CompDelta( ispin, omega, mu, GF, delta_d )
    use parameters
    use constants
    use numeric
    implicit none
    !
    ! INPUT parameters:
    !
    ! Spin
    integer, intent(in) :: ispin
    ! Energy and chemical potential
    real*8, intent(in) :: omega, mu
    ! Device Green's function
    complex*16, dimension(NDAO,NDAO), intent(in) :: GF
    ! Hybridization function for each block
    complex*16, dimension(NCorrBl,NMaxCorr,NMaxCorr), intent(out) :: delta_d

    integer :: iblock, iao, jao, i, j, iao_beg, iao_end, info

    delta_d = d_zero
    !
    ! Loop over correlated blocks
    !
    do iblock=1,NCorrBl
       iao_beg=CorrBeg(iblock)
       iao_end=CorrEnd(iblock)
       call CSetId(inv_g_d)
       !
       ! Loop over orbitals in correlated block to obtain local GF g_d
       !
       i = 0
       do iao=iao_beg,iao_end
          i = i + 1
          j = 0
          do jao=iao_beg,iao_end
             j = j + 1
             inv_g_d(i,j) = GF(iao,jao)
          end do
       end do
       !
       ! Compute inverse g0
       !
       info = CInv(inv_g_d) 
       if( info /=0 )then
          print '(A,I2,A,E20.10)', "Correlation/CompDelta/Warning: Could not invert Green's function for correlated subspace #", iblock, " at omega =", omega
       end if
       !
       ! Compute hybridization function for correlated block
       !
       do i = 1,ncorrao(iblock)
          do j = 1,ncorrao(iblock)
             delta_d(iblock,i,j) = (omega+mu)*s_d(iblock,i,j)-h_d(ispin,iblock,i,j)-inv_g_d(i,j)
          end do
       end do
    end do
    
  end subroutine CompDelta
    
  !
  ! Computes DFT+U potential for density matrix P 
  ! and adds it to Hamiltonian H
  !
  subroutine Add_DFT_plus_U_Pot( P, H )
    use constants
    implicit none
    ! Density matrix 
    real*8, dimension(NDSpin,NDAO,NDAO), intent(in) :: P
    ! Hamiltonian
    real*8, dimension(NDSpin,NDAO,NDAO), intent(inout) :: H
    integer :: ispin, iblock, iao
    real*8 :: Ndtot, Nd(2), VHart, EDCC(2)
    print *, "-------------------------------"
    print *, "Calculation of DFT+U potential:"
    print *, "-------------------------------"
    do iblock=1,NCorrBl
       print '(A,I2)', " Correlated Block #"
       print *, "Orbital occupations:"
       print '(A,14(F10.5))', " Spin-up:   ", (P(1,iao,iao),iao=CorrBeg(iblock),CorrEnd(iblock)) 
       print '(A,14(F10.5))', " Spin-down: ", (P(NDSpin,iao,iao),iao=CorrBeg(iblock),CorrEnd(iblock)) 
       ! Compute number of electrons in correlated subspace
       Nd(1) = d_zero; Nd(2) = d_zero
       do iao=CorrBeg(iblock),CorrEnd(iblock)
          Nd(1)=Nd(1)+P(1,iao,iao)
          Nd(2)=ND(2)+P(NDSpin,iao,iao)
       end do
       Ndtot = Nd(1)+Nd(2)
       ! Hartree term
       VHart = UCoul(iblock)*Ndtot
       ! Fully localized DCC
       EDCC(1) = UCoul(iblock)*(Ndtot-0.5) - JHund(iblock)*(Nd(1)-0.5)
       EDCC(2) = UCoul(iblock)*(Ndtot-0.5) - JHund(iblock)*(Nd(2)-0.5)
       print *, "Local DFT+U potential:"
       do ispin=1,NDSpin
          print '(A,14(F10.5))', " Spin-up:   ", &
            (UCoul(iblock)*(Ndtot-P(1,iao,iao))-JHund(iblock)*(Nd(1)-P(1,iao,iao))-EDCC(1), &
               iao=CorrBeg(iblock),CorrEnd(iblock))
          print '(A,14(F10.5))', " Spin-down: ", &
            (UCoul(iblock)*(Ndtot-P(NDSpin,iao,iao))-JHund(iblock)*(Nd(NDSpin)-P(NDSpin,iao,iao))-EDCC(NDSpin), &
               iao=CorrBeg(iblock),CorrEnd(iblock))
          do iao=CorrBeg(iblock),CorrEnd(iblock)
             ! Add on-site Hartree-Fock term minus DCC
             H(ispin,iao,iao) = H(ispin,iao,iao) &
                  ! HF term for direct repulsion U
                  + UCoul(iblock)*(Ndtot-P(ispin,iao,iao)) &
                  ! Fock term for Hund's rule coupling
                  - JHund(iblock)*(Nd(ispin)-P(ispin,iao,iao)) &
                  ! Double counting correction
                  - EDCC(ispin)
          end do
       end do
       print *
    end do
  end subroutine Add_DFT_plus_U_Pot

  !*************************************************
  ! Diagonalization of a correlated subspace       !
  ! prior to calculation of hybridization function !
  !*************************************************
  subroutine DiagCorrBlocks( HD, SD )
    use parameters
    use numeric, only: RMatPow, RSDiag, RSetId
    use util
    implicit none
    
    ! Orthogonalization and diagonalization matrix
    real*8, dimension(:,:,:), intent(inout) :: HD
    real*8, dimension(:,:), intent(inout) :: SD

    integer :: i, j, ibeg, iend, iblock, ispin, info, iunit, fstat
    character(len=100) :: istr

    real*8, dimension(NDAO,NDAO) :: temp, Trafo
    real*8, dimension(NDAO) :: E

    print *
    print *, "===================================="
    print *, "Diagonalization of correlated blocks"
    print *, "===================================="
    print *
    !
    ! 1. Individual Orthogonalization of each correlated block
    !
    call RSetID(temp)
    do iblock=1,NCorrBl
       ibeg = CorrBeg(iblock)
       iend = CorrEnd(iblock)
       print *, "  Block", iblock
       do ispin=1,NDSpin
          if( NDSpin == 2 .and. ispin == 1 ) print *, "Spin-up HD = "
          if( NDSpin == 2 .and. ispin == 2 ) print *, "Spin-down HD = "
          if( NDSpin == 1 ) print *, "HD = "
          call PrintRMatrix( HD(ispin,ibeg:iend,ibeg:iend) )
       end do
       temp(ibeg:iend,ibeg:iend) = SD(ibeg:iend,ibeg:iend)
    end do
    call RMatPow( temp, -0.5d0, Trafo)
    temp = matmul( Trafo, SD  )
    SD = matmul( temp, Trafo  )
    do ispin=1,NDSpin
       temp = matmul( HD(ispin,:,:), Trafo )
       HD(ispin,:,:) = matmul( Trafo, temp )
    end do
    !
    ! 2. Individual Diagonalization of each correlated block
    !
    call RSetID(Trafo)    
    do iblock=1,NCorrBl
       ibeg = CorrBeg(iblock)
       iend = CorrEnd(iblock)
       print '(A,I2,A,I4,A,I4)', " CorrBlock #", iblock, "  AOs: ", ibeg, "-", iend
       !
       ! Diagonalize Hamiltonian of correlated subspace
       !
       if( NDSpin == 1 ) Trafo(ibeg:iend,ibeg:iend) = HD(ibeg:iend,ibeg:iend,1)
       ! For spin-polarized systems: diagonalize average Hamiltonian
       if( NDSpin == 2 ) Trafo(ibeg:iend,ibeg:iend) = HD(ibeg:iend,ibeg:iend,1)+HD(ibeg:iend,ibeg:iend,2)
       call RSDiag( Trafo(ibeg:iend,ibeg:iend), E(ibeg:iend), info )          
       if( info /= 0 ) then
          print *, "Error: Correlation/DiagCorrBlocks: Diagonalization failure. Abort."
          stop
       end if
    end do
    !
    ! Apply unitary transformation
    !
    temp = matmul( SD, Trafo )
    SD = matmul( transpose(Trafo), temp )
    do ispin=1,NDSpin
       temp = matmul( HD(ispin,:,:), Trafo )
       HD(ispin,:,:) = matmul( transpose(Trafo), temp )
    end do
    !
    ! Print eigen energies and eigen states
    !
    do iblock=1,NCorrBl
       ibeg = CorrBeg(iblock)
       iend = CorrEnd(iblock)
       do ispin=1,NDSpin
          if(NDSpin == 1) print *, "Eigenenergies:"
          if(NDSpin == 2 .and. ispin == 1) print *, "Eigenenergies (spin-up):"
          if(NDSpin == 2 .and. ispin == 2) print *, "Eigenenergies (spin-down):"
          print '(21F14.8)', (HD(ispin,i,i),i=ibeg,iend)
       end do
       print *, "Eigenvectors (column-wise):"
       call PrintRMatrix(Trafo(ibeg:iend,ibeg:iend))
       print *
    end do

  end subroutine DiagCorrBlocks
  
     
end module Correlation
