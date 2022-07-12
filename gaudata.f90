!
! Module containing general system data obtained from Gaussian
!
module gaudata

  ! Molecule specification
  integer*4 :: NAtoms  
  integer*4, dimension(:), allocatable :: IAN
  real*8, dimension(:), allocatable :: AtmChg, C
  
  integer*4 :: NE, NAE, NBE

  ! Basis set
  integer*4 :: NBasis
  integer*4, dimension(:), allocatable :: ibfatm, ibftyp
  
contains

  ! Obtains data from Gaussian using gauopen interface
  subroutine gaudata_init( natoms_, ian_, c_, ne_, nae_, nbe_, atmchg_, nbasis_, ibfatm_, ibftyp_  )
    implicit none
    integer*4, intent(in) :: natoms_, ne_, nae_, nbe_, nbasis_
    integer*4, dimension(:), intent(in) :: ian_, ibfatm_, ibftyp_
    real*8, dimension(:), intent(in) :: c_, atmchg_
    integer*4 :: err
  
    NAtoms = natoms_
    NE = ne_
    NAE = nae_
    NBE = nbe_
    NBasis = nbasis_
    
    allocate( IAN(natoms), ibfatm(nbasis), ibftyp(nbasis), C(3*natoms), AtmChg(NAtoms),  stat=err)
    if( err /= 0 )then
       print *, "Problems allocating memory..."
       STOP
    end if
    IAN = ian_
    ibfatm = ibfatm_
    ibftyp = ibftyp_
    C = c_
    AtmChg = AtmChg_

    print *, "NBasis = ", NBasis

    print *, "IAN = ", IAN
    print *, "C = ", C
    
    
  end subroutine gaudata_init

  ! Cleaning up: deallocate arrays etc.
  subroutine gaudata_cleanup
    implicit none
    deallocate( IAN, ibfatm, ibftyp, C, AtmChg )
  end subroutine gaudata_cleanup
  
end module gaudata
