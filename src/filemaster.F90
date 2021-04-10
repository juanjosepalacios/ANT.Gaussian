!***************************************
!*                                     *
!*  ANT1D - FileMaster.f90             *
!*                                     *
!*  Administration of File units       *
!*                                     *
!***************************************
!*                                     *
!*  This source file is part of the    *
!*  ANT1D project.                     *
!*                                     *
!*  Copyright (c) 2006 - 2015 by       *
!*                                     *
!*  David Jacob                        *
!*                                     *
!*  MPI fuer Mikrostrukturphysik       *
!*  Weinberg 2                         *
!*  06120 Halle                        *
!*  Germany                            *
!*                                     *
!***************************************

!**************************************
! Module for administrating file units
!**************************************
! Warning: don't open/close files with 
! standard Fortran OPEN/CLOSE statement
! if using this module
module FileMaster
  implicit none
  save
  private

  public :: fopen, fclose

  !****************** 
  !*** File Units ***
  !******************

  integer, parameter :: min_unit=10, max_unit=1000
  integer, dimension(min_unit:max_unit) :: units = 0

!!$  ! Main input file
!!$  integer, parameter :: iunit_inpf = 10
!!$
!!$  ! Device Hamiltonian and Overlap files
!!$  integer, parameter :: iunit_hdf = 20, iunit_sdf = 21
!!$
!!$  ! Leads Hamiltonian and overlap
!!$  integer,dimension(2), parameter :: iunit_vlf = (/30,31/), iunit_slf = (/32,33/)
!!$ 
!!$  ! Lead bulk DOS
!!$  integer, parameter :: iunit_ldosf = 34
  
contains

  !***************************************
  ! Open a file connecting to a free unit
  !***************************************
  ! Finds a free unit and connect 
  ! the file given by fname to it
  ! Unit will be marked as used if
  ! file is succesfully connected
  integer function fopen( fname, mode, ios )
    use util
    implicit none
    character(len=*), intent(in) :: fname, mode
    integer, intent(out) :: ios
    character(len_trim(mode)) :: mode_lc
    integer :: iunit
    iunit = FindFreeUnit() 
    if( iunit > 0 )then
       mode_lc=mode
       call locase(mode_lc)
       if( mode_lc == 'scratch' )then
          open( unit=iunit, status=mode, iostat=ios )
       else
          open( unit=iunit, file=fname, status=mode, iostat=ios )
       end if
       if( ios == 0 )units(iunit) = 1
       if( ios /= 0 )iunit = -2
    end if
    fopen=iunit
  end function fopen

  !**************************************
  ! Close file by disconneting from unit
  !**************************************
  ! Closes a file and frees the 
  ! unit that was connected to it
  subroutine fclose( iunit )
    implicit none
    integer, intent(in) :: iunit
    close( iunit )
    units(iunit)=0
  end subroutine fclose

  !****************************************
  ! Look for a free unit to connect a file
  !****************************************
  integer function FindFreeUnit()
    implicit none
    integer :: iunit
    do iunit=min_unit,max_unit
       FindFreeUnit=iunit
       if( units(iunit) == 0 )exit
    end do
    if( iunit > max_unit )FindFreeUnit=-1
  end function FindFreeUnit

end module FileMaster
