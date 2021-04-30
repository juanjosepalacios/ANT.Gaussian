!***************************************
!*                                     *
!*  ANT1D - messages.f90               *
!*                                     *
!*  Standardized error messages        *
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
module Messages
  implicit none
  
contains
  !
  ! Aborts program safely
  !
  subroutine AbortProg
    implicit none
    stop
  end subroutine AbortProg   
  !
  ! Generic error message
  !
  subroutine ErrMessage( routine, message, abort )
    implicit none
    !
    ! Takes name of routine where error occured
    ! message string, and logical flag whether to abort
    ! program execution as dummy arguments
    !
    character(len=*), intent(in) :: routine, message
    logical,intent(in),optional :: abort
    print '(1000A)', routine, "/Error: ", message
    if( present(abort) .and. abort ) call AbortProg()
  end subroutine ErrMessage
  ! 
  ! Failure allocating memory
  ! 
  subroutine AllocErr( routine, array )
    implicit none
    character(len=*), intent(in) :: routine, array
    call ErrMessage( routine, "Failure allocating memory for "//array//".", .true.) 
  end subroutine AllocErr
  !
  ! Error opening a file
  !
  subroutine FileErr( routine, fname )
    implicit none
    character(len=*), intent(in) :: routine, fname
    call ErrMessage( routine, "Could not open file '"//fname//"'.", .true.)
  end subroutine FileErr
  ! 
  ! Failure parsing NAMELIST
  ! 
  subroutine NMLErr( routine, nml_name )
    implicit none
    character(len=*), intent(in) :: routine, nml_name
    call ErrMessage( routine, "Failure parsing namelist '"//nml_name//"'.",.true.)
  end subroutine NMLErr
  !
  ! Error reading a variable from file
  !
  subroutine ReadErr( routine, var_name )
    implicit none
    character(len=*), intent(in) :: routine, var_name
    call ErrMessage( routine, "Could not read '"//var_name//"'.", .true.)
  end subroutine ReadErr
  !
  ! Keyword was not found or read incorrectly
  !
  subroutine NotFoundErr( routine, keyword )
    implicit none
    character(len=*), intent(in) :: routine, keyword
    call ErrMessage( routine, "Keyword '"//keyword//"' not found.", .true.)
  end subroutine NotFoundErr
  !
  ! Missing END statement for some input Field
  !
  subroutine MissEndErr( routine, keyword )
    implicit none
    character(len=*), intent(in) :: routine, keyword
    call ErrMessage( routine, "Missing 'END' statement for keyword '"//keyword//"'.", .true.)
  end subroutine MissEndErr

end module Messages
