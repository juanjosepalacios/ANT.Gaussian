!**********************************************************
!*********************  ANT.G-2.4.1  **********************
!**********************************************************
!*                                                        *
!*  Copyright (c) by                                      *
!*                                                        *
!*  David Jacob                                           *
!*                                                        *
!*     Theory Department                                  *
!*     Max-Planck-Institute for Microstructure Physics    *
!*     Halle, 06120 (GERMANY)                             *
!*                                                        *
!**********************************************************
  MODULE util
!**********************************************************
! Module with useful subroutines
!**********************************************************
  implicit none
contains
  !
  ! Pretty print a real Matrix on standard output
  !
  subroutine PrintRMatrix( A )
    implicit none
    real*8,dimension(:,:),intent(in) :: A
    integer :: i,j,dim1,dim2
    dim1 = size(A,1) 
    dim2 = size(A,2) 
    do i=1,dim1
       print '(1000(ES14.4))', ( A(i,j), j=1,dim2 )
    end do
  end subroutine PrintRMatrix
  !
  ! Convert <name> to upper case (useful for parsing)
  !
  subroutine upcase( name )
    implicit none
    character(LEN=*), intent(inout) :: name
    integer :: i, ic
    do i = 1,len(name)
       ic = iachar( name(i:i) )
       if( ic >= 97 .and. ic <= 122 ) name(i:i) = achar( ic-32 )
    end do
  end subroutine upcase
  !
  ! Convert <name> to lower case
  !
  subroutine locase( name )
    implicit none
    character(LEN=*), intent(inout) :: name
    integer :: i, ic
    do i = 1,len(name)
       ic = iachar( name(i:i) )
       if( ic >= 65 .and. ic <= 90 ) name(i:i) = achar( ic+32 )
    end do
  end subroutine locase
  !
  ! Convert an integer value to a string of digits
  !
  subroutine int2str( ival, str )
    implicit none
    
    integer, intent(in) :: ival
    character(len=*), intent(out) :: str
    
    integer :: idig, iival, ndigs, ipos, pos1
    
    if( ival == 0 )then
       str(1:1)='0'
       str(2:)=' '
       return
    end if

    iival=ival
    pos1=1
    ndigs=0

    if( iival < 0 )then
       iival=-iival
       ndigs=1
       str(1:1)='-'
       pos1=2
       ndigs = 1
    end if

    do while( iival /= 0 .and. ndigs<len(str))
       ndigs=ndigs+1
       idig = mod( iival, 10 )
       do ipos=ndigs-1,pos1,-1
          str( ipos+1:ipos+1 ) = str(ipos:ipos)
       end do
       str( pos1:pos1 ) = achar( idig+iachar('0') )
       iival = iival / 10
       ipos = ipos + 1
    end do
    if( iival /= 0 )then
       str(:) = '*'
       return
    end if
    str(ndigs+1:len(str))=' '

  end subroutine int2str

end module util
