!**********************************************************
!*********************  ANT.G-2.5.2  **********************
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
  ! Read a complex matrix A from a file connected to iunit
  !
  ! with sparse matrix format: i j Aij
  ! End of matrix marker: 0 0 0
  !
  subroutine ReadSparseMatrix( iunit, A )
    use constants
    use messages
    implicit none

    integer, intent(in) :: iunit
    complex(double),dimension(:,:),intent(out) :: A

    integer :: i,j, ios, dim1, dim2
    complex(double) :: Aij
    character :: ch
    integer :: lnum

    dim1 = size(A,1) 
    dim2 = size(A,2) 

    lnum = 0
    do
       lnum=lnum+1
       ! Skip comment lines
       read (unit=iunit,fmt=*,iostat=ios), ch
       if( ch == '!' ) cycle
       backspace iunit
       read (unit=iunit,fmt=*,iostat=ios), i, j, Aij
       if( ios /= 0 )then
          print *, "matrix line number = ", lnum
          call ReadErr( "ReadSparseMatrix", "i, j, Aij" )
       end if
       if( i == 0 .and. j == 0 ) exit
       if( i > dim1 .or. i < 1 .or. j > dim2 .or. j < 1)&
            call ErrMessage( "Util/ReadSparseMatrix", "Matrix index out of range.", .true. )
       A(i,j) = Aij
    end do
    
  end subroutine ReadSparseMatrix

  !*******************************
  ! Print Matrix in Sparse format
  !*******************************
  subroutine PrintSparseMatrix( A )
    use constants
    implicit none
    real(double),dimension(:,:),intent(in) :: A
    integer :: i,j,dim1,dim2
    dim1 = size(A,1) 
    dim2 = size(A,2) 
    do i=1,dim1
       do j=1,dim2
          if(abs(A(i,j))>1.0d-6) print *, i, j, A(i,j)
       end do
    end do
    print *, 0, 0, 0.0d0
  end subroutine PrintSparseMatrix

  !******************************
  ! Read a complex matrix A from 
  ! a file connected to iunit
  ! with normal matrix format: 
  !
  !  A11 A12 ... A1N
  !  A21 A22 ... A2N
  !   :   :       :
  !  AM1 AM2 ... AMN
  !*****************************
  subroutine ReadMatrix( iunit, A )
    use constants
    use messages
    implicit none

    integer, intent(in) :: iunit
    complex(double),dimension(:,:),intent(out) :: A

    integer :: i,j, ios, dim1, dim2
    character :: ch

    dim1 = size(A,1) 
    dim2 = size(A,2) 
    ! Skip comment lines
    do 
       read (unit=iunit,fmt=*,iostat=ios), ch
       if( ch /= "!" )exit
    end do
    backspace iunit
    do i=1,dim1
       read(unit=iunit,fmt=*,iostat=ios), (A(i,j),j=1,dim2)
       if( ios /= 0 ) call ReadErr("Util/ReadMatrix", "A(i,j)")
    end do
  end subroutine ReadMatrix
  
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
  
  subroutine PrintCMatrix( A )
    implicit none
    complex*16,dimension(:,:),intent(in) :: A
    integer :: i,j,dim1,dim2
    dim1 = size(A,1)
    dim2 = size(A,2)
    do i=1,dim1
       !print '(1000(ES14.4))', ( A(i,j), j=1,dim2 )
       Write(*,'(100(g15.5,g15.5,2x))') ( real(A(i,j)),&
                                        AIMAG(A(i,j)), j=1,dim2 )
    end do
  end subroutine PrintCMatrix  
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
