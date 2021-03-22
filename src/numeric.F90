!*********************************************************!
!*********************  ANT.G-2.5.2  *********************!
!*********************************************************!
!                                                         !
!  Copyright (c) by                                       !
!                                                         !
!  Juan Jose Palacios (1)                                 !
!  David Jacob (2)                                        !
!                                                         !
! (1) Departamento de Fisica de la Materia Condensada     !
!     Universidad Autonoma de Madrid                      !    
!     28049 Madrid (SPAIN)                                !
! (2) Theory Department                                   !
!     Max-Planck-Institute for Microstructure Physics     !
!     Halle, 06120 (GERMANY)                              !
!                                                         !
!*********************************************************!
  MODULE numeric
!*********************************************************!
!  Module containing numerical routines                   !
!*********************************************************!
  IMPLICIT NONE
  CONTAINS

! *************************************************************
! Very simple sorting algorithm (not adequate for large arrays)
! *************************************************************
  SUBROUTINE sort(n,arr)
    IMPLICIT NONE

    INTEGER, INTENT(in) :: n
    REAL*8, DIMENSION(n), INTENT(inout) :: arr

    INTEGER :: i,j
    REAL*8 :: a

    DO j=2,n
       a=arr(j)
       DO i=j-1,1,-1
          IF(arr(i)<=a) EXIT 
          arr(i+1)=arr(i)
          arr(i) = a 
       END DO
    END DO
  END SUBROUTINE sort


  ! *******************************************************************
  ! Inversion of complex matrix based on LAPACK routines zgetrf, zgetri 
  ! *******************************************************************
  integer function CInv( Z )
#ifdef PGI
     USE lapack_blas, ONLY: zgetri,zgetrf
#endif
    IMPLICIT NONE
    external zgetri,zgetrf
    COMPLEX*16, DIMENSION(:,:),INTENT(inout) :: Z
    INTEGER :: n, ipiv(SIZE( Z, 1)), info
    COMPLEX*16, DIMENSION( 4*SIZE( Z, 1) ) :: work
    n = SIZE( Z, 1)
    CALL zgetrf(n,n,Z,n,ipiv,info)
    CALL zgetri(n,Z,n,ipiv,work,4*n,info)
    CInv=info
  END function CInv
  
  ! *******************************************************************
  ! Inversion of real matrix based on LAPACK routines dgetrf, dgetri 
  ! *******************************************************************
  integer function RInv( Z )
!    USE lapack_blas, ONLY: zgetri,zgetrf
!    use mkl_blas, ONLY: zgetri,zgetrf
    IMPLICIT NONE
    real, DIMENSION(:,:),INTENT(inout) :: Z
!    complex*16, DIMENSION(:,:) :: CZ
    complex*16, DIMENSION(SIZE( Z, 1),SIZE( Z, 2)) :: CZ
    integer :: n, ipiv(SIZE( CZ, 1)), info
    complex*16, DIMENSION( 4*SIZE( CZ, 1) ) :: work
    n = SIZE( Z, 1)
    CZ = dcmplx(Z)
    CALL zgetrf(n,n,CZ,n,ipiv,info)
!    CALL zgetrf(n,n,CZ,n+1,ipiv,info)
    CALL zgetri(n,CZ,n,ipiv,work,4*n,info)
!    CALL zgetri(n,CZ,n+1,ipiv,work,4*n+1,info)
    Z=real(CZ)
    RInv=info
  END function RInv

SUBROUTINE GAUSSJOLD(a,n,np,b,m,mp,ierr)

!   Purpose: Solution of the system of linear equations AX = B by
!      Gauss-Jordan elimination, where A is a matrix of order N and B is
!      an N x M matrix.  On output A is replaced by its matrix inverse
!      and B is preplaced by the corresponding set of solution vectors.

!   Source: W.H. Press et al, "Numerical Recipes," 1989, p. 28.

!   Modifications: 
!      1. Double  precision.
!      2. Error parameter IERR included.  0 = no error. 1 = singular 
!         matrix encountered; no inverse is returned.

!   Prepared by J. Applequist, 8/17/91.

!      implicit real(a-h,o-z)

!        Set largest anticipated value of N.
IMPLICIT NONE

real, DIMENSION(:,:), INTENT(INOUT) :: a,b
integer, INTENT(IN) :: n, np, m, mp
!integer, INTENT(IN) :: ierr
integer, INTENT(OUT) :: ierr
!Linear equation solution by Gauss-Jordan elimination, equation (2.1.1). a is an N ×N input
!coefficient matrix. b is an N × M input matrix containing M right-hand-side vectors. On
!output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of
!solution vectors.
integer*4, DIMENSION(size(a,1)) :: ipiv,indxr,indxc
!These arrays are used for bookkeeping on the pivoting.
logical, DIMENSION(size(a,1)) :: lpiv
real :: pivinv, big, dum
real, DIMENSION(size(a,1)) :: dumc
integer*4, TARGET :: irc(2)
integer :: i,l,nmax,j,k,ll
integer*4, POINTER :: irow,icol


      !parameter (nmax=500)
nmax=500
      !dimension a(np,np),b(np,mp),ipiv(nmax),indxr(nmax),indxc(nmax)
      ierr=0
do 11 j=1,n
  ipiv(j)=0
 11   continue
  do 22 i=1,n
    big=0.d0
    do 13 j=1,n
      if (ipiv(j).ne.1) then
        do 12 k=1,n
          if (ipiv(k).eq.0) then
            if (dabs(a(j,k)).ge.big) then
              big=dabs(a(j,k))
              irow=j
              icol=k
            endif
          else if (ipiv(k).gt.1) then
            ierr=1
            return
          endif
 12   continue
        endif
 13   continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
 14   continue
            do 15 l=1,m
              dum=b(irow,l)
              b(irow,l)=b(icol,l)
              b(icol,l)=dum
 15   continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.d0) then
          ierr=1
          return
        endif
        pivinv=1.d0/a(icol,icol)
        a(icol,icol)=1.d0
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
 16   continue
          do 17 l=1,m
            b(icol,l)=b(icol,l)*pivinv
 17   continue
            do 21 ll=1,n
              if (ll.ne.icol) then
                dum=a(ll,icol)
                a(ll,icol)=0.d0
                do 18 l=1,n
                  a(ll,l)=a(ll,l)-a(icol,l)*dum
 18   continue
                  do 19 l=1,m
                    b(ll,l)=b(ll,l)-b(icol,l)*dum
 19   continue
              endif
 21   continue
 22   continue
              do 24 l=n,1,-1
                if (indxr(l).ne.indxc(l)) then
                  do 23 k=1,n
                    dum=a(k,indxr(l))
                    a(k,indxr(l))=a(k,indxc(l))
                    a(k,indxc(l))=dum
 23   continue
                endif
 24   continue
                return
      END SUBROUTINE

! !ROUTINE: gaussj
!
! !INTERFACE:
      subroutine gaussj(a,n,np,b,m,mp)
      
! !DESCRIPTION:
!
! Linear equation solution by Gauss-Jordan elimination, equation
!{gaussj-1} below.
!{a(1:n,1:n)}
! is an imput matrix stored in an array of physical dimensions {np}
! by {np}.{b(1:n,1:m)} is an input matrix containing the
! {m} right-hand side vectors, stored in an array of physical
! dimensions {np} by {mp}. On output, {a(1:n,1:n)} is
! replaced by its matrix inverse, and {b(1:n,1:m)} is replaced by
! the corresponding set of solution vectors.
!
! !INPUT PARAMETERS:
      
      implicit none
      integer, intent(in) :: m
      integer, intent(in) :: mp
      integer, intent(in) :: n
      integer, intent(in) :: np
      
! !INPUT/OUTPUT PARAMETERS:

      real, intent(inout) :: a(1:np,1:np)
      real, intent(inout) :: b(1:np,1:mp)

! !LOCAL VARIABLES:
      
      integer :: i
      integer :: icol
      integer :: irow
      integer :: j
      integer :: k
      integer :: l
      integer :: ll
      integer :: indxc(1:np) ! Used for bookkeeping on the pivoting
      integer :: indxr(1:np) ! Used for bookkeeping on the pivoting
      integer :: ipiv(1:np) ! Used for bookkeeping on the pivoting
      
      real :: big
      real :: dum
      real :: pivinv

! !REVISION HISTORY:
!
! Original subroutine: gaussj.for (c) copr. 1986-92 numerical recipes
! software &30i..
! Last modified: 7th. Jul. 2005 by RGA
!
!EOP
!BOC
      do j=1,n
        ipiv(j)=0
      enddo ! j
!
! This is the main loop over the columns to be reduced
!
      do i=1,n
        big=0.0d0
!
! This is the outer loop of the search for a pivot element
!
        do j=1,n
          if(ipiv(j).ne.1)then
do k=1,n
              if (ipiv(k).eq.0) then
if (abs(a(j,k)).ge.big)then
big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
else if (ipiv(k).gt.1) then
stop 'singular matrix in gaussj'
              endif
enddo ! k
          endif
enddo ! j
        ipiv(icol)=ipiv(icol)+1
!
! We now have the pivot element, so we interchange rows, if needed,
! to put the pivot element on the diagonal. The columns are not
! physically interchanged, only relabeled: indxc(i), the column of
! the ith pivot element, is the ith column that is reduced, while
! indxr(i) is the row in which that pivot element was originally
! located. If indxr(i) neq indxc(i) there is an implied column
! interchange. With this form of bookkeeping, the solution b's will
! end up in the correct order, and the inverse matrix will be
! scrambled by columns.
!
        if (irow.ne.icol) then
do l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
          enddo ! l
          do l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
          enddo ! l
        endif
indxr(i)=irow
        indxc(i)=icol
!
! We are now ready to divide the pivot row by the pivot element,
! located at irow and icol.
!
        if (a(icol,icol).eq.0.0d0) stop 'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do l=1,n
          a(icol,l)=a(icol,l)*pivinv
        enddo ! l
        do l=1,m
          b(icol,l)=b(icol,l)*pivinv
        enddo ! l
!
! Next, we reduce the rows...
!
        do ll=1,n
          if(ll.ne.icol)then ! ... except for the pivor one, of course.
            dum=a(ll,icol)
            a(ll,icol)=0.0d0
            do l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
            enddo ! l
            do l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
            enddo ! l
          endif
enddo ! ll
      enddo ! i
!
! This is the end of the main loop over columns of the reduction. It
! only remains to unscramble the solution in view of the column
! interchanges. We do this by interchanging pairs of columns in the
! reverse order that the permutation was bilt up.
!
      do l=n,1,-1
        if(indxr(l).ne.indxc(l))then
do k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
          enddo ! k
        endif
enddo ! l

      return

end subroutine gaussj

SUBROUTINE GaussInverse(A)
IMPLICIT NONE
      real, dimension(:,:), intent(inout) :: A
      !real, dimension(:,:), intent(out) :: Ainv
      real, allocatable, dimension(:,:) :: Apiv
      integer :: i,n
      n=size(A,2)
      allocate(Apiv(n,n))
        Apiv(:,:) = 0.0
        do i=1,n
          Apiv(i,i) = 1.0
        end do
        call GAUSSJ(A,n,n,Apiv,n,n)
      return
      deallocate(Apiv)
END SUBROUTINE GaussInverse  
  

  ! *************************
  ! Trace of a complex matrix
  ! *************************
  REAL*8 FUNCTION RTrace( RMat )
    USE constants, ONLY: d_zero
    IMPLICIT NONE

    REAL*8, DIMENSION(:,:) :: RMat
    INTEGER :: i,n
    n=SIZE( RMAT, 1 )
    RTrace =d_zero
    DO i=1,n
       RTrace = RTrace + RMat(i,i)
    END DO
  END FUNCTION RTrace

  ! *************************
  ! Trace of a complex matrix
  ! *************************
  COMPLEX*16 FUNCTION CTrace( CMat )
    USE constants, ONLY: c_zero
    IMPLICIT NONE

    COMPLEX*16, DIMENSION(:,:) :: CMat
    INTEGER :: i,n
    n=SIZE( CMAT, 1 )
    CTrace =c_zero
    DO i=1,n
       CTrace = CTrace + CMat(i,i)
    END DO
  END FUNCTION CTrace

  !
  ! Set a matrix to Identity 
  !
  subroutine CSetId( A )
    implicit none
    complex*16,dimension(:,:),intent(out) :: A
    integer :: i
    A = 0
    do i=1,min(size(A,1),size(A,2))
       A(i,i)=1
    end do
  end subroutine CSetId

  !
  ! Set a matrix to Identity 
  !
  subroutine RSetId( A )
    implicit none
    real*8,dimension(:,:),intent(out) :: A
    integer :: i
    A = 0
    do i=1,min(size(A,1),size(A,2))
       A(i,i)=1
    end do
  end subroutine RSetId


  ! *************************************
  ! Gauss-Legendre quadrature
  ! - computes abscissas and weights for
  !   Gaussian quadrature using Legendre
  !   polynomials 
  ! x1,x2 : real integration interval
  ! x(n)  : array of abscissas            
  ! w(n)  : array of weights
  ! n     : number of abscissas
  ! *************************************
  SUBROUTINE gauleg(x1,x2,x,w,n)
    USE constants, ONLY: d_one, d_zero 
    IMPLICIT NONE

    REAL*8, INTENT(in) :: x1, x2
    REAL*8, DIMENSION(n), INTENT(out) :: x, w
    INTEGER, INTENT(in) :: n
    
    REAL*8, PARAMETER :: EPS = 3.0d-14, Pi=4.0d0*ATAN(1.0)!Pi=3.141592654d0
    INTEGER :: i,j,m,ncycle
    REAL*8 :: p1,p2,p3,pp,xl,xm,z,z1
    
    ! number of roots to find (symmetric in interval)
    m=(n+1)/2        
    xm=0.5d0*(x2+x1) 
    xl=0.5d0*(x2-x1)
    ! loop over roots to find
    DO i=1,m 
       z=DCOS(Pi*(i-0.25d0)/(n+0.5d0)) ! starting value for i-th root
       DO ! loop to find i-th root
          p1 = d_one  
          p2 = d_zero 
          ! p1 will contain n-th Legendre polynomial evaluated at z, 
          ! p2 the polynomial of one lower order after the next loop 
          ! down the recurrence relation for generating Legendre polynomials 
          DO j=1,n    
             p3=p2
             p2=p1
             p1=((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
          END DO
          ! derivative of polynomial n at z
          pp=n*(z*p1-p2)/(z*z-1.0d0)
          z1=z
          z=z1-p1/pp
          ! convergence criterion
          IF(ABS(z-z1)<EPS)EXIT
       END DO 
       ! scale to interval (x1,x2)
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.0d0*xl/((1.0d0-z*z)*pp*pp)
       w(n+1-i)=w(i)
    END DO
  END SUBROUTINE gauleg
  

  ! ************************
  ! Polynomial interpolation
  ! ************************
  SUBROUTINE polint(xa,ya,n,x,y,dy)
    IMPLICIT NONE 
    INTEGER, INTENT(in) :: n
    REAL*8, DIMENSION(n),INTENT(in) :: xa, ya
    REAL*8, INTENT(in)  :: x
    REAL*8, INTENT(out) :: y, dy

    INTEGER, PARAMETER :: nmax=10
    
    INTEGER :: i,m,ns
    REAL*8 :: den, dif, dift, ho, hp, w
    REAL*8, DIMENSION(nmax) :: c, d

    ns=1
    dif=ABS(x-xa(1))
    DO i=1,n
       dift=ABS(x-xa(i))
       IF( dift < dif ) THEN
          ns=i
          dif=dift
       END IF
       c(i) = ya(i)
       d(i) = ya(i)
    END DO
    y=ya(ns)
    ns=ns-1
    DO m=1,n-1
       DO i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          IF( den == 0.0d0 )THEN
             PRINT *, "POLINT/Failure in subroutine. Abort."
             STOP
          END IF
          den = w/den
          d(i)= hp*den
          c(i)= ho*den
       END DO
       IF( 2*ns < n-m )THEN
          dy=c(ns+1)
       ELSE
          dy=d(ns)
          ns=ns-1
       END IF
       y=y+dy
    END DO
  END SUBROUTINE polint


  ! *************************************
  ! Romberg integration on open interval 
  ! *************************************
  SUBROUTINE qromo(func,a,b,ss,rule)

    IMPLICIT NONE

    REAL*8, EXTERNAL :: func
    REAL*8, INTENT(in) :: a,b
    REAL*8, INTENT(inout) :: ss

    INTEGER, PARAMETER :: jmax=14, jmaxp=jmax+1, k=5, km=k-1
    REAL*8, PARAMETER :: eps = 1.0d-6
    
    INTEGER :: j
    REAL*8 :: dss
    REAL*8, DIMENSION(jmaxp) :: h, s

    h(1) = 1.0d0
    DO j=1,jmax
       CALL rule(func,a,b,s(j),j)
       IF( j >= k )THEN
          CALL polint(h(j-km),s(j-km),k,0.0d0,ss,dss)
          IF(ABS(dss)<=eps*ABS(ss))RETURN
       END IF
       s(j+1)=s(j)
       h(j+1)=h(j)/9.0d0
    END DO

    PRINT *, "QROMO/Too many quadrature steps."
    STOP
  END SUBROUTINE qromo
  
  ! *************************************
  ! Romberg integration on open interval
  ! *************************************
  SUBROUTINE qromo1DBL(func,a,b,ss,rule)

    IMPLICIT NONE

    real, EXTERNAL :: func
    real, INTENT(in) :: a,b
    real, INTENT(inout) :: ss

    integer, PARAMETER :: jmax=14, jmaxp=jmax+1, k=5, km=k-1 ! ORIGINAL jmax=14
!    integer, PARAMETER :: jmax=28, jmaxp=jmax+1, k=5, km=k-1 ! ORIGINAL jmax=14
!    real, PARAMETER :: eps = 1.0d-12 ! ORIGINAL CODE.
    real, PARAMETER :: eps = 1.0d-9 ! ORIGINAL THIS OR 1.0d-12
!    real, PARAMETER :: eps = 1.0d-4 ! ADDED BY C.SALGADO ON 2017-03-28.
    integer :: j
    real :: dss
    real, DIMENSION(jmaxp) :: h, s

    h(1) = 1.0d0
    DO j=1,jmax
       CALL rule(func,a,b,s(j),j)
       IF( j >= k )THEN
          CALL polint(h(j-km),s(j-km),k,0.0d0,ss,dss)
          Write(*,'(g15.5,A,g15.5,A)')ABS(dss),"  <",eps*ABS(ss)," ?"
          IF(ABS(dss)<=eps*ABS(ss))RETURN
       END IF
       s(j+1)=s(j)
       h(j+1)=h(j)/9.0d0
    END DO

    PRINT *, "QROMO1DBL/Too many quadrature steps."
    STOP
  END SUBROUTINE qromo1DBL

  ! *************************************
  ! Romberg integration on open interval 
  ! *************************************
  SUBROUTINE qromo1DBLOLD(func,a,b,ss,rule)

    IMPLICIT NONE

    REAL*8, EXTERNAL :: func
    REAL*8, INTENT(in) :: a,b
    REAL*8, INTENT(inout) :: ss

!    INTEGER, PARAMETER :: jmax=14, jmaxp=jmax+1, k=5, km=k-1 ! ORIGINAL jmax=14
    INTEGER, PARAMETER :: jmax=28, jmaxp=jmax+1, k=5, km=k-1
!    REAL*8, PARAMETER :: eps = 1.0d-6 ! ORIGINAL THIS OR 1.0d-12
    REAL*8, PARAMETER :: eps = 2.0d-4
    
    INTEGER :: j
    REAL*8 :: dss
    REAL*8, DIMENSION(jmaxp) :: h, s

    h(1) = 1.0d0
    DO j=1,jmax
       CALL rule(func,a,b,s(j),j)
       IF( j >= k )THEN
          CALL polint(h(j-km),s(j-km),k,0.0d0,ss,dss)
          Write(*,*)dss
          Write(*,*)ss
          Write(*,'(g15.5,A,g15.5,A)')ABS(dss),"  <",eps*ABS(ss)," ?"
          Write(*,'(F16.12,A,F16.12,A)')ABS(dss),"  <",eps*ABS(ss)," ?"
          IF(ABS(dss)<=eps*ABS(ss))RETURN
       END IF
       s(j+1)=s(j)
       h(j+1)=h(j)/9.0d0
    END DO

    PRINT *, "QROMO/Too many quadrature steps."
    STOP
  END SUBROUTINE qromo1DBLOLD

  SUBROUTINE qromoc(func,a,b,ss,rule)
    IMPLICIT NONE

    real, EXTERNAL :: func
    real, INTENT(in) :: a,b
    real, INTENT(inout) :: ss

    EXTERNAL rule

    integer, PARAMETER :: jmax=14, jmaxp=jmax+1, k=9, km=k-1
    !integer, PARAMETER :: jmax=140, jmaxp=jmax+1, k=5, km=k-1
    !real, PARAMETER :: eps = 1.0d-6
    real, PARAMETER :: eps = 1.0d-12
    
    integer :: j
    real :: dss
    real, DIMENSION(jmaxp) :: h, s

    !PRINT *, "I am in QROMOC"

    h(1) = 1.0d0
    DO j=1,jmax
       CALL rule(func,a,b,s(j),j)
       !PRINT *, j,s(j),h(j)
       IF( j >= k )THEN
          CALL polint(h(j-km),s(j-km),k,0.0d0,ss,dss)
          IF(ABS(dss)<=eps*ABS(ss))RETURN
       END IF
       s(j+1)=s(j)
       h(j+1)=h(j)/9.0d0
    END DO
    PRINT *, "QROMO/Too many quadrature steps."
    STOP
  END SUBROUTINE qromoc
  
  ! **********************************
  ! Midpoint rule for finite intervals
  ! **********************************
  SUBROUTINE midpnt(func,a,b,s,n)
    IMPLICIT NONE

    REAL*8, EXTERNAL :: func
    REAL*8, INTENT(in) :: a,b
    REAL*8, INTENT(inout) :: s
    INTEGER, INTENT(in) ::  n

    INTEGER :: it, j
    REAL*8 :: ddel, del, sum ,tnm, x
    
    IF( n == 1 ) THEN
        s=(b-a)*func(0.5d0*(a+b))
    ELSE
       it  = 3**(n-2)
       tnm = it
       del = (b-a)/(3.0d0*tnm)
       ddel = del + del
       x = a + 0.5d0*del
       sum = 0.0d0
       DO j=1,it
          sum = sum + func(x)
          x = x+ddel
          sum = sum + func(x)
          x = x+del
       END DO
       s = (s+(b-a)*sum/tnm)/3.0d0
    END IF 
  END SUBROUTINE midpnt
  

  ! *************************************
  ! Midpoint rule for infinite intervals
  ! *************************************
   SUBROUTINE midinf(funk, aa, bb, s, n )
    IMPLICIT NONE

    REAL*8, EXTERNAL :: funk
    REAL*8, INTENT(in)    :: aa, bb
    REAL*8, INTENT(inout)   :: s
    INTEGER, INTENT(in) :: n
    
    INTEGER :: it, j
    REAL*8 :: a,b, ddel, del, sum, tnm, func, x
    
    func(x) = funk(1.0d0/x)/x**2
    
    b = 1.0d0/aa
    a = 1.0d0/bb
    
    IF( n == 1 ) THEN
       s=(b-a)*func(0.5d0*(a+b))
    ELSE
       it  = 3**(n-2)
       tnm = it
       del = (b-a)/(3.0d0*tnm)
       ddel = del + del
       x = a + 0.5d0*del
       sum = 0.0d0
       DO j=1,it
          sum = sum + func(x)
          x = x+ddel
          sum = sum + func(x)
          x = x+del
       END DO
       s = (s+(b-a)*sum/tnm)/3.0d0
    END IF
   END SUBROUTINE midinf


  ! ****************************************************
  ! Compute generalized power of a complex matrix: 
  ! X = Z^a = U z^a conjg(U)
  ! where z is the diagonalized matrix corresponding to Z
  ! ****************************************************
  SUBROUTINE CMatPow( Z, a , X )
    USE constants, ONLY: c_zero, c_one
#ifdef PGI
    USE lapack_blas, ONLY: zheev,zgemm
#endif
    use ANTCommon
    IMPLICIT NONE
    external zheev,zgemm
    COMPLEX*16, DIMENSION(:,:), INTENT(in)  :: Z
    REAL*8, INTENT(in) :: a
    COMPLEX*16, DIMENSION(SIZE(Z,1),SIZE(Z,2)), INTENT(out) :: X
    COMPLEX*16, DIMENSION(SIZE(Z,1),SIZE(Z,2)) :: U,tmp
    REAL*8, DIMENSION(SIZE(Z,1)) :: zev
    COMPLEX*16, DIMENSION(2*SIZE(Z,1)-1) :: work
    REAL*8, DIMENSION(3*SIZE(Z,1)-2) :: rwork
    complex*16, DIMENSION(SIZE(Z,1)) :: zev2
    INTEGER :: N, INFO, i

    N = SIZE(Z,1)
    U=Z

    ! 1. Diagonalize Z

    CALL ZHEEV('V','U',N,U,N,zev,WORK,2*N-1,RWORK,INFO)

    IF( INFO /= 0 ) THEN
       WRITE(ifu_log,*)'CMatPow/Problems diagonalizing Z'
       WRITE(ifu_log,*)'INFO=',INFO
       STOP
    END IF

    ! 2. compute diagonal matrix z^a
    X=c_zero
    zev2=c_zero
    DO i=1, N
          X(i,i) = (zev2(i)+zev(i))**a
    ENDDO

    ! 3. X = U z^a conjg(U)
    CALL ZGEMM ('N','N',N,N,N,c_one,U,N,X,N,c_zero,tmp,N )
    CALL ZGEMM ('N','C',N,N,N,c_one,tmp,N,U,N,c_zero,X,N )
  END SUBROUTINE CMatPow

  ! ****************************************************
  ! Compute generalized power a real matrix: 
  ! X = Z^a = U z^a U^T
  ! where z is the diagonalied matrix corresponding to Z
  ! ****************************************************
  SUBROUTINE RMatPow( Z, a , X )
     USE constants, ONLY: d_zero, d_one
#ifdef PGI
    USE lapack_blas, ONLY: dsyev,dgemm
#endif
    use ANTCommon
    IMPLICIT NONE
    external dsyev,dgemm
    
    REAL*8, DIMENSION(:,:), INTENT(in)  :: Z
    REAL*8, INTENT(in) :: a
    REAL*8, DIMENSION(SIZE(Z,1),SIZE(Z,2)), INTENT(out) :: X

    REAL*8, DIMENSION(SIZE(Z,1),SIZE(Z,2)) :: U,tmp
    REAL*8, DIMENSION(SIZE(Z,1)) :: zev
    REAL*8, DIMENSION(3*SIZE(Z,1)) :: work

    INTEGER :: N, INFO, i
    
    N = SIZE(Z,1)
    U=Z
    ! 1. Diagonalize Z
    CALL DSYEV('V','U',N,U,N,zev,WORK,3*N,INFO)

    IF( INFO /= 0 ) THEN
       WRITE(ifu_log,*)'RMatPow/Problems diagonalizing Z'
       WRITE(ifu_log,*)'INFO=',INFO
       STOP
    END IF
    ! 2. compute diagonal matrix z^a
    X=d_zero
    DO i=1, N
      !IF ( ABS(zev(i)) .GT. 1.0D-10 ) THEN 
          X(i,i) = zev(i)**a
      !ELSE
      !   WRITE(ifu_log,*)'RMatPow/Eigenvalues of matrix Z too small'
      !   WRITE(ifu_log,*)'zev=', zev
      !   STOP
      !ENDIF
    ENDDO
    ! 3. X = U z^a conjg(U)
    CALL DGEMM ('N','N',N,N,N,d_one,U,N,X,N,d_zero,tmp,N )
    CALL DGEMM ('N','T',N,N,N,d_one,tmp,N,U,N,d_zero,X,N )
  END SUBROUTINE RMatPow

  !**************************************
  !* Diagonalize real symmetric matrix  *
  !* 'Wrapper' for LAPACK routine DSYEV *
  !**************************************
  SUBROUTINE RSDiag( A, w, info )
#ifdef PGI
    USE lapack_blas, ONLY: dsyev
#endif
    IMPLICIT NONE
    external dsyev
    !
    ! Matrix A 
    ! on input:  matrix to be diagonalized
    ! on output: transformation matrix
    ! 
    REAL*8, DIMENSION(:,:), INTENT(inout)  :: A
    !
    ! array w - on output: eigen values
    !
    REAL*8, DIMENSION(:), INTENT(out) :: w
    !
    ! info - on output: 
    ! succesful: 0
    ! failed: >0 - not converged
    !         <0 - ilegal value of argument -info 
    ! 
    INTEGER, INTENT(out) :: info
    !
    REAL*8, DIMENSION(3*SIZE(A,1)) :: work
    INTEGER :: N
    !
    N = SIZE( A, 1 )
    CALL DSYEV('V','U',N,A,N,w,WORK,3*N,INFO)
  END SUBROUTINE RSDiag

  !****************************************
  !* Diagonalize complex hermitian matrix *
  !* 'Wrapper' for LAPACK routine ZHEEV   *
  !****************************************
  SUBROUTINE CHDiag( A, w, info )
#ifdef PGI
    USE lapack_blas, ONLY: zheev
#endif
    IMPLICIT NONE
    external zheev
    !
    ! Matrix A 
    ! on input:  matrix to be diagonalized
    ! on output: transformation matrix
    ! 
    COMPLEX*16, DIMENSION(:,:), INTENT(inout)  :: A
    !
    ! array w - on output: eigen values
    !
    REAL*8, DIMENSION(:), INTENT(out) :: w
    !
    ! info - on output: 
    ! succesful: 0
    ! failed: >0 - not converged
    !         <0 - ilegal value of argument -info 
    ! 
    INTEGER, INTENT(out) :: info
    !
    REAL*8, DIMENSION(3*SIZE(A,1)-2) :: rwork
    COMPLEX*16, DIMENSION(3*SIZE(A,1)) :: work
    INTEGER :: N
    !
    N = SIZE( A, 1 )
    CALL ZHEEV('V','U',N,A,N,w,WORK,3*N,RWORK,INFO)
  END SUBROUTINE CHDiag


  !****************************************
  !* Diagonalize general complex matrix   *
  !* 'Wrapper' for LAPACK routine ZGEEV   *
  !****************************************
  SUBROUTINE CDiag( A, w, info )
#ifdef PGI
    USE lapack_blas, ONLY: zgeev
#endif
    IMPLICIT NONE
    external zgeev
    !
    ! Matrix A 
    ! on input:  matrix to be diagonalized
    ! on output: transformation matrix
    ! 
    COMPLEX*16, DIMENSION(:,:), INTENT(inout)  :: A
    !
    ! array w - on output: eigen values
    !
    COMPLEX*16, DIMENSION(:), INTENT(out) :: w
    !
    ! info - on output: 
    ! succesful: 0
    ! failed: >0 - not converged
    !         <0 - ilegal value of argument -info 
    ! 
    INTEGER, INTENT(out) :: info
    !
    REAL*8, DIMENSION(2*SIZE(A,1)) :: rwork
    COMPLEX*16, DIMENSION(4*SIZE(A,1)) :: work
    COMPLEX*16, DIMENSION(SIZE(A,1),SIZE(A,1)) :: VL, VR
    INTEGER :: N
    !
    N = SIZE( A, 1 )
    
    CALL ZGEEV('N','N',N,     A,N, w,   VL,1,VR,1,  WORK,4*N,  RWORK, INFO ) 
  END SUBROUTINE CDiag

!----------------------------------------------------------------------------------------------
!------- Main code to call jacobi diagonalization of a real symmetric matrix ------------------
!---------------------------------------------------------------------------------------------- 
!====================================================================
!  eigenvalues and eigenvectors of a real symmetric matrix
!  Method: calls Jacobi
!====================================================================
!
!  call Jacobi(a,x,abserr,n)
!
SUBROUTINE JacobiOld(a,x,abserr,n)
!===========================================================
! Evaluate eigenvalues and eigenvectors
! of a real symmetric matrix a(n,n): a*x = lambda*x 
! method: Jacoby method for symmetric matrices 
! Alex G. (December 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - number of equations
! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! output ...
! a(i,i) - eigenvalues
! x(i,j) - eigenvectors
! comments ...
!===========================================================
implicit none
integer i, j, k, n
double precision a(n,n),x(n,n)
double precision abserr, b2, bar
double precision beta, coeff, c, s, cs, sc

! initialize x(i,j)=0, x(i,i)=1
! *** the array operation x=0.0 is specific for Fortran 90/95
x = 0.0
do i=1,n
  x(i,i) = 1.0
end do

! find the sum of all off-diagonal elements (squared)
b2 = 0.0
do i=1,n
  do j=1,n
    if (i.ne.j) b2 = b2 + a(i,j)**2
  end do
end do

if (b2 <= abserr) return

! average for off-diagonal elements /2
!bar = 0.5*b2/float(n*n)
bar = 0.5*b2/dfloat(n*n)

do while (b2.gt.abserr)
  do i=1,n-1
    do j=i+1,n
      if (a(j,i)**2 <= bar) cycle  ! do not touch small elements
      b2 = b2 - 2.0*a(j,i)**2
      bar = 0.5*b2/dfloat(n*n)
! calculate coefficient c and s for Givens matrix
      beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
      coeff = 0.5*beta/sqrt(1.0+beta**2)
      s = sqrt(max(0.5+coeff,0.0))
      c = sqrt(max(0.5-coeff,0.0))
! recalculate rows i and j
      do k=1,n
        cs =  c*a(i,k)+s*a(j,k)
        sc = -s*a(i,k)+c*a(j,k)
        a(i,k) = cs
        a(j,k) = sc
      end do
! new matrix a_{k+1} from a_{k}, and eigenvectors 
      do k=1,n
        cs =  c*a(k,i)+s*a(k,j)
        sc = -s*a(k,i)+c*a(k,j)
        a(k,i) = cs
        a(k,j) = sc
        cs =  c*x(k,i)+s*x(k,j)
        sc = -s*x(k,i)+c*x(k,j)
        x(k,i) = cs
        x(k,j) = sc
      end do
    end do
  end do
end do
return
END SUBROUTINE
!----------------------------------------------------------------------------------------------
!---------- NUMERICAL RECIPES IN F90 JACOBI IMPLEMENTATION- -----------------------------------
!----------------------------------------------------------------------------------------------
!---------- Necessary to symmetrize for Jacobi ------------------------------------------------
SUBROUTINE balanc(a)
!USE nrtype; USE nrutil
IMPLICIT NONE
real, DIMENSION(:,:), INTENT(INOUT) :: a
real, PARAMETER :: RADX=radix(a),SQRADX=RADX**2
!Given an N × N matrix a, this routine replaces it by a balanced matrix with identical
!eigenvalues. A symmetric matrix is already balanced and is unaffected by this procedure.
!The parameter RADX is the machine’s floating-point radix.
integer :: i,last,ndum
real :: c,f,g,r,s
!ndum=assert_eq(size(a,1),size(a,2),’balanc’)
do
 last=1
 do i=1,size(a,1)
  c=sum(abs(a(:,i)))-a(i,i)
  r=sum(abs(a(i,:)))-a(i,i)
  if (c /= 0.0 .and. r /= 0.0) then
   g=r/RADX
   f=1.0
   s=c+r
   do
    if (c >= g) exit
    f=f*RADX
    c=c*SQRADX
   end do
   g=r*RADX
   do
    if (c <= g) exit
    f=f/RADX
    c=c/SQRADX
   end do
   if ((c+r)/f < 0.95*s) then
    last=0
    g=1.0/f
    a(i,:)=a(i,:)*g
    a(:,i)=a(:,i)*f
   end if
  end if
 end do
 if (last /= 0) exit
end do
END SUBROUTINE balanc

!*************************************************************
!* This subroutine computes all eigenvalues and eigenvectors *
!* of a real symmetric square matrix A(N,N). On output, ele- *
!* ments of A above the diagonal are destroyed. D(N) returns *
!* the eigenvalues of matrix A. V(N,N) contains, on output,  *
!* the eigenvectors of A by columns. THe normalization to    *
!* unity is made by main program before printing results.    *
!* NROT returns the number of Jacobi matrix rotations which  *
!* were required.                                            *
!* --------------------------------------------------------- *
!* Ref.:"NUMERICAL RECIPES, Cambridge University Press, 1986,*
!*       chap. 11, pages 346-348" [BIBLI 08].                *
!*************************************************************
Subroutine Jacobi(A,N,D,V,NROT)
integer N,NROT,ip,iq,i,j
real  A(1:N,1:N),D(1:N),V(1:N,1:N)
real, pointer :: B(:), Z(:)
real  c,g,h,s,sm,t,tau,theta,tresh

!integer*WORDSIZE malloc !added by Carlos Salgado

!allocate(B(1:100),stat=ialloc) !original
!allocate(Z(1:100),stat=ialloc) !original
!allocate(B(1:100))
!allocate(Z(1:100))
!allocate(B(1:100),stat=malloc)
!allocate(Z(1:100),stat=malloc)
allocate(B(1:N))
allocate(Z(1:N))

  do ip=1, N    !initialize V to identity matrix
    do iq=1, N
      V(ip,iq)=0.d0 
    end do
      V(ip,ip)=1.d0
  end do  
  do ip=1, N
    B(ip)=A(ip,ip)
    D(ip)=B(ip)
    Z(ip)=0.d0    
  end do
  NROT=0
  do i=1, 50
    sm=0.d0
    do ip=1, N-1     !sum off-diagonal elements
      do iq=ip+1, N
        sm=sm+DABS(A(ip,iq))
      end do
    end do
    if(sm==0.d0) return  !normal return
    if(i.lt.4) then
      tresh=0.2d0*sm**2
    else
      tresh=0.d0
    end if
    do ip=1, N-1
      do iq=ip+1, N
        g=100.d0*DABS(A(ip,iq))
!write(*,*)'g=100.d0*DABS(A(ip,iq))',i
! after 4 sweeps, skip the rotation if the off-diagonal element is small
        if((i.gt.4).and.(DABS(D(ip))+g.eq.DABS(D(ip))) &
  .and.(DABS(D(iq))+g.eq.DABS(D(iq)))) then
    A(ip,iq)=0.d0
        else if(DABS(A(ip,iq)).gt.tresh) then
   h=D(iq)-D(ip)
   if(DABS(h)+g.eq.DABS(h)) then
     t=A(ip,iq)/h
          else
     theta=0.5d0*h/A(ip,iq)  
            t=1.d0/(DABS(theta)+DSQRT(1.d0+theta**2))
     if(theta.lt.0.d0) t=-t
          end if
   c=1.d0/DSQRT(1.d0+t**2)
!write(*,*)'c=1.d0/DSQRT(1.d0+t**2)',i
   s=t*c
          tau=s/(1.d0+c)
   h=t*A(ip,iq)
   Z(ip)=Z(ip)-h
   Z(iq)=Z(iq)+h
   D(ip)=D(ip)-h
   D(iq)=D(iq)+h
   A(ip,iq)=0.d0
   do j=1, ip-1
     g=A(j,ip)
     h=A(j,iq)
     A(j,ip)=g-s*(h+g*tau)
     A(j,iq)=h+s*(g-h*tau)
          end do
   do j=ip+1, iq-1
     g=A(ip,j)
     h=A(j,iq)
     A(ip,j)=g-s*(h+g*tau)
     A(j,iq)=h+s*(g-h*tau)
          end do		      
   do j=iq+1, N
     g=A(ip,j)
     h=A(iq,j)
     A(ip,j)=g-s*(h+g*tau)
     A(iq,j)=h+s*(g-h*tau)
          end do		  
   do j=1, N
     g=V(j,ip)
     h=V(j,iq)
     V(j,ip)=g-s*(h+g*tau)
     V(j,iq)=h+s*(g-h*tau)
          end do		  
          NROT=NROT+1
        end if !if ((i.gt.4)...
      end do !main iq loop
    end do !main ip loop
    do ip=1, N
      B(ip)=B(ip)+Z(ip)
      D(ip)=B(ip)
      Z(ip)=0.d0
    end do
!write(*,*)'end do main i loop',i
  end do !main i loop
  !pause ' 50 iterations !'
  return
END SUBROUTINE

! end of file ujacobi.f90

SUBROUTINE eigsrt(d,v)
!USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,swap
IMPLICIT NONE
real, DIMENSION(:), INTENT(INOUT) :: d
real, DIMENSION(:,:), INTENT(INOUT) :: v
!Given the eigenvalues d and eigenvectors v as output from jacobi (§11.1) or tqli (§11.3),
!this routine sorts the eigenvalues into descending order, and rearranges the columns of v
!correspondingly. The method is straight insertion.
integer :: i,j,n
!n=assert_eq(size(d),size(v,1),size(v,2),’eigsrt’)
n = SIZE(d,1)
do i=1,n-1
j=imaxloc_r(d(i:n))+i-1
if (j /= i) then
call swap_r(d(i),d(j))
call swap_rv(v(:,i),v(:,j))
end if
end do
END SUBROUTINE eigsrt

SUBROUTINE ordks(d,v)
!USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,swap
IMPLICIT NONE
real, DIMENSION(:), INTENT(INOUT) :: d
real, DIMENSION(:), INTENT(INOUT) :: v
!Given the eigenvalues d and eigenvectors v as output from jacobi (§11.1) or tqli (§11.3),
!this routine sorts the eigenvalues into descending order, and rearranges the columns of v
!correspondingly. The method is straight insertion.
integer :: i,j,n
!n=assert_eq(size(d),size(v,1),size(v,2),’eigsrt’)
n = SIZE(d,1)
do i=1,n-1
j=imaxloc_r(d(i:n))+i-1
if (j /= i) then
call swap_r(d(i),d(j))
call swap_r(v(i),v(j))
end if
end do
END SUBROUTINE ordks

FUNCTION imaxloc_r(arr)
real, DIMENSION(:), INTENT(IN) :: arr
integer :: imaxloc_r
integer, DIMENSION(1) :: imax
imax=maxloc(arr(:))
imaxloc_r=imax(1)
END FUNCTION imaxloc_r
FUNCTION imaxloc_i(iarr)
integer, DIMENSION(:), INTENT(IN) :: iarr
integer, DIMENSION(1) :: imax
integer :: imaxloc_i
imax=maxloc(iarr(:))
imaxloc_i=imax(1)
END FUNCTION imaxloc_i
SUBROUTINE swap_i(a,b)
integer, INTENT(INOUT) :: a,b
integer :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_i
SUBROUTINE swap_r(a,b)
real, INTENT(INOUT) :: a,b
real :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_r
SUBROUTINE swap_rv(a,b)
real, DIMENSION(:), INTENT(INOUT) :: a,b
real, DIMENSION(SIZE(a)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_rv
!----------------------------------------------------------------------------------------------

  ! *********************************
  ! Bisection method for root finding 
  ! *********************************
  REAL *8 FUNCTION Bisec(func,x1,x2,xacc,jmax,j)
    IMPLICIT NONE
    ! func  : real function to find root for
    ! x1,x2 : interval bracketing root
    ! xacc  : accuracy for root
    REAL*8, EXTERNAL :: func
    REAL*8, INTENT(in) :: x1,x2, xacc
    INTEGER, INTENT(IN):: jmax 

    INTEGER, intent(out) :: j

    REAL*8 :: dx, f, fmid, xmid
 
    fmid = func(x2)
    f = func(x1)
    
    IF( f*fmid >= 0 ) THEN
       PRINT *, "Bisec/Error: Root is not bracketed (x1,x2). Abort."
       STOP
    END IF

    IF( f <  0 )THEN
       bisec = x1
       dx = x2-x1
    ELSE
       bisec = x2
       dx = x1-x2
    END IF
    DO j=0,jmax-1
       dx=dx*0.5d0
       xmid=bisec+dx
       fmid=func(xmid)
       IF(fmid <= 0 ) bisec=xmid
       IF(ABS(dx) < xacc .OR. fmid == 0 ) RETURN
       !PRINT '(I3,A,F11.5,A,F11.5)', j, ".  xmid=", xmid, "  fmid=", fmid
    END DO
    PRINT *, "Bisec/Warning: Too many bisections."
    
  END FUNCTION Bisec


  !**************************************************************************
  ! NUMERICAL METHODS: FORTRAN Programs, (c) John H. Mathews 1994
  ! To accompany the text:
  ! NUMERICAL METHODS for Mathematics, Science and Engineering, 2nd Ed, 1992
  ! Prentice Hall, Englewood Cliffs, New Jersey, 07632, U.S.A.
  ! This free software is complements of the author.
  !
  ! Algorithm 2.8 (Muller's Method).
  ! Section 2.5, Aitken's Process & Steffensen's & Muller's Methods, Page 97
  !**************************************************************************
  SUBROUTINE MULLER(F,P0,P1,P2,Delta,Epsilon,Max,P3,Z,K,Cond)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: Max
    INTEGER, INTENT(out) :: K,Cond
    REAL*8 :: A,B,C,Det,Disc,H0,H1,E0,E1
    REAL*8 :: Y0,Y1,Y2,RelErr,U,V
    REAL*8, INTENT(in) :: Delta,Epsilon
    REAL*8, INTENT(out) :: P3,Z
    REAL*8, INTENT(inout) :: P0,P1,P2
    REAL*8, PARAMETER :: Small=1E-10
    REAL*8, EXTERNAL :: F
    K=0
    Cond=0
    Y0=F(P0)
    Y1=F(P1)
    Y2=F(P2)
    DO WHILE ((K.LT.Max).AND.(Cond.EQ.0))
      !print *, "k =", k
      !print *, " P0 = ", P0, Y0
      !print *, " P1 = ", P1, Y1
      !print *, " P2 = ", P2, Y2
       H0=P0-P2
       H1=P1-P2
       C=Y2
       E0=Y0-C
       E1=Y1-C
       Det=H0*H1*(H0-H1)
       A=(E0*H1-H0*E1)/Det
       B=(H0*H0*E1-H1*H1*E0)/Det
       IF ((B*B).GT.(4*A*C)) THEN
          Disc=SQRT(B*B-4*A*C)
       ELSE
          DISC=0
       ENDIF
       IF (B.LT.0) Disc=-Disc
       Z=-2*C/(B+Disc)
       P3=P2+Z
       IF (ABS(P3-P1).LT.ABS(P3-P0)) THEN
          U=P1;P1=P0;P0=U;V=Y1;Y1=Y0;Y0=V
       ENDIF
       IF (ABS(P3-P2).LT.ABS(P3-P1)) THEN
          U=P2;P2=P1;P1=U;V=Y2;Y2=Y1;Y1=V
       ENDIF
       P2=P3
       Y2=F(P2)
       RelErr=ABS(Z)/(ABS(P2)+Small)
       IF ((RelErr.LT.Delta).AND.(ABS(Y2).LT.Epsilon)) Cond=1
       !IF (RelErr.LT.Delta) Cond=1
       !IF (ABS(Y2).LT.Epsilon) Cond=2
       !print *, "P3 = ", P3
       K=K+1
      !print *,'RelErr',RelErr,'Y2',Y2
    END DO
    RETURN
  END SUBROUTINE MULLER

  !**************************************************************************
  ! NUMERICAL METHODS: FORTRAN Programs, (c) John H. Mathews 1994
  ! To accompany the text:
  ! NUMERICAL METHODS for Mathematics, Science and Engineering, 2nd Ed, 1992
  ! Prentice Hall, Englewood Cliffs, New Jersey, 07632, U.S.A.
  ! This free software is complements of the author.
  !
  ! Algorithm 2.6 (Secant Method).
  ! Section 2.4, Newton-Raphson and Secant Methods, Page 85
  !**************************************************************************
  SUBROUTINE SECANT(F,P0,P1,Delta,Epsilon,Max,P2,Dp,Cond,K)
    
    IMPLICIT NONE
    
    REAL*8, PARAMETER :: Small=1E-10
    INTEGER :: Max, Cond, K
    REAL*8 :: Delta,Epsilon,Df,Dp,P0,P1,P2,Y0,Y1,Y2,RelErr
    REAL*8, EXTERNAL :: F
    
    K=0
    Cond=0
    Y0=F(P0)
    Y1=F(P1)
    DO WHILE ((K.LT.Max).AND.(Cond.EQ.0))
       Df=(Y1-Y0)/(P1-P0)
       IF (Df.EQ.0) THEN
          Cond=1
          Dp=P1-P0
          P2=P1
       ELSE
          Dp=Y1/Df
          P2=P1-Dp
       ENDIF
       Y2=F(P2)
       RelErr=ABS(Dp)/(ABS(P2)+Small)
       !IF (RelErr.LT.Delta) Cond=1
       !IF (ABS(Y2).LT.Epsilon) Cond=2
       IF ((RelErr.LT.Delta).AND.(ABS(Y2).LT.Epsilon)) Cond=1
       !IF ((Cond.EQ.2).AND.(ABS(Y2).LT.Epsilon)) Cond=3
       P0=P1
       P1=P2
       Y0=Y1
       Y1=Y2
       K=K+1
    END DO
    RETURN
  END SUBROUTINE SECANT
  

  SUBROUTINE SECANT_OMP(F,P0,P1,Delta,Epsilon,Max,P2,Dp,Cond,K)
    
    use omp_lib
!   USE IFLPORT
    IMPLICIT NONE
    
    REAL*8, PARAMETER :: Small=1E-10
    INTEGER :: Max, Cond, K
    REAL*8 :: Delta,Epsilon,Df,Dp,P0,P1,P2,Y0,Y1,Y2,RelErr
    REAL*8, EXTERNAL :: F
    
    K=0
    Cond=0
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS
!$OMP SECTION
    Y0=F(P0)
!$OMP SECTION
    Y1=F(P1)
!$OMP END SECTIONS NOWAIT
!$OMP END PARALLEL
    DO WHILE ((K.LT.Max).AND.(Cond.EQ.0))
       Df=(Y1-Y0)/(P1-P0)
       IF (Df.EQ.0) THEN
          Cond=1
          Dp=P1-P0
          P2=P1
       ELSE
          Dp=Y1/Df
          P2=P1-Dp
       ENDIF
       Y2=F(P2)
       RelErr=ABS(Dp)/(ABS(P2)+Small)
       IF (RelErr.LT.Delta) Cond=2
       IF (ABS(Y2).LT.Epsilon) Cond=3
       IF ((Cond.EQ.2).AND.(ABS(Y2).LT.Epsilon)) Cond=4
       P0=P1
       P1=P2
       Y0=Y1
       Y1=Y2
       K=K+1
    END DO
    RETURN
  END SUBROUTINE SECANT_OMP
  

  SUBROUTINE MULLER_OMP(F,P0,P1,P2,Delta,Epsilon,Max,P3,Z,K,Cond)
    use omp_lib
!   USE IFLPORT
    IMPLICIT NONE
    INTEGER, INTENT(in) :: Max
    INTEGER, INTENT(out) :: K,Cond
    REAL*8 :: A,B,C,Det,Disc,H0,H1,E0,E1
    REAL*8 :: Y0,Y1,Y2,RelErr,U,V
    REAL*8, INTENT(in) :: Delta,Epsilon
    REAL*8, INTENT(out) :: P3,Z
    REAL*8, INTENT(inout) :: P0,P1,P2
    REAL*8, PARAMETER :: Small=1E-10
    REAL*8, EXTERNAL :: F
    K=0
    Cond=0
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SECTIONS
!$OMP SECTION
    Y0=F(P0)
!$OMP SECTION
    Y1=F(P1)
!$OMP SECTION
    Y2=F(P2)
!$OMP END SECTIONS NOWAIT
!$OMP END PARALLEL
    DO WHILE ((K.LT.Max).AND.(Cond.EQ.0))
       !print *, "k =", k
       !print *, " P0 = ", P0, Y0
       !print *, " P1 = ", P1, Y1
       !print *, " P2 = ", P2, Y2
       H0=P0-P2
       H1=P1-P2
       C=Y2
       E0=Y0-C
       E1=Y1-C
       Det=H0*H1*(H0-H1)
       A=(E0*H1-H0*E1)/Det
       B=(H0*H0*E1-H1*H1*E0)/Det
       IF ((B*B).GT.(4*A*C)) THEN
          Disc=SQRT(B*B-4*A*C)
       ELSE
          DISC=0
       ENDIF
       IF (B.LT.0) Disc=-Disc
       Z=-2*C/(B+Disc)
       P3=P2+Z
       IF (ABS(P3-P1).LT.ABS(P3-P0)) THEN
          U=P1;P1=P0;P0=U;V=Y1;Y1=Y0;Y0=V
       ENDIF
       IF (ABS(P3-P2).LT.ABS(P3-P1)) THEN
          U=P2;P2=P1;P1=U;V=Y2;Y2=Y1;Y1=V
       ENDIF
       P2=P3
       Y2=F(P2)
       RelErr=ABS(Z)/(ABS(P2)+Small)
       !IF ((RelErr.LT.Delta).AND.(ABS(Y2).LT.Epsilon)) Cond=1
       IF ((RelErr.LT.Delta).OR.(ABS(Y2).LT.Epsilon)) Cond=1
      !print *, "P3 = ", P3
       K=K+1
    END DO
    RETURN
  END SUBROUTINE MULLER_OMP


END MODULE numeric
