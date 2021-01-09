!*********************************************************!
!*********************  ANT.G-2.5.0  *********************!
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
