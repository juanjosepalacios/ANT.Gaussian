!***************************************
!*                                     *
!*  ANT.G-2.5.0  - SpinOrbit.f90       *
!*                                     *
!*  Calculation of Spin-orbit coupling *
!*                                     *
!***************************************
!*                                     *
!*  Copyright (c) 2006-2018 by         *
!*                                     *
!*  Dr. David Jacob                    *
!*                                     *
!*  MPI fuer Mikrostrukturphysik       *
!*  Weinberg 2                         *
!*  06108 Halle                        *
!*  Germany                            *
!*				       *
!*  Wynand Dednam and                  *
!*  Prof. Andre E. Botha               *
!*                                     *
!*  University of South Africa         *
!*  28 Pioneer Ave                     *
!*  Florida Park                       *
!*  Roodepoort                         *
!*  1709                               *
!*  South Africa                       *
!*                                     *
!*  Sahar Pakdel and                   *
!*  Prof. Juan Jose Palacios Burgos    *
!*                                     *
!*  Departamento de Física             *
!*  de la Materia Condensada,          *
!*  Universidad Autónoma de Madrid     *
!*  28049 Madrid                       *
!*  Spain                              *
!*                                     *
!***************************************


!*********************************************************
!*                                                       *
!*  Module for calculation of Spin-orbit coupling based  *
!*  on subroutine by Joaquin Fernandez-Rossier           *
!*                                                       *
!*********************************************************
MODULE SpinOrbit
  IMPLICIT NONE

  PRIVATE
  
  REAL*8, DIMENSION(100000,2) :: matrix 
  COMMON /array/ matrix 
  
  PUBLIC :: CompHSO 
  
CONTAINS  
  
  !**************************************************************
  !*** Pauli matrices in basis of alpha- and beta-spin states ***
  !*** quantization axis z of |alpha>, |beta> parallel to     ***
  !*** global z axis                                          ***
  !**************************************************************
  SUBROUTINE PAULI_MATRICES(theta,phi,sigma_z,sigma_p,sigma_m,u,v,ustar,vstar)
    USE constants, ONLY : d_pi, ui
    IMPLICIT NONE

    !Input: direction of magnetization
    REAL*8, INTENT(in) :: theta, phi

    !Output: Pauli matrices 
    COMPLEX*16, DIMENSION(2,2), INTENT(out) :: sigma_z, sigma_p, sigma_m
    
    !Half of the magnetization direction angles
    REAL*8 :: halfphi, halftheta

    !Complex phases defining transformation matrix
    COMPLEX*16, INTENT(out) :: u, v, ustar,vstar

    INTEGER :: i,j
    
    sigma_z = 0.0d0
    sigma_p = 0.0d0
    sigma_m = 0.0d0

    u = dcmplx(0.0d0,0.0d0)
    v = dcmplx(0.0d0,0.0d0)
    ustar = dcmplx(0.0d0,0.0d0)
    vstar = dcmplx(0.0d0,0.0d0)
    
    halfphi=0.5d0*phi*d_pi/180.0d0
    halftheta=0.5d0*theta*d_pi/180.0d0
    
    u=dcos(halftheta)*dcmplx(dcos(halfphi),dsin(halfphi))
    v=dsin(halftheta)*dcmplx(dcos(halfphi),-dsin(halfphi))
    
    ustar=dconjg(u)                                                  
    vstar=dconjg(v)						
    
    sigma_z(1,1)=1.d0      !  u*ustar-v*vstar         
    sigma_z(1,2)=0.d0      !  -2.d0*u*v               
    sigma_z(2,1)=0.d0      !  -2.d0*ustar*vstar       
    sigma_z(2,2)=-1.d0     !  v*vstar-u*ustar         
                                                  
    sigma_p(1,1)=0.d0      !  2.d0*u*vstar        
    sigma_p(1,2)=2.d0      !  2.d0*u*u            
    sigma_p(2,1)=0.d0      !  -2.d0*vstar*vstar   
    sigma_p(2,2)=0.d0      !  -2.d0*vstar*u       
                                                  
    sigma_m(1,1)=0.d0      !  2.d0*v*ustar                        
    sigma_m(1,2)=0.d0      !  -2.d0*v*v                           
    sigma_m(2,1)=2.d0      !  2.d0*ustar*ustar                    
    sigma_m(2,2)=0.d0      !  -2.d0*ustar*v                       
    
    !PRINT *, " sigma_z = "
    !DO i=1,2
    !   PRINT '(2(F11.5))', ( REAL(sigma_z(i,j)), j=1,2 )
    !END DO
    !PRINT *, " sigma_p = "
    !DO i=1,2
    !   PRINT '(2(F11.5))', ( REAL(sigma_p(i,j)), j=1,2 )
    !END DO
    !PRINT *, " sigma_m = "
    !DO i=1,2
    !   PRINT '(2(F11.5))', ( REAL(sigma_m(i,j)), j=1,2 )
    !END DO

  END SUBROUTINE PAULI_MATRICES

  !*************************************************************************************
  !*** Angular momentum matrices for total angular momentum L = 0,1,2,... **************
  !*** Assumed basis: Lz eigen states, Order for L >= 2: m = 0, +1, -1, +2, -2, ...  ***
  !*************************************************************************************
  SUBROUTINE L_MATRICES( L, Lz, Lp, Lm )
    USE constants
    IMPLICIT NONE
    
    ! Input: total angular momentum number L = 0,1,2,...
    INTEGER, INTENT(in) :: L

    !Ouput: Angular momentum matrices
    COMPLEX*16, DIMENSION(2*L+1,2*L+1), INTENT(out) :: Lz,Lp,Lm
    COMPLEX*16, DIMENSION(2*L+1,2*L+1) :: temp
    COMPLEX*16 :: Us, c0, c1, c2
    COMPLEX*16, dimension(3,3) :: Up
    COMPLEX*16, dimension(5,5) :: Ud
    COMPLEX*16, dimension(7,7) :: Uf
    
    INTEGER  :: n1, n2, i, j
    REAL*8 :: m, m1, m2
    
    Lp = 0.0d0
    Lm = 0.0d0
    Lz = 0.0d0
    
    IF (L > 1) THEN
        DO n1=1,2*L+1
           ! m1 = 0, +1, -1, +2, -2, ... for L = 2, 3, ... shells in Gaussian
           m1 = n1 / 2 * ( 1 - 2*MOD(n1,2) )
           
           !<m1|Lz|m2> matrix is diagonal :
           Lz(n1,n1) = m1
            
           DO n2=1,2*L+1
              ! m2 = 0, +1, -1, +2, -2, ...
               m2 = n2 / 2 * ( 1 - 2*MOD(n2,2) )
              
              !< m1 | Lp | m2 > =
              IF( m1 == m2 + 1 ) THEN
                 Lp(n1,n2) = CMPLX( SQRT( DBLE(L*(L+1)-(m2*(m2+1))) ), 0.0d0 )
              !< m1 | Lm | m2 >=   
              ELSE IF( m1 == m2 - 1 ) THEN
                 Lm(n1,n2) = CMPLX( SQRT(DBLE(L*(L+1)-(m2*(m2-1)))), 0.0d0 )              
              ELSE
                 Lp(n1,n2) = c_zero
                 Lm(n1,n2) = c_zero
              END IF
              
           END DO       
        END DO
    ELSE IF (L == 1) THEN
         i = 0
         j = 0    
         DO n1=2*L+1,1,-1
           ! m1 = +1, -1, 0 for P shells in Gaussian
           m1 = n1 / 2 * ( 2*MOD(n1,2) - 1 )
           
           i = i + 1
           
           !<m1|Lz|m2> matrix is diagonal :
           Lz(i,i) = m1
            
           DO n2=2*L+1,1,-1
              ! m1 = +1, -1, 0 for P shells in Gaussian
               m2 = n2 / 2 * ( 2*MOD(n2,2) - 1 )
               
               j = j + 1
               
              !< m1 | Lp | m2 > =
              IF( m1 == m2 + 1 ) THEN
                 Lp(i,j) = CMPLX( SQRT( DBLE(L*(L+1)-(m2*(m2+1))) ), 0.0d0 )
              !< m1 | Lm | m2 >=   
              ELSE IF( m1 == m2 - 1 ) THEN
                 Lm(i,j) = CMPLX( SQRT(DBLE(L*(L+1)-(m2*(m2-1)))), 0.0d0 )              
              ELSE
                 Lp(i,j) = c_zero
                 Lm(i,j) = c_zero
              END IF
              !PRINT *, i, j
              IF (j == 3) THEN
                  j = 0
              END IF    
           END DO       
        END DO        
    END IF

    PRINT *, " Lz in spherical harmonic basis = "
    DO n1=1,2*L+1
       PRINT '(100(F11.5))', ( REAL(Lz(n1,n2)), n2=1,2*L+1 )  
    END DO
    PRINT *, " Lp in spherical harmonic basis = "
    DO n1=1,2*L+1
       PRINT '(100(F11.5))', ( REAL(Lp(n1,n2)), n2=1,2*L+1 )
    END DO
    PRINT *, " Lm in spherical harmonic basis = "
    DO n1=1,2*L+1
       PRINT '(100(F11.5))', ( REAL(Lm(n1,n2)), n2=1,2*L+1 )
    END DO    
    
    
    !
    ! Transformation matrices for transformation from spherical harmonic orbitals                               
    ! to cubic harmonic orbitals                                                                          
    !                                                                                                           
    ! Cubic harmonics are assumed to be ordered as in CRYSTAL:                                   
    ! l=0: s                                                                                     
    ! l=1: x, y, z                                                                                          
    ! l=2: z^2, xz, yz, x^2-y^2, xy                                                                    
    ! l=3: z^3, xz^2, yz^2, z(x^2-y^2), xyz, x(x^2-3y^2), y(3x^2-y^2)                                                                     
    !                                                                                           
    
    Us = c_one
    c0 = c_zero                                                     
    c1 = c_one
    c2 = 1.0/sqrt(2.0)
    Up = TRANSPOSE(reshape((/ -c2,     c2,    c0, &        
                    -c2/ui, -c2/ui, c0, &                  
                     c0,     c0,    c1  /), (/ 3, 3 /) ))  
    Ud = TRANSPOSE(reshape((/  c1,     c0,     c0,     c0,     c0,    &                                   
              c0,    -c2,     c2,     c0,     c0,    &                                   
              c0,    -c2/ui,  -c2/ui, c0,     c0,    &                                   
              c0,     c0,     c0,     c2,     c2,    &                                    
              c0,     c0,     c0,     c2/ui, -c2/ui  /), (/ 5, 5 /) ))                    
    Uf = TRANSPOSE(reshape((/  c1,     c0,     c0,     c0,     c0,     c0,     c0,   &         
              c0,    -c2,     c2,     c0,     c0,    c0,     c0,    &                
              c0,    -c2/ui,  -c2/ui, c0,     c0,    c0,     c0,    &                
              c0,     c0,     c0,     c2,     c2,    c0,     c0,    &                
              c0,     c0,     c0,     c2/ui, -c2/ui,  c0,     c0,    &
              c0,     c0,     c0,     c0,     c0,   -c2,     c2,    &
              c0,     c0,     c0,     c0,     c0,  -c2/ui,  -c2/ui  /), (/ 7, 7 /) )) 
              
    !
    ! transform matrices to cubic harmonic basis
    !    
    
    IF( L == 1 ) THEN
      PRINT *, " Up = "
      DO n1=1,2*L+1
         PRINT '(100(F11.5))', ( REAL(Up(n1,n2)), n2=1,2*L+1 )
      END DO      
      temp = MATMUL( Lz, CONJG(TRANSPOSE(Up)))
      Lz = MATMUL( Up, temp )
      
      temp = MATMUL( Lp, CONJG(TRANSPOSE(Up)))
      Lp = MATMUL( Up, temp )
      
      temp = MATMUL( Lm, CONJG(TRANSPOSE(Up)))
      Lm = MATMUL( Up, temp )
    ELSE IF (L == 2) THEN         
      PRINT *, " Ud = "
      DO n1=1,2*L+1
         PRINT '(100(F11.5))', ( REAL(Ud(n1,n2)), n2=1,2*L+1 ) 
      END DO         
      temp = MATMUL( Lz, CONJG(TRANSPOSE(Ud)))    
      Lz = MATMUL( Ud, temp )                     
                                                  
      temp = MATMUL( Lp, CONJG(TRANSPOSE(Ud)))    
      Lp = MATMUL( Ud, temp )                     
                                                  
      temp = MATMUL( Lm, CONJG(TRANSPOSE(Ud)))    
      Lm = MATMUL( Ud, temp )                     
    ELSE IF (L == 3) THEN         
      PRINT *, " Uf = "
      DO n1=1,2*L+1
         PRINT '(100(F11.5))', ( REAL(Uf(n1,n2)), n2=1,2*L+1 ) 
      END DO         
      temp = MATMUL( Lz, CONJG(TRANSPOSE(Uf)))    
      Lz = MATMUL( Uf, temp )                     
                                                  
      temp = MATMUL( Lp, CONJG(TRANSPOSE(Uf)))    
      Lp = MATMUL( Uf, temp )                     
                                                  
      temp = MATMUL( Lm, CONJG(TRANSPOSE(Uf)))    
      Lm = MATMUL( Uf, temp )                              
    END IF    
    
    
    PRINT *, " Lz in cartesian basis (real part) = "
    DO n1=1,2*L+1
       PRINT '(100(F11.5))', ( REAL(Lz(n1,n2)), n2=1,2*L+1 )
    END DO
    PRINT *, " Lz in cartesian basis (imaginary part) = "
    DO n1=1,2*L+1
       PRINT '(100(F11.5))', ( AIMAG(Lz(n1,n2)), n2=1,2*L+1 )
    END DO
    PRINT *, " Lp in cartesian basis (real part) = "
    DO n1=1,2*L+1
       PRINT '(100(F11.5))', ( REAL(Lp(n1,n2)), n2=1,2*L+1 )
    END DO
    PRINT *, " Lp in cartesian basis (imaginary part) = "
    DO n1=1,2*L+1
       PRINT '(100(F11.5))', ( AIMAG(Lp(n1,n2)), n2=1,2*L+1 )
    END DO
    PRINT *, " Lm in cartesian basis (real part) = "
    DO n1=1,2*L+1
       PRINT '(100(F11.5))', ( REAL(Lm(n1,n2)), n2=1,2*L+1 )
    END DO
    PRINT *, " Lm in cartesian basis (imaginary part) = "
    DO n1=1,2*L+1
       PRINT '(100(F11.5))', ( AIMAG(Lm(n1,n2)), n2=1,2*L+1 )
    END DO


  END SUBROUTINE L_MATRICES 
  
  !********************************************************************************************************************    
   
  
  FUNCTION R(x,ll,idx1,idx2)  ! For calculating N_0, the normalization constant of the CGTOs in Eq. (3) of Beilstein J. Nanotechnol. 2018, 9, 1015-1023
    IMPLICIT NONE
    
      REAL*8, INTENT(IN) :: x
      REAL*8 :: R 
      INTEGER :: j,ll,idx1,idx2
      
      R=0.0
      DO j=idx1,idx2 ! Runs over all the primitives in a given shell of a basis set
          R=R+matrix(j,2)*x**ll*exp(-matrix(j,1)*x**2) ! Here x^ll-1 with ll=l+1 depends on the shell type l = 0, 1, 2 !! 
      END DO
      R = R**2*x**2    

  END FUNCTION R  
    
  

  !********************************************************************************************************************    
 
  FUNCTION Rr(x,ll,idx1,idx2)  ! Corresponds to R_i(r) or R_j(r) in Eq. (5) of Beilstein J. Nanotechnol. 2018, 9, 1015-1023
    IMPLICIT NONE
    
      REAL*8, INTENT(IN) :: x
      REAL*8 :: Rr 
      INTEGER :: j,ll,idx1,idx2
      
      Rr=0.0
      DO j=idx1,idx2  ! Runs over all the primitives in a given shell of a basis set
          Rr=Rr+matrix(j,2)*x**ll*exp(-matrix(j,1)*x**2) ! Here x^ll-1 with ll=l+1 depends on the shell type l = 0, 1, 2 !!  
      END DO
          
  
  END FUNCTION Rr  
  
  !********************************************************************************************************************
     
! ---------------------------------------------------------
!    Subroutine to evaluate integrals using
!    10 point Newton-Cotes rule. See the        
!    Handbook of mathematical formulas  
!    and integrals by Alan Jeffrey      
! --------------------------------------------------------
  
  SUBROUTINE INTEGRATE(func,A,B,ll,result,idx1,idx2)  ! Calculates the normalizations N_0 over CGTOs in Eq. (3) of  Beilstein J. Nanotechnol. 2018, 9, 1015-1023
      IMPLICIT REAL*8 (A-H,O-Y)
      IMPLICIT INTEGER (I,N)
      DIMENSION XP(9999), WP(9999)
      INTEGER :: idx1,idx2,ll
      REAL*8 :: result, SUM1
     
      NQ = 730   ! NQ-1 MUST BE A MULTIPLE OF 9!!   
      NL = 1
      SUM1  = 0.0D0
      dx = (B-A)/DFLOAT(NQ)
      DO N = 1,NQ
      SUM1 = SUM1 + dx
      XP(N) = SUM1
      END DO

      !CALL SIMS(XP,NL,NQ,WP)
      CALL BODE10PT(A,B,XP,WP,NQ)

      result = 0.0d0
      DO 12 N = 1,NQ
       result = result + WP(N)*func(XP(N),ll,idx1,idx2)
 12   CONTINUE
         
      RETURN
      END SUBROUTINE INTEGRATE

 !********************************************************************************************************************       
! ----------------------------------------------------------      
!    Subroutine to evaluate integrals using   
!    10 point Newton-Cotes rule. See the      
!    Handbook of mathematical formulas        
!    and integrals by Alan Jeffrey            
! --------------------------------------------------------
  
  SUBROUTINE INTEGRATE1(ii,jj,A,B,result,idx11,idx12,idx21,idx22) ! Performs the integration in Eq. (5) of  Beilstein J. Nanotechnol. 2018, 9, 1015-1023
      IMPLICIT REAL*8 (A-H,O-Y)
      IMPLICIT INTEGER (I,N)
      DIMENSION XP(9999), WP(9999)
      INTEGER :: ii,jj,idx11,idx12,idx21,idx22
      REAL*8 :: result, SUM1, NormOrb1, NormOrb2
     
      INTERFACE AFunc
       FUNCTION func (y,ll,i1,i2)
          REAL*8 :: func
          REAL*8, INTENT(in) ::y
          INTEGER :: i1,i2,ll      
       END FUNCTION func                         
      END INTERFACE AFunc     
      
      procedure (func), pointer :: orb1 => null ()  
      procedure (func), pointer :: orb2 => null ()     
      
      orb1 => Rr
      CALL integrate(R,A,B,ii,NormOrb1,idx11,idx12)
      
      orb2 => Rr
      CALL integrate(R,A,B,jj,NormOrb2,idx21,idx22)
     
      NQ = 730 ! NQ-1 MUST BE A MULTIPLE OF 9!!
      NL = 1
      SUM1  = 0.0D0
      dx = (B-A)/DFLOAT(NQ)
      DO N = 1,NQ
      SUM1 = SUM1 + dx
      XP(N) = SUM1
      END DO

      !CALL SIMS(XP,NL,NQ,WP)
      CALL BODE10PT(A,B,XP,WP,NQ)

      result = 0.0d0
      DO 12 N = 1,NQ
         result = result + WP(N)*(1.0/XP(N))*orb1(XP(N),ii,idx11,idx12)*orb2(XP(N),jj,idx21,idx22)/sqrt(NormOrb1*NormOrb2) ! Eq. (5) Beilstein J. Nanotechnol. 2018, 9, 1015-1023
 12   CONTINUE
         
      RETURN
      END SUBROUTINE INTEGRATE1
      
! -------------------------------------------------
      SUBROUTINE BODE10PT(A,B,X,W,N)
! -------------------------------------------------
!     PREPARE POINTS X(N) AND WEIGHTS W(N) FOR
!     INTEGRATION BY 10 POINT NEWTON-COTES RULE.
!     SUBROUTINE WORKS WHEN N-1 IS A MULTIPLE OF 9
! -------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER (I,K,N)
      DIMENSION  W(N),X(N)

      IF ((N .LT. 10) .OR. (MOD((N-1),9) .NE. 0)) THEN
      WRITE(*,*) ' ERROR IN SUBROUTINE SIMPSON. '
      WRITE(*,*) ' N-1 MUST BE A MULTIPLE OF 9 '
      STOP
      END IF
      W(1) = 0.0D0
      X(1) = A
      X(N) = B
      H  = (B-A)/DFLOAT(N-1)
      H1 = 9.0D0*H/89600.0D0
      DO K = 2,N-1
       X(K) = X(K-1) + H
      END DO
      DO K= 1,N-1,9
      W(K)   =  2857.0D0*H1 + W(K)
      W(K+1) = 15741.0D0*H1
      W(K+2) =  1080.0D0*H1
      W(K+3) = 19344.0D0*H1
      W(K+4) =  5778.0D0*H1
      W(K+5) =  5778.0D0*H1
      W(K+6) = 19344.0D0*H1
      W(K+7) =  1080.0D0*H1
      W(K+8) = 15741.0D0*H1
      W(K+9) =  2857.0D0*H1
      END DO
      END SUBROUTINE                  

!! ===================================
!      SUBROUTINE SIMS(XP,NL,NR,T)
!! ===================================
!      IMPLICIT REAL*8(A-H,O-Z)
!      INTEGER, intent(IN) :: NL,NR
!      INTEGER :: NJ, NK, K, N9, M, N4, J
!      REAL*8, DIMENSION(NR) ::  T,XP
!
!      T(NL)=0.0D0
!      H=XP(NR)-XP(NR-1)
!      H1=9.D0*H/89600.D0                                              
!      NK=NL                                                           
!      NJ=NL                                                           
!      IF((NR-NL).LT.9)GO TO 25                                        
!      DO 10 K=NL,NR,9                                                 
!      N9=K+9                                                          
!      IF(N9.GT.NR)GO TO 25                                            
!      NK=N9                                                                                                                                     
!      T(K)   =  2857.D0*H1+T(K)
!      T(K+1) = 15741.D0*H1
!      T(K+2) =  1080.D0*H1
!      T(K+3) = 19344.D0*H1
!      T(K+4) =  5778.D0*H1
!      T(K+5) =  5778.D0*H1
!      T(K+6) = 19344.D0*H1
!      T(K+7) =  1080.D0*H1
!      T(K+8) = 15741.D0*H1
!      T(K+9) =  2857.D0*H1
! 10   CONTINUE
! 25   CONTINUE
!      H2 = 3.0D0/8.0D0*H
!      NJ = NK
!      IF(( NR-NK).LT.3)GO TO 35
!      DO 30 M = NK,NR,3
!      N4=M+3
!      IF(N4.GT.NR)GO TO 35
!      NJ=N4
!      T(M)   = T(M)+H2
!      T(M+1) = 3.0D0*H2
!      T(M+2) = 3.0D0*H2
! 30   T(M+3) = H2
! 35   CONTINUE
!      IF(NJ.EQ.NR)GO TO 60
!      DO 50 J=NJ,NR-1
!      T(J+1)= H/2.0D0
! 50   T(J)  = H/2.0D0 + T(J)
! 60   CONTINUE
!      RETURN
!      END SUBROUTINE SIMS   
  
  !**************************************************************
  !*** Compute matrix of SO Hamiltonian for a given basis set ***
  !**************************************************************
  SUBROUTINE CompHSO(hamil_SO,HD,hamil,SD,overlap_SO,NAOs,Nshell)
    USE parameters, ONLY: theta, phi, soc_cff_p, soc_cff_d, soc_cff_f, NSpinRot, SpinRotTheta, SpinRotPhi, socfac, NSOCFacAtom, SOCFacAtom, NSOCEdit, SOCEditP, SOCEditD, SOCEditF
    USE G09common, ONLY : GetNAtoms, GetShellT, GetShellC, GetAtm4Sh, GetShellN, GetShellA, GetShlADF, GetEXX, GetC1, GetC2, GetC3, GetC4, GetAN, GetAtmCo
    USE cluster, ONLY : LoAOrbNo, HiAOrbNo
    USE constants
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NAOs, NShell
    INTEGER :: L1, L2, L3, AtomID1, AtomID2, ShellT1, ShellT2, ShellC1, ShellC2, ShellNPrim1, ShellAindex1, ShlADFindex1, ShellNPrim2, ShellAindex2, ShlADFindex2
    INTEGER, DIMENSION(NAOs) :: AOT
    REAL*8, DIMENSION(NAOs) :: lambda

    !*********************************
    !Output: Atomic LS coupling matrix 
    !*********************************
    COMPLEX*16, DIMENSION(2,NAOs,2,NAOs) :: SROT, HROT, HSO  ! Need to reshape HSO to H_SOC(i,j) where i = 1, 2*Norb, j = 1, 2*Norb
    COMPLEX*16, DIMENSION(2*NAOs,2*NAOs), INTENT(OUT) :: hamil_SO
    COMPLEX*16, DIMENSION(2*NAOs,2*NAOs), INTENT(INOUT) :: hamil, overlap_SO
    
    REAL*8, DIMENSION(2,NAOs,NAOs), INTENT(IN) :: HD
    REAL*8, DIMENSION(NAOs,NAOs), INTENT(IN) :: SD

    INTEGER :: i, j, k, q, s1, s2, ish1, ish2, ispin , jspin, Z, acount
    REAL*8 :: theta_atom, phi_atom, result, A, B, x, zz, socfac_atom, soc_cff_p_atom, soc_cff_d_atom, soc_cff_f_atom

    COMPLEX*16 :: u, v, ustar, vstar, HROT1temp11, HROT1temp12, HROT1temp21, HROT1temp22, SROT1temp11, SROT1temp12, SROT1temp21, SROT1temp22
    COMPLEX*16 :: HROT2temp11, HROT2temp12, HROT2temp21, HROT2temp22, SROT2temp11, SROT2temp12, SROT2temp21, SROT2temp22
    COMPLEX*16, DIMENSION(2,2) :: sigma_z, sigma_p, sigma_m
    COMPLEX*16, DIMENSION(3,3) :: L_z1, L_p1, L_m1
    COMPLEX*16, DIMENSION(5,5) :: L_z2, L_p2, L_m2
    COMPLEX*16, DIMENSION(7,7) :: L_z3, L_p3, L_m3

    COMPLEX*16, DIMENSION(2,3,2,3) :: LS1
    COMPLEX*16, DIMENSION(2,5,2,5) :: LS2
    COMPLEX*16, DIMENSION(2,7,2,7) :: LS3    
    
    !**************************************************
    !Output: Spin-orbit constants X_ij within each atom 
    !**************************************************    
    REAL*8, DIMENSION(NShell,NShell) :: Xi
     
    REAL*8, PARAMETER :: e=1.602177e-19 !Coulomb
    REAL*8, PARAMETER :: me=9.1093897e-31 !kg
    REAL*8, PARAMETER :: c=2.99792458e8 !m/s
    REAL*8, PARAMETER :: hbar=1.054571800e-34 !J.s
    REAL*8, PARAMETER :: Ke=8.987552e9 ! 1/(4*pi*epsilon0)
    REAL*8, PARAMETER :: a0 = 5.29e-11 ! bohr radius    
    
    L1 = 1
    L2 = 2
    L3 = 3
    
    ! Lz, L+, L- Matrices
    PRINT *, "Calculating L matrices for p orbitals... "
    CALL L_Matrices( L1, L_z1, L_p1, L_m1 ) ! p-orbitals
    PRINT *, "Calculating L matrices for d orbitals... "
    CALL L_Matrices( L2, L_z2, L_p2, L_m2 ) ! d-orbitals
    PRINT *, "Calculating L matrices for f orbitals... "
    CALL L_Matrices( L3, L_z3, L_p3, L_m3 ) ! f-orbitals    
    
    PRINT *, "Total number of shells in the cluster:"
    PRINT *, NShell
      
    ish1 = 1
    
    DO i=1,NShell
    	  ShellT1 = GetShellT(i) 
    		DO k=1,2*ShellT1+1
    			AOT(ish1) = ShellT1
    			ish1 = ish1 + 1
            END DO
    END DO        
    
    PRINT *, "Shell type of each orbital:"
    DO i=1,NAOs
      PRINT *, AOT(i)
    END DO  
      
    PRINT *, "Basis set of each atom extracted from internal Gaussian variables:"
    DO i=1,NShell
       ShellT1 = GetShellT(i)
       ShellC1 = GetShellC(i)
       ShellNPrim1 = GetShellN(i)
       ShellAindex1 = GetShellA(i)
       ShlADFindex1 = GetShlADF(i)
       AtomID1 = GetAtm4Sh(i)
       !IF (i > 1 .and. AtomID1 > GetAtm4Sh(i-1)) THEN
       !   acount = acount + 1
       !END IF       
       DO j=1,ShellNPrim1
          !print *, ShellNPrim1 
          IF ( ShellT1 == 0) THEN
              matrix(ShellAindex1+j-1,1)=GetEXX(ShellAindex1+j-1)
              matrix(ShellAindex1+j-1,2)=GetC1(ShellAindex1+j-1)   
              PRINT *,matrix(ShellAindex1+j-1,1),matrix(ShellAindex1+j-1,2)              
          ELSE IF ( ShellT1 == 1 .and. (ShellC1 == 0 .or. ShellC1 == 1)) THEN
              matrix(ShellAindex1+j-1,1)=GetEXX(ShellAindex1+j-1)
              matrix(ShellAindex1+j-1,2)=GetC2(ShellAindex1+j-1)
              PRINT *,matrix(ShellAindex1+j-1,1),matrix(ShellAindex1+j-1,2)              
          ELSE IF (ShellT1 == 2 .and. ShellC1 == 2) THEN
              matrix(ShellAindex1+j-1,1)=GetEXX(ShellAindex1+j-1)                         
              matrix(ShellAindex1+j-1,2)=GetC3(ShlADFindex1+j-1)                          
              PRINT *,matrix(ShellAindex1+j-1,1),matrix(ShellAindex1+j-1,2)
          ELSE IF (ShellT1 == 3 .and. ShellC1 == 2) THEN       
              matrix(ShellAindex1+j-1,1)=GetEXX(ShellAindex1+j-1)                         
              matrix(ShellAindex1+j-1,2)=GetC4(ShlADFindex1+j-1)                          
              PRINT *,matrix(ShellAindex1+j-1,1),matrix(ShellAindex1+j-1,2)          
          END IF                 
       END DO                                                 
    END DO	 
    
    A=1.0e-12/a0
    B=1.0e-9/a0    
    
    zz = (1.0/a0**3)*(hbar**2)*(Ke/e)*((e**2)/(2.*(me*c)**2)) ! Divide by e to convert from Joules to eV     
                
    PRINT *, "Spin-orbit constants within each atom:"
    PRINT '(a10,2(a15),a20)', 'Atom No.','Shell type 1','Shell type 2','Spin-orbit constant'   
    DO i=1,NShell
        ShellT1 = GetShellT(i)
        ShellC1 = GetShellC(i)    
        ShellAindex1 = GetShellA(i)
        ShlADFindex1 = GetShlADF(i)
        ShellNPrim1 = GetShellN(i)
        AtomID1 = GetAtm4Sh(i)
        DO k=1,NShell
           ShellT2 = GetShellT(k)
           ShellC2 = GetShellC(k)          
           ShellAindex2 = GetShellA(k)
           ShlADFindex2 = GetShlADF(k)
           ShellNPrim2 = GetShellN(k)
           AtomID2 = GetAtm4Sh(k)
           IF (AtomID1 == AtomID2) THEN
              Z = GetAN(AtomID2)  ! Atomic number of atom AtomID2
              IF( NSOCFacAtom > 0) THEN  ! User-defined multiplicative SOC factor of atom AtomID2
                socfac_atom = SOCFacAtom(AtomID2)
              ELSE
                socfac_atom = socfac
              END IF              
              IF( NSOCEdit > 0) THEN  ! User-defined SOC coefficients of atom AtomID2
                soc_cff_p_atom = SOCEditP(AtomID2)
                soc_cff_d_atom = SOCEditD(AtomID2)
                soc_cff_f_atom = SOCEditF(AtomID2)
              ELSE
                soc_cff_p_atom = soc_cff_p
                soc_cff_d_atom = soc_cff_d
                soc_cff_f_atom = soc_cff_f
              END IF	               
              IF (socfac_atom > 0.0d0) THEN
                  CALL integrate1(ShellT1,ShellT2,A,B,result,ShellAindex1,ShellAindex1+ShellNPrim1-1,ShellAindex2,ShellAindex2+ShellNPrim2-1)                
                  Xi(i,k) = socfac_atom*zz*Z*result       
              ELSE IF (ShellT1 == 1 .and. ShellT2 == 1) THEN
                  Xi(i,k) = soc_cff_p_atom
              ELSE IF (ShellT1 == 2 .and. ShellT2 == 2) THEN                
                  Xi(i,k) = soc_cff_d_atom
              ELSE IF (ShellT1 == 3 .and. ShellT2 == 3) THEN 
                  Xi(i,k) = soc_cff_f_atom
              ELSE 
                  Xi(i,k) = 0.0 
              END IF
              PRINT '(3(i10),F25.10)',AtomID1,ShellT1,ShellT2,Xi(i,k)                                     
           END IF
        END DO                                                          
    END DO	               
       
!    PRINT *, "Spin-orbit couplings:"
!    DO i=1,NAOs
!      IF (AOT(i) == 0) THEN
!      	lambda(i) = 0.0d0
!      ELSE IF (AOT(i) == 1) THEN	
!        lambda(i) = soc_cff_p
!      ELSE IF (AOT(i) == 2) THEN
!        lambda(i) = soc_cff_d 
!      END IF  
!      PRINT *, lambda(i)
!    END DO
    
    CALL FLUSH(6)
    
    !*********************************************************************************************************************
    !Populate dummy matrix for rotation from local to global spin basis with elements from SELF-CONISTENT COLLINEAR matrix  
    !*********************************************************************************************************************  
    
    IF (NSpinRot > 0) then                                                                
    HROT=c_zero    
    SROT=c_zero                                                 
    
    do AtomID1=1,GetNAtoms()
       do AtomID2=1,GetNAtoms()     
           do i=LoAOrbNo(AtomID1),HiAOrbNo(AtomID1)
              do j=LoAOrbNo(AtomID2),HiAOrbNo(AtomID2)
                 HROT (1, i, 1, j) = dcmplx(HD(1,i,j),0.0d0)
                 HROT (2, i, 2, j) = dcmplx(HD(2,i,j),0.0d0)
                 SROT (1, i, 1, j) = dcmplx(SD(i,j),0.0d0)
                 SROT (2, i, 2, j) = dcmplx(SD(i,j),0.0d0)
              end do
           end do         
       end do
    end do
    
    end if
    

    !**************************************
    !Construct LS matrix for complete basis
    !**************************************
                
    HSO = c_zero
    LS1 = c_zero
    LS2 = c_zero
    LS3 = c_zero
    ish1 = 1
    ish2 = 1
    
    DO i = 1, NShell   
          ShellT1 = GetShellT(i) 
          ShellC1 = GetShellC(i)
          AtomID1 = GetAtm4Sh(i)
          DO j=1,2*ShellT1+1
             DO k = 1, NShell   
                ShellT2 = GetShellT(k) 
                ShellC2 = GetShellC(k)
                AtomID2 = GetAtm4Sh(k)   
          	DO q=1,2*ShellT2+1
          	  IF (AtomID2 == AtomID1) THEN 
                    IF( NSpinRot > 0) THEN
                      theta_atom = SpinRotTheta(AtomID2)
                      phi_atom = SpinRotPhi(AtomID2)
                    ELSE
                      theta_atom = theta
                      phi_atom = phi	
                    END IF                   	    
                    CALL Pauli_Matrices( theta_atom, phi_atom, sigma_z, sigma_p, sigma_m, u, v, ustar, vstar )                     	        
          	        DO s1 = 1,2        
          	            DO s2 = 1,2
                               IF( ShellT1 == 0 .and. ShellT2 == 0)THEN 
                                     HSO( s1, ish1, s2, ish2 ) = c_zero 
                               ELSE IF( ShellT1 == 1 .and. ShellT2 == 1 .and. (ShellC1 == 0 .or. ShellC1 == 1) .and. (ShellC2 == 0 .or. ShellC2 == 1)) THEN
                                     LS1( s1, j, s2, q ) &                       
                                          = 0.50d0 * L_z1(j,q) * sigma_z(s1,s2) &
                                          + 0.25d0 * L_p1(j,q) * sigma_m(s1,s2) &
                                          + 0.25d0 * L_m1(j,q) * sigma_p(s1,s2)  
                                     HSO( s1, ish1, s2, ish2) = Xi(i,k)*LS1(s1,j,s2,q)                              !lambda(ish1) * LS1(s1,j,s2,q)                                     
                               ELSE IF( ShellT1 == 2 .and. ShellT2 == 2 .and. ShellC1 == 2 .and. ShellC2 == 2) THEN
                                     LS2( s1, j, s2, q ) &                       
                                          = 0.50d0 * L_z2(j,q) * sigma_z(s1,s2) &
                                          + 0.25d0 * L_p2(j,q) * sigma_m(s1,s2) &
                                          + 0.25d0 * L_m2(j,q) * sigma_p(s1,s2)  
                                     HSO( s1, ish1, s2, ish2) = Xi(i,k)*LS2(s1,j,s2,q)                              !lambda(ish1) * LS2(s1,j,s2,q)
                               ELSE IF( ShellT1 == 3 .and. ShellT2 == 3 .and. ShellC1 == 2 .and. ShellC2 == 2) THEN
                                     LS3( s1, j, s2, q ) &                       
                                          = 0.50d0 * L_z3(j,q) * sigma_z(s1,s2) &
                                          + 0.25d0 * L_p3(j,q) * sigma_m(s1,s2) &
                                          + 0.25d0 * L_m3(j,q) * sigma_p(s1,s2)  
                                     HSO( s1, ish1, s2, ish2) = Xi(i,k)*LS3(s1,j,s2,q)                              !lambda(ish1) * LS3(s1,j,s2,q)                                
                               END IF                              
                        END DO
                    END DO
                    ! Rotate INTRA-ATOMIC SELF-CONSISTENT COLLINEAR and OVERLAP matrix elements from local to global spin space
                    IF (NSpinRot > 0) THEN 
                       !Hamiltonian                                                                
                       HROT (1, ish1, 1, ish2) =  u*ustar*dcmplx(HD(1,ish1,ish2),0.0d0)+v*vstar*dcmplx(HD(2,ish1,ish2),0.0d0) 
                       HROT (1, ish1, 2, ish2) =  v*ustar*dcmplx(HD(1,ish1,ish2),0.0d0)-v*ustar*dcmplx(HD(2,ish1,ish2),0.0d0) 
                       HROT (2, ish1, 1, ish2) =  u*vstar*dcmplx(HD(1,ish1,ish2),0.0d0)-u*vstar*dcmplx(HD(2,ish1,ish2),0.0d0) 
                       HROT (2, ish1, 2, ish2) =  v*vstar*dcmplx(HD(1,ish1,ish2),0.0d0)+u*ustar*dcmplx(HD(2,ish1,ish2),0.0d0)     
                       
                       !Overlap matrix
                       SROT (1, ish1, 1, ish2) =  u*ustar*dcmplx(SD(ish1,ish2),0.0d0)+v*vstar*dcmplx(SD(ish1,ish2),0.0d0)
                       SROT (1, ish1, 2, ish2) =  v*ustar*dcmplx(SD(ish1,ish2),0.0d0)-v*ustar*dcmplx(SD(ish1,ish2),0.0d0)
                       SROT (2, ish1, 1, ish2) =  u*vstar*dcmplx(SD(ish1,ish2),0.0d0)-u*vstar*dcmplx(SD(ish1,ish2),0.0d0)
                       SROT (2, ish1, 2, ish2) =  v*vstar*dcmplx(SD(ish1,ish2),0.0d0)+u*ustar*dcmplx(SD(ish1,ish2),0.0d0)
                       
                    END IF 
                  ! Rotate INTER-ATOMIC SELF-CONSISTENT COLLINEAR and OVERLAP matrix elements from local to global spin space    
                  ELSE IF (NSpinRot > 0) THEN !.and. distance > 0.0d0) THEN 
                      
                      CALL Pauli_Matrices( SpinRotTheta(AtomID1), SpinRotPhi(AtomID1), sigma_z, sigma_p, sigma_m, u, v, ustar, vstar )
                      ! Hamiltonian
                      HROT1temp11 = u*ustar*dcmplx(HD(1,ish1,ish2),0.0d0)+v*vstar*dcmplx(HD(2,ish1,ish2),0.0d0) 
                      HROT1temp12 = v*ustar*dcmplx(HD(1,ish1,ish2),0.0d0)-v*ustar*dcmplx(HD(2,ish1,ish2),0.0d0) 
                      HROT1temp21 = u*vstar*dcmplx(HD(1,ish1,ish2),0.0d0)-u*vstar*dcmplx(HD(2,ish1,ish2),0.0d0) 
                      HROT1temp22 = v*vstar*dcmplx(HD(1,ish1,ish2),0.0d0)+u*ustar*dcmplx(HD(2,ish1,ish2),0.0d0)     
                      ! Overlap matrix
                      SROT1temp11 = u*ustar*dcmplx(SD(ish1,ish2),0.0d0)+v*vstar*dcmplx(SD(ish1,ish2),0.0d0)
                      SROT1temp12 = v*ustar*dcmplx(SD(ish1,ish2),0.0d0)-v*ustar*dcmplx(SD(ish1,ish2),0.0d0)
                      SROT1temp21 = u*vstar*dcmplx(SD(ish1,ish2),0.0d0)-u*vstar*dcmplx(SD(ish1,ish2),0.0d0)
                      SROT1temp22 = v*vstar*dcmplx(SD(ish1,ish2),0.0d0)+u*ustar*dcmplx(SD(ish1,ish2),0.0d0)
                                                                                                                                      
                      CALL Pauli_Matrices( SpinRotTheta(AtomID2), SpinRotPhi(AtomID2), sigma_z, sigma_p, sigma_m, u, v, ustar, vstar )
                      ! Hamiltonian
                      HROT2temp11 = u*ustar*dcmplx(HD(1,ish1,ish2),0.0d0)+v*vstar*dcmplx(HD(2,ish1,ish2),0.0d0) 
                      HROT2temp12 = v*ustar*dcmplx(HD(1,ish1,ish2),0.0d0)-v*ustar*dcmplx(HD(2,ish1,ish2),0.0d0) 
                      HROT2temp21 = u*vstar*dcmplx(HD(1,ish1,ish2),0.0d0)-u*vstar*dcmplx(HD(2,ish1,ish2),0.0d0) 
                      HROT2temp22 = v*vstar*dcmplx(HD(1,ish1,ish2),0.0d0)+u*ustar*dcmplx(HD(2,ish1,ish2),0.0d0)     
                      ! Overlap matrix
                      SROT2temp11 = u*ustar*dcmplx(SD(ish1,ish2),0.0d0)+v*vstar*dcmplx(SD(ish1,ish2),0.0d0)
                      SROT2temp12 = v*ustar*dcmplx(SD(ish1,ish2),0.0d0)-v*ustar*dcmplx(SD(ish1,ish2),0.0d0)
                      SROT2temp21 = u*vstar*dcmplx(SD(ish1,ish2),0.0d0)-u*vstar*dcmplx(SD(ish1,ish2),0.0d0)
                      SROT2temp22 = v*vstar*dcmplx(SD(ish1,ish2),0.0d0)+u*ustar*dcmplx(SD(ish1,ish2),0.0d0)
				  
                      ! Average Hamiltonian elements  
                      HROT (1, ish1, 1, ish2) = (HROT1temp11 + HROT2temp11)/2.0d0
                      HROT (1, ish1, 2, ish2) = (HROT1temp12 + HROT2temp12)/2.0d0
                      HROT (2, ish1, 1, ish2) = (HROT1temp21 + HROT2temp21)/2.0d0
                      HROT (2, ish1, 2, ish2) = (HROT1temp22 + HROT2temp22)/2.0d0 
                      ! Average Overlap matrix elemnets
                      SROT (1, ish1, 1, ish2) = (SROT1temp11 + SROT2temp11)/2.0d0
                      SROT (1, ish1, 2, ish2) = (SROT1temp12 + SROT2temp12)/2.0d0
                      SROT (2, ish1, 1, ish2) = (SROT1temp21 + SROT2temp21)/2.0d0
                      SROT (2, ish1, 2, ish2) = (SROT1temp22 + SROT2temp22)/2.0d0 
                                          
                  END IF  
                ish2 = ish2 + 1
                END DO    
             END DO                      
             ish1 = ish1 + 1       
             !Print *, ish1, ish2, NAOs             
             ish2 = 1
           END DO
    END DO                
    
    !Print *, ish1, ish2      
    
    Print *, "LS1 (real part) = "                                                                          
    DO s1 = 1,2                                                                                
       DO i = 1,3                                                                              
          PRINT '(100(F11.5))', ( REAL(LS1(s1,i,1,j)), j=1,3 ), ( REAL(LS1(s1,i,2,j)), j=1,3 ) 
       END DO                                                                                  
    END DO   
    Print *, "LS1 (imaginary part) = "                 
    DO s1 = 1,2                                                                                
       DO i = 1,3                                                                              
          PRINT '(100(F11.5))', ( AIMAG(LS1(s1,i,1,j)), j=1,3 ), ( AIMAG(LS1(s1,i,2,j)), j=1,3 ) 
       END DO                                                                                  
    END DO                                                                                                                                                          
                                                                                               
    Print *, "LS2 (real part) = "                                                                          
    DO s1 = 1,2                                                                                
       DO i = 1,5                                                                              
          PRINT '(100(F11.5))', ( REAL(LS2(s1,i,1,j)), j=1,5 ), (REAL( LS2(s1,i,2,j)), j=1,5)  
       END DO                                                                                  
    END DO               
    Print *, "LS2 (imaginary part) = "
    DO s1 = 1,2                                                                                
       DO i = 1,5                                                                              
          PRINT '(100(F11.5))', ( AIMAG(LS2(s1,i,1,j)), j=1,5 ), (AIMAG( LS2(s1,i,2,j)), j=1,5)  
       END DO                                                                                  
    END DO                                                                              
                                                                                               
    Print *, "LS3 (real part) = "                                                                          
    DO s1 = 1,2                                                                                
       DO i = 1,7                                                                              
          PRINT '(100(F11.5))', ( REAL(LS3(s1,i,1,j)), j=1,7 ), (REAL( LS3(s1,i,2,j)), j=1,7)  
       END DO                                                                                  
    END DO   
    Print *, "LS3 (imaginary part) = "
    DO s1 = 1,2                                                                                
       DO i = 1,7                                                                              
          PRINT '(100(F11.5))', ( AIMAG(LS3(s1,i,1,j)), j=1,7 ), (AIMAG( LS3(s1,i,2,j)), j=1,7)  
       END DO                                                                                  
    END DO       
                                                                                                                                                  
    
    !*****************************************************
    !Construct hamil and hamil_SO matrices to return to ANT
    !*****************************************************
    
    DO i = 1,NAOs
       DO j = 1,NAOs
            !Up-Up
    		hamil_SO(i,j) = HSO(1,i,1,j) 
    		IF (NSPinROT > 0) THEN
    		   hamil(i,j) = HROT(1,i,1,j)
    		   overlap_SO(i,j) = SROT(1,i,1,j)
    		END IF   
    		!Up-Down
    		hamil_SO(i,j+NAOs) = HSO(1,i,2,j)
    		IF (NSPinROT > 0) THEN
    		   hamil(i,j+NAOs) = HROT(1,i,2,j)
    		   overlap_SO(i,j+NAOs) = SROT(1,i,2,j)
    		END IF   
            !Down-Up
    		hamil_SO(i+NAOs,j) = HSO(2,i,1,j)
    		IF (NSPinROT > 0) THEN
    		   hamil(i+NAOs,j) = HROT(2,i,1,j)
    		   overlap_SO(i+NAOs,j) = SROT(2,i,1,j)
    		END IF   
            !Down-Down
    		hamil_SO(i+NAOs,j+NAOs) = HSO(2,i,2,j)
    		IF (NSPinROT > 0) THEN
    		   hamil(i+NAOs,j+NAOs) = HROT(2,i,2,j) 
    		   overlap_SO(i+NAOs,j+NAOs) = SROT(2,i,2,j) 
    		END IF   
       END DO        
    END DO  
      
   CALL FLUSH(6)
  END SUBROUTINE CompHSO
 
END MODULE SpinOrbit
