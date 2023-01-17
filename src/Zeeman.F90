!***************************************
!*                                     *
!*  ANT.G-2.6.2  - Zeeman.f90          *
!*                                     *
!*  Calculation of Zeeman energy       *
!*                                     *
!***************************************
!*                                     *
!*  Copyright (c) 2006-2022 by         *
!*                                     *
!*  Dr. David Jacob                    *
!*                                     *
!*  MPI fuer Mikrostrukturphysik       *
!*  Weinberg 2                         *
!*  06108 Halle                        *
!*  Germany                            *
!*				                       *
!*  Dr. Wynand Dednam                  *
!*                                     *
!*  University of South Africa         *
!*  28 Pioneer Ave                     *
!*  Florida Park                       *
!*  Roodepoort                         *
!*  1709                               *
!*  South Africa                       *
!*                                     *
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
!*  Module for calculation of Zeeman energy based        *
!*  on subroutine by Joaquin Fernandez-Rossier           *
!*                                                       *
!*********************************************************
MODULE Zeeman
  IMPLICIT NONE
  
  PUBLIC :: CompHZM 
  
CONTAINS  
  
  !**************************************************************
  !*** Pauli z matrix in basis of alpha- and beta-spin states ***
  !*** quantization axis z of |alpha>, |beta> parallel to     ***
  !*** global z axis                                          ***
  !**************************************************************
  SUBROUTINE PAULI_MATRIX(sigma_z)
    IMPLICIT NONE

    !Output: Pauli z matrix 
    COMPLEX*16, DIMENSION(2,2), INTENT(out) :: sigma_z, sigma_p, sigma_m, sigma_x, sigma_y
    
    sigma_z = 0.0d0
    sigma_p = 0.0d0
    sigma_m = 0.0d0
    
    sigma_z(1,1)=1.d0        
    sigma_z(1,2)=0.d0        
    sigma_z(2,1)=0.d0        
    sigma_z(2,2)=-1.d0       
    
    sigma_p(1,1)=0.d0      
    sigma_p(1,2)=2.d0      
    sigma_p(2,1)=0.d0      
    sigma_p(2,2)=0.d0      
                           
    sigma_m(1,1)=0.d0                    
    sigma_m(1,2)=0.d0                    
    sigma_m(2,1)=2.d0                    
    sigma_m(2,2)=0.d0                       

  END SUBROUTINE PAULI_MATRIX

  !*************************************************************************************
  !*** Assumed basis: Lz eigen states, Order for L >= 2: m = 0, +1, -1, +2, -2, ...  ***
  !*************************************************************************************
  SUBROUTINE L_MATRIX(L, Lz)
    USE constants
    IMPLICIT NONE
    
    ! Input: total angular momentum number L = 0,1,2,...
    INTEGER, INTENT(in) :: L    

    !Ouput: Angular momentum matrices
    COMPLEX*16, DIMENSION(2*L+1,2*L+1), INTENT(out) :: Lz
    COMPLEX*16, DIMENSION(2*L+1,2*L+1) :: temp
    COMPLEX*16 :: Us, c0, c1, c2
    COMPLEX*16, dimension(3,3) :: Up
    COMPLEX*16, dimension(5,5) :: Ud
    COMPLEX*16, dimension(7,7) :: Uf
    
    INTEGER  :: n1, i
    REAL*8 :: m, m1
    
    Lz = 0.0d0
    
    IF (L > 1) THEN
        DO n1=1,2*L+1
           ! m1 = 0, +1, -1, +2, -2, ... for L = 2, 3, ... shells in Gaussian
           m1 = n1 / 2 * ( 1 - 2*MOD(n1,2) )
           
           !<m1|Lz|m2> matrix is diagonal :
           Lz(n1,n1) = m1
               
        END DO
    ELSE IF (L == 1) THEN
         i = 0    
         DO n1=2*L+1,1,-1
           ! m1 = +1, -1, 0 for P shells in Gaussian
           m1 = n1 / 2 * ( 2*MOD(n1,2) - 1 )
           
           i = i + 1
           
           !<m1|Lz|m2> matrix is diagonal :
           Lz(i,i) = m1
                 
        END DO        
    END IF

    !PRINT *, " Lz in spherical harmonic basis = "
    !DO n1=1,2*L+1
    !   PRINT '(100(F11.5))', ( REAL(Lz(n1,n2)), n2=1,2*L+1 )  
    !END DO 
    
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
                               c2/ui,  c2/ui, c0, &                  
                               c0,     c0,    c1  /), (/ 3, 3 /) ))  
    Ud = TRANSPOSE(reshape((/  c1,     c0,     c0,     c0,     c0,    &                                   
                               c0,    -c2,     c2,     c0,     c0,    &                                   
                               c0,     c2/ui,  c2/ui, c0,      c0,    &                                   
                               c0,     c0,     c0,     c2,     c2,    &                                    
                               c0,     c0,     c0,    -c2/ui,  c2/ui  /), (/ 5, 5 /) ))                    
    Uf = TRANSPOSE(reshape((/  c1,     c0,     c0,     c0,     c0,     c0,     c0,    &         
                               c0,    -c2,     c2,     c0,     c0,     c0,     c0,    &                
                               c0,     c2/ui,  c2/ui,  c0,     c0,     c0,     c0,    &                
                               c0,     c0,     c0,     c2,     c2,     c0,     c0,    &                
                               c0,     c0,     c0,    -c2/ui,  c2/ui,  c0,     c0,    &
                               c0,     c0,     c0,     c0,     c0,    -c2,     c2,    &
                               c0,     c0,     c0,     c0,     c0,     c2/ui,  c2/ui  /), (/ 7, 7 /) )) 
              
    !
    ! transform matrices to cubic harmonic basis
    !    
    
    IF( L == 1 ) THEN
      !PRINT *, " Up = "
      !DO n1=1,2*L+1
      !   PRINT '(100(F11.5))', ( REAL(Up(n1,n2)), n2=1,2*L+1 )
      !END DO      
      temp = MATMUL( Lz, CONJG(TRANSPOSE(Up)))
      Lz = MATMUL( Up, temp )
      
    ELSE IF (L == 2) THEN         
      !PRINT *, " Ud = "
      !DO n1=1,2*L+1
      !   PRINT '(100(F11.5))', ( REAL(Ud(n1,n2)), n2=1,2*L+1 ) 
      !END DO         
      temp = MATMUL( Lz, CONJG(TRANSPOSE(Ud)))    
      Lz = MATMUL( Ud, temp )                     
                                                                  
    ELSE IF (L == 3) THEN         
      !PRINT *, " Uf = "
      !DO n1=1,2*L+1
      !   PRINT '(100(F11.5))', ( REAL(Uf(n1,n2)), n2=1,2*L+1 ) 
      !END DO         
      temp = MATMUL( Lz, CONJG(TRANSPOSE(Uf)))    
      Lz = MATMUL( Uf, temp )                     
                                                                                
    END IF    
    
    
    !PRINT *, " Lz in cartesian basis (real part) = "
    !DO n1=1,2*L+1
    !   PRINT '(100(F11.5))', ( REAL(Lz(n1,n2)), n2=1,2*L+1 )
    !END DO
    !PRINT *, " Lz in cartesian basis (imaginary part) = "
    !DO n1=1,2*L+1
    !   PRINT '(100(F11.5))', ( AIMAG(Lz(n1,n2)), n2=1,2*L+1 )
    !END DO


  END SUBROUTINE L_MATRIX 
  
  !*****************************************************************
  !*** Compute matrix of Zeemn Hamiltonian for a given basis set ***
  !*****************************************************************
  SUBROUTINE CompHZM(hamil_ZM,NAOs,Nshell)
    USE parameters, ONLY: Bx, By, Bz, NZMAtom, ZMAtomBx, ZMAtomBy, ZMAtomBz
    USE G09common, ONLY : GetNAtoms, GetShellT, GetShellC, GetAtm4Sh
    USE cluster, ONLY : LoAOrbNo, HiAOrbNo
    USE constants
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NAOs, NShell
    INTEGER :: L1, L2, L3, AtomID1, AtomID2, ShellT1, ShellT2, ShellC1, ShellC2
    INTEGER, DIMENSION(NAOs) :: AOT

    !*********************************
    !Output: Atomic LS coupling matrix 
    !*********************************
    COMPLEX*16, DIMENSION(2,NAOs,2,NAOs) :: HZM  ! Need to reshape HZM to H_ZM(i,j) where i = 1, 2*Norb, j = 1, 2*Norb
    COMPLEX*16, DIMENSION(2*NAOs,2*NAOs), INTENT(OUT) :: hamil_ZM

    INTEGER :: i, j, k, q, s1, s2, ish1, ish2
    REAL*8 :: zz, Bx_atom, By_atom, Bz_atom

    COMPLEX*16, DIMENSION(2,2) :: sigma_z
    COMPLEX*16, DIMENSION(3,3) :: L_z1
    COMPLEX*16, DIMENSION(5,5) :: L_z2
    COMPLEX*16, DIMENSION(7,7) :: L_z3
    
    REAL*8, PARAMETER :: hbar=4.135667696e-15 !eV.s
    REAL*8, PARAMETER :: muB = 	5.7883818012e-5 ! Bohr Magneton eV/T
    
    L1 = 1
    L2 = 2
    L3 = 3
    
    ! Lz, L+, L- Matrices
    !PRINT *, "Calculating Lz matrix for p orbitals... "
    CALL L_Matrix( L1, L_z1) ! p-orbitals
    !PRINT *, "Calculating Lz matrix for d orbitals... "
    CALL L_Matrix( L2, L_z2) ! d-orbitals
    !PRINT *, "Calculating Lz matrix for f orbitals... "
    CALL L_Matrix( L3, L_z3) ! f-orbitals    
    
    !PRINT *, "Total number of shells in the cluster:"
    !PRINT *, NShell
      
    ish1 = 1
    
    DO i=1,NShell
    	  ShellT1 = GetShellT(i) 
    		DO k=1,2*ShellT1+1
    			AOT(ish1) = ShellT1
    			ish1 = ish1 + 1
            END DO
    END DO        
    
    !PRINT *, "Shell type of each orbital:"
    !DO i=1,NAOs
    !  PRINT *, AOT(i)
    !END DO  
      
    !PRINT *, "Basis set of each atom extracted from internal Gaussian variables:"
    
    zz = muB/hbar 
                

    !***************************************
    !Construct HZM matrix for complete basis
    !***************************************
                
    HZM = c_zero
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
                    IF( NZMAtom > 0) THEN  ! User-defined Zeeman field at atom AtomID2
                      Bx_atom = ZMAtomBx(AtomID2)
                      By_atom = ZMAtomBy(AtomID2)
                      Bz_atom = ZMAtomBz(AtomID2)
                    ELSE
                      Bx_atom = Bx
                      By_atom = By
                      Bz_atom = Bz
                    END IF	          	                   	    
                    CALL Pauli_Matrix( sigma_z, sigma_p, sigma_m)                     	        
          	        DO s1 = 1,2        
          	          DO s2 = 1,2
                         IF( ShellT1 == 0 .and. ShellT2 == 0)THEN   
                               ! Strong-field anomalous Zeeman effect. See Zettili, Quantum Mechanics, 2nd ed., section 9.2.3.4. pp 484-485     ! S = 0.5*hbar*sigma        
                               HZM( s1, ish1, s2, ish2 ) = zz*hbar*(0.5d0*(sigma_p(s1,s2)+sigma_m(s1,s2))*Bx_atom + 0.5d0*ui*(sigma_m(s1,s2)-sigma_p(s1,s2))*By_atom + sigma_z(s1,s2)*Bz_atom)    
                         ELSE IF( ShellT1 == 1 .and. ShellT2 == 1 .and. (ShellC1 == 0 .or. ShellC1 == 1) .and. (ShellC2 == 0 .or. ShellC2 == 1)) THEN
                               !HZM( s1, ish1, s2, ish2) = zz*(L_z1(j,q)+hbar*sigma_z(s1,s2))                                    
                               HZM( s1, ish1, s2, ish2) = zz*hbar*(0.5d0*(sigma_p(s1,s2)+sigma_m(s1,s2))*Bx_atom + 0.5d0*ui*(sigma_m(s1,s2)-sigma_p(s1,s2))*By_atom + sigma_z(s1,s2)*Bz_atom)
                         ELSE IF( ShellT1 == 2 .and. ShellT2 == 2 .and. ShellC1 == 2 .and. ShellC2 == 2) THEN
                               !HZM( s1, ish1, s2, ish2) = zz*(L_z2(j,q)+hbar*sigma_z(s1,s2))
                               HZM( s1, ish1, s2, ish2) = zz*hbar*(0.5d0*(sigma_p(s1,s2)+sigma_m(s1,s2))*Bx_atom + 0.5d0*ui*(sigma_m(s1,s2)-sigma_p(s1,s2))*By_atom + sigma_z(s1,s2)*Bz_atom)
                         ELSE IF( ShellT1 == 3 .and. ShellT2 == 3 .and. ShellC1 == 2 .and. ShellC2 == 2) THEN
                               !HZM( s1, ish1, s2, ish2) = zz*(L_z3(j,q)+hbar*sigma_z(s1,s2))
                               HZM( s1, ish1, s2, ish2) = zz*hbar*(0.5d0*(sigma_p(s1,s2)+sigma_m(s1,s2))*Bx_atom + 0.5d0*ui*(sigma_m(s1,s2)-sigma_p(s1,s2))*By_atom + sigma_z(s1,s2)*Bz_atom)
                         END IF                              
                      END DO
                    END DO                
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
                                                                                                                                                  
    
    !******************************************
    !Construct hamil_ZM matrix to return to ANT
    !******************************************
             
    
    DO i = 1,NAOs
       DO j = 1,NAOs
            !Up-Up
    		hamil_ZM(i,j) = HZM(1,i,1,j) 
    		!Up-Down
    		hamil_ZM(i,j+NAOs) = HZM(1,i,2,j)
            !Down-Up
    		hamil_ZM(i+NAOs,j) = HZM(2,i,1,j)
            !Down-Down
    		hamil_ZM(i+NAOs,j+NAOs) = HZM(2,i,2,j)
       END DO        
    END DO  
      
   CALL FLUSH(6)
  END SUBROUTINE CompHZM
 
END MODULE Zeeman
