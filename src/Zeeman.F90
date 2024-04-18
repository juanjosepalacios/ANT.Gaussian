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
  SUBROUTINE PAULI_MATRIX(sigma_z, sigma_p, sigma_m)
    IMPLICIT NONE

    !Output: Pauli z matrix 
    COMPLEX*16, DIMENSION(2,2), INTENT(out) :: sigma_z, sigma_p, sigma_m
    
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

   
  !*****************************************************************
  !*** Compute matrix of Zeemn Hamiltonian for a given basis set ***
  !*****************************************************************
  SUBROUTINE CompHZM(hamil_ZM,S,NAOs,Nshell)
    USE parameters, ONLY: Bx, By, Bz, NZMAtom, ZMAtomBx, ZMAtomBy, ZMAtomBz
    USE G09common, ONLY : GetNAtoms, GetShellT, GetShellC, GetAtm4Sh
    USE cluster, ONLY : LoAOrbNo, HiAOrbNo
    USE constants
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NAOs, NShell
    INTEGER :: AtomID1, AtomID2, ShellT1, ShellT2, ShellC1, ShellC2
    INTEGER, DIMENSION(NAOs) :: AOT

    !*********************************
    !Output: Atomic LS coupling matrix 
    !*********************************
    COMPLEX*16, DIMENSION(2,NAOs,2,NAOs) :: HZM  ! Need to reshape HZM to H_ZM(i,j) where i = 1, 2*Norb, j = 1, 2*Norb
    COMPLEX*16, DIMENSION(2*NAOs,2*NAOs), INTENT(OUT) :: hamil_ZM
    REAL*8, DIMENSION(NAOs,NAOs), INTENT(IN) :: S

    INTEGER :: i, j, k, q, s1, s2, ish1, ish2
    REAL*8 :: Bx_atom, By_atom, Bz_atom

    COMPLEX*16, DIMENSION(2,2) :: sigma_z, sigma_p, sigma_m
    
    REAL*8, PARAMETER :: hbar=4.135667696e-15 !eV.s
    REAL*8, PARAMETER :: muB = 	5.7883818012e-5 ! Bohr Magneton eV/T
    
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
                HZM( s1, ish1, s2, ish1) = muB*(0.5d0*(sigma_p(s1,s2)+sigma_m(s1,s2))*Bx_atom + 0.5d0*ui*(sigma_m(s1,s2)-sigma_p(s1,s2))*By_atom + sigma_z(s1,s2)*Bz_atom)
              END DO
            END DO                
          ish1 = ish1 + 1       
          !Print *, ish1, ish2, NAOs             
          END DO
    END DO                
    
    !Print *, ish1, ish2          
                                                                                                                                                  
    
    !******************************************
    !Construct hamil_ZM matrix to return to ANT
    !******************************************
             
    
    DO i = 1,NAOs
       DO j = 1,NAOs
            !Up-Up
    		hamil_ZM(i,j) = HZM(1,i,1,j)!* S(i,j)
    		!Up-Down
    		hamil_ZM(i,j+NAOs) = HZM(1,i,2,j)!* S(i,j)
            !Down-Up
    		hamil_ZM(i+NAOs,j) = HZM(2,i,1,j)!* S(i,j)
            !Down-Down
    		hamil_ZM(i+NAOs,j+NAOs) = HZM(2,i,2,j)!* S(i,j)
       END DO        
    END DO  
      
   CALL FLUSH(6)
  END SUBROUTINE CompHZM
 
END MODULE Zeeman
