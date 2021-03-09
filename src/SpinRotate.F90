!***************************************
!*                                     *
!*  ANT.G-2.5.2  - SpinRotate.F90      *
!*                                     *
!*  Calculation of Spin rotations      *
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
!*				                       *
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
!*  Module for calculation of Rotation matrix based      *
!*  on subroutine by Joaquin Fernandez-Rossier           *
!*                                                       *
!*********************************************************
MODULE SpinRotate
  IMPLICIT NONE
 
  PUBLIC :: CompHROT 
  
CONTAINS  
  
  !*****************************************************************
  !*** Rotation matrices in basis of alpha- and beta-spin states ***
  !*** quantization axis z of |alpha>, |beta> parallel to        ***
  !*** global z axis                                             ***
  !*****************************************************************
  SUBROUTINE ROT_MATRICES(theta,phi,u,v,ustar,vstar)
    USE constants, ONLY : d_pi, ui
    IMPLICIT NONE

    !Input: direction of magnetization
    REAL*8, INTENT(in) :: theta, phi
  
    !Half of the magnetization direction angles
    REAL*8 :: halfphi, halftheta, a, b

    !Complex phases defining transformation matrix
    COMPLEX*16, INTENT(out) :: u, v, ustar,vstar

    INTEGER :: i,j
    
    u = dcmplx(0.0d0,0.0d0)
    v = dcmplx(0.0d0,0.0d0)
    ustar = dcmplx(0.0d0,0.0d0)
    vstar = dcmplx(0.0d0,0.0d0)
    
    halfphi=0.5d0*phi*d_pi/180.0d0
    halftheta=0.5d0*theta*d_pi/180.0d0
    
    u=dcmplx(dcos(halfphi),-dsin(halfphi))
    v=dcmplx(dcos(halfphi),-dsin(halfphi))
    
    ustar=dconjg(u)                                                  
    vstar=dconjg(v)						
    
    a=dcos(halftheta)
    b=dsin(halftheta)
    
    u=a*u
    v=b*v
    
    ustar=a*ustar
    vstar=b*vstar
    
  END SUBROUTINE ROT_MATRICES
  
  !*******************************************************************
  !*** Compute matrix of Rotated Hamiltonian for a given basis set ***
  !*******************************************************************
  SUBROUTINE CompHROT(HD,hamilrot,SD,overlaprot,NAOs,Nshell)
    USE parameters, ONLY: theta, phi, NSpinRotAtom, SpinRotAtomTheta, SpinRotAtomPhi
    USE G09common, ONLY : GetNAtoms, GetShellT, GetAtm4Sh, GetShellN
    USE cluster, ONLY : LoAOrbNo, HiAOrbNo
    USE constants
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NAOs, NShell
    INTEGER :: AtomID1, AtomID2, ShellT1, ShellT2

    !*********************************
    !Output: Rotation matrix 
    !*********************************
    COMPLEX*16, DIMENSION(2*NAOs,2*NAOs), INTENT(INOUT) :: hamilrot, overlaprot
    
    REAL*8, DIMENSION(2,NAOs,NAOs), INTENT(IN) :: HD
    REAL*8, DIMENSION(NAOs,NAOs), INTENT(IN) :: SD

    INTEGER :: i, j, k, q, s1, s2, ish1, ish2, ispin , jspin
    REAL*8 :: theta_atom, phi_atom
    
    COMPLEX*16 :: u, v, ustar, vstar, HROT1temp11, HROT1temp12, HROT1temp21, HROT1temp22, SROT1temp11, SROT1temp22
    COMPLEX*16 :: HROT2temp11, HROT2temp12, HROT2temp21, HROT2temp22, SROT2temp11, SROT2temp22
    COMPLEX*16, DIMENSION(2,2) :: sigma_z, sigma_p, sigma_m
    
                                                                  
    !*******************************************
    !Construct Rotated matrix for complete basis
    !*******************************************                                              
                
    ish1 = 1
    ish2 = 1
    
    DO i = 1, NShell
      ShellT1 = GetShellT(i)        
      AtomID1 = GetAtm4Sh(i)
      DO j=1,2*ShellT1+1
         DO k = 1, NShell
            ShellT2 = GetShellT(k)                 
            AtomID2 = GetAtm4Sh(k)   
          	DO q=1,2*ShellT2+1
          	  IF (AtomID2 == AtomID1) THEN 
                    IF( NSpinRotAtom > 0) THEN
                      theta_atom = SpinRotAtomTheta(AtomID2)
                      phi_atom = SpinRotAtomPhi(AtomID2)
                    ELSE
                      theta_atom = theta
                      phi_atom = phi	
                    END IF                   	    
                    CALL Rot_Matrices( theta_atom, phi_atom, u, v, ustar, vstar )                     	        
                    ! Rotate INTRA-ATOMIC SELF-CONSISTENT COLLINEAR and OVERLAP matrix elements from local to global spin space
                    !Hamiltonian                                                                
                       ! Up-Up
                       hamilrot (ish1, ish2) =  u*ustar*dcmplx(HD(1,ish1,ish2),0.0d0)+v*vstar*dcmplx(HD(2,ish1,ish2),0.0d0) 
                       ! Up-Down
                       hamilrot (ish1, ish2+NAOs) =  ustar*vstar*(dcmplx(HD(1,ish1,ish2)-HD(2,ish1,ish2),0.0d0))
                       !Down-Up
                       hamilrot (ish1+NAOs, ish2) =  u*v*(dcmplx(HD(1,ish1,ish2)-HD(2,ish1,ish2),0.0d0)) 
                       !Down-Down
                       hamilrot (ish1+NAOs, ish2+NAOs) =  v*vstar*dcmplx(HD(1,ish1,ish2),0.0d0)+u*ustar*dcmplx(HD(2,ish1,ish2),0.0d0)     
                       
                    !Overlap matrix
                       ! Up-Up
                       overlaprot (ish1,ish2) =  u*ustar*dcmplx(SD(ish1,ish2),0.0d0)+v*vstar*dcmplx(SD(ish1,ish2),0.0d0)
                       !Down-Down
                       overlaprot (ish1+NAOs,ish2+NAOs) =  v*vstar*dcmplx(SD(ish1,ish2),0.0d0)+u*ustar*dcmplx(SD(ish1,ish2),0.0d0)
                       
                  ! Rotate INTER-ATOMIC SELF-CONSISTENT COLLINEAR and OVERLAP matrix elements from local to global spin space    
                  ELSE
                      
                      CALL Rot_Matrices( SpinRotAtomTheta(AtomID1), SpinRotAtomPhi(AtomID1), u, v, ustar, vstar )
                      ! Hamiltonian
                      HROT1temp11 = u*ustar*dcmplx(HD(1,ish1,ish2),0.0d0)+v*vstar*dcmplx(HD(2,ish1,ish2),0.0d0) 
                      HROT1temp12 = ustar*vstar*(dcmplx(HD(1,ish1,ish2)-HD(2,ish1,ish2),0.0d0))
                      HROT1temp21 = u*v*(dcmplx(HD(1,ish1,ish2)-HD(2,ish1,ish2),0.0d0)) 
                      HROT1temp22 = v*vstar*dcmplx(HD(1,ish1,ish2),0.0d0)+u*ustar*dcmplx(HD(2,ish1,ish2),0.0d0)     
                      ! Overlap matrix
                      SROT1temp11 = u*ustar*dcmplx(SD(ish1,ish2),0.0d0)+v*vstar*dcmplx(SD(ish1,ish2),0.0d0)
                      SROT1temp22 = v*vstar*dcmplx(SD(ish1,ish2),0.0d0)+u*ustar*dcmplx(SD(ish1,ish2),0.0d0)
                                                                                                                                      
                      CALL Rot_Matrices( SpinRotAtomTheta(AtomID2), SpinRotAtomPhi(AtomID2), u, v, ustar, vstar )
                      ! Hamiltonian
                      HROT2temp11 = u*ustar*dcmplx(HD(1,ish1,ish2),0.0d0)+v*vstar*dcmplx(HD(2,ish1,ish2),0.0d0) 
                      HROT2temp12 = ustar*vstar*(dcmplx(HD(1,ish1,ish2)-HD(2,ish1,ish2),0.0d0))
                      HROT2temp21 = u*v*(dcmplx(HD(1,ish1,ish2)-HD(2,ish1,ish2),0.0d0)) 
                      HROT2temp22 = v*vstar*dcmplx(HD(1,ish1,ish2),0.0d0)+u*ustar*dcmplx(HD(2,ish1,ish2),0.0d0)     
                      ! Overlap matrix
                      SROT2temp11 = u*ustar*dcmplx(SD(ish1,ish2),0.0d0)+v*vstar*dcmplx(SD(ish1,ish2),0.0d0)
                      SROT2temp22 = v*vstar*dcmplx(SD(ish1,ish2),0.0d0)+u*ustar*dcmplx(SD(ish1,ish2),0.0d0)
				  
                      ! Average Hamiltonian elements                      
                      hamilrot (ish1, ish2) = (HROT1temp11 + HROT2temp11)/2.0d0
                      hamilrot (ish1, ish2+NAOs) = (HROT1temp12 + HROT2temp12)/2.0d0
                      hamilrot (ish1+NAOs, ish2) = (HROT1temp21 + HROT2temp21)/2.0d0
                      hamilrot (ish1+NAOs, ish2+NAOs) = (HROT1temp22 + HROT2temp22)/2.0d0 
                      
                      ! Average Overlap matrix elemnets                     
                      overlaprot (ish1,ish2) = (SROT1temp11 + SROT2temp11)/2.0d0
                      overlaprot (ish1+NAOs,ish2+NAOs) = (SROT1temp22 + SROT2temp22)/2.0d0 
                                          
                  END IF  
                ish2 = ish2 + 1
                END DO    
             END DO                      
             ish1 = ish1 + 1       
             !Print *, ish1, ish2, NAOs             
             ish2 = 1
           END DO
    END DO                                                                                                                                          
         
   CALL FLUSH(6)
  END SUBROUTINE CompHROT
 
END MODULE SpinRotate
