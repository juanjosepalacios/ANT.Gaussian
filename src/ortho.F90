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
  MODULE Orthogonalization
!**********************************************************
! Module for Projective Orthogonalization Procedure       *
!---------------------------------------------------------*
!                                                         *
! Gram-Schmidth-like Orthogonalization with               *
! respect to some subspace C given by cix index array     *
!                                                         *
! Procedure:                                              *
! 1. Loewdin Orthogonalization of subspace C              *
! 2. Loewdin Orthogonalization of whole system            *
! 3. Unitary rotation of subspace back to                 *
!    original (orthogonalized) subspace C'                *
! 4. Orthogonalization of rest of system R = S/C          *
!                                                         *
!**********************************************************
  implicit none
contains
  
  subroutine ProjOrtho( cix, S, O )
    use constants
    use numeric, only: RSetId, RMatPow
    use util
    implicit none
    !----------------------------------------------------
    ! INPUT: 
    ! Index array specifying subspace C:
    ! .true.  = is part of subspace C
    ! .false. = not part of subspace C  
    logical, dimension(:), intent(in) :: cix
    !----------------------------------------------------
    ! INPUT/OUTPUT:
    ! Overlap matrix of entire system 
    real*8, dimension(:,:), intent(inout) :: S
    !----------------------------------------------------
    ! OUTPUT:
    ! Orthogonalization matrix
    real*8, dimension(:,:), intent(out) :: O
    !----------------------------------------------------

    ! temporary matrix, deorthogonalization matrix S^+1/2, and transformation matrix
    real*8, dimension(:,:),allocatable :: temp, SPH, Trafo
    
    ! Dimension of entore vector space, indixes
    integer :: ndim, i, j, k

    ndim = size(cix)
    if( .not. any(cix(:)) )then
       print *, "ORTHOGONALIZATION/ProjOrtho/Warning: Correlated block not specified. Projective Orthogonalization not possible." 
       stop
    end if

    if( ndim <= 0 .or. size(S,1) /= ndim .or. size(S,2) /=ndim .or. size(O,1) /= ndim .or. size(O,2) /= ndim )then
       print *, "Orthogonalization/ProjOrtho: Error. Abort."
       stop
    end if
    
    allocate( temp(ndim,ndim), SPH(ndim,ndim), Trafo(ndim,ndim) )

    print *
    print *, "============================"    
    print *, "Projective Orthogonalization"
    print *, "============================"    
    print *
    
    print *
    print *, "1. Orthogonalization of subspace C"
    print *
    ! Construct Overlap matrix within subspace C
    call RSetId(temp)
    do i=1,ndim
       do j=1,ndim
          if( cix(i) .and. cix(j) ) temp(i,j) = S(i,j)
       end do
    end do
    ! Loewdin Orthogonalization for subspace C
    call RMatPow( temp, -0.5d0, Trafo )
    temp = matmul( Trafo, S )
    S = matmul( temp, Trafo )

    O = Trafo

    print *
    print *, "2. Loewdin orthogonalization of entire system C+R"
    print *

    ! Loewdin orthogonalization matrix S^-1/2 for entire system
    call RSetID( Trafo )
    call RMatPow( S, -0.5d0, Trafo )
    ! Loewdin de-orthogonalization matrix S^+1/2
    call RSetID( SPH )
    call RMatPow( S, +0.5d0, SPH )
    ! Loewdin orthogonalization of entire system
    call RSetID( S )

    temp = matmul( O, Trafo )
    O = temp

    print *
    print *, "3. Rotation back to orthogonalized subspace C"
    print *
    
    call RSetId(Trafo)
    do i=1,ndim
       do j=1,ndim
          if( cix(j) ) Trafo(i,j) = SPH(i,j)
          if( .not. cix(j) )then
             do k=1,ndim
                if( cix(k) ) Trafo(i,j) = Trafo(i,j) - SPH(i,k)*SPH(k,j)
             end do
          end if
       end do
    end do

    temp = matmul( transpose(Trafo), S )
    S = matmul( temp, Trafo )

    temp = matmul( O, Trafo )
    O = temp

    print *
    print *, "4. Ortho-normalize subspace R"
    print *

    call RMatPow( S, -0.5d0, Trafo )
    call RSetId( S )

    temp = matmul( O, Trafo )
    O = temp

    deallocate( temp, SPH, Trafo )

    print *
    print *, "Projektive Orthogonalization finished."
    print *

  end subroutine ProjOrtho

end module Orthogonalization
