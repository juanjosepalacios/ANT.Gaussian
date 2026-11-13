!***********************************************************************
! Module lapack_blas for the NVIDIA build with the optimized BLAS/LAPACK
! of the NVIDIA HPC SDK (OpenBLAS, 64-bit integers, liblapack_ilp64 and
! libblas_ilp64).
!
! It replaces lapack_blas.f (the bundled reference routines) with the
! same names and arguments, so that the rest of ANT is unchanged. The
! NVIDIA libraries are linked under renamed symbols (antnv_ in front of
! every Fortran name, see the Makefile), because Gaussian's own BLAS and
! LAPACK (util.a, ATLAS) are linked into l502.exe with the same names.
!
! OpenBLAS here uses the same OpenMP runtime as ANT: called inside an ANT
! parallel region it runs on one thread, outside it uses all threads.
!***********************************************************************
module lapack_blas
  implicit none
  external antnv_dsyev, antnv_zgemm, antnv_zgetri, antnv_zgetrf, antnv_zgeev, antnv_zheev, antnv_dgemm
contains

  subroutine DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
    character          JOBZ, UPLO
    integer            INFO, LDA, LWORK, N
    double precision   A( LDA, * ), W( * ), WORK( * )
    call antnv_dsyev( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
  end subroutine DSYEV

  subroutine ZGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
    double complex     ALPHA, BETA
    integer            K, LDA, LDB, LDC, M, N
    character          TRANSA, TRANSB
    double complex     A( LDA, * ), B( LDB, * ), C( LDC, * )
    call antnv_zgemm( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
  end subroutine ZGEMM

  subroutine ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
    integer            INFO, LDA, LWORK, N
    integer            IPIV( * )
    complex*16         A( LDA, * ), WORK( * )
    call antnv_zgetri( N, A, LDA, IPIV, WORK, LWORK, INFO )
  end subroutine ZGETRI

  subroutine ZGETRF( M, N, A, LDA, IPIV, INFO )
    integer            INFO, LDA, M, N
    integer            IPIV( * )
    complex*16         A( LDA, * )
    call antnv_zgetrf( M, N, A, LDA, IPIV, INFO )
  end subroutine ZGETRF

  subroutine ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
    character          JOBVL, JOBVR
    integer            INFO, LDA, LDVL, LDVR, LWORK, N
    double precision   RWORK( * )
    complex*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), WORK( * )
    call antnv_zgeev( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
  end subroutine ZGEEV

  subroutine ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
    character          JOBZ, UPLO
    integer            INFO, LDA, LWORK, N
    double precision   RWORK( * ), W( * )
    complex*16         A( LDA, * ), WORK( * )
    call antnv_zheev( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
  end subroutine ZHEEV

  subroutine DGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
    double precision   ALPHA, BETA
    integer            K, LDA, LDB, LDC, M, N
    character          TRANSA, TRANSB
    double precision   A( LDA, * ), B( LDB, * ), C( LDC, * )
    call antnv_dgemm( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
  end subroutine DGEMM

end module lapack_blas
