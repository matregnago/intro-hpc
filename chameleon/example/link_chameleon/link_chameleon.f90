!!!
!
! @file link_chameleon.f90
!
!  CHAMELEON Fortran 90 interfaces using Fortran 2003 ISO C bindings
!  CHAMELEON is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 1.4.0
! @author Florent Pruvost
! @author Mathieu Faverge
! @date 2025-12-19
program fortran_example
    implicit none
    include 'chameleon_fortran.h'
    integer, parameter :: dp   = kind(0.d0) ! double precision
    integer, parameter :: NCPU = -1
    integer, parameter :: NGPU = 0
    integer, parameter :: N    = 500
    integer, parameter :: NRHS = 1
    double precision, dimension(N*N)    :: A, Acpy
    double precision, dimension(N*NRHS) :: B, X
    double precision :: anorm, bnorm, xnorm, res, eps=1.11022d-16
    integer :: info
    integer :: UPLO = ChamUpper
    logical :: hres


    ! Initialize CHAMELEON with main parameters
    call CHAMELEON_Init(NCPU, NGPU, info)

    ! generate A matrix with random values such that it is spd
    call CHAMELEON_dplgsy( dfloat(N), ChamUpperLower, N, A, N, 51, info )
    Acpy = A

    ! generate RHS
    call CHAMELEON_dplrnt( N, NRHS, B, N, 5673, info )
    X = B

    call CHAMELEON_dpotrf( UPLO, N, A, N, INFO )
    call CHAMELEON_dpotrs( UPLO, N, NRHS, A, N, X, N, info)

    ! compute norms to check the result
    call CHAMELEON_dlange( ChamInfNorm, N, N, Acpy, N, anorm)
    call CHAMELEON_dlange( ChamInfNorm, N, NRHS, B, N, bnorm)
    call CHAMELEON_dlange( ChamInfNorm, N, NRHS, X, N, xnorm)

    ! compute A*X-B, store the result in B
    call CHAMELEON_dgemm( ChamNoTrans, ChamNoTrans, N, NRHS, N, 1.d0, Acpy, N, X, N, -1.d0, B, N, info)
    call CHAMELEON_dlange( ChamInfNorm, N, NRHS, B, N, res)

    ! if hres = 0 then the test succeed
    ! else the test failed
    hres = .TRUE.
    hres = ( res / N / eps / (anorm * xnorm + bnorm ) > 100.0 )
    print *, "   ||Ax-b||       ||A||       ||x||       ||b|| ||Ax-b||/N/eps/(||A||||x||+||b||)"
    if (hres) then
        print *, res, anorm, xnorm, bnorm, res / N / eps / (anorm * xnorm + bnorm ), "FAILURE"
    else
        print *, res, anorm, xnorm, bnorm, res / N / eps / (anorm * xnorm + bnorm), "SUCCESS"
    endif

    ! Finalize CHAMELEON
    call CHAMELEON_Finalize(info)

end program fortran_example
