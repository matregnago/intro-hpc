!!!
!
! @file chameleon_zcf90.F90
!
!  CHAMELEON Fortran 90 interfaces using Fortran 2003 ISO C bindings
!  CHAMELEON is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 1.4.0
! @author Numerical Algorithm Group
! @author Mathieu Faverge
! @author Emmanuel Agullo
! @author Cedric Castagnede
! @author Florent Pruvost
! @date 2024-02-18
! @precisions mixed zc -> ds
!
! -- Inria
! -- (C) Copyright 2012
!
! This software is a computer program whose purpose is to process
! Matrices Over Runtime Systems @ Exascale (MORSE). More information
! can be found on the following website: http://www.inria.fr/en/teams/morse.
!
! This software is governed by the CeCILL-B license under French law and
! abiding by the rules of distribution of free software.  You can  use,
! modify and/ or redistribute the software under the terms of the CeCILL-B
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and,  more generally, to use and operate it in the
! same conditions as regards security.
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL-B license and that you accept its terms.
!
!!!

!
!     Copyright © 2011 The Numerical Algorithms Group Ltd. All rights reserved.
!
!     Redistribution and use in source and binary forms, with or without
!     modification, are permitted provided that the following conditions are
!     met:
!     - Redistributions of source code must retain the above copyright notice,
!       this list of conditions, and the following disclaimer.
!     - Redistributions in binary form must reproduce the above copyright
!       notice, this list of conditions and the following disclaimer listed in
!       this license in the documentation and/or other materials provided with
!       the distribution.
!     - Neither the name of the copyright holders nor the names of its
!       contributors may be used to endorse or promote products derived from
!       this software without specific prior written permission.
!
!     This software is provided by the copyright holders and contributors "as
!     is" and any express or implied warranties, including, but not limited
!     to, the implied warranties of merchantability and fitness for a
!     particular purpose are disclaimed. in no event shall the copyright owner
!     or contributors be liable for any direct, indirect, incidental, special,
!     exemplary, or consequential damages (including, but not limited to,
!     procurement of substitute goods or services; loss of use, data, or
!     profits; or business interruption) however caused and on any theory of
!     liability, whether in contract, strict liability, or tort (including
!     negligence or otherwise) arising in any way out of the use of this
!     software, even if advised of the possibility of such damage.
!
!
module chameleon_zc

      interface
         function CHAMELEON_zcgesv_c(N,NRHS,A,LDA,IPIV,B,LDB,X,LDX,ITER) &
          & bind(c, name='CHAMELEON_zcgesv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: CHAMELEON_zcgesv_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            type(c_ptr), value :: X
            integer(kind=c_int), value :: LDX
            type(c_ptr), value :: ITER
         end function CHAMELEON_zcgesv_c
      end interface

      interface
         function CHAMELEON_zcposv_c(uplo,N,NRHS,A,LDA,B,LDB,X,LDX,ITER) &
          & bind(c, name='CHAMELEON_zcposv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: CHAMELEON_zcposv_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            type(c_ptr), value :: X
            integer(kind=c_int), value :: LDX
            type(c_ptr), value :: ITER
         end function CHAMELEON_zcposv_c
      end interface

      interface
         function CHAMELEON_zcungesv_c(trans,N,NRHS,A,LDA,B,LDB,X,LDX,ITER) &
          & bind(c, name='CHAMELEON_zcungesv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: CHAMELEON_zcungesv_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            type(c_ptr), value :: X
            integer(kind=c_int), value :: LDX
            type(c_ptr), value :: ITER
         end function CHAMELEON_zcungesv_c
      end interface

      interface
         function CHAMELEON_zcgesv_Tile_c(A,IPIV,B,X,ITER) &
          & bind(c, name='CHAMELEON_zcgesv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: CHAMELEON_zcgesv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            type(c_ptr), value :: X
            type(c_ptr), value :: ITER
         end function CHAMELEON_zcgesv_Tile_c
      end interface

      interface
         function CHAMELEON_zcposv_Tile_c(uplo,A,B,X,ITER) &
          & bind(c, name='CHAMELEON_zcposv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: CHAMELEON_zcposv_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: X
            type(c_ptr), value :: ITER
         end function CHAMELEON_zcposv_Tile_c
      end interface

      interface
         function CHAMELEON_zcungesv_Tile_c(trans,A,T,B,X,ITER) &
          & bind(c, name='CHAMELEON_zcungesv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: CHAMELEON_zcungesv_Tile_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: X
            type(c_ptr), value :: ITER
         end function CHAMELEON_zcungesv_Tile_c
      end interface

      interface
         function CHAMELEON_zcgesv_Tile_Async_c(A,IPIV,B,X,ITER,sequence,request) &
          & bind(c, name='CHAMELEON_zcgesv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: CHAMELEON_zcgesv_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            type(c_ptr), value :: X
            type(c_ptr), value :: ITER
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
         end function CHAMELEON_zcgesv_Tile_Async_c
      end interface

      interface
         function CHAMELEON_zcposv_Tile_Async_c(uplo,A,B,X,ITER,sequence,request) &
          & bind(c, name='CHAMELEON_zcposv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: CHAMELEON_zcposv_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: X
            type(c_ptr), value :: ITER
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
         end function CHAMELEON_zcposv_Tile_Async_c
      end interface

      interface
         function CHAMELEON_zcungesv_Tile_Async_c(trans,A,T,B,X,ITER,sequence,request) &
          & bind(c, name='CHAMELEON_zcungesv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: CHAMELEON_zcungesv_Tile_Async_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: X
            type(c_ptr), value :: ITER
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
         end function CHAMELEON_zcungesv_Tile_Async_c
      end interface

  contains

      subroutine CHAMELEON_zcgesv(N,NRHS,A,LDA,IPIV,B,LDB,X,LDX,ITER,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: ITER
         integer(kind=c_int), intent(out), target :: IPIV(*)
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDX
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_double_complex), intent(out), target :: X(LDX,*)
         info = CHAMELEON_zcgesv_c(N,NRHS,c_loc(A),LDA,c_loc(IPIV),c_loc(B),LDB,c_loc(X),LDX,c_loc(ITER))
      end subroutine CHAMELEON_zcgesv

      subroutine CHAMELEON_zcposv(uplo,N,NRHS,A,LDA,B,LDB,X,LDX,ITER,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: ITER
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDX
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_double_complex), intent(out), target :: X(LDX,*)
         info = CHAMELEON_zcposv_c(uplo,N,NRHS,c_loc(A),LDA,c_loc(B),LDB,c_loc(X),LDX,c_loc(ITER))
      end subroutine CHAMELEON_zcposv

      subroutine CHAMELEON_zcungesv(trans,N,NRHS,A,LDA,B,LDB,X,LDX,ITER,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: ITER
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDX
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
         complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
         complex(kind=c_double_complex), intent(out), target :: X(LDX,*)
         info = CHAMELEON_zcungesv_c(trans,N,NRHS,c_loc(A),LDA,c_loc(B),LDB,c_loc(X),LDX,c_loc(ITER))
      end subroutine CHAMELEON_zcungesv

      subroutine CHAMELEON_zcgesv_Tile(A,IPIV,B,X,ITER,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
          integer(kind=c_int), intent(out), target :: ITER
         type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: X ! Arg managed by CHAMELEON: opaque to Fortran
         info = CHAMELEON_zcgesv_Tile_c(A,c_loc(IPIV),B,X,c_loc(ITER))
      end subroutine CHAMELEON_zcgesv_Tile

      subroutine CHAMELEON_zcposv_Tile(uplo,A,B,X,ITER,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: ITER
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: X ! Arg managed by CHAMELEON: opaque to Fortran
         info = CHAMELEON_zcposv_Tile_c(uplo,A,B,X,c_loc(ITER))
      end subroutine CHAMELEON_zcposv_Tile

      subroutine CHAMELEON_zcungesv_Tile(trans,A,T,B,X,ITER,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: ITER
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: X ! Arg managed by CHAMELEON: opaque to Fortran
         info = CHAMELEON_zcungesv_Tile_c(trans,A,T,B,X,c_loc(ITER))
      end subroutine CHAMELEON_zcungesv_Tile

      subroutine CHAMELEON_zcgesv_Tile_Async(A,IPIV,B,X,ITER,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
         integer(kind=c_int), intent(out), target :: ITER
         type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: X ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
         info = CHAMELEON_zcgesv_Tile_Async_c(A,c_loc(IPIV),B,X,c_loc(ITER),sequence,request)
      end subroutine CHAMELEON_zcgesv_Tile_Async

      subroutine CHAMELEON_zcposv_Tile_Async(uplo,A,B,X,ITER,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: ITER
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: X ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
         info = CHAMELEON_zcposv_Tile_Async_c(uplo,A,B,X,c_loc(ITER),sequence,request)
      end subroutine CHAMELEON_zcposv_Tile_Async

      subroutine CHAMELEON_zcungesv_Tile_Async(trans,A,T,B,X,ITER,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: ITER
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: X ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
         info = CHAMELEON_zcungesv_Tile_Async_c(trans,A,T,B,X,c_loc(ITER),sequence,request)
      end subroutine CHAMELEON_zcungesv_Tile_Async

end module chameleon_zc
