!!!
!
! @file chameleon_zf90_wrappers.F90
!
!  CHAMELEON fortran wrapper for BLAS and LAPACK subroutines.
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
! @precisions normal z -> c d s
!
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
! Wrappers to CHAMELEON functions are provided for the following BLAS
! subroutines since the CHAMELEON and BLAS interfaces match exactly:
! ZGEMM  CHAMELEON_zgemm
! ZHEMM  CHAMELEON_zhemm
! ZHER2K CHAMELEON_zher2k
! ZHERK  CHAMELEON_zherk
! ZSYMM  CHAMELEON_zsymm
! ZSYR2K CHAMELEON_zsyr2k
! ZSYRK  CHAMELEON_zsyrk
! ZTRMM  CHAMELEON_ztrmm
! ZTRSM  CHAMELEON_ztrsm
!
! Wrappers to CHAMELEON functions are provided for the following LAPACK
! subroutines since the CHAMELEON and LAPACK interfaces match exactly:
! ZGESV  CHAMELEON_zgesv
! ZGETRF CHAMELEON_zgetrf
! ZGETRS CHAMELEON_zgetrs
! ZHEGST CHAMELEON_zhegst
! ZLASWP CHAMELEON_zlaswp
! ZLAUUM CHAMELEON_zlauum
! ZPOSV  CHAMELEON_zposv
! ZSYSV  CHAMELEON_zsysv
! ZPOTRF CHAMELEON_zpotrf
! ZSYTRF CHAMELEON_zsytrf
! ZPOTRI CHAMELEON_zpotri
! ZPOTRS CHAMELEON_zpotrs
! ZSYTRS CHAMELEON_zsytrs
! ZTRTRI CHAMELEON_ztrtri
! ZLACPY CHAMELEON_zlacpy
! ZLASET CHAMELEON_zlaset
#define PRECISION_z

subroutine chameleon_wrap_ZGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: LDB
  integer, intent(in) :: N
  integer, intent(in) :: NRHS
  integer, intent(out) :: INFO
  integer, intent(out), target :: IPIV(*)
  complex(kind=wp), intent(inout), target :: A(LDA,*)
  complex(kind=wp), intent(inout), target :: B(LDB,*)
  if (.not. chameleon_initialized) call chameleon_init(24,INFO)
  ! write(*,*) " Calling CHAMELEON_ZGESV"
  call CHAMELEON_ZGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
end subroutine chameleon_wrap_ZGESV

subroutine chameleon_wrap_ZGETRF(M,N,A,LDA,IPIV,INFO)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: M
  integer, intent(in) :: N
  integer, intent(out) :: INFO
  integer, intent(out), target :: IPIV(*)
  complex(kind=wp), intent(inout), target :: A(LDA,*)
  if (.not. chameleon_initialized) call chameleon_init(24,INFO)
  ! write(*,*) " Calling CHAMELEON_ZGETRF"
  call CHAMELEON_ZGETRF(M,N,A,LDA,IPIV,INFO)
end subroutine chameleon_wrap_ZGETRF

subroutine chameleon_wrap_ZGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: LDB
  integer, intent(in) :: N
  integer, intent(in) :: NRHS
  integer, intent(in), target :: IPIV(*)
  integer, intent(out) :: INFO
  character, intent(in) :: TRANS
  complex(kind=wp), intent(in), target :: A(LDA,*)
  complex(kind=wp), intent(inout), target :: B(LDB,*)
  integer :: local_TRANS
  if(TRANS=='N' .or. TRANS=='n')then
     local_TRANS = ChamNoTrans
  else if(TRANS=='T' .or. TRANS=='t')then
     local_TRANS = ChamTrans
  else if(TRANS=='C' .or. TRANS=='c')then
     local_TRANS = ChamConjTrans
  else
     local_TRANS = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,INFO)
  ! write(*,*) " Calling CHAMELEON_ZGETRS"
  call CHAMELEON_ZGETRS(local_TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
end subroutine chameleon_wrap_ZGETRS

subroutine chameleon_wrap_ZHEGST(ITYPE,UPLO,N,A,LDA,B,LDB,INFO)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: ITYPE
  integer, intent(in) :: LDA
  integer, intent(in) :: LDB
  integer, intent(in) :: N
  integer, intent(out) :: INFO
  character, intent(in) :: UPLO
  complex(kind=wp), intent(in), target :: B(LDB,*)
  complex(kind=wp), intent(inout), target :: A(LDA,*)
  integer :: local_UPLO
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,INFO)
  ! write(*,*) " Calling CHAMELEON_ZHEGST"
  call CHAMELEON_ZHEGST(ITYPE,local_UPLO,N,A,LDA,B,LDB,INFO)
end subroutine chameleon_wrap_ZHEGST

subroutine chameleon_wrap_ZLASWP(N,A,LDA,K1,K2,IPIV,INCX)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: INCX
  integer, intent(in) :: K1
  integer, intent(in) :: K2
  integer, intent(in) :: LDA
  integer, intent(in) :: N
  integer, intent(in), target :: IPIV(*)
  complex(kind=wp), intent(inout), target :: A(LDA,*)
  integer :: local_ret
  if (.not. chameleon_initialized) call chameleon_init(24,local_ret)
  ! write(*,*) " Calling CHAMELEON_ZLASWP"
  call CHAMELEON_ZLASWP(N,A,LDA,K1,K2,IPIV,INCX,local_ret)
end subroutine chameleon_wrap_ZLASWP

subroutine chameleon_wrap_ZLAUUM(UPLO,N,A,LDA,INFO)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: N
  integer, intent(out) :: INFO
  character, intent(in) :: UPLO
  complex(kind=wp), intent(inout), target :: A(LDA,*)
  integer :: local_UPLO
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,INFO)
  ! write(*,*) " Calling CHAMELEON_ZLAUUM"
  call CHAMELEON_ZLAUUM(local_UPLO,N,A,LDA,INFO)
end subroutine chameleon_wrap_ZLAUUM

subroutine chameleon_wrap_ZPOSV(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: LDB
  integer, intent(in) :: N
  integer, intent(in) :: NRHS
  integer, intent(out) :: INFO
  character, intent(in) :: UPLO
  complex(kind=wp), intent(inout), target :: A(LDA,*)
  complex(kind=wp), intent(inout), target :: B(LDB,*)
  integer :: local_UPLO
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,INFO)
  ! write(*,*) " Calling CHAMELEON_ZPOSV"
  call CHAMELEON_ZPOSV(local_UPLO,N,NRHS,A,LDA,B,LDB,INFO)
end subroutine chameleon_wrap_ZPOSV

subroutine chameleon_wrap_ZSYSV(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: LDB
  integer, intent(in) :: N
  integer, intent(in) :: NRHS
  integer, intent(out) :: INFO
  character, intent(in) :: UPLO
  complex(kind=wp), intent(inout), target :: A(LDA,*)
  complex(kind=wp), intent(inout), target :: B(LDB,*)
  integer :: local_UPLO
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,INFO)
  ! write(*,*) " Calling CHAMELEON_ZSYSV"
  call CHAMELEON_ZSYSV(local_UPLO,N,NRHS,A,LDA,B,LDB,INFO)
end subroutine chameleon_wrap_ZSYSV

subroutine chameleon_wrap_ZPOTRF(UPLO,N,A,LDA,INFO)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: N
  integer, intent(out) :: INFO
  character, intent(in) :: UPLO
  complex(kind=wp), intent(inout), target :: A(LDA,*)
  integer :: local_UPLO
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,INFO)
  ! write(*,*) " Calling CHAMELEON_ZPOTRF"
  call CHAMELEON_ZPOTRF(local_UPLO,N,A,LDA,INFO)
end subroutine chameleon_wrap_ZPOTRF

subroutine chameleon_wrap_ZSYTRF(UPLO,N,A,LDA,INFO)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: N
  integer, intent(out) :: INFO
  character, intent(in) :: UPLO
  complex(kind=wp), intent(inout), target :: A(LDA,*)
  integer :: local_UPLO
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,INFO)
  ! write(*,*) " Calling CHAMELEON_ZSYTRF"
  call CHAMELEON_ZSYTRF(local_UPLO,N,A,LDA,INFO)
end subroutine chameleon_wrap_ZSYTRF

subroutine chameleon_wrap_ZPOTRI(UPLO,N,A,LDA,INFO)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: N
  integer, intent(out) :: INFO
  character, intent(in) :: UPLO
  complex(kind=wp), intent(inout), target :: A(LDA,*)
  integer :: local_UPLO
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,INFO)
  ! write(*,*) " Calling CHAMELEON_ZPOTRI"
  call CHAMELEON_ZPOTRI(local_UPLO,N,A,LDA,INFO)
end subroutine chameleon_wrap_ZPOTRI

subroutine chameleon_wrap_ZPOTRS(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: LDB
  integer, intent(in) :: N
  integer, intent(in) :: NRHS
  integer, intent(out) :: INFO
  character, intent(in) :: UPLO
  complex(kind=wp), intent(in), target :: A(LDA,*)
  complex(kind=wp), intent(inout), target :: B(LDB,*)
  integer :: local_UPLO
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,INFO)
  ! write(*,*) " Calling CHAMELEON_ZPOTRS"
  call CHAMELEON_ZPOTRS(local_UPLO,N,NRHS,A,LDA,B,LDB,INFO)
end subroutine chameleon_wrap_ZPOTRS

subroutine chameleon_wrap_ZSYTRS(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: LDB
  integer, intent(in) :: N
  integer, intent(in) :: NRHS
  integer, intent(out) :: INFO
  character, intent(in) :: UPLO
  complex(kind=wp), intent(in), target :: A(LDA,*)
  complex(kind=wp), intent(inout), target :: B(LDB,*)
  integer :: local_UPLO
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,INFO)
  ! write(*,*) " Calling CHAMELEON_ZSYTRS"
  call CHAMELEON_ZSYTRS(local_UPLO,N,NRHS,A,LDA,B,LDB,INFO)
end subroutine chameleon_wrap_ZSYTRS

subroutine chameleon_wrap_ZTRTRI(UPLO,DIAG,N,A,LDA,INFO)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: N
  integer, intent(out) :: INFO
  character, intent(in) :: DIAG
  character, intent(in) :: UPLO
  complex(kind=wp), intent(inout), target :: A(LDA,*)
  integer :: local_DIAG
  integer :: local_UPLO
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if(DIAG=='U' .or. DIAG=='u')then
     local_DIAG = ChamUnit
  else if(DIAG=='N' .or. DIAG=='n')then
     local_DIAG = ChamNonUnit
  else
     local_DIAG = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,INFO)
  ! write(*,*) " Calling CHAMELEON_ZTRTRI"
  call CHAMELEON_ZTRTRI(local_UPLO,local_DIAG,N,A,LDA,INFO)
end subroutine chameleon_wrap_ZTRTRI

subroutine chameleon_wrap_ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: K
  integer, intent(in) :: LDA
  integer, intent(in) :: LDB
  integer, intent(in) :: LDC
  integer, intent(in) :: M
  integer, intent(in) :: N
  character, intent(in) :: TRANSA
  character, intent(in) :: TRANSB
  complex(kind=wp), intent(in) :: ALPHA
  complex(kind=wp), intent(in) :: BETA
  complex(kind=wp), intent(in), target :: A(LDA,*)
  complex(kind=wp), intent(in), target :: B(LDB,*)
  complex(kind=wp), intent(inout), target :: C(LDC,*)
  integer :: local_TRANSA
  integer :: local_TRANSB
  integer :: local_ret
  if(TRANSA=='N' .or. TRANSA=='n')then
     local_TRANSA = ChamNoTrans
  else if(TRANSA=='T' .or. TRANSA=='t')then
     local_TRANSA = ChamTrans
  else if(TRANSA=='C' .or. TRANSA=='c')then
     local_TRANSA = ChamConjTrans
  else
     local_TRANSA = -1
  end if
  if(TRANSB=='N' .or. TRANSB=='n')then
     local_TRANSB = ChamNoTrans
  else if(TRANSB=='T' .or. TRANSB=='t')then
     local_TRANSB = ChamTrans
  else if(TRANSB=='C' .or. TRANSB=='c')then
     local_TRANSB = ChamConjTrans
  else
     local_TRANSB = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,local_ret)
  ! write(*,*) " Calling CHAMELEON_ZGEMM"
  call CHAMELEON_ZGEMM(local_TRANSA,local_TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
end subroutine chameleon_wrap_ZGEMM

#if defined(PRECISION_z) || defined(PRECISION_c)
subroutine chameleon_wrap_ZHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: LDB
  integer, intent(in) :: LDC
  integer, intent(in) :: M
  integer, intent(in) :: N
  character, intent(in) :: SIDE
  character, intent(in) :: UPLO
  complex(kind=wp), intent(in) :: ALPHA
  complex(kind=wp), intent(in) :: BETA
  complex(kind=wp), intent(in), target :: A(LDA,*)
  complex(kind=wp), intent(in), target :: B(LDB,*)
  complex(kind=wp), intent(inout), target :: C(LDC,*)
  integer :: local_SIDE
  integer :: local_UPLO
  integer :: local_ret
  if(SIDE=='L' .or. SIDE=='l')then
     local_SIDE = ChamLeft
  else if(SIDE=='R' .or. SIDE=='r')then
     local_SIDE = ChamRight
  else
     local_SIDE = -1
  end if
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,local_ret)
  ! write(*,*) " Calling CHAMELEON_ZHEMM"
  call CHAMELEON_ZHEMM(local_SIDE,local_UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
end subroutine chameleon_wrap_ZHEMM

subroutine chameleon_wrap_ZHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: K
  integer, intent(in) :: LDA
  integer, intent(in) :: LDB
  integer, intent(in) :: LDC
  integer, intent(in) :: N
  character, intent(in) :: TRANS
  character, intent(in) :: UPLO
  complex(kind=wp), intent(in) :: ALPHA
  complex(kind=wp), intent(in), target :: A(LDA,*)
  complex(kind=wp), intent(in), target :: B(LDB,*)
  complex(kind=wp), intent(inout), target :: C(LDC,*)
  double precision, intent(in) :: BETA
  integer :: local_TRANS
  integer :: local_UPLO
  integer :: local_ret
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if(TRANS=='N' .or. TRANS=='n')then
     local_TRANS = ChamNoTrans
  else if(TRANS=='T' .or. TRANS=='t')then
     local_TRANS = ChamTrans
  else if(TRANS=='C' .or. TRANS=='c')then
     local_TRANS = ChamConjTrans
  else
     local_TRANS = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,local_ret)
  ! write(*,*) " Calling CHAMELEON_ZHER2K"
  call CHAMELEON_ZHER2K(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
end subroutine chameleon_wrap_ZHER2K

subroutine chameleon_wrap_ZHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: K
  integer, intent(in) :: LDA
  integer, intent(in) :: LDC
  integer, intent(in) :: N
  character, intent(in) :: TRANS
  character, intent(in) :: UPLO
  complex(kind=wp), intent(in), target :: A(LDA,*)
  complex(kind=wp), intent(inout), target :: C(LDC,*)
  double precision, intent(in) :: ALPHA
  double precision, intent(in) :: BETA
  integer :: local_TRANS
  integer :: local_UPLO
  integer :: local_ret
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if(TRANS=='N' .or. TRANS=='n')then
     local_TRANS = ChamNoTrans
  else if(TRANS=='T' .or. TRANS=='t')then
     local_TRANS = ChamTrans
  else if(TRANS=='C' .or. TRANS=='c')then
     local_TRANS = ChamConjTrans
  else
     local_TRANS = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,local_ret)
  ! write(*,*) " Calling CHAMELEON_ZHERK"
  call CHAMELEON_ZHERK(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC,local_ret)
end subroutine chameleon_wrap_ZHERK
#endif

subroutine chameleon_wrap_ZSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: LDB
  integer, intent(in) :: LDC
  integer, intent(in) :: M
  integer, intent(in) :: N
  character, intent(in) :: SIDE
  character, intent(in) :: UPLO
  complex(kind=wp), intent(in) :: ALPHA
  complex(kind=wp), intent(in) :: BETA
  complex(kind=wp), intent(in), target :: A(LDA,*)
  complex(kind=wp), intent(in), target :: B(LDB,*)
  complex(kind=wp), intent(inout), target :: C(LDC,*)
  integer :: local_SIDE
  integer :: local_UPLO
  integer :: local_ret
  if(SIDE=='L' .or. SIDE=='l')then
     local_SIDE = ChamLeft
  else if(SIDE=='R' .or. SIDE=='r')then
     local_SIDE = ChamRight
  else
     local_SIDE = -1
  end if
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,local_ret)
  ! write(*,*) " Calling CHAMELEON_ZSYMM"
  call CHAMELEON_ZSYMM(local_SIDE,local_UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
end subroutine chameleon_wrap_ZSYMM

subroutine chameleon_wrap_ZSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: K
  integer, intent(in) :: LDA
  integer, intent(in) :: LDB
  integer, intent(in) :: LDC
  integer, intent(in) :: N
  character, intent(in) :: TRANS
  character, intent(in) :: UPLO
  complex(kind=wp), intent(in) :: ALPHA
  complex(kind=wp), intent(in) :: BETA
  complex(kind=wp), intent(in), target :: A(LDA,*)
  complex(kind=wp), intent(in), target :: B(LDB,*)
  complex(kind=wp), intent(inout), target :: C(LDC,*)
  integer :: local_TRANS
  integer :: local_UPLO
  integer :: local_ret
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if(TRANS=='N' .or. TRANS=='n')then
     local_TRANS = ChamNoTrans
  else if(TRANS=='T' .or. TRANS=='t')then
     local_TRANS = ChamTrans
  else if(TRANS=='C' .or. TRANS=='c')then
     local_TRANS = ChamConjTrans
  else
     local_TRANS = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,local_ret)
  ! write(*,*) " Calling CHAMELEON_ZSYR2K"
  call CHAMELEON_ZSYR2K(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
end subroutine chameleon_wrap_ZSYR2K

subroutine chameleon_wrap_ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: K
  integer, intent(in) :: LDA
  integer, intent(in) :: LDC
  integer, intent(in) :: N
  character, intent(in) :: TRANS
  character, intent(in) :: UPLO
  complex(kind=wp), intent(in) :: ALPHA
  complex(kind=wp), intent(in) :: BETA
  complex(kind=wp), intent(in), target :: A(LDA,*)
  complex(kind=wp), intent(inout), target :: C(LDC,*)
  integer :: local_TRANS
  integer :: local_UPLO
  integer :: local_ret
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if(TRANS=='N' .or. TRANS=='n')then
     local_TRANS = ChamNoTrans
  else if(TRANS=='T' .or. TRANS=='t')then
     local_TRANS = ChamTrans
  else if(TRANS=='C' .or. TRANS=='c')then
     local_TRANS = ChamConjTrans
  else
     local_TRANS = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,local_ret)
  ! write(*,*) " Calling CHAMELEON_ZSYRK"
  call CHAMELEON_ZSYRK(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC,local_ret)
end subroutine chameleon_wrap_ZSYRK

subroutine chameleon_wrap_ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: LDB
  integer, intent(in) :: M
  integer, intent(in) :: N
  character, intent(in) :: DIAG
  character, intent(in) :: SIDE
  character, intent(in) :: TRANSA
  character, intent(in) :: UPLO
  complex(kind=wp), intent(in) :: ALPHA
  complex(kind=wp), intent(in), target :: A(LDA,*)
  complex(kind=wp), intent(inout), target :: B(LDB,*)
  integer :: local_DIAG
  integer :: local_SIDE
  integer :: local_TRANSA
  integer :: local_UPLO
  integer :: local_ret
  if(SIDE=='L' .or. SIDE=='l')then
     local_SIDE = ChamLeft
  else if(SIDE=='R' .or. SIDE=='r')then
     local_SIDE = ChamRight
  else
     local_SIDE = -1
  end if
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if(TRANSA=='N' .or. TRANSA=='n')then
     local_TRANSA = ChamNoTrans
  else if(TRANSA=='T' .or. TRANSA=='t')then
     local_TRANSA = ChamTrans
  else if(TRANSA=='C' .or. TRANSA=='c')then
     local_TRANSA = ChamConjTrans
  else
     local_TRANSA = -1
  end if
  if(DIAG=='U' .or. DIAG=='u')then
     local_DIAG = ChamUnit
  else if(DIAG=='N' .or. DIAG=='n')then
     local_DIAG = ChamNonUnit
  else
     local_DIAG = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,local_ret)
  ! write(*,*) " Calling CHAMELEON_ZTRMM"
  call CHAMELEON_ZTRMM(local_SIDE,local_UPLO,local_TRANSA,local_DIAG,M,N,ALPHA,A,LDA,B,LDB,local_ret)
end subroutine chameleon_wrap_ZTRMM

subroutine chameleon_wrap_ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: LDB
  integer, intent(in) :: M
  integer, intent(in) :: N
  character, intent(in) :: DIAG
  character, intent(in) :: SIDE
  character, intent(in) :: TRANSA
  character, intent(in) :: UPLO
  complex(kind=wp), intent(in) :: ALPHA
  complex(kind=wp), intent(in), target :: A(LDA,*)
  complex(kind=wp), intent(inout), target :: B(LDB,*)
  integer :: local_DIAG
  integer :: local_SIDE
  integer :: local_TRANSA
  integer :: local_UPLO
  integer :: local_ret
  if(SIDE=='L' .or. SIDE=='l')then
     local_SIDE = ChamLeft
  else if(SIDE=='R' .or. SIDE=='r')then
     local_SIDE = ChamRight
  else
     local_SIDE = -1
  end if
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if(TRANSA=='N' .or. TRANSA=='n')then
     local_TRANSA = ChamNoTrans
  else if(TRANSA=='T' .or. TRANSA=='t')then
     local_TRANSA = ChamTrans
  else if(TRANSA=='C' .or. TRANSA=='c')then
     local_TRANSA = ChamConjTrans
  else
     local_TRANSA = -1
  end if
  if(DIAG=='U' .or. DIAG=='u')then
     local_DIAG = ChamUnit
  else if(DIAG=='N' .or. DIAG=='n')then
     local_DIAG = ChamNonUnit
  else
     local_DIAG = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,local_ret)
  ! write(*,*) " Calling CHAMELEON_ZTRSM"
  call CHAMELEON_ZTRSM(local_SIDE,local_UPLO,local_TRANSA,local_DIAG,M,N,ALPHA,A,LDA,B,LDB,local_ret)
end subroutine chameleon_wrap_ZTRSM

subroutine chameleon_wrap_ZLACPY(UPLO,M,N,A,LDA,B,LDB)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: LDB
  integer, intent(in) :: M
  integer, intent(in) :: N
  character, intent(in) :: UPLO
  complex(kind=wp), intent(in), target :: A(LDA,*)
  complex(kind=wp), intent(out), target :: B(LDB,*)
  integer :: local_UPLO
  integer :: local_ret
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,local_ret)
  ! write(*,*) " Calling CHAMELEON_ZLACPY"
  call CHAMELEON_ZLACPY(local_UPLO,M,N,A,LDA,B,LDB,local_ret)
end subroutine chameleon_wrap_ZLACPY

subroutine chameleon_wrap_ZLASET(UPLO,M,N,ALPHA,BETA,A,LDA)
  use iso_c_binding
  use chameleon
  implicit none
  integer, parameter :: wp = kind(0.0d0)
  integer, intent(in) :: LDA
  integer, intent(in) :: M
  integer, intent(in) :: N
  character, intent(in) :: UPLO
  complex(kind=wp), intent(in) :: ALPHA
  complex(kind=wp), intent(in) :: BETA
  complex(kind=wp), intent(inout), target :: A(LDA,*)
  integer :: local_UPLO
  integer :: local_ret
  if(UPLO=='U' .or. UPLO=='u')then
     local_UPLO = ChamUpper
  else if(UPLO=='L' .or. UPLO=='l')then
     local_UPLO = ChamLower
  else
     local_UPLO = -1
  end if
  if (.not. chameleon_initialized) call chameleon_init(24,local_ret)
  ! write(*,*) " Calling CHAMELEON_ZLASET"
  call CHAMELEON_ZLASET(local_UPLO,M,N,ALPHA,BETA,A,LDA,local_ret)
end subroutine chameleon_wrap_ZLASET
