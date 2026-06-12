!!!
!
! @file chameleon_zf90.F90
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
! @precisions normal z -> c d s
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
#define PRECISION_z

module chameleon_z
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !  FORTRAN API - math functions (simple interface)
  !
  interface
     function CHAMELEON_zgebrd_c(M,N,A,LDA,D,E,descT) &
          & bind(c, name='CHAMELEON_zgebrd')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgebrd_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: D
       type(c_ptr), value :: E
       type(c_ptr), value :: descT
     end function CHAMELEON_zgebrd_c
  end interface

  interface
     function CHAMELEON_zgecfi_c(m,n,A,fin,imb,inb,fout,omb,onb) &
          & bind(c, name='CHAMELEON_zgecfi')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgecfi_c
       integer(kind=c_int), value :: m
       integer(kind=c_int), value :: n
       type(c_ptr), value :: A
       integer(kind=c_int), value :: fin
       integer(kind=c_int), value :: imb
       integer(kind=c_int), value :: inb
       integer(kind=c_int), value :: fout
       integer(kind=c_int), value :: omb
       integer(kind=c_int), value :: onb
     end function CHAMELEON_zgecfi_c
  end interface

  interface
     function CHAMELEON_zgecfi_Async_c(m,n,A,f_in,imb,inb,f_out,omb,onb,sequence,request) &
          & bind(c, name='CHAMELEON_zgecfi_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgecfi_Async_c
       integer(kind=c_int), value :: m
       integer(kind=c_int), value :: n
       type(c_ptr), value :: A
       integer(kind=c_int), value :: f_in
       integer(kind=c_int), value :: imb
       integer(kind=c_int), value :: inb
       integer(kind=c_int), value :: f_out
       integer(kind=c_int), value :: omb
       integer(kind=c_int), value :: onb
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgecfi_Async_c
  end interface

  interface
     function CHAMELEON_zgelqf_c(M,N,A,LDA,descT) &
          & bind(c, name='CHAMELEON_zgelqf')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgelqf_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: descT
     end function CHAMELEON_zgelqf_c
  end interface

  interface
     function CHAMELEON_zgelqs_c(M,N,NRHS,A,LDA,descT,B,LDB) &
          & bind(c, name='CHAMELEON_zgelqs')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgelqs_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: NRHS
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: descT
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zgelqs_c
  end interface

  interface
     function CHAMELEON_zgels_c(trans,M,N,NRHS,A,LDA,descT,B,LDB) &
          & bind(c, name='CHAMELEON_zgels')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgels_c
       integer(kind=c_int), value :: trans
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: NRHS
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: descT
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zgels_c
  end interface

  interface
     function CHAMELEON_zgemm_c(transA,transB,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='CHAMELEON_zgemm')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgemm_c
       integer(kind=c_int), value :: transA
       integer(kind=c_int), value :: transB
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: K
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: C
       integer(kind=c_int), value :: LDC
     end function CHAMELEON_zgemm_c
  end interface

  interface
     function CHAMELEON_zgeqrf_c(M,N,A,LDA,descT) &
          & bind(c, name='CHAMELEON_zgeqrf')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgeqrf_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: descT
     end function CHAMELEON_zgeqrf_c
  end interface

  interface
     function CHAMELEON_zgeqrs_c(M,N,NRHS,A,LDA,descT,B,LDB) &
          & bind(c, name='CHAMELEON_zgeqrs')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgeqrs_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: NRHS
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: descT
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zgeqrs_c
  end interface

  interface
     function CHAMELEON_zgesv_c(N,NRHS,A,LDA,IPIV,B,LDB) &
          & bind(c, name='CHAMELEON_zgesv')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgesv_c
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: NRHS
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zgesv_c
  end interface

  interface
     function CHAMELEON_zgesv_incpiv_c(N,NRHS,A,LDA,descL,IPIV,B,LDB) &
          & bind(c, name='CHAMELEON_zgesv_incpiv')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgesv_incpiv_c
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: NRHS
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: descL
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zgesv_incpiv_c
  end interface

  interface
     function CHAMELEON_zgesvd_c(jobu,jobvt,M,N,A,LDA,S,U,LDU,VT,LDVT,descT) &
          & bind(c, name='CHAMELEON_zgesvd')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgesvd_c
       integer(kind=c_int), value :: jobu
       integer(kind=c_int), value :: jobvt
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: S
       type(c_ptr), value :: U
       integer(kind=c_int), value :: LDU
       type(c_ptr), value :: VT
       integer(kind=c_int), value :: LDVT
       type(c_ptr), value :: descT
     end function CHAMELEON_zgesvd_c
  end interface

  interface
     function CHAMELEON_zgetmi_c(m,n,A,fin,mb,nb) &
          & bind(c, name='CHAMELEON_zgetmi')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetmi_c
       integer(kind=c_int), value :: m
       integer(kind=c_int), value :: n
       type(c_ptr), value :: A
       integer(kind=c_int), value :: fin
       integer(kind=c_int), value :: mb
       integer(kind=c_int), value :: nb
     end function CHAMELEON_zgetmi_c
  end interface

  interface
     function CHAMELEON_zgetmi_Async_c(m,n,A,f_in,mb,inb,sequence,request) &
          & bind(c, name='CHAMELEON_zgetmi_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetmi_Async_c
       integer(kind=c_int), value :: m
       integer(kind=c_int), value :: n
       type(c_ptr), value :: A
       integer(kind=c_int), value :: f_in
       integer(kind=c_int), value :: mb
       integer(kind=c_int), value :: inb
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgetmi_Async_c
  end interface

  interface
     function CHAMELEON_zgetrf_c(M,N,A,LDA,IPIV) &
          & bind(c, name='CHAMELEON_zgetrf')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetrf_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: IPIV
     end function CHAMELEON_zgetrf_c
  end interface

  interface
     function CHAMELEON_zgetrf_incpiv_c(M,N,A,LDA,descL,IPIV) &
          & bind(c, name='CHAMELEON_zgetrf_incpiv')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetrf_incpiv_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: descL
       type(c_ptr), value :: IPIV
     end function CHAMELEON_zgetrf_incpiv_c
  end interface

  interface
     function CHAMELEON_zgetrf_nopiv_c(M,N,A,LDA) &
          & bind(c, name='CHAMELEON_zgetrf_nopiv')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetrf_nopiv_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
     end function CHAMELEON_zgetrf_nopiv_c
  end interface

  interface
     function CHAMELEON_zgetri_c(N,A,LDA,IPIV) &
          & bind(c, name='CHAMELEON_zgetri')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetri_c
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: IPIV
     end function CHAMELEON_zgetri_c
  end interface

  interface
     function CHAMELEON_zgetrs_c(trans,N,NRHS,A,LDA,IPIV,B,LDB) &
          & bind(c, name='CHAMELEON_zgetrs')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetrs_c
       integer(kind=c_int), value :: trans
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: NRHS
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zgetrs_c
  end interface

  interface
     function CHAMELEON_zgetrs_incpiv_c(trans,N,NRHS,A,LDA,descL,IPIV,B,LDB) &
          & bind(c, name='CHAMELEON_zgetrs_incpiv')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetrs_incpiv_c
       integer(kind=c_int), value :: trans
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: NRHS
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: descL
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zgetrs_incpiv_c
  end interface

  interface
     function CHAMELEON_zheev_c(jobz,uplo,N,A,LDA,W,descT,Q,LDQ) &
          & bind(c, name='CHAMELEON_zheev')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zheev_c
       integer(kind=c_int), value :: jobz
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: W
       type(c_ptr), value :: descT
       type(c_ptr), value :: Q
       integer(kind=c_int), value :: LDQ
     end function CHAMELEON_zheev_c
  end interface

  interface
     function CHAMELEON_zheevd_c(jobz,uplo,N,A,LDA,W,descT,Q,LDQ) &
          & bind(c, name='CHAMELEON_zheevd')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zheevd_c
       integer(kind=c_int), value :: jobz
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: W
       type(c_ptr), value :: descT
       type(c_ptr), value :: Q
       integer(kind=c_int), value :: LDQ
     end function CHAMELEON_zheevd_c
  end interface

  interface
     function CHAMELEON_zhegst_c(itype,uplo,N,A,LDA,B,LDB) &
          & bind(c, name='CHAMELEON_zhegst')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zhegst_c
       integer(kind=c_int), value :: itype
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zhegst_c
  end interface

  interface
     function CHAMELEON_zhegv_c(itype,jobz,uplo,N,A,LDA,B,LDB,W,descT,Q,LDQ) &
          & bind(c, name='CHAMELEON_zhegv')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zhegv_c
       integer(kind=c_int), value :: itype
       integer(kind=c_int), value :: jobz
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
       type(c_ptr), value :: W
       type(c_ptr), value :: descT
       type(c_ptr), value :: Q
       integer(kind=c_int), value :: LDQ
     end function CHAMELEON_zhegv_c
  end interface

  interface
     function CHAMELEON_zhegvd_c(itype,jobz,uplo,N,A,LDA,B,LDB,W,descT,Q,LDQ) &
          & bind(c, name='CHAMELEON_zhegvd')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zhegvd_c
       integer(kind=c_int), value :: itype
       integer(kind=c_int), value :: jobz
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
       type(c_ptr), value :: W
       type(c_ptr), value :: descT
       type(c_ptr), value :: Q
       integer(kind=c_int), value :: LDQ
     end function CHAMELEON_zhegvd_c
  end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
  interface
     function CHAMELEON_zhemm_c(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='CHAMELEON_zhemm')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zhemm_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: C
       integer(kind=c_int), value :: LDC
     end function CHAMELEON_zhemm_c
  end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
  interface
     function CHAMELEON_zher2k_c(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='CHAMELEON_zher2k')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zher2k_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: trans
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: K
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
       real(kind=c_double), value :: beta
       type(c_ptr), value :: C
       integer(kind=c_int), value :: LDC
     end function CHAMELEON_zher2k_c
  end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
  interface
     function CHAMELEON_zherk_c(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC) &
          & bind(c, name='CHAMELEON_zherk')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zherk_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: trans
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: K
       real(kind=c_double), value :: alpha
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       real(kind=c_double), value :: beta
       type(c_ptr), value :: C
       integer(kind=c_int), value :: LDC
     end function CHAMELEON_zherk_c
  end interface
#endif

  interface
     function CHAMELEON_zhetrd_c(jobz,uplo,N,A,LDA,D,E,descT,Q,LDQ) &
          & bind(c, name='CHAMELEON_zhetrd')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zhetrd_c
       integer(kind=c_int), value :: jobz
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: D
       type(c_ptr), value :: E
       type(c_ptr), value :: descT
       type(c_ptr), value :: Q
       integer(kind=c_int), value :: LDQ
     end function CHAMELEON_zhetrd_c
  end interface

  interface
     function CHAMELEON_zlacpy_c(uplo,M,N,A,LDA,B,LDB) &
          & bind(c, name='CHAMELEON_zlacpy')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlacpy_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zlacpy_c
  end interface

  interface
     function CHAMELEON_zlange_c(norm,M,N,A,LDA,work) &
          & bind(c, name='CHAMELEON_zlange')
       use iso_c_binding
       implicit none
       real(kind=c_double) :: CHAMELEON_zlange_c
       integer(kind=c_int), value :: norm
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: work
     end function CHAMELEON_zlange_c
  end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
  interface
     function CHAMELEON_zlanhe_c(norm,uplo,N,A,LDA,work) &
          & bind(c, name='CHAMELEON_zlanhe')
       use iso_c_binding
       implicit none
       real(kind=c_double) :: CHAMELEON_zlanhe_c
       integer(kind=c_int), value :: norm
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: work
     end function CHAMELEON_zlanhe_c
  end interface
#endif

  interface
     function CHAMELEON_zlansy_c(norm,uplo,N,A,LDA,work) &
          & bind(c, name='CHAMELEON_zlansy')
       use iso_c_binding
       implicit none
       real(kind=c_double) :: CHAMELEON_zlansy_c
       integer(kind=c_int), value :: norm
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: work
     end function CHAMELEON_zlansy_c
  end interface

  interface
     function CHAMELEON_zlaset_c(uplo,M,N,alpha,beta,A,LDA) &
          & bind(c, name='CHAMELEON_zlaset')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlaset_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       complex(kind=c_double_complex), value :: alpha
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
     end function CHAMELEON_zlaset_c
  end interface

  interface
     function CHAMELEON_zlaswp_c(N,A,LDA,K1,K2,IPIV,INCX) &
          & bind(c, name='CHAMELEON_zlaswp')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlaswp_c
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       integer(kind=c_int), value :: K1
       integer(kind=c_int), value :: K2
       type(c_ptr), value :: IPIV
       integer(kind=c_int), value :: INCX
     end function CHAMELEON_zlaswp_c
  end interface

  interface
     function CHAMELEON_zlaswpc_c(N,A,LDA,K1,K2,IPIV,INCX) &
          & bind(c, name='CHAMELEON_zlaswpc')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlaswpc_c
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       integer(kind=c_int), value :: K1
       integer(kind=c_int), value :: K2
       type(c_ptr), value :: IPIV
       integer(kind=c_int), value :: INCX
     end function CHAMELEON_zlaswpc_c
  end interface

  interface
     function CHAMELEON_zlauum_c(uplo,N,A,LDA) &
          & bind(c, name='CHAMELEON_zlauum')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlauum_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
     end function CHAMELEON_zlauum_c
  end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
  interface
     function CHAMELEON_zplghe_c(bump,uplo,N,A,LDA,seed) &
          & bind(c, name='CHAMELEON_zplghe')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zplghe_c
       real(kind=c_double), value :: bump
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       integer(kind=c_long_long), value :: seed
     end function CHAMELEON_zplghe_c
  end interface
#endif

  interface
     function CHAMELEON_zplgsy_c(bump,uplo,N,A,LDA,seed) &
          & bind(c, name='CHAMELEON_zplgsy')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zplgsy_c
       complex(kind=c_double_complex), value :: bump
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       integer(kind=c_long_long), value :: seed
     end function CHAMELEON_zplgsy_c
  end interface

  interface
     function CHAMELEON_zplrnt_c(M,N,A,LDA,seed) &
          & bind(c, name='CHAMELEON_zplrnt')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zplrnt_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       integer(kind=c_long_long), value :: seed
     end function CHAMELEON_zplrnt_c
  end interface

  interface
     function CHAMELEON_zplrnk_c(M,N,K,C,LDC,seedA,seedB) &
          & bind(c, name='CHAMELEON_zplrnk')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zplrnk_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: K
       type(c_ptr), value :: C
       integer(kind=c_int), value :: LDC
       integer(kind=c_long_long), value :: seedA
       integer(kind=c_long_long), value :: seedB
     end function CHAMELEON_zplrnk_c
  end interface

  interface
     function CHAMELEON_zposv_c(uplo,N,NRHS,A,LDA,B,LDB) &
          & bind(c, name='CHAMELEON_zposv')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zposv_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: NRHS
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zposv_c
  end interface

  interface
     function CHAMELEON_zsysv_c(uplo,N,NRHS,A,LDA,B,LDB) &
          & bind(c, name='CHAMELEON_zsysv')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsysv_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: NRHS
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zsysv_c
  end interface

  interface
     function CHAMELEON_zpotrf_c(uplo,N,A,LDA) &
          & bind(c, name='CHAMELEON_zpotrf')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zpotrf_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
     end function CHAMELEON_zpotrf_c
  end interface

  interface
     function CHAMELEON_zsytrf_c(uplo,N,A,LDA) &
          & bind(c, name='CHAMELEON_zsytrf')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsytrf_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
     end function CHAMELEON_zsytrf_c
  end interface

  interface
     function CHAMELEON_zpotri_c(uplo,N,A,LDA) &
          & bind(c, name='CHAMELEON_zpotri')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zpotri_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
     end function CHAMELEON_zpotri_c
  end interface

  interface
     function CHAMELEON_zpotrs_c(uplo,N,NRHS,A,LDA,B,LDB) &
          & bind(c, name='CHAMELEON_zpotrs')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zpotrs_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: NRHS
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zpotrs_c
  end interface

  interface
     function CHAMELEON_zsytrs_c(uplo,N,NRHS,A,LDA,B,LDB) &
          & bind(c, name='CHAMELEON_zsytrs')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsytrs_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: NRHS
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zsytrs_c
  end interface

  interface
     function CHAMELEON_zsymm_c(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='CHAMELEON_zsymm')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsymm_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: C
       integer(kind=c_int), value :: LDC
     end function CHAMELEON_zsymm_c
  end interface

  interface
     function CHAMELEON_zsyr2k_c(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) &
          & bind(c, name='CHAMELEON_zsyr2k')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsyr2k_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: trans
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: K
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: C
       integer(kind=c_int), value :: LDC
     end function CHAMELEON_zsyr2k_c
  end interface

  interface
     function CHAMELEON_zsyrk_c(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC) &
          & bind(c, name='CHAMELEON_zsyrk')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsyrk_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: trans
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: K
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: C
       integer(kind=c_int), value :: LDC
     end function CHAMELEON_zsyrk_c
  end interface

  interface
     function CHAMELEON_ztrmm_c(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB) &
          & bind(c, name='CHAMELEON_ztrmm')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_ztrmm_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: transA
       integer(kind=c_int), value :: diag
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: NRHS
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_ztrmm_c
  end interface

  interface
     function CHAMELEON_ztrsm_c(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB) &
          & bind(c, name='CHAMELEON_ztrsm')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_ztrsm_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: transA
       integer(kind=c_int), value :: diag
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: NRHS
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_ztrsm_c
  end interface

  interface
     function CHAMELEON_ztrsmpl_c(N,NRHS,A,LDA,descL,IPIV,B,LDB) &
          & bind(c, name='CHAMELEON_ztrsmpl')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_ztrsmpl_c
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: NRHS
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: descL
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_ztrsmpl_c
  end interface

  interface
     function CHAMELEON_ztrsmrv_c(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB) &
          & bind(c, name='CHAMELEON_ztrsmrv')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_ztrsmrv_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: transA
       integer(kind=c_int), value :: diag
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: NRHS
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_ztrsmrv_c
  end interface

  interface
     function CHAMELEON_ztrtri_c(uplo,diag,N,A,LDA) &
          & bind(c, name='CHAMELEON_ztrtri')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_ztrtri_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: diag
       integer(kind=c_int), value :: N
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
     end function CHAMELEON_ztrtri_c
  end interface

  interface
     function CHAMELEON_zunglq_c(M,N,K,A,LDA,descT,B,LDB) &
          & bind(c, name='CHAMELEON_zunglq')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zunglq_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: K
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: descT
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zunglq_c
  end interface

  interface
     function CHAMELEON_zungqr_c(M,N,K,A,LDA,descT,B,LDB) &
          & bind(c, name='CHAMELEON_zungqr')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zungqr_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: K
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: descT
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zungqr_c
  end interface

  interface
     function CHAMELEON_zunmlq_c(side,trans,M,N,K,A,LDA,descT,B,LDB) &
          & bind(c, name='CHAMELEON_zunmlq')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zunmlq_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: trans
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: K
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: descT
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zunmlq_c
  end interface

  interface
     function CHAMELEON_zunmqr_c(side,trans,M,N,K,A,LDA,descT,B,LDB) &
          & bind(c, name='CHAMELEON_zunmqr')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zunmqr_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: trans
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       integer(kind=c_int), value :: K
       type(c_ptr), value :: A
       integer(kind=c_int), value :: LDA
       type(c_ptr), value :: descT
       type(c_ptr), value :: B
       integer(kind=c_int), value :: LDB
     end function CHAMELEON_zunmqr_c
  end interface

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !  FORTRAN API - math functions (native interface)
  !
  interface
     function CHAMELEON_zgebrd_Tile_c(A,D,E,T) &
          & bind(c, name='CHAMELEON_zgebrd_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgebrd_Tile_c
       type(c_ptr), value :: A
       type(c_ptr), value :: D
       type(c_ptr), value :: E
       type(c_ptr), value :: T
     end function CHAMELEON_zgebrd_Tile_c
  end interface

  interface
     function CHAMELEON_zgelqf_Tile_c(A,T) &
          & bind(c, name='CHAMELEON_zgelqf_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgelqf_Tile_c
       type(c_ptr), value :: A
       type(c_ptr), value :: T
     end function CHAMELEON_zgelqf_Tile_c
  end interface

  interface
     function CHAMELEON_zgelqs_Tile_c(A,T,B) &
          & bind(c, name='CHAMELEON_zgelqs_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgelqs_Tile_c
       type(c_ptr), value :: A
       type(c_ptr), value :: T
       type(c_ptr), value :: B
     end function CHAMELEON_zgelqs_Tile_c
  end interface

  interface
     function CHAMELEON_zgels_Tile_c(trans,A,T,B) &
          & bind(c, name='CHAMELEON_zgels_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgels_Tile_c
       integer(kind=c_int), value :: trans
       type(c_ptr), value :: A
       type(c_ptr), value :: T
       type(c_ptr), value :: B
     end function CHAMELEON_zgels_Tile_c
  end interface

  interface
     function CHAMELEON_zgemm_Tile_c(transA,transB,alpha,A,B,beta,C) &
          & bind(c, name='CHAMELEON_zgemm_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgemm_Tile_c
       integer(kind=c_int), value :: transA
       integer(kind=c_int), value :: transB
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: C
     end function CHAMELEON_zgemm_Tile_c
  end interface

  interface
     function CHAMELEON_zgeqrf_Tile_c(A,T) &
          & bind(c, name='CHAMELEON_zgeqrf_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgeqrf_Tile_c
       type(c_ptr), value :: A
       type(c_ptr), value :: T
     end function CHAMELEON_zgeqrf_Tile_c
  end interface

  interface
     function CHAMELEON_zgeqrs_Tile_c(A,T,B) &
          & bind(c, name='CHAMELEON_zgeqrs_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgeqrs_Tile_c
       type(c_ptr), value :: A
       type(c_ptr), value :: T
       type(c_ptr), value :: B
     end function CHAMELEON_zgeqrs_Tile_c
  end interface

  interface
     function CHAMELEON_zgesv_Tile_c(A,IPIV,B) &
          & bind(c, name='CHAMELEON_zgesv_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgesv_Tile_c
       type(c_ptr), value :: A
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: B
     end function CHAMELEON_zgesv_Tile_c
  end interface

  interface
     function CHAMELEON_zgesv_incpiv_Tile_c(A,L,IPIV,B) &
          & bind(c, name='CHAMELEON_zgesv_incpiv_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgesv_incpiv_Tile_c
       type(c_ptr), value :: A
       type(c_ptr), value :: L
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: B
     end function CHAMELEON_zgesv_incpiv_Tile_c
  end interface

  interface
     function CHAMELEON_zgesvd_Tile_c(jobu,jobvt,A,S,U,VT,T) &
          & bind(c, name='CHAMELEON_zgesvd_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgesvd_Tile_c
       integer(kind=c_int), value :: jobu
       integer(kind=c_int), value :: jobvt
       type(c_ptr), value :: A
       type(c_ptr), value :: S
       type(c_ptr), value :: U
       type(c_ptr), value :: VT
       type(c_ptr), value :: T
     end function CHAMELEON_zgesvd_Tile_c
  end interface

  interface
     function CHAMELEON_zgetrf_Tile_c(A,IPIV) &
          & bind(c, name='CHAMELEON_zgetrf_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetrf_Tile_c
       type(c_ptr), value :: A
       type(c_ptr), value :: IPIV
     end function CHAMELEON_zgetrf_Tile_c
  end interface

  interface
     function CHAMELEON_zgetrf_incpiv_Tile_c(A,L,IPIV) &
          & bind(c, name='CHAMELEON_zgetrf_incpiv_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetrf_incpiv_Tile_c
       type(c_ptr), value :: A
       type(c_ptr), value :: L
       type(c_ptr), value :: IPIV
     end function CHAMELEON_zgetrf_incpiv_Tile_c
  end interface

  interface
     function CHAMELEON_zgetrf_nopiv_Tile_c(A) &
          & bind(c, name='CHAMELEON_zgetrf_nopiv_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetrf_nopiv_Tile_c
       type(c_ptr), value :: A
     end function CHAMELEON_zgetrf_nopiv_Tile_c
  end interface

  interface
     function CHAMELEON_zgetri_Tile_c(A,IPIV) &
          & bind(c, name='CHAMELEON_zgetri_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetri_Tile_c
       type(c_ptr), value :: A
       type(c_ptr), value :: IPIV
     end function CHAMELEON_zgetri_Tile_c
  end interface

  interface
     function CHAMELEON_zgetrs_Tile_c(trans,A,IPIV,B) &
          & bind(c, name='CHAMELEON_zgetrs_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetrs_Tile_c
       integer(kind=c_int), value :: trans
       type(c_ptr), value :: A
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: B
     end function CHAMELEON_zgetrs_Tile_c
  end interface

  interface
     function CHAMELEON_zgetrs_incpiv_Tile_c(A,L,IPIV,B) &
          & bind(c, name='CHAMELEON_zgetrs_incpiv_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetrs_incpiv_Tile_c
       type(c_ptr), value :: A
       type(c_ptr), value :: L
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: B
     end function CHAMELEON_zgetrs_incpiv_Tile_c
  end interface

  interface
     function CHAMELEON_zheev_Tile_c(jobz,uplo,A,W,T,Q,LDQ) &
          & bind(c, name='CHAMELEON_zheev_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zheev_Tile_c
       integer(kind=c_int), value :: jobz
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: W
       type(c_ptr), value :: T
       type(c_ptr), value :: Q
       integer(kind=c_int), value :: LDQ
     end function CHAMELEON_zheev_Tile_c
  end interface

  interface
     function CHAMELEON_zheevd_Tile_c(jobz,uplo,A,W,T,Q,LDQ) &
          & bind(c, name='CHAMELEON_zheevd_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zheevd_Tile_c
       integer(kind=c_int), value :: jobz
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: W
       type(c_ptr), value :: T
       type(c_ptr), value :: Q
       integer(kind=c_int), value :: LDQ
     end function CHAMELEON_zheevd_Tile_c
  end interface

  interface
     function CHAMELEON_zhegst_Tile_c(itype,uplo,A,B) &
          & bind(c, name='CHAMELEON_zhegst_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zhegst_Tile_c
       integer(kind=c_int), value :: itype
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: B
     end function CHAMELEON_zhegst_Tile_c
  end interface

  interface
     function CHAMELEON_zhegv_Tile_c(itype,jobz,uplo,A,B,W,T,Q) &
          & bind(c, name='CHAMELEON_zhegv_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zhegv_Tile_c
       integer(kind=c_int), value :: itype
       integer(kind=c_int), value :: jobz
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       type(c_ptr), value :: W
       type(c_ptr), value :: T
       type(c_ptr), value :: Q
     end function CHAMELEON_zhegv_Tile_c
  end interface

  interface
     function CHAMELEON_zhegvd_Tile_c(itype,jobz,uplo,A,B,W,T,Q) &
          & bind(c, name='CHAMELEON_zhegvd_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zhegvd_Tile_c
       integer(kind=c_int), value :: itype
       integer(kind=c_int), value :: jobz
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       type(c_ptr), value :: W
       type(c_ptr), value :: T
       type(c_ptr), value :: Q
     end function CHAMELEON_zhegvd_Tile_c
  end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
  interface
     function CHAMELEON_zhemm_Tile_c(side,uplo,alpha,A,B,beta,C) &
          & bind(c, name='CHAMELEON_zhemm_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zhemm_Tile_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: uplo
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: C
     end function CHAMELEON_zhemm_Tile_c
  end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
  interface
     function CHAMELEON_zher2k_Tile_c(uplo,trans,alpha,A,B,beta,C) &
          & bind(c, name='CHAMELEON_zher2k_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zher2k_Tile_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: trans
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       real(kind=c_double), value :: beta
       type(c_ptr), value :: C
     end function CHAMELEON_zher2k_Tile_c
  end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
  interface
     function CHAMELEON_zherk_Tile_c(uplo,trans,alpha,A,beta,C) &
          & bind(c, name='CHAMELEON_zherk_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zherk_Tile_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: trans
       real(kind=c_double), value :: alpha
       type(c_ptr), value :: A
       real(kind=c_double), value :: beta
       type(c_ptr), value :: C
     end function CHAMELEON_zherk_Tile_c
  end interface
#endif

  interface
     function CHAMELEON_zhetrd_Tile_c(jobz,uplo,A,D,E,T,Q,LDQ) &
          & bind(c, name='CHAMELEON_zhetrd_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zhetrd_Tile_c
       integer(kind=c_int), value :: jobz
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: D
       type(c_ptr), value :: E
       type(c_ptr), value :: T
       type(c_ptr), value :: Q
       integer(kind=c_int), value :: LDQ
     end function CHAMELEON_zhetrd_Tile_c
  end interface

  interface
     function CHAMELEON_zlacpy_Tile_c(uplo,A,B) &
          & bind(c, name='CHAMELEON_zlacpy_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlacpy_Tile_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: B
     end function CHAMELEON_zlacpy_Tile_c
  end interface

  interface
     function CHAMELEON_zlange_Tile_c(norm,A,work) &
          & bind(c, name='CHAMELEON_zlange_Tile')
       use iso_c_binding
       implicit none
       real(kind=c_double) :: CHAMELEON_zlange_Tile_c
       integer(kind=c_int), value :: norm
       type(c_ptr), value :: A
       type(c_ptr), value :: work
     end function CHAMELEON_zlange_Tile_c
  end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
  interface
     function CHAMELEON_zlanhe_Tile_c(norm,uplo,A,work) &
          & bind(c, name='CHAMELEON_zlanhe_Tile')
       use iso_c_binding
       implicit none
       real(kind=c_double) :: CHAMELEON_zlanhe_Tile_c
       integer(kind=c_int), value :: norm
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: work
     end function CHAMELEON_zlanhe_Tile_c
  end interface
#endif

  interface
     function CHAMELEON_zlansy_Tile_c(norm,uplo,A,work) &
          & bind(c, name='CHAMELEON_zlansy_Tile')
       use iso_c_binding
       implicit none
       real(kind=c_double) :: CHAMELEON_zlansy_Tile_c
       integer(kind=c_int), value :: norm
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: work
     end function CHAMELEON_zlansy_Tile_c
  end interface

  interface
     function CHAMELEON_zlaset_Tile_c(uplo,alpha,beta,A) &
          & bind(c, name='CHAMELEON_zlaset_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlaset_Tile_c
       integer(kind=c_int), value :: uplo
       complex(kind=c_double_complex), value :: alpha
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: A
     end function CHAMELEON_zlaset_Tile_c
  end interface

  interface
     function CHAMELEON_zlaswp_Tile_c(A,K1,K2,IPIV,INCX) &
          & bind(c, name='CHAMELEON_zlaswp_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlaswp_Tile_c
       type(c_ptr), value :: A
       integer(kind=c_int), value :: K1
       integer(kind=c_int), value :: K2
       type(c_ptr), value :: IPIV
       integer(kind=c_int), value :: INCX
     end function CHAMELEON_zlaswp_Tile_c
  end interface

  interface
     function CHAMELEON_zlaswpc_Tile_c(A,K1,K2,IPIV,INCX) &
          & bind(c, name='CHAMELEON_zlaswpc_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlaswpc_Tile_c
       type(c_ptr), value :: A
       integer(kind=c_int), value :: K1
       integer(kind=c_int), value :: K2
       type(c_ptr), value :: IPIV
       integer(kind=c_int), value :: INCX
     end function CHAMELEON_zlaswpc_Tile_c
  end interface

  interface
     function CHAMELEON_zlauum_Tile_c(uplo,A) &
          & bind(c, name='CHAMELEON_zlauum_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlauum_Tile_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
     end function CHAMELEON_zlauum_Tile_c
  end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
  interface
     function CHAMELEON_zplghe_Tile_c(bump,uplo,A,seed) &
          & bind(c, name='CHAMELEON_zplghe_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zplghe_Tile_c
       real(kind=c_double), value :: bump
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       integer(kind=c_long_long), value :: seed
     end function CHAMELEON_zplghe_Tile_c
  end interface
#endif

  interface
     function CHAMELEON_zplgsy_Tile_c(bump,uplo,A,seed) &
          & bind(c, name='CHAMELEON_zplgsy_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zplgsy_Tile_c
       complex(kind=c_double_complex), value :: bump
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       integer(kind=c_long_long), value :: seed
     end function CHAMELEON_zplgsy_Tile_c
  end interface

  interface
     function CHAMELEON_zplrnt_Tile_c(A,seed) &
          & bind(c, name='CHAMELEON_zplrnt_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zplrnt_Tile_c
       type(c_ptr), value :: A
       integer(kind=c_long_long), value :: seed
     end function CHAMELEON_zplrnt_Tile_c
  end interface

  interface
     function CHAMELEON_zplrnk_Tile_c(K,C,seedA,seedB) &
          & bind(c, name='CHAMELEON_zplrnk_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zplrnk_Tile_c
       integer(kind=c_int), value :: K
       type(c_ptr), value :: C
       integer(kind=c_long_long), value :: seedA
       integer(kind=c_long_long), value :: seedB
     end function CHAMELEON_zplrnk_Tile_c
  end interface

  interface
     function CHAMELEON_zposv_Tile_c(uplo,A,B) &
          & bind(c, name='CHAMELEON_zposv_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zposv_Tile_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: B
     end function CHAMELEON_zposv_Tile_c
  end interface

  interface
     function CHAMELEON_zsysv_Tile_c(uplo,A,B) &
          & bind(c, name='CHAMELEON_zsysv_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsysv_Tile_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: B
     end function CHAMELEON_zsysv_Tile_c
  end interface

  interface
     function CHAMELEON_zpotrf_Tile_c(uplo,A) &
          & bind(c, name='CHAMELEON_zpotrf_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zpotrf_Tile_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
     end function CHAMELEON_zpotrf_Tile_c
  end interface

  interface
     function CHAMELEON_zsytrf_Tile_c(uplo,A) &
          & bind(c, name='CHAMELEON_zsytrf_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsytrf_Tile_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
     end function CHAMELEON_zsytrf_Tile_c
  end interface

  interface
     function CHAMELEON_zpotri_Tile_c(uplo,A) &
          & bind(c, name='CHAMELEON_zpotri_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zpotri_Tile_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
     end function CHAMELEON_zpotri_Tile_c
  end interface

  interface
     function CHAMELEON_zpotrs_Tile_c(uplo,A,B) &
          & bind(c, name='CHAMELEON_zpotrs_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zpotrs_Tile_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: B
     end function CHAMELEON_zpotrs_Tile_c
  end interface

  interface
     function CHAMELEON_zsytrs_Tile_c(uplo,A,B) &
          & bind(c, name='CHAMELEON_zsytrs_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsytrs_Tile_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: B
     end function CHAMELEON_zsytrs_Tile_c
  end interface

  interface
     function CHAMELEON_zsymm_Tile_c(side,uplo,alpha,A,B,beta,C) &
          & bind(c, name='CHAMELEON_zsymm_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsymm_Tile_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: uplo
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: C
     end function CHAMELEON_zsymm_Tile_c
  end interface

  interface
     function CHAMELEON_zsyr2k_Tile_c(uplo,trans,alpha,A,B,beta,C) &
          & bind(c, name='CHAMELEON_zsyr2k_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsyr2k_Tile_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: trans
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: C
     end function CHAMELEON_zsyr2k_Tile_c
  end interface

  interface
     function CHAMELEON_zsyrk_Tile_c(uplo,trans,alpha,A,beta,C) &
          & bind(c, name='CHAMELEON_zsyrk_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsyrk_Tile_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: trans
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: C
     end function CHAMELEON_zsyrk_Tile_c
  end interface

  interface
     function CHAMELEON_ztrmm_Tile_c(side,uplo,transA,diag,alpha,A,B) &
          & bind(c, name='CHAMELEON_ztrmm_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_ztrmm_Tile_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: transA
       integer(kind=c_int), value :: diag
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       type(c_ptr), value :: B
     end function CHAMELEON_ztrmm_Tile_c
  end interface

  interface
     function CHAMELEON_ztrsm_Tile_c(side,uplo,transA,diag,alpha,A,B) &
          & bind(c, name='CHAMELEON_ztrsm_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_ztrsm_Tile_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: transA
       integer(kind=c_int), value :: diag
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       type(c_ptr), value :: B
     end function CHAMELEON_ztrsm_Tile_c
  end interface

  interface
     function CHAMELEON_ztrsmpl_Tile_c(A,L,IPIV,B) &
          & bind(c, name='CHAMELEON_ztrsmpl_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_ztrsmpl_Tile_c
       type(c_ptr), value :: A
       type(c_ptr), value :: L
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: B
     end function CHAMELEON_ztrsmpl_Tile_c
  end interface

  interface
     function CHAMELEON_ztrsmrv_Tile_c(side,uplo,transA,diag,alpha,A,B) &
          & bind(c, name='CHAMELEON_ztrsmrv_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_ztrsmrv_Tile_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: transA
       integer(kind=c_int), value :: diag
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       type(c_ptr), value :: B
     end function CHAMELEON_ztrsmrv_Tile_c
  end interface

  interface
     function CHAMELEON_ztrtri_Tile_c(uplo,diag,A) &
          & bind(c, name='CHAMELEON_ztrtri_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_ztrtri_Tile_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: diag
       type(c_ptr), value :: A
     end function CHAMELEON_ztrtri_Tile_c
  end interface

  interface
     function CHAMELEON_zunglq_Tile_c(A,T,B) &
          & bind(c, name='CHAMELEON_zunglq_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zunglq_Tile_c
       type(c_ptr), value :: A
       type(c_ptr), value :: T
       type(c_ptr), value :: B
     end function CHAMELEON_zunglq_Tile_c
  end interface

  interface
     function CHAMELEON_zungqr_Tile_c(A,T,B) &
          & bind(c, name='CHAMELEON_zungqr_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zungqr_Tile_c
       type(c_ptr), value :: A
       type(c_ptr), value :: T
       type(c_ptr), value :: B
     end function CHAMELEON_zungqr_Tile_c
  end interface

  interface
     function CHAMELEON_zunmlq_Tile_c(side,trans,A,T,B) &
          & bind(c, name='CHAMELEON_zunmlq_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zunmlq_Tile_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: trans
       type(c_ptr), value :: A
       type(c_ptr), value :: T
       type(c_ptr), value :: B
     end function CHAMELEON_zunmlq_Tile_c
  end interface

  interface
     function CHAMELEON_zunmqr_Tile_c(side,trans,A,T,B) &
          & bind(c, name='CHAMELEON_zunmqr_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zunmqr_Tile_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: trans
       type(c_ptr), value :: A
       type(c_ptr), value :: T
       type(c_ptr), value :: B
     end function CHAMELEON_zunmqr_Tile_c
  end interface

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !  FORTRAN API - math functions (asynchronous interface)
  !
  interface
     function CHAMELEON_zgebrd_Tile_Async_c(A,D,E,T,sequence,request) &
          & bind(c, name='CHAMELEON_zgebrd_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgebrd_Tile_Async_c
       type(c_ptr), value :: A
       type(c_ptr), value :: D
       type(c_ptr), value :: E
       type(c_ptr), value :: T
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgebrd_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zgelqf_Tile_Async_c(A,T,sequence,request) &
          & bind(c, name='CHAMELEON_zgelqf_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgelqf_Tile_Async_c
       type(c_ptr), value :: A
       type(c_ptr), value :: T
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgelqf_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zgelqs_Tile_Async_c(A,T,B,sequence,request) &
          & bind(c, name='CHAMELEON_zgelqs_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgelqs_Tile_Async_c
       type(c_ptr), value :: A
       type(c_ptr), value :: T
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgelqs_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zgels_Tile_Async_c(trans,A,T,B,sequence,request) &
          & bind(c, name='CHAMELEON_zgels_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgels_Tile_Async_c
       integer(kind=c_int), value :: trans
       type(c_ptr), value :: A
       type(c_ptr), value :: T
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgels_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zgemm_Tile_Async_c(transA,transB,alpha,A,B,beta,C,WS,sequence,request) &
          & bind(c, name='CHAMELEON_zgemm_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgemm_Tile_Async_c
       integer(kind=c_int), value :: transA
       integer(kind=c_int), value :: transB
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: C
       type(c_ptr), value :: WS
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgemm_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zgeqrf_Tile_Async_c(A,T,sequence,request) &
          & bind(c, name='CHAMELEON_zgeqrf_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgeqrf_Tile_Async_c
       type(c_ptr), value :: A
       type(c_ptr), value :: T
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgeqrf_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zgeqrs_Tile_Async_c(A,T,B,sequence,request) &
          & bind(c, name='CHAMELEON_zgeqrs_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgeqrs_Tile_Async_c
       type(c_ptr), value :: A
       type(c_ptr), value :: T
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgeqrs_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zgesv_Tile_Async_c(A,IPIV,B,sequence,request) &
          & bind(c, name='CHAMELEON_zgesv_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgesv_Tile_Async_c
       type(c_ptr), value :: A
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgesv_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zgesv_incpiv_Tile_Async_c(A,L,IPIV,B,sequence,request) &
          & bind(c, name='CHAMELEON_zgesv_incpiv_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgesv_incpiv_Tile_Async_c
       type(c_ptr), value :: A
       type(c_ptr), value :: L
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgesv_incpiv_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zgesvd_Tile_Async_c(jobu,jobvt,A,S,U,VT,T,sequence,request) &
          & bind(c, name='CHAMELEON_zgesvd_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgesvd_Tile_Async_c
       integer(kind=c_int), value :: jobu
       integer(kind=c_int), value :: jobvt
       type(c_ptr), value :: A
       type(c_ptr), value :: S
       type(c_ptr), value :: U
       type(c_ptr), value :: VT
       type(c_ptr), value :: T
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgesvd_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zgetrf_Tile_Async_c(A,IPIV,sequence,request) &
          & bind(c, name='CHAMELEON_zgetrf_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetrf_Tile_Async_c
       type(c_ptr), value :: A
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgetrf_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zgetrf_incpiv_Tile_Async_c(A,L,IPIV,sequence,request) &
          & bind(c, name='CHAMELEON_zgetrf_incpiv_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetrf_incpiv_Tile_Async_c
       type(c_ptr), value :: A
       type(c_ptr), value :: L
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgetrf_incpiv_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zgetrf_nopiv_Tile_Async_c(A,sequence,request) &
          & bind(c, name='CHAMELEON_zgetrf_nopiv_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetrf_nopiv_Tile_Async_c
       type(c_ptr), value :: A
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgetrf_nopiv_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zgetri_Tile_Async_c(A,IPIV,W,sequence,request) &
          & bind(c, name='CHAMELEON_zgetri_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetri_Tile_Async_c
       type(c_ptr), value :: A
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: W
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgetri_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zgetrs_Tile_Async_c(trans,A,IPIV,B,sequence,request) &
          & bind(c, name='CHAMELEON_zgetrs_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetrs_Tile_Async_c
       integer(kind=c_int), value :: trans
       type(c_ptr), value :: A
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgetrs_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zgetrs_incpiv_Tile_Async_c(A,L,IPIV,B,sequence,request) &
          & bind(c, name='CHAMELEON_zgetrs_incpiv_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zgetrs_incpiv_Tile_Async_c
       type(c_ptr), value :: A
       type(c_ptr), value :: L
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zgetrs_incpiv_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zheev_Tile_Async_c(jobz,uplo,A,W,T,Q,LDQ,sequence,request) &
          & bind(c, name='CHAMELEON_zheev_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zheev_Tile_Async_c
       integer(kind=c_int), value :: jobz
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: W
       type(c_ptr), value :: T
       type(c_ptr), value :: Q
       integer(kind=c_int), value :: LDQ
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zheev_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zheevd_Tile_Async_c(jobz,uplo,A,W,T,Q,LDQ,sequence,request) &
          & bind(c, name='CHAMELEON_zheevd_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zheevd_Tile_Async_c
       integer(kind=c_int), value :: jobz
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: W
       type(c_ptr), value :: T
       type(c_ptr), value :: Q
       integer(kind=c_int), value :: LDQ
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zheevd_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zhegst_Tile_Async_c(itype,uplo,A,B,sequence,request) &
          & bind(c, name='CHAMELEON_zhegst_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zhegst_Tile_Async_c
       integer(kind=c_int), value :: itype
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zhegst_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zhegv_Tile_Async_c(itype,jobz,uplo,A,B,W,T,Q,sequence,request) &
          & bind(c, name='CHAMELEON_zhegv_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zhegv_Tile_Async_c
       integer(kind=c_int), value :: itype
       integer(kind=c_int), value :: jobz
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       type(c_ptr), value :: W
       type(c_ptr), value :: T
       type(c_ptr), value :: Q
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zhegv_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zhegvd_Tile_Async_c(itype,jobz,uplo,A,B,W,T,Q,sequence,request) &
          & bind(c, name='CHAMELEON_zhegvd_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zhegvd_Tile_Async_c
       integer(kind=c_int), value :: itype
       integer(kind=c_int), value :: jobz
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       type(c_ptr), value :: W
       type(c_ptr), value :: T
       type(c_ptr), value :: Q
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zhegvd_Tile_Async_c
  end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
  interface
     function CHAMELEON_zhemm_Tile_Async_c(side,uplo,alpha,A,B,beta,C,sequence,request) &
          & bind(c, name='CHAMELEON_zhemm_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zhemm_Tile_Async_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: uplo
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: C
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zhemm_Tile_Async_c
  end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
  interface
     function CHAMELEON_zher2k_Tile_Async_c(uplo,trans,alpha,A,B,beta,C,sequence,request) &
          & bind(c, name='CHAMELEON_zher2k_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zher2k_Tile_Async_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: trans
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       real(kind=c_double), value :: beta
       type(c_ptr), value :: C
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zher2k_Tile_Async_c
  end interface
#endif

#if defined(PRECISION_z) || defined(PRECISION_c)
  interface
     function CHAMELEON_zherk_Tile_Async_c(uplo,trans,alpha,A,beta,C,sequence,request) &
          & bind(c, name='CHAMELEON_zherk_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zherk_Tile_Async_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: trans
       real(kind=c_double), value :: alpha
       type(c_ptr), value :: A
       real(kind=c_double), value :: beta
       type(c_ptr), value :: C
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zherk_Tile_Async_c
  end interface
#endif

  interface
     function CHAMELEON_zhetrd_Tile_Async_c(jobz,uplo,A,D,E,T,Q,LDQ,sequence,request) &
          & bind(c, name='CHAMELEON_zhetrd_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zhetrd_Tile_Async_c
       integer(kind=c_int), value :: jobz
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: D
       type(c_ptr), value :: E
       type(c_ptr), value :: T
       type(c_ptr), value :: Q
       integer(kind=c_int), value :: LDQ
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zhetrd_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zlacpy_Tile_Async_c(uplo,A,B,sequence,request) &
          & bind(c, name='CHAMELEON_zlacpy_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlacpy_Tile_Async_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zlacpy_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zlange_Tile_Async_c(norm,A,work,value,sequence,request) &
          & bind(c, name='CHAMELEON_zlange_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlange_Tile_Async_c
       integer(kind=c_int), value :: norm
       type(c_ptr), value :: A
       type(c_ptr), value :: work
       type(c_ptr), value :: value
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zlange_Tile_Async_c
  end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
  interface
     function CHAMELEON_zlanhe_Tile_Async_c(norm,uplo,A,work,value,sequence,request) &
          & bind(c, name='CHAMELEON_zlanhe_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlanhe_Tile_Async_c
       integer(kind=c_int), value :: norm
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: work
       type(c_ptr), value :: value
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zlanhe_Tile_Async_c
  end interface
#endif

  interface
     function CHAMELEON_zlansy_Tile_Async_c(norm,uplo,A,work,value,sequence,request) &
          & bind(c, name='CHAMELEON_zlansy_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlansy_Tile_Async_c
       integer(kind=c_int), value :: norm
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: work
       type(c_ptr), value :: value
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zlansy_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zlaset_Tile_Async_c(uplo,alpha,beta,A,sequence,request) &
          & bind(c, name='CHAMELEON_zlaset_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlaset_Tile_Async_c
       integer(kind=c_int), value :: uplo
       complex(kind=c_double_complex), value :: alpha
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: A
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zlaset_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zlaswp_Tile_Async_c(A,K1,K2,IPIV,INCX,sequence,request) &
          & bind(c, name='CHAMELEON_zlaswp_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlaswp_Tile_Async_c
       type(c_ptr), value :: A
       integer(kind=c_int), value :: K1
       integer(kind=c_int), value :: K2
       type(c_ptr), value :: IPIV
       integer(kind=c_int), value :: INCX
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zlaswp_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zlaswpc_Tile_Async_c(A,K1,K2,IPIV,INCX,sequence,request) &
          & bind(c, name='CHAMELEON_zlaswpc_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlaswpc_Tile_Async_c
       type(c_ptr), value :: A
       integer(kind=c_int), value :: K1
       integer(kind=c_int), value :: K2
       type(c_ptr), value :: IPIV
       integer(kind=c_int), value :: INCX
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zlaswpc_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zlauum_Tile_Async_c(uplo,A,sequence,request) &
          & bind(c, name='CHAMELEON_zlauum_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zlauum_Tile_Async_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zlauum_Tile_Async_c
  end interface

#if defined(PRECISION_z) || defined(PRECISION_c)
  interface
     function CHAMELEON_zplghe_Tile_Async_c(bump,uplo,A,seed,sequence,request) &
          & bind(c, name='CHAMELEON_zplghe_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zplghe_Tile_Async_c
       real(kind=c_double), value :: bump
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       integer(kind=c_long_long), value :: seed
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zplghe_Tile_Async_c
  end interface
#endif

  interface
     function CHAMELEON_zplgsy_Tile_Async_c(bump,uplo,A,seed,sequence,request) &
          & bind(c, name='CHAMELEON_zplgsy_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zplgsy_Tile_Async_c
       complex(kind=c_double_complex), value :: bump
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       integer(kind=c_long_long), value :: seed
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zplgsy_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zplrnt_Tile_Async_c(A,seed,sequence,request) &
          & bind(c, name='CHAMELEON_zplrnt_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zplrnt_Tile_Async_c
       type(c_ptr), value :: A
       integer(kind=c_long_long), value :: seed
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zplrnt_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zplrnk_Tile_Async_c(K,C,seedA,seedB,sequence,request) &
          & bind(c, name='CHAMELEON_zplrnk_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zplrnk_Tile_Async_c
       integer(kind=c_int), value :: K
       type(c_ptr), value :: C
       integer(kind=c_long_long), value :: seedA
       integer(kind=c_long_long), value :: seedB
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zplrnk_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zposv_Tile_Async_c(uplo,A,B,sequence,request) &
          & bind(c, name='CHAMELEON_zposv_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zposv_Tile_Async_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zposv_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zsysv_Tile_Async_c(uplo,A,B,sequence,request) &
          & bind(c, name='CHAMELEON_zsysv_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsysv_Tile_Async_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zsysv_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zpotrf_Tile_Async_c(uplo,A,sequence,request) &
          & bind(c, name='CHAMELEON_zpotrf_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zpotrf_Tile_Async_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zpotrf_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zsytrf_Tile_Async_c(uplo,A,sequence,request) &
          & bind(c, name='CHAMELEON_zsytrf_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsytrf_Tile_Async_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zsytrf_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zpotri_Tile_Async_c(uplo,A,sequence,request) &
          & bind(c, name='CHAMELEON_zpotri_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zpotri_Tile_Async_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zpotri_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zpotrs_Tile_Async_c(uplo,A,B,sequence,request) &
          & bind(c, name='CHAMELEON_zpotrs_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zpotrs_Tile_Async_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zpotrs_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zsytrs_Tile_Async_c(uplo,A,B,sequence,request) &
          & bind(c, name='CHAMELEON_zsytrs_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsytrs_Tile_Async_c
       integer(kind=c_int), value :: uplo
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zsytrs_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zsymm_Tile_Async_c(side,uplo,alpha,A,B,beta,C,sequence,request) &
          & bind(c, name='CHAMELEON_zsymm_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsymm_Tile_Async_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: uplo
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: C
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zsymm_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zsyr2k_Tile_Async_c(uplo,trans,alpha,A,B,beta,C,sequence,request) &
          & bind(c, name='CHAMELEON_zsyr2k_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsyr2k_Tile_Async_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: trans
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: C
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zsyr2k_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zsyrk_Tile_Async_c(uplo,trans,alpha,A,beta,C,sequence,request) &
          & bind(c, name='CHAMELEON_zsyrk_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zsyrk_Tile_Async_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: trans
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       complex(kind=c_double_complex), value :: beta
       type(c_ptr), value :: C
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zsyrk_Tile_Async_c
  end interface

  interface
     function CHAMELEON_ztrmm_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request) &
          & bind(c, name='CHAMELEON_ztrmm_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_ztrmm_Tile_Async_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: transA
       integer(kind=c_int), value :: diag
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_ztrmm_Tile_Async_c
  end interface

  interface
     function CHAMELEON_ztrsm_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request) &
          & bind(c, name='CHAMELEON_ztrsm_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_ztrsm_Tile_Async_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: transA
       integer(kind=c_int), value :: diag
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_ztrsm_Tile_Async_c
  end interface

  interface
     function CHAMELEON_ztrsmpl_Tile_Async_c(A,L,IPIV,B,sequence,request) &
          & bind(c, name='CHAMELEON_ztrsmpl_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_ztrsmpl_Tile_Async_c
       type(c_ptr), value :: A
       type(c_ptr), value :: L
       type(c_ptr), value :: IPIV
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_ztrsmpl_Tile_Async_c
  end interface

  interface
     function CHAMELEON_ztrsmrv_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request) &
          & bind(c, name='CHAMELEON_ztrsmrv_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_ztrsmrv_Tile_Async_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: transA
       integer(kind=c_int), value :: diag
       complex(kind=c_double_complex), value :: alpha
       type(c_ptr), value :: A
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_ztrsmrv_Tile_Async_c
  end interface

  interface
     function CHAMELEON_ztrtri_Tile_Async_c(uplo,diag,A,sequence,request) &
          & bind(c, name='CHAMELEON_ztrtri_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_ztrtri_Tile_Async_c
       integer(kind=c_int), value :: uplo
       integer(kind=c_int), value :: diag
       type(c_ptr), value :: A
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_ztrtri_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zunglq_Tile_Async_c(A,T,B,sequence,request) &
          & bind(c, name='CHAMELEON_zunglq_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zunglq_Tile_Async_c
       type(c_ptr), value :: A
       type(c_ptr), value :: T
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zunglq_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zungqr_Tile_Async_c(A,T,B,sequence,request) &
          & bind(c, name='CHAMELEON_zungqr_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zungqr_Tile_Async_c
       type(c_ptr), value :: A
       type(c_ptr), value :: T
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zungqr_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zunmlq_Tile_Async_c(side,trans,A,T,B,sequence,request) &
          & bind(c, name='CHAMELEON_zunmlq_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zunmlq_Tile_Async_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: trans
       type(c_ptr), value :: A
       type(c_ptr), value :: T
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zunmlq_Tile_Async_c
  end interface

  interface
     function CHAMELEON_zunmqr_Tile_Async_c(side,trans,A,T,B,sequence,request) &
          & bind(c, name='CHAMELEON_zunmqr_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_zunmqr_Tile_Async_c
       integer(kind=c_int), value :: side
       integer(kind=c_int), value :: trans
       type(c_ptr), value :: A
       type(c_ptr), value :: T
       type(c_ptr), value :: B
       type(c_ptr), value :: sequence
       type(c_ptr), value :: request
     end function CHAMELEON_zunmqr_Tile_Async_c
  end interface

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !  FORTRAN API - workspace allocation
  !
  interface
     function CHAMELEON_Alloc_Workspace_zgebrd_c(M,N,descT,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zgebrd')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zgebrd_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr) :: descT ! descT is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zgebrd_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zgeev_c(N,descT,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zgeev')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zgeev_c
       integer(kind=c_int), value :: N
       type(c_ptr) :: descT ! descT is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zgeev_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zgehrd_c(N,descT,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zgehrd')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zgehrd_c
       integer(kind=c_int), value :: N
       type(c_ptr) :: descT ! descT is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zgehrd_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zgelqf_c(M,N,T,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zgelqf')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zgelqf_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr) :: T ! T is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zgelqf_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zgelqf_Tile_c(M,N,descT,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zgelqf_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zgelqf_Tile_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr) :: descT ! descT is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zgelqf_Tile_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zgels_c(M,N,T,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zgels')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zgels_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr) :: T ! T is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zgels_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zgels_Tile_c(M,N,descT,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zgels_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zgels_Tile_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr) :: descT ! descT is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zgels_Tile_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zgeqrf_c(M,N,T,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zgeqrf')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zgeqrf_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr) :: T ! T is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zgeqrf_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zgeqrf_Tile_c(M,N,descT,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zgeqrf_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zgeqrf_Tile_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr) :: descT ! descT is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zgeqrf_Tile_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zgesv_incpiv_c(N,descL,IPIV,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zgesv_incpiv')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zgesv_incpiv_c
       integer(kind=c_int), value :: N
       type(c_ptr) :: descL ! descL is **, so pass by reference
       type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zgesv_incpiv_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zgesv_incpiv_Tile_c(N,descL,IPIV,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zgesv_incpiv_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zgesv_incpiv_Tile_c
       integer(kind=c_int), value :: N
       type(c_ptr) :: descL ! descL is **, so pass by reference
       type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zgesv_incpiv_Tile_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zgesvd_c(M,N,descT,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zgesvd')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zgesvd_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr) :: descT ! descT is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zgesvd_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zgetrf_incpiv_c(M,N,descL,IPIV,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zgetrf_incpiv')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zgetrf_incpiv_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr) :: descL ! descL is **, so pass by reference
       type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zgetrf_incpiv_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zgetrf_incpiv_Tile_c(N,descL,IPIV,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zgetrf_incpiv_Tile')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zgetrf_incpiv_Tile_c
       integer(kind=c_int), value :: N
       type(c_ptr) :: descL ! descL is **, so pass by reference
       type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zgetrf_incpiv_Tile_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zgetri_Tile_Async_c(A,W) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zgetri_Tile_Async')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zgetri_Tile_Async_c
       type(c_ptr), value :: A
       type(c_ptr), value :: W
     end function CHAMELEON_Alloc_Workspace_zgetri_Tile_Async_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zheev_c(M,N,descT,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zheev')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zheev_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr) :: descT ! descT is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zheev_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zheevd_c(M,N,descT,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zheevd')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zheevd_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr) :: descT ! descT is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zheevd_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zhegv_c(M,N,descT,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zhegv')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zhegv_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr) :: descT ! descT is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zhegv_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zhegvd_c(M,N,descT,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zhegvd')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zhegvd_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr) :: descT ! descT is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zhegvd_c
  end interface

  interface
     function CHAMELEON_Alloc_Workspace_zhetrd_c(M,N,descT,p,q) &
          & bind(c, name='CHAMELEON_Alloc_Workspace_zhetrd')
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: CHAMELEON_Alloc_Workspace_zhetrd_c
       integer(kind=c_int), value :: M
       integer(kind=c_int), value :: N
       type(c_ptr) :: descT ! descT is **, so pass by reference
       integer(kind=c_int), value :: p
       integer(kind=c_int), value :: q
     end function CHAMELEON_Alloc_Workspace_zhetrd_c
  end interface

contains

  subroutine CHAMELEON_zgebrd(M,N,A,LDA,D,E,T,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    real(kind=c_double), intent(out), target :: D(*)
    real(kind=c_double), intent(out), target :: E(*)
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgebrd_c(M,N,c_loc(A),LDA,c_loc(D),c_loc(E),T)
  end subroutine CHAMELEON_zgebrd

  subroutine CHAMELEON_zgelqf(M,N,A,LDA,T,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgelqf_c(M,N,c_loc(A),LDA,T)
  end subroutine CHAMELEON_zgelqf

  subroutine CHAMELEON_zgelqs(M,N,NRHS,A,LDA,T,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: NRHS
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgelqs_c(M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
  end subroutine CHAMELEON_zgelqs

  subroutine CHAMELEON_zgels(trans,M,N,NRHS,A,LDA,T,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: NRHS
    integer(kind=c_int), intent(in) :: trans
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgels_c(trans,M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
  end subroutine CHAMELEON_zgels

  subroutine CHAMELEON_zgemm(transA,transB,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: K
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: LDC
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: transA
    integer(kind=c_int), intent(in) :: transB
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
    complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
    info = CHAMELEON_zgemm_c(transA,transB,M,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
  end subroutine CHAMELEON_zgemm

  subroutine CHAMELEON_zgeqrf(M,N,A,LDA,T,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgeqrf_c(M,N,c_loc(A),LDA,T)
  end subroutine CHAMELEON_zgeqrf

  subroutine CHAMELEON_zgeqrs(M,N,NRHS,A,LDA,T,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: NRHS
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgeqrs_c(M,N,NRHS,c_loc(A),LDA,T,c_loc(B),LDB)
  end subroutine CHAMELEON_zgeqrs

  subroutine CHAMELEON_zgesv(N,NRHS,A,LDA,IPIV,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: NRHS
    integer(kind=c_int), intent(out), target :: IPIV(*)
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    info = CHAMELEON_zgesv_c(N,NRHS,c_loc(A),LDA,c_loc(IPIV),c_loc(B),LDB)
  end subroutine CHAMELEON_zgesv

  subroutine CHAMELEON_zgesv_incpiv(N,NRHS,A,LDA,L,IPIV,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: NRHS
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    type(c_ptr), value :: IPIV ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: L ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgesv_incpiv_c(N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
  end subroutine CHAMELEON_zgesv_incpiv

  subroutine CHAMELEON_zgesvd(jobu,jobvt,M,N,A,LDA,S,U,LDU,VT,LDVT,T,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDU
    integer(kind=c_int), intent(in) :: LDVT
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: jobu
    integer(kind=c_int), intent(in) :: jobvt
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(out), target :: U(LDU,*)
    complex(kind=c_double_complex), intent(out), target :: VT(LDVT,*)
    real(kind=c_double), intent(out), target :: S(*)
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgesvd_c(jobu,jobvt,M,N,c_loc(A),LDA,c_loc(S),c_loc(U),LDU,c_loc(VT),LDVT,T)
  end subroutine CHAMELEON_zgesvd

  subroutine CHAMELEON_zgetrf(M,N,A,LDA,IPIV,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(out), target :: IPIV(*)
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    info = CHAMELEON_zgetrf_c(M,N,c_loc(A),LDA,c_loc(IPIV))
  end subroutine CHAMELEON_zgetrf

  subroutine CHAMELEON_zgetrf_incpiv(M,N,A,LDA,L,IPIV,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    type(c_ptr), value :: IPIV ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: L    ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgetrf_incpiv_c(M,N,c_loc(A),LDA,L,IPIV)
  end subroutine CHAMELEON_zgetrf_incpiv

  subroutine CHAMELEON_zgetrf_nopiv(M,N,A,LDA,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    info = CHAMELEON_zgetrf_nopiv_c(M,N,c_loc(A),LDA)
  end subroutine CHAMELEON_zgetrf_nopiv

  subroutine CHAMELEON_zgetrs(trans,N,NRHS,A,LDA,IPIV,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: NRHS
    integer(kind=c_int), intent(in) :: trans
    integer(kind=c_int), intent(in), target :: IPIV(*)
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    info = CHAMELEON_zgetrs_c(trans,N,NRHS,c_loc(A),LDA,c_loc(IPIV),c_loc(B),LDB)
  end subroutine CHAMELEON_zgetrs

  subroutine CHAMELEON_zgetrs_incpiv(trans,N,NRHS,A,LDA,L,IPIV,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: NRHS
    integer(kind=c_int), intent(in) :: trans
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    type(c_ptr), value :: IPIV ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: L    ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgetrs_incpiv_c(trans,N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
  end subroutine CHAMELEON_zgetrs_incpiv

#if defined(PRECISION_z) || defined(PRECISION_c)
  subroutine CHAMELEON_zhemm(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: LDC
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
    complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
    info = CHAMELEON_zhemm_c(side,uplo,M,N,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
  end subroutine CHAMELEON_zhemm

  subroutine CHAMELEON_zherk(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: K
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDC
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: trans
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
    real(kind=c_double), intent(in) :: alpha
    real(kind=c_double), intent(in) :: beta
    info = CHAMELEON_zherk_c(uplo,trans,N,K,alpha,c_loc(A),LDA,beta,c_loc(C),LDC)
  end subroutine CHAMELEON_zherk

  subroutine CHAMELEON_zher2k(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: K
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: LDC
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: trans
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
    complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
    real(kind=c_double), intent(in) :: beta
    info = CHAMELEON_zher2k_c(uplo,trans,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
  end subroutine CHAMELEON_zher2k
#endif

  subroutine CHAMELEON_zheev(jobz,uplo,N,A,LDA,W,T,Q,LDQ,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDQ
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: jobz
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    real(kind=c_double), intent(out), target :: W(*)
    complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zheev_c(jobz,uplo,N,c_loc(A),LDA,c_loc(W),T,c_loc(Q),LDQ)
  end subroutine CHAMELEON_zheev

  subroutine CHAMELEON_zheevd(jobz,uplo,N,A,LDA,W,T,Q,LDQ,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDQ
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: jobz
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    real(kind=c_double), intent(out), target :: W(*)
    complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zheevd_c(jobz,uplo,N,c_loc(A),LDA,c_loc(W),T,c_loc(Q),LDQ)
  end subroutine CHAMELEON_zheevd

  subroutine CHAMELEON_zhegv(itype,jobz,uplo,N,A,LDA,B,LDB,W,T,Q,LDQ,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: LDQ
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: itype
    integer(kind=c_int), intent(in) :: jobz
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    real(kind=c_double), intent(out), target :: W(*)
    complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zhegv_c(itype,jobz,uplo,N,c_loc(A),LDA,c_loc(B),LDB,c_loc(W),T,c_loc(Q),LDQ)
  end subroutine CHAMELEON_zhegv

  subroutine CHAMELEON_zhegvd(itype,jobz,uplo,N,A,LDA,B,LDB,W,T,Q,LDQ,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: LDQ
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: itype
    integer(kind=c_int), intent(in) :: jobz
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    real(kind=c_double), intent(out), target :: W(*)
    complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zhegvd_c(itype,jobz,uplo,N,c_loc(A),LDA,c_loc(B),LDB,c_loc(W),T,c_loc(Q),LDQ)
  end subroutine CHAMELEON_zhegvd

  subroutine CHAMELEON_zhegst(itype,uplo,N,A,LDA,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: itype
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    info = CHAMELEON_zhegst_c(itype,uplo,N,c_loc(A),LDA,c_loc(B),LDB)
  end subroutine CHAMELEON_zhegst

  subroutine CHAMELEON_zhetrd(jobz,uplo,N,A,LDA,D,E,descT,Q,LDQ,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: jobz
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDQ
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    real(kind=c_double), intent(out), target :: D(*)
    real(kind=c_double), intent(out), target :: E(*)
    type(c_ptr), value :: descT ! Arg managed by CHAMELEON: opaque to Fortran
    complex(kind=c_double_complex), intent(inout), target :: Q(LDQ,*)
    info = CHAMELEON_zhetrd_c(jobz,uplo,N,c_loc(A),LDA,c_loc(D),c_loc(E),descT,c_loc(Q),LDQ)
  end subroutine CHAMELEON_zhetrd

  function CHAMELEON_zlange(norm,M,N,A,LDA,work)
    use iso_c_binding
    implicit none
    real(kind=c_double) :: CHAMELEON_zlange
    real(kind=c_double), intent(inout), target :: work(*)
    integer(kind=c_int), intent(in) :: norm
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    CHAMELEON_zlange = CHAMELEON_zlange_c(norm,M,N,c_loc(A),LDA,c_loc(work))
  end function CHAMELEON_zlange

#if defined(PRECISION_z) || defined(PRECISION_c)
  function CHAMELEON_zlanhe(norm,uplo,N,A,LDA,work)
    use iso_c_binding
    implicit none
    real(kind=c_double) :: CHAMELEON_zlanhe
    real(kind=c_double), intent(inout), target :: work(*)
    integer(kind=c_int), intent(in) :: norm
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: N
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    CHAMELEON_zlanhe = CHAMELEON_zlanhe_c(norm,uplo,N,c_loc(A),LDA,c_loc(work))
  end function CHAMELEON_zlanhe
#endif

  function CHAMELEON_zlansy(norm,uplo,N,A,LDA,work)
    use iso_c_binding
    implicit none
    real(kind=c_double) :: CHAMELEON_zlansy
    real(kind=c_double), intent(inout), target :: work(*)
    integer(kind=c_int), intent(in) :: norm
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: N
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    CHAMELEON_zlansy = CHAMELEON_zlansy_c(norm,uplo,N,c_loc(A),LDA,c_loc(work))
  end function CHAMELEON_zlansy

  subroutine CHAMELEON_zlaswp(N,A,LDA,K1,K2,IPIV,INCX,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: INCX
    integer(kind=c_int), intent(in) :: K1
    integer(kind=c_int), intent(in) :: K2
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in), target :: IPIV(*)
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    info = CHAMELEON_zlaswp_c(N,c_loc(A),LDA,K1,K2,c_loc(IPIV),INCX)
  end subroutine CHAMELEON_zlaswp

  subroutine CHAMELEON_zlauum(uplo,N,A,LDA,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    info = CHAMELEON_zlauum_c(uplo,N,c_loc(A),LDA)
  end subroutine CHAMELEON_zlauum

#if defined(PRECISION_z) || defined(PRECISION_c)
  subroutine CHAMELEON_zplghe(bump,uplo,N,A,LDA,seed,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    real(kind=c_double), intent(in) :: bump
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_long_long), intent(in) :: seed
    complex(kind=c_double_complex), intent(out), target :: A(LDA,*)
    info = CHAMELEON_zplghe_c(bump,N,c_loc(A),LDA,seed)
  end subroutine CHAMELEON_zplghe
#endif

  subroutine CHAMELEON_zplgsy(bump,uplo,N,A,LDA,seed,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    complex(kind=c_double_complex), intent(in) :: bump
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_long_long), intent(in) :: seed
    complex(kind=c_double_complex), intent(out), target :: A(LDA,*)
    info = CHAMELEON_zplgsy_c(bump,N,c_loc(A),LDA,seed)
  end subroutine CHAMELEON_zplgsy

  subroutine CHAMELEON_zplrnt(M,N,A,LDA,seed,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_long_long), intent(in) :: seed
    complex(kind=c_double_complex), intent(out), target :: A(LDA,*)
    info = CHAMELEON_zplrnt_c(M,N,c_loc(A),LDA,seed)
  end subroutine CHAMELEON_zplrnt

  subroutine CHAMELEON_zplrnk(M,N,K,C,LDC,seedA,seedB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: K
    integer(kind=c_int), intent(in) :: LDC
    integer(kind=c_long_long), intent(in) :: seedA
    integer(kind=c_long_long), intent(in) :: seedB
    complex(kind=c_double_complex), intent(out), target :: C(LDC,*)
    info = CHAMELEON_zplrnk_c(M,N,K,c_loc(C),LDC,seedA,seedB)
  end subroutine CHAMELEON_zplrnk

  subroutine CHAMELEON_zposv(uplo,N,NRHS,A,LDA,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: NRHS
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    info = CHAMELEON_zposv_c(uplo,N,NRHS,c_loc(A),LDA,c_loc(B),LDB)
  end subroutine CHAMELEON_zposv

  subroutine CHAMELEON_zpotrf(uplo,N,A,LDA,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    info = CHAMELEON_zpotrf_c(uplo,N,c_loc(A),LDA)
  end subroutine CHAMELEON_zpotrf

  subroutine CHAMELEON_zsytrf(uplo,N,A,LDA,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    info = CHAMELEON_zsytrf_c(uplo,N,c_loc(A),LDA)
  end subroutine CHAMELEON_zsytrf

  subroutine CHAMELEON_zpotri(uplo,N,A,LDA,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    info = CHAMELEON_zpotri_c(uplo,N,c_loc(A),LDA)
  end subroutine CHAMELEON_zpotri

  subroutine CHAMELEON_zpotrs(uplo,N,NRHS,A,LDA,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: NRHS
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    info = CHAMELEON_zpotrs_c(uplo,N,NRHS,c_loc(A),LDA,c_loc(B),LDB)
  end subroutine CHAMELEON_zpotrs

  subroutine CHAMELEON_zsytrs(uplo,N,NRHS,A,LDA,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: NRHS
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    info = CHAMELEON_zsytrs_c(uplo,N,NRHS,c_loc(A),LDA,c_loc(B),LDB)
  end subroutine CHAMELEON_zsytrs

  subroutine CHAMELEON_zsymm(side,uplo,M,N,alpha,A,LDA,B,LDB,beta,C,LDC,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: LDC
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
    complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
    info = CHAMELEON_zsymm_c(side,uplo,M,N,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
  end subroutine CHAMELEON_zsymm

  subroutine CHAMELEON_zsyrk(uplo,trans,N,K,alpha,A,LDA,beta,C,LDC,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: K
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDC
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: trans
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
    info = CHAMELEON_zsyrk_c(uplo,trans,N,K,alpha,c_loc(A),LDA,beta,c_loc(C),LDC)
  end subroutine CHAMELEON_zsyrk

  subroutine CHAMELEON_zsyr2k(uplo,trans,N,K,alpha,A,LDA,B,LDB,beta,C,LDC,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: K
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: LDC
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: trans
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(in), target :: B(LDB,*)
    complex(kind=c_double_complex), intent(inout), target :: C(LDC,*)
    info = CHAMELEON_zsyr2k_c(uplo,trans,N,K,alpha,c_loc(A),LDA,c_loc(B),LDB,beta,c_loc(C),LDC)
  end subroutine CHAMELEON_zsyr2k

  subroutine CHAMELEON_ztrmm(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: NRHS
    integer(kind=c_int), intent(in) :: diag
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: transA
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    info = CHAMELEON_ztrmm_c(side,uplo,transA,diag,N,NRHS,alpha,c_loc(A),LDA,c_loc(B),LDB)
  end subroutine CHAMELEON_ztrmm

  subroutine CHAMELEON_ztrsm(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: NRHS
    integer(kind=c_int), intent(in) :: diag
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: transA
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    info = CHAMELEON_ztrsm_c(side,uplo,transA,diag,N,NRHS,alpha,c_loc(A),LDA,c_loc(B),LDB)
  end subroutine CHAMELEON_ztrsm

  subroutine CHAMELEON_ztrsmpl(N,NRHS,A,LDA,L,IPIV,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: NRHS
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    type(c_ptr), value :: L    ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: IPIV ! Arg managed by CHAMELEON: opaque to Fortran
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    info = CHAMELEON_ztrsmpl_c(N,NRHS,c_loc(A),LDA,L,IPIV,c_loc(B),LDB)
  end subroutine CHAMELEON_ztrsmpl

  subroutine CHAMELEON_ztrtri(uplo,diag,N,A,LDA,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: diag
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    info = CHAMELEON_ztrtri_c(uplo,diag,N,c_loc(A),LDA)
  end subroutine CHAMELEON_ztrtri

  subroutine CHAMELEON_zunglq(M,N,K,A,LDA,T,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    integer(kind=c_int), intent(in) :: K
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(out), target :: B(LDB,*)
    info = CHAMELEON_zunglq_c(M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
  end subroutine CHAMELEON_zunglq

  subroutine CHAMELEON_zungqr(M,N,K,A,LDA,T,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    integer(kind=c_int), intent(in) :: K
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(out), target :: B(LDB,*)
    info = CHAMELEON_zungqr_c(M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
  end subroutine CHAMELEON_zungqr

  subroutine CHAMELEON_zunmlq(side,trans,M,N,K,A,LDA,T,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    integer(kind=c_int), intent(in) :: K
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: trans
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    info = CHAMELEON_zunmlq_c(side,trans,M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
  end subroutine CHAMELEON_zunmlq

  subroutine CHAMELEON_zunmqr(side,trans,M,N,K,A,LDA,T,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    integer(kind=c_int), intent(in) :: K
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: trans
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    info = CHAMELEON_zunmqr_c(side,trans,M,N,K,c_loc(A),LDA,T,c_loc(B),LDB)
  end subroutine CHAMELEON_zunmqr

  subroutine CHAMELEON_zgecfi(m,n,A,fin,imb,inb,fout,omb,onb,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    complex(kind=c_double_complex), intent(inout), target :: A(*)
    integer(kind=c_int), intent(in) :: fin
    integer(kind=c_int), intent(in) :: fout
    integer(kind=c_int), intent(in) :: imb
    integer(kind=c_int), intent(in) :: inb
    integer(kind=c_int), intent(in) :: omb
    integer(kind=c_int), intent(in) :: onb
    integer(kind=c_int), intent(in) :: m
    integer(kind=c_int), intent(in) :: n
    info = CHAMELEON_zgecfi_c(m,n,c_loc(A),fin,imb,inb,fout,omb,onb)
  end subroutine CHAMELEON_zgecfi

  subroutine CHAMELEON_zgetmi(m,n,A,fin,mb,nb,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    complex(kind=c_double_complex), intent(inout), target :: A(*)
    integer(kind=c_int), intent(in) :: fin
    integer(kind=c_int), intent(in) :: mb
    integer(kind=c_int), intent(in) :: nb
    integer(kind=c_int), intent(in) :: m
    integer(kind=c_int), intent(in) :: n
    info = CHAMELEON_zgetmi_c(m,n,c_loc(A),fin,mb,nb)
  end subroutine CHAMELEON_zgetmi

  subroutine CHAMELEON_zgetri(N,A,LDA,IPIV,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in), target :: IPIV(*)
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    info = CHAMELEON_zgetri_c(N,c_loc(A),LDA,c_loc(IPIV))
  end subroutine CHAMELEON_zgetri

  subroutine CHAMELEON_zlacpy(uplo,M,N,A,LDA,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(out), target :: B(LDB,*)
    info = CHAMELEON_zlacpy_c(uplo,M,N,c_loc(A),LDA,c_loc(B),LDB)
  end subroutine CHAMELEON_zlacpy

  subroutine CHAMELEON_zlaset(uplo,M,N,alpha,beta,A,LDA,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    info = CHAMELEON_zlaset_c(uplo,M,N,alpha,beta,c_loc(A),LDA)
  end subroutine CHAMELEON_zlaset

  subroutine CHAMELEON_zlaswpc(N,A,LDA,K1,K2,IPIV,INCX,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in), target :: IPIV(*)
    integer(kind=c_int), intent(in) :: INCX
    integer(kind=c_int), intent(in) :: K1
    integer(kind=c_int), intent(in) :: K2
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: N
    complex(kind=c_double_complex), intent(inout), target :: A(LDA,*)
    info = CHAMELEON_zlaswpc_c(N,c_loc(A),LDA,K1,K2,c_loc(IPIV),INCX)
  end subroutine CHAMELEON_zlaswpc

  subroutine CHAMELEON_ztrsmrv(side,uplo,transA,diag,N,NRHS,alpha,A,LDA,B,LDB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: diag
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: transA
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_int), intent(in) :: LDA
    integer(kind=c_int), intent(in) :: LDB
    integer(kind=c_int), intent(in) :: N
    integer(kind=c_int), intent(in) :: NRHS
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in), target :: A(LDA,*)
    complex(kind=c_double_complex), intent(inout), target :: B(LDB,*)
    info = CHAMELEON_ztrsmrv_c(side,uplo,transA,diag,N,NRHS,alpha,c_loc(A),LDA,c_loc(B),LDB)
  end subroutine CHAMELEON_ztrsmrv

  subroutine CHAMELEON_zgebrd_Tile(A,D,E,T,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    real(kind=c_double), intent(out), target :: D(*)
    real(kind=c_double), intent(out), target :: E(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgebrd_Tile_c(A,c_loc(D),c_loc(E),T)
  end subroutine CHAMELEON_zgebrd_Tile

  subroutine CHAMELEON_zgelqf_Tile(A,T,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgelqf_Tile_c(A,T)
  end subroutine CHAMELEON_zgelqf_Tile

  subroutine CHAMELEON_zgelqs_Tile(A,T,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgelqs_Tile_c(A,T,B)
  end subroutine CHAMELEON_zgelqs_Tile

  subroutine CHAMELEON_zgels_Tile(trans,A,T,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: trans
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgels_Tile_c(trans,A,T,B)
  end subroutine CHAMELEON_zgels_Tile

  subroutine CHAMELEON_zgemm_Tile(transA,transB,alpha,A,B,beta,C,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: transA
    integer(kind=c_int), intent(in) :: transB
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: C ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgemm_Tile_c(transA,transB,alpha,A,B,beta,C)
  end subroutine CHAMELEON_zgemm_Tile

  subroutine CHAMELEON_zgeqrf_Tile(A,T,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgeqrf_Tile_c(A,T)
  end subroutine CHAMELEON_zgeqrf_Tile

  subroutine CHAMELEON_zgeqrs_Tile(A,T,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgeqrs_Tile_c(A,T,B)
  end subroutine CHAMELEON_zgeqrs_Tile

  subroutine CHAMELEON_zgesv_Tile(A,IPIV,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(out), target :: IPIV(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgesv_Tile_c(A,c_loc(IPIV),B)
  end subroutine CHAMELEON_zgesv_Tile

  subroutine CHAMELEON_zgesv_incpiv_Tile(A,L,IPIV,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: IPIV ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: L ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgesv_incpiv_Tile_c(A,L,IPIV,B)
  end subroutine CHAMELEON_zgesv_incpiv_Tile

  subroutine CHAMELEON_zgesvd_Tile(jobu,jobvt,A,S,U,VT,T,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: jobu
    integer(kind=c_int), intent(in) :: jobvt
    real(kind=c_double), intent(out), target :: S(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: U ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: VT ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgesvd_Tile_c(jobu,jobvt,A,c_loc(S),U,VT,T)
  end subroutine CHAMELEON_zgesvd_Tile

  subroutine CHAMELEON_zgetrf_Tile(A,IPIV,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(out), target :: IPIV(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgetrf_Tile_c(A,c_loc(IPIV))
  end subroutine CHAMELEON_zgetrf_Tile

  subroutine CHAMELEON_zgetrf_incpiv_Tile(A,L,IPIV,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: IPIV ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: L ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgetrf_incpiv_Tile_c(A,L,IPIV)
  end subroutine CHAMELEON_zgetrf_incpiv_Tile

  subroutine CHAMELEON_zgetrf_nopiv_Tile(A,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgetrf_nopiv_Tile_c(A)
  end subroutine CHAMELEON_zgetrf_nopiv_Tile

  subroutine CHAMELEON_zgetrs_Tile(trans,A,IPIV,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: trans
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    integer(kind=c_int), intent(in), target :: IPIV(*)
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgetrs_Tile_c(trans,A,c_loc(IPIV),B)
  end subroutine CHAMELEON_zgetrs_Tile

  subroutine CHAMELEON_zgetrs_incpiv_Tile(A,L,IPIV,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: L ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: IPIV ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgetrs_incpiv_Tile_c(A,L,IPIV,B)
  end subroutine CHAMELEON_zgetrs_incpiv_Tile

#if defined(PRECISION_z) || defined(PRECISION_c)
  subroutine CHAMELEON_zhemm_Tile(side,uplo,alpha,A,B,beta,C,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: C ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zhemm_Tile_c(side,uplo,alpha,A,B,beta,C)
  end subroutine CHAMELEON_zhemm_Tile

  subroutine CHAMELEON_zherk_Tile(uplo,trans,alpha,A,beta,C,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: trans
    integer(kind=c_int), intent(in) :: uplo
    real(kind=c_double), intent(in) :: alpha
    real(kind=c_double), intent(in) :: beta
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: C ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zherk_Tile_c(uplo,trans,alpha,A,beta,C)
  end subroutine CHAMELEON_zherk_Tile

  subroutine CHAMELEON_zher2k_Tile(uplo,trans,alpha,A,B,beta,C,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: trans
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    real(kind=c_double), intent(in) :: beta
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: C ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zher2k_Tile_c(uplo,trans,alpha,A,B,beta,C)
  end subroutine CHAMELEON_zher2k_Tile
#endif

  subroutine CHAMELEON_zheev_Tile(jobz,uplo,A,W,T,Q,LDQ,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: jobz
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_int), intent(in) :: LDQ
    real(kind=c_double), intent(out), target :: W(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
    info = CHAMELEON_zheev_Tile_c(jobz,uplo,A,c_loc(W),T,c_loc(Q),LDQ)
  end subroutine CHAMELEON_zheev_Tile

  subroutine CHAMELEON_zheevd_Tile(jobz,uplo,A,W,T,Q,LDQ,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: jobz
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_int), intent(in) :: LDQ
    real(kind=c_double), intent(out), target :: W(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
    info = CHAMELEON_zheevd_Tile_c(jobz,uplo,A,c_loc(W),T,c_loc(Q),LDQ)
  end subroutine CHAMELEON_zheevd_Tile

  subroutine CHAMELEON_zhegv_Tile(itype,jobz,uplo,A,B,W,T,Q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: itype
    integer(kind=c_int), intent(in) :: jobz
    integer(kind=c_int), intent(in) :: uplo
    real(kind=c_double), intent(out), target :: W(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: Q ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zhegv_Tile_c(itype,jobz,uplo,A,B,c_loc(W),T,Q)
  end subroutine CHAMELEON_zhegv_Tile

  subroutine CHAMELEON_zhegvd_Tile(itype,jobz,uplo,A,B,W,T,Q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: itype
    integer(kind=c_int), intent(in) :: jobz
    integer(kind=c_int), intent(in) :: uplo
    real(kind=c_double), intent(out), target :: W(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: Q ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zhegvd_Tile_c(itype,jobz,uplo,A,B,c_loc(W),T,Q)
  end subroutine CHAMELEON_zhegvd_Tile

  subroutine CHAMELEON_zhegst_Tile(itype,uplo,A,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: itype
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zhegst_Tile_c(itype,uplo,A,B)
  end subroutine CHAMELEON_zhegst_Tile

  subroutine CHAMELEON_zhetrd_Tile(jobz,uplo,A,D,E,T,Q,LDQ,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: jobz
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_int), intent(in) :: LDQ
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    real(kind=c_double), intent(out), target :: D(*)
    real(kind=c_double), intent(out), target :: E(*)
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    complex(kind=c_double_complex), intent(out), target :: Q(LDQ,*)
    info = CHAMELEON_zhetrd_Tile_c(jobz,uplo,A,c_loc(D),c_loc(E),T,c_loc(Q),LDQ)
  end subroutine CHAMELEON_zhetrd_Tile

  function CHAMELEON_zlange_Tile(norm,A,work)
    use iso_c_binding
    implicit none
    real(kind=c_double) :: CHAMELEON_zlange_Tile
    integer(kind=c_int), intent(in) :: norm
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    real(kind=c_double), intent(inout), target :: work(*)
    CHAMELEON_zlange_Tile = CHAMELEON_zlange_Tile_c(norm,A,c_loc(work))
  end function CHAMELEON_zlange_Tile

#if defined(PRECISION_z) || defined(PRECISION_c)
  function CHAMELEON_zlanhe_Tile(norm,uplo,A,work)
    use iso_c_binding
    implicit none
    real(kind=c_double) :: CHAMELEON_zlanhe_Tile
    real(kind=c_double), intent(inout), target :: work(*)
    integer(kind=c_int), intent(in) :: norm
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    CHAMELEON_zlanhe_Tile = CHAMELEON_zlanhe_Tile_c(norm,uplo,A,c_loc(work))
  end function CHAMELEON_zlanhe_Tile
#endif

  function CHAMELEON_zlansy_Tile(norm,uplo,A,work)
    use iso_c_binding
    implicit none
    real(kind=c_double) :: CHAMELEON_zlansy_Tile
    real(kind=c_double), intent(inout), target :: work(*)
    integer(kind=c_int), intent(in) :: norm
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    CHAMELEON_zlansy_Tile = CHAMELEON_zlansy_Tile_c(norm,uplo,A,c_loc(work))
  end function CHAMELEON_zlansy_Tile

  subroutine CHAMELEON_zlaswp_Tile(A,K1,K2,IPIV,INCX,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: INCX
    integer(kind=c_int), intent(in) :: K1
    integer(kind=c_int), intent(in) :: K2
    integer(kind=c_int), intent(in), target :: IPIV(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zlaswp_Tile_c(A,K1,K2,c_loc(IPIV),INCX)
  end subroutine CHAMELEON_zlaswp_Tile

  subroutine CHAMELEON_zlauum_Tile(uplo,A,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zlauum_Tile_c(uplo,A)
  end subroutine CHAMELEON_zlauum_Tile

#if defined(PRECISION_z) || defined(PRECISION_c)
  subroutine CHAMELEON_zplghe_Tile(bump,uplo,A,seed,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    real(kind=c_double), intent(in) :: bump
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_long_long), intent(in) :: seed
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zplghe_Tile_c(bump,A,seed)
  end subroutine CHAMELEON_zplghe_Tile
#endif

  subroutine CHAMELEON_zplgsy_Tile(bump,uplo,A,seed,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    complex(kind=c_double_complex), intent(in) :: bump
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_long_long), intent(in) :: seed
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zplgsy_Tile_c(bump,A,seed)
  end subroutine CHAMELEON_zplgsy_Tile

  subroutine CHAMELEON_zplrnt_Tile(A,seed,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_long_long), intent(in) :: seed
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zplrnt_Tile_c(A,seed)
  end subroutine CHAMELEON_zplrnt_Tile

  subroutine CHAMELEON_zplrnk_Tile(K,C,seedA,seedB,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: K
    integer(kind=c_long_long), intent(in) :: seedA
    integer(kind=c_long_long), intent(in) :: seedB
    type(c_ptr), value :: C ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zplrnk_Tile_c(K,C,seedA,seedB)
  end subroutine CHAMELEON_zplrnk_Tile

  subroutine CHAMELEON_zposv_Tile(uplo,A,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zposv_Tile_c(uplo,A,B)
  end subroutine CHAMELEON_zposv_Tile

  subroutine CHAMELEON_zpotrf_Tile(uplo,A,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zpotrf_Tile_c(uplo,A)
  end subroutine CHAMELEON_zpotrf_Tile

  subroutine CHAMELEON_zsytrf_Tile(uplo,A,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zsytrf_Tile_c(uplo,A)
  end subroutine CHAMELEON_zsytrf_Tile

  subroutine CHAMELEON_zpotri_Tile(uplo,A,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zpotri_Tile_c(uplo,A)
  end subroutine CHAMELEON_zpotri_Tile

  subroutine CHAMELEON_zpotrs_Tile(uplo,A,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zpotrs_Tile_c(uplo,A,B)
  end subroutine CHAMELEON_zpotrs_Tile

  subroutine CHAMELEON_zsytrs_Tile(uplo,A,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zsytrs_Tile_c(uplo,A,B)
  end subroutine CHAMELEON_zsytrs_Tile

  subroutine CHAMELEON_zsymm_Tile(side,uplo,alpha,A,B,beta,C,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: C ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zsymm_Tile_c(side,uplo,alpha,A,B,beta,C)
  end subroutine CHAMELEON_zsymm_Tile

  subroutine CHAMELEON_zsyrk_Tile(uplo,trans,alpha,A,beta,C,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: trans
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: C ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zsyrk_Tile_c(uplo,trans,alpha,A,beta,C)
  end subroutine CHAMELEON_zsyrk_Tile

  subroutine CHAMELEON_zsyr2k_Tile(uplo,trans,alpha,A,B,beta,C,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: trans
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: C ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zsyr2k_Tile_c(uplo,trans,alpha,A,B,beta,C)
  end subroutine CHAMELEON_zsyr2k_Tile

  subroutine CHAMELEON_ztrmm_Tile(side,uplo,transA,diag,alpha,A,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: diag
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: transA
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_ztrmm_Tile_c(side,uplo,transA,diag,alpha,A,B)
  end subroutine CHAMELEON_ztrmm_Tile

  subroutine CHAMELEON_ztrsm_Tile(side,uplo,transA,diag,alpha,A,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: diag
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: transA
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_ztrsm_Tile_c(side,uplo,transA,diag,alpha,A,B)
  end subroutine CHAMELEON_ztrsm_Tile

  subroutine CHAMELEON_ztrsmpl_Tile(A,L,IPIV,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: L ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: IPIV ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_ztrsmpl_Tile_c(A,L,IPIV,B)
  end subroutine CHAMELEON_ztrsmpl_Tile

  subroutine CHAMELEON_ztrtri_Tile(uplo,diag,A,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: diag
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_ztrtri_Tile_c(uplo,diag,A)
  end subroutine CHAMELEON_ztrtri_Tile

  subroutine CHAMELEON_zunglq_Tile(A,T,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zunglq_Tile_c(A,T,B)
  end subroutine CHAMELEON_zunglq_Tile

  subroutine CHAMELEON_zungqr_Tile(A,T,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zungqr_Tile_c(A,T,B)
  end subroutine CHAMELEON_zungqr_Tile

  subroutine CHAMELEON_zunmlq_Tile(side,trans,A,T,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: trans
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zunmlq_Tile_c(side,trans,A,T,B)
  end subroutine CHAMELEON_zunmlq_Tile

  subroutine CHAMELEON_zunmqr_Tile(side,trans,A,T,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: trans
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zunmqr_Tile_c(side,trans,A,T,B)
  end subroutine CHAMELEON_zunmqr_Tile

  subroutine CHAMELEON_zgetri_Tile(A,IPIV,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in), target :: IPIV(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgetri_Tile_c(A,c_loc(IPIV))
  end subroutine CHAMELEON_zgetri_Tile

  subroutine CHAMELEON_zlacpy_Tile(uplo,A,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zlacpy_Tile_c(uplo,A,B)
  end subroutine CHAMELEON_zlacpy_Tile

  subroutine CHAMELEON_zlaset_Tile(uplo,alpha,beta,A,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zlaset_Tile_c(uplo,alpha,beta,A)
  end subroutine CHAMELEON_zlaset_Tile

  subroutine CHAMELEON_zlaswpc_Tile(A,K1,K2,IPIV,INCX,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in), target :: IPIV(*)
    integer(kind=c_int), intent(in) :: INCX
    integer(kind=c_int), intent(in) :: K1
    integer(kind=c_int), intent(in) :: K2
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zlaswpc_Tile_c(A,K1,K2,c_loc(IPIV),INCX)
  end subroutine CHAMELEON_zlaswpc_Tile

  subroutine CHAMELEON_ztrsmrv_Tile(side,uplo,transA,diag,alpha,A,B,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: diag
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: transA
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_ztrsmrv_Tile_c(side,uplo,transA,diag,alpha,A,B)
  end subroutine CHAMELEON_ztrsmrv_Tile

  subroutine CHAMELEON_zgetri_Tile_Async(A,IPIV,W,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in), target :: IPIV(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: W ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgetri_Tile_Async_c(A,c_loc(IPIV),W,sequence,request)
  end subroutine CHAMELEON_zgetri_Tile_Async

  subroutine CHAMELEON_zlange_Tile_Async(norm,A,work,value,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int) :: info
    real(kind=c_double), intent(out), target :: value
    real(kind=c_double), intent(inout), target :: work(*)
    integer(kind=c_int), intent(in) :: norm
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zlange_Tile_Async_c(norm,A,c_loc(work),c_loc(value),sequence,request)
  end subroutine CHAMELEON_zlange_Tile_Async

  subroutine CHAMELEON_zlansy_Tile_Async(norm,uplo,A,work,value,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int) :: info
    real(kind=c_double), intent(out), target :: value
    real(kind=c_double), intent(inout), target :: work(*)
    integer(kind=c_int), intent(in) :: norm
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zlansy_Tile_Async_c(norm,uplo,A,c_loc(work),c_loc(value),sequence,request)
  end subroutine CHAMELEON_zlansy_Tile_Async

  subroutine CHAMELEON_zgebrd_Tile_Async(A,D,E,T,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    real(kind=c_double), intent(out), target :: D(*)
    real(kind=c_double), intent(out), target :: E(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgebrd_Tile_Async_c(A,c_loc(D),c_loc(E),T,sequence,request)
  end subroutine CHAMELEON_zgebrd_Tile_Async

  subroutine CHAMELEON_zgecfi_Async(m,n,A,fin,imb,inb,fout,omb,onb,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    complex(kind=c_double_complex), intent(inout), target :: A(*)
    integer(kind=c_int), intent(in) :: fin
    integer(kind=c_int), intent(in) :: fout
    integer(kind=c_int), intent(in) :: imb
    integer(kind=c_int), intent(in) :: inb
    integer(kind=c_int), intent(in) :: omb
    integer(kind=c_int), intent(in) :: onb
    integer(kind=c_int), intent(in) :: m
    integer(kind=c_int), intent(in) :: n
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgecfi_Async_c(m,n,c_loc(A),fin,imb,inb,fout,omb,onb,sequence,request)
  end subroutine CHAMELEON_zgecfi_Async

  subroutine CHAMELEON_zgelqf_Tile_Async(A,T,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgelqf_Tile_Async_c(A,T,sequence,request)
  end subroutine CHAMELEON_zgelqf_Tile_Async

  subroutine CHAMELEON_zgelqs_Tile_Async(A,T,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgelqs_Tile_Async_c(A,T,B,sequence,request)
  end subroutine CHAMELEON_zgelqs_Tile_Async

  subroutine CHAMELEON_zgels_Tile_Async(trans,A,T,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: trans
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgels_Tile_Async_c(trans,A,T,B,sequence,request)
  end subroutine CHAMELEON_zgels_Tile_Async

  subroutine CHAMELEON_zgemm_Tile_Async(transA,transB,alpha,A,B,beta,C,WS,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: transA
    integer(kind=c_int), intent(in) :: transB
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: C ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: WS ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgemm_Tile_Async_c(transA,transB,alpha,A,B,beta,C,WS,sequence,request)
  end subroutine CHAMELEON_zgemm_Tile_Async

  subroutine CHAMELEON_zgeqrf_Tile_Async(A,T,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgeqrf_Tile_Async_c(A,T,sequence,request)
  end subroutine CHAMELEON_zgeqrf_Tile_Async

  subroutine CHAMELEON_zgeqrs_Tile_Async(A,T,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgeqrs_Tile_Async_c(A,T,B,sequence,request)
  end subroutine CHAMELEON_zgeqrs_Tile_Async

  subroutine CHAMELEON_zgesv_Tile_Async(A,IPIV,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(out), target :: IPIV(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgesv_Tile_Async_c(A,c_loc(IPIV),B,sequence,request)
  end subroutine CHAMELEON_zgesv_Tile_Async

  subroutine CHAMELEON_zgesv_incpiv_Tile_Async(A,L,IPIV,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: IPIV ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: L ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgesv_incpiv_Tile_Async_c(A,L,IPIV,B,sequence,request)
  end subroutine CHAMELEON_zgesv_incpiv_Tile_Async

  subroutine CHAMELEON_zgesvd_Tile_Async(jobu,jobvt,A,S,U,VT,T,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: jobu
    integer(kind=c_int), intent(in) :: jobvt
    real(kind=c_double), intent(out), target :: S(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: U ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: VT ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgesvd_Tile_Async_c(jobu,jobvt,A,c_loc(S),U,VT,T,sequence,request)
  end subroutine CHAMELEON_zgesvd_Tile_Async

  subroutine CHAMELEON_zgetmi_Async(m,n,A,fin,mb,nb,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    complex(kind=c_double_complex), intent(inout), target :: A(*)
    integer(kind=c_int), intent(in) :: fin
    integer(kind=c_int), intent(in) :: mb
    integer(kind=c_int), intent(in) :: nb
    integer(kind=c_int), intent(in) :: m
    integer(kind=c_int), intent(in) :: n
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgetmi_Async_c(m,n,c_loc(A),fin,mb,nb,sequence,request)
  end subroutine CHAMELEON_zgetmi_Async

  subroutine CHAMELEON_zgetrf_Tile_Async(A,IPIV,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(out), target :: IPIV(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgetrf_Tile_Async_c(A,c_loc(IPIV),sequence,request)
  end subroutine CHAMELEON_zgetrf_Tile_Async

  subroutine CHAMELEON_zgetrf_incpiv_Tile_Async(A,L,IPIV,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: IPIV ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: L ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgetrf_incpiv_Tile_Async_c(A,L,IPIV,sequence,request)
  end subroutine CHAMELEON_zgetrf_incpiv_Tile_Async

  subroutine CHAMELEON_zgetrf_nopiv_Tile_Async(A,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgetrf_nopiv_Tile_Async_c(A,sequence,request)
  end subroutine CHAMELEON_zgetrf_nopiv_Tile_Async

  subroutine CHAMELEON_zgetrs_Tile_Async(trans,A,IPIV,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: trans
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    integer(kind=c_int), intent(in), target :: IPIV(*)
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgetrs_Tile_Async_c(trans,A,c_loc(IPIV),B,sequence,request)
  end subroutine CHAMELEON_zgetrs_Tile_Async

  subroutine CHAMELEON_zgetrs_incpiv_Tile_Async(A,L,IPIV,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: L ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: IPIV ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zgetrs_incpiv_Tile_Async_c(A,L,IPIV,B,sequence,request)
  end subroutine CHAMELEON_zgetrs_incpiv_Tile_Async

#if defined(PRECISION_z) || defined(PRECISION_c)
  subroutine CHAMELEON_zhemm_Tile_Async(side,uplo,alpha,A,B,beta,C,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: C ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zhemm_Tile_Async_c(side,uplo,alpha,A,B,beta,C,sequence,request)
  end subroutine CHAMELEON_zhemm_Tile_Async

  subroutine CHAMELEON_zherk_Tile_Async(uplo,trans,alpha,A,beta,C,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: trans
    integer(kind=c_int), intent(in) :: uplo
    real(kind=c_double), intent(in) :: alpha
    real(kind=c_double), intent(in) :: beta
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: C ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zherk_Tile_Async_c(uplo,trans,alpha,A,beta,C,sequence,request)
  end subroutine CHAMELEON_zherk_Tile_Async

  subroutine CHAMELEON_zher2k_Tile_Async(uplo,trans,alpha,A,B,beta,C,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: trans
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    real(kind=c_double), intent(in) :: beta
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: C ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zher2k_Tile_Async_c(uplo,trans,alpha,A,B,beta,C,sequence,request)
  end subroutine CHAMELEON_zher2k_Tile_Async
#endif

  subroutine CHAMELEON_zheev_Tile_Async(jobz,uplo,A,W,T,Q,LDQ,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: jobz
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_int), intent(in) :: LDQ
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    real(kind=c_double), intent(out), target :: W(*)
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    complex(kind=c_double_complex), intent(out), target :: Q(LDQ, *)
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zheev_Tile_Async_c(jobz,uplo,A,c_loc(W),T,c_loc(Q),LDQ,sequence,request)
  end subroutine CHAMELEON_zheev_Tile_Async

  subroutine CHAMELEON_zheevd_Tile_Async(jobz,uplo,A,W,T,Q,LDQ,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: jobz
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_int), intent(in) :: LDQ
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    real(kind=c_double), intent(out), target :: W(*)
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    complex(kind=c_double_complex), intent(out), target :: Q(LDQ, *)
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zheevd_Tile_Async_c(jobz,uplo,A,c_loc(W),T,c_loc(Q),LDQ,sequence,request)
  end subroutine CHAMELEON_zheevd_Tile_Async

  subroutine CHAMELEON_zhegv_Tile_Async(itype,jobz,uplo,A,B,W,T,Q,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: itype
    integer(kind=c_int), intent(in) :: jobz
    integer(kind=c_int), intent(in) :: uplo
    real(kind=c_double), intent(out), target :: W(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: Q ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zhegv_Tile_Async_c(itype,jobz,uplo,A,B,c_loc(W),T,Q,sequence,request)
  end subroutine CHAMELEON_zhegv_Tile_Async

  subroutine CHAMELEON_zhegvd_Tile_Async(itype,jobz,uplo,A,B,W,T,Q,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: itype
    integer(kind=c_int), intent(in) :: jobz
    integer(kind=c_int), intent(in) :: uplo
    real(kind=c_double), intent(out), target :: W(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: Q ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zhegvd_Tile_Async_c(itype,jobz,uplo,A,B,c_loc(W),T,Q,sequence,request)
  end subroutine CHAMELEON_zhegvd_Tile_Async

  subroutine CHAMELEON_zhegst_Tile_Async(itype,uplo,A,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: itype
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zhegst_Tile_Async_c(itype,uplo,A,B,sequence,request)
  end subroutine CHAMELEON_zhegst_Tile_Async

  subroutine CHAMELEON_zhetrd_Tile_Async(jobz,uplo,A,D,E,T,Q,LDQ,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: jobz
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_int), intent(in) :: LDQ
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    real(kind=c_double), intent(out), target :: D(*)
    real(kind=c_double), intent(out), target :: E(*)
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    complex(kind=c_double_complex), intent(out), target :: Q(LDQ, *)
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zhetrd_Tile_Async_c(jobz,uplo,A,c_loc(D),c_loc(E),T,c_loc(Q),LDQ,sequence,request)
  end subroutine CHAMELEON_zhetrd_Tile_Async

#if defined(PRECISION_z) || defined(PRECISION_c)
  subroutine CHAMELEON_zlanhe_Tile_Async(norm,uplo,A,work,value,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    real(kind=c_double), intent(inout), target :: work(*)
    integer(kind=c_int), intent(in) :: norm
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_int), intent(out), target :: value
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zlanhe_Tile_Async_c(norm,uplo,A,c_loc(work),c_loc(value),sequence,request)
  end subroutine CHAMELEON_zlanhe_Tile_Async
#endif

  subroutine CHAMELEON_zlaswp_Tile_Async(A,K1,K2,IPIV,INCX,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: INCX
    integer(kind=c_int), intent(in) :: K1
    integer(kind=c_int), intent(in) :: K2
    integer(kind=c_int), intent(in), target :: IPIV(*)
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zlaswp_Tile_Async_c(A,K1,K2,c_loc(IPIV),INCX,sequence,request)
  end subroutine CHAMELEON_zlaswp_Tile_Async

  subroutine CHAMELEON_zlauum_Tile_Async(uplo,A,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zlauum_Tile_Async_c(uplo,A,sequence,request)
  end subroutine CHAMELEON_zlauum_Tile_Async

#if defined(PRECISION_z) || defined(PRECISION_c)
  subroutine CHAMELEON_zplghe_Tile_Async(bump,uplo,A,seed,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    real(kind=c_double), intent(in) :: bump
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_long_long), intent(in) :: seed
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zplghe_Tile_Async_c(bump,A,seed,sequence,request)
  end subroutine CHAMELEON_zplghe_Tile_Async
#endif

  subroutine CHAMELEON_zplgsy_Tile_Async(bump,uplo,A,seed,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    complex(kind=c_double_complex), intent(in) :: bump
    integer(kind=c_int), intent(in) :: uplo
    integer(kind=c_long_long), intent(in) :: seed
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zplgsy_Tile_Async_c(bump,A,seed,sequence,request)
  end subroutine CHAMELEON_zplgsy_Tile_Async

  subroutine CHAMELEON_zplrnt_Tile_Async(A,seed,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_long_long), intent(in) :: seed
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zplrnt_Tile_Async_c(A,seed,sequence,request)
  end subroutine CHAMELEON_zplrnt_Tile_Async

  subroutine CHAMELEON_zplrnk_Tile_Async(K,C,seedA,seedB,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: K
    integer(kind=c_long_long), intent(in) :: seedA
    integer(kind=c_long_long), intent(in) :: seedB
    type(c_ptr), value :: C ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zplrnk_Tile_Async_c(K,C,seedA,seedB,sequence,request)
  end subroutine CHAMELEON_zplrnk_Tile_Async

  subroutine CHAMELEON_zposv_Tile_Async(uplo,A,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zposv_Tile_Async_c(uplo,A,B,sequence,request)
  end subroutine CHAMELEON_zposv_Tile_Async

  subroutine CHAMELEON_zpotrf_Tile_Async(uplo,A,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zpotrf_Tile_Async_c(uplo,A,sequence,request)
  end subroutine CHAMELEON_zpotrf_Tile_Async

  subroutine CHAMELEON_zsytrf_Tile_Async(uplo,A,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zsytrf_Tile_Async_c(uplo,A,sequence,request)
  end subroutine CHAMELEON_zsytrf_Tile_Async

  subroutine CHAMELEON_zpotri_Tile_Async(uplo,A,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zpotri_Tile_Async_c(uplo,A,sequence,request)
  end subroutine CHAMELEON_zpotri_Tile_Async

  subroutine CHAMELEON_zpotrs_Tile_Async(uplo,A,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zpotrs_Tile_Async_c(uplo,A,B,sequence,request)
  end subroutine CHAMELEON_zpotrs_Tile_Async

  subroutine CHAMELEON_zsytrs_Tile_Async(uplo,A,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zsytrs_Tile_Async_c(uplo,A,B,sequence,request)
  end subroutine CHAMELEON_zsytrs_Tile_Async

  subroutine CHAMELEON_zsymm_Tile_Async(side,uplo,alpha,A,B,beta,C,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: C ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zsymm_Tile_Async_c(side,uplo,alpha,A,B,beta,C,sequence,request)
  end subroutine CHAMELEON_zsymm_Tile_Async

  subroutine CHAMELEON_zsyrk_Tile_Async(uplo,trans,alpha,A,beta,C,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: trans
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: C ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zsyrk_Tile_Async_c(uplo,trans,alpha,A,beta,C,sequence,request)
  end subroutine CHAMELEON_zsyrk_Tile_Async

  subroutine CHAMELEON_zsyr2k_Tile_Async(uplo,trans,alpha,A,B,beta,C,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: trans
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: C ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zsyr2k_Tile_Async_c(uplo,trans,alpha,A,B,beta,C,sequence,request)
  end subroutine CHAMELEON_zsyr2k_Tile_Async

  subroutine CHAMELEON_ztrmm_Tile_Async(side,uplo,transA,diag,alpha,A,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: diag
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: transA
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_ztrmm_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request)
  end subroutine CHAMELEON_ztrmm_Tile_Async

  subroutine CHAMELEON_ztrsm_Tile_Async(side,uplo,transA,diag,alpha,A,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: diag
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: transA
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_ztrsm_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request)
  end subroutine CHAMELEON_ztrsm_Tile_Async

  subroutine CHAMELEON_ztrsmpl_Tile_Async(A,L,IPIV,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: L ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: IPIV ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_ztrsmpl_Tile_Async_c(A,L,IPIV,B,sequence,request)
  end subroutine CHAMELEON_ztrsmpl_Tile_Async

  subroutine CHAMELEON_ztrtri_Tile_Async(uplo,diag,A,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: diag
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_ztrtri_Tile_Async_c(uplo,diag,A,sequence,request)
  end subroutine CHAMELEON_ztrtri_Tile_Async

  subroutine CHAMELEON_zunglq_Tile_Async(A,T,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zunglq_Tile_Async_c(A,T,B,sequence,request)
  end subroutine CHAMELEON_zunglq_Tile_Async

  subroutine CHAMELEON_zungqr_Tile_Async(A,T,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zungqr_Tile_Async_c(A,T,B,sequence,request)
  end subroutine CHAMELEON_zungqr_Tile_Async

  subroutine CHAMELEON_zunmlq_Tile_Async(side,trans,A,T,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: trans
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zunmlq_Tile_Async_c(side,trans,A,T,B,sequence,request)
  end subroutine CHAMELEON_zunmlq_Tile_Async

  subroutine CHAMELEON_zunmqr_Tile_Async(side,trans,A,T,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: trans
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: T ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zunmqr_Tile_Async_c(side,trans,A,T,B,sequence,request)
  end subroutine CHAMELEON_zunmqr_Tile_Async

  subroutine CHAMELEON_zlacpy_Tile_Async(uplo,A,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zlacpy_Tile_Async_c(uplo,A,B,sequence,request)
  end subroutine CHAMELEON_zlacpy_Tile_Async

  subroutine CHAMELEON_zlaset_Tile_Async(uplo,alpha,beta,A,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    complex(kind=c_double_complex), intent(in) :: beta
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zlaset_Tile_Async_c(uplo,alpha,beta,A,sequence,request)
  end subroutine CHAMELEON_zlaset_Tile_Async

  subroutine CHAMELEON_zlaswpc_Tile_Async(A,K1,K2,IPIV,INCX,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in), target :: IPIV(*)
    integer(kind=c_int), intent(in) :: INCX
    integer(kind=c_int), intent(in) :: K1
    integer(kind=c_int), intent(in) :: K2
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_zlaswpc_Tile_Async_c(A,K1,K2,c_loc(IPIV),INCX,sequence,request)
  end subroutine CHAMELEON_zlaswpc_Tile_Async

  subroutine CHAMELEON_ztrsmrv_Tile_Async(side,uplo,transA,diag,alpha,A,B,sequence,request,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: diag
    integer(kind=c_int), intent(in) :: side
    integer(kind=c_int), intent(in) :: transA
    integer(kind=c_int), intent(in) :: uplo
    complex(kind=c_double_complex), intent(in) :: alpha
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: B ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: sequence ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: request ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_ztrsmrv_Tile_Async_c(side,uplo,transA,diag,alpha,A,B,sequence,request)
  end subroutine CHAMELEON_ztrsmrv_Tile_Async

  subroutine CHAMELEON_Alloc_Workspace_zgelqf(M,N,T,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: T ! T is **, so pass by reference
    info = CHAMELEON_Alloc_Workspace_zgelqf_c(M,N,T)
  end subroutine CHAMELEON_Alloc_Workspace_zgelqf

  subroutine CHAMELEON_Alloc_Workspace_zgels(M,N,T,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: T ! T is **, so pass by reference
    info = CHAMELEON_Alloc_Workspace_zgels_c(M,N,T)
  end subroutine CHAMELEON_Alloc_Workspace_zgels

  subroutine CHAMELEON_Alloc_Workspace_zgeqrf(M,N,T,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: T ! T is **, so pass by reference
    info = CHAMELEON_Alloc_Workspace_zgeqrf_c(M,N,T)
  end subroutine CHAMELEON_Alloc_Workspace_zgeqrf

  subroutine CHAMELEON_Alloc_Workspace_zgesv_incpiv(N,L,IPIV,p,q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
    type(c_ptr) :: L ! L is **, so pass by reference
    integer(kind=c_int), value :: p
    integer(kind=c_int), value :: q
    info = CHAMELEON_Alloc_Workspace_zgesv_incpiv_c(N,L,IPIV,p,q)
  end subroutine CHAMELEON_Alloc_Workspace_zgesv_incpiv

  subroutine CHAMELEON_Alloc_Workspace_zgetrf_incpiv(M,N,L,IPIV,p,q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
    type(c_ptr) :: L ! L is **, so pass by reference
    integer(kind=c_int), value :: p
    integer(kind=c_int), value :: q
    info = CHAMELEON_Alloc_Workspace_zgetrf_incpiv_c(M,N,L,IPIV,p,q)
  end subroutine CHAMELEON_Alloc_Workspace_zgetrf_incpiv

  subroutine CHAMELEON_Alloc_Workspace_zgeev(N,T,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: T ! T is **, so pass by reference
    info = CHAMELEON_Alloc_Workspace_zgeev_c(N,T)
  end subroutine CHAMELEON_Alloc_Workspace_zgeev

  subroutine CHAMELEON_Alloc_Workspace_zgebrd(M,N,T,p,q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: T ! T is **, so pass by reference
    integer(kind=c_int), value :: p
    integer(kind=c_int), value :: q
    info = CHAMELEON_Alloc_Workspace_zgebrd_c(M,N,T,p,q)
  end subroutine CHAMELEON_Alloc_Workspace_zgebrd

  subroutine CHAMELEON_Alloc_Workspace_zgehrd(N,T,p,q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: T ! T is **, so pass by reference
    integer(kind=c_int), value :: p
    integer(kind=c_int), value :: q
    info = CHAMELEON_Alloc_Workspace_zgehrd_c(N,T,p,q)
  end subroutine CHAMELEON_Alloc_Workspace_zgehrd

  subroutine CHAMELEON_Alloc_Workspace_zgesvd(M,N,T,p,q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: T ! T is **, so pass by reference
    integer(kind=c_int), value :: p
    integer(kind=c_int), value :: q
    info = CHAMELEON_Alloc_Workspace_zgesvd_c(M,N,T,p,q)
  end subroutine CHAMELEON_Alloc_Workspace_zgesvd

  subroutine CHAMELEON_Alloc_Workspace_zheev(M,N,T,p,q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: T ! T is **, so pass by reference
    integer(kind=c_int), value :: p
    integer(kind=c_int), value :: q
    info = CHAMELEON_Alloc_Workspace_zheev_c(M,N,T,p,q)
  end subroutine CHAMELEON_Alloc_Workspace_zheev

  subroutine CHAMELEON_Alloc_Workspace_zheevd(M,N,T,p,q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: T ! T is **, so pass by reference
    integer(kind=c_int), value :: p
    integer(kind=c_int), value :: q
    info = CHAMELEON_Alloc_Workspace_zheevd_c(M,N,T,p,q)
  end subroutine CHAMELEON_Alloc_Workspace_zheevd

  subroutine CHAMELEON_Alloc_Workspace_zhegv(M,N,T,p,q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: T ! T is **, so pass by reference
    integer(kind=c_int), value :: p
    integer(kind=c_int), value :: q
    info = CHAMELEON_Alloc_Workspace_zhegv_c(M,N,T,p,q)
  end subroutine CHAMELEON_Alloc_Workspace_zhegv

  subroutine CHAMELEON_Alloc_Workspace_zhegvd(M,N,T,p,q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: T ! T is **, so pass by reference
    integer(kind=c_int), value :: p
    integer(kind=c_int), value :: q
    info = CHAMELEON_Alloc_Workspace_zhegvd_c(M,N,T,p,q)
  end subroutine CHAMELEON_Alloc_Workspace_zhegvd

  subroutine CHAMELEON_Alloc_Workspace_zhetrd(M,N,T,p,q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: T ! T is **, so pass by reference
    integer(kind=c_int), value :: p
    integer(kind=c_int), value :: q
    info = CHAMELEON_Alloc_Workspace_zhetrd_c(M,N,T,p,q)
  end subroutine CHAMELEON_Alloc_Workspace_zhetrd

  subroutine CHAMELEON_Alloc_Workspace_zgelqf_Tile(M,N,descT,p,q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: descT ! descT is **, so pass by reference
    integer(kind=c_int), value :: p
    integer(kind=c_int), value :: q
    info = CHAMELEON_Alloc_Workspace_zgelqf_Tile_c(M,N,descT,p,q)
  end subroutine CHAMELEON_Alloc_Workspace_zgelqf_Tile

  subroutine CHAMELEON_Alloc_Workspace_zgels_Tile(M,N,descT,p,q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: descT ! descT is **, so pass by reference
    integer(kind=c_int), value :: p
    integer(kind=c_int), value :: q
    info = CHAMELEON_Alloc_Workspace_zgels_Tile_c(M,N,descT,p,q)
  end subroutine CHAMELEON_Alloc_Workspace_zgels_Tile

  subroutine CHAMELEON_Alloc_Workspace_zgeqrf_Tile(M,N,descT,p,q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: M
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: descT ! descT is **, so pass by reference
    integer(kind=c_int), value :: p
    integer(kind=c_int), value :: q
    info = CHAMELEON_Alloc_Workspace_zgeqrf_Tile_c(M,N,descT,p,q)
  end subroutine CHAMELEON_Alloc_Workspace_zgeqrf_Tile

  subroutine CHAMELEON_Alloc_Workspace_zgesv_incpiv_Tile(N,descL,IPIV,p,q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
    type(c_ptr) :: descL ! descL is **, so pass by reference
    integer(kind=c_int), value :: p
    integer(kind=c_int), value :: q
    info = CHAMELEON_Alloc_Workspace_zgesv_incpiv_Tile_c(N,descL,IPIV,p,q)
  end subroutine CHAMELEON_Alloc_Workspace_zgesv_incpiv_Tile

  subroutine CHAMELEON_Alloc_Workspace_zgetrf_incpiv_Tile(N,descL,IPIV,p,q,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    integer(kind=c_int), intent(in) :: N
    type(c_ptr) :: IPIV ! IPIV is **, so pass by reference
    type(c_ptr) :: descL ! descL is **, so pass by reference
    integer(kind=c_int), value :: p
    integer(kind=c_int), value :: q
    info = CHAMELEON_Alloc_Workspace_zgetrf_incpiv_Tile_c(N,descL,IPIV,p,q)
  end subroutine CHAMELEON_Alloc_Workspace_zgetrf_incpiv_Tile

  subroutine CHAMELEON_Alloc_Workspace_zgetri_Tile_Async(A,W,info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(out) :: info
    type(c_ptr), value :: A ! Arg managed by CHAMELEON: opaque to Fortran
    type(c_ptr), value :: W ! Arg managed by CHAMELEON: opaque to Fortran
    info = CHAMELEON_Alloc_Workspace_zgetri_Tile_Async_c(A,W)
  end subroutine CHAMELEON_Alloc_Workspace_zgetri_Tile_Async

end module chameleon_z
