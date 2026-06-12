!**
!
! @file chameleon_fortran.h
! @brief Chameleon Fortran77 header
!
! @copyright 2009-2014 The University of Tennessee and The University of
!                      Tennessee Research Foundation. All rights reserved.
! @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 1.4.0
! @author Bilel Hadri
! @author Mathieu Faverge
! @author Cedric Castagnede
! @author Florent Pruvost
! @date 2024-02-18
!
!**
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
!**

!********************************************************************
!   CHAMELEON constants - precisions
!
      integer  ChamByte, ChamInteger, ChamRealFloat
      integer  ChamRealDouble, ChamComplexFloat, ChamComplexDouble
      parameter ( ChamByte          = 0 )
      parameter ( ChamInteger       = 1 )
      parameter ( ChamRealFloat     = 2 )
      parameter ( ChamRealDouble    = 3 )
      parameter ( ChamComplexFloat  = 4 )
      parameter ( ChamComplexDouble = 5 )

!********************************************************************
!   CHAMELEON constants - CBLAS & LAPACK
!
      integer ChamCM, ChamRM, ChamCCRB
      integer ChamCRRB, ChamRCRB, ChamRRRB
      parameter ( ChamCM         = 101 )
      parameter ( ChamRM         = 102 )
      parameter ( ChamCCRB       = 103 )
      parameter ( ChamCRRB       = 104 )
      parameter ( ChamRCRB       = 105 )
      parameter ( ChamRRRB       = 106 )

      integer  ChamNoTrans, ChamTrans, ChamConjTrans
      parameter ( ChamNoTrans    = 111 )
      parameter ( ChamTrans      = 112 )
      parameter ( ChamConjTrans  = 113 )

      integer ChamUpper, ChamLower
      integer ChamUpperLower
      parameter ( ChamUpper      = 121 )
      parameter ( ChamLower      = 122 )
      parameter ( ChamUpperLower = 123 )

      integer ChamNonUnit,ChamUnit
      parameter ( ChamNonUnit    = 131 )
      parameter ( ChamUnit       = 132 )

      integer ChamLeft,ChamRight
      parameter ( ChamLeft       = 141 )
      parameter ( ChamRight      = 142 )

      integer ChamOneNorm, ChamRealOneNorm
      integer ChamTwoNorm, ChamFrobeniusNorm
      integer ChamInfNorm, ChamRealInfNorm
      integer ChamMaxNorm, ChamRealMaxNorm
      parameter ( ChamOneNorm       = 171 )
      parameter ( ChamRealOneNorm   = 172 )
      parameter ( ChamTwoNorm       = 173 )
      parameter ( ChamFrobeniusNorm = 174 )
      parameter ( ChamInfNorm       = 175 )
      parameter ( ChamRealInfNorm   = 176 )
      parameter ( ChamMaxNorm       = 177 )
      parameter ( ChamRealMaxNorm   = 178 )

      integer ChamDistUniform
      integer ChamDistSymmetric
      integer ChamDistNormal
      parameter ( ChamDistUniform   = 201 )
      parameter ( ChamDistSymmetric = 202 )
      parameter ( ChamDistNormal    = 203 )

      integer ChamHermGeev
      integer ChamHermPoev
      integer ChamNonsymPosv
      integer ChamSymPosv
      parameter ( ChamHermGeev    = 241 )
      parameter ( ChamHermPoev    = 242 )
      parameter ( ChamNonsymPosv  = 243 )
      parameter ( ChamSymPosv     = 244 )

      integer ChamNoPacking
      integer ChamPackSubdiag
      integer ChamPackSupdiag
      integer ChamPackColumn
      integer ChamPackLowerBand
      integer ChamPackRow
      integer ChamPackUpeprBand
      integer ChamPackAll
      parameter ( ChamNoPacking     = 291 )
      parameter ( ChamPackSubdiag   = 292 )
      parameter ( ChamPackSupdiag   = 293 )
      parameter ( ChamPackColumn    = 294 )
      parameter ( ChamPackRow       = 295 )
      parameter ( ChamPackLowerBand = 296 )
      parameter ( ChamPackUpeprBand = 297 )
      parameter ( ChamPackAll       = 298 )

      integer ChamNoVec,ChamVec,ChamIvec
      parameter ( ChamNoVec = 301 )
      parameter ( ChamVec   = 302 )
      parameter ( ChamIvec  = 303 )

      integer ChamDirForward, ChamDirBackward
      parameter ( ChamDirForward    = 391 )
      parameter ( ChamDirBackward   = 392 )

      integer ChamColumnwise,ChamRowwise
      parameter ( ChamColumnwise = 401 )
      parameter ( ChamRowwise    = 402 )

!********************************************************************
!   CHAMELEON constants - boolean
!
      integer CHAMELEON_FALSE, CHAMELEON_TRUE
      parameter ( CHAMELEON_FALSE = 0 )
      parameter ( CHAMELEON_TRUE  = 1 )

!********************************************************************
!   State machine switches
!
      integer CHAMELEON_WARNINGS, CHAMELEON_ERRORS, CHAMELEON_AUTOTUNING
      integer CHAMELEON_DAG, CHAMELEON_GENERATE_TRACE, CHAMELEON_GENERATE_STATS
      integer CHAMELEON_PARALLEL_MODE, CHAMELEON_BOUND, CHAMELEON_PROGRESS
      integer CHAMELEON_GEMM3M, CHAMELEON_GENERIC
      parameter ( CHAMELEON_WARNINGS       = 1 )
      parameter ( CHAMELEON_ERRORS         = 2 )
      parameter ( CHAMELEON_AUTOTUNING     = 3 )
      parameter ( CHAMELEON_DAG            = 4 )
      parameter ( CHAMELEON_GENERATE_TRACE = 5 )
      parameter ( CHAMELEON_GENERATE_STATS = 6 )
      parameter ( CHAMELEON_PARALLEL_MODE  = 7 )
      parameter ( CHAMELEON_BOUND          = 8 )
      parameter ( CHAMELEON_PROGRESS       = 9 )
      parameter ( CHAMELEON_GEMM3M         = 10 )
      parameter ( CHAMELEON_GENERIC        = 11 )

!********************************************************************
!   CHAMELEON constants - configuration  parameters
!
      integer CHAMELEON_CONCURRENCY, CHAMELEON_TILE_SIZE
      integer CHAMELEON_INNER_BLOCK_SIZE, CHAMELEON_SCHEDULING_MODE
      integer CHAMELEON_HOUSEHOLDER_MODE, CHAMELEON_HOUSEHOLDER_SIZE
      integer CHAMELEON_TRANSLATION_MODE
      parameter ( CHAMELEON_CONCURRENCY      = 1 )
      parameter ( CHAMELEON_TILE_SIZE        = 2 )
      parameter ( CHAMELEON_INNER_BLOCK_SIZE = 3 )
      parameter ( CHAMELEON_SCHEDULING_MODE  = 4 )
      parameter ( CHAMELEON_HOUSEHOLDER_MODE = 5 )
      parameter ( CHAMELEON_HOUSEHOLDER_SIZE = 6 )
      parameter ( CHAMELEON_TRANSLATION_MODE = 7 )

!********************************************************************
!   CHAMELEON constants - householder mode
!
      integer ChamFlatHouseholder, ChamTreeHouseholder
      parameter ( ChamFlatHouseholder  = 1 )
      parameter ( ChamTreeHouseholder  = 2 )

!*********************************************************************
!   CHAMELEON constants - translation mode
!
      integer ChamInPlace, ChamOutOfPlace
      parameter ( ChamInPlace     = 1 )
      parameter ( ChamOutOfPlace  = 2 )

!********************************************************************
!   CHAMELEON constants - success & error codes
!
      integer CHAMELEON_SUCCESS, CHAMELEON_ERR_NOT_INITIALIZED
      integer CHAMELEON_ERR_REINITIALIZED, CHAMELEON_ERR_NOT_SUPPORTED
      integer CHAMELEON_ERR_ILLEGAL_VALUE, CHAMELEON_ERR_NOT_FOUND
      integer CHAMELEON_ERR_OUT_OF_MEMORY, CHAMELEON_ERR_INTERNAL_LIMIT
      integer CHAMELEON_ERR_UNALLOCATED, CHAMELEON_ERR_FILESYSTEM
      integer CHAMELEON_ERR_UNEXPECTED, CHAMELEON_ERR_SEQUENCE_FLUSHED
      parameter ( CHAMELEON_SUCCESS             =    0 )
      parameter ( CHAMELEON_ERR_NOT_INITIALIZED = -101 )
      parameter ( CHAMELEON_ERR_REINITIALIZED   = -102 )
      parameter ( CHAMELEON_ERR_NOT_SUPPORTED   = -103 )
      parameter ( CHAMELEON_ERR_ILLEGAL_VALUE   = -104 )
      parameter ( CHAMELEON_ERR_NOT_FOUND       = -105 )
      parameter ( CHAMELEON_ERR_OUT_OF_MEMORY   = -106 )
      parameter ( CHAMELEON_ERR_INTERNAL_LIMIT  = -107 )
      parameter ( CHAMELEON_ERR_UNALLOCATED     = -108 )
      parameter ( CHAMELEON_ERR_FILESYSTEM      = -109 )
      parameter ( CHAMELEON_ERR_UNEXPECTED      = -110 )
      parameter ( CHAMELEON_ERR_SEQUENCE_FLUSHED= -111 )

!********************************************************************
!   CHAMELEON constants - kernels options
!
      integer RUNTIME_PRIORITY_MIN, RUNTIME_PRIORITY_MAX
      parameter ( RUNTIME_PRIORITY_MIN = 0     )
      parameter ( RUNTIME_PRIORITY_MAX = 10000 )

!********************************************************************
!   CHAMELEON constants - scheduler properties
!
      integer PRIORITY, CALLBACK, REDUX
      parameter ( PRIORITY = 16 )
      parameter ( CALLBACK = 17 )
      parameter ( REDUX    = 18 )
