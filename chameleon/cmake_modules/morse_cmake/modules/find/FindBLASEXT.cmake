###
#
# @copyright (c) 2012-2024 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
# Copyright 2012-2024 Mathieu Faverge
# Copyright 2013-2024 Florent Pruvost
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file MORSE-Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of Morse, substitute the full
#  License text for the above reference.)
###
#
# - Extension of FindBLAS in order to provide both sequential and multi-threaded libraries when possible (e.g. Intel MKL)
#
# Result variables:
#
# - ``BLASEXT_FOUND`` if BLAS is found
# - ``BLAS_SEQ_LINKER_FLAGS`` sequential version of BLAS_LINKER_FLAGS
# - ``BLAS_MT_LINKER_FLAGS`` multi-threaded version of BLAS_LINKER_FLAGS
# - ``BLAS_SEQ_LIBRARIES`` sequential version of BLAS_LIBRARIES
# - ``BLAS_MT_LIBRARIES`` multi-threaded version of BLAS_LIBRARIES
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``BLAS::BLAS_SEQ``
#   The libraries of sequential blas, if found.
#
# ``BLAS::BLAS_MT``
#   The libraries of multi-threaded blas, if found.
#
#=============================================================================

if(NOT BLASEXT_FIND_QUIETLY)
  message(STATUS "FindBLASEXT: Try to find BLAS")
endif()

# BLASEXT first search for the main BLAS library available
if(BLASEXT_FIND_REQUIRED)
  find_package(BLAS QUIET REQUIRED)
else()
  find_package(BLAS QUIET)
endif()

set(BLA_VENDOR_COPY ${BLA_VENDOR})

macro(blasext_set_library VERSION)
  if (BLAS_FOUND)
    set(BLAS_${VERSION}_FOUND ${BLAS_FOUND})
    if(NOT TARGET BLAS::BLAS_${VERSION})
      add_library(BLAS::BLAS_${VERSION} INTERFACE IMPORTED)
    endif()
    if (BLAS_LINKER_FLAGS)
      set(BLAS_${VERSION}_LINKER_FLAGS ${BLAS_LINKER_FLAGS})
      set_target_properties(BLAS::BLAS_${VERSION} PROPERTIES
        INTERFACE_LINK_OPTIONS "${BLAS_LINKER_FLAGS}"
        )
    endif()
    if (BLAS_LIBRARIES)
      if(NOT BLASEXT_FIND_QUIETLY)
        message(STATUS "FindBLASEXT: Found BLAS ${BLA_VENDOR}")
        message(STATUS "FindBLASEXT: Store following libraries in BLAS_${VERSION}_LIBRARIES and target BLAS::BLAS_${VERSION} ${BLAS_LIBRARIES}")
      endif()
      set(BLAS_${VERSION}_LIBRARIES ${BLAS_LIBRARIES})
      set_target_properties(BLAS::BLAS_${VERSION} PROPERTIES
        INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}"
        )
    endif()
  endif()
endmacro()

if (BLAS_FOUND)

  if(BLAS_LIBRARIES MATCHES "mkl")
    if(NOT BLASEXT_FIND_QUIETLY)
      message(STATUS "FindBLASEXT: BLAS_LIBRARIES matches mkl")
    endif()

    # Look for the sequential MKL
    # ---------------------------
    set(BLA_VENDOR "Intel10_64lp_seq")
    unset(BLAS_FOUND)
    unset(BLAS_LINKER_FLAGS)
    unset(BLAS_LIBRARIES)
    find_package(BLAS QUIET)
    blasext_set_library( SEQ )

    # Look for the multi-threaded MKL
    # -------------------------------
    set(BLA_VENDOR "Intel10_64lp")
    unset(BLAS_FOUND)
    unset(BLAS_LINKER_FLAGS)
    unset(BLAS_LIBRARIES)
    find_package(BLAS QUIET)
    blasext_set_library( MT )

    # Restore the original library to make sure BLAS_LIBRARIES is set
    # ---------------------------------------------------------------
    set(BLA_VENDOR ${BLA_VENDOR_COPY})
    unset(BLAS_FOUND)
    unset(BLAS_LINKER_FLAGS)
    unset(BLAS_LIBRARIES)
    find_package(BLAS QUIET)

  else(BLAS_LIBRARIES MATCHES "mkl")

    blasext_set_library( SEQ )

  endif(BLAS_LIBRARIES MATCHES "mkl")

  # number of theads function detection
  if(BLAS_LIBRARIES MATCHES "blis|mkl|openblas")
    enable_language(Fortran)
    include(CheckFortranFunctionExists)
    set(CMAKE_REQUIRED_LIBRARIES "${BLAS_LIBRARIES}")
    if(BLAS_LIBRARIES MATCHES "blis")
      unset(HAVE_BLI_THREAD_SET_NUM_THREADS CACHE)
      check_fortran_function_exists(bli_thread_set_num_threads HAVE_BLI_THREAD_SET_NUM_THREADS)
    elseif(BLAS_LIBRARIES MATCHES "mkl")
      unset(HAVE_MKL_SET_NUM_THREADS CACHE)
      check_fortran_function_exists(mkl_set_num_threads HAVE_MKL_SET_NUM_THREADS)
    elseif(BLAS_LIBRARIES MATCHES "openblas")
      unset(HAVE_OPENBLAS_SET_NUM_THREADS CACHE)
      check_fortran_function_exists(openblas_set_num_threads HAVE_OPENBLAS_SET_NUM_THREADS)
    endif()
    set(CMAKE_REQUIRED_LIBRARIES)
    if (HAVE_BLI_THREAD_SET_NUM_THREADS OR HAVE_MKL_SET_NUM_THREADS OR HAVE_OPENBLAS_SET_NUM_THREADS)
      set(HAVE_BLAS_SET_NUM_THREADS 1)
    endif()
  endif()

else(BLAS_FOUND)
  if(NOT BLASEXT_FIND_QUIETLY)
    message(STATUS "FindBLASEXT: BLAS not found or BLAS_LIBRARIES does not match mkl")
  endif()
endif(BLAS_FOUND)

set(BLA_VENDOR ${BLA_VENDOR_COPY})

# check that BLASEXT has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BLASEXT DEFAULT_MSG
  BLAS_LIBRARIES)
