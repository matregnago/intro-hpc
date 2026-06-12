###
#
# @copyright (c) 2012-2020 Inria. All rights reserved.
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
# - Extension of FindLAPACK in order to provide both sequential and multi-threaded libraries when possible (e.g. Intel MKL)
#
# Result variables:
#
# - ``LAPACKEXT_FOUND`` if LAPACK is found
# - ``LAPACK_SEQ_LINKER_FLAGS`` sequential version of LAPACK_LINKER_FLAGS
# - ``LAPACK_MT_LINKER_FLAGS`` multi-threaded version of LAPACK_LINKER_FLAGS
# - ``LAPACK_SEQ_LIBRARIES`` sequential version of LAPACK_LIBRARIES
# - ``LAPACK_MT_LIBRARIES`` multi-threaded version of LAPACK_LIBRARIES
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``LAPACK::LAPACK_SEQ``
#   The libraries of sequential LAPACK, if found.
#
# ``LAPACK::LAPACK_MT``
#   The libraries of multi-threaded LAPACK, if found.
#
#=============================================================================

if(NOT LAPACKEXT_FIND_QUIETLY)
  message(STATUS "FindLAPACKEXT: Try to find LAPACK")
endif()

# LAPACKEXT first search LAPACK available
if(LAPACKEXT_FIND_REQUIRED)
  find_package(LAPACK QUIET REQUIRED)
else()
  find_package(LAPACK QUIET)
endif()

set(BLA_VENDOR_COPY ${BLA_VENDOR})

macro(lapackext_set_library VERSION)
  if (LAPACK_FOUND)
    set(LAPACK_${VERSION}_FOUND ${LAPACK_FOUND})
    if(NOT TARGET LAPACK::LAPACK_${VERSION})
      add_library(LAPACK::LAPACK_${VERSION} INTERFACE IMPORTED)
    endif()
    if (LAPACK_LINKER_FLAGS)
      set(LAPACK_${VERSION}_LINKER_FLAGS ${LAPACK_LINKER_FLAGS})
      set_target_properties(LAPACK::LAPACK_${VERSION} PROPERTIES
        INTERFACE_LINK_OPTIONS "${LAPACK_LINKER_FLAGS}"
        )
    endif()
    if (LAPACK_LIBRARIES)
      if(NOT LAPACKEXT_FIND_QUIETLY)
        message(STATUS "FindLAPACKEXT: Found LAPACK ${BLA_VENDOR}")
        message(STATUS "FindLAPACKEXT: Store following libraries in LAPACK_${VERSION}_LIBRARIES and target LAPACK::LAPACK_${VERSION} ${LAPACK_LIBRARIES}")
      endif()
      set(LAPACK_${VERSION}_LIBRARIES ${LAPACK_LIBRARIES})
      set_target_properties(LAPACK::LAPACK_${VERSION} PROPERTIES
        INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}"
        )
    endif()
  endif()
endmacro()

if (LAPACK_FOUND)

  if(LAPACK_LIBRARIES MATCHES "mkl")
    if(NOT LAPACKEXT_FIND_QUIETLY)
      message(STATUS "FindLAPACKEXT: LAPACK_LIBRARIES matches mkl")
    endif()

    # Look for the sequential MKL
    # ---------------------------
    message( DEBUG "FindLAPACKEXT: Search for sequential MKL")
    set(BLA_VENDOR "Intel10_64lp_seq")
    unset(LAPACK_FOUND)
    unset(LAPACK_LINKER_FLAGS)
    unset(LAPACK_LIBRARIES)
    find_package(LAPACK QUIET)
    lapackext_set_library( SEQ )

    # Look for the multi-threaded MKL
    # -------------------------------
    message( DEBUG "FindLAPACKEXT: Search for parallel MKL")
    set(BLA_VENDOR "Intel10_64lp")
    unset(LAPACK_FOUND)
    unset(LAPACK_LINKER_FLAGS)
    unset(LAPACK_LIBRARIES)
    find_package(LAPACK QUIET)
    lapackext_set_library( MT )

    # Restore the original library to make sure LAPACK_LIBRARIES is set
    # -----------------------------------------------------------------
    message( DEBUG "FindLAPACKEXT: Reset to default search")
    set(BLA_VENDOR ${BLA_VENDOR_COPY})
    unset(LAPACK_FOUND)
    unset(LAPACK_LINKER_FLAGS)
    unset(LAPACK_LIBRARIES)
    find_package(LAPACK QUIET)

  else(LAPACK_LIBRARIES MATCHES "mkl")

    lapackext_set_library( SEQ )

  endif(LAPACK_LIBRARIES MATCHES "mkl")

else(LAPACK_FOUND)
  if(NOT LAPACKEXT_FIND_QUIETLY)
    message(STATUS "FindLAPACKEXT: LAPACK not found or LAPACK_LIBRARIES does not match mkl")
  endif()
endif(LAPACK_FOUND)

set(BLA_VENDOR ${BLA_VENDOR_COPY})

# check that LAPACKEXT has been found
# -----------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACKEXT DEFAULT_MSG
  LAPACK_LIBRARIES)
