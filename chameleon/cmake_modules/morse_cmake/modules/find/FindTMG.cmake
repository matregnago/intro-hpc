###
#
# @copyright (c) 2012-2020 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013-2020 Florent Pruvost
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
# - Find TMG include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(TMG
#               [REQUIRED]             # Fail with error if tmg is not found
#              )
#
# This module finds headers and tmg library.
# Results are reported in variables:
#  TMG_FOUND             - True if headers and requested libraries were found
#  TMG_PREFIX            - installation path of the lib found
#  TMG_CFLAGS_OTHER      - tmglib compiler flags without headers paths
#  TMG_LDFLAGS_OTHER     - tmglib linker flags without libraries
#  TMG_INCLUDE_DIRS      - tmglib include directories
#  TMG_LIBRARY_DIRS      - tmglib link directories
#  TMG_LIBRARIES         - tmglib libraries to be linked (absolute path)
#
# Set TMG_STATIC to 1 to force using static libraries if exist.
# Set TMG_MT to TRUE  to force using multi-threaded lapack libraries if exists (Intel MKL).
# Set TMG_MT to FALSE to force using sequential lapack libraries if exists (Intel MKL).
# If TMG_MT is undefined, then TMG is linked to the default LAPACK.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::TMG``
#   The headers and libraries to use for TMG, if found.
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DTMG_DIR=path/to/tmg):
#  TMG_DIR              - Where to find the base directory of tmg
#  TMG_INCDIR           - Where to find the header files
#  TMG_LIBDIR           - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: TMG_DIR, TMG_INCDIR, TMG_LIBDIR
#
#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

# Set variables from environment if needed
# ----------------------------------------
morse_find_package_get_envdir(TMG)

# used to test a TMG function after
get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
if (NOT _LANGUAGES_ MATCHES Fortran)
  include(CheckFunctionExists)
else (NOT _LANGUAGES_ MATCHES Fortran)
  include(CheckFortranFunctionExists)
endif (NOT _LANGUAGES_ MATCHES Fortran)

# Duplicate lapack informations into tmg
# --------------------------------------
macro(tmg_init_variables LAPACK)
  # This macro can be called for initialization and/or
  # extension of tmg discovery

  message( DEBUG "[TMG] ${LAPACK}_LIBRARIES: ${${LAPACK}_LIBRARIES}")
  message( DEBUG "[TMG] ${LAPACK}_LIBRARY_DIRS: ${${LAPACK}_LIBRARY_DIRS}")
  message( DEBUG "[TMG] ${LAPACK}_INCLUDE_DIRS: ${${LAPACK}_INCLUDE_DIRS}")
  message( DEBUG "[TMG] ${LAPACK}_CFLAGS_OTHER: ${${LAPACK}_CFLAGS_OTHER}")
  message( DEBUG "[TMG] ${LAPACK}_LDFLAGS_OTHER: ${${LAPACK}_LDFLAGS_OTHER}")

  if (${LAPACK}_LIBRARIES)
    set(TMG_LAPACK ${LAPACK})
    list(APPEND TMG_LIBRARIES "${${LAPACK}_LIBRARIES}")
  else()
    set(TMG_LAPACK "TMG_LAPACK-NOTFOUND")
    set(TMG_LIBRARIES "TMG_LIBRARIES-NOTFOUND")
  endif()
  if (${LAPACK}_LINKER_FLAGS)
    list(APPEND TMG_LIBRARIES "${${LAPACK}_LINKER_FLAGS}")
  endif()
  if (${LAPACK}_INCLUDE_DIRS)
    list(APPEND TMG_INCLUDE_DIRS "${${LAPACK}_INCLUDE_DIRS}")
  endif()
  if (${LAPACK}_LIBRARY_DIRS)
    list(APPEND TMG_LIBRARY_DIRS "${${LAPACK}_LIBRARY_DIRS}")
  endif()
  if (${LAPACK}_CFLAGS_OTHER)
    list(APPEND TMG_CFLAGS_OTHER "${${LAPACK}_CFLAGS_OTHER}")
  endif()
  if (${LAPACK}_LDFLAGS_OTHER)
    list(APPEND TMG_LDFLAGS_OTHER "${${LAPACK}_LDFLAGS_OTHER}")
  endif()

  morse_cleanup_variables(TMG)

  message( DEBUG "[TMG] TMG_LAPACK: ${TMG_LAPACK}")
  message( DEBUG "[TMG] TMG_LIBRARIES: ${TMG_LIBRARIES}")
  message( DEBUG "[TMG] TMG_LIBRARY_DIRS: ${TMG_LIBRARY_DIRS}")
  message( DEBUG "[TMG] TMG_INCLUDE_DIRS: ${TMG_INCLUDE_DIRS}")
  message( DEBUG "[TMG] TMG_CFLAGS_OTHER: ${TMG_CFLAGS_OTHER}")
  message( DEBUG "[TMG] TMG_LDFLAGS_OTHER: ${TMG_LDFLAGS_OTHER}")

endmacro()

# Check if tmg functions exist in the lib
# ---------------------------------------
macro(tmg_check_library _verbose)
  morse_cmake_required_set(TMG)

  unset(TMG_WORKS CACHE)
  unset(TMG_dlatms_WORKS CACHE)
  unset(TMG_dlagsy_WORKS CACHE)

  if (NOT _LANGUAGES_ MATCHES Fortran)
    check_function_exists(dlatms TMG_dlatms_WORKS)
  else (NOT _LANGUAGES_ MATCHES Fortran)
    check_fortran_function_exists(dlatms TMG_dlatms_WORKS)
  endif (NOT _LANGUAGES_ MATCHES Fortran)

  if (NOT _LANGUAGES_ MATCHES Fortran)
    check_function_exists(dlagsy TMG_dlagsy_WORKS)
  else (NOT _LANGUAGES_ MATCHES Fortran)
    check_fortran_function_exists(dlagsy TMG_dlagsy_WORKS)
  endif (NOT _LANGUAGES_ MATCHES Fortran)

  if ( TMG_dlatms_WORKS AND TMG_dlagsy_WORKS )
    set(TMG_WORKS TRUE)
    mark_as_advanced(TMG_WORKS)
  endif()

  if(${_verbose})
    if((NOT TMG_WORKS) AND (NOT TMG_FIND_QUIETLY))
      message(STATUS "Looking for tmg: test of dlatms and dlagsy with tmg and lapack libraries fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "CMAKE_REQUIRED_DEFINITIONS: ${CMAKE_REQUIRED_DEFINITIONS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()

  morse_cmake_required_unset()
endmacro()

# TMG depends on LAPACK anyway, try to find it
if(TMG_FIND_REQUIRED)
  find_package(LAPACKEXT QUIET REQUIRED)
else()
  find_package(LAPACKEXT QUIET)
endif()

if(DEFINED TMG_MT)
  if (TMG_MT)
    tmg_init_variables("LAPACK_MT")
  else()
    tmg_init_variables("LAPACK_SEQ")
  endif()
else()
  tmg_init_variables("LAPACK")
endif()

# TMG depends on LAPACK
if (TMG_LAPACK)

  # Test link with lapack only
  tmg_check_library(FALSE)

  if(TMG_WORKS)

    if(NOT TMG_FIND_QUIETLY)
      message(STATUS "Looking for tmg: test with lapack succeeds")
    endif()

  else()

    if(NOT TMG_FIND_QUIETLY)
      message(STATUS "Looking for tmg : test with lapack fails")
      message(STATUS "Looking for tmg : try to find it elsewhere")
    endif()

    # Make sure TMG is invalidated
    unset(TMG_WORKS CACHE)
    unset(TMG_LIBRARIES)
    unset(TMG_INCLUDE_DIRS)
    unset(TMG_LIBRARY_DIRS)
    unset(TMG_CFLAGS_OTHER)
    unset(TMG_LDFLAGS_OTHER)

    # No include, let's set the the include_dir
    set(TMG_INCLUDE_DIRS "")

    # Looking for lib tmg
    # -------------------
    morse_find_library(TMG
      LIBRARIES tmglib tmg
      SUFFIXES  lib lib32 lib64
      OPTIONAL )

    # If found, add path to cmake variable
    # ------------------------------------
    if ((NOT TMG_tmglib_LIBRARY) AND (NOT TMG_tmg_LIBRARY))
      set(TMG_LIBRARIES    "TMG_LIBRARIES-NOTFOUND")
      set(TMG_LIBRARY_DIRS "TMG_LIBRARY_DIRS-NOTFOUND")
      if(NOT TMG_FIND_QUIETLY)
        message(STATUS "Looking for tmg -- lib tmg not found")
      endif()
    endif ()

    # We found a tmg library, let's check
    # -----------------------------------
    if(TMG_LIBRARIES)

      # check if static or dynamic lib
      morse_check_static_or_dynamic(TMG TMG_LIBRARIES)
      if(TMG_STATIC)
        set(STATIC "_STATIC")
      endif()

      # Extend the tmg variables with the lapack ones
      if (LAPACKE_STATIC OR CBLAS_STATIC OR BLA_STATIC)
        tmg_init_variables(${TMG_LAPACK})
      endif()

      # test link
      tmg_check_library(TRUE)
    endif()
  endif()

  # Additional settings for all cases (in or out lapack lib)
  # --------------------------------------------------------
  if ( TMG_LIBRARIES )

    list(GET TMG_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" DIRECTORY)
    if (NOT TMG_LIBRARY_DIRS)
      set(TMG_LIBRARY_DIRS "${first_lib_path}")
    endif()
    if (${first_lib_path} MATCHES "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)")
      string(REGEX REPLACE "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)" "" not_cached_dir "${first_lib_path}")
      set(TMG_PREFIX "${not_cached_dir}" CACHE PATH "Installation directory of TMG library" FORCE)
    else()
      set(TMG_PREFIX "${first_lib_path}" CACHE PATH "Installation directory of TMG library" FORCE)
    endif()
    mark_as_advanced(TMG_PREFIX)

  endif(TMG_LIBRARIES)

else()

  if(NOT TMG_FIND_QUIETLY)
    message(STATUS "TMG requires LAPACK but LAPACK has not been found. "
      "Please look for LAPACK first or change the MT support required (TMG_MT=${TMG_MT}).")
  endif()

endif()

# check that TMG has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TMG DEFAULT_MSG
  TMG_LIBRARIES
  TMG_WORKS)

# Add imported target
if (TMG_FOUND)
  morse_create_imported_target(TMG)
endif()
