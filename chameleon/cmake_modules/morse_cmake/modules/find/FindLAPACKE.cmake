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
# - Find LAPACKE include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(LAPACKE
#               [REQUIRED] # Fail with error if lapacke is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  LAPACKE depends on the following libraries:
#   - LAPACK
#
#  COMPONENTS are optional libraries LAPACKE could be linked with,
#  Use it to drive detection of a specific compilation chain
#  COMPONENTS can be some of the following:
#   - TMG: to check that LAPACKE provides the tmglib interface
#
#  LAPACKE_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found the following variables may be set
#  LAPACKE_PREFIX            - installation path of the lib found
#  <PREFIX>  = LAPACKE
#  <XPREFIX> = <PREFIX>        for common case
#  <XPREFIX> = <PREFIX>_STATIC for static linking
#  <XPREFIX>_FOUND          ... set to 1 if module(s) exist
#  <XPREFIX>_LIBRARIES      ... only the libraries (w/o the '-l')
#  <XPREFIX>_LIBRARY_DIRS   ... the paths of the libraries (w/o the '-L')
#  <XPREFIX>_LDFLAGS        ... all required linker flags
#  <XPREFIX>_LDFLAGS_OTHER  ... all other linker flags
#  <XPREFIX>_INCLUDE_DIRS   ... the '-I' preprocessor flags (w/o the '-I')
#  <XPREFIX>_CFLAGS         ... all required cflags
#  <XPREFIX>_CFLAGS_OTHER   ... the other compiler flags
#
# Set LAPACKE_STATIC to 1 to force using static libraries if exist.
# Set LAPACKE_MT to TRUE  to force using multi-threaded lapack libraries if exists (Intel MKL).
# Set LAPACKE_MT to FALSE to force using sequential lapack libraries if exists (Intel MKL).
# If LAPACKE_MT is undefined, then LAPACKE is linked to the default LAPACK.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::LAPACKE``
#   The headers and libraries to use for LAPACKE, if found.
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DLAPACKE_DIR=path/to/lapacke):
#  LAPACKE_DIR             - Where to find the base directory of lapacke
#  LAPACKE_INCDIR          - Where to find the header files
#  LAPACKE_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: LAPACKE_DIR, LAPACKE_INCDIR, LAPACKE_LIBDIR
#
# LAPACKE could be directly embedded in LAPACK library (ex: Intel MKL) so that
# we test a lapacke function with the lapack libraries found and set LAPACKE
# variables to LAPACK ones if test is successful. To skip this feature and
# look for a stand alone lapacke, please set LAPACKE_STANDALONE to TRUE

#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

# Set variables from environment if needed
# ----------------------------------------
morse_find_package_get_envdir(LAPACKE)

# Duplicate lapack informations into lapacke
# ------------------------------------------
macro(lapacke_init_variables LAPACK)
  # This macro can be called for initialization and/or
  # extension of lapacke discovery

  message( DEBUG "[LAPACKE] ${LAPACK}_LIBRARIES: ${${LAPACK}_LIBRARIES}")
  message( DEBUG "[LAPACKE] ${LAPACK}_LIBRARY_DIRS: ${${LAPACK}_LIBRARY_DIRS}")
  message( DEBUG "[LAPACKE] ${LAPACK}_INCLUDE_DIRS: ${${LAPACK}_INCLUDE_DIRS}")
  message( DEBUG "[LAPACKE] ${LAPACK}_CFLAGS_OTHER: ${${LAPACK}_CFLAGS_OTHER}")
  message( DEBUG "[LAPACKE] ${LAPACK}_LDFLAGS_OTHER: ${${LAPACK}_LDFLAGS_OTHER}")

  if (${LAPACK}_LIBRARIES)
    set(LAPACKE_LAPACK ${LAPACK})
    list(APPEND LAPACKE_LIBRARIES "${${LAPACK}_LIBRARIES}")
  else()
    set(LAPACKE_LAPACK "LAPACKE_LAPACK-NOTFOUND")
    set(LAPACKE_LIBRARIES "LAPACKE_LIBRARIES-NOTFOUND")
  endif()
  if (${LAPACK}_LINKER_FLAGS)
    list(APPEND LAPACKE_LIBRARIES "${${LAPACK}_LINKER_FLAGS}")
  endif()
  if (${LAPACK}_INCLUDE_DIRS)
    list(APPEND LAPACKE_INCLUDE_DIRS "${${LAPACK}_INCLUDE_DIRS}")
  endif()
  if (${LAPACK}_LIBRARY_DIRS)
    list(APPEND LAPACKE_LIBRARY_DIRS "${${LAPACK}_LIBRARY_DIRS}")
  endif()
  if (${LAPACK}_CFLAGS_OTHER)
    list(APPEND LAPACKE_CFLAGS_OTHER "${${LAPACK}_CFLAGS_OTHER}")
  endif()
  if (${LAPACK}_LDFLAGS_OTHER)
    list(APPEND LAPACKE_LDFLAGS_OTHER "${${LAPACK}_LDFLAGS_OTHER}")
  endif()

  morse_cleanup_variables(LAPACKE)

  message( DEBUG "[LAPACKE] LAPACKE_LAPACK: ${LAPACKE_LAPACK}")
  message( DEBUG "[LAPACKE] LAPACKE_LIBRARIES: ${LAPACKE_LIBRARIES}")
  message( DEBUG "[LAPACKE] LAPACKE_LIBRARY_DIRS: ${LAPACKE_LIBRARY_DIRS}")
  message( DEBUG "[LAPACKE] LAPACKE_INCLUDE_DIRS: ${LAPACKE_INCLUDE_DIRS}")
  message( DEBUG "[LAPACKE] LAPACKE_CFLAGS_OTHER: ${LAPACKE_CFLAGS_OTHER}")
  message( DEBUG "[LAPACKE] LAPACKE_LDFLAGS_OTHER: ${LAPACKE_LDFLAGS_OTHER}")

endmacro()

# Check if a lapacke function exists in the lib, and check if the
# advanced complex gemm functions are available
# -------------------------------------------------------------
macro(lapacke_check_library _verbose)
  morse_cmake_required_set(LAPACKE)

  unset(LAPACKE_WORKS CACHE)
  unset(LAPACKE_WITH_LASCL CACHE)
  unset(LAPACKE_dlatms_WORKS CACHE)

  check_function_exists(LAPACKE_dgeqrf LAPACKE_WORKS)
  mark_as_advanced(LAPACKE_WORKS)

  if (LAPACKE_WORKS)
    check_function_exists(LAPACKE_dlascl_work LAPACKE_WITH_LASCL)
    mark_as_advanced(LAPACKE_WITH_LASCL)

    if (LAPACKE_WITH_TMG)
      check_function_exists(LAPACKE_dlatms_work LAPACKE_dlatms_WORKS)
      mark_as_advanced(LAPACKE_dlatms_WORKS)
      if(NOT LAPACKE_dlatms_WORKS)
        if (NOT LAPACKE_FIND_QUIETLY)
          message(STATUS "Looking for lapacke tmg interface : test of LAPACKE_dlatms_work with lapacke and lapack libraries fails")
        endif()
        unset(LAPACKE_WORKS CACHE)
      endif()
    endif()
  endif()

  if(${_verbose})
    if((NOT LAPACKE_WORKS) AND (NOT LAPACKE_FIND_QUIETLY))
      message(STATUS "Looking for lapacke : test of lapacke_dscal with lapacke and lapack libraries fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "CMAKE_REQUIRED_DEFINITIONS: ${CMAKE_REQUIRED_DEFINITIONS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()

  morse_cmake_required_unset()
endmacro()

# Look  for the lapacke header files
# ---------------------------------
macro(lapacke_check_include)

  if ( NOT LAPACKE_INCLUDE_DIRS )

    # Try to find the lapacke header in the given paths
    # -------------------------------------------------
    if (LAPACKE_LIBRARIES MATCHES "mkl")
      set(LAPACKE_hdrs_to_find "mkl.h")
    else()
      set(LAPACKE_hdrs_to_find "lapacke.h")
    endif()

    # call cmake macro to find the header path
    # ----------------------------------------
    morse_find_path(LAPACKE
      HEADERS  ${LAPACKE_hdrs_to_find}
      SUFFIXES "include" "include/lapacke" "include/mkl")

  endif()

endmacro()

# to check LAPACKE components
set(LAPACKE_WITH_TMG OFF)
if( LAPACKE_FIND_COMPONENTS )
  foreach( component ${LAPACKE_FIND_COMPONENTS} )
    if (${component} STREQUAL "TMG")
      set(LAPACKE_WITH_TMG ON)
    endif()
  endforeach()
endif()

if (DEFINED LAPACKE_MT )
  set( LAPACKE_MT_internal ${LAPACKE_MT} )
endif()

# LAPACKE may depend on TMG (which depends on lapack), try to find it
if (LAPACKE_WITH_TMG)
  if (DEFINED LAPACKE_MT)
    set(TMG_MT ${LAPACKE_MT})
  endif()
  if(LAPACKE_FIND_REQUIRED)
    find_package(TMG REQUIRED)
  else()
    find_package(TMG)
  endif()

  # Lapacke depends only on TMG, which itself depends on the correct
  # lapack. No need to look for the MT/SEQ one anymore
  set( LAPACKE_dep TMG )
  unset( LAPACKE_MT_internal )
else()

  # LAPACKE depends on LAPACK, try to find it
  if(LAPACKE_FIND_REQUIRED)
    find_package(LAPACKEXT QUIET REQUIRED)
  else()
    find_package(LAPACKEXT QUIET)
  endif()

  set( LAPACKE_dep LAPACK )
endif()

if(DEFINED LAPACKE_MT_internal)
  if (LAPACKE_MT_internal)
    lapacke_init_variables("${LAPACKE_dep}_MT")
  else()
    lapacke_init_variables("${LAPACKE_dep}_SEQ")
  endif()
else()
  lapacke_init_variables("${LAPACKE_dep}")
endif()

# LAPACKE depends on LAPACK
if (LAPACKE_LAPACK)

  if (NOT LAPACKE_STANDALONE)

    # Check if the lapack library includes cblas
    lapacke_check_library(FALSE)

    # Lapack lib includes lapacke
    if(LAPACKE_WORKS)
      if(NOT LAPACKE_FIND_QUIETLY)
        message(STATUS "Looking for lapacke: test with lapack succeeded")
      endif()

      # Set the mkl library dirs for compatibility with former version
      # --------------------------------------------------------------
      if (LAPACKE_LIBRARIES MATCHES "mkl" AND DEFINED ENV{MKLROOT})
        set(LAPACKE_PREFIX "$ENV{MKLROOT}" CACHE PATH "Installation directory of LAPACKE library" FORCE)
        set(LAPACKE_LIBRARY_DIRS "${LAPACKE_PREFIX}/lib/intel64")
      endif()

      lapacke_check_include()

    endif()

  endif (NOT LAPACKE_STANDALONE)

  # test fails with lapack: try to find LAPACKE lib exterior to LAPACK
  if (LAPACKE_STANDALONE OR (NOT LAPACKE_WORKS))

    if((NOT LAPACKE_WORKS) AND (NOT LAPACKE_FIND_QUIETLY))
      message(STATUS "Looking for lapacke : test with lapack fails")
    endif()

    # Make sure LAPACKE is invalidated
    unset(LAPACKE_WORKS CACHE)
    unset(LAPACKE_LIBRARIES)
    unset(LAPACKE_INCLUDE_DIRS)
    unset(LAPACKE_LIBRARY_DIRS)
    unset(LAPACKE_CFLAGS_OTHER)
    unset(LAPACKE_LDFLAGS_OTHER)

    # try with pkg-config
    set(ENV_MKLROOT "$ENV{MKLROOT}")
    set(LAPACKE_GIVEN_BY_USER FALSE)
    if ( LAPACKE_DIR OR ( LAPACKE_INCDIR AND LAPACKE_LIBDIR ) OR ( ENV_MKLROOT ) )
      set(LAPACKE_GIVEN_BY_USER TRUE)
    endif()

    # Search for lapacke with pkg-config
    # --------------------------------
    find_package(PkgConfig QUIET)
    if( PKG_CONFIG_EXECUTABLE AND (NOT LAPACKE_GIVEN_BY_USER))

      if (BLA_STATIC)
        set(MKL_STR_BLA_STATIC "static")
      else()
        set(MKL_STR_BLA_STATIC "dynamic")
      endif()
      # try different blas
      if (BLA_VENDOR STREQUAL "Intel10_64lp")
        pkg_search_module(LAPACKE mkl-${MKL_STR_BLA_STATIC}-lp64-iomp)
      elseif(BLA_VENDOR STREQUAL "Intel10_64lp_seq")
        pkg_search_module(LAPACKE mkl-${MKL_STR_BLA_STATIC}-lp64-seq)
      elseif(BLA_VENDOR STREQUAL "Open")
        pkg_search_module(LAPACKE openblas)
      elseif(BLA_VENDOR STREQUAL "Generic")
        pkg_search_module(LAPACKE lapacke)
      else()
        pkg_search_module(LAPACKE lapacke)
        pkg_search_module(LAPACKE openblas)
        pkg_search_module(LAPACKE mkl-${MKL_STR_BLA_STATIC}-lp64-seq)
      endif()

      if (NOT LAPACKE_FIND_QUIETLY)
        if (LAPACKE_FOUND AND LAPACKE_LIBRARIES)
          message(STATUS "Looking for LAPACKE - found using PkgConfig")
        else()
          message(STATUS "${Magenta}Looking for LAPACKE - not found using PkgConfig."
            "\n   Perhaps you should add the directory containing lapacke.pc to"
            "\n   the PKG_CONFIG_PATH environment variable.${ColourReset}")
        endif()
      endif()

      if (LAPACKE_FOUND AND LAPACKE_LIBRARIES)
        set(LAPACKE_FOUND_WITH_PKGCONFIG "TRUE")
        morse_find_pkgconfig_libraries_absolute_path(LAPACKE)
      else()
        set(LAPACKE_FOUND_WITH_PKGCONFIG "FALSE")
      endif()

      if (LAPACKE_STATIC AND LAPACKE_STATIC_LIBRARIES)
        set(LAPACKE_DEPENDENCIES ${LAPACKE_STATIC_LIBRARIES})
        list(REMOVE_ITEM LAPACKE_DEPENDENCIES "lapacke")
        list(APPEND LAPACKE_LIBRARIES ${LAPACKE_DEPENDENCIES})
        set(LAPACKE_CFLAGS_OTHER ${LAPACKE_STATIC_CFLAGS_OTHER})
        set(LAPACKE_LDFLAGS_OTHER ${LAPACKE_STATIC_LDFLAGS_OTHER})
        if (NOT LAPACKE_FIND_QUIETLY)
          message(STATUS "LAPACKE_STATIC set to 1 by user, LAPACKE_LIBRARIES: ${LAPACKE_LIBRARIES}.")
        endif()
      endif()
    endif()

    # Search for lapacke the classical way
    # ------------------------------------
    if (NOT LAPACKE_FOUND_WITH_PKGCONFIG OR LAPACKE_GIVEN_BY_USER)
      # Looking for include
      # -------------------
      lapacke_check_include()

      # Looking for lib
      # ---------------
      morse_find_library(LAPACKE
        LIBRARIES lapacke
        SUFFIXES  lib lib32 lib64)

    endif (NOT LAPACKE_FOUND_WITH_PKGCONFIG OR LAPACKE_GIVEN_BY_USER)
  endif (LAPACKE_STANDALONE OR (NOT LAPACKE_WORKS))

  # check if static or dynamic lib
  morse_check_static_or_dynamic(LAPACKE LAPACKE_LIBRARIES)
  if(LAPACKE_STATIC)
    set(STATIC "_STATIC")
  endif()

  # We found a lapacke library with pkg-config or manually
  # ------------------------------------------------------
  if((NOT LAPACKE_WORKS) AND LAPACKE_LIBRARIES)

    # Need to add dependencies if not found with pkg-config
    # -----------------------------------------------------
    if (LAPACKE_STATIC OR CBLAS_STATIC OR BLA_STATIC)
      if (NOT LAPACKE_FOUND_WITH_PKGCONFIG)
        # Extend the discovered library with the lapack ones
        lapacke_init_variables(${LAPACKE_LAPACK})
      else()
        # save link with dependencies
        set(LAPACKE_LIBRARIES "${LAPACKE_STATIC_LIBRARIES}")
      endif()
    endif()

    # Test link
    lapacke_check_library(TRUE)

    list(GET LAPACKE_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" DIRECTORY)
    if (NOT LAPACKE_LIBRARY_DIRS)
      set(LAPACKE_LIBRARY_DIRS "${first_lib_path}")
    endif()
    if (NOT LAPACKE_PREFIX)
      if (${first_lib_path} MATCHES "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)")
        string(REGEX REPLACE "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)" "" not_cached_dir "${first_lib_path}")
        set(LAPACKE_PREFIX "${not_cached_dir}" CACHE PATH "Installation directory of LAPACKE library" FORCE)
      else()
        set(LAPACKE_PREFIX "${first_lib_path}" CACHE PATH "Installation directory of LAPACKE library" FORCE)
      endif()
    endif (NOT LAPACKE_PREFIX)
    if (NOT LAPACKE_INCLUDE_DIRS)
      if (EXISTS "${LAPACKE_PREFIX}/include")
        set(LAPACKE_INCLUDE_DIRS "${LAPACKE_PREFIX}/include")
      endif()
    endif()
    mark_as_advanced(LAPACKE_PREFIX)

  endif()

else(LAPACKE_LAPACK)

  if (NOT LAPACKE_FIND_QUIETLY)
    message(STATUS "LAPACKE requires LAPACK but LAPACK has not been found. "
      "Please look for LAPACK first or change the MT support required (LAPACKE_MT=${LAPACKE_MT}).")
  endif()

endif(LAPACKE_LAPACK)

# check that LAPACKE has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACKE DEFAULT_MSG
  LAPACKE_LIBRARIES
  LAPACKE_WORKS)

# Add imported target
if (LAPACKE_FOUND)
  morse_create_imported_target(LAPACKE)
endif()
