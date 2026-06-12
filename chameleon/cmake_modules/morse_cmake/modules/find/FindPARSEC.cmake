###
#
# @copyright (c) 2012-2020 Inria. All rights reserved.
# @copyright (c) 2012-2015 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
# Copyright 2012-2019 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013-2018 Florent Pruvost
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
#
###
#
# - Find PARSEC include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(PARSEC
#               [version] [EXACT]      # Minimum or EXACT version e.g. 1.1
#               [REQUIRED]             # Fail with error if parsec is not found
#              )
#
#  PARSEC depends on the following libraries:
#   - Threads, m, rt
#  Optional dependencies
#   - AYUDAME: ??
#   - CUDA
#   - HWLOC
#   - MPI
#
#  PARSEC_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = PARSEC
#  <PREFIX>_PREFIX          ... installation path of the lib found
#  <PREFIX>_VERSION         ... version of the lib
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
# Set PARSEC_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::PARSEC``
#   The headers and libraries to use for PARSEC if found.
#
#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

include(CheckSymbolExists)

if (NOT PARSEC_FIND_QUIETLY)
  message(STATUS "FindPARSEC needs pkg-config program and PKG_CONFIG_PATH set with parsec.pc file path.")
endif()

# use pkg-config to detect include/library dirs (if pkg-config is available)
# --------------------------------------------------------------------------
if(PKG_CONFIG_EXECUTABLE)
  unset(PARSEC_FOUND CACHE)
  pkg_search_module(PARSEC parsec)

  if (NOT PARSEC_FIND_QUIETLY)
    if (PARSEC_FOUND AND PARSEC_LIBRARIES)
      message(STATUS "Looking for PARSEC - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for PARSEC - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing libparsec.pc"
        "\n   to the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()

  # call cmake macro to find parsec bin programs
  if (PARSEC_PREFIX)
    # create list of binaries to find
    set(PARSEC_bins_to_find "parsec-ptgpp;parsec_ptgpp")
    foreach(parsec_bin ${PARSEC_bins_to_find})
      set(PARSEC_${parsec_bin} "PARSEC_${parsec_bin}-NOTFOUND")
      find_program(PARSEC_${parsec_bin}
        NAMES ${parsec_bin}
        HINTS ${PARSEC_PREFIX}
        PATH_SUFFIXES "bin"
        NO_PACKAGE_ROOT_PATH NO_CMAKE_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_FIND_ROOT_PATH)
    endforeach()

    if ( NOT PARSEC_parsec_ptgpp )
      set( PARSEC_parsec_ptgpp ${PARSEC_parsec-ptgpp} )
    endif()
  else()
    if (PARSEC_FIND_REQUIRED)
      message(FATAL_ERROR "PARSEC_PREFIX not defined by pkg_search_module")
    endif()
  endif()

  if (PARSEC_FOUND AND PARSEC_LIBRARIES)
    if (NOT PARSEC_INCLUDE_DIRS)
      pkg_get_variable(PARSEC_INCLUDE_DIRS parsec includedir)
    endif()
    set(PARSEC_FOUND_WITH_PKGCONFIG "TRUE")
    morse_find_pkgconfig_libraries_absolute_path(PARSEC)
  else()
    set(PARSEC_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

  if (PARSEC_STATIC AND PARSEC_STATIC_LIBRARIES)
    set (PARSEC_DEPENDENCIES ${PARSEC_STATIC_LIBRARIES})
    list (REMOVE_ITEM PARSEC_DEPENDENCIES "parsec")
    list (APPEND PARSEC_LIBRARIES ${PARSEC_DEPENDENCIES})
    set(PARSEC_CFLAGS_OTHER ${PARSEC_STATIC_CFLAGS_OTHER})
    set(PARSEC_LDFLAGS_OTHER ${PARSEC_STATIC_LDFLAGS_OTHER})
    if (NOT PARSEC_FIND_QUIETLY)
      message(STATUS "PARSEC_STATIC set to 1 by user, PARSEC_LIBRARIES: ${PARSEC_LIBRARIES}.")
    endif()
  endif()
endif()

# check a function to validate the find
if(PARSEC_FOUND AND PARSEC_LIBRARIES)

  # check if static or dynamic lib
  morse_check_static_or_dynamic(PARSEC PARSEC_LIBRARIES)
  if(PARSEC_STATIC)
    set(STATIC "_STATIC")
  else()
    set(STATIC "")
  endif()

  # set required libraries for link
  morse_set_required_test_lib_link(PARSEC)

  # test link
  unset(PARSEC_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(parsec_init PARSEC_WORKS)
  mark_as_advanced(PARSEC_WORKS)

  if(NOT PARSEC_WORKS)
    if(NOT PARSEC_FIND_QUIETLY)
      message(STATUS "Looking for parsec : test of parsec_init fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
      message(STATUS "Maybe PARSEC is linked with specific libraries. "
        "Have you tried with COMPONENTS (HWLOC, CUDA, MPI)? "
        "See the explanation in FindPARSEC.cmake.")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

endif(PARSEC_FOUND AND PARSEC_LIBRARIES)

# check that PARSEC has been found
# --------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PARSEC DEFAULT_MSG
  PARSEC_LIBRARIES
  PARSEC_parsec_ptgpp
  PARSEC_WORKS)

# Add imported targe
if (PARSEC_FOUND)
  morse_create_imported_target(PARSEC)
endif()
