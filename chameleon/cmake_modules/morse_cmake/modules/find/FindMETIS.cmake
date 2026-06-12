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
#
###
#
# - Find METIS include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(METIS
#               [REQUIRED]             # Fail with error if metis is not found
#              )
#
# This module finds headers and metis library.
# Results are reported in variables:
#  METIS_FOUND           - True if headers and requested libraries were found
#  METIS_PREFIX          - installation path of the lib found
#  METIS_INCLUDE_DIRS    - metis include directories
#  METIS_LIBRARY_DIRS    - Link directories for metis libraries
#  METIS_LIBRARIES       - metis component libraries to be linked
#  METIX_INTSIZE         - Number of octets occupied by a idx_t (IDXTYPEWIDTH)
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DMETIS_DIR=path/to/metis):
#  METIS_DIR             - Where to find the base directory of metis
#  METIS_INCDIR          - Where to find the header files
#  METIS_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: METIS_DIR, METIS_INCDIR, METIS_LIBDIR
#
# Set METIS_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::METIS``
#   The headers and libraries to use for METIS, if found.
#
#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

# Set variables from environment if needed
# ----------------------------------------
morse_find_package_get_envdir(METIS)

# Looking for include
# -------------------
morse_find_path(METIS
  HEADERS  metis.h
  SUFFIXES include include/metis)

# Looking for lib
# ---------------
morse_find_library(METIS
  LIBRARIES metis
  SUFFIXES lib lib32 lib64)

# check a function to validate the find
if(METIS_LIBRARIES)

  # Metis may depend on the M library and the static compilation does
  # not include the dependency, so we enforce it
  find_package(M QUIET)
  if ( M_FOUND )
    set( METIS_LIBRARIES "${METIS_LIBRARIES};${M_LIBRARIES}" )
  endif()

  morse_cmake_required_set(METIS)

  # test link
  unset(METIS_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(METIS_NodeND METIS_WORKS)
  mark_as_advanced(METIS_WORKS)

  if(NOT METIS_WORKS)
    if(NOT METIS_FIND_QUIETLY)
      message(STATUS "Looking for METIS : test of METIS_NodeND with METIS library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  morse_cmake_required_unset()

  list(GET METIS_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" DIRECTORY)
  if (NOT METIS_LIBRARY_DIRS)
    set(METIS_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(METIS_PREFIX "${not_cached_dir}" CACHE PATH "Installation directory of METIS library" FORCE)
  else()
    set(METIS_PREFIX "${first_lib_path}" CACHE PATH "Installation directory of METIS library" FORCE)
  endif()
  mark_as_advanced(METIS_DIR)
  mark_as_advanced(METIS_PREFIX)

endif(METIS_LIBRARIES)


if (METIS_WORKS)
  morse_cmake_required_set(METIS)

  # Check the size of METIS_Idx
  # ---------------------------------
  include(CheckCSourceRuns)
  #stdio.h and stdint.h should be included by metis.h directly
  set(METIS_C_TEST_METIS_Idx_4 "
#include <stdio.h>
#include <stdint.h>
#include <metis.h>
int main() {
  if (sizeof(idx_t) == 4)
    return 0;
  else
    return 1;
}
")

  set(METIS_C_TEST_METIS_Idx_8 "
#include <stdio.h>
#include <stdint.h>
#include <metis.h>
int main() {
  if (sizeof(idx_t) == 8)
    return 0;
  else
    return 1;
}
")
  unset(METIS_Idx_4 CACHE)
  unset(METIS_Idx_8 CACHE)
  check_c_source_runs("${METIS_C_TEST_METIS_Idx_4}" METIS_Idx_4)
  check_c_source_runs("${METIS_C_TEST_METIS_Idx_8}" METIS_Idx_8)
  if(NOT METIS_Idx_4)
    if(NOT METIS_Idx_8)
      set(METIS_INTSIZE -1)
    else()
      set(METIS_INTSIZE 8)
    endif()
  else()
    set(METIS_INTSIZE 4)
  endif()

  morse_cmake_required_unset()
endif()

# check that METIS has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS DEFAULT_MSG
  METIS_LIBRARIES
  METIS_WORKS)

# Add imported target
if (METIS_FOUND)
  morse_create_imported_target(METIS)
endif()
