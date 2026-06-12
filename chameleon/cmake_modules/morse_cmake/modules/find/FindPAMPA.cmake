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
# - Find PAMPA include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(PAMPA
#               [REQUIRED]             # Fail with error if pampa is not found
#              )
#
#  PAMPA depends on the following libraries:
#   - MPI
#   - PTSCOTCH
#
# This module finds headers and pampa library.
# Results are reported in variables:
#  PAMPA_FOUND             - True if headers and requested libraries were found
#  PAMPA_PREFIX            - installation path of the lib found
#  PAMPA_CFLAGS_OTHER      - pampa compiler flags without headers paths
#  PAMPA_LDFLAGS_OTHER     - pampa linker flags without libraries
#  PAMPA_INCLUDE_DIRS      - pampa include directories
#  PAMPA_LIBRARY_DIRS      - pampa link directories
#  PAMPA_LIBRARIES         - pampa libraries to be linked (absolute path)
#  PAMPA_INTSIZE         - Number of octets occupied by a PAMPA_Num
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DPAMPA=path/to/pampa):
#  PAMPA_DIR             - Where to find the base directory of pampa
#  PAMPA_INCDIR          - Where to find the header files
#  PAMPA_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: PAMPA_DIR, PAMPA_INCDIR, PAMPA_LIBDIR
#
# Set PAMPA_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::PAMPA``
#   The headers and libraries to use for PAMPA, if found.
#
#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

# Set variables from environment if needed
# ----------------------------------------
morse_find_package_get_envdir(PAMPA)

# PAMPA depends on MPI, try to find it
if (PAMPA_FIND_REQUIRED)
  find_package(MPI REQUIRED)
else()
  find_package(MPI)
endif()

# PAMPA depends on PAMPA, try to find it
if (PAMPA_FIND_REQUIRED)
  find_package(PTSCOTCH REQUIRED)
else()
  find_package(PTSCOTCH)
endif()

# Looking for include
# -------------------

# Try to find the pampa headers in the given paths
# ------------------------------------------------
morse_find_path(PAMPA
  HEADERS "pampa.h"
  SUFFIXES include include/pampa)

# Looking for lib
# ---------------
set(PAMPA_libs_to_find "pampa;pampaerr")
morse_find_library(PAMPA
  LIBRARIES ${PAMPA_libs_to_find}
  SUFFIXES lib lib32 lib64)

# check a function to validate the find
if(PAMPA_LIBRARIES)

  # check if static or dynamic lib
  morse_check_static_or_dynamic(PAMPA PAMPA_LIBRARIES)
  if(PAMPA_STATIC)
    set(STATIC "_STATIC")
  endif()

  set(REQUIRED_INCDIRS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)

  # PAMPA
  if (PAMPA_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS  "${PAMPA_INCLUDE_DIRS}")
  endif()
  if (PAMPA_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${PAMPA_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${PAMPA_LIBRARIES}")
  add_definitions(-Drestrict=__restrict) # pampa uses the restrict keyword
  list(APPEND REQUIRED_FLAGS "-Drestrict=__restrict")
  # MPI
  if (MPI_FOUND)
    if (MPI_C_INCLUDE_PATH)
      list(APPEND CMAKE_REQUIRED_INCLUDES "${MPI_C_INCLUDE_PATH}")
    endif()
    list(APPEND REQUIRED_LIBS "${MPI_C_LIBRARIES}")
  endif()
  # PTSCOTCH
  if (PTSCOTCH_INCLUDE_DIRS)
    list(APPEND REQUIRED_INCDIRS "${PTSCOTCH_INCLUDE_DIRS}")
  endif()
  if (PTSCOTCH_CFLAGS_OTHER)
    list(APPEND REQUIRED_FLAGS "${PTSCOTCH_CFLAGS_OTHER}")
  endif()
  if (PTSCOTCH_LDFLAGS_OTHER)
    list(APPEND REQUIRED_LDFLAGS "${PTSCOTCH_LDFLAGS_OTHER}")
  endif()
  foreach(libdir ${PTSCOTCH_LIBRARY_DIRS_DEP})
    if (libdir)
      list(APPEND REQUIRED_LIBDIRS "${libdir}")
    endif()
  endforeach()
  list(APPEND REQUIRED_LIBS "${PTSCOTCH_LIBRARIES}")
  # set required libraries for link
  set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
  if (REQUIRED_FLAGS)
    set(REQUIRED_FLAGS_COPY "${REQUIRED_FLAGS}")
    set(REQUIRED_FLAGS)
    set(REQUIRED_DEFINITIONS)
    foreach(_flag ${REQUIRED_FLAGS_COPY})
      if (_flag MATCHES "^-D")
       list(APPEND REQUIRED_DEFINITIONS "${_flag}")
      endif()
      string(REGEX REPLACE "^-D.*" "" _flag "${_flag}")
      list(APPEND REQUIRED_FLAGS "${_flag}")
    endforeach()
  endif()
  morse_finds_remove_duplicates()
  set(CMAKE_REQUIRED_DEFINITIONS "${REQUIRED_DEFINITIONS}")
  set(CMAKE_REQUIRED_FLAGS "${REQUIRED_LDFLAGS}")
  set(CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
  string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

  # test link
  unset(PAMPA_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(PAMPA_dmeshInit PAMPA_WORKS)
  mark_as_advanced(PAMPA_WORKS)

  if(PAMPA_WORKS)
    set(PAMPA_INCLUDE_DIRS "${REQUIRED_INCDIRS}")
    set(PAMPA_LIBRARY_DIRS "${REQUIRED_LIBDIRS}")
    set(PAMPA_CFLAGS_OTHER "${REQUIRED_FLAGS}")
    set(PAMPA_LDFLAGS_OTHER "${REQUIRED_LDFLAGS}")
    if (PAMPA_STATIC OR PTSCOTCH_STATIC)
      # save link with dependencies
      set(PAMPA_LIBRARIES "${REQUIRED_LIBS}")
    endif()
  else()
    if(NOT PAMPA_FIND_QUIETLY)
      message(STATUS "Looking for PAMPA : test of PAMPA_dmeshInit with PAMPA library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

  list(GET PAMPA_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" DIRECTORY)
  if (NOT PAMPA_LIBRARY_DIRS)
    set(PAMPA_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(PAMPA_PREFIX "${not_cached_dir}" CACHE PATH "Installation directory of PAMPA library" FORCE)
  else()
    set(PAMPA_PREFIX "${first_lib_path}" CACHE PATH "Installation directory of PAMPA library" FORCE)
  endif()
  mark_as_advanced(PAMPA_DIR)
  mark_as_advanced(PAMPA_PREFIX)

endif(PAMPA_LIBRARIES)

# Check the size of PAMPA_Num
# ---------------------------------
set(CMAKE_REQUIRED_INCLUDES ${PAMPA_INCLUDE_DIRS})

include(CheckCSourceRuns)
#stdio.h and stdint.h should be included by scotch.h directly
set(PAMPA_C_TEST_PAMPA_Num_4 "
#include <stdio.h>
#include <stdint.h>
#include <mpi.h>
#include <pampa.h>
int main(int argc, char **argv) {
  if (sizeof(PAMPA_Num) == 4)
    return 0;
  else
    return 1;
}
")

set(PAMPA_C_TEST_PAMPA_Num_8 "
#include <stdio.h>
#include <stdint.h>
#include <mpi.h>
#include <pampa.h>
int main(int argc, char **argv) {
  if (sizeof(PAMPA_Num) == 8)
    return 0;
  else
    return 1;
}
")
check_c_source_runs("${PAMPA_C_TEST_PAMPA_Num_4}" PAMPA_Num_4)
if(NOT PAMPA_Num_4)
  check_c_source_runs("${PAMPA_C_TEST_PAMPA_Num_8}" PAMPA_Num_8)
  if(NOT PAMPA_Num_8)
    set(PAMPA_INTSIZE -1)
  else()
    set(PAMPA_INTSIZE 8)
  endif()
else()
  set(PAMPA_INTSIZE 4)
endif()
set(CMAKE_REQUIRED_INCLUDES "")

# check that PAMPA has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PAMPA DEFAULT_MSG
  PAMPA_LIBRARIES
  PAMPA_WORKS)

# Add imported target
if (PAMPA_FOUND)
  morse_create_imported_target(PAMPA)
endif()
