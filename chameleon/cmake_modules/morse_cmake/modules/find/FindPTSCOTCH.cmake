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
# - Find PTSCOTCH include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(PTSCOTCH
#               [REQUIRED]             # Fail with error if ptscotch is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  PTSCOTCH depends on the following libraries:
#   - Threads
#   - MPI
#
#  COMPONENTS can be some of the following:
#   - ERREXIT: to choose libptscotcherrexit as error library instead of libptscotcherr (default)
#   - ESMUMPS: to activate detection of PT-Scotch with the esmumps interface
#
# This module finds headers and ptscotch library.
# Results are reported in variables:
#  PTSCOTCH_FOUND             - True if headers and requested libraries were found
#  PTSCOTCH_PREFIX            - installation path of the lib found
#  PTSCOTCH_CFLAGS_OTHER      - ptscotch compiler flags without headers paths
#  PTSCOTCH_LDFLAGS_OTHER     - ptscotch linker flags without libraries
#  PTSCOTCH_INCLUDE_DIRS      - ptscotch include directories
#  PTSCOTCH_LIBRARY_DIRS      - ptscotch link directories
#  PTSCOTCH_LIBRARIES         - ptscotch libraries to be linked (absolute path)
#  PTSCOTCH_INTSIZE          - Number of octets occupied by a SCOTCH_Num
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DPTSCOTCH=path/to/ptscotch):
#  PTSCOTCH_DIR              - Where to find the base directory of ptscotch
#  PTSCOTCH_INCDIR           - Where to find the header files
#  PTSCOTCH_LIBDIR           - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: PTSCOTCH_DIR, PTSCOTCH_INCDIR, PTSCOTCH_LIBDIR
#
# Set PTSCOTCH_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::PTSCOTCH``
#   The headers and libraries to use for PTSCOTCH, if found.
#
#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

# Set variables from environment if needed
# ----------------------------------------
morse_find_package_get_envdir(PTSCOTCH)

# Set the version to find
set(PTSCOTCH_LOOK_FOR_ERREXIT OFF)
set(PTSCOTCH_LOOK_FOR_ESMUMPS OFF)

if( PTSCOTCH_FIND_COMPONENTS )
  foreach( component ${PTSCOTCH_FIND_COMPONENTS} )
    if (${component} STREQUAL "ERREXIT")
      # means we look for ptscotcherrexit library
      set(PTSCOTCH_LOOK_FOR_ERREXIT ON)
    endif()
    if (${component} STREQUAL "ESMUMPS")
      # means we look for esmumps library
      set(PTSCOTCH_LOOK_FOR_ESMUMPS ON)
    endif()
  endforeach()
endif()

# PTSCOTCH depends on Threads, try to find it
if (PTSCOTCH_FIND_REQUIRED)
  find_package(Threads REQUIRED)
else()
  find_package(Threads)
endif()
if( THREADS_FOUND AND NOT THREADS_PREFER_PTHREAD_FLAG )
  libraries_absolute_path(CMAKE_THREAD_LIBS_INIT "")
endif ()

# PTSCOTCH depends on MPI, try to find it
if (PTSCOTCH_FIND_REQUIRED)
  find_package(MPI REQUIRED)
else()
  find_package(MPI)
endif()

# Looking for include
# -------------------
set(PTSCOTCH_hdrs_to_find "ptscotch.h;scotch.h")
morse_find_path(PTSCOTCH
  HEADERS  ${PTSCOTCH_hdrs_to_find}
  SUFFIXES include include/scotch include/ptscotch)

# Looking for lib
# ---------------
if (PTSCOTCH_LOOK_FOR_ERREXIT)
  list(APPEND PTSCOTCH_libs_to_find "ptscotch;ptscotcherrexit;scotch;scotcherrexit")
else()
  list(APPEND PTSCOTCH_libs_to_find "ptscotch;ptscotcherr;scotch;scotcherr")
endif()
if (PTSCOTCH_LOOK_FOR_ESMUMPS)
  list(INSERT PTSCOTCH_libs_to_find 0 "ptesmumps")
  list(APPEND PTSCOTCH_libs_to_find   "esmumps"  )
endif()

morse_find_library(PTSCOTCH
  LIBRARIES ${PTSCOTCH_libs_to_find}
  SUFFIXES  lib lib32 lib64)

# check a function to validate the find
if(PTSCOTCH_LIBRARIES)

  # check if static or dynamic lib
  morse_check_static_or_dynamic(PTSCOTCH PTSCOTCH_LIBRARIES)
  if(PTSCOTCH_STATIC)
    set(STATIC "_STATIC")
  endif()

  set(REQUIRED_INCDIRS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)

  # PTSCOTCH
  if (PTSCOTCH_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${PTSCOTCH_INCLUDE_DIRS}")
  endif()
  if (PTSCOTCH_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${PTSCOTCH_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${PTSCOTCH_LIBRARIES}")
  # MPI
  if (MPI_FOUND)
    if (MPI_C_INCLUDE_PATH)
      list(APPEND REQUIRED_INCDIRS "${MPI_C_INCLUDE_PATH}")
    endif()
    list(APPEND REQUIRED_LIBS "${MPI_C_LIBRARIES}")
  endif()
  # THREADS
  if (THREADS_PREFER_PTHREAD_FLAG)
    list(APPEND REQUIRED_FLAGS "${CMAKE_THREAD_LIBS_INIT}")
    list(APPEND REQUIRED_LDFLAGS "${CMAKE_THREAD_LIBS_INIT}")
  else()
    list(APPEND REQUIRED_LIBS "${CMAKE_THREAD_LIBS_INIT}")
  endif()
  set(Z_LIBRARY "Z_LIBRARY-NOTFOUND")
  find_library(Z_LIBRARY NAMES z)
  mark_as_advanced(Z_LIBRARY)
  if(Z_LIBRARY)
    list(APPEND REQUIRED_LIBS "${Z_LIBRARY}")
  endif()
  set(M_LIBRARY "M_LIBRARY-NOTFOUND")
  find_library(M_LIBRARY NAMES m)
  mark_as_advanced(M_LIBRARY)
  if(M_LIBRARY)
    list(APPEND REQUIRED_LIBS "${M_LIBRARY}")
  endif()
  set(RT_LIBRARY "RT_LIBRARY-NOTFOUND")
  find_library(RT_LIBRARY NAMES rt)
  mark_as_advanced(RT_LIBRARY)
  if(RT_LIBRARY)
    list(APPEND REQUIRED_LIBS "${RT_LIBRARY}")
  endif()

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
  unset(PTSCOTCH_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(SCOTCH_dgraphInit PTSCOTCH_WORKS)
  mark_as_advanced(PTSCOTCH_WORKS)

  if(PTSCOTCH_WORKS)
    set(PTSCOTCH_INCLUDE_DIRS "${REQUIRED_INCDIRS}")
    set(PTSCOTCH_LIBRARY_DIRS "${REQUIRED_LIBDIRS}")
    set(PTSCOTCH_CFLAGS_OTHER "${REQUIRED_FLAGS}")
    set(PTSCOTCH_LDFLAGS_OTHER "${REQUIRED_LDFLAGS}")
    # scotch is not giving its dependencies even in shared library mode ...
    #if (PTSCOTCH_STATIC)
      # save link with dependencies
      set(PTSCOTCH_LIBRARIES "${REQUIRED_LIBS}")
    #endif()
  else()
    if(NOT PTSCOTCH_FIND_QUIETLY)
      message(STATUS "Looking for PTSCOTCH : test of SCOTCH_dgraphInit with PTSCOTCH library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

  list(GET PTSCOTCH_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" DIRECTORY)
  if (NOT PTSCOTCH_LIBRARY_DIRS)
    set(PTSCOTCH_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(PTSCOTCH_PREFIX "${not_cached_dir}" CACHE PATH "Installation directory of PTSCOTCH library" FORCE)
  else()
    set(PTSCOTCH_PREFIX "${first_lib_path}" CACHE PATH "Installation directory of PTSCOTCH library" FORCE)
  endif()
  mark_as_advanced(PTSCOTCH_DIR)
  mark_as_advanced(PTSCOTCH_PREFIX)

endif(PTSCOTCH_LIBRARIES)

# Check the size of SCOTCH_Num
# ---------------------------------
set(CMAKE_REQUIRED_INCLUDES ${PTSCOTCH_INCLUDE_DIRS})

include(CheckCSourceRuns)
#stdio.h and stdint.h should be included by scotch.h directly
#mpi.h not included into ptscotch.h => MPI_comm undefined
set(PTSCOTCH_C_TEST_SCOTCH_Num_4 "
#include <stdio.h>
#include <stdint.h>
#include <mpi.h>
#include <ptscotch.h>
int main() {
  if (sizeof(SCOTCH_Num) == 4)
    return 0;
  else
    return 1;
}
")

set(PTSCOTCH_C_TEST_SCOTCH_Num_8 "
#include <stdio.h>
#include <stdint.h>
#include <mpi.h>
#include <ptscotch.h>
int main() {
  if (sizeof(SCOTCH_Num) == 8)
    return 0;
  else
    return 1;
}
")

unset(PTSCOTCH_Num_4 CACHE)
unset(PTSCOTCH_Num_8 CACHE)
check_c_source_runs("${PTSCOTCH_C_TEST_SCOTCH_Num_4}" PTSCOTCH_Num_4)
check_c_source_runs("${PTSCOTCH_C_TEST_SCOTCH_Num_8}" PTSCOTCH_Num_8)
if(NOT PTSCOTCH_Num_4)
  if(NOT PTSCOTCH_Num_8)
    set(PTSCOTCH_INTSIZE -1)
  else()
    set(PTSCOTCH_INTSIZE 8)
  endif()
else()
  set(PTSCOTCH_INTSIZE 4)
endif()
set(CMAKE_REQUIRED_INCLUDES "")

# check that PTSCOTCH has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PTSCOTCH DEFAULT_MSG
  PTSCOTCH_LIBRARIES
  PTSCOTCH_WORKS)

# Add imported target
if (PTSCOTCH_FOUND)
  morse_create_imported_target(PTSCOTCH)
endif()
