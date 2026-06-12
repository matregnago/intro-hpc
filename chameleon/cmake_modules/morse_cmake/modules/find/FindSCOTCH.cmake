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
# - Find SCOTCH include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(SCOTCH
#               [REQUIRED]             # Fail with error if scotch is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  COMPONENTS can be some of the following:
#   - ERREXIT: to choose libscotcherrexit as error library instead of libscotcherr (default)
#   - ESMUMPS: to activate detection of Scotch with the esmumps interface
#
# This module finds headers and scotch library.
# Results are reported in variables:
#  SCOTCH_FOUND           - True if headers and requested libraries were found
#  SCOTCH_PREFIX          - installation path of the lib found
#  SCOTCH_INCLUDE_DIRS    - scotch include directories
#  SCOTCH_LIBRARY_DIRS    - Link directories for scotch libraries
#  SCOTCH_LIBRARIES       - scotch component libraries to be linked
#  SCOTCH_INTSIZE         - Number of octets occupied by a SCOTCH_Num
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DSCOTCH=path/to/scotch):
#  SCOTCH_DIR             - Where to find the base directory of scotch
#  SCOTCH_INCDIR          - Where to find the header files
#  SCOTCH_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: SCOTCH_DIR, SCOTCH_INCDIR, SCOTCH_LIBDIR
#
# Set SCOTCH_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::SCOTCH``
#   The headers and libraries to use for SCOTCH, if found.
#
#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

# Set variables from environment if needed
# ----------------------------------------
morse_find_package_get_envdir(SCOTCH)

# Set the version to find
set(SCOTCH_LOOK_FOR_ERREXIT OFF)
set(SCOTCH_LOOK_FOR_ESMUMPS OFF)

if( SCOTCH_FIND_COMPONENTS )
  foreach( component ${SCOTCH_FIND_COMPONENTS} )
    if (${component} STREQUAL "ERREXIT")
      # means we look for scotcherrexit library
      set(SCOTCH_LOOK_FOR_ERREXIT ON)
    endif()
    if (${component} STREQUAL "ESMUMPS")
      # means we look for esmumps library
      set(SCOTCH_LOOK_FOR_ESMUMPS ON)
    endif()
  endforeach()
endif()

# SCOTCH may depend on Threads, try to find it
if (SCOTCH_FIND_REQUIRED)
  find_package(Threads REQUIRED)
else()
  find_package(Threads)
endif()
if( THREADS_FOUND AND NOT THREADS_PREFER_PTHREAD_FLAG)
  libraries_absolute_path(CMAKE_THREAD_LIBS_INIT "")
endif ()

# Looking for include
# -------------------
set(SCOTCH_hdrs_to_find "scotch.h")
morse_find_path(SCOTCH
  HEADERS  ${SCOTCH_hdrs_to_find}
  SUFFIXES include include/scotch)

# Looking for lib
# ---------------
if (SCOTCH_LOOK_FOR_ERREXIT)
  list(APPEND SCOTCH_libs_to_find "scotch;scotcherrexit")
else()
  list(APPEND SCOTCH_libs_to_find "scotch;scotcherr")
endif()
if (SCOTCH_LOOK_FOR_ESMUMPS)
  list(INSERT SCOTCH_libs_to_find 0 "esmumps")
endif()

morse_find_library(SCOTCH
  LIBRARIES ${SCOTCH_libs_to_find}
  SUFFIXES  lib lib32 lib64)

# check a function to validate the find
if(SCOTCH_LIBRARIES)

  # check if static or dynamic lib
  morse_check_static_or_dynamic(SCOTCH SCOTCH_LIBRARIES)
  if(SCOTCH_STATIC)
    set(STATIC "_STATIC")
  endif()

  set(REQUIRED_INCDIRS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)

  # SCOTCH
  if (SCOTCH_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS  "${SCOTCH_INCLUDE_DIRS}")
  endif()
  if (SCOTCH_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${SCOTCH_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${SCOTCH_LIBRARIES}")
  # THREADS
  if (THREADS_PREFER_PTHREAD_FLAG)
    list(APPEND REQUIRED_FLAGS "${CMAKE_THREAD_LIBS_INIT}")
    list(APPEND REQUIRED_LDFLAGS "${CMAKE_THREAD_LIBS_INIT}")
  else()
    if (CMAKE_THREAD_LIBS_INIT)
      list(APPEND REQUIRED_LIBS "${CMAKE_THREAD_LIBS_INIT}")
    endif()
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
  morse_finds_remove_duplicates()
  # set required libraries for link
  set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
  set(CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
  string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

  # test link
  unset(SCOTCH_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(SCOTCH_graphInit SCOTCH_WORKS)
  mark_as_advanced(SCOTCH_WORKS)

  # test scotch version
  unset(HAVE_SCOTCH_CONTEXT_INIT CACHE)
  check_function_exists(SCOTCH_contextInit HAVE_SCOTCH_CONTEXT_INIT)
  mark_as_advanced(HAVE_SCOTCH_CONTEXT_INIT)

  if(SCOTCH_WORKS)
    set(SCOTCH_LIBRARY_DIRS "${REQUIRED_LIBDIRS}")
    set(SCOTCH_INCLUDE_DIRS "${REQUIRED_INCDIRS}")
    set(SCOTCH_CFLAGS_OTHER "${REQUIRED_FLAGS}")
    set(SCOTCH_LDFLAGS_OTHER "${REQUIRED_LDFLAGS}")
    # scotch is not giving its dependencies even in shared library mode ...
    #if (SCOTCH_STATIC)
      # save link with dependencies
      set(SCOTCH_LIBRARIES "${REQUIRED_LIBS}")
    #endif()
  else()
    if(NOT SCOTCH_FIND_QUIETLY)
      message(STATUS "Looking for SCOTCH : test of SCOTCH_graphInit with SCOTCH library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

  list(GET SCOTCH_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" DIRECTORY)
  if (NOT SCOTCH_LIBRARY_DIRS)
    set(SCOTCH_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(SCOTCH_PREFIX "${not_cached_dir}" CACHE PATH "Installation directory of SCOTCH library" FORCE)
   else()
    set(SCOTCH_PREFIX "${first_lib_path}" CACHE PATH "Installation directory of SCOTCH library" FORCE)
  endif()
  mark_as_advanced(SCOTCH_DIR)
  mark_as_advanced(SCOTCH_PREFIX)

endif(SCOTCH_LIBRARIES)

# Check the size of SCOTCH_Num
# ---------------------------------
set(CMAKE_REQUIRED_INCLUDES ${SCOTCH_INCLUDE_DIRS})

include(CheckCSourceRuns)
#stdio.h and stdint.h should be included by scotch.h directly
set(SCOTCH_C_TEST_SCOTCH_Num_4 "
#include <stdio.h>
#include <stdint.h>
#include <scotch.h>
int main() {
  if (sizeof(SCOTCH_Num) == 4)
    return 0;
  else
    return 1;
}
")

set(SCOTCH_C_TEST_SCOTCH_Num_8 "
#include <stdio.h>
#include <stdint.h>
#include <scotch.h>
int main() {
  if (sizeof(SCOTCH_Num) == 8)
    return 0;
  else
    return 1;
}
")

unset(SCOTCH_Num_4 CACHE)
unset(SCOTCH_Num_8 CACHE)
check_c_source_runs("${SCOTCH_C_TEST_SCOTCH_Num_4}" SCOTCH_Num_4)
check_c_source_runs("${SCOTCH_C_TEST_SCOTCH_Num_8}" SCOTCH_Num_8)
if(NOT SCOTCH_Num_4)
  if(NOT SCOTCH_Num_8)
    set(SCOTCH_INTSIZE -1)
  else()
    set(SCOTCH_INTSIZE 8)
  endif()
else()
  set(SCOTCH_INTSIZE 4)
endif()
set(CMAKE_REQUIRED_INCLUDES "")

# check that SCOTCH has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCOTCH DEFAULT_MSG
  SCOTCH_LIBRARIES
  SCOTCH_WORKS)

# Add imported target
if (SCOTCH_FOUND)
  morse_create_imported_target(SCOTCH)
endif()
