###
# WARNING: only HYPRE lib is searched for now
# it is surely too simple, must be completed
###
#
# @copyright (c) 2016-2020 Inria. All rights reserved.
#
# Copyright 2016-2020 Florent Pruvost
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
# - Find HYPRE include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(HYPRE
#               [REQUIRED]) # Fail with error if hypre is not found
#
#  HYPRE depends on the following libraries:
#   - MPI
#
# This module finds headers and hypre library.
# Results are reported in variables:
#  HYPRE_FOUND             - True if headers and requested libraries were found
#  HYPRE_INCLUDE_DIRS      - hypre include directories
#  HYPRE_LIBRARY_DIRS      - Link directories for hypre libraries
#  HYPRE_LIBRARIES         - hypre component libraries to be linked
#
# Set HYPRE_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::HYPRE``
#   The headers and libraries to use for HYPRE, if found.
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DHYPRE_DIR=path/to/hypre):
#  HYPRE_DIR             - Where to find the base directory of hypre
#  HYPRE_INCDIR          - Where to find the header files
#  HYPRE_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: HYPRE_DIR, HYPRE_INCDIR, HYPRE_LIBDIR

#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

# Set variables from environment if needed
# ----------------------------------------
morse_find_package_get_envdir(HYPRE)

set(HYPRE_GIVEN_BY_USER "FALSE")
if ( HYPRE_DIR OR ( HYPRE_INCDIR AND HYPRE_LIBDIR ) )
  set(HYPRE_GIVEN_BY_USER "TRUE")
endif()

# HYPRE depends on MPI, try to find it
if (HYPRE_FIND_REQUIRED)
  find_package(MPI REQUIRED)
else()
  find_package(MPI)
endif()

# Looking for include
# -------------------
morse_find_path(HYPRE
  HEADERS HYPRE.h
  SUFFIXES include include/hypre)

# Looking for lib
# ---------------
morse_find_library( HYPRE
  LIBRARIES HYPRE
  SUFFIXES lib lib32 lib64)

# check a function to validate the find
if(HYPRE_LIBRARIES)

  # check if static or dynamic lib
  morse_check_static_or_dynamic(HYPRE HYPRE_LIBRARIES)
  if(HYPRE_STATIC)
    set(STATIC "_STATIC")
  endif()

  set(REQUIRED_INCDIRS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)

  # HYPRE
  if (HYPRE_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${HYPRE_INCLUDE_DIRS}")
  endif()
  if (HYPRE_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${HYPRE_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${HYPRE_LIBRARIES}")
  # MPI
  if (MPI_FOUND)
    if (MPI_C_INCLUDE_PATH)
      list(APPEND REQUIRED_INCDIRS "${MPI_C_INCLUDE_PATH}")
    endif()
    list(APPEND REQUIRED_LIBS "${MPI_C_LIBRARIES}")
  endif()
  # libm
  set(M_LIBRARY "M_LIBRARY-NOTFOUND")
  find_library(M_LIBRARY NAMES m)
  mark_as_advanced(M_LIBRARY)
  if(M_LIBRARY)
    list(APPEND REQUIRED_LIBS "${M_LIBRARY}")
  endif()
  # libstdc++
  set(stdcpp_LIBRARY "stdcpp_LIBRARY-NOTFOUND")
  find_library(stdcpp_LIBRARY NAMES stdc++ NO_SYSTEM_ENVIRONMENT_PATH)
  mark_as_advanced(stdcpp_LIBRARY)
  if(stdcpp_LIBRARY)
    list(APPEND REQUIRED_LIBS "${stdcpp_LIBRARY}")
  endif()

  # set required libraries for link
  morse_finds_remove_duplicates()
  set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
  set(CMAKE_REQUIRED_LIBRARIES)
  foreach(lib_dir ${REQUIRED_LIBDIRS})
    list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${lib_dir}")
  endforeach()
  list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
  string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

  # test link
  unset(HYPRE_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(HYPRE_StructGridCreate HYPRE_WORKS)
  mark_as_advanced(HYPRE_WORKS)

  if(HYPRE_WORKS)
    if (HYPRE_STATIC)
      # save link with dependencies
      set(HYPRE_LIBRARIES "${REQUIRED_LIBS}")
      set(HYPRE_LIBRARY_DIRS "${REQUIRED_LIBDIRS}")
      set(HYPRE_INCLUDE_DIRS "${REQUIRED_INCDIRS}")
    endif()
  else()
    if(NOT HYPRE_FIND_QUIETLY)
      message(STATUS "Looking for HYPRE : test of HYPRE_StructGridCreate with HYPRE library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

  if (HYPRE_LIBRARY_DIRS)
    list(GET HYPRE_LIBRARY_DIRS 0 first_lib_path)
  else()
    list(GET HYPRE_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" DIRECTORY)
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(HYPRE_PREFIX "${not_cached_dir}" CACHE PATH "Installation directory of HYPRE library" FORCE)
  else()
    set(HYPRE_PREFIX "${first_lib_path}" CACHE PATH "Installation directory of HYPRE library" FORCE)
  endif()
  mark_as_advanced(HYPRE_DIR)
  mark_as_advanced(HYPRE_PREFIX)

endif(HYPRE_LIBRARIES)

# check that HYPRE has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HYPRE DEFAULT_MSG
  HYPRE_LIBRARIES
  HYPRE_WORKS)

# Add imported target
if (HYPRE_FOUND)
  morse_create_imported_target(HYPRE)
endif()
