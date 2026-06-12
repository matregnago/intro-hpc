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
# - Find PARMETIS include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(PARMETIS
#               [REQUIRED]             # Fail with error if parmetis is not found
#              )
#
#  PARMETIS depends on the following libraries:
#   - METIS
#   - MPI
#
# This module finds headers and parmetis library.
# Results are reported in variables:
#  PARMETIS_FOUND             - True if headers and requested libraries were found
#  PARMETIS_PREFIX            - installation path of the lib found
#  PARMETIS_CFLAGS_OTHER      - parmetis compiler flags without headers paths
#  PARMETIS_LDFLAGS_OTHER     - parmetis linker flags without libraries
#  PARMETIS_INCLUDE_DIRS      - parmetis include directories
#  PARMETIS_LIBRARY_DIRS      - parmetis link directories
#  PARMETIS_LIBRARIES         - parmetis libraries to be linked (absolute path)
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DPARMETIS_DIR=path/to/parmetis):
#  PARMETIS_DIR              - Where to find the base directory of parmetis
#  PARMETIS_INCDIR           - Where to find the header files
#  PARMETIS_LIBDIR           - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: PARMETIS_DIR, PARMETIS_INCDIR, PARMETIS_LIBDIR
#
# Set METIS_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::PARMETIS``
#   The headers and libraries to use for PARMETIS, if found.
#
#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

# Set variables from environment if needed
# ----------------------------------------
morse_find_package_get_envdir(PARMETIS)

# PARMETIS depends on METIS, try to find it
if(PARMETIS_FIND_REQUIRED)
  find_package(METIS REQUIRED)
else()
  find_package(METIS)
endif()

# PARMETIS depends on MPI, try to find it
if(PARMETIS_FIND_REQUIRED)
  find_package(MPI REQUIRED)
else()
  find_package(MPI)
endif()

# Looking for include
# -------------------
morse_find_path(PARMETIS
  HEADERS parmetis.h
  SUFFIXES include include/parmetis)

# Looking for lib
# ---------------
morse_find_library(PARMETIS
  LIBRARIES parmetis
  SUFFIXES lib lib32 lib64)

# check a function to validate the find
if(PARMETIS_LIBRARIES)

  # check if static or dynamic lib
  morse_check_static_or_dynamic(PARMETIS PARMETIS_LIBRARIES)
  if(PARMETIS_STATIC)
    set(STATIC "_STATIC")
  endif()

  set(REQUIRED_INCDIRS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)

  # PARMETIS
  if (PARMETIS_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${PARMETIS_INCLUDE_DIRS}")
  endif()
  if (PARMETIS_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${PARMETIS_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${PARMETIS_LIBRARIES}")
  # METIS
  if (METIS_FOUND)
    if (METIS_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${METIS_INCLUDE_DIRS}")
    endif()
    foreach(libdir ${METIS_LIBRARY_DIRS})
      if (libdir)
        list(APPEND REQUIRED_LIBDIRS "${libdir}")
      endif()
    endforeach()
    list(APPEND REQUIRED_LIBS "${METIS_LIBRARIES}")
  endif()
  # MPI
  if (MPI_FOUND)
    if (MPI_C_INCLUDE_PATH)
      list(APPEND REQUIRED_INCDIRS "${MPI_C_INCLUDE_PATH}")
    endif()
    list(APPEND REQUIRED_LIBS "${MPI_C_LIBRARIES}")
  endif()
  # m
  find_library(M_LIBRARY NAMES m)
  mark_as_advanced(M_LIBRARY)
  if(M_LIBRARY)
    list(APPEND REQUIRED_LIBS "${M_LIBRARY}")
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
  unset(PARMETIS_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(ParMETIS_V3_NodeND PARMETIS_WORKS)
  mark_as_advanced(PARMETIS_WORKS)

  if(PARMETIS_WORKS)
    set(PARMETIS_INCLUDE_DIRS "${REQUIRED_INCDIRS}")
    set(PARMETIS_LIBRARY_DIRS "${REQUIRED_LIBDIRS}")
    set(PARMETIS_CFLAGS_OTHER "${REQUIRED_FLAGS}")
    set(PARMETIS_LDFLAGS_OTHER "${REQUIRED_LDFLAGS}")
    if (PARMETIS_STATIC OR METIS_STATIC)
      # save link with dependencies
      set(PARMETIS_LIBRARIES "${REQUIRED_LIBS}")
    endif()
  else()
    if(NOT PARMETIS_FIND_QUIETLY)
      message(STATUS "Looking for PARMETIS : test of ParMETIS_V3_NodeND with PARMETIS library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

  list(GET PARMETIS_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" DIRECTORY)
  if (NOT PARMETIS_LIBRARY_DIRS)
    set(PARMETIS_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(PARMETIS_PREFIX "${not_cached_dir}" CACHE PATH "Installation directory of PARMETIS library" FORCE)
  else()
    set(PARMETIS_PREFIX "${first_lib_path}" CACHE PATH "Installation directory of PARMETIS library" FORCE)
  endif()
  mark_as_advanced(PARMETIS_DIR)
  mark_as_advanced(PARMETIS_PREFIX)

endif(PARMETIS_LIBRARIES)

# check that PARMETIS has been found
# ----------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PARMETIS DEFAULT_MSG
  PARMETIS_LIBRARIES
  PARMETIS_WORKS)

# Add imported target
if (PARMETIS_FOUND)
  morse_create_imported_target(PARMETIS)
endif()
