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
# - Find QUARK include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(QUARK
#               [REQUIRED]             # Fail with error if quark is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  QUARK depends on the following libraries:
#   - Threads
#
#  COMPONENTS are optional libraries QUARK could be linked with,
#  Use it to drive detection of a specific compilation chain
#  COMPONENTS can be some of the following:
#   - HWLOC: to activate detection of QUARK linked with hwloc
#
# This module finds headers and quark library.
# Results are reported in variables:
#  QUARK_FOUND             - True if headers and requested libraries were found
#  QUARK_PREFIX            - installation path of the lib found
#  QUARK_CFLAGS_OTHER      - quark compiler flags without headers paths
#  QUARK_LDFLAGS_OTHER     - quark linker flags without libraries
#  QUARK_INCLUDE_DIRS      - quark include directories
#  QUARK_LIBRARY_DIRS      - quark link directories
#  QUARK_LIBRARIES         - quark libraries to be linked (absolute path)
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DQUARK=path/to/quark):
#  QUARK_DIR              - Where to find the base directory of quark
#  QUARK_INCDIR           - Where to find the header files
#  QUARK_LIBDIR           - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: QUARK_DIR, QUARK_INCDIR, QUARK_LIBDIR
#
# Set QUARK_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::QUARK``
#   The headers and libraries to use for QUARK, if found.
#
#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

# Set variables from environment if needed
# ----------------------------------------
morse_find_package_get_envdir(QUARK)

# QUARK may depend on HWLOC
# try to find it specified as COMPONENTS during the call
set(QUARK_LOOK_FOR_HWLOC FALSE)

if( QUARK_FIND_COMPONENTS )
  foreach( component ${QUARK_FIND_COMPONENTS} )
    if(${component} STREQUAL "HWLOC")
      set(QUARK_LOOK_FOR_HWLOC TRUE)
    endif()
  endforeach()
endif()

# QUARK may depend on Threads, try to find it
if (QUARK_FIND_REQUIRED)
  find_package(Threads REQUIRED)
else()
  find_package(Threads)
endif()
if( THREADS_FOUND AND NOT THREADS_PREFER_PTHREAD_FLAG )
  libraries_absolute_path(CMAKE_THREAD_LIBS_INIT "")
endif ()

# QUARK may depend on HWLOC, try to find it
if (QUARK_LOOK_FOR_HWLOC)
  if (QUARK_FIND_REQUIRED AND QUARK_FIND_REQUIRED_HWLOC)
    find_package(HWLOC REQUIRED)
  else()
    find_package(HWLOC)
  endif()
endif()

# Looking for include
# -------------------
morse_find_path(QUARK
  HEADERS  quark.h
  SUFFIXES include include/quark include/plasma)

# Looking for lib
# ---------------
morse_find_library(QUARK
  LIBRARIES quark
  SUFFIXES lib lib32 lib64)

# check a function to validate the find
if(QUARK_LIBRARIES)

  # check if static or dynamic lib
  morse_check_static_or_dynamic(QUARK QUARK_LIBRARIES)
  if(QUARK_STATIC)
    set(STATIC "_STATIC")
  endif()

  set(REQUIRED_INCDIRS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)

  # QUARK
  if (QUARK_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${QUARK_INCLUDE_DIRS}")
  endif()
  if (QUARK_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${QUARK_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${QUARK_LIBRARIES}")
  # HWLOC
  if (HWLOC_FOUND AND QUARK_LOOK_FOR_HWLOC)
    if (HWLOC_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${HWLOC_INCLUDE_DIRS}")
    endif()
    if (HWLOC_CFLAGS_OTHER)
      list(APPEND REQUIRED_FLAGS "${HWLOC_CFLAGS_OTHER}")
    endif()
    if (HWLOC_LDFLAGS_OTHER)
      list(APPEND REQUIRED_LDFLAGS "${HWLOC_LDFLAGS_OTHER}")
    endif()
    if (HWLOC_LIBRARY_DIRS)
      list(APPEND REQUIRED_LIBDIRS "${HWLOC_LIBRARY_DIRS}")
    endif()
    list(APPEND REQUIRED_LIBS "${HWLOC_LIBRARIES}")
  endif()
  # THREADS
  if (THREADS_PREFER_PTHREAD_FLAG)
    list(APPEND REQUIRED_FLAGS "${CMAKE_THREAD_LIBS_INIT}")
    list(APPEND REQUIRED_LDFLAGS "${CMAKE_THREAD_LIBS_INIT}")
  else()
    list(APPEND REQUIRED_LIBS "${CMAKE_THREAD_LIBS_INIT}")
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
  unset(QUARK_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(QUARK_New QUARK_WORKS)
  mark_as_advanced(QUARK_WORKS)

  if(QUARK_WORKS)
    set(QUARK_INCLUDE_DIRS "${REQUIRED_INCDIRS}")
    set(QUARK_LIBRARY_DIRS "${REQUIRED_LIBDIRS}")
    set(QUARK_CFLAGS_OTHER "${REQUIRED_FLAGS}")
    set(QUARK_LDFLAGS_OTHER "${REQUIRED_LDFLAGS}")
    if (QUARK_STATIC OR HWLOC_STATIC)
      # save link with dependencies
      set(QUARK_LIBRARIES "${REQUIRED_LIBS}")
    endif()
  else()
    if(NOT QUARK_FIND_QUIETLY)
      message(STATUS "Looking for QUARK : test of QUARK_New with QUARK library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
      message(STATUS "Maybe QUARK is linked with specific libraries like. "
        "Have you tried with COMPONENTS (HWLOC)? See the explanation in FindQUARK.cmake.")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

  list(GET QUARK_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" DIRECTORY)
  if (NOT QUARK_LIBRARY_DIRS)
    set(QUARK_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(QUARK_PREFIX "${not_cached_dir}" CACHE PATH "Installation directory of QUARK library" FORCE)
  else()
    set(QUARK_PREFIX "${first_lib_path}" CACHE PATH "Installation directory of QUARK library" FORCE)
  endif()
  mark_as_advanced(QUARK_PREFIX)

endif(QUARK_LIBRARIES)

# check that QUARK has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(QUARK DEFAULT_MSG
  QUARK_LIBRARIES
  QUARK_WORKS)

# Add imported target
if (QUARK_FOUND)
  morse_create_imported_target(QUARK)
endif()
