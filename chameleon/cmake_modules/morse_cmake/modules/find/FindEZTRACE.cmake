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
# - Find EZTRACE include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(EZTRACE
#               [REQUIRED] # Fail with error if eztrace is not found
#               [COMPONENTS <comp1> <comp2> ...]) # dependencies
#
#  EZTRACE depends on the following libraries:
#   - threads (e.g. pthreads)
#   - bfd
#   - cuda (optional)
#   - iberty (optional)
#   - mpi (optional)
#
#  COMPONENTS can be some of the following:
#   - MPI:  to detect the MPI version of EZTRACE
#
# This module finds headers and eztrace library using pkg-config file.
#  if found with pkg-config the following variables are set
#  <PREFIX>  = EZTRACE
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
# Set EZTRACE_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::EZTRACE``
#   The headers and libraries to use for EZTRACE, if found.
#
#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

if (NOT EZTRACE_FIND_QUIETLY)
  message(STATUS "FindEZTRACE needs pkg-config program and PKG_CONFIG_PATH set with eztrace.pc file path.")
endif()

# Set the version to find
set(EZTRACE_LOOK_FOR_MPI  OFF)

if( EZTRACE_FIND_COMPONENTS )
  foreach( component ${EZTRACE_FIND_COMPONENTS} )
    if (${component} STREQUAL "MPI")
      # means we look for eztrace-mpi library
      set(EZTRACE_LOOK_FOR_MPI ON)
    endif()
  endforeach()
endif()

# Use pkg-config to detect include/library dirs
# ---------------------------------------------
if (PKG_CONFIG_EXECUTABLE)
  unset(EZTRACE_FOUND CACHE)
  pkg_search_module(EZTRACE eztrace)

  if (NOT EZTRACE_FIND_QUIETLY)
    if (EZTRACE_FOUND AND EZTRACE_LIBRARIES)
      message(STATUS "Looking for EZTRACE - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for EZTRACE - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing eztrace.pc to"
        "\n   the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()

  if (EZTRACE_FOUND AND EZTRACE_LIBRARIES)

    # Set include dirs and libraries absolute path.
    if (NOT EZTRACE_INCLUDE_DIRS)
      pkg_get_variable(EZTRACE_INCLUDE_DIRS eztrace includedir)
    endif()

    # We have found EZTrace, now we look for the components.
    if(EZTRACE_LOOK_FOR_MPI)
      morse_find_library(EZTRACE_MPI
          LIBRARIES eztrace-mpi
          SUFFIXES  lib lib32 lib64)
      if (EZTRACE_MPI_LIBRARIES)
        message(STATUS "Looking for EZTRACE-MPI - found")
      else()
        message(STATUS "${Magenta}Looking for EZTRACE-MPI - not found.${ColourReset}")
      endif()

      list(APPEND EZTRACE_LINK_LIBRARIES ${EZTRACE_MPI_LIBRARIES})
      if (MPI_FOUND)
        if (MPI_C_INCLUDE_PATH)
          list(APPEND EZTRACE_INCLUDE_DIRS ${MPI_C_INCLUDE_PATH})
        endif()
        list(APPEND EZTRACE_LINK_LIBRARIES ${MPI_C_LIBRARIES})
      endif()
    endif(EZTRACE_LOOK_FOR_MPI)

    set(EZTRACE_FOUND_WITH_PKGCONFIG "TRUE")
    morse_find_pkgconfig_libraries_absolute_path(EZTRACE)
  else()
      set(EZTRACE_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

  if (EZTRACE_STATIC AND EZTRACE_STATIC_LIBRARIES)
    set (EZTRACE_DEPENDENCIES ${EZTRACE_STATIC_LIBRARIES})
    list (REMOVE_ITEM EZTRACE_DEPENDENCIES "eztrace")
    list (APPEND EZTRACE_LIBRARIES ${EZTRACE_DEPENDENCIES})
    set(EZTRACE_CFLAGS_OTHER ${EZTRACE_STATIC_CFLAGS_OTHER})
    set(EZTRACE_LDFLAGS_OTHER ${EZTRACE_STATIC_LDFLAGS_OTHER})
    if (NOT EZTRACE_FIND_QUIETLY)
      message(STATUS "EZTRACE_STATIC set to 1 by user, EZTRACE_LIBRARIES: ${EZTRACE_LIBRARIES}.")
    endif()
  endif()
endif()

# check a function to validate the find
if(EZTRACE_FOUND AND EZTRACE_LIBRARIES)

  # check if static or dynamic lib
  morse_check_static_or_dynamic(EZTRACE EZTRACE_LIBRARIES)
  if(EZTRACE_STATIC)
    set(STATIC "_STATIC")
  else()
    set(STATIC "")
  endif()

  # set required libraries for link
  morse_set_required_test_lib_link(EZTRACE)

  # test link
  unset(EZTRACE_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(eztrace_start EZTRACE_WORKS)
  mark_as_advanced(EZTRACE_WORKS)

  if(NOT EZTRACE_WORKS)
    if(NOT EZTRACE_FIND_QUIETLY)
      message(STATUS "Looking for eztrace : test of eztrace_topology_init with eztrace library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

endif(EZTRACE_FOUND AND EZTRACE_LIBRARIES)

# check that EZTRACE has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EZTRACE DEFAULT_MSG
  EZTRACE_LIBRARIES
  EZTRACE_WORKS)

# Add imported targe
if (EZTRACE_FOUND)
  morse_create_imported_target(EZTRACE)
endif()
