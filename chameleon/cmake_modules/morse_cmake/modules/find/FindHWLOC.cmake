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
# - Find HWLOC include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(HWLOC
#               [REQUIRED]) # Fail with error if hwloc is not found
#
# This module finds headers and hwloc library using pkg-config file.
#  if found with pkg-config the following variables are set
#  <PREFIX>  = HWLOC
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
# Set HWLOC_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::HWLOC``
#   The headers and libraries to use for HWLOC, if found.
#
#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

include(CheckStructHasMember)
include(CheckCSourceCompiles)

if (NOT HWLOC_FIND_QUIETLY)
  message(STATUS "FindHWLOC needs pkg-config program and PKG_CONFIG_PATH set HWLOC.pc file path.")
endif()

# use pkg-config to detect include/library dirs (if pkg-config is available)
# --------------------------------------------------------------------------
if (PKG_CONFIG_EXECUTABLE)
  unset(HWLOC_FOUND CACHE)
  pkg_search_module(HWLOC hwloc)

  if (NOT HWLOC_FIND_QUIETLY)
    if (HWLOC_FOUND AND HWLOC_LIBRARIES)
      message(STATUS "Looking for HWLOC - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for HWLOC - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing hwloc.pc to"
        "\n   the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()
  if (HWLOC_FOUND AND HWLOC_LIBRARIES)
    if (NOT HWLOC_INCLUDE_DIRS)
      pkg_get_variable(HWLOC_INCLUDE_DIRS hwloc includedir)
    endif()
    set(HWLOC_FOUND_WITH_PKGCONFIG "TRUE")
    morse_find_pkgconfig_libraries_absolute_path(HWLOC)
  else()
    set(HWLOC_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

  if (HWLOC_STATIC AND HWLOC_STATIC_LIBRARIES)
    set (HWLOC_DEPENDENCIES ${HWLOC_STATIC_LIBRARIES})
    list (REMOVE_ITEM HWLOC_DEPENDENCIES "hwloc")
    list (APPEND HWLOC_LIBRARIES ${HWLOC_DEPENDENCIES})
    set(HWLOC_CFLAGS_OTHER ${HWLOC_STATIC_CFLAGS_OTHER})
    set(HWLOC_LDFLAGS_OTHER ${HWLOC_STATIC_LDFLAGS_OTHER})
    if (NOT HWLOC_FIND_QUIETLY)
      message(STATUS "HWLOC_STATIC set to 1 by user, HWLOC_LIBRARIES: ${HWLOC_LIBRARIES}.")
    endif()
  endif()
endif()

# check a function to validate the find
if(HWLOC_FOUND AND HWLOC_LIBRARIES)

  # check if static or dynamic lib
  morse_check_static_or_dynamic(HWLOC HWLOC_LIBRARIES)
  if(HWLOC_STATIC)
    set(STATIC "_STATIC")
  else()
    set(STATIC "")
  endif()

  # set required libraries for link
  morse_set_required_test_lib_link(HWLOC)

  # test link
  unset(HWLOC_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(hwloc_topology_init HWLOC_WORKS)
  mark_as_advanced(HWLOC_WORKS)

  if(NOT HWLOC_WORKS)
    if(NOT HWLOC_FIND_QUIETLY)
      message(STATUS "Looking for hwloc : test of hwloc_topology_init with hwloc library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

  # Add imported target
  morse_create_imported_target(HWLOC)

endif(HWLOC_FOUND AND HWLOC_LIBRARIES)

# check that HWLOC has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HWLOC DEFAULT_MSG
  HWLOC_LIBRARIES
  HWLOC_WORKS)

# Add imported targe
if (HWLOC_FOUND)
  morse_create_imported_target(HWLOC)
endif()