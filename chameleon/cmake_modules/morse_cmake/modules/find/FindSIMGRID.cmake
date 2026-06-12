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
# - Find SIMGRID include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(SIMGRID
#               [REQUIRED]) # Fail with error if simgrid is not found
#
# This module finds headers and simgrid library using pkg-config file.
#  if found with pkg-config the following variables are set
#  <PREFIX>  = SIMGRID
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
# Set SIMGRID_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::SIMGRID``
#   The headers and libraries to use for SIMGRID, if found.
#
#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

if (NOT SIMGRID_FIND_QUIETLY)
  message(STATUS "FindSIMGRID needs pkg-config program and PKG_CONFIG_PATH set with simgrid.pc file path.")
endif()

# use pkg-config to detect include/library dirs (if pkg-config is available)
# --------------------------------------------------------------------------
if (PKG_CONFIG_EXECUTABLE)
  unset(SIMGRID_FOUND CACHE)
  pkg_search_module(SIMGRID simgrid)

  if (NOT SIMGRID_FIND_QUIETLY)
    if (SIMGRID_FOUND AND SIMGRID_LIBRARIES)
      message(STATUS "Looking for SIMGRID - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for SIMGRID - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing simgrid.pc to the"
        "\n   PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()
  if (SIMGRID_FOUND AND SIMGRID_LIBRARIES)
    if (NOT SIMGRID_INCLUDE_DIRS)
      pkg_get_variable(SIMGRID_INCLUDE_DIRS simgrid includedir)
    endif()
    set(SIMGRID_FOUND_WITH_PKGCONFIG "TRUE")
    morse_find_pkgconfig_libraries_absolute_path(SIMGRID)
  else()
    set(SIMGRID_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

  if (SIMGRID_STATIC AND SIMGRID_STATIC_LIBRARIES)
    set (SIMGRID_DEPENDENCIES ${SIMGRID_STATIC_LIBRARIES})
    list (REMOVE_ITEM SIMGRID_DEPENDENCIES "simgrid")
    list (APPEND SIMGRID_LIBRARIES ${SIMGRID_DEPENDENCIES})
    set(SIMGRID_CFLAGS_OTHER ${SIMGRID_STATIC_CFLAGS_OTHER})
    set(SIMGRID_LDFLAGS_OTHER ${SIMGRID_STATIC_LDFLAGS_OTHER})
    if (NOT SIMGRID_FIND_QUIETLY)
      message(STATUS "SIMGRID_STATIC set to 1 by user, SIMGRID_LIBRARIES: ${SIMGRID_LIBRARIES}.")
    endif()
  endif()
endif()

# check a function to validate the find
if(SIMGRID_FOUND AND SIMGRID_LIBRARIES)

  # check if static or dynamic lib
  morse_check_static_or_dynamic(SIMGRID SIMGRID_LIBRARIES)
  if(SIMGRID_STATIC)
    set(STATIC "_STATIC")
  else()
    set(STATIC "")
  endif()

  # set required libraries for link
  morse_set_required_test_lib_link(SIMGRID)

  # test link
  unset(SIMGRID_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(simgrid_init SIMGRID_WORKS)
  mark_as_advanced(SIMGRID_WORKS)

  if(SIMGRID_WORKS)
    if (SIMGRID_STATIC)
      # save link with dependencies
      set(SIMGRID_STATIC_LIBRARIES "${REQUIRED_LIBS}")
      set(SIMGRID_STATIC_LIBRARY_DIRS "${REQUIRED_LIBDIRS}")
      set(SIMGRID_STATIC_INCLUDE_DIRS "${REQUIRED_INCDIRS}")
      set(SIMGRID_STATIC_CFLAGS_OTHER "${REQUIRED_FLAGS}")
      set(SIMGRID_STATIC_LDFLAGS_OTHER "${REQUIRED_LDFLAGS}")
    endif()
  else()
    if(NOT SIMGRID_FIND_QUIETLY)
      message(STATUS "Looking for simgrid : test of fut_keychange with simgrid library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

endif(SIMGRID_FOUND AND SIMGRID_LIBRARIES)

# check that SIMGRID has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SIMGRID DEFAULT_MSG
  SIMGRID_LIBRARIES
  SIMGRID_WORKS)

  # Add imported targe
if (SIMGRID_FOUND)
  morse_create_imported_target(SIMGRID)
endif()
