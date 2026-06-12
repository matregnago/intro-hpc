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
# - Find PAPI include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(PAPI
#               [REQUIRED]) # Fail with error if papi is not found
#
# This module finds headers and papi library using pkg-config file.
#  if found with pkg-config the following variables are set
#  <PREFIX>  = PAPI
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
# Set PAPI_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::PAPI``
#   The headers and libraries to use for PAPI, if found.
#
#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

if (NOT PAPI_FIND_QUIETLY)
  message(STATUS "FindPAPI needs pkg-config program and PKG_CONFIG_PATH set with papi.pc file path.")
endif()

# use pkg-config to detect include/library dirs (if pkg-config is available)
# --------------------------------------------------------------------------
if (PKG_CONFIG_EXECUTABLE)
  unset(PAPI_FOUND CACHE)
  pkg_search_module(PAPI papi)

  if (NOT PAPI_FIND_QUIETLY)
    if (PAPI_FOUND AND PAPI_LIBRARIES)
      message(STATUS "Looking for PAPI - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for PAPI - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing papi.pc to the"
        "\n   PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()
  if (PAPI_FOUND AND PAPI_LIBRARIES)
    if (NOT PAPI_INCLUDE_DIRS)
      pkg_get_variable(PAPI_INCLUDE_DIRS papi includedir)
    endif()
    set(PAPI_FOUND_WITH_PKGCONFIG "TRUE")
    morse_find_pkgconfig_libraries_absolute_path(PAPI)
  else()
    set(PAPI_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

  if (PAPI_STATIC AND PAPI_STATIC_LIBRARIES)
    set (PAPI_DEPENDENCIES ${PAPI_STATIC_LIBRARIES})
    list (REMOVE_ITEM PAPI_DEPENDENCIES "papi")
    list (APPEND PAPI_LIBRARIES ${PAPI_DEPENDENCIES})
    set(PAPI_CFLAGS_OTHER ${PAPI_STATIC_CFLAGS_OTHER})
    set(PAPI_LDFLAGS_OTHER ${PAPI_STATIC_LDFLAGS_OTHER})
    if (NOT PAPI_FIND_QUIETLY)
      message(STATUS "PAPI_STATIC set to 1 by user, PAPI_LIBRARIES: ${PAPI_LIBRARIES}.")
    endif()
  endif()
endif()

# check a function to validate the find
if(PAPI_FOUND AND PAPI_LIBRARIES)

  # check if static or dynamic lib
  morse_check_static_or_dynamic(PAPI PAPI_LIBRARIES)
  if(PAPI_STATIC)
    set(STATIC "_STATIC")
  else()
    set(STATIC "")
  endif()

  # set required libraries for link
  morse_set_required_test_lib_link(PAPI)

  # test link
  unset(PAPI_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(PAPI_start PAPI_WORKS)
  mark_as_advanced(PAPI_WORKS)

  if(PAPI_WORKS)
    if (PAPI_STATIC)
      # save link with dependencies
      set(PAPI_STATIC_LIBRARIES "${REQUIRED_LIBS}")
      set(PAPI_STATIC_LIBRARY_DIRS "${REQUIRED_LIBDIRS}")
      set(PAPI_STATIC_INCLUDE_DIRS "${REQUIRED_INCDIRS}")
      set(PAPI_STATIC_CFLAGS_OTHER "${REQUIRED_FLAGS}")
      set(PAPI_STATIC_LDFLAGS_OTHER "${REQUIRED_LDFLAGS}")
    endif()
  else()
    if(NOT PAPI_FIND_QUIETLY)
      message(STATUS "Looking for papi : test of PAPI_start with papi library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

endif(PAPI_FOUND AND PAPI_LIBRARIES)

# check that PAPI has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PAPI DEFAULT_MSG
  PAPI_LIBRARIES
  PAPI_WORKS)

  # Add imported targe
if (PAPI_FOUND)
  morse_create_imported_target(PAPI)
endif()
