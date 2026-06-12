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
# - Find GTG include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(GTG
#               [REQUIRED]) # Fail with error if gtg is not found
#
# This module finds headers and gtg library using pkg-config file.
#  if found with pkg-config the following variables are set
#  <PREFIX>  = GTG
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
# Set GTG_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::GTG``
#   The headers and libraries to use for GTG, if found.
#
#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

if (NOT GTG_FIND_QUIETLY)
  message(STATUS "FindGTG needs pkg-config program and PKG_CONFIG_PATH set with gtg.pc file path.")
endif()

# use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
if (PKG_CONFIG_EXECUTABLE)
  unset(GTG_FOUND CACHE)
  pkg_search_module(GTG gtg)

  if (NOT GTG_FIND_QUIETLY)
    if (GTG_FOUND AND GTG_LIBRARIES)
      message(STATUS "Looking for GTG - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for GTG - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing gtg.pc to the"
        "\n   PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()
  if (GTG_FOUND AND GTG_LIBRARIES)
    if (NOT GTG_INCLUDE_DIRS)
      pkg_get_variable(GTG_INCLUDE_DIRS gtg includedir)
    endif()
    set(GTG_FOUND_WITH_PKGCONFIG "TRUE")
    morse_find_pkgconfig_libraries_absolute_path(GTG)
  else()
    set(GTG_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

  if (GTG_STATIC AND GTG_STATIC_LIBRARIES)
    set (GTG_DEPENDENCIES ${GTG_STATIC_LIBRARIES})
    list (REMOVE_ITEM GTG_DEPENDENCIES "gtg")
    list (APPEND GTG_LIBRARIES ${GTG_DEPENDENCIES})
    set(GTG_CFLAGS_OTHER ${GTG_STATIC_CFLAGS_OTHER})
    set(GTG_LDFLAGS_OTHER ${GTG_STATIC_LDFLAGS_OTHER})
    if (NOT GTG_FIND_QUIETLY)
      message(STATUS "GTG_STATIC set to 1 by user, GTG_LIBRARIES: ${GTG_LIBRARIES}.")
    endif()
  endif()
endif()

# check a function to validate the find
if(GTG_FOUND AND GTG_LIBRARIES)

  # check if static or dynamic lib
  morse_check_static_or_dynamic(GTG GTG_LIBRARIES)
  if(GTG_STATIC)
    set(STATIC "_STATIC")
  else()
    set(STATIC "")
   endif()

  # set required libraries for link
  morse_set_required_test_lib_link(GTG)

  # test link
  unset(GTG_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(initTrace GTG_WORKS)
  mark_as_advanced(GTG_WORKS)

  if(NOT GTG_WORKS)
    if(NOT GTG_FIND_QUIETLY)
      message(STATUS "Looking for gtg : test of GTG_start with gtg library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

endif(GTG_FOUND AND GTG_LIBRARIES)

# check that GTG has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GTG DEFAULT_MSG
  GTG_LIBRARIES
  GTG_WORKS)

# Add imported targe
if (GTG_FOUND)
  morse_create_imported_target(GTG)
endif()