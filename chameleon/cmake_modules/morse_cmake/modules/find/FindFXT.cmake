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
# - Find FXT include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(FXT
#               [REQUIRED]) # Fail with error if fxt is not found
#
# This module finds headers and FXT library using pkg-config file.
#
#  if found with pkg-config the following variables are set
#  <PREFIX>  = FXT
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
# Set FXT_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::FXT``
#   The headers and libraries to use for FXT, if found.
#
#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

if (NOT FXT_FIND_QUIETLY)
  message(STATUS "FindFXT needs pkg-config program and PKG_CONFIG_PATH set with fxt.pc file path.")
endif()

# use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
if (PKG_CONFIG_EXECUTABLE)
  unset(FXT_FOUND CACHE)
  pkg_search_module(FXT fxt)

  message(STATUS "FXT_LIBRARY_DIRS ${FXT_LIBRARY_DIRS}")
  message(STATUS "FXT_LIBDIR ${FXT_LIBDIR}")

  if (NOT FXT_FIND_QUIETLY)
    if (FXT_FOUND AND FXT_LIBRARIES)
      message(STATUS "Looking for FXT - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for FXT - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing fxt.pc to the"
        "\n   PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()
  if (FXT_FOUND AND FXT_LIBRARIES)
    if (NOT FXT_INCLUDE_DIRS)
      pkg_get_variable(FXT_INCLUDE_DIRS fxt includedir)
    endif()
    set(FXT_FOUND_WITH_PKGCONFIG "TRUE")
    morse_find_pkgconfig_libraries_absolute_path(FXT)
  else()
    set(FXT_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

  if (FXT_STATIC AND FXT_STATIC_LIBRARIES)
    set (FXT_DEPENDENCIES ${FXT_STATIC_LIBRARIES})
    list (REMOVE_ITEM FXT_DEPENDENCIES "fxt")
    list (APPEND FXT_LIBRARIES ${FXT_DEPENDENCIES})
    set(FXT_CFLAGS_OTHER ${FXT_STATIC_CFLAGS_OTHER})
    set(FXT_LDFLAGS_OTHER ${FXT_STATIC_LDFLAGS_OTHER})
    if (NOT FXT_FIND_QUIETLY)
      message(STATUS "FXT_STATIC set to 1 by user, FXT_LIBRARIES: ${FXT_LIBRARIES}.")
    endif()
  endif()
endif()


# check a function to validate the find
if(FXT_FOUND AND FXT_LIBRARIES)

  # check if static or dynamic lib
  morse_check_static_or_dynamic(FXT FXT_LIBRARIES)
  if(FXT_STATIC)
    set(STATIC "_STATIC")
  else()
    set(STATIC "")
  endif()

  # set required libraries for link
  morse_set_required_test_lib_link(FXT)

  # test link
  unset(FXT_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(fut_keychange FXT_WORKS)
  mark_as_advanced(FXT_WORKS)

  if(NOT FXT_WORKS)
    if(NOT FXT_FIND_QUIETLY)
      message(STATUS "Looking for fxt : test of fut_keychange with fxt library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

endif(FXT_FOUND AND FXT_LIBRARIES)

# check that FXT has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FXT DEFAULT_MSG
  FXT_LIBRARIES
  FXT_WORKS)

# Add imported targe
if (FXT_FOUND)
  morse_create_imported_target(FXT)
endif()