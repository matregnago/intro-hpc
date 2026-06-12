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
# - Find LITL include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(LITL
#               [REQUIRED]) # Fail with error if litl is not found
#
#  LITL depends on the following libraries:
#   - threads (e.g. pthreads)
#
# This module finds headers and litl library using pkg-config file.
#  if found with pkg-config the following variables are set
#  <PREFIX>  = LITL
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
# Set LITL_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::LITL``
#   The headers and libraries to use for LITL, if found.
#
#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

if (NOT LITL_FIND_QUIETLY)
  message(STATUS "FindLITL needs pkg-config program and PKG_CONFIG_PATH set with litl.pc file path.")
endif()

# Use pkg-config to detect include/library dirs
# ---------------------------------------------
if (PKG_CONFIG_EXECUTABLE)
  unset(LITL_FOUND CACHE)
  pkg_search_module(LITL litl)

  if (NOT LITL_FIND_QUIETLY)
    if (LITL_FOUND AND LITL_LIBRARIES)
      message(STATUS "Looking for LITL - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for LITL - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing litl.pc to"
        "\n   the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()

  if (LITL_FOUND AND LITL_LIBRARIES)
    # Set include dirs and libraries absolute path.
    if (NOT LITL_INCLUDE_DIRS)
      pkg_get_variable(LITL_INCLUDE_DIRS litl includedir)
    endif()
    set(LITL_FOUND_WITH_PKGCONFIG "TRUE")
    morse_find_pkgconfig_libraries_absolute_path(LITL)
  else()
      set(LITL_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

  if (LITL_STATIC AND LITL_STATIC_LIBRARIES)
    set (LITL_DEPENDENCIES ${LITL_STATIC_LIBRARIES})
    list (REMOVE_ITEM LITL_DEPENDENCIES "litl")
    list (APPEND LITL_LIBRARIES ${LITL_DEPENDENCIES})
    set(LITL_CFLAGS_OTHER ${LITL_STATIC_CFLAGS_OTHER})
    set(LITL_LDFLAGS_OTHER ${LITL_STATIC_LDFLAGS_OTHER})
    if (NOT LITL_FIND_QUIETLY)
      message(STATUS "LITL_STATIC set to 1 by user, LITL_LIBRARIES: ${LITL_LIBRARIES}.")
    endif()
  endif()
endif()

# check a function to validate the find
if(LITL_FOUND AND LITL_LIBRARIES)

  # check if static or dynamic lib
  morse_check_static_or_dynamic(LITL LITL_LIBRARIES)
  if(LITL_STATIC)
    set(STATIC "_STATIC")
  else()
    set(STATIC "")
  endif()

  # set required libraries for link
  morse_set_required_test_lib_link(LITL)

  # test link
  unset(LITL_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(litl_time_initialize LITL_WORKS)
  mark_as_advanced(LITL_WORKS)

  if(NOT LITL_WORKS)
    if(NOT LITL_FIND_QUIETLY)
      message(STATUS "Looking for litl : test of litl_time_initialize with litl library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

endif(LITL_FOUND AND LITL_LIBRARIES)

# check that LITL has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LITL DEFAULT_MSG
  LITL_LIBRARIES
  LITL_WORKS)

# Add imported targe
if (LITL_FOUND)
  morse_create_imported_target(LITL)
endif()