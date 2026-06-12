###
#
# @copyright (c) 2012-2020 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
# Copyright 2020 Florent Pruvost
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
# - Find M (Math library)
# Use this module by invoking find_package with the form:
#  find_package(M
#               [REQUIRED]) # Fail with error if m is not found
#
# The following variables are set if found:
#
#   M_LIBRARY gives the path to the standard library "m"
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::M``
#   The library to use for M, if found.
#
#=============================================================================

include(FindPackageHandleStandardArgs)

# tests used in this script is not compliant with -Werror or -Werror=...
# remove it temporarily from C flags
set(CMAKE_C_FLAGS_COPY "${CMAKE_C_FLAGS}" CACHE STRING "" FORCE)
string(REGEX REPLACE "-Werror[^ ]*" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

# check if we can call math directly without linking explicitly to libm
include(CheckFunctionExists)
check_function_exists(sqrt HAVE_MATH)

# if works with the compiler we do not need anything else, variables are empty
if(HAVE_MATH)

  set(M_LIBRARY)

  if(NOT TARGET MORSE::M)
    add_library(MORSE::M INTERFACE IMPORTED)
  endif()

  find_package_handle_standard_args(M DEFAULT_MSG)

else()

  # look for libm
  find_library(M_LIBRARY m)

  # check call to math
  set(CMAKE_REQUIRED_LIBRARIES ${M_LIBRARY})
  check_function_exists(sqrt LIBM_MATH_WORKS)
  unset(CMAKE_REQUIRED_LIBRARIES)

  # check and set M_FOUND
  find_package_handle_standard_args(M DEFAULT_MSG M_LIBRARY LIBM_MATH_WORKS)
  mark_as_advanced(M_LIBRARY LIBM_MATH_WORKS)

  # add imported target
  if(M_FOUND)
    if(NOT TARGET MORSE::M)
      add_library(MORSE::M INTERFACE IMPORTED)
      set_target_properties(MORSE::M PROPERTIES
        INTERFACE_LINK_LIBRARIES "${M_LIBRARY}")
    endif()
  endif()

endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_COPY}" CACHE STRING "" FORCE)
