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
# - Find RT (librt, libposix4 - POSIX.1b Realtime Extensions library)
# Use this module by invoking find_package with the form:
#  find_package(RT
#               [REQUIRED]) # Fail with error if rt is not found
#
# The following variables are set if found:
#
#   RT_LIBRARIES gives the library librt
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::RT``
#   The headers and libraries to use for RT, if found.
#
#=============================================================================

find_library(RT_LIBRARY rt)
set(RT_LIBRARIES ${RT_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(RT DEFAULT_MSG RT_LIBRARY)
mark_as_advanced(RT_LIBRARY)

# add imported target
if(RT_FOUND)
  if(NOT TARGET MORSE::RT)
    add_library(MORSE::RT INTERFACE IMPORTED)
    set_target_properties(MORSE::RT PROPERTIES
      INTERFACE_LINK_LIBRARIES "${RT_LIBRARIES}")
  endif()
endif()
