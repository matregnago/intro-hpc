###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2019 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
#  @file MorseInit.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 1.0.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 13-07-2012
#
###

# Path to Morse modules
get_filename_component(MORSE_CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_FILE} DIRECTORY CACHE)

# Global Morse options
option(MORSE_ENABLE_WARNING         "Enable warning messages from the compiler" OFF)
option(MORSE_ENABLE_STATIC_ANALYSIS "Enable warning messages from external static analyzers such as cppcheck " OFF)
option(MORSE_ENABLE_COVERAGE        "Enable flags for coverage test" OFF)
option(MORSE_ENABLE_COLOR_MESSAGE   "Enable colors in messages" OFF)
option(MORSE_ENABLE_ENV_INCLUDE     "Enable the use of INCLUDE, INCLUDE_PATH to populate CMAKE_INCLUDE_PATH in that order" ON)
option(MORSE_ENABLE_ENV_LIBRARY     "Enable the use of LIBRARY_PATH, LD_LIBRARY_PATH to populate CMAKE_LIBRARY_PATH in that order" ON)
#option(MORSE_VERBOSE_FIND_PACKAGE "Add additional messages concerning packages not found" OFF)
#message(STATUS "MORSE_VERBOSE_FIND_PACKAGE is set to OFF, turn it ON to get"
#        "   information about packages not found")


# This include is required to check symbols of libs in the main CMakeLists.txt
include(CheckFunctionExists)

# This include is required to check defines in headers
include(CheckIncludeFiles)

if (MORSE_ENABLE_COLOR_MESSAGE)
  # colorize messages
  include(ColorizeMessage)
endif()

# Define some auxilary flags
include(AuxilaryFlags)

# Define some variables to set info about ressources
include(Ressources)

# Add the path where we handle our FindFOO.cmake to seek for liraries
list(APPEND CMAKE_MODULE_PATH ${MORSE_CMAKE_MODULE_PATH}/find)

# To load some macros used in Finds (could be useful for other projects)
include(FindMorseInit)

### Build type
set( CMAKE_BUILD_TYPE_DROP_LIST "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
include(Sanitizer)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are None, Debug, Release, RelWithDebInfo, MinSizeRel, ..." FORCE)
endif(NOT CMAKE_BUILD_TYPE)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${CMAKE_BUILD_TYPE_DROP_LIST})

if ( MORSE_ENABLE_ENV_LIBRARY )
  set(_libdir "$ENV{LIBRARY_PATH}")
  if (WIN32)
    set(_libdir "${_libdir}:$ENV{LIB}")
  elseif (APPLE)
    set(_libdir "${_libdir}:$ENV{DYLD_LIBRARY_PATH}")
  else ()
    set(_libdir "${_libdir}:$ENV{LD_LIBRARY_PATH}")
  endif ()
  string(REPLACE ":" ";" _libdir "${_libdir}")
  string(REPLACE ";;" ";" _libdir "${_libdir}")

  list(APPEND _libdir "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
  list(REMOVE_DUPLICATES _libdir)

  list(APPEND CMAKE_LIBRARY_PATH "${_libdir}" )
endif()

if ( MORSE_ENABLE_ENV_INCLUDE )
  if(WIN32)
    set( _incdir "$ENV{INCLUDE}" )
  else()
    set( _incdir "$ENV{INCLUDE}:" )
    set( _incdir "${_incdir}:$ENV{C_INCLUDE_PATH}" )
    set( _incdir "${_incdir}:$ENV{CPATH}" )
    set( _incdir "${_incdir}:$ENV{INCLUDE_PATH}" )
  endif()

  string( REPLACE ":"  ";" _incdir "${_incdir}" )
  string( REPLACE ";;" ";" _incdir "${_incdir}" )

  list(APPEND _incdir "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
  list(REMOVE_DUPLICATES _incdir)

  list(APPEND CMAKE_INCLUDE_PATH "${_incdir}" )
endif()

##
## @end file MorseInit.cmake
##
