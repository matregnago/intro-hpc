###
#
# @copyright (c) 2019 Inria. All rights reserved.
#
###
#
# - Find CPPCHECK bin directory
# Use this module by invoking find_package with the form:
#  find_package(CPPCHECK
#               [REQUIRED]) # Fail with error if cppcheck is not found
#
# This module finds headers and cppcheck library.
# Results are reported in variables:
#  CPPCHECK_FOUND           - True if cppcheck executable is found
#  CPPCHECK_EXECUTABLE      - cppcheck executable directory found
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DCPPCHECK_DIR=path/to/cppcheck):
#  CPPCHECK_DIR             - Where to find the base directory of cppcheck
#  CPPCHECK_BINDIR          - Where to find the executable cppcheck
# The module can also look for the following environment variables if paths
# are not given as cmake variable: CPPCHECK_DIR, CPPCHECK_BINDIR

# Common macros to use in finds
include(FindMorseInit)

# Set variables from environment if needed
# ----------------------------------------
morse_find_package_get_envdir(CPPCHECK)

set(CPPCHECK_GIVEN_BY_USER "FALSE")
if ( CPPCHECK_DIR OR CPPCHECK_BINDIR )
  set(CPPCHECK_GIVEN_BY_USER "TRUE")
endif()

if (NOT CPPCHECK_FIND_QUIETLY)
  message(STATUS "Looking for CPPCHECK")
endif()

# Looking for executable
# ----------------------

# Add system include paths to search include
# ------------------------------------------
unset(where_to_look)
if(CPPCHECK_BINDIR)
  list(APPEND where_to_look "${CPPCHECK_BINDIR}")
elseif(CPPCHECK_DIR)
  list(APPEND where_to_look "${CPPCHECK_DIR}")
else()
  string(REPLACE ":" ";" _path_env "$ENV{PATH}")
  list(APPEND where_to_look "${_path_env}")
endif()
list(REMOVE_DUPLICATES where_to_look)

find_program(CPPCHECK_EXECUTABLE
  NAMES cppcheck
  HINTS ${where_to_look}
)
mark_as_advanced(CPPCHECK_EXECUTABLE)

if(CPPCHECK_EXECUTABLE)
  execute_process(COMMAND ${CPPCHECK_EXECUTABLE} "--version" OUTPUT_VARIABLE CPPCHECK_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX REPLACE "Cppcheck " "" CPPCHECK_VERSION "${CPPCHECK_VERSION}")
endif()

# check that CPPCHECK has been found
# ----------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CPPCHECK DEFAULT_MSG
                                  CPPCHECK_EXECUTABLE
                                  CPPCHECK_VERSION)
