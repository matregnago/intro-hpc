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
###
#
# - Find SUITESPARSE include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(SUITESPARSE
#               [REQUIRED] # Fail with error if suitesparse is not found
#              )
#
#  SUITESPARSE depends on the following libraries:
#   - BLAS
#   - LAPACK
#   - METIS
##
# This module finds headers and suitesparse library.
# Results are reported in variables:
#  SUITESPARSE_FOUND             - True if headers and requested libraries were found
#  SUITESPARSE_PREFIX            - installation path of the lib found
#  SUITESPARSE_CFLAGS_OTHER      - suitesparse compiler flags without headers paths
#  SUITESPARSE_LDFLAGS_OTHER     - suitesparse linker flags without libraries
#  SUITESPARSE_INCLUDE_DIRS      - suitesparse include directories
#  SUITESPARSE_LIBRARY_DIRS      - suitesparse link directories
#  SUITESPARSE_LIBRARIES         - suitesparse libraries to be linked (absolute path)
#
# Set SUITESPARSE_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::SUITESPARSE``
#   The headers and libraries to use for SUITESPARSE, if found.
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DSUITESPARSE_DIR=path/to/suitesparse):
#  SUITESPARSE_DIR              - Where to find the base directory of suitesparse
# The module can also look for the following environment variables if paths
# are not given as cmake variable: SUITESPARSE_DIR

#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

# Set variables from environment if needed
# ----------------------------------------
morse_find_package_get_envdir(SUITESPARSE)

if (NOT SUITESPARSE_FIND_QUIETLY)
  message(STATUS "Looking for SUITESPARSE")
endif()

if (NOT SUITESPARSE_FIND_QUIETLY)
  message(STATUS "Looking for SUITESPARSE - PkgConfig not used")
endif()

# Required dependencies
# ---------------------
if (NOT SUITESPARSE_FIND_QUIETLY)
  message(STATUS "Looking for SUITESPARSE - Try to detect metis")
endif()
if (SUITESPARSE_FIND_REQUIRED)
  find_package(METIS REQUIRED)
else()
  find_package(METIS)
endif()

# SUITESPARSE depends on LAPACK
#------------------------------
if (NOT SUITESPARSE_FIND_QUIETLY)
  message(STATUS "Looking for SUITESPARSE - Try to detect LAPACK")
endif()
if (SUITESPARSE_FIND_REQUIRED)
  find_package(LAPACK REQUIRED)
else()
  find_package(LAPACK)
endif()

# Looking for SUITESPARSE
# -----------------

# Looking for include
# -------------------
set(SUITESPARSE_hdrs_to_find
  amd.h
  btf.h
  ccolamd.h
  colamd.h
  cs.h
  klu.h
  ldl.h
  #RBio.h
  spqr.hpp
  SuiteSparse_config.h
  umfpack.h)

morse_find_path(SUITESPARSE
  HEADERS  ${SUITESPARSE_hdrs_to_find}
  SUFFIXES include include/suitesparse
  OPTIONAL )

# Make sure we found at least SuiteSparse_config.h
# ------------------------------------------------
if (NOT SUITESPARSE_SuiteSparse_config.h_DIRS)
  set(SUITESPARSE_INCLUDE_DIRS "SUITESPARSE_INCLUDE_DIRS-NOTFOUND")
  if (NOT SUITESPARSE_FIND_QUIETLY)
    message(STATUS "Looking for suitesparse -- SuiteSparse_config.h not found")
  endif()
endif()

# Looking for lib
# ---------------
set(SUITESPARSE_libs_to_find
  cholmod
  cxsparse
  klu
  ldl
  spqr
  umfpack
  amd
  btf
  camd
  ccolamd
  colamd
  #rbio
  suitesparseconfig
  )

morse_find_library(SUITESPARSE
  LIBRARIES ${SUITESPARSE_libs_to_find}
  SUFFIXES lib lib32 lib64
  OPTIONAL)

# Make sure we found at least suitesparseconfig
# ---------------------------------------------
if (NOT SUITESPARSE_suitesparseconfig_LIBRARY)
  set(SUITESPARSE_LIBRARIES "SUITESPARSE_LIBRARIES-NOTFOUND")
  if(NOT SUITESPARSE_FIND_QUIETLY)
    message(STATUS "Looking for suitesparse -- libsuitesparseconfig.a/so not found")
  endif()
endif()

# check a function to validate the find
if(SUITESPARSE_LIBRARIES)

  # check if static or dynamic lib
  morse_check_static_or_dynamic(SUITESPARSE SUITESPARSE_LIBRARIES)
  if(SUITESPARSE_STATIC)
    set(STATIC "_STATIC")
  endif()

  set(REQUIRED_INCDIRS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)

  # SUITESPARSE
  if (SUITESPARSE_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${SUITESPARSE_INCLUDE_DIRS}")
  endif()
  foreach(libdir ${SUITESPARSE_LIBRARY_DIRS})
    if (libdir)
      list(APPEND REQUIRED_LIBDIRS "${libdir}")
    endif()
  endforeach()
  set(REQUIRED_LIBS "${SUITESPARSE_LIBRARIES}")
  # METIS
  if (METIS_FOUND)
    if (METIS_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${METIS_INCLUDE_DIRS}")
    endif()
    foreach(libdir ${METIS_LIBRARY_DIRS})
      if (libdir)
        list(APPEND REQUIRED_LIBDIRS "${libdir}")
      endif()
    endforeach()
    list(APPEND REQUIRED_LIBS "${METIS_LIBRARIES}")
  endif()
  # LAPACK
  if (LAPACK_FOUND)
    if (LAPACK_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${LAPACK_INCLUDE_DIRS}")
    endif()
    if (LAPACK_CFLAGS_OTHER)
      list(APPEND REQUIRED_FLAGS "${LAPACK_CFLAGS_OTHER}")
    endif()
    if (LAPACK_LDFLAGS_OTHER)
      list(APPEND REQUIRED_LDFLAGS "${LAPACK_LDFLAGS_OTHER}")
    endif()
    if (LAPACK_LIBRARY_DIRS)
      list(APPEND REQUIRED_LIBDIRS "${LAPACK_LIBRARY_DIRS}")
    endif()
    list(APPEND REQUIRED_LIBS "${LAPACK_LIBRARIES}")
  endif()
  # others
  set(M_LIBRARY "M_LIBRARY-NOTFOUND")
  find_library(M_LIBRARY NAMES m)
  mark_as_advanced(M_LIBRARY)
  if(M_LIBRARY)
    list(APPEND REQUIRED_LIBS "${M_LIBRARY}")
  endif()

  # set required libraries for link
  set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
  if (REQUIRED_FLAGS)
    set(REQUIRED_FLAGS_COPY "${REQUIRED_FLAGS}")
    set(REQUIRED_FLAGS)
    set(REQUIRED_DEFINITIONS)
    foreach(_flag ${REQUIRED_FLAGS_COPY})
      if (_flag MATCHES "^-D")
       list(APPEND REQUIRED_DEFINITIONS "${_flag}")
      endif()
      string(REGEX REPLACE "^-D.*" "" _flag "${_flag}")
      list(APPEND REQUIRED_FLAGS "${_flag}")
    endforeach()
  endif()
  morse_finds_remove_duplicates()
  set(CMAKE_REQUIRED_DEFINITIONS "${REQUIRED_DEFINITIONS}")
  set(CMAKE_REQUIRED_FLAGS "${REQUIRED_LDFLAGS}")
  set(CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
  string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

  # test link
  include(CheckFortranFunctionExists)
  unset(SUITESPARSE_WORKS CACHE)
  check_function_exists(SuiteSparse_start SUITESPARSE_WORKS)
  mark_as_advanced(SUITESPARSE_WORKS)

  if(SUITESPARSE_WORKS)
    set(SUITESPARSE_LIBRARY_DIRS "${REQUIRED_LIBDIRS}")
    set(SUITESPARSE_INCLUDE_DIRS "${REQUIRED_INCDIRS}")
    set(SUITESPARSE_CFLAGS_OTHER "${REQUIRED_FLAGS}")
    set(SUITESPARSE_LDFLAGS_OTHER "${REQUIRED_LDFLAGS}")
    if (SUITESPARSE_STATIC OR METIS_STATIC OR BLA_STATIC)
      # save link with dependencies
      set(SUITESPARSE_LIBRARIES "${REQUIRED_LIBS}")
    endif()
  else()
    if(NOT SUITESPARSE_FIND_QUIETLY)
      message(STATUS "Looking for SUITESPARSE : test of symbol SuiteSparse_start fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
      message(STATUS "Maybe SUITESPARSE is linked with specific libraries. ")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

  list(GET SUITESPARSE_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" DIRECTORY)
  if (NOT SUITESPARSE_LIBRARY_DIRS)
    set(SUITESPARSE_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(SUITESPARSE_PREFIX "${not_cached_dir}" CACHE PATH "Installation directory of SUITESPARSE library" FORCE)
  else()
    set(SUITESPARSE_PREFIX "${first_lib_path}" CACHE PATH "Installation directory of SUITESPARSE library" FORCE)
  endif()
  mark_as_advanced(SUITESPARSE_DIR)
  mark_as_advanced(SUITESPARSE_PREFIX)

endif(SUITESPARSE_LIBRARIES)

# check that SUITESPARSE has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUITESPARSE DEFAULT_MSG
  SUITESPARSE_LIBRARIES
  SUITESPARSE_WORKS)

# Add imported target
if (SUITESPARSE_FOUND)
  morse_create_imported_target(SUITESPARSE)
endif()
