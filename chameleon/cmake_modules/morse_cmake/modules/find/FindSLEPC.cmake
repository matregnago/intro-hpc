#
# @copyright (c) 2018-2020 Inria. All rights reserved.
#
# Marc Fuentes
# Vincent Perrier
# Florent Pruvost
#
# https://gitlab.inria.fr/sed-bso/findPetsc
#
# INPUT
# -----
#
# SLEPC is not installed in a standard way on Unix systems so that
# this module requires hints to know where PETSC is installed. Please
# give the installation directory (contains ./include/petsc.h, ./lib/, etc):
# 1. by setting the SLEPC_DIR variable
#    a. as an environment variable, e.g.
#       $ export SLEPC_DIR=/usr/lib/slepcdir/3.6.2/x86_64-linux-gnu-real
#    b. or as an CMake variable, e.g.
#       $ cmake .. -DSLEPC_DIR==/usr/lib/slepcdir/3.6.2/x86_64-linux-gnu-real
# 2. or by using the pkg-config mechanism, e.g.
#    $ export PKG_CONFIG_PATH=/usr/lib/slepcdir/3.6.2/x86_64-linux-gnu-real/lib/pkgconfig:$PKG_CONFIG_PATH
#
# OUTPUT
# -------
# SLEPC_INCLUDE_DIRS - the PETSC include directories
# SLEPC_LIBRARIES    - Link these to use PETSC
# SLEPC_LIBRARY_DIRS - Link these to use PETSC
# PETSC_MPIEXEC - Executable for running MPI programs ????
#
# if pkg-config is used i.e. pkgconfig installed, PETSC.pc file path
# in the PKG_CONFIG_PATH environment variable and SLEPC_DIR not set
# then the following variables are set (or empty)
#
# SLEPC_VERSION           ... version of the module
# SLEPC_PREFIX            ... prefix-directory of the module
# SLEPC_INCLUDEDIR        ... include-dir of the module
# SLEPC_LIBDIR            ... lib-dir of the module
#
# <PREFIX>_FOUND          ... set to 1 if module(s) exist
# <PREFIX>_LIBRARIES      ... only the libraries (w/o the '-l')
# <PREFIX>_LIBRARY_DIRS   ... the paths of the libraries (w/o the '-L')
# <PREFIX>_LDFLAGS        ... all required linker flags
# <PREFIX>_LDFLAGS_OTHER  ... all other linker flags
# <PREFIX>_INCLUDE_DIRS   ... the '-I' preprocessor flags (w/o the '-I')
# <PREFIX>_CFLAGS         ... all required cflags
# <PREFIX>_CFLAGS_OTHER   ... the other compiler flags
#
# <PREFIX> = SLEPC        for common case
# <PREFIX> = SLEPC_STATIC for static linking
#
# SLEPC_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::SLEPC``
#   The headers and libraries to use for SLEPC, if found.
#
# find_package(SLEPC [QUIET] [REQUIRED])
#
# Setting these changes the behavior of the search
# SLEPC_DIR - directory in which SLEPC is installed

# create a cmake cache variable
set(SLEPC_DIR "" CACHE PATH "Installation directory of SLEPC library")
if (NOT SLEPC_FIND_QUIETLY AND NOT SLEPC_DIR)
  message(STATUS "A cache variable, namely SLEPC_DIR, has been set
  to specify a custom installation directory of SLEPC")
endif()

# Use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
if( PKG_CONFIG_EXECUTABLE AND NOT SLEPC_DIR )
  pkg_search_module(SLEPC SLEPc)
  if (NOT SLEPC_FIND_QUIETLY)
    if (SLEPC_FOUND AND SLEPC_LIBRARIES)
      message(STATUS "Looking for SLEPC - found using PkgConfig")
    else()
      message(STATUS "Looking for SLEPC - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing SLEPC.pc to"
        "\n   the PKG_CONFIG_PATH environment variable.")
    endif()
  endif()
  set(SLEPC_DIR "${SLEPC_PREFIX}")
  if (SLEPC_FOUND AND SLEPC_LIBRARIES)
    if (NOT SLEPC_INCLUDE_DIRS)
      pkg_get_variable(SLEPC_INCLUDE_DIRS SLEPc includedir)
    endif()
    set(SLEPC_FOUND_WITH_PKGCONFIG "TRUE")
    morse_find_pkgconfig_libraries_absolute_path(SLEPC)
  else()
    set(SLEPC_FOUND_WITH_PKGCONFIG "FALSE")
  endif()
  if (SLEPC_STATIC AND SLEPC_STATIC_LIBRARIES)
    set (SLEPC_DEPENDENCIES ${SLEPC_STATIC_LIBRARIES})
    list (REMOVE_ITEM SLEPC_DEPENDENCIES "slepc")
    list (APPEND SLEPC_LIBRARIES ${SLEPC_DEPENDENCIES})
    if (NOT SLEPC_FIND_QUIETLY)
      message(STATUS "SLEPC_STATIC set to 1 by user, SLEPC_LIBRARIES: ${SLEPC_LIBRARIES}.")
    endif()
  endif()
endif()

# consider using the env. var. SLEPC_DIR if not directly given through the CMake cache var.
if (NOT SLEPC_DIR AND DEFINED ENV{SLEPC_DIR})
  set(SLEPC_DIR "$ENV{SLEPC_DIR}")
endif()

if (SLEPC_DIR)
  if (EXISTS ${SLEPC_DIR})
    if (EXISTS ${SLEPC_DIR}/include/slepc.h)
        if (NOT SLEPC_FIND_QUIETLY)
          message(STATUS "SLEPC_DIR = ${SLEPC_DIR} contains include/slepc.h")
        endif()
        set(SLEPC_INCLUDE_DIRS "${SLEPC_DIR}/include")
    else()
      if (SLEPC_FIND_REQUIRED)
        message(FATAL_ERROR "include/slepc.h not found in SLEPC_DIR = ${SLEPC_DIR}")
      endif()
    endif()
    if (SLEPC_STATIC)
      if(EXISTS ${SLEPC_DIR}/lib/libslepc.a)
         set(SLEPC_LIBRARIES "${SLEPC_DIR}/lib/libslepc.a")
      else()
        if (SLEPC_FIND_REQUIRED)
          message(FATAL_ERROR "lib/libslepc.a not found in SLEPC_DIR = ${SLEPC_DIR}")
        endif()
      endif()
    else()
      if (EXISTS ${SLEPC_DIR}/lib/libslepc.so)
        set(SLEPC_LIBRARIES "${SLEPC_DIR}/lib/libslepc.so")
      else()
        if (SLEPC_FIND_REQUIRED)
          message(FATAL_ERROR "lib/libslepc.so not found in SLEPC_DIR = ${SLEPC_DIR}")
        endif()
      endif()
    endif()
    if (SLEPC_DIR AND SLEPC_LIBRARIES)
      set(SLEPC_PREFIX ${SLEPC_DIR})
    endif()
  else()
    if (SLEPC_FIND_REQUIRED)
      message(FATAL_ERROR "SLEPC_DIR defined, but ${SLEPC_DIR} does not exist")
    endif()
  endif()
else()
  if (SLEPC_FIND_REQUIRED)
    message(FATAL_ERROR "\
SLEPC is not installed in a standard way on Unix systems so that
this module requires hints to know where SLEPC is installed. Please
give the installation directory (contains ./include/slepc.h, ./lib/, etc):
1. by setting the SLEPC_DIR variable
   a. as an environment variable, e.g.
      $ export SLEPC_DIR=/usr/lib/slepcdir/3.6.2/x86_64-linux-gnu-real
   b. or as an CMake variable, e.g.
      $ cmake .. -DSLEPC_DIR==/usr/lib/slepcdir/3.6.2/x86_64-linux-gnu-real
2. or by using the pkg-config mechanism, e.g.
   $ export PKG_CONFIG_PATH=/usr/lib/slepcdir/3.6.2/x86_64-linux-gnu-real/lib/pkgconfig:$PKG_CONFIG_PATH\
    ")
  endif()
endif()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (SLEPC
  "SLEPC could not be found. Be sure to set SLEPC_DIR or use with pkg-config (slepc.pc in PKG_CONFIG_PATH)."
  SLEPC_INCLUDE_DIRS SLEPC_LIBRARIES)

# Add imported target
if (SLEPC_FOUND)
  morse_create_imported_target(SLEPC)
endif()