###
#
# @copyright (c) 2012-2020 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
# Copyright 2012-2019 Inria
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
# - Find FFTW Version 3 include dirs and libraries
# Default configuration will find the real double precision fftw library version
# without THREADS|OMP.
# Use this module by invoking find_package with the form:
#  find_package(FFTW
#               [REQUIRED] # Fail with error if fftw is not found
#               [COMPONENTS MKL]
#
#  COMPONENTS can be some of the following:
#   - MKL:     to detect the FFTW from Intel MKL
#   - ESSL:    to detect the FFTW from IBM ESSL
#   - THREADS: to detect the Threads version of FFTW
#   - OMP:     to detect the OpenMP version of FFTW
#   - SIMPLE:  to detect the FFTW simple precision fftw3f
#   - LONG:    to detect the FFTW long double precision fftw3l
#   - QUAD:    to detect the FFTW quadruple precision fftw3q
#
# This module finds headers and fftw library.
# Results are reported in variables:
#  FFTW_FOUND             - True if headers and requested libraries were found
#  FFTW_PREFIX            - installation path of the lib found
#  FFTW_CFLAGS_OTHER      - fftw compiler flags without headers paths
#  FFTW_LDFLAGS_OTHER     - fftw linker flags without libraries
#  FFTW_INCLUDE_DIRS      - fftw include directories
#  FFTW_LIBRARY_DIRS      - fftw link directories
#  FFTW_LIBRARIES         - fftw libraries to be linked (absolute path)
#
#  FFTW_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = FFTW3F or FFTW3 or FFTW3L or FFTW3Q
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
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::FFTW``
#   The headers and libraries to use for FFTW, if found.
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DFFTW_DIR=path/to/fftw):
#  FFTW_DIR             - Where to find the base directory of fftw
#  FFTW_INCDIR          - Where to find the header files
#  FFTW_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: FFTW_DIR, FFTW_INCDIR, FFTW_LIBDIR
# For MKL case and if no paths are given as hints, we will try to use the MKLROOT
# environment variable

#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

# Set the version to find
set(FFTW_LOOK_FOR_MKL OFF)
set(FFTW_LOOK_FOR_ESSL OFF)
set(FFTW_LOOK_FOR_THREADS OFF)
set(FFTW_LOOK_FOR_OMP OFF)
set(FFTW_LOOK_FOR_FFTW_SIMPLE OFF)
set(FFTW_LOOK_FOR_FFTW_LONG OFF)
set(FFTW_LOOK_FOR_FFTW_QUAD OFF)

if( FFTW_FIND_COMPONENTS )
  foreach( component ${FFTW_FIND_COMPONENTS} )
    if (${component} STREQUAL "THREADS")
      # means we look for the Threads version of FFTW
      set(FFTW_LOOK_FOR_THREADS ON)
    endif()
    if (${component} STREQUAL "OMP")
      # means we look for the OpenMP version of FFTW
      set(FFTW_LOOK_FOR_OMP ON)
    endif()
    if (${component} STREQUAL "SIMPLE")
      # means we look for FFTW simple precision (fftw3f)
      set(FFTW_LOOK_FOR_FFTW_SIMPLE ON)
      set(FFTW_LOOK_FOR_FFTW_LONG OFF)
      set(FFTW_LOOK_FOR_FFTW_QUAD OFF)
    endif()
    if (${component} STREQUAL "LONG")
      # means we look for FFTW long double precision (fftw3l)
      set(FFTW_LOOK_FOR_FFTW_SIMPLE OFF)
      set(FFTW_LOOK_FOR_FFTW_LONG ON)
      set(FFTW_LOOK_FOR_FFTW_QUAD OFF)
    endif()
    if (${component} STREQUAL "QUAD")
      # means we look for FFTW quad precision (fftw3q)
      set(FFTW_LOOK_FOR_FFTW_SIMPLE OFF)
      set(FFTW_LOOK_FOR_FFTW_LONG OFF)
      set(FFTW_LOOK_FOR_FFTW_QUAD ON)
    endif()
    if (${component} STREQUAL "MKL")
      # means we look for the Intel MKL version of FFTW
      set(FFTW_LOOK_FOR_MKL ON)
      if (FFTW_LOOK_FOR_FFTW_LONG)
        if (NOT FFTW_FIND_QUIETLY)
          message(WARNING "Looking for FFTW -- long precision functions do not exist in MKL FFTW")
        endif()
        set(FFTW_LOOK_FOR_FFTW_LONG OFF)
      endif()
      if (FFTW_LOOK_FOR_FFTW_QUAD)
        if (NOT FFTW_FIND_QUIETLY)
          message(WARNING "Looking for FFTW -- quadruple functions do not exist in MKL FFTW")
        endif()
        set(FFTW_LOOK_FOR_FFTW_QUAD OFF)
      endif()
    endif()
    if (${component} STREQUAL "ESSL")
      # means we look for the Intel MKL version of FFTW
      set(FFTW_LOOK_FOR_ESSL ON)
      if (FFTW_LOOK_FOR_FFTW_LONG)
        if (NOT FFTW_FIND_QUIETLY)
          message(WARNING "Looking for FFTW -- long precision functions do not exist in FFTW_ESSL")
        endif()
        set(FFTW_LOOK_FOR_FFTW_LONG OFF)
      endif()
      if (FFTW_LOOK_FOR_FFTW_QUAD)
        if (NOT FFTW_FIND_QUIETLY)
          message(WARNING "Looking for FFTW -- quadruple functions do not exist in FFTW_ESSL")
        endif()
        set(FFTW_LOOK_FOR_FFTW_QUAD OFF)
      endif()
      if (FFTW_LOOK_FOR_OMP)
        if (NOT FFTW_FIND_QUIETLY)
          message(WARNING "Looking for FFTW -- FFTW_ESSL does not use OpenMP")
        endif()
        set(FFTW_LOOK_FOR_OMP OFF)
      endif()
    endif()
  endforeach()
endif()

if (FFTW_LOOK_FOR_THREADS)
  if (NOT FFTW_FIND_QUIETLY)
    message(STATUS "FFTW looks for threads")
  endif()
  if (FFTW_FIND_REQUIRED AND FFTW_FIND_REQUIRED_THREADS)
    find_package(Threads REQUIRED)
  else()
    find_package(Threads)
  endif()
endif()

if (FFTW_LOOK_FOR_OMP)
  if (NOT FFTW_FIND_QUIETLY)
    message(STATUS "FFTW looks for openmp")
  endif()
  if (FFTW_FIND_REQUIRED AND FFTW_FIND_REQUIRED_OMP)
    find_package(OpenMP REQUIRED)
  else()
    find_package(OpenMP)
  endif()
endif()

if (FFTW_LOOK_FOR_MKL)
  if (NOT FFTW_FIND_QUIETLY)
    message(STATUS "FFTW looks for threads and Intel MKL")
  endif()
  if (FFTW_LOOK_FOR_THREADS)
    set(BLA_VENDOR "Intel10_64lp")
  else()
    set(BLA_VENDOR "Intel10_64lp_seq")
  endif()
  if (FFTW_FIND_REQUIRED AND FFTW_FIND_REQUIRED_MKL)
    find_package(Threads REQUIRED)
    find_package(BLASEXT REQUIRED)
  else()
    find_package(Threads)
    find_package(BLASEXT)
  endif()
endif()

if (FFTW_LOOK_FOR_ESSL)
  if (NOT FFTW_FIND_QUIETLY)
    message(STATUS "FFTW looks for IBM ESSL")
  endif()
  if (FFTW_LOOK_FOR_THREADS)
    set(BLA_VENDOR "IBMESSLMT")
  else()
    set(BLA_VENDOR "IBMESSL")
  endif()
  if (FFTW_FIND_REQUIRED AND FFTW_FIND_REQUIRED_ESSL)
    find_package(BLAS REQUIRED)
  else()
    find_package(BLAS)
  endif()
endif()

if( THREADS_FOUND AND NOT THREADS_PREFER_PTHREAD_FLAG)
  libraries_absolute_path(CMAKE_THREAD_LIBS_INIT "")
endif ()

# Set variables from environment if needed
# ----------------------------------------
morse_find_package_get_envdir(FFTW)

set(FFTW_GIVEN_BY_USER "FALSE")
if ( FFTW_DIR OR ( FFTW_INCDIR AND FFTW_LIBDIR ) )
  set(FFTW_GIVEN_BY_USER "TRUE")
endif()


# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
if (NOT FFTW_LOOK_FOR_MKL AND NOT FFTW_LOOK_FOR_ESSL)
  find_package(PkgConfig QUIET)
  if( PKG_CONFIG_EXECUTABLE AND NOT FFTW_GIVEN_BY_USER )

    set(FFTW_INCLUDE_DIRS)
    set(FFTW_LIBRARY_DIRS)
    set(FFTW_LIBRARIES)

    if(FFTW_LOOK_FOR_FFTW_SIMPLE)
      pkg_search_module(FFTW3F fftw3f)
      pkg_search_module(FFTW3 fftw3)
      if (FFTW3F_FOUND)
        if (NOT FFTW3F_INCLUDE_DIRS)
          pkg_get_variable(FFTW3F_INCLUDE_DIRS fftw3f includedir)
        endif()
        if (NOT FFTW_FIND_QUIETLY)
          message(STATUS "Looking for FFTW3F - found using PkgConfig")
        endif()
        if (FFTW3F_LIBRARIES)
          morse_find_pkgconfig_libraries_absolute_path(FFTW3F)
          list(APPEND FFTW_LIBRARIES "${FFTW3F_LIBRARIES}")
        endif()
        if(FFTW3F_INCLUDE_DIRS)
          list(APPEND FFTW_INCLUDE_DIRS "${FFTW3F_INCLUDE_DIRS}")
        else()
          if (NOT FFTW_FIND_QUIETLY)
            message(WARNING "FFTW3F_INCLUDE_DIRS is empty using PkgConfig."
              "Perhaps the path to fftw3f headers is already present in your"
              "CPATH/C(PLUS)_INCLUDE_PATH environment variables.")
          endif()
        endif()
        if(FFTW3F_LIBRARY_DIRS)
          list(APPEND FFTW_LIBRARY_DIRS "${FFTW3F_LIBRARY_DIRS}")
        endif()
      else(FFTW3F_FOUND)
        if (NOT FFTW_FIND_QUIETLY)
          message(STATUS "Looking for FFTW3F - not found using PkgConfig."
            "\n   Perhaps you should add the directory containing fftw3f.pc to"
            "\n   the PKG_CONFIG_PATH environment variable.")
        endif()
      endif(FFTW3F_FOUND)
    elseif(FFTW_LOOK_FOR_FFTW_LONG)
      pkg_search_module(FFTW3L fftw3l)
      pkg_search_module(FFTW3 fftw3)
      if (FFTW3L_FOUND)
        if (NOT FFTW3L_INCLUDE_DIRS)
          pkg_get_variable(FFTW3L_INCLUDE_DIRS fftw3l includedir)
        endif()
        if (NOT FFTW_FIND_QUIETLY)
          message(STATUS "Looking for FFTW3L - found using PkgConfig")
        endif()
        if (FFTW3L_LIBRARIES)
          morse_find_pkgconfig_libraries_absolute_path(FFTW3L)
          list(APPEND FFTW_LIBRARIES "${FFTW3L_LIBRARIES}")
        endif()
        if(FFTW3L_INCLUDE_DIRS)
          list(APPEND FFTW_INCLUDE_DIRS "${FFTW3L_INCLUDE_DIRS}")
        else()
          if (NOT FFTW_FIND_QUIETLY)
            message(WARNING "FFTW3L_INCLUDE_DIRS is empty using PkgConfig."
              "Perhaps the path to fftw3l headers is already present in your"
              "CPATH/C(PLUS)_INCLUDE_PATH environment variables.")
          endif()
        endif()
        if(FFTW3L_LIBRARY_DIRS)
          list(APPEND FFTW_LIBRARY_DIRS "${FFTW3L_LIBRARY_DIRS}")
        endif()
      else(FFTW3L_FOUND)
        if (NOT FFTW_FIND_QUIETLY)
          message(STATUS "Looking for FFTW3L - not found using PkgConfig."
            "\n   Perhaps you should add the directory containing fftw3l.pc to"
            "\n   the PKG_CONFIG_PATH environment variable.")
        endif()
      endif(FFTW3L_FOUND)
    elseif(FFTW_LOOK_FOR_FFTW_QUAD)
      pkg_search_module(FFTW3Q fftw3q)
      pkg_search_module(FFTW3 fftw3)
      if (FFTW3Q_FOUND)
        if (NOT FFTW3Q_INCLUDE_DIRS)
          pkg_get_variable(FFTW3Q_INCLUDE_DIRS fftw3q includedir)
        endif()
        if (NOT FFTW_FIND_QUIETLY)
          message(STATUS "Looking for FFTW3Q - found using PkgConfig")
        endif()
        if (FFTW3Q_LIBRARIES)
          morse_find_pkgconfig_libraries_absolute_path(FFTW3Q)
          list(APPEND FFTW_LIBRARIES "${FFTW3Q_LIBRARIES}")
        endif()
        if(FFTW3Q_INCLUDE_DIRS)
          list(APPEND FFTW_INCLUDE_DIRS "${FFTW3Q_INCLUDE_DIRS}")
        else()
          if (NOT FFTW_FIND_QUIETLY)
            message(WARNING "FFTW3Q_INCLUDE_DIRS is empty using PkgConfig."
              "Perhaps the path to fftw3q headers is already present in your"
              "CPATH/C(PLUS)_INCLUDE_PATH environment variables.")
          endif()
        endif()
        if(FFTW3Q_LIBRARY_DIRS)
          list(APPEND FFTW_LIBRARY_DIRS "${FFTW3Q_LIBRARY_DIRS}")
        endif()
      else(FFTW3Q_FOUND)
        if (NOT FFTW_FIND_QUIETLY)
          message(STATUS "Looking for FFTW3Q - not found using PkgConfig."
            "\n   Perhaps you should add the directory containing fftw3q.pc to"
            "\n   the PKG_CONFIG_PATH environment variable.")
        endif()
      endif(FFTW3Q_FOUND)
    else()
      pkg_search_module(FFTW3 fftw3)
      if (FFTW3_FOUND AND FFTW3_LIBRARIES)
        morse_find_pkgconfig_libraries_absolute_path(FFTW3)
      endif()
    endif()
    if (FFTW3_FOUND)
      if (NOT FFTW_FIND_QUIETLY)
        message(STATUS "Looking for FFTW3 - found using PkgConfig")
      endif()
      if (NOT FFTW3_INCLUDE_DIRS)
        pkg_get_variable(FFTW3_INCLUDE_DIRS fftw3 includedir)
      endif()
      if (FFTW3_LIBRARIES)
        morse_find_pkgconfig_libraries_absolute_path(FFTW3)
        list(APPEND FFTW_LIBRARIES "${FFTW3_LIBRARIES}")
      endif()
      if(FFTW3_INCLUDE_DIRS)
            list(APPEND FFTW_INCLUDE_DIRS "${FFTW3_INCLUDE_DIRS}")
      else()
        if (NOT FFTW_FIND_QUIETLY)
          message(WARNING "FFTW3_INCLUDE_DIRS is empty using PkgConfig."
            "Perhaps the path to fftw3 headers is already present in your"
            "CPATH/C(PLUS)_INCLUDE_PATH environment variables.")
        endif()
      endif()
      if(FFTW3_LIBRARY_DIRS)
            list(APPEND FFTW_LIBRARY_DIRS "${FFTW3_LIBRARY_DIRS}")
      endif()
    else(FFTW3_FOUND)
      if (NOT FFTW_FIND_QUIETLY)
        message(STATUS "Looking for FFTW3 - not found using PkgConfig."
          "\n   Perhaps you should add the directory containing fftw3.pc to"
          "\n   the PKG_CONFIG_PATH environment variable.")
      endif()
    endif(FFTW3_FOUND)

    if (FFTW_FOUND AND FFTW_LIBRARIES)
      set(FFTW_FOUND_WITH_PKGCONFIG "TRUE")
    else()
      set(FFTW_FOUND_WITH_PKGCONFIG "FALSE")
    endif()

  endif( PKG_CONFIG_EXECUTABLE AND NOT FFTW_GIVEN_BY_USER )

endif(NOT FFTW_LOOK_FOR_MKL AND NOT FFTW_LOOK_FOR_ESSL)

if( (NOT PKG_CONFIG_EXECUTABLE) OR
    (PKG_CONFIG_EXECUTABLE AND NOT FFTW_FOUND) OR
    FFTW_GIVEN_BY_USER OR
    FFTW_LOOK_FOR_MKL  OR
    FFTW_LOOK_FOR_ESSL
    )

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
  set(PATH_TO_LOOK_FOR)
  if (DEFINED ENV{MKLROOT})
    list(APPEND PATH_TO_LOOK_FOR "$ENV{MKLROOT}/include/fftw")
  endif()

  if (FFTW_LOOK_FOR_ESSL)
    set(FFTW3_HEADER_TO_FIND "fftw3_essl.h")
  else()
    set(FFTW3_HEADER_TO_FIND "fftw3.h")
  endif()

  morse_find_path(FFTW
    HEADERS  ${FFTW3_HEADER_TO_FIND}
    SUFFIXES "include" "include/fftw" "fftw"
    HINTS    ${PATH_TO_LOOK_FOR})

  # Looking for lib
  # ---------------

  # Add system library paths to search lib
  # --------------------------------------
  set(PATH_TO_LOOK_FOR)
  if (DEFINED $ENV{MKLROOT})
    list(APPEND PATH_TO_LOOK_FOR "$ENV{MKLROOT}/lib")
  endif()

  set( _lib_suffixes lib )
  if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
    list(APPEND _lib_suffixes "lib64" "lib/intel64" )
  else()
    list(APPEND _lib_suffixes "lib32" "lib/ia32" )
  endif()

  if(FFTW_LOOK_FOR_FFTW_SIMPLE)
    set(FFTW_PREC "f")
    set(FFTW_PREC_TESTFUNC "s")
  elseif(FFTW_LOOK_FOR_FFTW_LONG)
    set(FFTW_PREC "l")
    set(FFTW_PREC_TESTFUNC "l")
  elseif(FFTW_LOOK_FOR_FFTW_QUAD)
    set(FFTW_PREC "q")
    set(FFTW_PREC_TESTFUNC "q")
  else()
    set(FFTW_PREC "")
    set(FFTW_PREC_TESTFUNC "d")
  endif()

  set(FFTW_LIBRARIES "")
  set(FFTW_LIBRARY_DIRS "")

  if(NOT FFTW_LOOK_FOR_MKL)
    if (FFTW_LOOK_FOR_THREADS)
      set(FFTW_libs_to_find "fftw3${FFTW_PREC}_threads;fftw3${FFTW_PREC};fftw3")
    elseif (FFTW_LOOK_FOR_OMP)
      set(FFTW_libs_to_find "fftw3${FFTW_PREC}_omp;fftw3${FFTW_PREC};fftw3")
    else()
      set(FFTW_libs_to_find "fftw3${FFTW_PREC};fftw3")
    endif()
    if (FFTW_LOOK_FOR_FFTW_QUAD)
      if (NOT FFTW_LOOK_FOR_MKL AND NOT FFTW_LOOK_FOR_ESSL)
        list(APPEND FFTW_libs_to_find "quadmath")
      endif()
    endif()

    if (FFTW_LOOK_FOR_ESSL)
      set(FFTW_libs_to_find "fftw3_essl")
    endif()

    # Try to find the fftw lib in the given paths
    # ----------------------------------------------
    morse_find_library(FFTW
      LIBRARIES ${FFTW_libs_to_find}
      SUFFIXES  ${_lib_suffixes}
      HINTS     ${PATH_TO_LOOK_FOR})

  endif(NOT FFTW_LOOK_FOR_MKL)

  if (FFTW_LOOK_FOR_MKL OR FFTW_LOOK_FOR_ESSL)

    # FFTW relies on blas libs
    if (FFTW_LOOK_FOR_THREADS)
      if (FFTW_LOOK_FOR_MKL)
          if (BLAS_MT_LIBRARIES)
            list(APPEND FFTW_LIBRARIES "${BLAS_MT_LIBRARIES}")
          if (NOT FFTW_FIND_QUIETLY)
            message(STATUS "Multithreaded FFTW has been found: ${FFTW_LIBRARIES}")
          endif()
        else()
          if (NOT FFTW_FIND_QUIETLY)
            if (FFTW_FIND_REQUIRED AND FFTW_FIND_REQUIRED_MKL)
              message(FATAL_ERROR "FFTW is required but not found.")
            else()
              message(STATUS "Multithreaded FFTW not found.")
            endif()
          endif()
      endif(BLAS_MT_LIBRARIES)
      elseif (FFTW_LOOK_FOR_ESSL)
        if (FFTW_LIBRARIES AND BLAS_LIBRARIES)
          list(APPEND FFTW_LIBRARIES "${BLAS_LIBRARIES}")
          if (NOT FFTW_FIND_QUIETLY)
            message(STATUS "Multithreaded FFTW has been found: ${FFTW_LIBRARIES}")
          endif()
        else()
          if (NOT FFTW_FIND_QUIETLY)
            if (FFTW_FIND_REQUIRED AND FFTW_FIND_REQUIRED_MKL)
              message(FATAL_ERROR "FFTW is required but not found.")
            else()
              message(STATUS "Multithreaded FFTW not found.")
            endif()
          endif()
        endif(FFTW_LIBRARIES AND BLAS_LIBRARIES)
      endif()
    else(FFTW_LOOK_FOR_THREADS)
      if (FFTW_LOOK_FOR_MKL)
          if (BLAS_SEQ_LIBRARIES)
          list(APPEND FFTW_LIBRARIES "${BLAS_SEQ_LIBRARIES}")
          if (NOT FFTW_FIND_QUIETLY)
            message(STATUS "FFTW has been found: ${FFTW_LIBRARIES}")
          endif()
        else()
          if (NOT FFTW_FIND_QUIETLY)
            if (FFTW_FIND_REQUIRED AND FFTW_FIND_REQUIRED_MKL)
              message(FATAL_ERROR "FFTW is required but not found.")
            else()
              message(STATUS "FFTW not found.")
            endif()
          endif()
      endif(BLAS_SEQ_LIBRARIES)
      elseif (FFTW_LOOK_FOR_ESSL)
        if (FFTW_LIBRARIES AND BLAS_LIBRARIES)
          list(APPEND FFTW_LIBRARIES "${BLAS_LIBRARIES}")
          if (NOT FFTW_FIND_QUIETLY)
            message(STATUS "FFTW has been found: ${FFTW_LIBRARIES}")
          endif()
        else()
          if (NOT FFTW_FIND_QUIETLY)
            if (FFTW_FIND_REQUIRED AND FFTW_FIND_REQUIRED_MKL)
              message(FATAL_ERROR "FFTW is required but not found.")
            else()
              message(STATUS "FFTW not found.")
            endif()
          endif()
        endif(FFTW_LIBRARIES AND BLAS_LIBRARIES)
      endif()
    endif(FFTW_LOOK_FOR_THREADS)

    if (BLAS_LIBRARY_DIRS)
      list(APPEND FFTW_LIBRARY_DIRS "${BLAS_LIBRARY_DIRS}")
    else()
      if (NOT FFTW_FIND_QUIETLY)
        message(WARNING "FFTW_LIBRARY_DIRS may not be complete because BLAS_LIBRARY_DIRS is empty.")
      endif()
    endif()

  endif(FFTW_LOOK_FOR_MKL OR FFTW_LOOK_FOR_ESSL)

  list(REMOVE_DUPLICATES FFTW_INCLUDE_DIRS)
  list(REMOVE_DUPLICATES FFTW_LIBRARY_DIRS)

  # check if one lib is NOTFOUND
  foreach(lib ${FFTW_LIBRARIES})
    if (NOT lib)
      set(FFTW_LIBRARIES "FFTW_LIBRARIES-NOTFOUND")
    endif()
  endforeach()

endif( (NOT PKG_CONFIG_EXECUTABLE) OR
  (PKG_CONFIG_EXECUTABLE AND NOT FFTW_FOUND) OR
  FFTW_GIVEN_BY_USER OR
  FFTW_LOOK_FOR_MKL  OR
  FFTW_LOOK_FOR_ESSL
  )

# check a function to validate the find
if(FFTW_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)

  # FFTW
  if (FFTW_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${FFTW_INCLUDE_DIRS}")
  endif()
  if (FFTW_CFLAGS_OTHER)
    set(REQUIRED_FLAGS "${FFTW_CFLAGS_OTHER}")
  endif()
  if (FFTW_LDFLAGS_OTHER)
    set(REQUIRED_LDFLAGS "${FFTW_LDFLAGS_OTHER}")
  endif()
  if (FFTW_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${FFTW_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${FFTW_LIBRARIES}")
  # THREADS
  if (FFTW_LOOK_FOR_THREADS)
    if (THREADS_PREFER_PTHREAD_FLAG)
      list(APPEND REQUIRED_FLAGS "${CMAKE_THREAD_LIBS_INIT}")
    else()
      list(APPEND REQUIRED_LIBS "${CMAKE_THREAD_LIBS_INIT}")
    endif()
  endif()
  # OMP
  if(FFTW_LOOK_FOR_OMP)
    list(APPEND REQUIRED_FLAGS "${OPENMP_C_FLAGS}")
  endif()
  # MKL
  if(FFTW_LOOK_FOR_MKL)
    if (THREADS_PREFER_PTHREAD_FLAG)
      list(APPEND REQUIRED_FLAGS "${CMAKE_THREAD_LIBS_INIT}")
    else()
      list(APPEND REQUIRED_LIBS "${CMAKE_THREAD_LIBS_INIT}")
    endif()
    if (CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
      list(APPEND REQUIRED_LDFLAGS "-Wl,--no-as-needed")
    endif()
  endif()
  # m
  find_library(M_LIBRARY NAMES m)
  mark_as_advanced(M_LIBRARY)
  if(M_LIBRARY)
    list(APPEND REQUIRED_LIBS "-lm")
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
  unset(FFTW_WORKS CACHE)
  include(CheckFunctionExists)
  if (FFTW_LOOK_FOR_ESSL)
    check_function_exists(${FFTW_PREC_TESTFUNC}fftw_execute FFTW_WORKS)
  else()
    check_function_exists(${FFTW_PREC_TESTFUNC}fftw_execute_ FFTW_WORKS)
  endif()
  mark_as_advanced(FFTW_WORKS)

  if(FFTW_WORKS)
    # save link with dependencies
    set(FFTW_LIBRARIES "${REQUIRED_LIBS}")
    set(FFTW_LIBRARY_DIRS "${REQUIRED_LIBDIRS}")
    set(FFTW_INCLUDE_DIRS "${REQUIRED_INCDIRS}")
    set(FFTW_CFLAGS_OTHER "${REQUIRED_FLAGS}")
    set(FFTW_LDFLAGS_OTHER "${REQUIRED_LDFLAGS}")
  else()
    if(NOT FFTW_FIND_QUIETLY)
      message(STATUS "Looking for FFTW : test of ${FFTW_PREC_TESTFUNC}fftw_execute_ with fftw library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

  list(GET FFTW_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" DIRECTORY)
  if (NOT FFTW_LIBRARY_DIRS)
    set(FFTW_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)")
    string(REGEX REPLACE "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)" "" not_cached_dir "${first_lib_path}")
    set(FFTW_PREFIX "${not_cached_dir}" CACHE PATH "Installation directory of FFTW library" FORCE)
  else()
    set(FFTW_PREFIX "${first_lib_path}" CACHE PATH "Installation directory of FFTW library" FORCE)
  endif()
  mark_as_advanced(FFTW_DIR)
  mark_as_advanced(FFTW_PREFIX)

endif(FFTW_LIBRARIES)

# check that FFTW has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG
  FFTW_LIBRARIES
  FFTW_WORKS)

# Add imported target
if (FFTW_FOUND)
  if(NOT TARGET FFTW::FFTW)
    add_library(FFTW::FFTW INTERFACE IMPORTED)
    set_property(TARGET FFTW::FFTW APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}")
    set_property(TARGET FFTW::FFTW APPEND PROPERTY INTERFACE_LINK_DIRECTORIES "${FFTW_LIBRARY_DIRS}")
    set_property(TARGET FFTW::FFTW APPEND PROPERTY INTERFACE_LINK_LIBRARIES "${FFTW_LIBRARIES}")
    set_property(TARGET FFTW::FFTW APPEND PROPERTY INTERFACE_COMPILE_OPTIONS "${FFTW_CFLAGS_OTHER}")
    set_property(TARGET FFTW::FFTW APPEND PROPERTY INTERFACE_LINNK_OPTIONS "${FFTW_LDFLAGS_OTHER}")
  endif()
  morse_create_imported_target(FFTW)
endif()
