###
#
# @copyright (c) 2012-2020 Inria. All rights reserved.
# @copyright (c) 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
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
# - Find STARPU include dirs and libraries  using pkg-config file.
# Use this module by invoking find_package with the form:
#  find_package(STARPU
#               [version] [EXACT] # Minimum or EXACT version e.g. 1.4
#               [REQUIRED] # Fail with error if starpu is not found
#               [COMPONENTS components...]) # to search specific variants
#
#   Optional dependencies
#   - BLAS
#   - CUDA
#   - HIP
#   - FORTRAN interface
#   - FXT
#   - HDF5
#   - MPI
#   - OPENCL
#   - OPENMP
#   - SIMGRID
#
#  COMPONENTS can be some of the following:
#   - MPI: to search STARPU configured with MPI
#   - CUDA: to search STARPU configured with CUDA
#   - HIP: to search STARPU configured with HIP
#   - OPENCL: to search STARPU configured with OpenCL
#   - FXT: to search STARPU configured with FXT
#   - SIMGRID: to search STARPU configured with SimGrid
#   - FORTRAN: to search STARPU Fortran interface
#
#  STARPU_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = STARPU
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
#  STARPU_FORTRAN_MOD            - Points to the StarPU Fortran interface module fstarpu_mod.f90
#
# Set STARPU_STATIC to 1 to force using static libraries if exist.
#
# This module defines the following :prop_tgt:`IMPORTED` target:
#
# ``MORSE::STARPU``
#   The CMake target to use for STARPU if found.
#
#=============================================================================

# Common macros to use in finds
include(FindMorseInit)

if (NOT STARPU_FIND_QUIETLY)
  message(STATUS "FindSTARPU needs pkg-config program and PKG_CONFIG_PATH set with starpu-x.y.pc (x.y the version) file path.")
endif()

if (NOT STARPU_FIND_VERSION)
  set(STARPU_FIND_VERSION "1.0")
  message(STATUS "FindSTARPU needs a version to check minimal API requirement. As it is not given we set to 1.0 by default.")
endif()

# use pkg-config to detect include/library dirs (if pkg-config is available)
# --------------------------------------------------------------------------
if (PKG_CONFIG_EXECUTABLE)
  unset(STARPU_FOUND CACHE)
  if ("${STARPU_FIND_COMPONENTS}" MATCHES "MPI")
    pkg_search_module(STARPU starpumpi-${STARPU_FIND_VERSION})
  else()
    pkg_search_module(STARPU starpu-${STARPU_FIND_VERSION})
  endif()
  if (NOT STARPU_FIND_QUIETLY)
    if (STARPU_FOUND AND STARPU_LIBRARIES)
      message(STATUS "Looking for STARPU - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for STARPU - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing libstarpumpi.pc, libstarpu.pc to"
        "\n   the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()
  if (STARPU_FOUND AND STARPU_LIBRARIES)
    if (NOT STARPU_INCLUDE_DIRS)
      pkg_get_variable(STARPU_INCLUDE_DIRS libstarpu includedir)
    endif()
    set(STARPU_FOUND_WITH_PKGCONFIG "TRUE")
    morse_find_pkgconfig_libraries_absolute_path(STARPU)
  else()
    set(STARPU_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

  if (STARPU_STATIC AND STARPU_STATIC_LIBRARIES)
    set (STARPU_DEPENDENCIES ${STARPU_STATIC_LIBRARIES})
    list (REMOVE_ITEM STARPU_DEPENDENCIES "starpu")
    list (REMOVE_ITEM STARPU_DEPENDENCIES "starpumpi")
    list (APPEND STARPU_LIBRARIES ${STARPU_DEPENDENCIES})
    set(STARPU_CFLAGS_OTHER ${STARPU_STATIC_CFLAGS_OTHER})
    set(STARPU_LDFLAGS_OTHER ${STARPU_STATIC_LDFLAGS_OTHER})
    if (NOT STARPU_FIND_QUIETLY)
      message(STATUS "STARPU_STATIC set to 1 by user, STARPU_LIBRARIES: ${STARPU_LIBRARIES}.")
    endif()
  endif()
endif()

# checks to validate the find
if (STARPU_FOUND AND STARPU_LIBRARIES)

  # check if static or dynamic lib
  morse_check_static_or_dynamic(STARPU STARPU_LIBRARIES)
  if (STARPU_STATIC)
    set(STATIC "_STATIC")
  else()
    set(STATIC "")
  endif()

  # tests used in this script is not compliant with -Werror or -Werror=...
  # remove it temporarily from C flags
  set(CMAKE_C_FLAGS_COPY "${CMAKE_C_FLAGS}" CACHE STRING "" FORCE)
  string(REGEX REPLACE "-Werror[^ ]*" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

  # set required libraries for link
  morse_set_required_test_lib_link(STARPU)

  # test link
  unset(STARPU_INIT CACHE)
  include(CheckFunctionExists)
  check_function_exists(starpu_init STARPU_INIT)
  mark_as_advanced(STARPU_INIT)
  if (NOT STARPU_INIT AND NOT STARPU_FIND_QUIETLY)
      message(STATUS "Looking for starpu : link test with starpu_init call fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Please check that your starpu pkg-config files are fine. "
        "To debug this use cmake option --debug-trycompile and execute 'make' "
        "command in the given directory.")
  endif()

  # add additional checks such as compilation and user's required components:
  # MPI, CUDA, HIP, OPENCL, FXT, SIMGRID, FORTRAN
  if (STARPU_INIT)

    set(STARPU_CHECK TRUE)

    # check basic starpu init call compiles
    morse_cmake_required_set(STARPU)
    unset(STARPU_COMPILES CACHE)
    include(CheckCSourceCompiles)
    if ("${STARPU_LIBRARIES}" MATCHES "starpumpi" )
      # there is a problem with the pkg-config file of starpumpi
      # it does not provide the MPI include dirs and libraries while required
      # for the compilation
      if (NOT MPI_C_FOUND)
        find_package(MPI REQUIRED C)
      endif()
      list(APPEND CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_DIRS})
      list(APPEND CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
      check_c_source_compiles("
#include <starpu.h>
#include <starpu_mpi.h>
int main(void) {
  return starpu_mpi_init_conf(NULL, NULL, 1, MPI_COMM_WORLD, NULL);
}
"       STARPU_COMPILES)
    else ()
      check_c_source_compiles("
#include <starpu.h>
int main(void) {
  return starpu_init(NULL);
}
"       STARPU_COMPILES)
    endif ()
    mark_as_advanced(STARPU_COMPILES)
    if (NOT STARPU_COMPILES AND STARPU_FIND_REQUIRED )
      message(WARNING "We were not able to compile a C source file calling "
      "starpu_init() or starpu_mpi_init_conf(). Please check that your starpu "
      "pkg-config files are fine. To debug this use cmake option "
      "--debug-trycompile and execute 'make' command in the given directory.")
      set(STARPU_CHECK FALSE)
    endif()

    # check for required CUDA component
    if ("${STARPU_FIND_COMPONENTS}" MATCHES "CUDA")
      unset(STARPU_CUBLAS_GET_LOCAL_HANDLE CACHE)
      check_function_exists(starpu_cublas_get_local_handle STARPU_CUBLAS_GET_LOCAL_HANDLE)
      mark_as_advanced(STARPU_CUBLAS_GET_LOCAL_HANDLE)
      if (NOT STARPU_CUBLAS_GET_LOCAL_HANDLE AND STARPU_FIND_REQUIRED )
        if (NOT STARPU_FIND_QUIETLY)
          message(WARNING "Required starpu component CUDA not found")
        endif()
        set(STARPU_CHECK FALSE)
      endif()
    endif()

    # check for required HIP component
    if ("${STARPU_FIND_COMPONENTS}" MATCHES "HIP")
      unset(STARPU_HIPBLAS_GET_LOCAL_HANDLE CACHE)
      check_function_exists(starpu_hipblas_get_local_handle STARPU_HIPBLAS_GET_LOCAL_HANDLE)
      mark_as_advanced(STARPU_HIPBLAS_GET_LOCAL_HANDLE)
      if (NOT STARPU_HIPBLAS_GET_LOCAL_HANDLE AND STARPU_FIND_REQUIRED )
        if (NOT STARPU_FIND_QUIETLY)
          message(WARNING "Required starpu component HIP not found")
        endif()
        set(STARPU_CHECK FALSE)
      endif()
    endif()

    # check for required OPENCL component
    if ("${STARPU_FIND_COMPONENTS}" MATCHES "OPENCL")
      unset(STARPU_OPENCL_GET_CONTEXT CACHE)
      check_function_exists(starpu_opencl_get_context STARPU_OPENCL_GET_CONTEXT)
      mark_as_advanced(STARPU_OPENCL_GET_CONTEXT)
      if (NOT STARPU_OPENCL_GET_CONTEXT AND STARPU_FIND_REQUIRED )
        if (NOT STARPU_FIND_QUIETLY)
          message(WARNING "Required starpu component OPENCL not found")
        endif()
        set(STARPU_CHECK FALSE)
      endif()
    endif()

    # check for required FXT component
    if ("${STARPU_FIND_COMPONENTS}" MATCHES "FXT")
      unset(STARPU_FXT_START_PROFILING CACHE)
      check_function_exists(starpu_fxt_start_profiling STARPU_FXT_START_PROFILING)
      mark_as_advanced(STARPU_FXT_START_PROFILING)
      if (NOT STARPU_FXT_START_PROFILING AND STARPU_FIND_REQUIRED)
        if (NOT STARPU_FIND_QUIETLY)
          message(WARNING "Required starpu component FXT not found")
        endif()
        set(STARPU_CHECK FALSE)
      endif()
    endif()

    # check for required SIMGRID component
    if ("${STARPU_FIND_COMPONENTS}" MATCHES "SIMGRID")
      check_c_source_compiles("
#include <starpu.h>
#ifdef STARPU_SIMGRID
int main() { return 0; }
#else
#error \"STARPU_SIMGRID is not defined\"
#endif
"       STARPU_SIMGRID_DEFINED)
      mark_as_advanced(STARPU_SIMGRID_DEFINED)
      if (NOT STARPU_SIMGRID_DEFINED AND STARPU_FIND_REQUIRED)
        if (NOT STARPU_FIND_QUIETLY)
          message(WARNING "Required starpu component SIMGRID not found")
        endif()
        set(STARPU_CHECK FALSE)
      endif()
    endif()

    # check for required Fortran interface
    if ("${STARPU_FIND_COMPONENTS}" MATCHES "FORTRAN")
      unset(FSTARPU_TASK_INSERT CACHE)
      check_function_exists(fstarpu_task_insert FSTARPU_TASK_INSERT)
      mark_as_advanced(FSTARPU_TASK_INSERT)
      if (FSTARPU_TASK_INSERT)
        set(STARPU_FMOD_DIRS "STARPU_FMOD_DIRS-NOTFOUND")
        find_path(STARPU_FMOD_DIRS
          NAMES "fstarpu_mod.f90"
          HINTS ${STARPU_INCLUDE_DIRS})
        if (STARPU_FMOD_DIRS)
          set(STARPU_FORTRAN_MOD "${STARPU_FMOD_DIRS}/fstarpu_mod.f90")
          mark_as_advanced(STARPU_FORTRAN_MOD)
          if (NOT STARPU_FIND_QUIETLY)
            message(STATUS "Found ${STARPU_FORTRAN_MOD}")
          endif()
        elseif(STARPU_FIND_REQUIRED)
          if (NOT STARPU_FIND_QUIETLY)
            message(WARNING "Looking for starpu fstarpu_mod.f90: not found")
          endif()
          set(STARPU_CHECK FALSE)
        endif()
      elseif(STARPU_FIND_REQUIRED)
        if (NOT STARPU_FIND_QUIETLY)
          message(WARNING "Looking for starpu fstarpu_task_insert: not found")
        endif()
        set(STARPU_CHECK FALSE)
      endif()
    endif()

    morse_cmake_required_unset()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_COPY}" CACHE STRING "" FORCE)

  endif(STARPU_INIT)

endif(STARPU_FOUND AND STARPU_LIBRARIES)

# check that STARPU has been found
# --------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(STARPU
  REQUIRED_VARS STARPU_INCLUDE_DIRS
                STARPU_LIBRARIES
                STARPU_INIT
                STARPU_COMPILES
                STARPU_CHECK
  VERSION_VAR   STARPU_VERSION)

 # Add imported targe
if (STARPU_FOUND)
  morse_create_imported_target(STARPU)
endif()