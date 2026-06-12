###
#
# @copyright (c) 2012-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                          Univ. Bordeaux. All rights reserved.
# @copyright (c) 2009-2023 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
#
###
#
#  @file SetCMakeCudaArchitectures.cmake
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
#  @author Mathieu Faverge
#  @date 04-07-2023
#
###
if (NOT CUDAToolkit_FOUND)
  message( STATUS "SetCMakeCudaArchitectures: This file should not be included if CUDA Toolkit has not bee found" )
  return()
endif()

set(CMAKE_CUDA_STANDARD 11)

# Define the architectures (inspired from the MAGMA project)
set( ${CMAKE_PROJECT_NAME}_CUDA_TARGETS "Maxwell Pascal Volta Ampere" CACHE STRING "CUDA architectures to compile for; one or more of Fermi, Kepler, Maxwell, Pascal, Volta, Turing, Ampere, Hopper, or valid sm_[0-9][0-9]" )

# NVCC options for the different cards
# sm_xx is binary, compute_xx is PTX for forward compatibility
# MIN_ARCH is the lowest requested version

if(WIN32)
  # Disable separable compilation on Windows because object linking list
  # becomes too long when building multiple archs and MSVC throws errors
  set(CUDA_SEPARABLE_COMPILATION OFF)
else()
  set(CUDA_SEPARABLE_COMPILATION ON)
endif()

set(__cuda_architectures)

# Architectures by names
# ----------------------
if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES Fermi)
  set( ${CMAKE_PROJECT_NAME}_CUDA_TARGETS "${${CMAKE_PROJECT_NAME}_CUDA_TARGETS} sm_20" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES Kepler)
  set( ${CMAKE_PROJECT_NAME}_CUDA_TARGETS "${${CMAKE_PROJECT_NAME}_CUDA_TARGETS} sm_30 sm_35 sm_37" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES Maxwell)
  set( ${CMAKE_PROJECT_NAME}_CUDA_TARGETS "${${CMAKE_PROJECT_NAME}_CUDA_TARGETS} sm_50" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES Pascal)
  set( ${CMAKE_PROJECT_NAME}_CUDA_TARGETS "${${CMAKE_PROJECT_NAME}_CUDA_TARGETS} sm_60" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES Volta)
  set( ${CMAKE_PROJECT_NAME}_CUDA_TARGETS "${${CMAKE_PROJECT_NAME}_CUDA_TARGETS} sm_70" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES Turing)
  set( ${CMAKE_PROJECT_NAME}_CUDA_TARGETS "${${CMAKE_PROJECT_NAME}_CUDA_TARGETS} sm_75" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES Ampere)
  set( ${CMAKE_PROJECT_NAME}_CUDA_TARGETS "${${CMAKE_PROJECT_NAME}_CUDA_TARGETS} sm_80" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES Hopper)
  set( ${CMAKE_PROJECT_NAME}_CUDA_TARGETS "${${CMAKE_PROJECT_NAME}_CUDA_TARGETS} sm_90" )
endif()

# Architectures versions
# ----------------------
if ( (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES sm_20) AND (CUDA_VERSION VERSION_LESS "8.0") )
  if (NOT MIN_ARCH)
    set( MIN_ARCH 200 )
  endif()
  list(APPEND __cuda_architectures 20)
  message( STATUS "    compile for CUDA arch 2.0 (Fermi)" )
endif()

if ( (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES sm_30) AND (CUDA_VERSION VERSION_LESS "10.0") )
  if (NOT MIN_ARCH)
    set( MIN_ARCH 300 )
  endif()
  list(APPEND __cuda_architectures 30)
  message( STATUS "    compile for CUDA arch 3.0 (Kepler)" )
endif()

if ( (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES sm_35) AND (CUDA_VERSION VERSION_LESS "11.0") )
  if (NOT MIN_ARCH)
    set( MIN_ARCH 300 )
  endif()
  list(APPEND __cuda_architectures 35)
  message( STATUS "    compile for CUDA arch 3.5 (Kepler)" )
endif()

if ( (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES sm_50) AND (CUDA_VERSION VERSION_LESS "11.0") )
  if (NOT MIN_ARCH)
    set( MIN_ARCH 500 )
  endif()
  list(APPEND __cuda_architectures 50)
  message( STATUS "    compile for CUDA arch 5.0 (Maxwell)" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES sm_52)
  if (NOT MIN_ARCH)
    set( MIN_ARCH 520 )
  endif()
  list(APPEND __cuda_architectures 52)
  message( STATUS "    compile for CUDA arch 5.2 (Maxwell)" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES sm_53)
  if (NOT MIN_ARCH)
    set( MIN_ARCH 530 )
  endif()
  list(APPEND __cuda_architectures 53)
  message( STATUS "    compile for CUDA arch 5.3 (Maxwell)" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES sm_60)
  if (NOT MIN_ARCH)
    set( MIN_ARCH 600 )
  endif()
  list(APPEND __cuda_architectures 60)
  message( STATUS "    compile for CUDA arch 6.0 (Pascal)" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES sm_61)
  if (NOT MIN_ARCH)
    set( MIN_ARCH 610 )
  endif()
  list(APPEND __cuda_architectures 61)
  message( STATUS "    compile for CUDA arch 6.1 (Pascal)" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES sm_62)
  if (NOT MIN_ARCH)
    set( MIN_ARCH 620 )
  endif()
  list(APPEND __cuda_architectures 62)
  message( STATUS "    compile for CUDA arch 6.2 (Pascal)" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES sm_70)
  if (NOT MIN_ARCH)
    set( MIN_ARCH 700 )
  endif()
  list(APPEND __cuda_architectures 70)
  message( STATUS "    compile for CUDA arch 7.0 (Volta)" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES sm_71)
  if (NOT MIN_ARCH)
    set( MIN_ARCH 710 )
  endif()
  list(APPEND __cuda_architectures 71)
  message( STATUS "    compile for CUDA arch 7.1 (Volta)" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES sm_75)
  if (NOT MIN_ARCH)
    set( MIN_ARCH 750 )
  endif()
  list(APPEND __cuda_architectures 75)
  message( STATUS "    compile for CUDA arch 7.5 (Turing)" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES sm_80)
  if (NOT MIN_ARCH)
    set( MIN_ARCH 800 )
  endif()
  list(APPEND __cuda_architectures 80)
  message( STATUS "    compile for CUDA arch 8.0 (Ampere)" )
endif()

if (${CMAKE_PROJECT_NAME}_CUDA_TARGETS MATCHES sm_90)
  if (NOT MIN_ARCH)
    set( MIN_ARCH 900 )
  endif()
  list(APPEND __cuda_architectures 90)
  message( STATUS "    compile for CUDA arch 9.0 (Hopper)" )
endif()

if (NOT MIN_ARCH)
  message( FATAL_ERROR "${CMAKE_PROJECT_NAME}_CUDA_TARGETS must contain one or more of Fermi, Kepler, Maxwell, Pascal, Volta, Turing, Ampere, or valid sm_[0-9][0-9]" )
endif()

# Remove extra
# ------------
if(CUDA_VERSION VERSION_GREATER_EQUAL "8.0")
  list(REMOVE_ITEM __cuda_architectures "20" "21")
endif()

if(CUDA_VERSION VERSION_GREATER_EQUAL "9.0")
  list(REMOVE_ITEM __cuda_architectures "20" "21")
endif()

if(CUDA_VERSION VERSION_GREATER_EQUAL "10.0")
  list(REMOVE_ITEM __cuda_architectures "30" "32")
endif()

if(CUDA_VERSION VERSION_GREATER_EQUAL "11.0")
  list(REMOVE_ITEM __cuda_architectures "35" "50")
endif()

set(CMAKE_CUDA_ARCHITECTURES "${__cuda_architectures}")
