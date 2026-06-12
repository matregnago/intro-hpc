###
#
# @file CMakeLists.txt
#
# @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
###
#
#  @project CHAMELEON
#  CHAMELEON is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 1.4.0
#  @author Florent Pruvost
#  @date 2025-12-19
#
###
cmake_minimum_required(VERSION 3.12)
project(LINK_CHAMELEON C)

#
# By default the Fortran example is not enabled as the fortran wrapper is no longer maintained
#
option(FORTRAN "Enable Fortran example" OFF)
option(CXX     "Enable C++ example"     ON )

find_package(CHAMELEON REQUIRED)

if (TARGET CHAMELEON::chameleon)
  message(STATUS "CHAMELEON_BIN_DIR ${CHAMELEON_BIN_DIR}")
  message(STATUS "CHAMELEON_INC_DIR ${CHAMELEON_INC_DIR}")
  message(STATUS "CHAMELEON_LIB_DIR ${CHAMELEON_LIB_DIR}")
else()
  message(FATAL_ERROR "Target CHAMELEON::chameleon is not found, check your CHAMELEONConfig.cmake.")
endif()

include(CheckSymbolExists)
set(CMAKE_REQUIRED_INCLUDES "${CHAMELEON_INC_DIR}")
check_symbol_exists( CHAMELEON_PREC_S  "chameleon/config.h" HAVE_CHAMELEON_PREC_S  )
check_symbol_exists( CHAMELEON_PREC_D  "chameleon/config.h" HAVE_CHAMELEON_PREC_D  )
check_symbol_exists( CHAMELEON_PREC_C  "chameleon/config.h" HAVE_CHAMELEON_PREC_C  )
check_symbol_exists( CHAMELEON_USE_MPI "chameleon/config.h" HAVE_CHAMELEON_USE_MPI )
unset(CMAKE_REQUIRED_INCLUDES)

set(CMAKE_BUILD_RPATH ${CHAMELEON_LIB_DIR})
enable_testing()
include(CTest)

include(CheckLanguage)

#
# Check link in C
#
if (HAVE_CHAMELEON_PREC_S)
  add_executable(link_chameleon_c link_chameleon.c)
  target_link_libraries(link_chameleon_c PRIVATE CHAMELEON::chameleon)
  add_test(link_chameleon_c link_chameleon_c)
endif()

#
# Check link in Fortran
#
if (HAVE_CHAMELEON_PREC_D AND FORTRAN)
  check_language(Fortran)
  if( NOT CMAKE_Fortran_COMPILER)
    message(STATUS "No Fortran support")
  else()
    enable_language(Fortran)

    add_executable(link_chameleon_f link_chameleon.f90)
    target_link_libraries(link_chameleon_f PRIVATE CHAMELEON::chameleon)
    add_test(link_chameleon_f link_chameleon_f)
  endif()
endif()

#
# Check link in C++
#
if (HAVE_CHAMELEON_PREC_C AND CXX)
  check_language(CXX)
  if( NOT CMAKE_CXX_COMPILER)
    message(STATUS "No C++ support")
  else()
    enable_language(CXX)

    add_executable(link_chameleon_cpp link_chameleon.cpp)
    target_link_libraries(link_chameleon_cpp PRIVATE CHAMELEON::chameleon)

    # Add the CXX MPI library if needed
    if (HAVE_CHAMELEON_USE_MPI)
      find_package(MPI REQUIRED COMPONENTS CXX)
      target_link_libraries(link_chameleon_cpp PRIVATE MPI::MPI_CXX)
    endif()

    add_test(link_chameleon_cpp link_chameleon_cpp)
  endif()
endif()
###
### END CMakeLists.txt
###
