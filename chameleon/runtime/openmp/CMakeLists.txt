###
#
# @file openmp/CMakeLists.txt
#
# @copyright 2009-2015 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
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
# @version 1.4.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @author Philippe Virouleau
#  @date 2025-12-19
#
###
cmake_minimum_required(VERSION 3.12)

find_package(OpenMP REQUIRED)
morse_export_imported_target(OpenMP::OpenMP_C openmp chameleon)

if ( OPENMP_FOUND )
  message("-- ${Blue}Add definition CHAMELEON_SCHED_OPENMP"
    " - Activate OpenMP in Chameleon${ColourReset}")
  # TODO: check for compiler support for "_Pragma"
else ( OPENMP_FOUND )
  message(FATAL_ERROR "Something went wrong when finding OpenMP")
  return()
endif ( OPENMP_FOUND )

# Define the list of headers
# --------------------------
set(RUNTIME_HDRS
  include/chameleon_openmp.h
  )

# Force generation of headers
# ---------------------------
add_custom_target(
  runtime_openmp_include
  ALL SOURCES ${RUNTIME_HDRS})

# Installation
# ------------
install(
  FILES ${RUNTIME_HDRS}
  DESTINATION include/runtime/openmp )

# Generate the Chameleon common for all possible precisions
# ---------------------------------------------------------
set(RUNTIME_COMMON_GENERATED "")
set(ZSRC ${CONTROL_ZSRC} )

precisions_rules_py(
  RUNTIME_COMMON_GENERATED "${ZSRC}"
  PRECISIONS "${CHAMELEON_PRECISION}"
  TARGETDIR "control")

set(RUNTIME_COMMON
  ${CONTROL_SRC}
  ${RUNTIME_COMMON_GENERATED}
  )

# Generate the Chameleon sources for all possible precisions
# ----------------------------------------------------------
set(RUNTIME_SRCS_GENERATED "")
set(ZSRC
  ${CODELETS_ZSRC}
  )

precisions_rules_py(RUNTIME_SRCS_GENERATED "${ZSRC}"
  PRECISIONS "${CHAMELEON_PRECISION}"
  TARGETDIR "codelets")

set(RUNTIME_SRCS
  ${RUNTIME_COMMON}
  ${RUNTIME_SRCS_GENERATED}
  ${CODELETS_SRC}
  )

# Force generation of sources
# ---------------------------
add_custom_target(openmp_sources ALL SOURCES ${RUNTIME_SRCS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};runtime_openmp_include;openmp_sources" CACHE INTERNAL "List of targets of sources")

# Add library
# -----------
add_library(chameleon_openmp ${RUNTIME_SRCS})
set_target_properties(chameleon_openmp PROPERTIES VERSION ${CHAMELEON_VERSION})
set_target_properties(chameleon_openmp PROPERTIES SOVERSION ${CHAMELEON_VERSION_MAJOR})
set_property(TARGET chameleon_openmp PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

add_dependencies(chameleon_openmp
  chameleon_include
  control_include
  coreblas_include
  runtime_openmp_include
  openmp_sources
  )

target_include_directories(chameleon_openmp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)
target_include_directories(chameleon_openmp PRIVATE
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/include>
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/coreblas/include>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/coreblas/include>)

target_link_libraries(chameleon_openmp PUBLIC OpenMP::OpenMP_C)
target_link_libraries(chameleon_openmp PRIVATE coreblas)
target_link_libraries(chameleon_openmp PRIVATE MORSE::M)

# necessary cause of include control/common.h which depends on hqr, mpi, cuda, ...
# TODO: remove this dependency as soon as headers are clean
target_link_libraries(chameleon_openmp PUBLIC hqr)

# export target
install(EXPORT chameleon_openmpTargets
        NAMESPACE CHAMELEON::
        DESTINATION lib/cmake/chameleon
        )

# installation
# ------------
install(TARGETS chameleon_openmp
        EXPORT chameleon_openmpTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )

###
### END CMakeLists.txt
###
