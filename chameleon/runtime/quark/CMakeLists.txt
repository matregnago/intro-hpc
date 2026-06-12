###
#
# @file quark/CMakeLists.txt
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
#  @author Guillaume Sylvand
#  @date 2025-12-19
#
###
cmake_minimum_required(VERSION 3.12)

find_package( QUARK REQUIRED COMPONENTS HWLOC )
morse_export_imported_target(MORSE::QUARK quark chameleon)

# Add definition and include_dir if found
if(QUARK_FOUND)
  message("-- ${Blue}Add definition CHAMELEON_SCHED_QUARK"
    " - Activate QUARK in Chameleon${ColourReset}")
endif(QUARK_FOUND)

# # Generate headers for all possible precisions
# # --------------------------------------------
# set(RUNTIME_HDRS_GENERATED "")
# set(ZHDR
# )
#
# precisions_rules_py(RUNTIME_HDRS_GENERATED "${ZHDR}"
#   PRECISIONS "s;d;c;z;ds;zc"
#   TARGETDIR  "include")

# Define the list of headers
# --------------------------
set(RUNTIME_HDRS
  include/chameleon_quark.h
  include/core_blas_dag.h
  )

# # Add generated headers
# # ---------------------
# foreach( hdr_file ${RUNTIME_HDRS_GENERATED} )
#   list(APPEND RUNTIME_HDRS ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
# endforeach()

# Force generation of headers
# ---------------------------
add_custom_target(
  runtime_quark_include
  ALL SOURCES ${RUNTIME_HDRS})

# Installation
# ------------
install(
  FILES ${RUNTIME_HDRS}
  DESTINATION include/runtime/quark )

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
add_custom_target(quark_sources ALL SOURCES ${RUNTIME_SRCS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};runtime_quark_include;quark_sources" CACHE INTERNAL "List of targets of sources")

# Add library
# -----------
add_library(chameleon_quark ${RUNTIME_SRCS})
set_target_properties(chameleon_quark PROPERTIES VERSION ${CHAMELEON_VERSION})
set_property(TARGET chameleon_quark PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

add_dependencies(chameleon_quark
  chameleon_include
  control_include
  coreblas_include
  runtime_quark_include
  quark_sources
)

target_include_directories(chameleon_quark PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)
target_include_directories(chameleon_quark PRIVATE
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/include>
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/coreblas/include>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/coreblas/include>)

target_link_libraries(chameleon_quark PUBLIC MORSE::QUARK)
target_link_libraries(chameleon_quark PUBLIC coreblas)
target_link_libraries(chameleon_quark PRIVATE MORSE::M)

# necessary cause of include control/common.h which depends on hqr, mpi, cuda, ...
# TODO: remove this dependency as soon as headers are clean
target_link_libraries(chameleon_quark PUBLIC hqr)

  # export target
install(EXPORT chameleon_quarkTargets
        NAMESPACE CHAMELEON::
        DESTINATION lib/cmake/chameleon
        )

# installation
# ------------
install(TARGETS chameleon_quark
        EXPORT chameleon_quarkTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )


###
### END CMakeLists.txt
###
