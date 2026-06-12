###
#
# @file CMakeLists.txt
#
# @copyright 2009-2014 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
#  @author Florent Pruvost
#  @author Guillaume Sylvand
#  @author Mathieu Faverge
#  @author Loris Lucido
#  @date 2025-12-19
#
###

# Generate the chameleon sources for all possible precisions
# ------------------------------------------------------
set(GPUHIPBLAS_SRCS_GENERATED "")
set(ZSRC
    hip_zgemm.c
    hip_zhemm.c
    hip_zher2k.c
    hip_zherk.c
    hip_zpotrf.c
    hip_zsymm.c
    hip_zsyr2k.c
    hip_zsyrk.c
    hip_ztrmm.c
    hip_ztrsm.c
    )

precisions_rules_py(
  GPUHIPBLAS_SRCS_GENERATED "${ZSRC}"
  PRECISIONS "${CHAMELEON_PRECISION}")

set(GPUHIPBLAS_SRCS
  ${GPUHIPBLAS_SRCS_GENERATED}
  hip_hgemm.c
  hipglobal.c
  )

# Force generation of sources
# ---------------------------
add_custom_target(gpuhipblas_sources ALL SOURCES ${GPUHIPBLAS_SRCS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};gpuhipblas_sources" CACHE INTERNAL "List of targets of sources")

# Compile step
# ------------
add_library(gpuhipblas ${GPUHIPBLAS_SRCS})
set_target_properties(gpuhipblas PROPERTIES VERSION ${CHAMELEON_VERSION})
set_target_properties(gpuhipblas PROPERTIES SOVERSION ${CHAMELEON_VERSION_MAJOR})
add_dependencies(gpuhipblas gpuhipblas_include gpuhipblas_sources)
target_include_directories(gpuhipblas PUBLIC
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/gpuhipblas/include>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/gpuhipblas/include>
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)
set_property(TARGET gpuhipblas PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

target_link_libraries(gpuhipblas PRIVATE coreblas)
target_link_libraries(gpuhipblas PUBLIC roc::hipsolver)
target_link_libraries(gpuhipblas PUBLIC roc::hipblas)
target_link_libraries(gpuhipblas PUBLIC MORSE::M)

# export target coreblas
install(EXPORT gpuhipblasTargets
        NAMESPACE CHAMELEON::
        DESTINATION lib/cmake/chameleon
        )

# installation
# ------------
install(TARGETS gpuhipblas
        EXPORT gpuhipblasTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )

###
### END CMakeLists.txt
###
