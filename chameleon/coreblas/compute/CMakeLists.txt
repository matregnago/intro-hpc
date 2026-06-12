###
#
# @file CMakeLists.txt
#
# @copyright 2009-2014 The University of Tennessee and The University of
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
#  @author Matthieu Kuhn
#  @author Abel Calluaud
#  @date 2025-12-19
#
###

option( CHAMELEON_NO_KERNELS "Do not call the numerical kernels" OFF)

# Generate the chameleon sources for all possible precisions
# ------------------------------------------------------

set(COREBLAS_SRCS_GENERATED "")
set(ZSRC
    core_dlag2z.c
    core_dzasum.c
    core_zaxpy.c
    core_zcesca.c
    core_zgeadd.c
    core_zlascal.c
    core_zlatm1.c
    core_zgelqt.c
    core_zgemv.c
    core_zgemm.c
    core_zgeqrt.c
    core_zgesplit.c
    core_zgesum.c
    core_zgessm.c
    core_zgessq.c
    core_zgetf2_nopiv.c
    core_zgetrf.c
    core_zgetrf_incpiv.c
    core_zgetrf_nopiv.c
    core_zgetrf_panel.c
    core_zgram.c
    core_zhe2ge.c
    core_zherfb.c
    core_zhemm.c
    core_zher2k.c
    core_zherk.c
    core_zhessq.c
    core_zlacpy.c
    core_zlag2c.c
    core_zlange.c
    core_zlanhe.c
    core_zlansy.c
    core_zlantr.c
    core_zlaset.c
    core_zlaset2.c
    core_zlaswp.c
    core_zlatro.c
    core_zlauum.c
    core_zpamm.c
    core_zparfb.c
    core_zpemv.c
    core_zplghe.c
    core_zplgsy.c
    core_zplrnt.c
    core_zplssq.c
    core_zpotrf.c
    core_zssssm.c
    core_zsymm.c
    core_zsyr2k.c
    core_zsyrk.c
    core_zsyssq.c
    core_zsytf2_nopiv.c
    core_ztplqt.c
    core_ztpqrt.c
    core_ztpmlqt.c
    core_ztpmqrt.c
    core_ztradd.c
    core_ztrasm.c
    core_ztrmm.c
    core_ztrsm.c
    core_ztrssq.c
    core_ztrtri.c
    core_ztslqt.c
    core_ztsmlq.c
    core_ztsmqr.c
    core_ztsmlq_hetra1.c
    core_ztsmqr_hetra1.c
    core_ztsqrt.c
    core_ztstrf.c
    core_zttlqt.c
    core_zttmlq.c
    core_zttmqr.c
    core_zttqrt.c
    core_zunmlq.c
    core_zunmqr.c
    core_zprint.c
    )
if( CHAMELEON_USE_HMATOSS )
  list( APPEND ZSRC
    hmat_z.c )
endif()

if ( CHAMELEON_NO_KERNELS )
  list( APPEND ZSRC
    core_ztile_empty.c
    core_zctile_empty.c
    )
else()
  list( APPEND ZSRC
    core_ztile.c
    core_zctile.c
    )
endif()

precisions_rules_py(COREBLAS_SRCS_GENERATED "${ZSRC}"
                    PRECISIONS "${CHAMELEON_PRECISION}")

set(COREBLAS_SRCS
  global.c
  core_ipiv_to_perm.c
  ${COREBLAS_SRCS_GENERATED}
)

# Force generation of sources
# ---------------------------
add_custom_target(coreblas_sources ALL SOURCES ${COREBLAS_SRCS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};coreblas_sources" CACHE INTERNAL "List of targets of sources")

# Compile step
# ------------
add_library(coreblas ${COREBLAS_SRCS})
set_target_properties(coreblas PROPERTIES VERSION ${CHAMELEON_VERSION})
set_target_properties(coreblas PROPERTIES SOVERSION ${CHAMELEON_VERSION_MAJOR})
add_dependencies(coreblas coreblas_include coreblas_sources)
target_include_directories(coreblas PUBLIC
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/coreblas/include>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/coreblas/include>
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)
set_property(TARGET coreblas PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

if( CHAMELEON_USE_HMATOSS )
  target_link_libraries(coreblas PUBLIC HMAT::hmat )
endif()

target_link_libraries(coreblas PRIVATE MORSE::LAPACKE)
target_link_libraries(coreblas PRIVATE MORSE::CBLAS)
target_link_libraries(coreblas PUBLIC MORSE::M)

# export target coreblas
install(EXPORT coreblasTargets
        NAMESPACE CHAMELEON::
        DESTINATION lib/cmake/chameleon
        )

# Installation
# ------------
install(TARGETS coreblas
        EXPORT coreblasTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )

###
### END CMakeLists.txt
###
