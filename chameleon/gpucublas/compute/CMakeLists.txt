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
#  @author Brieuc Nicolas
#  @date 2025-12-19
#
###

# To define CMAKE_CUDA_COMPILER
cmake_minimum_required(VERSION 3.18)

# Generate the chameleon sources for all possible precisions
# ------------------------------------------------------
set(GPUCUBLAS_SRCS_GENERATED "")
set(ZSRC
    cuda_zgeadd.c
    cuda_zgemerge.c
    cuda_zgemm.c
    cuda_zhemm.c
    cuda_zher2k.c
    cuda_zherfb.c
    cuda_zherk.c
    cuda_zlarfb.c
    cuda_zparfb.c
    cuda_zpotrf.c
    cuda_zsymm.c
    cuda_zsyr2k.c
    cuda_zsyrk.c
    cuda_ztpmlqt.c
    cuda_ztpmqrt.c
    cuda_ztrmm.c
    cuda_ztrsm.c
    cuda_ztsmlq.c
    cuda_ztsmqr.c
    cuda_zttmlq.c
    cuda_zttmqr.c
    cuda_zunmlqt.c
    cuda_zunmqrt.c
    )

if ( CHAMELEON_NO_KERNELS )
  list( APPEND ZSRC
    cuda_ztile_empty.c
    cuda_zctile_empty.c
    )
else()
  list( APPEND ZSRC
    cuda_ztile.c
    cuda_zctile.c
    )
endif()

if ( GPUCUBLAS_HAVE_CUDA_TOOLKIT )
  set(ZSRC
    ${ZSRC}
    cuda_zlag2c.cu
    cuda_zlatro.cu
    cuda_zlacpy.cu
    cuda_zlaset.cu
  )
endif()
if ( GPUCUBLAS_HAVE_CUDA_HALF )
  set(ZSRC
    ${ZSRC}
    cuda_dlag2h.cu
  )
endif()

# Former MAGMA files that are no longer supported
# if( CHAMELEON_USE_MAGMA )
#   set(ZSRC
#     ${ZSRC}
#     cuda_zgelqt.c
#     cuda_zgeqrt.c
#     cuda_zgessm.c
#     cuda_zgetrf.c
#     cuda_zlauum.c
#     cuda_zpotrf.c
#     cuda_zssssm.c
#     cuda_ztrtri.c
#     cuda_ztslqt.c
#     cuda_ztsqrt.c
#     cuda_ztstrf.c
#     )
# endif()

precisions_rules_py(
  GPUCUBLAS_SRCS_GENERATED "${ZSRC}"
  PRECISIONS "${CHAMELEON_PRECISION}")

set(GPUCUBLAS_SRCS
  ${GPUCUBLAS_SRCS_GENERATED}
  cudaglobal.c
)

if (GPUCUBLAS_HAVE_CUBLASHGEMM)
  set(GPUCUBLAS_SRCS
    ${GPUCUBLAS_SRCS}
    cuda_hgemm.c
  )
  # Need to use CXX compiler to have the __half support and access to cublasHgemm()
  set_source_files_properties( cuda_hgemm.c PROPERTIES LANGUAGE CXX )
endif()

if (GPUCUBLAS_HAVE_CUBLASGEMMEX)
  set(GPUCUBLAS_SRCS
    ${GPUCUBLAS_SRCS}
    cuda_gemmex.c
  )
endif()

if ( GPUCUBLAS_HAVE_CUBLASHGEMM OR GPUCUBLAS_HAVE_CUBLASGEMMEX )
  if ( CHAMELEON_NO_KERNELS )
    list( APPEND GPUCUBLAS_SRCS
      cuda_htile_empty.c
      )
  else()
    list( APPEND GPUCUBLAS_SRCS
      cuda_htile.c
      )
  endif()
endif()

# Force generation of sources
# ---------------------------
add_custom_target(gpucublas_sources ALL SOURCES ${GPUCUBLAS_SRCS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};gpucublas_sources" CACHE INTERNAL "List of targets of sources")

# Compile step
# ------------
add_library(gpucublas ${GPUCUBLAS_SRCS})
set_target_properties(gpucublas PROPERTIES VERSION ${CHAMELEON_VERSION})
set_target_properties(gpucublas PROPERTIES SOVERSION ${CHAMELEON_VERSION_MAJOR})
add_dependencies(gpucublas gpucublas_include gpucublas_sources)
target_include_directories(gpucublas PUBLIC
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/gpucublas/include>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/gpucublas/include>
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>)
set_property(TARGET gpucublas PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

target_link_libraries(gpucublas PRIVATE coreblas CUDA::cusolver CUDA::cublas)
target_link_libraries(gpucublas PUBLIC MORSE::M)

set_target_properties(gpucublas PROPERTIES
  CUDA_SEPARABLE_COMPILATION OFF)

#target_include_directories( gpucublas PRIVATE ${CUDAToolkit_INCLUDE_DIRS})
#target_link_libraries( gpucublas PRIVATE CUDA::cusolver CUDA::cublas CUDA::cudart )

# export target coreblas
install(EXPORT gpucublasTargets
        NAMESPACE CHAMELEON::
        DESTINATION lib/cmake/chameleon
        )

# installation
# ------------
install(TARGETS gpucublas
        EXPORT gpucublasTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )

###
### END CMakeLists.txt
###
