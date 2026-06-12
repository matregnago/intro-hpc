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
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @author Loris Lucido
#  @date 2024-02-18
#
###

# Generate header files
# ---------------------
set(GPUHIPBLAS_HDRS_GENERATED "")
set(ZHDR
    gpuhipblas/gpuhipblas_z.h
)
precisions_rules_py(
  GPUHIPBLAS_HDRS_GENERATED "${ZHDR}"
  TARGETDIR gpuhipblas
  PRECISIONS "s;d;c;z;zc;ds" )

# Define the list of headers
# --------------------------
set(GPUHIPBLAS_HDRS
    gpuhipblas.h
    )

# Add generated headers
# ---------------------
foreach( hdr_file ${GPUHIPBLAS_HDRS_GENERATED} )
  list(APPEND GPUHIPBLAS_HDRS ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
endforeach()

# Force generation of headers
# ---------------------------
add_custom_target(gpuhipblas_include ALL SOURCES ${GPUHIPBLAS_HDRS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};gpuhipblas_include" CACHE INTERNAL "List of targets of sources")

# Installation
# ------------
install( FILES gpuhipblas.h
         DESTINATION include )

install( FILES ${GPUHIPBLAS_HDRS}
         DESTINATION include/gpuhipblas )

###
### END CMakeLists.txt
###
