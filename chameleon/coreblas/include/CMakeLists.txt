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
#  @author Abel Calluaud
#  @date 2024-03-27
#
###

# Generate header files
# ---------------------
set(COREBLAS_HDRS_GENERATED "")
set(ZHDR
  coreblas/coreblas_z.h
  coreblas/coreblas_zc.h
  coreblas/coreblas_ztile.h
  coreblas/coreblas_zctile.h
)
if( CHAMELEON_USE_HMATOSS )
  list( APPEND ZHDR
    coreblas/hmat_z.h )
endif()

precisions_rules_py(
  COREBLAS_HDRS_GENERATED "${ZHDR}"
  TARGETDIR  coreblas
  PRECISIONS "s;d;c;z;zc;ds" )

# Define the list of headers
# --------------------------
set(COREBLAS_HDRS
    coreblas.h
    coreblas/cblas_wrapper.h
    coreblas/cblas.h
    coreblas/lapacke.h
    coreblas/lapacke_config.h
    coreblas/lapacke_mangling.h
    coreblas/random.h
    )
if( CHAMELEON_USE_HMATOSS )
  list( APPEND COREBLAS_HDRS
    coreblas/hmat.h )
endif()

# Add generated headers
# ---------------------
foreach( hdr_file ${COREBLAS_HDRS_GENERATED} )
  list(APPEND COREBLAS_HDRS ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
endforeach()

# Force generation of headers
# ---------------------------
add_custom_target(coreblas_include ALL SOURCES ${COREBLAS_HDRS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};coreblas_include" CACHE INTERNAL "List of targets of sources")

# Installation
# ------------
install( FILES coreblas.h
         DESTINATION include )

install( FILES ${COREBLAS_HDRS}
         DESTINATION include/coreblas )

###
### END CMakeLists.txt
###
