###
#
# @file CMakeLists.txt
#
# @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 2025-01-22
#
###

# Generate the chameleon headers for all possible precisions
# ----------------------------------------------------------
set(CHAMELEON_LAPACK_HDRS_GENERATED "")
set(ZHDR
  include/chameleon/chameleon_zlapack.h
)

precisions_rules_py(CHAMELEON_LAPACK_HDRS_GENERATED "${ZHDR}"
  TARGETDIR  "include/chameleon"
  PRECISIONS "z;c;d;s" )

set(CHAMELEON_LAPACK_HDRS "")
foreach( hdr_file ${CHAMELEON_LAPACK_HDRS_GENERATED} )
  list(APPEND CHAMELEON_LAPACK_HDRS ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
endforeach()

# Force generation of headers
# ---------------------------
add_custom_target(chameleon_lapack_include
  ALL SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/include/chameleon_lapack.h
              ${CHAMELEON_LAPACK_HDRS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};chameleon_lapack_include" CACHE INTERNAL "List of targets of sources")

# Installation of headers
# -----------------------
install( FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/chameleon_lapack.h
         DESTINATION include )
install( FILES ${CHAMELEON_LAPACK_HDRS}
         DESTINATION include/chameleon )

# Define the list of sources
# --------------------------
set(LAPACK_SRCS_GENERATED "")
set(ZSRC
  src/lapack_zgemm.c
  src/lapack_zgetrf.c
  src/lapack_zhemm.c
  src/lapack_zher2k.c
  src/lapack_zherk.c
  src/lapack_zlacpy.c
  src/lapack_zlange.c
  src/lapack_zlanhe.c
  src/lapack_zlansy.c
  src/lapack_zlantr.c
  src/lapack_zlaset.c
  src/lapack_zlauum.c
  src/lapack_zposv.c
  src/lapack_zpotrf.c
  src/lapack_zpotri.c
  src/lapack_zpotrs.c
  src/lapack_zsymm.c
  src/lapack_zsyr2k.c
  src/lapack_zsyrk.c
  src/lapack_ztrmm.c
  src/lapack_ztrsm.c
  src/lapack_ztrtri.c
  )
precisions_rules_py(LAPACK_SRCS_GENERATED "${ZSRC}"
  PRECISIONS "${CHAMELEON_PRECISION}")
set(CHAMELEON_LAPACK_SRCS
  ${LAPACK_SRCS_GENERATED}
  src/lapack_api_common.c
  )

# Force generation of sources
# ---------------------------
add_custom_target(chameleon_lapack_sources ALL SOURCES "${CHAMELEON_LAPACK_SRCS}")
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};chameleon_lapack_sources" CACHE INTERNAL "List of targets of sources")

# Define library
# --------------
add_library(chameleon_lapack ${CHAMELEON_LAPACK_SRCS})
set_target_properties(chameleon_lapack PROPERTIES VERSION ${CHAMELEON_VERSION})
set_target_properties(chameleon_lapack PROPERTIES SOVERSION ${CHAMELEON_VERSION_MAJOR})

# Define dependencies
# -------------------
add_dependencies(chameleon_lapack chameleon_lapack_include)
target_include_directories(chameleon_lapack PUBLIC
  $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/lapack_api/include>
  $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/lapack_api/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(chameleon_lapack PUBLIC chameleon)
if (CHAMELEON_SIMULATION)
  target_include_directories(chameleon_lapack PUBLIC $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/coreblas/include>)
else()
  target_link_libraries(chameleon_lapack PUBLIC coreblas)
endif()
set_property(TARGET chameleon_lapack PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

# Export target
# -------------
install(EXPORT chameleon_lapackTargets
        NAMESPACE CHAMELEON::
        DESTINATION lib/cmake/chameleon
        )

# Installation of library and cmake target
# ----------------------------------------
install(TARGETS chameleon_lapack
        EXPORT chameleon_lapackTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )

###
### END CMakeLists.txt
###
