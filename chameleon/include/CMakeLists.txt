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
#  @date 2024-03-16
#
###

# Generate the chameleon headers for all possible precisions
# ----------------------------------------------------------
set(CHAMELEON_HDRS_GENERATED "")
set(ZHDR
  chameleon/chameleon_z.h
  chameleon/chameleon_zc.h
  chameleon/tasks_z.h
  chameleon/tasks_zc.h
)

precisions_rules_py(CHAMELEON_HDRS_GENERATED "${ZHDR}"
  TARGETDIR  "chameleon"
  PRECISIONS "z;c;d;s;zc;ds" )

configure_file("chameleon/config.h.in"
               "chameleon/config.h"
               @ONLY)

# Define the list of headers
# --------------------------
set(CHAMELEON_HDRS
  ${CMAKE_CURRENT_BINARY_DIR}/chameleon/config.h
  ${CMAKE_CURRENT_BINARY_DIR}/chameleon/mangling.h
  chameleon/constants.h
  chameleon/runtime.h
  chameleon/runtime_struct.h
  chameleon/simulate.h
  chameleon/struct.h
  chameleon/struct_context.h
  chameleon/tasks.h
  chameleon/timer.h
  chameleon/types.h
  chameleon/flops.h
  chameleon/descriptor_helpers.h
  )

# Add generated headers
# ---------------------
foreach( hdr_file ${CHAMELEON_HDRS_GENERATED} )
  list(APPEND CHAMELEON_HDRS ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
endforeach()

# Force generation of headers
# ---------------------------
add_custom_target(chameleon_include
  ALL SOURCES chameleon.h ${CHAMELEON_HDRS} )
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};chameleon_include" CACHE INTERNAL "List of targets of sources")

# Installation
# ------------
install( FILES chameleon.h
         DESTINATION include )

install( FILES ${CHAMELEON_HDRS}
         DESTINATION include/chameleon )

###
### END CMakeLists.txt
###
