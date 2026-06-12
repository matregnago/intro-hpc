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
#  @date 2025-12-19
#
###

cmake_minimum_required(VERSION 3.12)

# Looking for doxygen
# -------------------
find_package(Doxygen REQUIRED)

# Generate documentation
# ----------------------
if(DOXYGEN_EXECUTABLE)

    # Create file chameleon.dox
    # -------------------------
    configure_file("conf.dox.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/chameleon.dox"
                   @ONLY)

    # Message during configure
    #-------------------------
    message(STATUS "Looking for doxygen - found")

    # Doxygen documentation
    # ---------------------
    add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/index.html
                       COMMAND ${DOXYGEN_EXECUTABLE}
                       ARGS    ${CMAKE_CURRENT_BINARY_DIR}/chameleon.dox
                       DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/chameleon.dox
                               ${CMAKE_CURRENT_SOURCE_DIR}/main.dox
                      )
    add_custom_target(doxygen-out ALL
                      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/index.html
                              chameleon_all_sources
                              )

    # Installation for html version
    # -----------------------------
    install(DIRECTORY   ${CMAKE_CURRENT_BINARY_DIR}/html
            DESTINATION share/chameleon/doc/dev)

else()
    message(STATUS "Looking for doxygen - not found")

endif()

###
### END CMakeLists.txt
###
