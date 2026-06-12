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
#  @date 2025-12-19
#
###

if (NOT EZTRACE_FOUND)
    find_package(EZTRACE)
endif()

if (EZTRACE_FOUND AND EZTRACE_DIR_FOUND AND CHAMELEON_USE_EZTRACE)

    set(EZTRACE_eztrace_create_plugin_DIR "EZTRACE_eztrace_create_plugin_DIR-NOTFOUND")
    find_path(EZTRACE_eztrace_create_plugin_DIR
      NAMES eztrace_create_plugin
      HINTS ${EZTRACE_DIR_FOUND}/bin)
    mark_as_advanced(EZTRACE_eztrace_create_plugin_DIR)

    if (EZTRACE_eztrace_create_plugin_DIR)

        set(EZTRACE_CREATE_PLUGIN "${EZTRACE_eztrace_create_plugin_DIR}/eztrace_create_plugin")

        add_custom_command(
            OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output
            COMMAND ${EZTRACE_CREATE_PLUGIN}
            ARGS ${CMAKE_CURRENT_SOURCE_DIR}/gpuhipblas_eztrace_module
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gpuhipblas_eztrace_module
            )
        add_custom_target(
            eztrace-module-chameleon_hip-dir ALL
            DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output
            )
        add_custom_command(
            OUTPUT libeztrace-convert-chameleon_hip.so
            COMMAND make
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output
            DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output
            )
        add_custom_target(
            eztrace-module-chameleon_hip-libs ALL
            DEPENDS libeztrace-convert-chameleon_hip.so
            )
        # installation
        # ------------
        install(
            FILES
            ${CMAKE_CURRENT_BINARY_DIR}/output/libeztrace-autostart-chameleon_hip.so
            ${CMAKE_CURRENT_BINARY_DIR}/output/libeztrace-chameleon_hip.so
            ${CMAKE_CURRENT_BINARY_DIR}/output/libeztrace-convert-chameleon_hip.so
            DESTINATION lib
            )

    endif (EZTRACE_eztrace_create_plugin_DIR)

endif ()

###
### END CMakeLists.txt
###
