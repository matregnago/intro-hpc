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
#  @author Brieuc Nicolas
#  @date 2025-12-19
#
###
if (EZTRACE_FOUND AND CHAMELEON_USE_EZTRACE)
    if( EZTRACE_VERSION VERSION_GREATER "2.0.0" )

        set( COREBLAS_EZTRACE_MODULES core tcore)
        foreach( _module ${COREBLAS_EZTRACE_MODULES} )

            set( MODULE_SRCS ${_module}_z.c ${_module}_zc.c)
            precisions_rules_py(MODULE_SRCS_GENERATED "${MODULE_SRCS}"
                    PRECISIONS "${CHAMELEON_PRECISION}")

            add_library(eztrace-coreblas_${_module} SHARED ${MODULE_SRCS_GENERATED} ${_module}.c)

            target_link_libraries(eztrace-coreblas_${_module} PUBLIC MORSE::EZTRACE)

            add_dependencies(eztrace-coreblas_${_module} coreblas_include)

            target_include_directories(eztrace-coreblas_${_module} PUBLIC
                $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/coreblas/include>
                $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/coreblas/include>
                $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/include>
                $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/include>
            )

            # Install the shared library
            install(
                TARGETS eztrace-coreblas_${_module}
                LIBRARY DESTINATION lib
            )
            unset(MODULE_SRCS_GENERATED)
        endforeach()

    else()
        set(EZTRACE_eztrace_create_plugin_DIR "EZTRACE_eztrace_create_plugin_DIR-NOTFOUND")
        find_path(EZTRACE_eztrace_create_plugin_DIR
        NAMES eztrace_create_plugin
        HINTS ${EZTRACE_DIR_FOUND}/bin)
        mark_as_advanced(EZTRACE_eztrace_create_plugin_DIR)

        if (EZTRACE_eztrace_create_plugin_DIR)

            set(EZTRACE_CREATE_PLUGIN "${EZTRACE_eztrace_create_plugin_DIR}/eztrace_create_plugin")

            set( COREBLAS_EZTRACE_MODULES core tcore )

            foreach( _module ${COREBLAS_EZTRACE_MODULES} )
                add_custom_command(
                    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_module}
                    COMMAND ${EZTRACE_CREATE_PLUGIN} -o ${CMAKE_CURRENT_BINARY_DIR}/${_module} > /dev/null
                    ARGS ${CMAKE_CURRENT_SOURCE_DIR}/coreblas_${_module}_eztrace_module
                    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/coreblas_${_module}_eztrace_module
                    )
                add_custom_target(
                    eztrace-module-chameleon_${_module}-dir ALL
                    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${_module}
                    )
                add_custom_command(
                    OUTPUT libeztrace-convert-chameleon_${_module}.so
                    COMMAND make
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_module}
                    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${_module}
                    )
                add_custom_target(
                    eztrace-module-chameleon_${_module}-libs ALL
                    DEPENDS libeztrace-convert-chameleon_${_module}.so
                    )
                # installation
                # ------------
                install(
                    FILES
                    ${CMAKE_CURRENT_BINARY_DIR}/${_module}/libeztrace-autostart-chameleon_${_module}.so
                    ${CMAKE_CURRENT_BINARY_DIR}/${_module}/libeztrace-chameleon_${_module}.so
                    ${CMAKE_CURRENT_BINARY_DIR}/${_module}/libeztrace-convert-chameleon_${_module}.so
                    LIBRARY DESTINATION lib
                    )
            endforeach()
        endif (EZTRACE_eztrace_create_plugin_DIR)
    endif()
endif (EZTRACE_FOUND AND CHAMELEON_USE_EZTRACE)

###
### END CMakeLists.txt
###
