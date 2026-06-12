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
        set( MODULE_SRCS gpucublas_z.c)

        if ( GPUCUBLAS_HAVE_CUDA_TOOLKIT )
            list(APPEND MODULE_SRCS gpucublas_zc.c)
        endif()

        precisions_rules_py(MODULE_SRCS_GENERATED "${MODULE_SRCS}"
        PRECISIONS "${CHAMELEON_PRECISION}")

        if ( GPUCUBLAS_HAVE_CUDA_HALF OR GPUCUBLAS_HAVE_CUBLASGEMMEX )
            list(APPEND MODULE_SRCS_GENERATED gpucublas_h.c)
        endif()

        add_library(eztrace-gpucublas SHARED ${MODULE_SRCS_GENERATED} gpucublas.c)
        add_dependencies(eztrace-gpucublas
            gpucublas_include)

        target_link_libraries(eztrace-gpucublas PUBLIC MORSE::EZTRACE)
        target_link_libraries(eztrace-gpucublas PRIVATE coreblas CUDA::cusolver CUDA::cublas)

        target_include_directories(eztrace-gpucublas PUBLIC
            $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/gpucublas/include>
            $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/gpucublas/include>
            $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}/include>
            $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}/include>
            $<INSTALL_INTERFACE:include>)

        # Install the shared library
        install(
            TARGETS eztrace-gpucublas
            LIBRARY DESTINATION lib
        )
    else()
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
                ARGS ${CMAKE_CURRENT_SOURCE_DIR}/gpucublas_eztrace_module
                DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gpucublas_eztrace_module
                )
            add_custom_target(
                eztrace-module-chameleon_cuda-dir ALL
                DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output
                )
            add_custom_command(
                OUTPUT libeztrace-convert-chameleon_cuda.so
                COMMAND make
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output
                DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output
                )
            add_custom_target(
                eztrace-module-chameleon_cuda-libs ALL
                DEPENDS libeztrace-convert-chameleon_cuda.so
                )
            # installation
            # ------------
            install(
                FILES
                ${CMAKE_CURRENT_BINARY_DIR}/output/libeztrace-autostart-chameleon_cuda.so
                ${CMAKE_CURRENT_BINARY_DIR}/output/libeztrace-chameleon_cuda.so
                ${CMAKE_CURRENT_BINARY_DIR}/output/libeztrace-convert-chameleon_cuda.so
                DESTINATION lib
                )
        endif (EZTRACE_eztrace_create_plugin_DIR)
    endif()
endif ()

###
### END CMakeLists.txt
###
