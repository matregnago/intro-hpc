###
#
# @file CMakeLists.txt
#
# @copyright (c) 2017 Inria. All rights reserved.
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
#  @date 2025-12-19
#
###

cmake_minimum_required(VERSION 3.12)

# Create files in binary dir
# --------------------------
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/homepage.org
               ${CMAKE_CURRENT_BINARY_DIR}/homepage.org
               COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/news.org
               ${CMAKE_CURRENT_BINARY_DIR}/news.org
               COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/publish.el
               ${CMAKE_CURRENT_BINARY_DIR}/publish.el
               COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.org.in
               ${CMAKE_CURRENT_BINARY_DIR}/version.org
               @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/users_guide.org.in
               ${CMAKE_CURRENT_BINARY_DIR}/users_guide.org
               @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/CONTRIBUTING.org
               ${CMAKE_CURRENT_BINARY_DIR}/CONTRIBUTING.org
               COPYONLY)
set(FIGURES_HOMEPAGE
    chameleon.svg
    )
foreach(_fig ${FIGURES_HOMEPAGE})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${_fig}
                   ${CMAKE_CURRENT_BINARY_DIR}/${_fig}
                   COPYONLY)
endforeach()

set(FIGURES_USERGUIDE
    tile_lu.jpg
    tile_layout.jpg
    trace_qr.jpg
    potri_async.png
    chameleon_header.png
    lstopo-sirocco24.png
    )
foreach(_fig ${FIGURES_USERGUIDE})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/chapters/${_fig}
                   ${CMAKE_CURRENT_BINARY_DIR}/${_fig}
                   COPYONLY)
endforeach()

# Looking for emacs
# -----------------
find_program(EMACS_COMPILER emacs REQUIRED)

if(EMACS_COMPILER)
    message(STATUS "Looking for emacs - found")

    # Add target for homepage
    # -----------------------
    add_custom_command(OUTPUT  index.html
                       COMMAND ${EMACS_COMPILER}
                       ARGS    --batch --no-init-file --load publish.el --funcall org-publish-all
                       DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/homepage.org
                               ${CMAKE_CURRENT_BINARY_DIR}/news.org
                               ${CMAKE_CURRENT_BINARY_DIR}/version.org
                               ${CMAKE_CURRENT_BINARY_DIR}/users_guide.org
                               ${CMAKE_CURRENT_SOURCE_DIR}/chapters/introduction.org
                               ${CMAKE_CURRENT_SOURCE_DIR}/chapters/installing.org
                               ${CMAKE_CURRENT_SOURCE_DIR}/chapters/using.org
                               ${CMAKE_CURRENT_SOURCE_DIR}/chapters/parallel_worker.org
                               ${CMAKE_CURRENT_SOURCE_DIR}/chapters/performances.org
                               ${CMAKE_CURRENT_BINARY_DIR}/CONTRIBUTING.org
                               ${CMAKE_CURRENT_BINARY_DIR}/publish.el
                               ${CMAKE_CURRENT_BINARY_DIR}/chameleon.svg)

    add_custom_target(doc-homepage ALL DEPENDS index.html)

    # Installation
    # ------------
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/index.html
            DESTINATION share/chameleon/doc/)
    foreach(_fig ${FIGURES_HOMEPAGE} ${FIGURES_USERGUIDE})
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${_fig}
                DESTINATION share/chameleon/doc)
    endforeach()

else(EMACS_COMPILER)
    message(STATUS "Looking for emacs - not found")
endif(EMACS_COMPILER)

###
### END CMakeLists.txt
###
