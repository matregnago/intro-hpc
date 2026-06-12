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
#  CHAMELEON example routines
#  CHAMELEON is a software package provided by Inria Bordeaux - Sud-Ouest, LaBRI,
#  University of Bordeaux, Bordeaux INP
#
# @version 1.4.0
#  @author Florent Pruvost
#  @author Guillaume Sylvand
#  @author Mathieu Faverge
#  @date 2024-02-18
#
###

# sources depending on CBLAS/LAPACKE
set(LAPACK_SOURCES
    step0.c
   )
foreach(_src ${LAPACK_SOURCES})
    get_filename_component(_name_exe ${_src} NAME_WE)
    add_executable(${_name_exe} ${_src})
    target_link_libraries(${_name_exe} PRIVATE
      MORSE::LAPACKE
      MORSE::CBLAS)
    install(TARGETS ${_name_exe}
            DESTINATION bin/chameleon/lapack_to_chameleon)
    # for chameleon/include/chameleon/timer.h
    target_include_directories(${_name_exe} PRIVATE
      ${CHAMELEON_SOURCE_DIR}/include)
    # for chameleon/coreblas/include/coreblas/cblas.h
    target_include_directories(${_name_exe} PRIVATE
      ${CHAMELEON_SOURCE_DIR}/coreblas/include)
endforeach()

# sources depending on chameleon
set(CHAMELEON_SOURCES
    step1.c
    step2.c
    step3.c
    step4.c
    step5.c
    step6.c
    step7.c
)

foreach(_src ${CHAMELEON_SOURCES})
    get_filename_component(_name_exe ${_src} NAME_WE)
    add_executable(${_name_exe} ${_src})
    target_link_libraries(${_name_exe} PRIVATE
      chameleon
      coreblas
      MORSE::LAPACKE)
    # for chameleon/include/chameleon/timer.h
    target_include_directories(${_name_exe} PRIVATE
      ${CHAMELEON_SOURCE_DIR}/include)
    # for chameleon/coreblas/include/coreblas/cblas.h
    target_include_directories(${_name_exe} PRIVATE
      ${CHAMELEON_SOURCE_DIR}/coreblas/include
      ${CHAMELEON_BINARY_DIR}/coreblas/include)
    install(TARGETS ${_name_exe}
            DESTINATION bin/chameleon/lapack_to_chameleon)
endforeach()

#-------- Tests ---------
include(CTestLists.cmake)

###
### END CMakeLists.txt
###
