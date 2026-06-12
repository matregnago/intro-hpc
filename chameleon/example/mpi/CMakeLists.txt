###
#
# @file CMakeLists.txt
#
# @copyright 2024-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
#  @author Mathieu Faverge
#  @date 2024-05-06
#
###

set(MPICMD mpiexec --bind-to none -n 4)

set(SOURCES
    comm_split.c
   )

foreach(_src ${SOURCES})
  get_filename_component(_name_exe ${_src} NAME_WE)
  add_executable(${_name_exe} ${_src})
  target_link_libraries(${_name_exe} PRIVATE chameleon coreblas MORSE::LAPACKE)
  install(TARGETS ${_name_exe} DESTINATION bin/chameleon/mpi)
  add_test(example_mpi_${_name_exe} ${MPICMD} ./${_name_exe})
endforeach()

###
### END CMakeLists.txt
###
