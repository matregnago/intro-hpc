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

if (CHAMELEON_USE_HMATOSS AND NOT CHAMELEON_HMAT_EXTERNAL)
  add_subdirectory(hmat-oss)
endif()
add_subdirectory(include)
if (NOT CHAMELEON_SIMULATION)
  add_subdirectory(compute)
  add_subdirectory(eztrace_module)
endif()

###
### END CMakeLists.txt
###
