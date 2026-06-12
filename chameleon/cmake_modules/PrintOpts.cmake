###
#
# @file PrintOpts.cmake
#
# @copyright 2009-2014 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2026 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
#  @author Abel Calluaud
#  @author Brieuc Nicolas
#  @date 2025-12-19
#
###

if (BLAS_LIBRARIES MATCHES "mkl_core")
  set(BLAS_VENDOR_FOUND "Intel")
elseif(BLAS_LIBRARIES MATCHES "openblas")
  set(BLAS_VENDOR_FOUND "OpenBLAS")
elseif(BLAS_LIBRARIES MATCHES "blis")
  set(BLAS_VENDOR_FOUND "BLIS")
elseif(BLAS_LIBRARIES MATCHES "Accelerate")
  set(BLAS_VENDOR_FOUND "Apple")
elseif(BLAS_LIBRARIES MATCHES "essl")
  set(BLAS_VENDOR_FOUND "IBMESSL")
elseif(BLAS_LIBRARIES MATCHES "blas")
  set(BLAS_VENDOR_FOUND "Generic")
endif()
if (LAPACK_LIBRARIES MATCHES "mkl_core")
  set(LAPACK_VENDOR_FOUND "Intel")
elseif(LAPACK_LIBRARIES MATCHES "openblas")
  set(LAPACK_VENDOR_FOUND "OpenBLAS")
elseif(LAPACK_LIBRARIES MATCHES "flame")
  set(LAPACK_VENDOR_FOUND "FLAME")
elseif(LAPACK_LIBRARIES MATCHES "Accelerate")
  set(LAPACK_VENDOR_FOUND "Apple")
elseif(LAPACK_LIBRARIES MATCHES "lapack")
  set(LAPACK_VENDOR_FOUND "Generic")
endif()

set(dep_message "\nConfiguration of Chameleon:\n"
        "       BUILDNAME ...........: ${BUILDNAME}\n"
        "       SITE ................: ${SITE}\n"
        "\n"
        "       Compiler: C .........: ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})\n"
        "       Compiler: Fortran ...: ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})\n")
if(CHAMELEON_USE_MPI)
  set(dep_message "${dep_message}"
  "       Compiler: MPI .......: ${MPI_C_COMPILER}\n"
  "       compiler flags ......: ${MPI_C_COMPILE_FLAGS}\n")
endif()
set(dep_message "${dep_message}"
"       Linker: .............: ${CMAKE_LINKER}\n"
"\n"
"       Build type ..........: ${CMAKE_BUILD_TYPE}\n"
"       Build shared ........: ${BUILD_SHARED_LIBS}\n"
"       CFlags ..............: ${CMAKE_C_FLAGS}\n"
"       LDFlags .............: ${CMAKE_C_LINK_FLAGS}\n"
"       EXE LDFlags .........: ${CMAKE_EXE_LINKER_FLAGS}\n"
"\n"
"       Floating points precisions enabled\n"
"       CHAMELEON_PRECISIONS : ${CHAMELEON_PRECISION}\n"
"\n"
"       Implementation paradigm\n"
"       CUDA ................: ${CHAMELEON_USE_CUDA}\n"
"       HIP .................: ${CHAMELEON_USE_HIP}\n"
"       MPI .................: ${CHAMELEON_USE_MPI}\n"
"\n"
"       Runtime specific\n"
"       OPENMP ..............: ${CHAMELEON_SCHED_OPENMP}\n"
"       PARSEC ..............: ${CHAMELEON_SCHED_PARSEC}\n"
"       QUARK ...............: ${CHAMELEON_SCHED_QUARK}\n"
"       STARPU ..............: ${CHAMELEON_SCHED_STARPU}\n"
"\n"
"       Kernels specific\n"
"       BLAS ................: ${BLAS_VENDOR_FOUND}\n"
"       LAPACK...............: ${LAPACK_VENDOR_FOUND}\n"
"       HMAT-OSS.............: ${CHAMELEON_USE_HMATOSS}\n"
"\n"
"       Simulation mode .....: ${CHAMELEON_SIMULATION}\n"
"       EZTrace Module.......: ${CHAMELEON_USE_EZTRACE}\n"
"\n"
"       Binaries to build\n"
"       documentation ........: ${CHAMELEON_ENABLE_DOC}\n"
"       example ..............: ${CHAMELEON_ENABLE_EXAMPLE}\n"
"       testing ..............: ${CHAMELEON_ENABLE_TESTING}\n"
"\n")

set(dep_message "${dep_message}"
"\n"
"       INSTALL_PREFIX ......: ${CMAKE_INSTALL_PREFIX}\n\n")

string(REPLACE ";" " " dep_message_wsc "${dep_message}")
message(${dep_message})
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/config.log "${dep_message_wsc}")
message(STATUS "Configuration is done - A summary of the current configuration"
"\n   has been written in ${CMAKE_CURRENT_BINARY_DIR}/config.log")
# installation
# ------------
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.log DESTINATION share/chameleon)
