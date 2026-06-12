###
#
# @file ChameleonPkgConfig.cmake
#
# @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
###
#
#  @project Chameleon
#  @version 1.4.0
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @author Abel Calluaud
#  @date 2024-12-16
#
###
include(GenPkgConfig)

###
#
# GENERATE_PKGCONFIG_FILE: generate files chameleon_lapack.pc, chameleon.pc, coreblas.pc and gpucublas.pc or gpuhipblas.pc
#
###
macro(chameleon_generate_pkgconfig_files)

  # The link flags specific to this package and any required libraries
  # that don't support PkgConfig
  set(CHAMELEON_LAPACK_PKGCONFIG_LIBS chameleon_lapack)
  set(CHAMELEON_PKGCONFIG_LIBS        chameleon       )
  set(COREBLAS_PKGCONFIG_LIBS         coreblas        )
  set(GPUCUBLAS_PKGCONFIG_LIBS        gpucublas       )
  set(GPUHIPBLAS_PKGCONFIG_LIBS       gpuhipblas      )

  # The link flags for private libraries required by this package but not
  # exposed to applications
  set(CHAMELEON_LAPACK_PKGCONFIG_LIBS_PRIVATE "")
  set(CHAMELEON_PKGCONFIG_LIBS_PRIVATE        "")
  set(COREBLAS_PKGCONFIG_LIBS_PRIVATE         "")
  set(GPUCUBLAS_PKGCONFIG_LIBS_PRIVATE        "")
  set(GPUHIPBLAS_PKGCONFIG_LIBS_PRIVATE       "")

  # A list of packages required by this package
  set(CHAMELEON_LAPACK_PKGCONFIG_REQUIRED "chameleon")
  set(CHAMELEON_PKGCONFIG_REQUIRED        "hqr")
  set(COREBLAS_PKGCONFIG_REQUIRED         "")
  set(GPUCUBLAS_PKGCONFIG_REQUIRED        "")
  set(GPUHIPBLAS_PKGCONFIG_REQUIRED       "")

  # A list of private packages required by this package but not exposed to
  # applications
  set(CHAMELEON_LAPACK_PKGCONFIG_REQUIRED_PRIVATE "")
  set(CHAMELEON_PKGCONFIG_REQUIRED_PRIVATE        "")
  set(COREBLAS_PKGCONFIG_REQUIRED_PRIVATE         "")
  set(GPUCUBLAS_PKGCONFIG_REQUIRED_PRIVATE        "")
  set(GPUHIPBLAS_PKGCONFIG_REQUIRED_PRIVATE       "")

  if(CHAMELEON_SCHED_OPENMP)
    list(APPEND CHAMELEON_PKGCONFIG_LIBS -lchameleon_openmp)
    list(APPEND CHAMELEON_PKGCONFIG_LIBS_PRIVATE "${OpenMP_C_LIBRARIES}")
  elseif(CHAMELEON_SCHED_PARSEC)
    list(APPEND CHAMELEON_PKGCONFIG_LIBS -lchameleon_parsec)
    list(APPEND CHAMELEON_PKGCONFIG_LIBS_PRIVATE "${PARSEC_LIBRARIES}")
  elseif(CHAMELEON_SCHED_QUARK)
    list(APPEND CHAMELEON_PKGCONFIG_LIBS -lchameleon_quark)
    list(APPEND CHAMELEON_PKGCONFIG_LIBS_PRIVATE "${QUARK_LIBRARIES_DEP}")
  elseif(CHAMELEON_SCHED_STARPU)
    list(APPEND CHAMELEON_PKGCONFIG_LIBS -lchameleon_starpu)
    if ( CHAMELEON_SIMULATION )
      # need PUBLIC here because simulation mode with starpu requires to replace the
      # main function by a starpu main one
      if ( CHAMELEON_USE_MPI )
        list(APPEND CHAMELEON_PKGCONFIG_REQUIRED starpumpi-${CHAMELEON_STARPU_VERSION})
      else()
        list(APPEND CHAMELEON_PKGCONFIG_REQUIRED starpu-${CHAMELEON_STARPU_VERSION})
      endif()
    else()
      if ( CHAMELEON_USE_MPI )
        list(APPEND CHAMELEON_PKGCONFIG_REQUIRED_PRIVATE starpumpi-${CHAMELEON_STARPU_VERSION})
      else()
        list(APPEND CHAMELEON_PKGCONFIG_REQUIRED_PRIVATE starpu-${CHAMELEON_STARPU_VERSION})
      endif()
    endif()
  endif()

  if(NOT CHAMELEON_SIMULATION)

    list(APPEND COREBLAS_PKGCONFIG_LIBS_PRIVATE
      ${LAPACKE_LIBRARIES}
      ${CBLAS_LIBRARIES}
    )
    list(APPEND CHAMELEON_PKGCONFIG_REQUIRED "coreblas")

    if(CHAMELEON_USE_CUDA)
      list(APPEND GPUCUBLAS_PKGCONFIG_LIBS_PRIVATE ${CUDA_CUSOLVER_LIBRARY} ${CUDA_CUBLAS_LIBRARIES})
      list(APPEND CHAMELEON_PKGCONFIG_REQUIRED "gpucublas")
    endif()

    if(CHAMELEON_USE_HIP)
      list(APPEND GPUHIPBLAS_PKGCONFIG_LIBS_PRIVATE ${HIPBLAS_LIBRARIES})
      list(APPEND GPUHIPBLAS_PKGCONFIG_LIBS_PRIVATE ${HIP_LIBRARIES})
      list(APPEND CHAMELEON_PKGCONFIG_REQUIRED "gpuhipblas")
    endif()

  endif(NOT CHAMELEON_SIMULATION)

  list(APPEND CHAMELEON_PKGCONFIG_LIBS_PRIVATE
    ${M_LIBRARIES}
  )

  if(CHAMELEON_USE_MPI)
    if("${MPI_C_LIBRARIES}" MATCHES "mpich")
      set(MPI_NAME "mpich")
    elseif("${MPI_C_LIBRARIES}" MATCHES "mvapich2")
      set(MPI_NAME "mvapich2")
    elseif("${MPI_C_LIBRARIES}" MATCHES "madmpi" OR "${MPI_C_LIBRARIES}" MATCHES "nmad")
      set(MPI_NAME "nmad")
    elseif("${MPI_C_LIBRARIES}" MATCHES "openmpi")
      set(MPI_NAME "ompi")
    else()
      set(MPI_NAME "mpi")
    endif()
    find_file(MPIPC_PATH "${MPI_NAME}.pc"
              HINTS ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}
              ENV PKG_CONFIG_PATH
              PATH_SUFFIXES pkgconfig
              NO_CACHE)
    if(MPIPC_PATH)
      list(APPEND CHAMELEON_PKGCONFIG_REQUIRED ${MPI_NAME})
    else()
      list(APPEND CHAMELEON_PKGCONFIG_INCS "${MPI_C_INCLUDE_DIRS}")
      list(APPEND CHAMELEON_PKGCONFIG_LIBS "${MPI_C_LIBRARIES}")
    endif()
  endif()

  generate_pkgconfig_files(
    ${CMAKE_SOURCE_DIR}/lib/pkgconfig/chameleon.pc.in
    PROJECTNAME CHAMELEON )

  generate_pkgconfig_files(
    ${CMAKE_SOURCE_DIR}/lib/pkgconfig/chameleon_lapack.pc.in
    PROJECTNAME CHAMELEON_LAPACK )

  generate_pkgconfig_files(
    ${CMAKE_SOURCE_DIR}/lib/pkgconfig/coreblas.pc.in
    PROJECTNAME COREBLAS )

  if(NOT CHAMELEON_SIMULATION)
    if(CHAMELEON_USE_CUDA)
      generate_pkgconfig_files(
        ${CMAKE_SOURCE_DIR}/lib/pkgconfig/gpucublas.pc.in
        PROJECTNAME GPUCUBLAS )
    endif()

    if(CHAMELEON_USE_HIP)
      generate_pkgconfig_files(
        ${CMAKE_SOURCE_DIR}/lib/pkgconfig/gpuhipblas.pc.in
        PROJECTNAME GPUHIPBLAS )
    endif()

  endif(NOT CHAMELEON_SIMULATION)
endmacro()

##
## @end file ChameleonPkgConfig.cmake
##
