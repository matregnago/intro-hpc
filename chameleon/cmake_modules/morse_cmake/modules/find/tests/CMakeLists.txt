cmake_minimum_required(VERSION 3.12)
project(TEST_MORSE_CMAKE_MODULES_FIND Fortran C CXX)

# location of Morse modules to get some specific macros
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../")
include(MorseInit)

# main variable: control the list of libraries to find thanks to find_package
# pay attention that package names must be given with capital letters
set(PACKAGES "" CACHE STRING "List of packages to find, ex: BLAS;STARPU;PASTIX")

# specific components to look for with packages
set(LAPACKE_COMPONENTS "" CACHE STRING "List of specific dependencies to look for with LAPACKE, ex: TMG")
set(QUARK_COMPONENTS "" CACHE STRING "List of specific dependencies to look for with QUARK, ex: HWLOC")
set(PETSC_COMPONENTS "" CACHE STRING "List of specific dependencies to look for with PETSc, ex: CXX|C")
set(FFTW_COMPONENTS "" CACHE STRING "List of specific dependencies to look for with FFTW, ex: MKL;ESSL;THREADS;OMP;SIMPLE;LONG;QUAD")
set(MUMPS_COMPONENTS "" CACHE STRING "List of specific dependencies to look for with MUMPS, ex: MPI;SEQ;SCOTCH;PTSCOTCH;METIS;PARMETIS;OPENMP")
set(STARPU_COMPONENTS "" CACHE STRING "List of specific starpu variants, ex: MPI")

# to enable ctests
option(ENABLE_CTEST "Enable Testing: will test all supported packages" OFF)

foreach(_library ${PACKAGES})

  # in case _library is unset inside find_package
  set (_library_cpy ${_library})

  if (${_library}_COMPONENTS)
    find_package(${_library} COMPONENTS ${${_library}_COMPONENTS})
  else()
    find_package(${_library})
  endif()
  set (_library ${_library_cpy})

  if (${_library}_FOUND)
    if ( ${_library} MATCHES "LAPACKEXT" )
      set( _library LAPACK )
    elseif( ${_library} MATCHES "BLASEXT" )
      set( _library BLAS )
    endif()

    if (${_library}_PREFIX)
      message(STATUS "${_library}_PREFIX ${${_library}_PREFIX}")
    endif()
    if (${_library}_VERSION)
      message(STATUS "${_library}_VERSION ${${_library}_VERSION}")
    endif()
    if (${_library}_STATIC)
      message(STATUS "${_library}_STATIC ${${_library}_STATIC}")
    endif()
    if (${_library}_LIBRARIES)
      message(STATUS "${_library}_LIBRARIES found: ${${_library}_LIBRARIES}")
    else()
      message(WARNING "${_library}_LIBRARIES not found: ${${_library}_LIBRARIES}")
    endif()
    if (${_library}_STATIC_LIBRARIES)
      message(STATUS "${_library}_STATIC_LIBRARIES found: ${${_library}_STATIC_LIBRARIES}")
    endif()
    if (${_library}_LIBRARY_DIRS)
      message(STATUS "${_library}_LIBRARY_DIRS found: ${${_library}_LIBRARY_DIRS}")
    else()
      message(WARNING "${_library}_LIBRARY_DIRS not found: ${${_library}_LIBRARY_DIRS}")
    endif()
    if (${_library}_INCLUDE_DIRS)
      message(STATUS "${_library}_INCLUDE_DIRS found: ${${_library}_INCLUDE_DIRS}")
    else()
      message(WARNING "${_library}_INCLUDE_DIRS not found: ${${_library}_INCLUDE_DIRS}")
    endif()
    if (${_library}_CFLAGS_OTHER)
      message(STATUS "${_library}_CFLAGS_OTHER found: ${${_library}_CFLAGS_OTHER}")
    endif()
    if (${_library}_LDFLAGS_OTHER)
      message(STATUS "${_library}_LDFLAGS_OTHER found: ${${_library}_LDFLAGS_OTHER}")
    endif()
    if (TARGET MORSE::${_library})
      get_target_property(_INCLUDES MORSE::${_library} INTERFACE_INCLUDE_DIRECTORIES)
      get_target_property(_DIRECTORIES MORSE::${_library} INTERFACE_LINK_DIRECTORIES)
      get_target_property(_LIBRARIES MORSE::${_library} INTERFACE_LINK_LIBRARIES)
      get_target_property(_CFLAGS MORSE::${_library} INTERFACE_COMPILE_OPTIONS)
      get_target_property(_LDFLAGS MORSE::${_library} INTERFACE_LINK_OPTIONS)
      message(STATUS "IMPORTED TARGET MORSE::${_library} INTERFACE_INCLUDE_DIRECTORIES ${_INCLUDES}")
      message(STATUS "IMPORTED TARGET MORSE::${_library} INTERFACE_LINK_DIRECTORIES ${_DIRECTORIES}")
      message(STATUS "IMPORTED TARGET MORSE::${_library} INTERFACE_LINK_LIBRARIES ${_LIBRARIES}")
      message(STATUS "IMPORTED TARGET MORSE::${_library} INTERFACE_COMPILE_OPTIONS ${_CFLAGS}")
      message(STATUS "IMPORTED TARGET MORSE::${_library} INTERFACE_LINK_OPTIONS ${_LDFLAGS}")
    endif()
    # for BLAS and LAPACK or other external/official Find that we check
    if (TARGET ${_library}::${_library})
      get_target_property(_INCLUDES ${_library}::${_library} INTERFACE_INCLUDE_DIRECTORIES)
      get_target_property(_DIRECTORIES ${_library}::${_library} INTERFACE_LINK_DIRECTORIES)
      get_target_property(_LIBRARIES ${_library}::${_library} INTERFACE_LINK_LIBRARIES)
      get_target_property(_CFLAGS  ${_library}::${_library} INTERFACE_COMPILE_OPTIONS)
      get_target_property(_LDFLAGS ${_library}::${_library} INTERFACE_LINK_OPTIONS)
      message(STATUS "IMPORTED TARGET ${_library}::${_library} INTERFACE_INCLUDE_DIRECTORIES ${_INCLUDES}")
      message(STATUS "IMPORTED TARGET ${_library}::${_library} INTERFACE_LINK_DIRECTORIES ${_DIRECTORIES}")
      message(STATUS "IMPORTED TARGET ${_library}::${_library} INTERFACE_LINK_LIBRARIES ${_LIBRARIES}")
      message(STATUS "IMPORTED TARGET ${_library}::${_library} INTERFACE_COMPILE_OPTIONS ${_CFLAGS}")
      message(STATUS "IMPORTED TARGET ${_library}::${_library} INTERFACE_LINK_OPTIONS ${_LDFLAGS}")
    endif()
  else()
    message(FATAL_ERROR "${_library} NOT FOUND !!")
  endif()

  set( ORDERINGS SCOTCH PTSCOTCH METIS PAMPA )
  if ( ${_library} IN_LIST ORDERINGS )
    if ( ${_library}_INTSIZE EQUAL -1 )
      message(FATAL_ERROR "${_library} could not find integer size !!")
    endif()
  endif()

endforeach()

# Add CTest rules
if (ENABLE_CTEST)

  enable_testing()
  include(CTest)

  set(CTEST_PACKAGES_LIST
      # AL4SAN # not installed in the common image hpclib/hiepacs, see hpclib/al4san to get it
      BLAS
      BLASEXT
      CBLAS
      CPPCHECK
      # EZTRACE # not installed yet, we should do it
      FFTW
      FXT
      # GTG # parsec build fails with gtg : undefined reference to pajeGetProcFile
      HWLOC
      # HYPRE # not installed cf. https://gitlab.inria.fr/solverstack/docker/-/blob/master/dockerfile-distrib
      LAPACK
      LAPACKEXT
      LAPACKE
      M
      METIS
      MPIEXT
      MUMPS
      # PAMPA # not compatible with openmpi 4
      PAPI
      PARMETIS
      PARSEC
      PETSC
      PTSCOTCH
      QUARK
      SCALAPACK
      SCOTCH
      SIMGRID
      SLEPC
      STARPU
      SUITESPARSE
      RT
      TMG)

  foreach(_package ${CTEST_PACKAGES_LIST})
    add_test(FIND${_package} cmake ${CMAKE_SOURCE_DIR} -DPACKAGES=${_package})
  endforeach()

endif()

###
### END CMakeLists.txt
###
