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
#  @date 2024-07-15
#
###

# Add CUDA kernel if compiler and toolkit are available
# -----------------------------------------------------
include(CheckLanguage)
check_language(CUDA)

if(CMAKE_CUDA_COMPILER)
  enable_language(CUDA)
  find_package(CUDAToolkit)
else()
  message(STATUS "CUDA language is not supported")
endif()

if (CUDAToolkit_FOUND)
  set(GPUCUBLAS_HAVE_CUDA_TOOLKIT ON CACHE INTERNAL "Indicate if cuda kernels are enabled or not" FORCE)

  include(SetCMakeCudaArchitectures)

  # Check compatibility between the CUDA version and targets
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 11.1)
    # Convert the space-separated string into a proper CMake list
    separate_arguments(CHAMELEON_CUDA_TARGET_LIST UNIX_COMMAND "${CHAMELEON_CUDA_TARGETS}")
    if("Fermi" IN_LIST CHAMELEON_CUDA_TARGET_LIST OR "Kepler" IN_LIST CHAMELEON_CUDA_TARGET_LIST)
      message(STATUS "Old architecture support requested (Fermi or Kepler found)")
      message(FATAL_ERROR "CMAKE_CUDA_COMPILER_VERSION is ${CMAKE_CUDA_COMPILER_VERSION}."
        " CUDA version > 11.0 has removed support for Fermi/Kepler,"
        " adapt CHAMELEON_CUDA_TARGETS list to remove those.")
    endif()
  endif()

  if (CUDAToolkit_VERSION VERSION_GREATER_EQUAL "7.5")
    set(GPUCUBLAS_HAVE_CUDA_HALF ON CACHE INTERNAL "Indicate if half precision support is enabled or not" FORCE)
  else()
    set(GPUCUBLAS_HAVE_CUDA_HALF OFF CACHE INTERNAL "Indicate if half precision support is enabled or not" FORCE)
  endif()
else()
  set(GPUCUBLAS_HAVE_CUDA_TOOLKIT OFF CACHE INTERNAL "Indicate if cuda kernels are enabled or not" FORCE)
  set(GPUCUBLAS_HAVE_CUDA_HALF OFF CACHE INTERNAL "Indicate if half precision support is enabled or not" FORCE)
endif()

if ( GPUCUBLAS_HAVE_CUDA_HALF )
  ##morse_cmake_required_set( CUBLAS )
  set(CMAKE_REQUIRED_LIBRARIES CUDA::cublas)

  check_function_exists(cublasHgemm GPUCUBLAS_HAVE_CUBLASHGEMM)
  if ( GPUCUBLAS_HAVE_CUBLASHGEMM )
    message("-- ${Blue}Add definition HAVE_CUBLASHGEMM${ColourReset}")
  endif()
  check_function_exists(cublasGemmEx GPUCUBLAS_HAVE_CUBLASGEMMEX)
  if ( GPUCUBLAS_HAVE_CUBLASGEMMEX )
    message("-- ${Blue}Add definition HAVE_CUBLASGEMMEX${ColourReset}")
  endif()
  morse_cmake_required_unset()
endif()

add_subdirectory(include)
add_subdirectory(compute)
add_subdirectory(eztrace_module)

###
### END CMakeLists.txt
###
