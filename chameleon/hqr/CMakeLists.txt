###
#
#  @copyright 2017-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 0.1.0
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 2017-04-26
#
###
cmake_minimum_required(VERSION 3.12)
project (HQR C)

# Check if compiled independently or within another project
if ( ${HQR_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set( BUILD_SUBPROJECT OFF )
else()
  set( BUILD_SUBPROJECT ON )
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")
include(GenHQRPkgConfig)
include (CheckCCompilerFlag)

# The current version number
# --------------------------
set (HQR_VERSION_MAJOR 0)
set (HQR_VERSION_MINOR 1)
set (HQR_VERSION_MICRO 0)

set( HQR_VERSION "${HQR_VERSION_MAJOR}.${HQR_VERSION_MINOR}.${HQR_VERSION_MICRO}" )

# Set warnings for debug builds
# -----------------------------
check_c_compiler_flag( "-Wall" HAVE_WALL )
if( HAVE_WALL )
    set( C_WFLAGS "${C_WFLAGS} -Wall" )
endif( HAVE_WALL )
check_c_compiler_flag( "-Wextra" HAVE_WEXTRA )
if( HAVE_WEXTRA )
    set( C_WFLAGS "${C_WFLAGS} -Wextra" )
endif( HAVE_WEXTRA )

# Check for attribute fallthrough
# -------------------------------
check_c_source_compiles("
#include <stdarg.h>
int main(void) {
  int a = 2;
  switch( a ){
  case 0:
     __attribute__((fallthrough));
  default:
      ;
  }
  return 0;
}
"
  HAVE_FALLTHROUGH
  FAIL_REGEX ".*fallthrough.*" )

### Misc options
if ( NOT BUILD_SUBPROJECT )
  option(BUILD_SHARED_LIBS
    "Build shared libraries" OFF)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are None, Debug, Release, RelWithDebInfo and MinSizeRel." FORCE)
  endif(NOT CMAKE_BUILD_TYPE)
endif()

# The compiled library code is here
add_subdirectory(src)

###############################################################################
# Executables and Tests #
#########################
enable_testing()
include(CTest)

add_subdirectory(testings)

#-- Add a custom target to generate tags
add_custom_target (tags
  COMMAND "git ls-files | xargs etags" )
