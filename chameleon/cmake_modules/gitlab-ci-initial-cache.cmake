set(BUILD_SHARED_LIBS "ON" CACHE BOOL "")
set(CMAKE_INSTALL_PREFIX "$ENV{PWD}/install-$ENV{VERSION}" CACHE PATH "")
set(CMAKE_VERBOSE_MAKEFILE "ON" CACHE BOOL "")

if(CMAKE_C_FLAGS)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror -fprofile-update=atomic" CACHE STRING "")
else()
  set(CMAKE_C_FLAGS "-Werror -fprofile-update=atomic" CACHE STRING "")
endif()

option(MORSE_ENABLE_WARNING  "Enable warning messages"        ON)
option(MORSE_ENABLE_COVERAGE "Enable flags for coverage test" ON)

set(CHAMELEON_TESTINGS_VENDOR ON CACHE BOOL "")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "")

set(BLA_VENDOR "Generic" CACHE STRING "Set to one of the BLAS/LAPACK Vendors to search for BLAS only from the specified vendor.")
