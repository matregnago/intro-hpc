###
#
#  @copyright 2017-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 0.1.0
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 2020-12-14
#
###
cmake_minimum_required(VERSION 3.20)

# Note that headers are optional, and do not affect add_library, but they will not
# show up in IDEs unless they are listed in add_library.
set(HEADER_LIST
  ${HQR_SOURCE_DIR}/include/libhqr.h
  ${HQR_SOURCE_DIR}/include/libhqr_internal.h
  ${HQR_SOURCE_DIR}/include/libhqr_queue.h
)

set(SOURCE_LIST
  # Low level tree functions
  low_flat.c
  low_binary.c
  low_fibonacci.c
  low_greedy.c
  low_greedy1p.c
  low_adaptiv.c
  # High level tree functions
  high_flat.c
  high_binary.c
  high_fibonacci.c
  high_greedy.c
  # Direct access tree functions
  systolic.c
  svd.c
  hqr.c
  tphqr.c
  tshqr.c
  mtxtree.c
  # Others
  check.c
  gendot.c
  gensvg.c
  print.c
  queue.c
)

# Make an automatic library - will be static or dynamic based on user setting
add_library(hqr ${SOURCE_LIST} ${HEADER_LIST})

set_target_properties(hqr PROPERTIES VERSION ${HQR_VERSION})

# We need this directory, to build.
target_include_directories(hqr PUBLIC
  $<BUILD_INTERFACE:${HQR_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
set_property(TARGET hqr PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

# Users of our library will need the header defining the API
set_target_properties(hqr PROPERTIES PUBLIC_HEADER ${HQR_SOURCE_DIR}/include/libhqr.h)

# Depends on m at link
target_link_libraries(hqr PRIVATE m)

# Add definitions
if ( HAVE_FALLTHROUGH )
  target_compile_definitions(hqr PRIVATE LIBHQR_HAVE_FALLTHROUGH)
endif()

# Target export
install(EXPORT hqrTargets
        NAMESPACE HQR::
        DESTINATION lib/cmake/hqr)

# Install executables, headers, lib, targets file
install(TARGETS hqr
        EXPORT hqrTargets
        RUNTIME DESTINATION bin
        PUBLIC_HEADER DESTINATION include
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

# see https://cmake.org/cmake/help/latest/module/CMakePackageConfigHelpers.html
set(BIN_INSTALL_DIR "bin/"     CACHE PATH "Location of binary files" )
set(INC_INSTALL_DIR "include/" CACHE PATH "Location of header files" )
set(LIB_INSTALL_DIR "lib/"     CACHE PATH "Location of library files" )

include(CMakePackageConfigHelpers)
configure_package_config_file(HQRConfig.cmake.in
                              ${CMAKE_CURRENT_BINARY_DIR}/HQRConfig.cmake
                              INSTALL_DESTINATION ${LIB_INSTALL_DIR}/cmake/hqr
                              PATH_VARS BIN_INSTALL_DIR INC_INSTALL_DIR LIB_INSTALL_DIR)
write_basic_package_version_file(HQRConfigVersion.cmake
                                 VERSION ${HQR_VERSION}
                                 COMPATIBILITY AnyNewerVersion)

# Install config files
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/HQRConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/HQRConfigVersion.cmake
        DESTINATION lib/cmake/hqr)

# Config files (pkg_config, hqr_env.sh)
generate_hqr_pkgconfig_file()
generate_env_file()
