###
#
# @copyright (c) 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
#  @file GenHQRPkgConfig.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 0.1.0
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 10-11-2014
#
###

###
#
# GENERATE_PKGCONFIG_FILE: generate files hqr.pc
#
###
macro(GENERATE_HQR_PKGCONFIG_FILE)

  set(HQR_PKGCONFIG_LIBS "-lhqr")
  set(HQR_PKGCONFIG_LIBS_PRIVATE "-lm")
  set(HQR_PKGCONFIG_REQUIRED "")
  set(HQR_PKGCONFIG_REQUIRED_PRIVATE "")

  #clean_lib_list(HQR)

  set(_output_hqr_file "${CMAKE_BINARY_DIR}/hqr.pc")
  configure_file(
    "${HQR_SOURCE_DIR}/lib/pkgconfig/hqr.pc.in"
    "${_output_hqr_file}"
    @ONLY
    )
  install(
    FILES ${_output_hqr_file}
    DESTINATION lib/pkgconfig
    )

endmacro(GENERATE_HQR_PKGCONFIG_FILE)

###
#
# generate_env_file: generate files pastix.pc
#
###
macro(generate_env_file)

    # Create .sh file
    # ---------------
    configure_file(
      "${HQR_SOURCE_DIR}/hqr_env.sh.in"
      "${CMAKE_BINARY_DIR}/bin/hqr_env.sh" @ONLY)

    # installation
    # ------------
    install(FILES "${CMAKE_BINARY_DIR}/bin/hqr_env.sh"
      DESTINATION bin)

endmacro(generate_env_file)

##
## @end file GenPkgConfig.cmake
##
