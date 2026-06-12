###
#
#  @copyright 2017-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 0.1.0
#  @author Mathieu Faverge
#  @date 2017-04-26
#
###

if ( BUILD_SUBPROJECT )
  set( HQR_INSTALL_DIR bin/hqr )
else()
  set( HQR_INSTALL_DIR bin )
endif()

set(TESTINGS
  draw_hqr.c
  draw_systolic.c
  draw_tphqr.c
  testing_hqr.c
  testing_svd.c
  testing_systolic.c
  testing_tphqr.c
  testing_print.c
  )

foreach (_file ${TESTINGS})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we}
    ${_file} common.c)
  target_link_libraries(${_name_we} PRIVATE hqr)

  install(
    TARGETS ${_name_we}
    RUNTIME DESTINATION ${HQR_INSTALL_DIR} )
  install(
    FILES   ${_file}
    DESTINATION ${HQR_INSTALL_DIR} )

endforeach()

if ( NOT BUILD_SUBPROJECT )
  add_test( draw_hqr      ./draw_hqr      -M 27 -N 7 -P 2 -d -a 3 -p 2 )
  add_test( draw_systolic ./draw_systolic -M 27 -N 7 -P 3 -L 3 -l 2 )
  add_test( draw_tphqr    ./draw_tphqr    -M 27 -P 2)
  add_test( test_hqr      ./testing_hqr      -X )
  add_test( test_svd      ./testing_svd      -X )
  add_test( test_systolic ./testing_systolic -X )
  add_test( test_tphqr    ./testing_tphqr    -X )
  add_test( single_hqr      ./testing_hqr      -v 2 )
  add_test( single_svd      ./testing_svd      -v 2 )
  add_test( single_systolic ./testing_systolic -v 2 )
  add_test( single_tphqr    ./testing_tphqr    -v 2 )
  add_test( usage         ./testing_tphqr -h )
  add_test( test_print    ./testing_print -M 27 -N 7 -P 3 )
endif()
