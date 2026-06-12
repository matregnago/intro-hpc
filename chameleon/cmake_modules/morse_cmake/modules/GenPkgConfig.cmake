###
#
# @copyright (c) 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                          Univ. Bordeaux. All rights reserved.
#
###
#
#  @file GenPkgConfig.cmake
#
#  @project MORSE
#
#  @version 0.9.0
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @date 2018-06-21
#
#.rst:
# GenPkgConfig
# ------------
#
#   The GenPkgConfig module intend to provide functions to generate a
#   pkg-config file with the right format, and a environment file that
#   the user can source to setup its environment with the library
#   installed.
#
# generate_pkgconfig_files(
#   file1 file2 ...
#   [PROJECTNAME name]
#   [LIBS lib1 lib2 ...]
#   [LIBS_PRIVATE lib1 lib2 ...]
#   [REQUIRED pkg1 pkg2 ...]
#   [REQUIRED_PRIVATE pkg1 pkg2 ...]
# )
# Where:
#   - file1, file2, .. are the different input files for the
#     pkg-config configuration files.
#   - PROJECTNAME defines the prefix of the variables to use for the
#     default following values ${PROJECTNAME}_PKGCONFIG_XXX with XXX
#     being part of (LIBS, LIBS_PRIVATE, REQUIRED and
#     REQUIRED_PRIVATE)
#   - LIBS defines the libs of the current project for dynamic
#     linking.
#   - LIBS_PRIVATE defines the libs of the current project and its
#     required dependencies not handled by pkg-config for static
#     linking.
#   - LIBS defines the additional pkg-config packages needed by the
#     project for dynamic linking.
#   - LIBS_PRIVATE defines the additional pkg-config packages needed
#     by the project for static linking.
#
# The generated files are then installed in the subdirectory
# lib/pkgconfig of the prefix directory.
#
# generate_env_file(
#   [PROJECTNAME name]
# )
# Where:
#   - PROJECTNAME defines the project name of the library if different
#     from CMAKE_PROJECT_NAME.
#
# The generated file is installed in the subdirectory bin of the
# prefix directory.
#
# generate_newenv_file(
#   [PACKAGE name]
#   [PACKAGE_TEST test]
# )
# Where:
#   - PACKAGE defines the project name of the library if different
#     from CMAKE_PROJECT_NAME.
#   - PACKAGE_TEST defines a test string to check if the installation
#   - is correct within the generated file.
#
# While the first version generates a file named $PACKAGE_env.sh, the
# new one generates a $PACKAGE_env file to avoid conflict
#
# The generated file is installed in the subdirectory bin of the
# prefix directory.
#
###

###
#
# gpc_convert_incstyle_to_pkgconfig(): convert a libraries list to
# follow the pkg-config style used in gpc_clean_lib_list
#
###
macro(gpc_convert_incstyle_to_pkgconfig _inclist)
  set(${_inclist}_CPY "${${_inclist}}")
  set(${_inclist} "")
  foreach(_dep ${${_inclist}_CPY})
    if (${_dep} MATCHES "^-")
      # any compiler flag (e.g. -pthread, -DUSE_FOO): keep as it is
      list(APPEND ${_inclist} ${_dep})
    elseif (EXISTS ${_dep})
      # an absolute path to includes: add the -I
      list(APPEND ${_inclist} "-I${_dep}")
    endif()
  endforeach()
endmacro(gpc_convert_incstyle_to_pkgconfig)

###
#
# gpc_convert_libstyle_to_pkgconfig(): convert a libraries list to
# follow the pkg-config style used in gpc_clean_lib_list
#
###
macro(gpc_convert_libstyle_to_pkgconfig _liblist)

  # Start by removing duplicates
  list(REVERSE ${_liblist})
  list(REMOVE_DUPLICATES ${_liblist})
  list(REVERSE ${_liblist})

  # Convert to pkg-config file format
  set(${_liblist}_CPY "${${_liblist}}")
  set(${_liblist} "")
  # we consider a list of libraries as input
  foreach(_dep ${${_liblist}_CPY})
    if (${_dep} MATCHES "^/|^[A-Z]:/")
      # case library is an absolute path (e.g. /opt/foo/lib/libfoo.so)
      get_filename_component(dep_libname ${_dep} NAME)
      get_filename_component(dep_libdir  ${_dep} PATH)
      foreach( _ext ${CMAKE_FIND_LIBRARY_SUFFIXES} )
        string(REGEX REPLACE "${_ext}$" "" dep_libname "${dep_libname}")
      endforeach()
      string(REGEX REPLACE "^lib" "" dep_libname "${dep_libname}")
      list(APPEND ${_liblist} -L${dep_libdir} -l${dep_libname})
    elseif(${_dep} MATCHES "^-l|^-L")
      # case library given with -l (e.g. -lfoo) or -L (e.g. -L/opt/foo/lib -lfoo)
      list(APPEND ${_liblist} ${_dep})
    else()
      # case library given by name (e.g. foo)
      list(APPEND ${_liblist} "-l${_dep}")
    endif()
  endforeach()
endmacro(gpc_convert_libstyle_to_pkgconfig)

###
#
# gpc_clean_lib_list(): clean libraries lists to follow the pkg-config
# style used in GENERATE_PKGCONFIG_FILE
#
###
macro(gpc_clean_lib_list _package)
  if ( ${_package}_PKGCONFIG_CFLAGS )
    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_CFLAGS)
    gpc_convert_incstyle_to_pkgconfig(${_package}_PKGCONFIG_CFLAGS)
    string(REPLACE ";" " " ${_package}_PKGCONFIG_CFLAGS "${${_package}_PKGCONFIG_CFLAGS}")
  endif()
  if ( ${_package}_PKGCONFIG_INCS )
    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_INCS)
    gpc_convert_incstyle_to_pkgconfig(${_package}_PKGCONFIG_INCS)
    string(REPLACE ";" " " ${_package}_PKGCONFIG_INCS "${${_package}_PKGCONFIG_INCS}")
  endif()
  if ( ${_package}_PKGCONFIG_LIBS )
    gpc_convert_libstyle_to_pkgconfig(${_package}_PKGCONFIG_LIBS)
    string(REPLACE ";" " " ${_package}_PKGCONFIG_LIBS "${${_package}_PKGCONFIG_LIBS}")
  endif()
  if ( ${_package}_PKGCONFIG_LIBS_PRIVATE )
    gpc_convert_libstyle_to_pkgconfig(${_package}_PKGCONFIG_LIBS_PRIVATE)
    string(REPLACE ";" " " ${_package}_PKGCONFIG_LIBS_PRIVATE "${${_package}_PKGCONFIG_LIBS_PRIVATE}")
  endif()
  if ( ${_package}_PKGCONFIG_REQUIRED )
    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_REQUIRED)
    string(REPLACE ";" " " ${_package}_PKGCONFIG_REQUIRED "${${_package}_PKGCONFIG_REQUIRED}")
  endif()
  if ( ${_package}_PKGCONFIG_REQUIRED_PRIVATE )
    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_REQUIRED_PRIVATE)
    string(REPLACE ";" " " ${_package}_PKGCONFIG_REQUIRED_PRIVATE "${${_package}_PKGCONFIG_REQUIRED_PRIVATE}")
  endif()
endmacro(gpc_clean_lib_list)

###
#
# generate_pkgconfig_file: generate pkg-config file
#
###
macro(generate_pkgconfig_files)
  set(_options )
  set(_oneValueArgs PROJECTNAME)
  set(_multiValueArgs LIBS LIBS_PRIVATE REQUIRED REQUIRED_PRIVATE)
  set(MN "generate_pkgconfig_files")
  cmake_parse_arguments(${MN}
    "${_options}" "${_oneValueArgs}"
    "${_multiValueArgs}" ${ARGN} )

  if ( NOT DEFINED ${MN}_PROJECTNAME )
    set( ${MN}_PROJECTNAME ${CMAKE_PROJECT_NAME} )
  endif()

  set(ARGN ${${MN}_UNPARSED_ARGUMENTS})

  # The link flags specific to this package and any required libraries
  # that don't support PkgConfig
  list(APPEND ${${MN}_PROJECTNAME}_PKGCONFIG_LIBS ${${MN}_LIBS})

  # The link flags for private libraries required by this package but not
  # exposed to applications
  list(APPEND ${${MN}_PROJECTNAME}_PKGCONFIG_LIBS_PRIVATE ${${MN}_LIBS_PRIVATE})

  # A list of packages required by this package
  list(APPEND ${${MN}_PROJECTNAME}_PKGCONFIG_REQUIRED ${${MN}_REQUIRED})

  # A list of private packages required by this package but not exposed to
  # applications
  list(APPEND ${${MN}_PROJECTNAME}_PKGCONFIG_REQUIRED_PRIVATE ${${MN}_REQUIRED_PRIVATE})

  # Define required package
  # -----------------------
  gpc_clean_lib_list(${${MN}_PROJECTNAME})

  foreach(f IN LISTS ARGN)
    get_filename_component(fname "${f}" NAME_WE)

    # Create .pc files
    # ----------------
    configure_file(
      ${f}
      "${CMAKE_BINARY_DIR}/lib/pkgconfig/${fname}.pc" @ONLY)

    # installation
    # ------------
    install(FILES
      "${CMAKE_BINARY_DIR}/lib/pkgconfig/${fname}.pc"
      DESTINATION lib/pkgconfig )

  endforeach()

endmacro(generate_pkgconfig_files)

###
#
# generate_env_file: generate file to setup environment variables
#
###
macro(generate_env_file)

  set(_options )
  set(_oneValueArgs PROJECTNAME)
  set(_multiValueArgs )
  set(MN "generate_env_file")
  cmake_parse_arguments(generate_env_file
    "${_options}" "${_oneValueArgs}"
    "${_multiValueArgs}" ${ARGN} )

  if ( NOT DEFINED ${MN}_PROJECTNAME )
    set( ${MN}_PROJECTNAME ${CMAKE_PROJECT_NAME} )
  endif()

  string(TOLOWER ${${MN}_PROJECTNAME} LONAME)
  string(TOUPPER ${${MN}_PROJECTNAME} UPNAME)

  # Create .sh file
  # ---------------
  configure_file(
    "${MORSE_CMAKE_MODULE_PATH}/env.sh.in"
    "${CMAKE_CURRENT_BINARY_DIR}/bin/${LONAME}_env.sh" @ONLY)

  # installation
  # ------------
  # Check if the PREFIX is in the PATH or not
  # If yes, no need to install this file
  string(REPLACE ":" ";" pathlist "$ENV{PATH}")
  if ( "${CMAKE_INSTALL_PREFIX}/bin" IN_LIST pathlist )
    message( STATUS "No installation of ${LONAME}_env.sh - already in the default environment" )
  else ()
    message( STATUS "Install ${LONAME}_env.sh in ${CMAKE_INSTALL_PREFIX}/bin" )
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/bin/${LONAME}_env.sh"
      DESTINATION bin)
  endif()

endmacro(generate_env_file)

###
#
# generate_env_file2: generate file to setup environment variables (StarPU-like)
#
###
macro(generate_newenv_file)

  set(_options )
  set(_oneValueArgs PACKAGE PACKAGE_TEST)
  set(_multiValueArgs )
  set(MN "generate_env_file2")
  cmake_parse_arguments(generate_env_file2
    "${_options}" "${_oneValueArgs}"
    "${_multiValueArgs}" ${ARGN} )

  set(MN "generate_pkgconfig_files")
  cmake_parse_arguments(${MN}
    "${_options}" "${_oneValueArgs}"
    "${_multiValueArgs}" ${ARGN} )

  if ( NOT DEFINED ${MN}_PACKAGE )
    set( PACKAGE ${CMAKE_PROJECT_NAME} )
  else()
    set( PACKAGE ${${MN}_PACKAGE} )
  endif()

  if ( NOT DEFINED ${MN}_PACKAGE_TEST )
    set( PACKAGE_TEST "0 = 0" )
  else()
    set( PACKAGE_TEST ${${MN}_PACKAGE_TEST} )
  endif()

  string(TOLOWER ${PACKAGE} PACKAGE_LC)
  string(TOUPPER ${PACKAGE} PACKAGE_UC)

  # Create .sh file
  # ---------------
  configure_file(
    "${MORSE_CMAKE_MODULE_PATH}/env.in"
    "${CMAKE_CURRENT_BINARY_DIR}/bin/${PACKAGE_LC}_env" @ONLY)

  # installation
  # ------------
  # Check if the PREFIX is in the PATH or not
  # If yes, no need to install this file
  string(REPLACE ":" ";" pathlist "$ENV{PATH}")
  if ( "${CMAKE_INSTALL_PREFIX}/bin" IN_LIST pathlist )
    message( STATUS "No installation of ${PACKAGE_LC}_env - already in the default environment" )
  else ()
    message( STATUS "Install ${PACKAGE_LC}_env in ${CMAKE_INSTALL_PREFIX}/bin" )
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/bin/${PACKAGE_LC}_env"
      DESTINATION bin)
  endif()

endmacro(generate_newenv_file)

##
## @end file GenPkgConfig.cmake
##
