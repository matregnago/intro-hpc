###
#
# @copyright (c) 2019 Inria. All rights reserved.
#
###
#
#  @file FindMorseCommon.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 1.0.0
#  @author Florent Pruvost
#  @date 13-04-2018
#
###
cmake_minimum_required(VERSION 3.12)

# Macro to cleanup packages variables and avoid duplications
# ----------------------------------------------------------
macro(morse_cleanup_variables _prefix)
  #  <PREFIX>_LIBRARIES      ... only the libraries (w/o the '-l')
  #  <PREFIX>_LIBRARY_DIRS   ... the paths of the libraries (w/o the '-L')
  #  <PREFIX>_LDFLAGS_OTHER  ... all other linker flags
  #  <PREFIX>_INCLUDE_DIRS   ... the '-I' preprocessor flags (w/o the '-I')
  #  <PREFIX>_CFLAGS_OTHER   ... the other compiler flags
  if (${_prefix}_LIBRARIES)
    list(REVERSE ${_prefix}_LIBRARIES)
    list(REMOVE_DUPLICATES ${_prefix}_LIBRARIES)
    list(REVERSE ${_prefix}_LIBRARIES)
  endif()
  if (${_prefix}_LIBRARY_DIRS)
    list(REMOVE_DUPLICATES ${_prefix}_LIBRARY_DIRS)
  endif()
  if (${_prefix}_LDFLAGS_OTHER)
    list(REMOVE_DUPLICATES ${_prefix}_LDFLAGS_OTHER)
  endif()
  if (${_prefix}_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES ${_prefix}_INCLUDE_DIRS)
  endif()
  if (${_prefix}_CFLAGS_OTHER)
    list(REMOVE_DUPLICATES ${_prefix}_CFLAGS_OTHER)
  endif()
endmacro()

# clean these variables before using them in CMAKE_REQUIRED_* variables in
# check_function_exists
macro(morse_finds_remove_duplicates)
  if (REQUIRED_DEFINITIONS)
    list(REMOVE_DUPLICATES REQUIRED_DEFINITIONS)
  endif()
  if (REQUIRED_INCDIRS)
    list(REMOVE_DUPLICATES REQUIRED_INCDIRS)
  endif()
  if (REQUIRED_FLAGS)
    list(REMOVE_DUPLICATES REQUIRED_FLAGS)
  endif()
  if (REQUIRED_LDFLAGS)
    list(REMOVE_DUPLICATES REQUIRED_LDFLAGS)
  endif()
  if (REQUIRED_LIBS)
    list(REVERSE REQUIRED_LIBS)
    list(REMOVE_DUPLICATES REQUIRED_LIBS)
    list(REVERSE REQUIRED_LIBS)
  endif()
endmacro()

# try to guess libraries extension, if dyanmic or static
# if dynamic we set ${package}_STATIC to 0, else to 1
# if we do not find the extension ${package}_STATIC is not set at all
macro(morse_check_static_or_dynamic package libraries)
  list(GET ${libraries} 0 _first_lib)
  get_filename_component(_suffix ${_first_lib} EXT)

  message( DEBUG "[mcstod] package ${package}")
  message( DEBUG "[mcstod] libraries ${libraries} ${${libraries}}")
  message( DEBUG "[mcstod] _suffix ${_suffix} ${_first_lib}")

  if (NOT _suffix)
    unset (_lib_path CACHE)
    find_library(_lib_path ${_first_lib} HINTS ${${package}_LIBDIR} ${${package}_LIBRARY_DIRS} NO_DEFAULT_PATH)

    message( DEBUG "[mcstod] Could not find suffix, try to find the library again" )
    message( DEBUG "[mcstod] _first_lib ${_first_lib}"   )
    message( DEBUG "[mcstod] ${${package}_LIBRARY_DIRS}" )
    message( DEBUG "[mcstod] _lib_path ${_lib_path}"     )

    get_filename_component(_suffix ${_lib_path} EXT)

    message( DEBUG "[mcstod] _suffix ${_suffix}")
  endif()
  if (_suffix)
    # some libraries provide the version number so that the suffix becomes
    # something like .3.so for example. Thus only keep .so
    if (_suffix MATCHES "\\.so$" OR _suffix MATCHES "\\.so\\.")
      set(_suffix ".so")
    endif()
    if (_suffix MATCHES "\\.dylib$" OR _suffix MATCHES "\\.dylib\\.")
      set(_suffix ".dylib")
    endif()
    set(${package}_STATIC 0)
    if (WIN32)
      if(${_suffix} MATCHES "\\.lib$")
        set(${package}_STATIC 1)
      endif()
    endif ()
    if (APPLE)
      if(${_suffix} MATCHES "\\.lib$")
        set(${package}_STATIC 1)
      endif()
    endif ()
    if(${_suffix} MATCHES "\\.a$")
      set(${package}_STATIC 1)
    endif()

    # Check that the extension is known
    if(NOT ${package}_STATIC)
      if ( NOT ${_suffix} IN_LIST CMAKE_FIND_LIBRARY_SUFFIXES )
        message( WARNING "${package} library has an unknown extension (${_suffix})")
      endif()
    endif()
  else()
    message( WARNING "${package} could not detect library extension")
  endif()
endmacro()

# add imported target for non-cmake projects or projects which do not provide
# "PROJECT"Config.cmake file at installation
macro(morse_create_imported_target name)

  if(NOT TARGET MORSE::${name})

    # initialize imported target
    add_library(MORSE::${name} INTERFACE IMPORTED)

    if (TARGET PkgConfig::${name})
      get_target_property(_INCLUDES  PkgConfig::${name} INTERFACE_INCLUDE_DIRECTORIES)
      get_target_property(_LIBDIRS   PkgConfig::${name} INTERFACE_LINK_DIRECTORIES)
      get_target_property(_LIBRARIES PkgConfig::${name} INTERFACE_LINK_LIBRARIES)
      get_target_property(_CFLAGS    PkgConfig::${name} INTERFACE_COMPILE_OPTIONS)
      get_target_property(_LDFLAGS   PkgConfig::${name} INTERFACE_LINK_OPTIONS)

      string(REGEX REPLACE "\"" "" _CFLAGS "${_CFLAGS}")
      set_target_properties(MORSE::${name} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${_INCLUDES}")
      set_target_properties(MORSE::${name} PROPERTIES INTERFACE_LINK_DIRECTORIES    "${_LIBDIRS}")
      set_target_properties(MORSE::${name} PROPERTIES INTERFACE_LINK_LIBRARIES      "${_LIBRARIES}")
      set_target_properties(MORSE::${name} PROPERTIES INTERFACE_COMPILE_OPTIONS     "${_CFLAGS}")
      set_target_properties(MORSE::${name} PROPERTIES INTERFACE_LINK_OPTIONS        "${_LDFLAGS}")
    else ()
      if (${name}_INCLUDE_DIRS)
        set_target_properties(MORSE::${name} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${${name}_INCLUDE_DIRS}")
      endif()
      if (${name}_LIBRARY_DIRS)
        set_target_properties(MORSE::${name} PROPERTIES INTERFACE_LINK_DIRECTORIES "${${name}_LIBRARY_DIRS}")
      elseif (${name}_LIBDIR)
        set_target_properties(MORSE::${name} PROPERTIES INTERFACE_LINK_DIRECTORIES "${${name}_LIBDIR}")
        set (${name}_LIBRARY_DIRS ${${name}_LIBDIR})
      endif()
      if (${name}_LIBRARIES)
        set_target_properties(MORSE::${name} PROPERTIES INTERFACE_LINK_LIBRARIES "${${name}_LIBRARIES}")
      endif()
      if (${name}_CFLAGS_OTHER)
        string(REGEX REPLACE "\"" "" ${name}_CFLAGS_OTHER "${${name}_CFLAGS_OTHER}")
        set_target_properties(MORSE::${name} PROPERTIES INTERFACE_COMPILE_OPTIONS "${${name}_CFLAGS_OTHER}")
      endif()
      if (${name}_LDFLAGS_OTHER)
        set_target_properties(MORSE::${name} PROPERTIES INTERFACE_LINK_OPTIONS "${${name}_LDFLAGS_OTHER}")
      endif()
    endif()

  endif (NOT TARGET MORSE::${name})

  set(debug_morse_create_imported_target "FALSE")
  if (debug_morse_create_imported_target)
    if (TARGET MORSE::${name})
      get_target_property(_INCLUDES MORSE::${name} INTERFACE_INCLUDE_DIRECTORIES)
      get_target_property(_DIRECTORIES MORSE::${name} INTERFACE_LINK_DIRECTORIES)
      get_target_property(_LIBRARIES MORSE::${name} INTERFACE_LINK_LIBRARIES)
      get_target_property(_CFLAGS MORSE::${name} INTERFACE_COMPILE_OPTIONS)
      get_target_property(_LDFLAGS MORSE::${name} INTERFACE_LINK_OPTIONS)
      message(STATUS "IMPORTED TARGET ${name}:
                      _INCLUDES ${_INCLUDES}
                      _DIRECTORIES ${_DIRECTORIES}
                      _LIBRARIES ${_LIBRARIES}
                      _CFLAGS ${_CFLAGS}
                      _LDFLAGS ${_LDFLAGS}")
    endif()
  endif()

endmacro()

# export existing imported target ${targetname} in cmake file
# ${targetfilename}Targets.cmake and installed in lib/cmake/${cmakesubdir}/. To be used
# by "PROJECT"Config.cmake file at installation
macro(morse_export_imported_target targetname targetfilename cmakesubdir)

  if (TARGET ${targetname})
    get_target_property(_INCLUDES  ${targetname} INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(_LIBDIRS   ${targetname} INTERFACE_LINK_DIRECTORIES)
    get_target_property(_LIBRARIES ${targetname} INTERFACE_LINK_LIBRARIES)
    get_target_property(_CFLAGS    ${targetname} INTERFACE_COMPILE_OPTIONS)
    get_target_property(_DEFS      ${targetname} INTERFACE_COMPILE_DEFINITIONS)
    get_target_property(_LDFLAGS   ${targetname} INTERFACE_LINK_OPTIONS)

    if (NOT _INCLUDES)
      set(_INCLUDES "")
    endif()
    if (NOT _LIBDIRS)
      set(_LIBDIRS "")
    endif()
    if (NOT _LIBRARIES)
      set(_LIBRARIES "")
    endif()
    if (NOT _CFLAGS)
      set(_CFLAGS "")
    endif()
    if (NOT _DEFS)
      set(_DEFS "")
    endif()
    if (NOT _LDFLAGS)
      set(_LDFLAGS "")
    endif()

    get_target_property(_IMPORTED_CONF ${targetname} IMPORTED_CONFIGURATIONS)
    if (NOT _IMPORTED_CONF)
      set(_IMPORTED_CONF "NOCONFIG")
    endif()
    get_target_property(_IMPORTED_DEP ${targetname} IMPORTED_LINK_DEPENDENT_LIBRARIES_${_IMPORTED_CONF})
    get_target_property(_IMPORTED_LOCATION ${targetname} IMPORTED_LOCATION_${_IMPORTED_CONF})
    get_target_property(_IMPORTED_SONAME ${targetname} IMPORTED_SONAME_${_IMPORTED_CONF})
    if (_IMPORTED_DEP)
      list(PREPEND _LIBRARIES ${_IMPORTED_DEP})
    endif()
    if (_IMPORTED_LOCATION)
      list(PREPEND _LIBRARIES ${_IMPORTED_LOCATION})
    endif()

    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${targetfilename}Targets.cmake"
"if (NOT TARGET ${targetname})
add_library(${targetname} INTERFACE IMPORTED)
set_target_properties(${targetname} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES \"${_INCLUDES}\")
set_target_properties(${targetname} PROPERTIES INTERFACE_LINK_DIRECTORIES    \"${_LIBDIRS}\")
set_target_properties(${targetname} PROPERTIES INTERFACE_LINK_LIBRARIES      \"${_LIBRARIES}\")
set_target_properties(${targetname} PROPERTIES INTERFACE_COMPILE_OPTIONS     \"${_CFLAGS}\")
set_target_properties(${targetname} PROPERTIES INTERFACE_COMPILE_DEFINITIONS \"${_DEFS}\")
set_target_properties(${targetname} PROPERTIES INTERFACE_LINK_OPTIONS        \"${_LDFLAGS}\")
endif()
")

    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${targetfilename}Targets.cmake"
            DESTINATION "lib/cmake/${cmakesubdir}")
  else()
    set(list_var "${ARGV}")
    if(NOT "QUIET" IN_LIST list_var)
      message(WARNING "morse_export_imported_target_without_namespace: try to export target "
        "${targetname} but it does not exist.")
    endif()
  endif()

endmacro()

# Set the CMAKE_REQUIRED_... porperties to check libraries
# --------------------------------------------------------
macro(morse_set_required_test_lib_link name)
  set(CMAKE_REQUIRED_INCLUDES "${${name}${STATIC}_INCLUDE_DIRS}")
  if (${name}${STATIC}_CFLAGS_OTHER)
    set(REQUIRED_FLAGS_COPY "${${name}${STATIC}_CFLAGS_OTHER}")
    set(REQUIRED_FLAGS)
    set(REQUIRED_DEFINITIONS)
    foreach(_flag ${REQUIRED_FLAGS_COPY})
      if (_flag MATCHES "^-D")
       list(APPEND REQUIRED_DEFINITIONS "${_flag}")
      endif()
      string(REGEX REPLACE "^-D.*" "" _flag "${_flag}")
      list(APPEND REQUIRED_FLAGS "${_flag}")
    endforeach()
  endif()
  morse_finds_remove_duplicates()
  set(CMAKE_REQUIRED_DEFINITIONS "${REQUIRED_DEFINITIONS}")
  set(CMAKE_REQUIRED_FLAGS "${REQUIRED_LDFLAGS}")
  set(CMAKE_REQUIRED_LIBRARIES)
  list(APPEND CMAKE_REQUIRED_LIBRARIES "${${name}${STATIC}_LDFLAGS_OTHER}")
  if (${name}${STATIC}_LIBRARY_DIRS)
    foreach(_dir ${${name}${STATIC}_LIBRARY_DIRS})
      list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${_dir}")
    endforeach()
  endif()
  list(APPEND CMAKE_REQUIRED_LIBRARIES "${${name}${STATIC}_LIBRARIES}")
  string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")
endmacro()

# Set the CMAKE_REQUIRED_... porperties to check libraries
# --------------------------------------------------------
macro(morse_cmake_required_set prefix)
  if (${prefix}_CFLAGS_OTHER)
    set(REQUIRED_FLAGS_COPY "${${prefix}_CFLAGS_OTHER}")
    set(REQUIRED_FLAGS)
    set(REQUIRED_DEFINITIONS)
    foreach(_flag ${REQUIRED_FLAGS_COPY})
      if (_flag MATCHES "^-D")
       list(APPEND REQUIRED_DEFINITIONS "${_flag}")
      endif()
      string(REGEX REPLACE "^-D.*" "" _flag "${_flag}")
      list(APPEND REQUIRED_FLAGS "${_flag}")
    endforeach()
  endif()
  set(CMAKE_REQUIRED_DEFINITIONS "${REQUIRED_DEFINITIONS}")
  set(CMAKE_REQUIRED_FLAGS       "${REQUIRED_FLAGS}")
  set(CMAKE_REQUIRED_INCLUDES    "${${prefix}_INCLUDE_DIRS}")
  set(CMAKE_REQUIRED_LIBRARIES   "${${prefix}_LDFLAGS_OTHER}")
  foreach(_dir ${${prefix}_LIBRARY_DIRS})
    list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${_dir}")
  endforeach()
  list(APPEND CMAKE_REQUIRED_LIBRARIES "${${prefix}_LIBRARIES}")
  string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")
endmacro()

# Unset the CMAKE_REQUIRED_... properties
# ---------------------------------------
macro(morse_cmake_required_unset)
  set(CMAKE_REQUIRED_DEFINITIONS)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_LIBRARIES)
endmacro()

# Transform relative path into absolute path for libraries found with the
# pkg_search_module cmake macro
# _prefix: the name of the CMake variable used when pkg_search_module was called
# e.g. for pkg_search_module(BLAS blas) _prefix would be BLAS
macro(morse_find_pkgconfig_libraries_absolute_path _prefix)

  set(${_prefix}_LIBRARIES_COPY "${${_prefix}_LINK_LIBRARIES}")
  set(${_prefix}_LIBRARIES "")
  foreach(_library ${${_prefix}_LIBRARIES_COPY})
    # The full path is given, let's store it and move to the next one
    if(EXISTS "${_library}")
      list(APPEND ${_prefix}_LIBRARIES ${_library})
      continue()
    endif()

    set (CMAKE_FIND_LIBRARY_SUFFIXES_COPY ${CMAKE_FIND_LIBRARY_SUFFIXES})
    get_filename_component(_ext "${_library}" EXT)
    list(FIND CMAKE_FIND_LIBRARY_SUFFIXES "${_ext}" _index)

    # Define the extension to look for
    if (${_index} GREATER -1)
      get_filename_component(_library "${_library}" NAME_WE)
      set(CMAKE_FIND_LIBRARY_SUFFIXES ${_ext} ${CMAKE_FIND_LIBRARY_SUFFIXES})
    else()
      if (${_prefix}_STATIC)
        if (WIN32)
          set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
        endif ()
        if (APPLE)
          set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
        else ()
          set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
        endif ()
      endif()
    endif()

    find_library( _library_path
      NAMES ${_library}
      HINTS ${${_prefix}_LIBDIR} ${${_prefix}_LIBRARY_DIRS}
      )

    if (_library_path)
      list(APPEND ${_prefix}_LIBRARIES ${_library_path})
    else()
      if (${_prefix}_STATIC)
        message(STATUS "${_prefix}_STATIC ${${_prefix}_STATIC}")
      endif()
      message(WARNING "Dependency of ${_prefix} '${_library}' NOT FOUND with suffixes ${CMAKE_FIND_LIBRARY_SUFFIXES} in ${${_prefix}_LIBDIR} ${${_prefix}_LIBRARY_DIRS}")
    endif()

    unset(_library_path CACHE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_COPY})

  endforeach()
  set (${_prefix}_LIBRARIES "${${_prefix}_LIBRARIES}" CACHE INTERNAL "" FORCE)

endmacro()

# install necessary morse modules files (mods), in dest, when distribute a cmake
# lib depending on it
macro(morse_install_finds mods dest )
  # install specific dependencies of the caller, given in mods
  foreach(_mod ${${mods}})
    install(FILES ${MORSE_CMAKE_MODULE_PATH}/find/Find${_mod}.cmake
            DESTINATION ${dest})
  endforeach()

  # install other necessary morse files containing macros
  set(morse_find_core "../ParseArguments.cmake;FindHeadersAndLibs.cmake;FindMorseCommon.cmake;FindMorseInit.cmake;LibrariesAbsolutePath.cmake;PrintFindStatus.cmake;MORSE-Copyright.txt")
  foreach(_file ${morse_find_core})
     install(FILES ${MORSE_CMAKE_MODULE_PATH}/find/${_file}
           DESTINATION ${dest})
  endforeach()
endmacro()

macro(morse_find_package_get_envdir _package)
  # This macro checks for the existence of environment variables
  # that would help to find packages.  For each variable
  # ${_package}_DIR, ${_package}_INCDIR, and ${_package}_LIBDIR is the
  # environment variable exists and the cmake cache one is not
  # defined, then it defines the cmake variable to the environment
  # one.

  if(NOT ${_package}_DIR AND DEFINED ENV{${_package}_DIR})
    set(${_package}_DIR "$ENV{${_package}_DIR}"
      CACHE PATH "Installation directory of the ${_package} package given by the user")
    mark_as_advanced(${_package}_DIR)
  endif()

  if(NOT ${_package}_INCDIR AND DEFINED ENV{${_package}_INCDIR})
    set(${_package}_INCDIR "$ENV{${_package}_INCDIR}"
      CACHE PATH "Header installation directory the ${_package} package given by the user")
    mark_as_advanced(${_package}_INCDIR)
  endif()

  if(NOT ${_package}_LIBDIR AND DEFINED ENV{${_package}_LIBDIR})
    set(${_package}_LIBDIR "$ENV{${_package}_LIBDIR}"
      CACHE PATH "Libraries installation directory the ${_package} package given by the user")
    mark_as_advanced(${_package}_LIBDIR)
  endif()

  if(NOT ${_package}_BINDIR AND DEFINED ENV{${_package}_BINDIR})
    set(${_package}_BINDIR "$ENV{${_package}_BINDIR}"
      CACHE PATH "Binaries installation directory the ${_package} package given by the user")
    mark_as_advanced(${_package}_BINDIR)
  endif()
endmacro()

macro(morse_find_path _package )
  # This macro is an extension of the find_path function to search
  # explictly for files in the user defined directory first before
  # looking for it in system directories.
  parse_arguments(mfp "HEADERS;SUFFIXES;HINTS" "OPTIONAL" ${ARGN})

  message( DEBUG "[MFP/${_package}] HEADERS  : ${mfp_HEADERS}"  )
  message( DEBUG "[MFP/${_package}] SUFFIXES : ${mfp_SUFFIXES}" )
  message( DEBUG "[MFP/${_package}] HINTS    : ${mfp_HINTS}"    )
  message( DEBUG "[MFP/${_package}] OPTIONAL : ${mfp_OPTIONAL}" )

  # Adjust verbosity
  # ----------------
  if(${_package}_FIND_QUIETLY)
    set(_outlvl VERBOSE)
  else()
    set(_outlvl STATUS)
  endif()

  # Prepare the package variables to fill
  # -------------------------------------
  set(${_package}_INCLUDE_DIRS)

  # Call cmake macro to find all the headers path
  # ---------------------------------------------
  foreach(_hdrfile ${mfp_HEADERS})
    set(${_package}_${_hdrfile}_DIRS "${_package}_${_hdrfile}_DIRS-NOTFOUND")

    string(REGEX REPLACE "include/" "" _suffixes "${mfp_SUFFIXES}")
    if(${_package}_INCDIR)
      find_path(${_package}_${_hdrfile}_DIRS
        NAMES ${_hdrfile}
        HINTS ${${_package}_INCDIR}
        PATH_SUFFIXES ${_suffixes}
        NO_PACKAGE_ROOT_PATH NO_CMAKE_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_FIND_ROOT_PATH)
    else()
      if(${_package}_DIR)
        find_path(${_package}_${_hdrfile}_DIRS
          NAMES ${_hdrfile}
          HINTS ${${_package}_DIR}
          PATH_SUFFIXES ${mfp_SUFFIXES}
          NO_PACKAGE_ROOT_PATH NO_CMAKE_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_FIND_ROOT_PATH)
      else()
        find_path(${_package}_${_hdrfile}_DIRS
          NAMES ${_hdrfile}
          HINTS ${mfp_HINTS}
          PATH_SUFFIXES ${_suffixes})
      endif()
    endif()
    mark_as_advanced(${_package}_${_hdrfile}_DIRS)

    if( ${_package}_${_hdrfile}_DIRS )
      list(APPEND ${_package}_INCLUDE_DIRS ${${_package}_${_hdrfile}_DIRS} )
      message( VERBOSE    "[MFP/${_package}] Found ${_hdrfile} in ${${_package}_${_hdrfile}_DIRS}" )
    else()
      if ( NOT ${mfp_OPTIONAL} )
        set( ${_package}_INCLUDE_DIRS "${_package}_INCLUDE_DIRS-NOTFOUND" )
        break()
      endif()
      message( ${_outlvl} "[MFP/${_package}] ${_hdrfile} not found" )
    endif()
  endforeach()

  list(REMOVE_DUPLICATES ${_package}_INCLUDE_DIRS)
endmacro()

macro(morse_find_library _package )
  # This macro is an extension of the find_library function to search
  # explicitly for files in the user defined directory first before
  # looking for it in system directories.
  parse_arguments(mfl "LIBRARIES;SUFFIXES;HINTS" "OPTIONAL" ${ARGN})

  message( DEBUG "[MFL/${_package}] LIBRARIES : ${mfl_LIBRARIES}" )
  message( DEBUG "[MFL/${_package}] SUFFIXES  : ${mfl_SUFFIXES}"  )
  message( DEBUG "[MFL/${_package}] HINTS     : ${mfl_HINTS}"     )
  message( DEBUG "[MFL/${_package}] OPTIONAL  : ${mfl_OPTIONAL}"  )

  # Adjust verbosity
  # ----------------
  if(${_package}_FIND_QUIETLY)
    set(_outlvl VERBOSE)
  else()
    set(_outlvl STATUS)
  endif()

  # Set the suffixes for static libraries
  # -------------------------------------
  set (CMAKE_FIND_LIBRARY_SUFFIXES_COPY ${CMAKE_FIND_LIBRARY_SUFFIXES})
  if (${_package}_STATIC)
    if (WIN32)
      set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
    endif ( WIN32 )
    if (APPLE)
      set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
    endif()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else()
    # for ubuntu's libblas3gf and libscalapack3gf packages
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} .so.3gf)
  endif()

  # Prepare the package variables to fill
  # -------------------------------------
  set(${_package}_LIBRARIES)
  set(${_package}_LIBRARY_DIRS)

  # Call cmake macro to find each library
  # -------------------------------------
  foreach(_library ${mfl_LIBRARIES})
    set(${_package}_${_library}_LIBRARY "${_package}_${_library}_LIBRARY-NOTFOUND")

    if(${_package}_LIBDIR)
      find_library(${_package}_${_library}_LIBRARY
        NAMES ${_library}
        HINTS ${${_package}_LIBDIR}
        NO_PACKAGE_ROOT_PATH NO_CMAKE_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_FIND_ROOT_PATH)
    else()
      if(${_package}_DIR)
        find_library(${_package}_${_library}_LIBRARY
          NAMES ${_library}
          HINTS ${${_package}_DIR}
          PATH_SUFFIXES ${mfl_SUFFIXES}
          NO_PACKAGE_ROOT_PATH NO_CMAKE_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_FIND_ROOT_PATH)
      else()
        find_library(${_package}_${_library}_LIBRARY
          NAMES ${_library}
          HINTS ${mfl_HINTS})
      endif()
    endif()
    mark_as_advanced(${_package}_${_library}_LIBRARY)

    if( ${_package}_${_library}_LIBRARY )
      # Register the found library to the list
      # --------------------------------------
      list(APPEND ${_package}_LIBRARIES ${${_package}_${_library}_LIBRARY})

      get_filename_component(_path ${${_package}_${_library}_LIBRARY} PATH)
      list(APPEND ${_package}_LIBRARY_DIRS ${_path})

      message( VERBOSE    "[MFL/${_package}] Found ${_library} - ${${_package}_${_library}_LIBRARY}" )
    else()
      message( ${_outlvl} "[MFL/${_package}] ${_library} not found" )
      if ( NOT ${mfl_OPTIONAL} )
        set( ${_package}_LIBRARIES    "${_package}_LIBRARIES-NOTFOUND"    )
        set( ${_package}_LIBRARY_DIRS "${_package}_LIBRARY_DIRS-NOTFOUND" )
        break()
      endif()
    endif()
  endforeach()
  list(REMOVE_DUPLICATES ${_package}_LIBRARY_DIRS)

  mark_as_advanced(${_package}_LIBRARIES)
  mark_as_advanced(${_package}_LIBRARY_DIRS)

  # Restore suffixes
  # ----------------
  set (CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_COPY})
endmacro()

##
## @end file FindMorseCommon
##
