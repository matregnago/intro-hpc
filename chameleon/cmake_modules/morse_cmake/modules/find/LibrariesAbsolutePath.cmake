###
#
# @copyright (c) 2019 Inria. All rights reserved.
#
###
#
#  @file LibrariesAbsolutePath.cmake
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

# Transform relative path into absolute path for libraries
# lib_list (input/output): the name of the CMake variable containing libraries, e.g. BLAS_LIBRARIES
# hints_paths (input): additional paths to add when looking for libraries
macro(LIBRARIES_ABSOLUTE_PATH lib_list hints_paths)

  # copy the lib list
  set (${lib_list}_COPY "${${lib_list}}")

  # reset the lib list to populate
  set(${lib_list} "")
  foreach(_library ${${lib_list}_COPY})
    if(EXISTS "${_library}")
      # if already an absolute path, nothing special to do
      list(APPEND ${lib_list} ${_library})
    else()
      # replace pattern -lfoo -> foo
      string(REGEX REPLACE "^-l" "" _library "${_library}")

      # remove extensions if exist
      if( ${CMAKE_VERSION} VERSION_GREATER "3.14.0" )
        get_filename_component(_lext "${_library}" LAST_EXT)
        get_filename_component(_ext "${_library}" EXT)
        if ( "${_lext}" IN_LIST CMAKE_FIND_LIBRARY_SUFFIXES )
          get_filename_component(_library "${_library}" NAME_WLE)
        elseif ( "${_ext}" IN_LIST CMAKE_FIND_LIBRARY_SUFFIXES )
          get_filename_component(_library "${_library}" NAME_WE)
        endif()
      else()
        get_filename_component(_ext "${_library}" EXT)
        if ( "${_ext}" IN_LIST CMAKE_FIND_LIBRARY_SUFFIXES )
          get_filename_component(_library "${_library}" NAME_WE)
        endif()
      endif()

      # try to find the lib
      find_library(_library_path
        NAMES ${_library}
        HINTS ${hints_paths}
        )
      if (_library_path)
          list(APPEND ${lib_list} ${_library_path})
      else()
          message(FATAL_ERROR "Dependency of ${lib_list} '${_library}' NOT FOUND")
      endif()
      unset(_library_path CACHE)
    endif()
  endforeach()
endmacro()

##
## @end file LibrariesAbsolutePath.cmake
##
