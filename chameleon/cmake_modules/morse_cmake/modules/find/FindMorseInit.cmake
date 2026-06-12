###
#
# @copyright (c) 2019 Inria. All rights reserved.
#
###
#
#  @file FindMorseInit.cmake
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
#  @date 24-04-2018
#
###


# This include is required to check symbols of libs
include(CheckFunctionExists)

# This include is required to check defines in headers
include(CheckIncludeFiles)

# Factorize some piece of code
include(ParseArguments)

# Factorize some piece of code
include(FindMorseCommon)

# To find headers and libs
include(FindHeadersAndLibs)

# To transform relative path into absolute for a list of libraries
include(LibrariesAbsolutePath)

# Some macros to print status when search for headers and libs
include(PrintFindStatus)

# To use pkg_search_module macro
set(FPHSA_NAME_MISMATCHED 1) # Suppress warnings, see https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
include(FindPkgConfig)
unset(FPHSA_NAME_MISMATCHED)

##
## @end file FindMorseInit.cmake
##
