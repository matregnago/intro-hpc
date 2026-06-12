#!/usr/bin/env bash
###
#
#  @file link_cmake.sh
#  @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.4.0
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @date 2025-12-19
#
# Check that linking with the project is ok when using cmake.
#
###
set -ex

if [ $# -lt 1 ]
then
    echo "usage: ./link_cmake.sh path_to_install"
    exit 1
fi

cmake -S ../../example/link_chameleon -B build -DCMAKE_PREFIX_PATH=$1
cmake --build build --verbose

if [[ "$VERSION" != "starpu_simgrid" ]]; then
    ctest --test-dir build --verbose
fi

rm -rf build
