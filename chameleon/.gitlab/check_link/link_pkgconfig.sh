#!/usr/bin/env bash
###
#
#  @file link_pkgconfig.sh
#  @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.4.0
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @date 2024-06-17
#
# Check that linking with the project is ok when using pkg-config.
#
###
set -ex

static=""
chameleon_path=""

while [ $# -gt 0 ]
do
    case $1 in
        --static | -static | -s )
            static="--static"
            ;;
        * )
            chameleon_path=$1
            ;;
    esac
    shift
done

if [ -z $chameleon_path ]
then
    echo """
usage: ./link_pkgconfig.sh path_to_chameleon_install [--static|-static|-s]
   use the --static parameter if chameleon is static (.a)
   env. var. CC and FC must be defined to C and Fortran90 compilers
"""
    exit 1
fi

export PKG_CONFIG_PATH=$chameleon_path/lib/pkgconfig:$PKG_CONFIG_PATH

mkdir -p build
cd build

FLAGS=`pkg-config $static --cflags chameleon`
if [[ "$SYSTEM" == "macosx" ]]; then
    FLAGS="-Wl,-rpath,$chameleon_path/lib $FLAGS"
fi
LIBS=`pkg-config $static --libs chameleon`
$CC $FLAGS ../../../example/link_chameleon/link_chameleon.c $LIBS -o link_chameleon_c

if [[ "$VERSION" != "starpu_simgrid" ]]; then
    ./link_chameleon_c
fi

cd ..
rm -rf build
