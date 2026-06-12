#!/usr/bin/env bash
set -e
set -x

if [[ "$SYSTEM" == "linux" ]]; then

  # on linux we perform a coverage analysis
  cmake -B build -S . -DCMAKE_INSTALL_PREFIX=$PWD/install \
        -DBUILD_SHARED_LIBS=ON -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DCMAKE_C_FLAGS="-O0 -g -fPIC --coverage -Werror -Wall -fdiagnostics-show-option -fno-inline" \
        -DCMAKE_EXE_LINKER_FLAGS="--coverage" -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

else

  # no coverage analysis on other platforms
  cmake -B build -S . -DCMAKE_INSTALL_PREFIX=$PWD/install \
        -DBUILD_SHARED_LIBS=ON -DCMAKE_VERBOSE_MAKEFILE=ON

fi
cmake --build build -j 4
cmake --install build
