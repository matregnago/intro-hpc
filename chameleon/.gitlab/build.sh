#!/usr/bin/env bash
set -ex

SCAN=""

# Configure with CMake
case $SYSTEM in

  guix)
    echo "build on guix"
    cmake -B build-${VERSION} -S . -C cmake_modules/gitlab-ci-initial-cache.cmake $BUILD_OPTIONS
    ;;

  linux)
    echo "build on linux"
    source .gitlab-ci-env.sh $CHAM_CI_ENV_ARG
    if [[ $CI_COMMIT_REF_NAME == $CI_DEFAULT_BRANCH ]]; then
      SCAN="scan-build -plist --intercept-first --exclude CMakeFiles --analyze-headers -o analyzer_reports "
    fi
    eval '${SCAN}cmake -B build-${VERSION} -S . -C cmake_modules/gitlab-ci-initial-cache.cmake $BUILD_OPTIONS'
    ;;

  macosx)
    echo "build on macosx"
    # check starpu is already installed and install it if necessary
    DEP_INSTALLED=`brew ls --versions starpu | cut -d " " -f 2`
    if [[ -z "${DEP_INSTALLED}" ]]; then
      # dep is not installed, we have to install it
      brew install starpu
    fi
    # clang is used on macosx and it is not compatible with MORSE_ENABLE_COVERAGE=ON
    # to avoid the Accelerate framework and get Openblas we use BLA_PREFER_PKGCONFIG
    cmake -B build-${VERSION} -S . \
          -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$PWD/install-${VERSION} \
          -DMORSE_ENABLE_COVERAGE=OFF -DBLA_PREFER_PKGCONFIG=ON
    ;;

  windows)
    echo "build on windows"
    # on windows the mpi_f08 interface is missing, see https://www.scivision.dev/windows-mpi-msys2/
    # do not use static libraries because executables are too large and the build
    # directory can reach more than 10Go
    cmake -GNinja -B build-${VERSION} -S . \
          -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$PWD/install-${VERSION} \
          -DCHAMELEON_USE_MPI=OFF
    ;;
  *)
    echo "The SYSTEM environment variable is $SYSTEM. It is not one of : guix, linux, macosx, windows -> exit 1."
    exit 1
    ;;
esac

# Compile
eval '${SCAN}cmake --build build-${VERSION} -j ${CMAKE_BUILD_PARALLEL_LEVEL} > /dev/null'

# Install
cmake --install build-${VERSION}


#
# Check link to chameleon
#

# Set the path variables
if [[ "$SYSTEM" == "windows" ]]; then
  export PATH="/c/Windows/WinSxS/x86_microsoft-windows-m..namespace-downlevel_31bf3856ad364e35_10.0.19041.1_none_21374cb0681a6320":$PATH
fi
source install-${VERSION}/bin/chameleon_env

# Set the compiler
if [[ "$SYSTEM" == "macosx" ]]; then
  export CC=clang
  if brew ls --versions chameleon > /dev/null; then brew remove --force --ignore-dependencies chameleon; fi
else
  export CC=gcc
fi
export FC=gfortran

cd .gitlab/check_link/
# 1) using cmake:
./link_cmake.sh $PWD/../../install-${VERSION}
# 2) using pkg-config:
./link_pkgconfig.sh $PWD/../../install-${VERSION}

cd ../..
# Remove the install directory to save space, except for starpu_cuda and starpu_hip because on
# plafrim with gpu (shell runner, not docker) the install directory is used for the tests
if [[ "$VERSION" != "starpu_cuda" && "$VERSION" != "starpu_hip" ]]; then
  rm -r install-${VERSION}
fi
