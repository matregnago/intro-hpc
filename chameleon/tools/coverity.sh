#!/usr/bin/env bash
#
# @file coverity.sh
#
# @copyright 2024-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
# @version 1.4.0
# @author Florent Pruvost
# @date 2024-12-03
#
set -e
set -x
CHAMELEON_SRC_DIR=${CHAMELEON_SRC_DIR:-$PWD}

# install plotting dependencies in the testing docker image
CURRENTUSER=`whoami`
if [[ -f /.dockerenv && $CURRENTUSER == "gitlab" ]]; then
  export PKG_CONFIG_PATH=/home/gitlab/install/starpu/lib/pkgconfig:$PKG_CONFIG_PATH
fi

cmake -B build-$VERSION -C cmake_modules/gitlab-ci-initial-cache.cmake -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_USE_MPI=ON
cov-build --dir cov-int/ cmake --build build-$VERSION -j4
tar czvf chameleon.tgz cov-int/
curl --form token=$COVERITY_TOKEN \
     --form email=florent.pruvost@inria.fr \
     --form file=@chameleon.tgz \
     --form version="`git rev-parse --short HEAD`" \
     --form description="" \
     https://scan.coverity.com/builds?project=Chameleon
