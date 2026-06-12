#!/usr/bin/env bash
###
#
#  @file chameleon_guix.sh
#  @copyright 2018-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.4.0
#  @author Florent Pruvost
#  @date 2025-12-19
#
###
set -ex

# Configure and Build Chameleon
if [ -d build-$NODE-$MPI ]; then
  rm build-$NODE-$MPI -r
fi
cmake -B build-$NODE-$MPI $CHAMELEON_BUILD_OPTIONS
cmake --build build-$NODE-$MPI -j20  > /dev/null
export CHAMELEON_BUILD=$PWD/build-$NODE-$MPI
export HFI_NO_CPUAFFINITY=1
if [[ "${NODE}" == "sirocco" ]]; then
    export LD_PRELOAD="/usr/lib64/libcuda.so"
fi


# clean old benchmarks
if [ -d tools/bench/$PLATFORM/results ]; then
  rm tools/bench/$PLATFORM/results -r
fi
# Execute jube benchmarks
jube run tools/bench/$PLATFORM/chameleon.xml --tag gemm potrf geqrf --include-path tools/bench/$PLATFORM/parameters/$NODE --id $JUBE_ID
#jube run tools/bench/$PLATFORM/chameleon-test.xml --tag gemm potrf geqrf --include-path tools/bench/$PLATFORM/parameters/$NODE --id $JUBE_ID
# jube analysis
jube analyse tools/bench/$PLATFORM/results --id $JUBE_ID
# jube report
jube result tools/bench/$PLATFORM/results --id $JUBE_ID > chameleon\-$NODE\-$MPI.csv
cat chameleon\-$NODE\-$MPI.csv

# send results to the elasticsearch server
#ls guix.json
# python3 tools/bench/jube/add_result.py -e https://elasticsearch.bordeaux.inria.fr -t hiepacs -p "chameleon" -m $MPI chameleon\-$NODE\-$MPI.csv
#python3 tools/bench/jube/add_result.py -e https://elasticsearch.bordeaux.inria.fr -t hiepacs -p "chameleon-test" -m $MPI chameleon\-$NODE\-$MPI.csv
