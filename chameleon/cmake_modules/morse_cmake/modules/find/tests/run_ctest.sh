#!/bin/bash
set -ex
cmake -B build -S ./modules/find/tests \
               -DENABLE_CTEST=ON \
               -DLAPACKE_COMPONENTS="TMG" \
               -DQUARK_COMPONENTS="HWLOC" \
               -DSTARPU_COMPONENTS="MPI"
ctest --test-dir build --no-compress-output --verbose --output-junit report.xml
