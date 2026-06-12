#!/usr/bin/env bash
set -e
set -x

if [[ "$SYSTEM" == "linux" ]]; then
  source install/bin/hqr_env.sh
fi
cd build
if [[ "$SYSTEM" == "windows" ]]; then
  # this is required with BUILD_SHARED_LIBS=ON
  export PATH="/c/Windows/WinSxS/x86_microsoft-windows-m..namespace-downlevel_31bf3856ad364e35_10.0.19041.1_none_21374cb0681a6320":$PATH
  export PATH=$PWD/src:$PATH
fi
ctest --output-on-failure --no-compress-output -T Test --output-junit ../junit.xml
if [[ "$SYSTEM" == "linux" ]]; then
  # clang is used on macosx and it is not compatible with --coverage option
  # so that we can only make the coverage report on the linux runner with gcc
  cd ..
  lcov --capture --directory build -q --output-file hqr.lcov
  lcov --summary hqr.lcov
  lcov_cobertura hqr.lcov --output coverage.xml
fi
