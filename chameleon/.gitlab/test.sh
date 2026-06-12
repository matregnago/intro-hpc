#!/usr/bin/env bash
set -e
set -x

if [[ "$VERSION" == "starpu_cuda" || "$VERSION" == "starpu_hip" ]]; then
  sed -i "s#realpath \(.*\)/install-starpu_cuda#realpath ${PWD}/install-starpu_cuda#g" ./install-${VERSION}/bin/chameleon_env
  sed -i "s#realpath \(.*\)/install-starpu_hip#realpath ${PWD}/install-starpu_hip#g" ./install-${VERSION}/bin/chameleon_env
  source ./install-${VERSION}/bin/chameleon_env
else
  source ./.gitlab-ci-env.sh $CHAM_CI_ENV_ARG
fi
cd build-$VERSION
ctest --output-on-failure --no-compress-output $TESTS_RESTRICTION --output-junit ../${LOGNAME}-junit.xml
if [[ "$SYSTEM" == "linux" ]]; then
  # clang is used on macosx and it is not compatible with MORSE_ENABLE_COVERAGE=ON
  # so that we can only make the coverage report on the linux runner with gcc
  cd ..
  lcov --directory build-$VERSION --capture --output-file ${LOGNAME}.lcov
fi
