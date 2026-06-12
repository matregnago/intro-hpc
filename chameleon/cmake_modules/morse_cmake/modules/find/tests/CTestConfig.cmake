## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "MorseCmake")
set(CTEST_NIGHTLY_START_TIME "00:00:00 GMT")

set(CTEST_DROP_METHOD "https")
set(CTEST_DROP_SITE "cdash-ci.inria.fr")
set(CTEST_DROP_LOCATION "/submit.php?project=MorseCmake")
set(CTEST_CURL_OPTIONS "CURLOPT_SSL_VERIFYPEER_OFF")

#--------------------------------------------------------------------
# BUILDNAME variable construction
# This variable will be used to set the build name which will appear
# on the Chameleon dashboard http://cdash.inria.fr/CDash/
#--------------------------------------------------------------------
# Start with the short system name, e.g. "Linux", "FreeBSD" or "Windows"
if(NOT BUILDNAME)

  set(BUILDNAME "${CMAKE_SYSTEM_NAME}")

endif()
