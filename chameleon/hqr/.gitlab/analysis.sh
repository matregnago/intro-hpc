#!/usr/bin/env bash
###
#
#  @file analysis.sh
#  @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.0.0
#  @author Mathieu Faverge
#  @date 2021-03-25
#
###

# Performs an analysis of HQR source code:
# - we consider to be in HQR's source code root
# - we consider having the junit.xml file and the coverage file hqr-coverage.xml in the root directory
# - we consider having cppcheck, rats, sonar-scanner programs available in the environment
set -e
set -x
if [ $# -gt 0 ]
then
    BUILDDIR=$1
fi
BUILDDIR=${BUILDDIR:-build}

# List source files:
rm -f filelist.txt
git ls-files | grep "\.[ch]$" > filelist.txt

# Undefine this because not relevant in our configuration
export UNDEFINITIONS="-UWIN32 -UWIN64 -U_MSC_EXTENSIONS -U_MSC_VER -U__SUNPRO_C -U__SUNPRO_CC -U__sun -Usun -U__cplusplus"

# run cppcheck analysis
cppcheck -v -f --project=build/compile_commands.json --language=c --platform=unix64 --enable=all --xml --xml-version=2 --suppress=checkersReport --suppress=missingIncludeSystem ${UNDEFINITIONS} ${DEFINITIONS} 2> hqr-cppcheck.xml

# run rats analysis
rats -w 3 --xml  `cat filelist.txt` > hqr-rats.xml

# create the sonarqube config file
cat > sonar-project.properties << EOF
sonar.host.url=https://sonarqube.inria.fr/sonarqube

sonar.projectKey=solverstack_hqr_AZJSz1v4sbMNg1jXgm3B
sonar.qualitygate.wait=true

sonar.links.homepage=$CI_PROJECT_URL
sonar.links.scm=$CI_REPOSITORY_URL
sonar.links.ci=$CI_PROJECT_URL/pipelines
sonar.links.issue=$CI_PROJECT_URL/issues

sonar.projectDescription=Library for hierarchical QR/LQ reduction trees
sonar.projectVersion=0.1

sonar.scm.disabled=false
sonar.scm.provider=git
sonar.scm.exclusions.disabled=true

sonar.sourceEncoding=UTF-8
sonar.sources=include, src, testings
sonar.inclusions=`cat filelist.txt | xargs echo | sed 's/ /, /g'`

sonar.cxx.jsonCompilationDatabase=${BUILDDIR}/compile_commands.json
sonar.cxx.file.suffixes=.h,.c
sonar.cxx.errorRecoveryEnabled=true
sonar.cxx.gcc.encoding=UTF-8
sonar.cxx.gcc.regex=(?<file>.*):(?<line>[0-9]+):[0-9]+:\\\x20warning:\\\x20(?<message>.*)\\\x20\\\[(?<id>.*)\\\]
sonar.cxx.gcc.reportPaths=hqr-build*.log
sonar.cxx.xunit.reportPaths=junit.xml
sonar.cxx.cobertura.reportPaths=coverage.xml
sonar.cxx.cppcheck.reportPaths=hqr-cppcheck.xml
sonar.cxx.rats.reportPaths=hqr-rats.xml
EOF

# run sonar analysis + publish on sonarqube
sonar-scanner -X > sonar.log
