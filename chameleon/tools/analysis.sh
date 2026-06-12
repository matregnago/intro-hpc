#!/usr/bin/env bash
###
#
#  @file analysis.sh
#  @copyright 2013-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.4.0
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @date 2025-12-19
#
###

# Performs an analysis of Chameleon source code:
# - we consider to be in Chameleon's source code root
# - we consider having the build *.log files in the root directory
# - we consider having the test reports *.junit files in the root directory
# - we consider having the coverage file chameleon_coverage.xml in the root directory
# - we consider having cppcheck, rats, sonar-scanner programs available in the environment

# filter sources:
# - consider generated files in ${BUILDDIR}
# - exclude base *z* files to avoid duplication
# - exclude cblas.h and lapacke-.h because not really part of chameleon and make cppcheck analysis too long
set -e
set -x
if [ $# -gt 0 ]
then
    BUILDDIR=$1
fi
BUILDDIR=${BUILDDIR:-build}

TOOLSDIR=$(dirname $0)
$TOOLSDIR/find_sources.sh

# Undefine this because not relevant in our configuration
export UNDEFINITIONS="-UCHAMELEON_USE_OPENCL -UWIN32 -UWIN64 -U_MSC_EXTENSIONS -U_MSC_VER -U__SUNPRO_C -U__SUNPRO_CC -U__sun -Usun -U__cplusplus"

# run cppcheck analysis
CPPCHECK_OPT=" -v -f --language=c --platform=unix64 --enable=all --xml --xml-version=2  --suppress=checkersReport --suppress=missingIncludeSystem ${UNDEFINITIONS}"
cppcheck $CPPCHECK_OPT --file-list=./filelist_none.txt 2> chameleon_cppcheck.xml
cppcheck $CPPCHECK_OPT -DPRECISION_s -UPRECISION_d -UPRECISION_c -UPRECISION_z --file-list=./filelist_s.txt 2>> chameleon_cppcheck.xml
cppcheck $CPPCHECK_OPT -DPRECISION_d -UPRECISION_s -UPRECISION_c -UPRECISION_z --file-list=./filelist_d.txt 2>> chameleon_cppcheck.xml
cppcheck $CPPCHECK_OPT -DPRECISION_c -UPRECISION_s -UPRECISION_d -UPRECISION_z --file-list=./filelist_c.txt 2>> chameleon_cppcheck.xml
cppcheck $CPPCHECK_OPT -DPRECISION_z -UPRECISION_s -UPRECISION_d -UPRECISION_c --file-list=./filelist_z.txt 2>> chameleon_cppcheck.xml

# merge all different compile_commands.json files from build-* directories into one single file for sonarqube
jq -s 'map(.[])' $PWD/build-*/compile_commands.json > compile_commands.json

# create the sonarqube config file
cat > sonar-project.properties << EOF
sonar.host.url=https://sonarqube.inria.fr/sonarqube

sonar.projectKey=solverstack_chameleon_AZJTCfl1sbMNg1jXgm3k
sonar.qualitygate.wait=true

sonar.links.homepage=$CI_PROJECT_URL
sonar.links.scm=$CI_REPOSITORY_URL
sonar.links.ci=$CI_PROJECT_URL/pipelines
sonar.links.issue=$CI_PROJECT_URL/issues

sonar.projectDescription=Dense linear algebra subroutines for heterogeneous and distributed architectures
sonar.projectVersion=1.4.1

sonar.scm.disabled=false
sonar.scm.provider=git
sonar.scm.exclusions.disabled=true

sonar.sources=build-openmp/runtime/openmp, build-parsec/runtime/parsec, build-quark/runtime/quark, build-starpu, compute, control, coreblas, example, include, runtime, testing
sonar.inclusions=`cat filelist.txt | sed ':a;N;$!ba;s/\n/, /g'`
sonar.coverage.exclusions=build-starpu/gpucublas/**/*,build-starpu/gpuhipblas/**/*
sonar.sourceEncoding=UTF-8
sonar.cxx.file.suffixes=.h,.c
sonar.cxx.errorRecoveryEnabled=true
sonar.cxx.gcc.encoding=UTF-8
sonar.cxx.gcc.regex=(?<file>.*):(?<line>[0-9]+):[0-9]+:\\\x20warning:\\\x20(?<message>.*)\\\x20\\\[(?<id>.*)\\\]
sonar.cxx.gcc.reportPaths=chameleon_build.log
sonar.cxx.xunit.reportPaths=*junit.xml
sonar.cxx.cobertura.reportPaths=chameleon_coverage.xml
sonar.cxx.cppcheck.reportPaths=chameleon_cppcheck.xml
sonar.cxx.clangsa.reportPaths=build-openmp/analyzer_reports/*/*.plist, build-parsec/analyzer_reports/*/*.plist, build-quark/analyzer_reports/*/*.plist, build-starpu/analyzer_reports/*/*.plist, build-starpu_simgrid/analyzer_reports/*/*.plist
sonar.cxx.jsonCompilationDatabase=compile_commands.json
EOF
echo "====== sonar-project.properties ============"
cat sonar-project.properties
echo "============================================"

# run sonar analysis + publish on sonarqube
sonar-scanner -X > sonar.log
