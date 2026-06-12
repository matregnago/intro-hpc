#!/usr/bin/env bash
###
#
#  @file coverage.sh
#  @copyright 2013-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.4.0
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 2024-02-18
#
###
#
# Group all the coverage files together and display the summary
# to produce the final chameleon.lcov file.
#
###
INPUT_FILES=""
for name in $( ls -1 chameleon*.lcov | grep -v simgrid)
do
    INPUT_FILES="$INPUT_FILES -a $name";
done
lcov $INPUT_FILES -o chameleon.lcov
lcov --summary chameleon.lcov
genhtml -o coverage chameleon.lcov
lcov_cobertura chameleon.lcov --output chameleon_coverage.xml
