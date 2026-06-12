#!/bin/sh
#
# @file find_sources.sh
#
# @copyright 2020-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
# @version 1.4.0
# @author Mathieu Faverge
# @author Florent Pruvost
# @date 2024-12-03
#
#set -x

SRCDIR_TO_ANALYZE="build-openmp/runtime/openmp build-parsec/runtime/parsec build-quark/runtime/quark build-starpu compute control coreblas example include runtime testing"

echo $PWD
rm -f filelist.txt
for dir in ${SRCDIR_TO_ANALYZE}
do
    find $dir -name '*\.[ch]' >> filelist.txt
done

# Remove all CMakeFiles generated file
sed -i '/CMakeFiles/d' filelist.txt

# Remove installed files
sed -i '/build.*\/install.*/d' filelist.txt

# Remove all hmat files (external package)
sed -i '/coreblas\/hmat-oss/d' filelist.txt
sed -i '/testing\/test_fembem/d' filelist.txt

# Remove original files used for precision generation
for file in `git grep "@precisions" | awk -F ":" '{ print $1 }'`
do
    sed -i "\:^$file.*:d" filelist.txt
done

# Remove fortran header file that is incorrectly parsed as a C file
sed -i '/fortran\.h/d' filelist.txt

# Remove external header files
for file in coreblas/include/coreblas/cblas.h coreblas/include/coreblas/lapacke.h coreblas/include/coreblas/lapacke_config.h coreblas/include/coreblas/lapacke_mangling.h
do
    sed -i "\:^$file.*:d" filelist.txt
done

rm -f filelist_*.txt
for name in $(cat filelist.txt)
do
    test=$(grep "@generated" $name | wc -l)
    if [ $test -gt 0 ]
    then
        prec=$(grep "@generated" $name | sed 's/^.*[scdz] -> \([sdcz]\).*$/\1/')
        echo $name >> filelist_${prec}.txt
    else
        echo $name >> filelist_none.txt
    fi
done
