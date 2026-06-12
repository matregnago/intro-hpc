#!/usr/bin/env bash
#
# @file pages.sh
#
# @copyright 2020-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
# @version 1.4.0
# @author Florent Pruvost
# @date 2024-12-03
#
set -e
set -x
CHAMELEON_SRC_DIR=${CHAMELEON_SRC_DIR:-$PWD}

mkdir -p tmp_fig
cd tmp_fig

## need to generate figures from last benchmarks
# get the most recent chameleon commit date for which benchs have been performed
commit_sha=`curl -X GET "https://elasticsearch.bordeaux.inria.fr/hiepacs-chameleon_perf/_search?pretty" -H 'Content-Type: application/json' -d'
{
  "query": { "match_all": {} },
  "sort": [
    { "Commit_date_chameleon": "desc" }
  ],
  "size": 1
}
' | grep "Commit_sha_chameleon" | awk '{ print $3" "$4}' |sed -e "s#,##g" |sed -e "s#\"##g" |sed -e "s# ##g"`
echo $commit_sha

# generate the csv file from elasticsearch for the given chameleon commit
python3 ${CHAMELEON_SRC_DIR}/tools/bench/jube/get_result.py -e https://elasticsearch.bordeaux.inria.fr -t hiepacs -p chameleon -c $commit_sha

# generate the figures
Rscript ${CHAMELEON_SRC_DIR}/tools/bench/jube/GenFigures.R

# add the performance chapter to the doc. we need the performances files to add
# this chapter that are not commited in the sources so that this chapter is not
# here by default
cat >> ${CHAMELEON_SRC_DIR}/doc/user/users_guide.org.in <<EOF
** Chameleon Performances on PlaFRIM
:PROPERTIES:
:CUSTOM_ID: doc-perf
:END:
Chameleon commit: *$commit_sha*.
#+INCLUDE: @CMAKE_CURRENT_SOURCE_DIR@/chapters/performances.org
EOF

cd ..

## Build the doc
VERSION=${VERSION:-pages}
cmake -S ${CHAMELEON_SRC_DIR} -B build-$VERSION -DCHAMELEON_ENABLE_DOC=ON
cmake --build build-$VERSION --target doc --verbose -j5

## Copy files in public/ used as an artefact (zip archive) to upload on gitlab pages, see
mkdir -p public/dev/

## Homepage (user's guide): https://solverstack.gitlabpages.inria.fr/chameleon/index.html
cp build-$VERSION/doc/user/*.html public/
cp build-$VERSION/doc/user/*.png public/
cp build-$VERSION/doc/user/*.jpg public/
cp build-$VERSION/doc/user/*.svg public/

## API (doxygen): https://solverstack.gitlabpages.inria.fr/chameleon/dev/index.html
cp -r build-$VERSION/doc/dev/html/* public/dev/

## Images such as last performances on plafrim
cp tmp_fig/* public/

## lcov code coverage: https://solverstack.gitlabpages.inria.fr/chameleon/coverage/
cp -r coverage public/
