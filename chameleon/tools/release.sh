#!/usr/bin/env bash
###
#
#  @file release.sh
#  @copyright 2013-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @brief Script to generate the release when pushing a tag
#
#  @version 1.4.0
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @date 2025-12-19
#
###

#
# Steps to do the release:
#    - Update information in the code (see update_release.sh)
#    - Update the ChangeLog
#    - Create a tag named vx.x.x and push it on solverstack (will trigger the CI to generate the release)
#
changelog=""
function gen_changelog()
{
    local firstline=$( grep -n "^chameleon-" ChangeLog | head -n 1 | cut -d ':' -f 1 )
    firstline=$(( firstline + 2 ))
    #echo $firstline
    local lastline=$( grep -n "^chameleon-" ChangeLog | head -n 2 | tail -n 1 | cut -d ':' -f 1 )
    lastline=$(( lastline - 1 ))
    #echo $lastline

    changelog="Changes:\n"
    for i in `seq $firstline $lastline`
    do
        local line=$( head -n $i ChangeLog | tail -n 1 )
        changelog="$changelog$line\\n"
        #echo $line
    done

    changelog="$changelog\nWARNING: Download the source archive by clicking on the link __Download release__ above, please do not consider the automatic Source code links as they are missing the submodules.\nVisit the [documentation](https://solverstack.gitlabpages.inria.fr/chameleon/#quickstart-install) to see how to install Chameleon."
}

release=""
function get_release()
{
    local firstline=$( grep -n "^chameleon-" ChangeLog | head -n 1 | cut -d ':' -f 1 )
    release=$( head -n $firstline ChangeLog | tail -n 1 | sed 's/chameleon\-//' )
}

set -x

# Get the release name through the branch name, and through the ChangeLog file.
# Both have to match to be correct
RELEASE_NAME=`echo $CI_COMMIT_TAG | cut -d v -f 2`
get_release

if [ -z "$RELEASE_NAME" -o -z "$release" -o "$RELEASE_NAME" != "$release" ]
then
    echo "Commit name $RELEASE_NAME is different from ChangeLog name $release"
    exit 1
fi

# generate the archive
wget https://raw.githubusercontent.com/Kentzo/git-archive-all/master/git_archive_all.py
python3 git_archive_all.py --force-submodules chameleon-$RELEASE_NAME.tar.gz

# upload the source archive to the Gitlab's Package registry
curl --header "JOB-TOKEN: $CI_JOB_TOKEN" --upload-file ./chameleon-$RELEASE_NAME.tar.gz "https://gitlab.inria.fr/api/v4/projects/$CI_PROJECT_ID/packages/generic/source/$CI_COMMIT_TAG/chameleon-$RELEASE_NAME.tar.gz"

# extract the change log from ChangeLog
gen_changelog
echo $changelog

# Try to remove the release if it already exists
curl --request DELETE --header "JOB-TOKEN: $CI_JOB_TOKEN" https://gitlab.inria.fr/api/v4/projects/$CI_PROJECT_ID/releases/v$RELEASE_NAME

# create the release associated to the tag
COMMAND=`echo curl --header \"Content-Type: application/json\" --header \"JOB-TOKEN: $CI_JOB_TOKEN\" \
  --data \'{ \"name\": \"v$RELEASE_NAME\", \
            \"tag_name\": \"v$RELEASE_NAME\", \
            \"ref\": \"$CI_COMMIT_REF_NAME\", \
            \"description\": \"$changelog\", \
            \"assets\": { \"links\": [{ \"name\": \"Download release chameleon-$RELEASE_NAME.tar.gz\", \"url\": \"https://gitlab.inria.fr/api/v4/projects/$CI_PROJECT_ID/packages/generic/source/$CI_COMMIT_TAG/chameleon-$RELEASE_NAME.tar.gz\" }] } }\' \
  --request POST https://gitlab.inria.fr/api/v4/projects/$CI_PROJECT_ID/releases`
eval $COMMAND
