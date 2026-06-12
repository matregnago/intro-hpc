#!/usr/bin/env bash

success=1

check_rebase()
{
    hash_master=$(git show-ref -s origin/master)
    hash_common=$(git merge-base origin/master ${CI_COMMIT_SHA})
    if [ "${hash_master}" = "${hash_common}" ]
    then
        echo "check_rebase: OK"
        return 0
    else
        echo "check_rebase: Rebase is required"
        success=0
        return 1
    fi
}

check_draft()
{
    if [ "${CI_PIPELINE_SOURCE}" = "merge_request_event" ]
    then
        draft=$( echo ${CI_MERGE_REQUEST_TITLE} | sed "s/^Draft.*$/Draft/" )
        wip=$( echo ${CI_MERGE_REQUEST_TITLE} | sed "s/^WIP.*$/WIP/" )

        if [ "$draft" = "Draft" ]
        then
            echo "check_draft: Merge request is in draft mode"
            success=0
            return 1
        fi

        if [ "$wip" = "WIP" ]
        then
            echo "check_draft: Merge request is in WIP"
            success=0
            return 1
        fi

        # if [ "${CI_MERGE_REQUEST_APPROVED}" != "true" ]
        # then
        #     echo "check_approval: Merge request not yet approved"
        #     success=0
        #     return 1
        # fi

        echo "check_draft: Merge request is ok"
    fi

    return 0
}

check_header()
{
    echo " Checking file headers: "
    TOOLSDIR=$(dirname $0)/../tools

    $TOOLSDIR/check_header.sh
    rc=$?
    if [ $rc -eq 0 ]
    then
        echo "Check header: SUCCESS"
    else
        echo "Check header: FAILED"
        success=0
    fi
}

if [ $# -lt 1 ]
then
    echo "Usage: $0 [rebase|draft|header]"
    exit 1
fi

echo ""
echo "----------------------------------------------------"
case $1 in
    rebase)
	check_rebase
	;;
    draft)
	check_draft
	;;
    header)
	check_header
	;;
    *)
        echo "Usage: $0 [rebase|draft|header]"
	exit 1
	;;
esac

if [ $success -eq 0 ]
then
    exit 1
    # We could cancel the job, but then the log is not pushed in time to the web interface
    #curl --request POST --header "PRIVATE-TOKEN: ${PIPELINE_TOKEN}" "https://gitlab.inria.fr/api/v4/projects/$CI_PROJECT_ID/pipelines/$CI_PIPELINE_ID/cancel"
fi
