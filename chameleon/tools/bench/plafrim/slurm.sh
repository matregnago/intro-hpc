#!/usr/bin/env bash
#
# @file slurm.sh
#
# @copyright 2020-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
# @version 1.4.0
# @author Florent Pruvost
# @date 2024-04-30
#

# Check the environment
echo $PLATFORM
echo $NODE
env |grep ^CI
env |grep ^SLURM
env |grep ^JUBE_ID
env |grep ^MPI
env |grep ^STARPU
env |grep ^CHAMELEON

set -x

# Parameters of the Slurm jobs
TIME=02:00:00
PART=routage
NP=$SLURM_NP
CONS=$SLURM_CONSTRAINTS
EXCL=

# Submit jobs
export JOB_NAME=chameleon\-$NODE\-$MPI\-$NP
sbatch --wait --job-name="$JOB_NAME" --output="$JOB_NAME.out" --error="$JOB_NAME.err" --nodes=$NP --time=$TIME --partition=$PART --constraint=$CONS --exclude=$EXCL --exclusive --ntasks-per-node=1 --threads-per-core=1 ./tools/bench/chameleon_guix.sh
err=$?

cat chameleon\-$NODE\-$MPI\-$NP.err
cat chameleon\-$NODE\-$MPI\-$NP.out

# exit with error code from the sbatch command
exit $err
