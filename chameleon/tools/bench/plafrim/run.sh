#!/usr/bin/env bash
#
# @file run.sh
#
# @copyright 2020-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
# @version 1.4.0
# @author Florent Pruvost
# @date 2025-12-19
#
set -x

# Unset the binding environment of the CI for this specific case
unset STARPU_MPI_NOBIND
unset STARPU_WORKERS_NOBIND

# to avoid a lock during fetching chameleon branch in parallel
export XDG_CACHE_HOME=/tmp/guix-$$

# take advantage of the existing local Guix cache
mkdir -p $XDG_CACHE_HOME
cp -ar $HOME/.cache/guix $XDG_CACHE_HOME

# save guix commits
#guix describe --format=json > guix.json
guix time-machine \
  --channels=./tools/bench/guix-channels.scm -- describe --format=json > guix.json

# define env var depending on the node type
if [[ "$NODE" == "bora" ]]; then
  export SLURM_CONSTRAINTS="bora,omnipath"
  export CHAMELEON_BUILD_OPTIONS="-DCHAMELEON_USE_MPI=ON -DCMAKE_BUILD_TYPE=Release"
  export STARPU_HOSTNAME="bora"
elif [[ "$NODE" == "sirocco" ]]; then
  export SLURM_CONSTRAINTS="sirocco,v100"
  export CHAMELEON_BUILD_OPTIONS="-DCHAMELEON_USE_MPI=ON -DCHAMELEON_USE_CUDA=ON -DCMAKE_BUILD_TYPE=Release"
  export STARPU_HOSTNAME="sirocco"
else
  echo "$0: Please set the NODE environment variable to bora or sirocco."
  exit -1
fi

# define env var and guix rule to use depending on the mpi vendor
GUIX_ENV="chameleon-mkl"
export MPI_OPTIONS=""
if [[ "$MPI" == "openmpi" ]]; then
  export MPI_OPTIONS="--bind-to board"
  GUIX_ENV_MPI=""
  GUIX_ADHOC_MPI=""
  if [[ "$NODE" == "sirocco" ]]; then
    GUIX_ENV="chameleon-cuda-mkl"
  fi
elif [[ "$MPI" == "nmad" ]]; then
  export MPI_OPTIONS="-DPIOM_DEDICATED=1 -DPIOM_DEDICATED_WAIT=1 hwloc-bind --cpubind machine:0"
  GUIX_ENV_MPI="--with-input=openmpi=nmad"
  GUIX_ADHOC_MPI="which gzip zlib tar inetutils util-linux procps nmad"
else
  echo "$0: Please set the MPI environnement variable to openmpi or nmad."
  exit -1
fi
GUIX_ADHOC="coreutils gawk grep hwloc jube perl slurm@23"
GUIX_RULE="-D $GUIX_ENV $GUIX_ENV_MPI $GUIX_ADHOC $GUIX_ADHOC_MPI"

# 1. Run benchmark
guix time-machine \
  --channels=./tools/bench/guix-channels.scm \
  -- shell --pure \
  --preserve="PLATFORM|NODE|^CI|^SLURM|^JUBE|^MPI|^STARPU|^CHAMELEON" \
  $GUIX_RULE \
  -- /bin/bash --norc ./tools/bench/plafrim/slurm.sh
err1=$?

# 2. Upload results
guix time-machine \
  --channels=./tools/bench/guix-channels-elk.scm \
  -- shell --pure --preserve="proxy$" \
  curl nss-certs python python-click python-certifi python-elasticsearch python-gitpython \
  -- python3 tools/bench/jube/add_result.py -e https://elasticsearch.bordeaux.inria.fr \
  -t hiepacs -p chameleon -m $MPI chameleon\-$NODE\-$MPI.csv
err2=$?

# 3. Clean up
rm -rf /tmp/guix-$$

# Final exit code = any non-zero status
final_err=$(( err1 || err2 ))
exit "$final_err"
