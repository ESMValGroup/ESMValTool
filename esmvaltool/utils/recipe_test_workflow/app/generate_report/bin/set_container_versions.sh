#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
set -eux

echo Running new bash script

export CAPTURE_OUTPUT=True
export QUIET_MODE=True
export YESTERDAYS_CONTAINER_PATH="${ROSE_DATACPT10M}/container/esmvaltool.sif"

export ESMVAL_VERSIONS_CURRENT=$(env-file esmvaltool version)

if [[ -f ${YESTERDAYS_CONTAINER_PATH} ]]; then
    export ESMVAL_VERSIONS_PREVIOUS=$(CONTAINER_PATH=${YESTERDAYS_CONTAINER_PATH} env-file esmvaltool version)
fi

echo ESMVal Current Versions:  $ESMVAL_VERSIONS_CURRENT
echo ESMVal Previous Versions: $ESMVAL_VERSIONS_PREVIOUS
