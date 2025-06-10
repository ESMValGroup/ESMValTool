#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
set -eu

YESTERDAYS_CONTAINER_PATH="${ROSE_DATACPT10M}/container/esmvaltool.sif"

ESMVAL_VERSIONS_CURRENT=$(QUIET_MODE=True CAPTURE_OUTPUT=True env-file esmvaltool version)
export ESMVAL_VERSIONS_CURRENT

if [[ -f ${YESTERDAYS_CONTAINER_PATH} ]]; then
    ESMVAL_VERSIONS_PREVIOUS=$(QUIET_MODE=True CAPTURE_OUTPUT=True CONTAINER_PATH=${YESTERDAYS_CONTAINER_PATH} env-file esmvaltool version)
    export ESMVAL_VERSIONS_PREVIOUS
fi
