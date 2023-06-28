#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
set -eux
source /etc/profile.d/conda.sh

export CONDA_PKGS_DIRS=/tmp/conda/pkgs/
export CONDA_ENVS_PATH=${CYLC_WORKFLOW_RUN_DIR}/conda/envs/

conda create -y --name rtw-env mamba

conda activate ${CYLC_WORKFLOW_RUN_DIR}/conda/envs/rtw-env

mamba env update -y --file ${ESMVALTOOL_DIR}/environment.yml