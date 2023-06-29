#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
# set -eux

export CONDA_PKGS_DIRS=${CYLC_WORKFLOW_SHARE_DIR}/rtw/conda/pkgs/
export CONDA_ENVS_PATH=${CYLC_WORKFLOW_SHARE_DIR}/rtw/conda/envs/

conda create -y --name rtw-env --override-channels -c conda-forge mamba

eval "$(conda shell.bash hook)"
conda activate rtw-env

mamba env update --file ${ESMVALTOOL_DIR}/environment.yml
