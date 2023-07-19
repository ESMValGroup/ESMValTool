#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
# set -eux

conda create -y --name rtw-env --override-channels -c conda-forge mamba

eval "$(conda shell.bash hook)"
conda activate rtw-env

mamba env update --file ${ESMVALTOOL_DIR}/environment.yml

conda activate esmvaltool
cd $ESMVALCORE_DIR
pip install .
cd $ESMVALTOOL_DIR
pip install .
