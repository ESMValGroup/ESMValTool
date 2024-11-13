#!/bin/sh
set -eo pipefail
unset PYTHONPATH

. ~/conda/etc/profile.d/conda.sh
conda activate esmvaltool_v2.3

esmvaltool run "$1"
