#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
set -eux

FOLDER_NAME="${OUTPUT_DIR}/recipe_${KGO_METRIC}*"
KGO="${KGO_ROOT_PATH}/recipe_${KGO_METRIC}*"

# Create a variable that defines the path to the compare script.
COMPARE_SCRIPT="${ESMVALTOOL_DIR}/esmvaltool/utils/testing/regression/compare.py"

# Pick up local installation of imagehash. Slightly fragile: see #3436
export PYTHONPATH=${PYTHONPATH}:${HOME}/.local/lib/python3.11/site-packages

# Run the compare script.
python "${COMPARE_SCRIPT}" ${KGO} ${FOLDER_NAME}
