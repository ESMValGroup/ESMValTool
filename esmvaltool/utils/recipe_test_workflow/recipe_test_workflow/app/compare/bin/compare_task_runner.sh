#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
set -eux

FOLDER_NAME="${OUTPUT_DIR}/recipe_${KGO_METRIC}*"
KGO="${KGO_ROOT_PATH}/recipe_${KGO_METRIC}*"

# Create a variable that defines the path to the compare script.
COMPARE_SCRIPT=${SSS_TAG_DIR}/lib/python3.10/site-packages/esmvaltool/utils/testing/regression/compare.py

# Run the compare script.
python "${COMPARE_SCRIPT}" ${FOLDER_NAME} ${KGO}
