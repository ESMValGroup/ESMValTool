#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
#set -eux

FOLDER_NAME=$(find "${OUTPUT_DIR}" -type d -name "recipe_${KGO_METRIC}*")

# Add cd to python installation here


KGO=$(find "${KGO_ROOT_PATH}" -type d -name "recipe_${KGO_METRIC}*")
# Create a variable that defines the path to the compare script.
COMPARE_SCRIPT=${SSS_TAG_DIR}/lib/python3.10/site-packages/esmvaltool/utils/testing/regression/compare.py

# Run the compare script.
python "${COMPARE_SCRIPT}" "${FOLDER_NAME}" "${KGO}"
