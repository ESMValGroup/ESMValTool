#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
set -eux

FOLDER_NAME="${OUTPUT_DIR}/${RECIPE_NAME}*"
KGO="${KGO_ROOT_PATH}/${RECIPE_NAME}*"

# Create a variable that defines the path to the compare script.
COMPARE_SCRIPT="${ESMVALTOOL_DIR}/esmvaltool/utils/testing/regression/compare.py"

# Run the compare script.
python "${COMPARE_SCRIPT}" ${KGO} ${FOLDER_NAME}
