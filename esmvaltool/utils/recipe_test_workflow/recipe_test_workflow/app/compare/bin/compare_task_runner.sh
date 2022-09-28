#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
#set -eux

FOLDER_NAME=$(find "${OUTPUT_DIR}" -type d -name "recipe_${KGO_METRIC}*")

# Fix this line later
cd /net/project/ukmo/scitools/opt_scitools/conda/environments/esmvaltool-2.6.0/lib/python3.10/site-packages/esmvaltool


KGO=$(find "${KGO_PATH}" -type d -name "recipe_${KGO_METRIC}*")
python utils/testing/regression/compare.py "${FOLDER_NAME}" "${KGO}"
