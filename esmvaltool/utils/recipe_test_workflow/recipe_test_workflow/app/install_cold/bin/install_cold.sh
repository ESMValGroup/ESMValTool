#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
set -eux

SOURCE_PATH="${CYLC_WORKFLOW_RUN_DIR}/site/${SITE}-env"
TARGET_DIR="${CYLC_WORKFLOW_RUN_DIR}/bin"
ENV_FILE="rtw-env"

# Create the 'bin' directory in the installed workflow.
mkdir "${TARGET_DIR}"

# Copy the environment file to the 'bin' directory.
cp "${SOURCE_PATH}" "${TARGET_DIR}/${ENV_FILE}"

# Manually and temporarily install imagehash (it will be included in the
# next ESMValTool community environment).
rtw-env pip install imagehash -t "${CYLC_WORKFLOW_RUN_DIR}/lib"
