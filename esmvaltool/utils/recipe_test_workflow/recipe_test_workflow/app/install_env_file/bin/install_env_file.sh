#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
set -eux

# Copy the site specific environment file to the 'share/bin' directory in the
# installed Cylc workflow (this directory is automatically added to the
# ${PATH} by Cylc).
SOURCE_PATH="${CYLC_WORKFLOW_RUN_DIR}/site/${SITE}-env"
TARGET_DIR="${CYLC_WORKFLOW_RUN_DIR}/share/bin"
ENV_FILE="rtw-env"

# Create the 'share/bin' directory in the installed workflow.
mkdir "${TARGET_DIR}"

# Copy the environment file to the 'share/bin' directory.
cp "${SOURCE_PATH}" "${TARGET_DIR}/${ENV_FILE}"
