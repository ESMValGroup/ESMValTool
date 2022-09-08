#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
set -eux

# Remove source directories, if they exists.
if [[ -d ${ESMVALTOOL_DIR} ]]; then
    rm -rf "${ESMVALTOOL_DIR}"
fi

if [[ -d ${ESMVALCORE_DIR} ]]; then
    rm -rf "${ESMVALCORE_DIR}"
fi

SOURCE_PATH="${CYLC_WORKFLOW_RUN_DIR}/site/${SITE}-env"
TARGET_DIR="${CYLC_WORKFLOW_RUN_DIR}/bin"
ENV_FILE="rrw-env"

# Create the 'bin' directory in the installed workflow.
mkdir "${TARGET_DIR}"

# Copy the environment file to the 'bin' directory.
cp "${SOURCE_PATH}" "${TARGET_DIR}/${ENV_FILE}"

# Checkout ESMValTool.
git clone -b "${BRANCH}" "${ESMVALTOOL_URL}" "${ESMVALTOOL_DIR}"
git clone -b "${BRANCH}" "${ESMVALCORE_URL}" "${ESMVALCORE_DIR}"

# Install ESMValTool.
cd "${ESMVALTOOL_DIR}"
# pip install . --prefix="${CYLC_WORKFLOW_RUN_DIR}"
