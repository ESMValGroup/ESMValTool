#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
set -eux

# Remove source directory, if it exists.
if [[ -d ${SOURCE_DIR} ]]; then
    rm -rf "${SOURCE_DIR}"
fi

# Checkout ESMValTool.
#git clone -b "${BRANCH}" "${ESMVALTOOL_URL}" "${SOURCE_DIR}"

# Install ESMValTool.
#cd "${SOURCE_DIR}"
#pip install . --prefix="${RRW_INSTALL_ROOT}"

SOURCE_PATH="${CYLC_WORKFLOW_RUN_DIR}/site/${SITE}-env"
TARGET_DIR="${CYLC_WORKFLOW_RUN_DIR}/bin"
ENV_FILE="rrw-env"

# Create the 'bin' directory in the installed workflow.
mkdir "${TARGET_DIR}"

# Copy the environment file to the 'bin' directory.
cp "${SOURCE_PATH}" "${TARGET_DIR}/${ENV_FILE}"
