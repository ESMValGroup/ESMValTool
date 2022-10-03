#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
set -eux


# Copy the site specific environment launch file to the cylc run /bin
# Where it will be identified by rose task-run
SOURCE_PATH="${CYLC_WORKFLOW_RUN_DIR}/site/${SITE}-env"
TARGET_DIR="${CYLC_WORKFLOW_RUN_DIR}/bin"
ENV_FILE="rtw-env"

# Create the 'bin' directory in the installed workflow.
mkdir "${TARGET_DIR}"

# Copy the environment file to the 'bin' directory.
cp "${SOURCE_PATH}" "${TARGET_DIR}/${ENV_FILE}"

# Remove esmvaltool and esmvalcore directories, if they exist.
if [[ -d ${ESMVALTOOL_DIR} ]]; then
    rm -rf "${ESMVALTOOL_DIR}"
fi

if [[ -d ${ESMVALCORE_DIR} ]]; then
    rm -rf "${ESMVALCORE_DIR}"
fi

# Checkout main branch for ESMValTool and ESMValCore from github.
git clone -q -b "${BRANCH}" "${ESMVALTOOL_URL}" "${ESMVALTOOL_DIR}"
git clone -q -b "${BRANCH}" "${ESMVALCORE_URL}" "${ESMVALCORE_DIR}"
