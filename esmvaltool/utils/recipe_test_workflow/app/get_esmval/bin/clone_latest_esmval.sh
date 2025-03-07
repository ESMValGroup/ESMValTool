#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
set -eux

# Remove the ESMValTool and ESMValCore directories, if they exist.
if [[ -d ${ESMVALTOOL_DIR} ]]; then
    rm -rf "${ESMVALTOOL_DIR}"
fi

if [[ -d ${ESMVALCORE_DIR} ]]; then
    rm -rf "${ESMVALCORE_DIR}"
fi

# Checkout the specified branch for ESMValCore and ESMValTool. Use the
# quiet ('-q') option to prevent the progress status from being reported
# (which is done via done via 'stderr').
git clone -q -b "${BRANCH}" "${ESMVALTOOL_URL}" "${ESMVALTOOL_DIR}"
git clone -q -b "${BRANCH}" "${ESMVALCORE_URL}" "${ESMVALCORE_DIR}"
