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

# Checkout the specified branch for ESMValCore and ESMValTool to a `src` dir. Use the
# quiet ('-q') option to prevent the progress status from being reported
# (which is done via done via 'stderr').
git clone -q -b "${BRANCH}" "${ESMVALTOOL_URL}" "${ESMVALTOOL_SOURCE_DIR}"
git clone -q -b "${BRANCH}" "${ESMVALCORE_URL}" "${ESMVALCORE_SOURCE_DIR}"

# Copy ESMValCore and ESMValTool into 'lib/python', adding them to the PYTHONPATH.
# 'pip install' does this and also triggers `setuptools-scm` to increment the
# package version number. Otherwise `version` returns the last tagged
# release, not the current development version.
rtw-env pip install "${ESMVALTOOL_SOURCE_DIR}" --target="${ESMVALTOOL_DIR}" --no-deps
rtw-env pip install "${ESMVALCORE_SOURCE_DIR}" --target="${ESMVALCORE_DIR}" --no-deps
