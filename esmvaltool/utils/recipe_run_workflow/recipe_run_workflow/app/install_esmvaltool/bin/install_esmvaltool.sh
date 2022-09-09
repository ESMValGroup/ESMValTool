#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
set -eux

# install esmvaltool from main branch
cd "${ESMVALTOOL_DIR}"
pip install . --prefix="${CYLC_WORKFLOW_RUN_DIR}"

# install esmvaltool from main branch
cd "${ESMVALCORE_DIR}"
pip install . --prefix="${CYLC_WORKFLOW_RUN_DIR}"

# this currently fails due to permissions issues when trying to install.
# currently assuming it is trying to overwrite something in the installed
# scitools module? (the launch file from task 1)
