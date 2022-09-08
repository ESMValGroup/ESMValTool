#!/bin/bash
set -eux

# Remove source directory, if it exists.
if [[ -d ${SOURCE_DIR} ]]; then
    rm -rf "${SOURCE_DIR}"
fi
printf 'This is running.'
# Checkout ESMValTool.
#git clone -b "${BRANCH}" "${ESMVALTOOL_URL}" "${SOURCE_DIR}"

# Install ESMValTool.
#cd "${SOURCE_DIR}"
#pip install . --prefix="${RRW_INSTALL_ROOT}"
