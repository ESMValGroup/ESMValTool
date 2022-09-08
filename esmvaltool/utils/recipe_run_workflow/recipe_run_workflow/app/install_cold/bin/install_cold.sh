#!/bin/bash
set -eux

# Remove source directory, if it exists.
if [[ -d ${SOURCE_DIR} ]]; then
    rm -rf "${SOURCE_DIR}"
fi

# Checkout uthpy.
git clone -b "${BRANCH}" "${URL}" "${SOURCE_DIR}"

# Install uthpy.
cd "${SOURCE_DIR}"
pip install . --prefix="${UTH_INSTALL_ROOT}"
