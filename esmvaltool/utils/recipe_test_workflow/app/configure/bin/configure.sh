#!/bin/bash
# Create the ESMValTool configuration directory for this cycle.
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
set -eux

cp "${SITE_CONFIG_DIR}"/*.yml "${ESMVALTOOL_CONFIG_DIR}/"

if [[ -n "${DATA_SOURCES_FILE}" ]]; then
    env-file esmvaltool config copy "${DATA_SOURCES_FILE}" --overwrite
fi

# Fail here rather than in every 'process' task if the configuration is invalid.
env-file esmvaltool config show > /dev/null
