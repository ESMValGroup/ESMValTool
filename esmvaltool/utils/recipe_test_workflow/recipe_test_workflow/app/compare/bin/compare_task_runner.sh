#!/bin/bash
# Send the output from 'set -x' to 'stdout' rather than 'stderr'.
BASH_XTRACEFD=1
#set -eux

echo output dir is $OUTPUT_DIR

FOLDER_NAME=$(ls -d ${OUTPUT_DIR}/recipe_${KGO_NAME}*)


echo this is the folder name $folder_name
# Fix this line later
cd /net/project/ukmo/scitools/opt_scitools/conda/environments/esmvaltool-2.6.0/lib/python3.10/site-packages/esmvaltool

$KGO_PATH
$KGO_NAME
$FOLDER_NAME

KGO=$(ls -d ${KGO_PATH}/recipe_${KGO_NAME}*)
python utils/testing/regression/compare.py ${KGO} "${FOLDER_NAME}"
