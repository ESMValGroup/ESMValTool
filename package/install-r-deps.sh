#!/bin/bash

if ! command -v Rscript > /dev/null 2>&1
then
    echo "Executable Rscript not found! Exiting..." >> "$PREFIX/.messages.txt"
    exit 1
else
    Rscript "$PREFIX"/lib/python*/site-packages/esmvaltool/install/R/setup.R >> "$PREFIX/.messages.txt" 2>&1
fi
