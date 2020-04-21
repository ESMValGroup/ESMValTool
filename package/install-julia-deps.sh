#!/bin/bash

if ! command -v julia > /dev/null 2>&1
then
    echo "Executable julia not found! Exiting..." >> "$PREFIX/.messages.txt"
    exit 1
else
    julia "$PREFIX"/lib/python*/site-packages/esmvaltool/install/Julia/setup.jl >> "$PREFIX/.messages.txt" 2>&1
fi
