#!/bin/bash

tools=(R Rscript julia)

for tool in "${tools[@]}"; do
    if ! command -v "$tool" > /dev/null 2>&1; then
        echo "Executable $tool not found! Exiting..." >> $PREFIX/.messages.txt
        exit 1
    fi
done

Rscript $PREFIX/lib/python*/site-packages/esmvaltool/install/R/setup.R >> $PREFIX/.messages.txt
julia $PREFIX/lib/python*/site-packages/esmvaltool/install/Julia/setup.jl >> $PREFIX/.messages.txt
