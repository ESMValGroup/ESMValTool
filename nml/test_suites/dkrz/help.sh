#!/bin/bash
cmd="/pf/b/b309070/util/translate2ESGF/translate2ESGF.py"
for i in namelist*Diur*xml
do
    echo $i
    echo python  $cmd '-n $i' $(grep -r '<model>' $i |sed -e "s/<model>/'/g" -e "s/<\/model>/'/g")
    ./args "$(grep -r '<model>' $i |sed -e "s/<model>/'/g" -e "s/<\/model>/'/g")"
done

