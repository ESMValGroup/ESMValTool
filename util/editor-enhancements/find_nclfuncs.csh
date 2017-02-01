#!/bin/csh

# Interface
set list1 = `grep "undef(" ../../interface_scripts/*.ncl | awk -F '"' '{print $2}'`
# Libraries
set list2 = `grep "undef(" ../../diag_scripts/lib/ncl/*.ncl | awk -F '"' '{print $2}'`
# Diag aux
set list3 = `grep -ir --include="*.ncl" "undef(" ../../diag_scripts/aux/ | awk -F '"' '{print $2}'`
# Plot scripts
set list4 = `grep "undef(" ../../plot_scripts/ncl/*.ncl | awk -F '"' '{print $2}'`
# Reformat
set list5 = `grep "undef(" ../../reformat_scripts/default/reformat_default_func.ncl | awk -F '"' '{print $2}'`
set list6 = `grep "undef(" ../../reformat_scripts/EMAC/reformat_EMAC_func.ncl | awk -F '"' '{print $2}'`
set list7 = `grep "undef(" ../../reformat_scripts/obs/reformat_obs_func.ncl | awk -F '"' '{print $2}'`
set list8 = `grep "undef(" ../../reformat_scripts/constants.ncl | awk -F '"' '{print $2}'`

set list = ($list1 $list2 $list3 $list4 $list5 $list6 $list7 $list8 "calculate")

set tM = '\\\|'

set str = ""

foreach ll ($list)
    set str = $str$ll$tM
end

echo $str
