#!/bin/csh

set list = `find ../../interface_scripts/ ../../diag_scripts/ ../cmorizers/obs/ -type f -name "*.ncl" -exec grep 'undef(' {} \; | awk -F '"' '{print $2}'`

echo $#list

set tM = '\\\|'
set str = ""

foreach ll ($list)
  set str = $str$ll$tM
end

echo $str
