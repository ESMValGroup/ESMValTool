#!/bin/csh
###############################################################################
## REFORMAT SCRIPT FOR THE TAKAHASHI14 OBSERVATIONAL DATA
###############################################################################
##
## Tier
##    Tier 1: ops4mips dataset.
##
## Source
##    
##
## Download and processing instructions
##    
##    
##
## Caveats
##    - Does not jet use the ${ESMAVallTool_OBSPATH}!!!
##    - Script should be modified to process xco2 and xch4 in a loop!
##
## Modification history
##    20160714-A_wenz_sa: written.
##
###############################################################################

set inpath  = "/export/pa_data02/ESMVal/obs/Tier1/ESACCI-GHG"
set outpath = "/export/pa_data02/ESMVal/obs/Tier1/ESACCI-GHG"

if (! -d $outpath) then
    mkdir -p $outpath
endif

##set in and out file:
set infile   = $inpath/xco2_ghgcci_l3_v100_200301_201412.nc
set outfile  = $outpath/xco2_ESACCI-GHG_l3_v100_200301-201412.nc

##select variables from infile:
cdo selvar,xco2        $infile $outpath/temp_xco2_ESACCI-GHG.nc
cdo selvar,xco2_nobs   $infile $outpath/temp_xco2_nobs_ESACCI-GHG.nc
cdo selvar,xco2_stderr $infile $outpath/temp_xco2_stderr_ESACCI-GHG.nc
cdo selvar,xco2_stddev $infile $outpath/temp_xco2_stddev_ESACCI-GHG.nc

##change units of Stddev and Stderr:
ncap2 -s xco2_stderr="xco2_stddev*1.e6" $outpath/temp_xco2_stderr_ESACCI-GHG.nc $outpath/temp_xco2_stderr.nc
cdo chunit,"1","ppmv" $outpath/temp_xco2_stderr.nc $outpath/xco2_stderr_ESACCI-GHG_l3_v100_200301-201412.nc

ncap2 -s xco2_stddev="xco2_stddev*1.e6" $outpath/temp_xco2_stddev_ESACCI-GHG.nc $outpath/temp_xco2_stddev.nc
cdo chunit,"1","ppmv" $outpath/temp_xco2_stddev.nc $outpath/xco2_stddev_ESACCI-GHG_l3_v100_200301-201412.nc

##change var name to standard:
cdo chvar,xco2_nobs,xco2Nobs     $outpath/xco2_nobs_ESACCI-GHG_l3_v100_200301-201412.nc   $outpath/xco2Nobs_ESACCI-GHG_l3_v100_200301-201412.nc 
cdo chvar,xco2_stderr,xco2Stderr $outpath/xco2_stderr_ESACCI-GHG_l3_v100_200301-201412.nc $outpath/xco2Stderr_ESACCI-GHG_l3_v100_200301-201412.nc
cdo chvar,xco2_stddev,xco2Stddev $outpath/xco2_stddev_ESACCI-GHG_l3_v100_200301-201412.nc $outpath/xco2Stddev_ESACCI-GHG_l3_v100_200301-201412.nc

##remove temporary files:
if (-e $outfile) then
    rm $outpath/xco2_nobs_ESACCI-GHG_l3_v100_200301-201412.nc
    rm $outpath/xco2_stderr_ESACCI-GHG_l3_v100_200301-201412.nc
    rm $outpath/xco2_stddev_ESACCI-GHG_l3_v100_200301-201412.nc
    rm $outpath/temp_xco2_*.nv
    echo 'created '$outfile
else
    echo 'error: no output written'
endif


##do the same again for xh4:
set infile   = $inpath/xch4_ghgcci_l3_v100_200301_201412.nc
set outfile  = $outpath/xch4_ESACCI-GHG_l3_v100_200301-201412.nc

cdo selvar,xch4        $infile $outpath/xch4_ESACCI-GHG.nc
cdo selvar,xch4_nobs   $infile $outpath/xch4_nobs_ESACCI-GHG.nc
cdo selvar,xch4_stderr $infile $outpath/xch4_stderr_ESACCI-GHG.nc
cdo selvar,xch4_stddev $infile $outpath/xch4_stddev_ESACCI-GHG.nc

ncap2 -s xco2_stderr="xco2_stddev*1.e9" $outpath/temp_xco2_stderr_ESACCI-GHG.nc $outpath/temp_xco2_stderr.nc
cdo chunit,"1","ppmb" $outpath/temp_xco2_stderr.nc $outpath/xco2_stderr_ESACCI-GHG_l3_v100_200301-201412.nc

ncap2 -s xco2_stddev="xco2_stddev*1.e9" $outpath/temp_xco2_stddev_ESACCI-GHG.nc $outpath/temp_xco2_stddev.nc
cdo chunit,"1","ppmb" $outpath/temp_xco2_stddev.nc $outpath/xco2_stddev_ESACCI-GHG_l3_v100_200301-201412.nc

cdo chvar,xch4_nobs,xch4Nobs     $outpath/xch4_nobs_ESACCI-GHG_l3_v100_200301-201412.nc   $outpath/xch4Nobs_ESACCI-GHG_l3_v100_200301-201412.nc 
cdo chvar,xch4_stderr,xch4Stderr $outpath/xch4_stderr_ESACCI-GHG_l3_v100_200301-201412.nc $outpath/xch4Stderr_ESACCI-GHG_l3_v100_200301-201412.nc
cdo chvar,xch4_stddev,xch4Stddev $outpath/xch4_stddev_ESACCI-GHG_l3_v100_200301-201412.nc $outpath/xch4Stddev_ESACCI-GHG_l3_v100_200301-201412.nc


if (-e $outfile) then
    rm $outpath/xch4_nobs_ESACCI-GHG_l3_v100_200301-201412.nc
    rm $outpath/xch4_stderr_ESACCI-GHG_l3_v100_200301-201412.nc
    rm $outpath/xch4_stddev_ESACCI-GHG_l3_v100_200301-201412.nc
    rm $outpath/temp_xch4_*.nv
    echo 'created '$outfile
else
    echo 'error: no output written'
endif
