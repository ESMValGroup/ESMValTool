#!/usr/bin/env bash 
###############################################################################
## REFORMAT SCRIPT FOR THE GlobAlbedo OBSERVATIONAL DATA
###############################################################################
##
## Tier
##    Tier 2: other freely-available datasets.
##
## Source
##    Download at
##    
##
## Download and processing instructions
##    
##
## Caveats
##
## Modification history
##    20170131-A_muel_bn: written.
##
###############################################################################

cdo chunit,"1","W m-2 / W m-2" -settunits,days -setreftime,1970-01-01,12:00 -settaxis,1998-01-15,12:00,1months -chname,DHR_SW,alb -selname,DHR_SW OBS_GlobAlbedo_sat_online_T2Ms_ORIG_199801-200512.nc OBS_GlobAlbedo_sat_black-online_T2Ms_alb_199801-200512.nc
