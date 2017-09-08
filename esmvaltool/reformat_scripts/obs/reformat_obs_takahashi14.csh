#!/usr/bin/env csh -eu
###############################################################################
## REFORMAT SCRIPT FOR THE TAKAHASHI14 OBSERVATIONAL DATA
###############################################################################
##
## Tier
##    Tier 2: other freely-available dataset.
##
## Source
##    http://cdiac.ornl.gov/ftp/oceans/NDP_094/
##
## Download and processing instructions
##    Select "TALK_TCO2_pCO2_GLOB_Grid_Dat.nc" 
##    Run this script (requires NCO, http://nco.sourceforge.net/)
##
## Caveats
##
## Modification history
##    20151111-A_laue_ax: written.
##
###############################################################################

set inpath="${ESMValTool_RAWOBSPATH}/Tier2/takahashi14"
set outpath="${ESMValTool_OBSPATH}/Tier2/takahashi14"

if (! -d $outpath) then
    mkdir -p $outpath
endif

set infile=$inpath/TALK_TCO2_pCO2_GLOB_Grid_Dat.nc
set outfile=$outpath/OBS_takahashi14_ocean_1_TO2Ms_talk_200501-200512.nc

if (-e $infile) then
    echo 'input file = '$infile
else
    echo 'error: input file not found - '$infile
    exit
endif

cdo remapbil,r360x180 $infile talk_monthly_ref_takahashi14_reg_2005-2005.nc
cdo chname,alkali,talk talk_monthly_ref_takahashi14_reg_2005-2005.nc  talk_monthly_ref_takahashi14_reg_2005-2005_chname.nc
ncap2 -s talk="talk*1.023e-3" talk_monthly_ref_takahashi14_reg_2005-2005_chname.nc talk_monthly_ref_takahashi14_reg_2005-2005_unit.nc
cdo chunit,"micromol/kg","mol m-3" talk_monthly_ref_takahashi14_reg_2005-2005_unit.nc $outfile

ncatted -O -a title,global,a,c,"Monthly surface climatological Total Alkalinity" $outfile 
ncatted -O -a source,global,a,c,"http://cdiac.ornl.gov/ftp/oceans/NDP_094/" $outfile
ncatted -O -a tier,global,a,c,"2" $outfile 
ncatted -O -a period,global,a,c,"2005" $outfile
ncatted -O -a reference,global,a,c,"Takahashi, T., Sutherland, S. C., Chipman, D. W., Goddard, J. G., Ho, C., Newberger, T., Sweeney, C., and Munro, D. R.: Climatological distributions of pH, pCO2, total CO2, alkalinity, and CaCO3 saturation in the global surface ocean, and temporal changes at selected locations, Mar. Chem., 164, 95-125, doi:10.1016/j.marchem.2014.06.004, 2014." $outfile

if (-e $outfile) then
    rm talk_monthly_ref_takahashi14_reg_2005-2005.nc
    rm talk_monthly_ref_takahashi14_reg_2005-2005_chname.nc
    rm talk_monthly_ref_takahashi14_reg_2005-2005_unit.nc
    echo 'created '$outfile
else
    echo 'error: no output written'
endif

