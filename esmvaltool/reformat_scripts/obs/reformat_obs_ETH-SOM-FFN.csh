#!/usr/bin/env csh -eu
###############################################################################
## REFORMAT SCRIPT FOR ETH OBSERVATIONAL DATA
###############################################################################
##
## Tier
##    Tier 2: other freely-available dataset.
##
## Source
##    http://cdiac.ornl.gov/ftp/oceans/spco2_1998_2011_ETH_SOM-FFN/
##
## Last access
##
## Download and processing instructions
##    Download the file spco2_1998-2011_ETH_SOM-FFN_CDIAC_G05.nc.zip
##    unzip spco2_1998-2011_ETH_SOM-FFN_CDIAC_G05.nc.zip
##    Run this script (requires requires NCO, http://nco.sourceforge.net/ and
##    CDO, http://https://code.zmaw.de/projects/cdo/)
##
## Caveats
##
## Modification history
##    20151111-A_laue_ax: written.
##
###############################################################################

set inpath="${ESMValTool_RAWOBSPATH}/Tier2/ETH-SOM-FFN"
set outpath="${ESMValTool_OBSPATH}/Tier2/ETH-SOM-FFN"

if (! -d $outpath) then
    mkdir -p $outpath
endif

set infile=$inpath/spco2_1998-2011_ETH_SOM-FFN_CDIAC_G05.nc
set outfile="$outpath/OBS_ETH-SOM-FFN_ocean_1_TO2Ms_spco2_199801-201112.nc"

if (-e $infile) then
    echo 'input file = '$infile
else
    echo 'error: input file not found - '$infile
    exit
endif

ncks -v lon,lat,time,spco2_raw $infile tmp.nc
ncrename -v spco2_raw,spco2 tmp.nc
cdo chunit,muatm,uatm tmp.nc spco2_monthly_ref_ETH-SOM-FFN_reg_1998-2011_varname_unit.nc
cdo selvar,spco2 spco2_monthly_ref_ETH-SOM-FFN_reg_1998-2011_varname_unit.nc spco2_monthly_ref_ETH-SOM-FFN_reg_1998-2011_selvar.nc
ncatted -a origin,time,c,c,"seconds since 2000-01-01 00:00:00" spco2_monthly_ref_ETH-SOM-FFN_reg_1998-2011_selvar.nc spco2_monthly_ref_ETH-SOM-FFN_reg_1998-2011_origin.nc
ncks -O -x -v nb2,time_bnds spco2_monthly_ref_ETH-SOM-FFN_reg_1998-2011_origin.nc $outfile
###ncwa -a nb2 out.nc $outfile

ncatted -O -a title,global,a,c,"One Surface Ocean pCO2 Mapping Intercomparison (SOCOM) product: ETH-SOM-FFN" $outfile
ncatted -O -a source,global,a,c,"http://cdiac.ornl.gov/oceans/SPCO2_1998_2011_ETH_SOM_FFN.html" $outfile
ncatted -O -a tier,global,a,c,"2" $outfile
ncatted -O -a period,global,a,c,"1998-2011" $outfile
ncatted -O -a reference,global,a,c,"Landschützer, P., Gruber, N., Bakker, D.C.E., Schuster, U.: Recent variability of the global ocean carbon sink, Global Biogeochemical Cycles, 28, doi: 10.1002/2014GB004853, 2014" $outfile

if (-e $outfile) then
    rm -f tmp.nc
###    rm -f out.nc
    rm -f spco2_monthly_ref_ETH-SOM-FFN_reg_1998-2011_varname_unit.nc
    rm -f spco2_monthly_ref_ETH-SOM-FFN_reg_1998-2011_selvar.nc
    rm -f spco2_monthly_ref_ETH-SOM-FFN_reg_1998-2011_origin.nc
    echo 'created '$outfile
else
    echo 'error: no output written'
endif

