#!/usr/bin/env csh -eu
###############################################################################
## REFORMAT SCRIPT FOR THE SeaWIFS OBSERVATIONAL DATA
###############################################################################
##
## Tier
##    Tier 2: other freely-available dataset.
##
## Source
##    http://cdiac.ornl.gov/ftp/oceans/SOCATv2/SOCATv2_Gridded_Dat/SOCAT_tracks_gridded_monthly_v2.nc.zip
##
## Last access
##
## Download and processing instructions
##    Download data and unpack zip file 
##    (unzip SOCAT_tracks_gridded_monthly_v2.nc.zip) 
##    Run this script (requires NCO, http://nco.sourceforge.net/)
##
## Caveats
##
## Modification history
##    20151111-A_laue_ax: written.
##
###############################################################################

set inpath="${ESMValTool_RAWOBSPATH}/Tier2/SOCAT"
set outpath="${ESMValTool_OBSPATH}/Tier2/SOCAT"

if (! -d $outpath) then
    mkdir -p $outpath
endif

set infile=$inpath/SOCAT_tracks_gridded_monthly_v2.nc
set outfile=$outpath/OBS_SOCAT_ocean_1_TO2M_spco2_197001-201112.nc

if (-e $infile) then
    echo 'input file = '$infile
else
    echo 'error: input file not found - '$infile
    exit
endif

ncks -v XLON,YLAT,TMNTH,FCO2_AVE_UNWTD $infile SOCAT_tracks_gridded_fco2_ave_unwtd.nc
ncrename -v FCO2_AVE_UNWTD,spco2 SOCAT_tracks_gridded_fco2_ave_unwtd.nc
ncrename -v TMNTH,time SOCAT_tracks_gridded_fco2_ave_unwtd.nc
ncrename -d TMNTH,time SOCAT_tracks_gridded_fco2_ave_unwtd.nc
ncrename -d YLAT,lat SOCAT_tracks_gridded_fco2_ave_unwtd.nc
ncrename -v YLAT,lat SOCAT_tracks_gridded_fco2_ave_unwtd.nc
ncrename -v XLON,lon SOCAT_tracks_gridded_fco2_ave_unwtd.nc
ncrename -d XLON,lon SOCAT_tracks_gridded_fco2_ave_unwtd.nc
ncatted -a calendar,time,c,c,noleap SOCAT_tracks_gridded_fco2_ave_unwtd.nc
ncatted -O -a title,global,a,c,"Surface Ocean CO2 Atlas (SOCAT version 2)" SOCAT_tracks_gridded_fco2_ave_unwtd.nc 
ncatted -O -a source,global,a,c,"http://cdiac.ornl.gov/ftp/oceans/SOCATv2/SOCATv2_Gridded_Dat/SOCAT_tracks_gridded_monthly_v2.nc.zip" SOCAT_tracks_gridded_fco2_ave_unwtd.nc 
ncatted -O -a tier,global,a,c,"2" SOCAT_tracks_gridded_fco2_ave_unwtd.nc 
ncatted -O -a period,global,a,c,"1970-2011" SOCAT_tracks_gridded_fco2_ave_unwtd.nc
ncatted -O -a reference,global,a,c,"Bakker, D.C.E., et al., An update to the Surface Ocean CO2 Atlas (SOCAT version 2), Earth Syst. Sci. Data, 6, 69-90, doi:10.5194/essd-6-69-2014, 2014." SOCAT_tracks_gridded_fco2_ave_unwtd.nc 
mv SOCAT_tracks_gridded_fco2_ave_unwtd.nc $outfile

if (-e $outfile) then
    echo 'created '$outfile
else
    echo 'error: no output written'
endif

