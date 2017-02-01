#!/bin/csh
#;;#############################################################################
#;;
#;; Tier
#;;    Tier 2: freely available non-obs4mips datasets
#;;
#;; Source
#;;    Download at
#;;    http://cdiac.ornl.gov/ftp/oceans/SOCATv2/SOCATv2_Gridded_Dat/SOCAT_tracks_gridded_monthly_v2.nc.zip
#;;
#;; Download and processing instructions
#;;    *) Download data and unpack zip file (unzip SOCAT_tracks_gridded_monthly_v2.nc.zip) 
#;;    *) Run this script (requires NCO, http://nco.sourceforge.net/)
#;;
#;; Caveats
#;;
#;; Modification history
#;;    20151111-A_laue_ax: written.
#;;
#;;#############################################################################

set inpath=/export/pa_data01/ESMVal/obs/RAW/Tier2/SOCAT
set outpath=/export/pa_data02/ESMVal/obs/Tier2/SOCAT

if (! -d $outpath) then
    mkdir $outpath
endif

set infile=$inpath/SOCAT_tracks_gridded_monthly_v2.nc
set outfile=$outpath/spco2_monthly_ref_SOCAT_reg_1970-2011.nc

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

