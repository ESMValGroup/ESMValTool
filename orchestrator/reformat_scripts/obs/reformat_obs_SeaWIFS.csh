#!/usr/bin/env csh -eu
###############################################################################
## REFORMAT SCRIPT FOR THE SeaWIFS OBSERVATIONAL DATA
###############################################################################
##
## Tier
##    Tier 2: other freely-available dataset.
##
## Source
##    Download original files at http://oceandata.sci.gsfc.nasa.gov/SeaWiFS/Mapped/Monthly/
##    Alternatively, download preprocessed SeaWiFS data at 
##    http://lgmacweb.env.uea.ac.uk/green_ocean/data/biogeochemistry/seawifs_chl.nc 
##
## Download and processing instructions (original SeaWiFS files)
##    Select "9 km" and "chl", then download all files (e.g., S19972441997273.L3m_MO_CHL_chl_ocx_9km.nc)
##    Data is missing for the following 3 time slots: 2008032-2008060, 2008061-2008091, 2009121-2009151
##    --> create symbolic links pointing to the empty file S2008001-2008031.L3m_MO_CHL_chl_ocx_9km.nc:
##      ln -s S20080012008031.L3m_MO_CHL_chl_ocx_9km.nc S20080322008060.L3m_MO_CHL_chl_ocx_9km.nc
##      ln -s S20080012008031.L3m_MO_CHL_chl_ocx_9km.nc S20080612008091.L3m_MO_CHL_chl_ocx_9km.nc
##      ln -s S20080012008031.L3m_MO_CHL_chl_ocx_9km.nc S20091212009151.L3m_MO_CHL_chl_ocx_9km.nc
##    Alternatively, the preprocessed SeaWiFS data can be downloaded from
##    http://lgmacweb.env.uea.ac.uk/green_ocean/data/biogeochemistry/seawifs_chl.nc
##    Then set variable "use_original_data" to "0" (in this script)
##    Run this script (requires NCO, http://nco.sourceforge.net/)
##
## Caveats
##    When processing the original SeaWiFs files, ncrename might produce the error message
##    "NetCDF: Operation not allowed in data mode". The reason for this is unknown but the
##    error message does not seem to have any relevant effect on the output.
##
## Modification history
##    20151111-A_laue_ax: written.
##
###############################################################################

set inpath="${ESMValTool_RAWOBSPATH}/Tier2/SeaWIFS"
set outpath="${ESMValTool_OBSPATH}/Tier2/SeaWIFS"

# ------------------------------
# user switch: use_original_data
# ------------------------------
# 0 = use preprocessed SeaWiFS data downloaded from
#     http://lgmacweb.env.uea.ac.uk/green_ocean/data/biogeochemistry/seawifs_chl.nc
# 1 = use original SeaWiFS data downloaded from
#     http://oceandata.sci.gsfc.nasa.gov/SeaWiFS/Mapped/Monthly/

set use_original_data=0

if (! -d $outpath) then
    mkdir -p $outpath
endif

set outfile=$outpath/OBS_SeaWIFS_ocean_1_TO2Ms_chl_199701-201012.nc

if ($use_original_data == "1") then
    # In case of using the original SeaWiFS files, a time dimension
    # has to be added before concatenating the files.

    set test=`ls -l $inpath/S??????????????.L3m_MO_CHL_chl_ocx_9km.nc | wc | awk '{print $1}'`
    if ($test != "160") then
        echo 'error: incorrect number of SeaWiFS input files (expected 160, found '$test')'
        echo 'check if symbolic links are set (see instructions in the header section of this script)'
        echo 'or if other input files are missing'
        exit
    endif

    foreach file ($inpath/S??????????????.L3m_MO_CHL_chl_ocx_9km.nc)
        echo $file

        # get time range from file name

        set df=`echo $file:t | cut -c6-8`     # first day
        set dl=`echo $file:t | cut -c13-15`   # last year
        set yf=`echo $file:t | cut -c2-5`     # first year
        set yl=`echo $file:t | cut -c9-12`    # last year

        # calculate time as yyyy.xxx (with days as fraction of a year)
        # assuming calendar = "no leap" (see below)

        set date=`echo $df $yf $dl $yl | awk '{print $2 + (($4-$2)*365+$1+$3-1)/2/365}'`
        echo $date

        ncecat -O -u time -v chl_ocx $file tmp.nc
        ncap2 -O -s 'defdim("time",-1);time[time]='${date}'' tmp.nc  $file:t:r_t.nc
    end

    ncrcat -O *_t.nc seawifs_chl_ncks.nc
    ncap -O -c -s "time=float(time-1997.71)*31536000.)" seawifs_chl_ncks.nc tmp.nc
    ncrename -v chl_ocx,chl tmp.nc

else # use preprocessed SeaWiFS data
    # remove unneeded "DEPTH" dimension
    ncwa -y rms -a DEPTH $inpath/seawifs_chl.nc seawifs_chl_ncks.nc
    ncrename -d LONGITUDE,lon seawifs_chl_ncks.nc
    ncrename -d LATITUDE,lat seawifs_chl_ncks.nc
    ncrename -d TIME,time seawifs_chl_ncks.nc
    ncrename -v TIME,time seawifs_chl_ncks.nc
    ncrename -v LONGITUDE,lon seawifs_chl_ncks.nc
    ncrename -v LATITUDE,lat seawifs_chl_ncks.nc
    ncap -O -c -s "time=float(time-1997.668)*31536000.)" seawifs_chl_ncks.nc tmp.nc
endif

ncatted -O -a units,time,c,c,"seconds since 1997-09-15 00:00:00" tmp.nc seawifs_chl_ncatted.nc
ncatted -O -a time_origin,time,c,c,"1997-SEP-15 00:00:00" seawifs_chl_ncatted.nc $outfile
ncatted -O -a calendar,time,c,c,noleap $outfile

ncatted -O -a title,global,a,c,"Chlorophyll concentration data from Sea-viewing WIde Field-of-view Sensor (SeaWIFS)" $outfile 
ncatted -O -a source,global,a,c,"http://oceancolor.gsfc.nasa.gov/SeaWiFS/" $outfile 
ncatted -O -a tier,global,a,c,"2" $outfile 
ncatted -O -a period,global,a,c,"1997-2010" $outfile
ncatted -O -a reference,global,a,c,"SeaWiFS Project, NASA Goddard Space Flight Center" $outfile 

if (-e $outfile) then
    if (-e tmp.nc) rm -f tmp.nc
    rm -f seawifs_chl_ncks.nc
    rm -f seawifs_chl_ncatted.nc
    echo 'created '$outfile
else
    echo 'error: no output written'
endif
