#!/bin/bash

# This script downloads MERRA-2 data
#
# Instructions:
#  - Create an Earthdata account: https://urs.earthdata.nasa.gov/
#  - Link your Earthdata account to 'GES DISC' by following instructions given here:
#    https://disc.gsfc.nasa.gov/earthdata-login
#  - Setup download with wget: https://disc.gsfc.nasa.gov/data-access#mac_linux_wget
#  - Change the destination directory specified below
#  - Check the years
#  - Run this script
# 
# Additional resources:
#    - file specification document is available here: https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf)



# The destination for the files to be downloaded
destination=/net/exo/landclim/PROJECTS/C3S/datadir/rawobsdir/Tier3/MERRA-2/

# The years that should be downloaded.
yrs=({1980..2019})

for yr in ${yrs[@]}; do
    wget -r -nd -nH --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --directory-prefix=${destination} --auth-no-challenge=on --keep-session-cookies -A '*.nc4' https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/M2TMNXLND.5.12.4/${yr}/ -R "*.tmp"
done
