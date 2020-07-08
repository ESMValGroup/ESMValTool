#!/bin/bash

# This script downloads MERRA-2 data
#
# Instructions:
#  - Create an Earthdata account: https://urs.earthdata.nasa.gov/
#  - Link your Earthdata account to 'GES DISC' by following instructions given here:
#    https://disc.gsfc.nasa.gov/earthdata-login
#  - Setup download with wget: https://disc.gsfc.nasa.gov/data-access#mac_linux_wget
#  - Check the years
#  - Run this script as follow ./download_merra2.sh path_to_save_the_files
#
#
# Additional resources:
#    - file specification document is available here: https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf)

if [[ $# -eq 0 ]] ; then
    echo 'Please provide the download directory as first argument to this script.'
    exit 1
fi

# The years that should be downloaded.
yrs=({1980..2020})


# Monthly data
# for yr in ${yrs[@]}; do
#    echo $yr
#    wget -r -nd -nH --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --directory-prefix=$1 --auth-no-challenge=on --keep-session-cookies -A '*_lnd_*.nc4' https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/M2TMNXLND.5.12.4/${yr}/ -R "*.tmp"
# done

# Daily data
for yr in ${yrs[@]}; do
   echo $yr
   wget -r -nd -nH --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --directory-prefix=$1 --auth-no-challenge=on --keep-session-cookies -A '*_lnd_*.nc4' https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_DIURNAL/M2TUNXLND.5.12.4/${yr}/ -R "*.tmp"
done

# Also download constants file
# wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies directory-prefix=$1 --auth-no-challenge=on --keep-session-cookies https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/M2C0NXLND.5.12.4/1980/MERRA2_100.const_2d_lnd_Nx.00000000.nc4
