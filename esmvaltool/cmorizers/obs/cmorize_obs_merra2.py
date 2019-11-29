"""ESMValTool CMORizer for MERRA2 data.

Tier
    Tier 3: restricted datasets (i.e., dataset which requires a registration
 to be retrieved or provided upon request to the respective contact or PI).

Source
    https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels

Last access
    20191129

Download and processing instructions
    - Create an Earthdata account: https://urs.earthdata.nasa.gov/
    - Link your Earthdata account to 'GES DISC' by following instructions given here:
      https://disc.gsfc.nasa.gov/earthdata-login
    - Setup download with wget: https://disc.gsfc.nasa.gov/data-access#mac_linux_wget
    - Download with: 
        wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies --content-disposition https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_DIURNAL/M2TUNXLND.5.12.4/1980/MERRA2_100.tavgU_2d_lnd_Nx.198001.nc4

"""

