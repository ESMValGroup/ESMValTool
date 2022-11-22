"""ESMValTool CMORizer for NOAA-CIRES-20CR data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://psl.noaa.gov/data/gridded/data.ncep.reanalysis2.html

Last access
    20220906

Download and processing instructions
    To facilitate the download, the links to the ftp server are provided.
    Since the filenames are sometimes identical across different
    save the data in three subdirectories in input_dir_path.

    ftp://ftp.cdc.noaa.gov/Projects/20thC_ReanV2/Monthlies/

    Subdirectory surface/:
        pr_wtr.eatm.mon.mean.nc
        cldwtr.eatm.mon.mean.nc

    Subdirectory surface_gauss/:
        tcdc.eatm.mon.mean.nc
        ulwrf.ntat.mon.mean.nc
        uswrf.ntat.mon.mean.nc

    Subdirectory pressure/:
        shum.mon.mean.nc

Caveats

"""
from .ncep_ncar_r1 import cmorization

# The following line makes it clear that the above import is not an error
cmorization
