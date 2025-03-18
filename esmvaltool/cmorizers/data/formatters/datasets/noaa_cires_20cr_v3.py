"""ESMValTool CMORizer for NOAA-CIRES-20CR-V3 data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://psl.noaa.gov/data/gridded/data.ncep.reanalysis2.html

Last access
    20230327

Download and processing instructions
    To facilitate the download, the links to the ftp server are provided.

    ftp://ftp.cdc.noaa.gov/Datasets/20thC_ReanV3/Monthlies/

        pr_wtr.eatm.mon.mean.nc
        cldwtr.eatm.mon.mean.nc
        tcdc.eatm.mon.mean.nc
        ulwrf.ntat.mon.mean.nc
        uswrf.ntat.mon.mean.nc
        csulf.ntat.mon.mean.nc
        csusf.ntat.mon.mean.nc
        shum.mon.mean.nc

Caveats

"""

from .ncep_ncar_r1 import cmorization

# The following line makes it clear that the above import is not an error
__all__ = ["cmorization"]
