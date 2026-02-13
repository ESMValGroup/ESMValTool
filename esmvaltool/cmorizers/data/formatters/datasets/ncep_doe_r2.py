"""ESMValTool CMORizer for NCEP-DOE-R2 data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://psl.noaa.gov/data/gridded/data.ncep.reanalysis2.html

Last access
    20220906

Download and processing instructions
    To facilitate the download, the links to the https server are provided.

      https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/Monthlies/
      pressure/
        rhum.mon.mean.nc
        air.mon.mean.nc
      https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/Monthlies/
      gaussian_grid/
        tcdc.eatm.mon.mean.nc
      https://downlooads.psl.noaa.gov/Datasets/ncep.reanalysis2/Monthlies/
      surface/
        pr_wtr.eatm.mon.mean.nc

Caveats

"""

from .ncep_ncar_r1 import cmorization

# The following line makes it clear that the above import is not an error
__all__ = ["cmorization"]
