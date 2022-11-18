"""ESMValTool CMORizer for NCEP-DOE-R2 data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://psl.noaa.gov/data/gridded/data.ncep.reanalysis2.html

Last access
    20220906

Download and processing instructions
    To facilitate the download, the links to the ftp server are provided.
    Since the filenames are sometimes identical across different
    save the data in two subdirectories in input_dir_path.
    Subdirectory pressure/:
      ftp://ftp.cdc.noaa.gov/Projects/Datasets/ncep.reanalysis2/Monthlies/pressure/
        rhum.mon.mean.nc
        air.mon.mean.nc
        
    Subdirectory surface/:
      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis2/Monthlies/gaussian_grid
        tcdc.eatm.mon.mean.nc
      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis2/Monthlies/surface
        pr_wtr.eatm.mon.mean.nc

    #Select the section "Pressure" and "Surface" and download the variables
    #listed below. Since raw data on pressure levels and for surface have the
    #same file and variable name, save the data in two different subdirectories
    #"press" and "surf" in input_dir_path.

Caveats

"""
from .ncep_ncar_r1 import cmorization

# The following line makes it clear that the above import is not an error
cmorization