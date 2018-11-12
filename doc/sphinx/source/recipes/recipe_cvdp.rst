CVDP Package
============

About
-----

The NCAR Climate Variability Package (CVDP) is a stand alone ncl based application for the analysis of climate variability in models and observations, see [1].

Requirements
------------

+ cvdp (https://github.com/NCAR/CVDP-ncl)
+ nco (optional for creating netcdf files)

Configuration
-------------

To use it within the ESMValTool set the environment variable *CVDP_ROOT* to the top level directory of the CVDP package.

Also adapt the *config-developer.yml* file at the keyword: *output_file* to let the preprosessor output end in
*YYYYMM-YYYYMM.nc*
e.g.:
    -  output_file: '[project]_[dataset]_[mip]_[exp]_[ensemble]_[field]_[short_name]_[start_year]-[end_year]'
    +  output_file: '[project]_[dataset]_[mip]_[exp]_[ensemble]_[field]_[short_name]_[start_year]01-[end_year]12'


References
----------
[1] http://www.cesm.ucar.edu/working_groups/CVC/cvdp/

