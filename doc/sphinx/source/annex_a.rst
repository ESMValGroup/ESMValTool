More tables
***********

:numref:`tab_direc_struc` Directory structure of the ESMValTool sorted by file type.

.. tabularcolumns:: |p{5.2cm}|p{10.3cm}|

.. _tab_direc_struc:

+-------------------------------------------------------------------------------+
| *Namelists*                                                                   |
+-------------------------------+-----------------------------------------------+
| nml/namelist_XyZ.xml          | Namelists for specifying general parameters,  |
|                               | input data and diagnostics to run.            |
+-------------------------------+-----------------------------------------------+

.. tabularcolumns:: |p{5.2cm}|p{10.3cm}|

+-------------------------------------------------------------------------------+
| *Configuration files*                                                         |
+-------------------------------+-----------------------------------------------+
| nml/cfg_XyZ/cfg_XyZ_*.typ     | Configuration files for diagnostic scripts.   |
|                               | The suffix ".typ" specifies the language the  |
|                               | routine is written in. Note: there is usually |
|                               | than one configuration script per diagnostic  |
|                               | set.                                          |
+-------------------------------+-----------------------------------------------+

.. tabularcolumns:: |p{5.2cm}|p{10.3cm}|

+-------------------------------------------------------------------------------+
| *Scripts*                                                                     |
+-------------------------------+-----------------------------------------------+
| main.py                       | Driver script controlling the overall program |
|                               | flow                                          |
+-------------------------------+-----------------------------------------------+
| diag_scripts/                 | Directory containing all diagnostics called   |
|                               | by the namelists. Supporting routines are     |
| - MyDiag.ncl                  | placed in "diag_scripts/lib" under the        |
| - SeaIce_polcon.ncl           | subdirectory corresponding to the programming |
| - SAMonsoon.ncl               | language used (NCL, Python, R).               |
| - etc.                        |                                               |
+-------------------------------+-----------------------------------------------+
| aux/                          | Functions and procedures specific to          |
|                               | a given diagnostic are stored in the          |
| - <diagnostic>/               | subdirectory                                  |
|                               | diag_scripts/aux/<diagnostic>.                |
+-------------------------------+-----------------------------------------------+
| lib/                          | Functions that are called by the              |
|                               | diag_scripts, for example statistics.typ      |
| - ncl/                        | collects all statistical functions in a       |
|                               | single file. When adding a new function,      |
|   - ensemble.ncl              | it must be added to list in the header.       |
|   - regrid.ncl                |                                               |
|   - statistics.ncl            |                                               |
|   - style.ncl                 |                                               |
|   - ...                       |                                               |
|                               |                                               |
| - python/                     |                                               |
|                               |                                               |
|   - ensemble.py               |                                               |
|   - regrid.py                 |                                               |
|   - statistics.py             |                                               |
|   - style.py                  |                                               |
|   - ...                       |                                               |
|                               |                                               |
| - ... for other languages     |                                               |
|   (e.g., R)                   |                                               |
+-------------------------------+-----------------------------------------------+

.. tabularcolumns:: |p{5.2cm}|p{10.3cm}|

+-------------------------------+-----------------------------------------------+
| plot_scripts/                 | Plotting routines; files should have an       |
|                               | intuitive name for their purpose. Data to be  |
| - ncl/                        | plotted may be passed to them directly or via |
|                               | netCDF files.                                 |
|   - plotting1.ncl             |                                               |
|   - plotting2.ncl             |                                               |
|   - ...                       |                                               |
|                               |                                               |
| - python/                     |                                               |
|                               |                                               |
|   - plotting1.py              |                                               |
|   - plotting2.py              |                                               |
|                               |                                               |
| - ... for other languages     |                                               |
|   (e.g., R)                   |                                               |
+-------------------------------+-----------------------------------------------+
| interface_data/               | Inter-process communication, e.g., between    |
|                               | Python and NCL/R, is done by sourcing NCL/R   |
|                               | specific files updated on the fly in this     |
|                               | folder. These intermediate files are based on |
|                               | templates files for different languages.      |
+-------------------------------+-----------------------------------------------+
| interface_scripts/            | Routines called from the workflow manager     |
|                               | script "main.py", mainly used to handle the   |
|                               | control flow of the tool, e.g., parsing       |
|                               | namelists, updating temporary files in the    |
|                               | folder interface_data/, etc).                 |
+-------------------------------+-----------------------------------------------+

.. tabularcolumns:: |p{5.2cm}|p{10.3cm}|

+-------------------------------+-----------------------------------------------+
| reformat_scripts/             | Routines for checking or reformatting raw     |
|                               | input data.                                   |
+-------------------------------+-----------------------------------------------+
| reformat_scripts/cmor/        | Contains the CMOR tables, defined as          |
|                               | plain-text files CMOR_<var>.dat, where <var>  |
| - CMOR_<var>.dat              | is the variable's standard name. Additional   |
|                               | tables can be added by the users, e.g.,       |
|                               | from http://www2-pcmdi.llnl.gov/cmor/tables/. |
+-------------------------------+-----------------------------------------------+
| reformat_scripts/default/     | Contains the reformat routines, in two NCL    |
|                               | scripts.                                      |
| - reformat_default_main.ncl   |                                               |
| - reformat_default_func.ncl   |                                               |
+-------------------------------+-----------------------------------------------+
| reformat_scripts/ECEARTH/     | Contains the EC-Earth/NEMO-specific reformat  |
|                               | routines.                                     |
| - reformat_ECEARTH_main.ncl   |                                               |
| - reformat_ECEARTH_func.ncl   |                                               |
| - names_ECEARTH.dat           |                                               |
| - make_lsm3d.sc               |                                               |
+-------------------------------+-----------------------------------------------+
| reformat_scripts/EMAC         | Contains the EMAC-specific reformat routines. |
|                               |                                               |
| - reformat_EMAC_main.ncl      |                                               |
| - reformat_EMAC_func.ncl      |                                               |
| - names_EMAC.dat              |                                               |
+-------------------------------+-----------------------------------------------+
| reformat_scripts/fixes/       | Contains the user-defined, project- and       |
|                               | model-specific fixes, defined as NCL          |
| - <project>_<model>_fixes.ncl | scripts <project>_<model>_fixes.ncl.          |
|                               | A template is also provided for the user      |
|                               | to add more fixes.                            |
+-------------------------------+-----------------------------------------------+
| reformat_scripts/obs/         | Contains specific reformat routines for       |
|                               | "cmorizing" observational data.               |
+-------------------------------+-----------------------------------------------+
| reformat_scripts/constants.ncl| Contains general-purpose functions and        |
|                               | procedure, called by the default, the         |
|                               | ECEARTH- and the EMAC-specific routines.      |
+-------------------------------+-----------------------------------------------+
| reformat_scripts/recognized\_ | Provides a list of possible alternative units |
| units.dat                     | to the CMOR standard and the corresponding    |
|                               | conversion factor. Can be extended by the     |
|                               | user.                                         |
+-------------------------------+-----------------------------------------------+
| reformat_scripts/recognized\_ | Provides a list of possible alternative       |
| vars.dat                      | variable names to the CMOR standard           |
|                               | names. Can be extended by the user.           |
+-------------------------------+-----------------------------------------------+
| reformat_scripts/variable\_   | Declaration of variables, variable specific   |
| defs/                         | attributes and calculation of derived         |
|                               | variables                                     |
+-------------------------------+-----------------------------------------------+

.. tabularcolumns:: |p{5.2cm}|p{10.3cm}|

+-------------------------------------------------------------------------------+
| *Data folders*                                                                |
+-------------------------------+-----------------------------------------------+
|                               | The data folders are specified in             |
|                               | nml/namelist_*, and thus may be different     |
|                               | from the defaults given here. These folders   |
|                               | contain the output generated by the ESMValTool|
|                               | and are created on the fly if needed. Note    |
|                               | that these folders do not need to be in       |
|                               | the same directory as the source code. They   |
|                               | can be arbitrarily specified  in the namelist |
|                               | as path relative to the root path. Using      |
|                               | symbolic links is another option  to separate |
|                               | the actual data from the code.                |
+-------------------------------+-----------------------------------------------+
| climo/                        | Quality checked and derived netCDF files,     |
|                               | reformatted from the original data.           |
+-------------------------------+-----------------------------------------------+
| plots/                        | Destination directory for graphics files.     |
+-------------------------------+-----------------------------------------------+
| work/                         | Miscellaneous files produced during run-time, |
|                               | e.g., optional netCDF output and              |
|                               | references/acknowledgements.                  |
+-------------------------------+-----------------------------------------------+


.. _workflow_reformat:

Workflow of reformat routines
*****************************

**Control flow of reformat_default**

The reformat_default_main.ncl script sets the global variables as defined in reformat.py (input and output paths, variable name and field, model name and ensemble, etc.) and then performs a list of operations calling various functions and procedures defined in reformat_default_func.ncl. The workflow is as follows:

* find grid type: the data can be defined on a standard rectilinear grid or on an irregular grid. In the latter case, the script does not modify the grid properties and additionally attaches the area field (the area weights) for the irregular grid to the output file. The location of the area file is typically defined as an entry in the namelist, for example by using the project class CMIP5_gridfile where the final entry is the full path to the area file, see :numref:`tab_proj_spec`.
* read variable: the selected variable is read from the input file. If the variable is not found, the reading function checks for possible alternative variable names (as specified in recognized_vars.dat), before issuing an error message.
* apply project- and model-specific fixes: if a fixing procedure is found in the fixes/ directory for the selected project and model, it is called at this point in order to apply the user-defined corrections to the data.
* create time-series: the variable is read for the selected time range (start_year-end_year) and a time-series is created.
* rank/field consistency: the consistency of variable's rank with the given field (T3M, T2Mz, T2Ms, etc.). A simple calculation of the zonal mean is performed in case a rank 4 variable is provided with T2?z field.
* check fill values: a default missing values is assigned if the variable does not have one. The function then looks for data values that might represent undefined missing values. Currently it considers: -999., -9999., -99999., -999999., 1.e15, -1.e34. Finally, the ESMValTool default missing values (1.e20) is assigned as a standard _FillValue to the variable.
* reformat time coordinate: the time coordinate is reformatted according to the CMOR standard. If a calendar attribute is not assigned, the standard is used. The consistency of the time-series with the selected time range is checked.
* reformat vertical coordinate (applies only to certain fields and to rectilinear grids): the vertical coordinate is assigned "Pa" units, converting from the most common alternative units (mbar, bar, hPa) if required. The ordering is set from top to bottom (monotonically decreasing).
* reformat latitude coordinate (applies only to certain fields and to rectilinear grids): the ordering is set from South to North (monotonically increasing).
* reformat longitude coordinate (applies only to certain fields and to rectilinear grids): the ordering is set from 0 to 360 degrees.
* check units: consistency of the variable's units with the CMOR standard is checked. The CMOR table for the selected variable must be available in the CMOR/ directory (an error message is issued otherwise). Units renaming and conversion can also be performed, if the corresponding information is given in recognized_units.dat.
* set variable attributes: the CMOR standard attributes are assigned to the selected variable. The corresponding CMOR table must be available in the CMOR/ directory (an error message is issued otherwise).
* write output file: the variable reformatted according to the CMOR standard is written in the selected output file.
* add info for irregular grids (applies only to irregular grids): the area file of the irregular grid is added, this file may later be used for averaging. 


**Control flow of reformat_ECEARTH**

This reformat procedure can be used to convert raw EC-Earth/NEMO files to a format that complies with the ESMValTool requirements. It performs the following additional operations compared with the default workflow:

1. find_name: translate the EC-Earth/NEMO name to a CMOR name using the table names_ECEARTH.dat.
2. sub_staggergrid: determine grid type (T, U, V) and add that information to the filename.
3. mask_land: land points have the value 0 in the raw files, not a fill value (missing value). This routine sets land points (as in the landmask file) to fill values.
4. rename_time: rename time variable from EC-Earth name to standard name and remove the attribute _FillValue.
5. rename_lev: vertical coordinate name in raw files depends on grid, rename it to lev. Requires the external input table names_ECEARTH.dat.
6. add_ijcoordvars: add i and j variables and assign them as coordinate variables.
7. convert_units: unit conversions that cannot be handled by check_units.
8. add_ECEARTHgrid_info: add ECEARTH grid info (lat, lon, areacello and grid sizes) to the output. 

**Control flow of reformat_EMAC**

The workflow is similar to the default case, but some additional operations specific to the EMAC model are performed in addition:

9. find messy version: the MESSy version with which the EMAC output has been generated is read from the data.
10. find EMAC name: the EMAC name of the selected variable is found from the table in names_EMAC.dat (an error message is issued if not defined). For complex variables (i.e., variables not directly available as EMAC output but derivable from other EMAC variables), a user-defined recipe can be provided in reformat_scripts/EMAC/recipes/EMAC_recipes_xxx.ncl to derive it.
11. check field consistency: reads from names_EMAC.dat file the list of allowed fields for the selected variable (for example is not possible to select total column ozone toz as T3M field).
12. check vertical integration type (only for T2?s types): reads from names_EMAC.dat the option for the vertical coordinate (C for column integration, S for surface value).
13. start the time loop: the EMAC output is assumed to be monthly-aggregated (monthly averages are optional). The data are read starting from January of the start_year to December of the end_year.
14. extract variable: the selected variable is searched in the EMAC output. If multiple files for a given month/year combination contain the selected variable, the following priority list is applied: time coordinate matching the field type (monthly mean or daily output), data from tracer_gp and tr_* streams/channels, first file in the list. For complex variables, the corresponding user-defined recipe is called (reformat_scripts/EMAC/recipes/EMAC_recipes_xxx.ncl). For T2?z types, the data are interpolated on constant pressure levels (defined in reformat_scripts/constants.ncl).
15. create time series: within the time loop, a time series start_year-end_year is created.
16. reformat coordinates, check units, set variable attributes and write output: these operations are applied exactly as in the default case.

The user can extend the reformat_scripts/EMAC/recipes/EMAC_recipes_xxx.ncl in order to calculate additional (derived) 
variables not directly available in EMAC.

