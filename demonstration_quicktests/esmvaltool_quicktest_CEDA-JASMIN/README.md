EXPLAINING THE QUICK TEST
==========================
This file will tell you how to run a very quick test of ESMValTool on the CEDA-Jasmin cluster;
you don't need to install ESMValTool nor do you need to get any data, all these are done for you
automatically. Just get the

```
run_esmvaltool_test.sh
```

runner and source it from wherver you are on Jasmin. Follow the on-screen messages and when it's done,
come back here and read the explanations below. You should have a fair understanding of ESMValTool after this.
When you feel brave enough check out the documentation pages at:

http://esmvaltool.readthedocs.io/en/refactoring_backend/

Have fun!

1. Explaining config and namelist files
========================================

Get into the directory that you ran esmvaltool and have a first look:
(in my case the dates are different than yours, esmvaltool creates the subdirs based on timestamps so 
you can run multiple instances at different times):

(after installing like a big boy, conda envirinment will be called esmvaltool)

```
(esmvaltool) [valeriu@jasmin-sci1 valeriu]$ cd ESMValTool_QuickTest/
(esmvaltool) [valeriu@jasmin-sci1 ESMValTool_QuickTest]$ ls -la
total 320
drwxr-sr-x  3 valeriu gws_ncas_cms  4096 May  1 12:16 .
drwxr-sr-x 78 valeriu gws_ncas_cms 16384 May  1 12:16 ..
drwxr-sr-x  3 valeriu gws_ncas_cms  4096 May  1 12:16 PREPROCESSOR
-rw-r--r--  1 valeriu gws_ncas_cms  1109 May  1 12:16 config.yml
-rw-r--r--  1 valeriu gws_ncas_cms  1910 May  1 12:16 namelist.yml
```

 - config.yml contains the user specific run settings:

```
###############################################################################
# User's configuration file for the ESMValTool
###############################################################################
---

# Diagnostics create plots? [true]/false
write_plots: true
# Diagnositcs write NetCDF files? [true]/false
write_netcdf: true
# Set the console log level debug, [info], warning, error
log_level: info
# verbosity is deprecated and will be removed in the future
# verbosity: 1
# Exit on warning? true/[false]
exit_on_warning: false
# Plot file format? [ps]/pdf/png/eps/epsi
output_file_type: pdf
# Destination directory
output_dir: ./PREPROCESSOR
# Save intermediary cubes in the preprocessor true/[false]
save_intermediary_cubes: true
#max_parallel_tasks: 4

rootpath:
  # Rootpath to CMIP5 data
  CMIP5: /home/users/valeriu/ESMValTool_KIT
  # Rootpath to OBS data
  OBS: /home/users/valeriu/ESMValTool_KIT
  # Default
  default: /home/users/valeriu/ESMValTool_KIT

# Directory structure for input data: [default]/BADC/DKRZ/ETHZ/etc
# See config-developer.yml for definitions.
drs:
  CMIP5: default
```

The namelist.yml file contains directives and parameters that tell the code which simulation models to use, which observation data to use, what are the preprocessing steps and their respective workflow parameters (see below for a better explanation of the preprocessors) and ultimately what diagnostics to run and what are the specific (external) scripts for diagnostics:

```
models:
  - {model: bcc-csm1-1,       project: CMIP5, exp: historical, ensemble: r1i1p1, start_year: 2000, end_year: 2002}
  - {model: GFDL-ESM2G,       project: CMIP5, exp: historical, ensemble: r1i1p1, start_year: 2000, end_year: 2002}
  - {model: MPI-ESM-LR,       project: CMIP5, exp: historical, ensemble: r1i1p1, start_year: 2000, end_year: 2002}
  - {model: MPI-ESM-MR,       project: CMIP5, exp: historical, ensemble: r1i1p1, start_year: 2000, end_year: 2002}

preprocessors:
  pp850:
    extract_levels:
      levels: 85000
      scheme: linear
    regrid:
      target_grid: ERA-Interim
      scheme: linear
    mask_fillvalues:
      threshold_fraction: 0.95
    multi_model_statistics:
      span: overlap
      statistics: [mean, median]
      #exclude: [GFDL-ESM2G, reference_model, NCEP]

diagnostics:
  ta850:
    description: Air temperature at 850 hPa global.
    variables:
      ta:
        preprocessor: pp850
        reference_model: ERA-Interim
        alternative_model: NCEP
        mip: Amon
        field: T3M
    additional_models:
      - {model: ERA-Interim,  project: OBS,  type: reanaly,  version: 1,  start_year: 2000,  end_year: 2002,  tier: 3}
      - {model: NCEP,         project: OBS,  type: reanaly,  version: 1,  start_year: 2000,  end_year: 2002,  tier: 2}
    scripts:
      cycle:
        script: perfmetrics/main.ncl
        plot_type: cycle       # Plot type ('cycle' [time], 'zonal' [plev, lat], 'latlon' [lat, lon], 'cycle_latlon' [time, lat, lon])
        time_avg: monthlyclim  # Time average ('opt' argument of time_operations.ncl)
        region: Global         # Selected region ('Global', 'Tropics', 'NH extratropics', 'SH extratropics')
        plot_stddev: ref_model # Plot standard deviation ('all', 'none', 'ref_model' or given model name)
        legend_outside: true   # Plot legend in a separate file
```

Let's have a look at the output now: PREPROCESSOR contains the output from the run:

```
[valeriu@jasmin-sci2 ESMValTool_QuickTest]$ ls -la PREPROCESSOR/namelist_20180501_111616/
total 384
drwxr-sr-x 6 valeriu gws_ncas_cms 4096 May  1 12:18 .
drwxr-sr-x 3 valeriu gws_ncas_cms 4096 May  1 12:16 ..
drwxr-sr-x 3 valeriu gws_ncas_cms 4096 May  1 12:18 plots -- where the plots live
drwxr-sr-x 3 valeriu gws_ncas_cms 4096 May  1 12:16 preproc -- where all the preproc files live
drwxr-sr-x 3 valeriu gws_ncas_cms 4096 May  1 12:18 run -- where the runfiles live
drwxr-sr-x 3 valeriu gws_ncas_cms 4096 May  1 12:18 work -- where the diagnostic specific files live
```

Let's have a walk around these subdirectories! We start with /preproc:

2. The ```/preproc``` directory and "preprocessing" concept
======================================================

The preproc directory contains the files that result from running the Preprocessor module. To better understand the contents of
this directory and the meaning the files have for you, we need to look at the namelist.yml file, and specifically, at the "preprocessors"
section of this file:

This section looks like this (in this example case):

```
preprocessors:
  pp850:
    extract_levels:
      levels: 85000
      scheme: linear
    regrid:
      target_grid: ERA-Interim
      scheme: linear
    mask_fillvalues:
      threshold_fraction: 0.95
    multi_model_statistics:
      span: overlap
      statistics: [mean, median]
      #exclude: [GFDL-ESM2G, reference_model, NCEP]
```

Things to note and be aware of:

 - the section is actually called "preprocessors" since a namelist may contain more than one preprocessor section; they are listed in order
and labelled by the sub-section name, in this case "pp850" (you will understand why pp850 when we grt to the vertical level selection "extract_levels"); these
preprocessor instances may bear any name, but usually one tens to call them by a significant label like "pp" (preprocessor) + whatever it is specific to the preprocessor, like here...
 - ...extract_level: and that means we want to extract 850hPa vertical level pr pressure level (hence pp850) from all the models and observational data; the scheme
specifies the type of interpolation that will be performed when level extracting (for the cases where the exact 850hPa level is not present in the data, but rather two neighnoring levels e.g. 800hPa and 900hPa);
 - regrid: sets the parameters for horizontal regridding of each of the models and obs data; in this case we are regridding on the native grid of the ERA-Interim obs data; again, "scheme" specifies the type of interpolation used for regridding;
 - mask_fillvalues: this module creates a common mask for all models and obs data that combines all the missing values masks from each model and obs data; the masking accepts a statistical parameter ("threshold_fraction") which sets the threshold of accepting/rejecting a certain grid cell based on the number of missing values in the data (here, if 95% or less of values are missing then we reject that particular grid cell);
 - multi_model_statistics combines all the already vertical level selected, horizontally regridded and masked data sets and computes a number of statistics (here mean and median) from all these models and obs data; the "overlap" parameter tells the code that these statistics will be computed only on the overlap in time from all the used models and obs data sets; these are written to netCDF files and can be further treated as new models for further operations (not the case here); you can exclude certain combinations of models or models and obs sets from the multi model computations by setting "exclude";

These preprocessor functions and their working parameters are set in the namelist; there are other preprocessing functionalities that could be called in the "preprocessors/name of preprocessor" sections e.g. area select, land/ocean masks etc. It is to be noted that these data processing actions that comprise the preprocessor are run BEFORE the actual diagnostic is run, and are direct requirements of the actual diagnostic. These steps are a result of functional standardization of common data processing steps that were previously present in a large number of diagnostics; taking them out of the diagnostic and standardizing them via iris cubes manipulation brings forth both a more centralized and more flexible approach to data and significant increases in performance. They are, however, not the only steps that are applied to all models and obs data - a nice way of observing the workflow of the generalized preprocessor engine is to inspect the data directories created in /preproc/name_of_diagnostic/any_given_model_and_its_parameters/ e.g.:

```
[valeriu@jasmin-sci2 ESMValTool_QuickTest]$ ls -la PREPROCESSOR/namelist_20180501_111616/preproc/ta850_pp850_ta/CMIP5_GFDL-ESM2G_Amon_historical_r1i1p1_T3M_ta_2000-2002
total 613568
drwxr-sr-x 2 valeriu gws_ncas_cms      4096 May  1 12:18 .
drwxr-sr-x 8 valeriu gws_ncas_cms      4096 May  1 12:18 ..
-rw-r--r-- 1 valeriu gws_ncas_cms 105797202 May  1 12:17 00_load_cubes.nc
-rw-r--r-- 1 valeriu gws_ncas_cms 105802154 May  1 12:17 01_fix_metadata.nc
-rw-r--r-- 1 valeriu gws_ncas_cms 105799052 May  1 12:17 02_concatenate.nc
-rw-r--r-- 1 valeriu gws_ncas_cms 105799052 May  1 12:17 03_cmor_check_metadata.nc
-rw-r--r-- 1 valeriu gws_ncas_cms  31765371 May  1 12:17 04_extract_time.nc
-rw-r--r-- 1 valeriu gws_ncas_cms  31765371 May  1 12:17 05_fix_data.nc
-rw-r--r-- 1 valeriu gws_ncas_cms   1903877 May  1 12:17 06_extract_levels.nc
-rw-r--r-- 1 valeriu gws_ncas_cms  16688990 May  1 12:17 07_regrid.nc
-rw-r--r-- 1 valeriu gws_ncas_cms  16688990 May  1 12:18 08_mask_fillvalues.nc
-rw-r--r-- 1 valeriu gws_ncas_cms  16688990 May  1 12:18 09_multi_model_statistics.nc
-rw-r--r-- 1 valeriu gws_ncas_cms  16688990 May  1 12:18 10_cmor_check_data.nc
```

These types of directories are created automatically when setting "save_intermediary_cubes" to "true" in the config.yml file; they reflect the order in which each input data file gets preprocessed so to prep it for the diagnostic run (and use the ordered prefix to show you whci files get created when), so in our case:

- ```00_load_cubes.nc``` -- contains the loaded data from netCDF files in rootpath: in config.yml (raw loading, no checks or fixes applied);
- ```01_fix_metadata.nc``` -- fixes to metadata are applied;
- ```02_concatenate.nc``` -- stitch together multiple files from a given model to cover the required time span;
- ```03_cmor_check_metadata.nc``` -- specific CMOR checks are applied;
- ```04_extract_time.nc``` -- time gate data so you select and keep only the data for the required time span;
- ```05_fix_data.nc``` -- specific data fixes;
- ```06_extract_levels.nc``` -- vertical level selection, discard all other unneeded levels;
- ```07_regrid.nc``` -- horizontal regridding on the vertical level extracted data;
- ```08_mask_fillvalues.nc``` -- apply the missing values mask on the regridded data;
- ```09_multi_model_statistics.nc``` -- compute the multi model statistics on the regridded data;
- ```10_cmor_check_data.nc``` -- further final CMOR checks;

These files are not usually written to disk unless the user really wants to (usually for debug purposes); instead they are passed from one module to another in the form of iris cubes in memory; not writing to disk is much more efficient in terms of run time and obviously, in terms of disk space; for an example number 5 CMIP models and two OBS models, writing to file the intermediary files may increase the preprocessing time up to five times.

Once all these preprocessing steps are completed succesfully, the diagnostic can now start running and will be using the preprocessed files. These preprocessed files, irrespective of the writing of intermediary files or not, are located in the /preproc/name_of_diagnostic/ directory e.g.:

```
[valeriu@jasmin-sci2 ESMValTool_QuickTest]$ ls -la PREPROCESSOR/namelist_20180501_111616/preproc/ta850_pp850_ta/*nc*
-rw-r--r-- 1 valeriu gws_ncas_cms 16688888 May  1 12:18 PREPROCESSOR/namelist_20180501_111616/preproc/ta850_pp850_ta/CMIP5_GFDL-ESM2G_Amon_historical_r1i1p1_T3M_ta_2000-2002.nc
-rw-r--r-- 1 valeriu gws_ncas_cms 16685281 May  1 12:18 PREPROCESSOR/namelist_20180501_111616/preproc/ta850_pp850_ta/CMIP5_MPI-ESM-LR_Amon_historical_r1i1p1_T3M_ta_2000-2002.nc
-rw-r--r-- 1 valeriu gws_ncas_cms 16685281 May  1 12:18 PREPROCESSOR/namelist_20180501_111616/preproc/ta850_pp850_ta/CMIP5_MPI-ESM-MR_Amon_historical_r1i1p1_T3M_ta_2000-2002.nc
-rw-r--r-- 1 valeriu gws_ncas_cms 16685275 May  1 12:18 PREPROCESSOR/namelist_20180501_111616/preproc/ta850_pp850_ta/CMIP5_bcc-csm1-1_Amon_historical_r1i1p1_T3M_ta_2000-2002.nc
-rw-r--r-- 1 valeriu gws_ncas_cms 16677926 May  1 12:18 PREPROCESSOR/namelist_20180501_111616/preproc/ta850_pp850_ta/MultiModelMean_T3M_ta_2000-2002.nc
-rw-r--r-- 1 valeriu gws_ncas_cms 16677926 May  1 12:18 PREPROCESSOR/namelist_20180501_111616/preproc/ta850_pp850_ta/MultiModelMedian_T3M_ta_2000-2002.nc
-rw-r--r-- 1 valeriu gws_ncas_cms 16683489 May  1 12:18 PREPROCESSOR/namelist_20180501_111616/preproc/ta850_pp850_ta/OBS_ERA-Interim_reanaly_1_T3M_ta_2000-2002.nc
-rw-r--r-- 1 valeriu gws_ncas_cms 16683517 May  1 12:18 PREPROCESSOR/namelist_20180501_111616/preproc/ta850_pp850_ta/OBS_NCEP_reanaly_1_T3M_ta_2000-2002.nc
-rw-r--r-- 1 valeriu gws_ncas_cms     6611 May  1 12:18 PREPROCESSOR/namelist_20180501_111616/preproc/ta850_pp850_ta/ta_info.ncl
```

Note that ALL files are there, both from CMIP simulations and OBS data; each of these will have contained the data processed according to the preprocessor steps and after all its steps have been applied ie if the last preprocessor setting was mask_fillvalues, the file e.g. ```PREPROCESSOR/namelist_20180501_111616/preproc/ta850_pp850_ta/CMIP5_MPI-ESM-MR_Amon_historical_r1i1p1_T3M_ta_2000-2002.nc``` contains ```MPI-ESM-MR``` data that has been raw loaded, metadata fixed, concatenated for time completion, time gated between 2000 and 2002 (inclusive of these years, so 3 years worth of data), ..., masked for missing values of variable ta, atmospheric montly MIP, historical experiment. This file together with the others obtained in the same manner is now ready to be passed to the diagnostic:

3. Diagnostic and its settings
===============================

You can see the diagnostic settings in the namelist.yml:

```
diagnostics:
  ta850:
    description: Air temperature at 850 hPa global.
    variables:
      ta:
        preprocessor: pp850
        reference_model: ERA-Interim
        alternative_model: NCEP
        mip: Amon
        field: T3M
    additional_models:
      - {model: ERA-Interim,  project: OBS,  type: reanaly,  version: 1,  start_year: 2000,  end_year: 2002,  tier: 3}
      - {model: NCEP,         project: OBS,  type: reanaly,  version: 1,  start_year: 2000,  end_year: 2002,  tier: 2}
    scripts:
      cycle:
        script: perfmetrics/main.ncl
        plot_type: cycle       # Plot type ('cycle' [time], 'zonal' [plev, lat], 'latlon' [lat, lon], 'cycle_latlon' [time, lat, lon])
        time_avg: monthlyclim  # Time average ('opt' argument of time_operations.ncl)
        region: Global         # Selected region ('Global', 'Tropics', 'NH extratropics', 'SH extratropics')
        plot_stddev: ref_model # Plot standard deviation ('all', 'none', 'ref_model' or given model name)
        legend_outside: true   # Plot legend in a separate file
        styleset: CMIP5        # Plot style
```

"diagnostics" may include a large number of diagnostics, each of them requiring data preprocessing identified by "preprocessor:"; in our case we run a single diagnostic called "ta850" (variable is Air Temperature, ta, and we need just the 850hPa vertical level). Each diagnostic may be run on a multitude of variables: "variables:" contains the name of the variables and the preprocessor settings and the reference and alternative models. These models are used for both the diagnostic per se (e.g. when plotting, comparisons are made with the reference models) and for preprocessing via the "additional_models" list of dictionaries - here you see that ERA-Interim is used for regridding all the simulation CMIP models onto it.

"scripts" contains the names and directives and parameters for the diagnostic scripts: here we are running perfmetrics/main.ncl NCL diagnostic script, that wants to know what plot type, what temporal averaging, region etc. are needed.

4. The diagnostic output
=========================

The single plot this diagnostic outputs is in: ```PREPROCESSOR/namelist_20180501_111616/plots/ta850/cycle/``` and it's really basic -- the montly global mean of temperatures with OBS models as reference. Quite a bit of hoolabaloo for just such a plot :)
