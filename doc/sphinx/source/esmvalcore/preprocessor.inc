.. _preprocessor:

************
Preprocessor
************
The ESMValTool preprocessor can be used to perform all types of climate data pre-processing needed before indices or diagnostics can be calculated. It is a base component for many other diagnostics and metrics shown on this portal. It can be applied to tailor the climate model data to the need of the user for its own calculations.

Features of the ESMValTool Climate data pre-processor are:

* Regridding
* Geographical area selection
* Aggregation of data
* Provenance tracking of the calculations
* Model statistics
* Multi-model mean
* and many more

Variable derivation
===================
Documentation of _derive.py


Time manipulation
=================
The _time.py module contains the following preprocessor functions:

* extract_time: Extract a time range from a cube.
* extract_season: Extract only the times that occur within a specific season.
* extract_month: Extract only the times that occur within a specific month.
* time_average: Take the weighted average over the time dimension.
* seasonal_mean: Produces a mean for each season (DJF, MAM, JJA, SON)
* annual_mean: Produces an annual or decadal mean.
* regrid_time: Aligns the time axis of each dataset to have common time points and calendars.

1. extract_time
---------------

This function subsets a dataset between two points in times. It removes all
times in the dataset before the first time and after the last time point.
The required arguments are relatively self explanatory:

* start_year
* start_month
* start_day
* end_year
* end_month
* end_day

These start and end points are set using the datasets native calendar.
All six arguments should be given as integers - the named month string
will not be accepted.

See also :func:`esmvalcore.preprocessor.extract_time`.


2. extract_season
-----------------

Extract only the times that occur within a specific season.

This function only has one argument: `season`. This is the named season to
extract. ie: DJF, MAM, JJA, SON.

Note that this function does not change the time resolution. If your original
data is in monthly time resolution, then this function will return three
monthly datapoints per year.

If you want the seasonal average, then this function needs to be combined with
the seasonal_mean function, below.

See also :func:`esmvalcore.preprocessor.extract_season`.


3. extract_month
----------------

The function extracts the times that occur within a specific month.
This function only has one argument: `month`. This value should be an integer
between 1 and 12 as the named month string will not be accepted.

See also :func:`esmvalcore.preprocessor.extract_month`.


4. time_average
---------------

This functions takes the weighted average over the time dimension. This
function requires no arguments and removes the time dimension of the cube.

See also :func:`esmvalcore.preprocessor.time_average`.


5. seasonal_mean
----------------

This function produces a seasonal mean for each season (DJF, MAM, JJA, SON).
Note that this function will not check for missing time points. For instance,
if you are looking at the DJF field, but your datasets starts on January 1st,
the first DJF field will only contain data from January and February.

We recommend using the extract_time to start the dataset from the following
December and remove such biased initial datapoints.

See also :func:`esmvalcore.preprocessor.seasonal_mean`.


6. annual_mean
--------------

This function produces an annual or a decadal mean. The only argument is the
decadal boolean switch. When this switch is set to true, this function
will output the decadal averages.

See also :func:`esmvalcore.preprocessor.annual_mean`.


7. regrid_time
--------------

This function aligns the time points of each component dataset so that the dataset
iris cubes can be subtracted. The operation makes the datasets time points common and
sets common calendars; it also resets the time bounds and auxiliary coordinates to
reflect the artifically shifted time points. Current implementation for monthly
and daily data; the frequency is set automatically from the variable CMOR table
unless a custom frequency is set manually by the user in recipe.


Area manipulation
=================
The _area.py module contains the following preprocessor functions:

* extract_region: Extract a region from a cube based on lat/lon corners.
* zonal_means: Calculates the zonal or meridional means.
* area_statistics: Calculates the average value over a region.
* extract_named_regions: Extract a specific region from in the region cooordinate.


1. extract_region
-----------------

This function masks data outside a rectagular region requested. The boundairies
of the region are provided as latitude and longitude coordinates in the
arguments:

* start_longitude
* end_longitude
* start_latitude
* end_latitude

Note that this function can only be used to extract a rectangular region.

See also :func:`esmvalcore.preprocessor.extract_region`.


2. zonal_means
--------------

The function calculates the zonal or meridional means. While this function is
named `zonal_mean`, it can be used to apply several different operations in
an zonal or meridional direction.
This function takes two arguments:

* coordinate: Which direction to apply the operation: latitude or longitude
* mean_type: Which operation to apply: mean, std_dev, variance, median, min or max

See also :func:`esmvalcore.preprocessor.zonal_means`.


3. area_statistics
------------------

This function calculates the average value over a region - weighted by the
cell areas of the region.

This function takes the argument,
operator: the name of the operation to apply.

This function can be used to apply several
different operations in the horizonal plane: mean, standard deviation, median
variance, minimum and maximum.

Note that this function is applied over the entire dataset. If only a specific
region, depth layer or time period is required, then those regions need to be
removed using other preprocessor operations in advance.

See also :func:`esmvalcore.preprocessor.area_statistics`.


4. extract_named_regions
------------------------

This function extract a specific named region from the data. This function
takes the following argument: `regions` which is either a string or a list
of strings of named regions. Note that the dataset must have a `region`
cooordinate which includes a list of strings as values. This function then
matches the named regions against the requested string.

See also :func:`esmvalcore.preprocessor.extract_named_regions`.


Volume manipulation
===================
The _volume.py module contains the following preprocessor functions:

* extract_volume: Extract a specific depth range from a cube.
* volume_statistics: Calculate the volume-weighted average.
* depth_integration: Integrate over the depth dimension.
* extract_transect: Extract data along a line of constant latitude or longitude.
* extract_trajectory: Extract data along a specified trajectory.


1. extract_volume
-----------------

Extract a specific range in the z-direction from a cube.  This function
takes two arguments, a minimum and a maximum (`z_min` and `z_max`,
respectively) in the z direction.

Note that this requires the requested z-coordinate range to be the
same sign as the iris cube. ie, if the cube has z-coordinate as
negative, then z_min and z_max need to be negative numbers.

See also :func:`esmvalcore.preprocessor.extract_volume`.


2. volume_statistics
--------------------

This function calculates the volume-weighted average across three dimensions,
but maintains the time dimension. The following arguments are required:

This function takes the argument: operator, which defines the
operation to apply over the volume.

No depth coordinate is required as this is determined by iris. This
function works best when the fx_files provide the cell volume.

See also :func:`esmvalcore.preprocessor.volume_statistics`.


3. depth_integration
--------------------

This function integrate over the depth dimension. This function does a
weighted sum along the z-coordinate, and removes the z direction of the output
cube. This preprocessor takes no arguments.

See also :func:`esmvalcore.preprocessor.depth_integration`.


4. extract_transect
-------------------

This function extract data along a line of constant latitude or longitude.
This function takes two arguments, although only one is strictly required.
The two arguments are `latitude` and `longitude`. One of these arguments
needs to be set to a float, and the other can then be either ignored or set to
a minimum or maximum value.
Ie: If we set latitude to 0 N and leave longitude blank, it would produce a
cube along the equator. On the other hand, if we set latitude to 0 and then
set longitude to `[40., 100.]` this will produce a transect of the equator
in the indian ocean.

See also :func:`esmvalcore.preprocessor.extract_transect`.


5. extract_trajectory
---------------------

This function extract data along a specified trajectory.
The three areguments are: latitudes and longitudes are the coordinates of the
trajectory.

If two points are provided, the `number_points` argument is used to set a
the number of places to extract between the two end points.

If more than two points are provided, then
extract_trajectory will produce a cube which has extrapolated the data
of the cube to those points, and `number_points` is not needed.

Note that this function uses the expensive interpolate method, but it may be
necceasiry for irregular grids.

See also :func:`esmvalcore.preprocessor.extract_trajectory`.


CMORization and dataset-specific fixes
======================================
Documentation of _reformat.py, check.py and fix.py

Vertical interpolation
======================
Documentation of _regrid.py (part 1)

Masking
=======
Documentation of _mask.py (part 1)

1. Introduction to masking
--------------------------

Certain metrics and diagnostics need to be computed and performed on restricted regions of the Globe; ESMValTool supports subsetting the input data on land mass, oceans and seas, ice. This is achived by masking the model data and keeping only the values associated with grid points that correspond to e.g. land mass
or oceans and seas; masking is done either by using standard mask files that have the same grid resolution as the model data (these files are usually produced
at the same time with the model data and are called fx files) or, in the absence of these files, by using Natural Earth masks. Natural Earth masks, even if they are not model-specific, represent a good approximation since their grid resolution is almost always much higher than the model data, and they are constantly updated with changing
geographical features.

2. Land-sea masking
-------------------

In ESMValTool v2 land-seas-ice masking can be done in two places: in the preprocessor, to apply a mask on the data before any subsequent preprocessing step, and before
running the diagnostic, or in the disgnostic phase. We present both these implementations below.

To mask out seas in the preprocessor step, simply add `mask_landsea:` as a preprocessor step in the `preprocessor` of your choice section of the recipe, example:

.. code-block:: bash

    preprocessors:
      my_masking_preprocessor:
        mask_landsea:
          mask_out: sea

The tool will retrieve the corresponding `fx: stfof` type of mask for each of the used variables and apply the mask so that only the land mass points are
kept in the data after applying the mask; conversely, it will retrieve the `fx: sftlf` files when land needs to be masked out.
`mask_out` accepts: land or sea as values. If the corresponding fx file is not found (some models are missing these
type of files; observational data is missing them altogether), then the tool attempts to mask using Natural Earth mask files (that are vectorized rasters).
Note that the resolutions for the Natural Earth masks are much higher than any usual CMIP model: 10m for land and 50m for ocean masks.

3. Ice masking
--------------

Note that for masking out ice the preprocessor is using a different function, this so that both land and sea or ice can be masked out without
losing generality. To mask ice out one needs to add the preprocessing step much as above:

.. code-block:: bash

    preprocessors:
      my_masking_preprocessor:
        mask_landseaice:
          mask_out: ice

To keep only the ice, one needs to mask out landsea, so use that as value for mask_out. As in the case of mask_landsea, the tool will automatically
retrieve the `fx: sftgif` file corresponding the the used variable and extract the ice mask from it.

4. mask files
-------------

At the core of the land/sea/ice masking in the preprocessor are the mask files (whether it be fx type or Natural Earth type of files); these files (bar Natural Earth)
can be retrived and used in the diagnostic phase as well or solely. By specifying the `fx_files:` key in the variable in diagnostic in the recipe, and populating it
with a list of desired files e.g.:

.. code-block:: bash

    variables:
      ta:
        preprocessor: my_masking_preprocessor
          fx_files: [sftlf, sftof, sftgif, areacello, areacella]

Such a recipe will automatically retrieve all the `[sftlf, sftof, sftgif, areacello, areacella]`-type fx files for each of the variables that are needed for
and then, in the diagnostic phase, these mask files will be available for the developer to use them as they need to. They `fx_files` attribute of the big `variable`
nested dictionary that gets passed to the diagnostic is, in turn, a dictionary on its own, and members of it can be accessed in the diagnostic through a simple loop over
the `config` diagnostic variable items e.g.:

.. code-block:: bash

    for filename, attributes in config['input_data'].items():
        sftlf_file = attributes['fx_files']['sftlf']
        areacello_file = attributes['fx_files']['areacello']

5. Missing values masks
-----------------------

Missing (masked) values can be a nuisance especially when dealing with multimodel ensembles and having to compute
multimodel statistics; different numbers of missing data from dataset to datest may introduce biases and artifically
assign more weight to the datasets that have less missing data. This is handled in ESMValTool via the missing values
masks: two types of such masks are available: one for the multimodel case and another for the single model case.

The multimodel missing values mask (mask_fillvalues) is a preprocessor step that usually comes after all the single-model
steps (regridding, area selection etc) have been performed; in a nutshell, it combines missing values masks from individual
models into a multimodel missing values mask; the individual model masks are built according to common criteria:
the user chooses a time window in which missing data points are counted, and if the number of missing data points relative
to the number of total data points in a window is less than a chosen fractional theshold, the window is discarded i.e.
all the points in the window are masked (set to missing).

.. code-block:: bash

    preprocessors:
      missing_values_preprocessor:
        mask_fillvalues:
          threshold_fraction: 0.95
          min_value: 19.0
          time_window: 10.0

In the example above, the fractional threshold for missing data vs. total data is set to 95% and the time window is set to
10.0 (units of the time coordinate units). Optionally, a minimum value threshold can be applied, in this case it is set
to 19.0 (in units of the variable units).

A similar preprocessor step exists for the single-dataset: mask_window_threshold (with the same arguments as mask_fillvalues).

6. Min, max and interval masking
--------------------------------

Thresholding on minimum and maximum accepted data values can also be performed: masks are constructed based on the
results of thresholding; inside and outside interval thresholding and masking can also be performed. These functions
are mask_above_threshold, mask_below_threshold, mask_inside_range, and mask_outside_range.

Horizontal regridding
=====================
Documentation of _regrid.py (part 2)

Masking of missing values
=========================
Documentation of _mask.py (part 2)

Multi-model statistics
======================
Documentation of_multimodel.py

Time-area statistics
====================
Documentation of _area_pp.py and _volume_pp.py

Information on maximum memory required
======================================
In the most general case, we can set upper limits on the maximum memory the anlysis will require:


Ms = (R + N) x F_eff - F_eff - when no multimodel analysis is performed;
Mm = (2R + N) x F_eff - 2F_eff - when multimodel analysis is performed;

where

Ms: maximum memory for non-multimodel module
Mm: maximum memory for multimodel module
R: computational efficiency of module; R is typically 2-3
N: number of datasets
F_eff: average size of data per dataset where F_eff = e x f x F
where e is the factor that describes how lazy the data is (e = 1 for fully realized data)
and f describes how much the data was shrunk by the immediately previous module eg
time extraction, area selection or level extraction; note that for fix_data f relates only to the time extraction, if data is exact in time (no time selection) f = 1 for fix_data

so for cases when we deal with a lot of datasets (R + N = N), data is fully realized, assuming an average size of 1.5GB for 10 years of 3D netCDF data, N datasets will require


Ms = 1.5 x (N - 1) GB
Mm = 1.5 x (N - 2) GB
=======
