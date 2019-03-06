:: _preprocessor:

***********************
ESMValTool Preprocessor
***********************
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


3. extract_month
----------------

The function extracts the times that occur within a specific month.
This function only has one argument: `month`. This value should be an integer
between 1 and 12 as the named month string will not be accepted.


4. time_average
---------------

This functions takes the weighted average over the time dimension. This
function requires no arguments and removes the time dimension of the cube.


5. seasonal_mean
----------------

This function produces a seasonal mean for each season (DJF, MAM, JJA, SON).
Note that this function will not check for missing time points. For instance,
if you are looking at the DJF field, but your datasets starts on January 1st,
the first DJF field will only contain data from January and February.

We recommend using the extract_time to start the dataset from the following
December and remove such biased initial datapoints.


6. annual_mean
--------------

This function produces an annual or a decadal mean. The only argument is the
decadal boolean switch. When this switch is set to True, this function
will output the decadal averages.


Area manipulation
=================
The _area.py module contains the following preprocessor functions:

* extract_region: Extract a region from a cube based on lat/lon corners.
* zonal_means: Calculates the zonal or meridional means.
* average_region: Calculates the average value over a region.
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


2. zonal_means
--------------

The function calculates the zonal or meridional means. While this function is
named `zonal_mean`, it can be used to apply several different operations in
an zonal or meridional direction.
 This function takes two arguments:

* coordinate: Which direction to apply the operation: latitude or longitude
* mean_type: Which operation to apply: mean, stdev, variance, median, min or max


3. average_region
-----------------

This function calculates the average value over a region - weighted by the
cell areas of the region.

This function takes three arguments:
coord1: the name of the coordinate in the first direction.
coord2: the name of the coordinate in the second dimension.
operator: the name of the operation to apply (default: mean).

While this function is named `average_region`, it can be used to apply several
different operations in the horizonal plane: mean, standard deviation, median
variance, minimum and maximum.

Note that this function is applied over the entire dataset. If only a specific
region, depth layer or time period is required, then those regions need to be
removed using other preprocessor operations in advance.


4. extract_named_regions
------------------------

This function extract a specific named region from the data. This function
takes the following argument: `regions` which is either a string or a list
of strings of named regions. Note that the dataset must have a `region`
cooordinate which includes a list of strings as values. This function then
matches the named regions against the requested string.


Volume manipulation
===================
The _volume.py module contains the following preprocessor functions:

* extract_volume: Extract a specific depth range from a cube.
* average_volume: Calculate the volume-weighted average.
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


2. average_volume
-----------------

This function calculates the volume-weighted average across three dimensions,
but maintains the time dimension. The following arguments are required:

coord1: the name of the coordinate in the first direction.
coord2: the name of the coordinate in the second dimension.

No depth coordinate is required as this is determined by iris. This
function works best when the fx_files provide the cell volume.


3. depth_integration
--------------------

This function integrate over the depth dimension. This function does a
weighted sum along the z-coordinate, and removes the z direction of the output
cube. This preprocessor takes no arguments.

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


CMORization and dataset-specific fixes
======================================
Documentation of _reformat.py, check.py and fix.py

Vertical interpolation
======================
Documentation of _regrid.py (part 1)

Land/Sea/Ice Masking
====================
Documentation of _mask.py (part 1)

Certain metrics and diagnostics need to be computed and performed on restricted regions of the Globe; ESMValTool supports subsetting the input data on land mass, oceans and seas, ice. This is achived by masking the model data and keeping only the values associated with grid points that correspond to e.g. land mass
or oceans and seas; masking is done either by using standard mask files that have the same grid resolution as the model data (these files are usually produced
at the same time with the model data and are called fx files) or, in the absence of these files, by using Natural Earth masks. Natural Earth masks, even if they are not model-specific, represent a good approximation since their grid resolution is almost always much higher than the model data, and they are constantly updated with changing
geographical features.

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

Note that for masking out ice the preprocessor is using a different function, this so that both land and sea or ice can be masked out without
losing generality. To mask ice out one needs to add the preprocessing step much as above:

.. code-block:: bash

    preprocessors:
      my_masking_preprocessor:
        mask_landseaice:
          mask_out: ice

To keep only the ice, one needs to mask out landsea, so use that as value for mask_out. As in the case of mask_landsea, the tool will automatically
retrieve the `fx: sftgif` file corresponding the the used variable and extract the ice mask from it.

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
the 'config' diagnostic variable items e.g.:

.. code-block:: bash

    for filename, attributes in config['input_data'].items():
        sftlf_file = attributes['fx_files']['sftlf']
        areacello_file = attributes['fx_files']['areacello']


Horizontal regridding
=====================
Documentation of _regrid.py (part 2)

Masking of missing values
=========================
Documentation of _mask.py (part 2)

Multi-model statistics
======================
Documentation of_multimodel.py
