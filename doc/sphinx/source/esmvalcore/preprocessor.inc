.. _preprocessor:

************
Preprocessor
************
The ESMValTool preprocessor can be used to perform a broad range of operations
on the input data before diagnostics or metrics are applied. The
preprocessor performs these operations in a centralized, documented and
efficient way, thus reducing the data processing load on the diagnostics side.

Each of the preprocessor operations is written in a dedicated python module and
all of them receive and return an Iris cube, working sequentially on the data
with no interactions between them.  The order
in which the preprocessor operations is applied is set by default in order to
minimize the loss of information due to, for example, temporal and spatial
subsetting or multi-model averaging. Nevertheless, the user is free to change
such order to address specific scientific requirements, but keeping in mind
that some operations must be necessarily performed in a specific order. This is
the case, for instance, for multi-model statistics, which required the model to
be on a common grid and therefore has to be called after the regridding module.

In this section, each of the preprocessor modules is described in detail
following the default order in which they are applied:

* `Variable derivation`_.
* `CMOR check and dataset-specific fixes`_.
* `Vertical interpolation`_.
* `Land/Sea/Ice masking`_.
* `Horizontal regridding`_.
* `Masking of missing values`_.
* `Multi-model statistics`_.
* `Time operations`_.
* `Area operations`_.
* `Volume operations`_.
* `Unit conversion`_.


Variable derivation
===================

The variable derivation module allows to derive variables which are not in the
CMIP standard data request using standard variables as input. The typical use
case of this operation is the evaluation of a variable which is only available
in an observational dataset but not in the models. In this case a derivation
function is provided by the ESMValTool in order to calculate the variable and
perform the comparison. For example, several observational datasets deliver
total column ozone as observed variable (`toz`), but CMIP models only provide
the ozone 3D field. In this case, a derivation function is provided to
vertically integrate the ozone and obtain total column ozone for direct
comparison with the observations.

To contribute a new derived variable, it is also necessary to define a name for
it and to provide the corresponding CMOR table. This is to guarantee the proper
metadata definition is attached to the derived data. Such custom CMOR tables
are collected as part of the `ESMValTool core package
<https://github.com/ESMValGroup/ESMValCore>`_. By default, the variable
derivation will be applied only if not already available in the input data, but 
the derivation can be forced by setting the appropriate flag.

.. code-block:: yaml

  variables:
    toz:
      derive: true
      force_derivation: false

The required arguments for this module are two boolean switches:
* derive: activate variable derivation
* force_derivation: force variable derivation even if the variable is
directly available in the input data.

See also :func:`esmvalcore.preprocessor.derive`.


CMOR check and dataset-specific fixes
======================================
.. warning::
   Documentation of _reformat.py, check.py and fix.py to be added


Vertical interpolation
======================
.. warning::
   Documentation of  _regrid.py (part 1) to be added


Land/Sea/Ice masking
====================

Certain metrics and diagnostics need to be computed and performed on specific
domains on the globe. The ESMValTool preprocessor supports filtering
the input data on continents, oceans/seas and ice. This is achived by masking
the model data and keeping only the values associated with grid points that
correspond to, e.g., land, ocean or ice surfaces, as specified by the
user. Where possible, the masking is realized using the standard mask files
provided together with the model data as part of the CMIP data request (the
so-called fx variable). In the absence of these files, the Natural Earth masks
are used: although these are not model-specific, they represent a good
approximation since they have a much higher resolution than most of the models
and they are regularly updated with changing geographical features.

In ESMValTool, land-sea-ice masking can be done in two places: in the
preprocessor, to apply a mask on the data before any subsequent preprocessing
step and before running the diagnostic, or in the diagnostic scripts
themselves. We present both these implementations below.

To mask out a certain domain (e.g., sea) in the preprocessor,
`mask_landsea` can be used:

.. code-block:: yaml

    preprocessors:
      preproc_mask:
        mask_landsea:
          mask_out: sea

and requires only one argument:
* mask_out: either land or sea.

The preprocessor automatically retrieves the corresponding mask (`fx: stfof` in
this case) and applies it so that sea-covered grid cells are set to
missing. Conversely, it retrieves the `fx: sftlf` mask when land need to be
masked out, respectively. If the corresponding fx file is not found (which is
the case for some models and almost all observational datasets), the
preprocessor attempts to mask the data using Natural Earth mask files (that are
vectorized rasters). As mentioned above, the spatial resolution of the the
Natural Earth masks are much higher than any typical global model (10m for
land and 50m for ocean masks).

Note that for masking out ice sheets, the preprocessor uses a different
function, to ensure that both land and sea or ice can be masked out without
losing generality. To mask ice out, `mask_landseaice` can be used:

.. code-block:: yaml

  preprocessors:
    preproc_mask:
      mask_landseaice:
        mask_out: ice

and requires only one argument:
* mask_out: either landsea or ice.

As in the case of `mask_landsea`, the preprocessor automatically retrieves the
`fx: sftgif` mask.

Another option is to just read the fx masks as any other CMOR variable and use
it within a diagnostic script. This can be done in the variable dictionary by
specifiying the desired fx variables (masks):

.. warning::
  Code snippet, text and link to function to be added (after #1037 and #1075
  are closed).


Horizontal regridding
=====================
.. warning::
   Documentation of  _regrid.py (part 2) to be added


Masking of missing values
=========================
.. warning::
   Documentation of  _mask.py (part 2) to be added


Multi-model statistics
======================

.. warning::
   Documentation of _multimodel.py to be added.

Information on maximum memory required: In the most general case, we can set
upper limits on the maximum memory the analysis will require: 

Ms = (R + N) x F_eff - F_eff - when no multimodel analysis is performed;
Mm = (2R + N) x F_eff - 2F_eff - when multimodel analysis is performed;

where

Ms: maximum memory for non-multimodel module
Mm: maximum memory for multimodel module
R: computational efficiency of module (typically 2-3)
N: number of datasets
F_eff: average size of data per dataset where F_eff = e x f x F
where e is the factor that describes how lazy the data is (e = 1 for fully
realized data) and f describes how much the data was shrunk by the immediately
previous module e.g. time extraction, area selection or level extraction; note
that for fix_data f relates  only to the time extraction, if data is exact in
time (no time selection) f = 1 for fix_data.

So for cases when we deal with a lot of datasets (R + N = N), data is fully
realized, assuming an average size of 1.5GB for 10 years of 3D netCDF data, N
datasets will require: 

Ms = 1.5 x (N - 1) GB
Mm = 1.5 x (N - 2) GB


Time operations
===============

The time operations module contains a broad set of functions to subset data and apply
statistical operators along the temporal coordinate of the input data:

| `1. extract_time`_: extract a specified time range from a cube.
| `2. extract_season`_: extract only the times that occur within a specific
      season.
| `3. extract_month`_: extract only the times that occur within a specific
      month. 
| `4. time_average`_: take the weighted average over the entire time dimension.
| `5. seasonal_mean`_: produce a mean for each season (DJF, MAM, JJA, SON)
| `6. annual_mean`_: produce an annual or decadal mean.
| `7. regrid_time`_: align the time axis of each dataset to have common time
      points and calendars.

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

These start and end points are set using the datasets native calendar. All six
arguments should be given as integers, named month strings (e.g., March) will
not be accepted. Note that start_year and end_year can be omitted, as they are
filled in automatically from the dataset definition if not specified
here (end_year will be the value in the dataset definition + 1).

See also :func:`esmvalcore.preprocessor.extract_time`.

2. extract_season
-----------------

Extract only the times that occur within a specific season.

This function only has one argument:

* season: DJF, MAM, JJA, or SON

Note that this function does not change the time resolution. If your original
data is in monthly time resolution, then this function will return three
monthly datapoints per year.

To calculate a seasonal average, this function needs to be combined with the
seasonal_mean function, below.

See also :func:`esmvalcore.preprocessor.extract_season`.

3. extract_month
----------------

The function extracts the times that occur within a specific month.
This function only has one argument:

* month: [1-12]

Note that named month strings will not be accepted.

See also :func:`esmvalcore.preprocessor.extract_month`.

4. time_average
---------------

This function takes the weighted average over the time dimension. This
function requires no arguments and removes the time dimension of the cube.

See also :func:`esmvalcore.preprocessor.time_average`.

5. seasonal_mean
----------------

This function produces a seasonal mean for each season (DJF, MAM, JJA, SON).
Note that this function will not check for missing time points. For instance,
if the DJF field is selected, but the input datasets starts on January 1st,
the first DJF field will only contain data from January and February.

We recommend using the extract_time to start the dataset from the following
December and remove such biased initial datapoints.

See also :func:`esmvalcore.preprocessor.seasonal_mean`.

6. annual_mean
--------------

This function produces an annual or a decadal mean. It takes a single boolean
switch as argument:
* decadal: set this to true to calculate decadal averages instead of annual
averages.

See also :func:`esmvalcore.preprocessor.annual_mean`.

7. regrid_time
--------------

This function aligns the time points of each component dataset to allow the
subtraction of two Iris cubes from different datasets. The operation makes the
datasets time points common and sets common calendars; it also resets the time
bounds and auxiliary coordinates to reflect the artifically shifted time
points. The current implementation works only for monthly and daily data.
        
See also :func:`esmvalcore.preprocessor.regrid_time`.


Area operations
===============

.. warning::
   Need to be adapted after renaming action in #1123

The area manipulation module contains the following preprocessor functions:

| `1. extract_region`_: extract a region from a cube based on lat/lon corners.
| `2. zonal_means`_: calculate the zonal or meridional means.
| `3. area_statistics`_: calculate the average value over a region.
| `4. extract_named_regions`_: extract a region from a cube given its name.

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

* coordinate: Which direction to apply the operation: latitude or longitude.
* mean_type: Which operation to apply: mean, std_dev, variance, median, min or
* max. 

See also :func:`esmvalcore.preprocessor.zonal_means`.

3. area_statistics
------------------

This function calculates the average value over a region - weighted by the
cell areas of the region.

This function takes one argument:
* operator: the name of the operation to apply.

This function can be used to apply several different operations in the
horizonal plane: mean, standard deviation, median variance, minimum and
maximum.

Note that this function is applied over the entire dataset. If only a specific
region, depth layer or time period is required, then those regions need to be
removed using other preprocessor operations in advance.

See also :func:`esmvalcore.preprocessor.area_statistics`.

4. extract_named_regions
------------------------

This function extract a specific named region from the data. 
This function takes onw argument: 

* regions: either a string or a list of strings of named regions. 

Note that the dataset must have a `region` cooordinate which includes a list of
strings as values. This function then matches the named regions against the
requested string.

See also :func:`esmvalcore.preprocessor.extract_named_regions`.


Volume operations
=================

The volume operations module contains the following preprocessor functions:

| `1. extract_volume`_: extract a specific depth range from a cube.
| `2. volume_statistics`_: calculate the volume-weighted average.
| `3. depth_integration`_: integrate over the depth dimension.
| `4. extract_transect`_: extract data along a line of constant latitude or
      longitude. 
| `5. extract_trajectory`_: extract data along a specified trajectory.

1. extract_volume
-----------------

This function extracts a specific range in the z-direction from a cube.
This function takes two arguments:

* z_min: minimum in the z direction
* z_max: maximum in the z direction

Note that this requires the requested z-coordinate range to be the same sign as
the Iris cube, i.e. if the cube has z-coordinate as negative, then z_min and
z_max need to be negative numbers.

See also :func:`esmvalcore.preprocessor.extract_volume`.

2. volume_statistics
--------------------

This function calculates the volume-weighted average across three dimensions,
but maintains the time dimension.

This function takes one argument:
* operator: operation to apply over the volume (at the moment only mean is implemented)

No depth coordinate is required as this is determined by Iris. This
function works best when the fx files provide the cell volume.

See also :func:`esmvalcore.preprocessor.volume_statistics`.


3. depth_integration
--------------------

This function integrates over the depth dimension. It performs a weighted sum
along the z-coordinate, and removes the z direction of the output cube. It takes no arguments.

See also :func:`esmvalcore.preprocessor.depth_integration`.

4. extract_transect
-------------------

This function extracts data along a line of constant latitude or longitude.
This function takes two arguments, although only one is strictly required:
* latitude
* longitude

One of these arguments needs to be set to a float, and the other can then be
either ignored or set to a minimum or maximum value. For example, if latitude
is set to 0 and longitude is left blank, the function would produce a cube
along the equator. If latitude is set to to 0 and longitude to `[40., 100.]` it
will produce a transect of the equator in the Indian Ocean.

See also :func:`esmvalcore.preprocessor.extract_transect`.

5. extract_trajectory
---------------------

This function extracts data along a specified trajectory. It requires three
arguments:
* latitude_points: list of latitude coordinates
* longitude_points: list of longiute coordinates
* number_points: if two points are provided, the `number_points` argument is
used to set the number of places to extract between the two end points.

If more than two points are provided, then extract_trajectory will produce a
cube which has extrapolated the data of the cube to those points, and
`number_points` is not needed. Note that this function uses the expensive
interpolate method, but it may be necessary for irregular grids.

See also :func:`esmvalcore.preprocessor.extract_trajectory`.


Unit conversion
===============
