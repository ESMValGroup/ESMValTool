:: _preprocessor:

***********************
ESMValTool Preprocessor
***********************

Variable derivation
===================
Documentation of _derive.py

Time and area subsetting
========================
Documentation of _time_area.py

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

Time-area statistics
====================
Documentation of _area_pp.py and _volume_pp.py
