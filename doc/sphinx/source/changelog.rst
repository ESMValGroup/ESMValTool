.. _changelog:

Changelog
=========

.. _changelog-v2-13-0:

v2.13.0
-------

Highlights
~~~~~~~~~~

Constellation of CMORizers and diagnostics updated/implemented as part of the
`Rapid Evaluation Framework (REF) <https://wcrp-cmip.org/cmip-phases/cmip7/rapid-evaluation-framework/>`__:

- CMORizers:

   - Update of the ESACCI-OZONE CMORizer
   - ESACCI-SEAICE CMORizer
   - Update of the OSI-450 (OSISAF/CCI) CMORizer

- :ref:`Diagnostics <recipes_REF>`:

   - Seasonal cycle and time series of Arctic/Antarctic sea ice area
   - :ref:`Zero Emissions Commitment (ZEC) <recipes_zec>`
   - Cloud radiative effects
   - :ref:`Transient Climate Response to Emissions (TCRE) <recipes_tcre>`
   - :ref:`ENSO metrics <recipes_enso_metrics>`
   - 2D histograms and scatterplots for clouds
   - Regional historical changes
   - :ref:`Evaluation of climate variables at Global warming levels <recipes_recipe_calculate_gwl_exceedance_stats>`
   - :ref:`Sea surface temperature bias (IPCC AR6 Ch3 Fig. 3.24) <recipes_ipccwg1ar6ch3>`
   - :ref:`Climate drivers for fire <recipe_ref_fire>`
   - :ref:`Rate of sea ice area loss per degree of warming <recipe_seaice_sensitivity>`
   - Ozone metrics and climatology

This release includes

Backwards incompatible changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Refactor diagnostic `monitor/multi_datasets.py` (:pull:`4002`) by :user:`schlunma`
   Plot option `annual_mean_kwargs`: This will be ignored now. Previously,
   this could lead to errors in rare cases. Instead, explicitly specify a
   dataset that uses the annual_statistics preprocessor in the recipe,
   for example `esmvaltool/recipes/monitor/recipe_monitor_with_refs.yml#L132-L139`.
-  Retire recipe schlund20jgr (:pull:`4022`) by :user:`axel-lauer`
-  Retire julia recipes (:pull:`4045`) by :user:`bettina-gier`
-  Retire psyplot recipe and diagnostic (:pull:`4063`) by :user:`valeriupredoi`

Community
~~~~~~~~~

-  Added stale action (:pull:`4100`) by :user:`schlunma`

Documentation
~~~~~~~~~~~~~

-  Update citation and changelog for version 2.12 (:pull:`3945`) by :user:`sloosvel`
-  Use Python-only environment for ReadTheDocs builds (:pull:`3955`) by :user:`valeriupredoi`
-  Add explanatory note in documentation advice on Codacy (:pull:`4020`) by :user:`alistairsellar`
-  Extend broken recipe policy (:pull:`3913`) by :user:`axel-lauer`
-  Remove some more Julia references from documentation (but keep some) (:pull:`4073`) by :user:`valeriupredoi`
-  Add footer with legal info (:pull:`4072`) by :user:`axel-lauer`
-  Add mention of NSC HPC servers (:pull:`4084`) by :user:`jlenh`
-  Update figure selection and layout of the gallery (:pull:`4060`) by :user:`bettina-gier`
-  Update timeline for release 2.13.0 (:pull:`4096`) by :user:`jlenh`
-  Refer to esmvaltool.org webpage for development team information (:pull:`3892`) by :user:`bouweandela`
-  Remove mention of Julia in documentation ahead of discontinued support in v2.13 (:pull:`4156`) by :user:`jlenh`
-  Update `cartopy` documentation link to fix documentation build fail (:pull:`4179`) by :user:`jlenh`
-  Add debugging help (:pull:`2638`) by :user:`zklaus`

Diagnostics
~~~~~~~~~~~

-  Add recipe for seasonal cycle and time series of Arctic/Antarctic sea ice area (REF) (:pull:`3891`) by :user:`axel-lauer`
-  Add recipe and diagnostic for ZEC (Zero Emissions Commitment) (:pull:`3897`) by :user:`bettina-gier`
-  Create scatterplots for the REF (:pull:`3923`) by :user:`LisaBock`
-  Evaluate key climate variables at Global warming levels (:pull:`3927`) by :user:`rswamina`
-  Add SST bias diagnostic for IPCC AR6 Ch3 Fig. 3.24. Can be used for REF. (:pull:`3944`) by :user:`malininae`
-  Update recipe_combined_indices.yml for R versions larger than 4.2 (:pull:`3967`) by :user:`bouweandela`
-  Added option to plot horizontal lines in monitoring line plots (:pull:`3977`) by :user:`schlunma`
-  Fix R recipes for R > 4.2 (:pull:`3985`) by :user:`sloosvel`
-  REF ozone recipe (:pull:`3988`) by :user:`hb326`
-  Faster and lower memory version of recipe_ref_scatterplot.yml (:pull:`3992`) by :user:`bouweandela`
-  Update GPCP and ERA5 obs4MIPs datasets in recipe_ref_scatterplot.yml (:pull:`4011`) by :user:`bouweandela`
-  Add option to specify caption in recipe for `monitor/multi_datasets.py` (:pull:`4028`) by :user:`schlunma`
-  Optimize caption produced by `tcr.py` (:pull:`4042`) by :user:`schlunma`
-  Optimize caption produced by `ecs.py` (:pull:`4040`) by :user:`schlunma`
-  Fix tcre metric diagnostic (:pull:`4044`) by :user:`valeriupredoi`
-  Optimize `tcre.py` diagnostic (caption and usability via Python API) (:pull:`4056`) by :user:`schlunma`
-  Add additional diagnostics to recipe_ocean_quadmap.yml (:pull:`3953`) by :user:`SophieHall2024`
-  New CH4 lifetime diagnostic (:pull:`3507`) by :user:`FranziskaWinterstein`
-  New recipe: rate of sea ice area loss per degree warming (:pull:`3983`) by :user:`alistairsellar`
-  Remove unused option `annual_mean_kwargs` from all recipes (:pull:`4133`) by :user:`schlunma`
-  Aeronet format updates and diagnostic changes (:pull:`4131`) by :user:`nchawang`
-  Use obs4MIPs consistently in all recipes (:pull:`4142`) by :user:`bouweandela`
-  Add missing method in `seaice_drift.py` (:pull:`4161`) by :user:`sloosvel`
-  Update the name of the `remapscon2` operator in R recipes (:pull:`4164`) by :user:`jlenh`
-  Trace gases (CO2, CH4, N2O) surface concentration diagnostic (:pull:`4014`) by :user:`jlenh`
-  Climate drivers for fire for REF (:pull:`3975`) by :user:`jlenh`
-  Exchange toz reference data in ozone REF recipe (:pull:`4008`) by :user:`axel-lauer`
-  Add captions to REF recipe for basic sea ice area diags (:pull:`4033`) by :user:`axel-lauer`
-  CLIVAR ENSO metrics for REF (:pull:`3972`) by :user:`flicj191`
-  Update and combine drought diagnostics (:pull:`3907`) by :user:`lukruh`
-  Regional historical changes for REF (:pull:`3926`) by :user:`katjaweigel`
-  Add more datasets to recipe_easy_ipcc.yml (:pull:`3920`) by :user:`bouweandela`
-  ACCESS-NRI ENSO diagnostics CLIVAR porting (:pull:`3905`) by :user:`flicj191`
-  Correct typo in "pearsonr" (:pull:`4140`) by :user:`katjaweigel`
-  Update time ranges in REF 'cloud radiative effect' recipe (:pull:`4057`) by :user:`LisaBock`
-  Fix caption and time range for REF diagnostic "Cloud scatterplots" (:pull:`4066`) by :user:`LisaBock`
-  Portrait plot example: remove sea surface temperature (:pull:`4099`) by :user:`lukruh`
-  Recipes ipccwg1ar5ch9: Update CERES and GPCP versions to match ESGF (:pull:`3440`) by :user:`rbeucher`
-  Recipe_schlund20esd Remove CERES-EBAF version to fix ESGF search (:pull:`3441`) by :user:`rbeucher`
-  Adding basic sanity checks (:pull:`4120`) by :user:`LisaBock`
-  Adding optional data output for seaborn_diag.py (:pull:`4144`) by :user:`katjaweigel`
-  Fix `recipe_williams09climdyn_CREM.yml` by replacing `np.NAN` with `np.nan` (:pull:`4160`) by :user:`jlenh`
-  Fix `recipe_carvalhais14nat.yml` to avoid the usage of `vars()` inside a function (:pull:`4170`) by :user:`jlenh`
-  Fix small matplotlib errors in `recipe_zmnam.yml` (:pull:`4177`) by :user:`jlenh`
-  Fix `_fix_lons` function in `recipe_flato13ipcc_figure_914.yml` due to `numpy` error [#4034] (:pull:`4178`) by :user:`jlenh`
-  Improve the handling of options in recipe_ref_fire.yml (:pull:`4184`) by :user:`bouweandela`
-  Fix small matplotlib errors in `recipe_zmnam.yml` (:pull:`4177`) by :user:`jlenh`
-  Fix minor NCL errors for `recipe_russel18jgr.yml` (:pull:`4183`) by :user:`jlenh`
-  Fix specification for WOA dataset in recipe_ocean_quadmap (:pull:`4176`) by :user:`TomasTorsvik`
-  Fix autoassess normalise error (:pull:`4187`) by :user:`alistairsellar`

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Update ESACCI-OZONE CMORizer(formatter and downloader) for REF (:pull:`3899`) by :user:`diegokam`
-  Update CMORizer for OSI-450 aka OSISAF/CCI (sea ice concentration) (:pull:`3963`) by :user:`axel-lauer`
-  Adding HadEX3 CMORizer (:pull:`3636`) by :user:`malininae`
-  Downloader & cmorizer for NOAA GML surface flask trace gases  (:pull:`3847`) by :user:`jlenh`
-  Use nested update to parse configuration options given to CMORizer interface (:pull:`4113`) by :user:`schlunma`
-  Aeronet-header-fix to handle V3 data (:pull:`4101`) by :user:`nchawang`
-  Cmorizer for the Wetland Area and Dynamics for Methane Modeling (WAD2M) dataset (:pull:`4110`) by :user:`jlenh`
-  Cmorizer for the Global Lakes and Wetlands Database (GLWD) version 2.0 dataset (:pull:`4112`) by :user:`jlenh`
-  Update ESACCI Cloud CMORizer (daily and monthly data)  (:pull:`3756`) by :user:`diegokam`
-  Fix coords in ESACCI-CLOUD cmorizer (:pull:`4129`) by :user:`LisaBock`
-  Added CMORizer for Yang2020 data (:pull:`4090`) by :user:`schlunma`
-  Add CMORized NOAA GML surface flask concentration measurement of N2O (:pull:`4059`) by :user:`jlenh`
-  CMORizer for ESACCI-BIOMASS (:pull:`4121`) by :user:`axel-lauer`
-  Faster OSI-450-nh and OSI-450-sh CMORizer and respect Dask settings (:pull:`3990`) by :user:`bouweandela`
-  CMORiser for the IAPv4.2 global ocean temperature dataset (:pull:`3887`) by :user:`rbeucher`
-  Update recipe_check_obs.yml to match the latest AERONET data (:pull:`4171`) by :user:`hb326`
-  Update version facet for ESACCI-CLOUD in recipe_check_obs.yml (:pull:`4166`) by :user:`jlenh`
-  Remove `end_year` facet for NOAA-GML datasets in `recipe_check_obs.yml` (:pull:`4167`) by :user:`jlenh`
-  New version of JRA-25 CMORizer (:pull:`3470`) by :user:`axel-lauer`
-  Add CAMS cmorizer (:pull:`3749`) by :user:`bettina-gier`

Automatic testing
~~~~~~~~~~~~~~~~~

-  Use mamba>=2 for Circle CI upstream development test (:pull:`3949`) by :user:`valeriupredoi`
-  Use python-only environment file for Circle build documentation (:pull:`3981`) by :user:`valeriupredoi`
-  Use conda exec instead of mamba for circleCI `build_documentation` test (:pull:`4003`) by :user:`valeriupredoi`
-  Run pre-commit on main (:pull:`4026`) by :user:`schlunma`
-  Try set a longer no output time limit for Circle tests (:pull:`4035`) by :user:`valeriupredoi`
-  Only run imagehash tests on Linux machines (:pull:`4039`) by :user:`schlunma`
-  Do not run conda install package Github Action with Python 3.13 (yet) (:pull:`4071`) by :user:`valeriupredoi`
-  Fix Cicle CI tests (prepend Julia to PATH for the containerized run tests) (:pull:`4086`) by :user:`valeriupredoi`
-  Automatically open pull requests to update GitHub Actions (:pull:`4143`) by :user:`valeriupredoi`
-  Run esmvalbot test PR twice a month (:pull:`3976`) by :user:`bouweandela`
-  Add DKRZ resources for RTW to `dkrz-recipes.jinja` (:pull:`3986`) by :user:`ehogan`
-  Add JASMIN resources for RTW to `jasmin-recipes.jinja` (:pull:`3987`) by :user:`ehogan`
-  Use site directories in the RTW (:pull:`3984`) by :user:`ehogan`
-  Update MO paths in the RTW (:pull:`4006`) by :user:`ehogan`
-  Select directories for comparison carefully in the RTW (:pull:`4007`) by :user:`ehogan`
-  Implement additional user configuration options in RTW (:pull:`4009`) by :user:`ehogan`
-  Generate "recipe status" HTML after each Recipe Test Workflow cycle (:pull:`4013`) by :user:`chrisbillowsMO`
-  Update RTW to use the `git` scheme (:pull:`4030`) by :user:`ehogan`
-  Add housekeeping to RTW (:pull:`4019`) by :user:`ehogan`
-  Revert "Update RTW to use the `git` scheme" (:pull:`4041`) by :user:`ehogan`
-  Update location of `recipe_consecdrydays.yml` in the RTW (:pull:`4052`) by :user:`ehogan`
-  Monitor Tests workflow fails (:pull:`4065`) by :user:`chrisbillowsMO`
-  RTW status report shows missing compare tasks as succeeded (:pull:`4064`) by :user:`chrisbillowsMO`
-  Use compare command in recipe test workflow (:pull:`4087`) by :user:`alistairsellar`
-  Separate `BRANCH` variable for recipe_test_workflow cloning Core and Tool (:pull:`4070`) by :user:`alistairsellar`
-  Update RTW to use new user configuration (:pull:`4017`) by :user:`ehogan`

Installation
~~~~~~~~~~~~

-  Remove Julia from conda dependency environment (but keep tests for `esmvaltool install Julia`)
   and add docs how to install from source (:pull:`3921`) by :user:`valeriupredoi`
-  Unpin upper limit `r-base` (:pull:`3971`) by :user:`valeriupredoi`
-  Drop support for Python 3.10 (:pull:`4027`) by :user:`schlunma`
-  Support Python 3.13 (:pull:`3805`) by :user:`valeriupredoi`
-  Tackle Julia in DockerHub container builds: remove nstallation of packages
   (and remove dependencies in Project.jl bar :beer: YAML) (:pull:`4083`) by :user:`valeriupredoi`
-  Use `pyproject.toml` and retire `setup.py` with NO Ruff configuration changes (:pull:`4088`) by :user:`valeriupredoi`
-  Update environment to use ESMValCore v2.13.0rc1 (:pull:`4169`) by :user:`bouweandela`
-  Update ESMValCore version also for OSX (:pull:`4174`) by :user:`bouweandela`

Improvements
~~~~~~~~~~~~

-  Replace flake8, yapf, and isort with ruff (:pull:`3893`) by :user:`bouweandela`
-  Remove `test_recipe` command (:pull:`4031`) by :user:`alistairsellar`
-  Replace deprecated horizontal regridding scheme `linear_extrapolate` with generic scheme (:pull:`4081`) by :user:`schlunma`
-  Add compare command, via new develop command group (:pull:`4054`) by :user:`alistairsellar`
-  Add extra Ruff rules and switch to `requires-python = ">=3.11"` in pyproject.toml (:pull:`4085`) by :user:`valeriupredoi`
-  Move dependabot configuration to correct path (:pull:`4145`) by :user:`bouweandela`
-  Use trusted publisher routine for PyPI and retire token-based type (:pull:`4150`) by :user:`valeriupredoi`
-  Fix pypa Github Action version (:pull:`4151`) by :user:`valeriupredoi`

.. _changelog-v2-12-0:

v2.12.0
-------

Highlights
~~~~~~~~~~

- Recipe :ref:`recipe_portrait_CMIP.yml <recipe_portrait>` allows the creation
  of portrait plots to visualize performance metrics.

- A set of new :ref:`recipes <recipe_benchmarking>` allow plotting
  arbitrary preprocessor output to benchmark a model simulation with
  other datasets.

This release includes

Backwards incompatible changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Remove recipe filler utility (:pull:`3777`) by :user:`schlunma`.
   ESMValTool accepts the use of wildcards in the recipe, and a filled
   version of the recipe is automatically created and stored in the
   `run` directory.

Deprecations
~~~~~~~~~~~~

-  Adapt ESMValTool to new configuration (:pull:`3761`) by :user:`schlunma`

   - The single configuration file ``config-user.yml`` has been deprecated in favour of configuration directories.
     By default, the directory ``~/.config/esmvaltool`` will be considered.
     To switch to the new format run:

       .. code-block:: bash

          mkdir -p ~/.config/esmvaltool && mv ~/.esmvaltool/config-user.yml ~/.config/esmvaltool

     You can also specify the location of the configuration directory with the ``--config_dir`` flag.
     Please refer to :ref:`esmvalcore:config` for a detailed description on how to configure the tool.

Bug fixes
~~~~~~~~~

-  Fix for setting global attributes in cmorizers (:pull:`3717`) by :user:`LisaBock`

Community
~~~~~~~~~

-  Add a security policy (:pull:`3456`) by :user:`ehogan`

Documentation
~~~~~~~~~~~~~

-  Ignore autosummary warning in documentation build (:pull:`3718`) by :user:`valeriupredoi`
-  Retire mambaforge - addendum to 3774 (:pull:`3778`) by :user:`valeriupredoi`
-  Update Docker builds badge in README (:pull:`3783`) by :user:`valeriupredoi`
-  Readthedocs configuration/builds: revert to miniconda before miniforge is available (:pull:`3785`) by :user:`valeriupredoi`
-  Add info on obs tiers to docu (:pull:`3624`) by :user:`axel-lauer`
-  Fix wrong version in title in Readthedocs documentation (:pull:`3815`) by :user:`valeriupredoi`
-  Removed wrong description of option in portrait plot docstring (:pull:`3908`) by :user:`schlunma`
-  Update broken recipes (:pull:`3940`) by :user:`sloosvel`
-  Add next release schedule (:pull:`3794`) by :user:`sloosvel`

Diagnostics
~~~~~~~~~~~

-  More flexible file loading in `monitor/multi_datasets.py` (:pull:`3728`) by :user:`schlunma`
-  Add option to plot time on x-axis in monitoring Hovmoeller plots (:pull:`3732`) by :user:`FranziskaWinterstein`
-  Use `transform_first=True` for contourf plots with Robinson projection to avoid cartopy bug (:pull:`3789`) by :user:`schlunma`
-  Fix contourf plots for masked data (:pull:`3797`) by :user:`schlunma`
-  Fix issue related to removal/change of private function imported in `diag_scripts/shared/_supermeans.py` (deprecation in iris=3.11) (:pull:`3810`) by :user:`valeriupredoi`
-  Allow setting `matplotlib.rcParams` in `monitor/multi_datasets.py` (:pull:`3844`) by :user:`schlunma`
-  Adding figures of Bock and Lauer (2024) (:pull:`3526`) by :user:`LisaBock`
-  Added Python portrait plot diagnostic (:pull:`3551`) by :user:`lukruh`
-  Fix Autoassess diagnostic for new matplotlib 3.10 api change (:pull:`3917`) by :user:`valeriupredoi`
-  Fix diagnostic mpqb (:pull:`3918`) by :user:`sloosvel`
-  Add diagnostic to calculate Transient Climate Response to Emissions (TCRE) (:pull:`3904`) by :user:`schlunma`
-  Benchmarking recipes (Lauer et al.) (:pull:`3598`) by :user:`axel-lauer`
-  Add Cloud Radiative Effects recipe for REF (:pull:`3903`) by :user:`LisaBock`

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Update esacci-soilmoisture(v08.1) downloader and CMORizer (Python version) (:pull:`3676`) by :user:`diegokam`
-  Update ESACCI Landcover CMORizer (python version) and downloader (pft yearly data, v2.0.8) (:pull:`3727`) by :user:`diegokam`
-  CMORizer for JRA-55 (:pull:`3141`) by :user:`axel-lauer`
-  Replace underscores in facets of observational datasets with minuses (:pull:`3840`) by :user:`schlunma`
-  Add TropFlux CMORiser (:pull:`3863`) by :user:`rbeucher`
-  CMORizer for ESACCI-SEAICE (:pull:`3821`) by :user:`axel-lauer`
-  Update ESACCI-SST cmorizer to v3.0 (:pull:`3697`) by :user:`LisaBock`
-  Adding pr, tauu, tauv NOAA-CIRES-20CR-V2 CMORISER (:pull:`3763`) by :user:`max-anu`
-  Adding pr, tauu, tauv, tos to NCEP2 CMORISer (:pull:`3765`) by :user:`max-anu`
-  Adding a CMORiser for CMAP data for pr (:pull:`3766`) by :user:`max-anu`
-  Update NSIDC_G02202_sh CMORiser to add bounds for lat,lon and time (:pull:`3744`) by :user:`flicj191`

Automatic testing
~~~~~~~~~~~~~~~~~

-  Remove obsolete utility `esmvt_rose_wrapper` and its documentation and very obsolete `mip_convert` cmorizer (:pull:`3759`) by :user:`valeriupredoi`
-  Remove obsolete and inactive `tests/system` tests (:pull:`3760`) by :user:`valeriupredoi`
-  Retire Mambaforge (:pull:`3774`) by :user:`valeriupredoi`
-  Pin mamba<2 for conda-lock: solution by Ben Mares @maresb (:pull:`3771`) by :user:`valeriupredoi`
-  Update comment in conda lock creation Github action (:pull:`3788`) by :user:`valeriupredoi`
-  Run a periodic esmvalbot test via automated PR (:pull:`3465`) by :user:`valeriupredoi`
-  Add init files in package-like diagnostic directories to comply with `pytest >=8.3.4` (:pull:`3834`) by :user:`valeriupredoi`
-  Fix failing sklearn test (:pull:`3850`) by :user:`schlunma`
-  Fix circle ci nightly test installation from source development mode (and test mode) that times out (:pull:`3864`) by :user:`valeriupredoi`
-  Remove portrait plot test recipe (:pull:`3871`) by :user:`schlunma`
-  Restrict cron Github Actions tests to run from forks (:pull:`3894`) by :user:`valeriupredoi`
-  Force mamba 2 in Github Actions that solve environment (:pull:`3929`) by :user:`valeriupredoi`
-  Force install mamba2 in CircleCI and add ipython and pin it to <9.0 (:pull:`3947`) by :user:`valeriupredoi`
-  Use mamba>=2 for Circle CI upstream development test  (:pull:`3949`) by :user:`valeriupredoi`

Installation
~~~~~~~~~~~~

-  Retire support for Python 3.9 (:pull:`3683`) by :user:`valeriupredoi`
-  [Julia] pin `curl <8.10` to restrict `libcurl <8.10` so Julia installs packages correctly (:pull:`3755`) by :user:`valeriupredoi`
-  Add support for Python=3.12 (:pull:`3501`) by :user:`valeriupredoi`
-  Pin cartopy to `cartopy<0.24` (:pull:`3768`) by :user:`valeriupredoi`
-  Pin pys2index >=0.1.5 in osx environment (:pull:`3792`) by :user:`valeriupredoi`
-  Update environment: pin `iris>=3.11`, unpin `cartopy` and allow for `numpy >=2` (:pull:`3811`) by :user:`valeriupredoi`
-  Unpin pandas  (:pull:`3924`) by :user:`valeriupredoi`

Improvements
~~~~~~~~~~~~

-  Update ERA5 renaming script for hourly (:pull:`3630`) by :user:`malininae`
-  Avoid masking issues in Dask 2024.8.0 (:pull:`3736`) by :user:`bouweandela`
-  Dark mode compatible transparent background logo (:pull:`3751`) by :user:`lukruh`
-  Change authors name (:pull:`3806`) by :user:`lukruh`
-  Recipe Test Workflow (RTW) prototype (:pull:`3210`) by :user:`ehogan`
-  Update obs4MIPs names to match ESGF see #2974 (:pull:`3439`) by :user:`rbeucher`
-  Recipe_wenzel16jclim: Remove CERES-EBAF version to fix ESGF (:pull:`3442`) by :user:`rbeucher`
-  Update RTW checks (:pull:`3823`) by :user:`ehogan`
-  Add cccma collaborators (:pull:`3831`) by :user:`malininae`
-  Add DKRZ as a site to RTW (:pull:`3837`) by :user:`ehogan`
-  Add an option to disable using the distributed scheduler from the diagnostic script (:pull:`3787`) by :user:`bouweandela`
-  Use site-specific lists of recipes in the RTW (:pull:`3856`) by :user:`chrisbillowsMO`
-  Add links to summary page that are required for hosting at DKRZ (:pull:`3866`) by :user:`bouweandela`
-  Add `filename` to list of ignored attributes in `compare.py` (:pull:`3919`) by :user:`sloosvel`
-  Add unit conversion to `DU` to recipes that use `toz` (:pull:`3784`) by :user:`schlunma`

.. _changelog-v2-11-0:

v2.11.0
-------
Highlights

- Two new recipes have been added:

  - Recipe :ref:`recipe_aod_aeronet_assess.yml <recipe_aod_aeronet_assess>`
    evaluates model aerosol optical depth (AOD) climatologies against ground
    based observations from the AeroNET measurement network.
  - Recipe :ref:`recipe_climate_patterns.yml <recipes_climate_patterns>`
    generates climate patterns from CMIP6 model datasets.

- The ESACCI-WATERVAPOUR CMORizer now includes daily data and uses the
  officially released CDR2 data.
- Support for 5 new datasets have been added:

  - AeroNET
  - ANU Climate 2.0 Australian data
  - Australian Gridded Climate Data(AGCD) precipitation
  - NOAA-ERSST
  - NSIDC-G02202-sh sea ice fraction

- NEW TREND: First time release manager shout-outs!

  - This is the first ESMValTool release managed by the Met Office! We want to
    shout this out - and for all future first time release managers to
    shout-out - to celebrate the growing, thriving ESMValTool community.

This release includes

Bug fixes
~~~~~~~~~

-  Recipe_ocean_quadmap: Update ATSR to match ESGF name (:pull:`3443`) by :user:`rbeucher`
-  Fix recipe_bock20jgr_fig_8-10.yml (:pull:`3665`) by :user:`LisaBock`
-  Update the list of datasets used in ``recipe_easy_ipcc.yml`` (:pull:`3710`) by :user:`bouweandela`

Documentation
~~~~~~~~~~~~~

-  Improve release tools and documentation (:pull:`3462`) by :user:`bouweandela`
-  Fix a typo in the references file (:pull:`3499`) by :user:`bouweandela`
-  Fix recipe path in ``recipe_perfmetrics.rst`` (:pull:`3532`) by :user:`TomasTorsvik`
-  Improved description of model evaluation recipes (:pull:`3541`) by :user:`schlunma`
-  Remove double word in cmorizer documentation (:pull:`3553`) by :user:`bettina-gier`
-  Fix Codacy badge (:pull:`3558`) by :user:`bouweandela`
-  Update the release schedule for v2.11.0 (:pull:`3573`) by :user:`ehogan`
-  Improve the formatting of the recipe documentation template (:pull:`3652`) by :user:`mo-gill`
-  Add introduction material on the main documentation page (:pull:`3628`) by :user:`bouweandela`
-  Avoid warning in documentation build (:pull:`3675`) by :user:`bouweandela`
-  Update the list of broken recipes for ``v2.11.0`` (:pull:`3706`) by :user:`ehogan`

Diagnostics
~~~~~~~~~~~

-  ``monitor/multi_dataset.py`` improvements: allow data w/o ``timerange`` and improve text formatting (:pull:`3528`) by :user:`schlunma`
-  Allow datasets without ``project`` in multi_datasets.py (:pull:`3552`) by :user:`schlunma`
-  Prevent overlapping time axis tick labels in monitoring recipe (:pull:`3682`) by :user:`schlunma`

New recipe
~~~~~~~~~~

-  Add support for aerosol optical depth climatology metrics to the AutoAssess replacement (:pull:`3048`) by :user:`catherinehardacre`
-  CMIP6 climate patterns  (:pull:`2785`) by :user:`mo-gregmunday`

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Add cmorizer scripts for NOAA-ERSST. (:pull:`1799`) by :user:`bjoernbroetz`
-  Update OceanSODA-ETHZ CMORizer with new source file (:pull:`3535`) by :user:`TomasTorsvik`
-  Add CMORizer script for NSIDC-G02202-sh sea ice fraction (:pull:`3512`) by :user:`flicj191`
-  CMORizer Australian Gridded Climate Data(AGCD) precipitation (:pull:`3445`) by :user:`flicj191`
-  Extend CMORizer NCEP-DOE-R2 (:pull:`3469`) by :user:`axel-lauer`
-  Add comment to recipe_lauer13jclim regarding UWisc being superseded by MAC-LWP (:pull:`3537`) by :user:`rbeucher`
-  Recipe_autoassess_landsurface_surfrad: Remove CERES-EBAF version to fix ESGF search (:pull:`3438`) by :user:`rbeucher`
-  Updating ESACCI-WATERVAPOUR cmorizer (:pull:`3282`) by :user:`malininae`
-  CMORiser for ANU Climate 2.0 Australian data (:pull:`3511`) by :user:`flicj191`
-  Add AERONET cmorizer (:pull:`3227`) by :user:`zklaus`
-  Update CRU CMORizer (:pull:`3381`) by :user:`lukruh`
-  Fix recipe_check_obs to be aligned with DKRZ (:pull:`3673`) by :user:`LisaBock`
-  Update AERONET data version (:pull:`3692`) by :user:`ehogan`

Automatic testing
~~~~~~~~~~~~~~~~~

-  Move code into function in batch job generation script (:pull:`3491`) by :user:`bouweandela`
-  Fix sklearn tests (:pull:`3506`) by :user:`schlunma`
-  Pinned sklearn>=1.4.0 (:pull:`3508`) by :user:`schlunma`
-  Update sklearn tests to be compatible with current pytest version (pytest >=8.0.0) (:pull:`3517`) by :user:`schlunma`
-  Update sklearn tests to be compatible with current pytest version (pytest >=8.0.0) Part 2 (:pull:`3518`) by :user:`schlunma`
-  [Circle/CI]Fix `test_installation_from_conda` Circle CI tests (:pull:`3538`) by :user:`valeriupredoi`
-  [Github Actions] install git in OSX and add environment inspection (:pull:`3581`) by :user:`valeriupredoi`
-  [CI Github Actions] Update (outdated) actions versions that produce Node.js warnings (:pull:`3586`) by :user:`valeriupredoi`
-  Fix ``flake8==7`` linting issues (:pull:`3634`) by :user:`valeriupredoi`
-  Use ``importlib`` as the import mode for ``pytest`` (:pull:`3672`) by :user:`ehogan`

Installation
~~~~~~~~~~~~

-  Update dependencies (:pull:`3487`) by :user:`bouweandela`
-  Merge v2.10.x into main (:pull:`3489`) by :user:`schlunma`
-  Add imagehash package as an ESMValTool dependency (:pull:`3557`) by :user:`alistairsellar`
-  Unpin ``r-akima`` (:pull:`3564`) by :user:`valeriupredoi`
-  Adding pys2index dependency (:pull:`3577`) by :user:`ljoakim`
-  Pin esmpy <8.6.0 (:pull:`3585`) by :user:`valeriupredoi`
-  Pin R <4.3.0 (:pull:`3689`) by :user:`ehogan`
-  Pin importlib_metadata <8 (:pull:`3700`) by :user:`ehogan`
-  Pin matplotlib <3.9.0 on ESMValTool release branch (:pull:`3712`) by :user:`ehogan`

Dependency updates
~~~~~~~~~~~~~~~~~~

-  Fix for ``recipe_seaice_drift.yml``: fix CRS transformer for "North Pole Stereographic" (:pull:`3531`) by :user:`flicj191`
-  Fixed attribute handling in austral_jet/main.ncl for iris>=3.8 (:pull:`3603`) by :user:`schlunma`
-  Fixed attribute handling in emergent constraint diagnostic for iris>=3.8 (:pull:`3605`) by :user:`schlunma`
-  Update the name of the remapcon2 operator in R recipes (:pull:`3611`) by :user:`ehogan`
-  Use ``iris.FUTURE.save_split_attrs = True`` to remove iris warning in many diagnostics (:pull:`3651`) by :user:`schlunma`
-  Avoid concatenation error in recipe_pcrglobwb.yml (:pull:`3645`) by :user:`bouweandela`
-  Update `scipy.integrate.simps` import (:pull:`3704`) by :user:`ehogan`

Improvements
~~~~~~~~~~~~

-  Add native6, OBS6 and RAWOBS rootpaths to metoffice template in config-user-example.yml and remove temporary dir (:pull:`3613`) by :user:`alistairsellar`

.. _changelog-v2-10-0:

v2.10.0
-------
Highlights

-  Add a realistic IPCC example recipe that reproduces figure 9.3 from AR6. It
   computes the mean sea-surface temperature anomaly between 1850-2100 over all
   available CMIP6 models. See the :ref:`recipe documentation <recipe_examples>`
   or read the `blog post <https://blog.esciencecenter.nl/easy-ipcc-powered-by-esmvalcore-19a0b6366ea7>`__
   for more information.

-  Added more plot types to monitoring diagnostic: Hovmoeller Z vs. time,
   Hovmoeller time vs latlon, variable vs. latitude are now available. See the
   :ref:`recipe documentation <recipe_monitor>` for more information.

-  Add support for 4 new datasets:

   - NOAA-CIRES-20CR v3 reanalysis
   - NASA MERRA reanalysis
   - NOAA marine boundary layer data for CH4
   - MOBO-DIC2004-2019

   See :ref:`supported_datasets` and :ref:`inputdata_observations` for more
   information.

-  Many recipes now have up-to-date obs4MIPs dataset names so required data can
   automatically be downloaded from ESGF.

This release includes

Bug fixes
~~~~~~~~~

-  Update recipe shapeselect to work with shapely v2 (:pull:`3283`) :user:`lukruh`
-  Correctly handle ``~`` when reading ``plot_folder`` option of monitoring diagnostic (:pull:`3449`) :user:`schlunma`
-  Fixed provenance tracking for NCL multipanel PNGs (:pull:`3332`) :user:`schlunma`
-  Fixed plot paths in NCL provenance tracking (:pull:`3422`) :user:`schlunma`
-  Fix erroneous file_type handling in certain NCL diagnostics (:pull:`3474`) :user:`zklaus`
-  Fix NCL provenance tracking (:pull:`3477`) :user:`schlunma`
-  Fix plots and provenance in Russell diagnostics (:pull:`3479`) :user:`schlunma`

Documentation
~~~~~~~~~~~~~

-  Add merge instructions to release instructions (:pull:`3292`) :user:`remi-kazeroni`
-  Update release schedule after release of v2.9.0 (:pull:`3289`) :user:`remi-kazeroni`
-  Add list of failing recipes for v2.9.0 release (:pull:`3294`) :user:`remi-kazeroni`
-  Update ``mamba`` version in readthedocs configuration docs builds (:pull:`3310`) :user:`valeriupredoi`
-  Add Romain Beucher to citation file as contributor (:pull:`3318`) :user:`valeriupredoi`
-  Removed recipe_carvalhais14nat from list of broken recipes (:pull:`3319`) :user:`remi-kazeroni`
-  Add `OBS-maintainers <https://github.com/orgs/ESMValGroup/teams/obs-maintainers>`__ team to documentation on OBS data maintenance and CMORizer reviews (:pull:`3335`) :user:`remi-kazeroni`
-  Add Pauline Bonnet to citation file (:pull:`3347`) :user:`Paulinebonnet111`
-  Ensure compatible zstandard and zstd in readthedocs builds (:pull:`3362`) :user:`zklaus`
-  Fix documentation build (:pull:`3397`) :user:`bouweandela`
-  Minor updates to release tools (:pull:`3216`) :user:`bouweandela`
-  Enhance provenance documentation (:pull:`3305`) :user:`alistairsellar`
-  Re-add communities and grants in zenodo file (:pull:`3416`) :user:`valeriupredoi`
-  Update Anconda badge in README (:pull:`3375`, :pull:`3453`) :user:`valeriupredoi`

Diagnostics
~~~~~~~~~~~

-  Slight refactoring of diagnostic script ``galytska23/select_variables_for_tigramite.py`` for generality and portability (:pull:`3298`) :user:`valeriupredoi` and :user:`egalytska`
-  Allow custom variable grouping in diagnostic script ``monitor/multi_datasets.py`` (:pull:`3343`) :user:`schlunma`
-  Extended monitor diagnostic with plot type variable vs. latitude (:pull:`3340`) :user:`ellensarauer`
-  Add Hovmoeller Z vs. time plot to monitoring diagnostic (:pull:`3345`) :user:`cubeme` and :user:`helgehr`
-  Adding Hovmoeller time vs latlon plots to monitoring recipes (:pull:`3341`) :user:`lukruh` and :user:`jeremykraftdlr`
-  Implied heat transport new diagnostic (:pull:`3177`) :user:`mo-abodas`
-  Recipe changes for new statistics preprocessors (percentiles) (:pull:`3351`) :user:`schlunma`
-  Add a realistic example recipe (:pull:`3356`) :user:`Peter9191` and :user:`bouweandela`
-  Support ``CenteredNorm`` in diagnostic monitor/multidatasets.py (:pull:`3415`) :user:`schlunma`
-  Use new preprocessor statistics calling convention for recipe_easy_ipcc.yml (:pull:`3418`) :user:`bouweandela`
-  Adapt to changed style scheme name in matplotlib (:pull:`3475`) :user:`zklaus`
-  Add version to dataset in python example recipe to avoid "Unknown file format" issue on JASMIN (:pull:`3322`) :user:`ehogan`
-  Add the dataset version in the heatwaves_coldwaves recipe to avoid the "Unknown file format" issue on JASMIN (:pull:`3373`) :user:`ehogan`

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Cmorizer for NOAA-CIRES-20CR v3 reanalysis (clt, clwvi, hus, prw, rlut, rlutcs, rsut, rsutcs) (:pull:`3137`) :user:`LisaBock`
-  CMORizer for NASA MERRA reanalysis (:pull:`3039`) :user:`axel-lauer`
-  Download and formatting of NOAA marine boundary layer data for CH4 (NOAA-MBL-CH4) (:pull:`3301`) :user:`FranziskaWinterstein`
-  Added CMORizer for MOBO-DIC2004-2019 (:pull:`3297`) :user:`schlunma`
-  Update obs4MIPs dataset names in quantilebias recipe (:pull:`3330`) :user:`rbeucher`
-  Update obs4MIPs dataset names in Schlund20esd recipe (:pull:`3329`) :user:`rbeucher`
-  Update obs4MIPs dataset names in flatoipcc recipes (:pull:`3328`) :user:`rbeucher`
-  Update obs4mips dataset names in clouds recipes (:pull:`3326`) :user:`rbeucher`
-  Update Obs4MIPs dataset names in ECS recipes (:pull:`3327`) :user:`rbeucher`
-  Update obs4mips dataset names in Bock et al recipes (:pull:`3324`, :pull:`3389` and :pull:`3473`) :user:`rbeucher` and :user:`bouweandela`
-  Update obs4mips dataset names in radiation budget recipe (:pull:`3323`) :user:`rbeucher`
-  Update Obs4MIPs dataset names in perfmetrics CMIP5 recipe (:pull:`3325`) :user:`rbeucher`

Automatic testing
~~~~~~~~~~~~~~~~~

-  Made sklearn test backwards-compatible with sklearn < 1.3 (:pull:`3285`) :user:`schlunma`
-  Update conda lock creation Github Action workflow and ship updated conda-lock file (:pull:`3307`, :pull:`3407`) :user:`valeriupredoi`
-  Compress all bash shell setters into one default option per GitHub Action workflow (:pull:`3315`) :user:`valeriupredoi`
-  Remove deprecated option ``offline`` from CI configuration (:pull:`3367`) :user:`schlunma`

Installation
~~~~~~~~~~~~

-  Use ESMValCore v2.10 (:pull:`3486`) :user:`bouweandela`

Improvements
~~~~~~~~~~~~

-  Merge v2.9.x into main  (:pull:`3286`) :user:`schlunma`
-  Allow NCL unit conversion `kg s-1` -> `GtC y-1` (:pull:`3300`) :user:`schlunma`

.. _changelog-v2-9-0:

v2.9.0
------

Highlights
~~~~~~~~~~

-  A new :ref:`diagnostic <api.esmvaltool.diag_scripts.seaborn_diag>` has been
   added to provide a high-level interface to
   `seaborn <https://seaborn.pydata.org/>`__,
   a Python data visualization library based on
   `matplotlib <https://matplotlib.org/>`__.
   See the :ref:`recipe documentation <recipes_seaborn_diag>` for more
   information.

-  We have included a new recipe and diagnostic that represent the major
   physical processes that describe Arctic-midlatitude teleconnections and
   provide the basis for the CMIP6 model evaluation for the further application
   of causal discovery.
   The results are discussed in the article
   `"Causal model evaluation of Arctic-midlatitude teleconnections in CMIP6" <https://doi.org/10.1002/essoar.10512569.1>`__
   by Galytska et al. (in review in Journal of Geophysical Research: Atmospheres).

-  It is now possible to use the
   `Dask distributed scheduler <https://docs.dask.org/en/latest/deploying.html>`__,
   which can
   `significantly reduce the run-time of recipes <https://github.com/ESMValGroup/ESMValCore/pull/2049#pullrequestreview-1446279391>`__.
   Configuration examples and advice are available in the
   :ref:`ESMValCore documentation <esmvalcore:config-dask>`.
   If configured, the Dask distributed scheduler will also be used by diagnostic
   scripts written in Python, so make sure to use
   `lazy data <https://scitools-iris.readthedocs.io/en/latest/userguide/real_and_lazy_data.html#real-and-lazy-data>`__
   wherever it is possible in your (new) diagnostics.
   More work on improving the computational performance is planned, so please
   share your experiences, good and bad, with this new feature in
   `ESMValGroup/ESMValCore#1763 <https://github.com/ESMValGroup/ESMValCore/discussions/1763>`__.

This release includes

Bug fixes
~~~~~~~~~

-  Fixed usage of ``work_dir`` in some CMORizer scripts (:pull:`3192`) :user:`remi-kazeroni`
-  Realize data for scalar cube in `recipe_carvalhais14nat` to avert issue from dask latest (2023.6.0) (:pull:`3265`) :user:`valeriupredoi`
-  Fix failing ``mlr`` diagnostic test by adding new scikit-learn default tag (:pull:`3273`) :user:`remi-kazeroni`
-  Fix ordering of models in perfmetrics diagnostic script (:pull:`3275`) :user:`LisaBock`

Documentation
~~~~~~~~~~~~~

-  Update release schedule after v2.8.0 (:pull:`3138`) :user:`remi-kazeroni`
-  Added reference entry for Winterstein (:pull:`3154`) :user:`FranziskaWinterstein`
-  Show logo on PyPI (:pull:`3185`) :user:`valeriupredoi`
-  Add Release Managers for v2.9.0 and v2.10.0 (:pull:`3184`) :user:`remi-kazeroni`
-  Fix readthedocs build with esmpy>=8.4.0 and missing ESMFMKFILE variable (:pull:`3205`) :user:`valeriupredoi`
-  Add ESMValCore release v2.8.1 into the documentation (:pull:`3235`) :user:`remi-kazeroni`
-  Modified links to the tutorial (:pull:`3236`) :user:`remi-kazeroni`
-  Fix gitter badge in README (:pull:`3258`) :user:`remi-kazeroni`
-  Add release notes for v2.9.0 (:pull:`3266`) :user:`bouweandela`

Diagnostics
~~~~~~~~~~~

-  New plot_type 1d_profile in monitor  (:pull:`3178`) :user:`FranziskaWinterstein`
-  Add Seaborn diagnostic (:pull:`3155`) :user:`schlunma`
-  New recipe and diagnostic for Arctic-midlatitude research (:pull:`3021`) :user:`egalytska`
-  Generate climatology on the fly for AutoAssess soil moisture (:pull:`3197`) :user:`alistairsellar`
-  Remove "fx_variables" from recipe_tebaldi21esd.yml (:pull:`3211`) :user:`hb326`
-  Remove "fx_variables" from ipccwg1ar5ch9 recipes (:pull:`3215`) :user:`katjaweigel`
-  Remove "fx_variables" from recipe_wenzel14jgr.yml (:pull:`3212`) :user:`hb326`
-  Update obs4MIPs dataset to the current naming scheme in recipe_smpi.yml (:pull:`2991`) :user:`bouweandela`
-  Fixed pandas diagnostics for pandas>=2.0.0 (:pull:`3209`) :user:`schlunma`
-  Update recipe_impact.yml to work with newer versions of `pandas` (:pull:`3220`) :user:`bouweandela`
-  Add variable long names to provenance record in monitoring diagnostics (:pull:`3222`) :user:`bsolino`

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Add CMORizer for GPCP-SG (pr) (:pull:`3150`) :user:`FranziskaWinterstein`
-  Extension of NASA MERRA2 CMORizer (cl, cli, clivi, clw, clwvi) (:pull:`3167`) :user:`axel-lauer`

Automatic testing
~~~~~~~~~~~~~~~~~

-  Add a CircleCI-testing-specific ``recipe_python_for_CI.yml`` to avoid calling geolocator/Nominatim over CI (:pull:`3159`) :user:`valeriupredoi`
-  Check if Python minor version changed after Julia install in development installation test (:pull:`3213`) :user:`valeriupredoi`
-  Fix tests using deprecated ``esmvalcore._config`` module that has been removed in ESMValCore v2.9 (:pull:`3204`) :user:`valeriupredoi`

Installation
~~~~~~~~~~~~

-  Add support for Python=3.11 (:pull:`3173`) :user:`valeriupredoi`
-  Drop python=3.8 support (:pull:`3193`) :user:`valeriupredoi`
-  Repair generation of conda lock files (:pull:`3148`) :user:`valeriupredoi`
-  Modernize lock creation script and repair lock generation (:pull:`3174`) :user:`valeriupredoi`
-  Pin numpy !=1.24.3 due to severe masking bug (:pull:`3182`) :user:`valeriupredoi`
-  Update xesmf to versions >= 0.4.0 (:pull:`2728`) :user:`zklaus`
-  Update esmpy import for ESMF version 8.4.0 or larger (:pull:`3188`) :user:`valeriupredoi`
-  Relax the pin on iris to allow the use of older versions for performance reasons (:pull:`3270`) :user:`bouweandela`
-  Use ESMValCore v2.9.0 (:pull:`3274`) :user:`bouweandela`

Improvements
~~~~~~~~~~~~

-  Update pre-commit hooks (:pull:`3189`) :user:`bouweandela`
-  Add support for using a dask distributed scheduler (:pull:`3151`) :user:`bouweandela`

.. _changelog-v2-8-0:

v2.8.0
------

Highlights
~~~~~~~~~~

-  This release includes the diagnostics for reproducing figures 3.9, 3.19,
   3.42 and 3.43 of the IPCC AR6 WG1 report.
   See :ref:`recipe documentation <recipes_ipccwg1ar6ch3>` about added recipes.
-  A new set of recipes and diagnostics has been included to evaluate cloud
   climatologies from CMIP models as used in `Lauer et al. (2023), J. Climate
   <https://doi.org/10.1175/JCLI-D-22-0181.1>`__.
   See :ref:`recipe documentation <recipes_clouds>` about added recipes.
-  Addition of a set of recipes for extreme events, regional and impact
   evaluation as used in `Weigel et al. (2021), GMD
   <https://doi.org/10.5194/gmd-14-3159-2021>`__ and in IPCC AR5.
   See :ref:`recipe documentation <recipes_ipccwg1ar5ch9>` about added recipes.

Highlights from ESMValCore v2.8.0 :ref:`here<esmvalcore:changelog-v2-8-0>`:

- ESMValCore now supports wildcards in recipes and offers improved support
  for ancillary variables and dataset versioning.
- Support for CORDEX datasets in a rotated pole coordinate system has been added.
- Native :ref:`ICON <esmvalcore:read_icon>` output is now made UGRID-compliant
  on-the-fly.
- The Python API has been extended with the addition of three modules:
  :mod:`esmvalcore.config`, :mod:`esmvalcore.dataset`, and
  :mod:`esmvalcore.local`
- The preprocessor :func:`~esmvalcore.preprocessor.multi_model_statistics`
  has been extended to support more use-cases.

This release includes:

Backwards incompatible changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please read the descriptions of the linked pull requests for detailed upgrade instructions.

-  Deprecated features scheduled for removal in v2.8.0 or earlier have now been removed
   (:pull:`2941`)
   :user:`schlunma`.
   Removed ``esmvaltool.iris_helpers.var_name_constraint`` (has been deprecated
   in v2.6.0; please use :class:`iris.NameConstraint` with the keyword argument
   ``var_name`` instead).
   Removed `write_netcdf` and `write_plots` from `recipe_filer.py`.
-  No files from the ``native6`` project will be found if a non-existent version
   of a dataset is specified (`#3041 <https://github.com/ESMValGroup/ESMValTool/pull/3041>`_)
   :user:`remi-kazeroni`.
   The tool now searches for exact ``version`` of ``native6`` datasets.
   Therefore, it is necessary to make sure that the version number in the
   directory tree matches with the version number in the recipe to find the files.
-  The conversion of precipitation units from monitoring diagnostic is now done
   at the preprocessor stage
   (`#3049 <https://github.com/ESMValGroup/ESMValTool/pull/3049>`_)
   :user:`schlunma`.
   To use the unit conversion for precipitation in the new version of this
   diagnostic, add it as a preprocessor for the precipitation dataset to the
   recipe.

Bug fixes
~~~~~~~~~

-  Fix for provenance records from `seaice_tsline.ncl` (:pull:`2938`) :user:`axel-lauer`
-  Fix in `validation.py` for resolving datasets with identical names by using distinct aliases (:pull:`2955`) :user:`FranziskaWinterstein`
-  Bugfix: masking of non-significant differences in `zonal.ncl` (perfmetrics) (:pull:`2957`) :user:`axel-lauer`
-  Fix typo in `perfmetrics/main.ncl` to add tropopause (:pull:`2966`) :user:`FranziskaWinterstein`
-  Fix .png bug in `wenzel16nat` diagnostics (:pull:`2976`) :user:`axel-lauer`
-  `Recipe_ocean_Landschuetzer2016`: Fix typo in filename to run model vs OBS diagnostics (:pull:`2997`) :user:`TomasTorsvik`
-  Fix read_cmor in NCL utilities (:pull:`3007`) :user:`axel-lauer`
-  Removed usages of deprecated features that cause diagnostic crashes (:pull:`3009`) :user:`schlunma`
-  Replace removed `matplotlib.pyplot.savefig` option `additional_artists` (:pull:`3075`) :user:`schlunma`
-  Added missing comma to `sommer17joss.bibtex` (:pull:`3078`) :user:`schlunma`
-  Fix call of output_type in `aux_plotting.ncl` (:pull:`3083`) :user:`LisaBock`
-  Remove colorbar from `bbox_extra_artists` (:pull:`3087`) :user:`schlunma`
-  Fix `MPI-ESM1-2-HR` entries in `recipe_tebaldi21esd` (:pull:`3093`) :user:`remi-kazeroni`
-  Fix bug in provenance writing of `perfmetrics` recipes v2.8.0 (:pull:`3098`) :user:`axel-lauer`
-  Fix `recipe_sea_surface_salinity` for v2.8 (:pull:`3102`) :user:`sloosvel`
-  Fix variable `short_name` and metadata for ESACCI-LST CMORizer (:pull:`3104`) :user:`remi-kazeroni`
-  Fix `recipe_carvalhais14`: replace outline patch with splines (:pull:`3111`) :user:`valeriupredoi`
-  Replace deprecated function `cm.register_cmap` with `mpl.colormaps.register` for `recipe_ arctic_ocean` (:pull:`3112`) :user:`TomasTorsvik`
-  Fix `recipe_extract_shape.yml` (lacking caption for provenance) (:pull:`3126`) :user:`valeriupredoi`

Community
~~~~~~~~~

-  Update documentation on pre-installed versions on HPC clusters (:pull:`2934`) :user:`remi-kazeroni`

Deprecations
~~~~~~~~~~~~

-  Remove radiation recipes that have been superseded by :ref:`recipe_radiation_budget <recipes_radiation_budget>` along with associated diagnostic scripts (`#3115 <https://github.com/ESMValGroup/ESMValTool/pull/3115>`_) :user:`alistairsellar`

Documentation
~~~~~~~~~~~~~

-  Backward compatibility policy (:pull:`2879`) :user:`alistairsellar`
-  Suppress installing and reinstalling dependencies with pip during readthedocs builds (:pull:`2913`) :user:`valeriupredoi`
-  Update installation instructions (:pull:`2939`) :user:`bouweandela`
-  Update documentation for `recipe_extreme_index` (:pull:`2951`) :user:`katjaweigel`
-  Update documentation and `recipe_check_obs` (ERA5) (:pull:`2952`) :user:`axel-lauer`
-  Updated ICON dataset entry in documentation (:pull:`2954`) :user:`schlunma`
-  Add Franziska Winterstein as collaborator in CITATION file (:pull:`3001`) :user:`valeriupredoi`
-  Update release schedule for v2.7.0 and v2.8.0 (:pull:`3010`) :user:`remi-kazeroni`
-  Add ESMValCore Bugfix release v2.7.1 to the release overview table (:pull:`3028`) :user:`valeriupredoi`
-  Detailed instructions for release procedure: running recipes and analyzing the output (:pull:`3032`) :user:`valeriupredoi`
-  Link backward compatibility policy to top level of ESMValCore changelog  (:pull:`3052`) :user:`alistairsellar`
-  Update release instructions (:pull:`3066`) :user:`remi-kazeroni`
-  Updated docs and tests regarding new `search_esgf` option (:pull:`3069`) :user:`schlunma`
-  Update script to draft release notes (:pull:`3070`) :user:`remi-kazeroni`
-  Synchronize documentation table of contents with ESMValCore (:pull:`3073`) :user:`bouweandela`
-  Update environment handling in release documentation (:pull:`3096`) :user:`remi-kazeroni`
-  Clarify use (or not) of Jasmin climatology files by soil moisture & permafrost recipes (:pull:`3103`) :user:`alistairsellar`
-  Add link to recipe portal in the gallery page (:pull:`3113`) :user:`remi-kazeroni`
-  Improve stratosphere documentation (:pull:`3114`) :user:`alistairsellar`
-  Added note to documentation that not all datasets used in `schlund20jgr` recipes are available on ESGF (:pull:`3121`) :user:`schlunma`
-  Draft changelog for `v2.8.0` (:pull:`3124`) :user:`remi-kazeroni`
-  Documenting broken recipes after recipe testing for releases (:pull:`3129`) :user:`remi-kazeroni`
-  Increase ESMValTool version to 2.8.0 and update release dates (:pull:`3136`) :user:`remi-kazeroni`

Diagnostics
~~~~~~~~~~~

-  Cloud diagnostics for Lauer et al. (2023) (:pull:`2750`) :user:`axel-lauer`
-  Splitting of `flato13ipcc.yml` into separate recipes and adding recipes for regional Figures (:pull:`2156`) :user:`katjaweigel`
-  Adding IPCC AR6 Chapter 3 Figure  3.43 - Pattern Correlation (:pull:`2772`) :user:`LisaBock`
-  Adding IPCC AR6 Chapter 3 Fig. 3.42 - Perfmetrics (:pull:`2856`) :user:`LisaBock`
-  Comment missing datasets and remove deprecated argument in `recipe_climate_change_hotspot` (:pull:`2920`) :user:`sloosvel`
-  Add plot type `annual_cycle` to multi-dataset monitoring diagnostic (:pull:`2922`) :user:`schlunma`
-  Adding IPCC AR6 Chapter 3 Fig. 3.19 - Speed-Up Of Zonal Mean Wind (:pull:`2984`) :user:`LisaBock`
-  Adding IPCC AR6 Chapter 3 Fig. 3.9 - Attribution (:pull:`2986`) :user:`LisaBock`
-  Obs4mips CERES-EBAF: update version to latest available through esgf in `recipe_validation.yml` (:pull:`3002`) :user:`valeriupredoi`
-  Improve flexibility of cloud diagnostics (:pull:`3016`) :user:`axel-lauer`
-  Let `recipe_impact.yml` write a CSV file that can directly be used in C4I portal (:pull:`2258`) :user:`Peter9192`
-  Fix version numbers of native6 datasets in recipes (`#3041`_) :user:`remi-kazeroni`
-  Removed automatic conversion of precipitation units from monitoring diagnostic (`#3049`_) :user:`schlunma`.
-  Updated recipes for ESMValCore v2.8 (:pull:`3064`) :user:`schlunma`
-  Fix `cos22esd` for release of 2.8 (:pull:`3097`) :user:`sloosvel`
-  Diagnostic for `recipe_autoassess_stratosphere.yml`: remove unused feature incompatible with Matplotlib=3.7.1 (:pull:`3089`) :user:`valeriupredoi`
-  Fix numpy deprecation in `hype` diagnostic (:pull:`3101`) :user:`Peter9192`
-  Remove superseded radiation recipes (`#3115`_) :user:`alistairsellar`
-  Removed `fx_variables` in `recipe_mpqb_xch4` and `recipe_lauer22jclim_fig8` (:pull:`3117`) :user:`axel-lauer`
-  Update Python example recipe (:pull:`3119`) :user:`bouweandela`
-  Updated figure settings to account for newer matplotlib version (:pull:`3133`) :user:`katjaweigel`

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Earth System Data Cube (ESDC) cmorizer (:pull:`2799`) :user:`bsolino`
-  Added CMORizer for Landschützer2020 (spco2) (:pull:`2908`) :user:`schlunma`
-  Added CMORizer for MOBO-DIC_MPIM (dissic) (:pull:`2909`) :user:`schlunma`
-  Added CMORizer for OceanSODA-ETHZ (areacello, co3os, dissicos, fgco2, phos, spco2, talkos) (:pull:`2915`) :user:`schlunma`
-  Extension of ERA-Interim CMORizer (cl, cli, clw, lwp, rlut, rlutcs, rsut, rsutcs) (:pull:`2923`) :user:`axel-lauer`
-  Add JRA-25 cmorizer (clt, hus, prw, rlut, rlutcs, rsut, rsutcs) (:pull:`2927`) :user:`LisaBock`
-  New CMORizers for datasets from the NCEP family (NCEP-DOE-R2, NCEP-NCAR-R1, NOAA-CIRES-20CR) (:pull:`2931`) :user:`hb326`
-  Updates to the recipes that use the NCEP reanalysis dataset (:pull:`2932`) :user:`hb326`
-  MERRA2 cmorizer convert vertical level coordinate units from hPa to Pa (:pull:`3003`) :user:`valeriupredoi`
-  MERRA2 cmorizer set UNLIMITED time coordinate (:pull:`3006`) :user:`valeriupredoi`
-  Added CMORizers for TCOM-CH4 (CH4) and TCOM-N2O (N2O) (:pull:`3014`) :user:`schlunma`
-  Update HadISST cmorizer to include recent years (:pull:`3027`) :user:`remi-kazeroni`

Automatic testing
~~~~~~~~~~~~~~~~~

-  Add DKRZ/Levante batch scripts for release recipe running (:pull:`2883`) :user:`valeriupredoi`
-  Remove `pytest-flake8` and call the use of `flake8` straight (:pull:`2904`) :user:`valeriupredoi`
-  Unpin `flake8` (:pull:`2937`) :user:`valeriupredoi`
-  Fix failing tests that use deprecated feature of `sklearn` (:pull:`2961`) :user:`schlunma`
-  Fix recipe loading tests for esmvalcore before and after version 2.8 (:pull:`3020`) :user:`valeriupredoi`
-  Update recipe load test for v2.8 (:pull:`3040`) :user:`bouweandela`
-  Test running recipes with the development version of ESMValCore (:pull:`3072`) :user:`bouweandela`
-  Fix `test_naming.py` so it doesn't let through directories that need be ignored (:pull:`3082`) :user:`valeriupredoi`
-  Conda environment files for interim use of `esmvalcore=2.8.0rc1` (:pull:`3090`) :user:`valeriupredoi`
-  Move `flake8` check to a step separate from installation on CircleCI (:pull:`3105`) :user:`bouweandela`
-  Recreate conda lock file to harpoon esmvalcore=2.8.0rc1 (:pull:`3108`) :user:`valeriupredoi`
-  Update batch script generation to run all recipes in one command (:pull:`3130`) :user:`remi-kazeroni`

Installation
~~~~~~~~~~~~

-  Merge release branch `release_270stable` in main so we pick up unsquashed commits and set the correct version 2.7.0 for main (and up version in CITATION.cff) (:pull:`2896`) :user:`valeriupredoi`
-  Unpin `NetCDF4` (:pull:`2929`) :user:`valeriupredoi`
-  Unpin `cf-units` (:pull:`2930`) :user:`bouweandela`
-  Set the version number on the development branches to one minor version more than the last release  (:pull:`2964`) :user:`bouweandela`
-  Pin `shapely<2.0.0` for linux64 (:pull:`2970`) :user:`valeriupredoi`
-  Unpin `matplotlib` (:pull:`3068`) :user:`valeriupredoi`
-  Add `packaging` as direct dependency to ESMValTool (:pull:`3099`) :user:`valeriupredoi`
-  Re-pin sphinx to latest (6.1.3) and add nbsphinx to the environment (:pull:`3118`) :user:`valeriupredoi`
-  Conda environment files for esmvalcore=2.8.0rc2 (:pull:`3120`) :user:`remi-kazeroni`
-  Remove rc (release candidates) conda channel and re-pin esmvalcore to new stable 2.8 (:pull:`3131`) :user:`valeriupredoi`

Improvements
~~~~~~~~~~~~

-  Read `config-user.yml` using `esmvalcore.config` module (:pull:`2736`) :user:`bouweandela`
-  Make results of recipes `schlund20jgr_*.yml` deterministic (:pull:`2900`) :user:`schlunma`
-  `Recipe_gier2020bg.yml`: add sorting to SA barplot (:pull:`2905`) :user:`bettina-gier`
-  Add the outline of a climatological tropopause to the zonalmean_profile plots (:pull:`2947`) :user:`FranziskaWinterstein`
-  Update data finder imports (:pull:`2958`) :user:`bouweandela`
-  Add support for the upcoming ESMValCore v2.8 release to the recipe filler tool (:pull:`2995`) :user:`bouweandela`
-  Updated monitoring diagnostics with netCDF output and additional logging (:pull:`3029`) :user:`schlunma`
-  Use aliases in perfmetrics (:pull:`3058`) :user:`FranziskaWinterstein`


.. _changelog-v2-7-0:

v2.7.0
------

Highlights
~~~~~~~~~~

-  This release has seen the inclusion of the code for figures 3.3, 3.4, 3.5, 3,13 and 3.15 of the IPCC AR6 WG1 report, see them in the `new documentation <https://docs.esmvaltool.org/en/latest/recipes/recipe_ipccwg1ar6ch3.html>`__
-  We have also included new diagnostics and recipe necessary to produce the plots and tables for the journal article "Climate model projections from the Scenario Model Intercomparison Project (ScenarioMIP) of CMIP6" by `Tebaldi et al. in ESD 2020-68 <https://doi.org/10.5194/esd-2020-68>`__ from 2021; also see the `recipe entry <https://docs.esmvaltool.org/en/latest/recipes/recipe_tebaldi21esd.html>`__
-  We have also extended the support for MERRA2 observational dataset, by adding support for a large number of variables, including 3D variables, see the `table of supported obs datasets <https://docs.esmvaltool.org/en/latest/input.html#supported-datasets-for-which-a-cmorizer-script-is-available>`__

Backwards incompatible changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Remove installation of R dependencies from the help message (:pull:`2761`) :user:`remi-kazeroni`

Bug fixes
~~~~~~~~~

-  Fix misplaced provenance records from IPCC AR5 Ch.12 diags (:pull:`2758`) :user:`axel-lauer`
-  Fix `esmvaltool.utils.testing.regression.compare` module to run with Python<3.10 too (:pull:`2778`) :user:`valeriupredoi`
-  Fixed small bug that could lead to wrong pr units in `monitor/multi_datasets.py` (:pull:`2788`) :user:`schlunma`
-  Pin `xgboost>1.6.1` so we avert documentation failing to build with `1.6.1` (:pull:`2780`) :user:`valeriupredoi`
-  Pin `matplotlib-base<3.6.0` to avoid conflict from `mapgenerator` that fails doc builds (:pull:`2830`) :user:`valeriupredoi`
-  Fixed wrong latitudes in NDP CMORizer (:pull:`2832`) :user:`schlunma`
-  Fix indexer in Autoassess supermeans module use a tuple of `(slice(), idx, idx)` (:pull:`2838`) :user:`valeriupredoi`
-  Replace xarray ufuncs with bogstandard numpy in weighting/climwip/calibrate_sigmas.py (:pull:`2848`) :user:`valeriupredoi`
-  Fix units MERRA2 CMORizer (:pull:`2850`) :user:`axel-lauer`
-  Fix bug when using log-scale y-axis for ocean transects. (:pull:`2862`) :user:`TomasTorsvik`

Community
~~~~~~~~~

-  Add MO-paths to config file (:pull:`2784`) `mo-tgeddes <https://github.com/mo-tgeddes>`__

Deprecations
~~~~~~~~~~~~

-  Recipe `recipe_esacci_oc.yml` replace with new regrid scheme `nearest_extrapolate` (:pull:`2841`) :user:`valeriupredoi`

Documentation
~~~~~~~~~~~~~

-  Update release schedule for v2.7 (:pull:`2747`) :user:`bouweandela`
-  Add Met Office installation method (:pull:`2751`) `mo-tgeddes <https://github.com/mo-tgeddes>`__
-  Add release dates for 2023 (:pull:`2769`) :user:`remi-kazeroni`
-  Made `maintainer` entry mandatory for published recipes (:pull:`2703`) :user:`schlunma`
-  Use command with current command line opts for `cffconvert` in documentation (:pull:`2791`) :user:`valeriupredoi`
-  Update CMORizer documentation with command options (:pull:`2795`) :user:`remi-kazeroni`
-  Fixed broken link for monthly meetings (:pull:`2806`) :user:`remi-kazeroni`
-  Update MO obs4MIPs paths in the user configuration file (:pull:`2813`) `mo-tgeddes <https://github.com/mo-tgeddes>`__
-  Fix Windows incompatible file names in documentation of recipe_climate_change_hotspot.yml (:pull:`2823`) :user:`ledm`
-  Update documentation for the Landschuetzer 2016 recipe. (:pull:`2801`) :user:`TomasTorsvik`
-  Fixed anaconda badge in README (:pull:`2866`) :user:`valeriupredoi`
-  Update release strategy notes (:pull:`2734`) :user:`sloosvel`
-  Add documentation on how to handle CMORizers for multiple dataset versions (:pull:`2730`) :user:`remi-kazeroni`
-  Extending documentation: recipe maintainer + broken recipe policy (:pull:`2719`) :user:`axel-lauer`

Diagnostics
~~~~~~~~~~~

-  Recipe and diagnostics for : Tebaldi et al.,ESD, 2021 (:pull:`2052`) `debe-kevin <https://github.com/debe-kevin>`__
-  Figures for IPCC AR6 WG1 Chapter 3 (Atmosphere) (:pull:`2533`) :user:`LisaBock`

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Update CERES-EBAF to Ed4.1 (:pull:`2752`) :user:`axel-lauer`
-  New CMORizer for CALIPSO-ICECLOUD (:pull:`2753`) :user:`axel-lauer`
-  New CMORizer for CLOUDSAT-L2 (:pull:`2754`) :user:`axel-lauer`
-  Update MERRA2 cmorizer with extra 2D and 3D variables (:pull:`2774`) :user:`valeriupredoi`

Automatic testing
~~~~~~~~~~~~~~~~~

-  Pin `netcdf4 != 1.6.1` since that is spitting large numbers of SegFaults (:pull:`2796`) :user:`valeriupredoi`

Installation
~~~~~~~~~~~~

-  Increase esmvalcore version to 2.7.0 in environment files (:pull:`2860`) :user:`valeriupredoi`
-  Add iris-esmf-regrid as a dependency (:pull:`2880`) :user:`zklaus`

Improvements
~~~~~~~~~~~~

-  Fix tebaldi21esd (:pull:`2749`) :user:`axel-lauer`
-  Added option to show basic statistics in plots of `monitor/multi_datasets.py` (:pull:`2790`) :user:`schlunma`
-  Remove retracted datasets from `recipe_climate_change_hotspot` (:pull:`2854`) :user:`sloosvel`


.. _changelog-v2-6-0:

v2.6.0
------

Highlights
~~~~~~~~~~

- A new monitoring diagnostic has been added to allow the comparison of model runs against reference datasets. For details, see :ref:`Monitoring diagnostic to show multiple datasets in one plot (incl. biases) <api.esmvaltool.diag_scripts.monitor.multi_datasets>`.
- A tool has been developed to compare the output of recipe runs against previous runs, in order to detect in an automated way breaking changes between releases. Find more information in :ref:`Comparing recipe runs <compare_recipe_runs>`.
- The recipe :ref:`Climate Change Hotspot <recipe_climate_change_hotspot.rst>` allows to compute hotspots in any rectangular region.

Please also note the highlights from the corresponding ESMValCore release :ref:`here<esmvalcore:changelog-v2-6-0>`.
Thanks to that ESMValTool has gained the following features:

- A new set of CMOR fixes is now available in order to load native EMAC model output and CMORize it on the fly.
- The version number of ESMValCore is now automatically generated using `setuptools_scm <https://github.com/pypa/setuptools_scm/#default-versioning-scheme>`__, which extracts Python package versions from git metadata.

This release includes

Bug fixes
~~~~~~~~~

-  Fix dtype for Marrmot recipe results (:pull:`2646`) :user:`SarahAlidoost`
-  Adapt test_fix_coords to new version of cf-units (:pull:`2707`) :user:`zklaus`
-  Fix nested axes in `recipe_martin18_grl` and `recipe_li17natcc` (:pull:`2712`) :user:`lukruh`
-  Update common_climdex_preprocessing_for_plots.R (:pull:`2727`) :user:`earnone`

Community
~~~~~~~~~

-  Collecting github user names for config-references (:pull:`2677`) :user:`lukruh`

Deprecations
~~~~~~~~~~~~

-  Deprecate the function `esmvaltool.diag_scripts.shared.var_name_constraint`. This function is scheduled for removal in v2.8.0. Please use :class:`iris.NameConstraint` with the keyword argument var_name instead: this is an exact replacement. (:pull:`2655`) :user:`schlunma`

Documentation
~~~~~~~~~~~~~

-  Documentation Improvements (:pull:`2580`) :user:`stacristo`
-  Fixed broken label in the documentation (:pull:`2616`) :user:`remi-kazeroni`
-  Add readthedocs configuration file (:pull:`2627`) :user:`bouweandela`
-  Update the command for building the documentation (:pull:`2622`) :user:`bouweandela`
-  Added DKRZ-Levante to `config-user-example.yml` (:pull:`2632`) :user:`remi-kazeroni`
-  Improved documentation on native dataset support (:pull:`2635`) :user:`schlunma`
-  Add documentation on building and uploading Docker images (:pull:`2662`) :user:`bouweandela`
-  Remove support for Mistral in `config-user-example.yml` (:pull:`2667`) :user:`remi-kazeroni`
-  Add note to clarify that CORDEX support is work in progress (:pull:`2682`) :user:`bouweandela`
-  Restore accidentally deleted text from input data docs (:pull:`2683`) :user:`bouweandela`
-  Add running settings note in `recipe_wenzel16nat.yml` documentation (:pull:`2692`) :user:`sloosvel`
-  Add a note on transferring permissions to the release manager (:pull:`2688`) :user:`bouweandela`
-  Update documentation on ESMValTool module at DKRZ (:pull:`2696`) :user:`remi-kazeroni`
-  Add note on how to run recipe_wenzel14jgr.yml (:pull:`2717`) :user:`sloosvel`
-  Added conda forge feedstock repo link in README (:pull:`2555`) :user:`valeriupredoi`

Diagnostics
~~~~~~~~~~~

-  Compute bias instead of correlation in `compare_salinity.py` (:pull:`2642`) :user:`sloosvel`
-  Update monitor diagnostics (:pull:`2608`) :user:`schlunma`
-  Add new Psyplot diagnostic (:pull:`2653`) :user:`schlunma`
-  Reduce memory usage of lisflood recipe (:pull:`2634`) :user:`sverhoeven`
-  Provenance in ocean diagnostics (:pull:`2651`) :user:`tomaslovato`
-  Extend monitor diagnostics with multi-dataset plots (:pull:`2657`) :user:`schlunma`
-  Recipe and diagnostics to plot climate change hotspots: Cos et al., ESD 2022 (:pull:`2614`) :user:`pepcos`
-  Update plots of consecutive dry days recipe (:pull:`2671`) :user:`bouweandela`
-  Fix the format of ids in Hype forcing files (:pull:`2679`) :user:`SarahAlidoost`
-  WFlow diagnostic script: remove manual rechunking (:pull:`2680`) :user:`Peter9192`

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Extending the HadCRUT5 cmorizer (:pull:`2509`) :user:`LisaBock`
-  Cmorize Kadow2020 dataset (:pull:`2513`) :user:`LisaBock`
-  Cmorize NOAAGlobalTemp dataset (:pull:`2515`) :user:`LisaBock`
-  Add option to CMORize ts as tos in ESACCI data (:pull:`2731`) :user:`sloosvel`

Automatic testing
~~~~~~~~~~~~~~~~~

-  Add a tool for comparing recipe runs to previous runs (:pull:`2613`) :user:`bouweandela`
-  Ignore NCL interface files when comparing recipe runs (:pull:`2673`) :user:`bouweandela`
-  Add a short version of recipe deangelis15nat for testing (:pull:`2685`) :user:`katjaweigel`
-  Expanded recipe output comparison tool to better handle absolute paths in output (:pull:`2709`) :user:`schlunma`
-  Update development infrastructure (:pull:`2663`) :user:`bouweandela`

Installation
~~~~~~~~~~~~

-  Removed `package/meta.yaml` and all references to it (:pull:`2612`) :user:`schlunma`

Improvements
~~~~~~~~~~~~

-  Improved handling of weights in MLR diagnostics (:pull:`2625`) :user:`schlunma`
-  Fixed order of variables in perfemetrics plot of Anav13jclim recipe (:pull:`2706`) :user:`schlunma`
-  Added input file sorting to many diagnostic to make output exactly reproducible (:pull:`2710`) :user:`schlunma`
-  Removed 'ancestors' attributes before saving netcdf files in emergent constraints diagnostics (:pull:`2713`) :user:`schlunma`

.. _changelog-v2-5-0:

v2.5.0
------

Highlights
~~~~~~~~~~

- A new recipe to plot generic preprocessor output is now available. For details, see :ref:`recipe_monitor`.
- The CMORization of observational and other datasets has been overhauled. For many datasets, an automatic download script is now available. For details, see :ref:`inputdata_observations` and :ref:`new-cmorizer`.

Please also note the highlights from the corresponding ESMValCore release :ref:`here<esmvalcore:changelog-v2-5-0>`.
Thanks to that ESMValTool has gained the following features:

- The new preprocessor ``extract_location`` can extract arbitrary locations on the Earth.
- Time ranges can now be extracted using the `ISO 8601 format <https://en.wikipedia.org/wiki/ISO_8601>`_.
- The new preprocessor ``ensemble_statistics`` can calculate arbitrary statistics over all ensemble members of a simulation.


This release includes

Backwards incompatible changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Streamline observations download (:pull:`1657`) `Javier Vegas-Regidor <https://github.com/jvegreg>`__. This change removes the ``cmorize_obs`` command which has previously been used to CMORize observations and other datasets. The new command ``esmvaltool data`` provides many new features apart from the CMORization (``esmvaltool data format``), for example, automatic downloading of observational datasets (``esmvaltool data download``). More details on this can be found :ref:`here<inputdata_observations>` and :ref:`here<new-cmorizer>`.
-  Dropped Python 3.7 (:pull:`2585`) :user:`schlunma`. ESMValTool v2.5.0 dropped support for Python 3.7. From now on Python >=3.8 is required to install ESMValTool. The main reason for this is that conda-forge dropped support for Python 3.7 for OSX and arm64 (more details are given `here <https://github.com/ESMValGroup/ESMValTool/issues/2584#issuecomment-1063853630>`__).

Bug fixes
~~~~~~~~~

-  Remove the use of `esmvalgroup` channel from the conda install Github Action workflow (:pull:`2420`) :user:`valeriupredoi`
-  Ignore .pymon-journal file in test discovery (:pull:`2491`) :user:`zklaus`
-  Relocate pytest-monitor outputted database `.pymon` so `.pymon-journal` file should not be looked for by `pytest` (:pull:`2501`) :user:`valeriupredoi`
-  Re-establish Python 3.7 compatibility (:pull:`2506`) :user:`zklaus`
-  Update intersphinx mapping (:pull:`2531`) :user:`zklaus`
-  Fixed `KeyError` in `recipe_ocean_bgc.yml` (:pull:`2540`) :user:`schlunma`
-  Corrected ESACCI-SEA-SURFACE-SALINITY from OBS to OBS6 (:pull:`2542`) :user:`axel-lauer`
-  Fixed `recipe_kcs.yml` (:pull:`2541`) :user:`schlunma`
-  Fix MDER diagnostic regression_stepwise (:pull:`2545`) :user:`axel-lauer`
-  Fix for recipe_wenzel16nat (:pull:`2547`) :user:`axel-lauer`
-  Fixed `recipe_carvalhais14nat` and removed deprecated use of np.float (:pull:`2558`) :user:`schlunma`
-  Fix `recipe_wenzel14jgr` (:pull:`2577`) :user:`remi-kazeroni`
-  Fixed various recipes by removing faulty or non-available datasets (:pull:`2563`) :user:`schlunma`
-  Remove missing CMIP5 data from 2 recipes (:pull:`2579`) :user:`remi-kazeroni`
-  Fix `recipe_seaice` (:pull:`2578`) :user:`remi-kazeroni`
-  Fix `recipe_climwip_brunner20esd` (:pull:`2581`) :user:`remi-kazeroni`

Deprecations
~~~~~~~~~~~~

-  Remove `--use-feature=2020-resolver` command line option for obsolete pip 2020 solver (:pull:`2493`) :user:`valeriupredoi`
-  Renamed vertical regridding schemes in affected recipes (:pull:`2487`) :user:`schlunma`

Documentation
~~~~~~~~~~~~~

-  Update release manager for v2.5 (:pull:`2429`) :user:`axel-lauer`
-  Mention ENES Climate Analytics service (:pull:`2438`) :user:`bouweandela`
-  Add recipe overview page (:pull:`2439`) :user:`bouweandela`
-  Fix pointer to Tutorial lesson on preprocessor from 05 to 06 (:pull:`2473`) :user:`valeriupredoi`
-  Removed obsolete option `synda-download` from documentation (:pull:`2485`) :user:`schlunma`
-  Update CMUG XCH4 docu figure (:pull:`2502`) :user:`axel-lauer`
-  Add Python=3.10 to package info, update Circle CI auto install and documentation for Python=3.10 (:pull:`2503`) :user:`schlunma`
-  Unify user configuration file (:pull:`2507`) :user:`schlunma`
-  Synchronized `config-user.yml` with version from ESMValCore (:pull:`2516`) :user:`schlunma`
-  CITATION.cff fix and automatic validation of your citation metadata (:pull:`2517`) :user:`abelsiqueira`
-  Add backwards incompatible changes at the top of the release notes draft (:pull:`2431`) :user:`bouweandela`
-  Fixed intersphinx mapping of `scipy` (:pull:`2523`) :user:`schlunma`
-  Add authors to citation cff (:pull:`2525`) :user:`SarahAlidoost`
-  Update documentation on running a recipe (:pull:`2432`) :user:`bouweandela`
-  Fix recipe `hydrology/recipe_wflow.yml` (:pull:`2549`) :user:`remi-kazeroni`
-  Update `draft_release_notes.py` for new release (:pull:`2553`) :user:`schlunma`
-  Added stand with Ukraine badge (:pull:`2565`) :user:`valeriupredoi`
-  Updated CREM docu (recipe_williams09climdyn.yml) (:pull:`2567`) :user:`axel-lauer`
-  First draft for v2.5.0 changelog (:pull:`2554`) :user:`schlunma`
-  Replace nonfunctional Github Actions badge with cool one in README (:pull:`2582`) :user:`valeriupredoi`
-  Updated changelog (:pull:`2589`) :user:`schlunma`
-  Updated release strategy with current release and upcoming release (:pull:`2597`) :user:`schlunma`
-  Increased ESMValTool version to 2.5.0 (:pull:`2600`) :user:`schlunma`

Diagnostics
~~~~~~~~~~~

-  AutoAssess: Add new diagnostic for radiation budget (:pull:`2282`) :user:`Jon-Lillis`
-  CMUG Sea Surface Salinity dataset and diagnostic (:pull:`1832`) `Javier Vegas-Regidor <https://github.com/jvegreg>`__
-  Recipe with new diagnostics for ESA-CMUG H2O (:pull:`1834`) :user:`katjaweigel`
-  Cleaned Schlund et al. (2020) recipe and fixed small bugs in corresponding diagnostic (:pull:`2484`) :user:`schlunma`
-  Add ESA CCI LST cmorizer and diagnostic (:pull:`1897`) :user:`morobking`
-  XCH4 ESA CMUG diagnostics (subset of the MPQB diagnostics) (:pull:`1960`) :user:`hb326`
-  Add support for ESACCI Ocean Color (Chlorophyll) observations (:pull:`2055`) `ulrikaw-cloud <https://github.com/ulrikaw-cloud>`__
-  Updated `recipe_zmnam.yml` with hemisphere selection (:pull:`2230`) :user:`fserva`
-  Add recipe and diagnostic scripts to compute figures of D9.4 of ISENES3 (:pull:`2441`) :user:`sloosvel`
-  Save resampled climates from KCS diagnostic local_resampling.py (:pull:`2221`) :user:`Emmadd`
-  Use years from KCS recipe (:pull:`2223`) :user:`Emmadd`
-  Recipe to plot generic output from the preprocessor (:pull:`2184`) `Javier Vegas-Regidor <https://github.com/jvegreg>`__
-  Fixed provenance tracking for emergent constraint diagnostics (:pull:`2573`) :user:`schlunma`

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Ensure dummy data for cmorize_obs_woa test are written to the correct directory (:pull:`2451`) :user:`ehogan`
-  Add ESA CCI LST cmorizer and diagnostic (see previous section `Diagnostics`)

Automatic testing
~~~~~~~~~~~~~~~~~

-  Run a nightly Github Actions workflow to monitor tests memory per test (configurable for other metrics too) and lists the slowest 100 tests (:pull:`2449`) :user:`valeriupredoi`
-  Fix individual pytest runs broken due to missing explicit imports from `iris` and adding a couple missing package markers (:pull:`2455`) :user:`valeriupredoi`
-  Add Python=3.10 to Github Actions and switch to Python=3.10 for the Github Action that builds the PyPi package (:pull:`2488`) :user:`valeriupredoi`
-  Switch all github actions from miniconda to mambaforge (:pull:`2498`) :user:`zklaus`
-  Pin `flake8<4` to have actual FLAKE8 error printed if tests fail and not garbage (:pull:`2492`) :user:`valeriupredoi`
-  Implementing conda lock (:pull:`2193`) :user:`valeriupredoi`
-  [Docker] Update Docker container builds with correct installations of Julia (:pull:`2530`) :user:`valeriupredoi`
- Update Linux condalock file (various pull requests) github-actions[bot]

Installation
~~~~~~~~~~~~

-  Comment out release candidate channel in environment.yml (:pull:`2417`) :user:`zklaus`
-  Comment out rc channel in osx environment file (:pull:`2421`) :user:`valeriupredoi`
-  Add `python-cdo` as conda-forge dependency in environment files to ensure `cdo` gets used from conda-forge and not pip (:pull:`2469`) :user:`valeriupredoi`
-  Install rasterio from conda-forge and avoid issues from python=3.10 (:pull:`2479`) :user:`valeriupredoi`
-  Updated dependencies with new ESMValCore version (:pull:`2599`) :user:`schlunma`

Improvements
~~~~~~~~~~~~

-  Remove use of OBS and use CMIP instead in `examples/recipe_ncl.yml` (:pull:`2494`) :user:`valeriupredoi`
-  Expanded `recipe_preprocessor_test.yml` to account for new `multi_model_statistics` features (:pull:`2519`) :user:`schlunma`
-  Updated piControl periods for recipes that use KACE-1-0-G (:pull:`2537`) :user:`schlunma`
-  Reduced time range in `recipe_globwat.yml` (:pull:`2548`) :user:`schlunma`
-  Removed models with missing data from recipe_williams09climdyn.yml (:pull:`2566`) :user:`axel-lauer`
-  Restored original versions of `recipe_schlund20esd.yml` and `recipe_meehl20sciadv.yml` (:pull:`2583`) :user:`schlunma`


.. _changelog-v2-4-0:

v2.4.0
------

Highlights
~~~~~~~~~~

- ESMValTool is moving from Conda to Mamba as the preferred installation method. This will speed up the
  installation and comes with some improvements behind the scenes.
  Read more about it at :ref:`Move to Mamba<move-to-mamba>` and in :ref:`the installation guide<install>`.

Please also note the highlights from the corresponding ESMValCore release :ref:`here<esmvalcore:changelog-v2-4-0>`.
Thanks to that ESMValTool has gained the following features:

- Download any missing data that is available on the ESGF automatically.
- Resume previous runs, reusing expensive pre-processing results.


This release includes

Bug fixes
~~~~~~~~~

-  Fixed `recipe_meehl20sciadv.yml` for ESMValCore 2.3 (:pull:`2253`) :user:`schlunma`
-  Fix provenance of NCL figures created using the log_provenance function (:pull:`2279`) :user:`bouweandela`
-  Fix bug in ClimWIP brunner19 recipe when plotting (:pull:`2226`) :user:`lukasbrunner`
-  Pin docutils <0.17 to fix sphinx build with rtd theme (:pull:`2312`) :user:`zklaus`
-  Fix example recipes (:pull:`2338`) :user:`valeriupredoi`
-  Do not add bounds to plev (plev19) in era interim cmorizer (:pull:`2328`) :user:`valeriupredoi`
-  Fix problem with pip 21.3 that prevents installation from source (:pull:`2344`) :user:`zklaus`
-  Add title to recipe embedded in test_diagnostic_run.py (:pull:`2353`) :user:`zklaus`
-  Fix capitalization of obs4MIPs (:pull:`2368`) :user:`bouweandela`
-  Specify that areacella is needed for area statistics in the Python example recipe (:pull:`2371`) :user:`bouweandela`
-  Enabling variable `obs550lt1aer` in recipes (:pull:`2388`) :user:`remi-kazeroni`
-  Update a diagnostic to new Iris version (:pull:`2390`) :user:`katjaweigel`
-  Fixed bug in provenance tracking of ecs_scatter.ncl (:pull:`2391`) :user:`schlunma`
-  Fix provenance issue in pv_capacity_factor.R (:pull:`2392`) :user:`katjaweigel`
-  Remove obsolete write_plots option from R diagnostics (:pull:`2395`) :user:`zklaus`
-  Fix arctic ocean diagnostic (:pull:`2397`) :user:`zklaus`
-  Fix sea ice drift recipe and script (:pull:`2404`) :user:`sloosvel`
-  Adapt diagnostic script to new version of iris (:pull:`2403`) :user:`zklaus`
-  Fix ocean multimap (:pull:`2406`) :user:`zklaus`
-  Fix diagnostic that uses `xarray`: `dtype` correctly set and harmonize `xarray` and `matplotlib` (:pull:`2409`) :user:`zklaus`
-  Deactivate provenance logging for plots in thermodyn toolbox (:pull:`2414`) :user:`zklaus`

Deprecations
~~~~~~~~~~~~

-  Removed write_plots and write_netcdf from some NCL diagnostics (:pull:`2293`) :user:`schlunma`
-  Fixed provenance logging of all python diagnostics by removing 'plot_file' entry (:pull:`2296`) :user:`schlunma`
-  Do not deprecate classes Variable, Variables and Datasets on a specific version (:pull:`2286`) :user:`schlunma`
-  Remove obsolete write_netcdf option from ncl diagnostic scripts (:pull:`2387`) :user:`zklaus`
-  Remove write plots from ocean diagnostics (:pull:`2393`) :user:`valeriupredoi`
-  More removals of instances of `write_plots` from Python diagnostics (appears to be the final removal from Py diags) (:pull:`2394`) :user:`valeriupredoi`

Documentation
~~~~~~~~~~~~~

-  List Manuel Schlund as release manager for v2.5 (:pull:`2268`) :user:`bouweandela`
-  GlobWat fix download links and gdal command (:pull:`2334`) :user:`babdollahi`
-  Add titles to recipes authored by `predoi_valeriu` (:pull:`2333`) :user:`valeriupredoi`
-  Added titles to recipes maintained by lauer_axel (:pull:`2332`) :user:`axel-lauer`
-  Update the documentation of the GRACE CMORizer (:pull:`2349`) :user:`remi-kazeroni`
-  Add titles in BSC recipes (:pull:`2351`) :user:`sloosvel`
-  Update esmvalcore dependency to 2.4.0rc1 (:pull:`2348`) :user:`zklaus`
-  Add titles to recipes maintained by Peter Kalverla (:pull:`2356`) :user:`Peter9192`
-  Adding titles to the recipes with maintainer hb326 (:pull:`2358`) :user:`hb326`
-  Add title for zmnam as for #2354 (:pull:`2363`) :user:`fserva`
-  Added recipe titles the the ocean recipes.  (:pull:`2364`) :user:`ledm`
-  Update recipe_thermodyn_diagtool.yml - add title (:pull:`2365`) :user:`ValerioLembo`
-  Fix provenance of figures of several R diagnostics (:pull:`2300`) :user:`bouweandela`
-  Adding titles to Mattia's recipes (:pull:`2367`) :user:`remi-kazeroni`
-  Adding titles to wenzel recipes (:pull:`2366`) :user:`hb326`
-  Fix formatting of some recipe titles merged from PR 2364 (:pull:`2372`) :user:`zklaus`
-  Adding titles to Bjoern's recipes (:pull:`2369`) :user:`remi-kazeroni`
-  Add titles to ocean recipes (maintainer Lovato) (:pull:`2375`) :user:`tomaslovato`
-  Add titles for three c3s-magic recipes (:pull:`2378`) :user:`zklaus`
-  Add title for recipe maintained by Ruth Lorenz (:pull:`2379`) :user:`zklaus`
-  Fix toymodel recipe (:pull:`2381`) Javier Vegas-Regidor
-  Added titles for recipes of maintainer `schlund_manuel` (:pull:`2377`) :user:`schlunma`
-  Write_plots and titles for deangelis15nat, li17natcc, martin18grl, pv_capacity_factor (:pull:`2382`) :user:`katjaweigel`
-  Add titles for some recipes (:pull:`2383`) :user:`zklaus`
-  Adding titles for recipes by von Hardenberg and Arnone (:pull:`2384`) :user:`zklaus`
-  Last two missing titles (:pull:`2386`) :user:`valeriupredoi`
-  Update documentation on downloading data (:pull:`2370`) :user:`bouweandela`
-  Fix installation instructions for Julia (:pull:`2335`) :user:`zklaus`
-  Fix provenance of Julia example diagnostic (:pull:`2289`) :user:`bouweandela`
-  Added notes on use of mamba in the installation documentation chapter (:pull:`2236`) :user:`valeriupredoi`
-  Update version number for 2.4.0 release (:pull:`2410`) :user:`zklaus`
-  Update release schedule for 2.4.0 (:pull:`2412`) :user:`zklaus`
-  Update changelog for 2.4.0 release (:pull:`2411`) :user:`zklaus`

Diagnostics
~~~~~~~~~~~

-  Add all available CMIP5 and CMIP6 models to recipe_impact.yml (:pull:`2251`) :user:`bouweandela`
-  Add Fig. 6, 7 and 9 of Bock20jgr (:pull:`2252`) :user:`LisaBock`
-  Generalize `recipe_validation*` diagnostic to work with identical control and experiment dataset names (:pull:`2284`) :user:`valeriupredoi`
-  Add missing preprocessor to recipe_gier2020bg and adapt to available data (:pull:`2399`) :user:`bettina-gier`
-  Removed custom version of `AtmosphereSigmaFactory` in diagnostics (:pull:`2405`) :user:`schlunma`

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Replace recipe_era5.yml with recipe_daily_era5.yml (:pull:`2182`) :user:`SarahAlidoost`
-  Update WOA cmorizer for WOA18 and WOA13v2 (:pull:`1812`) :user:`LisaBock`
-  GLODAP v2.2016 ocean data cmorizer (:pull:`2185`) :user:`tomaslovato`
-  Updated GCP CMORizer (:pull:`2295`) :user:`schlunma`

Automatic testing
~~~~~~~~~~~~~~~~~

-  Add a cylc suite to run all recipes (:pull:`2219`) :user:`bouweandela`
-  Retire test with Python 3.6 from full development Github Actions test (:pull:`2229`) :user:`valeriupredoi`
-  Remove Python 3.6 tests from GitHub Actions (:pull:`2264`) :user:`valeriupredoi`
-  Unpin upper bound for iris (previously was at <3.0.4) (:pull:`2266`) :user:`valeriupredoi`
-  Pin latest esmvalcore to allow use of the bugfix release 2.3.1 always (:pull:`2269`) :user:`valeriupredoi`
-  Add apt update so Julia gets found and installed by Docker (:pull:`2290`) :user:`valeriupredoi`
-  Use mamba for environment update and creation in the Docker container build on DockerHub (:pull:`2297`) :user:`valeriupredoi`
-  Docker container experimental - run a full env solve with mamba instead of a conda update (:pull:`2306`) :user:`valeriupredoi`
-  Full use of mamba in Github Actions source install test and use generic Python 3.7 (removing the very specific 3.7.10) (:pull:`2287`) :user:`valeriupredoi`
-  Replace use of conda with mamba for conda_install test on Circle CI (:pull:`2237`) :user:`valeriupredoi`
-  Update circleci configuration (:pull:`2357`) :user:`zklaus`

Installation
~~~~~~~~~~~~

-  Remove `mpich` from conda dependencies list (:pull:`2343`) :user:`valeriupredoi`

Improvements
~~~~~~~~~~~~

-  Add script for extracting a list of input files from the provenance (:pull:`2278`) :user:`bouweandela`
-  Update github actions (:pull:`2360`) :user:`zklaus`
-  Removed 'write_plots' from all NCL diagnostics (:pull:`2331`) :user:`axel-lauer`
-  Update and modernize `config-user-example.yml` (:pull:`2374`) :user:`valeriupredoi`


.. _changelog-v2-3-0:

v2.3.0
------

This release includes

Bug fixes
~~~~~~~~~

-  Indent block to pick up and raise exception if cmorizer data not found (TierX dir is not there) (:pull:`1877`) :user:`valeriupredoi`
-  Skip recipe filler tests until we have a new release since GA tests are failing (:pull:`2089`) :user:`valeriupredoi`
-  Fixed broken link to contributions in README (:pull:`2102`) :user:`schlunma`
-  Fix recipe filler for the case the variable doesn't contain short_name (:pull:`2104`) :user:`valeriupredoi`
-  Add fix for iris longitude bug to ClimWIP (:pull:`2107`) :user:`lukasbrunner`
-  Update for outdated link to reference Déandreis et al. (2014). (:pull:`2076`) :user:`katjaweigel`
-  Fixed recipes for ESMValCore 2.3.0 (:pull:`2203`) :user:`schlunma`
-  Fix the WFDE5 cmorizer (:pull:`2211`) :user:`remi-kazeroni`
-  Fix broken CMORizer log message if no Tier directory exists (:pull:`2207`) :user:`jmrgonza`
-  Fix bug in ClimWIP basic test recipe when plotting (:pull:`2225`) :user:`lukasbrunner`
-  Fix bug in ClimWIP advanced test recipe when plotting (:pull:`2227`) :user:`lukasbrunner`
-  Adjust time range for the `WDFE5` dataset in the `recipe_check_obs.yml` (:pull:`2232`) :user:`remi-kazeroni`
-  Fix plot and provenance of recipe_consecdrydays (:pull:`2244`) :user:`bouweandela`

Documentation
~~~~~~~~~~~~~

-  Improving the README.md file with a more appealing look and bit more info (:pull:`2065`) :user:`valeriupredoi`
-  Update plot title martin18grl (:pull:`2080`) :user:`katjaweigel`
-  Update contribution guidelines (:pull:`2031`) :user:`bouweandela`
-  Update links in pull request template to point to latest documentation (:pull:`2083`) :user:`bouweandela`
-  Update release schedule (:pull:`2081`) :user:`bouweandela`
-  Updates to contribution guidelines (:pull:`2092`) :user:`bouweandela`
-  Update documentation for ERA5 with new variables (:pull:`2111`) :user:`lukasbrunner`
-  Add OSX installation instructions to docs (:pull:`2115`) :user:`bvreede`
-  Instructions to use pre-installed versions on HPC clusters (:pull:`2197`) :user:`remi-kazeroni`
-  Add functional Autoassess diagnostics: land surface metrics: permafrost, soil moisture, surface radiation (:pull:`2170`) :user:`valeriupredoi`
-  Add citation info in `recipe_eady_growth_rate.yml` (:pull:`2188`) :user:`sloosvel`
-  Update version number to 2.3.0 (:pull:`2213`) :user:`zklaus`
-  Update release schedule for 2.3.0 (:pull:`2247`) :user:`zklaus`
-  Changelog update to v2.3.0 (:pull:`2214`) :user:`zklaus`

Diagnostics
~~~~~~~~~~~

-  Added figures 8 and 10 to recipe_bock20jgr.yml (:pull:`2074`) :user:`schlunma`
-  Add hydrological forcing comparison recipe (:pull:`2013`) :user:`stefsmeets`
-  Added recipe for Meehl et al., Sci. Adv. (2020) (:pull:`2094`) :user:`schlunma`
-  Add GlobWat recipe and diagnostic  (:pull:`1808`) :user:`babdollahi`
-  Add ClimWIP recipe to reproduce Brunner et al. 2019 (:pull:`2109`) :user:`lukasbrunner`
-  Update Climwip recipe to reproduce brunner2020esd (:pull:`1859`) :user:`ruthlorenz`
-  Update recipe_thermodyn_diagtool.yml: code improvements and more user options (:pull:`1391`) :user:`ValerioLembo`
-  Remove model AWI-CM-1-1-MR from recipe_impact.yml (:pull:`2238`) :user:`bouweandela`
-  PV capacity factor for ESMValTool GMD paper  (:pull:`2153`) :user:`katjaweigel`

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Cmorize wfde5 (:pull:`1991`) :user:`mwjury`
-  Make cmorizer utils funcs public in utilities.py and add some numpy style docstrings (:pull:`2206`) :user:`valeriupredoi`
-  CMORizer for CLARA-AVHRR cloud data (:pull:`2101`) :user:`axel-lauer`
-  Update of ESACCI-CLOUD CMORizer (:pull:`2144`) :user:`axel-lauer`

Automatic testing
~~~~~~~~~~~~~~~~~

-  Force latest Python in empty environment in conda install CI test (:pull:`2069`) :user:`valeriupredoi`
-  Removed imports from private sklearn modules and improved test coverage of custom_sklearn.py (:pull:`2078`) :user:`schlunma`
-  Move private _(global)_stock_cube from esmvacore.preprocessor._regrid to cmorizer (:pull:`2087`) :user:`valeriupredoi`
-  Try mamba install esmvaltool (:pull:`2125`) :user:`valeriupredoi`
-  Reinstate OSX Github Action tests (:pull:`2110`) :user:`valeriupredoi`
-  Pin mpich to avoid default install of 3.4.1 and 3.4.2 with external_0 builds (:pull:`2220`) :user:`valeriupredoi`
-  Include test sources in distribution (:pull:`2234`) :user:`zklaus`
-  Pin `iris<3.0.4` to ensure we still (sort of) support Python 3.6 (:pull:`2246`) :user:`valeriupredoi`

Installation
~~~~~~~~~~~~

-  Fix conda build by skipping documentation test (:pull:`2058`) Javier Vegas-Regidor
-  Update pin on esmvalcore pick up esmvalcore=2.3.0 (:pull:`2200`) :user:`valeriupredoi`
-  Pin Python to 3.9 for development installation (:pull:`2208`) :user:`bouweandela`

Improvements
~~~~~~~~~~~~

-  Add EUCP and IS-ENES3 projects to config-references (:pull:`2066`) :user:`Peter9192`
-  Fix flake8 tests on CircleCI (:pull:`2070`) :user:`bouweandela`
-  Added recipe filler. (:pull:`1707`) :user:`ledm`
-  Update use of fx vars to new syntax  (:pull:`2145`) :user:`sloosvel`
-  Add recipe for climate impact research (:pull:`2072`) :user:`Peter9192`
-  Update references "master" to "main" (:pull:`2172`) :user:`axel-lauer`
-  Force git to ignore VSCode workspace files (:pull:`2186`) Javier Vegas-Regidor
-  Update to new ESMValTool logo (:pull:`2168`) :user:`axel-lauer`
-  Python cmorizers for CDR1 and CDR2 ESACCI H2O (TCWV=prw) data. (:pull:`2152`) :user:`katjaweigel`
-  Remove obsolete conda package (closes #2100) (:pull:`2103`) :user:`zklaus`

.. _changelog-v2-2-0:

v2.2.0
------

Highlights
~~~~~~~~~~

ESMValTool is now using the recently released `Iris 3 <https://scitools-iris.readthedocs.io/en/latest/whatsnew/3.0.html>`__.
We acknowledge that this change may impact your work, as Iris 3 introduces
several changes that are not backward-compatible, but we think that moving forward is the best
decision for the tool in the long term.


This release includes

Bug fixes
~~~~~~~~~

-  Bugfix: time weights in time_operations (:pull:`1956`) :user:`axel-lauer`
-  Fix issues with bibtex references (:pull:`1955`) :user:`stefsmeets`
-  Fix ImportError for `configure_logging` (:pull:`1976`) :user:`stefsmeets`
-  Add required functional parameters for extract time in recipe_er5.yml (:pull:`1978`) :user:`valeriupredoi`
-  Revert "Fix ImportError for `configure_logging`" (:pull:`1992`) :user:`bouweandela`
-  Fix import of esmvalcore _logging module in cmorize_obs.py (:pull:`2020`) :user:`valeriupredoi`
-  Fix logging import in cmorize_obs again since last merge was nulled by pre-commit hooks (:pull:`2022`) :user:`valeriupredoi`
-  Refactor the functions in derive_evspsblpot due to new iris (:pull:`2023`) :user:`SarahAlidoost`
-  Avoid importing private ESMValCore functions in CMORizer (:pull:`2027`) :user:`bouweandela`
-  Fix extract_seasons in validation recipe  (:pull:`2054`) Javier Vegas-Regidor

Deprecations
~~~~~~~~~~~~

-  Deprecate classes Variable, Variables and Datasets (:pull:`1944`) :user:`schlunma`
-  Python 3.9: remove pynio as dependency and replace with rasterio and pin Matplotlib>3.3.1 and pin cartopy>=0.18 (:pull:`1997`) :user:`valeriupredoi`
-  Removed write_plots and write_netcdf in some python diagnostics (:pull:`2036`) :user:`schlunma`

Documentation
~~~~~~~~~~~~~

-  Update instructions on making a release (:pull:`1867`) :user:`bouweandela`
-  Update review.rst (:pull:`1917`) :user:`axel-lauer`
-  Add guidance on how to review a pull request (:pull:`1872`) :user:`bouweandela`
-  Adding tutorial links to documentation (:pull:`1927`) :user:`hb326`
-  Added bibtex file for schlund20jgr (:pull:`1928`) :user:`schlunma`
-  Documentation contact added the actual email for the mailing list (:pull:`1938`) :user:`valeriupredoi`
-  Make CircleCI badge specific to main branch (:pull:`1831`) :user:`bouweandela`
-  Documentation on how to move code from a private repository to a public repository (:pull:`1920`) :user:`hb326`
-  Refine pull request review guidelines (:pull:`1924`) :user:`stefsmeets`
-  Update release schedule (:pull:`1948`) :user:`zklaus`
-  Improve contact info and move to more prominent location (:pull:`1950`) :user:`bouweandela`
-  Add some maintainers to some recipes that are missing them (:pull:`1970`) :user:`valeriupredoi`
-  Update core team info (:pull:`1973`) :user:`axel-lauer`
-  Combine installation from source instructions and add common issues (:pull:`1971`) :user:`bouweandela`
-  Update iris documentation URL for sphinx (:pull:`2003`) :user:`bouweandela`
-  Fix iris documentation link(s) with new iris3 location on readthedocs (:pull:`2012`) :user:`valeriupredoi`
-  Document how to run tests for installation verification  (:pull:`1847`) :user:`valeriupredoi`
-  List Remi Kazeroni as a code owner and sole merger of CMORizers (:pull:`2017`) :user:`bouweandela`
-  Install documentation: mention that we build conda package with python>=3.7 (:pull:`2030`) :user:`valeriupredoi`
-  Recipe and documentation update for ERA5-Land. (:pull:`1906`) :user:`katjaweigel`
-  Update changelog and changelog tool for v2.2.0 (:pull:`2043`) Javier Vegas-Regidor
-  Final update to the changelog for v2.2.0 (:pull:`2056`) Javier Vegas-Regidor

Diagnostics
~~~~~~~~~~~

-  Add mapplot diagnostic to ClimWIP (:pull:`1864`) :user:`lukasbrunner`
-  Add the option to weight variable groups in ClimWIP (:pull:`1856`) :user:`lukasbrunner`
-  Implementation of ensemble member recognition to the ClimWIP diagnostic (:pull:`1852`) :user:`lukasbrunner`
-  Restructure ClimWIP (:pull:`1919`) :user:`lukasbrunner`
-  Diagnostic for recipe_eyring13jgr.yml Fig. 12 (:pull:`1922`) :user:`LisaBock`
-  Added changes in shared functions necessary for schlund20esd (:pull:`1967`) :user:`schlunma`
-  Adding recipe and diagnostics for Gier et al 2020 (:pull:`1914`) :user:`bettina-gier`
-  Added recipe, diagnostics and documentation for Schlund et al., ESD (2020) (:pull:`2015`) :user:`schlunma`
-  Add PRIMAVERA Eady Growth Rate diagnostic (:pull:`1285`) :user:`sloosvel`
-  Implement shape parameter calibration for ClimWIP (:pull:`1905`) :user:`lukasbrunner`

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Extended ESRL cmorizer (:pull:`1937`) :user:`bettina-gier`
-  Cmorizer for GRACE data (:pull:`1694`) :user:`bascrezee`
-  Cmorizer for latest ESACCI-SST data (:pull:`1895`) :user:`valeriupredoi`
-  Fix longitude in ESRL cmorizer (:pull:`1988`) :user:`bettina-gier`
-  Selectively turn off fixing bounds for coordinates during cmorization with utilities.py (:pull:`2014`) :user:`valeriupredoi`
-  Cmorize hadcrut5 (:pull:`1977`) :user:`mwjury`
-  Cmorize gpcc masking (:pull:`1995`) :user:`mwjury`
-  Cmorize_utils_save_1mon_Amon (:pull:`1990`) :user:`mwjury`
-  Cmorize gpcc fix (:pull:`1982`) :user:`mwjury`
-  Fix flake8 raised by develop test in cmorize_obs_gpcc.py (:pull:`2038`) :user:`valeriupredoi`

Automatic testing
~~~~~~~~~~~~~~~~~

-  Switched miniconda conda setup hooks for Github Actions workflows (:pull:`1913`) :user:`valeriupredoi`
-  Fix style issue (:pull:`1929`) :user:`bouweandela`
-  Fix mlr test with solution that works for CentOS too (:pull:`1936`) :user:`valeriupredoi`
-  Temporary deactivation Github Actions on OSX (:pull:`1939`) :user:`valeriupredoi`
-  Fix conda installation test on CircleCI (:pull:`1952`) :user:`bouweandela`
-  Github Actions: change time for cron job that installs from conda (:pull:`1969`) :user:`valeriupredoi`
-  CI upload relevant artifacts for test job (:pull:`1999`) :user:`valeriupredoi`
-  Github Actions test that runs with the latest ESMValCore main (:pull:`1989`) :user:`valeriupredoi`
-  Introduce python 39 in Github Actions tests (:pull:`2029`) :user:`valeriupredoi`
-  Remove test for conda package installation on Python 3.6 (:pull:`2033`) :user:`valeriupredoi`
-  Update codacy coverage reporter to fix coverage (:pull:`2039`) :user:`bouweandela`

Installation
~~~~~~~~~~~~

-  Simplify installation of R development dependencies (:pull:`1930`) :user:`bouweandela`
-  Fix docker build (:pull:`1934`) :user:`bouweandela`
-  Use new conda environment for installing ESMValTool in Docker containers (:pull:`1993`) :user:`bouweandela`
-  Fix conda build (:pull:`2026`) :user:`bouweandela`

Improvements
~~~~~~~~~~~~

-  Allow multiple references for a cmorizer script (:pull:`1953`) :user:`SarahAlidoost`
-  Add GRACE to the recipe check_obs (:pull:`1963`) :user:`remi-kazeroni`
-  Align ESMValTool to ESMValCore=2.2.0 (adopt iris3, fix environment for new Core release) (:pull:`1874`) :user:`stefsmeets`
-  Make it possible to use write_plots and write_netcdf from recipe instead of config-user.yml (:pull:`2018`) :user:`bouweandela`
-  Revise lisflood and hype recipes (:pull:`2035`) :user:`SarahAlidoost`
-  Set version to 2.2.0 (:pull:`2042`) Javier Vegas-Regidor

.. _changelog-v2-1-1:

v2.1.1
------

This release includes

Improvements
~~~~~~~~~~~~

- Fix the conda build on CircleCI (:pull:`1883`) :user:`bouweandela`
- Pin matplotlib to <3.3 and add compilers (:pull:`1898`) :user:`bouweandela`
- Pin esmvaltool subpackages to the same version and build as the esmvaltool conda package (:pull:`1899`) :user:`bouweandela`

Documentation
~~~~~~~~~~~~~

- Release notes v2.1.1 (:pull:`1932`) :user:`valeriupredoi`

.. _changelog-v2-1-0:

v2.1.0
------

This release includes

Diagnostics
~~~~~~~~~~~

-  Add extra steps to diagnostic to make output of hydrology/recipe_lisflood.yml usable by the LISFLOOD model (:pull:`1737`) :user:`JaroCamphuijsen`
-  Recipe to reproduce the 2014 KNMI Climate Scenarios (kcs). (:pull:`1667`) :user:`Peter9192`
-  Implement the climwip weighting scheme in a recipe and diagnostic (:pull:`1648`) :user:`JaroCamphuijsen`
-  Remove unreviewed autoassess recipes (:pull:`1840`) :user:`valeriupredoi`
-  Changes in shared scripts for Schlund et al., JGR: Biogeosciences, 2020 (:pull:`1845`) :user:`schlunma`
-  Updated derivation test recipe (:pull:`1790`) :user:`schlunma`
-  Support for multiple model occurrence in perf main (:pull:`1649`) :user:`bettina-gier`
-  Add recipe and diagnostics for Schlund et al., JGR: Biogeosciences, 2020 (:pull:`1860`) :user:`schlunma`
-  Adjust recipe_extract_shape.yml to recent changes in the example diagnostic.py (:pull:`1880`) :user:`bouweandela`

Documentation
~~~~~~~~~~~~~

-  Add pip installation instructions (:pull:`1783`) :user:`bouweandela`
-  Add installation instruction for R and Julia dependencies tot pip install (:pull:`1787`) :user:`bouweandela`
-  Avoid autodocsumm 0.2.0 and update documentation build dependencies (:pull:`1794`) :user:`bouweandela`
-  Add more information on working on cluster attached to ESGF node (:pull:`1821`) :user:`bouweandela`
-  Add release strategy to community documentation (:pull:`1809`) :user:`zklaus`
-  Update esmvaltool run command everywhere in documentation (:pull:`1820`) :user:`bouweandela`
-  Add more info on documenting a recipe (:pull:`1795`) :user:`bouweandela`
-  Improve the Python example diagnostic and documentation (:pull:`1827`) :user:`bouweandela`
-  Improve description of how to use draft_release_notes.py (:pull:`1848`) :user:`bouweandela`
-  Update changelog for release 2.1 (:pull:`1886`) :user:`valeriupredoi`

Improvements
~~~~~~~~~~~~

-  Fix R installation in WSL (:pull:`1789`) Javier Vegas-Regidor
-  Add pre-commit for linting/formatting (:pull:`1796`) :user:`stefsmeets`
-  Speed up tests on CircleCI and use pytest to run them (:pull:`1804`) :user:`bouweandela`
-  Move pre-commit excludes to top-level and correct order of lintr and styler (:pull:`1805`) :user:`stefsmeets`
-  Remove isort setup to fix formatting conflict with yapf (:pull:`1815`) :user:`stefsmeets`
-  GitHub Actions (:pull:`1806`) :user:`valeriupredoi`
-  Fix yapf-isort import formatting conflict (:pull:`1822`) :user:`stefsmeets`
-  Replace vmprof with vprof as the default profiler (:pull:`1829`) :user:`bouweandela`
-  Update ESMValCore v2.1.0 requirement (:pull:`1839`) Javier Vegas-Regidor
-  Pin iris to version 2 (:pull:`1881`) :user:`bouweandela`
-  Pin eccodes to not use eccodes=2.19.0 for cdo to work fine (:pull:`1869`) :user:`valeriupredoi`
-  Increase version to 2.1.0 and add release notes (:pull:`1868`) :user:`valeriupredoi`
-  Github Actions Build Packages and Deploy tests (conda and PyPi) (:pull:`1858`) :user:`valeriupredoi`

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Added CMORizer for Scripps-CO2-KUM (:pull:`1857`) :user:`schlunma`

.. _changelog-v2-0-0:

v2.0.0
------

This release includes

Bug fixes
~~~~~~~~~

-  Fix pep8-naming errors and fix zmnam diagnostic (:pull:`1702`) :user:`bouweandela`
-  Fix keyword argument in cmorize_obs (:pull:`1721`) :user:`mattiarighi`
-  Fixed JMA-TRANSCOM CMORizer (:pull:`1735`) :user:`schlunma`
-  Fix bug in extract_doi_value (:pull:`1734`) :user:`bascrezee`
-  Fix small errors in the arctic_ocean diagnostic (:pull:`1722`) :user:`koldunovn`
-  Flatten ancestor lists for diag_spei.R and diag_spi.R. (:pull:`1745`) :user:`katjaweigel`
-  Fix for recipe_ocean_ice_extent.yml (:pull:`1744`) :user:`mattiarighi`
-  Fix recipe_combined_indices.yml provenance (:pull:`1746`) Javier Vegas-Regidor
-  Fix provenance in recipe_multimodel_products (:pull:`1747`) Javier Vegas-Regidor
-  Exclude FGOALS-g2 due to ESMValCore issue #728 (:pull:`1749`) :user:`mattiarighi`
-  Fix recipe_modes_of_variability (:pull:`1753`) Javier Vegas-Regidor
-  Flatten lists for ancestors for hyint to prevent nested lists. (:pull:`1752`) :user:`katjaweigel`
-  Fix bug in cmorize_obs_eppley_vgpm_modis.py (#1729) (:pull:`1759`) :user:`tomaslovato`
-  Correct mip for clltkisccp in example derive preprocessor recipe (:pull:`1768`) :user:`bouweandela`
-  Update date conversion in recipe_hype.yml (:pull:`1769`) :user:`bouweandela`
-  Fix recipe_correlation.yml (:pull:`1767`) :user:`bouweandela`
-  Add attribute positive: down to plev coordinate in ERA-Interim CMORizer (:pull:`1771`) :user:`bouweandela`
-  Fix sispeed in recipe_preprocessor_derive_test (:pull:`1772`) Javier Vegas-Regidor
-  Fix extreme events and extreme index ancestors (:pull:`1774`) :user:`katjaweigel`
-  Correct date in output filenames of ERA5 CMORizer recipe (:pull:`1773`) :user:`bouweandela`
-  Exclude WOA from multi-model stats in recipe_ocean_bgc (:pull:`1778`) :user:`mattiarighi`

Diagnostics
~~~~~~~~~~~

-  Enhancement of the hyint recipe to include etccdi indices (:pull:`1133`) :user:`earnone`
-  Add lazy regridding for wflow diagnostic (:pull:`1630`) :user:`bouweandela`
-  Miles default domains to include lat=0 (:pull:`1626`) :user:`jhardenberg`
-  Miles: selection of reference dataset based on experiment (:pull:`1632`) :user:`jhardenberg`
-  New recipe/diagnostic:  recipe_li17natcc.yml for Axels GMD Paper (:pull:`1567`) :user:`katjaweigel`
-  New recipe/diagnostics: recipe_deangelis_for_gmdpart4.yml for Axels GMD Paper (:pull:`1576`) :user:`katjaweigel`
-  EWaterCycle: Add recipe to prepare input for LISFLOOD (:pull:`1298`) :user:`sverhoeven`
-  Use area weighted regridding in wflow diagnostic (:pull:`1643`) :user:`bouweandela`
-  Workaround for permetrics recipe until Iris3 (:pull:`1674`) :user:`mattiarighi`
-  C3S_511_MPQB_bas-features (:pull:`1465`) :user:`bascrezee`
-  Additional Land perfmetrics (:pull:`1641`) :user:`bettina-gier`
-  Necessary diagnostic from eyring06jgr for the release of version2 (:pull:`1686`) :user:`hb326`
-  Drought characteristics based on Martin2018 and SPI for gmd paper (:pull:`1689`) :user:`katjaweigel`
-  Additional features and bugfixes for recipe anav13clim (:pull:`1723`) :user:`bettina-gier`
-  Gmd laueretal2020 revisions (:pull:`1725`) :user:`axel-lauer`
-  Wenzel16nature (:pull:`1692`) :user:`zechlau`
-  Add mask albedolandcover (:pull:`1673`) :user:`bascrezee`
-  IPCC AR5 fig. 9.3 (seasonality) (:pull:`1726`) :user:`axel-lauer`
-  Added additional emergent constraints on ECS (:pull:`1585`) :user:`schlunma`
-  A diagnostic to evaluate the turnover times of land ecosystem carbon (:pull:`1395`) koir-su
-  Removed multi_model_statistics step in recipe_oceans_example.yml as a workaround (:pull:`1779`) :user:`valeriupredoi`

Documentation
~~~~~~~~~~~~~

-  Extend getting started instructions to obtain config-user.yml (:pull:`1642`) :user:`Peter9192`
-  Extend information about native6 support on RTD (:pull:`1652`) :user:`Peter9192`
-  Update citation of ESMValTool paper in the doc (:pull:`1664`) :user:`mattiarighi`
-  Updated references to documentation (now docs.esmvaltool.org) (:pull:`1679`) :user:`axel-lauer`
-  Replace dead link with ESGF link. (:pull:`1681`) :user:`mattiarighi`
-  Add all European grants to Zenodo (:pull:`1682`) :user:`bouweandela`
-  Update Sphinx to v3 or later (:pull:`1685`) :user:`bouweandela`
-  Small fix to number of models in ensclus documentation (:pull:`1691`) :user:`jhardenberg`
-  Move draft_release_notes.py from ESMValCore to here and update (:pull:`1701`) :user:`bouweandela`
-  Improve the installation instructions (:pull:`1634`) :user:`valeriupredoi`
-  Improve description of how to implement provenance in diagnostic (:pull:`1750`) :user:`SarahAlidoost`
-  Update command line interface documentation and add links to ESMValCore configuration documentation (:pull:`1776`) :user:`bouweandela`
-  Documentation on how to find shapefiles for hydrology recipes (:pull:`1777`) :user:`JaroCamphuijsen`

Improvements
~~~~~~~~~~~~

-  Pin flake8<3.8.0 (:pull:`1635`) :user:`valeriupredoi`
-  Update conda package path in more places (:pull:`1636`) :user:`bouweandela`
-  Remove curly brackets around issue number in pull request template (:pull:`1637`) :user:`bouweandela`
-  Fix style issue in test (:pull:`1639`) :user:`bouweandela`
-  Update Codacy badges (:pull:`1662`) :user:`bouweandela`
-  Support extra installation methods in R (:pull:`1360`) Javier Vegas-Regidor
-  Add ncdf4.helpers package as a dependency again (:pull:`1678`) :user:`bouweandela`
-  Speed up conda installation (:pull:`1677`) :user:`bouweandela`
-  Update CMORizers and recipes for ESMValCore v2.0.0 (:pull:`1699`) :user:`SarahAlidoost`
-  Update setup.py for PyPI package (:pull:`1700`) :user:`bouweandela`
-  Cleanup recipe headers before the release (:pull:`1740`) :user:`mattiarighi`
-    Add colortables as esmvaltool subcommand (:pull:`1666`) Javier Vegas-Regidor
-  Increase version to v2.0.0 (:pull:`1756`) :user:`bouweandela`
-  Update job script (:pull:`1757`) :user:`mattiarighi`
-  Read authors and description from .zenodo.json (:pull:`1758`) :user:`bouweandela`
-  Update docker recipe to install from source (:pull:`1651`) Javier Vegas-Regidor

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Cmorize aphro ma (:pull:`1555`) :user:`mwjury`
-  Respectable testing for cmorizers/obs/utilities.py and cmorizers/obs/cmorize_obs.py (:pull:`1517`) :user:`valeriupredoi`
-  Fix start year in recipe_check_obs (:pull:`1638`) :user:`mattiarighi`
-  Cmorizer for the PERSIANN-CDR precipitation data (:pull:`1633`) :user:`hb326`
-  Cmorize eobs (:pull:`1554`) :user:`mwjury`
-  Update download cds satellite lai fapar (:pull:`1654`) :user:`bascrezee`
-  Added monthly mean vars (ta, va, zg) to era5 cmorizer via recipe (:pull:`1644`) :user:`egalytska`
-  Make format time check more flexible (:pull:`1661`) :user:`mattiarighi`
-  Exclude od550lt1aer from recipe_check_obs.yml (:pull:`1720`) :user:`mattiarighi`
-  PERSIANN-CDR cmorizer update: adding the capability to save monthly mean files (:pull:`1728`) :user:`hb326`
-  Add standard_name attribute to lon and lat in cmorize_obs_esacci_oc.py (:pull:`1760`) :user:`tomaslovato`
-  Allow for incomplete months on daily frequency in cmorizer ncl utilities (:pull:`1754`) :user:`mattiarighi`
-  Fix AURA-TES cmorizer (:pull:`1766`) :user:`mattiarighi`

.. _changelog-v2-0-0b4:

v2.0.0b4
--------

This release includes

Bug fixes
~~~~~~~~~

-  Fix HALOE plev coordinate (:pull:`1590`) :user:`mattiarighi`
-  Fix tro3 units in HALOE (:pull:`1591`) :user:`mattiarighi`

Diagnostics
~~~~~~~~~~~

-  Applicate sea ice negative feedback (:pull:`1299`) Javier Vegas-Regidor
-  Add Russell18jgr ocean diagnostics (:pull:`1592`) :user:`bouweandela`
-  Refactor marrmot recipe and diagnostic to use ERA5 daily data made by new cmorizer (:pull:`1600`) :user:`SarahAlidoost`
-  In recipe_wflow, use daily ERA5 data from the new cmorizer. (:pull:`1599`) :user:`Peter9192`
-  In wflow diagnostic, calculate PET after(!) interpolation and lapse rate correction (:pull:`1618`) :user:`jeromaerts`
-  Fixed wenz14jgr (:pull:`1562`) :user:`zechlau`
-  Update portrait_plot.ncl (:pull:`1625`) :user:`bettina-gier`

Documentation
~~~~~~~~~~~~~

-  Restructure documentation (:pull:`1587`) :user:`bouweandela`
-  Add more links to documentation (:pull:`1595`) :user:`bouweandela`
-  Update links in readme (:pull:`1598`) :user:`bouweandela`
-  Minor improvements to installation documentation (:pull:`1608`) :user:`bouweandela`
-  Add info for new mailing list to documentation. (:pull:`1607`) :user:`bjoernbroetz`
-  Update making a release documentation (:pull:`1627`) :user:`bouweandela`

Improvements
~~~~~~~~~~~~

-  Avoid broken pytest-html plugin (:pull:`1583`) :user:`bouweandela`
-  Remove reference section in config-references.yml (:pull:`1545`) :user:`SarahAlidoost`
-  Various improvements to development infrastructure (:pull:`1570`) :user:`bouweandela`
-  Install scikit-learn from conda, remove libunwind as a direct dependency (:pull:`1611`) :user:`valeriupredoi`
-  Create conda subpackages and enable tests (:pull:`1624`) :user:`bouweandela`

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Cmorizer for HALOE (:pull:`1581`) :user:`mattiarighi`
-  Add CMORizer for CT2019 (:pull:`1604`) :user:`schlunma`

For older releases, see the release notes on https://github.com/ESMValGroup/ESMValTool/releases.
