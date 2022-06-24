.. _changelog:

Changelog
=========


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
- The new preprocessor ``ensemble_statistics`` can calculate arbitrary statitics over all ensemble members of a simulation.


This release includes

Backwards incompatible changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Streamline observations download (`#1657 <https://github.com/ESMValGroup/ESMValTool/pull/1657>`__) `Javier Vegas-Regidor <https://github.com/jvegreg>`__. This change removes the ``cmorize_obs`` command which has previously been used to CMORize observations and other datasets. The new command ``esmvaltool data`` provides many new features apart from the CMORization (``esmvaltool data format``), for example, automatic downloading of observational datasets (``esmvaltool data download``). More details on this can be found :ref:`here<inputdata_observations>` and :ref:`here<new-cmorizer>`.
-  Dropped Python 3.7 (`#2585 <https://github.com/ESMValGroup/ESMValTool/pull/2585>`__) `Manuel Schlund <https://github.com/schlunma>`__. ESMValTool v2.5.0 dropped support for Python 3.7. From now on Python >=3.8 is required to install ESMValTool. The main reason for this is that conda-forge dropped support for Python 3.7 for OSX and arm64 (more details are given `here <https://github.com/ESMValGroup/ESMValTool/issues/2584#issuecomment-1063853630>`__).

Bug fixes
~~~~~~~~~

-  Remove the use of `esmvalgroup` channel from the conda install Github Action workflow (`#2420 <https://github.com/ESMValGroup/ESMValTool/pull/2420>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Ignore .pymon-journal file in test discovery (`#2491 <https://github.com/ESMValGroup/ESMValTool/pull/2491>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Relocate pytest-monitor outputted database `.pymon` so `.pymon-journal` file should not be looked for by `pytest` (`#2501 <https://github.com/ESMValGroup/ESMValTool/pull/2501>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Re-establish Python 3.7 compatibility (`#2506 <https://github.com/ESMValGroup/ESMValTool/pull/2506>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Update intersphinx mapping (`#2531 <https://github.com/ESMValGroup/ESMValTool/pull/2531>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fixed `KeyError` in `recipe_ocean_bgc.yml` (`#2540 <https://github.com/ESMValGroup/ESMValTool/pull/2540>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Corrected ESACCI-SEA-SURFACE-SALINITY from OBS to OBS6 (`#2542 <https://github.com/ESMValGroup/ESMValTool/pull/2542>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Fixed `recipe_kcs.yml` (`#2541 <https://github.com/ESMValGroup/ESMValTool/pull/2541>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix MDER diagnostic regression_stepwise (`#2545 <https://github.com/ESMValGroup/ESMValTool/pull/2545>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Fix for recipe_wenzel16nat (`#2547 <https://github.com/ESMValGroup/ESMValTool/pull/2547>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Fixed `recipe_carvalhais14nat` and removed deprecated use of np.float (`#2558 <https://github.com/ESMValGroup/ESMValTool/pull/2558>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix `recipe_wenzel14jgr` (`#2577 <https://github.com/ESMValGroup/ESMValTool/pull/2577>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Fixed various recipes by removing faulty or non-available datasets (`#2563 <https://github.com/ESMValGroup/ESMValTool/pull/2563>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Remove missing CMIP5 data from 2 recipes (`#2579 <https://github.com/ESMValGroup/ESMValTool/pull/2579>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Fix `recipe_seaice` (`#2578 <https://github.com/ESMValGroup/ESMValTool/pull/2578>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Fix `recipe_climwip_brunner20esd` (`#2581 <https://github.com/ESMValGroup/ESMValTool/pull/2581>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__

Deprecations
~~~~~~~~~~~~

-  Remove `--use-feature=2020-resolver` command line option for obsolete pip 2020 solver (`#2493 <https://github.com/ESMValGroup/ESMValTool/pull/2493>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Renamed vertical regridding schemes in affected recipes (`#2487 <https://github.com/ESMValGroup/ESMValTool/pull/2487>`__) `Manuel Schlund <https://github.com/schlunma>`__

Documentation
~~~~~~~~~~~~~

-  Update release manager for v2.5 (`#2429 <https://github.com/ESMValGroup/ESMValTool/pull/2429>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Mention ENES Climate Analytics service (`#2438 <https://github.com/ESMValGroup/ESMValTool/pull/2438>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add recipe overview page (`#2439 <https://github.com/ESMValGroup/ESMValTool/pull/2439>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix pointer to Tutorial lesson on preprocessor from 05 to 06 (`#2473 <https://github.com/ESMValGroup/ESMValTool/pull/2473>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Removed obsolete option `synda-download` from documentation (`#2485 <https://github.com/ESMValGroup/ESMValTool/pull/2485>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Update CMUG XCH4 docu figure (`#2502 <https://github.com/ESMValGroup/ESMValTool/pull/2502>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Add Python=3.10 to package info, update Circle CI auto install and documentation for Python=3.10 (`#2503 <https://github.com/ESMValGroup/ESMValTool/pull/2503>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Unify user configuration file (`#2507 <https://github.com/ESMValGroup/ESMValTool/pull/2507>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Synchronized `config-user.yml` with version from ESMValCore (`#2516 <https://github.com/ESMValGroup/ESMValTool/pull/2516>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  CITATION.cff fix and automatic validation of your citation metadata (`#2517 <https://github.com/ESMValGroup/ESMValTool/pull/2517>`__) `Abel Siqueira <https://github.com/abelsiqueira>`__
-  Add backwards incompatible changes at the top of the release notes draft (`#2431 <https://github.com/ESMValGroup/ESMValTool/pull/2431>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fixed intersphinx mapping of `scipy` (`#2523 <https://github.com/ESMValGroup/ESMValTool/pull/2523>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Add authors to citation cff (`#2525 <https://github.com/ESMValGroup/ESMValTool/pull/2525>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Update documentation on running a recipe (`#2432 <https://github.com/ESMValGroup/ESMValTool/pull/2432>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix recipe `hydrology/recipe_wflow.yml` (`#2549 <https://github.com/ESMValGroup/ESMValTool/pull/2549>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Update `draft_release_notes.py` for new release (`#2553 <https://github.com/ESMValGroup/ESMValTool/pull/2553>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Added stand with Ukraine badge (`#2565 <https://github.com/ESMValGroup/ESMValTool/pull/2565>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Updated CREM docu (recipe_williams09climdyn.yml) (`#2567 <https://github.com/ESMValGroup/ESMValTool/pull/2567>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  First draft for v2.5.0 changelog (`#2554 <https://github.com/ESMValGroup/ESMValTool/pull/2554>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Replace nonfunctional Github Actions badge with cool one in README (`#2582 <https://github.com/ESMValGroup/ESMValTool/pull/2582>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Updated changelog (`#2589 <https://github.com/ESMValGroup/ESMValTool/pull/2589>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Updated release strategy with current release and upcoming release (`#2597 <https://github.com/ESMValGroup/ESMValTool/pull/2597>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Increased ESMValTool version to 2.5.0 (`#2600 <https://github.com/ESMValGroup/ESMValTool/pull/2600>`__) `Manuel Schlund <https://github.com/schlunma>`__

Diagnostics
~~~~~~~~~~~

-  AutoAssess: Add new diagnostic for radiation budget (`#2282 <https://github.com/ESMValGroup/ESMValTool/pull/2282>`__) `Jon Lillis <https://github.com/Jon-Lillis>`__
-  CMUG Sea Surface Salinity dataset and diagnostic (`#1832 <https://github.com/ESMValGroup/ESMValTool/pull/1832>`__) `Javier Vegas-Regidor <https://github.com/jvegreg>`__
-  Recipe with new diagnostics for ESA-CMUG H2O (`#1834 <https://github.com/ESMValGroup/ESMValTool/pull/1834>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Cleaned Schlund et al. (2020) recipe and fixed small bugs in corresponding diagnostic (`#2484 <https://github.com/ESMValGroup/ESMValTool/pull/2484>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Add ESA CCI LST cmorizer and diagnostic (`#1897 <https://github.com/ESMValGroup/ESMValTool/pull/1897>`__) `morobking <https://github.com/morobking>`__
-  XCH4 ESA CMUG diagnostics (subset of the MPQB diagnostics) (`#1960 <https://github.com/ESMValGroup/ESMValTool/pull/1960>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Add support for ESACCI Ocean Color (Chlorophyll) observations (`#2055 <https://github.com/ESMValGroup/ESMValTool/pull/2055>`__) `ulrikaw-cloud <https://github.com/ulrikaw-cloud>`__
-  Updated `recipe_zmnam.yml` with hemisphere selection (`#2230 <https://github.com/ESMValGroup/ESMValTool/pull/2230>`__) `fserva <https://github.com/fserva>`__
-  Add recipe and diagnostic scripts to compute figures of D9.4 of ISENES3 (`#2441 <https://github.com/ESMValGroup/ESMValTool/pull/2441>`__) `sloosvel <https://github.com/sloosvel>`__
-  Save resampled climates from KCS diagnostic local_resampling.py (`#2221 <https://github.com/ESMValGroup/ESMValTool/pull/2221>`__) `Emma Daniels <https://github.com/Emmadd>`__
-  Use years from KCS recipe (`#2223 <https://github.com/ESMValGroup/ESMValTool/pull/2223>`__) `Emma Daniels <https://github.com/Emmadd>`__
-  Recipe to plot generic output from the preprocessor (`#2184 <https://github.com/ESMValGroup/ESMValTool/pull/2184>`__) `Javier Vegas-Regidor <https://github.com/jvegreg>`__
-  Fixed provenance tracking for emergent constraint diagnostics (`#2573 <https://github.com/ESMValGroup/ESMValTool/pull/2573>`__) `Manuel Schlund <https://github.com/schlunma>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Ensure dummy data for cmorize_obs_woa test are written to the correct directory (`#2451 <https://github.com/ESMValGroup/ESMValTool/pull/2451>`__) `Emma Hogan <https://github.com/ehogan>`__
-  Add ESA CCI LST cmorizer and diagnostic (see previous section `Diagnostics`)

Automatic testing
~~~~~~~~~~~~~~~~~

-  Run a nightly Github Actions workflow to monitor tests memory per test (configurable for other metrics too) and lists the slowest 100 tests (`#2449 <https://github.com/ESMValGroup/ESMValTool/pull/2449>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix individual pytest runs broken due to missing explicit imports from `iris` and adding a couple missing package markers (`#2455 <https://github.com/ESMValGroup/ESMValTool/pull/2455>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Add Python=3.10 to Github Actions and switch to Python=3.10 for the Github Action that builds the PyPi package (`#2488 <https://github.com/ESMValGroup/ESMValTool/pull/2488>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Switch all github actions from miniconda to mambaforge (`#2498 <https://github.com/ESMValGroup/ESMValTool/pull/2498>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Pin `flake8<4` to have actual FLAKE8 error printed if tests fail and not garbage (`#2492 <https://github.com/ESMValGroup/ESMValTool/pull/2492>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Implementing conda lock (`#2193 <https://github.com/ESMValGroup/ESMValTool/pull/2193>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  [Docker] Update Docker container builds with correct installations of Julia (`#2530 <https://github.com/ESMValGroup/ESMValTool/pull/2530>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
- Update Linux condalock file (various pull requests) github-actions[bot]

Installation
~~~~~~~~~~~~

-  Comment out release candidate channel in environment.yml (`#2417 <https://github.com/ESMValGroup/ESMValTool/pull/2417>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Comment out rc channel in osx environment file (`#2421 <https://github.com/ESMValGroup/ESMValTool/pull/2421>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Add `python-cdo` as conda-forge dependency in environment files to ensure `cdo` gets used from conda-forge and not pip (`#2469 <https://github.com/ESMValGroup/ESMValTool/pull/2469>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Install rasterio from conda-forge and avoid issues from python=3.10 (`#2479 <https://github.com/ESMValGroup/ESMValTool/pull/2479>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Updated dependencies with new ESMValCore version (`#2599 <https://github.com/ESMValGroup/ESMValTool/pull/2599>`__) `Manuel Schlund <https://github.com/schlunma>`__

Improvements
~~~~~~~~~~~~

-  Remove use of OBS and use CMIP instead in `examples/recipe_ncl.yml` (`#2494 <https://github.com/ESMValGroup/ESMValTool/pull/2494>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Expanded `recipe_preprocessor_test.yml` to account for new `multi_model_statistics` features (`#2519 <https://github.com/ESMValGroup/ESMValTool/pull/2519>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Updated piControl periods for recipes that use KACE-1-0-G (`#2537 <https://github.com/ESMValGroup/ESMValTool/pull/2537>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Reduced time range in `recipe_globwat.yml` (`#2548 <https://github.com/ESMValGroup/ESMValTool/pull/2548>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Removed models with missing data from recipe_williams09climdyn.yml (`#2566 <https://github.com/ESMValGroup/ESMValTool/pull/2566>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Restored original versions of `recipe_schlund20esd.yml` and `recipe_meehl20sciadv.yml` (`#2583 <https://github.com/ESMValGroup/ESMValTool/pull/2583>`__) `Manuel Schlund <https://github.com/schlunma>`__


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

-  Fixed `recipe_meehl20sciadv.yml` for ESMValCore 2.3 (`#2253 <https://github.com/ESMValGroup/ESMValTool/pull/2253>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix provenance of NCL figures created using the log_provenance function (`#2279 <https://github.com/ESMValGroup/ESMValTool/pull/2279>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix bug in ClimWIP brunner19 recipe when plotting (`#2226 <https://github.com/ESMValGroup/ESMValTool/pull/2226>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Pin docutils <0.17 to fix sphinx build with rtd theme (`#2312 <https://github.com/ESMValGroup/ESMValTool/pull/2312>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix example recipes (`#2338 <https://github.com/ESMValGroup/ESMValTool/pull/2338>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Do not add bounds to plev (plev19) in era interim cmorizer (`#2328 <https://github.com/ESMValGroup/ESMValTool/pull/2328>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix problem with pip 21.3 that prevents installation from source (`#2344 <https://github.com/ESMValGroup/ESMValTool/pull/2344>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Add title to recipe embedded in test_diagnostic_run.py (`#2353 <https://github.com/ESMValGroup/ESMValTool/pull/2353>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix capitalization of obs4MIPs (`#2368 <https://github.com/ESMValGroup/ESMValTool/pull/2368>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Specify that areacella is needed for area statistics in the Python example recipe (`#2371 <https://github.com/ESMValGroup/ESMValTool/pull/2371>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Enabling variable `obs550lt1aer` in recipes (`#2388 <https://github.com/ESMValGroup/ESMValTool/pull/2388>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Update a diagnostic to new Iris version (`#2390 <https://github.com/ESMValGroup/ESMValTool/pull/2390>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Fixed bug in provenance tracking of ecs_scatter.ncl (`#2391 <https://github.com/ESMValGroup/ESMValTool/pull/2391>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix provenance issue in pv_capacity_factor.R (`#2392 <https://github.com/ESMValGroup/ESMValTool/pull/2392>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Remove obsolete write_plots option from R diagnostics (`#2395 <https://github.com/ESMValGroup/ESMValTool/pull/2395>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix arctic ocean diagnostic (`#2397 <https://github.com/ESMValGroup/ESMValTool/pull/2397>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix sea ice drift recipe and script (`#2404 <https://github.com/ESMValGroup/ESMValTool/pull/2404>`__) `sloosvel <https://github.com/sloosvel>`__
-  Adapt diagnostic script to new version of iris (`#2403 <https://github.com/ESMValGroup/ESMValTool/pull/2403>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix ocean multimap (`#2406 <https://github.com/ESMValGroup/ESMValTool/pull/2406>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix diagnostic that uses `xarray`: `dtype` correctly set and harmonize `xarray` and `matplotlib` (`#2409 <https://github.com/ESMValGroup/ESMValTool/pull/2409>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Deactivate provenance logging for plots in thermodyn toolbox (`#2414 <https://github.com/ESMValGroup/ESMValTool/pull/2414>`__) `Klaus Zimmermann <https://github.com/zklaus>`__

Deprecations
~~~~~~~~~~~~

-  Removed write_plots and write_netcdf from some NCL diagnostics (`#2293 <https://github.com/ESMValGroup/ESMValTool/pull/2293>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fixed provenance logging of all python diagnostics by removing 'plot_file' entry (`#2296 <https://github.com/ESMValGroup/ESMValTool/pull/2296>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Do not deprecate classes Variable, Variables and Datasets on a specific version (`#2286 <https://github.com/ESMValGroup/ESMValTool/pull/2286>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Remove obsolete write_netcdf option from ncl diagnostic scripts (`#2387 <https://github.com/ESMValGroup/ESMValTool/pull/2387>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Remove write plots from ocean diagnostics (`#2393 <https://github.com/ESMValGroup/ESMValTool/pull/2393>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  More removals of instances of `write_plots` from Python diagnostics (appears to be the final removal from Py diags) (`#2394 <https://github.com/ESMValGroup/ESMValTool/pull/2394>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Documentation
~~~~~~~~~~~~~

-  List Manuel Schlund as release manager for v2.5 (`#2268 <https://github.com/ESMValGroup/ESMValTool/pull/2268>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  GlobWat fix download links and gdal command (`#2334 <https://github.com/ESMValGroup/ESMValTool/pull/2334>`__) `Banafsheh Abdollahi <https://github.com/babdollahi>`__
-  Add titles to recipes authored by `predoi_valeriu` (`#2333 <https://github.com/ESMValGroup/ESMValTool/pull/2333>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Added titles to recipes maintained by lauer_axel (`#2332 <https://github.com/ESMValGroup/ESMValTool/pull/2332>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Update the documentation of the GRACE CMORizer (`#2349 <https://github.com/ESMValGroup/ESMValTool/pull/2349>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Add titles in BSC recipes (`#2351 <https://github.com/ESMValGroup/ESMValTool/pull/2351>`__) `sloosvel <https://github.com/sloosvel>`__
-  Update esmvalcore dependency to 2.4.0rc1 (`#2348 <https://github.com/ESMValGroup/ESMValTool/pull/2348>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Add titles to recipes maintained by Peter Kalverla (`#2356 <https://github.com/ESMValGroup/ESMValTool/pull/2356>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Adding titles to the recipes with maintainer hb326 (`#2358 <https://github.com/ESMValGroup/ESMValTool/pull/2358>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Add title for zmnam as for #2354 (`#2363 <https://github.com/ESMValGroup/ESMValTool/pull/2363>`__) `fserva <https://github.com/fserva>`__
-  Added recipe titles the the ocean recipes.  (`#2364 <https://github.com/ESMValGroup/ESMValTool/pull/2364>`__) `Lee de Mora <https://github.com/ledm>`__
-  Update recipe_thermodyn_diagtool.yml - add title (`#2365 <https://github.com/ESMValGroup/ESMValTool/pull/2365>`__) `ValerioLembo <https://github.com/ValerioLembo>`__
-  Fix provenance of figures of several R diagnostics (`#2300 <https://github.com/ESMValGroup/ESMValTool/pull/2300>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Adding titles to Mattia's recipes (`#2367 <https://github.com/ESMValGroup/ESMValTool/pull/2367>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Adding titles to wenzel recipes (`#2366 <https://github.com/ESMValGroup/ESMValTool/pull/2366>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Fix formatting of some recipe titles merged from PR 2364 (`#2372 <https://github.com/ESMValGroup/ESMValTool/pull/2372>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Adding titles to Bjoern's recipes (`#2369 <https://github.com/ESMValGroup/ESMValTool/pull/2369>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Add titles to ocean recipes (maintainer Lovato) (`#2375 <https://github.com/ESMValGroup/ESMValTool/pull/2375>`__) `Tomas Lovato <https://github.com/tomaslovato>`__
-  Add titles for three c3s-magic recipes (`#2378 <https://github.com/ESMValGroup/ESMValTool/pull/2378>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Add title for recipe maintained by Ruth Lorenz (`#2379 <https://github.com/ESMValGroup/ESMValTool/pull/2379>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix toymodel recipe (`#2381 <https://github.com/ESMValGroup/ESMValTool/pull/2381>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Added titles for recipes of maintainer `schlund_manuel` (`#2377 <https://github.com/ESMValGroup/ESMValTool/pull/2377>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Write_plots and titles for deangelis15nat, li17natcc, martin18grl, pv_capacity_factor (`#2382 <https://github.com/ESMValGroup/ESMValTool/pull/2382>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Add titles for some recipes (`#2383 <https://github.com/ESMValGroup/ESMValTool/pull/2383>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Adding titles for recipes by von Hardenberg and Arnone (`#2384 <https://github.com/ESMValGroup/ESMValTool/pull/2384>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Last two missing titles (`#2386 <https://github.com/ESMValGroup/ESMValTool/pull/2386>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update documentation on downloading data (`#2370 <https://github.com/ESMValGroup/ESMValTool/pull/2370>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix installation instructions for Julia (`#2335 <https://github.com/ESMValGroup/ESMValTool/pull/2335>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Fix provenance of Julia example diagnostic (`#2289 <https://github.com/ESMValGroup/ESMValTool/pull/2289>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Added notes on use of mamba in the installation documentation chapter (`#2236 <https://github.com/ESMValGroup/ESMValTool/pull/2236>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update version number for 2.4.0 release (`#2410 <https://github.com/ESMValGroup/ESMValTool/pull/2410>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Update release schedule for 2.4.0 (`#2412 <https://github.com/ESMValGroup/ESMValTool/pull/2412>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Update changelog for 2.4.0 release (`#2411 <https://github.com/ESMValGroup/ESMValTool/pull/2411>`__) `Klaus Zimmermann <https://github.com/zklaus>`__

Diagnostics
~~~~~~~~~~~

-  Add all available CMIP5 and CMIP6 models to recipe_impact.yml (`#2251 <https://github.com/ESMValGroup/ESMValTool/pull/2251>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add Fig. 6, 7 and 9 of Bock20jgr (`#2252 <https://github.com/ESMValGroup/ESMValTool/pull/2252>`__) `Lisa Bock <https://github.com/LisaBock>`__
-  Generalize `recipe_validation*` diagnostic to work with identical control and experiment dataset names (`#2284 <https://github.com/ESMValGroup/ESMValTool/pull/2284>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Add missing preprocessor to recipe_gier2020bg and adapt to available data (`#2399 <https://github.com/ESMValGroup/ESMValTool/pull/2399>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Removed custom version of `AtmosphereSigmaFactory` in diagnostics (`#2405 <https://github.com/ESMValGroup/ESMValTool/pull/2405>`__) `Manuel Schlund <https://github.com/schlunma>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Replace recipe_era5.yml with recipe_daily_era5.yml (`#2182 <https://github.com/ESMValGroup/ESMValTool/pull/2182>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Update WOA cmorizer for WOA18 and WOA13v2 (`#1812 <https://github.com/ESMValGroup/ESMValTool/pull/1812>`__) `Lisa Bock <https://github.com/LisaBock>`__
-  GLODAP v2.2016 ocean data cmorizer (`#2185 <https://github.com/ESMValGroup/ESMValTool/pull/2185>`__) `Tomas Lovato <https://github.com/tomaslovato>`__
-  Updated GCP CMORizer (`#2295 <https://github.com/ESMValGroup/ESMValTool/pull/2295>`__) `Manuel Schlund <https://github.com/schlunma>`__

Automatic testing
~~~~~~~~~~~~~~~~~

-  Add a cylc suite to run all recipes (`#2219 <https://github.com/ESMValGroup/ESMValTool/pull/2219>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Retire test with Python 3.6 from full development Github Actions test (`#2229 <https://github.com/ESMValGroup/ESMValTool/pull/2229>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Remove Python 3.6 tests from GitHub Actions (`#2264 <https://github.com/ESMValGroup/ESMValTool/pull/2264>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Unpin upper bound for iris (previously was at <3.0.4) (`#2266 <https://github.com/ESMValGroup/ESMValTool/pull/2266>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Pin latest esmvalcore to allow use of the bugfix release 2.3.1 always (`#2269 <https://github.com/ESMValGroup/ESMValTool/pull/2269>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Add apt update so Julia gets found and installed by Docker (`#2290 <https://github.com/ESMValGroup/ESMValTool/pull/2290>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Use mamba for environment update and creation in the Docker container build on DockerHub (`#2297 <https://github.com/ESMValGroup/ESMValTool/pull/2297>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Docker container experimental - run a full env solve with mamba instead of a conda update (`#2306 <https://github.com/ESMValGroup/ESMValTool/pull/2306>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Full use of mamba in Github Actions source install test and use generic Python 3.7 (removing the very specific 3.7.10) (`#2287 <https://github.com/ESMValGroup/ESMValTool/pull/2287>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Replace use of conda with mamba for conda_install test on Circle CI (`#2237 <https://github.com/ESMValGroup/ESMValTool/pull/2237>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update circleci configuration (`#2357 <https://github.com/ESMValGroup/ESMValTool/pull/2357>`__) `Klaus Zimmermann <https://github.com/zklaus>`__

Installation
~~~~~~~~~~~~

-  Remove `mpich` from conda dependencies list (`#2343 <https://github.com/ESMValGroup/ESMValTool/pull/2343>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Improvements
~~~~~~~~~~~~

-  Add script for extracting a list of input files from the provenance (`#2278 <https://github.com/ESMValGroup/ESMValTool/pull/2278>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update github actions (`#2360 <https://github.com/ESMValGroup/ESMValTool/pull/2360>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Removed 'write_plots' from all NCL diagnostics (`#2331 <https://github.com/ESMValGroup/ESMValTool/pull/2331>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Update and modernize `config-user-example.yml` (`#2374 <https://github.com/ESMValGroup/ESMValTool/pull/2374>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__


.. _changelog-v2-3-0:

v2.3.0
------

This release includes

Bug fixes
~~~~~~~~~

-  Indent block to pick up and raise exception if cmorizer data not found (TierX dir is not there) (`#1877 <https://github.com/ESMValGroup/ESMValTool/pull/1877>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Skip recipe filler tests until we have a new release since GA tests are failing (`#2089 <https://github.com/ESMValGroup/ESMValTool/pull/2089>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fixed broken link to contributions in README (`#2102 <https://github.com/ESMValGroup/ESMValTool/pull/2102>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix recipe filler for the case the variable doesn't contain short_name (`#2104 <https://github.com/ESMValGroup/ESMValTool/pull/2104>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Add fix for iris longitude bug to ClimWIP (`#2107 <https://github.com/ESMValGroup/ESMValTool/pull/2107>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Update for outdated link to reference Déandreis et al. (2014). (`#2076 <https://github.com/ESMValGroup/ESMValTool/pull/2076>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Fixed recipes for ESMValCore 2.3.0 (`#2203 <https://github.com/ESMValGroup/ESMValTool/pull/2203>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix the WFDE5 cmorizer (`#2211 <https://github.com/ESMValGroup/ESMValTool/pull/2211>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Fix broken CMORizer log message if no Tier directory exists (`#2207 <https://github.com/ESMValGroup/ESMValTool/pull/2207>`__) `jmrgonza <https://github.com/jmrgonza>`__
-  Fix bug in ClimWIP basic test recipe when plotting (`#2225 <https://github.com/ESMValGroup/ESMValTool/pull/2225>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Fix bug in ClimWIP advanced test recipe when plotting (`#2227 <https://github.com/ESMValGroup/ESMValTool/pull/2227>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Adjust time range for the `WDFE5` dataset in the `recipe_check_obs.yml` (`#2232 <https://github.com/ESMValGroup/ESMValTool/pull/2232>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Fix plot and provenance of recipe_consecdrydays (`#2244 <https://github.com/ESMValGroup/ESMValTool/pull/2244>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Documentation
~~~~~~~~~~~~~

-  Improving the README.md file with a more appealing look and bit more info (`#2065 <https://github.com/ESMValGroup/ESMValTool/pull/2065>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update plot title martin18grl (`#2080 <https://github.com/ESMValGroup/ESMValTool/pull/2080>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Update contribution guidelines (`#2031 <https://github.com/ESMValGroup/ESMValTool/pull/2031>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update links in pull request template to point to latest documentation (`#2083 <https://github.com/ESMValGroup/ESMValTool/pull/2083>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update release schedule (`#2081 <https://github.com/ESMValGroup/ESMValTool/pull/2081>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Updates to contribution guidelines (`#2092 <https://github.com/ESMValGroup/ESMValTool/pull/2092>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update documentation for ERA5 with new variables (`#2111 <https://github.com/ESMValGroup/ESMValTool/pull/2111>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Add OSX installation instructions to docs (`#2115 <https://github.com/ESMValGroup/ESMValTool/pull/2115>`__) `Barbara Vreede <https://github.com/bvreede>`__
-  Instructions to use pre-installed versions on HPC clusters (`#2197 <https://github.com/ESMValGroup/ESMValTool/pull/2197>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Add functional Autoassess diagnostics: land surface metrics: permafrost, soil moisture, surface radiation (`#2170 <https://github.com/ESMValGroup/ESMValTool/pull/2170>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Add citation info in `recipe_eady_growth_rate.yml` (`#2188 <https://github.com/ESMValGroup/ESMValTool/pull/2188>`__) `sloosvel <https://github.com/sloosvel>`__
-  Update version number to 2.3.0 (`#2213 <https://github.com/ESMValGroup/ESMValTool/pull/2213>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Update release schedule for 2.3.0 (`#2247 <https://github.com/ESMValGroup/ESMValTool/pull/2247>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Changelog update to v2.3.0 (`#2214 <https://github.com/ESMValGroup/ESMValTool/pull/2214>`__) `Klaus Zimmermann <https://github.com/zklaus>`__

Diagnostics
~~~~~~~~~~~

-  Added figures 8 and 10 to recipe_bock20jgr.yml (`#2074 <https://github.com/ESMValGroup/ESMValTool/pull/2074>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Add hydrological forcing comparison recipe (`#2013 <https://github.com/ESMValGroup/ESMValTool/pull/2013>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Added recipe for Meehl et al., Sci. Adv. (2020) (`#2094 <https://github.com/ESMValGroup/ESMValTool/pull/2094>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Add GlobWat recipe and diagnostic  (`#1808 <https://github.com/ESMValGroup/ESMValTool/pull/1808>`__) `Banafsheh Abdollahi <https://github.com/babdollahi>`__
-  Add ClimWIP recipe to reproduce Brunner et al. 2019 (`#2109 <https://github.com/ESMValGroup/ESMValTool/pull/2109>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Update Climwip recipe to reproduce brunner2020esd (`#1859 <https://github.com/ESMValGroup/ESMValTool/pull/1859>`__) `Ruth Lorenz <https://github.com/ruthlorenz>`__
-  Update recipe_thermodyn_diagtool.yml: code improvements and more user options (`#1391 <https://github.com/ESMValGroup/ESMValTool/pull/1391>`__) `ValerioLembo <https://github.com/ValerioLembo>`__
-  Remove model AWI-CM-1-1-MR from recipe_impact.yml (`#2238 <https://github.com/ESMValGroup/ESMValTool/pull/2238>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  PV capacity factor for ESMValTool GMD paper  (`#2153 <https://github.com/ESMValGroup/ESMValTool/pull/2153>`__) `katjaweigel <https://github.com/katjaweigel>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Cmorize wfde5 (`#1991 <https://github.com/ESMValGroup/ESMValTool/pull/1991>`__) `mwjury <https://github.com/mwjury>`__
-  Make cmorizer utils funcs public in utilities.py and add some numpy style docstrings (`#2206 <https://github.com/ESMValGroup/ESMValTool/pull/2206>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  CMORizer for CLARA-AVHRR cloud data (`#2101 <https://github.com/ESMValGroup/ESMValTool/pull/2101>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Update of ESACCI-CLOUD CMORizer (`#2144 <https://github.com/ESMValGroup/ESMValTool/pull/2144>`__) `Axel Lauer <https://github.com/axel-lauer>`__

Automatic testing
~~~~~~~~~~~~~~~~~

-  Force latest Python in empty environment in conda install CI test (`#2069 <https://github.com/ESMValGroup/ESMValTool/pull/2069>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Removed imports from private sklearn modules and improved test coverage of custom_sklearn.py (`#2078 <https://github.com/ESMValGroup/ESMValTool/pull/2078>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Move private _(global)_stock_cube from esmvacore.preprocessor._regrid to cmorizer (`#2087 <https://github.com/ESMValGroup/ESMValTool/pull/2087>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Try mamba install esmvaltool (`#2125 <https://github.com/ESMValGroup/ESMValTool/pull/2125>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Reinstate OSX Github Action tests (`#2110 <https://github.com/ESMValGroup/ESMValTool/pull/2110>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Pin mpich to avoid default install of 3.4.1 and 3.4.2 with external_0 builds (`#2220 <https://github.com/ESMValGroup/ESMValTool/pull/2220>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Include test sources in distribution (`#2234 <https://github.com/ESMValGroup/ESMValTool/pull/2234>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Pin `iris<3.0.4` to ensure we still (sort of) support Python 3.6 (`#2246 <https://github.com/ESMValGroup/ESMValTool/pull/2246>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Installation
~~~~~~~~~~~~

-  Fix conda build by skipping documentation test (`#2058 <https://github.com/ESMValGroup/ESMValTool/pull/2058>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Update pin on esmvalcore pick up esmvalcore=2.3.0 (`#2200 <https://github.com/ESMValGroup/ESMValTool/pull/2200>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Pin Python to 3.9 for development installation (`#2208 <https://github.com/ESMValGroup/ESMValTool/pull/2208>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Improvements
~~~~~~~~~~~~

-  Add EUCP and IS-ENES3 projects to config-references (`#2066 <https://github.com/ESMValGroup/ESMValTool/pull/2066>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Fix flake8 tests on CircleCI (`#2070 <https://github.com/ESMValGroup/ESMValTool/pull/2070>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Added recipe filler. (`#1707 <https://github.com/ESMValGroup/ESMValTool/pull/1707>`__) `Lee de Mora <https://github.com/ledm>`__
-  Update use of fx vars to new syntax  (`#2145 <https://github.com/ESMValGroup/ESMValTool/pull/2145>`__) `sloosvel <https://github.com/sloosvel>`__
-  Add recipe for climate impact research (`#2072 <https://github.com/ESMValGroup/ESMValTool/pull/2072>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Update references "master" to "main" (`#2172 <https://github.com/ESMValGroup/ESMValTool/pull/2172>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Force git to ignore VSCode workspace files (`#2186 <https://github.com/ESMValGroup/ESMValTool/pull/2186>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Update to new ESMValTool logo (`#2168 <https://github.com/ESMValGroup/ESMValTool/pull/2168>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Python cmorizers for CDR1 and CDR2 ESACCI H2O (TCWV=prw) data. (`#2152 <https://github.com/ESMValGroup/ESMValTool/pull/2152>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Remove obsolete conda package (closes #2100) (`#2103 <https://github.com/ESMValGroup/ESMValTool/pull/2103>`__) `Klaus Zimmermann <https://github.com/zklaus>`__

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

-  Bugfix: time weights in time_operations (`#1956 <https://github.com/ESMValGroup/ESMValTool/pull/1956>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Fix issues with bibtex references (`#1955 <https://github.com/ESMValGroup/ESMValTool/pull/1955>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Fix ImportError for `configure_logging` (`#1976 <https://github.com/ESMValGroup/ESMValTool/pull/1976>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Add required functional parameters for extract time in recipe_er5.yml (`#1978 <https://github.com/ESMValGroup/ESMValTool/pull/1978>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Revert "Fix ImportError for `configure_logging`" (`#1992 <https://github.com/ESMValGroup/ESMValTool/pull/1992>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix import of esmvalcore _logging module in cmorize_obs.py (`#2020 <https://github.com/ESMValGroup/ESMValTool/pull/2020>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix logging import in cmorize_obs again since last merge was nulled by pre-commit hooks (`#2022 <https://github.com/ESMValGroup/ESMValTool/pull/2022>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Refactor the functions in derive_evspsblpot due to new iris (`#2023 <https://github.com/ESMValGroup/ESMValTool/pull/2023>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Avoid importing private ESMValCore functions in CMORizer (`#2027 <https://github.com/ESMValGroup/ESMValTool/pull/2027>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix extract_seasons in validation recipe  (`#2054 <https://github.com/ESMValGroup/ESMValTool/pull/2054>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__

Deprecations
~~~~~~~~~~~~

-  Deprecate classes Variable, Variables and Datasets (`#1944 <https://github.com/ESMValGroup/ESMValTool/pull/1944>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Python 3.9: remove pynio as dependency and replace with rasterio and pin Matplotlib>3.3.1 and pin cartopy>=0.18 (`#1997 <https://github.com/ESMValGroup/ESMValTool/pull/1997>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Removed write_plots and write_netcdf in some python diagnostics (`#2036 <https://github.com/ESMValGroup/ESMValTool/pull/2036>`__) `Manuel Schlund <https://github.com/schlunma>`__

Documentation
~~~~~~~~~~~~~

-  Update instructions on making a release (`#1867 <https://github.com/ESMValGroup/ESMValTool/pull/1867>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update review.rst (`#1917 <https://github.com/ESMValGroup/ESMValTool/pull/1917>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Add guidance on how to review a pull request (`#1872 <https://github.com/ESMValGroup/ESMValTool/pull/1872>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Adding tutorial links to documentation (`#1927 <https://github.com/ESMValGroup/ESMValTool/pull/1927>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Added bibtex file for schlund20jgr (`#1928 <https://github.com/ESMValGroup/ESMValTool/pull/1928>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Documentation contact added the actual email for the mailing list (`#1938 <https://github.com/ESMValGroup/ESMValTool/pull/1938>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Make CircleCI badge specific to main branch (`#1831 <https://github.com/ESMValGroup/ESMValTool/pull/1831>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Documentation on how to move code from a private repository to a public repository (`#1920 <https://github.com/ESMValGroup/ESMValTool/pull/1920>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Refine pull request review guidelines (`#1924 <https://github.com/ESMValGroup/ESMValTool/pull/1924>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Update release schedule (`#1948 <https://github.com/ESMValGroup/ESMValTool/pull/1948>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Improve contact info and move to more prominent location (`#1950 <https://github.com/ESMValGroup/ESMValTool/pull/1950>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add some maintainers to some recipes that are missing them (`#1970 <https://github.com/ESMValGroup/ESMValTool/pull/1970>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update core team info (`#1973 <https://github.com/ESMValGroup/ESMValTool/pull/1973>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Combine installation from source instructions and add common issues (`#1971 <https://github.com/ESMValGroup/ESMValTool/pull/1971>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update iris documentation URL for sphinx (`#2003 <https://github.com/ESMValGroup/ESMValTool/pull/2003>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix iris documentation link(s) with new iris3 location on readthedocs (`#2012 <https://github.com/ESMValGroup/ESMValTool/pull/2012>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Document how to run tests for installation verification  (`#1847 <https://github.com/ESMValGroup/ESMValTool/pull/1847>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  List Remi Kazeroni as a code owner and sole merger of CMORizers (`#2017 <https://github.com/ESMValGroup/ESMValTool/pull/2017>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Install documentation: mention that we build conda package with python>=3.7 (`#2030 <https://github.com/ESMValGroup/ESMValTool/pull/2030>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Recipe and documentation update for ERA5-Land. (`#1906 <https://github.com/ESMValGroup/ESMValTool/pull/1906>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Update changelog and changelog tool for v2.2.0 (`#2043 <https://github.com/ESMValGroup/ESMValTool/pull/2043>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Final update to the changelog for v2.2.0 (`#2056 <https://github.com/ESMValGroup/ESMValTool/pull/2056>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__

Diagnostics
~~~~~~~~~~~

-  Add mapplot diagnostic to ClimWIP (`#1864 <https://github.com/ESMValGroup/ESMValTool/pull/1864>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Add the option to weight variable groups in ClimWIP (`#1856 <https://github.com/ESMValGroup/ESMValTool/pull/1856>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Implementation of ensemble member recognition to the ClimWIP diagnostic (`#1852 <https://github.com/ESMValGroup/ESMValTool/pull/1852>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Restructure ClimWIP (`#1919 <https://github.com/ESMValGroup/ESMValTool/pull/1919>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Diagnostic for recipe_eyring13jgr.yml Fig. 12 (`#1922 <https://github.com/ESMValGroup/ESMValTool/pull/1922>`__) `Lisa Bock <https://github.com/LisaBock>`__
-  Added changes in shared functions necessary for schlund20esd (`#1967 <https://github.com/ESMValGroup/ESMValTool/pull/1967>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Adding recipe and diagnostics for Gier et al 2020 (`#1914 <https://github.com/ESMValGroup/ESMValTool/pull/1914>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Added recipe, diagnostics and documentation for Schlund et al., ESD (2020) (`#2015 <https://github.com/ESMValGroup/ESMValTool/pull/2015>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Add PRIMAVERA Eady Growth Rate diagnostic (`#1285 <https://github.com/ESMValGroup/ESMValTool/pull/1285>`__) `sloosvel <https://github.com/sloosvel>`__
-  Implement shape parameter calibration for ClimWIP (`#1905 <https://github.com/ESMValGroup/ESMValTool/pull/1905>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Extended ESRL cmorizer (`#1937 <https://github.com/ESMValGroup/ESMValTool/pull/1937>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Cmorizer for GRACE data (`#1694 <https://github.com/ESMValGroup/ESMValTool/pull/1694>`__) `bascrezee <https://github.com/bascrezee>`__
-  Cmorizer for latest ESACCI-SST data (`#1895 <https://github.com/ESMValGroup/ESMValTool/pull/1895>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix longitude in ESRL cmorizer (`#1988 <https://github.com/ESMValGroup/ESMValTool/pull/1988>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Selectively turn off fixing bounds for coordinates during cmorization with utilities.py (`#2014 <https://github.com/ESMValGroup/ESMValTool/pull/2014>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Cmorize hadcrut5 (`#1977 <https://github.com/ESMValGroup/ESMValTool/pull/1977>`__) `mwjury <https://github.com/mwjury>`__
-  Cmorize gpcc masking (`#1995 <https://github.com/ESMValGroup/ESMValTool/pull/1995>`__) `mwjury <https://github.com/mwjury>`__
-  Cmorize_utils_save_1mon_Amon (`#1990 <https://github.com/ESMValGroup/ESMValTool/pull/1990>`__) `mwjury <https://github.com/mwjury>`__
-  Cmorize gpcc fix (`#1982 <https://github.com/ESMValGroup/ESMValTool/pull/1982>`__) `mwjury <https://github.com/mwjury>`__
-  Fix flake8 raised by develop test in cmorize_obs_gpcc.py (`#2038 <https://github.com/ESMValGroup/ESMValTool/pull/2038>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Automatic testing
~~~~~~~~~~~~~~~~~

-  Switched miniconda conda setup hooks for Github Actions workflows (`#1913 <https://github.com/ESMValGroup/ESMValTool/pull/1913>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix style issue (`#1929 <https://github.com/ESMValGroup/ESMValTool/pull/1929>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix mlr test with solution that works for CentOS too (`#1936 <https://github.com/ESMValGroup/ESMValTool/pull/1936>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Temporary deactivation Github Actions on OSX (`#1939 <https://github.com/ESMValGroup/ESMValTool/pull/1939>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix conda installation test on CircleCI (`#1952 <https://github.com/ESMValGroup/ESMValTool/pull/1952>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Github Actions: change time for cron job that installs from conda (`#1969 <https://github.com/ESMValGroup/ESMValTool/pull/1969>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  CI upload relevant artifacts for test job (`#1999 <https://github.com/ESMValGroup/ESMValTool/pull/1999>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Github Actions test that runs with the latest ESMValCore main (`#1989 <https://github.com/ESMValGroup/ESMValTool/pull/1989>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Introduce python 39 in Github Actions tests (`#2029 <https://github.com/ESMValGroup/ESMValTool/pull/2029>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Remove test for conda package installation on Python 3.6 (`#2033 <https://github.com/ESMValGroup/ESMValTool/pull/2033>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update codacy coverage reporter to fix coverage (`#2039 <https://github.com/ESMValGroup/ESMValTool/pull/2039>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Installation
~~~~~~~~~~~~

-  Simplify installation of R development dependencies (`#1930 <https://github.com/ESMValGroup/ESMValTool/pull/1930>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix docker build (`#1934 <https://github.com/ESMValGroup/ESMValTool/pull/1934>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Use new conda environment for installing ESMValTool in Docker containers (`#1993 <https://github.com/ESMValGroup/ESMValTool/pull/1993>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix conda build (`#2026 <https://github.com/ESMValGroup/ESMValTool/pull/2026>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Improvements
~~~~~~~~~~~~

-  Allow multiple references for a cmorizer script (`#1953 <https://github.com/ESMValGroup/ESMValTool/pull/1953>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Add GRACE to the recipe check_obs (`#1963 <https://github.com/ESMValGroup/ESMValTool/pull/1963>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Align ESMValTool to ESMValCore=2.2.0 (adopt iris3, fix environment for new Core release) (`#1874 <https://github.com/ESMValGroup/ESMValTool/pull/1874>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Make it possible to use write_plots and write_netcdf from recipe instead of config-user.yml (`#2018 <https://github.com/ESMValGroup/ESMValTool/pull/2018>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Revise lisflood and hype recipes (`#2035 <https://github.com/ESMValGroup/ESMValTool/pull/2035>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Set version to 2.2.0 (`#2042 <https://github.com/ESMValGroup/ESMValTool/pull/2042>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__

.. _changelog-v2-1-1:

v2.1.1
------

This release includes

Improvements
~~~~~~~~~~~~

- Fix the conda build on CircleCI (`#1883 <https://github.com/ESMValGroup/ESMValTool/pull/1883>`__) `Bouwe Andela <https://github.com/bouweandela>`__
- Pin matplotlib to <3.3 and add compilers (`#1898 <https://github.com/ESMValGroup/ESMValTool/pull/1898>`__) `Bouwe Andela <https://github.com/bouweandela>`__
- Pin esmvaltool subpackages to the same version and build as the esmvaltool conda package (`#1899 <https://github.com/ESMValGroup/ESMValTool/pull/1899>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Documentation
~~~~~~~~~~~~~

- Release notes v2.1.1 (`#1932 <https://github.com/ESMValGroup/ESMValTool/pull/1932>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

.. _changelog-v2-1-0:

v2.1.0
------

This release includes

Diagnostics
~~~~~~~~~~~

-  Add extra steps to diagnostic to make output of hydrology/recipe_lisflood.yml usable by the LISFLOOD model (`#1737 <https://github.com/ESMValGroup/ESMValTool/pull/1737>`__) `Jaro Camphuijsen <https://github.com/JaroCamphuijsen>`__
-  Recipe to reproduce the 2014 KNMI Climate Scenarios (kcs). (`#1667 <https://github.com/ESMValGroup/ESMValTool/pull/1667>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Implement the climwip weighting scheme in a recipe and diagnostic (`#1648 <https://github.com/ESMValGroup/ESMValTool/pull/1648>`__) `Jaro Camphuijsen <https://github.com/JaroCamphuijsen>`__
-  Remove unreviewed autoassess recipes (`#1840 <https://github.com/ESMValGroup/ESMValTool/pull/1840>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Changes in shared scripts for Schlund et al., JGR: Biogeosciences, 2020 (`#1845 <https://github.com/ESMValGroup/ESMValTool/pull/1845>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Updated derivation test recipe (`#1790 <https://github.com/ESMValGroup/ESMValTool/pull/1790>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Support for multiple model occurrence in perf main (`#1649 <https://github.com/ESMValGroup/ESMValTool/pull/1649>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Add recipe and diagnostics for Schlund et al., JGR: Biogeosciences, 2020 (`#1860 <https://github.com/ESMValGroup/ESMValTool/pull/1860>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Adjust recipe_extract_shape.yml to recent changes in the example diagnostic.py (`#1880 <https://github.com/ESMValGroup/ESMValTool/pull/1880>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Documentation
~~~~~~~~~~~~~

-  Add pip installation instructions (`#1783 <https://github.com/ESMValGroup/ESMValTool/pull/1783>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add installation instruction for R and Julia dependencies tot pip install (`#1787 <https://github.com/ESMValGroup/ESMValTool/pull/1787>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Avoid autodocsumm 0.2.0 and update documentation build dependencies (`#1794 <https://github.com/ESMValGroup/ESMValTool/pull/1794>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add more information on working on cluster attached to ESGF node (`#1821 <https://github.com/ESMValGroup/ESMValTool/pull/1821>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add release strategy to community documentation (`#1809 <https://github.com/ESMValGroup/ESMValTool/pull/1809>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Update esmvaltool run command everywhere in documentation (`#1820 <https://github.com/ESMValGroup/ESMValTool/pull/1820>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add more info on documenting a recipe (`#1795 <https://github.com/ESMValGroup/ESMValTool/pull/1795>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Improve the Python example diagnostic and documentation (`#1827 <https://github.com/ESMValGroup/ESMValTool/pull/1827>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Improve description of how to use draft_release_notes.py (`#1848 <https://github.com/ESMValGroup/ESMValTool/pull/1848>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update changelog for release 2.1 (`#1886 <https://github.com/ESMValGroup/ESMValTool/pull/1886>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Improvements
~~~~~~~~~~~~

-  Fix R installation in WSL (`#1789 <https://github.com/ESMValGroup/ESMValTool/pull/1789>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Add pre-commit for linting/formatting (`#1796 <https://github.com/ESMValGroup/ESMValTool/pull/1796>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Speed up tests on CircleCI and use pytest to run them (`#1804 <https://github.com/ESMValGroup/ESMValTool/pull/1804>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Move pre-commit excludes to top-level and correct order of lintr and styler (`#1805 <https://github.com/ESMValGroup/ESMValTool/pull/1805>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Remove isort setup to fix formatting conflict with yapf (`#1815 <https://github.com/ESMValGroup/ESMValTool/pull/1815>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  GitHub Actions (`#1806 <https://github.com/ESMValGroup/ESMValTool/pull/1806>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix yapf-isort import formatting conflict (`#1822 <https://github.com/ESMValGroup/ESMValTool/pull/1822>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Replace vmprof with vprof as the default profiler (`#1829 <https://github.com/ESMValGroup/ESMValTool/pull/1829>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update ESMValCore v2.1.0 requirement (`#1839 <https://github.com/ESMValGroup/ESMValTool/pull/1839>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Pin iris to version 2 (`#1881 <https://github.com/ESMValGroup/ESMValTool/pull/1881>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Pin eccodes to not use eccodes=2.19.0 for cdo to work fine (`#1869 <https://github.com/ESMValGroup/ESMValTool/pull/1869>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Increase version to 2.1.0 and add release notes (`#1868 <https://github.com/ESMValGroup/ESMValTool/pull/1868>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Github Actions Build Packages and Deploy tests (conda and PyPi) (`#1858 <https://github.com/ESMValGroup/ESMValTool/pull/1858>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Added CMORizer for Scripps-CO2-KUM (`#1857 <https://github.com/ESMValGroup/ESMValTool/pull/1857>`__) `Manuel Schlund <https://github.com/schlunma>`__

.. _changelog-v2-0-0:

v2.0.0
------

This release includes

Bug fixes
~~~~~~~~~

-  Fix pep8-naming errors and fix zmnam diagnostic (`#1702 <https://github.com/ESMValGroup/ESMValTool/pull/1702>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix keyword argument in cmorize_obs (`#1721 <https://github.com/ESMValGroup/ESMValTool/pull/1721>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Fixed JMA-TRANSCOM CMORizer (`#1735 <https://github.com/ESMValGroup/ESMValTool/pull/1735>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix bug in extract_doi_value (`#1734 <https://github.com/ESMValGroup/ESMValTool/pull/1734>`__) `bascrezee <https://github.com/bascrezee>`__
-  Fix small errors in the arctic_ocean diagnostic (`#1722 <https://github.com/ESMValGroup/ESMValTool/pull/1722>`__) `Nikolay Koldunov <https://github.com/koldunovn>`__
-  Flatten ancestor lists for diag_spei.R and diag_spi.R. (`#1745 <https://github.com/ESMValGroup/ESMValTool/pull/1745>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Fix for recipe_ocean_ice_extent.yml (`#1744 <https://github.com/ESMValGroup/ESMValTool/pull/1744>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Fix recipe_combined_indices.yml provenance (`#1746 <https://github.com/ESMValGroup/ESMValTool/pull/1746>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Fix provenance in recipe_multimodel_products (`#1747 <https://github.com/ESMValGroup/ESMValTool/pull/1747>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Exclude FGOALS-g2 due to ESMValCore issue #728 (`#1749 <https://github.com/ESMValGroup/ESMValTool/pull/1749>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Fix recipe_modes_of_variability (`#1753 <https://github.com/ESMValGroup/ESMValTool/pull/1753>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Flatten lists for ancestors for hyint to prevent nested lists. (`#1752 <https://github.com/ESMValGroup/ESMValTool/pull/1752>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Fix bug in cmorize_obs_eppley_vgpm_modis.py (#1729) (`#1759 <https://github.com/ESMValGroup/ESMValTool/pull/1759>`__) `Tomas Lovato <https://github.com/tomaslovato>`__
-  Correct mip for clltkisccp in example derive preprocessor recipe (`#1768 <https://github.com/ESMValGroup/ESMValTool/pull/1768>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update date conversion in recipe_hype.yml (`#1769 <https://github.com/ESMValGroup/ESMValTool/pull/1769>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix recipe_correlation.yml (`#1767 <https://github.com/ESMValGroup/ESMValTool/pull/1767>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add attribute positive: down to plev coordinate in ERA-Interim CMORizer (`#1771 <https://github.com/ESMValGroup/ESMValTool/pull/1771>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix sispeed in recipe_preprocessor_derive_test (`#1772 <https://github.com/ESMValGroup/ESMValTool/pull/1772>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Fix extreme events and extreme index ancestors (`#1774 <https://github.com/ESMValGroup/ESMValTool/pull/1774>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Correct date in output filenames of ERA5 CMORizer recipe (`#1773 <https://github.com/ESMValGroup/ESMValTool/pull/1773>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Exclude WOA from multi-model stats in recipe_ocean_bgc (`#1778 <https://github.com/ESMValGroup/ESMValTool/pull/1778>`__) `Mattia Righi <https://github.com/mattiarighi>`__

Diagnostics
~~~~~~~~~~~

-  Enhancement of the hyint recipe to include etccdi indices (`#1133 <https://github.com/ESMValGroup/ESMValTool/pull/1133>`__) `Enrico Arnone <https://github.com/earnone>`__
-  Add lazy regridding for wflow diagnostic (`#1630 <https://github.com/ESMValGroup/ESMValTool/pull/1630>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Miles default domains to include lat=0 (`#1626 <https://github.com/ESMValGroup/ESMValTool/pull/1626>`__) `Jost von Hardenberg <https://github.com/jhardenberg>`__
-  Miles: selection of reference dataset based on experiment (`#1632 <https://github.com/ESMValGroup/ESMValTool/pull/1632>`__) `Jost von Hardenberg <https://github.com/jhardenberg>`__
-  New recipe/diagnostic:  recipe_li17natcc.yml for Axels GMD Paper (`#1567 <https://github.com/ESMValGroup/ESMValTool/pull/1567>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  New recipe/diagnostics: recipe_deangelis_for_gmdpart4.yml for Axels GMD Paper (`#1576 <https://github.com/ESMValGroup/ESMValTool/pull/1576>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  EWaterCycle: Add recipe to prepare input for LISFLOOD (`#1298 <https://github.com/ESMValGroup/ESMValTool/pull/1298>`__) `Stefan Verhoeven <https://github.com/sverhoeven>`__
-  Use area weighted regridding in wflow diagnostic (`#1643 <https://github.com/ESMValGroup/ESMValTool/pull/1643>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Workaround for permetrics recipe until Iris3 (`#1674 <https://github.com/ESMValGroup/ESMValTool/pull/1674>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  C3S_511_MPQB_bas-features (`#1465 <https://github.com/ESMValGroup/ESMValTool/pull/1465>`__) `bascrezee <https://github.com/bascrezee>`__
-  Additional Land perfmetrics (`#1641 <https://github.com/ESMValGroup/ESMValTool/pull/1641>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Necessary diagnostic from eyring06jgr for the release of version2 (`#1686 <https://github.com/ESMValGroup/ESMValTool/pull/1686>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Drought characteristics based on Martin2018 and SPI for gmd paper (`#1689 <https://github.com/ESMValGroup/ESMValTool/pull/1689>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Additional features and bugfixes for recipe anav13clim (`#1723 <https://github.com/ESMValGroup/ESMValTool/pull/1723>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Gmd laueretal2020 revisions (`#1725 <https://github.com/ESMValGroup/ESMValTool/pull/1725>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Wenzel16nature (`#1692 <https://github.com/ESMValGroup/ESMValTool/pull/1692>`__) `zechlau <https://github.com/zechlau>`__
-  Add mask albedolandcover (`#1673 <https://github.com/ESMValGroup/ESMValTool/pull/1673>`__) `bascrezee <https://github.com/bascrezee>`__
-  IPCC AR5 fig. 9.3 (seasonality) (`#1726 <https://github.com/ESMValGroup/ESMValTool/pull/1726>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Added additional emergent constraints on ECS (`#1585 <https://github.com/ESMValGroup/ESMValTool/pull/1585>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  A diagnostic to evaluate the turnover times of land ecosystem carbon (`#1395 <https://github.com/ESMValGroup/ESMValTool/pull/1395>`__) `koir-su <https://github.com/koir-su>`__
-  Removed multi_model_statistics step in recipe_oceans_example.yml as a workaround (`#1779 <https://github.com/ESMValGroup/ESMValTool/pull/1779>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Documentation
~~~~~~~~~~~~~

-  Extend getting started instructions to obtain config-user.yml (`#1642 <https://github.com/ESMValGroup/ESMValTool/pull/1642>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Extend information about native6 support on RTD (`#1652 <https://github.com/ESMValGroup/ESMValTool/pull/1652>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Update citation of ESMValTool paper in the doc (`#1664 <https://github.com/ESMValGroup/ESMValTool/pull/1664>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Updated references to documentation (now docs.esmvaltool.org) (`#1679 <https://github.com/ESMValGroup/ESMValTool/pull/1679>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Replace dead link with ESGF link. (`#1681 <https://github.com/ESMValGroup/ESMValTool/pull/1681>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Add all European grants to Zenodo (`#1682 <https://github.com/ESMValGroup/ESMValTool/pull/1682>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update Sphinx to v3 or later (`#1685 <https://github.com/ESMValGroup/ESMValTool/pull/1685>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Small fix to number of models in ensclus documentation (`#1691 <https://github.com/ESMValGroup/ESMValTool/pull/1691>`__) `Jost von Hardenberg <https://github.com/jhardenberg>`__
-  Move draft_release_notes.py from ESMValCore to here and update (`#1701 <https://github.com/ESMValGroup/ESMValTool/pull/1701>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Improve the installation instructions (`#1634 <https://github.com/ESMValGroup/ESMValTool/pull/1634>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Improve description of how to implement provenance in diagnostic (`#1750 <https://github.com/ESMValGroup/ESMValTool/pull/1750>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Update command line interface documentation and add links to ESMValCore configuration documentation (`#1776 <https://github.com/ESMValGroup/ESMValTool/pull/1776>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Documentation on how to find shapefiles for hydrology recipes (`#1777 <https://github.com/ESMValGroup/ESMValTool/pull/1777>`__) `Jaro Camphuijsen <https://github.com/JaroCamphuijsen>`__

Improvements
~~~~~~~~~~~~

-  Pin flake8<3.8.0 (`#1635 <https://github.com/ESMValGroup/ESMValTool/pull/1635>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update conda package path in more places (`#1636 <https://github.com/ESMValGroup/ESMValTool/pull/1636>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Remove curly brackets around issue number in pull request template (`#1637 <https://github.com/ESMValGroup/ESMValTool/pull/1637>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix style issue in test (`#1639 <https://github.com/ESMValGroup/ESMValTool/pull/1639>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update Codacy badges (`#1662 <https://github.com/ESMValGroup/ESMValTool/pull/1662>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Support extra installation methods in R (`#1360 <https://github.com/ESMValGroup/ESMValTool/pull/1360>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Add ncdf4.helpers package as a dependency again (`#1678 <https://github.com/ESMValGroup/ESMValTool/pull/1678>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Speed up conda installation (`#1677 <https://github.com/ESMValGroup/ESMValTool/pull/1677>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update CMORizers and recipes for ESMValCore v2.0.0 (`#1699 <https://github.com/ESMValGroup/ESMValTool/pull/1699>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Update setup.py for PyPI package (`#1700 <https://github.com/ESMValGroup/ESMValTool/pull/1700>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Cleanup recipe headers before the release (`#1740 <https://github.com/ESMValGroup/ESMValTool/pull/1740>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-    Add colortables as esmvaltool subcommand (`#1666 <https://github.com/ESMValGroup/ESMValTool/pull/1666>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Increase version to v2.0.0 (`#1756 <https://github.com/ESMValGroup/ESMValTool/pull/1756>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update job script (`#1757 <https://github.com/ESMValGroup/ESMValTool/pull/1757>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Read authors and description from .zenodo.json (`#1758 <https://github.com/ESMValGroup/ESMValTool/pull/1758>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update docker recipe to install from source (`#1651 <https://github.com/ESMValGroup/ESMValTool/pull/1651>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Cmorize aphro ma (`#1555 <https://github.com/ESMValGroup/ESMValTool/pull/1555>`__) `mwjury <https://github.com/mwjury>`__
-  Respectable testing for cmorizers/obs/utilities.py and cmorizers/obs/cmorize_obs.py (`#1517 <https://github.com/ESMValGroup/ESMValTool/pull/1517>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix start year in recipe_check_obs (`#1638 <https://github.com/ESMValGroup/ESMValTool/pull/1638>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Cmorizer for the PERSIANN-CDR precipitation data (`#1633 <https://github.com/ESMValGroup/ESMValTool/pull/1633>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Cmorize eobs (`#1554 <https://github.com/ESMValGroup/ESMValTool/pull/1554>`__) `mwjury <https://github.com/mwjury>`__
-  Update download cds satellite lai fapar (`#1654 <https://github.com/ESMValGroup/ESMValTool/pull/1654>`__) `bascrezee <https://github.com/bascrezee>`__
-  Added monthly mean vars (ta, va, zg) to era5 cmorizer via recipe (`#1644 <https://github.com/ESMValGroup/ESMValTool/pull/1644>`__) `Evgenia Galytska <https://github.com/egalytska>`__
-  Make format time check more flexible (`#1661 <https://github.com/ESMValGroup/ESMValTool/pull/1661>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Exclude od550lt1aer from recipe_check_obs.yml (`#1720 <https://github.com/ESMValGroup/ESMValTool/pull/1720>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  PERSIANN-CDR cmorizer update: adding the capability to save monthly mean files (`#1728 <https://github.com/ESMValGroup/ESMValTool/pull/1728>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Add standard_name attribute to lon and lat in cmorize_obs_esacci_oc.py (`#1760 <https://github.com/ESMValGroup/ESMValTool/pull/1760>`__) `Tomas Lovato <https://github.com/tomaslovato>`__
-  Allow for incomplete months on daily frequency in cmorizer ncl utilities (`#1754 <https://github.com/ESMValGroup/ESMValTool/pull/1754>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Fix AURA-TES cmorizer (`#1766 <https://github.com/ESMValGroup/ESMValTool/pull/1766>`__) `Mattia Righi <https://github.com/mattiarighi>`__

.. _changelog-v2-0-0b4:

v2.0.0b4
--------

This release includes

Bug fixes
~~~~~~~~~

-  Fix HALOE plev coordinate (`#1590 <https://github.com/ESMValGroup/ESMValTool/pull/1590>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Fix tro3 units in HALOE (`#1591 <https://github.com/ESMValGroup/ESMValTool/pull/1591>`__) `Mattia Righi <https://github.com/mattiarighi>`__

Diagnostics
~~~~~~~~~~~

-  Applicate sea ice negative feedback (`#1299 <https://github.com/ESMValGroup/ESMValTool/pull/1299>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Add Russell18jgr ocean diagnostics (`#1592 <https://github.com/ESMValGroup/ESMValTool/pull/1592>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Refactor marrmot recipe and diagnostic to use ERA5 daily data made by new cmorizer (`#1600 <https://github.com/ESMValGroup/ESMValTool/pull/1600>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  In recipe_wflow, use daily ERA5 data from the new cmorizer. (`#1599 <https://github.com/ESMValGroup/ESMValTool/pull/1599>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  In wflow diagnostic, calculate PET after(!) interpolation and lapse rate correction (`#1618 <https://github.com/ESMValGroup/ESMValTool/pull/1618>`__) `Jerom Aerts <https://github.com/jeromaerts>`__
-  Fixed wenz14jgr (`#1562 <https://github.com/ESMValGroup/ESMValTool/pull/1562>`__) `zechlau <https://github.com/zechlau>`__
-  Update portrait_plot.ncl (`#1625 <https://github.com/ESMValGroup/ESMValTool/pull/1625>`__) `Bettina Gier <https://github.com/bettina-gier>`__

Documentation
~~~~~~~~~~~~~

-  Restructure documentation (`#1587 <https://github.com/ESMValGroup/ESMValTool/pull/1587>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add more links to documentation (`#1595 <https://github.com/ESMValGroup/ESMValTool/pull/1595>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update links in readme (`#1598 <https://github.com/ESMValGroup/ESMValTool/pull/1598>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Minor improvements to installation documentation (`#1608 <https://github.com/ESMValGroup/ESMValTool/pull/1608>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add info for new mailing list to documentation. (`#1607 <https://github.com/ESMValGroup/ESMValTool/pull/1607>`__) `Björn Brötz <https://github.com/bjoernbroetz>`__
-  Update making a release documentation (`#1627 <https://github.com/ESMValGroup/ESMValTool/pull/1627>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Improvements
~~~~~~~~~~~~

-  Avoid broken pytest-html plugin (`#1583 <https://github.com/ESMValGroup/ESMValTool/pull/1583>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Remove reference section in config-references.yml (`#1545 <https://github.com/ESMValGroup/ESMValTool/pull/1545>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Various improvements to development infrastructure (`#1570 <https://github.com/ESMValGroup/ESMValTool/pull/1570>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Install scikit-learn from conda, remove libunwind as a direct dependency (`#1611 <https://github.com/ESMValGroup/ESMValTool/pull/1611>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Create conda subpackages and enable tests (`#1624 <https://github.com/ESMValGroup/ESMValTool/pull/1624>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Cmorizer for HALOE (`#1581 <https://github.com/ESMValGroup/ESMValTool/pull/1581>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Add CMORizer for CT2019 (`#1604 <https://github.com/ESMValGroup/ESMValTool/pull/1604>`__) `Manuel Schlund <https://github.com/schlunma>`__

For older releases, see the release notes on https://github.com/ESMValGroup/ESMValTool/releases.
