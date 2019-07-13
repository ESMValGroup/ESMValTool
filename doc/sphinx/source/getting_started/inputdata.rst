.. _inputdata:

**********
Input data
**********

Models
======

Observations
============

Observational and reanalysis products in the standard CF/CMOR format used in CMIP and required by the ESMValTool are available via the obs4mips (https://esgf-node.llnl.gov/projects/obs4mips/) and ana4mips (https://esgf.nccs.nasa.gov/projects/ana4mips/) proejcts, respectively. Their use is strongly recommended, when possible.

Other datasets not available in these archives can be obtained by the user from the respective sources and reformatted to the CF/CMOR standard using the cmorizers included in the ESMValTool. The cmorizers are dataset-specific scripts that can be run once to generate a local pool of observational datasets for usage with the ESMValTool. The necessary information to download and process the data is provided in the header of each cmorizing script. These scripts also serve as template to create new cmorizers for datasets not yet included. Note that dataset cmorized for ESMValTool v1 may not be working with v2, due to the much stronger constraints on metadata set by the Iris library.

To cmorize one or more datasets, run:

.. code-block:: bash

    cmorize_obs -c [CONFIG_FILE] -o [DATASET_LIST]

The path to the raw data to be cmorized must be specified in the CONFIG_FILE as RAWOBS. Within this path, the data are expected to be organized in subdirectories corresponding to the data tier: Tier2 for freely-available datasets (other than obs4mips and ana4mips) and Tier3 for restricted datasets (i.e., dataset which requires a registration to be retrieved or provided upon request to the respective contact or PI). The cmorization follows the CMIP5 CMOR tables. The resulting output is saved in the output_dir, again following the Tier structure. The output file names follow the definition given in ``config-developer.yml`` for the ``OBS`` project: ``OBS_[dataset]_[type]_[version]_[mip]_[short_name]_YYYYMM_YYYYMM.nc``, where ``type`` may be ``sat`` (satellite data), ``reanaly`` (reanalysis data), ``ground`` (ground observations), ``clim`` (derived climatologies), ``campaign`` (aircraft campaign).


At the moment, cmorize_obs supports Python and NCL scripts.

A list of the datasets for which a cmorizers is available is provided in the following table.

.. tabularcolumns:: |p{3cm}|p{6cm}|p{3cm}|p{3cm}|

+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Dataset                      | Variables (MIP)                                                                                      | Tier | Script language |
+==============================+======================================================================================================+======+=================+
| AURA-TES                     | tro3 (Amon)                                                                                          |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CDS-SATELLITE-SOIL-MOISTURE  | sm (Lmon), smStderr (Lmon)                                                                           |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CDS-XCH4                     | xch4 (Amon)                                                                                          |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CDS-XCO2                     | xco2 (Amon)                                                                                          |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CERES-SYN1deg                | rlds, rldscs, rlus, rluscs, rlut, rlutcs, rsds, rsdscs, rsus, rsuscs, rsut, rsutcs (3hr)             |   3  | NCL             |
|                              | rlds, rldscs, rlus, rlut, rlutcs, rsds, rsdt, rsus, rsut, rsutcs (Amon)                              |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CRU                          | tas, pr (Amon)                                                                                       |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ERA-Interim                  | clivi, clt, clwvi, hfds, hur, hus, pr, prw, ps, psl, ta, tas, tauu, tauv, ts, ua, va, wap, zg (Amon) |   3  | NCL             |
|                              | pr, psl, tas, tasmin, tasmax, zg (day), sftlf (fx), tos (Omon)                                       |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-AEROSOL               | abs550aer, od550aer, od550aerStderr, od550lt1aer, od870aer, od870aerStderr (aero)                    |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-CLOUD                 | clivi, clt, cltStderr, clwvi (Amon)                                                                  |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-FIRE                  | burntArea (Lmon)                                                                                     |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-LANDCOVER             | baresoilFrac, cropFrac, grassFrac, shrubFrac, treeFrac (Lmon)                                        |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-OC                    | chl                                                                                                  |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-OZONE                 | toz, tozStderr, tro3prof, tro3profStderr (Amon)                                                      |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-SOILMOISTURE          | dos, dosStderr, sm, smStderr (Lmon)                                                                  |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-SST                   | ts, tsStderr (Amon)                                                                                  |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GHCN                         | pr (Amon)                                                                                            |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HadCRUT3                     | tas, tasa (Amon)                                                                                     |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HadCRUT4                     | tas, tasa (Amon)                                                                                     |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HadISST                      | sic (OImon), tos (Omon), ts (Amon)                                                                   |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| LandFlux-EVAL                | et, etStderr (Lmon)                                                                                  |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Landschuetzer2016            | fgco2 (Omon), spco2 (Omon), dpco2 (Omon)                                                             |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MODIS                        | cliwi, clt, clwvi, iwpStderr, lwpStderr (Amon), od550aer (aero)                                      |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MTE                          | gpp, gppStderr (Lmon)                                                                                |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NCEP                         | hur, hus, pr, ta, tas, ua, va, wap, zg (Amon)                                                        |   2  | NCL             |
|                              | pr, rlut, ua, va (day)                                                                               |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NIWA-BS                      | toz, tozStderr (Amon)                                                                                |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| PATMOS-x                     | clt (Amon)                                                                                           |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| UWisc                        | clwvi, lwpStderr (Amon)                                                                              |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| WOA                          | no3, o2, po4, si (Oyr), so, thetao (Omon)                                                            |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
