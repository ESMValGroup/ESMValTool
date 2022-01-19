.. _inputdata:

********************
Obtaining input data
********************

ESMValTool supports input data from climate models participating in
`CMIP6 <https://www.wcrp-climate.org/wgcm-cmip/wgcm-cmip6>`__,
`CMIP5 <https://www.wcrp-climate.org/wgcm-cmip/wgcm-cmip5>`__,
`CMIP3 <https://www.wcrp-climate.org/wgcm-cmip/wgcm-cmip3>`__, and
`CORDEX <https://cordex.org/>`__
as well as observations, reanalysis, and any other data, provided that it
adheres to the
`CF conventions <https://cfconventions.org/>`__
and the data is described in a
`CMOR table <http://pcmdi.github.io/software/cmorTable/index.html>`__
as used in the various
`Climate Model Intercomparison Projects <http://pcmdi.github.io/mips/>`__.
This section provides some guidelines for unfamiliar users.

Because the amount of data required by ESMValTool is typically large, it is
recommended that you use the tool on a compute cluster where the data is
already available, for example because it is connected to an
`ESGF node <https://esgf.llnl.gov/index.html>`__.
Examples of such compute clusters are
`Mistral <https://www.dkrz.de/up/systems/mistral>`__
and
`Jasmin <https://www.jasmin.ac.uk/>`__,
but many more exist around the world.

If you do not have access to such a facility through your institute or the
project you are working on, you can request access by applying for the
`ENES Climate Analytics Service <https://portal.enes.org/data/data-metadata-service/climate-analytics-service>`__
or, if you need longer term access or more computational resources, the
`IS-ENES3 Trans-national Access call <https://portal.enes.org/data/data-metadata-service/analysis-platforms>`__.

If the options above are not available to you, ESMValTool also offers a feature
to make it easy to download CMIP6, CMIP5, CMIP3, CORDEX, and obs4MIPs from ESGF.

The chapter in the ESMValCore documentation on
:ref:`finding data <esmvalcore:findingdata>` explains how to
configure the ESMValTool so it can find locally available data and/or
download it from ESGF if it isn't available locally yet.

Models
======

If you do not have access to a compute cluster with the data already mounted,
the ESMValTool can automatically download any required data that is available on ESGF.
This is the recommended approach for first-time users to obtain some data for
running ESMValTool.
For example, run

.. code-block:: bash

    esmvaltool run --offline=False examples/recipe_python.yml

to run the default example recipe and automatically download the required data
to the directory ``~/climate_data``.
The data only needs to be downloaded once, every following run will re-use
previously downloaded data stored in this directory.
See :ref:`esmvalcore:config-esgf` for a more in depth explanation and the
available configuration options.

Alternatively, you can use an external tool called
`Synda <http://prodiguer.github.io/synda/index.html>`__
to maintain your own collection of ESGF data.

Observations
============

Observational and reanalysis products in the standard CF/CMOR format used in CMIP and required by the ESMValTool are available via the obs4MIPs and ana4mips projects at the ESGF (e.g., https://esgf-data.dkrz.de/projects/esgf-dkrz/). Their use is strongly recommended, when possible.

Other datasets not available in these archives can be obtained by the user from the respective sources and reformatted to the CF/CMOR standard. ESMValTool currently support two ways to perform this reformatting (aka 'CMORization'). The first is to use a CMORizer script to generate a local pool of reformatted data that can readily be used by the ESMValTool. The second way is to implement specific 'fixes' for your dataset. In that case, the reformatting is performed 'on the fly' during the execution of an ESMValTool recipe (note that one of the first preprocessor tasks is 'CMOR checks and fixes'). Below, both methods are explained in more detail.

Using a CMORizer script
-----------------------

ESMValTool comes with a set of CMORizers readily available.
The CMORizers are dataset-specific scripts that can be run once to generate
a local pool of CMOR-compliant data. The necessary information to download
and process the data is provided in the header of each CMORizing script.
These scripts also serve as template to create new CMORizers for datasets not
yet included.
Note that datasets CMORized for ESMValTool v1 may not be working with v2, due
to the much stronger constraints on metadata set by the iris library.

To CMORize one or more datasets, run:

.. code-block:: bash

    cmorize_obs -c [CONFIG_FILE] -o [DATASET_LIST]

The path to the raw data to be CMORized must be specified in the
:ref:`user configuration file<config-user>` as RAWOBS.
Within this path, the data are expected to be organized in subdirectories
corresponding to the data tier: Tier2 for freely-available datasets (other
than obs4MIPs and ana4mips) and Tier3 for restricted datasets (i.e., dataset
which requires a registration to be retrieved or provided upon request to
the respective contact or PI).
The CMORization follows the
`CMIP5 CMOR tables <https://github.com/PCMDI/cmip5-cmor-tables>`_ or
`CMIP6 CMOR tables <https://github.com/PCMDI/cmip6-cmor-tables>`_ for the
OBS and OBS6 projects respectively.
The resulting output is saved in the output_dir, again following the Tier
structure.
The output file names follow the definition given in
:ref:`config-developer file <esmvalcore:config-developer>` for the ``OBS``
project:

.. code-block::

    [project]_[dataset]_[type]_[version]_[mip]_[short_name]_YYYYMM_YYYYMM.nc

where ``project`` may be OBS (CMIP5 format) or OBS6 (CMIP6 format), ``type``
may be ``sat`` (satellite data), ``reanaly`` (reanalysis data),
``ground`` (ground observations), ``clim`` (derived climatologies),
``campaign`` (aircraft campaign).

At the moment, cmorize_obs supports Python and NCL scripts.

.. _cmorization_as_fix:

CMORization as a fix
--------------------
ESMValCore also provides support for some datasets in their native format.
In this case, the steps needed to reformat the data are executed as datasets
fixes during the execution of an ESMValTool recipe, as one of the first
preprocessor steps, see :ref:`fixing data <esmvalcore:fixing_data>`.
Compared to the workflow described above, this has the advantage that the user
does not need to store a duplicate (CMORized) copy of the data.
Instead, the CMORization is performed 'on the fly' when running a recipe.
The native6 project supports files named according to the format defined in
the :ref:`config-developer file <esmvalcore:config-developer>`.
Some of ERA5, ERA5-Land and MSWEP data are currently supported, see
:ref:`supported datasets <supported_datasets>`.

To use this functionality, users need to provide a path for the ``native6``
project data in the :ref:`user configuration file<config-user>`.
Then, in the recipe, they can refer to the native6 project.
For example:

.. code-block:: yaml

    datasets:
    - {dataset: ERA5, project: native6, type: reanaly, version: '1', tier: 3, start_year: 1990, end_year: 1990}

More examples can be found in the diagnostics ``ERA5_native6`` in the recipe
`examples/recipe_check_obs.yml <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/examples/recipe_check_obs.yml>`_.

.. _supported_datasets:

Supported datasets
------------------
A list of the datasets for which a CMORizers is available is provided in the following table.

.. tabularcolumns:: |p{3cm}|p{6cm}|p{3cm}|p{3cm}|

+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Dataset                      | Variables (MIP)                                                                                      | Tier | Script language |
+==============================+======================================================================================================+======+=================+
| APHRO-MA                     | pr, tas (day), pr, tas (Amon)                                                                        |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| AURA-TES                     | tro3 (Amon)                                                                                          |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| BerkelyEarth                 | tas, tasa (Amon), sftlf (fx)                                                                         |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CALIPSO-GOCCP                | clcalipso (cfMon)                                                                                    |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CDS-SATELLITE-ALBEDO         | bdalb (Lmon), bhalb (Lmon)                                                                           |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CDS-SATELLITE-LAI-FAPAR      | fapar (Lmon), lai (Lmon)                                                                             |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CDS-SATELLITE-SOIL-MOISTURE  | sm (day), sm (Lmon)                                                                                  |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CDS-UERRA                    | sm (E6hr)                                                                                            |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CDS-XCH4                     | xch4 (Amon)                                                                                          |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CDS-XCO2                     | xco2 (Amon)                                                                                          |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CERES-EBAF                   | rlut, rlutcs, rsut, rsutcs (Amon)                                                                    |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CERES-SYN1deg                | rlds, rldscs, rlus, rluscs, rlut, rlutcs, rsds, rsdscs, rsus, rsuscs, rsut, rsutcs (3hr)             |   3  | NCL             |
|                              | rlds, rldscs, rlus, rlut, rlutcs, rsds, rsdt, rsus, rsut, rsutcs (Amon)                              |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CLARA-AVHRR                  | clt, clivi, lwp (Amon)                                                                               |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CowtanWay                    | tasa (Amon)                                                                                          |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CRU                          | tas, pr (Amon)                                                                                       |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CT2019                       | co2s (Amon)                                                                                          |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Duveiller2018                | albDiffiTr13                                                                                         |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| E-OBS                        | tas, tasmin, tasmax, pr, psl (day, Amon)                                                             |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Eppley-VGPM-MODIS            | intpp (Omon)                                                                                         |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ERA5 [#note1]_               | clt, evspsbl, evspsblpot, mrro, pr, prsn, ps, psl, ptype, rls, rlds, rlns, rlus [#note2]_, rsds,     |   3  | n/a             |
|                              | rsns, rsus [#note2]_, rsdt, rss, uas, vas, tas, tasmax, tasmin, tdps, ts, tsn (E1hr/Amon), orog (fx) |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ERA5-Land [#note1]_          | pr                                                                                                   |   3  | n/a             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ERA-Interim                  | clivi, clt, clwvi, evspsbl, hur, hus, pr, prsn, prw, ps, psl, rlds, rsds, rsdt, ta, tas, tauu, tauv, |   3  | Python          |
|                              | ts, ua, uas, va, vas, wap, zg (Amon), ps, rsdt (CFday), clt, pr, prsn, psl, rsds, rss, ta, tas,      |      |                 |
|                              | tasmax, tasmin, uas, va, vas, zg (day), evspsbl, tdps, ts, tsn, rss, tdps (Eday), tsn (LImon), hfds, |      |                 |
|                              | tos (Omon), orog, sftlf (fx)                                                                         |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ERA-Interim-Land             | sm (Lmon)                                                                                            |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-AEROSOL               | abs550aer, od550aer, od550aerStderr, od550lt1aer, od870aer, od870aerStderr (aero)                    |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-CLOUD                 | clivi, clt, cltStderr, lwp, rlut, rlutcs, rsut, rsutcs, rsdt, rlus, rsus, rsuscs (Amon)              |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-FIRE                  | burntArea (Lmon)                                                                                     |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-LANDCOVER             | baresoilFrac, cropFrac, grassFrac, shrubFrac, treeFrac (Lmon)                                        |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-OC                    | chl (Omon)                                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-OZONE                 | toz, tozStderr, tro3prof, tro3profStderr (Amon)                                                      |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-SEA-SURFACE-SALINITY  | sos (Omon)                                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-SOILMOISTURE          | dos, dosStderr, sm, smStderr (Lmon)                                                                  |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-SST                   | ts, tsStderr (Amon)                                                                                  |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-WATERVAPOUR           | prw (Amon)                                                                                           |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESRL                         | co2s (Amon)                                                                                          |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| FLUXCOM                      | gpp (Lmon)                                                                                           |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GCP2018                      | fgco2 (Omon), nbp (Lmon)                                                                             |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GCP2020                      | fgco2 (Omon), nbp (Lmon)                                                                             |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GHCN                         | pr (Amon)                                                                                            |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GHCN-CAMS                    | tas (Amon)                                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GISTEMP                      | tasa (Amon)                                                                                          |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GLODAP                       | dissic, ph, talk (Oyr)                                                                               |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GPCC                         | pr (Amon)                                                                                            |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GRACE                        | lweGrace (Lmon)                                                                                      |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HadCRUT3                     | tas, tasa (Amon)                                                                                     |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HadCRUT4                     | tas, tasa (Amon), tasConf5, tasConf95                                                                |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HadCRUT5                     | tas (Amon)                                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HadISST                      | sic (OImon), tos (Omon), ts (Amon)                                                                   |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HALOE                        | tro3, hus (Amon)                                                                                     |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HWSD                         | cSoil (Lmon), areacella (fx), sftlf (fx)                                                             |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ISCCP-FH                     | alb, prw, ps, rlds, rlus, rlut, rlutcs, rsds, rsdt, rsus, rsut, rsutcs, tas, ts (Amon)               |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| JMA-TRANSCOM                 | nbp (Lmon), fgco2 (Omon)                                                                             |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| LAI3g                        | lai (Lmon)                                                                                           |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| LandFlux-EVAL                | et, etStderr (Lmon)                                                                                  |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Landschuetzer2016            | dpco2, fgco2, spco2 (Omon)                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MAC-LWP                      | lwp, lwpStderr (Amon)                                                                                |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MERRA2                       | sm (Lmon)                                                                                            |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MLS-AURA                     | hur, hurStderr (day)                                                                                 |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MODIS                        | cliwi, clt, clwvi, iwpStderr, lwpStderr (Amon), od550aer (aero)                                      |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MSWEP [#note1]_              | pr                                                                                                   |   3  | n/a             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MTE                          | gpp, gppStderr (Lmon)                                                                                |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NCEP                         | hur, hus, pr, ta, tas, ua, va, wap, zg (Amon)                                                        |   2  | NCL             |
|                              | pr, rlut, ua, va (day)                                                                               |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NDP                          | cVeg (Lmon)                                                                                          |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NIWA-BS                      | toz, tozStderr (Amon)                                                                                |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NSIDC-0116-[nh|sh]           | usi, vsi (day)                                                                                       |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| OSI-450-[nh|sh]              | sic (OImon), sic (day)                                                                               |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| PATMOS-x                     | clt (Amon)                                                                                           |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| PERSIANN-CDR                 | pr (Amon), pr (day)                                                                                  |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| PHC                          | thetao, so                                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| PIOMAS                       | sit (day)                                                                                            |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| REGEN                        | pr (day, Amon)                                                                                       |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Scripps-CO2-KUM              | co2s (Amon)                                                                                          |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| UWisc                        | clwvi, lwpStderr (Amon)                                                                              |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| WFDE5                        | tas, pr (Amon, day)                                                                                  |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| WOA                          | thetao, so, tos, sos (Omon)                                                                          |   2  | Python          |
|                              | no3, o2, po4, si (Oyr)                                                                               |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+

.. [#note1] CMORization is built into ESMValTool through the native6 project, so there is no separate CMORizer script.

.. [#note2] Derived on the fly from down & net radiation.
