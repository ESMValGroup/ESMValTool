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

.. _cordex_note:

.. note::

    CORDEX support is still
    `work in progress <https://github.com/orgs/ESMValGroup/projects/11>`__.
    Contributions, in the form of
    :ref:`pull request reviews <reviewing>` or
    :ref:`pull requests <esmvalcore:contributing>`
    are most welcome. We are particularly interested in contributions from
    people with good understanding of the CORDEX project and its standards.

This section provides an introduction to getting (access to) climate data
for use with ESMValTool.

Because the amount of data required by ESMValTool is typically large, it is
recommended that you use the tool on a compute cluster where the data is
already available, for example because it is connected to an
`ESGF node <https://esgf.llnl.gov/index.html>`__.
Examples of such compute clusters are
`Levante <https://docs.dkrz.de/doc/levante/index.html>`__
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
ESMValTool also provides support to download some observational dataset from source.

The chapter in the ESMValCore documentation on
:ref:`finding data <esmvalcore:findingdata>` explains how to
configure ESMValTool so it can find locally available data and/or
download it from ESGF if it isn't available locally yet.


.. _inputdata_models:

Models
======

If you do not have access to a compute cluster with the data already mounted,
ESMValTool can automatically download any required data that is available on
ESGF.
This is the recommended approach for first-time users to obtain some data for
running ESMValTool.
For example, run

.. code-block:: bash

    esmvaltool run --search_esgf=when_missing examples/recipe_python.yml

to run the default example recipe and automatically download the required data
to the directory ``~/climate_data``.
The data only needs to be downloaded once, every following run will reuse
previously downloaded data stored in this directory.
See :ref:`esmvalcore:config-esgf` for a more in depth explanation and the
available configuration options.

Alternatively, you can use an external tool called
`Synda <http://prodiguer.github.io/synda/index.html>`__
to maintain your own collection of ESGF data.


.. _inputdata_observations:

Observations
============

Observational and reanalysis products in the standard CF/CMOR format used in
CMIP and required by ESMValTool are available via the obs4MIPs and ana4mips
projects at the ESGF (e.g., https://esgf-data.dkrz.de/projects/esgf-dkrz/).
Their use is strongly recommended, when possible.

Other datasets not available in these archives can be obtained by the user from
the respective sources and reformatted to the CF/CMOR standard.
ESMValTool currently supports two ways to perform this reformatting (aka
'CMORization'):

#. Using a CMORizer script: The first is to use a CMORizer script to generate a
   local pool of reformatted data that can readily be used by ESMValTool.  This
   method is described in detail below.

#. Using fixes for on-the-fly CMORization: The second way is to implement
   specific :ref:`'fixes' <esmvalcore:fixing_data>` for your dataset.  In that
   case, the reformatting is performed 'on the fly' during the execution of an
   ESMValTool recipe (note that one of the first preprocessor tasks is 'CMOR
   checks and fixes').  Details on this second method are given at the
   :ref:`end of this chapter <inputdata_native_datasets>`.

Tiers
-----

All observational datasets are grouped into in three tiers:

* **Tier 1**: obs4mips and ana4mips datasets. These datasets are publicly and freely available without any license restrictions. These datasets do not need any reformatting and can be used as is with ESMValTool.
* **Tier 2** other freely available datasets that are not obs4mips. There are no license restrictions. These datasets need to be reformatted to be used with ESMValTool ('CMORization', see above). 
* **Tier 3** restricted datasets. Datasets which require registration to be downloaded or that can only be obtained upon request from the respective authors. License restrictions do not allow us to redistribute Tier 3 datasets. The data have to be obtained and reformatted by the user ('CMORization', see above).

[!NOTE]
.. _tier3_note:
For some of the Tier 3 datasets, we obtained permission from the dataset providers to share the data among ESMValTool users on HPC systems. These Tier 3 datasets are marked with an asterisk in the table in section :ref:`supported datasets below<supported_datasets>`.

An overview of the Tier 2 and Tier 3 datasets for which a CMORizing script is available in ESMValTool v2.0 is given in section :ref:`supported datasets below<supported_datasets>`.

A collection of readily CMORized OBS and OBS6 datasets can be accessed directly on CEDA/JASMIN and DKRZ. At CEDA/JASMIN
OBS and OBS6 data is stored in the `esmeval` Group Workspace (GWS), and to be granted read (and execute) permissions to the
GWS, one must apply at https://accounts.jasmin.ac.uk/services/group_workspaces/esmeval/ ; after permission has been granted, the user
is encouraged to use the data locally, and not move it elsewhere, to minimize both data transfers and
stale disk usage; to note that Tier 3 data is subject to data protection restrictions; for further inquiries,
the GWS is administered by [Valeriu Predoi](mailto:valeriu.predoi@ncas.ac.uk).

Using a CMORizer script
-----------------------

ESMValTool comes with a set of CMORizers readily available.
The CMORizers are dataset-specific scripts that can be run once to generate a
local pool of CMOR-compliant data.
The necessary information to download and process the data is provided in the
header of each CMORizing script.
These scripts also serve as template to create new CMORizers for datasets not
yet included.
Note that datasets CMORized for ESMValTool v1 may not be working with v2, due
to the much stronger constraints on metadata set by the iris library.

ESMValTool provides the ``esmvaltool data`` command line tool, which can be
used to download and format datasets.

To list the available commands, run

.. code-block:: bash

    esmvaltool data --help

It is also possible to get help on specific commands, e.g.

.. code-block:: bash

    esmvaltool data download --help

The list of datasets supported by ESMValTool through a CMORizer script can be
obtained with:

.. code-block:: bash

    esmvaltool data list

Datasets for which auto-download is supported can be downloaded with:

.. code-block:: bash

    esmvaltool data download --config_file [CONFIG_FILE] [DATASET_LIST]

Note that all Tier3 and some Tier2 datasets for which auto-download is supported
will require an authentication. In such cases enter your credentials in your
``~/.netrc`` file as explained
`here <https://www.gnu.org/software/inetutils/manual/html_node/The-_002enetrc-file.html>`_.

An entry to the ``~/.netrc`` should look like:

.. code-block:: bash

    machine [server_name] login [user_name] password [password]

Make sure that the permissions of the ``~/.netrc`` file are set so only you and administrators
can read it, i.e.

.. code-block:: bash

    chmod 600 ~/.netrc
    ls -l ~/.netrc

The latter command should show ``-rw-------``.

For other datasets, downloading instructions can be obtained with:

.. code-block:: bash

    esmvaltool data info [DATASET]

To CMORize one or more datasets, run:

.. code-block:: bash

    esmvaltool data format --config_file [CONFIG_FILE] [DATASET_LIST]

The ``rootpath`` to the raw data to be CMORized must be specified in the
:ref:`configuration <esmvalcore:config_options>` as ``RAWOBS``.
Within this path, the data are expected to be organized in subdirectories
corresponding to the data tier: Tier2 for freely-available datasets (other than
obs4MIPs and ana4mips) and Tier3 for restricted datasets (i.e., dataset which
requires a registration to be retrieved or provided upon request to the
respective contact or PI).
The CMORization follows the `CMIP5 CMOR tables
<https://github.com/PCMDI/cmip5-cmor-tables>`_ or `CMIP6 CMOR tables
<https://github.com/PCMDI/cmip6-cmor-tables>`_ for the OBS and OBS6 projects
respectively.
The resulting output is saved in the output_dir, again following the Tier
structure.
The output file names follow the definition given in :ref:`config-developer
file <esmvalcore:config-developer>` for the ``OBS`` project:

.. code-block::

    [project]_[dataset]_[type]_[version]_[mip]_[short_name]_YYYYMM_YYYYMM.nc

where ``project`` may be OBS (CMIP5 format) or OBS6 (CMIP6 format), ``type``
may be ``sat`` (satellite data), ``reanaly`` (reanalysis data),
``ground`` (ground observations), ``clim`` (derived climatologies),
``campaign`` (aircraft campaign).

At the moment, ``esmvaltool data format`` supports Python and NCL scripts.

.. _supported_datasets:

Supported datasets for which a CMORizer script is available
-----------------------------------------------------------

A list of the datasets for which a CMORizers is available is provided in the following table.

.. tabularcolumns:: |p{3cm}|p{6cm}|p{3cm}|p{3cm}|

+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Dataset                      | Variables (MIP)                                                                                      | Tier | Script language |
+==============================+======================================================================================================+======+=================+
| AERONET                      | od440aer, od550aer, od870aer (AERmon)                                                                |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| AGCD                         | pr (Amon)                                                                                            |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ANU Climate                  | pr, tas, tasmin, tasmax (Amon)                                                                       |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| APHRO-MA                     | pr, tas (day), pr, tas (Amon)                                                                        |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| AURA-TES                     | tro3 (Amon)                                                                                          |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| BerkelyEarth                 | tas, tasa (Amon), sftlf (fx)                                                                         |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CALIPSO-GOCCP                | clcalipso (cfMon)                                                                                    |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CALIPSO-ICECLOUD* [#t3]_     | cli (AMon)                                                                                           |   3  | NCL             |
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
| CLARA-AVHRR                  | clt, clivi, clwvi, lwp (Amon)                                                                        |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CLOUDSAT-L2                  | clw, clivi, clwvi, lwp (Amon)                                                                        |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CMAP                         | pr (Amon)                                                                                            |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CowtanWay                    | tasa (Amon)                                                                                          |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CRU                          | tas, tasmin, tasmax,  pr, clt (Amon), evspsblpot (Emon)                                              |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| CT2019                       | co2s (Amon)                                                                                          |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Duveiller2018                | albDiffiTr13                                                                                         |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| E-OBS                        | tas, tasmin, tasmax, pr, psl (day, Amon)                                                             |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Eppley-VGPM-MODIS            | intpp (Omon)                                                                                         |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ERA5 [#note1]_               | cl, clt, evspsbl, evspsblpot, mrro, pr, prsn, ps, psl, ptype, rls, rlds, rlns, rlus [#note2]_, rsds, |   3  | n/a             |
|                              | rsns, rsus [#note2]_, rsdt, rss, uas, vas, tas, tasmax, tasmin, tdps, ts, tsn (E1hr/Amon), orog (fx) |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ERA5-Land [#note1]_          | pr                                                                                                   |   3  | n/a             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ERA-Interim                  | cl, cli, clivi, clt, clw, clwvi, evspsbl, hfds, hur, hus, lwp, orog, pr, prsn, prw, ps, psl, rlds,   |   3  | Python          |
|                              | rlut, rlutcs, rsds, rsdt, rss, rsut, rsutcs, sftlf, ta, tas, tasmax, tasmin, tauu, tauv, tdps, tos,  |      |                 |
|                              | ts, tsn, ua, uas, va, vas, wap, zg                                                                   |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ERA-Interim-Land             | sm (Lmon)                                                                                            |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-AEROSOL               | abs550aer, od550aer, od550aerStderr, od550lt1aer, od870aer, od870aerStderr (aero)                    |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-CLOUD                 | clivi, clt, cltStderr, clwvi, lwp, rlut, rlutcs, rsut, rsutcs, rsdt, rlus, rsus, rsuscs (Amon)       |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-FIRE                  | burntArea (Lmon)                                                                                     |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-LANDCOVER v1.6.1      | baresoilFrac, cropFrac, grassFrac, shrubFrac, treeFrac (Lmon)                                        |   2  | NCL             |
|                              |                                                                                                      |      | (CMORizer       |
|                              |                                                                                                      |      | available until |
|                              |                                                                                                      |      | ESMValTool      |
|                              |                                                                                                      |      | v2.11.0)        |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-LANDCOVER v2.0.8      | baresoilFrac, cropFrac, grassFrac, shrubFrac, treeFrac (Lmon, frequency=yr)                          |   2  | Python          |
|                              |                                                                                                      |      | (CMORizer       |
|                              |                                                                                                      |      | available since |
|                              |                                                                                                      |      | ESMValTool      |
|                              |                                                                                                      |      | v2.12.0)        |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-LST                   | ts (Amon)                                                                                            |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-OC                    | chl (Omon)                                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-OZONE                 | toz, tozStderr, tro3prof, tro3profStderr (Amon)                                                      |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-SEA-SURFACE-SALINITY  | sos (Omon)                                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-SOILMOISTURE          | sm (Eday, Lmon), smStderr (Eday)                                                                     |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-SST                   | ts, tsStderr (Amon)                                                                                  |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-WATERVAPOUR           | prw (Amon)                                                                                           |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESDC                         | tas, tasmax, tasmin (Amon)                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESRL                         | co2s (Amon)                                                                                          |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| FLUXCOM* [#t3]_              | gpp (Lmon)                                                                                           |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GCP2018                      | fgco2 (Omon [#note3]_), nbp (Lmon [#note3]_)                                                         |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GCP2020                      | fgco2 (Omon [#note3]_), nbp (Lmon [#note3]_)                                                         |   2  | Python          |
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
| GPCP-SG                      | pr (Amon)                                                                                            |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GRACE                        | lweGrace (Lmon)                                                                                      |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HadCRUT3                     | tas, tasa (Amon)                                                                                     |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HadCRUT4                     | tas, tasa (Amon), tasConf5, tasConf95                                                                |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HadCRUT5                     | tas, tasa (Amon)                                                                                     |   2  | Python          |
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
| JRA-25                       | clt, hus, prw, rlut, rlutcs, rsut, rsutcs (Amon)                                                     |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| JRA-55                       | cli, clivi, clw, clwvi, clt, prw, rlus, rlut, rlutcs, rsus, rsuscs, rsut, rsutcs, ta, tas, wap (Amon)|   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Kadow2020                    | tasa (Amon)                                                                                          |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| LAI3g                        | lai (Lmon)                                                                                           |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| LandFlux-EVAL                | et, etStderr (Lmon)                                                                                  |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Landschuetzer2016            | dpco2, fgco2, spco2 (Omon)                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Landschuetzer2020            | spco2 (Omon)                                                                                         |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MAC-LWP* [#t3]_              | lwp, lwpStderr (Amon)                                                                                |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MERRA                        | cli, clivi, clt, clw, clwvi, hur, hus, lwp, pr, prw, ps, psl, rlut, rlutcs, rsdt, rsut, rsutcs, ta,  |   3  | NCL             |
|                              | tas, ts, ua, va, wap, zg (Amon)                                                                      |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MERRA2* [#t3]_               | sm (Lmon)                                                                                            |   3  | Python          |
|                              | clt, pr, evspsbl, hfss, hfls, huss, prc, prsn, prw, ps, psl, rlds, rldscs, rlus, rlut, rlutcs, rsds, |      |                 |
|                              | rsdscs, rsdt, tas, tasmin, tasmax, tauu, tauv, ts, uas, vas, rsus, rsuscs, rsut, rsutcs, ta, ua, va, |      |                 |
|                              | tro3, zg, hus, wap, hur, cl, clw, cli, clwvi, clivi (Amon)                                           |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MLS-AURA* [#t3]_             | hur, hurStderr (day)                                                                                 |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MOBO-DIC-MPIM                | dissic (Omon)                                                                                        |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MOBO-DIC2004-2019            | dissic (Omon)                                                                                        |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MODIS                        | cliwi, clt, clwvi, iwpStderr, lwpStderr (Amon), od550aer (aero)                                      |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MSWEP [#note1]_              | pr                                                                                                   |   3  | n/a             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| MTE* [#t3]_                  | gpp, gppStderr (Lmon)                                                                                |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NCEP-NCAR-R1                 | clt, hur, hurs, hus, pr, prw, psl, rlut, rlutcs, rsut, rsutcs, sfcWind, ta, tas,                     |   2  | Python          |
|                              | tasmax, tasmin, ts, ua, va, wap, zg (Amon)                                                           |      |                 |
|                              | pr, rlut, ua, va (day)                                                                               |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NCEP-DOE-R2                  | clt, hur, prw, ta, wap, pr, tauu, tauv, tos (Amon)                                                   |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NDP                          | cVeg (Lmon)                                                                                          |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NIWA-BS* [#t3]_              | toz, tozStderr (Amon)                                                                                |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NOAA-CIRES-20CR-V2           | clt, clwvi, hus, prw, rlut, rsut, pr, tauu, tauv (Amon)                                              |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NOAA-CIRES-20CR-V3           | clt, clwvi, hus, prw, rlut, rlutcs, rsut, rsutcs (Amon)                                              |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NOAA-ERSSTv3b                | tos (Omon)                                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NOAA-ERSSTv5                 | tos (Omon)                                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NOAA-MBL-CH4                 | ch4s (Amon)                                                                                          |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NOAAGlobalTemp               | tasa (Amon)                                                                                          |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NSIDC-0116-[nh|sh] [#note4]_ | usi, vsi (day)                                                                                       |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| NSIDC-g02202-[sh]            | siconc (SImon)                                                                                       |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| OceanSODA-ETHZ               | areacello (Ofx), co3os, dissicos, fgco2, phos, spco2, talkos (Omon)                                  |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| OSI-450-[nh|sh]              | sic (OImon), sic (day)                                                                               |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| PATMOS-x                     | clt (Amon)                                                                                           |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| PERSIANN-CDR                 | pr (Amon), pr (day)                                                                                  |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| PHC                          | thetao, so (Omon [#note3]_)                                                                          |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| PIOMAS                       | sit (day)                                                                                            |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| REGEN                        | pr (day, Amon)                                                                                       |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| Scripps-CO2-KUM              | co2s (Amon)                                                                                          |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| TCOM-CH4                     | ch4 (Amon [#note3]_)                                                                                 |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| TCOM-N2O                     | n2o (Amon [#note3]_)                                                                                 |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| TROPFLUX                     | ts, rlut, rsds, tauu, tauv, hfds (Amon)                                                              |   2  | Python          |
|                              | tos (Omon)                                                                                           |      |                 |
|                              | hfls (Lmon, Amon)                                                                                    |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| UWisc* [#t3]_                | clwvi, lwpStderr (Amon)                                                                              |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| WFDE5                        | tas, pr (Amon, day)                                                                                  |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| WOA                          | thetao, so, tos, sos (Omon)                                                                          |   2  | Python          |
|                              | no3, o2, po4, si (Oyr)                                                                               |      |                 |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+

.. [#t3] We obtained permission from the dataset provider to share this dataset
         among ESMValTool users on HPC systems.

.. [#note1] CMORization is built into ESMValTool through the native6 project,
            so there is no separate CMORizer script.

.. [#note2] Derived on the fly from down & net radiation.

.. [#note3] The frequency of this variable differs from the one specified in
            the table. The correct entry that needs to be used in the recipe
            can be found in the corresponding section of `recipe_check_obs.yml
            <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/examples/recipe_check_obs.yml>`__.

.. [#note4] The cmoriser requires PROJ>=9.3. Previous version of PROJ will return an error:
            ``Internal Proj Error: proj_create: unhandled axis direction: UNKNOWN)``
            You can check the version of PROJ in your conda environment by running:
            ``conda list PROJ``.

.. _inputdata_native_datasets:

Datasets in native format
=========================

ESMValCore also provides support for some datasets in their native format.
In this case, the steps needed to reformat the data are executed as dataset
fixes during the execution of an ESMValTool recipe, as one of the first
preprocessor steps, see :ref:`fixing data <esmvalcore:fixing_data>`.
Compared to the workflow described above, this has the advantage that the user
does not need to store a duplicate (CMORized) copy of the data.
Instead, the CMORization is performed 'on the fly' when running a recipe.
Native datasets can be hosted either under a dedicated project (usually done
for native model output) or under project ``native6`` (usually done for native
reanalysis/observational products).
These projects are configured in the :ref:`config-developer file
<esmvalcore:configure_native_models>`.

A list of all currently supported native datasets is :ref:`provided here
<esmvalcore:read_native_datasets>`.
A detailed description of how to include new native datasets is given
:ref:`here <esmvalcore:add_new_fix_native_datasets>`.

To use this functionality, users need to provide a ``rootpath`` in the
:ref:`configuration <config_option_rootpath>` for the ``native6`` project data
and/or the dedicated project used for the native dataset, e.g., ``ICON``.
Then, in the recipe, they can refer to those projects.
For example:

.. code-block:: yaml

    datasets:
    - {project: native6, dataset: ERA5, type: reanaly, version: v1, tier: 3, start_year: 1990, end_year: 1990}
    - {project: ICON, dataset: ICON, exp: icon-2.6.1_atm_amip_R2B5_r1i1p1f1, mip: Amon, short_name: tas, start_year: 2000, end_year: 2014}

For project ``native6``, more examples can be found in the diagnostics
``ERA5_native6`` in the recipe `examples/recipe_check_obs.yml
<https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/examples/recipe_check_obs.yml>`_.
