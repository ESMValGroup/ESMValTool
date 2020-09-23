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
`IS-ENES3 Trans-national Access call <https://portal.enes.org/data/data-metadata-service/analysis-platforms>`__.

If the options above are not available to you, ESMValTool also offers features
to make it easier to download the data.

Models
======

ESMValTool will look for existing data in the directories specified in the
user configuration file. Alternatively, it can use an external
tool called `Synda <http://prodiguer.github.io/synda/index.html>`__. If you
do not have access to a compute cluster with the data already mounted, this is
the recommended approach for first-time users to obtain some data for
running ESMValTool.

Installing Synda for use from ESMValTool
----------------------------------------
Here, we describe the basic steps to configure EMSValTool so it can use Synda
to download CMIP6 or CMIP5 model data.

To install Synda, follow the steps listed in the
`Synda installation documentation <http://prodiguer.github.io/synda/sdt/conda_install.html>`__.
(This description assumes that Synda is installed using Conda.)
As the last step, Synda will ask to set your openID credentials.
Therefore, you'll need to create an account on an ESGF node, e.g.
`the ESGF node at Lawrence Livermore National Laboratory <https://esgf-node.llnl.gov/projects/esgf-llnl/>`__
and join a Data Access Control Group, e.g. 'CMIP5 Research'. For more information, see
`the ESGF user guide <https://esgf.github.io/esgf-user-support/user_guide.html>`__.

Once you have set up Synda, you'll need to configure ESMValTool to find
your Synda installation. Note that it is not possible to combine the two in a
single
`conda environment <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments>`__,
because Synda requires python 2 and ESMValTool requires Python 3.
Running

.. code-block:: bash

    which synda

on the command line, while your synda environment is active, will print its location.
To make the ``synda`` program usable from ESMValTool we suggest
creating a directory

.. code-block:: bash

    mkdir ~/bin

and appending that folder to your ``PATH`` environment variable,
e.g. by adding the following line to your ``~/.bashrc`` file:

.. code-block:: bash

    PATH=$PATH:$HOME/bin

Finally, in the new bin folder, make a link to synda:

.. code-block:: bash

    ln -s /path/to/conda/envs/synda/bin/synda ~/bin/synda

Now, ESMValTool should be able to find your Synda installation. First time
users can now continue with :ref:`Running ESMValTool <running>`.

Observations
============

Observational and reanalysis products in the standard CF/CMOR format used in CMIP and required by the ESMValTool are available via the obs4mips and ana4mips projects at the ESGF (e.g., https://esgf-data.dkrz.de/projects/esgf-dkrz/). Their use is strongly recommended, when possible.

Other datasets not available in these archives can be obtained by the user from the respective sources and reformatted to the CF/CMOR standard. ESMValTool currently support two ways to perform this reformatting (aka 'cmorization'). The first is to use a cmorizer script to generate a local pool of reformatted data that can readily be used by the ESMValTool. The second way is to implement specific 'fixes' for your dataset. In that case, the reformatting is performed 'on the fly' during the execution of an ESMValTool recipe (note that one of the first preprocessor tasks is 'cmor checks and fixes'). Below, both methods are explained in more detail.

Using a cmorizer script
-----------------------

ESMValTool comes with a set of cmorizers readily available. The cmorizers are dataset-specific scripts that can be run once to generate a local pool of CMOR-compliant data. The necessary information to download and process the data is provided in the header of each cmorizing script. These scripts also serve as template to create new cmorizers for datasets not yet included. Note that datasets cmorized for ESMValTool v1 may not be working with v2, due to the much stronger constraints on metadata set by the iris library.

To cmorize one or more datasets, run:

.. code-block:: bash

    cmorize_obs -c [CONFIG_FILE] -o [DATASET_LIST]

The path to the raw data to be cmorized must be specified in the CONFIG_FILE as RAWOBS. Within this path, the data are expected to be organized in subdirectories corresponding to the data tier: Tier2 for freely-available datasets (other than obs4mips and ana4mips) and Tier3 for restricted datasets (i.e., dataset which requires a registration to be retrieved or provided upon request to the respective contact or PI). The cmorization follows the CMIP5 CMOR tables. The resulting output is saved in the output_dir, again following the Tier structure. The output file names follow the definition given in ``config-developer.yml`` for the ``OBS`` project: ``OBS_[dataset]_[type]_[version]_[mip]_[short_name]_YYYYMM_YYYYMM.nc``, where ``type`` may be ``sat`` (satellite data), ``reanaly`` (reanalysis data), ``ground`` (ground observations), ``clim`` (derived climatologies), ``campaign`` (aircraft campaign).

At the moment, cmorize_obs supports Python and NCL scripts.

.. _cmorization_as_fix:

Cmorization as a fix
--------------------
As of early 2020, ESMValTool also provides (limited) support for data in their native format. In this case, the steps needed to reformat the data are executed as datasets fixes during the execution of an ESMValTool recipe, as one of the first preprocessor steps. Compared to the workflow described above, this has the advantage that the user does not need to store a duplicate (cmorized) copy of the data. Instead, the cmorization is performed 'on the fly' when running a recipe. ERA5 is the first dataset for which this 'cmorization on the fly' is supported.

To use this functionality, users need to provide a path for the ``native6`` project data in the :ref:`user configuration file<config-user>`. Then, in the recipe, they can refer to the native6 project, like so:

.. code-block:: yaml

    datasets:
    - {dataset: ERA5, project: native6, type: reanaly, version: '1', tier: 3, start_year: 1990, end_year: 1990}

Currently, the native6 project only supports ERA5 data in the format defined in the `config-developer file <https://github.com/ESMValGroup/ESMValCore/blob/a9312a7d5be4fa3aac55c0b2ef089c6b4e1a61a9/esmvalcore/config-developer.yml#L191-L201>`_. The filenames correspond to the default filenames from `era5cli <https://era5cli.readthedocs.io>`_ To support other datasets as well, we need to make it possible to have a dataset specific DRS. This is still on the horizon.

While it is not strictly necessary, it may still be useful in some cases to create a local pool of cmorized observations. This can be achieved by using a cmorizer *recipe*. For an example, see `recipe_era5.yml <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/recipes/cmorizers/recipe_era5.yml>`_. This recipe reads native, hourly ERA5 data, performs a daily aggregation preprocessor, and then calls a diagnostic that operates on the data. In this example, the diagnostic renames the data to the standard OBS6 format. The output are thus daily, cmorized ERA5 data, that can be used through the OBS6 project. As such, this example recipe does exactly the same as the cmorizer scripts described above: create a local pool of cmorized data. The advantage, in this case, is that the daily aggregation is performed only once, which can save a lot of time and compute if it is used often.

The example cmorizer recipe can be run like any other ESMValTool recipe:

.. code-block:: bash

    esmvaltool -c [CONFIG_FILE] cmorizers/recipe_era5.yml

(Note that the ``recipe_era5.yml`` adds the next day of the new year to the input data. This is because one of the fixes needed for the ERA5 data is to shift (some of) the data half an hour back in time, resulting in a missing record on the last day of the year.)

To add support for new variables using this method, one needs to add dataset-specific fixes to the ESMValCore. For more information about fixes, see: `fixing data <https://docs.esmvaltool.org/projects/esmvalcore/en/latest/develop/fixing_data.html#fixing-data>`_.


Supported datasets
------------------
A list of the datasets for which a cmorizers is available is provided in the following table.

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
| ERA5 [*]_                    | clt, evspsbl, evspsblpot, mrro, pr, prsn, ps, psl, ptype, rls, rlds, rsds, rsdt, rss, uas, vas, tas, |   3  | n/a             |
|                              | tasmax, tasmin, tdps, ts, tsn (E1hr/Amon), orog (fx)                                                 |      |                 |
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
| ESACCI-CLOUD                 | clivi, clt, cltStderr, clwvi (Amon)                                                                  |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-FIRE                  | burntArea (Lmon)                                                                                     |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-LANDCOVER             | baresoilFrac, cropFrac, grassFrac, shrubFrac, treeFrac (Lmon)                                        |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-OC                    | chl (Omon)                                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-OZONE                 | toz, tozStderr, tro3prof, tro3profStderr (Amon)                                                      |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-SOILMOISTURE          | dos, dosStderr, sm, smStderr (Lmon)                                                                  |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESACCI-SST                   | ts, tsStderr (Amon)                                                                                  |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| ESRL                         | co2s (Amon)                                                                                          |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| FLUXCOM                      | gpp (Lmon)                                                                                           |   3  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GCP                          | fgco2 (Omon), nbp (Lmon)                                                                             |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GHCN                         | pr (Amon)                                                                                            |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GHCN-CAMS                    | tas (Amon)                                                                                           |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GISTEMP                      | tasa (Amon)                                                                                          |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| GPCC                         | pr (Amon)                                                                                            |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HadCRUT3                     | tas, tasa (Amon)                                                                                     |   2  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| HadCRUT4                     | tas, tasa (Amon)                                                                                     |   2  | NCL             |
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
| UWisc                        | clwvi, lwpStderr (Amon)                                                                              |   3  | NCL             |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+
| WOA                          | no3, o2, po4, si (Oyr), so, thetao (Omon)                                                            |   2  | Python          |
+------------------------------+------------------------------------------------------------------------------------------------------+------+-----------------+

.. [*] ERA5 cmorization is built into ESMValTool through the native6 project, so there is no separate cmorizer script.
