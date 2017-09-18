.. _documentation:

Scientific documentation of a diagnostic script or metrics set
**************************************************************

An important part of the implementation of a new diagnostic script is the documentation of the diagnostic itself as well as the documentation of the observational data sets used on the ESMValTool development team wiki (see :numref:`wiki` for details on the wiki, see also :numref:`writing` for general guidelines on documentation).
The former should comply with the standard template for new diagnostics (see section :numref:`std_diag` below) and the latter should include instructions how to download the observational data and, if necessary, scripts to convert it to the format required in ESMValTool, see section :numref:`mod_obs_data` below.

.. _std_diag:

Standard template
=================

When implementing a new diagnostic script or metrics set, it should be documented on the ESMValTool development team wiki in the OPEN or PRIVATE Github repository (see :numref:`wiki`) starting from the standard template given below:

.. tabularcolumns:: |p{4cm}|p{11cm}|

+---------------------------+--------------------------------------------------------------------------+
| **Title of diagnostic/performance metrict**                                                          |
+===========================+==========================================================================+
| **Developers**            | first name surname 1 (DLR tag 1), first name surname 2 (DLR tag 2), etc. |
+---------------------------+--------------------------------------------------------------------------+
| **Contributors**          | first name surname 1 (DLR tag 1), first name surname 2 (DLR tag 2), etc. |
+---------------------------+--------------------------------------------------------------------------+
| **Date of documentation** | yyyy-mm-dd                                                               |
+---------------------------+--------------------------------------------------------------------------+
| **Name of standard**      | For Xyz the following naming conventions is used:                        |
| **namelist (XyZ)**        |                                                                          |
|                           | For papers:                                                              |
|                           |                                                                          |
|                           | XyZ=SurnameYearJournalabbreviation                                       |
|                           |                                                                          |
|                           | (e.g., stocker12jgr, stocker12sci1, stocker12sci2).                      |
|                           |                                                                          |
|                           | For copies of reports that are not publicly available:                   |
|                           |                                                                          |
|                           | XyZ=OrgYearTitleabbrev                                                   |
|                           |                                                                          |
|                           | (e.g., unep10water, unep11gap,roysoc09geoengineering).                   |
|                           |                                                                          |
|                           |                                                                          |
|                           | For grouped set of diagnostics and performance metrics that do not follow|
|                           | a published paper or report:                                             |
|                           |                                                                          |
|                           | an intuitive name that describes the science theme                       |
|                           |                                                                          |
|                           | (e.g., XyZ=aerosol, MyDiag, SAMonsoon, SeaIce).                          |
+---------------------------+--------------------------------------------------------------------------+
| **User settings**         | list of all settings that have to be checked/changed by a user in order  |
|                           | to run the diagnostic (e.g., pathnames, configuration files, color       |
|                           | tables, supported model names, etc.)                                     |
+---------------------------+--------------------------------------------------------------------------+
| **Brief summary**         | 1-3 sentence summary                                                     |
+---------------------------+--------------------------------------------------------------------------+
| **Status**                | Planned/work in progress/finished and integrated into development branch |
+---------------------------+--------------------------------------------------------------------------+
| **CMOR variable name**    | e.g., tas (atmos, monthly mean, longitude latitude plevs time)           |
| **(realm, frequency,**    |                                                                          |
| **dimension)**            |                                                                          |
+---------------------------+--------------------------------------------------------------------------+
| **Link to git repository**| e.g., https://github.com/axel-lauer/ESMValTool/tree/cloud                |
| **(feature branch)**      |                                                                          |
+---------------------------+--------------------------------------------------------------------------+
| ``1.`` Overview                                                                                      |
| ... Insert text here                                                                                 |
| ``2.`` Available Diagnostics                                                                         |
| ... Insert text here                                                                                 |
| ``3.`` Specific Routines                                                                             |
| ... Contains a description of specific routines being developed for the given diagnostic that helps  |
| ... to identify common code (which should then go in the ``lib/``)                                   |
| ``4.`` Observations and Scripts (also see Model and observational data below)                        |
| ... Insert text here                                                                                 |
| ``5.`` Test Cases (see also Automated testing, section :numref:`auto_test`)                          |
| ... Insert text here                                                                                 |
| ``6.`` References                                                                                    |
| ... REF1                                                                                             |
| ... REF2                                                                                             |
| ... etc.                                                                                             |
| ``7.`` Sample Plots                                                                                  |
| ... Please insert sample plots for all plot types produced by the namelist                           |
+---------------------------+--------------------------------------------------------------------------+


.. _mod_obs_data:

Model and observational data
============================

Overview
--------

When possible, observations from the obs4MIPs/ana4MIPs archives are used in the model evaluation (see :numref:`diag_avail`).
These data are freely available from the ESGF in the same format as the CMIP simulations and can be directly used in the ESMValTool using the obs4mips or ana4mips class in the namelist (see also section :numref:`mod_obs_run`).

Important links

https://www.earthsystemcog.org/projects/obs4mips/satellite_data_products

Nightly scan across nodes

https://www.earthsystemcog.org/search/obs4mips/?template=obs4mips&limit=200

Observational data sets not available in these archives need to be reformatted according to the CF/CMOR standard before they can be used.
In this case a reference to the official URL is provided such that a user can get the latest version of the data set as well as a description and a script how to convert the data set to the format required by the ESMValTool. These conversion scripts are collected in *reformat_scripts/obs/reformat_obs_<NAME>.ncl*.
The reformatting routines must be documented with a standard header providing all information required to retrieve and process the data, as well as their availability (Tier 1, Tier 2, or Tier 3).

All observations are tiered as follows:

    * Tier 1: data sets from the obs4MIPs and ana4MIPs archives
    * Tier 2: other freely available data sets
    * Tier 3: restricted data sets (e.g., license agreement required)

For Tier 2 and 3 data, the developer shall also provide links and helper scripts through the reformatting routines, following the template for the standard header described in section for the reformatting routines.
An example can be found here:

.. centered::
    *reformat_scripts/obs/reformat_obs_AURA-MLS-OMI.ncl*.

An overview on the available reformatting scripts for Tier 2 and 3 data is given in :numref:`tab_obs_data`.
The reformatted observational data (Tier 2 and Tier 3) must be named according to the OBS class defintion, which considers the following naming convention:

.. centered::
    OBS_<name>_<case>_<ensemble>_<field>_<variable>_<YYY1M1>-<YYY2M2>.nc

where:

<name> is the name of the satellite, instrument, campaign, network, model, etc. (e.g., ERA-Interim, AERONET, AURA-MLS-OMI, etc.)

<case> is the observation type (insitu, ground, sat, reanaly, campaign, etc.)

<ensemble> is the version number, processing level or station code (for ground-based networks), use 1 if not available.

It is also possible to split the output in multiple files, like in the CMIP5 class, e.g. _200101-200512.nc, 200601_201012.nc, 201101-201512.nc, etc. This is particularly useful for daily data, which are usually too large to be collected in a single file covering the whole time period.

Standard header for the reformatting routines for observational data
--------------------------------------------------------------------

This is a template of the standard header for the reformat_obs routines.
The parts in red are the ones to be modified by the author.
The modification history is given in reverse chronological order (i.e., most recent on top) and the last entry always contains the written statement.
The author of each entry in the modification history shall be indicated with the author tag, as given in the master reference file (*doc/MASTER_authors-refs-acknow.txt*), e.g., A_surn_na = surname, name.
All lines should be limited to a maximum of 79 characters.

.. code-block:: ncl

    ;;#############################################################################
    ;; REFORMAT SCRIPT FOR THE [OBSERVATION NAME] OBSERVATIONAL DATA
    ;;#############################################################################
    ;;
    ;; Tier
    ;;    [Information on data availability, possible options are:]
    ;;    Tier 1: obs4MIPs or ana4MIPs
    ;;    Tier 2: other freely-available data set
    ;;    Tier 3: restricted data set
    ;;
    ;; Source
    ;;    [URL to the data source or the reference]
    ;;
    ;; Last access
    ;;    [YYYYMMDD]
    ;;
    ;; Download and processing instructions
    ;;    [Short explanation on how to download and process the data]
    ;;
    ;; Caveats
    ;;    [List possible caveats or limitations of this script]
    ;;    [Features to-be-implemented shall also be mentioned here]
    ;;
    ;; Modification history
    ;;    [YYYYMMDD-A_xxxx_yy: extended...]
    ;;    [YYYYMMDD-A_xxxx_yy: written.]
    ;;
    ;; #############################################################################

    load ...
    load ...

