.. _documentation:

Scientific documentation of a diagnostic script or metrics set
**************************************************************

An important part of the implementation of a new diagnostic script is the documentation of the diagnostic itself as well as the documentation of the observational data sets used on the ESMValTool development team wiki (see section TODO:ref 12.4 for details on the wiki, see also section TODO:ref 7 for general guidelines on documentation).
The former should comply with the standard template for new diagnostics (see section TODO:ref 8.1 below) and the latter should include instructions how to download the observational data and, if necessary, scripts to convert it to the format required in ESMValTool, see section TODO:ref 8.2 below.

Standard template
=================

When implementing a new diagnostic script or metrics set, it should be documented on the ESMValTool development team wiki in the OPEN or PRIVATE Github repository (see section TODO:ref 12.4) starting from the standard template given below:

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
|                           | (e.g., stocker12jgr,stocker12sci1,stocker12sci2).                        |
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
| . . . Insert text here                                                                               |
| ``2.`` Available Diagnostics                                                                         |
| . . . Insert text here                                                                               |
| ``3.`` Specific Routines                                                                             |
| . . . Contains a description of specific routines being developed for the given diagnostic that helps|
| . . . to identify common code (which should then go in the ``lib/``)                                 |
| ``4.`` Observations and Scripts (also see Model and observational data below)                        |
| . . . Insert text here                                                                               |
| ``5.`` Test Cases (see also Automated testing, section TODO:ref 7.9)                                 |
| . . . Insert text here                                                                               |
| ``6.`` References                                                                                    |
| . . . REF1                                                                                           |
| . . . REF2                                                                                           |
| . . . etc.                                                                                           |
| ``7.`` Sample Plots                                                                                  |
+---------------------------+--------------------------------------------------------------------------+

Model and observational data
============================

Overview
--------

Standard header for the reformatting routines for observational data
--------------------------------------------------------------------

