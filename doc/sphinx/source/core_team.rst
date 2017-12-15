.. _tasks_and_responsibilities:

Tasks and responsibilities
**************************

Mailing list
============

A mailing list has been set up for all general and technical questions on the ESMValTool such as, for instance,
questions on installation, application or development. You are encouraged to subscribe to the ESMValTool user
mailing list by sending an email to Listserv@dlr.de with the following text:

.. centered:: *subscribe ESMValTool-Usr*

.. _core_dev_team:

The ESMValTool core development team
====================================

* Deutsches Zentrum für Luft- und Raumfahrt (DLR), Institut für Physik der Atmosphäre, Germany (PI)

  ESMValTool Core PI and Developer: contact for requests to use the ESMValTool and for collaboration with the
  development team, access to the PRIVATE GitHub repository (see Part :numref:`annex_b`)

* Alfred-Wegener-Institute Bremerhaven (AWI), Germany (overseeing EU Horizon 2020 APPLICATE and TRR181 ESMValTool work)
* Barcelona Computing Center (BSC), Spain (overseeing EU Horizon 2020 PRIMAVERA ESMValTool work)
* Ludwig Maximilian University of Munich, Germany (overseeing EU Horizon 2020 CRESCENDO ESMValTool work)
* Netherlands e-Science Center (NLeSC), The Netherlands
* University of Reading, United Kingdom

Contacts for specific diagnostic sets are the respective authors, as listed in the corresponding diagnostic
documentation and in the source code.

Pull requests
=============

The development of the ESMValTool is a community effort and involves a number of tasks to allow for a smooth integration of new code into the ESMValTool (*for technical details on how to work with the version controlled repository for the ESMValTool source code see Part* :numref:`annex_b`). These tasks are divided into:

1) responsibilities for the ESMValTool Core Development Team that oversees the development and prepares releases of the ESMValTool and
2) responsibilities for the individual developers contributing code for new diagnostics or metrics.

ESMValTool Core Development Team
--------------------------------

#. Responsible for the *MASTER BRANCH* and *DEVELOPMENT BRANCH* including pull requests by developers (see below)
#. Oversee collection of observations that are required to actually run all diagnostics that are in the *DEVELOPMENT BRANCH*
#. Responsible for general quality control, maintenance, and documentation (including compliance with coding rules)
#. Oversee collection of testing data (observations/reference model/idealized data) that are required for automated testing (host a local copy)
#. Responsible for new releases of the ESMValTool

ESMValTool developers
---------------------

#. Accept the **ESMValTool license and ESMValTool Development terms of use** (see http://www.esmvaltool.org/license.html).
#. Provide documentation that is compliant with **documentation templates for diagnostics and metrics sets** (see Sections :numref:`std_namelist` and :numref:`std_diag`).
#. Provide well documented code that follows the coding rules and standards. It is recommended to use an existing diagnostic in the target language as a template for the development. Check with the Core Development Team if unsure which diagnostic to use as a template.
#. For each pull request to implement a diagnostic set into the *DEVELOPMENT BRANCH*:

   **Scientific analysis**

   * Provide the code for **all diagnostics and metrics** that are called; all source code files should include a standard header as defined in Section :numref:`std_diag`.
   * Standard namelist (with settings so that it runs on, if possible, all CMIP5 models) and corresponding plots (for documentation), the namelist should include a standard header as defined in Section :numref:`std_namelist`.
   * Provide documentation of the diagnostic following the structure defined in the template doc/sphinx/source/namelists/namelist_template.rst (reStructuredText (RST) syntax) including example images (copy to doc/sphinx/source/namelists/figures/<namelist>/).
   * Provide the **full set of observations** that allows a scientific application of the full standard namelist list (indicate source and if applicable license issues).
   * Provide that **observations are documented**. If new observations are introduced, add an entry to the table listing the available observations (:numref:`tab_obs_data`) in doc/sphinx/source/running.rst) and that a reformat routine is available if the original source does not follow the CMOR standard.

   **Automated testing (see Section :numref:`auto_test`)**

   * Provide the code for automated testing for the diagnostic set that should be integrated into the *DEVELOPMENT BRANCH* (see Part :numref:`annex_b`).
   * Provide a **namelist for automated testing**.
   * Provide a **reduced and small set of observations/reference model/idealized data** for each diagnostic that is called by the testing namelist.
   * Provide **netCDF output + example plots for automated testing** based on the reduced dataset and the standard namelist as a reference.

#. **Name a contact person** providing (scientific) support for your diagnostics.

