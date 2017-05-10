.. _core_team:

The ESMValTool core development team
************************************

.. _core_dev_team:

Main contacts
=============

A mailing list has been set up for all general and technical questions on the ESMValTool such as, for instance,
questions on installation, application or development. You are encouraged to subscribe to the ESMValTool user
mailing list by sending an email to Listserv@dlr.de with the following text:

.. centered:: *subscribe ESMValTool-Usr*

Core development team
=====================

* **Veronika Eyring (DLR, Germany)**, veronika.eyring at dlr.de

  ESMValTool Core PI and Developer: contact for requests to use the ESMValTool and for collaboration with the
  development team, access to the PRIVATE GitHub repository (see section 12 TODOref)

* **Axel Lauer (DLR, Germany)**, axel.lauer at dlr.de

  ESMValTool Core Developer: contact for technical and scientific questions, access to the PRIVATE GitHub
  repository (see section 12 TODOref)

* **Mattia Righi (DLR, Germany)**, mattia.righi at dlr.de

  ESMValTool Core Developer: contact for technical and scientific questions

* **Björn Brötz (DLR, Germany)**, bjoern.broetz at dlr.de

  ESMValTool Core Developer: contact for technical questions

Contacts for specific diagnostic sets are the respective authors, as listed in the corresponding diagnostic
documentation (including the ESMValTool development team wiki) and in the source code.

Pull requests
=============

This section describes the general workflow of how new diagnostics are integrated into the ESMValTool and the
responsibilities of the developer contribution to the ESMValTool. *For technical details on how to work with the
version controlled repository for the ESMValTool source code see section 12 TODOref.*

Workflow core development team
------------------------------

The following workflow followed by the ESMValTool core development team takes place whenever a developer
requests integration of a diagnostics set into the *development branch* (see section 12 TODOref):

#. Check that the developer submits a **standard namelist** that calls a set of diagnostics / metrics

#. Check that the related documentation on the ESMValTool development team wiki is compliant with **documentation templates for diagnostics and metrics sets** (see section 7.1 TODOref).

#. Check that the code follows **coding rules and standard** (see section 7.5 TODOref).

#. Check that a **namelist** is provided for **automated testing** that runs on 2-3 models and a small set of observations/reference model/idealized data.

   * Verify that such a **reduced and small set of observations/reference model/idealized data** is delivered for each diagnostic that is called by the standard namelist.
   * Verify that an **example plot + netCDF for automated testing** created with this reduced data set is provided for each diagnostic that is called by the standard namelist as a reference.

#. Check that also the **full set of observations** is provided that allows a sophisticated scientific application of the (full) standard namelist.

#. Check that the **observations are documented** on the ESMValTool development team wiki and that a reformat routine is available in case the original source is not in CMOR standard.

#. Run the **automated testing** with all available diagnostics.

#. Iterate with developer(s) on points 1-7 until the above items are fulfilled and the reference plots for all standard namelists included in the *DEVELOPMENT BRANCH* (see section 12 TODOref) can be reproduced.

Responsibilities of ESMValTool developers
-----------------------------------------

1. Accept the **ESMValTool license agreement / terms of use** (see http://www.esmvaltool.org/license.html).

2. Provide documentation on the ESMValTool developmet team wiki that is compliant with **documentation templates for diagnostics and metrics sets**.

3. Provide well documented code that follows the **coding rules and standards**.

4. **For each pull request** to implement a diagnostic set into the *DEVELOPMENT BRANCH* (see section 12 TODOref).

   **Scientific analysis**

   * Provide the code **for all diagnostics and metrics** that are called.
   * Standard namelist running on (if possible) all CMIP5 models and corresponding plots that are produced (for the
   * wiki and the user's guide).
   * Provide the **full set of observations** that allows a sophisticated scientific application of the full standard namelist list (indicate source and if applicable license issues).
   * Provide **documentation for the observations** on the ESMValTool development team wiki and a reformat routine if the original source does not follow the CMOR standard.

   **Automated testing (see section 7.9 TODOref)**

   * Provide the **code for automated testing** for the diagnostic set that should be integrated into the
   * *DEVELOPMENT BRANCH* (see section 12 TODOref).
   * Provide a **namelist for automated testing**.
   * Provide a **reduced and small set of observations/reference model/idealized data** for each diagnostic that is called by the testing namelist.
   * Provide **NetCDF + example plots for automated testing** based on the reduced data set and the standard namelist as a reference.

5. **Name a contact person** providing (scientific) support for your diagnostics.

