Guidance document to prepare the Single Product Quality Brief (SPQB) in ESMValTool v2
=====================================================================================

0. Preparation of the Deliverable Report with the beta version of the SPQB in ESMValTool
----------------------------------------------------------------------------------------
The scope of this document is to provide a step-by-step guidance on the preparation of the Single Products Quality Briefs (SPQB) with the ESMValTool. This document is evolutive. The current beta version of the SPQB namelist does not implement the full functionalities. Therefore, this guidance will be updated for future reporting deliverables. The reports need to be produced in the following way:

The reports will be produced with the following procedure:
1.	Download the dataset from the CDS and ensure that the format is ESMValTool compliant (CMOR) (Section 2).
2.	Prepare the (semi-)automatic sub-reports with the ESMValTool (Section 3).
3.	Follow the report description (Section 4) and make use of the template that can be found in the file C3S_D511_template_SPQB_12102018 with the specific guidelines on how to compile it. **Note**: This is an example for SSTs. Adjust these templates to your analyzed data set following the yellow guidance attaching the subreports.
4.	Put all your report in the official deliverable document.

**Please note**: The content of the deliverable is specific to each ECV and Evaluators can make their choices in terms of plots / content of their summary, keeping in mind that reports should be readable and contain synthetic information following the basic content described in the C3S_D511_template_SPQB_12102018 template.

**Please**, if already not done, start your work with the ESMValTool and refer to the general template to prepare first version of your report. Final templates with proper deliverable numbers are sent to ECV responsible to finalize the report.


1. Which ECVs shall be evaluated, when and by who?
--------------------------------------------------
Below is a table of prioritized variables as reminder. Green: some datasets are available on CDS. Yellow: no datasets available in CDS, please check regularly (latest 15/10/2018 for red marked cells). Red ECVs prioritized for October 2018.

+------------+------------+-----+---------------+-------------+
| ECV        | Products   | WP  | Who evaluates | Release     |
+============+============+=====+===============+=============+
| Surface    | ERA5,      | 3   | NUIM/CNR      | October 2018|
| Temperature| E-OBS      |     |               |             |
+------------+------------+-----+---------------+-------------+
| Wind       | ERA5       | 3   | NUIM          | April 2019  |
+------------+------------+-----+---------------+-------------+
| Humidity   | ERA5       | 3   | NUIM          | April 2019  |
+------------+------------+-----+---------------+-------------+
| Pressure   | ERA5,      | 3   | NUIM          | April 2018  |
|            | E-OBS      | 3   |               |             |
+------------+------------+-----+---------------+-------------+
