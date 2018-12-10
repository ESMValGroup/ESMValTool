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

**Please**, if already not done, start your work with the ESMValTool and refer to the general template to prepare a first version of your report. Final templates with proper deliverable numbers are sent to ECV responsible to finalize the report.


1. Which ECVs will be evaluated, when and by who?
--------------------------------------------------
Below is a table of prioritized variables as reminder. *: some datasets are available on CDS. ^: no datasets available in CDS, please check regularly (latest 15/10/2018 for bold marked deadlines). Bold marked ECVs prioritized for October 2018.

+------------+-----------------------+-----+---------------+-----------------+
| ECV        | Products              | WP  | Who evaluates | Release         |
+============+=======================+=====+===============+=================+
| Surface    | ERA5*,                | 3   | NUIM/CNR      | **October**     |
| Temperature| E-OBS*                |     |               | **2018**        |
+------------+-----------------------+-----+---------------+-----------------+
| Wind       | ERA5*                 | 3   | NUIM          | April 2019      |
+------------+-----------------------+-----+---------------+-----------------+
| Humidity   | ERA5*                 | 3   | NUIM          | April 2019      |
+------------+-----------------------+-----+---------------+-----------------+
| Pressure   | ERA5*,                | 3   | NUIM          | April 2019      |
|            | E-OBS*                |     |               |                 |
+------------+-----------------------+-----+---------------+-----------------+
| Air        | ERA5*                 | 4   | CNR           | **October**     |
| Temperature|                       |     |               | **2018**        |
+------------+-----------------------+-----+---------------+-----------------+
| Ozone      | ESA-CCI Total column*,| 4   | ENEA          | April 2019      |
|            | ESA-CCI Tropospheric*,|     |               |                 |
|            | ESA-CCI Profiles*,    |     |               |                 |
|            | ERA5*                 |     |               |                 |
+------------+-----------------------+-----+---------------+-----------------+
| CO2        | ESA-CCI*              | 4   | DLR           | April 2019      |
+------------+-----------------------+-----+---------------+-----------------+
| Winds      | ERA5*                 | 4   | CNR           | April 2019      |
+------------+-----------------------+-----+---------------+-----------------+
| SST        | ESA-CCI OSTIA*,       | 5   | CNR           | **October 2018**|
|            | (ERA5)*               | 5   | CNR           | **October 2018**|
+------------+-----------------------+-----+---------------+-----------------+
| Sea Ice    | ESA-CCI OSTIA*,       | 5   | ULB           | April 2019      |
|            | (ERA5)*               | 5   | ULB           | April 2019      |
+------------+-----------------------+-----+---------------+-----------------+
| Sea Level, | ^                     | 5   | CNR, ENEA,    | April 2019      |
| Surface    |                       |     | CSIC          |                 |
| Properties |                       |     |               |                 |
+------------+-----------------------+-----+---------------+-----------------+
| Soil       | ERA5**,               | 6   | ETHz          | **October 2018**|
| Moisture   | (no ESA-CCI)**        |     |               |                 |
+------------+-----------------------+-----+---------------+-----------------+
| Albedo     | SPOT-VGT*             | 6   | LMU           | April 2019      |
+------------+-----------------------+-----+---------------+-----------------+
| Snow       | ^                     | 6   |               | April 2019      |
+------------+-----------------------+-----+---------------+-----------------+
| Glaciers   | Randolph*,            | 6   | VUB           | April 2019      |
|            | Elevation and mass*   |     |               |                 |
+------------+-----------------------+-----+---------------+-----------------+


2. How to access the datasets and make them compatible to ESMValTool
--------------------------------------------------------------------

For the first phase (until the ESMValTool is available through the CDS) it is necessary to perform the reporting off-line and download the datasets on a local machine or cluster where ESMValTool is run. 

The extraction of the datasets is done:
 1. Through the CDS itself 
 2. Using the CDS-Toolbox
 3. Using the CDS API (`<https://pypi.org/project/cdsapi/>`_ with an example script `here <https://github.com/bascrezee/c3s_tools/blob/master/retrieve_era5.py>`_)

We advise to make use of the second option since the toolbox treats the variables following the CMOR tables and allows the necessary sub-setting and resolution degradation. 

The procedure is as follows:
  * Go to the CDS and choose the dataset you want: https://climate.copernicus.eu/climate-data-store
  *	Select a small subset example
  *	Perform a download after accepting the agreement
  *	On the bottom of the page tick: Show Toolbox Request
  *	Open the toolbox from the options at the top of the page
  *	In the toolbox window open the sample code 01 retrieve data 
  *	Load it
  *	Make a copy of it
  *	Replace the command after *"data=“* with the Toolbox request script
  *	Run it
  *	Download the data from the right column of the screen

Below are two examples of the resulting script for 1D and 3D data.

In order to have data compliant with the CDM of the Toolbox, an operation on the data needs to be performed. The simple retrieval does not require the data are CDM compliant. Whatever is the initial format of the data (grib, netcdf, zip...) the operation assures the retrieved dataset is netcdf and CDM compliant, e.g. ct.climate.daily_mean.


