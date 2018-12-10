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

*  Go to the CDS and choose the dataset you want: `<https://climate.copernicus.eu/climate-data-store>`_
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

2.1 Extract 1D from ERA-5
^^^^^^^^^^^^^^^^^^^^^^^^^^

See example code on how to extract 1D ERA5 data from the CDS via the Toolbox `here <https://github.com/ESMValGroup/ESMValTool-private/blob/C3S_511_v2/esmvaltool/toolbox_1D_ERA5.py>`_


**2.2 Extract 1-month daily 3D temperature from ERA-5**

See example code on how to extract 3D ERA5 data from the CDS via the Toolbox `here <https://github.com/ESMValGroup/ESMValTool-private/blob/C3S_511_v2/esmvaltool/toolbox_3D_ERA5.py>`_

**Please note:** The last two lines perform a daily mean and produce a netcdf output file.


**2.3 Data Limitations**

Several datasets may need to be reduced in size due to limitations of memory space for operations of ESMValTool. The maximum size of the dataset is determined by the machine which is used to run the ESMValTool. Rough guidelines for planning the dataset size would be:

*	Machine storage available: ~dataset size x 3
*	Machine memory available: ~dataset size x 3
*	Machine minimum CPU requirements: single core

We strongly suggest not exceeding these limits.
This might require downscaling the datasets in temporal or spatial resolution for the reports. We advise to download parts of the required dataset from CDS to estimate the full size your data (e.g. size(one month) x 12 x number of years = full size). Then estimate which spatial aggregation (0.25x0.25 => 0.5x0.5 makes 2x2=4 times less space in memory or storage) or temporal aggregation (from daily to monthly data this makes ~30 times less space) is needed to a) make it possible for you to produce the reports and b) have a minimum of alteration (e.g. if spatial aggregation is done with averaging or nearest neighbor depends on your dataset, or if you have not enough space by a factor of 5, monthly means are not necessary). The exact extend of the downscaling (which coordinates are downscaled) is ECV dependent, and therefore has to be decided by the ECV expert user.
We suggest the following 3 subsets as an approach for reporting 4D variables, if processing time is available:

*	Full time resolution on one chosen level and reduced spatial resolution
*	Full vertical resolution and reduced spatial and temporal resolution
*	Full spatial resolution and reduced vertical and temporal resolution

**Please note:** There is an ongoing discussion for homogenization of this approach for ERA5 data.


**2.4 If your dataset is not CMORized…**

Data is required to be adherent to CMOR tables to be treated by the ESMValTool. The tool will crash if this is not the case. 
Reference to CMOR can be found in (`<https://cmor.llnl.gov/>`_). Please note: we are not using the CMOR program (CMOR = Climate Model Output Rewriter) itself, only the definitions provided and described by it!
CMOR tables reporting the definitions are available at: `<https://github.com/PCMDI/cmip5-cmor-tables>`_
Either, you perform any adjustments with the widely known tools (e.g. `cdo, nco <https://www.unidata.ucar.edu/software/netcdf/software.html>`_) or you make use of the CDS toolbox procedures (**recommended**) as described in Section 2.

3. Preparation of the SPQB reports with the ESMValTool
------------------------------------------------------

**3.1 Installation of the ESMValTool on the local servers**

For the installation of the necessary python modules and ncl to be able to run the ESMValTool, please follow the steps outlined below:

1. GitHub:

*	Open a GitHub account (`<http://www.github.com>`_).
*	Send your GitHub user name to "Axel.Lauer@dlr.de" to request access to the private branch of the ESMValTool with the note that you work for the C3S_511 service. 
*	After Axel adds you to the ESMVal group on GitHub, you should have access to `<https://github.com/ESMValGroup/ESMValTool-private/tree/development>`_
*	Familiarize yourself with GitHub and the ESMValTool workflow. An introduction can be found here: `<http://esmvaltool.readthedocs.io/en/latest/annex_b.html>`_

2. Read the installation instructions that are given in the ESMValTool manual.

*	You should start with installing the ESMValTool on your (Linux) computer or your institute’s computing facilities (e.g. a cluster). A step-by-step installation guide is given in the User Manual: `<http://esmvaltool.readthedocs.io/en/latest/install.html>`_
*	Additionally, please install sphinx (‘conda install sphinx’), and make sure that you have Latex installed on your machine.	
*	If further support is needed for the installation or the recommended test, please contact your IT people.

If there are still technical issues after you followed the outlined steps, please contact the C3S_511 service support ("C3S_511_Support@dlr.de") with a proper problem description including configuration details.


**3.2 Running the namelist for the SPQB**

When you are familiar with the ESMValTool after following Section 3.1, git checkout the branch C3S_511_beta. Follow the description below in addition to the general guidance from the ESMValTool.

Before you run the SPQB namelist, you should check and update the following files:

*	Diagnostic specific cfg-file: This file is called “cfg_C3S_511.py” and is located in the directory “ESMValTool-private/nml/cfg_C3S_511/”. Here you can specify your preferences about 3D variable levels, your ECV specific color scheme, and your preferred output. More detailed instructions on how to do this are given below.

*	namelist: You will have to adjust the namelist to specify the data set that you want to produce the SPQB for. The namelist is called “namelist_C3S_511_SPQB_beta_wpp.xml”, and it is located in the directory “ESMValTool-private/nml/”. There are three parts in the namelist that need adjustments:

1.	Adjust the file path/name to your specific environment cfg-file (line 2 of the namelist) that includes the file paths for your specific working environment. The file is a xml-file, and is probably called something like “config_private.xml”
2.	In the diagnostics part of the namelist (this starts with the keyword <DIAGNOSTICS>), adjust all necessary parts (e.g. diagnostic specific cfg-file, <variable ref_model="??">, <field_type>, <model> …), so that your specific ECV can be read and processed.

After you have adjusted these files, you can run the SPQB namelist as described in the ESMValTool manual. Please be aware that you will have to run the namelist twice to produce the final reports with all additional input! Between the first and the second run you will have to finalize some files (these will be described in the sections below), so that this information can be added to the report during the second run of the namelist.


**3.3 C3S_511 SPQB Configuration (cfg) file options**

*Definition of levels for 3D variables*

If you have to provide reports for a 3D variable with the SPQB namelist, you have the option to specify the levels in a list you want to provide figures for in the reports in a configuration file (cfg-file, specified in the namelist). Your selection should be based on your expert opinion on which levels need to be shown to characterize the specific ECV. Please keep in mind that the number of figures shown in the reports for a 2D variable is multiplied by the number of levels you specify in the cfg-file (e.g. 3 levels selected -> 3 x number of trend plots for a 2D variable), so please select your levels carefully to avoid too many figures in the reports!

The levels that you specify have to be given in the respective unit, and they have to be available in the dataset that you assess. There is no level interpolation available (since this would provide information in the QB that is not available in the dataset)! If the level you specify is not available in the dataset, the ESMValTool will crash while running the SPQB namelist. If you are unsure, which levels are available, you can run the tool once before and you will get information from the first run.

*Definition of data color map*

Please specify a custom color map in the cfg file. This color map will then be used in the graphs for the mean and variability. Possible color maps are available here: `<https://matplotlib.org/examples/color/colormaps_reference.html>`_ 

*Definition of latex output*

*	For debugging purpose, you can put the latex option to True (“show_latex=True”). If you have installed ‘sphinx’ and ‘latex’ correctly, you should get the output from producing the pdf-files of the different reports. 
*	Recommended setting: the option False (“show_latex=False”). This allows you to avoid the production of the pdf-files every time you run the SPQB namelist, as the output is lengthy. 
*	If you have problems with producing latex output (latex is not running smoothly on your machine where you run the ESMValTool) you can only produce the figures and the latex file in sphinx compatible format (“show_latex=None”), and port these files to another latex compatible machine to compile them separately to a pdf file outside of your ESMValTool environment. With this option, the SPQB namelist automatically copies the latex script for the creation of the reports to reporting directories. Their structure is self-explaining.

**3.4 Input/Output structure for the SPQB**

When you set up your general ESMValTool configuration (in your ESMValTool-private directory), you defined your work directory. Within this directory, you will have two relevant and SPQB related subdirectories, one called “c3s_511” and one called “reporting”. The reporting directory contains your pdf output, or, if the latex-option was set to “None”, the respective built and source directories. This is the output you need for the reports. The directory c3s_511 contains all editable files. The first run produces all files according to your data set name in the namelist. Therefore, if you run the ESMValTool with the same namelist a second time, it will read in these files and check if information was added or, with some specific files, corrected. The ESMValTool also reports the set up files in the terminal output. Please check these files for:

*	Adding additional text to the reports.
*	Filling out information needed for SMM, APM, etc.
*	Adjusting information like the original resolution for the gcos requirements checks.

The following subsections explain the needs for the single reports in more detail.

**3.5 How to add customary text for the individual reports**

You will have to run the SPQB namelist for each ECV twice to be able to display customary text.

The first time you run the SPQB namelist, you produce the figures for each report with their respective figure numbers and figure captions. During that first run, an empty text file for each report that you want to produce is created. These will be located in individual report directories (e.g. a directory named ‘overview_input’). You can then add your text with the figure interpretation and comments in the text file that is available in each of these directories. **Please note:** There is no txt-file in the folder ‘smm_input’ but a csv-file instead! You can add the custom text there. (This is differing to the other text files to future upward compatibility.)

If you want to add a specific format to your customized text (more than having it appear as plain text), you will have to add the text in the 'reStructuredText'-format (see `<http://docutils.sourceforge.net/docs/user/rst/quickref.html>`_ for a brief documentation about the rst-format). The following webpages show examples on how to add references and footnotes to an rst-style file:

*	`<http://docutils.sourceforge.net/docs/user/rst/quickref.html#hyperlink-targets>`_
*	`<http://docutils.sourceforge.net/docs/user/rst/quickref.html#footnotes>`_

A useful online tool for reStructuredText can be found here: `<http://rst.ninjs.org/#>`_. It allows you to format your text without having to run the ESMValTool over and over.

After you have added your customized text, you have to run the SPQB namelist a second time to get the text added to the report(s). If you do not edit the text file for a report, no individual text will be added to that final report (which means that the contents of the text file produced with the first run of the ESMValTool is not displayed if it was not altered!). If you do add individual text, it will appear in the final report after the Table of Contents, any lists and before the figures to follow a paper draft design. **Please note:** for each report you will have to add the customary text individually!

**3.6 What to provide for the System Maturity Matrix (SMM) and how**

The C3S_511 SMM is derived from the Core Climax and adapted to the specific service needs. In order to produce a SMM in the report it is necessary to fill the “SMM_CORE_CLIMAX_c3s_Adapted_v5.0.xlsx” xlsx-file provided with the code of the SPQB namelist. This file contains the different specifications for each SMM category and subcategory. Additional help can be found in “SMM_Guide_for_USERS_C3S_511_v1.pdf” and the “CORE_CLIMAX_MANUAL.pdf”.

Guidelines and help for understanding the different categories and subcategories, as well as an export option for the required csv-file (asked for by the ESMValTool in your c3s_511 directory) are provided in the above mentioned xlsx-file, as well as in the rest of this section where the SMM fields are mapped to the fields expected from the EQCO (Evaluation and Quality Control for Observations) Service of the CDS reported in *CAPITAL ITALIC*. **Please note:** Some categories are left blank on purpose as they are currently not relevant for C3S_511 purposes but left in for completeness of the approach or for eventually later automatic filling when coupled to the QATs.

After filling, export the respective sheet then to the file that is requested by the first run of the ESMValTool (c3s_511 directory). Run the ESMValTool a second time for your chosen ECV. Afterwards, the SMM table cells should be colored. The colors represent the different subcategories and are based on the numbers that you have added to the csv-file.

**Please note:** if the SMM csv-file, which is requested by the ESMValTool after the first run, is empty; if not changed, the table cells will not contain colors after the second time you run the ESMValTool!

Guidance on the different categories of the SMM to be integrated with the xlsx-file:

**Software readiness:** The section is left blank on purpose.

**Metadata:** Metadata information has to be tested from the datafile itself. Access to the metadata may be done for instance using basic instructions (e.g. *ncdump –h [filename]*). Please note: this has to be applied to the original data (subset).

* Standard: *Is there any standard used?* 
Check the used metadata convention (original file) and whether the convention is CF_Convention or if there is any tool to translate the used standard to the CF_Convention. (ESMValTool does not run without this Convention.)

* Collection Level: *Is there the possibility to read in metadata?*
Sufficient for use – basic geolocation and sensor/platform identification
Enhanced detailed metadata (see as example the necessary fields in `<https://data.noaa.gov/datasetsearch/>`_)

* File level
The section is currently left blank on purpose. 

**Documentation:** All information on documentation are currently gathered from the documentation available on the CDS itself and / or from the EQCO framework.

*Formal description of scientific methodology: which level of description? -> see PRODUCT GENERATION: DOCUMENTATION & REFERENCES*

*Formal Product User Guide: is it available and updated? -> see QUALITY INDICATORS: DOCUMENTATION & REFERENCES*

*Formal Validation report: is it available and updated and reports uncertainties? -> see PRODUCT VALIDATION: DOCUMENTATION & REFERENCES*

*Formal Description of operations concepts* -> The section is currently left blank on purpose.

**Uncertainty:** All information on uncertainty are currently gathered from the documentation available on the CDS itself and / or from the EQCO framework.

*Standard: level of standard used for uncertainty* ->	see *UNCERTAINTY CHARACTERISATION*: Metrologically Assessed

*Uncertainty Validation* -> see *PRODUCT VALIDATION: DOCUMENTATION & REFERENCES*

*Uncertainty Quantification* -> see *PRODUCT VALIDATION: DOCUMENTATION & REFERENCES*

*Quality Monitoring* ->	see *QUALITY INDICATORS: QUALITY CONTROL*

**Public access, feedback, update:** The section is currently left blank on purpose.

**Usage:** The section is currently left blank on purpose.

In addition to these guidelines, Core Climax heritage material is available through `<https://drive.google.com/open?id=1hm5IHx-Nxl3ouVjwGwwuPT2tVmsKeL1g>`_.


**3.7 What to provide for GCOS requirement and how [optional]**

The GCOS requirements are in principle checked automatically based on the internal table, reporting these and the scan of the data performed by the ESMValTool. Nevertheless, several datasets are larger than the recommended size (see Section 2.3) and it might be necessary to reduce their size via resolution degradation. In this case it is necessary to adjust the real dataset resolution/temporal coverage in a correctional step. The calculated values can be found in the file “[ECV name]_gcos_values_editable.csv” and can be edited therein. 
Please note: only change values if really necessary.

**3.8 What to provide for the ESM evaluation and how**

The ESM evaluation report consists of two parts. The first part is a graphical/tabular display about the suitability of the temporal and spatial resolution of the ECV for ESM evaluation (based on the assessment of the expert users), the second part is a list of references where the respective ECV, in the product or a similar product, has been used previously for ESM evaluation.

For the first part, the expert user will have to provide an estimate about the necessary length of an ECV to be useful for the following applications related to ESM evaluation: mean/climatology, trends, and variability. The estimates have to be added to the csv-file “[ECV name]_esmeval_expert.csv”, that will be created when you run the ESMValTool the first time. The file is located in the folder “work/c3s_511/esmeval_input/”. The information given in this csv-file is then compared to the respective ECV’s temporal and spatial resolution, which will ultimately result in colored table cells for this first part of the ESM evaluation.

For the second part, the expert user will have to provide information about references of the recent literature about the usage of the respective ECV in ESM evaluation studies. The information should be added to the file “esmeval_expert.csv” in the directory “diag_scripts/aux/C3S_511/lib/predef/”. The following pieces of information about the references are needed:

* ECV
*	Product
*	Dataset(s): name of the dataset(s) that are used within the reference for ESM evaluation
*	Title of reference
*	Author(s) of reference: only provide the first authors last name, initial of the first name, and then add ‘et al.’ (e.g. Lauer, A., et al.)
*	Year of publication
*	DOI (reference): doi of the reference
*	Keywords: only provide key words here that describe for what purpose the ECV has been used in the respective ESM evaluation (e.g. climatology, trends, etc.). Please do not provide the keywords of the study here!

There are already plenty of examples in the file. Please follow their example, and add more references if necessary. The information from this csv-file will then be added to the report as bullet point list. If you update this file, please make sure that the content is made available to all the ESMValTool users, either by adding them to the branch and informing WP2 to check it, or sending the library file to WP2. 

**3.9 What to provide for the Sectoral Information System section and how**

Sectoral will come in the form of GCOS-like requirements based on the (SIS) User’s feedback. For the beta version they are not implemented. Demonstration will be provided on demand as “place holder” in the reports following the example in the template.

**3.10 What to provide for the Application Performance Matrix (APM) and how**

APM is not implemented fully in the beta version. Demonstration will be provided on demand as “place holder” in the reports following the example in the template.


4. Deliverable template
-----------------------

The template for the requested reports can be found in the file C3S_D511_template_SPQB_12102018.docx. This template requires the following specific actions:

*	Adjusting any ECV/dataset/producer information on the cover page and the file name.
*	Adding additional information on a possible preprocessing of the data (spatial/temporal aggregation, CMORization).
*	Adding the executive summary of the report with selected figures, text and references.
*	Produce the full reports and add them to the word template

Now, you are potentially done, and if satisfied can upload it to the shared repository for internal evaluation.

**Please note:** all SPQB report specific texts can be added with the procedure described in section 3.5.


5. How to report any problems with the SPQB
-------------------------------------------

The C3S_511 service is implementing a quality assurance system that requires to track the development and the problems of the service components. 

a) It is then needed to fill a report feedback form to report any problem/difficulty/missing functionality or information parts encountered found during the reporting activity. These are collected in a specific file Report_Feedback_ECV_Products.doc where ECV products may be specified to each report.

b) In addition, specific issues concerning the ESMValTool and eventual developments are tracked on GitHub in the following way:

The ESMValTool code is hosted on GitHub (`<www.github.com>`_). In the private branch of the ESMValTool (which you should have all access to; if not, please refer to the Workshop Agenda for the Workshop in May 2018; there is a detailed description on how to be added to the hub!), there is a project ‘C3S_511’ (`<https://github.com/ESMValGroup/ESMValTool-private/projects/1>`_, only accessible if you are logged in and part of the hub). The project page contains four different cards which are called ‘To Do’, ‘In Progress’, ‘Done’ and ‘Questions’. 

If something is developed for a specific ECV that might be of interest for other ECVs as well, it is recommended to add a note to the column ‘In Progress’ to let other service members know what is developed and how. This might prevent duplication of code development and produce synergies.

If there are any questions or comments about the most recent release of the SPQB namelist, it is recommended to add these as a note or an issue to the card ‘Questions’. In doing so, the questions and comments are available for all service members to see, it is possible to trace back who had posted the question/comment, and it can be made sure that all comments/questions are answered and dealt with. Please do not send any questions/comments directly to DLR or LMU, but post them on GitHub to ensure that all comments and questions can be dealt with, and that we can trace our efforts/work!


