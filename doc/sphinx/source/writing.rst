.. _writing:

Writing a diagnostic script or a metrics set
********************************************

Basics
======

The development of a new diagnostic (or set of diagnostics) requires the
following steps before getting started:

	* Creating a *FEATURE BRANCH* (see Part :numref:`annex_b`) in the applicable project subdirectory of the Git repository (via Git, see Section :numref:`git_repository`). Developers are encouraged to work actively through the Git repository. Regular "commits" to the repository help to document changes introduced to the ESMValTool and allow for efficient sharing of code with other developers.
	* Creating a standard namelist following the template described in Section :numref:`std_namelist`.

**General coding rules and conventions**

	* Regular updates of the *FEATURE BRANCH* (see Part :numref:`annex_b`) are strongly recommended in order to keep it synchronized with the *DEVELOPMENT BRANCH* (see Part :numref:`annex_b`).
	* Modularizing all diagnostic scripts as much as possible, using the general-purpose code in *lib/* and separating the diagnostic calculations from the plotting routines.
	* Before creating new functions or procedures, it should be considered to use or extend the existing routines within *lib/*. Each header (see Section :numref:`std_diag`) provides an overview of the already implemented functions and procedures.
	* Functions and procedures specific to a given diagnostic shall go in the subdirectory diag_scripts/*aux/<diagnostic>* (see Section :numref:`tab_direc_struc`).
	* Main namelist, diag_scripts, functions and procedures shall be documented within the respective file using the templates (see Section :numref:`std_namelist` and Section :numref:`std_diag`).
	* Each diag_script shall contain a call to the function *write_reference* (see also Section :numref:`ack_log`) in order to generate a respective acknowledgements log file (Section :numref:`ack_log`).

The reintegration of the feature branch into the *DEVELOPMENT BRANCH* (see Part 
:numref:`annex_b`) can only be done by the core development team (see Section :numref:`core_dev_team`) who
shall be contacted as soon as the *FEATURE BRANCH* is ready for integration into
the *DEVELOPMENT BRANCH*. Before contacting the core development team the
following items should be checked:

	* The new *FEATURE BRANCH* runs with different configuration options.
	* If the *lib/* routines have been modified, all the diagnostics using these routines have to be tested (see automated testing, Section :numref:`auto_test`).
	* The new code complies with the coding rules and standards (see Section :numref:`doc_soft`) and follows the ESMValTool directory structure (see :numref:`tab_direc_struc`).
	* All authors, contributors and data are properly acknowledged and referenced in the acknowledgements log file (see Section :numref:`ack_log`).
	* If the new observational data are used, the scripts to "cmorize" these data shall also be made available and placed as *reformat_obs_<name>* into the folder *reformat_scripts/obs/*. Once the *FEATURE BRANCH* has been integrated into the *DEVELOPMENT BRANCH* (see Part :numref:`annex_b`), it shall be deleted from the repository.


.. _std_diag:

Standard template
=================

All (diagnostic) scripts and namelists in the ESMValTool are documented
following the standards defined by templates (see Section :numref:`std_namelist` for the namelist
template). The following describes the standard header for diagnostics
scripts. The parts marked as *[text]* are the ones to be modified by the author.

	* The modification history is in reverse chronological order (i.e., most recent on top) and the last entry always contains the "written" statement (optionally with a statement such as "based on" if derived from existing code).
	* The author of each entry in the modification history is indicated with the author id as given in the author list in the master reference file (*doc/MASTER_authors-refs-acknow.txt*, e.g., A_surn_na = surname, name).
	* All lines should be limited to a maximum of 79 characters (see Section :numref:`rules`). Exceptions can be made to improve the readability of the code.


.. code-block:: none

   ;;#############################################################################
   ;; TITLE OF THE DIAGNOSTIC
   ;; Author: [Name Surname (Affiliation, Country)]
   ;; [PROJECT-NAME]
   ;;#############################################################################
   ;; Description
   ;;    [A short description of the diagnostic]
   ;;    [Additional description of the diagnostic]
   ;;    [Add more bullets if required]
   ;;
   ;; Required diag_script_info attributes (diagnostics specific)
   ;;    [att1]: [short description]
   ;;        [keep the indentation if more lines are needed]
   ;;    [att2]: [short description]
   ;;
   ;; Optional diag_script_info attributes (diagnostic specific)
   ;;    [att1]: [short description]
   ;;    [att2]: [short description]
   ;;
   ;; Required variable_info attributes (variable specific)
   ;;    [att1]: [short description]
   ;;    [att2]: [short description]
   ;;
   ;; Optional variable_info attributes (variable specific)
   ;;    [att1]: [short description]
   ;;    [att2]: [short description]
   ;;
   ;; Caveats
   ;;    [List possible caveats or limitations of this diagnostic]
   ;;    [Features to-be-implemented shall also be mentioned here]
   ;;
   ;; Modification history
   ;;    [YYYYMMDD-A_X4Y4: extended...]
   ;;    [YYYYMMDD-A_X3Y3: bug-fixed...]
   ;;    [YYYYMMDD-A_X2Y2: adapted to...]
   ;;    [YYYYMMDD-A-X1Y1: written.]
   ;;
   ;; #############################################################################
   
   load ...
   load ...
   
   begin
       ...
       ...
   end
	   
	   
.. _lib:

Library functions
=================

The folder *diag_scripts/lib/* contains general purpose routines used by several diagnostic scripts, these library routines are grouped in subfolders by language, i.e.,

.. centered::
   *diag_scripts/lib/ncl*

.. centered::
   *diag_scripts/lib/python*

.. centered::
   *diag_scripts/lib/R*

Library routines are grouped into individual files by topic, some examples for the NCL library routines are:

	* *diag_scripts/lib/ncl/latlon.ncl*: routines to compute grid cell areas, weighted area averages, etc...

	* *diag_scripts/lib/ncl/regridding.ncl*: routines interfacing the ESMF regridding functions in NCL

	* *diag_scripts/lib/ncl/statistics.ncl*: statistical routines not (yet) implemented in the standard distribution of NCL 

	* *diag_scripts/lib/ncl/style.ncl*: centralized control of NCL plot styles, e.g., defines line colors/dashes/thickness for each model name in CMIP5, based on the style files in *diag_scripts/lib/ncl/styles/*.

For further details on the library functions, see the documentation given in
the header of the functions themselves (see Section :numref:`std_diag` for a template).



.. _plot_func:

Plotting functions
==================

The folder *plot_scripts/* contains general purpose routines used for plotting
by the diagnostic scripts. The plotting functions should facilitate the
separation of computing the diagnostic and displaying the result. To this end
they should handle both the case when called directly from the diagnostic
script (with data to visualize as an argument), and the case when the computed
diagnostic is passed along as a netCDF file. These plotting routines are
grouped in subfolders by language,

	* *plot_scripts/ncl*

	* *plot_scripts/python*

	* *plot_scripts/R*

Each subfolder further groups the plotting routines into files by topic, e.g.,
for the NCL library routines:

	* *plot_scripts/ncl/contour_maps.ncl*: interfaces NCL plotting routines for contour map plots, contour polar maps and adding markers to contour maps

	* *plot_scripts/nc/scatterplot.ncl*: interfaces NCL plotting routines for of scatter plots

For further details on the plotting functions, see the inline documentation in the functions themselves.



.. _new_vars:

Adding new variables
====================


Adding new variables requires changes to *reformat_scripts/recognized_vars.dat*
(Section :numref:`rec_vars`) and possibly also to *reformat_scripts/recognized_units.dat*
(Section :numref:`rec_units`). In addition, a new definition file
*variable_defs/<varname>.ncl* is needed (Section :numref:`var_def`; see :numref:`tab_var_def` for a list
of currently available variable definition scripts). If the variable is a
**non-derived** variable (explained in Section :numref:`var_def`) it also needs to be defined
in a file named *reformat_scripts/cmor/CMOR_<variable>.dat* (see Section :numref:`cmor`).

.. note: New variables have to be added to :numref:`tab_var_def` (in doc/sphinx/source/namelists.rst).

.. _rec_vars:

reformat_scripts/recognized_vars.dat
------------------------------------

New variables have to be added to *reformat_scripts/recognized_vars.dat*. Two
lines are added per variable:

	* |  std_name = varname
	  |  standard CMOR variable name

	* |  alt_name = alternative name 1, alternative name 2, ...
	  |  comma separated list of alternative variable names

**Example (surface pressure)**

	* std_name = ps
	* alt_name = aps,PS,psurf

The ESMValTool reformat scripts will look for variable "varname" in the input
files. If not found, the alternative variable names "alternative name 1",
"alternative name 2", etc. are tried before an error message is issued that
the variable could not be found.


.. _rec_units:

reformat_scripts/recognized_units.dat
-------------------------------------

The file *reformat_scripts/recognized_units.dat* contains a list of known
units. If needed, the unit of the newly added variable can be added. There are
two lines per unit:


	* |  std_name = unit
	  |  standard CMOR unit

	* |  alt_name = alternative unit
	  |  comma separated list of possible alternative units and corresponding conversion factor, defined as units[cmor] = units[alternative] * factor

**Example (dobson units)**

	* std_unit = DU
	* alt_unit = g m-2, 4.6707e-5, kg m-2, mol m-2, 2.2414e-3


.. _var_def:

variable_defs/varname.ncl
-------------------------

The file *variable_defs/<varname>.ncl* is a NCL script containing the
declaration of the variable "varname" including its specific attributes. In
case of derived variables, a function "calculate" calculating the derived
variable must be defined in the script *<varname>.ncl* (see :numref:`tab_var_def` for a list
of currently available variable definition scripts).

|

**Remarks**

    #. For derived variables, a statement specifying the (standard, non-derived) variables required to calculate the derived variable is needed. In the example given below, this statement in the beginning of the NCL script looks like

	.. centered::
	   *;  Requires: rsut:T2*s,rsutcs:T2*s*

       In this example, the two standard variables "rsut" and "rsutcs" are needed to calculate the shortwave cloud forcing.

    #. Variable attributes are specified as attributes of the variable "variable_info" (see examples below). In order to activate the variable attributes, "variable_info" must be set to "True". Some examples for frequently used attributes are:

        * variable_info\@derived = False (True)
        * variable_info\@long_name = "..."
        * variable_info\@units = "..."
        * variable_info\@standard_name = "..."
        * variable_info\@short_name = "..."


**Example (precipitation, standard variable)**

.. code-block:: none

   ; Requires: none
   variable_info = True
   variable_info@derived = False

**Example (shortwave cloud forcing, derived variable)**

.. code-block:: none

   ; Requires: rsut:T2*s,rsutcs:T2*s

   [...]

   variable_info = True
   variable_info@derived = True
   variable_info@long_name = "CS Shortwave cloud radiation effect"
   variable_info@units = "W m-2"

   undef("calculate")
   function calculate( index [1] : integer, \
                       variable [1] : string, \
                       field_type [1] : string )
   ;;                  return_val [1] : logical
   ;; Arguments:
   ;;    index - index to current infile defined in the 
   ;;            'interface_data/ncl.interface'-file
   ;;    variable - Current variable as string
   ;;    field_type - string with field type classification
   ;; Return value:
   ;;    data_new  logical

   local tmp, tmp1, tmp2, dum1, dum2, dum, i, verbosity
   begin
       data_new = True
       tmp1 = read_data(index, "rsut", "T2Ms")
       tmp2 = read_data(index, "rsutcs", "T2Ms")
       dum1 = extract_data(index, tmp1, -1, 0, 0)
       dum2 = extract_data(index, tmp2, -1, 0, 0)

       dum = dum1
       dum = dum2 - dum1
       dum@long_name = variable_info@long_name
       dum@units = variable_info@units
       add_data_var(index, data_new, dum, variable)

       return(data_new)
   end



.. _cmor:

reformat_scripts/cmor/CMOR_variable.dat
---------------------------------------

Each standard variable (non-derived) also needs a configuration file indicating the expected units of the variable. The expected units are read from the file *reformat_scripts/cmor/CMOR_variable.dat* which follows the definitions in the official CMOR tables for CMIP5. If this file is missing for a specific variable, it can be downloaded from http://pcmdi.github.io/cmor-site/tables.html. If a CMOR table for the new variable is not available, the user can create a new one based on the existing tables (e.g., following the example in *reformat_scripts/cmor/CMOR_mmrbcfree.dat* based on *reformat_scripts/cmor/CMOR_mmrbc.dat*).

**Example, reformat_scripts/cmor/CMOR_pr.dat**

.. code-block:: none

   SOURCE: CMIP5   
   !============
   variable_entry:    pr  
   !============
   modeling_realm:    atmos
   !----------------------------------
   ! Variable attributes:
   !----------------------------------
   standard_name:     precipitation_flux
   units:             kg m-2 s-1 
   cell_methods:      time: mean
   cell_measures:     area: areacella
   long_name:         Precipitation
   comment:           at surface; includes both liquid and solid phases from  all types
                      of clouds (both large-scale and convective)
   !----------------------------------
   ! Additional variable information:
   !----------------------------------
   dimensions:        longitude latitude time
   out_name:          pr
   type:              real
   valid_min:         0   
   valid_max:         0.001254
   ok_min_mean_abs:   2.156e-05
   ok_max_mean_abs:   3.215e-05
   !----------------------------------


.. _rules:

Coding rules and standards
==========================

The purpose of the code conventions used in ESMValTool is to ensure a high
degree of consistency in the code layout. Consistently structured code
increases readability and understanding of the code making it easier for
developers and users work with a given piece of the code base. It is important
to emphasize two points:

	* Checking the code consistency should be done by software as this allows the check to be done automatically.
	* Code checkers are available at *util/ncl-checker/pep8.py* (NCL) and *util/pep8-checker/pep8.py* (Python).

The code conventions are guidelines and should be treated as such. There are circumstances when it is advisable, for various reasons such as improved readability, to ignore some of the guidelines.

**Code conventions used for Python**

Python code should conform to the PEP-8 style guide [PEP8 2001]. Recommended
tools to check Python code is the official PEP8-checker that is provided with
the ESMValTool distribution (*util/pep8-checker/pep8.py*) and PyFlakes.

To use it on a python file, cd into util/pep8-checker/, and run,

.. code:: bash

   cd util/pep8-checker
   python pep8.py <path-to-python-file>

Python: Pyflakes

Besides the PEP8-checker also the use of the 'pyflake'-tool is recommended (see the pyflakes homepage https://pypi.python.org/pypi/pyflakes for details). For a local install of pyflakes, try virtualenv, e.g., if the virtualenv already is installed, run

.. code:: bash

	source sandbox-pybot/bin/activate 
	pip install --upgrade pyflakes 
	pyflakes <python-file>

**Code conventions for NCL**

NCL code in ESMValTool should follow the PEP-8 style guides. An NCL adapted version of the Python PEP-8 checker is available in the ESMValTool repository (*util/ncl-checker/pep8.py*). Please note that the NCL checker may report some false-positive (e.g., the reading symbol -> is not recognized as such).

To use the NCL version of the PEP8-checker provided with the ESMValTool distribution, run

.. code:: bash

   cd util/ncl-checker
   python pep8.py <path-to-NCL-file> 
 
The NCL-version is adaption of the Pyhton checker and works satisfactorily as
long as one keeps in mind the false positives it finds due to language
differences between Python and NCL. These false positives may be addressed in
the future depending on priorities.

**Code conventions for R**

The code conventions for R should conform to the formatting produced by the R parser tree. This method is further described at "Tidying-R-code" (https://cran.r-project.org/doc/manuals/R-exts.html#Tidying-R-code). Note that this method can only be considered semi-automatic since it does preserve comments (they need to be repatched) and does not produce very nice line breaks.


.. _doc_soft:

Documentation of software
=========================


In order to ensure that all code can be maintained, all diagnostic packages must be well documented. It is the responsibility of the software developers to embed their documentation into the namelist and source code. For details see Sections :numref:`source_doc` and :numref:`documentation`. Documentation systems exist to organize embedded documentation into well structured, linked documents.

	* **R:** documentation should follow CRAN guidance.
	* **Python:** the Sphinx package allows embedded documentation to be assembled into indexed web pages (see Section :numref:`source_doc`)
	* **NCL and namelists:** a Sphinx extension has been developed to extract code documentation for NCL and namelists (see Section :numref:`source_doc`)


.. _ack_log:

The acknowledgements log file
=============================


The acknowledgements log file automatically created by each diagnostic (see also Section :numref:`diag_avail`) is written by the function *write_references* (*interface_scripts/messaging.ncl*, see below), which uses the tags defined in the master reference/acknowledgements file (*doc/MASTER_authors-refs-acknow.txt*) as input. This master file lists all authors and contributors (tags starting with A\_), the diagnostic references (tags with D\_), references for observational data (tags E\_) and projects (tags P\_).

**The function write_references**

The function write_references (defined in *interface_scripts/messaging.ncl*) should be called at the end of each diagnostic script in order to write the acknowledgements log file. The function has the arguments "author(s)", "contributors", "diagnostics", "observations", "projects" which are arrays of strings. All strings ("tags") used must be defined in the master reference file *doc/MASTER_authors-refs-acknow.txt*. The tags are then replaced by the function *write_references* with their definition when writing the acknowledgements log file. All tags in the master reference file are sorted by category of which there are four in total:

.. code-block:: none

	A_xxx = authors, contributors (xxx = author name)
	e.g., A_###

	D_xxx = diagnostics
	e.g., D_righi15gmd = Righi et al., Geosci. Model Dev., 8, 733-768 doi:10.5194/gmd-8-
          733-2015, 2015.
	
	E_xxx = observational data
	e.g., E_era40 = ERA40

	P_xxx = project
	e.g., P_embrace = EU FP7 project EMBRACE

	write_references(diag_script, \
	        "A_###", \
 		(/"D_righi15gmd", "D_gleckler08jgr"/), \
        	(/"E_kalnay96bams", "E_erainterim", "E_airs", "E_ceresebaf", "E_srb"/), \
		(/"P_embrace", "P_esmval"/))


.. _source_doc:

Documentation of source code
============================

The Sphinx documentation generator (http://sphinx-doc.org) is used to organize
and format ESMValTool documentation, including text which has been extracted
from source code. Sphinx can help to create documentation in a variety of
formats, including HTML, LaTeX (and hence printable PDF), manual pages and
plain text.

Sphinx may be obtained from http://sphinx-doc.org/install.html; an overview of
its workings is available at http://sphinx-doc.org/tutorial.html. In
ESMValTool, Sphinx has been used to set up the files in *doc/sphinx*. Running
*make <target>* in that directory will cause the documentation to be built, and
its output placed in the *build/<target>* subdirectory. Here, *<target>* is the
format required  for example, *html, latexpdf, man* or *text* for the four example
formats mentioned above. Running *make* by itself will generate a complete list
of output formats.

Sphinx was originally developed for documenting Python code, and one of its
features is that it is able  using the so-called autodoc extension  to extract
documentation strings from Python source files and use them in the
documentation it generates. This feature apparently does not exist for NCL
source files (such as those which are used in ESMValTool), but it has been
mimicked (or  more-or-less  reverse-engineered) here via the Python script
*doc/sphinx/scripts/process_ncl_docs.py*, which walks through a subset of the
ESMValTool NCL scripts, extracts function names, argument lists and
descriptions (from the comments immediately following the function
definition), and assembles them in a subdirectory of *doc/sphinx/source*. These
output files are in the so-called reStructuredText format (see, e.g.,
http://docutils.sourceforge.net/rst.html), which is the markup language used
by Sphinx; running make in *doc/sphinx* builds the ESMValTool documentation from
them, as noted above.

.. note:: See Section :numref:`std_sphinx` for more details on how to document a new diagnostic.

.. _auto_test:

Automated testing
=================

Any changes to a programming code have the risk of introducing unwanted side effects on some other parts of a code and introduce bugs. Routine and automated testing is therefore essential to maximize the code quality and ensure integrity of all diagnostics implemented within ESMValTool.



Setup and general workflow
--------------------------

Automated testing within the ESMValTool is implemented on two complementary
levels:

	* **unittests** are used to verify that small code units (e.g. functions/subroutines) provide the expected results
	* **integration** testing is used to verify that a diagnostic integrates well into the ESMValTool framework and that a diagnostic provides expected results. This is verified by comparison of the results against a set of reference data generated during the implementation of the diagnostic.


**Installation of the test environment**

All scripts required to run the test environment are provided together with
the ESMValTool code. Two external python packages are required which can be
installed using the python package manager (pip;
https://pypi.python.org/pypi/pip) as follows in a linux environment:

.. code:: bash

   # install nosetests (https://nose.readthedocs.org/en/latest/)``
   pip install nose
   # install easytest
   pip install easytest


**General functionality of testing framework**

Each diagnostic is expected to produce a set of well-defined results. These are files in a variety of formats and types (e.g. graphics, data files, ASCII files ...). While testing results of a diagnostic, a special namelist file is executed by ESMValTool which runs a diagnostic on a limited set of test data only. A small test data set is chosen to minimize executing time for testing while ensuring on the other hand that the diagnostic produces the correct results. The following general tests are implemented at the moment for diagnostics with available test data:

	* **Check for file availability:** a check is performed that all required output data have been successfully generated by the diagnostic. A missing file is always an indicator for a failure of the program.
	* **File checksum:** While the previous test only checks if a file is available, the checksum verifies if the content of a file is similar. Currently the MD5 checksum is used to verify that contents of a file are the same. The MD5 checksum is a good proxy for the similarity of two files and is used regularly to ensure integrity between files when transferring files between different computers.
	* **Graphics check:** For graphic files an additional test is therefore implemented which verifies that two graphical outputs are identical. This is in particular useful to verify that outputs of a diagnostic remain the same after code changes.


**Testing the ESMValTool diagnostics**

Unittests are implemented for each diagnostic independently. Details on
running unittests using **nose** is as simple as going to the ESMValTool root
directory and then execute the following shell command:

.. code:: bash

   # run nosetests
   nosetests

This will search recursively for test files and execute these tests. A
statistic on success and failures is provided at the end of execution. More
details on using nose can be found in the package’s documentation
(https://nose.readthedocs.org/en/latest/).

To run integration tests for each diagnostic, a small script needs to be
written once. An example for a file named esmvaltooltest.py is provided in
Section :numref:`test_imp`. To run all tests for diagnostics implemented in this
file the following command needs to be executed:

.. code:: bash

   # run integration tests
   python esmvaltooltest.py

A summary of success and failures is provided as output.


.. _test_imp:

Example test implementation for a diagnostic
--------------------------------------------

In the following an example is given how to implement a test environment for a new diagnostic with just a few lines of code.
File: esmvaltooltest.py

.. code-block:: python

   """
   sample script for ESMValTool testing
   """

   from esmvaltool import ESMValToolTest

   """
   Define a new class for testing a particular diagnostic
   """

   class PerfMetricCMIP5Test(ESMValToolTest):
       def __init__(self):
           # 1) define here the name of the test namelist
           nml_name = 'namelist_perfmetrics_CMIP5_test.xml'

           # 2) specify here the full path of the namelist
           # (relative to ESMValTool root)
           nml = 'nml/test_suites/dlr/' + nml_name

           # 3) define here the location of the reference data directory
           #    note that it is expected that the directory has the same
           #    name as the namelist
           refdir = esmval_dir + os.sep + os.path.splitext(nml_name)[0] + \
                    '/output/plots/'

           # initialize the parent class
           super(PerfMetricCMIP5Test,self).__init__(nml=nml,
                 refdirectory=refdir, esmval_dir=esmval_dir)

   # --------------------------------------------

   # This is how you run a test
   PT = PerfMetricCMIP5Test()  # create instance of test class
   PT.run_nml()  # run the testing namelist
   PT.run_tests(execute=False, graphics=None,
                checksum_files='all',files='all')  # perform tests

Checklist
=========

The following table can be used as a list of items to be done/checked when writing a new diagnostic.

  :numref:`tab_checklist` Example checklist for implementing new diagnostics and new observational datasets.

.. tabularcolumns:: |p{4.5cm}|p{4.5cm}|p{4.5cm}|

.. _tab_checklist:

+---------------------------+-------------------------------------------------------+-------------------------------------+
| diagnostic / namelist     | model data                                            | observational data                  |
+===========================+=======================================================+=====================================+
| header in diagnostic code | preprocessing (reformatting routines - if applicable) | preprocessing (reformatting script) |
+---------------------------+-------------------------------------------------------+-------------------------------------+
| header in namelist        | list of tools (if applicable)                         | list of tools (if applicable)       |
+---------------------------+-------------------------------------------------------+-------------------------------------+
| documentation of          | references                                            | header in reformatting script       |
| diagnostic (.rst file +   |                                                       | (if applicable)                     |
| example images, see       |                                                       |                                     |
| Section                   |                                                       |                                     |
| :numref:`documentation`)  |                                                       |                                     |
+---------------------------+-------------------------------------------------------+-------------------------------------+
| provenance (tagging)      | test data                                             | references                          |
| implemented in diagnostic |                                                       |                                     |
| code (see Section         |                                                       |                                     |
| :numref:`tagging`)        |                                                       |                                     |
+---------------------------+-------------------------------------------------------+-------------------------------------+
| testing namelist (runs    |                                                       | test data                           |
| diagnostic with reduced   |                                                       |                                     |
| input data)               |                                                       |                                     |
+---------------------------+-------------------------------------------------------+-------------------------------------+
| coding rules and          |                                                       |                                     |
| standards (see Section    |                                                       |                                     |
| :numref:`rules`) have     |                                                       |                                     |
| been followed             |                                                       |                                     |
+---------------------------+-------------------------------------------------------+-------------------------------------+
| *FEATURE BRANCH* updated  |                                                       |                                     |
| with latest               |                                                       |                                     |
| *DEVELOPMENT BRANCH*      |                                                       |                                     |
+---------------------------+-------------------------------------------------------+-------------------------------------+
| list of tools, libraries, |                                                       |                                     |
| etc.                      |                                                       |                                     |
+---------------------------+-------------------------------------------------------+-------------------------------------+
| references                |                                                       |                                     |
+---------------------------+-------------------------------------------------------+-------------------------------------+
| contact person for        |                                                       |                                     |
| scientific support        |                                                       |                                     |
+---------------------------+-------------------------------------------------------+-------------------------------------+

