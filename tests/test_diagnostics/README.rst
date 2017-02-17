Testing in ESMValTool
=====================

This documentation briefly describes how testing of ESMValTool works. This document is supposed to become part of the SPHINX based documentation later and should be moved there and modified accordingly.

Why testing?
------------

TBD

Getting started
---------------

Prerequesite #1: install `easytest`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ESMValTool testing facility is based on the python package `easytest <https://github.com/pygeo/easytest>`_ . If `easytest` is not installed yet on your machine, install it within seconds using the instructions provided `here <http://easytest.readthedocs.org/en/latest/>`_. There you find also the information on the types of tests supported so far. To use the `easytest` functionality for testing entire ESMValTool diagnostics only a few lines of code and reference data is needed. The implementation of a test for a new diagnostic will be described in the following.

Prerequesite #2: nosetests
~~~~~~~~~~~~~~~~~~~~~~~~~~

`nosetests <https://nose.readthedocs.org/en/latest/>`_ should have been already installed. If not, do so as it is needed for the entire testing of ESMValTool.


Prerequesite #3: provide reference data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You should collect all reference data (= expected results) for a diagnostic into a directory. This is easiest achieved by running your diagnostic script after finalizing the development and then put the output (typically ./work/plots ./work/climo) into some directory which we will call the *reference directory* subsequently.


Implement a test for a new diagnostic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here you go ... To implement the test for your diagnostic only a few more steps are needed:

1. create your test file. The filename should include your diagnostic name::

    cd ./tests/test_diagnostics
    cp test_namelist_TEMPLATE.py test_namelist_<your_diagnostic_name>.py

2. edit your test file as follows (only two changes are needed)::

    [...]
    def __init__(self):

        # 1) specify here the full path of the namelist to be tested (relative to ESMValTool root)
        nml = 'nml/test_suites/dlr/namelist_NAMEOFDIAGNOSTIC.xml'  # <<<<<<<<<  put here the name of the diagnostic to execute for testing

        # 2) define here the location of the reference directory
        #    note that it is expeced that the directory has the same name as the namelist
        refdir = '/reference/data/directory'  # <<<<<<<<<<<<<<<  put here your reference data directory
        super(MyDiagnosticTest,self).__init__(nml=nml, refdirectory=refdir, esmval_dir=esmval_dir)


That's it! To run the tests you simply do::

    nosetests test_namelist_<your_diagnostic_name>.py

This will run your diagnostics and after this was sucessfully completed the tests are performed. In the end you should get an::

    O.K!

In any other case, failures (e.g. missing files, files with different content) will be reported. Further interested, then keep reading ...


How does the testing work in principle?
---------------------------------------

The general concept of the testing framework is that it compares results of a diagnostic with reference data generated once with the same inputs for the same diagnostic. The following checks can be performed:

1. check if all output files are available (Filecheck)
2. check that content of output files is the same (MD5 checksum check)
3. check that graphic files look the same (graphic check)  - TBD in development

What is needed?

* a namelist for your diagnostic tailored for your tests
* test data
* a script that implements your test

For more detailed information about the different test implementations refer to
the `easytest` `documentation <http://easytest.readthedocs.org/en/latest/>`_.


Best practice
~~~~~~~~~~~~~

* small test data packages: During the development, the tests will be frequently executed. You should
therefore use a small set of testdata to test the functionality of your
diagnostic
* to obtain reliable test results it is recommended to clean up all output
  directories (work, plots, temp) before running the tests. Othwise you don't
  know if old files were used for the testing.


Known issues
~~~~~~~~~~~~

The following issues are known:

* postscript file content can not be tested: As postscript output has always
  different header information, the MD5 checksum will always differ. The
  similarity of two postscript files can therefore currently not be checked.
* portability of tests across different user machines is currently not ensured.
  When tests are run on a different machine, the test data package needs to be
  available and filenames in the test namelists need to be adapted. In the mid
  term future this should be used by using environment variables.
