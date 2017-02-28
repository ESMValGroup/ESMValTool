Testing in ESMValTool
=====================

This documentation briefly describes how testing of ESMValTool works. 

Why testing?
------------

In a nutshell: To make sure that the software does what it is supposed to do.

The overall objective of automated testing in the ESMValTool is to ensure that

* the software was properly installed on the system
* ensure that the different software components (backend, diagnostics,
  plotting) results in correct results
  
The following figure shows how testing integrates in the overall software development workflow. After a new feature has been implemented, tests should ensure that the altered code of the feature branch does not cause conflicts with the original version of the code. Once tests have passed sucessfully, the feature branch can be merged into the development or master branch.
  

.. figure:: testing_workflow.png
   :scale: 50 %
   :alt: Testing in the continous integration (CI) workflow

   The role of testing in the continous integration (CI) lifecyle of a software project


Testing levels in ESMValTool
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Testing is done on different levels within the ESMValTool

 * `unittests <https://en.wikipedia.org/wiki/Unit_testing>`_ are used to verify that small pieces of the ESMValTool software work correctly. This can comprise the test of specific functions, modules or classes
 * Testing of entire diagnostics are done to verify that a diagnostic produces the right output
 
 When developing a new diagnostic for the ESMValTool or improve some other components of the ESMValTool framework, you should always consider to implement appropriate tests. You might implement tests after having implemented the code functionality, it might be however also worth to consider a more agile `test driven development cycle <https://en.wikipedia.org/wiki/Test-driven_development>`_.
 
 
Unittesting in the ESMValTool
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Test for individual pieces of code are currantly implemented largely based on the python core libraray `unittest <https://docs.python.org/2/library/unittest.html>`_. Please look there or in the `tests`directory for examples how to implement tests for python code.


How does the testing of diagnostics in ESMValTool work in principle?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The general concept of the testing framework is that it compares results of a diagnostic with reference data generated once with the same inputs for the same diagnostic. The following checks can be currently performed:

1. check if all output files are available (Filecheck)
2. check that content of output files is the same (MD5 checksum check)
3. check that output file sizes are greater than zero bytes
4. check that acknowledgments are provided for each diagnostic with the processing

The overall testing approach is currently based on the philosophy that

a) testing should be fast to execture
b) should not require the user to download larger data volumes

For that reason the testing is largely based on synthetic data as an input which allows to control the input and thus also expected output of a specific diagnostic.

**What is needed?**

* a namelist for your diagnostic tailored for your tests
* a script that implements your test


Getting started
---------------

The steps required to sucessfully implement testing for the ESMValTool is described in the folliwing.

To be able to use the testing a number of prerequesites need to be fulfilled.
These are typically already installed when the ESMValTool is installed using
`conda <https://conda.io/>`_.

Prerequesite #1: install `easytest`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ESMValTool testing facility is based on the python package `easytest <https://github.com/pygeo/easytest>`_ . If `easytest` is not installed yet on your machine yet, install it within seconds using the instructions provided `here <http://easytest.readthedocs.org/en/latest/>`_. There you find also the information on the types of tests supported so far. To use the `easytest` functionality for testing entire ESMValTool diagnostics only a few lines of code and reference data is needed. The implementation of a test for a new diagnostic will be described in the following.

Prerequesite #2: nosetests
~~~~~~~~~~~~~~~~~~~~~~~~~~

`nosetests <https://nose.readthedocs.org/en/latest/>`_ should have been already installed. If not, do so as it is needed for a convenient testing experience. Test could also be exectuted without `nosetests``, but this makes life so much easier. 

Prerequesite #3: install `dummydata`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `dummydata <https://github.com/pygeo/dummydata>`_ package is required for the generation of synthetic datasets to be used within the testing. When the ESMValTool is installed via `conda` this will have been installed already automatically. Otherwise install the package like described in its documentation.


How to implement a test for a new diagnostic? THIS SECTION NEEDS COMPLETE REVISION
---------------------------------------------

To implement a test for a new diagnostic only a few steps are required.

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






