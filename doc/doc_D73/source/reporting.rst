Reporting Service for ESMValTool
=====================

Why reporting?
------------

The ESMValTool is supposed to provide automated reporting on the calculated diagnostics. The automated reporting service is supposed to be web-based, so the server preparing the data, can also provide the resulting report without forcing the user to download data or install additional software. The automated reporting service should care for gathering the output from the ESMValTool and present it in a flexible manner.

Reporting in ESMValTool
~~~~~~~~~~~~~~~~~~~~~~~

The current reporting service composes a functional and flexible html-based website that is “printer-friendly” to also provide results in a downloadable pdf-format.

Two different approaches were primarily considered for realization:

* The first approach is, to implement a serial reporting service that is running a **collector (and organizer)** after the ESMValTool produced the output for the reporting. Therefore, the ESMValTool must deliver collectable information on the structure of the outputs and their pattern for the report. This approach allows to prevents redundant and time consuming calculations.

* The other approach is, to implement the collector (and organizer) as part of the ESMValTool runtime. An ESMValTool **run-time environment** is needed for this. This approach basically is preferred, as the reporting service is in charge of producing and managing the output of the ESMValTool, including reading or setting up namelists, etc. The advantage is, that results from multiple namelist can be easily incorporated. Therefore, information on the reporting structure is not mandatory. It is also possible to easily report on diagnostics that are currently not supporting such kind of information.

The current version is a hybrid form of the aforementioned approaches. If the reporting service recognizes an ESMValTool namelist as input, the tool acts as a run-time environment for the tool and collects multiple diagnostic blocks' output into seperately reported parts. If the reporting service receives a specific report namelist, former results are gathered from predefined search directories and are prepared based on specific grouping instructions.

*put images here*

*missing: new hybrid*

.. figure:: reporting_post_workflow.png
   :scale: 100 %
   :alt: Reporting service as a collector

   The reporting service implemented as collector for ESMValTool output based on specific reporting namelists

.. figure:: reporting_envi_workflow.png
   :scale: 100 %
   :alt: Reporting service as an environment

   The reporting service implemented as environment for ESMValTool output distributing original namelists


Requirements
------------

MetaData for Files

V1: time synchronous setup

V2: tags version in nmls, tagged images




1) Specify MetaData
~~~~~~~~~~~~~~~~~~~

xml object

current limitation


2) Specify namelist tags
~~~~~~~~~~~~~~~~~~~~~~~~

Global

Diagnostic


3) Specify report namelist 
~~~~~~~~~~~~~~~~~~~~~~~~~~

Tags, Folders


Examples
--------






