Reporting Service for ESMValTool
=====================

Why reporting?
------------

The ESMValTool is supposed to provide automated reporting on the calculated diagnostics. The automated reporting service is supposed to be web-based, so the server preparing the data, can also provide the resulting report without forcing the user to download data or install additional software. The automated reporting service should care for gathering the output from the ESMValTool and present it in a flexible manner.

Reporting in ESMValTool
~~~~~~~~~~~~~~~~~~~~~~~

The current draft composes a functional and flexible html-based website that is “printer-friendly” to also provide results in a downloadable pdf-format.

Two different approaches were considered for realization:

* The first approach is, to implement a serial reporting service that is running a collector (and organizer) after the ESMValTool produced the output for the reporting. Therefore the ESMValTool must deliver collectable information on the structure of the outputs and their pattern for the report.

* The other approach is, to implement the collector (and organizer) as part of the ESMValTool runtime. An ESMValTool a run-time environment would be needed for this. Currently, this option is preferred, as the reporting service is in charge of producing and managing the output of the ESMValTool, including reading or setting up namelists, etc. The advantage is, that results from multiple namelist can be easily incorporated. Therefore, information on the reporting structure is not mandatory. It is also possible to easily report on diagnostics that are currently not supporting such kind of information.

The delivered version can be seen as a hybrid form of reporting service. If the reporting service recognizes an ESMValTool namelist as input, the tool acts as a processing environment for the tool and collects multiple diagnostic block outputs into seperately reported parts. If the reporting service receives a specific report namelist, former results are gathered from committed search directories and prepared based on specific grouping instructions.

*put images here*

*V1, V2, new hybrid*


Requirements
------------

V1: time synchronous setup

V2: tags version in nmls, tagged images

MetaData for Files


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






