Introduction
============

An automated testing and reporting framework for the ESMValTool has been developed within CRESCENDO. The technical solution developed will be briefly outlined in this report which is based on the technical documentation also provided to the ESMValTool user. As the main effort has been however the development of the actual technical solution itself, it is summarized in the following first which software components have been developed for the implementation of the technical solution.

Software package documentation
------------------------------

The automated testing and reporting framework will be officially released to the climate research community in a forthcoming release of the ESMValTool. To ensure traceability of the different software components developed and contributing to the ESMValTool testing and reporting framework, we document in the following table the various necessariy software packages.

+---------------+------------+--------------------------------------------------------+
| Name          | version    | source                                                 |
+===============+============+========================================================+
| ESMValTool    | development | https://github.com/ESMValGroup/ESMValTool-private.git |
+---------------+------------|--------------------------------------------------------+
| easytest      | 0.1.5      | PyPi, conda                                            |
+---------------+------------+--------------------------------------------------------+
| dummydata     | 0.1.2      | Pypi, conda                                            |
+---------------+------------+--------------------------------------------------------+
| ???reporting  |            |                                                        |
+---------------+------------+--------------------------------------------------------+

The detailed documentation of the testing and reporting is provided in the following sections. It is consistent with the technical documentation which is also made available to the users of the ESMValTool.
