.. _config-user:

*************
Configuration
*************

The ``esmvaltool`` command is provided by the ESMValCore package, the
documentation on configuring ESMValCore can be found
`here <https://docs.esmvaltool.org/projects/esmvalcore/en/latest/quickstart/configure.html>`_.
In particular, it is recommended to read the section on the
`User configuration file <https://docs.esmvaltool.org/projects/esmvalcore/en/latest/quickstart/configure.html#user-configuration-file>`_
and the section on
`Finding data <https://docs.esmvaltool.org/projects/esmvalcore/en/latest/quickstart/find_data.html>`_.

To install the default configuration file in the default location, run

 .. code:: bash

	  esmvaltool config get_config_user

Note that this file needs to be customized using the instructions above, so
the ``esmvaltool`` command can find the data on your system, before it can run
a recipe.

There is a lesson available in the 
`ESMValTool tutorial <https://esmvalgroup.github.io/ESMValTool_Tutorial/>`_
that describes how to personalize the configuration file. It can be found
`at this site <https://esmvalgroup.github.io/ESMValTool_Tutorial/03-configuration/index.html>`_.
