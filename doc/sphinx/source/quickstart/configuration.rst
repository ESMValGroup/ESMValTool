.. _config:

*************
Configuration
*************

The ``esmvaltool`` command is provided by the ESMValCore package, the
documentation on configuring ESMValCore can be found
:ref:`here <esmvalcore:config>`.
An overview of all configuration options can be found
:ref:`here <esmvalcore:config_options>`.
In particular, it is recommended to read the section on how to :ref:`specify
configuration options  <esmvalcore:config_for_cli>` and the section on
:ref:`data sources <esmvalcore:config-data-sources>`.

To install the default configuration in the default location, run

.. code:: bash

   esmvaltool config copy defaults/config-user.yml

Note that this needs to be customized using the instructions above, so
the ``esmvaltool`` command can find the data on your system, before it can run
a recipe.

There is a lesson available in the
`ESMValTool tutorial <https://tutorial.esmvaltool.org/>`_
that describes how to personalize the configuration. It can be found
`at this site <https://tutorial.esmvaltool.org/03-configuration/index.html>`_.
