Quick Start
===========

.. include:: ../common.txt

* Checkout the |RTW|::

    git clone git@github.com:ESMValGroup/ESMValTool.git

* Configure the |RTW|::

    cd ESMValTool/esmvaltool/utils/recipe_test_workflow/recipe_test_workflow
    rose edit

* Run the |RTW| at the Met Office, where ``<run-name>`` is a unique run name
  relevant to the current configuration::

    cylc install --run-name=<run-name> -O metoffice
    cylc play recipe_test_workflow/<run-name>

* Browse the logs using `Cylc Review`_, a web service for browsing logs via an
  HTTP interface.
