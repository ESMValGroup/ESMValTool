Quick Start
===========

.. include:: ../common.txt

* Checkout the |RTW|::

    git clone git@github.com:MetOffice/climate-assessment-workflow.git

* Configure the |RTW|::

    cd climate-assessment-workflow/climate-assessment-workflow
    rose edit

* Run the |RTW| at the Met Office, where ``<run-name>`` is a unique run name
  relevant to the current configuration::

    cylc install --run-name=<run-name> -O metoffice
    cylc play climate-assessment-workflow/<run-name>

* Browse the logs using `Cylc Review`_, a web service for browsing logs via an
  HTTP interface.
