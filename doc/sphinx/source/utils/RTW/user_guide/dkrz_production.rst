.. _dkrz_production:

Running the |RTW| on DKRZ in "production mode"
==============================================

.. include:: ../common.txt

"Production" is used to mean an intention to run the |RTW|
on a continual, live basis
to ensure recipe issues are captured when they occur.

In practice,
the production option generates a web page
with test results at:
https://esmvaltool.dkrz.de/shared/esmvaltool/rtw/status_report.html

To run the |RTW| in production mode:

* add the following line to your ``~/.bashrc`` file
  to ensure the Cylc and Rose executables can be found::

    export PATH=/work/bd0854/metomi/bin:$PATH

* run the RTW on DKRZ,
  installing the RTW into a ``recipe_test_workflow_production`` directory::

    cd ESMValTool/esmvaltool/utils/recipe_test_workflow
    cylc vip -O dkrz -O production -n recipe_test_workflow_production
