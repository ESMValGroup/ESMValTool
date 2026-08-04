.. _quick_start_guide:

Quick Start Guide
=================

.. include:: ../common.txt

* Clone the ESMValTool repository (which contains the |RTW|)::

    git clone git@github.com:ESMValGroup/ESMValTool.git

* Configure the |RTW|::

    cd ESMValTool/esmvaltool/utils/recipe_test_workflow
    rose edit

* Run the |RTW|:

  * at the Met Office::

      export CYLC_VERSION=8
      cylc vip -O metoffice

  * on JASMIN:

    * add the following line to your ``~/.bashrc`` file to ensure the Cylc and
      Rose executables can be found::

        export PATH=/apps/jasmin/metomi/bin:$PATH

    * SSH to the Rose and Cylc server::

        ssh -X cylc

    * run the RTW on JASMIN::

        cd ESMValTool/esmvaltool/utils/recipe_test_workflow
        export CYLC_VERSION=8
        cylc vip -O jasmin

  * on DKRZ:

    * add the following line to your ``~/.bashrc`` file to ensure the Cylc and
      Rose executables can be found::

        export PATH=/work/bd0854/metomi/bin:$PATH

    * run the RTW on DKRZ::

        cd ESMValTool/esmvaltool/utils/recipe_test_workflow
        cylc vip -O dkrz

* Optionally browse the logs using `Cylc Review`_,
  a web service for browsing logs via an HTTP interface.
