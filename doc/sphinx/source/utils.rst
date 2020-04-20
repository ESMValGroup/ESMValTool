.. _utils:

Utilities
*********

This section provides information on small tools that are available in the ``esmvaltool/utils`` directory.


nclcodestyle
============

A tool for checking the style of NCL code, based on pycodestyle.
Install ESMValTool in development mode (``pip install -e '.[develop]'``) to make it available.
To use it, run

.. code-block:: bash


    nclcodestyle /path/to/file.ncl


xml2yml
=======

A tool for converting version 1 recipes to version 2 recipes.
See the README.md file in the directory esmvaltool/utils/xml2yml for detailed usage instructions.


testing
=======

Tools for testing recipes.

test recipe settings
--------------------

A tool for generating recipes with various diagnostic settings, to test of those work.
