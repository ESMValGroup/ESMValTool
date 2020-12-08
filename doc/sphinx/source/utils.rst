.. _utils:

Utilities
*********

This section provides information on small tools that are available in the
`esmvaltool/utils <https://github.com/ESMValGroup/ESMValTool/tree/master/esmvaltool/utils>`_
directory.

.. _draft_release_notes.py:

draft_release_notes.py
======================

This is a script for drafting release notes based on the titles and labels of
the GitHub pull requests that have been merged since the previous release.

To use the tool, install the package pygithub_:

.. code-block:: bash

   pip install pygithub

Create a `GitHub access token`_ (leave all boxes for additional
permissions unchecked) and store it in the file ``~/.github_api_key``.

Edit the script and update the date and time of the previous release. If needed,
change it so it uses the correct repository.

Run the script:

.. code-block:: bash

   python esmvaltool/utils/draft_release_notes.py

Review the resulting output (in ``.rst`` format) and if anything needs changing,
change it on GitHub and re-run the script until the changelog looks acceptable.
In particular, make sure that pull requests have the correct label, so they are
listed in the correct category.
Finally, copy and paste the generated content at the top of the changelog.

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
Install ESMValTool in development mode (``pip install -e '.[develop]'``) to make it available.
To use it, run

.. code-block:: bash

    test_recipe --help

.. _GitHub access token: https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line
.. _pygithub: https://pygithub.readthedocs.io/en/latest/introduction.html
