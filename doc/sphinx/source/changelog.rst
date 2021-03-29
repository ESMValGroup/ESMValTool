.. _changelog:

Changelog
=========

.. _changelog-v2-2-0:

v2.2.0
------

Highlights
~~~~~~~~~~

ESMValTool is now using the recently released `Iris 3 <https://scitools-iris.readthedocs.io/en/latest/whatsnew/3.0.html>`__.
We acknowledge that this change may impact your work, as Iris 3 introduces
several changes that are not backward-compatible, but we think that moving forward is the best
decision for the tool in the long term.


This release includes

Bug fixes
~~~~~~~~~

-  Bugfix: time weights in time_operations (`#1956 <https://github.com/ESMValGroup/ESMValTool/pull/1956>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Fix issues with bibtex references (`#1955 <https://github.com/ESMValGroup/ESMValTool/pull/1955>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Fix ImportError for `configure_logging` (`#1976 <https://github.com/ESMValGroup/ESMValTool/pull/1976>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Add required functional parameters for extract time in recipe_er5.yml (`#1978 <https://github.com/ESMValGroup/ESMValTool/pull/1978>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Revert "Fix ImportError for `configure_logging`" (`#1992 <https://github.com/ESMValGroup/ESMValTool/pull/1992>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix import of esmvalcore _logging module in cmorize_obs.py (`#2020 <https://github.com/ESMValGroup/ESMValTool/pull/2020>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix logging import in cmorize_obs again since last merge was nulled by pre-commit hooks (`#2022 <https://github.com/ESMValGroup/ESMValTool/pull/2022>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Refactor the functions in derive_evspsblpot due to new iris (`#2023 <https://github.com/ESMValGroup/ESMValTool/pull/2023>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Avoid importing private ESMValCore functions in CMORizer (`#2027 <https://github.com/ESMValGroup/ESMValTool/pull/2027>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix extract_seasons in validation recipe  (`#2054 <https://github.com/ESMValGroup/ESMValTool/pull/2054>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__

Deprecations
~~~~~~~~~~~~

-  Deprecate classes Variable, Variables and Datasets (`#1944 <https://github.com/ESMValGroup/ESMValTool/pull/1944>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Python 3.9: remove pynio as dependency and replace with rasterio and pin Matplotlib>3.3.1 and pin cartopy>=0.18 (`#1997 <https://github.com/ESMValGroup/ESMValTool/pull/1997>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Removed write_plots and write_netcdf in some python diagnostics (`#2036 <https://github.com/ESMValGroup/ESMValTool/pull/2036>`__) `Manuel Schlund <https://github.com/schlunma>`__

Documentation
~~~~~~~~~~~~~

-  Update instructions on making a release (`#1867 <https://github.com/ESMValGroup/ESMValTool/pull/1867>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update review.rst (`#1917 <https://github.com/ESMValGroup/ESMValTool/pull/1917>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Add guidance on how to review a pull request (`#1872 <https://github.com/ESMValGroup/ESMValTool/pull/1872>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Adding tutorial links to documentation (`#1927 <https://github.com/ESMValGroup/ESMValTool/pull/1927>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Added bibtex file for schlund20jgr (`#1928 <https://github.com/ESMValGroup/ESMValTool/pull/1928>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Documentation contact added the actual email for the mailing list (`#1938 <https://github.com/ESMValGroup/ESMValTool/pull/1938>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Make CircleCI badge specific to master branch (`#1831 <https://github.com/ESMValGroup/ESMValTool/pull/1831>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Documentation on how to move code from a private repository to a public repository (`#1920 <https://github.com/ESMValGroup/ESMValTool/pull/1920>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Refine pull request review guidelines (`#1924 <https://github.com/ESMValGroup/ESMValTool/pull/1924>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Update release schedule (`#1948 <https://github.com/ESMValGroup/ESMValTool/pull/1948>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Improve contact info and move to more prominent location (`#1950 <https://github.com/ESMValGroup/ESMValTool/pull/1950>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add some maintainers to some recipes that are missing them (`#1970 <https://github.com/ESMValGroup/ESMValTool/pull/1970>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update core team info (`#1973 <https://github.com/ESMValGroup/ESMValTool/pull/1973>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Combine installation from source instructions and add common issues (`#1971 <https://github.com/ESMValGroup/ESMValTool/pull/1971>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update iris documentation URL for sphinx (`#2003 <https://github.com/ESMValGroup/ESMValTool/pull/2003>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix iris documentation link(s) with new iris3 location on readthedocs (`#2012 <https://github.com/ESMValGroup/ESMValTool/pull/2012>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Document how to run tests for installation verification  (`#1847 <https://github.com/ESMValGroup/ESMValTool/pull/1847>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  List Remi Kazeroni as a code owner and sole merger of CMORizers (`#2017 <https://github.com/ESMValGroup/ESMValTool/pull/2017>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Install documentation: mention that we build conda package with python>=3.7 (`#2030 <https://github.com/ESMValGroup/ESMValTool/pull/2030>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Recipe and documentation update for ERA5-Land. (`#1906 <https://github.com/ESMValGroup/ESMValTool/pull/1906>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Update changelog and changelog tool for v2.2.0 (`#2043 <https://github.com/ESMValGroup/ESMValTool/pull/2043>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Final update to the changelog for v2.2.0 (`#2056 <https://github.com/ESMValGroup/ESMValTool/pull/2056>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__

Diagnostics
~~~~~~~~~~~

-  Add mapplot diagnostic to ClimWIP (`#1864 <https://github.com/ESMValGroup/ESMValTool/pull/1864>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Add the option to weight variable groups in ClimWIP (`#1856 <https://github.com/ESMValGroup/ESMValTool/pull/1856>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Implementation of ensemble member recognition to the ClimWIP diagnostic (`#1852 <https://github.com/ESMValGroup/ESMValTool/pull/1852>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Restructure ClimWIP (`#1919 <https://github.com/ESMValGroup/ESMValTool/pull/1919>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__
-  Diagnostic for recipe_eyring13jgr.yml Fig. 12 (`#1922 <https://github.com/ESMValGroup/ESMValTool/pull/1922>`__) `Lisa Bock <https://github.com/LisaBock>`__
-  Added changes in shared functions necessary for schlund20esd (`#1967 <https://github.com/ESMValGroup/ESMValTool/pull/1967>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Adding recipe and diagnostics for Gier et al 2020 (`#1914 <https://github.com/ESMValGroup/ESMValTool/pull/1914>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Added recipe, diagnostics and documentation for Schlund et al., ESD (2020) (`#2015 <https://github.com/ESMValGroup/ESMValTool/pull/2015>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Add PRIMAVERA Eady Growth Rate diagnostic (`#1285 <https://github.com/ESMValGroup/ESMValTool/pull/1285>`__) `sloosvel <https://github.com/sloosvel>`__
-  Implement shape parameter calibration for ClimWIP (`#1905 <https://github.com/ESMValGroup/ESMValTool/pull/1905>`__) `Lukas Brunner <https://github.com/lukasbrunner>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Extended ESRL cmorizer (`#1937 <https://github.com/ESMValGroup/ESMValTool/pull/1937>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Cmorizer for GRACE data (`#1694 <https://github.com/ESMValGroup/ESMValTool/pull/1694>`__) `bascrezee <https://github.com/bascrezee>`__
-  Cmorizer for latest ESACCI-SST data (`#1895 <https://github.com/ESMValGroup/ESMValTool/pull/1895>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix longitude in ESRL cmorizer (`#1988 <https://github.com/ESMValGroup/ESMValTool/pull/1988>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Selectively turn off fixing bounds for coordinates during cmorization with utilities.py (`#2014 <https://github.com/ESMValGroup/ESMValTool/pull/2014>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Cmorize hadcrut5 (`#1977 <https://github.com/ESMValGroup/ESMValTool/pull/1977>`__) `mwjury <https://github.com/mwjury>`__
-  Cmorize gpcc masking (`#1995 <https://github.com/ESMValGroup/ESMValTool/pull/1995>`__) `mwjury <https://github.com/mwjury>`__
-  Cmorize_utils_save_1mon_Amon (`#1990 <https://github.com/ESMValGroup/ESMValTool/pull/1990>`__) `mwjury <https://github.com/mwjury>`__
-  Cmorize gpcc fix (`#1982 <https://github.com/ESMValGroup/ESMValTool/pull/1982>`__) `mwjury <https://github.com/mwjury>`__
-  Fix flake8 raised by develop test in cmorize_obs_gpcc.py (`#2038 <https://github.com/ESMValGroup/ESMValTool/pull/2038>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Automatic testing
~~~~~~~~~~~~~~~~~

-  Switched miniconda conda setup hooks for Github Actions workflows (`#1913 <https://github.com/ESMValGroup/ESMValTool/pull/1913>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix style issue (`#1929 <https://github.com/ESMValGroup/ESMValTool/pull/1929>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix mlr test with solution that works for CentOS too (`#1936 <https://github.com/ESMValGroup/ESMValTool/pull/1936>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Temporary deactivation Github Actions on OSX (`#1939 <https://github.com/ESMValGroup/ESMValTool/pull/1939>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix conda installation test on CircleCI (`#1952 <https://github.com/ESMValGroup/ESMValTool/pull/1952>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Github Actions: change time for cron job that installs from conda (`#1969 <https://github.com/ESMValGroup/ESMValTool/pull/1969>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  CI upload relevant artifacts for test job (`#1999 <https://github.com/ESMValGroup/ESMValTool/pull/1999>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Github Actions test that runs with the latest ESMValCore master (`#1989 <https://github.com/ESMValGroup/ESMValTool/pull/1989>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Introduce python 39 in Github Actions tests (`#2029 <https://github.com/ESMValGroup/ESMValTool/pull/2029>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Remove test for conda package installation on Python 3.6 (`#2033 <https://github.com/ESMValGroup/ESMValTool/pull/2033>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update codacy coverage reporter to fix coverage (`#2039 <https://github.com/ESMValGroup/ESMValTool/pull/2039>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Installation
~~~~~~~~~~~~

-  Simplify installation of R development dependencies (`#1930 <https://github.com/ESMValGroup/ESMValTool/pull/1930>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix docker build (`#1934 <https://github.com/ESMValGroup/ESMValTool/pull/1934>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Use new conda environment for installing ESMValTool in Docker containers (`#1993 <https://github.com/ESMValGroup/ESMValTool/pull/1993>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix conda build (`#2026 <https://github.com/ESMValGroup/ESMValTool/pull/2026>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Improvements
~~~~~~~~~~~~

-  Allow multiple references for a cmorizer script (`#1953 <https://github.com/ESMValGroup/ESMValTool/pull/1953>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Add GRACE to the recipe check_obs (`#1963 <https://github.com/ESMValGroup/ESMValTool/pull/1963>`__) `Rémi Kazeroni <https://github.com/remi-kazeroni>`__
-  Align ESMValTool to ESMValCore=2.2.0 (adopt iris3, fix environment for new Core release) (`#1874 <https://github.com/ESMValGroup/ESMValTool/pull/1874>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Make it possible to use write_plots and write_netcdf from recipe instead of config-user.yml (`#2018 <https://github.com/ESMValGroup/ESMValTool/pull/2018>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Revise lisflood and hype recipes (`#2035 <https://github.com/ESMValGroup/ESMValTool/pull/2035>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Set version to 2.2.0 (`#2042 <https://github.com/ESMValGroup/ESMValTool/pull/2042>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__

.. _changelog-v2-1-1:

v2.1.1
------

This release includes

Improvements
~~~~~~~~~~~~

- Fix the conda build on CircleCI (`#1883 <https://github.com/ESMValGroup/ESMValTool/pull/1883>`__) `Bouwe Andela <https://github.com/bouweandela>`__
- Pin matplotlib to <3.3 and add compilers (`#1898 <https://github.com/ESMValGroup/ESMValTool/pull/1898>`__) `Bouwe Andela <https://github.com/bouweandela>`__
- Pin esmvaltool subpackages to the same version and build as the esmvaltool conda package (`#1899 <https://github.com/ESMValGroup/ESMValTool/pull/1899>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Documentation
~~~~~~~~~~~~~

- Release notes v2.1.1 (`#1932 <https://github.com/ESMValGroup/ESMValTool/pull/1932>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

.. _changelog-v2-1-0:

v2.1.0
------

This release includes

Diagnostics
~~~~~~~~~~~

-  Add extra steps to diagnostic to make output of hydrology/recipe_lisflood.yml usable by the LISFLOOD model (`#1737 <https://github.com/ESMValGroup/ESMValTool/pull/1737>`__) `Jaro Camphuijsen <https://github.com/JaroCamphuijsen>`__
-  Recipe to reproduce the 2014 KNMI Climate Scenarios (kcs). (`#1667 <https://github.com/ESMValGroup/ESMValTool/pull/1667>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Implement the climwip weighting scheme in a recipe and diagnostic (`#1648 <https://github.com/ESMValGroup/ESMValTool/pull/1648>`__) `Jaro Camphuijsen <https://github.com/JaroCamphuijsen>`__
-  Remove unreviewed autoassess recipes (`#1840 <https://github.com/ESMValGroup/ESMValTool/pull/1840>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Changes in shared scripts for Schlund et al., JGR: Biogeosciences, 2020 (`#1845 <https://github.com/ESMValGroup/ESMValTool/pull/1845>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Updated derivation test recipe (`#1790 <https://github.com/ESMValGroup/ESMValTool/pull/1790>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Support for multiple model occurrence in perf main (`#1649 <https://github.com/ESMValGroup/ESMValTool/pull/1649>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Add recipe and diagnostics for Schlund et al., JGR: Biogeosciences, 2020 (`#1860 <https://github.com/ESMValGroup/ESMValTool/pull/1860>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Adjust recipe_extract_shape.yml to recent changes in the example diagnostic.py (`#1880 <https://github.com/ESMValGroup/ESMValTool/pull/1880>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Documentation
~~~~~~~~~~~~~

-  Add pip installation instructions (`#1783 <https://github.com/ESMValGroup/ESMValTool/pull/1783>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add installation instruction for R and Julia dependencies tot pip install (`#1787 <https://github.com/ESMValGroup/ESMValTool/pull/1787>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Avoid autodocsumm 0.2.0 and update documentation build dependencies (`#1794 <https://github.com/ESMValGroup/ESMValTool/pull/1794>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add more information on working on cluster attached to ESGF node (`#1821 <https://github.com/ESMValGroup/ESMValTool/pull/1821>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add release strategy to community documentation (`#1809 <https://github.com/ESMValGroup/ESMValTool/pull/1809>`__) `Klaus Zimmermann <https://github.com/zklaus>`__
-  Update esmvaltool run command everywhere in documentation (`#1820 <https://github.com/ESMValGroup/ESMValTool/pull/1820>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add more info on documenting a recipe (`#1795 <https://github.com/ESMValGroup/ESMValTool/pull/1795>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Improve the Python example diagnostic and documentation (`#1827 <https://github.com/ESMValGroup/ESMValTool/pull/1827>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Improve description of how to use draft_release_notes.py (`#1848 <https://github.com/ESMValGroup/ESMValTool/pull/1848>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update changelog for release 2.1 (`#1886 <https://github.com/ESMValGroup/ESMValTool/pull/1886>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Improvements
~~~~~~~~~~~~

-  Fix R installation in WSL (`#1789 <https://github.com/ESMValGroup/ESMValTool/pull/1789>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Add pre-commit for linting/formatting (`#1796 <https://github.com/ESMValGroup/ESMValTool/pull/1796>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Speed up tests on CircleCI and use pytest to run them (`#1804 <https://github.com/ESMValGroup/ESMValTool/pull/1804>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Move pre-commit excludes to top-level and correct order of lintr and styler (`#1805 <https://github.com/ESMValGroup/ESMValTool/pull/1805>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Remove isort setup to fix formatting conflict with yapf (`#1815 <https://github.com/ESMValGroup/ESMValTool/pull/1815>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  GitHub Actions (`#1806 <https://github.com/ESMValGroup/ESMValTool/pull/1806>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix yapf-isort import formatting conflict (`#1822 <https://github.com/ESMValGroup/ESMValTool/pull/1822>`__) `Stef Smeets <https://github.com/stefsmeets>`__
-  Replace vmprof with vprof as the default profiler (`#1829 <https://github.com/ESMValGroup/ESMValTool/pull/1829>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update ESMValCore v2.1.0 requirement (`#1839 <https://github.com/ESMValGroup/ESMValTool/pull/1839>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Pin iris to version 2 (`#1881 <https://github.com/ESMValGroup/ESMValTool/pull/1881>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Pin eccodes to not use eccodes=2.19.0 for cdo to work fine (`#1869 <https://github.com/ESMValGroup/ESMValTool/pull/1869>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Increase version to 2.1.0 and add release notes (`#1868 <https://github.com/ESMValGroup/ESMValTool/pull/1868>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Github Actions Build Packages and Deploy tests (conda and PyPi) (`#1858 <https://github.com/ESMValGroup/ESMValTool/pull/1858>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Added CMORizer for Scripps-CO2-KUM (`#1857 <https://github.com/ESMValGroup/ESMValTool/pull/1857>`__) `Manuel Schlund <https://github.com/schlunma>`__

.. _changelog-v2-0-0:

v2.0.0
------

This release includes

Bug fixes
~~~~~~~~~

-  Fix pep8-naming errors and fix zmnam diagnostic (`#1702 <https://github.com/ESMValGroup/ESMValTool/pull/1702>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix keyword argument in cmorize_obs (`#1721 <https://github.com/ESMValGroup/ESMValTool/pull/1721>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Fixed JMA-TRANSCOM CMORizer (`#1735 <https://github.com/ESMValGroup/ESMValTool/pull/1735>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  Fix bug in extract_doi_value (`#1734 <https://github.com/ESMValGroup/ESMValTool/pull/1734>`__) `bascrezee <https://github.com/bascrezee>`__
-  Fix small errors in the arctic_ocean diagnostic (`#1722 <https://github.com/ESMValGroup/ESMValTool/pull/1722>`__) `Nikolay Koldunov <https://github.com/koldunovn>`__
-  Flatten ancestor lists for diag_spei.R and diag_spi.R. (`#1745 <https://github.com/ESMValGroup/ESMValTool/pull/1745>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Fix for recipe_ocean_ice_extent.yml (`#1744 <https://github.com/ESMValGroup/ESMValTool/pull/1744>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Fix recipe_combined_indices.yml provenance (`#1746 <https://github.com/ESMValGroup/ESMValTool/pull/1746>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Fix provenance in recipe_multimodel_products (`#1747 <https://github.com/ESMValGroup/ESMValTool/pull/1747>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Exclude FGOALS-g2 due to ESMValCore issue #728 (`#1749 <https://github.com/ESMValGroup/ESMValTool/pull/1749>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Fix recipe_modes_of_variability (`#1753 <https://github.com/ESMValGroup/ESMValTool/pull/1753>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Flatten lists for ancestors for hyint to prevent nested lists. (`#1752 <https://github.com/ESMValGroup/ESMValTool/pull/1752>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Fix bug in cmorize_obs_eppley_vgpm_modis.py (#1729) (`#1759 <https://github.com/ESMValGroup/ESMValTool/pull/1759>`__) `Tomas Lovato <https://github.com/tomaslovato>`__
-  Correct mip for clltkisccp in example derive preprocessor recipe (`#1768 <https://github.com/ESMValGroup/ESMValTool/pull/1768>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update date conversion in recipe_hype.yml (`#1769 <https://github.com/ESMValGroup/ESMValTool/pull/1769>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix recipe_correlation.yml (`#1767 <https://github.com/ESMValGroup/ESMValTool/pull/1767>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add attribute positive: down to plev coordinate in ERA-Interim CMORizer (`#1771 <https://github.com/ESMValGroup/ESMValTool/pull/1771>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix sispeed in recipe_preprocessor_derive_test (`#1772 <https://github.com/ESMValGroup/ESMValTool/pull/1772>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Fix extreme events and extreme index ancestors (`#1774 <https://github.com/ESMValGroup/ESMValTool/pull/1774>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Correct date in output filenames of ERA5 CMORizer recipe (`#1773 <https://github.com/ESMValGroup/ESMValTool/pull/1773>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Exclude WOA from multi-model stats in recipe_ocean_bgc (`#1778 <https://github.com/ESMValGroup/ESMValTool/pull/1778>`__) `Mattia Righi <https://github.com/mattiarighi>`__

Diagnostics
~~~~~~~~~~~

-  Enhancement of the hyint recipe to include etccdi indices (`#1133 <https://github.com/ESMValGroup/ESMValTool/pull/1133>`__) `Enrico Arnone <https://github.com/earnone>`__
-  Add lazy regridding for wflow diagnostic (`#1630 <https://github.com/ESMValGroup/ESMValTool/pull/1630>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Miles default domains to include lat=0 (`#1626 <https://github.com/ESMValGroup/ESMValTool/pull/1626>`__) `Jost von Hardenberg <https://github.com/jhardenberg>`__
-  Miles: selection of reference dataset based on experiment (`#1632 <https://github.com/ESMValGroup/ESMValTool/pull/1632>`__) `Jost von Hardenberg <https://github.com/jhardenberg>`__
-  New recipe/diagnostic:  recipe_li17natcc.yml for Axels GMD Paper (`#1567 <https://github.com/ESMValGroup/ESMValTool/pull/1567>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  New recipe/diagnostics: recipe_deangelis_for_gmdpart4.yml for Axels GMD Paper (`#1576 <https://github.com/ESMValGroup/ESMValTool/pull/1576>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  EWaterCycle: Add recipe to prepare input for LISFLOOD (`#1298 <https://github.com/ESMValGroup/ESMValTool/pull/1298>`__) `Stefan Verhoeven <https://github.com/sverhoeven>`__
-  Use area weighted regridding in wflow diagnostic (`#1643 <https://github.com/ESMValGroup/ESMValTool/pull/1643>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Workaround for permetrics recipe until Iris3 (`#1674 <https://github.com/ESMValGroup/ESMValTool/pull/1674>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  C3S_511_MPQB_bas-features (`#1465 <https://github.com/ESMValGroup/ESMValTool/pull/1465>`__) `bascrezee <https://github.com/bascrezee>`__
-  Additional Land perfmetrics (`#1641 <https://github.com/ESMValGroup/ESMValTool/pull/1641>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Necessary diagnostic from eyring06jgr for the release of version2 (`#1686 <https://github.com/ESMValGroup/ESMValTool/pull/1686>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Drought characteristics based on Martin2018 and SPI for gmd paper (`#1689 <https://github.com/ESMValGroup/ESMValTool/pull/1689>`__) `katjaweigel <https://github.com/katjaweigel>`__
-  Additional features and bugfixes for recipe anav13clim (`#1723 <https://github.com/ESMValGroup/ESMValTool/pull/1723>`__) `Bettina Gier <https://github.com/bettina-gier>`__
-  Gmd laueretal2020 revisions (`#1725 <https://github.com/ESMValGroup/ESMValTool/pull/1725>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Wenzel16nature (`#1692 <https://github.com/ESMValGroup/ESMValTool/pull/1692>`__) `zechlau <https://github.com/zechlau>`__
-  Add mask albedolandcover (`#1673 <https://github.com/ESMValGroup/ESMValTool/pull/1673>`__) `bascrezee <https://github.com/bascrezee>`__
-  IPCC AR5 fig. 9.3 (seasonality) (`#1726 <https://github.com/ESMValGroup/ESMValTool/pull/1726>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Added additional emergent constraints on ECS (`#1585 <https://github.com/ESMValGroup/ESMValTool/pull/1585>`__) `Manuel Schlund <https://github.com/schlunma>`__
-  A diagnostic to evaluate the turnover times of land ecosystem carbon (`#1395 <https://github.com/ESMValGroup/ESMValTool/pull/1395>`__) `koir-su <https://github.com/koir-su>`__
-  Removed multi_model_statistics step in recipe_oceans_example.yml as a workaround (`#1779 <https://github.com/ESMValGroup/ESMValTool/pull/1779>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__

Documentation
~~~~~~~~~~~~~

-  Extend getting started instructions to obtain config-user.yml (`#1642 <https://github.com/ESMValGroup/ESMValTool/pull/1642>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Extend information about native6 support on RTD (`#1652 <https://github.com/ESMValGroup/ESMValTool/pull/1652>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  Update citation of ESMValTool paper in the doc (`#1664 <https://github.com/ESMValGroup/ESMValTool/pull/1664>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Updated references to documentation (now docs.esmvaltool.org) (`#1679 <https://github.com/ESMValGroup/ESMValTool/pull/1679>`__) `Axel Lauer <https://github.com/axel-lauer>`__
-  Replace dead link with ESGF link. (`#1681 <https://github.com/ESMValGroup/ESMValTool/pull/1681>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Add all European grants to Zenodo (`#1682 <https://github.com/ESMValGroup/ESMValTool/pull/1682>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update Sphinx to v3 or later (`#1685 <https://github.com/ESMValGroup/ESMValTool/pull/1685>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Small fix to number of models in ensclus documentation (`#1691 <https://github.com/ESMValGroup/ESMValTool/pull/1691>`__) `Jost von Hardenberg <https://github.com/jhardenberg>`__
-  Move draft_release_notes.py from ESMValCore to here and update (`#1701 <https://github.com/ESMValGroup/ESMValTool/pull/1701>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Improve the installation instructions (`#1634 <https://github.com/ESMValGroup/ESMValTool/pull/1634>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Improve description of how to implement provenance in diagnostic (`#1750 <https://github.com/ESMValGroup/ESMValTool/pull/1750>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Update command line interface documentation and add links to ESMValCore configuration documentation (`#1776 <https://github.com/ESMValGroup/ESMValTool/pull/1776>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Documentation on how to find shapefiles for hydrology recipes (`#1777 <https://github.com/ESMValGroup/ESMValTool/pull/1777>`__) `Jaro Camphuijsen <https://github.com/JaroCamphuijsen>`__

Improvements
~~~~~~~~~~~~

-  Pin flake8<3.8.0 (`#1635 <https://github.com/ESMValGroup/ESMValTool/pull/1635>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Update conda package path in more places (`#1636 <https://github.com/ESMValGroup/ESMValTool/pull/1636>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Remove curly brackets around issue number in pull request template (`#1637 <https://github.com/ESMValGroup/ESMValTool/pull/1637>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Fix style issue in test (`#1639 <https://github.com/ESMValGroup/ESMValTool/pull/1639>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update Codacy badges (`#1662 <https://github.com/ESMValGroup/ESMValTool/pull/1662>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Support extra installation methods in R (`#1360 <https://github.com/ESMValGroup/ESMValTool/pull/1360>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Add ncdf4.helpers package as a dependency again (`#1678 <https://github.com/ESMValGroup/ESMValTool/pull/1678>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Speed up conda installation (`#1677 <https://github.com/ESMValGroup/ESMValTool/pull/1677>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update CMORizers and recipes for ESMValCore v2.0.0 (`#1699 <https://github.com/ESMValGroup/ESMValTool/pull/1699>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Update setup.py for PyPI package (`#1700 <https://github.com/ESMValGroup/ESMValTool/pull/1700>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Cleanup recipe headers before the release (`#1740 <https://github.com/ESMValGroup/ESMValTool/pull/1740>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-    Add colortables as esmvaltool subcommand (`#1666 <https://github.com/ESMValGroup/ESMValTool/pull/1666>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Increase version to v2.0.0 (`#1756 <https://github.com/ESMValGroup/ESMValTool/pull/1756>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update job script (`#1757 <https://github.com/ESMValGroup/ESMValTool/pull/1757>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Read authors and description from .zenodo.json (`#1758 <https://github.com/ESMValGroup/ESMValTool/pull/1758>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update docker recipe to install from source (`#1651 <https://github.com/ESMValGroup/ESMValTool/pull/1651>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Cmorize aphro ma (`#1555 <https://github.com/ESMValGroup/ESMValTool/pull/1555>`__) `mwjury <https://github.com/mwjury>`__
-  Respectable testing for cmorizers/obs/utilities.py and cmorizers/obs/cmorize_obs.py (`#1517 <https://github.com/ESMValGroup/ESMValTool/pull/1517>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Fix start year in recipe_check_obs (`#1638 <https://github.com/ESMValGroup/ESMValTool/pull/1638>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Cmorizer for the PERSIANN-CDR precipitation data (`#1633 <https://github.com/ESMValGroup/ESMValTool/pull/1633>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Cmorize eobs (`#1554 <https://github.com/ESMValGroup/ESMValTool/pull/1554>`__) `mwjury <https://github.com/mwjury>`__
-  Update download cds satellite lai fapar (`#1654 <https://github.com/ESMValGroup/ESMValTool/pull/1654>`__) `bascrezee <https://github.com/bascrezee>`__
-  Added monthly mean vars (ta, va, zg) to era5 cmorizer via recipe (`#1644 <https://github.com/ESMValGroup/ESMValTool/pull/1644>`__) `Evgenia Galytska <https://github.com/egalytska>`__
-  Make format time check more flexible (`#1661 <https://github.com/ESMValGroup/ESMValTool/pull/1661>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Exclude od550lt1aer from recipe_check_obs.yml (`#1720 <https://github.com/ESMValGroup/ESMValTool/pull/1720>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  PERSIANN-CDR cmorizer update: adding the capability to save monthly mean files (`#1728 <https://github.com/ESMValGroup/ESMValTool/pull/1728>`__) `Birgit Hassler <https://github.com/hb326>`__
-  Add standard_name attribute to lon and lat in cmorize_obs_esacci_oc.py (`#1760 <https://github.com/ESMValGroup/ESMValTool/pull/1760>`__) `Tomas Lovato <https://github.com/tomaslovato>`__
-  Allow for incomplete months on daily frequency in cmorizer ncl utilities (`#1754 <https://github.com/ESMValGroup/ESMValTool/pull/1754>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Fix AURA-TES cmorizer (`#1766 <https://github.com/ESMValGroup/ESMValTool/pull/1766>`__) `Mattia Righi <https://github.com/mattiarighi>`__

.. _changelog-v2-0-0b4:

v2.0.0b4
--------

This release includes

Bug fixes
~~~~~~~~~

-  Fix HALOE plev coordinate (`#1590 <https://github.com/ESMValGroup/ESMValTool/pull/1590>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Fix tro3 units in HALOE (`#1591 <https://github.com/ESMValGroup/ESMValTool/pull/1591>`__) `Mattia Righi <https://github.com/mattiarighi>`__

Diagnostics
~~~~~~~~~~~

-  Applicate sea ice negative feedback (`#1299 <https://github.com/ESMValGroup/ESMValTool/pull/1299>`__) `Javier Vegas-Regidor <https://github.com/jvegasbsc>`__
-  Add Russell18jgr ocean diagnostics (`#1592 <https://github.com/ESMValGroup/ESMValTool/pull/1592>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Refactor marrmot recipe and diagnostic to use ERA5 daily data made by new cmorizer (`#1600 <https://github.com/ESMValGroup/ESMValTool/pull/1600>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  In recipe_wflow, use daily ERA5 data from the new cmorizer. (`#1599 <https://github.com/ESMValGroup/ESMValTool/pull/1599>`__) `Peter Kalverla <https://github.com/Peter9192>`__
-  In wflow diagnostic, calculate PET after(!) interpolation and lapse rate correction (`#1618 <https://github.com/ESMValGroup/ESMValTool/pull/1618>`__) `Jerom Aerts <https://github.com/jeromaerts>`__
-  Fixed wenz14jgr (`#1562 <https://github.com/ESMValGroup/ESMValTool/pull/1562>`__) `zechlau <https://github.com/zechlau>`__
-  Update portrait_plot.ncl (`#1625 <https://github.com/ESMValGroup/ESMValTool/pull/1625>`__) `Bettina Gier <https://github.com/bettina-gier>`__

Documentation
~~~~~~~~~~~~~

-  Restructure documentation (`#1587 <https://github.com/ESMValGroup/ESMValTool/pull/1587>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add more links to documentation (`#1595 <https://github.com/ESMValGroup/ESMValTool/pull/1595>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Update links in readme (`#1598 <https://github.com/ESMValGroup/ESMValTool/pull/1598>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Minor improvements to installation documentation (`#1608 <https://github.com/ESMValGroup/ESMValTool/pull/1608>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Add info for new mailing list to documentation. (`#1607 <https://github.com/ESMValGroup/ESMValTool/pull/1607>`__) `Björn Brötz <https://github.com/bjoernbroetz>`__
-  Update making a release documentation (`#1627 <https://github.com/ESMValGroup/ESMValTool/pull/1627>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Improvements
~~~~~~~~~~~~

-  Avoid broken pytest-html plugin (`#1583 <https://github.com/ESMValGroup/ESMValTool/pull/1583>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Remove reference section in config-references.yml (`#1545 <https://github.com/ESMValGroup/ESMValTool/pull/1545>`__) `SarahAlidoost <https://github.com/SarahAlidoost>`__
-  Various improvements to development infrastructure (`#1570 <https://github.com/ESMValGroup/ESMValTool/pull/1570>`__) `Bouwe Andela <https://github.com/bouweandela>`__
-  Install scikit-learn from conda, remove libunwind as a direct dependency (`#1611 <https://github.com/ESMValGroup/ESMValTool/pull/1611>`__) `Valeriu Predoi <https://github.com/valeriupredoi>`__
-  Create conda subpackages and enable tests (`#1624 <https://github.com/ESMValGroup/ESMValTool/pull/1624>`__) `Bouwe Andela <https://github.com/bouweandela>`__

Observational and re-analysis dataset support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Cmorizer for HALOE (`#1581 <https://github.com/ESMValGroup/ESMValTool/pull/1581>`__) `Mattia Righi <https://github.com/mattiarighi>`__
-  Add CMORizer for CT2019 (`#1604 <https://github.com/ESMValGroup/ESMValTool/pull/1604>`__) `Manuel Schlund <https://github.com/schlunma>`__

For older releases, see the release notes on https://github.com/ESMValGroup/ESMValTool/releases.
