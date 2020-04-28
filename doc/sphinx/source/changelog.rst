Changelog
=========

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
