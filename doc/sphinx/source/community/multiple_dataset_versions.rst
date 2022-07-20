Support for several dataset versions
************************************

If you plan to update a CMORizer script to account for a newer version of a dataset,
indicate in the issue or pull request if support for previous versions should be kept.
If the dataset is used in recipes, please also indicate if the recipes should be updated with
the newest dataset version.

Policy to add and retire support
================================

Support for older versions should preferably be kept as long as the data are publicly available.
This ensures reproducibility and eases comparison of results of recipes using this dataset.

If previous dataset versions are no longer available or data issues have been fixed in a newer
dataset version, it is preferable to only keep support for the latest version.

Naming conventions
==================

If the data structure is rather similar between versions, a single CMORizer script (e.g. woa.py)
and config file (e.g. WOA.yml) should be favored to handle multiple verions and avoid code duplication.
Version-dependent data fixes can be applied based on the ``version`` keys defined in the config file.

In some cases, it can be simpler to use different names for different dataset versions (e.g. GCP2018 and
GCP2020). CMORizer scripts and config files should be named accordingly.
