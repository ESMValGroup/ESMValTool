.. _dataset-versions:

Support for multiple versions of a dataset
******************************************

If you plan to update a CMORizer script to support a newer version of an
existing dataset, indicate in the issue or pull request if support for previous
versions should be kept.
If the dataset is used in recipes, please also indicate if the recipes should
be updated with the newest dataset version.

Policy for dropping support for older dataset versions
======================================================

Support for older versions should preferably be kept as long as the data are
publicly available.
This ensures reproducibility and eases comparison of results of recipes using
this dataset.

Even when previous dataset versions are no longer available or data issues have
been fixed in a newer dataset version, it is preferable to keep support for the
previous version in addition to supporting the newer version.
In such cases, it is recommended to ask the recipe maintainers of recipes using
the older version of the dataset to update to the newer version if possible so
that support for the old version can be dropped in the future.

Naming conventions
==================

If the data structure is rather similar between versions, a single CMORizer
script (e.g. woa.py) and config file (e.g. WOA.yml) should be favored to handle
multiple versions and avoid code duplication.
Version-dependent data fixes can be applied based on the ``version`` keys
defined in the config file.

In some cases, it can be simpler to use different names for different dataset
versions (e.g. GCP2018 and GCP2020).
CMORizer scripts and config files should be named accordingly.
