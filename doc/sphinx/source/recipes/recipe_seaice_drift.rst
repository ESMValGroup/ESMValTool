.. _recipes_seaice_drift:

Seaice drift
============

Overview
--------
[Please, fill this]


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_seaice_drift.yml


Diagnostics are stored in diag_scripts/seaice_drift/

    * seaice_drift.py: [Please, fill this]


User settings in recipe
-----------------------

#. Script diag_shapeselect.py

   *Required settings (scripts)*

    One of the following two combinations is required:

    1. Latitude threshold:

        * latitude_threshold: metric will be computed north of this latitude value

    2. Polygon:

        * polygon: metric will be computed inside the give polygon. Polygon is defined as a list of (lon, lat) tuple

        * polygon_name: name of the region defined by the polygon


Variables
---------

* sispeed, sithick, siconc (daily)
