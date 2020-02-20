.. _recipes_seaice_drift:

Seaice drift
============

Overview
--------
This recipe allows to quantify the relationships between Arctic sea-ice drift
speed, concentration and thickness (Docquier et al., 2017). A decrease in
concentration or thickness, as observed in recent decades in the Arctic Ocean
(Kwok, 2018; Stroeve and Notz, 2018), leads to reduced sea-ice strength and
internal stress, and thus larger sea-ice drift speed (Rampal et al., 2011).
This in turn could provide higher export of sea ice out of the Arctic Basin,
resulting in lower sea-ice concentration and further thinning. Olason and
Notz (2014) investigate the relationships between Arctic sea-ice drift speed,
concentration and thickness using satellite and buoy observations.
They show that both seasonal and recent long-term changes in sea ice drift are
primarily correlated to changes in sea ice concentration and thickness.
This recipe allows to quantify these relationships in climate models.

In this recipe, four process-based metrics are computed based on the multi-year
monthly mean sea-ice drift speed, concentration and thickness, averaged over
the Central Arctic.

The first metric is the ratio between the modelled drift-concentration slope
and the observed drift-concentration slope. The second metric is similar to the
first one, except that sea-ice thickness is involved instead of sea-ice
concentration. The third metric is the normalised distance between the model
and observations in the drift-concentration space. The fourth metric is similar
to the third one, except that sea-ice thickness is involved instead of sea-ice
concentration.

Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_seaice_drift.yml


Diagnostics are stored in diag_scripts/seaice_drift/

    * seaice_drift.py: Compute the fourth metrics and generate the associated plots


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
