.. raw:: latex

   \section*{Known issues}

.. raw:: html

   <h1> Known issues </h1>

**NCL 6.3.0**

* WAMonsoon: There is a problem with the routine mreg_part_corr (local opt should be removed and undef("mreg_part_corr") should be positioned just before.
 * A missing "undef" before the declaration of the function "mreg_part_corr" in the NCL library "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" causes the ESMValTool to crash when this library is loaded more than once within the same execution. This unfortunately happens quite often within the tool, since libraries can be loaded multiple times within the various scripts (diag_scripts, lib/ncl, plot_scripts). *At this time, only versions of NCL older than v6.3.0 can be used with the ESMValTool.* This problem will hopefully be fixed with the next NCL release.

**Namelist "SouthernHemisphere"**

* Models with "dash = 0" in their style definition (diag_scripts/lib/python/style.cfg) will not be plotted in the "fraction occurrence histograms of binned cloud cover" plots created by the diagnostic script diag_scripts/SoutherHemisphere_scatter.py.

