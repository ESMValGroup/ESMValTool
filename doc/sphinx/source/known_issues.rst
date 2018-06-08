ESMValTool known issues
***********************

**NCL 6.2.1+**

* The function "isfilepresent" has been updated in NCL 6.2.1 in a way that is backwards-incompatible. A new function "fileexists" has been introduced with NCL v6.2.1:

  https://www.ncl.ucar.edu/prev_releases.shtml#BackwardsIncompatibleChanges6.2.1 

  This issue can be addressed by using the ESMVal-function "isfilepresent_esmval" instead of "isfilepresent" and "fileexists".

**NCL 6.3.0**

* WAMonsoon: There is a problem with the routine mreg_part_corr (local opt should be removed and undef("mreg_part_corr") should be positioned just before.
 * A missing "undef" before the declaration of the function "mreg_part_corr" in the NCL library "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" causes the ESMValTool to crash when this library is loaded more than once within the same execution. This unfortunately happens quite often within the tool, since libraries can be loaded multiple times within the various scripts (diag_scripts, lib/ncl, plot_scripts). *At this time, only versions of NCL older than v6.3.0 can be used with the ESMValTool.* This problem will hopefully be fixed with the next NCL release.

**Namelist "SouthernHemisphere"**

* Models with "dash = 0" in their style definition (diag_scripts/lib/python/style.cfg) will not be plotted in the "fraction occurrence histograms of binned cloud cover" plots created by the diagnostic script diag_scripts/SoutherHemisphere_scatter.py.

