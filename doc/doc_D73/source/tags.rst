.. _Lotags:

List of tags
============

The following list provides a number for possible tags and suggestions where to put them.
The format is: "**tag (where to put them)** description"

namelist's name (namelist: Global; automatically added):
    For an easy comparison of different result sets from the same namelist (e.g. after changing a diagnostic).

Auto_Diag_[000-999] (namelist: Diagnostic; automatically added):
    For environment based reports, diagnostic blocks within one namelist are automatically numbered. These tags should not be used within collector based reports.

land (namelist: Global, Diagnostic):
	This namelist/diagnostic considers land surface information.

ocean (namelist: Global, Diagnostic):
	This namelist/diagnostic considers ocean information.

atmosphere (namelist: Global, Diagnostic):
	This namelist/diagnostic considers atmosphere information.

temporal/spatial resolution (namelist: Global, Diagnostic; e.g. monthly, daily, 1x1°):
    The final temporal/spatial resolution of the results.

variable name (namelist: Diagnostic, within code for plots/tables)
    An (alternative) variable short name can be defined within the Diagnostic or dynamically added to plots in the diagnostic script.

variable long name (namelist: Diagnostic)
    A common variable long name is added.

model/observation name (namelist: Diagnostic, within code for plots/tables)
    If the diagnostic is called for a specific model only, or dynamically for plots that only show a summary of a specific dataset.

Multi (namelist: Global, Diagnostic; within code for plots/tables)
    Multiple variables *and* models are analyzed, e.g. perfmetrics.

MultiMod (namelist: Global, Diagnostic; within code for plots/tables)
    Multiple models are analyzed.

MultiVar (namelist: Global, Diagnostic; within code for plots/tables)
    Multiple variables are analyzed.

basic (within code for plots/tables)
    This tag is shared by multiple diagnostics in python using the "./diag_scripts/aux/LMU_ESACCI-diagnostics/diagnostic.py" library.

gmd (within code for plots)
    A spatial plot of the temporal mean difference (= global mean difference).

gmt (within code for plots)
    A spatial plot of the temporal mean time series.

TimeS (within code for plots/tables)
    Time series in plot or table.

4plot (within code for plots)
    Short for a set of four plots that show global mean time series plots of two data sets and their absolute and relative difference.

stattab (within code for tables)
    Tag for statstical information within tables.

reg (within code for plots/tables)
    A subset plot or table for specific predefined regions.

TCorr (within code for plots)
    A spatial plot of the temporal correlation of two time series.

trend (within code for plots)
    A spatial plot of the temporal trend of a time series.
   
Note: Tables are not yet implemented in reporting service!





