.. _config_files:

Configuration files
*******************

nml/cfg_diag/cfg_diag*.typ
==========================

Diagnostic-specific settings are passed via configuration files. These are
collected in the nml directory under subdirectories named like the
corresponding diagnostic (e.g., *cfg_aerosol/, cfg_perfmetrics/*). The suffix
".typ" specifies the language the routine is written in.

There might be more than one configuration script per diagnostic set. All
cfg_* files for a diagnostic set need to be in the same folder specified by
the <diag_script_cfg_dir> entry of the namelist (:emphasis:`nml/namelist_*.xml`) (see
Section :numref:`diag_tag` and :numref:`tab_diag_tags`).

The configuration settings are specified as attributes of the variable
"diag_script_info" (in NCL via "diag_script_info\@attribute = ...", see example
below). In order to activate these attributes, "diag_script_info" must be set
to "True".

**Example (NCL)**::

	  diag_script_info = True

	  diag_script_info@projection = "CylindricalEquidistant" ; map projection,
					                         ; e.g. Mollweide, Mercator
	  diag_script_info@styleset = "CMIP5" ; "CMIP5", "DEFAULT"
	  diag_script_info@colormap = "WhiteBlueGreenYellowRed" ; e.g., WhiteBlueGreenYellowRed,
								; rainbow
	  diag_script_info@ncdf = "default" ; enable to output to netCDF; either use "default"
					    ; or give a full file name

**Example (Python)**::

	  class diag_script:
	    def __init__(self):
	      self.info = True
	      self.color = 1
	      self.factor = 2.0e-3
	      self.name = "example"

	  diag_script_info = diag_script()

**Example (R)**::

	  diag_script_info<-new()
	  diag_script_info[["begin_ref_year"]]<-1970
	  diag_script_info[["end_ref_year"]]<-2000
	  diag_script_info[["timescale"]]<-3
	  diag_script_info[["seasons"]]<-c("ann", "djf", "mam", "jja", "son")
