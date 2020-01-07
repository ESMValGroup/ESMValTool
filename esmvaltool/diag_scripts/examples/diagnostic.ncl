; #############################################################################
; EXAMPLE SCRIPT FOR NCL DIAGNOSTICS
; #############################################################################
; Author: Mattia Righi (DLR, Germany), Klaus-Dirk Gottschaldt (DLR, Germany)
; #############################################################################
;
; Description
;    Simple script to demonstrate the implementation of a NCL diagnostic.
;    It shows how the select the metadata for specific variable(s) and
;    dataset(s) and how to read the corresponding fields.
;    The input is averaged over time to create a global contour plot.
;    Output to NetCDF files is optional (and controlled via diag_script_info)
;    Time averaging and plotting is done in separate routines.
;    All information necessary for plotting is added as attributes to the
;    variable containing the data. Attribute names "diag_script", "var"
;    and those starting with "res_" are reserved for that purpose.
;    The data variable (incl. metadata) can be passed to the plot script
;    directly (most efficient), or via a netCDF file (might be useful for
;    non-NCL routines).
;    A workstation ("wks" = file to plot to) will by default be generated
;    in the plotting routine. However, the default will be overridden, if
;    a "graphic" object is passed to the plotting routine. The latter
;    option is intended to be used in special cases, e.g. with
;    res_gsnDraw = False and res_gsnFrame = False.
;
; Required diag_script_info attributes
;    projection: map projection to use
;    colormap: color table to use
;    ncdf: full path and name for *.nc output
;
; Optional diag_script_info attributes
;    MyDiag_title: string for annotation
;
; Caveats
;
; Modification history
;    20181107-righi_mattia: revised and extended to include new metadata
;                           handling
;    20151027-lauer_axel: moved call to 'write_references' to the beginning
;                         of the code
;    20150325-lauer_axel: modified reference tags used for acknowledgements
;                         (projects, observations, etc)
;    20140312-righi_mattia: adapted to new time_operations structure
;    20121130-gottschaldt_klaus-dirk: written
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"  ; load metadata

load "$diag_scripts/shared/statistics.ncl"  ; load statistics functions

load "$diag_scripts/shared/plot/style.ncl"  ; load plot style functions
load "$diag_scripts/shared/plot/contour_maps.ncl"  ; load plot function

begin

  ; DIAG_SCRIPT is a global variables defined in interface.ncl
  enter_msg(DIAG_SCRIPT, "")

  ; Shortcuts to current variable and field
  ; For multiple variables, use different indices of variable_info
  var0 = variable_info[0]@short_name
  info0 = select_metadata_by_name(input_file_info, var0)

  ; The number of variables is the size of the variable_info list
  n_variables = ListCount(variable_info)

end


begin

  ; The metadata for this diagnostic are defined in the input_file_info list
  ; in interface.ncl for all variables and datasets: use the convenience
  ; function to extract the metadata for a specific variable or dataset

  ; Retrieve metadata items by attribute (short_name and field)
  atts = True
  atts@short_name = var0
  info = select_metadata_by_atts(input_file_info, atts)
  delete(atts)

  ; Simple version using only the variable's short name (works for most cases)
  info_simple = select_metadata_by_name(input_file_info, var0)

  ; The number of datasets for var0 is the size of the returned list
  n_datasets = ListCount(info)

  ; Read data for the first dataset in the list
  ; Loop over all n_datasets info items to read all datasets in a loop
  dataset = read_data(info[0])
  printVarSummary(dataset)

  ; Retrieve metadata item for the reference dataset
  atts = True
  atts@short_name = var0
  atts@dataset = variable_info[0]@reference_dataset
  info_ref = select_metadata_by_atts(input_file_info, atts)
  delete(atts)

  ; Read data for the reference dataset
  reference = read_data(info_ref[0])
  printVarSummary(reference)

  ; Calculate time-average of the reference dataset
  ; A function in esmaltool/diag_scripts/shared/statistics.ncl is used
  timeavg = time_operations(reference, -1, -1, "average", "annualclim", True)

  ; Set annotations
  annots = project_style(info_ref, diag_script_info, "annots")

  ; Plot resources
  ; Add to timeavg, as attributes with prefix "res_".
  ; Provide attributes only that differ from the defaults given at
  ; http://www.ncl.ucar.edu/Document/Graphics/Resources/list_alpha_res.shtml
  timeavg@res_gsnMaximize       = True   ; Use full page for the plot
  timeavg@res_cnFillOn          = True   ; Color plot desired
  timeavg@res_cnLineLabelsOn    = False  ; Contour lines
  timeavg@res_cnLinesOn         = False
  ; Colors
  ; http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
  timeavg@res_cnFillPalette     = diag_script_info@colormap  ; Set colormap
  ; Annotation
  timeavg@res_tiMainString      = "NCL example diagnostic"
  timeavg@res_gsnLeftString     = annots
  timeavg@res_gsnRightString    = variable_info[0]@long_name + " in " + \
    variable_info[0]@units
  ; Specified in ./nml/cfg_MyDiag/cfg_MyDiag.ncl
  timeavg@res_mpProjection      = diag_script_info@projection
  ; Set explicit contour levels
  if(isatt(diag_script_info, "explicit_cn_levels")) then
    timeavg@res_cnLevelSelectionMode = "ExplicitLevels"
    timeavg@res_cnLevels = diag_script_info@explicit_cn_levels
  end if

  ; Other Metadata: diag_script, var
  ; Add to timeavg, as attributes without prefix
  if(isatt(timeavg, "diag_script")) then  ; Add to existing entries
    temp = timeavg@diag_script
    delete(timeavg@diag_script)
    timeavg@diag_script = array_append_record(temp, (/diag_script/), 0)
    delete(temp)
  else  ; Add as new attribute
    timeavg@diag_script = (/DIAG_SCRIPT/)
  end if
  timeavg@var = var0  ; Overwrite existing entry
  if(isatt(variable_info, "long_name")) then
    timeavg@var_long_name = variable_info[0]@long_name
  end if
  if(isatt(variable_info, "units")) then
    timeavg@var_units = variable_info[0]@units
  end if

  ; Optional output to netCDF
  ; (controlled by diag_script_info)
  if(isatt(diag_script_info, "ncdf")) then
    timeavg@ncdf = diag_script_info@ncdf
    ; Function in ~/interface_scripts/auxiliary.ncl
    ncdf_outfile = ncdf_write(timeavg, diag_script_info@ncdf)
  end if

  ; Create the plot
  ; There are two options to pass data to the plotting routine: directly
  ; or via a netCDF file. All information needed to create a graphics file
  ; are passed with the data. It is possible to provide "dummy_for_wks" and
  ; "dummy_for_var" as arguments in the call.
  ; You may invoke a wks here, and/or provide the variable name (see example
  ; in next section). The routine "contour_map" is a function in
  ; ~/diag_scripts/shared/plot/aux_plotting.ncl
  ;
  ; Option 1: pass data directly
  ; -> may only work with NCL (i.e. ok here)
  ; Dummies trigger defaults
  map = contour_map("dummy_for_wks", timeavg, "dummy_for_var")
  ; Option 2: pass data via netCDF file
  ; -> useful for non-NCL routines & metrics
  ; map = contour_map("dummy_for_wks", ncdf_outfile,"dummy_for_var")

  ; Call provenance logger
  log_provenance(ncdf_outfile, \
                 map@outfile, \
                 "Mean of variable: " + var0, \
                 "mean", \
                 "global", \
                 "geo", \
                 (/"righi_mattia", "gottschaldt_klaus-dirk"/), \
                 (/"acknow_author"/), \
                 metadata_att_as_array(info0, "filename"))

end
