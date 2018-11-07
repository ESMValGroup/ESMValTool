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
;    ncdf: full path and name for *.nc output
;
; Optional diag_script_info attributes
;    MyDiag_title: string for annotation
;
; Caveats
;
; Modification history
;    20181107-A_righ_ma: revised and extended to include new metadata handling
;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;                        of the code
;    20150325-A_laue_ax: modified reference tags used for acknowledgements
;                        (projects, observations, etc)
;    20140312-A_righ_ma: adapted to new time_operations structure
;    20130618-A_gott_kl: adjusted to r738 of
;                        https://svn.dlr.de/ESM-Diagnostic/sources/trunk
;    20121130-A_gott_kl: written
;
; #############################################################################

load "interface_scripts/interface.ncl"  ; load metadata

load "./diag_scripts/shared/statistics.ncl"  ; load statistics functions

load "./diag_scripts/shared/plot/style.ncl"  ; load plot style functions
load "./diag_scripts/shared/plot/contour_maps.ncl"  ; load plot function


begin

  ; _DIAG_SCRIPT is a global variables defined in interface.ncl
  enter_msg(_DIAG_SCRIPT, "")

  ; Shortcuts to current variable and field
  var0 = variable_info[0]@short_name
  field_type0 = variable_info[0]@field

end


begin

  ; The metadata for this diagnostic are defined in the input_file_info list
  ; in interface.ncl for all variables and datasets: use the convenience
  ; function to extract the metadata for a specific variable or dataset

  ; Retrieve metadata items by attribute (short_name and field)
  ; For multiple variables, use different indices of variable_info
  info = select_metadata_by_atts(\
    input_file_info, \
    [/[/"short_name", variable_info[0]@short_name/], \
      [/"field", variable_info[0]@field/] \
      /])

  ; Simple version using only variable's short name (work for most cases)
  ; info = select_metadata_by_name(input_file_info, var0)

  ; Read data for the first dataset in the list
  ; Loop over all ListCount(info) items to read all datasets in a loop
  dataset = read_data(info[0])
  printVarSummary(dataset)

  ; Read fx-variable for the second dataset in the list
  ; fx-variable are specified in the recipe, the field is returned as a list
  ; where also auxiliary coordinates are stored (if available).
  ; If no fx-variable is found a missing value is returned
  fx_var = read_fx_data(info[1], "sftlf")
  printVarSummary(fx_var[0])
  
  ; Retrieve metadata item for the reference dataset
  info_ref = select_metadata_by_atts( \
    input_file_info, \
    [/[/"short_name", variable_info[0]@short_name/], \
      [/"field", variable_info[0]@field/], \
      [/"dataset", variable_info[0]@reference_dataset/] \
      /])

  ; Read data for the reference dataset
  reference = read_data(info_ref[0])
  printVarSummary(reference)

  ; Calculate time-average of the reference dataset
  ; A function in esmaltool/diag_scripts/shared/statistics.ncl is used
  timeavg = time_operations(reference, -1, -1, "average", "annualclim", True)

  ; Set style dependent annotations
  annots = project_style(info, diag_script_info, "annots")

  print(annots)

end
