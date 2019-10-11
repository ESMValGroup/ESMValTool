; #############################################################################
; PLOT SCRIPTS FOR CONTOUR PLOTS
; #############################################################################
; Please consider using or extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;
;    function contour_map
;    function contour_map_polar
;    function contour_map_ce
;    procedure add_markers_to_map
;
; #############################################################################

load "$diag_scripts/../interface_scripts/auxiliary.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

load "$diag_scripts/shared/set_operators.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"

; #############################################################################
undef("contour_map")
function contour_map(wks_in[1],
                     source,
                     varname[1]: string)
;
; Arguments
;    wks_in: workstations (graphic object or default will be used).
;    source: data to be plotted or a NetCDF filename with data.
;    varname: variable name in the file.
;
; Source prototype
;    source[*,*]
;    source!0 = lat
;    source!1 = lon
;
; Return value
;    A graphic variable.
;
; Description
;    Wrapper for gsn_csm_contour_map.
;    Opens default wks, if not provided as argument of type "graphic".
;    Defines default resources, which are overridden by argument res.
;    Creates plot, according to wks & res.
;
; Caveats
;    Selection of defaults for res almost arbitrary
;    Please check results of all scripts that use this routine if modifying
;    the defaults!
;    Input via netCDF not yet implemented
;
; Modification history
;    20131104-gottschaldt_klaus-dirk: written.
;
local funcname, data, defaults, diag_script, var, wks, res, res_new
begin

  funcname = "contour_map"
  scriptname = "diag_scripts/shared/plot/contour_maps.ncl"
  enter_msg(scriptname, funcname)

  ; Get data, either directly or via netCDF file
  if (typeof(source) .eq. "string") then
    ; Function in ~/interface_scripts/auxiliary.ncl
    data = ncdf_read(source, varname)
  else
    data = source
    copy_VarMeta(source, data)
  end if

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname .eq.  defaults)) then  ; Use default
    ; Function in ~/interface_scripts/auxiliary.ncl
    var = att2var(data, "var")
  else  ; Use explicitly given name
    var = varname
    diag_script = att2var(data, "diag_script")
  end if
  ; Function in ~/interface_scripts/auxiliary.ncl
  diag_script = att2var(data, "diag_script")

  ; Check if a valid wks has been provided, otherwise invoke default
  wks = get_wks(wks_in, diag_script, var)  ; Function in aux_plotting.ncl

  ; Set default resources
  res = True
  res@gsnDraw           = True  ; Draw
  res@gsnFrame          = True  ; Advance frame
  res@tiMainString      = "Output of " + diag_script   ; Main title
  res@gsnLeftString     = ""               ; Annotation
  if (isatt(var, "long_name")) then
    res@gsnRightString    = var@long_name
  else
    res@gsnRightString    = var
  end if
  if (isatt(var, "units")) then
    res@gsnRightString    = res@gsnRightString + " in " + var@units
  end if

  ; Override defaults with "res_" attributes of "data"
  ; Function in ~/interface_scripts/auxiliary.ncl
  res_new = att2var(data, "res_")
  copy_VarMeta(res_new, res)

  ; Create plot
  map = gsn_csm_contour_map(wks, data, res)

  ; outfile name
  if (isatt(wks, "fullname")) then
    map@outfile = wks@fullname
  else
    map@outfile = wks@name
    error_msg("w", scriptname, funcname, "wks@fullname " + \
              "missing, consider to use get_wks to open wks")
  end if

  if (isatt(res, "gsnDraw")) then
    if(res@gsnDraw .eq. True) then
      log_info(" Wrote " + wks@name)
    end if
  else
    log_info(" Wrote " + wks@name)
  end if

  leave_msg(scriptname, funcname)
  return(map)

end

; #############################################################################
undef("contour_map_polar")
function contour_map_polar(wks_in[1],
                           source,
                           varname[1]: string)
;
; Arguments
;    wks_in: workstations (graphic object or default will be used).
;    source: data to be plotted or a NetCDF filename with data.
;    varname: variable name in the file.
;
; Source prototype
;    source[*,*]
;    source!0 = lat
;    source!1 = lon
;
; Return value
;    A graphic variable.
;
; Description
;    Wrapper for gsn_csm_contour_map_polar.
;    Opens default wks, if not provided as argument of type "graphic".
;    Defines default resources, which are overridden by argument res_in.
;    Creates plot according to wks & res, unless res_in@gsnDraw is set
;      to 'False' in your diag_script (to enable panelling, for example).
;
; Caveats
;    Selection of defaults for res almost arbitrary
;    Please check results of all scripts that use this routine if modifying
;    the defaults!
;    Input via netCDF not yet implemented
;
; Modification history
;    20140623-senftleben_daniel: now takes res as attributes of source.
;    20131218-senftleben_daniel: written.
;
local funcname, data, defaults, var, diag_script, wks, res, res_in, res1
begin

  funcname = "contour_map_polar"
  scriptname = "diag_scripts/shared/plot/contour_maps.ncl"
  enter_msg(scriptname, funcname)

  ; Get data, either directly or via netCDF file
  if (typeof(source) .eq. "string") then
    ; Function in ./interface_scripts/auxiliary.ncl
    data = ncdf_read(source, varname)
  else
    data = source
    copy_VarMeta(source, data)
  end if

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname .eq. defaults)) then  ; Use default
    ; Function in ./interface_scripts/auxiliary.ncl
    var = att2var(data, "var")
  else  ; Use explicitely given name
    var = varname
  end if
  if (isatt(data, "diag_script")) then
    diag_script = att2var(data, "diag_script")
  else
    error_msg("f", scriptname, funcname, "name of diag_script " + \
              "has to be passed as attribute of source")
  end if

  ; Check if a valid wks has been provided, otherwise invoke default
  wks = get_wks(wks_in, diag_script, var)

  ; Set default resources
  res              = True
  res@gsnDraw      = True  ; Draw
  res@gsnFrame     = True  ; Advance frame
  res@gsnLeftStringOrthogonalPosF = 0.1  ; shift the string up a bit
  res@tiMainString = "Output of " + diag_script  ; Main title
  res@cnFillOn     = True
  res@lbLabelBarOn = True
  if (isatt(var, "long_name")) then
    res@gsnRightString = var@long_name
  else
    res@gsnRightString = var
  end if
  if (isatt(var, "units")) then
    res@gsnRightString = res@gsnRightString + " in " + var@units
  end if

  ; Specify the hemisphere
  if (diag_script_info@region .eq. "Antarctic") then
    res@gsnPolar          = "SH"
    res@mpMaxLatF         = diag_script_info@max_lat  ; Default is 0
  else if (diag_script_info@region .eq. "Arctic")
    res@gsnPolarNH        = True
    res@mpMinLatF         = diag_script_info@min_lat  ; Default is 0
  else if (diag_script_info@region .eq. "global")
    error_msg("f", scriptname, funcname, "for global data use " + \
              "contour_map in ./diag_scripts/shared/plot/contour_maps.ncl")
  else
    error_msg("f", scriptname, funcname, "no (valid) region " + \
              "selection: " + diag_script_info@region)
  end if
  end if
  end if

  ; Override defaults with res_in set in diag_script
  ; Function in ~/interface_scripts/auxiliary.ncl
  res_in = att2var(data, "res_")
  copy_VarMeta(res_in, res)

  ; For sea ice: draw sea ice extent as a red line at sic = 0.15
  if (isatt(diag_script_info, "contour_extent")) then
    if (diag_script_info@contour_extent .eq. "y") then
      res1                       = True
      res1@gsnDraw               = False
      res1@gsnFrame              = False
      res1@cnLevelSelectionMode  = "ExplicitLevels"
      res1@cnLinesOn             = True
      res1@cnLineColor           = "red"
      res1@cnLineLabelsOn        = False
      res1@cnLevelFlag           = "LineOnly"
      res1@cnLevels              = (/0.15/)
      res1@cnFillOn              = False
      res1@cnInfoLabelOn         = False
      res1@cnLineThicknessF = 5.  ; Thicker contours
      areaplot = \
        gsn_csm_contour_map_overlay(wks_in, data, data, res, res1)
    else
      areaplot = gsn_csm_contour_map_polar(wks_in, data, res)
    end if
  else
    areaplot = gsn_csm_contour_map_polar(wks_in, data, res)
  end if

  leave_msg(scriptname, funcname)
  return(areaplot)

end

; #############################################################################
undef("contour_map_ce")
function contour_map_ce(wks_in[1],
                        source,
                        varname[1]: string)
;
; Arguments
;    wks_in: workstations (graphic object or default will be used).
;    source: data to be plotted or a NetCDF filename with data.
;    varname: variable name in the file.
;
; Source prototype
;    source[*,*]
;    source!0 = lat
;    source!1 = lon
;
; Return value
;    A graphic variable.
;
; Description
;    Wrapper for gsn_csm_contour_map_ce.
;    Opens default wks, if not provided as argument of type "graphic".
;    Defines default resources, which are overridden by argument res.
;    Creates plot, according to wks & res.
;
; Caveats
;    Selection of defaults for res almost arbitrary
;    Please check results of all scripts that use this routine if modifying
;    the defaults!
;    Input via netCDF not yet implemented
;
; Modification history
;    20161220-lauer_axel: added option to specificy plotting projection
;                         (via attribute data@projection)
;    20140228-righi_mattia: written.
;
local data, var, wks, res
begin

  funcname = "contour_map_ce"
  scriptname = "diag_scripts/shared/plot/contour_maps.ncl"
  enter_msg(scriptname, funcname)

  ; Get data, either directly or via netCDF
  if (typeof(source).eq."string") then
    data = ncdf_read(source, varname)
  else
    data = source
  end if

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname.eq.defaults)) then
    var = att2var(data, "var")
  else
    var = varname
  end if

  ; Check if a valid wks has been provided, otherwise invoke default
  wks = get_wks(wks_in, DIAG_SCRIPT, var)

  ; Default resources
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@cnFillOn = True
  res@tiMainFont = 22
  res@tiMainString = var
  if (isatt(data, "long_name")) then
    res@gsnLeftString = data@long_name
  else
    res@gsnLeftString = var
  end if
  if (isatt(data, "units")) then
    res@gsnLeftString  = res@gsnLeftString + \
      " [" + format_units(data@units) + "]"
  end if
  res@gsnCenterString = ""
  res@gsnRightString = ""

  ; Set right label
  if (isatt(data, "time_avg")) then
    res@gsnRightString = data@time_avg
  else
    res@gsnRightString = ""
  end if

  ; Override defaults with "res_" attributes of "data"
  res_new = att2var(data, "res_")
  copy_VarMeta(res_new, res)

  ; Create plot
  if (isatt(data, "projection")) then
    res@mpProjection = data@projection
    if (res@mpProjection.eq."Robinson") then
      res@mpPerimOn = False
      res@mpGridLineColor = -1
      res@mpGridAndLimbOn = True
      res@lbAutoManage = False
      res@lbLabelFontHeightF = 0.0125
    end if
    plot = gsn_csm_contour_map(wks_in, data, res)
  else
    plot = gsn_csm_contour_map_ce(wks_in, data, res)
  end if

  leave_msg(scriptname, funcname)
  return(plot)

end

; #############################################################################
undef("add_markers_to_map")
procedure add_markers_to_map(wks_in[1], plot[1]: graphic, lat[*]: numeric, \
                             lon[*]: numeric, data[*]: numeric)
;
; Arguments
;    wks_in: input workstation.
;    plot: a graphic object representing a contour plot.
;    lat: an array of latitude positions for the marker.
;    lon: an array of longitude positions for the marker.
;    data: an array of values to be overlaid as colored markers.
;
; Description
;    Overlays markers to an existing map plot, using the same color coding of
;    the associated labelbar.
;
; Caveats
;
; Modification history:
;    20140214-righi_mattia: written.
;
local funcname, class_name, views, cn_levels, lb_colors, cmap, colors, \
  loc_res, ndata, mfill, mbord, ii, idx, str
begin

  funcname = "add_markers_to_map"
  scriptname = "diag_scripts/shared/plot/contour_maps.ncl"
  enter_msg(scriptname, funcname)

  ; Check data consistency
  if (dimsizes(lat).ne.dimsizes(lon) .or. \
      dimsizes(data).ne.dimsizes(lon) .or. \
      dimsizes(data).ne.dimsizes(lat)) then
    error_msg("f", scriptname, funcname, "size of input not " + \
              "consistent")
  end if

  ; Check that the input plot is a contour plot and extract attributes
  class_name = NhlClassName(plot)
  if(class_name.eq."mapPlotClass") then
    if(isatt(plot, "contour")) then
      getvalues plot@contour
        "pmAnnoViews" : views
        "cnLevels" : cn_levels
      end getvalues
    end if
  else
    error_msg("f", scriptname, funcname, "input plot must be " + \
              "a contour plot")
  end if

  ; Extract labelbar colors
  labelbar = views(ind(NhlClassName(views).eq."labelBarClass"))
  getvalues labelbar
    "lbFillColors" : lb_colors
  end getvalues

  ; Extract associated colormap
  cmap = gsn_retrieve_colormap(wks_in)
  colors = cmap(lb_colors, :)

  ; Set marker resources
  loc_res = True
  loc_res@gsMarkerThicknessF = 1.5
  loc_res@gsMarkerSizeF      = 0.015
  ndata = dimsizes(data)
  mfill = new(ndata, graphic)
  mbord = new(ndata, graphic)

  ; Loop over input data and add marker to the plot
  do ii = 0, ndata - 1
    if (.not.ismissing(data(ii))) then

      ; Find value-color correspondance
      idx = min(ind(data(ii).lt.cn_levels))
      if (ismissing(idx)) then
        idx = dimsizes(cn_levels)
      end if

      ; Filled dot
      loc_res@gsMarkerIndex = 16
      loc_res@gsMarkerColor = colors(idx, :)
      mfill = gsn_add_polymarker(wks_in, plot, lon(ii), lat(ii), loc_res)
      delete(loc_res@gsMarkerColor)

      ; Hollow dot (as border)
      loc_res@gsMarkerIndex = 4
      loc_res@gsMarkerColor = "black"
      mbord = gsn_add_polymarker(wks_in, plot, lon(ii), lat(ii), loc_res)
      delete(loc_res@gsMarkerColor)

      ; Append graphics as plot attributes
      str = unique_string("mfill")
      plot@$str$ = mfill
      str = unique_string("mbord")
      plot@$str$ = mbord

    end if
  end do

end
