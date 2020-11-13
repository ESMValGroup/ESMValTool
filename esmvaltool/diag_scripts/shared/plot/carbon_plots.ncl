; #######################################################################
; carbon_plots.ncl
; #######################################################################
; Plots:
; * single_barchart         : one value per model
;
; Written by Sabrina Wenzel (sabrina.zechlau@dlr.de)
; #######################################################################
load "$diag_scripts/../interface_scripts/auxiliary.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
; #######################################################################

undef("barchart")
function barchart(wks_in[1], \
                  source, \
                  varname, \
                  inlist[1])
; Arguments:
;  * wks     : workstation for were to draw the plot
;  * data    : array for plotting on Y-axis
;  * varname : variable to be plotted
;  * inlist  : input_file_info

local plot, idx_obs, idx_mod, sres, dim_MOD, colors, data
begin

  ; Check for valid data
  funcname = "barchart"
  scriptname = "diag_scripts/shared/plot/carbon_plots.ncl"
  enter_msg(scriptname, funcname)

  ; Get data, either directly or via netCDF
  if (typeof(source).eq."string") then
    data = ncdf_read(source, varname)
  else
    data = source
  end if

  ; Style options
  colors  = project_style(inlist, diag_script_info, "colors")

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname.eq.defaults)) then
    var = att2var(data, "var")
  else
    var = varname
  end if

  ; Check if a valid wks has been provided, otherwise invoke default
  wks = get_wks(wks_in, DIAG_SCRIPT, var)
  gsn_define_colormap(wks, array_append_record((/"white", "black"/), \
                      colors, 0))

  ; Define names for x-ticks and dimsizes
  dim_MOD  = dimsizes(data(:, 0))

  XStg = data&datasets(0)
  xValue = ispan(1, dim_MOD, 1)

  ; Resources
  res = True
  res@gsnDraw       = False
  res@gsnFrame      = False
  res@gsnMaximize   = True
  res@gsnXYBarChart = True
  res@gsnXYBarChartBarWidth = 0.30

  res@trXMinF   = 1-0.5
  res@trXMaxF   = dim_MOD+0.5

  res@vpWidthF  = 0.5
  res@vpHeightF = 0.4

  res@tmXBMode             = "Explicit"  ; explicit labels
  res@tmXBValues           = xValue
  res@tmXBLabels           = data&datasets
  res@tmXBLabelFontHeightF = 0.025
  res@tmXBLabelAngleF      = 30
  res@tmXBLabelJust        = "TopRight"
  res@tmXTLabelFontHeightF = 0.020
  res@tmYLLabelFontHeightF = 0.025

  res@tiMainFontHeightF  = 0.025
  res@tiXAxisFontHeightF = 0.025
  res@tiMainString       = var
  res@tiYAxisString      = data&yaxis
  res@gsnYRefLine        = 0.
  res@gsnXYBarChartColors2  = ispan(2, dimsizes(colors) + 1, 1)

  ; Override defaults with "res_" attributes of "data"
  res_new = att2var(data, "res_")
  copy_VarMeta(res_new, res)

  plot = gsn_csm_xy(wks, xValue, data(:, 0), res)

  return(plot)

end
