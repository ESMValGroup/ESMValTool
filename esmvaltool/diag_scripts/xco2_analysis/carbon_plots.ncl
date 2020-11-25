; #############################################################################
; PLOT SCRIPTS FOR SOME PLOTS OF XCO2 analysis recipe (gier20bg)
; #############################################################################
; Please consider using or extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;
;    function panel_plots
;
; #############################################################################

load "$diag_scripts/shared/set_operators.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/style.ncl"

; #############################################################################

undef("panel_plots")
function panel_plots(wks_in[1],
                     source1,
                     source2,
                     varname1[1]: string,
                     varname2[1]: string,
                     obs_pres,
                     items: list)
;
; Arguments;
;    wks_in: workstations (graphic object or default will be used).
;    source1: data to be plotted on x-axis or a NetCDF filename
;    source2: data to be plotted on y-axis or a NetCDF filename
;    varname: variable name, needed for netCDF files with multiple variables
;    obs_pres: True/False: is obs present (and first index)
;    items: list of input_file_info items
;
; Return value
;    A graphic variable.
;
; Description
;    Creates panels of x-y plots vertically arranging different models
;
; Caveats
;
; Modification history
;    20200227_gier_bettina: ported to v2
;    20180328_gier_bettina: written
;
local funcname, scriptname, verbosity, data
begin

    funcname = "panel_plots"
    scriptname = "diag_scripts/xco2_analysis/carbon_plots.ncl"
    enter_msg(scriptname, funcname)

    start_year = source1&time(0)
    end_year = source1&time(dimsizes(source1&time)-1)
    nyear = end_year - start_year + 1

    if isatt(source1, "trange") then
      trange = data@trange
    else if iscoord(source1, "time") then
      trange = tostring(start_year) + " - " + tostring(end_year)
    else
      trange = ""
    end if
    end if
    if isatt(source1, "region") then
      region = ", " + source1@region + ", "
    else
      region = ", "
    end if

    if isatt(source1, "long_name") then
      varname1_f = source1@long_name
    else
      varname1_f = varname1
    end if

    if isatt(source2, "long_name") then
      varname2_f = source2@long_name
    else
      varname2_f = varname2
    end if

    var_yaxis = False
    if isatt(source1, "var_yaxis") then
      var_yaxis = source1@var_yaxis
    elseif ((isatt(source2, "var_yaxis")) .and. (.not. var_yaxis)) then
      var_yaxis = source2@var_yaxis
    end if

    yaxisstr =  varname2_f + " [" + source2@units + "]"

    ; Check if a valid wks has been provided, otherwise invoke default
    wks = get_wks(wks_in, DIAG_SCRIPT, varname1)

    ; Plot histogram
    res = True
    res@tiYAxisString = ""
    res@tiYAxisFontHeightF = 0.04
    res@tmYLLabelFontHeightF = 0.03
    res@tiXAxisString = varname1_f + " [" + source1@units + "]"
    res@tiXAxisFontHeightF = 0.04
    res@tmXBLabelFontHeightF = 0.03
    res@pmTickMarkDisplayMode = "Conditional"
    res@tiMainFontHeightF = 0.02
    res@tiMainString = ""

    colors = project_style(items, diag_script_info, "colors")
    if obs_pres then
      colors := array_append_record("black", colors(1:), 0)
    else
      colors := colors(1:)
    end if
    ; Add multi-model mean as it's not in diag_script_info
    colors := array_append_record(colors, "red", 0)

    ; calc min and max dT,
    min_X = min(source1)
    max_X = max(source1)
    res@gsnMaximize = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@trXMinF = min_X - 0.05*(max_X - min_X)
    res@trXMaxF = max_X + 0.05*(max_X - min_X)
    res@trYMinF = min(source2) - 0.5 * (max(source2) - min(source2))
    res@trYMaxF = max(source2) + 0.5 * (max(source2) - min(source2))
    res@vpHeightF = 0.2
    res@vpWidthF = 0.8
    res@tmLabelAutoStride = True

    resM = True
    resM@gsMarkerSizeF = 0.01
    resM@gsMarkerThicknessF = 4.

    regression = new((/dimsizes(source2&model), 2/), float)
    cor = new((/dimsizes(source2&model), 2/), float)

    lres = True
    lres@gsLineThicknessF = 3.
    txres = True
    txres@txFontHeightF = 0.02

    modnames = source1&model
    if obs_pres then
      modnames(0) = "OBS"
    end if

    plots = new(dimsizes(source1&model), graphic)

    ; Override defaults with "res_" attributes of "data"
    res_new = att2var(source1, "res_")
    copy_VarMeta(res_new, res)

    ; Determine scale for variable y-axis
    if var_yaxis then
      my_y_scale = max(abs(dim_max_n(source2, 1) - dim_min_n(source2, 1)))
    end if

    ; Make each separate plot
    do imod = 0, dimsizes(source1&model)-1

      if imod.eq.toint(dimsizes(source1&model)/2) then
        res@tiYAxisString = yaxisstr
      else
        res@tiYAxisString = " "
      end if

      ; Variable y-Axis if wanted
      if var_yaxis then
        res@trYMinF = min(source2(imod, :)) - 0.3 * my_y_scale
        res@trYMaxF = min(source2(imod, :)) + 1.5 * my_y_scale
      end if

      ; Add Regression lines

      x_fit = fspan(res@trXMinF - 1, res@trXMaxF + 1, 100)

      rc = regline_stats(source1(imod, :), source2(imod, :))
      regression(imod, 0) = (/rc/)
      regression(imod, 1) = (/rc@stderr(1)/)
      cor(imod, 0) = rc@r  ; pearson-correlation coefficient
      cor(imod, 1) = rc@pval(1)  ; p-value

      y_fit = x_fit*rc + rc@yintercept

      lres@gsLineColor = colors(imod)
      res@xyLineColor = colors(imod)
      txres@txFontColor = colors(imod)

      plots(imod) = gsn_csm_xy(wks, x_fit, y_fit, res)

      if isatt(varname2, "mean") .and. varname2@mean then
        v2_mean = rc@yave
        v2_mean = decimalPlaces(v2_mean, 2, True)
        txadd = ", mean = " + sprintf("%1.2f", v2_mean)
        delete(v2_mean)
      else
        txadd = ""
      end if

      txres@txJust = "centerLeft"
      plots@$unique_string("mod")$ = \
        gsn_add_text(wks, plots(imod), modnames(imod) + ": " \
                     + sprintf("%1.3f", regression(imod, 0)) \
                     + " ~F18~S~F~ " \
                     + sprintf("%1.3f", regression(imod, 1)) + txadd, \
                     res@trXMinF + 0.02*(res@trXMaxF - res@trXMinF), \
                     res@trYMaxF - 0.13*(res@trYMaxF - res@trYMinF), txres)
      txres@txJust = "centerRight"
      plots@$unique_string("cor")$ = \
        gsn_add_text(wks, plots(imod), "r: " \
                     + sprintf("%1.2f", cor(imod, 0)) + ", p-value: " \
                     + sprintf("%1.2f", cor(imod, 1)), \
                     res@trXMaxF - 0.02*(res@trXMaxF - res@trXMinF),\
                     res@trYMinF + 0.13*(res@trYMaxF - res@trYMinF), txres)

      do iyr = 0, nyear-1
        if nyear.lt.15 then
          resM@gsMarkerIndex = iyr + 2
          stride = 1
        else
          if nyear.lt.50 then
            stride = 5
          elseif nyear.lt.100 then
            stride = 10
          else
            stride = 20
          end if
          resM@gsMarkerIndex = 2 + toint(iyr/stride)
        end if
        resM@gsMarkerColor = colors(imod)
        plots@$unique_string("markers")$ = \
          gsn_add_polymarker(wks, plots(imod), source1(imod, iyr), \
                             source2(imod, iyr), resM)
      end do

      delete(x_fit)
      delete(y_fit)

    end do

    ; Add Marker Legend
    lgres = True
    lgres@lgMonoMarkerIndex = False
    lgres@vpWidthF = 0.3
    lgres@vpHeightF = 0.5
    lgres@lgPerimOn = False
    lgres@lgMarkerSizeF = 0.01
    lgres@lgMonoMarkerColor = True
    lgres@lgMarkerColor = "black"
    lgres@lgAutoManage = False
    lgres@lgLabelFontHeightF = 0.06
    lgres@lgItemType = "Markers"

    if stride.eq.1 then
      nlabels = nyear
      labels = tostring(ispan(start_year, end_year, 1))
      lgres@lgMarkerIndexes = ispan(2, nyear+1, 1)
      lgres@lgItemOrder = ispan(nyear-1, 0, 1)
    else
      nlabels = toint((nyear)/stride) + 1
      labels = new(nlabels, string)
      labels(nlabels-1) = tostring(start_year + (nlabels-1)*stride) \
        + " - " + tostring(end_year)
      do lb = 0, nlabels - 2
        labels(lb) = tostring(start_year + lb*stride) + " - " + \
                     tostring(start_year + (lb+1)*stride)
      end do
      lgres@lgMarkerIndexes = ispan(2, 2+nlabels, 1)
      lgres@lgItemOrder = ispan(nlabels-1, 0, 1)
    end if

    legend = gsn_create_legend(wks, nlabels, labels, lgres)

    anres = True
    anres@amJust = "TopLeft"
    ; Parallel -> x-direction
    anres@amParallelPosF = 0.45
    anres@amOrthogonalPosF = -0.5
    annoid1 = gsn_add_annotation(plots(0), legend, anres)

    res1                     = True
    res2                     = True
    res1@gsnMaximize         = True
    res2@gsnAttachPlotsXAxis = True

    newplot = gsn_attach_plots(plots(0), plots(1:), res1, res2)

    draw(plots(0))

    frame(wks)

    leave_msg(scriptname, funcname)

    return(plots(0))
end
