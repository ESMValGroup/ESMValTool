;;#############################################################################
;; GENERAL ROUTINES FOR MANAGING PLOT RESOURCES
;;#############################################################################
;; Please consider using or extending existing routines before adding new ones.
;; Check the header of each routine for documentation.      
;;
;; Contents:
;;     function panel_two_by_one
;;     function panel_three_by_one
;;     function panel_n_by_cols
;;     function add_blank_plot_title
;;     procedure add_explicit_title_font_size_if_defined
;;     procedure three_by_one_labelbar
;;     procedure adjust_lbstride_if_fonth_lt
;;     procedure two_by_one_labelbar
;;     procedure n_by_cols_labelbar
;;     procedure inset_labelbar
;;
;;#############################################################################
begin
    loadscript(exact_panel_position_file)
end
load "diag_scripts/aux/SAMonsoon/exact_panel_positions_generic.ncl"
load "./interface_scripts/messaging.ncl"

undef("add_explicit_title_font_size_if_defined")
procedure add_explicit_title_font_size_if_defined(resource [1]:logical,\
                                                  res_string [1]: string)
;; res  -- Input resource
;; id   -- String defining title font height
;; Description:
;;   Checks if panel_position file also defines title font height, if so applies it
begin
    if (isdefined(res_string)) then
        if (res_string .eq. "n_by_cols_1_tiMainFontHeightF") then
            resource@tiMainFontHeightF = n_by_cols_1_tiMainFontHeightF
            resource@tiMainOffsetYF = n_by_cols_1_tiMainOffsetYF
        else if (res_string .eq. "n_by_cols_2_tiMainFontHeightF") then
            resource@tiMainFontHeightF = n_by_cols_2_tiMainFontHeightF
            resource@tiMainOffsetYF = n_by_cols_2_tiMainOffsetYF
        else if (res_string .eq. "n_by_cols_3_tiMainFontHeightF") then
            resource@tiMainFontHeightF = n_by_cols_3_tiMainFontHeightF
            resource@tiMainOffsetYF = n_by_cols_3_tiMainOffsetYF
        else if (res_string .eq. "n_by_cols_3_tmFontHeightF") then
            resource@tmXBLabelFontHeightF = n_by_cols_3_tmFontHeightF
            resource@tmYLLabelFontHeightF = n_by_cols_3_tmFontHeightF
        end if
        end if
        end if
        end if
    end if
end

undef("panel_two_by_one")
function panel_two_by_one(res [1] : logical,
                          id  [1] : integer)
;;                 return_val [1] : logical
;; res  -- Input resource
;; id   -- Index to plot
;; Description:
;;   Updates the current plot resource with settings specific for plot 'id'
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering panel_two_by_one.ncl", verbosity, 6)
    res@vpKeepAspect = True
    res@gsnStringFontHeightF = 0.011
    if (id .eq. 0) then
        res@tmXBLabelsOn = True
        res@lbLabelBarOn = False
        res@vpXF      = two_by_one_first_vpXF
        res@vpYF      = two_by_one_first_vpYF
        res@vpWidthF  = two_by_one_first_vpWidthF
        res@vpHeightF = two_by_one_first_vpHeightF
    else if (id .eq. 1) then
        res@tmXBLabelsOn = True
        res@lbLabelBarOn = True
        res@vpXF      = two_by_one_second_vpXF
        res@vpYF      = two_by_one_second_vpYF
        res@vpWidthF  = two_by_one_second_vpWidthF
        res@vpHeightF = two_by_one_second_vpHeightF
    else
        print("fatal: panel_two_by_one() - Invalid value for id=" + id)
        status_exit(1)
    end if
    end if
    info_output(">>>>>>>> Leaving panel_two_by_one.ncl", verbosity, 6)
    return res
end

undef("panel_three_by_one")
function panel_three_by_one(res [1] : logical,
                            id  [1] : integer)
;;                   return_val [1] : logical
;; res  -- Input resource
;; id   -- Index to plot
;; Description:
;;   Updates the current plot resource with settings specific for plot 'id'
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering panel_three_by_one.ncl", verbosity, 6)
    res@vpKeepAspect = True
    res@gsnStringFontHeightF = 0.011
    if (id .eq. 0) then
        res@tmXBLabelsOn = False
        res@lbLabelBarOn = False
        res@vpXF      = three_by_one_first_vpXF
        res@vpYF      = three_by_one_first_vpYF
        res@vpWidthF  = three_by_one_first_vpWidthF
        res@vpHeightF = three_by_one_first_vpHeightF
    else if (id .eq. 1) then
        res@tmXBLabelsOn = False
        res@lbLabelBarOn = False
        res@vpXF      = three_by_one_second_vpXF
        res@vpYF      = three_by_one_second_vpYF
        res@vpWidthF  = three_by_one_second_vpWidthF
        res@vpHeightF = three_by_one_second_vpHeightF
    else if (id .eq. 2) then
        res@tmXBLabelsOn = True
        res@lbLabelBarOn = True
        res@lbTopMarginF = 0.3
        res@vpXF      = three_by_one_third_vpXF
        res@vpYF      = three_by_one_third_vpYF
        res@vpWidthF  = three_by_one_third_vpWidthF
        res@vpHeightF = three_by_one_third_vpHeightF
    else
        error_msg("fatal", "SAMonsoon_panel_positions.ncl", "panel_three_by_one", \
                  "Invalid value for id=" + id)
    end if
    end if
    end if
    info_output(">>>>>>>> Leaving panel_three_by_one.ncl", verbosity, 6)
    return res
end

undef("three_by_one_labelbar")
procedure three_by_one_labelbar(wks[1]:graphic,
                                plot[1]:graphic,
                                settings[1]:logical)
local amres, annoid, contour, lbid,
    lblabels, lbres,  nboxes, resource_handle
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering three_by_one_labelbar.ncl", verbosity, 6)

    if (isatt(plot, "contour")) then
        resource_handle = plot@contour
    else
        resource_handle = plot
    end if

    getvalues resource_handle
        "cnLevels" : cnLevels
        "cnFillColors" : cnFillColors
    end getvalues
    lbres                   = True
    lbres@lbAutoManage      = True
    lbres@lbFillColors      = cnFillColors
    lbres@lbMonoFillPattern = True
    lbres@lbOrientation     = "Vertical"
    lbres@lbPerimOn         = False
    lbres@lbTitleOffsetF    = 0.0
    lbres@lbTitleExtentF    = 0.1
    lbres@gsnDraw           = True
    lbres@gsnMaximize       = True
    lbres@lbJustification   = "topleft"
    lbres@lbLeftMarginF     = 0.0
    lbres@lbTopMarginF      = 0.0
    lbres@lbLabelPosition   = "Right"
    lbres@lbBoxMinorExtentF = 0.3
    lbres@lbTitleString     = settings@lb_title
    lbres@lbTitlePosition   = "Top"
    lbres@lbTitleAngleF     = 0
    lbres@vpWidthF          = 0.09
    lbres@lbLabelAlignment  = "InteriorEdges"
    lbres@vpHeightF         = three_by_one_first_vpYF - three_by_one_second_vpYF + three_by_one_second_vpHeightF

    nboxes = dimsizes(cnLevels)
    lblabels = sprintf("%5.2f", cnLevels)
    lbres@gsnDraw = True
    lbres@gsnFrame = False
    lbid = gsn_create_labelbar(wks, nboxes, lblabels, lbres)

    hor_displacement = 0.0
    if(isatt(settings, "lb_hor_displacement")) then
        hor_displacement = hor_displacement + settings@lb_hor_displacement
    end if
    amres = True
    amres@amZone           = 2
    amres@amSide           = "Right"
    amres@amParallelPosF   = three_by_one_lb_amres_adjust  ; Move down
    amres@amOrthogonalPosF = 0.05 + hor_displacement ; Center labelbar
    annoid = gsn_add_annotation(plot, lbid, amres)
    info_output(">>>>>>>> Leaving three_by_one_labelbar.ncl", verbosity, 6)
end

undef("adjust_lbstride_if_fonth_lt")
procedure adjust_lbstride_if_fonth_lt(lbid[1]:graphic,
                                      fontheight_min[1]:float,
                                      new_lbstride[1]:integer)
; Adjust lbLabelStride if font is too small (with lbAutoManage=True
; this will increase the font height
begin
    getvalues lbid
        "lbLabelFontHeightF" : lblabelfontheightf
    end getvalues
    if (lblabelfontheightf .lt. 0.003) then
        setvalues lbid
            "lbLabelStride" : 3
        end setvalues
    else if (lblabelfontheightf .lt. 0.007) then
        setvalues lbid
            "lbLabelStride" : 2
        end setvalues
    end if
    end if
end

undef("two_by_one_labelbar")
procedure two_by_one_labelbar(wks, plot, title)
local amres, annoid, contour, lbid,
    lblabels, lbres,  nboxes, resource_handle,
    fontheight_min, new_lbstride
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering two_by_one_labelbar.ncl", verbosity, 6)

    if (isatt(plot, "contour")) then
        resource_handle = plot@contour
    else
        resource_handle = plot
    end if

    getvalues resource_handle
        "cnLevels" : cnLevels
        "cnFillColors" : cnFillColors
    end getvalues
    lbres                   = True
    lbres@gsnDraw           = True
    lbres@gsnMaximize       = True
    lbres@lbAutoManage      = True
    lbres@lbBoxMinorExtentF = 0.3
    lbres@lbFillColors      = cnFillColors
    lbres@lbJustification   = "BottomLeft"
    lbres@lbLabelAlignment  = "InteriorEdges"
    lbres@lbLabelPosition   = "Right"
    lbres@lbLeftMarginF     = 0.0
    lbres@lbMonoFillPattern = True
    lbres@lbOrientation     = "Horizontal"
    lbres@lbPerimOn         = False
    lbres@lbTitleAngleF     = 0
    lbres@lbTitleExtentF    = 0.15
    lbres@lbTitleOffsetF    = -0.05
    lbres@lbTitlePosition   = "Bottom"
    lbres@lbTitleString     = title
    lbres@lbTopMarginF      = 0.0
    lbres@vpHeightF         = 0.1
    lbres@vpWidthF          = two_by_one_second_vpWidthF

    nboxes = dimsizes(cnLevels)
    lblabels = sprintf("%5.2f", cnLevels)
    lbres@gsnDraw = True
    lbres@gsnFrame = False
    lbid = gsn_create_labelbar(wks, nboxes, lblabels, lbres)

    ;; If font height seems small, increase skip some labels
    fontheight_min = 0.007
    new_lbstride = 2
    adjust_lbstride_if_fonth_lt(lbid, fontheight_min, new_lbstride)

    amres = True
    amres@amZone           = 2
    amres@amSide           = "Bottom"
    amres@amJust           = "CenterLeft"
    amres@amParallelPosF   = 0.0  ; Move right
    amres@amOrthogonalPosF = two_by_one_lb_amres_adjust  ; Move down
    annoid = gsn_add_annotation(plot, lbid, amres)
    info_output(">>>>>>>> Leaving two_by_one_labelbar.ncl", verbosity, 6)
end


undef("panel_n_by_cols")
function panel_n_by_cols(res [1] : logical,
                         id  [1] : integer,
                         rows  [1] : float,
                         cols [1] : integer,
                         no_figs [1] : integer)
;;                return_val [1] : logical
;; res  -- Input resource
;; id   -- Index to plot
;; rows -- Number of models/cols (=a float)
;; cols -- How many columns plot have
;; no_figs -- Number of figures on this page
;; Description:
;;   Updates the current plot resource with settings specific for plot 'id'
local verbosity, font_height_string, rows_int
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering panel_n_by_cols.ncl", verbosity, 6)

    font_title_height = "n_by_cols_" + cols + "_tiMainFontHeightF"
    font_label_height = "n_by_cols_" + cols + "_tmFontHeightF"
    add_explicit_title_font_size_if_defined(res, font_title_height)
    add_explicit_title_font_size_if_defined(res, font_label_height)

    ;; Add xlabels or not...
    if (id + 1 .gt. no_figs - cols) then
        res@tmXBLabelsOn = True
    else
        res@tmXBLabelsOn = False
    end if
    res@gsnStringFontHeightF = 0.011

    rows_int = tointeger(ceil(rows))
    res@vpKeepAspect = True
    if (cols .eq. 1) then
        res@tmYLLabelsOn = True
        res@lbLabelBarOn = False
        if (id .eq. 0) then
            res@vpXF      = n_by_cols_1_first_vpXF
            res@vpYF      = n_by_cols_1_first_vpYF
            res@vpWidthF  = n_by_cols_1_first_vpWidthF
            res@vpHeightF = n_by_cols_1_first_vpHeightF
        else if (id .eq. 1) then
            res@vpXF      = n_by_cols_1_second_vpXF
            res@vpYF      = n_by_cols_1_second_vpYF
            res@vpWidthF  = n_by_cols_1_second_vpWidthF
            res@vpHeightF = n_by_cols_1_second_vpHeightF
        else if (id .eq. 2) then
            res@vpXF      = n_by_cols_1_third_vpXF
            res@vpYF      = n_by_cols_1_third_vpYF
            res@vpWidthF  = n_by_cols_1_third_vpWidthF
            res@vpHeightF = n_by_cols_1_third_vpHeightF
        else if (id .eq. 3) then
            res@vpXF      = n_by_cols_1_fourth_vpXF
            res@vpYF      = n_by_cols_1_fourth_vpYF
            res@vpWidthF  = n_by_cols_1_fourth_vpWidthF
            res@vpHeightF = n_by_cols_1_fourth_vpHeightF
        else
            error_msg("fatal", "SAMonsoon_panel_positions.ncl", \
                      "panel_n_by_cols", "Invalid value for id=" + id)
        end if
        end if
        end if
        end if
    else if (cols .eq. 2) then
        if (id .eq. 0) then
            res@tmYLLabelsOn = True
            res@tiYAxisOn    = True
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_2_first_vpXF
            res@vpYF      = n_by_cols_2_first_vpYF
            res@vpWidthF  = n_by_cols_2_first_vpWidthF
            res@vpHeightF = n_by_cols_2_first_vpHeightF
        else if (id .eq. 1) then
            res@tmYLLabelsOn = False
            res@tiYAxisOn    = False
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_2_second_vpXF
            res@vpYF      = n_by_cols_2_second_vpYF
            res@vpWidthF  = n_by_cols_2_second_vpWidthF
            res@vpHeightF = n_by_cols_2_second_vpHeightF
        else if (id .eq. 2) then
            res@tmYLLabelsOn = True
            res@tiYAxisOn    = True
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_2_third_vpXF
            res@vpYF      = n_by_cols_2_third_vpYF
            res@vpWidthF  = n_by_cols_2_third_vpWidthF
            res@vpHeightF = n_by_cols_2_third_vpHeightF
        else if (id .eq. 3) then
            res@tmYLLabelsOn = False
            res@tiYAxisOn    = False
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_2_fourth_vpXF
            res@vpYF      = n_by_cols_2_fourth_vpYF
            res@vpWidthF  = n_by_cols_2_fourth_vpWidthF
            res@vpHeightF = n_by_cols_2_fourth_vpHeightF
        else if (id .eq. 4) then
            res@tmYLLabelsOn = True
            res@tiYAxisOn    = True
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_2_fifth_vpXF
            res@vpYF      = n_by_cols_2_fifth_vpYF
            res@vpWidthF  = n_by_cols_2_fifth_vpWidthF
            res@vpHeightF = n_by_cols_2_fifth_vpHeightF
        else if (id .eq. 5) then
            res@tmYLLabelsOn = False
            res@tiYAxisOn    = False
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_2_sixth_vpXF
            res@vpYF      = n_by_cols_2_sixth_vpYF
            res@vpWidthF  = n_by_cols_2_sixth_vpWidthF
            res@vpHeightF = n_by_cols_2_sixth_vpHeightF
        else if (id .eq. 6) then
            res@tmYLLabelsOn = True
            res@tiYAxisOn    = True
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_2_seventh_vpXF
            res@vpYF      = n_by_cols_2_seventh_vpYF
            res@vpWidthF  = n_by_cols_2_seventh_vpWidthF
            res@vpHeightF = n_by_cols_2_seventh_vpHeightF
        else if (id .eq. 7) then
            res@tmYLLabelsOn = False
            res@tiYAxisOn    = False
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_2_eighth_vpXF
            res@vpYF      = n_by_cols_2_eighth_vpYF
            res@vpWidthF  = n_by_cols_2_eighth_vpWidthF
            res@vpHeightF = n_by_cols_2_eighth_vpHeightF
        else if (id .eq. 8) then
            res@tmYLLabelsOn = True
            res@tiYAxisOn    = True
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_2_ninth_vpXF
            res@vpYF      = n_by_cols_2_ninth_vpYF
            res@vpWidthF  = n_by_cols_2_ninth_vpWidthF
            res@vpHeightF = n_by_cols_2_ninth_vpHeightF
        else if (id .eq. 9) then
            res@tmYLLabelsOn = False
            res@tiYAxisOn    = False
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_2_tenth_vpXF
            res@vpYF      = n_by_cols_2_tenth_vpYF
            res@vpWidthF  = n_by_cols_2_tenth_vpWidthF
            res@vpHeightF = n_by_cols_2_tenth_vpHeightF
        else
            error_msg("fatal", "SAMonsoon_panel_positions.ncl", \
                      "panel_n_by_cols", "Invalid value for id=" + id)
        end if
        end if
        end if
        end if
        end if
        end if
        end if
        end if
        end if
        end if

    else if (cols .eq. 3) then
        res@gsnStringFontHeightF = 0.008
        if (id .eq. 0) then
            res@tmYLLabelsOn = True
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_3_first_vpXF
            res@vpYF      = n_by_cols_3_first_vpYF
            res@vpWidthF  = n_by_cols_3_first_vpWidthF
            res@vpHeightF = n_by_cols_3_first_vpHeightF
        else if (id .eq. 1) then
            res@tmYLLabelsOn = False
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_3_second_vpXF
            res@vpYF      = n_by_cols_3_second_vpYF
            res@vpWidthF  = n_by_cols_3_second_vpWidthF
            res@vpHeightF = n_by_cols_3_second_vpHeightF
        else if (id .eq. 2) then
            res@tmYLLabelsOn = False
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_3_third_vpXF
            res@vpYF      = n_by_cols_3_third_vpYF
            res@vpWidthF  = n_by_cols_3_third_vpWidthF
            res@vpHeightF = n_by_cols_3_third_vpHeightF
        else if (id .eq. 3) then
            res@tmYLLabelsOn = True
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_3_fourth_vpXF
            res@vpYF      = n_by_cols_3_fourth_vpYF
            res@vpWidthF  = n_by_cols_3_fourth_vpWidthF
            res@vpHeightF = n_by_cols_3_fourth_vpHeightF
        else if (id .eq. 4) then
            res@tmYLLabelsOn = False
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_3_fifth_vpXF
            res@vpYF      = n_by_cols_3_fifth_vpYF
            res@vpWidthF  = n_by_cols_3_fifth_vpWidthF
            res@vpHeightF = n_by_cols_3_fifth_vpHeightF
        else if (id .eq. 5) then
            res@tmYLLabelsOn = False
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_3_sixth_vpXF
            res@vpYF      = n_by_cols_3_sixth_vpYF
            res@vpWidthF  = n_by_cols_3_sixth_vpWidthF
            res@vpHeightF = n_by_cols_3_sixth_vpHeightF
        else if (id .eq. 6) then
            res@tmYLLabelsOn = True
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_3_seventh_vpXF
            res@vpYF      = n_by_cols_3_seventh_vpYF
            res@vpWidthF  = n_by_cols_3_seventh_vpWidthF
            res@vpHeightF = n_by_cols_3_seventh_vpHeightF
        else if (id .eq. 7) then
            res@tmYLLabelsOn = False
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_3_eigth_vpXF
            res@vpYF      = n_by_cols_3_eigth_vpYF
            res@vpWidthF  = n_by_cols_3_eigth_vpWidthF
            res@vpHeightF = n_by_cols_3_eigth_vpHeightF
        else if (id .eq. 8) then
            res@tmYLLabelsOn = False
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_3_ninth_vpXF
            res@vpYF      = n_by_cols_3_ninth_vpYF
            res@vpWidthF  = n_by_cols_3_ninth_vpWidthF
            res@vpHeightF = n_by_cols_3_ninth_vpHeightF
        else if (id .eq. 9) then
            res@tmYLLabelsOn = True
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_3_tenth_vpXF
            res@vpYF      = n_by_cols_3_tenth_vpYF
            res@vpWidthF  = n_by_cols_3_tenth_vpWidthF
            res@vpHeightF = n_by_cols_3_tenth_vpHeightF
        else if (id .eq. 10) then
            res@tmYLLabelsOn = False
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_3_eleventh_vpXF
            res@vpYF      = n_by_cols_3_eleventh_vpYF
            res@vpWidthF  = n_by_cols_3_eleventh_vpWidthF
            res@vpHeightF = n_by_cols_3_eleventh_vpHeightF
        else if (id .eq. 11) then
            res@tmYLLabelsOn = False
            res@lbLabelBarOn = False
            res@vpXF      = n_by_cols_3_twelveth_vpXF
            res@vpYF      = n_by_cols_3_twelveth_vpYF
            res@vpWidthF  = n_by_cols_3_twelveth_vpWidthF
            res@vpHeightF = n_by_cols_3_twelveth_vpHeightF
        else
            error_msg("fatal", "SAMonsoon_panel_positions.ncl", \
                      "panel_n_by_cols", "Invalid value for id=" + id)
        end if
        end if
        end if
        end if
        end if
        end if
        end if
        end if
        end if
        end if
        end if
        end if
    else
        error_msg("fatal", "SAMonsoon_panel_positions.ncl", \
                  "panel_n_by_cols", "Invalid value for cols=" + cols)
    end if
    end if
    end if
    info_output(">>>>>>>> Leaving panel_n_by_cols.ncl", verbosity, 6)
    return res
end

undef("n_by_cols_labelbar")
procedure n_by_cols_labelbar(wks, plot, plot_res, rows, cols, title, plottype_lbres)
local verbosity, resource_handle, lbres, offset, nboxes, max_cnfill_idx,\
    plotsettings_prefix
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering n_by_cols_labelbar.ncl", verbosity, 6)
    offset = 0.0

    ;; Labelbar resource
    lbres                   = True

    if (isatt(plot_res, "contour")) then
        resource_handle = plot_res@contour
    else
        resource_handle = plot_res
    end if

    getvalues resource_handle
        "cnLevels" : cnLevels
        "cnFillColors" : cnFillColors
    end getvalues

    ;; Copy the end cnFillColor to end-1
    if (isatt(plottype_lbres, "shift_down_maxcolor")) then
        if (plottype_lbres@shift_down_maxcolor) then
            max_cnfill_idx = dimsizes(cnFillColors)
            cnFillColors(max_cnfill_idx - 2) = cnFillColors(max_cnfill_idx - 1)
        end if
    end if

    ;; Manual update of lb-labels?
    if (isatt(plottype_lbres, "lblabels")) then
        lbres@lbLabelStrings = plottype_lbres@lblabels
    end if
    ;; Default labelbar resource settings
    lbres@gsnMaximize       = True
    lbres@lbAutoManage      = True
    lbres@lbBoxMinorExtentF = 0.2
    lbres@lbFillColors      = cnFillColors
    lbres@lbJustification   = "TopRight"
    lbres@lbLabelAlignment  = "InteriorEdges"
    lbres@lbLabelPosition   = "Bottom"
    lbres@lbLeftMarginF     = 0.0
    lbres@lbMonoFillPattern = True
    lbres@lbOrientation     = "Horizontal"
    lbres@lbPerimOn         = False
    lbres@lbTitleAngleF     = 0
    lbres@lbTitleExtentF    = 0.3
    lbres@lbTitleOffsetF    = 0.0
    lbres@lbTitlePosition   = "Bottom"
    lbres@lbTitleString     = title
    lbres@lbTopMarginF      = 0.0
    lbres@vpHeightF          = 0.06
    if (cols .eq. 1) then
        lbres@vpWidthF     = n_by_cols_1_first_vpWidthF
        offset = -0.01
    else if (cols .eq. 2) then
        lbres@vpWidthF     = n_by_cols_2_second_vpXF - n_by_cols_2_first_vpXF + n_by_cols_2_first_vpWidthF
        offset = offset + n_by_cols_2_lb_offset
    else if(cols .eq. 3) then
        lbres@vpWidthF     = n_by_cols_3_third_vpXF + n_by_cols_3_third_vpWidthF - n_by_cols_3_first_vpXF
        offset = offset + n_by_cols_3_lb_offset
    else
        error_msg("fatal", "SAMonsoon_panel_positions.ncl", \
                  "n_by_cols_labelbar", "Invalid value for cols=" + cols)
    end if
    end if
    end if

    ;; Update defaults with plottype specific labelbar resources
    if (plottype_lbres) then
        update_res = getvaratts(plottype_lbres)
        do idx_res = 0, dimsizes(update_res) - 1
            curr_res = update_res(idx_res)
            if(isatt(lbres, curr_res)) then
                delete(lbres@$curr_res$)
                lbres@$curr_res$ = plottype_lbres@$curr_res$
            end if
        end do
    end if

    nboxes = dimsizes(cnLevels)
    lblabels = sprintf("%5.2f", cnLevels)
    lbres@gsnFrame = False
    lbid = gsn_create_labelbar(wks, nboxes, lblabels, lbres)
    if (cols .eq. 1) then
        setvalues lbid
            "lbAutoManage" : False
            "lbLabelFontHeightF" : 0.01
        end setvalues
    end if
    getvalues lbid
        "lbLabelFontHeightF" : lblabelfh
        "lbTitleFontHeightF" : lbtitlefh
    end getvalues
    if (lblabelfh .gt. lbtitlefh) then
        setvalues lbid
            "lbAutoManage" : False
            "lbLabelFontHeightF" : lbtitlefh
        end setvalues
    end if
    amres = True
    amres@amZone           = 2
    amres@amSide           = "Bottom"
    amres@amParallelPosF   = 0.5            ; Center labelbar.
    amres@amOrthogonalPosF = 0.0 + offset   ; Move down, away from plot
    annoid = gsn_add_annotation(plot, lbid, amres)

    info_output(">>>>>>>> Leaving n_by_cols_labelbar.ncl", verbosity, 6)

end

undef("add_blank_plot_title")
function add_blank_plot_title(wks [1] : graphic,
                              title [1] : string,
                              rows [1] : float,
                              cols [1] : integer)
local adjust_blank_plot_title_position,  blank_plot,                bres,\
    bvpHeightF,                        bvpWidthF,                 bvpXF,\
    bvpYF,                             extra_blank_plot_spacing,  rows_int,\
    text,                              txres,                     verbosity
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering add_blank_plot_title.ncl", verbosity, 6)
    rows_int = tointeger(ceil(rows))
    extra_blank_plot_spacing = 0.095
    adjust_blank_plot_title_position = 0.03
    if (cols .eq. 1) then
        ;; Blank plot settings
        bvpXF = n_by_cols_1_first_vpXF
        bvpYF = n_by_cols_1_first_vpYF
        bvpWidthF = n_by_cols_1_first_vpWidthF
        if (rows_int .eq. 1) then
            bvpHeightF = n_by_cols_1_first_vpYF\
                       - n_by_cols_1_second_vpYF\
                       + extra_blank_plot_spacing / 3.0
        else if (rows_int .eq. 2) then
            bvpHeightF = n_by_cols_1_first_vpYF\
                       - n_by_cols_1_second_vpYF\
                       + n_by_cols_1_first_vpHeightF\
                       + extra_blank_plot_spacing
        else if (rows_int .eq. 3) then
            bvpHeightF = n_by_cols_1_first_vpYF\
                       - n_by_cols_1_third_vpYF\
                       + n_by_cols_1_first_vpHeightF\
                       + extra_blank_plot_spacing
        else if (rows_int .eq. 4) then
            bvpHeightF = n_by_cols_1_first_vpYF\
                       - n_by_cols_1_fourth_vpYF\
                       + n_by_cols_1_first_vpHeightF\
                       + extra_blank_plot_spacing
        else
            error_msg("fatal", "SAMonsoon_panel_positions.ncl", \
                      "add_blank_plot_title", \
                      "Invalid value for rows_int=" + rows_int)
        end if
        end if
        end if
        end if
    else if (cols .eq. 2) then
        ;; Blank plot settings
        bvpXF = n_by_cols_2_first_vpXF
        bvpYF = n_by_cols_2_first_vpYF
        bvpWidthF = n_by_cols_2_second_vpXF\
                  + n_by_cols_2_second_vpWidthF\
                  - n_by_cols_2_first_vpXF
        if (rows_int .eq. 2) then
            bvpHeightF = n_by_cols_2_first_vpYF\
                  - n_by_cols_2_third_vpYF\
                  + n_by_cols_2_first_vpHeightF\
                  + extra_blank_plot_spacing
        else if (rows_int .eq. 3) then
            bvpHeightF = n_by_cols_2_first_vpYF\
                       - n_by_cols_2_fifth_vpYF\
                       + n_by_cols_2_first_vpHeightF\
                       + extra_blank_plot_spacing
        else if (rows_int .eq. 4) then
            bvpHeightF = n_by_cols_2_first_vpYF\
                       - n_by_cols_2_seventh_vpYF\
                       + n_by_cols_2_first_vpHeightF\
                       + extra_blank_plot_spacing
        else if (rows_int .eq. 5) then
            bvpHeightF = n_by_cols_2_first_vpYF\
                       - n_by_cols_2_ninth_vpYF\
                       + n_by_cols_2_first_vpHeightF\
                       + extra_blank_plot_spacing
        else
            error_msg("fatal", "SAMonsoon_panel_positions.ncl", \
                      "add_blank_plot_title", \
                      "Invalid value for rows_int=" + rows_int)
        end if
        end if
        end if
        end if
    else if (cols .eq. 3) then
        bvpXF = n_by_cols_3_first_vpXF
        bvpYF = n_by_cols_3_first_vpYF
        bvpWidthF = n_by_cols_3_third_vpXF\
                  + n_by_cols_3_third_vpWidthF\
                  - n_by_cols_3_first_vpXF
        if (rows_int .eq. 3) then
            bvpHeightF = n_by_cols_3_first_vpYF\
                       - n_by_cols_3_seventh_vpYF\
                       + n_by_cols_3_first_vpHeightF + extra_blank_plot_spacing
        else if (rows_int .eq. 4) then
            bvpHeightF = n_by_cols_3_first_vpYF\
                       - n_by_cols_3_tenth_vpYF\
                       + n_by_cols_3_first_vpHeightF + extra_blank_plot_spacing
        else
            error_msg("fatal", "SAMonsoon_panel_positions.ncl", \
                      "add_blank_plot_title", \
                      "Invalid value for rows_int=" + rows_int)
        end if
        end if
    else
        error_msg("fatal", "SAMonsoon_panel_positions.ncl", \
                  "add_blank_plot_title", \
                  "Invalid value for rows_int=" + rows_int)
    end if
    end if
    end if

    bres           = True
    bres@vpXF      = bvpXF  ; Global variables set in the functions above
    bres@vpYF      = bvpYF + extra_blank_plot_spacing - adjust_blank_plot_title_position
    bres@vpWidthF  = bvpWidthF
    bres@vpHeightF = bvpHeightF
    bres@gsnDraw   = False
    bres@gsnFrame  = False

    ;; Turn off all tickmarks, labels, and borders
    bres@tmXBOn       = False
    bres@tmYLOn       = False
    bres@tmXTOn       = False
    bres@tmYROn       = False
    bres@tmXBBorderOn = False
    bres@tmXTBorderOn = False
    bres@tmYLBorderOn = False
    bres@tmYRBorderOn = False

    blank_plot = gsn_blank_plot(wks, bres)

    ;; Resources for adding text
    txres               = True
    txres@txFontHeightF = 0.02
    txres@txJust        = "CenterCenter"

    text = gsn_add_text(wks, blank_plot, title, 0.5, 0.99, txres)
    info_output(">>>>>>>> Leaving add_blank_plot_title.ncl", verbosity, 6)
    return blank_plot
end

undef("inset_labelbar")
procedure inset_labelbar(wks, plot, plot_res, title, plottype_lbres)
local verbosity, contour, lbres, nboxes, vpHeight
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering inset_labelbar.ncl", verbosity, 6)

    lbres = True
    contour = plot@contour
    getvalues contour
        "cnLevels" : cnLevels
        "cnFillColors" : cnFillColors
    end getvalues
    getvalues plot
        "vpHeightF" : vpHeightF
    end getvalues
    ;; Default labelbar resource settings
    lbres                   = True
    lbres@lbAutoManage      = True
    lbres@lbFillColors      = cnFillColors
    lbres@lbMonoFillPattern = True
    lbres@lbOrientation     = "Vertical"
    lbres@lbPerimOn         = False
    lbres@lbTitleOffsetF    = 0.0
    lbres@lbTitleExtentF    = 0.0
    lbres@gsnDraw           = True
    lbres@gsnMaximize       = True
    lbres@lbLeftMarginF     = 0.0
    lbres@lbTopMarginF      = 0.0
    lbres@lbLabelPosition   = "Right"
    lbres@lbBoxMinorExtentF = 0.2
    lbres@lbTitleOn         = False
    lbres@lbLabelAlignment  = "InteriorEdges"
    nboxes = dimsizes(cnLevels)
    lblabels = sprintf("%5.2f", cnLevels)
    lbres@vpHeightF = vpHeightF
    lbres@vpWidthF = 0.05
    lbres@lbPerimOn = True
    lbres@lbPerimFill = "SolidFill"
    lbres@lbPerimFillColor = "Background"

    lbres@gsnDraw = True
    lbres@gsnFrame = False
    lbid = gsn_create_labelbar(wks, nboxes, lblabels, lbres)

    amres = True
    amres@amZone           = 2
    amres@amSide           = "Top"
    amres@amParallelPosF   = lbres@vpWidthF + 0.01    ; Center labelbar.
    amres@amOrthogonalPosF = -1.0            ; Move down, away from plot
    annoid = gsn_add_annotation(plot, lbid, amres)
    info_output(">>>>>>>> Leaving inset_labelbar.ncl", verbosity, 6)
end
