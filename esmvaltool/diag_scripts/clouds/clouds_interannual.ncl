; CLOUDS_INTERANNUAL
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; PROJECT-NAME EMBRACE
; ############################################################################
; Description
;   Calculates the interannual variability estimated as the temporal standard
;   deviation calculated from monthly mean anomalies after subtracting the
;   climatological mean seasonal cycle.
;
; Required diag_script_info attributes (diagnostic specific)
;   none
;
; Optional diag_script_info attributes (diagnostic specific)
;   colormap:        e.g., WhiteBlueGreenYellowRed, rainbow
;   explicit_cn_levels: use these contour levels for plotting
;   extrafiles:      write plots for individual models to separate files
;                    (True, False)
;   projection:      map projection, e.g., Mollweide, Mercator
;
; Required variable_info attributes (variable specific)
;   none
;
; Optional variable_info attributes (variable specific)
;   long_name: description of variable
;   reference_dataset: name of reference datatset
;
; Caveats
;   none
;
; Modification history
;   20190220-lauer_axel: added provenance to output (v2.0)
;   20181120-lauer_axel: adapted code to multi-variable capable framework
;   20180923-lauer_axel: added writing of results to netcdf
;   20180611-lauer_axel: code rewritten for ESMValTool v2.0
;   20170620-lauer_axel: added tags for reporting
;   20160901-lauer_axel: added regridding option 1 deg x 1 deg
;   20151027-lauer_axel: moved call to 'write_references' to the beginning
;                        of the code
;   20150415-lauer_axel: written.
;
; ############################################################################

; #####################################
; # load external NCL code, if needed #
; #####################################

; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/contour_maps.ncl"

begin
  enter_msg(DIAG_SCRIPT, "")

  var0 = variable_info[0]@short_name
  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info0)
  if (isatt(variable_info[0], "reference_dataset")) then
    refname = variable_info[0]@reference_dataset
  end if
  names = metadata_att_as_array(info0, "dataset")
  infiles = metadata_att_as_array(info0, "filename")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  set_default_att(diag_script_info, "colormap", "WhiteBlueGreenYellowRed")
  set_default_att(diag_script_info, "extrafiles", False)
  set_default_att(diag_script_info, "projection", "CylindricalEquidistant")

  extrafiles = diag_script_info@extrafiles

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

  if (diag_script_info@write_plots.eq."True") then
    write_plots = True
  else
    write_plots = False
  end if

  ; get multi-model mean index (if present)

  mm_ind = ind(names .eq. "MultiModelMean")

  if (ismissing(mm_ind)) then
    mm_ind = -1
  end if

  ref_ind = -1  ; set to invalid value

  ; if reference dataset has been defined, use it so plots can be sorted
  if (isvar("refname")) then
    ref_ind = ind(names .eq. refname)
  end if

end

begin
  ind_all_sorted = ispan(0, dim_MOD - 1, 1)  ; create array

  if (ref_ind .ge. 0) then
    ind_wo_ref = ind(names .ne. refname)
    ind_all_sorted(0) = ref_ind
    ind_all_sorted(1:dim_MOD - 1) = ind_wo_ref
  end if

  maps = new(dim_MOD, graphic)  ; collect individual maps in a graphic array

  ; ###########################################
  ; # get data and average time               #
  ; ###########################################

  do ii = 0, dim_MOD - 1

    imod = ind_all_sorted(ii)

    if (isvar("data1")) then
      delete(data1)
    end if

    log_info("processing " + names(imod))

    if (isvar("A0")) then
      delete(A0)
    end if

    A0 = read_data(info0[imod])

    ; check dimensions

    dims = getvardims(A0)
    if (dimsizes(dims) .lt. 2) then
      error_msg("f", DIAG_SCRIPT, "", dimsizes(dims) + \
                " dimensions, need 2 or 3")
    end if
    idx = ind(dims .eq. "lat")
    if (ismissing(idx)) then
      error_msg("f", DIAG_SCRIPT, "", "no lat dimension")
    end if
    idx = ind(dims .eq. "lon")
      if (ismissing(idx)) then
      error_msg("f", DIAG_SCRIPT, "", "no lon dimension")
    end if

    if (var0.eq."pr") then
      ; convert from kg m-2 s-1 to mm day-1
      A0 = A0 * 86400.0
      A0@units = "mm day-1"
    end if

    ; subtract climatological seasonal cycle from time series

    if (isvar("timeseries")) then
      delete(timeseries)
    end if

    timeseries = calcMonAnomTLL(A0, time_operations(A0, -1, -1, \
                                "average", "monthlyclim", True))

    ; calculate temporal standard deviation for each grid cell

    data1 = dim_stddev_n_Wrap(timeseries, 0)

    ; normalize standard deviation and convert to percent

    if (isvar("mean")) then
      delete(mean)
    end if
    mean = time_operations(A0, -1, -1, "average", "annualclim", True)
    ; replace "epsilon" values with missing value
    mean = where(abs(mean).lt.1.0e-4, mean@_FillValue, mean)
    data1 = 100.0 * data1 / abs(mean)

    ; ###########################################
    ; # Style dependent annotation              #
    ; ###########################################
    ; retrieve unique strings describing the data
    ; function in ./diag_scripts/lib/ncl/style.ncl

    annots = project_style(info0, diag_script_info, "annots")

    ; ###########################################
    ; # plot ressources                         #
    ; ###########################################

    data1@res_gsnMaximize    = True      ; use full page for the plot
    data1@res_cnFillOn       = True      ; color plot desired
    data1@res_cnLineLabelsOn = False     ; contour lines
    data1@res_cnLinesOn      = False

    ; colors
    ; http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml

    if (isdefined("pal")) then
      delete(pal)
    end if
    pal = read_colormap_file(diag_script_info@colormap)

    ; annotation

    data1@res_tiMainOn             = False
    data1@res_gsnLeftStringFontHeightF = 0.015
    data1@res_cnLevelSelectionMode = "ExplicitLevels"

    if (diag_script_info@projection.eq."Robinson") then
      data1@res_mpPerimOn = False  ; turn off perimeter around map
      data1@res_mpGridLineColor = -1
      data1@res_mpGridAndLimbOn = True
    end if

    data1@res_mpOutlineOn      = True
    data1@res_mpFillOn         = False

    ; variable specific plotting settings

    if (any((/"clt"/).eq.var0)) then
      data1@res_cnLevels            = ispan(5, 50, 5)
    else
      data1@res_cnLevels            = ispan(5, 100, 5)
    end if

    if (var0.eq."lwp") then
      data1@res_mpOutlineOn        = False
      data1@res_mpFillOn           = True
      data1@res_mpLandFillColor    = "Black"
;      delete(pal)
;      pal = read_colormap_file("$diag_scripts/shared/plot/rgb/qcm3.rgb")
    end if

    nboxes = dimsizes(data1@res_cnLevels)
    clen = dimsizes(pal)
    stride = max((/1, ((clen(0)-1) - 2) / nboxes /))
    fill_colors = ispan(2, clen(0) - 1, stride)
    data1@res_cnFillColors = fill_colors

    data1@res_lbLabelBarOn         = False
    data1@res_gsnRightString       = ""

    ; map attributes

    data1@res_mpFillDrawOrder         = "PostDraw"    ; draw map last
    data1@res_cnMissingValFillColor   = "Gray"

    ; no tickmarks and no labels

    data1@res_tmYLLabelsOn         = False
    data1@res_tmYLOn               = False
    data1@res_tmYRLabelsOn         = False
    data1@res_tmYROn               = False
    data1@res_tmXBLabelsOn         = False
    data1@res_tmXBOn               = False
    data1@res_tmXTLabelsOn         = False
    data1@res_tmXTOn               = False
    data1@res_cnInfoLabelOn        = False    ; turn off cn info label

    data1@res_mpProjection         = diag_script_info@projection

    ; set explicit contour levels

    if (isatt(diag_script_info, "explicit_cn_levels")) then
      data1@res_cnLevelSelectionMode = "ExplicitLevels"
      data1@res_cnLevels = diag_script_info@explicit_cn_levels
    end if

    ; ###########################################
    ; # other Metadata: diag_script, var        #
    ; ###########################################
    ; add to data1, as attributes without prefix

    if (isatt(data1, "diag_script")) then  ; add to existing entries
      temp = data1@diag_script
      delete(data1@diag_script)
      data1@diag_script = array_append_record(temp, (/DIAG_SCRIPT/), 0)
      delete(temp)
    else  ; add as new attribute
      data1@diag_script = (/DIAG_SCRIPT/)
    end if
    data1@var = var0  ; Overwrite existing entry
    if (isatt(variable_info[0], "long_name")) then
      data1@var_long_name = variable_info[0]@long_name
    else
      data1@var_long_name = var0
    end if
    data1@var_units    = "%"

    ; copy attributes for netCDF output

    data1@long_name = "interannual variability " + data1@var_long_name
    data1@units = data1@var_units

    ; ###########################################
    ; # create the plot                         #
    ; ###########################################

    data1@res_gsnFrame       = False  ; don't advance frame
    data1@res_gsnDraw        = False

    ; function in aux_plotting.ncl

    if (ii.eq.0) then
      if (.not.extrafiles) then
        wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_interannual_" \
                      + var0)
      end if
;      drawNDCGrid(wks)  ; debugging option
    end if

    if (extrafiles) then
      if (isvar("wks")) then
        delete(wks)
      end if
      wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_interannual_" \
                    + var0 + "_" + annots(imod))
    end if

    maps(ii) = contour_map(wks, data1, var0)

    if (extrafiles) then
      if (write_plots) then  ; add labels
        txres = True
        txres@txFontHeightF = 0.03
        txres@txJust = "BottomRight"
        txres@txPerimOn = True
        txres@txBackgroundFillColor = "white"
        text = gsn_add_text(wks, maps(ii), annots(imod), 170, -80, txres)
        draw(maps(ii))
        frame(wks)
        plotfile = maps@outfile
      else
        plotfile = ""
      end if

      ; ##########################################
      ; # output each dataset to separate netCDF #
      ; ##########################################

      nc_filename = work_dir + "clouds_interannual_" + var0 + "_" \
        + annots(imod) + ".nc"
      nc_outfile = ncdf_write(data1, nc_filename)

      ; -------------------------------------------------------------------
      ; write provenance info
      ; -------------------------------------------------------------------

      statistics = (/"clim", "var"/)
      domain = ("glob")
      plottype = ("geo")
      climofile = infiles(imod)
      caption = "Interannual variability of variable " + var0 + \
                " from dataset " + annots(imod) + "."

      log_provenance(nc_outfile, plotfile, caption, statistics, domain, \
                     plottype, "", "", climofile)

    else  ; extrafiles .eq. false

      ; #########################################
      ; # output all datasets to common netCDF  #
      ; #########################################

      nc_filename = work_dir + "clouds_interannual_" + var0 + ".nc"
      nc_filename@existing = "append"
      data1@var = var0 + "_var_" + annots(imod)
      nc_outfile = ncdf_write(data1, nc_filename)

    end if  ; if extrafiles
  end do  ; ii-loop (datasets)

  if (write_plots) then
    pres                            = True    ; needed to override
                                              ; panelling defaults
    pres@gsnPanelLabelBar           = True    ; add common colorbar
    ; print dataset name on each panel
    pres@gsnPanelFigureStrings      = annots(ind_all_sorted)
    pres@gsnPanelFigureStringsFontHeightF = 0.007
    pres@lbLabelFontHeightF         = 0.01
    pres@lbAutoManage               = False
    pres@lbTopMarginF               = 0.1
    pres@lbTitleOn                  = True
    pres@lbTitleFontHeightF         = 0.009
    pres@lbTitlePosition            = "Bottom"
    pres@lbTitleString              = "~F8~s~F21~" + var0 + " (%)"
    pres@lbPerimOn                  = False     ; draw line around label
                                                ; bar area
    pres@gsnPanelCenter             = False
    pres@pmLabelBarOrthogonalPosF   = -0.01     ; shift label bar a bit to
                                                ; the bottom
;    pres@gsnPanelDebug              = True

    if (extrafiles) then
      if (isvar("wks")) then
        delete(wks)
      end if
      ; plot legend
      wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_interannual_" \
                    + var0 + "_legend")
      pres@lbMonoFillPattern = True
      pres@lbOrientation = "Horizontal"
      pres@vpWidthF = 0.7
      pres@vpHeightF = 0.1
      pres@lbLabelFontHeightF = 0.015
      pres@lbLabelAlignment = "InteriorEdges"
      pres@lbTitleFontHeightF = 0.015
      labels = tostring(data1@res_cnLevels)
      pres@lbFillColors = fill_colors
      gsn_labelbar_ndc(wks, nboxes, labels, 0.1, 0.9, pres)
    else
      outfile = panelling(wks, maps, (dim_MOD + 3) / 4, 4, pres)
      log_info(" Wrote " + outfile)
    end if
  else
    outfile = ""
  end if  ; if write_plots

  ; ------------------------------------------------------------------------
  ; write provenance to common netcdf and plot file
  ; ------------------------------------------------------------------------

  if (.not. extrafiles) then
    statistics = (/"clim", "var"/)
    domain = "global"
    plottype = "geo"
    caption = "Interannual variability of variable " + var0 + "."
    log_provenance(nc_outfile, outfile, caption, statistics, domain, \
                   plottype, "", "", infiles)
  end if

  leave_msg(DIAG_SCRIPT, "")

end
