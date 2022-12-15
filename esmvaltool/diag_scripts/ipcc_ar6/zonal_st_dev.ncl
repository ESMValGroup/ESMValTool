; #############################################################################
;  PLOTS ZONAL MEAN STANDARD DEVIATION
;  Author: Lisa Bock (lisa.bock@dlr.de)
;  IPCC AR6
; #############################################################################
;
;  Description
;     Global climate variability of surface temperature as represented by the
;     standard deviation of zonal-mean values
;     This script has been developed to reproduce Fig. 9.33 in
;     IPCC AR5 Chapter 9.
;
;  Optional variable_info attributes (variable specific)
;     reference_dataset: reference dataset; REQUIRED when calculating
;                        anomalies
;
;  Required diag_script_info attributes (diagnostics specific)
;     styleset: as in diag_scripts/shared/plot/style.ncl functions
;
;  Optional diag_script_info attributes (diagnostic specific)
;     plot_legend: if true, plot legend will be plotted
;     plot_units: variable unit for plotting
;     multi_model_mean: if true, multi-model mean and uncertaintiy will be
;                       plotted
;
;  Modification history
;  20190304-bock_lisa: written
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/ensemble.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/xy_line.ncl"


begin

  enter_msg(DIAG_SCRIPT, "")

  ; Check for required settings
  exit_if_missing_atts(diag_script_info, "styleset")

  var0 = "tas"
  info_items = select_metadata_by_name(input_file_info, var0)
  datasetnames = metadata_att_as_array(info_items, "dataset")
  dim_MOD = ListCount(info_items)
  dim_VAR = ListCount(variable_info)

  if (dim_VAR .gt. 1) then
    var1 = "tasa"
    info_items1 = select_metadata_by_name(input_file_info, var1)
    datasetnames1 = metadata_att_as_array(info_items1, "dataset")
    dim_MOD1 = ListCount(info_items1)
    do i = 0, dim_MOD1-1
      ListAppend(info_items, info_items1[i])
    end do
    datasetnames := array_append_record(datasetnames, datasetnames1, 0)
    dim_MOD = dim_MOD + dim_MOD1
  end if

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "plot_legend", True)
  set_default_att(diag_script_info, "multi_model_mean", True)

  ; Create output plot directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Set index of the reference (and alternative) dataset
  if (isatt(variable_info[0], "reference_dataset")) then
    ref_ind = ind(datasetnames.eq.variable_info[0]@reference_dataset)
    A_ref = read_data(info_items[ref_ind])
  else
    A_ref = read_data(info_items[0])
  end if

  ; Determine time range
  start_year = min(metadata_att_as_array(info_items, "start_year"))
  end_year = max(metadata_att_as_array(info_items, "end_year"))
  all_years = ispan(start_year, end_year, 1)

  ; ------------------------------------------------------------------
  ; Calculation of standard deviation
  ; ------------------------------------------------------------------

  ; Create model array
  model_arr = new((/dim_MOD, dimsizes(A_ref&lat)/), float)
  model_arr!0 = "model"
  model_arr!1 = "lat"
  model_arr&model = datasetnames
  model_arr&lat = A_ref&lat

  ; Loop over models
  do imod = 0, dim_MOD - 1

    ; Read data
    A0 = read_data(info_items[imod])

    ; Convert units for plotting (if required)
    if (isatt(diag_script_info, "plot_units")) then
      if (info_items[imod]@short_name .ne. "tasa" .and. \
          info_items[imod]@units .ne. diag_script_info@plot_units) then
        A0 = convert_units(A0, diag_script_info@plot_units)
      end if
    end if

    ; detrending all timeseries
    A0_detrend = dtrend_quadratic_msg_n(A0, False, False, 0)

    A0 = (/A0_detrend/)

    procmod = dim_avg_Wrap(A0)
    delete(A0)
    delete(A0_detrend)

    A0_timavg = time_operations(procmod, -1, -1, "average", "yearly", True)
    delete(procmod)

    ; calculate zonal standard deviation
    A0_stddev = dim_stddev_n_Wrap(A0_timavg, 0)
    delete(A0_timavg)
    if(any(A0_stddev.eq."-nan")) then
      A0_stddev = where(A0_stddev.eq."-nan", A0_stddev@_FillValue, A0_stddev)
    end if

    if (dimsizes(A0_stddev&lat).ne.dimsizes(A_ref&lat)) then
      dnames = getVarDimNames(A0_stddev)
      dx = ind(dnames.eq."lat")
      var_hreg = linint1_n_Wrap(A0_stddev&lat, A0_stddev, False, \
                                A_ref&lat, 0, dx)
    else
      var_hreg = A0_stddev
      copy_VarAtts(A0_stddev, var_hreg)
    end if

    model_arr(imod, :) = (/var_hreg/)
    if (imod.eq.0) then
        copy_VarAtts(var_hreg, model_arr)
    end if
    delete(A0_stddev)

  end do

  ; ------------------------------------------------------------------
  ; Calculate multi model mean and stddev
  ; ------------------------------------------------------------------

  model_arr_stddev = new((/4, dimsizes(A_ref&lat)/), float)
  model_arr_stddev!1 = "lat"
  model_arr_stddev&lat = A_ref&lat

  model_arr_stddev(0, :) = dim_avg_n_Wrap(model_arr, 0)
  model_arr_stddev(1, :) = dim_stddev_n_Wrap(model_arr, 0)
  model_arr_stddev(2, :) = model_arr_stddev(0, :) - model_arr_stddev(1, :)
  model_arr_stddev(3, :) = model_arr_stddev(0, :) + model_arr_stddev(1, :)

  ; ------------------------------------------------------------------
  ; Output to NetCDF
  ; ------------------------------------------------------------------

  out_path = config_user_info@work_dir
  system("mkdir -p " + out_path)
  out_path = out_path + var0 + "_std_dev_zonmean.nc"
  model_arr@ncdf = out_path
  model_arr@diag_script = DIAG_SCRIPT
  model_arr@var = var0
  ncdf_outfile = ncdf_write(model_arr, out_path)

  ; ------------------------------------------------------------------
  ; Plotting
  ; ------------------------------------------------------------------

  ; Define workstation
  outfile = config_user_info@plot_dir + var0 + "_std_dev_zonmean"
  wks = gsn_open_wks(file_type, outfile)

  ; Select colors and other plotting attributes
  ; (see ./diag_scripts/shared/plot/style.ncl)
  colors := project_style(info_items, diag_script_info, "colors")
  dashes := project_style(info_items, diag_script_info, "dashes")
  thicks := project_style(info_items, diag_script_info, "thicks")
  annots := project_style(info_items, diag_script_info, "annots")

  ; Select colors and other plotting attributes for multi-model mean
  if (diag_script_info@multi_model_mean) then
    items_mmm = NewList("fifo")
    mmm_var = True
    mmm_var@dataset = "MultiModelMean"
    ListPush(items_mmm, mmm_var)
    colors_mmm = project_style(items_mmm, diag_script_info, "colors")
    dashes_mmm = project_style(items_mmm, diag_script_info, "dashes")
    thicks_mmm = project_style(items_mmm, diag_script_info, "thicks")
    annots_mmm = project_style(items_mmm, diag_script_info, "annots")
  end if

  ; Set/copy resources
  res                   = True
  res@xyLineColors      = colors  ; change line color
  res@xyLineThicknesses = thicks  ; make 2nd lines thicker
  res@xyMonoDashPattern = False
  res@xyDashPatterns    = dashes
  res@gsnDraw           = False   ; don't draw yet
  res@gsnFrame          = False   ; don't advance frame yet
  res@tiMainFontHeightF = .02    ; font height
  res@tiXAxisFontHeightF = .015
  res@tiYAxisFontHeightF = .015
  res@gsnMaximize       = True
  res@vpHeightF = 0.4
  res@vpWidthF = 0.8

  res@tmYROn       = False
  res@tmXTOn       = False
  res@tmYRBorderOn = False
  res@tmXTBorderOn = False

  res@tmXBMode = "Explicit"
  d2r  = get_d2r("f")
  res@tmXBValues = sin(d2r * ispan(-60, 60, 30))
  res@tmXBLabels = (/"60~F34~0~F~S", "30~F34~0~F~S", "0~F34~0~F~", \
                     "30~F34~0~F~N", "60~F34~0~F~N"/)
  res@tmXBMinorValues = sin(d2r * ispan(-90, 90, 10))
  res@trXMinF = -1.
  res@trXMaxF = 1.

  res@trYMinF  = min(model_arr) - 0.05 * (max(model_arr) - min(model_arr))
  res@trYMaxF  = max(model_arr) + 0.05 * (max(model_arr) - min(model_arr))

  res@tiMainString = "Temporal variability of near-surface air temperature"
  res@tiXAxisString = "Latitude"

  if (model_arr@units .eq. "degC") then
    units = "~F34~0~F~ C"
  else
    units = model_arr@units
  end if
  res@tiYAxisString = "Standard deviation of temperature (" + units + ")"

  model_arr&lat = sin(d2r * model_arr&lat)

  ; Plot
  plot = gsn_csm_xy(wks, model_arr&lat, model_arr, res)

  ; Add multi model mean and stddev
  if (diag_script_info@multi_model_mean) then
    ; Stddev
    res_stddev = True
    copy_VarMeta(res, res_stddev)
    res_stddev@gsnXYFillColors = "LightGrey"
    delete(res_stddev@xyLineColors)
    res_stddev@xyLineColor = -1  ; Make lines transparent
    shading_plot = gsn_csm_xy(wks, model_arr&lat, model_arr_stddev(2:3, :), \
                              res_stddev)
    overlay(plot, shading_plot)
    ; MMM
    delete([/res@xyLineThicknesses, res@xyLineColors, res@xyDashPatterns/])
    res@xyLineThicknesses = thicks_mmm
    res@xyLineColors = colors_mmm
    res@xyDashPatterns = dashes_mmm
    mmm = gsn_csm_xy(wks, model_arr&lat, model_arr_stddev(0, :), res)
    overlay(plot, mmm)
  end if

  ; ***********************************************
  ; legend resources
  ; ***********************************************

  if (diag_script_info@plot_legend) then

    lgres                    = True
    lgres@lgItemType         = "Lines"   ; show lines only (default)
    lgres@lgLabelFontHeightF = 1.0       ; set the legend label font thickness
    lgres@vpWidthF           = 0.15      ; width of legend (NDC)
    lgres@vpHeightF          = 0.7       ; height of legend (NDC)
    lgres@lgPerimColor       = "gray"    ; draw the box perimeter in orange
    lgres@lgPerimThicknessF  = 1.0       ; thicken the box perimeter

    lgres@lgLineColors      = colors
    lgres@lgDashIndexes     = dashes
    lgres@lgLineThicknesses = thicks + 3
    labels = annots
    nitems = dimsizes(labels)
    lgres@lgItemOrder = ispan(nitems - 1, 0, 1)

    ; Create legend
    lbid = gsn_create_legend(wks, nitems, labels, lgres)

    amres = True
    amres@amParallelPosF   = 0.65
    amres@amOrthogonalPosF = 0.0
    annoid1 = gsn_add_annotation(plot, lbid, amres)

  else

    colors := project_style(info_items1, diag_script_info, "colors")
    dashes := project_style(info_items1, diag_script_info, "dashes")
    annots := project_style(info_items1, diag_script_info, "annots")
    colors := array_append_record(colors, (/"red"/), 0)
    dashes := array_append_record(dashes, (/0/), 0)
    annots := array_append_record(annots, (/"MultiModelMean"/), 0)

    lgres                    = True
    lgres@lgItemType         = "Lines"   ; show lines only (default)
    lgres@lgLabelFontHeightF = .06       ; set the legend label font thickness
    lgres@vpWidthF           = 0.2      ; width of legend (NDC)
    lgres@vpHeightF          = 0.15      ; height of legend (NDC)
    lgres@lgPerimColor       = "gray"    ; draw the box perimeter in orange
    lgres@lgPerimThicknessF  = 1.0       ; thicken the box perimeter
    lgres@lgPerimFill        = 0
    lgres@lgPerimFillColor   = "white"

    lgres@lgLineColors      = colors
    lgres@lgDashIndexes     = dashes
    lgres@lgLineThicknessF = 3

    ; Create legend
    lbid = gsn_create_legend(wks, dimsizes(annots), annots, lgres)

    amres = True
    amres@amParallelPosF   = -0.25
    amres@amOrthogonalPosF = - 0.25
    annoid1 = gsn_add_annotation(plot, lbid, amres)

  end if

  resP = True
  resP@gsnMaximize = True

  gsn_panel(wks, plot, (/1, 1/), resP)

  log_info(" wrote " + outfile + "." + file_type)

  ; ------------------------------------------------------------------------
  ; write provenance to netcdf output and plot file(s) (mean)
  ; ------------------------------------------------------------------------

  authors = (/"bock_lisa"/)
  statistics = (/"stddev"/)
  domain = ("global")
  plottype = ("line")
  references = (/"eyring21ipcc"/)
  prov_caption = "Temporal variability of Near-Surface Air Temperature"
  infiles = metadata_att_as_array(info_items, "filename")

  log_provenance(out_path, outfile + "." + file_type, prov_caption, \
                 statistics, domain, plottype, authors, references, infiles)

  leave_msg(DIAG_SCRIPT, "")

end
