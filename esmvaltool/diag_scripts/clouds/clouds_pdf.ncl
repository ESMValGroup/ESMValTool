; CLOUDS_PDF
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; PROJECT-NAME ESA-CCI CMUG
; ############################################################################
; Description
;   Calculates temporal/spatial pdf's of 2-d (cloud) parameters for selected
;   regions.
;
; Required diag_script_info attributes (diagnostic specific)
;   none
;
; Optional diag_script_info attributes (diagnostic specific)
;   filename_add:       optionally add this string to output filenames
;   region:             show only selected geographic region given as latmin,
;                       latmax, lonmin, lonmax
;   styleset:        "CMIP5", "DEFAULT"
;
; Required variable attributes (variable specific)
;   none
;
; Optional variable_info attributes (variable specific)
;   reference_dataset: reference dataset
;
; Caveats
;   none
;
; Modification history
;   20190916-lauer_axel: written.
;
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/xy_line.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  var0 = variable_info[0]@short_name
  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info0)
  if (isatt(variable_info[0], "reference_dataset")) then
    refname = variable_info[0]@reference_dataset
  end if
  names = metadata_att_as_array(info0, "dataset")
  projects = metadata_att_as_array(info0, "project")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "filename_add", "")
  set_default_att(diag_script_info, "styleset", "DEFAULT")
  set_default_att(diag_script_info, "multi_model_mean", False)

  if (diag_script_info@filename_add .ne. "") then
    filename_add = "_" + diag_script_info@filename_add
  else
    filename_add = ""
  end if

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

  ref_ind = -1  ; set to invalid value

  ; if attribute is present, use it so datasets can be sorted
  if (isvar("refname")) then
    ; set reference model
    ref_ind = ind(names .eq. refname)
    if (ismissing(ref_ind)) then
      log_info("warning: reference dataset (" + refname + ") not found.")
      ref_ind = -1
    end if
  end if

  climofiles = metadata_att_as_array(info0, "filename")

  outfile = ""

end

begin
  ; ############
  ; # get data #
  ; ############

  plot = new((/dim_MOD/), graphic)

  ind_all_sorted = ispan(0, dim_MOD - 1, 1)  ; create array

  if (ref_ind .ge. 0) then
    ind_wo_ref = ind(names .ne. refname)
    ind_all_sorted(0) = ref_ind
    ind_all_sorted(1:dim_MOD - 1) = ind_wo_ref
  end if

  bins = 20

  data_all = new((/dim_MOD, bins/), float)
  data_all!0 = "models"
  data_all&models = names

  ; filenames for netcdf output

  nc_filename = work_dir + "clouds_pdf_" + var0 + filename_add + ".nc"
  nc_filename@existing = "append"

  do ii = 0, dim_MOD - 1

    imod = ind_all_sorted(ii)
    log_info("processing " + names(imod))

    if (isvar("data1")) then
      delete(data1)
    end if

    data1 = read_data(info0[imod])

    ; check dimensions

    dims = getvardims(data1)
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

    ; if requested, select geographical region

    if (isatt(diag_script_info, "region")) then
      region = diag_script_info@region
      data1 := area_operations(data1, region(0), region(1), region(2), \
                               region(3), "extract", False)
      delete(region)
    end if

    ; variable specific plotting settings

    res = True

    if (var0.eq."pr") then
      ; convert from kg m-2 s-1 to mm day-1
      data1 = data1 * 86400.0
      data1@units = "mm day-1"
    end if

    if (var0.eq."lwp") then
    end if

    if (var0.eq."tas") then
      ; convert from K to degC
      data1 = data1 - 273.15
      data1@units = "degC"
    end if

    if (var0.eq."clt") then
      res@gsnHistogramClassIntervals = ispan(0, 100, 10)
    end if

    if (var0.eq."clivi") then
    end if

    if (var0.eq."clwvi") then
    end if

    if (var0.eq."swcre") then
    end if

    if (var0.eq."lwcre") then
    end if

    if (var0.eq."netcre") then
    end if

    ; ###########################################
    ; # other Metadata: diag_script, var        #
    ; ###########################################
    ; add to data1 as attributes without prefix

    if (isatt(data1, "diag_script")) then  ; add to existing entries
      temp = data1@diag_script
      delete(data1@diag_script)
      data1@diag_script = array_append_record(temp, (/DIAG_SCRIPT/), 0)
      delete(temp)
    else  ; add as new attribute
      data1@diag_script = (/DIAG_SCRIPT/)
    end if

    if (isatt(variable_info[0], "long_name")) then
      data1@var_long_name = variable_info[0]@long_name
    end if

    data1@var = var0

    if (isatt(variable_info[0], "units")) then
      data1@var_units = variable_info[0]@units
    else
      data1@var_units = ""
    end if

    ; ###########################################
    ; # create the plot                         #
    ; ###########################################

    ; function in aux_plotting.ncl

    if (ii.eq.0) then
      wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_pdf_" + var0 \
                    + filename_add)
      wks_line = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_pdf_" + var0 \
                         + filename_add + "_line")
    end if

;    res@gsnDraw                  = False
;    res@gsnFrame                 = False
    res@gsnHistogramComputePercentages = True
    res@gsnHistogramPercentSign = True
    res@tiXAxisString = data1@long_name + " (" + data1@units + ")"
    res@tmYLOn = False
    res@tiYAxisOn = False
    res@tiMainString             = names(imod)

    plot(ii) = gsn_histogram(wks, ndtooned(data1), res)

    ; mandatory netcdf output

    if (ii .eq. 0) then
      opt = False
    else
      opt = True
      opt@bin_min = bin_min
      opt@bin_max = bin_max
    end if

    data_out = pdfx(data1, bins, opt)
    data_out@diag_script = (/DIAG_SCRIPT/)
    data_out@var = var0 + "_pdf_" + names(imod)
    nc_outfile = ncdf_write(data_out, nc_filename)

    if (ii .eq. 0) then
      bin_min = data_out@bin_bound_min
      bin_max = data_out@bin_bound_max
      data_all!1 = var0
      data_all&$var0$ = data_out@bin_center
      data_all@units = data_out@units
    end if

    data_all(imod, :) = tofloat(data_out)

    delete(data_out)

    ; =======================================================================

  end do  ; ii-loop (models)

  ; sort plots if needed (observations go first)

  plottmp = ispan(0, dim_MOD - 1, 1)
  plotind = plottmp

  ; move plots of observational datasets (if present) into the first line(s)
  ; of the panel plot

  j = 0
  do i = 0, dimsizes(plottmp) - 1
    if (i.eq.ref_ind) then
      plotind(j) = plottmp(i)
      j = j + 1
    else if (plottmp(i) .lt. dimsizes(projects)) then
      if (isStrSubset(str_lower(projects(plottmp(i))), "obs")) then
        plotind(j) = plottmp(i)
        j = j + 1
      end if
    end if
    end if
  end do

  do i = 0, dimsizes(plottmp) - 1
    if ((isStrSubset(str_lower(projects(plottmp(i))), \
                     "obs")).or.(i.eq.ref_ind)) then
    else
      plotind(j) = plottmp(i)
      j = j + 1
    end if
  end do

  pres                      = True    ; needed to override
                                      ; panelling defaults

  pres@gsnPanelFigureStrings = names(plotind)

  pres@gsnPanelFigureStringsFontHeightF = min((/0.01, 0.01 * 6.0 \
                                              / tofloat((dim_MOD + 1) / 2)/))
  pres@lbLabelFontHeightF               = min((/0.015, 0.01 * 6.0 \
                                              / tofloat((dim_MOD + 1) / 2)/))
  pres@lbAutoManage                     = False
  pres@lbTopMarginF                     = 0.1
  pres@lbTitleOn                        = True
  pres@lbTitleFontHeightF               = min((/0.015, 0.01 * 6.0 \
                                              / tofloat((dim_MOD + 1) / 2)/))
  pres@lbTitlePosition                  = "Bottom"
  pres@lbTitleString                    = data1@long_name + " (" \
    + data1@units + ")"
  pres@lbPerimOn                        = False   ; draw line around label
                                                  ; bar area
  pres@gsnPanelCenter                   = False
  if (dim_MOD.le.8) then
    pres@pmLabelBarOrthogonalPosF       = -0.03
  else
    pres@pmLabelBarOrthogonalPosF       = -0.01   ; shift label bar a bit to
                                                  ; the bottom
  end if

  pres@gsnPanelRowSpec = True             ; tell panel what order to plt
  pres@gsnPanelYWhiteSpacePercent = 5
  pres@gsnPanelXWhiteSpacePercent = 5
  if (isatt(diag_script_info, "PanelTop")) then
    top = tofloat(diag_script_info@PanelTop)
  else
    top = 0.99  ; default
  end if
  pres@gsnPanelTop = top

  if (isvar("plotsperline")) then
    delete(plotsperline)
  end if

  n = toint((dim_MOD + 1) / 4)
  if (n .lt. 1) then
    n = 1
  end if
  plotsperline = new(n, integer)
  plotsperline = 4

  if ((isStrSubset(str_lower(projects(plotind(0))), \
                   "obs")).and. \
     .not.(isStrSubset(str_lower(projects(plotind(1))), \
                       "obs"))) then
    plotsperline(0) = 1
  end if

  if (sum(plotsperline).gt.dimsizes(plotind)) then
    plotsperline(dimsizes(plotsperline) - 1) = 1
  end if

  if (sum(plotsperline).lt.dimsizes(plotind)) then
    xadd = 1
    xtmp = array_append_record(plotsperline, xadd, 0)
    delete(plotsperline)
    plotsperline = xtmp
    delete(xtmp)
  end if

  gsn_panel(wks, plot(plotind), plotsperline, pres)
  outfile = wks@fullname

  ; ========================================================================

  delete(res)
  res = True
  res@tiXAxisString = data1@long_name + " (" + data1@units + ")"
  res@tiYAxisString = "Frequency (" + data_all@units + ")"
  diag_script_info@xy_line_legend = False
  data_all@legend_outside = True
  wks_line@legendfile = get_plot_dir() + "clouds_pdf_" + var0 \
    + filename_add + "_legend"
  xy_line(wks_line, data_all, data_all&$var0$, dim_stddev_n(data_all, 0), \
          res, info0)
  outfile_line = wks_line@fullname

  ; ========================================================================

  log_info("Wrote " + outfile)
  log_info("Wrote " + outfile_line)

  ; ------------------------------------------------------------------------
  ; write provenance to netcdf output and plot file(s) (mean)
  ; ------------------------------------------------------------------------

  statistics = (/"clim", "pdf"/)
  if (isatt(diag_script_info, "region")) then
    domain = "reg"
  else
    domain = "global"
  end if
  plottype = "geo"

  caption = "Pdf of variable " + var0 + ")."
  log_provenance(nc_outfile, outfile, caption, statistics, \
                 domain, plottype, "", "", climofiles)

  ; ========================================================================

  leave_msg(DIAG_SCRIPT, "")

end
