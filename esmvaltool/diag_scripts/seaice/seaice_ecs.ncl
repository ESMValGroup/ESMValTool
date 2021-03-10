; #############################################################################
; EMERGENT CONSTRAINTS FOR SEA ICE
; Author: Senftleben, Daniel (DLR, Germany)
; CRESCENDO project
; #############################################################################
; Description
;   - Produces a scatter plot following Massonnet et al., 2012:
;     historical Sept sea ice extent (SSIE) trend vs hist SSIE long-term mean
;   - Produces two scatter plots following IPCC AR5 Fig. 12.31 a), d):
;     - historical SSIE mean vs 1st year of disappearance (YOD) RCP8.5
;     - historical SSIE trend vs YOD RCP8.5
;
; Required diag_script_info attributes
;   hist_exp: name of historical experiment (string)
;   month: selected month (1, 2, ..., 12) or annual mean ("A")
;   rcp_exp: name of RCP experiment (string)
;   region: region to be analyzed ( "Arctic" or "Antarctic")
;
; Optional diag_script_info attributes
;   fill_pole_hole: fill observational hole at North pole (default: False)
;   styleset: color style (e.g. "CMIP5")
;
; Optional variable_info attributes (variable specific)
;   reference_dataset: reference dataset
;
; Caveats
;   - CMOR units expected: areacello in m2, sic in %
;   - only 1 reference dataset supported
;
; Modification history
;    20191028-lauer_axel: rewritten for v2.0
;    20161212-senftleben_daniel: implemented.
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

load "$diag_scripts/seaice/seaice_aux.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/scatterplot.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  ; Get metadata items
  att = True
  att@mip = "OImon"
  info = select_metadata_by_atts(input_file_info, att)  ; variable
  var0 = info[0]@short_name
  datasets = metadata_att_as_array(info, "dataset")
  ensembles = metadata_att_as_array(info, "ensemble")
  tmpexps = metadata_att_as_array(info, "exp")
  exps = where(.not.ismissing(tmpexps), tmpexps, "")
  delete(tmpexps)
  dim_MOD = dimsizes(datasets)
  att@mip = "fx"
  areainfo = select_metadata_by_atts(input_file_info, att)  ; area
  areadatasets = metadata_att_as_array(areainfo, "dataset")
  delete(att)

  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"hist_exp", "rcp_exp", "month", \
                                           "region"/))

  hist_exp = diag_script_info@hist_exp
  rcp_exp = diag_script_info@rcp_exp

  ; Check for reference dataset
  if (isatt(info[0], "reference_dataset")) then
    refname = info[0]@reference_dataset
  else
    refname = ""
  end if

  ref_ind = -1  ; set to invalid value

  ; if attribute is present, use it so correlations can be calculated
  if (refname.ne."") then
    ; set reference model
    ref_ind = ind(datasets .eq. refname)
    if (ismissing(ref_ind)) then
      log_info("warning: reference dataset (" + refname + ") not found.")
      ref_ind = -1
    end if
  end if

  ; Check fill_pole_hole option
  if (.not.isatt(diag_script_info, "fill_pole_hole")) then
    diag_script_info@fill_pole_hole = False
  else
    if (diag_script_info@fill_pole_hole) then
      if (diag_script_info@region.eq."Antarctic") then
        diag_script_info@fill_pole_hole = False
        error_msg("f", DIAG_SCRIPT, "", "No pole hole at Antarctic " + \
                  "- setting diag_script_info@fill_pole_hole to False")
      end if
    end if
  end if

  ; Create output plot directory
  if (diag_script_info@write_plots) then
    plot_dir = config_user_info@plot_dir
    system("mkdir -p " + plot_dir)
  end if

  ; Create output directory (work)
  system("mkdir -p " + config_user_info@work_dir)

  ; Plot file type
  file_type = config_user_info@output_file_type

  ; Determine time range
  start_years = metadata_att_as_array(info, "start_year")
  end_years   = metadata_att_as_array(info, "end_year")

  ; Get experiment indices

  ; historical experiments (with no reference dataset - if specified)
  idcs_hist_tmp = ind(exps.eq.hist_exp)
  itmp = ind(idcs_hist_tmp.ne.ref_ind)
  idcs_hist = idcs_hist_tmp(itmp)

  delete([/itmp, idcs_hist_tmp/])
  idcs_rcp_tmp = ind(exps.eq.rcp_exp)
  idcs_rcp = new(dimsizes(idcs_rcp_tmp), integer)

  if (dimsizes(idcs_hist).ne.dimsizes(idcs_rcp_tmp)) then
        error_msg("fatal", DIAG_SCRIPT, "", "need same number of " \
                  + hist_exp + " and " + rcp_exp + " simulations")
  end if

  ; find indices of RCP experiments corresponding to historical experiments
  ; (in case both lists of datasets are not in the same order)

  datasets_rcp = datasets(idcs_rcp_tmp)  ; names of RCP datasets
  datasets_hist = datasets(idcs_hist)  ; names of historical datasets
  do imod = 0, dimsizes(idcs_hist) - 1
    idcs_rcp(imod) = idcs_rcp_tmp(ind(datasets_rcp.eq.datasets_hist(imod)))
  end do

  delete([/datasets_rcp, datasets_hist, idcs_rcp_tmp/])

  if (any(ismissing(idcs_rcp))) then
        error_msg("fatal", DIAG_SCRIPT, "", "list of " \
                  + hist_exp + " simulations does not contain all datasets " \
                  + "specified in list of " + rcp_exp + " simulations")
  end if

  hist_ind = 0  ; counter for historical simulations
  rcp_ind  = 0  ; counter for RCP simulations

  ; Set up time arrays
  time_hist = ispan(max(toint(start_years(idcs_hist))), \
                    min(toint(end_years(idcs_hist))), 1)
  time_rcp  = ispan(max(toint(start_years(idcs_rcp))), \
                    min(toint(end_years(idcs_rcp))), 1)

  ; Initialize arrays collecting data
  SSIE_MEAN = new((/dimsizes(idcs_hist)/), float)  ; long-term mean
  SSIE_MEAN!0 = "datasets"
  SSIE_MEAN&datasets = datasets(idcs_hist) + "_" + ensembles(idcs_hist)
  SSIE_TRND = SSIE_MEAN  ; trend of September Arctic sea ice extent
  SSIE_TRND@units = "10~S~3~N~ km~S~2~N~ / 10 yr"
  YOD       = new((/dimsizes(idcs_rcp)/), float)  ; year of disappearance
  YOD!0     = "datasets"
  YOD&datasets = datasets(idcs_rcp) + "_" + ensembles(idcs_rcp)

  ; Loop over datasets
  do imod = 0, dim_MOD - 1

    log_info("Processing dataset " + datasets(imod))

    sic = read_data(info[imod])

    ; Check grid type
    l_irreg = True
    if (sic!1.eq."lat" .and. sic!2.eq."lon") then
      if (dimsizes(dimsizes(sic&lat)).eq.1 .and. \
          dimsizes(dimsizes(sic&lon)).eq.1) then
        l_irreg = False
      end if
    end if

    ; Look for area variable
    if (any(datasets(imod).eq.areadatasets)) then

      ; Read areacello
      log_info("Using variable areacello for the area")
      area = read_data(areainfo[imod])
      locinfo = areainfo[imod]
      locinfo@short_name = "lat"  ; read auxiliary latitude coordinate
      arealat = read_data(locinfo)
      delete(locinfo)

      ; Check grid consistency
      if (.not.isconform(sic(0, :, :), area)) then
        if (l_irreg) then
          error_msg("f", DIAG_SCRIPT, "", "areacello not conformal")
        else
          delete(area)
          delete(arealat)
          log_info("areacello has inconsistent dimensions, " + \
                   "calculating area for a regular grid instead")
          area = map_area(sic&lat, sic&lon)
          arealat = area&lat
        end if
      end if

    else

      if (.not.l_irreg) then
        log_info("Calculating area for a regular grid")
        area = map_area(sic&lat, sic&lon)
        arealat = area&lat
      else
        error_msg("f", DIAG_SCRIPT, "", "cannot calculate area")
      end if

    end if

    ; Calcluate sea ice extent
    SSIE  = sea_ice_area(sic, area, arealat, l_irreg, \
                         diag_script_info@region, \
                         start_years(imod), end_years(imod), \
                         diag_script_info@month, "extent", \
                         diag_script_info@fill_pole_hole)

    if ((exps(imod).eq.hist_exp) .or. (imod.eq.ref_ind)) then
      ; long-term mean of sea ice extent [10^6 km^2]
      mean = dim_avg_Wrap(SSIE)
      ; trend of sea ice extent [10^3 km^2 decade^-1]
      trend = regline(time_hist, SSIE * 1.e3) * 10.
      if (imod.eq.ref_ind) then
        ssie_mean_ref = mean
        ssie_trnd_ref = trend
      else
        SSIE_MEAN(hist_ind) = mean
        SSIE_TRND(hist_ind) = trend
        hist_ind = hist_ind + 1
      end if
    end if

    ; calculate first year of disappearance (YOD)
    ;  = 1st year when 5-yr smoothed SSIE drops below 10^6 km^2
    ;    for more than 5 consecutive years
    if (exps(imod).eq.rcp_exp) then
      YOD(rcp_ind) = calc_yod(SSIE, 1., 5)
      rcp_ind = rcp_ind + 1
    end if

    ; Clean up
    delete([/sic, area, arealat, SSIE/])
  end do  ; imod

  ; ===========================================================================
  ; plotting section
  ; ===========================================================================

  smonths = (/"January", "February", "March", "April", \
              "May", "June", "July", "August", \
              "September", "October", "November", "December"/)

  if (diag_script_info@month .eq. "A") then
    i2month = "annual"
  else
    i2month = smonths(toint(diag_script_info@month) - 1)
  end if

  ; -----------------------------------------------
  ; scatter plot historical MEAN vs historical TRND
  ; -----------------------------------------------

  outf = "SSIE-MEAN_vs_SSIE-TRND_" + var0 + "_extend_" \
         + diag_script_info@region + "_" + i2month + "_" + min(time_hist) \
         + "-" + max(time_hist)

  toplot = (/SSIE_MEAN, SSIE_TRND/)
  toplot!0 = "datasets"  ; Required by plot routine, becomes axes labels
  toplot&datasets = (/i2month + " " + diag_script_info@region \
                     + " SSIE averaged~C~over " \
                     + min(time_hist) + "-" + max(time_hist) \
                     + " [" + SSIE_MEAN@units + "]", \
                     "Trend in " + i2month + " " + diag_script_info@region \
                     + " SSIE~C~over " \
                     + min(time_hist) + "-" + max(time_hist) \
                     + " [" + SSIE_TRND@units + "]"/)
  caption1 = toplot&datasets(0) + " vs " + toplot&datasets(1) + "."
  wks_MEAN_TRND = get_wks(0, DIAG_SCRIPT, outf)

  min_mean = min((/min(SSIE_MEAN), ssie_mean_ref/))
  max_mean = max((/max(SSIE_MEAN), ssie_mean_ref/))
  min_trnd = min((/min(SSIE_TRND), ssie_trnd_ref/))
  max_trnd = max((/max(SSIE_TRND), ssie_trnd_ref/))

  toplot@res_tiMainString = diag_script_info@region + " Sea Ice Extent (SSIE)"
  toplot@res_gsnCenterString = ""
  toplot@res_trXMinF  = decimalPlaces((min_mean - 1), 0, False)
  toplot@res_trXMaxF  = decimalPlaces((max_mean + 1), 0, False)
  toplot@res_trYMinF  = decimalPlaces((min_trnd - 100), 0, False)
  toplot@res_trYMaxF  = decimalPlaces((max_trnd + 100), 0, False)
  toplot@res_xyMarker = 16
  colors = project_style(info, diag_script_info, "colors")
  toplot@res_xyMarkerColors = colors
  plot_MEAN_TRND = scatterplot(wks_MEAN_TRND, toplot, var0, False, False, info)

  res = True
  res@gsLineColor = (/0.75, 0.75, 0.75/)
  res@gsLineThicknessF = 4.0
  res@tfPolyDrawOrder = "Draw"

  x = (/ssie_mean_ref, ssie_mean_ref/)
  y = (/toplot@res_trYMinF, toplot@res_trYMaxF/)
  ref0 = gsn_add_polyline(wks_MEAN_TRND, plot_MEAN_TRND, x, y, res)

  x = (/toplot@res_trXMinF, toplot@res_trXMaxF/)
  y = (/ssie_trnd_ref, ssie_trnd_ref/)
  ref1 = gsn_add_polyline(wks_MEAN_TRND, plot_MEAN_TRND, x, y, res)

  delete([/x, y/])

  draw(plot_MEAN_TRND)
  frame(wks_MEAN_TRND)

  ; ---------------------------------------
  ; scatter plot historical MEAN vs RCP YOD
  ; ---------------------------------------

  outf2 = "SSIE-MEAN_vs_YOD_" + var0 + "_extend_" \
    + diag_script_info@region + "_" + i2month + "_" + min(time_hist) \
    + "-" + max(time_rcp)

  wks_MEAN_YOD = get_wks(0, DIAG_SCRIPT, outf2)

  ; overwrite array --> retain most attributes
  toplot = (/SSIE_MEAN, YOD/)
  toplot&datasets = (/i2month + " " + diag_script_info@region \
                     + " SSIE averaged~C~over " \
                     + min(time_hist) + "-" + max(time_hist) \
                     + " [" + SSIE_MEAN@units + "]", \
                     "First year of near disappearance"/)
  caption2 = toplot&datasets(0) + " vs " + toplot&datasets(1) + "."
  toplot@res_gsnCenterString = ""
  toplot@res_trYMinF = 2000.
  toplot@res_trYMaxF = 2100.
  plot_MEAN_YOD = scatterplot(wks_MEAN_YOD, toplot, var0, False, False, info)

  ; reference datasets (observations)
  x = (/ssie_mean_ref, ssie_mean_ref/)
  y = (/toplot@res_trYMinF, toplot@res_trYMaxF/)
  ref2 = gsn_add_polyline(wks_MEAN_YOD, plot_MEAN_YOD, x, y, res)
  delete(x)

  draw(plot_MEAN_YOD)
  frame(wks_MEAN_YOD)

  ; ---------------------------------------
  ; scatter plot historical TRND vs RCP YOD
  ; ---------------------------------------

  toplot = (/SSIE_TRND, YOD/)
  toplot&datasets = (/"Trend in " + i2month + " " + diag_script_info@region \
                     + " SSIE~C~over " \
                     + min(time_hist) + "-" + max(time_hist) \
                     + " [" + SSIE_TRND@units + "]", \
                     "First year of near disappearance"/)
  caption3 = toplot&datasets(0) + " vs " + toplot&datasets(1) + "."

  outf3 = "SSIE-TRND_vs_YOD_" + var0 + "_extend_" \
    + diag_script_info@region + "_" + i2month + "_" \
    + min(time_hist) + "-" + max(time_rcp)

  wks_TRND_YOD = get_wks(0, DIAG_SCRIPT, outf3)
  toplot@res_trXMinF = decimalPlaces((min_trnd - 100), 0, False)
  toplot@res_trXMaxF = decimalPlaces((max_trnd + 100), 0, False)
  plot_TRND_YOD = scatterplot(wks_TRND_YOD, toplot, var0, False, False, info)

  ; reference datasets (observations)
  x = (/ssie_trnd_ref, ssie_trnd_ref/)
  ref3 = gsn_add_polyline(wks_TRND_YOD, plot_TRND_YOD, x, y, res)
  delete([/x, y, res/])

  draw(plot_TRND_YOD)
  frame(wks_TRND_YOD)

  ; ===========================================================================
  ; Provenance

  if (diag_script_info@region.eq."Arctic") then
    domain = "nh"
  else if (diag_script_info@region.eq."Antarctic") then
    domain = "sh"
  end if
  end if

  infiles = metadata_att_as_array(info, "filename")
  infiles_hist = infiles(idcs_hist)
  infiles_rcp = infiles(idcs_rcp)

  ; --------------------------------------------------------------

  SSIE_MEAN@var = "SSIE_mean"
  SSIE_MEAN@diag_script = DIAG_SCRIPT
  SSIE_TRND@var = "SSIE_trnd"
  SSIE_TRND@diag_script = DIAG_SCRIPT
  delete(SSIE_TRND&datasets@_FillValue)
  outf = config_user_info@work_dir + outf + ".nc"

  nc_ext = ncdf_write(SSIE_MEAN, outf)
  outf@existing = "append"
  nc_ext = ncdf_write(SSIE_TRND, outf)

  log_provenance(nc_ext, wks_MEAN_TRND@fullname, caption1, \
                 (/"mean", "trend"/), domain, "scatter", \
                 (/"senftleben_daniel", "lauer_axel"/), \
                 "massonnet12tc", infiles_hist)

  ; --------------------------------------------------------------

  YOD@var = "yod"
  YOD@diag_script = DIAG_SCRIPT
  delete(YOD&datasets@_FillValue)
  outf2 = config_user_info@work_dir + outf2 + ".nc"

  nc_ext2 = ncdf_write(SSIE_MEAN, outf2)
  outf2@existing = "append"
  nc_ext2 = ncdf_write(YOD, outf2)

  log_provenance(nc_ext2, wks_MEAN_YOD@fullname, caption2, \
                 (/"mean", "other"/), domain, "scatter", \
                 (/"senftleben_daniel", "lauer_axel"/), \
                 "", infiles)

  ; --------------------------------------------------------------

  outf3 = config_user_info@work_dir + outf3 + ".nc"

  nc_ext3 = ncdf_write(SSIE_TRND, outf3)
  outf3@existing = "append"
  nc_ext3 = ncdf_write(YOD, outf3)

  log_provenance(nc_ext3, wks_TRND_YOD@fullname, caption3, \
                 (/"trend", "other"/), domain, "scatter", \
                 (/"senftleben_daniel", "lauer_axel"/), \
                 "", infiles)

  leave_msg(DIAG_SCRIPT, "")

end
