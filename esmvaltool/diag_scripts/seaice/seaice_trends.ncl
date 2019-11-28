; #############################################################################
; Histograms of sea ice trend distribution (CMIP5): Figure 9.24 of IPPC AR5
; Author: Senftleben, Daniel (DLR, Germany)
; ESMVal project
; #############################################################################
; Description
;
;
; Required diag_script_info attributes (diagnostics specific)
;   month: selected month (1, 2, ..., 12) or annual mean ("A")
;   region: region to be analyzed ( "Arctic" or "Antarctic")
;
; Optional diag_script_info attributes (diagnostic specific)
;   fill_pole_hole: fill observational hole at North pole (default: False)
;   ref_model: array of references plotted as vertical lines
;
; Required variable_info attributes (variable specific)
;
; Caveats
;   CMOR units expected: areacello in m2, sic in %
;
; Modification history
;    20191025-lauer_axel: rewritten for v2.0
;    20170619-senftleben_daniel: written.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

load "$diag_scripts/seaice/seaice_aux.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  ; Get metadata items
  att = True
  att@mip = "OImon"
  info = select_metadata_by_atts(input_file_info, att)  ; variable
  var0 = info[0]@short_name
  datasets = metadata_att_as_array(info, "dataset")
  ensembles = metadata_att_as_array(info, "ensemble")
  exps = metadata_att_as_array(info, "exp")
  dim_MOD = dimsizes(datasets)
  att@mip = "fx"
  areainfo = select_metadata_by_atts(input_file_info, att)  ; area
  areadatasets = metadata_att_as_array(areainfo, "dataset")
  delete(att)

  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"month", "region"/))

  ; Check for reference dataset
  if (isatt(info[0], "reference_dataset")) then
    ref_models = info[0]@reference_dataset
    rcp_tag = ""
    dim_REF = dimsizes(ref_models)
  else
    ref_models = ""
    rcp_tag = ", " + exps(0)
    dim_REF = 0
  end if

  ; Check fill_pole_hole option
  if (.not.isatt(diag_script_info, "fill_pole_hole")) then
    diag_script_info@fill_pole_hole = False
  else
    if (diag_script_info@fill_pole_hole) then
      if (diag_script_info@region.eq."Antarctic") then
        diag_script_info@fill_pole_hole = False
        error_msg("w", DIAG_SCRIPT, "", "No pole hole at Antarctic " + \
                  "- setting diag_script_info@fill_pole_hole to False")
      end if
    end if
  end if

  ; Create output plot directory
  if (config_user_info@write_plots) then
    plot_dir = config_user_info@plot_dir
    system("mkdir -p " + plot_dir)
  end if

  ; Create output directory (work)
  system("mkdir -p " + config_user_info@work_dir)

  ; Plot file type
  file_type = config_user_info@output_file_type

  ; Create arrays to collect SSIE mean and trend from each model
  SSIE_trnd = new(dim_MOD - dim_REF, float)
  SSIE_trnd@units = "million sq km / decade"
  if (dim_REF .gt. 0) then
    SSIE_trnd_ref = new(dim_REF, float)
    SSIE_trnd_ref@units = "million sq km / decade"
  end if

  ; Determine time range
  start_years = metadata_att_as_array(info, "start_year")
  start_year  = min(start_years)
  end_years   = metadata_att_as_array(info, "end_year")
  end_year    = max(end_years)
  all_years   = ispan(start_year, end_year, 1)

  ; strings for caption and plot title

  year0 = tostring(start_year)
  if (start_year .ne. max(start_years)) then
    year0 = year0 + "..." + max(start_years)
  end if

  year1 = tostring(end_year)
  if (end_year .ne. min(end_years)) then
    year1 = min(end_years) + "..." + year1
  end if

  ; Set up counter(s)
  ref_count = 0

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

    ; Calculate SSIE trend
    trend_tmp = dtrend(SSIE, True)
    if (dim_REF .gt. 0) .and. \
      any(ref_models.eq.datasets(imod)) then
      SSIE_trnd_ref(ref_count) = trend_tmp@slope * 10.  ; per decade
      ref_count = ref_count + 1
    else
      SSIE_trnd(imod - ref_count) = trend_tmp@slope * 10.  ; per decade
      print(SSIE_trnd(imod - ref_count))
    end if

    ; Clean up
    delete([/sic, area, arealat, SSIE/])

  end do  ; imod

  ; ===========================================================================
  ; Plot histogram

  smonths = (/"January", "February", "March", "April", \
              "May", "June", "July", "August", \
              "September", "October", "November", "December"/)

  if (diag_script_info@month .eq. "A") then
    i2month = "annual"
  else
    i2month = smonths(toint(diag_script_info@month) - 1)
  end if

  out_ext = "trend_" + var0 + "_extend_" \
            + diag_script_info@region + "_" + i2month + "_histogram"

  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, out_ext)

  res = True
  res@tiMainString = i2month + " " + diag_script_info@region \
    + " sea ice extent trends " + year0 + "-" + year1 + rcp_tag
  res@tiXAxisString = "10~S~6~N~ km~S~2~N~ per decade"
  res@tiYAxisString = "Number of models"
  xAxis = fspan(-1.6, 0.2, 19)
  res@gsnHistogramBinIntervals = xAxis
  res@tmXBLabelAngleF = 265
  res@gsnHistogramBarWidthPercent = 100.
  res@gsFillColor = "dodgerblue3"
  res@gsnFrame = False
  res@gsnDraw = False
  plot = gsn_histogram(wks, SSIE_trnd, res)

  ; Add reference data as vertical lines
  if (dim_REF .gt. 0) then
    getvalues plot  ; Get plot parameters
      "trYMinF" : ymin
      "trYMaxF" : ymax
    end getvalues
    xmin = min(xAxis)
    xmax = max(xAxis)
    res_ref = True
    res_ref@gsLineThicknessF = 6.
    ref_colors = (/"firebrick1", "black", "green2", "gold"/)
    do iref = 0, dim_REF - 1
      res_ref@gsLineColor = ref_colors(iref)
      ref_val = SSIE_trnd_ref(iref)
      x = (ref_val - xmin) / (xmax - xmin)
      if (iref.eq.0) then
        ref_line0 = gsn_add_polyline(wks, plot, (/x, x/), \
                                     (/ymin, ymax/), res_ref)
      else if (iref.eq.1) then
        ref_line1 = gsn_add_polyline(wks, plot, (/x, x/), \
                                     (/ymin, ymax/), res_ref)
      else if (iref.eq.2) then
        ref_line2 = gsn_add_polyline(wks, plot, (/x, x/), \
                                     (/ymin, ymax/), res_ref)
      else if (iref.eq.3) then
        ref_line3 = gsn_add_polyline(wks, plot, (/x, x/), \
                                     (/ymin, ymax/), res_ref)
      else
        break  ; four ref lines should be enough ;)
      end if
      end if
      end if
      end if

      delete([/ref_val, x/])
    end do
  end if

  draw(plot)
  frame(wks)

  ; ===========================================================================
  ; Provenance

  if (diag_script_info@region.eq."Arctic") then
    domain = "nh"
  else if (diag_script_info@region.eq."Antarctic") then
    domain = "sh"
  end if
  end if

  infiles = metadata_att_as_array(info, "filename")

  alltrends = new(dim_MOD, float)
  alltrends(0:dim_MOD - dim_REF - 1) = SSIE_trnd
  if (dim_REF .gt. 0) then
    alltrends(dim_MOD - dim_REF - 1:dim_MOD - 1) = SSIE_trnd_ref
  end if
  alltrends!0 = "dataset"
  alltrends&dataset = datasets
  alltrends@var = "SSIE_trnd"
  alltrends@diag_script = DIAG_SCRIPT
  caption = res@tiMainString

  nc_ext = ncdf_write(alltrends, config_user_info@work_dir + out_ext + ".nc")

  log_provenance(nc_ext, wks@fullname, caption, \
                 (/"trend"/), domain, "histogram", \
                 (/"senftleben_daniel", "lauer_axel"/), \
                 "", infiles)

  leave_msg(DIAG_SCRIPT, "")

end
