; CLOUDS_LIFRAC_SCATTER
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; ############################################################################
; Description
;   Calculates average liquid / ice fraction vs. temperature.
;
; Required diag_script_info attributes (diagnostic specific)
;   none
;
; Optional diag_script_info attributes (diagnostic specific)
;   filename_add:       optionally add this string to plot filesnames
;   min_mass:           minimum cloud condensate (same units as clw, cli)
;   mm_mean_median:     calculate multi-model mean and meadian
;   nbins: number of equally spaced bins (ta (x-axis)), default = 20
;   panel_labels:       label individual panels (true, false)
;   PanelTop:           manual override for "@gnsPanelTop" used by panel
;                       plot(s)
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
;   20230117-lauer_axel: added support for ICON (code from Manuel)
;   20210302-lauer_axel: written.
;
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/dataset_selection.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  diag = "clouds_lifrac_scatter.ncl"
  variables = get_unique_values(metadata_att_as_array(variable_info, \
                                                      "short_name"))
  dim_VAR = dimsizes(variables)
  refname = new(dim_VAR, string)

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + variables + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; make sure required variables are available

  varlist = (/"clw", "cli", "ta"/)
  idx = new(dimsizes(varlist), integer)

  do i = 0, dimsizes(varlist) - 1
    idx(i) = ind(variables .eq. varlist(i))
  end do

  if (any(ismissing(idx))) then
    errstr = "diagnostic " + diag + " requires the following variable(s): " \
             + str_join(varlist, ", ")
    error_msg("f", DIAG_SCRIPT, "", errstr)
  end if

  ; save input files for writing provenance

  infiles = metadata_att_as_array(input_file_info, "filename")

  ; get reference datasets (if present) and check that number of datasets
  ; is equal for each variable

  do i = 0, dim_VAR - 1
    var = variables(idx(i))
    var_info = select_metadata_by_name(variable_info, var)
    var_info := var_info[0]
    if (isatt(var_info, "reference_dataset")) then
      refname(i) = var_info@reference_dataset
    end if
    info = select_metadata_by_name(input_file_info, var)
    if (i .eq. 0) then
      dim_MOD = ListCount(info)
    else
      dim_test = ListCount(info)
      if (dim_test .ne. dim_MOD) then
        error_msg("f", DIAG_SCRIPT, "", "number of datasets for variable " \
                  + variables(i) + " does not match number of datasets for " \
                  + variables(0))
      end if
    end if
    delete(info)
    delete(var)
    delete(var_info)
  end do

  delete(idx)

  ; Set default values for non-required diag_script_info attributes

  set_default_att(diag_script_info, "filename_add", "")
  set_default_att(diag_script_info, "panel_labels", True)
  set_default_att(diag_script_info, "min_mass", 1.0e-15)
  set_default_att(diag_script_info, "mm_mean_median", True)
  set_default_att(diag_script_info, "nbins", 20)

  if (diag_script_info@filename_add .ne. "") then
    filename_add = "_" + diag_script_info@filename_add
  else
    filename_add = ""
  end if

  nbins = toint(diag_script_info@nbins)

  panel_labels = diag_script_info@panel_labels
  min_mass = diag_script_info@min_mass
  mm_mean_median = diag_script_info@mm_mean_median

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

end

begin
  ; ############
  ; # get data #
  ; ############

  info_clw = select_metadata_by_name(input_file_info, "clw")
  names_clw = metadata_att_as_array(info_clw, "dataset")
  projects_clw = metadata_att_as_array(info_clw, "project")
  info_cli = select_metadata_by_name(input_file_info, "cli")
  names_cli = metadata_att_as_array(info_cli, "dataset")
  projects_cli = metadata_att_as_array(info_cli, "project")
  info_ta = select_metadata_by_name(input_file_info, "ta")
  names_ta = metadata_att_as_array(info_ta, "dataset")
  projects_ta = metadata_att_as_array(info_ta, "project")

  refidx_clw = ind(names_clw .eq. refname(0))
  refidx_cli = ind(names_cli .eq. refname(1))
  refidx_ta = ind(names_ta .eq. refname(2))

  if (ismissing(refidx_clw) .or. ismissing(refidx_cli) .or. \
      ismissing(refidx_ta)) then
    if (ismissing(refidx_clw)) then
      str = refname(0)
    end if
    if (ismissing(refidx_cli)) then
      if (isdefined("str")) then
        str = str + " + " + refname(1)
      else
        str = refname(1)
      end if
    end if
    if (ismissing(refidx_ta)) then
      if (isdefined("str")) then
        str = str + " + " + refname(2)
      else
        str = refname(2)
      end if
    end if

    error_msg("f", DIAG_SCRIPT, "", "the following reference dataset(s) " \
              + "are not available: " + str)
  end if

  add_dim_MOD = 0

  if (mm_mean_median) then
    ; check if enough model datasets are available to calculate multi-model
    ; mean and median ice/liquid fractions

    ; find all indices of models w/o MultiModelMean/MultiModelMedian
    ; (if present)
    idxmod = get_mod(names_clw, projects_clw)

    if ((idxmod(0) .eq. -1) .or. (dimsizes(idxmod) .le. 2)) then
      log_info("Not enough model datasets to calculate multi-model " \
               + "mean and median. Setting of 'mm_mean_median' will" \
               + " be ignored.")
      mm_mean_median = False
    else  ; more than one model found
      add_dim_MOD = 2
    end if
    delete(idxmod)
  end if

  resulti = new((/dim_MOD + add_dim_MOD, nbins/), float)
  resultl = new((/dim_MOD + add_dim_MOD, nbins/), float)
  bincenter = new((/nbins/), float)

  xmax = 300.  ; Kelvin (27°C)
  xmin = 230.  ; Kelvin (-43°C)
  binsize = tofloat(xmax - xmin) / nbins

  do n = 0, nbins - 1
    x0 = n * binsize
    x1 = x0 + binsize
    bincenter(n) = xmin + 0.5 * (x0 + x1)
  end do

  ; -------------------------------------------------------------------------

  ; create index vector with reference dataset as first entry (index 0)
  ; so reference icefrac is calculated first

  idataset = ispan(0,  dim_MOD - 1, 1)
  idataset0 = idataset(0)
  i1 = ind(idataset .eq. refidx_clw)
  idataset(0) = idataset(i1)
  idataset(i1) = idataset0

  do ii = 0, dim_MOD - 1

    i = idataset(ii)

    atts_cli = True
    atts_cli@short_name = "cli"

    atts_clw = True
    atts_clw@short_name = "clw"

    atts_ta = True
    atts_ta@short_name = "ta"

    ; reference datasets may have different names
    if (i .eq. refidx_clw) then
      atts_clw@dataset = refname(0)
      atts_cli@dataset = refname(1)
      atts_ta@dataset = refname(2)
    ; all other datasets: force same dataset name for clw, cli and ta
    else
      atts_clw@dataset = names_clw(i)
      atts_cli@dataset = names_clw(i)
      atts_ta@dataset = names_clw(i)
    end if

    ; read cli

    info = select_metadata_by_atts(input_file_info, atts_cli)
    cli = read_data(info[0])
    delete(info)

    ; read clw

    info = select_metadata_by_atts(input_file_info, atts_clw)
    clw = read_data(info[0])
    delete(info)

    ; read ta

    info = select_metadata_by_atts(input_file_info, atts_ta)
    ta = read_data(info[0])
    delete(info)

    ; check dimensions

    dims_clw = dimsizes(clw)
    dims_cli = dimsizes(cli)
    dims_ta = dimsizes(ta)

    dimerror = False

    if (dimsizes(dims_cli) .eq. dimsizes(dims_clw)) then
      if (any(dims_cli - dims_clw .ne. 0)) then
        dimerror = True
      end if
    else
      dimerror = True
    end if

    if (dimsizes(dims_cli) .eq. dimsizes(dims_ta)) then
      if (any(dims_cli - dims_ta .ne. 0)) then
        dimerror = True
      end if
    else
      dimerror = True
    end if

    if (dimerror) then
      error_msg("f", DIAG_SCRIPT, "", "dimensions of datasets " \
                + atts_cli@dataset + " (variable cli), " \
                + atts_clw@dataset + " (variable clw), " \
                + atts_ta@dataset + " (variable ta) do not match.")
    end if

    delete(dims_cli)
    delete(dims_clw)
    delete(dims_ta)

    ; calculate ice fraction

    ; filter valid values (needed for some models)

    cli = where(cli .lt. 0.0, cli@_FillValue, cli)
    cli = where(isnan_ieee(cli), cli@_FillValue, cli)
    clw = where(clw .lt. 0.0, clw@_FillValue, clw)
    clw = where(isnan_ieee(clw), clw@_FillValue, clw)

    mass = cli + clw
    mass = where(mass .lt. min_mass, mass@_FillValue, mass)

    icefrac = 100.0 * cli / mass
    liqfrac = 100.0 * clw / mass

    ; ========================================================================

    ; output zonal means of cli with max icefrac limited to observed icefrac
    ; ----------------------------------------------------------------------
    ; as modeled and observed timeseries of icefrac may differ in length and
    ; years covered, we use average seasonal cycles instead of the full time
    ; series

    copy_VarMeta(cli, icefrac)
    icedims = dimsizes(cli)
    icedimnames = getvardims(cli)
    lonidx = ind(icedimnames .eq. "lon")
    if (ismissing(lonidx)) then
      error_msg("f", DIAG_SCRIPT, "", "no lon dimension")
    end if

    if (ii .eq. 0) then
      system("rm " + config_user_info@work_dir + "zonal.nc")
      zonalfile = addfile(config_user_info@work_dir + "zonal.nc", "c")
    end if

    cli_avg = time_operations(cli, -1, -1, "average", "monthlyclim", True)
    zcli = dim_avg_n_Wrap(cli_avg, lonidx)
    delete(zcli&month)
    zcli!0 = "time"
    zcli&time = ispan(0, 11, 1)

    var = "f_" + tostring(atts_cli@dataset)
    zonalfile->$var$ = zcli

    delete(cli_avg)
    delete(zcli)

    ; ========================================================================

    cli1d = ndtooned(cli)

    delete(clw)
    delete(mass)

    icefrac_units = "%"
    liqfrac_units = "%"
    icefrac_long_name = "cloud ice fraction"
    liqfrac_long_name = "cloud liquid fraction"
    icefrac_var = "icefrac"
    liqfrac_var = "liqfrac"
    ref_ind = refidx_cli
    if (ismissing(ref_ind)) then
      ref_ind = -1
    end if
    names = names_cli
    projects = projects_cli

    ; if reference datasets for clw, cli, ta are from different sources
    if (refname(0) .ne. refname(1)) then
      names(refidx_cli) = refname(0) + "/" + refname(1)
    end if
    if (refname(0) .ne. refname(2) .and. refname(1) .ne. refname(2)) then
      names(refidx_cli) = names(refidx_cli) + "/" + refname(2)
    end if

    ; ========================================================================

    ; bin data

    ice1d = ndtooned(icefrac)
    delete(icefrac)
    liq1d = ndtooned(liqfrac)
    delete(liqfrac)
    ta1d = ndtooned(ta)
    ta_units = ta@units
    delete(ta)

    resulti@_FillValue = ice1d@_FillValue
    resultl@_FillValue = ice1d@_FillValue

    cli_limited1d = cli1d

    do n = 0, nbins - 1
      resulti(i, n) = 0.0
      resultl(i, n) = resultl@_FillValue
      x0 = xmin + n * binsize
      x1 = x0 + binsize
      idx0 = ind((ta1d .gt. x0) .and. (ta1d .le. x1))
      if (.not.all(ismissing(idx0))) then
        nv = num(.not.ismissing(ice1d(idx0)))
        if (nv .gt. 0) then
          resulti(i, n) = avg(ice1d(idx0))
          ; ------------------------------------------------------------------

          iref = idataset(0)
          ifrac_ref = resulti(iref, n)
          ifrac_mod = resulti(i, n)

          if (.not.ismissing(ifrac_ref) .and. .not.ismissing(ifrac_mod)) then
            if (ifrac_mod .gt. 0.0) then
              corr_fac = ifrac_ref / ifrac_mod
              if (corr_fac .lt. 1.0) then
                cli_limited1d(idx0) = cli1d(idx0) * corr_fac
              end if
            end if
          end if
          ; ------------------------------------------------------------------
        end if
        nv = num(.not.ismissing(liq1d(idx0)))
        if (nv .gt. 0) then
          resultl(i, n) = avg(liq1d(idx0))
        end if
      end if
      delete(idx0)
    end do

    cli_limited = onedtond(cli_limited1d, icedims)
    copy_VarMeta(cli, cli_limited)

    cli_avg = time_operations(cli_limited, -1, -1, "average", \
                              "monthlyclim", True)
    zcli = dim_avg_n_Wrap(cli_avg, lonidx)
    delete(zcli&month)
    zcli!0 = "time"
    zcli&time = ispan(0, 11, 1)

    var = "l_" + tostring(atts_cli@dataset)
    zonalfile->$var$ = zcli

    delete(cli_avg)
    delete(zcli)
    delete(cli_limited)
    delete(cli_limited1d)

    delete(ice1d)
    delete(liq1d)
    delete(ta1d)
    delete(cli1d)
    delete(icedims)
    delete(icedimnames)

    delete(cli)

  ; ==========================================================================

  end do  ; loop over all datasets (ii)

  ; calculate average and median over all models excluding possibly present
  ; MultiModelMean, MultiModelMedian and reference dataset

  if (mm_mean_median) then
    index_my_mm = dim_MOD
    index_my_med = dim_MOD + 1

    xnames = array_append_record(names, (/"Average", "Median"/), 0)
    delete(names)
    names = xnames
    delete(xnames)

    mm_ind = ind(names .eq. "MultiModelMean")
    mmed_ind = ind(names .eq. "MultiModelMedian")

    if (ismissing(mm_ind)) then
      mm_ind = -1
    end if
    if (ismissing(mmed_ind)) then
      mmed_ind = -1
    end if

    imod = ispan(0, dim_MOD - 1, 1)
    idx0 = ind((imod .ne. ref_ind) .and. (imod .ne. mm_ind) .and. \
               (imod .ne. mmed_ind))

    mmed25i = new(nbins, float)
    mmed75i = new(nbins, float)
    mmed25l = new(nbins, float)
    mmed75l = new(nbins, float)

    mmstdi = new(nbins, float)
    mmstdl = new(nbins, float)

    if (.not.all(ismissing(idx0))) then
      do n = 0, nbins - 1
        nv = num(.not.ismissing(resulti(idx0, n)))
        if (nv .gt. 2) then
          stat = stat_dispersion(resulti(idx0, n), False)
          resulti(index_my_mm, n) = stat(0)   ; mean
          mmstdi(n) = stat(1)                 ; standard deviation
          mmed25i(n) = stat(6)                ; lower quartile
          resulti(index_my_med, n) = stat(8)  ; median
          mmed75i(n) = stat(10)               ; upper quartile
          ; if number of datasets < 4, quartiles cannot be calculated
          if (ismissing(mmed25i(n)) .or. ismissing(mmed75i(n))) then
            mmed25i(n) = 0.0
            mmed75i(n) = 0.0
          end if
        else
          resulti(index_my_mm, n) = 0.0
          mmstdi(n) = 0.0
          mmed25i(n) = 0.0
          resulti(index_my_med, n) = 0.0
          mmed75i(n) = 0.0
        end if
        nv = num(.not.ismissing(resultl(idx0, n)))
        if (nv .gt. 2) then
          stat = stat_dispersion(resultl(idx0, n), False)
          resultl(index_my_mm, n) = stat(0)   ; mean
          mmstdl(n) = stat(1)                 ; standard deviation
          mmed25l(n) = stat(6)                ; lower quartile
          resultl(index_my_med, n) = stat(8)  ; median
          mmed75l(n) = stat(10)               ; upper quartile
          ; if number of datasets < 4, quartiles cannot be calculated
          if (ismissing(mmed25l(n)) .or. ismissing(mmed75l(n))) then
            mmed25l(n) = 0.0
            mmed75l(n) = 0.0
          end if
        else
          resultl(index_my_mm, n) = 0.0
          mmstdl(n) = 0.0
          mmed25l(n) = 0.0
          resultl(index_my_med, n) = 0.0
          mmed75l(n) = 0.0
        end if
      end do
    end if

    if (isdefined("stat")) then
      delete(stat)
    end if
    delete(idx0)
    delete(imod)
  else
    index_my_mm = -1
    index_my_med = -1
  end if

;  if (mm_ind .ge. 0) then
;    mmstdi = new(nbins, float)
;    mmstdl = new(nbins, float)
;    mmstdi = 0.0
;    mmstdl = 0.0
;    imod = ispan(0, dim_MOD - 1, 1)
;    idx0 = ind((imod .ne. ref_ind) .and. (imod .ne. mm_ind) .and. \
;               (imod .ne. mmed_ind))
;    if (.not.all(ismissing(idx0))) then
;      do n = 0, nbins - 1
;        nv = num(.not.ismissing(resulti(idx0, n)))
;        if (nv .gt. 2) then
;          mmstdi(n) = stddev(resulti(idx0, n))
;        end if
;        nv = num(.not.ismissing(resultl(idx0, n)))
;        if (nv .gt. 2) then
;          mmstdl(n) = stddev(resultl(idx0, n))
;        end if
;      end do
;    end if
;    delete(idx0)
;    delete(imod)
;  end if
;
;  ; if MultiModelMedian is present, calculate 25 and 75 pecentiles
;
;  if (mmed_ind .ge. 0) then
;    mmed25i = new(nbins, float)
;    mmed75i = new(nbins, float)
;    mmed25l = new(nbins, float)
;    mmed75l = new(nbins, float)
;    mmed25i = 0.0
;    mmed75i = 0.0
;    mmed25l = 0.0
;    mmed75l = 0.0
;    imod = ispan(0, dim_MOD - 1, 1)
;    idx0 = ind((imod .ne. ref_ind) .and. (imod .ne. mm_ind) .and. \
;               (imod .ne. mmed_ind))
;    if (.not.all(ismissing(idx0))) then
;      do n = 0, nbins - 1
;        nv = num(.not.ismissing(resulti(idx0, n)))
;        if (nv .gt. 2) then
;          stat = stat_dispersion(resulti(idx0, n), False)
;          mmed25i(n) = stat(6)
;          mmed75i(n) = stat(10)
;        end if
;        nv = num(.not.ismissing(resultl(idx0, n)))
;        if (nv .gt. 2) then
;          stat = stat_dispersion(resultl(idx0, n), False)
;          mmed25l(n) = stat(6)
;          mmed75l(n) = stat(10)
;        end if
;      end do
;    end if
;    delete(idx0)
;    delete(imod)
;    delete(stat)
;  end if

  ; ###########################################
  ; # netCDF output                           #
  ; ###########################################

  nc_filename_i = work_dir + "clouds_lifrac_scatter_ice" + \
    filename_add + ".nc"
  nc_filename_l = work_dir + "clouds_lifrac_scatter_liq" + \
    filename_add + ".nc"

  resulti!0 = "model"
  resulti!1 = "bin"
  resulti&model = str_sub_str(names, "/", "-")
  resulti&bin = bincenter
  resulti@diag_script = (/DIAG_SCRIPT/)

  copy_VarMeta(resulti, resultl)

  resulti@var = icefrac_var
  resulti@var_long_name = icefrac_long_name
  resulti@var_units = icefrac_units

  resultl@var = liqfrac_var
  resultl@var_long_name = liqfrac_long_name
  resultl@var_units = liqfrac_units

  nc_outfile_i = ncdf_write(resulti, nc_filename_i)
  nc_outfile_l = ncdf_write(resultl, nc_filename_l)

  ; ###########################################
  ; # create the plots                        #
  ; ###########################################

  plots_i = new(dim_MOD + add_dim_MOD, graphic)
  plots_l = new(dim_MOD + add_dim_MOD, graphic)
  res = True

  wks_i = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_lifrac_scatter_" + \
                  "ice" + filename_add)
  wks_l = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_lifrac_scatter_" + \
                  "liq" + filename_add)

  data_i = new((/2, nbins/), float)
  data_l = new((/2, nbins/), float)
  if (ref_ind .gt. 0) then
    data_i(0, :) = resulti(ref_ind, :)
    data_l(0, :) = resultl(ref_ind, :)
  else
    data_i(0, :) = data@_FillValue
    data_l(0, :) = data@_FillValue
  end if

  do ii = 0, dim_MOD + add_dim_MOD - 1
    if (ii .eq. refidx_ta) then
      continue
    end if

    res@gsnDraw        = False  ; do not draw yet
    res@gsnFrame       = False  ; don't advance frame

    res@xyMarkLineMode    = "Markers"
    res@xyMarkers         = (/4, 16/)
    res@xyMonoMarkerColor = False
    res@xyMarkerColors    = (/"black", "red"/)
    res@xyMarkerSizeF     = 0.01
    res@tmLabelAutoStride = True

    res@tiMainFontHeightF = 0.025
    res@tiYAxisFontHeightF = 0.025
    res@tiXAxisFontHeightF = 0.025
    res@tiXAxisString = varlist(2) + " (" + ta_units + ")"

    data_i(1, :) = resulti(ii, :)
    data_l(1, :) = resultl(ii, :)

    res@tiMainString = names(ii)

    res@tiYAxisString = icefrac_var + " (" + icefrac_units + ")"
    plots_i(ii) = gsn_csm_xy(wks_i, resulti&bin, data_i, res)
    res@tiYAxisString = liqfrac_var + " (" + liqfrac_units + ")"
    plots_l(ii) = gsn_csm_xy(wks_l, resultl&bin, data_l, res)

    if (ii .eq. index_my_mm) then
      res_std = True
      res_std@gsnDraw        = False  ; do not draw yet
      res_std@gsnFrame       = False  ; don't advance frame
      res_std@gsnXYFillColors = "lightpink"
      res_std@xyLineColor = -1  ; Make lines transparent
      mmstddev = new((/2, nbins/), float)
      mmstddev(0, :) = resulti(ii, :) - mmstdi(:)
      mmstddev(1, :) = resulti(ii, :) + mmstdi(:)
      plotstd_i = gsn_csm_xy(wks_i, resulti&bin, mmstddev, res_std)
      overlay(plots_i(ii), plotstd_i)
      mmstddev(0, :) = resultl(ii, :) - mmstdl(:)
      mmstddev(1, :) = resultl(ii, :) + mmstdl(:)
      plotstd_l = gsn_csm_xy(wks_l, resultl&bin, mmstddev, res_std)
      delete(mmstddev)
      overlay(plots_l(ii), plotstd_l)
    end if

    if (ii .eq. index_my_med) then
      res_med = True
      res_med@gsnDraw        = False  ; do not draw yet
      res_med@gsnFrame       = False  ; don't advance frame
      res_med@gsnXYFillColors = "lightpink"
      res_med@xyLineColor = -1  ; Make lines transparent
      mmed = new((/2, nbins/), float)
      mmed(0, :) = mmed25i(:)
      mmed(1, :) = mmed75i(:)
      plotmed_i = gsn_csm_xy(wks_i, resulti&bin, mmed, res_med)
      overlay(plots_i(ii), plotmed_i)
      mmed(0, :) = mmed25l(:)
      mmed(1, :) = mmed75l(:)
      plotmed_l = gsn_csm_xy(wks_l, resultl&bin, mmed, res_med)
      delete(mmed)
      overlay(plots_l(ii), plotmed_l)
    end if

    draw(plots_i(ii))
    frame(wks_i)
    draw(plots_l(ii))
    frame(wks_l)
  end do

  pres                      = True    ; needed to override
                                      ; panelling defaults
  pres@gsnPanelCenter                   = False

  idx0 = ind(.not.ismissing(plots_i))
  n = dimsizes(idx0)
  pres@gsnPanelFigureStrings = names(idx0)
  pres@gsnPanelFigureStringsFontHeightF = min((/0.008, 0.008 * 6.0 \
                                              / tofloat((dim_MOD + 1) / 2)/))
  pres@lbLabelFontHeightF               = min((/0.01, 0.01 * 6.0 \
                                              / tofloat((dim_MOD + 1) / 2)/))

  outfile_i = panelling(wks_i, plots_i(idx0), (n + 3) / 4, 4, pres)
  log_info("Wrote " + outfile_i)

  outfile_l = panelling(wks_l, plots_l(idx0), (n + 3) / 4, 4, pres)
  log_info("Wrote " + outfile_l)
  delete(idx0)

  ; ==========================================================================

  ; ----------------------------------------------------------------------
  ; write provenance to netcdf output (and plot file)
  ; ----------------------------------------------------------------------

  statistics = (/"clim", "mean", "pdf"/)
  domain = "reg"
  plottype = "scatter"
  caption_i = "Scatterplot of air temperature (x) vs. cloud ice fraction (y)."
  caption_l = "Scatterplot of air temperature (x) vs. " + \
    "cloud liquid fraction (y)."
  log_provenance(nc_outfile_i, outfile_i, caption_i, statistics, \
                 domain, plottype, "", "", infiles)
  log_provenance(nc_outfile_l, outfile_l, caption_l, statistics, \
                 domain, plottype, "", "", infiles)

  leave_msg(DIAG_SCRIPT, "")

end
