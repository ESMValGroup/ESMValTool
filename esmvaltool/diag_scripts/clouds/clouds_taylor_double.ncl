; ############################################################################
; CLOUDS_TAYLOR_DOUBLE
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; PROJECT-NAME EMBRACE
; ############################################################################
; Description
;   Calculates the performance of models in reproducing 2-d fields of annual
;   mean or seasonal (DFJ, MAM, JJA, SON) mean cloud properties
;   (Taylor diagramms). The code is based on
;     1) clouds_taylor.ncl
;     2) perfmetrics_main.ncl
;     3) perfmetrics_taylor.ncl
;     4) perfmetrics_taylor_collect.ncl
;   Note: This code requires that all data are provided on the same grid.
;
; Required diag_script_info attributes (diagnostic specific)
;   none
;
; Optional diag_script_info attributes (diagnostic specific)
;   filename_add:    legacy feature: arbitrary string to be added to all
;                    filenames of plots and netcdf output produced
;                    (default = "")
;   multiobs_exclude: list of *observational* datasets to be excluded when
;                    calculating uncertainty estimates from multiple
;                    observational datasets (see also multiobs_uncertainty)
;   multiobs_uncertainty: calculate uncertainty estimates from multiple
;                       observational datasets (true, false); by default,
;                    all "obs", "obs6", "obs4mips" and "native6" datasets
;                    are used; any of such datasets can be explicitely
;                    excluded when also specifying "multiobs_exclude"
;   projectcolors:   colors for each projectgroups
;                    (e.g. (/"(/0.0, 0.0, 1.0/)", "(/1.0, 0.0, 0.0/)"/)
;   projectgroups:   calculated mmm per "projectgroup"
;                    (e.g. (/"cmip5", "cmip6")/)
;   styleset:        "CMIP5", "DEFAULT"
;                    (if not set, CLOUDS_TAYLOR_DOUBLE will create a color
;                    table and symbols for plotting)
;   timemean:        time averaging
;                    - annualclim (default) = 1 plot annual mean
;                    - seasonalclim = 4 plots (DJF, MAM, JJA, SON)
;   var:             short_name of variable to process (default = "" - use
;                    first variable in variable list)
;
; Required variable_info attributes (variable specific)
;   reference_dataset: name of reference data set
;
; Optional variable attributes (variable specific)
;   none
;
; Caveats
;   KNOWN LIMITATIONS
;     1) only 2-dim variables are currently supported
;
; Modification history
;   20211108-lauer_axel: written.
;
; ###########################################################################


load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/taylor_plot.ncl"
load "$diag_scripts/shared/dataset_selection.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  set_default_att(diag_script_info, "multiobs_exclude", "")
  set_default_att(diag_script_info, "multiobs_uncertainty", False)
  set_default_att(diag_script_info, "projectcolors", "(/0.0, 0.0, 1.0/)")
  set_default_att(diag_script_info, "projectgroups", "")
  set_default_att(diag_script_info, "timemean", "annualclim")
  set_default_att(diag_script_info, "valid_fraction", 0.5)
  set_default_att(diag_script_info, "var", "")

  projectgroups = diag_script_info@projectgroups
  if (projectgroups(0) .ne. "") then
    dim_GROUPS = dimsizes(projectgroups)
    dim_COLORS = dimsizes(diag_script_info@projectcolors)
    projectcolors_mmm = new(dim_COLORS, string)
    projectcolors = new(dim_GROUPS, string)
    newrgb = new(3, float)
    do i = 0, dim_GROUPS - 1
      if (i .lt. dim_COLORS) then
        projectcolors(i) = diag_script_info@projectcolors(i)
      else
        projectcolors(i) = "(/0.0, 0.0, 1.0/)"
      end if

      do j = 0, 2
        rgb = str_get_field(projectcolors(i), j + 1, ",")
        rgb = str_sub_str(rgb, "(", "")
        rgb = str_sub_str(rgb, ")", "")
        rgb = str_sub_str(rgb, "/", "")
        newrgb(j) = 0.75 * tofloat(rgb)
      end do

      projectcolors_mmm(i) = sprintf("(/%f, ", newrgb(0)) + \
        sprintf("%f, ", newrgb(1)) + sprintf("%f/)", newrgb(2))
    end do
    groupnames = projectgroups + "_mmm"

  else
    dim_GROUPS = 0
    groupnames = ""
  end if

  variables = metadata_att_as_array(variable_info, "short_name")
  numvars = dimsizes(variables)

  if (diag_script_info@var .eq. "") then
    var0 = variable_info[0]@short_name
  else
    var0 = diag_script_info@var
  end if

  varidx = ind(variables .eq. var0)
  if (ismissing(varidx)) then
    errstr = "diagnostic " + diag + " requires the following variable: var0"
    error_msg("f", DIAG_SCRIPT, "", errstr)
  end if

  flag_multiobs_unc = diag_script_info@multiobs_uncertainty
  multiobs_exclude = diag_script_info@multiobs_exclude

  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info0)
  if (isatt(variable_info[varidx], "reference_dataset")) then
    refname = variable_info[varidx]@reference_dataset
  end if
  if (isatt(variable_info[varidx], "alternative_dataset")) then
    refname2 = variable_info[varidx]@alternative_dataset
  end if
  names = metadata_att_as_array(info0, "dataset")
  projects = metadata_att_as_array(info0, "project")
  infiles = metadata_att_as_array(input_file_info, "filename")

  if (.not. flag_multiobs_unc) then
    ; check for reference model definition
    if (.not.isvar("refname")) then
      error_msg("f", DIAG_SCRIPT, "", "no reference dataset defined in recipe")
    end if
  end if

  ; print info on variable(s) to log file

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT)
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info("variable: " + variables(varidx))
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  timemean = diag_script_info@timemean
  numseas  = 1          ; default
  season   = (/"annual"/)

  ; time averaging: at the moment, only "annualclim" and "seasonalclim"
  ; are supported

  if (flag_multiobs_unc .and. timemean .ne. "annualclim") then
    log_info("multiobs_uncertainty = True is currently supported for annual" \
             + " means only (timemean = annualclim). Setting " \
             + " multiobs_uncertainty to False.")
    flag_multiobs_unc = False
  end if

  if (timemean.eq."seasonalclim") then
    numseas = 4
    delete(season)
    season = (/"DJF", "MAM", "JJA", "SON"/)
  end if

  ; create string for caption (netcdf provenance)

  allseas = season(0)
  do is = 1, numseas - 1
    allseas = allseas + "/" + season(i)
  end do

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

end

begin
  ; ========================================================================
  ; ========================== initialization ==============================
  ; ========================================================================

  ; get reference model

  ; find indices of all OBS and obs4mips datasets
  ; (treat ERA5 as special case)

  if (flag_multiobs_unc) then
    ; find indices of all OBS and obs4mips datasets (including "native6" ERA5)

    idxobs = get_obs(names, projects, multiobs_exclude)

    if (idxobs(0) .eq. -1) then
      flag_multiobs_unc = False
      log_info("No OBS or obs4mips datasets found. Setting " \
               + " multiobs_uncertainty to False.")
    else
      refname = "REF"
      ref_ind = dimsizes(names)
      names := array_append_record(names, (/refname/), 0)
      dim_MOD = dim_MOD + 1
    end if
  end if

  if (.not. flag_multiobs_unc) then
    ; if attribute is present, use it so correlations can be calculated
    if (isvar("refname")) then
      ; set reference model
      ref_ind = ind(names .eq. refname)
      if (ismissing(ref_ind)) then
        error_msg("f", DIAG_SCRIPT, "", "reference dataset (" + refname + \
                  ") not found.")
      end if
    else
      error_msg("f", DIAG_SCRIPT, "", "no reference dataset defined in recipe")
    end if
  end if

  ; get multi-model mean index (if present)

  mm_ind = ind(names .eq. "MultiModelMean")

  if (ismissing(mm_ind)) then
    mm_ind = -1
  end if

  if (isatt(diag_script_info, "filename_add")) then
    filename_add = "_" + diag_script_info@filename_add
  else
    filename_add = ""
  end if

  if (dim_GROUPS .gt. 0) then
    names := array_append_record(names, groupnames, 0)
    projects := array_append_record(projects, "group_" + projectgroups, 0)
    dim_MOD = dim_MOD + dim_GROUPS
  end if

  ; ========================================================================
  ; ============================ statistics ================================
  ; ========================================================================

  ; ---------------------------------------------------------
  ; if requested, calculate multi-observational mean and standard deviation

  if (flag_multiobs_unc) then
    nobs = dimsizes(idxobs)

    ; step 1: calculate multi-obs mean

    do i = 0, nobs - 1
      A0 = read_data(info0[idxobs(i)])

      ; calculate time average
      mean = time_operations(A0, -1, -1, "average", "annualclim", True)
      delete(A0)

      if (i .eq. 0) then
        dims = dimsizes(mean)
        newdims = new(dimsizes(dims) + 1, integer)
        newdims(0) = nobs
        newdims(1:dimsizes(newdims) - 1) = dims
        ref_tmp = new(newdims, float)
        delete(dims)
      end if
      ref_tmp(i, :, :) = mean
    end do
    delete(mean)
    ; note: we are using dim_avg_n_Warp so missing values are ignored
    ;       when averaging
    ref_avg = dim_avg_n_Wrap(ref_tmp, 0)
    delete(ref_tmp)
  end if

  ; read data and calculate time average

  ivar = varidx

  info = select_metadata_by_name(input_file_info, variables(ivar))

  ; loop over all datasets not including groups means (to be calculated
  ; separately)

  do imod = 0, dim_MOD - 1 - dim_GROUPS

    if (imod .ne. ref_ind .or. .not.flag_multiobs_unc) then
      data_temp = read_data(info[imod])
      data1 = time_operations(data_temp, -1, -1, "average", \
                              timemean, True)

      delete(data_temp)
    else
      data1 = ref_avg
      delete(ref_avg)
    end if

    if (isdefined("global_mask")) then
      if (numseas.eq.1) then
        data1 = data1 + global_mask
      else
        do is = 0, numseas - 1
          data1(is, :, :) = data1(is, :, :) + global_mask
        end do
      end if
    end if

    if (.not.isdefined("data")) then
      dim_data = array_append_record((/dim_MOD/), dimsizes(data1), 0)
      data = new(dim_data, float)
      data!0 = "model"
      data&model = names
      if (numseas.eq.1) then
        copy_VarCoords(data1, data(imod, :, :))
      else
        copy_VarCoords(data1, data(imod, :, :, :))
      end if
    end if

    dim_data = dimsizes(data)
    rank = dimsizes(dim_data)
    if (numseas.eq.1) then
      data(imod, :, :) = data1
    else
      data(imod, :, :, :) = data1
    end if
    delete(data1)

    ; debug output

    ; debugfile->$input_file_info@dataset(imod)$ = data1
  end do  ; imod loop

  ; ------------------------------------
  ; calculate group means (if requested)
  ; ------------------------------------

  j = 0
  do imod = dim_MOD - dim_GROUPS, dim_MOD - 1
    idxgroup = ind(projects .eq. projectgroups(j))
    if (.not.all(ismissing(idxgroup))) then
      isize = dimsizes(idxgroup)
      if (isize .gt. 1) then
        if (numseas.eq.1) then
          data(imod, :, :) = dim_avg_n(data(idxgroup, :, :), 0)
        else
          data(imod, :, :, :) = dim_avg_n(data(idxgroup, :, :, :), 0)
        end if
      else
        if (numseas.eq.1) then
          data(imod, :, :) = data(idxgroup, :, :)
        else
          data(imod, :, :, :) = data(idxgroup, :, :, :)
        end if
      end if
    end if
    delete(idxgroup)
    j = j + 1
  end do

  ; ------------------------------------

  ; number of used reference data sets

  dim_REF = 1
  if (isvar("refname2")) then
    dim_REF = dim_REF + 1
  end if

  ; define result variable

  val = new((/dim_MOD - dim_REF, 2, numseas/), float)
  val!0 = "models"
  val!1 = "statistic"
  val!2 = "time"
  val&models = names(ind((names .ne. refname)))
  val&statistic = (/"stddev_ratio", "correlation"/)
  val&time = ispan(0, numseas - 1, 1)

  ; extract data

  do is = 0, numseas - 1

    if (numseas.eq.1) then
      ref = data(ref_ind, :, :)
    else
      ref = data(ref_ind, is, :, :)
    end if

    ; loop over models (skip reference)

    mod_i = 0
    do imod = 0, dim_MOD - 1

      ; skip reference

      if (imod.eq.ref_ind) then
        log_info("Ref_dataset = " + data&model(imod))
        continue
      end if

      if (numseas.eq.1) then
        var = data(imod, :, :)
      else
        var = data(imod, is, :, :)
      end if

      ; common mask
      mask1 = where(.not.ismissing(var), 0., var@_FillValue)
      mask2 = where(.not.ismissing(ref), 0., ref@_FillValue)
      amask = mask1 + mask2
      delete(mask1)
      delete(mask2)
      refmasked = ref
      refmasked = refmasked + amask
      var = var + amask
      delete(amask)

      ; calculate rate of standard deviations

      ; val(mod_i, 0, is) = calculate_metric(var, ref, "stddev_ratio_nowgt")
      val(mod_i, 0, is) = calculate_metric(var, refmasked, "stddev_ratio")

      ; calculate pattern correlation

      ; val(mod_i, 1, is) = calculate_metric(var, ref, "correlation_nowgt")
      val(mod_i, 1, is) = calculate_metric(var, refmasked, "correlation")

      delete(var)
      delete(refmasked)
      mod_i = mod_i + 1

    end do  ; loop over datasets
  end do  ; loop over seasons

  ; attach attributes to the results

  val@title = "taylor metrics"
  val@long_name = "Taylor diagram"
  val@diag_script = (/DIAG_SCRIPT/)
  val@var = "var"

  ; ------------------- calculate/estimate RMS of observations ---------------

  rmsobs = new((/numseas/), float)

  do is = 0, numseas - 1
    if (numseas.eq.1) then
      ref = data(ref_ind, :, :)
    else
      ref = data(ref_ind, is, :, :)
    end if

    ; This code is equivalent to the function "calculate_metric"
    ; (statistics.ncl)
    weights = map_area(ref&lat, ref&lon)
    ; optional: no weights --> reset weights
    ; wgt1d = 1.0

    ; convert to 1-D arrays
    ref1d = ndtooned(ref)
    wgt1d = ndtooned(weights)
    wgt1d = wgt1d / dim_avg_n_Wrap(wgt1d, 0)

    rmsobs(is) = 0.0

    ; normalize RMS by dividing by standard deviation of obs (ref)
;    avg_ref = dim_avg_wgt_Wrap(ref1d, wgt1d, 1)
;    rmsobs(is) =  rmsobs(is) / sqrt(dim_avg_wgt_Wrap((ref1d - avg_ref) ^ 2, \
;                                                     wgt1d, 1))

  end do  ; loop over seasons

  if (isvar("var")) then
    delete(var)
  end if
  if (isvar("ref")) then
    delete(ref)
  end if
  if (isvar("weights")) then
    delete(weights)
  end if
  if (isvar("var1d")) then
    delete(var1d)
  end if
  if (isvar("ref1d")) then
    delete(ref1d)
  end if
  if (isvar("wgt1d")) then
    delete(wgt1d)
  end if

  ; ========================================================================
  ; ============================= plotting =================================
  ; ========================================================================

  nummods = dim_MOD - dim_REF

  colors  = new(nummods, string)
  markers = new(nummods, integer)

  ratio = new((/nummods, 1, numseas/), float)  ; standard deviation ratios
  ratio = val(:, 0, :)
  cc = new((/nummods, 1, numseas/), float)   ; correlations
  cc = val(:, 1, :)

  if (isatt(diag_script_info, "styleset")) then
    colortab = project_style(info0, diag_script_info, "colors")
    markertab = project_style(info0, diag_script_info, "markers")
  else
    colortab = (/"(/0.00, 0.00, 0.59/)", "(/0.00, 0.39, 1.00/)", \
                 "(/0.20, 1.00, 1.00/)", "(/0.20, 0.88, 0.00/)", \
                 "(/1.00, 0.88, 0.00/)", "(/1.00, 0.59, 0.00/)", \
                 "(/1.00, 0.20, 0.00/)", "(/0.59, 0.00, 0.00/)", \
                 "(/0.78, 0.00, 0.78/)", "(/0.59, 0.00, 0.59/)", \
                 "(/0.90, 0.90, 0.90/)", "(/0.70, 0.70, 0.70/)", \
                 "(/0.50, 0.50, 0.50/)", "(/0.30, 0.30, 0.30/)"/)
    markertab = (/16, 4, 5, 0/)
  end if

  plotfile = new((/numseas/), string)

  do is = 0, numseas - 1
    if (isvar("wks")) then
      delete(wks)
    end if

    if (isvar("legendwks")) then
      delete(legendwks)
    end if

    if (isvar("plot")) then
      delete(plot)
    end if

    wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_taylor_" + var0 \
                  + "_" + season(is) + filename_add)

    ; create new marker: filled star

    mstring = "z"
    fontnum = 35
    size   = 1.75  ; 1.5
    angle  = 0.0

    new_index = NhlNewMarker(wks, mstring, fontnum, 0.0, 0.0, \
                             1.0, size, angle)

    ; create new marker: filled dot

    mstring = "m"
    fontnum = 37
    size   = 1.0
    angle  = 0.0

    new_index_dot = NhlNewMarker(wks, mstring, fontnum, 0.0, 0.0, \
                                 1.0, size, angle)

    legendwks = wks

    i = 0
    idx = 0

    if (isatt(diag_script_info, "styleset")) then
      do ii = 0, dim_MOD - 1

        if (ii.eq.ref_ind) then
          continue
        end if

        ; define group means (if present)
        testidx = ind(groupnames .eq. names(ii))
        if (.not.ismissing(testidx)) then
          colors(i) = projectcolors_mmm(testidx)
          markers(i) = new_index
          i = i + 1
          continue
        end if

        colors(i)  = colortab(i)
        markers(i) = markertab(i)

        ; override multi-model mean (if present)
        if (ii .eq. mm_ind) then
          colors(i)  = "(/0.00, 0.00, 0.00/)"  ; black
          markers(i) = new_index
        end if

        ; override colors and markers of all group members
        ; (if groups are defined)
        testidx = ind(projectgroups .eq. projects(ii))
        if (.not.ismissing(testidx)) then
          colors(i) = projectcolors(testidx)
          markers(i) = new_index_dot  ; 16
        end if

        i = i + 1
      end do
    else
      do ii = 0, dim_MOD - 1

        if (ii.eq.ref_ind) then
          continue
        end if

        ; define group means (if present)
        testidx = ind(groupnames .eq. names(ii))
        if (.not.ismissing(testidx)) then
          colors(i) = projectcolors_mmm(testidx)
          markers(i) = new_index
          i = i + 1
          continue
        end if

        if (ii .eq. mm_ind) then
          colors(i)  = "(/0.00, 0.00, 0.00/)"  ; black
          markers(i) = new_index
          i = i + 1
        else
          do n = 0, dim_MOD / dimsizes(colortab)
            colors(i)  = colortab(idx)
            markers(i) = markertab(n)
            i = i + 1
            if (i.ge.dim_MOD) then
              break
            end if
          end do
          idx = idx + 1
          if (idx.ge.dimsizes(colortab)) then
            idx = dimsizes(colortab) - 1
          end if
        end if

        ; override colors and markers of all group members
        ; (if groups are defined)
        testidx = ind(projectgroups .eq. projects(ii))
        if (.not.ismissing(testidx)) then
          colors(i) = projectcolors(testidx)
          markers(i) = new_index_dot  ; 16
        end if

      end do
    end if

    ropts = True
    ropts@gsMarkerSizeF = 0.0125
    ropts@Colors     = colors
    ropts@Markers    = markers  ; marker styles
    ropts@caseLabels = val&models

    if (dim_GROUPS .gt. 0) then
      ropts@legend_filter = projectgroups
      ropts@projects = projects
    end if

    ropts@stnRad        = (/ 0.5, 1.5, 2.0, 2.5 /)
    ropts@OneX          = "1.00"
    ropts@ccRays        = (/ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, \
                           0.9, 0.95, 0.99 /)
    ropts@ccRays_color  = "Black"
    ropts@centerDiffRMS = True
    ropts@centerDiffRMS_color = "PaleGreen1"  ; "LightGray"

    ropts@printDataLabels = False  ; print numeric label for each data pt.
    ropts@reverseOrder  = True  ; reverse order of items in legend
    ropts@legendBox     = True  ; draw box around legend

    ; legend parameters

    ropts@legendXpos        = 0.6
    ropts@legendYpos        = -0.5
    ropts@legendWidth       = 0.2

    ropts@caseLabelsFontHeightF = 0.1
    ropts@plotSize              = 0.6

    ropts@legendExtraFile = False
    if (dim_MOD.ge.20) then
      ropts@caseLabelsFontHeightF = 0.06
      ropts@plotSize              = 0.45
    end if

    ; plot
    ropts@tiMainString   = var0 + " (" + season(is) + ")"  ; title
    ropts@tiMainOffsetYF = 0.025
    ; observational uncertainties as normalized RMS error
    ropts@rmsobs = rmsobs(is)

    plot = taylor_diagram(wks, legendwks, ratio(:, :, is), cc(:, :, is), \
                          ropts)

    plotfile(is) = wks@fullname
    log_info("Wrote " + plotfile)
  end do  ; is-loop (seasons)

  ; write netCDF output

  if (any(rmsobs .gt. 0.0)) then
    val@RMSE_observations = rmsobs
  end if

  nc_filename = work_dir + "clouds_taylor_" + var0 + filename_add + ".nc"
  nc_outfile = ncdf_write(val, nc_filename)

  ; ------------------------------------------------------------------------
  ; write provenance to netcdf output and plot file(s)
  ; ------------------------------------------------------------------------

  statistics = "clim"
  domain = "global"
  plottype = "taylor"

  do is = 0, numseas - 1
    ; note: because function log_provenance does not yet support to attach
    ;       different captions to netcdf (contains all seasons) and plots
    ;       (contain one season each), the caption cannot specifiy the
    ;       season plotted; using "annual" or "DJF/MAM/JJA/SON" instead.

    caption = "Taylor diagram for variable " + var0 + " (" + allseas \
              + "), reference = " + refname + "."

    log_provenance(nc_outfile, plotfile, caption, statistics, domain, \
                   plottype, "", "", infiles)
  end do

  leave_msg(DIAG_SCRIPT, "")

end
