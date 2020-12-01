; #############################################################################
; MAIN SCRIPT FOR PERFORMANCE METRICS
; Authors: Mattia Righi (DLR, Germany) and Franziska Frank (DLR, Germany)
; ESMVal project
; #############################################################################
;
; Description
;    Calculates and (optionally) plots annual/seasonal cycles, zonal means,
;    lat-lon fields and time-lat-lon fields from input 2D/3D monthly data.
;    The calculated fields can be plotted as difference w.r.t. a given
;    reference dataset. It also calculates grading and taylor metrics.
;    Input data have to be regridded to a common grid in the preprocessor.
;
; Required diag_script_info attributes
;    plot_type: cycle (time), zonal (plev, lat), latlon (lat, lon) or
;               cycle_latlon (time, lat, lon)
;    time_avg: type of time average (see time_operations in
;              diag_scripts/shared/statistics.ncl)
;    region: selected region (see select_region in
;            diag_scripts/shared/latlon.ncl)
;
; Optional diag_script_info attributes
;    styleset (for cycle): as in diag_scripts/shared/plot/style.ncl functions
;    plot_stddev (for cycle): plot standard deviation
;    legend_outside (for cycle): save legend in a separate file
;    t_test (for zonal and latlon): calculate t-test in difference plots
;                                   (default: False)
;    conf_level (for zonal and latlon): confidence level for the t-test
;                                       (default: False)
;    projection: map projection for lat-lon plots (default:
;                CylindricalEquidistant)
;    plot_diff: draw difference plots (default: False)
;    calc_grading: calculate grading (default: False)
;    stippling: use stippling to mark stat. significant differences (default:
;               False = mask out non-significant differences in gray)
;    show_global_avg: diplay global avaerage as right string on lat-lon plots
;                     (default: False)
;    metric: grading metric (if calc_grading is True)
;    normalization: metric normalization (for RMSD and BIAS metrics)
;    abs_levs: (min, max, spacing) contour levels for absolute plot
;    diff_levs: (min, max, spacing) contour levels for difference plot
;    zonal_cmap (for zonal): color table (default: "amwg_blueyellowred")
;    zonal_ymin (for zonal): minimum pressure on the plots (default: 5. hPa)
;    latlon_cmap (for latlon): color table (default: "amwg_blueyellowred")
;    plot_units: plotting units (if different from standard CMOR units)
;
; Required variable_info attributes:
;    reference_dataset: reference dataset to compare with (usualy observations)
;
; Optional variable_info attributes:
;    alternative_dataset: a second dataset to compare with
;
; Caveats
;
; Modification history
;    20190405-righi_mattia: added provenance logging
;    20190315-hassler_birgit: extended to smpi metric
;    20180503-righi_mattia: completely rewritten and modularized
;    20171215-righi_mattia: merged with perfmetrics_grading and
;                           permetrics_taylor.ncl
;    20171124-righi_mattia: completely revised to adapt it to the new backend
;                           (level selection, regridding and masking now done
;                           by the python preprocessor)
;    20161220-lauer_axel: added option to set map projection for lat-lon plots
;                         (diag_script_info@projection)
;                         added option to choose how to plot t-test results:
;                         stippling or masking out in gray (lat-lon plots only)
;    20161019-lauer_axel: changed plotting of t-test results:
;                         now stippling significant grid cells (old version:
;                         masking out non-significant values in gray)
;    20160628-righi_mattia: moving ref_model specification from cfg- files to
;                        recipe file
;    20160628-senftleben_daniel: added regridding for irregular grids
;                                (ESMF_regrid).
;    20151027-lauer_axel: moved call to 'write_references' to the beginning
;                         of the code.
;    20151013-righi_mattia: fixed t-test mask in lat-lon difference plots.
;    20150325-lauer_axel: modified reference tags used for acknowledgements
;                         (projects, observations, etc.).
;    20150119-gottschaldt_klaus-dirk: removed "grid", "region" from req_atts
;                                     (for T2Ms vmrco).
;    20150113-gottschaldt_klaus-dirk: reconciled generalised regridding with
;                                     T1* & T0*
;    20140905-righi_mattia: consistent regridding and missing values mask.
;    20140701-gottschaldt_klaus-dirk: Adapted for T1M.
;    20140630-gottschaldt_klaus-dirk: Adapted for T0Ms.
;    20131203-frank_franziska: written.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/regridding.ncl"
load "$diag_scripts/shared/ensemble.ncl"
load "$diag_scripts/shared/scaling.ncl"

load "$diag_scripts/shared/plot/style.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  ; Get variables and datasets
  var0 = variable_info[0]@short_name
  info_items = select_metadata_by_name(input_file_info, var0)
  nDatasets = ListCount(info_items)

  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"plot_type", "time_avg", "region"/))

  ; Define region
  region = select_region(diag_script_info@region)

  ; Store required attributes
  ptype = diag_script_info@plot_type

  if (all(ptype.ne. \
          (/"cycle", "zonal", "latlon", "cycle_latlon", "cycle_zonal"/))) then
    error_msg("f", DIAG_SCRIPT, "", "plot_type " + ptype + " is not a " + \
              "supported plot_type in this diagnostic")
  end if

  ; Check for plot-type specific settings
  if (ptype.eq."cycle") then
    exit_if_missing_atts(diag_script_info, \
                         (/"legend_outside", "styleset", "plot_stddev"/))
  end if

  if ((ptype.eq."zonal" .or. ptype.eq."latlon") .and. \
      diag_script_info@region.ne."global") then
    error_msg("f", DIAG_SCRIPT, "", "plot_type " + ptype + \
              " implemented only for region 'global'")
  end if

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "projection", "CylindricalEquidistant")
  set_default_att(diag_script_info, "plot_diff", False)
  set_default_att(diag_script_info, "calc_grading", False)
  set_default_att(diag_script_info, "stippling", False)
  set_default_att(diag_script_info, "t_test", False)
  set_default_att(diag_script_info, "show_global_avg", False)
  set_default_att(diag_script_info, "zonal_ymin", 5.)
  set_default_att(diag_script_info, "zonal_cmap", "amwg_blueyellowred")
  set_default_att(diag_script_info, "latlon_cmap", "amwg_blueyellowred")

  ; Check consistency of diff plots settings
  if (diag_script_info@t_test .and. .not.diag_script_info@plot_diff) then
    error_msg("f", DIAG_SCRIPT, "", "plot_diff must be True to apply t-test")
  end if
  if (diag_script_info@t_test .and. .not.diag_script_info@conf_level) then
    error_msg("f", DIAG_SCRIPT, "", \
              "conf_level must be specified to apply t-test")
  end if

  ; Check metric
  if (diag_script_info@calc_grading) then
    exit_if_missing_atts(diag_script_info, (/"metric", "normalization"/))
    if (dimsizes(diag_script_info@metric).ne.\
        dimsizes(diag_script_info@normalization)) then
      error_msg("f", DIAG_SCRIPT, "", "normalization must be " + \
                "provided for each requested metric")
    end if
  end if

  ; Set dataset names
  datasetnames = metadata_att_as_array(info_items, "dataset")

  ; Save list of preproc files for provenance in collect.ncl
  preproc_files = metadata_att_as_array(info_items, "filename")

  ; Check for reference dataset definition
  if (variable_info[0]@reference_dataset.eq."None") then
    error_msg("f", DIAG_SCRIPT, "", "no reference dataset is specified")
  end if

  ; Set index of the reference (and alternative) dataset
  l_altern = False
  nobs = 1
  ref_ind = ind(datasetnames.eq.variable_info[0]@reference_dataset)
  ref_inds = ref_ind
  if (isatt(variable_info[0], "alternative_dataset")) then
    l_altern = True
    nobs = 2
    alt_ind = ind(datasetnames.eq.variable_info[0]@alternative_dataset)
    ref_inds := (/ref_ind, alt_ind/)
  end if

  ; Create output plot directory
  if (config_user_info@write_plots) then
    plot_dir = config_user_info@plot_dir
    system("mkdir -p " + plot_dir)
  end if

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Grading settings
  if (diag_script_info@calc_grading) then

    ; Define grading arrays
    nmetrics = dimsizes(diag_script_info@metric)
    ncdf_dir = new(nmetrics, string)
    nModels = dimsizes(datasetnames) - nobs

    grading = new((/nmetrics, 1, nModels, nobs/), float)
    grading!0 = "metric"
    grading!1 = "diagnostics"  ; dummy coord. to facilitate appending
    grading!2 = "models"
    grading!3 = "reference"
    grading&diagnostics = \
      variable_info[0]@diagnostic + "-" + diag_script_info@region
    grading&models = remove_index(datasetnames, ref_inds)
    if (isdim(grading, "reference")) then
      grading&reference = datasetnames(ref_inds)
    end if

    ; Special case Taylor
    if (any(diag_script_info@metric.eq."taylor")) then
      nModels = dimsizes(datasetnames) - 1  ; always 1 reference dataset
      taylor = new((/1, nModels, 2/), float)
      taylor!0 = "diagnostics"  ; dummy coord. to facilitate appending
      taylor!1 = "models"
      taylor!2 = "statistic"
      taylor&diagnostics = \
        variable_info[0]@diagnostic + "-" + diag_script_info@region
      taylor&statistic = (/"stddev_ratio", "correlation"/)
      taylor&models = remove_index(datasetnames, ref_ind)
    end if

    ; Special case SMPI
    if (any(diag_script_info@metric.eq."SMPI")) then
      nModels = dimsizes(datasetnames) - 1  ; always 1 reference model
      smpi = new((/diag_script_info@smpi_n_bootstrap + 1, nModels/), float)
      smpi!0 = "bootstrap_member"
      smpi!1 = "models"
      smpi&bootstrap_member = ispan(0, diag_script_info@smpi_n_bootstrap, 1)
      smpi&models = remove_index(datasetnames, ref_ind)
    end if

    ; Define grading filename
    do met = 0, nmetrics - 1
      ncdf_dir(met) = config_user_info@work_dir + \
        diag_script_info@metric(met) + "_" + \
        variable_info[0]@diagnostic + "-" + diag_script_info@region + ".nc"
    end do

  end if

  ; Load plot-type-specific script
  loadscript("$diag_scripts/perfmetrics/" + ptype + ".ncl")

end

begin

  ; Call plot-type-specific script
  perfmetrics_ptype_script()

  ; Finalize grading calculations
  if (diag_script_info@calc_grading) then

    do met = 0, nmetrics - 1

      if (diag_script_info@metric(met).eq."RMSD") then

        metric = grading(met, :, :, :)
        metric@title = diag_script_info@metric(met) + " metric"
        metric@long_name = \
          "Grading table of metric " + diag_script_info@metric(met)
        metric@var = "grade"

        ; Normalization
        do iobs = 0, nobs - 1
          metric(:, :, iobs) = \
            normalize_metric(metric(:, :, iobs), \
                             diag_script_info@normalization(met))
        end do

        ; Reduce dimensionality if no alternative dataset
        if (.not.l_altern) then
          metric := metric(:, :, 0)
          delete(metric@reference)
        end if

        ; Provenance information
        statistics = (/"rmsd"/)
        authors = (/"frank_franziska", "righi_mattia", "eyring_veronika"/)
        references = (/"righi15gmd", "gleckler08jgr"/)

      elseif (diag_script_info@metric(met).eq."BIAS") then

        metric = grading(met, :, :, :)
        metric@title = diag_script_info@metric(met) + " metric"
        metric@long_name = \
          "Grading table of metric " + diag_script_info@metric(met)
        metric@var = "grade"

        ; Normalization
        do iobs = 0, nobs - 1
          metric(:, :, iobs) = \
            normalize_metric(metric(:, :, iobs), \
                             diag_script_info@normalization(met))
        end do

        ; Reduce dimensionality if no alternative dataset
        if (.not.l_altern) then
          metric := metric(:, :, 0)
          delete(metric@reference)
        end if

        ; Provenance information
        statistics = (/"diff"/)
        authors = (/"frank_franziska", "righi_mattia", "eyring_veronika"/)
        references = (/"righi15gmd", "gleckler08jgr"/)

      elseif (diag_script_info@metric(met).eq."taylor") then

        metric = taylor
        metric@title = diag_script_info@metric(met) + " metric"
        metric@long_name = \
          "Grading table of metric " + diag_script_info@metric(met)
        metric@var = "grade"

        ; Provenance information
        statistics = (/"rmsd", "corr"/)
        authors = (/"frank_franziska", "righi_mattia", "eyring_veronika"/)
        references = (/"righi15gmd", "gleckler08jgr"/)

      elseif (diag_script_info@metric(met).eq."SMPI") then

        metric = smpi
        metric@title = "metrics"
        metric@long_name = "1 variable's Performance Index for " + \
          "the Single Model Performance Index"
        metric@var = "performance_index"
        metric@invar = var0
        metric@ensemble_name = diag_script_info@normalization(met)

        ; Normalization
        ens_idx = new(dimsizes(metric&models), integer)
        atts = True
        atts@project = diag_script_info@normalization(met)
        info = select_metadata_by_atts(input_file_info, atts)
        delete(atts)

        do ii = 0, dimsizes(ens_idx) - 1
          if (dimsizes(info).ne.0) then
            ens_idx(ii) = ii
          end if
        end do
        if (all(ismissing(ens_idx))) then
          error_msg("f", DIAG_SCRIPT, "", "No datasets for the selected " + \
                    "normalization (" + diag_script_info@normalization(met) + \
                    ") found")
        end if
        ens_idx := ens_idx(ind(.not.ismissing(ens_idx)))
        do iboot = 0, dimsizes(metric&bootstrap_member)-1
          metric(iboot, :) = metric(iboot, :) / avg(metric(iboot, ens_idx))
        end do

        ; Provenance information
        statistics = "smpi"
        authors = (/"gier_bettina", "hassler_birgit"/)
        references = (/"rk2008bams"/)

      else
        error_msg("f", DIAG_SCRIPT, "", "unrecognized metric " + \
                  diag_script_info@metric(met))
      end if

      ; Common attributes
      metric@metric = diag_script_info@metric(met)
      metric@diag_script = (/DIAG_SCRIPT/)
      metric@region = region@name
      metric@num_preproc_files = dimsizes(preproc_files)
      do imod = 0, metric@num_preproc_files - 1
        num_preproc = "preproc_file_" + imod
        metric@$num_preproc$ = preproc_files(imod)
      end do
      metric@ncdf_dir = ncdf_dir(met)
      ncdf_outfile = ncdf_write(metric, metric@ncdf_dir)

      ; Call provenance logger
      log_provenance(ncdf_outfile, "n/a", "n/a", statistics, \
                     diag_script_info@region, "other", authors, references, \
                     preproc_files)
      delete([/statistics, authors, references/])

      ; Write results of temporary grading list
      temp_dir = config_user_info@work_dir + "/" + \
        diag_script_info@metric(met) + ".nc"

      if (fileexists(temp_dir)) then
        temp_file = addfile(temp_dir, "r")
        temp_list = temp_file->temp_list
        temp_list := tostring(temp_list)
        temp_list := array_append_record(temp_list, ncdf_dir(met), 0)
        temp_list := tochar(temp_list)
        system("rm -f " + temp_dir)
      else
        ncdf_char = tochar(ncdf_dir(met))
        temp_list = new((/1, dimsizes(ncdf_char)/), character)
        temp_list(0, :) = ncdf_char
      end if

      ; Create new file and add list
      temp = addfile(temp_dir, "c")
      temp->temp_list = temp_list

      delete([/metric, temp_dir, temp_list, ncdf_char/])

    end do

  end if

  leave_msg(DIAG_SCRIPT, "")

end
