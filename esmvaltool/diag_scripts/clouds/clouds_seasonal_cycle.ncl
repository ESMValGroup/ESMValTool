; CLOUDS_SEASONAL_CYCLE
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; ############################################################################
; Description
;   Calculates the amplitude and month of maximum values from the
;   climatological mean seasonal cycle.
;
; Required diag_script_info attributes (diagnostic specific)
;   none
;
; Optional diag_script_info attributes (diagnostic specific)
;   colormap:        e.g., WhiteBlueGreenYellowRed, rainbow
;   epsilon:         "epsilon" value to be replaced with missing values
;   explicit_cn_levels: use these contour levels for plotting
;   filename_add:    optionally add this string to plot filesnames
;   projection:      map projection, e.g., Mollweide, Mercator
;   var:             short_name of variable to process (default = "" i.e. use
;                    first variable in variable list)
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
;   20230117-lauer_axel: added support for ICON (code from Manuel)
;   20210415-lauer_axel: written.
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
load "$diag_scripts/shared/dataset_selection.ncl"

begin
  enter_msg(DIAG_SCRIPT, "")

  set_default_att(diag_script_info, "colormap", "BkBlAqGrYeOrReViWh200")
  set_default_att(diag_script_info, "epsilon", 1.0e-4)
  set_default_att(diag_script_info, "filename_add", "")
  set_default_att(diag_script_info, "projection", "CylindricalEquidistant")
  set_default_att(diag_script_info, "var", "")

  if (diag_script_info@var .eq. "") then
    var0 = variable_info[0]@short_name
  else
    var0 = diag_script_info@var
  end if

  variables = metadata_att_as_array(variable_info, "short_name")
  if (.not. any(variables .eq. var0)) then
    errstr = "diagnostic " + diag + " requires the following variable: " + var0
    error_msg("f", DIAG_SCRIPT, "", errstr)
  end if

  var0_info = select_metadata_by_name(variable_info, var0)
  var0_info := var0_info[0]
  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info0)

  if (isatt(var0_info, "reference_dataset")) then
    refname = var0_info@reference_dataset
  end if

  names = metadata_att_as_array(info0, "dataset")
  projects = metadata_att_as_array(info0, "project")
  infiles = metadata_att_as_array(info0, "filename")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

  ref_ind = -1  ; set to invalid value

  ; if reference dataset has been defined, use it so plots can be sorted

  if (isvar("refname")) then
    ref_ind = ind(names .eq. refname)
  end if

  if (diag_script_info@filename_add .ne. "") then
    filename_add = "_" + diag_script_info@filename_add
  else
    filename_add = ""
  end if

  ; find indices of all OBS and obs4mips datasets (including "native6" ERA5)

  idxobs = get_obs(names, projects, "")

  if (idxobs(0) .eq. -1) then
    flag_multiobs = False
  else
    flag_multiobs = True
  end if

  ; find all indices of models w/o MultiModelMean/MultiModelMedian (if present)

  idxmod = get_mod(names, projects)

  if (idxmod(0) .eq. -1) then  ; no model found
    flag_multimod = False
  elseif (dimsizes(idxmod) .eq. 1) then  ; one model found
    flag_multimod = False
  else  ; more than one model found
    flag_multimod = True
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

    ; calculate climatological seasonal cycle from time series

    if (isvar("timeseries")) then
      delete(timeseries)
    end if

    timeseries = time_operations(A0, -1, -1, "average", "monthlyclim", True)

    ; calculate amplitude at each grid cell

    monmean = dim_avg_n_Wrap(timeseries, 0)
    monmean = where(abs(monmean) .le. diag_script_info@epsilon, \
                    monmean@_FillValue, monmean)

    monmin = monmean
    monmax = monmin

    do i = 0, 11
      monmin = where(timeseries(i, :, :) .lt. monmin, timeseries(i, :, :), \
                     monmin)
      monmax = where(timeseries(i, :, :) .gt. monmax, timeseries(i, :, :), \
                     monmax)
    end do

;    monind = dim_maxind(timeseries, 0)

;    data1 = 1.0 + monind
    data1 = 100.0 * abs((monmax - monmin) / monmean)
    copy_VarMeta(monmin, data1)

    ; create arrays for multi-obs and multi-model averages (if requested)

    if (ii .eq. 0) then
      multidim = dimsizes(data1)
      if (flag_multiobs) then
        newdims = array_append_record(dimsizes(idxobs), multidim, 0)
        multiobs_all = new(newdims, float)
      end if
      if (flag_multimod) then
        newdims = array_append_record(dimsizes(idxmod), multidim, 0)
        multimod_all = new(newdims, float)
      end if
    end if

    same_grid = False

    ; calculate multi-obs and multi-model averages (if requested)

    if (flag_multiobs) then
      iidx = ind(idxobs .eq. imod)
      if (.not.ismissing(iidx)) then
        dims1 = dimsizes(data1)
        dimerror = False
        if (dimsizes(multidim) .eq. dimsizes(dims1)) then
          if (any(multidim - dims1 .ne. 0)) then
            dimerror = True
          end if
        else
          dimerror = True
        end if
        if (dimerror) then
          error_msg("f", DIAG_SCRIPT, "", "dimensions of datasets " \
                    + "do not match. Use preprocessor to regrid data to " \
                    + "common grid.")
        end if
        multiobs_all(iidx, :, :) = data1
        same_grid = True
      end if
    end if

    if (flag_multimod) then
      iidx = ind(idxmod .eq. imod)
      if (.not.ismissing(iidx)) then
        dims1 = dimsizes(data1)
        dimerror = False
        if (dimsizes(multidim) .eq. dimsizes(dims1)) then
          if (any(multidim - dims1 .ne. 0)) then
            dimerror = True
          end if
        else
          dimerror = True
        end if
        if (dimerror) then
          error_msg("f", DIAG_SCRIPT, "", "dimensions of datasets " \
                    + "do not match. Use preprocessor to regrid data to " \
                    + "common grid.")
        end if
        multimod_all(iidx, :, :) = data1
        same_grid = True
      end if
    end if

    ; ###########################################
    ; # Style dependent annotation              #
    ; ###########################################
    ; retrieve unique strings describing the data
    ; function in ./diag_scripts/lib/ncl/style.ncl

    annots = project_style(info0, diag_script_info, "annots")

    ; ###########################################
    ; # plot ressources                         #
    ; ###########################################

    res = True

    res@cnFillOn       = True      ; color plot desired
    res@cnLineLabelsOn = False     ; contour lines
    res@cnLinesOn      = False

    ; colors
    ; http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml

    if (isdefined("pal")) then
      delete(pal)
    end if
    pal = read_colormap_file(diag_script_info@colormap)
;    pal = read_colormap_file("Cat12")
    ; annotation

    res@tiMainString = names(imod)
    res@tiMainFontHeightF = 0.025
    res@gsnStringFontHeightF = 0.02
    res@cnLevelSelectionMode = "ExplicitLevels"

    if (diag_script_info@projection.eq."Robinson") then
      res@mpPerimOn = False  ; turn off perimeter around map
      res@mpGridLineColor = -1
      res@mpGridAndLimbOn = True
    end if

    res@mpOutlineOn      = True
    res@mpFillOn         = False

    ; variable specific plotting settings

    res@cnLevels            = ispan(10, 200, 10)
;    res@cnLevels            = ispan(2, 12, 1)
;    res@lbLabelStrings = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", \
;                           "O", "N", "D"/)
;    res@lbLabelAlignment = "BoxCenters"

    nboxes = dimsizes(res@cnLevels)
    clen = dimsizes(pal)
    stride = max((/1, ((clen(0) - 1) - 2) / nboxes /))
    fill_colors = ispan(1, clen(0) - 1, stride)
    res@cnFillColors = pal(fill_colors, :)

    res@gsnRightString  = ""
    res@gsnLeftString   = ""
    res@gsnCenterString = ""

    gavg = area_operations(data1, -90., 90., 0., 360., "average", True)

    if (.not.ismissing(gavg)) then
      res@gsnLeftString = "mean = " + sprintf("%6.3f", gavg)
    else
      res@gsnLeftString = ""
    end if

    ; map attributes

    res@mpFillDrawOrder         = "PostDraw"    ; draw map last
    res@cnMissingValFillColor   = "Gray"
;    res@cnFillMode     = "RasterFill"       ; Raster Mode

    ; no tickmarks and no labels

    res@tmYLLabelsOn         = False
    res@tmYLOn               = False
    res@tmYRLabelsOn         = False
    res@tmYROn               = False
    res@tmXBLabelsOn         = False
    res@tmXBOn               = False
    res@tmXTLabelsOn         = False
    res@tmXTOn               = False
    res@cnInfoLabelOn        = False    ; turn off cn info label

    res@mpProjection         = diag_script_info@projection

    ; set explicit contour levels

    if (isatt(diag_script_info, "explicit_cn_levels")) then
      res@cnLevelSelectionMode = "ExplicitLevels"
      res@cnLevels = diag_script_info@explicit_cn_levels
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
    if (isatt(var0_info, "long_name")) then
      data1@long_name = var0_info@long_name
    else
      data1@long_name = var0
    end if
    data1@units    = "%"
;    data1@units    = "month"

    ; copy attributes for netCDF output

    data1@long_name = "interannual variability " + data1@long_name

    ; ###########################################
    ; # create the plot                         #
    ; ###########################################

    res@lbTitleString      = data1@units
    res@lbTitlePosition    = "Bottom"
    res@lbTitleFontHeightF = 0.02
    res@lbLabelFontHeightF = 0.02

    ; function in aux_plotting.ncl

    if (ii.eq.0) then
      wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_seasonal_cycle_" \
                    + var0 + filename_add)
    end if

    maps(ii) = gsn_csm_contour_map(wks, data1, res)

    ; #########################################
    ; # output all datasets to common netCDF  #
    ; #########################################

    nc_filename = work_dir + "clouds_seasonal_cycle_" + var0 + ".nc"
    nc_filename@existing = "append"
    data1@var = var0 + "_var_" + annots(imod)
    nc_outfile = ncdf_write(data1, nc_filename)

  end do  ; ii-loop (datasets)

  ; create panel plot

  pres                            = True    ; needed to override
                                            ; panelling defaults
  ; print dataset name on each panel
  pres@gsnPanelFigureStrings      = annots(ind_all_sorted)
  pres@gsnPanelFigureStringsFontHeightF = 0.007
  pres@lbLabelFontHeightF         = 0.01
  pres@lbAutoManage               = False
  pres@lbTopMarginF               = 0.1
  pres@lbPerimOn                  = False     ; draw line around label
                                              ; bar area
  pres@gsnPanelCenter             = False
  pres@pmLabelBarOrthogonalPosF   = -0.01     ; shift label bar a bit to
                                              ; the bottom
  outfile = panelling(wks, maps, (dim_MOD + 3) / 4, 4, pres)

  ; plot multi-obs and multi-model average (if requested)

  if (flag_multiobs) then
    multiobs = dim_avg_n(multiobs_all, 0)
    delete(multiobs_all)
    copy_VarMeta(data1, multiobs)
    gavg = area_operations(multiobs, -90., 90., 0., 360., "average", True)
    res@gsnLeftString = "mean = " + sprintf("%6.3f", gavg)
    res@tiMainString = "Multi-obs average"
    map_multiobs = gsn_csm_contour_map(wks, multiobs, res)
  end if

  if (flag_multimod) then
    multimod = dim_avg_n(multimod_all, 0)
    delete(multimod_all)
    copy_VarMeta(data1, multimod)

    mask1 = multiobs
    mask2 = multimod
    mask1 = where(.not.ismissing(mask1), 0., mask1@_FillValue)
    mask2 = where(.not.ismissing(mask2), 0., mask2@_FillValue)
    amask = mask1 + mask2
    delete(mask1)
    delete(mask2)
    refmasked = multiobs
    refmasked = refmasked + amask
    datmasked = multimod
    datmasked = datmasked + amask
    delete(amask)
    corr = calculate_metric(refmasked, datmasked, "correlation")
    gavg = area_operations(datmasked, -90., 90., 0., 360., "average", True)
    rmsd = calculate_metric(refmasked, datmasked, "RMSD")

    delete(refmasked)
    delete(datmasked)

;    gavg = area_operations(multimod, -90., 90., 0., 360., "average", True)
    res@gsnLeftString = "mean = " + sprintf("%6.3f", gavg)
    res@gsnCenterString = "corr = " + sprintf("%6.3f", corr)
    res@gsnRightString = "rmsd = " + sprintf("%6.3f", rmsd)
    res@tiMainString = "Multi-model average"
    map_multimod = gsn_csm_contour_map(wks, multimod, res)
  end if

  log_info(" Wrote " + outfile)

  ; ------------------------------------------------------------------------
  ; write provenance to common netcdf and plot file
  ; ------------------------------------------------------------------------

  statistics = (/"clim", "var"/)
  domain = "global"
  plottype = "geo"
  caption = "Relative amplitude of the climatological mean seasonal cycle " \
            + "of variable " + var0 + "."
  log_provenance(nc_outfile, outfile, caption, statistics, domain, \
                 plottype, "", "", infiles)

  leave_msg(DIAG_SCRIPT, "")

end
