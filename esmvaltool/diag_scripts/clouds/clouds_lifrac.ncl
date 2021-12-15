; CLOUDS_LIFRAC
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; ############################################################################
; Description
;   Calculates zonal average liquid / ice fraction as a function of height.
;
; Required diag_script_info attributes (diagnostic specific)
;   none
;
; Optional diag_script_info attributes (diagnostic specific)
;   embracesetup:       True = 2 plots per line, False = 4 plots per line
;                       (default)
;   extralegend:        create additional plot with legend only
;   filename_add:       optionally add this string to plot filesnames
;   min_mass:           minimum cloud condensate (same units as clw, cli)
;   panel_labels:       label individual panels (true, false)
;   PanelTop:           manual override for "@gnsPanelTop" used by panel
;                       plot(s)
;
; Required variable attributes (variable specific)
;   none
;
; Optional variable_info attributes (variable specific)
;   long_name:         variable description
;   reference_dataset: reference dataset
;   units:             variable units (for labeling plot only)
;
; Caveats
;   none
;
; Modification history
;   20200327-lauer_axel: written.
;
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  diag = "clouds_lifrac.ncl"
  variables = get_unique_values(metadata_att_as_array(variable_info, \
                                                      "short_name"))
  dim_VAR = dimsizes(variables)

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + variables + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; make sure required variables are available

  varlist = (/"clw", "cli"/)
  idx = new(dimsizes(varlist), integer)

  do i = 0, dimsizes(varlist) - 1
    idx(i) = ind(variables .eq. varlist(i))
  end do

  if (any(ismissing(idx))) then
    errstr = "diagnostic " + diag + " requires the following variable(s): " \
             + str_join(varlist, ", ")
    error_msg("f", DIAG_SCRIPT, "", errstr)
  end if

  delete(varlist)
  delete(idx)

  ; make sure that only the variables "clw" and "cli" are used and that their
  ; order is correct
  variables := (/"clw", "cli"/)
  dim_VAR = dimsizes(variables)
  refname = new(dim_VAR, string)

  ; save input files for writing provenance

  infiles = metadata_att_as_array(input_file_info, "filename")

  ; get reference datasets (if present) and check that number of datasets
  ; is equal for each variable

  do i = 0, dim_VAR - 1
    var = variables(i)
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

  ; Set default values for non-required diag_script_info attributes

  set_default_att(diag_script_info, "embrace_setup", False)
  set_default_att(diag_script_info, "extralegend", False)
  set_default_att(diag_script_info, "filename_add", "")
  set_default_att(diag_script_info, "panel_labels", True)
  set_default_att(diag_script_info, "min_mass", 1.0e-15)

  if (diag_script_info@filename_add .ne. "") then
    filename_add = "_" + diag_script_info@filename_add
  else
    filename_add = ""
  end if

  embracesetup = diag_script_info@embrace_setup
  panel_labels = diag_script_info@panel_labels
  extralegend = diag_script_info@extralegend
  min_mass = diag_script_info@min_mass

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

end

begin
  ; ############
  ; # get data #
  ; ############

  maps_ice = new(dim_MOD, graphic)
  maps_liq = new(dim_MOD, graphic)

  info_clw = select_metadata_by_name(input_file_info, "clw")
  names_clw = metadata_att_as_array(info_clw, "dataset")
  projects_clw = metadata_att_as_array(info_clw, "project")
  info_cli = select_metadata_by_name(input_file_info, "cli")
  names_cli = metadata_att_as_array(info_cli, "dataset")
  projects_cli = metadata_att_as_array(info_cli, "project")

  refidx_clw = ind(names_clw .eq. refname(0))
  refidx_cli = ind(names_cli .eq. refname(1))

  if (ismissing(refidx_clw) .or. ismissing(refidx_cli)) then
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

    error_msg("f", DIAG_SCRIPT, "", "the following reference dataset(s) " \
              + "are not available: " + str)
  end if

  do ii = 0, dim_MOD - 1
    atts_cli = True
    atts_cli@short_name = "cli"

    atts_clw = True
    atts_clw@short_name = "clw"

    ; reference datasets may have different names
    if (ii .eq. refidx_clw) then
      atts_cli@dataset = refname(1)
      atts_clw@dataset = refname(0)
    ; all other datasets: force same dataset name for clw and cli
    else
      atts_cli@dataset = names_clw(ii)
      atts_clw@dataset = names_clw(ii)
    end if

    ; read cli

    info = select_metadata_by_atts(input_file_info, atts_cli)
    cli = read_data(info[0])
    delete(info)

    ; read clw

    info = select_metadata_by_atts(input_file_info, atts_clw)
    clw = read_data(info[0])
    delete(info)

    ; check dimensions

    dims_cli = dimsizes(cli)
    dims_clw = dimsizes(clw)

    dimerror = False

    if (dimsizes(dims_cli) .eq. dimsizes(dims_clw)) then
      if (any(dims_cli - dims_clw .ne. 0)) then
        dimerror = True
      end if
    else
      dimerror = True
    end if

    if (dimerror) then
      error_msg("f", DIAG_SCRIPT, "", "dimensions of datasets " \
                + atts_cli@dataset + " (variable cli) and " \
                + atts_clw@dataset + " (variable clw) do not match.")
    end if

    ; check dimensions

    if (dimsizes(dims_cli) .lt. 2 .or. dimsizes(dims_clw) .lt. 2) then
      error_msg("f", DIAG_SCRIPT, "", "all variables need to have 3 " + \
                " dimensions")
    end if

    do i = 0, dim_VAR - 1
      if (variables(i) .eq. "cli") then
        dims = getvardims(cli)
        vcoordi = dims(1)  ; save name of vertical coordinate variable
        if (ii .eq. 0) then
          vcoord0_var = cli&$vcoordi$
          vcoord0 = vcoordi
        end if
      else
        dims = getvardims(clw)
        vcoordl = dims(1)  ; save name of vertical coordinate variable
      end if
      idx = ind(dims .eq. "lat")
      if (ismissing(idx)) then
        error_msg("f", DIAG_SCRIPT, "", variables(i) + ": no lat dimension")
      end if
      idx = ind(dims .eq. "time")
      if (ismissing(idx)) then
        error_msg("f", DIAG_SCRIPT, "", variables(i) + ": no time dimension")
      end if
      delete(dims)
    end do

    delete(dims_cli)
    delete(dims_clw)

    ; average over time

    ice = time_operations(cli, -1, -1, "average", "annualclim", True)
    delete(cli)
    liq = time_operations(clw, -1, -1, "average", "annualclim", True)
    delete(clw)

    ; calculate ice fraction

    mass = ice + liq
    mass = where(mass .lt. min_mass, mass@_FillValue, mass)

    icefrac = 100.0 * ice / mass
    liqfrac = 100.0 * liq / mass

    copy_VarMeta(ice, icefrac)
    copy_VarMeta(liq, liqfrac)

    delete(ice)
    delete(liq)
    delete(mass)

    icefrac@units = "%"
    liqfrac@units = "%"
    icefrac@long_name = "cloud ice fraction"
    liqfrac@long_name = "cloud liquid fraction"
    icefrac@var = "icefrac"
    liqfrac@var = "liqfrac"
    ref_ind = refidx_cli
    if (ismissing(ref_ind)) then
      ref_ind = -1
    end if
    imod = ii
    names = names_cli
    projects = projects_cli

    ; if reference datasets for clw and cli are from different sources
    if (refname(0) .ne. refname(1)) then
      names(refidx_cli) = refname(0) + "/" + refname(1)
    end if

    ; ###########################################
    ; # Style dependent annotation              #
    ; ###########################################
    ; retrieve unique strings describing the data
    ; function in ./diag_scripts/shared/plot/style.ncl

    ; ###########################################
    ; # plot ressources                         #
    ; ###########################################

    res = True

    res@cnFillOn       = True      ; color plot desired
    res@cnLineLabelsOn = False     ; contour lines

    ; colors
    ; http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml

    ; annotation

;    res@tiMainOn             = False
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLinesOn            = False
;    res@lbLabelBarOn         = False
    res@gsnRightString       = ""
    res@cnMissingValFillColor = "Gray"
    res@cnInfoLabelOn      = False    ; turn off cn info label

    nc_filename_ice = work_dir + "clouds_" + icefrac@var + "_mean.nc"
    nc_filename_liq = work_dir + "clouds_" + liqfrac@var + "_mean.nc"
    nc_filename_ice@existing = "append"
    nc_filename_liq@existing = "append"

    if (isatt(icefrac&$vcoordi$, "standard_name")) then
      res@tiYAxisString = icefrac&$vcoordi$@standard_name
      if (isatt(icefrac&$vcoordi$, "units")) then
          res@tiYAxisString = res@tiYAxisString + " (" \
            + icefrac&$vcoordi$@units + ")"
      end if
    end if

    if (isatt(icefrac&lat, "standard_name")) then
      res@tiXAxisString = icefrac&lat@standard_name
    end if

    ; set explicit contour levels

    if (isatt(diag_script_info, "explicit_cn_levels")) then
      res@cnLevelSelectionMode = "ExplicitLevels"
      res@cnLevels = diag_script_info@explicit_cn_levels
    end if

    if (.not.isatt(res, "cnLevels")) then
      res@cnLevels = fspan(5, 95, 10)
    end if

    ; ###########################################
    ; # other Metadata: diag_script, var        #
    ; ###########################################

    icefrac@diag_script = (/DIAG_SCRIPT/)
    icefrac@var_long_name = icefrac@long_name
    icefrac@var_units = icefrac@units

    liqfrac@diag_script = (/DIAG_SCRIPT/)
    liqfrac@var_long_name = liqfrac@long_name
    liqfrac@var_units = liqfrac@units

    res@tiMainString             = names(ii)

    ; ###########################################
    ; # create the plots                        #
    ; ###########################################

;    res@gsnDraw        = False  ; do not draw yet
;    res@gsnFrame       = False  ; don't advance frame

    ; function in aux_plotting.ncl

    if (ii.eq.0) then
      ; note: an array of workspaces (i.e. wks(numseas)) does not work as
      ;       attributes cannot be assigned to each array element
      ;       individually
      wksice = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_zonal_" + \
                       icefrac@var + filename_add)
      wksliq = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_zonal_" + \
                       liqfrac@var + filename_add)
    end if

    maps_ice(imod) = gsn_csm_pres_hgt(wksice, icefrac, res)
    maps_liq(imod) = gsn_csm_pres_hgt(wksliq, liqfrac, res)

    ; mandatory netcdf output

    ; *** icefrac ***

    data1 = icefrac
    data1@var = str_sub_str(names(imod), "/", "-")

    ; vertical coordinate might have a different name, which prevents
    ; all data to be written to the same netcdf file
    ; --> rename vertical coordinate to match the one of the first dataset
    ;     written to the netcdf file
    delete(data1&$vcoordi$)
    data1!0 = vcoord0
    data1&$vcoord0$ = vcoord0_var

    nc_outfile_ice = ncdf_write(data1, nc_filename_ice)

    delete(data1)

    ; *** liqfrac ***

    data1 = liqfrac
    data1@var = str_sub_str(names(imod), "/", "-")

    ; vertical coordinate might have a different name, which prevents
    ; all data to be written to the same netcdf file
    ; --> rename vertical coordinate to match the one of the first dataset
    ;     written to the netcdf file
    delete(data1&$vcoordl$)
    data1!0 = vcoord0
    data1&$vcoord0$ = vcoord0_var

    nc_outfile_liq = ncdf_write(data1, nc_filename_liq)

    delete(data1)

    ; =====================================================================

  end do  ; ii-loop (models)

  ; save default color map in case it is needed later for optionally
  ; plotting color bar to a separate file

  tmp_colors = gsn_retrieve_colormap(wksice)
  cdims = dimsizes(tmp_colors)
  nboxes = dimsizes(res@cnLevels)
  clen = cdims(0)
  stride = max((/1, ((clen(0)-1) - 2) / nboxes /))
  fill_colors = ispan(2, clen(0) - 1, stride)
  mean_colors = tmp_colors(fill_colors, :)
  delete(tmp_colors)
  delete(fill_colors)
  delete(cdims)

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
;    pres@gsnPanelLabelBar     = True    ; add common colorbar
  if (panel_labels) then
    ; print dataset name on each panel
    pres@gsnPanelFigureStrings = names(plotind)
  end if
  pres@gsnPanelFigureStringsFontHeightF = min((/0.008, 0.008 * 6.0 \
                                              / tofloat((dim_MOD + 1) / 2)/))
  pres@lbLabelFontHeightF               = min((/0.01, 0.01 * 6.0 \
                                              / tofloat((dim_MOD + 1) / 2)/))
  pres@lbAutoManage                     = False
  pres@lbTopMarginF                     = 0.1
  pres@lbTitleOn                        = True
  pres@lbTitleFontHeightF               = min((/0.015, 0.01 * 6.0 \
                                              / tofloat((dim_MOD + 1) / 2)/))
  pres@lbTitlePosition                  = "Bottom"
  pres@gsnPanelCenter                   = False
  if (dim_MOD.le.8) then
    pres@pmLabelBarOrthogonalPosF       = -0.03
  else
    pres@pmLabelBarOrthogonalPosF       = -0.01   ; shift label bar a bit to
                                                  ; the bottom
  end if

  if (embracesetup) then
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

    plotsperline = new((dim_MOD + 1) / 2, integer)
    plotsperline = 2

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

    pres_ice = pres
    pres_liq = pres
    pres_ice@lbTitleString = icefrac@long_name + " (" + icefrac@units + ")"
    pres_liq@lbTitleString = liqfrac@long_name + " (" + liqfrac@units + ")"

    gsn_panel(wksice, maps_ice(plotind), plotsperline, pres_ice)
    gsn_panel(wksliq, maps_liq(plotind), plotsperline, pres_liq)
    outfile_ice = wksice@fullname
    outfile_liq = wksliq@fullname

  else  ; if embracesetup
    pres_ice = pres
    pres_liq = pres
    pres_ice@lbTitleString = icefrac@long_name + " (" + icefrac@units + ")"
    pres_liq@lbTitleString = liqfrac@long_name + " (" + liqfrac@units + ")"

    outfile_ice = panelling(wksice, maps_ice(plotind), (dim_MOD + 3) / 4, \
                            4, pres_ice)
    outfile_liq = panelling(wksliq, maps_liq(plotind), (dim_MOD + 3) / 4, \
                            4, pres_liq)
  end if  ; if embracesetup

  log_info("Wrote " + outfile_ice + ", " + outfile_liq)

  ; ----------------------------------------------------------------------
  ; write provenance to netcdf output and plot file(s) (mean)
  ; ----------------------------------------------------------------------

  statistics = (/"clim", "mean"/)
  domain = "global"
  plottype = "zonal"

  caption = "Zonal mean values for variable " + icefrac@var + "."
  log_provenance(nc_outfile_ice, outfile_ice, caption, statistics, \
                 domain, plottype, "", "", infiles)

  caption = "Zonal mean values for variable " + liqfrac@var + "."
  log_provenance(nc_outfile_liq, outfile_liq, caption, statistics, \
                 domain, plottype, "", "", infiles)

  ; ======================================================================

  ; optionally save legends to extra files

  if (extralegend) then
    nboxes = dimsizes(res@cnLevels) + 1
    wksleg = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_lifrac_legend")
    pres@lbMonoFillPattern = True
    pres@lbOrientation = "Horizontal"
    pres@vpWidthF = 0.7
    pres@vpHeightF = 0.1
    pres@lbLabelFontHeightF = 0.015
    pres@lbLabelAlignment = "InteriorEdges"
    pres@lbTitleFontHeightF = 0.015
    pres@lbTitleString = icefrac@long_name + " (" + icefrac@units + ")"

    labels = tostring(res@cnLevels)

    ; remove trailing zeros from strings

    do i = 0, dimsizes(labels) - 1
      i1  = str_index_of_substr(labels(i), ".", -1)
      if (.not.ismissing(i1)) then
        tmp = stringtochar(labels(i))
        do j = dimsizes(tmp) - 2, i1, 1
          if ((tmp(j).ne.".").and.(tmp(j).ne."0")) then
            break
          end if
        end do
        labels(i) = chartostring(tmp(0:j))
        delete(tmp)
      end if
    end do

    if (isatt(icefrac, "res_cnFillColors")) then
      pres@lbFillColors = res@cnFillColors
    else if (isatt(icefrac, "res_cnFillPalette")) then
      pres@lbFillColors = res@cnFillPalette
    else
      pres@lbFillColors = mean_colors  ; default colors
    end if
    end if

    gsn_labelbar_ndc(wksleg, nboxes, labels, 0.1, 0.9, pres)

    delete(wksleg)
    delete(labels)
    delete(pres@lbFillColors)

  end if  ; if (extralegend)

  ; ==========================================================================

  leave_msg(DIAG_SCRIPT, "")

end
