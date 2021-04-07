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
;   filename_add:    optionally add this string to plot filesnames
;   projection:      map projection, e.g., Mollweide, Mercator
;   var:             short_name of variable to process (default = "" - use
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
;   20210407-lauer_axel: added option to speficfy variable if more than one
;                        variable is present
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

  set_default_att(diag_script_info, "colormap", "WhiteBlueGreenYellowRed")
  set_default_att(diag_script_info, "filename_add", "")
  set_default_att(diag_script_info, "projection", "CylindricalEquidistant")
  set_default_att(diag_script_info, "var", "")

  if (diag_script_info@var .eq. "") then
    var0 = variable_info[0]@short_name
  else
    var0 = diag_script_info@var
  end if

  variables = metadata_att_as_array(variable_info, "short_name")
  varidx = ind(variables .eq. var0)
  if (ismissing(varidx)) then
    errstr = "diagnostic " + diag + " requires the following variable: var0"
    error_msg("f", DIAG_SCRIPT, "", errstr)
  end if

  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info0)

  if (isatt(variable_info[varidx], "reference_dataset")) then
    refname = variable_info[varidx]@reference_dataset
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

  if (diag_script_info@filename_add .ne. "") then
    filename_add = "_" + diag_script_info@filename_add
  else
    filename_add = ""
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

    ; annotation

    res@tiMainString = names(imod)
    res@gsnLeftStringFontHeightF = 0.015
    res@cnLevelSelectionMode = "ExplicitLevels"

    if (diag_script_info@projection.eq."Robinson") then
      res@mpPerimOn = False  ; turn off perimeter around map
      res@mpGridLineColor = -1
      res@mpGridAndLimbOn = True
    end if

    res@mpOutlineOn      = True
    res@mpFillOn         = False

    ; variable specific plotting settings

    if (any((/"clt"/).eq.var0)) then
      res@cnLevels            = ispan(5, 50, 5)
    else
      res@cnLevels            = ispan(5, 100, 5)
    end if

    nboxes = dimsizes(res@cnLevels)
    clen = dimsizes(pal)
    stride = max((/1, ((clen(0)-1) - 2) / nboxes /))
    fill_colors = ispan(2, clen(0) - 1, stride)
    res@cnFillColors = fill_colors

    res@gsnRightString       = ""

    gavg = area_operations(data1, -90., 90., 0., 360., "average", True)
    if (.not.ismissing(gavg)) then
      res@gsnLeftString = "mean = " + sprintf("%6.3f", gavg)
    else
      res@gsnLeftString = ""
    end if

    ; map attributes

    res@mpFillDrawOrder         = "PostDraw"    ; draw map last
    res@cnMissingValFillColor   = "Gray"

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
    if (isatt(variable_info[varidx], "long_name")) then
      data1@long_name = variable_info[varidx]@long_name
    else
      data1@long_name = var0
    end if
    data1@units    = "%"

    ; copy attributes for netCDF output

    data1@long_name = "interannual variability " + data1@long_name

    ; ###########################################
    ; # create the plot                         #
    ; ###########################################

    res@lbTitleString      = data1@units
    res@lbTitlePosition    = "Bottom"
    res@lbTitleFontHeightF = 0.015

    ; function in aux_plotting.ncl

    if (ii.eq.0) then
      wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_interannual_" \
                    + var0 + filename_add)
    end if

    maps(ii) = gsn_csm_contour_map(wks, data1, res)

    ; #########################################
    ; # output all datasets to common netCDF  #
    ; #########################################

    nc_filename = work_dir + "clouds_interannual_" + var0 + ".nc"
    nc_filename@existing = "append"
    data1@var = var0 + "_var_" + annots(imod)
    nc_outfile = ncdf_write(data1, nc_filename)

  end do  ; ii-loop (datasets)

  if (write_plots) then
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
    log_info(" Wrote " + outfile)
  else
    outfile = ""
  end if  ; if write_plots

  ; ------------------------------------------------------------------------
  ; write provenance to common netcdf and plot file
  ; ------------------------------------------------------------------------

  statistics = (/"clim", "var"/)
  domain = "global"
  plottype = "geo"
  caption = "Interannual variability of variable " + var0 + "."
  log_provenance(nc_outfile, outfile, caption, statistics, domain, \
                 plottype, "", "", infiles)

  leave_msg(DIAG_SCRIPT, "")

end
