; CLOUDS_LIFRAC_SCATTER_POSTPROC
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; ############################################################################
; Description
;   Read ice fraction from clouds_lifrac_scatter.ncl and plot results for
;   cmip5/cmip6 model pairs defined below.
;
; Modification history
;   20211113-lauer_axel: written.
;
; ############################################################################

begin

  inpath = "/pf/b/b380103/workesm/esmvaltool_output/" \
           + "recipe_clouds_lifrac_scatter_20211110_123241/work"

  diag = "clouds_lifrac_scatter_postproc.ncl"

  models = (/(/"bcc-csm1-1", "BCC-CSM2-MR"/), \
             (/"CESM1-WACCM", "CESM2-WACCM"/), \
             (/"GFDL-CM3", "GFDL-CM4"/), \
             (/"GISS-E2-H", "GISS-E2-1-H"/), \
             (/"HadGEM2-CC", "HadGEM3-GC31-LL"/), \
             (/"MIROC-ESM", "MIROC-ES2L"/), \
             (/"MPI-ESM-LR", "MPI-ESM1-2-LR"/), \
             (/"MRI-ESM1", "MRI-ESM2-0"/), \
             (/"NorESM1-M", "NorESM2-LM"/)/)

  refname = "CLOUDSAT-L2-CALIPSO-ICECLOUD-ERA5"

  dims = dimsizes(models)
  nmod = dims(0)

  nbins = 20

  results = new((/nmod, nbins, 2/), float)
  refdata = new(nbins, float)
  mmmedian = new((/nbins, 2/), float)
  mmp10 = new((/nbins, 2/), float)
  mmp90 = new((/nbins, 2/), float)

  ; CMIP5
  fname5 = inpath + "/lifrac_scatter_cmip5/lifrac_scatter/" \
    + "clouds_lifrac_scatter_ice_cmip5.nc"
  infile5 = addfile(fname5, "r")
  icefrac5 = infile5->icefrac
  ; CMIP6
  fname6 = inpath + "/lifrac_scatter_cmip6/lifrac_scatter/" \
    + "clouds_lifrac_scatter_ice_cmip6.nc"
  infile6 = addfile(fname6, "r")
  icefrac6 = infile6->icefrac

  results!1 = "bin"
  results&bin = icefrac5&bin

  tmp = infile5->model
  dims = dimsizes(tmp)
  idxtmp = ispan(0, dims(0) - 1, 1)
  models5 = new(dims(0), string)
  do i = 0, dims(0) - 1
    models5(i) = tostring(tmp(i, :))
    if (isStrSubset(models5(i), "Median") .or. \
        isStrSubset(models5(i), "Mean") .or. \
        isStrSubset(models5(i), "Average") .or. \
        isStrSubset(models5(i), refname)) then
      idxtmp(i) = -1
    end if
  end do
  delete(tmp)
  delete(dims)

  idx5 = ind(idxtmp .ge. 0)
  delete(idxtmp)

  tmp = infile6->model
  dims = dimsizes(tmp)
  idxtmp = ispan(0, dims(0) - 1, 1)
  models6 = new(dims(0), string)
  do i = 0, dims(0) - 1
    models6(i) = tostring(tmp(i, :))
    if (isStrSubset(models6(i), "Median") .or. \
        isStrSubset(models6(i), "Mean") .or. \
        isStrSubset(models6(i), "Average") .or. \
        isStrSubset(models6(i), refname)) then
      idxtmp(i) = -1
    end if
  end do
  delete(tmp)
  delete(dims)

  idx6 = ind(idxtmp .ge. 0)
  delete(idxtmp)

  do imod = 0, nmod - 1
    idx = ind(models5 .eq. models(imod, 0))
    if (.not. ismissing(idx)) then
      results(imod, :, 0) = icefrac5(idx, :)
    end if
    idx = ind(models6 .eq. models(imod, 1))
    if (.not. ismissing(idx)) then
      results(imod, :, 1) = icefrac6(idx, :)
    end if
    idx = ind(models6 .eq. refname)
    if (.not. ismissing(idx)) then
      refdata(:) = icefrac6(idx, :)
    end if
  end do

  ; calculate multi-model median and p10 and p90

  do n = 0, nbins - 1
    selection = icefrac5(idx5, n)
    itmp = ind(.not.ismissing(selection))
    if (.not. ismissing(itmp(0))) then
      sorted = selection(itmp)
      qsort(sorted)
      i10 = toint(dimsizes(sorted) * 0.1 + 0.5)
      i50 = toint(dimsizes(sorted) * 0.5 + 0.5)
      i90 = toint(dimsizes(sorted) * 0.9 - 0.5)
      mmp10(n, 0) = sorted(i10)
      mmmedian(n, 0) = sorted(i50)
      mmp90(n, 0) = sorted(i90)
      delete(sorted)
    else
      mmp10(n, 0) = mmp10@_FillValue
      mmedian(n, 0) = mmmedian@_FillValue
      mmp90(n, 0) = mmp90@_FillValue
    end if
    delete(selection)
    delete(itmp)

    selection = icefrac6(idx6, n)
    itmp = ind(.not.ismissing(selection))
    if (.not. ismissing(itmp(0))) then
      sorted = selection(itmp)
      qsort(sorted)
      i10 = toint(dimsizes(sorted) * 0.1 + 0.5)
      i50 = toint(dimsizes(sorted) * 0.5 + 0.5)
      i90 = toint(dimsizes(sorted) * 0.9 - 0.5)
      mmp10(n, 1) = sorted(i10)
      mmmedian(n, 1) = sorted(i50)
      mmp90(n, 1) = sorted(i90)
      delete(sorted)
    else
      mmp10(n, 1) = mmp10@_FillValue
      mmedian(n, 0) = mmmedian@_FillValue
      mmp90(n, 1) = mmp90@_FillValue
    end if
    delete(selection)
    delete(itmp)
  end do

  ; ###########################################
  ; # create the plots                        #
  ; ###########################################

  plots = new(nmod + 1, graphic)

  res = True
  res@gsnDraw        = False  ; do not draw yet
  res@gsnFrame       = False  ; don't advance frame
  res@xyMarkLineMode    = "MarkLines"
  res@xyMarkers         = (/16, 16, 16/)
  res@xyMonoMarkerColor = False
  res@xyMarkerColors    = (/"blue", "red", "black"/)
  res@xyDashPatterns    = (/0, 0, 0/)
  res@xyMarkerSizeF     = 0.01
  res@tmLabelAutoStride = True
  res@tiMainFontHeightF = 0.025
  res@tiYAxisFontHeightF = 0.025
  res@tiXAxisFontHeightF = 0.025
  res@tiXAxisString = "temperature (K)"
  res@tiYAxisString = "icefrac (%)"

  lgres                    = True
  lgres@vpWidthF           = 0.2                   ; width of legend (NDC)
  lgres@vpHeightF          = 0.08                   ; height of legend (NDC)
  lgres@lgAutoManage       = False                 ; Necessary to set font hgt
  lgres@lgMonoLineColor    = True
  lgMonoMarkerIndex        = False
  lgMonoMarkerColor        = False
  lgres@lgLineColors       = "black"
  lgres@lgMonoItemType     = True                  ; more than one type
  lgres@lgItemType         = res@xyMarkLineMode    ; line/marker
  lgres@lgDashIndexes      = res@xyDashPatterns    ; dash indexes
  lgres@lgMarkerIndexes    = res@xyMarkers         ; marker indexes
  lgres@lgMarkerColors     = res@xyMarkerColors    ; marker colors
  lgres@lgLabelFontHeightF = 0.12
  lgres@lgLineLabelsOn     = False
  lgres@lgMarkerSizeF      = 0.01

  wks = gsn_open_wks("pdf", "clouds_lifrac_scatter_post_ice")

  data = new((/3, nbins/), float)

  do ii = 0, nmod
    if (ii .lt. nmod) then
      data(0, :) = results(ii, :, 0)
      data(1, :) = results(ii, :, 1)
      data(2, :) = refdata(:)
      labels = (/models(ii, 0), models(ii, 1), "observations"/)
    else
      data(0, :) = (/mmmedian(:, 0)/)
      data(1, :) = (/mmmedian(:, 1)/)
      data(2, :) = refdata(:)
      labels = (/"CMIP5 MMMedian", "CMIP6 MMMedian", "observations"/)
    end if

    plots(ii) = gsn_csm_xy(wks, results&bin, data, res)

    gsn_legend_ndc(wks, 3, labels, 0.58, 0.78, lgres)

    if (ii .eq. nmod) then
      res2 = True
      res2@gsnDraw         = False   ; do not draw yet
      res2@gsnFrame        = False  ; don't advance frame
      res2@xyLineColor     = -1  ; Make lines transparent
      mmvar = new((/2, nbins/), float)
      mmvar(0, :) = mmp10(:, 0)
      mmvar(1, :) = mmp90(:, 0)
      res2@gsnXYFillColors = (/1.0, 0.9, 0.9, 1.0/)  ; "lightpink"
      plotvar5 = gsn_csm_xy(wks, results&bin, mmvar, res2)
      overlay(plots(ii), plotvar5)
      mmvar(0, :) = mmp10(:, 1)
      mmvar(1, :) = mmp90(:, 1)
      res2@gsnXYFillColors = (/0.9, 0.9, 1.0, 0.3/)  ; "lightpink"
      plotvar6 = gsn_csm_xy(wks, results&bin, mmvar, res2)
      overlay(plots(ii), plotvar6)
      delete(mmvar)
    end if

    draw(plots(ii))
    frame(wks)
  end do

end
