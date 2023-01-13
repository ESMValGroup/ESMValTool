; CLOUDS_LIFRAC_SCATTER_POSTPROC
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; ############################################################################
; Description
;   Read ice fraction from clouds_lifrac_scatter.ncl and plot results for
;   cmip5/cmip6 model pairs defined below.
;
; Required diag_script_info attributes (diagnostic specific)
;   models: array of CMIP5/CMIP6 model pairs to be compared
;   refname: name of reference dataset
;
; Optional diag_script_info attributes (diagnostic specific)
;   nbins: number of bins used by clouds_lifrac_scatter.ncl (default = 20)
;   reg: region (string) (default = "")
;   t_int: array of temperatures for printing additional diagnostics
;
; Required variable attributes (variable specific)
;   none
;
; Optional variable_info attributes (variable specific)
;   none
;
; Caveats
;   none
;
; Modification history
;   20211113-lauer_axel: written.
;
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT)
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  set_default_att(diag_script_info, "nbins", 20)
  set_default_att(diag_script_info, "reg", "")
  set_default_att(diag_script_info, "t_int", (/240.0, 250.0, 260.0, 270.0/))

  models = diag_script_info@models
  t_int = diag_script_info@t_int

  if (isatt(diag_script_info, "refname")) then
    refname = diag_script_info@refname
  else
    error_msg("f", DIAG_SCRIPT, "", "no reference dataset defined (refname).")
  end if

  dims = dimsizes(models)
  nmod = dims(0)

  nbins = diag_script_info@nbins
  reg = diag_script_info@reg

  test = str_lower(diag_script_info@input_files)

  i = dimsizes(test)
  if (i .ne. 2) then
    error_msg("f", DIAG_SCRIPT, "", "number of datasets (ancestors) is " \
              + tostring(i) + ", expected 2")
  end if

  if (reg .ne. "") then
    reg = "_" + reg
  end if

  fname5 = ""
  fname6 = ""

  do i = 0, dimsizes(test) - 1
    if (isStrSubset(test(i), "cmip5")) then
      fname5 = diag_script_info@input_files(i) + \
               "/clouds_lifrac_scatter_ice_cmip5" + reg + ".nc"
    end if
    if (isStrSubset(test(i), "cmip6")) then
      fname6 = diag_script_info@input_files(i) + \
               "/clouds_lifrac_scatter_ice_cmip6" + reg + ".nc"
    end if
  end do

  if (fname5 .eq. "") then
    error_msg("f", DIAG_SCRIPT, "", "no dataset (ancestor) for CMIP5 found.")
  end if
  if (fname6 .eq. "") then
    error_msg("f", DIAG_SCRIPT, "", "no dataset (ancestor) for CMIP6 found.")
  end if

  results = new((/nmod, nbins, 2/), float)
  refdata = new(nbins, float)
  mmmedian = new((/nbins, 2/), float)
  mmmean = new((/nbins, 2/), float)
  mmp10 = new((/nbins, 2/), float)
  mmp90 = new((/nbins, 2/), float)
  mmstd = new((/nbins, 2/), float)

  ; CMIP5
  infile5 = addfile(fname5, "r")
  icefrac5 = infile5->icefrac
  ; CMIP6
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

  model_ok = new((/nmod, 2/), logical)
  model_ok = False
  ref_ok = False

  do imod = 0, nmod - 1
    idx = ind(models5 .eq. models(imod, 0))
    if (.not. ismissing(idx)) then
      results(imod, :, 0) = icefrac5(idx, :)
      model_ok(imod, 0) = True
    end if
    idx = ind(models6 .eq. models(imod, 1))
    if (.not. ismissing(idx)) then
      results(imod, :, 1) = icefrac6(idx, :)
      model_ok(imod, 1) = True
    end if
    idx = ind(models6 .eq. refname)
    if (.not. ismissing(idx)) then
      refdata(:) = icefrac6(idx, :)
      ref_ok = True
    end if
  end do

  if (any(.not. model_ok(:, 0))) then
    inotok = ind(.not. model_ok(:, 0))
    log_info("warning: could not find data for the following CMIP5 model(s) " \
             + str_join(models(inotok, 0), ","))
  end if
  if (any(.not. model_ok(:, 1))) then
    inotok = ind(.not. model_ok(:, 1))
    log_info("warning: could not find data for the following CMIP6 model(s) " \
             + str_join(models(inotok, 1), ","))
  end if

  ; calculate multi-model, p10 and p90, mean and standard deviation

  print("T (K), CMIP5-CMIP6 (abs), CMIP5-CMIP6 (rel), CMIP5, CMIP6, ref")
  print("--------------------------------------------------------------")

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
      mmmean(n, 0) = avg(selection)
      mmstd(n, 0) = stddev(selection)
    else
      mmp10(n, 0) = mmp10@_FillValue
      mmedian(n, 0) = mmmedian@_FillValue
      mmp90(n, 0) = mmp90@_FillValue
      mmmean(n, 0) = mmmean@_FillValue
      mmstd(n, 0) = mmstd@_FillValue
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
      mmmean(n, 1) = avg(selection)
      mmstd(n, 1) = stddev(selection)
    else
      mmp10(n, 1) = mmp10@_FillValue
      mmedian(n, 0) = mmmedian@_FillValue
      mmp90(n, 1) = mmp90@_FillValue
      mmmean(n, 1) = mmmean@_FillValue
      mmstd(n, 1) = mmstd@_FillValue
    end if

    delete(selection)
    delete(itmp)

    delta = mmmean(n, 0) - mmmean(n, 1)

    if (n .ge. 1) then
      do it = 0, dimsizes(t_int) - 1
        if ((results&bin(n - 1) .lt. t_int(it)) .and. \
            (results&bin(n) .gt. t_int(it))) then
          dx = results&bin(n) - results&bin(n - 1)
          dy5 = mmmean(n, 0) - mmmean(n - 1, 0)
          dy6 = mmmean(n, 1) - mmmean(n - 1, 1)
          dyr = refdata(n) - refdata(n - 1)
          r5 = mmmean(n - 1, 0) + dy5 / dx * (t_int(it) - results&bin(n - 1))
          r6 = mmmean(n - 1, 1) + dy6 / dx * (t_int(it) - results&bin(n - 1))
          rr = refdata(n - 1) + dyr / dx * (t_int(it) - results&bin(n - 1))

          print(tostring(t_int(it)) + ", " + tostring(r5) + ", " \
                + tostring(r6) + ", " + tostring(rr))
        end if
      end do
    end if

    mmmean = where(mmmean .eq. 0, mmmean@_FillValue, mmmean)

    print(tostring(results&bin(n)) + ", " + tostring(delta) + " (" \
          + tostring(delta / mmmean(n, 0) * 100.0 + "%)") \
          + ", " + tostring(mmmean(n, 0)) + ", " + tostring(mmmean(n, 1)) \
          + ", " + tostring(refdata(n)))
  end do

  ; ###########################################
  ; # create the plots                        #
  ; ###########################################

  plots = new(nmod + 2, graphic)

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

  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, \
                "clouds_lifrac_scatter_post_ice" + reg)

  data = new((/3, nbins/), float)

  do ii = 0, nmod + 1
    if (ii .lt. nmod) then
      data(0, :) = results(ii, :, 0)
      data(1, :) = results(ii, :, 1)
      data(2, :) = refdata(:)
      labels = (/models(ii, 0), models(ii, 1), "observations"/)
    else if (ii .eq. nmod) then
      data(0, :) = (/mmmedian(:, 0)/)
      data(1, :) = (/mmmedian(:, 1)/)
      data(2, :) = refdata(:)
      labels = (/"CMIP5 MMMedian", "CMIP6 MMMedian", "observations"/)
    else if (ii .eq. nmod + 1) then
      data(0, :) = (/mmmean(:, 0)/)
      data(1, :) = (/mmmean(:, 1)/)
      data(2, :) = refdata(:)
      labels = (/"CMIP5 MMMean", "CMIP6 MMMean", "observations"/)
    end if
    end if
    end if

    plots(ii) = gsn_csm_xy(wks, results&bin, data, res)

    gsn_legend_ndc(wks, 3, labels, 0.58, 0.78, lgres)

    if (ii .ge. nmod) then
      res2 = True
      res2@gsnDraw         = False   ; do not draw yet
      res2@gsnFrame        = False  ; don't advance frame
      res2@xyLineColor     = -1  ; Make lines transparent
      mmvar = new((/2, nbins/), float)
      if (ii .eq. nmod) then  ; median + percentiles
        mmvar(0, :) = mmp10(:, 0)
        mmvar(1, :) = mmp90(:, 0)
      else if (ii .eq. nmod + 1) then  ; mean + stddev
        mmvar(0, :) = mmmean(:, 0) - mmstd(:, 0)
        mmvar(1, :) = mmmean(:, 0) + mmstd(:, 0)
      end if
      end if
      res2@gsnXYFillColors = (/0.9, 0.9, 1.0, 1.0/)  ; "lightblue"
      plotvar5 = gsn_csm_xy(wks, results&bin, mmvar, res2)
      overlay(plots(ii), plotvar5)
      if (ii .eq. nmod) then  ; median + percentiles
        mmvar(0, :) = mmp10(:, 1)
        mmvar(1, :) = mmp90(:, 1)
      else if (ii .eq. nmod + 1) then
        mmvar(0, :) = mmmean(:, 1) - mmstd(:, 1)
        mmvar(1, :) = mmmean(:, 1) + mmstd(:, 1)
      end if
      end if
      res2@gsnXYFillColors = (/1.0, 0.9, 0.9, 0.3/)  ; "lightpink"
      plotvar6 = gsn_csm_xy(wks, results&bin, mmvar, res2)
      overlay(plots(ii), plotvar6)
      delete(mmvar)
    end if

    draw(plots(ii))
    frame(wks)
  end do

end
