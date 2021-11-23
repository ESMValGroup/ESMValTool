; CLOUDS_PDF_POSTPROC
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; ############################################################################
; Description
;
; Modification history
;   20211122-lauer_axel: written.
;
; ############################################################################

begin

  inpath = "/pf/b/b380103/workesm/esmvaltool_output/" \
           + "recipe_clouds_pdf_20211105_185734/work"

  diag = "clouds_pdf_postproc.ncl"

  refname = "ESACCI-CLOUD"

  ; CMIP5
  fname5 = inpath + "/pdf_clt_djf_cmip5/southernocean/" \
    + "clouds_pdf_clt_so_cmip5.nc"
  infile5 = addfile(fname5, "r")

varnames = getfilevarnames(infile5)
print(varnames)
exit


  count5 = infile5->count
  ; CMIP6
  fname6 = inpath + "/pdf_clt_djf_cmip6/southernocean/" \
    + "clouds_pdf_clt_so_cmip6.nc"
  infile6 = addfile(fname6, "r")
  count6 = infile6->count

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

  idxref5 = ind(models5 .eq. refname)
  idxref6 = ind(models6 .eq. refname)

  ref5 = where(count5(idxref5, :, :) .gt. 1.0e-3, count5(idxref5, :, :), count5@_FillValue)
  ref6 = where(count6(idxref6, :, :) .gt. 1.0e-3, count6(idxref6, :, :), count6@_FillValue)

  c5 = dim_avg_n_Wrap(count5, 0)
  c5 = where(c5 .gt. 1.0e-3, c5, c5@_FillValue)
  c6 = dim_avg_n_Wrap(count6, 0)
  c6 = where(c6 .gt. 1.0e-3, c6, c6@_FillValue)

  diff5 = c5
  diff5 = (diff5 - ref5); / ref5 * 100.0
  diff6 = c6
  diff6 = (diff6 - ref6); / ref6 * 100.0
  delta56 = c6
  delta56 = (delta56 - c5) / c5 * 100.0

  system("rm diff.nc")
  debugfile = addfile("diff.nc", "c")
  debugfile->diff5 = diff5
  debugfile->diff6 = diff6
  debugfile->delta = delta56

end
