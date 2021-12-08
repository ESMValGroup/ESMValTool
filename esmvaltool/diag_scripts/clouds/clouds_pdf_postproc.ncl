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

;  reg = "so"
;  reg = "sep"
  reg = "itcz"
;  season = "djf"
  season = "ann"

  ; CMIP5
  fname5 = inpath + "/pdf_clt_" + season + "_cmip5/" + reg + "/" \
    + "clouds_pdf_clt_" + reg + "_cmip5.nc"
  infile5 = addfile(fname5, "r")
  varnames5 = getfilevarnames(infile5)
  dims = dimsizes(varnames5)
  idxtmp = ispan(0, dims(0) - 1, 1)
  do i = 0, dims(0) - 1
    if (isStrSubset(varnames5(i), "_pdf_")) then
      if (isStrSubset(varnames5(i), refname)) then
        idxtmp(i) = -1
        refidx5 = i
      end if
    else
      idxtmp(i) = -1
    end if
  end do
  delete(dims)

  idx5 = ind(idxtmp .ge. 0)
  delete(idxtmp)

  ; CMIP6
  fname6 = inpath + "/pdf_clt_" + season + "_cmip6/" + reg + "/" \
    + "clouds_pdf_clt_" + reg + "_cmip6.nc"
  infile6 = addfile(fname6, "r")
  varnames6 = getfilevarnames(infile6)
  dims = dimsizes(varnames6)
  idxtmp = ispan(0, dims(0) - 1, 1)
  do i = 0, dims(0) - 1
    if (isStrSubset(varnames6(i), "_pdf_")) then
      if (isStrSubset(varnames6(i), refname)) then
        idxtmp(i) = -1
        refidx6 = i
      end if
    else
      idxtmp(i) = -1
    end if
  end do
  delete(dims)

  idx6 = ind(idxtmp .ge. 0)
  delete(idxtmp)

  ; read and average data

  ref5 = infile5->$varnames5(refidx5)$
  ref6 = infile6->$varnames6(refidx6)$

  do i = 0, dimsizes(idx5) - 1
    data = infile5->$varnames5(idx5(i))$
    if (i .eq. 0) then
      mmm5 = data
    else
      mmm5 = mmm5 + data
    end if
  end do

  mmm5 = mmm5 / dimsizes(idx5)
  delete(data)

  do i = 0, dimsizes(idx6) - 1
    data = infile6->$varnames6(idx6(i))$
    if (i .eq. 0) then
      mmm6 = data
    else
      mmm6 = mmm6 + data
    end if
  end do
  mmm6 = mmm6 / dimsizes(idx6)
  delete(data)

  x = infile6->x

  idx60 = ind(x .ge. 60.)
  idx80 = ind(x .ge. 80.)

  print("cumulative frequencies clt >= 60%")
  print("cmip5: " + sum(mmm5(idx60)))
  print("cmip6: " + sum(mmm6(idx60)))
  print("ref: " + sum(ref5(idx60)))
  print("")
  print("cumulative frequencies clt >= 80%")
  print("cmip5: " + sum(mmm5(idx80)))
  print("cmip6: " + sum(mmm6(idx80)))
  print("ref: " + sum(ref5(idx80)))

  system("rm postpdf.nc")
  debugfile = addfile("postpdf.nc", "c")
  debugfile->mmm5 = mmm5
  debugfile->mmm6 = mmm6
  debugfile->ref5 = ref5
  debugfile->ref6 = ref6

end
