; CLOUDS_SCATTER_POSTPROC
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; ############################################################################
; Description
;
; Modification history
;   20211129-lauer_axel: written.
;
; ############################################################################

begin

  inpath = "/home/b/b380103/workesm/esmvaltool_output/" \
           + "recipe_lauer22jclim_fig9-11ab_scatter_20220708_090817/work"

  diag = "clouds_scatter_postproc.ncl"

  refname = "ESACCI-CLOUD"
  mmmname = "Multi-model average"
;  reg_short = "so"
  reg_short = "sep"
;  reg_short = "itcz"
;  reg_long = "southernocean"
  reg_long = "sepacific"
;  reg_long = "itcz"
  var = "totalcwp"
;  var = "lwp"
;  var = "swcre"

  ; CMIP5
  fname5 = inpath + "/scatter_" + reg_long + "_cmip5/scatter_" + var + "/" \
    + "clouds_scatter_clt_" + var + "_" + reg_short + "_cmip5.nc"
  infile5 = addfile(fname5, "r")

  tmp = infile5->model
  dims = dimsizes(tmp)
  models5 = new(dims(0), string)
  do i = 0, dims(0) - 1
    models5(i) = tostring(tmp(i, :))
    if (isStrSubset(models5(i), refname)) then
      refidx5 = i
    end if
    if (isStrSubset(models5(i), mmmname)) then
      mmmidx5 = i
    end if
  end do
  delete(tmp)
  delete(dims)
  data5 = infile5->$var$

  ; CMIP6
  fname6 = inpath + "/scatter_" + reg_long + "_cmip6/scatter_" + var + "/" \
    + "clouds_scatter_clt_" + var + "_" + reg_short + "_cmip6.nc"
  infile6 = addfile(fname6, "r")

  tmp = infile6->model
  dims = dimsizes(tmp)
  models6 = new(dims(0), string)
  do i = 0, dims(0) - 1
    models6(i) = tostring(tmp(i, :))
    if (isStrSubset(models6(i), refname)) then
      refidx6 = i
    end if
    if (isStrSubset(models6(i), mmmname)) then
      mmmidx6 = i
    end if
  end do
  delete(tmp)
  delete(dims)
  data6 = infile6->$var$

  ratio5 = data5(mmmidx5, :) / data5(refidx5, :)
  ratio6 = data6(mmmidx6, :) / data6(refidx6, :)

  avg5 = avg(ratio5)
  avg6 = avg(ratio6)

  print("cmip5 average ratio = " + avg5)
  print("cmip6 average ratio = " + avg6)

end
