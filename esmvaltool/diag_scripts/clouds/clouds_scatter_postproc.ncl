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

  inpath = "/pf/b/b380103/workesm/esmvaltool_output/" \
           + "recipe_clouds_scatter_20211108_095716/work"

  diag = "clouds_scatter_postproc.ncl"

  refname = "ESACCI-CLOUD"
  mmmname = "Multi-model average"

  ; CMIP5
  fname5 = inpath + "/scatter_southernocean_cmip5/scatter_totalcwp/" \
    + "clouds_scatter_clt_totalcwp_so_cmip5.nc"
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
  totalcwp5 = infile5->totalcwp

  ; CMIP6
  fname6 = inpath + "/scatter_southernocean_cmip6/scatter_totalcwp/" \
    + "clouds_scatter_clt_totalcwp_so_cmip6.nc"
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
  totalcwp6 = infile6->totalcwp

  ratio5 = totalcwp5(mmmidx5, :) / totalcwp5(refidx5, :)
  ratio6 = totalcwp6(mmmidx6, :) / totalcwp6(refidx6, :)

  avg5 = avg(ratio5)
  avg6 = avg(ratio6)

  print("cmip5 average ratio = " + avg5)
  print("cmip6 average ratio = " + avg6)

end
