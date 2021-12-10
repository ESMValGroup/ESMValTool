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

  inpath = "/pf/b/b380103/workesm/esmvaltool_output/" + \
           "recipe_clouds_taylor2_20211109_085428/work/clouds_tay"
  var = "clivi"

  diag = "clouds_taylor_show.ncl"

  fname = inpath + "/tay" + var + "/clouds_taylor_" + var + ".nc"
  infile = addfile(fname, "r")
  invar = infile->var

  tmp = infile->models
  dims = dimsizes(tmp)
  model = new(dims(0), string)
  do i = 0, dims(0) - 1
    model(i) = tostring(tmp(i, :))
  end do
  delete(tmp)

  do i = 0, dims(0) - 1
    print(model(i) + " ratio=" + invar(i, 0, 0) + " corr=" + invar(i, 1, 0))
  end do

end
