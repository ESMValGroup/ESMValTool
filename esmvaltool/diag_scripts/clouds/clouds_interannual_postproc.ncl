; CLOUDS_INTERANNUAL_POSTPROC
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; ############################################################################
; Description
;
; Modification history
;   20211210-lauer_axel: written.
;
; ############################################################################

begin

  inpath = "/pf/b/b380103/workesm/esmvaltool_output/" \
           + "recipe_clouds_interannual_20211105_172010/work"

  diag = "clouds_interannual_postproc.ncl"

  refnames = (/"ESACCI-CLOUD", "CLARA-AVHRR", "CERES-EBAF", \
               "ESACCI-WATERVAPOUR", "MODIS", "PATMOS-x", "ISCCP-FH", \
               "MAC-LWP", "SSMI-MERIS", "CLOUDSAT-L2"/)

  vars = (/"clivi", "clt", "lwcre", "lwp", "prw", "swcre"/)

  do iv = 0, dimsizes(vars) - 1
    var = vars(iv)

    print("*** " + var + " ***")

    ; CMIP5
    fname5 = inpath + "/clouds_intera_cmip5/var" + var + "/" \
        + "clouds_interannual_" + var + ".nc"
    infile5 = addfile(fname5, "r")
    varnames5 = getfilevarnames(infile5)
    dims = dimsizes(varnames5)
    idxtmp_o = ispan(0, dims(0) - 1, 1)
    idxtmp_m = ispan(0, dims(0) - 1, 1)
    do i = 0, dims(0) - 1
      if (isStrSubset(varnames5(i), var + "_var_")) then
        imark = 0
        do ir = 0, dimsizes(refnames) - 1
          if (isStrSubset(varnames5(i), refnames(ir))) then
            imark = 1
          end if
        end do
        if (imark .eq. 0) then
          idxtmp_o(i) = -1
        else
          idxtmp_m(i) = -1
        end if
        if (isStrSubset(varnames5(i), "MultiModelMean")) then
          idxtmp_m(i) = -1
          idxtmp_o(i) = -1
        end if
      else
        idxtmp_o(i) = -1
        idxtmp_m(i) = -1
      end if
    end do
    idx5_o = ind(idxtmp_o .ge. 0)
    idx5_m = ind(idxtmp_m .ge. 0)
    delete(dims)
    delete(idxtmp_o)
    delete(idxtmp_m)

    ; CMIP6
    fname6 = inpath + "/clouds_intera_cmip6/var" + var + "/" \
        + "clouds_interannual_" + var + ".nc"
    infile6 = addfile(fname6, "r")
    varnames6 = getfilevarnames(infile6)
    dims = dimsizes(varnames6)
    idxtmp_o = ispan(0, dims(0) - 1, 1)
    idxtmp_m = ispan(0, dims(0) - 1, 1)
    do i = 0, dims(0) - 1
      if (isStrSubset(varnames6(i), var + "_var_")) then
        imark = 0
        do ir = 0, dimsizes(refnames) - 1
          if (isStrSubset(varnames6(i), refnames(ir))) then
            imark = 1
          end if
        end do
        if (imark .eq. 0) then
          idxtmp_o(i) = -1
        else
          idxtmp_m(i) = -1
        end if
        if (isStrSubset(varnames6(i), "MultiModelMean")) then
          idxtmp_m(i) = -1
          idxtmp_o(i) = -1
        end if
      else
        idxtmp_o(i) = -1
        idxtmp_m(i) = -1
      end if
    end do
    idx6_o = ind(idxtmp_o .ge. 0)
    idx6_m = ind(idxtmp_m .ge. 0)
    delete(dims)
    delete(idxtmp_o)
    delete(idxtmp_m)

    ; read and average data

    do i = 0, dimsizes(idx5_o) - 1
      data = infile5->$varnames5(idx5_o(i))$
      if (i .eq. 0) then
        multidim = dimsizes(data)
        newdims = array_append_record(dimsizes(idx5_o), multidim, 0)
        multiobs5 = new(newdims, float)
      end if
      multiobs5(i, :, :) = data
    end do
    ref5 = dim_avg_n(multiobs5, 0)
    delete(data)
    delete(multiobs5)
    do i = 0, dimsizes(idx5_m) - 1
      data = infile5->$varnames5(idx5_m(i))$
      if (i .eq. 0) then
        multidim = dimsizes(data)
        newdims = array_append_record(dimsizes(idx5_m), multidim, 0)
        multimod5 = new(newdims, float)
      end if
      multimod5(i, :, :) = data
    end do
    mmm5 = dim_avg_n(multimod5, 0)
    delete(data)
    delete(multimod5)

    do i = 0, dimsizes(idx6_o) - 1
      data = infile6->$varnames6(idx6_o(i))$
      if (i .eq. 0) then
        multidim = dimsizes(data)
        newdims = array_append_record(dimsizes(idx6_o), multidim, 0)
        multiobs6 = new(newdims, float)
      end if
      multiobs6(i, :, :) = data
    end do
    ref6 = dim_avg_n(multiobs6, 0)
    delete(data)
    delete(multiobs6)
    do i = 0, dimsizes(idx6_m) - 1
      data = infile6->$varnames6(idx6_m(i))$
      if (i .eq. 0) then
        multidim = dimsizes(data)
        newdims = array_append_record(dimsizes(idx6_m), multidim, 0)
        multimod6 = new(newdims, float)
      end if
      multimod6(i, :, :) = data
    end do
    mmm6 = dim_avg_n(multimod6, 0)
    delete(data)
    delete(multimod6)

    outfile = "postproc_" + var + ".nc"
    system("rm " + outfile)
    debugfile = addfile(outfile, "c")
    debugfile->mmm5 = mmm5
    debugfile->mmm6 = mmm6
    debugfile->ref5 = ref5
    debugfile->ref6 = ref6
    delta5 = mmm5
    delta5 = delta5 - ref5
    delta6 = mmm6
    delta6 = delta6 - ref6
    delta56 = mmm5
    delta56 = delta56 - mmm6
    debugfile->delta5 = delta5
    debugfile->delta6 = delta6
    debugfile->delta56 = delta56

    delete(mmm5)
    delete(mmm6)
    delete(delta5)
    delete(delta6)
    delete(delta56)
    delete(ref5)
    delete(ref6)
    delete(idx5_o)
    delete(idx5_m)
    delete(idx6_o)
    delete(idx6_m)
    delete(varnames5)
    delete(varnames6)

  end do  ; variable loop (iv)
end
