load "./interface_data/ncl.interface"

undef("set_ref_mask")
function set_ref_mask(plot_typ_info, imod, models, refModel,  A0)
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering set_ref_mask", verbosity, 4)
    if (isatt(diag_script_info, "refModel")) then
        info_output(" imod: " + imod\
                              + " model " + models@name(imod)\
                              + " refModel " + refModel, verbosity, 4)
        if (imod .eq. 0  .and. models@name(imod) .eq. refModel) then
            ref_data = conform(A0, 1.0, -1)
            ref_mask = mask(ref_data, A0 .lt. abs(A0@_FillValue), True)
        else if (imod .eq. 0) then
            print("fatal: any reference model must be input first " + refModel + models@name(imod))
            status_exit(1)
        end if
        end if
    else
        ref_mask = conform(A0, 1.0, -1)
    end if
    if(.not. all(iscoord(A0, (/"lat", "lon"/)))) then
        print("fatal: can't use masking by reference model for irregular grids")
        status_exit(1)
    end if
    ref_mask@_FillValue = A0@_FillValue
    ref_mask!1 = "lat"
    ref_mask&lat = A0&lat
    ref_mask!2 = "lon"
    ref_mask&lon = A0&lon
    ref_mask!0 = "time"
    ref_mask&time = A0&time

    ;; output mask for debugging purposes only
    filn =  "../../plots/GO_tsline" + "mask" + imod +".nc"
    if (isfilepresent_esmval(filn)) then
        system("/bin/rm " + filn)
    end if
    fout_mask = addfile(filn, "c")
    fout_mask->lon = A0&lon
    fout_mask->lat = A0&lat
    fout_mask->time = A0&time
    fout_mask->mask = ref_mask

    info_output(">>>>>>>> Leaving set_ref_mask", verbosity, 4)
    return(ref_mask)
end
