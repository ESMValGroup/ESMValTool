load "./interface_data/ncl.interface"

undef("remove_anom")
procedure remove_anom(imod, aYdata, ns, sYdata)
begin
    ;; annual
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering remove_anom", verbosity, 4)
    tmp = avg(aYdata(imod, :))
    aYdata(imod, :) = aYdata(imod, :) - tmp
    delete(tmp)
    ;; seasonal
    do s = 0, ns - 1
        tmp = avg(sYdata(imod, :, s))
        sYdata(imod, :, s) = sYdata(imod, :, s) - tmp
        delete(tmp)
    end do
    info_output(">>>>>>>> Leaving remove_anom", verbosity, 4)
end
