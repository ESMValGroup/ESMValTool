load "./interface_data/ncl.interface"

undef("output_area")
procedure output_area(diri_out: string, imod: integer, area_lon, area_lat, area_temp)
begin
    filn = diri_out + "/" + "ext_area" + imod + ".nc"
    if (isfilepresent_esmval(filn)) then
        system("/bin/rm " + filn)
    end if
    fout_area = addfile(filn, "c")
    fout_area->lons = area_lon
    fout_area->lats = area_lat
    fout_area->area = area_temp
    delete(fout_area)
end
