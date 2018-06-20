;;#############################################################################
;;# EXTERNAL SEA ICE FUNCTIONS
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    function sea_ice_area
;;
;;#############################################################################

load "./interface_scripts/interface.ncl"
load "./diag_scripts/shared/statistics.ncl"

;;#############################################################################
undef("sea_ice_area")
function sea_ice_area(sic:numeric,
                      imod[1]:integer,
                      data_pointer[1]:logical,
                      region[1]:string,
                      month[1]:string,
                      operation[1]:string,
                      fill_pole_hole[1]:logical)
;;
;; Arguments:
;;    sic            : numeric array of rank 3 (time, lon, lat) containing
;;                     sea ice concentration (sic)
;;    imod           : index of current model in data_pointer
;;    data_pointer   : see interface_scripts/data_handling.ncl
;;    region         : "Arctic" or "Antarctic"
;;    month          : [1, 12] = select month (e.g., "3" = Mar, "9" = Sep)
;;                     "A" = annual mean
;;                     "all" = no average, take all input
;;    operation      : "area" or "extent"
;;    fill_pole_hole : True = fill northern polar hole with sic=100%
;;
;; Return value:
;;    1-d array (time) containing sea ice area or sea ice extent
;;
;; Description:
;;    Calculates sea ice area or extent from sea ice concentration
;;
;; Caveats
;;
;; Modification history:
;;    20180226-A_senf_da: written.
;;
local funcname, scriptname, verbosity, cnt, modelfile, f, vnames, l_irreg, \
      date, month_ind, sic_month, check_month, smonths, month_string, latstep, \
      max_lat, lat1D, actual_lat, data1D, i1D, i2D, dim_i2D, npts, ndim, \
      ring, area_lat, area_lon, area_tmp, SIA_tmp, SIA_tmp2
begin

    funcname = "sea_ice_area"
    scriptname = "diag_scripts/lib/ncl/latlon.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)

    ;; Check region and units
    if (region.ne."Arctic" .and. region.ne."Antarctic") then
        error_msg("f", scriptname, funcname, "Argument 3 (region) invalid.")
    end if
    if (sic@units.eq."%") then
        sic = sic / 100.  ; Subsequent code expects %/100
        sic@units = "%/100"
    end if

    ;; Check for irregular grid (i.e., if "areacello" exists in file)
    cnt = interface_get_dictkeys("sic", imod) + "_filepath"
    modelfile = data_pointer@$cnt$
    f = addfile(modelfile, "r")
    vnames = getfilevarnames(f)
    if (any(vnames.eq."areacello")) then
        l_irreg = True
    else
        l_irreg = False
    end if

    ;; Calculate time coordinate
    date = cd_calendar(sic&time, 0)  ; Create calendar array
    if (month.eq."A") then
        sic_month = time_operations(sic, \
                                    toint(models@start_year(imod)), \
                                    toint(models@end_year(imod)), \
                                    "average", "yearly", True)
        sic_month!0 = "time"
        month_string = "Annual mean "
    else if (month.eq."all") then
        sic_month = sic
        month_string = ""
    else  ; "month" must be a number between 1 and 12
        check_month = tostring(ispan(1,12,1))
        if (.not.any(check_month.eq.month)) then
            error_msg("f", scriptname, funcname, "Argument 4 (month) invalid.")
        end if
        smonths = (/"January", "February", "March", "April", "May", \
                    "June", "July", "August", "September", "October", \
                    "November", "December"/)
        month_string = smonths(ind(check_month.eq.month)) + " "
        month_ind = ind(date(:,1).eq.toint(month))
        sic_month = sic(month_ind, :, :)
    end if
    end if

    ;; Fill pole hole
    if (fill_pole_hole .and. region.eq."Arctic") then
        if (l_irreg) then
            area_lat = extract_area(imod, data_pointer, "lat", "sic")
        else  ; Make 2D grid compatible with irregular grid handling below
            area_lat = conform_dims(dimsizes(map_area(sic&lat, sic&lon)), \
                                    sic&lat, 0)
        end if
        latstep = 0.1
        max_lat = 90.
        lat1D = ndtooned(area_lat)
        do itime=0, dimsizes(sic&time)-1  ; Hole depends on time
            actual_lat = max_lat
            data1D = ndtooned(sic(itime,:,:))
            if (all(ismissing(data1D))) then  ; Time step missing
                continue
            end if
            do while(True)  ; Will break when hole is filled
                i1D = ind(lat1D.ge.actual_lat .and. \
                lat1D.lt.actual_lat+latstep)
                if (.not. all(ismissing(i1D))) then
                    i2D     = ind_resolve(ind(lat1D.ge.actual_lat),\
                    dimsizes(area_lat))
                    dim_i2D = dimsizes(i2D)
                    npts    = dim_i2D(0)  ; Number of elements
                    ndim    = dim_i2D(1)  ; Rank of "lat" 
                    if (any(ismissing(data1D(i1D)))) then
                        do n=0, npts-1
                            ring = sic(itime,i2D(n,0),i2D(n,1))
                            sic(itime,i2D(n,0),i2D(n,1)) = \
                            where(ismissing(ring), 1., ring)
                            delete(ring)
                        end do
                    else  ; Out of the hole!
                        delete([/i1D, i2D, dim_i2D, npts, ndim/])
                        break
                    end if
                    delete([/i2D, dim_i2D, npts, ndim/])
                end if
                delete([/i1D/])  
                actual_lat = actual_lat - latstep
            end do  ; while
            delete(data1D)
        end do  ; itime
        delete([/area_lat, lat1D/])
    end if

    ;; Calculate sea ice area or extent
    if (operation.eq."extent") then
        sic_month = where(sic_month.ge.0.15, 1., 0.)  ; Make sic binary for SIE
    end if
    if (l_irreg) then  ; Irregular grid
        area_lat = extract_area(imod, data_pointer, "lat", "sic")
        area_tmp = extract_area(imod, data_pointer, "areacello", "sic")
        if (region.eq."Arctic") then
            area_tmp = where(area_lat .gt. 0.0, area_tmp / 1.e12, 0.0)
        else
            area_tmp = where(area_lat .lt. 0.0, area_tmp / 1.e12, 0.0)
        end if
        SIA_tmp  = conform(sic_month, area_tmp, (/1, 2/)) * sic_month
        copy_VarMeta(sic_month, SIA_tmp)
        SIA_tmp2 = dim_sum_Wrap(SIA_tmp)
        SIA      = dim_sum_Wrap(SIA_tmp2)
    else  ; Regular grid
        area_tmp = map_area(sic_month&lat, sic_month&lon)
        area_tmp = area_tmp / 1.e12  ; Million square km
        SIA_tmp  = conform(sic_month, area_tmp, (/1, 2/)) * sic_month
        copy_VarMeta(sic_month, SIA_tmp)
        if (region.eq."Arctic") then
            SIA = area_operations(SIA_tmp, 0, 90, 0, 360, "sum", False)
        else
            SIA = area_operations(SIA_tmp, -90, 0, 0, 360, "sum", False)
        end if
    end if
    SIA@units = "10~S~6~N~ km~S~2~N~"
    if (operation.eq."extent") then
        SIA@long_name = month_string + region + " Sea Ice Extent"
        SIA@standard_name = "sea_ice_extent"
    else if (operation.eq."area") then
        SIA@long_name = month_string + region + " Sea Ice Area"
        SIA@standard_name = "sea_ice_area"
    end if
    end if

    return(SIA)

end
