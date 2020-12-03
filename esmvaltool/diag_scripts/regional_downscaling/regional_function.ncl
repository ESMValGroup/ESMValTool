undef("select_region1")
function select_region1(region:string)
;;
;; Arguments
;;    region: a string specifying the region to be selected.
;;
;; Return value
;;    An array with the region boundary as (latmin, latmax, lonmin, lonmax)
;;    with the name of the region as a string attribute @name.
;;
;; Description
;;    Translates a region specification into lat/lon boundaries and a region
;;    name as an attribute.
;;
;; Modification history
;;    20141205-A_gott_kl: adjusted names to Righi et al. (2015).
;;    20140410-A_fran_fr: extended to midlat, equatorial and polar regions.
;;    20140129-A_fran_fr: written.
;;
local funcname, scriptname, verbosity, region
begin

    funcname = "select_region1"
    scriptname = "diag_scripts/lib/ncl/latlon.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
   ; enter_msg(scriptname, funcname, 10)

    if (region.eq."Global") then
        region_array = (/-90., 90., 0., 360./)
        region_array@name = "Glob"
        ;;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."Tropics") then
        region_array = (/-20., 20., 0., 360./)
        region_array@name = "Trop"
        ;;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."NH extratropics") then
        region_array = (/20., 90., 0., 360./)
        region_array@name = "NHext"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."SH extratropics") then
        region_array = (/-90., -20., 0., 360./)
        region_array@name = "SHext"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."NH equatorial") then
        region_array = (/0., 30., 0., 360./)
        region_array@name = "NHtrop"
        ;;;;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."SH equatorial") then
        region_array = (/-30., 0., 0., 360./)
        region_array@name = "SHtrop"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."Northern Hemisphere") then
        region_array = (/20., 90., 0., 360./)
        region_array@name = "NH"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."Southern Hemisphere") then
        region_array = (/-90., -20., 0., 360./)
        region_array@name = "SH"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."NH midlatitudes") then
        region_array = (/35., 60., 0., 360./)
        region_array@name = "NHmidlat"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."SH midlatitudes") then
        region_array = (/-60., -35., 0., 360./)
        region_array@name = "SHmidlat"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."Arctic") then
        region_array = (/60., 90., 0., 360./)
        region_array@name = "NHpolar"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."Antarctic") then
        region_array = (/-90., -60., 0., 360./)
        region_array@name = "SHpolar"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

    if (region.eq."Equatorial") then
        region_array = (/-10., 10., 0., 360./)
        region_array@name = "EQ"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

; Region latitudes and longitudes following Figure 9.38 Chapter 9 AR5
; (as Comunication from Chou Sin Chan)
;1- ENA Eastern North America
; Lon = 250 to 350
;Lat = 30 to 67

    if (region.eq."ENA") then
        region_array = (/30., 67., 250., 350./)
        region_array@name = "ENA"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;2- WNA Western North America
;Lon = 190 to 250
;Lat = 30 to 67

    if (region.eq."WNA") then
        region_array = (/30., 67., 190., 250./)
        region_array@name = "WNA"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

;3- CAM Central America
;Lon = 240 to 300
;Lat = 10 to 30
    if (region.eq."CAM") then
        region_array = (/10., 30., 240., 300./)
        region_array@name = "CAM"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;4- AMZ Amazon
;Lon = 277 to 329
;Lat = -20 to 10
    if (region.eq."TSA") then
        region_array = (/-20., 10., 277., 329./)
        region_array@name = "TSA"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if

;5- SSA Southern South America
;Lon = 277 to 329
;Lat = -56 to -20
    if (region.eq."SSA") then
        region_array = (/-56., -20., 277., 329./)
        region_array@name = "SSA"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;6- EUM Europe and Mediterranean
;Lon = 350 to 40
;Lat = 30 to 67
    if (region.eq."EUM") then
        region_array = (/30., 67., 350., 40./)
        region_array@name = "ENA"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;7- NAF Northern Africa
;Lon = 342 to 65
;Lat = 16 to 30
    if (region.eq."NAF") then
        region_array = (/16., 30., 342., 65./)
        region_array@name = "NAF"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;8- SAF South Africa
;Lon = 10 to 51
;Lat = -35 to -10
    if (region.eq."SAF") then
        region_array = (/-35., -10., 10., 51./)
        region_array@name = "SAF"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;9- CAF central Africa
;Lon = 342 to 60
;Lat = -10 to 16
    if (region.eq."CAF") then
        region_array = (/-10., 16., 342., 60./)
        region_array@name = "CAF"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;10- NAS Nothern Asia
;Lon = 40 to 167
;Lat = 48 to 67
    if (region.eq."NAS") then
        region_array = (/48., 67., 40., 167./)
        region_array@name = "NAS"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;11- EAS Eastern Asia
;Lon = 100 to 167
;Lat = 20 to 48
    if (region.eq."EAS") then
        region_array = (/20., 48., 100., 167./)
        region_array@name = "EAS"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;12- CAS Central Asia
;Lon = 40 to 100
;Lat = 30 to 48
    if (region.eq."CAS") then
        region_array = (/30., 48., 40., 100./)
        region_array@name = "CAS"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;13- SEA Southeast Asia
;Lon = 95 to 151
;Lat = -11 to 20
    if (region.eq."SEA") then
        region_array = (/-11.,20.,95., 151./)
        region_array@name = "SEA"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;14- SAS South Asia
;Lon = 65 to 95
;Lat = 5 to 30
    if (region.eq."SAS") then
        region_array = (/5., 30., 65., 95./)
        region_array@name = "SAS"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
;15- AUS Oceania
;Lon = 112 to 180
;Lat = -48 to -11
     if (region.eq."AUS") then
        region_array = (/-48., -11., 112., 180./)
        region_array@name = "AUS"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
    end if
     if (region.eq."Arctic") then
        region_array = (/67.5,90.,0., 360./)
        region_array@name = "Arctic"
        ;leave _msg(scriptname, funcname, 10)
        return(region_array)
     end if


end
;***************************
