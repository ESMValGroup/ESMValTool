load "./diag_scripts/lib/ncl/set_operators.ncl"

undef("write_infobox")
function write_infobox(info [1]       : logical,
                       validItems [*] : string)
;;                     return_val [1] : string
;; Arguments:
;;     info -- boolean with attributes defined in variable_defs/*-files
;;     validItems -- valid items to put in the infobox
;; Return value:
;;     infobox -- string to write to plot
;; Description:
;;     Read the info-attributes and returns them in a formatted string
;;     Items not in the 'validItems' list will be excluded from the return
;;     string.
local infobox
begin
    infobox = ""
    if (isatt(info, "masks") .and. inlist("masks", validItems)) then
        infobox = "Mask(s) used:  " + info@masks(0)
        do imasks = 1, dimsizes(info@masks) - 1
            infobox = infobox + "~C~                 " + info@masks(imasks)  ; ~C~ = new line
        end do

    end if

    if (isatt(info, "cf_latrange") .and. inlist("cf_latrange", validItems)) then
        if (strlen(infobox).gt.0) then
            infobox = infobox + "~C~"
        end if
        infobox = infobox + "latrange: " + info@cf_latrange(0) + "-" + info@cf_latrange(1)
    end if

    if (isatt(info, "cf_lonrange") .and. inlist("cf_lonrange", validItems)) then
        if (strlen(infobox).gt.0) then
            infobox = infobox + "~C~"
        end if
        infobox = infobox + "lonrange: " + info@cf_lonrange(0) + "-" + info@cf_lonrange(1)
    end if

    if (isatt(info, "scd_latrange") .and. inlist("scd_latrange", validItems)) then
        if (strlen(infobox).gt.0) then
            infobox = infobox + "~C~"
        end if
        infobox = infobox + "latrange: " + info@cf_latrange(0) + "-" + info@cf_latrange(1)
    end if

    if (isatt(info, "scd_lonrange") .and. inlist("scd_lonrange", validItems)) then
        if (strlen(infobox).gt.0) then
            infobox = infobox + "~C~"
        end if
        infobox = infobox + "lonrange: " + info@cf_lonrange(0) + "-" + info@cf_lonrange(1)
    end if

    if (isatt(info, "season") .and. inlist("season", validItems)) then
        if (strlen(infobox).gt.0) then
            infobox = infobox + "~C~"
        end if
        infobox = infobox + "season: " + info@season
    end if

    if (isatt(info, "years") .and. inlist("years", validItems)) then
        if (strlen(infobox).gt.0) then
            infobox = infobox + "~C~"
        end if
        infobox = infobox + "years: " + info@years
    end if
    return(infobox)
end
