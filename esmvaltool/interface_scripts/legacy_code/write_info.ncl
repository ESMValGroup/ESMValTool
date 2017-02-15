load "./interface_scripts/info_output.ncl"

undef("write_info")
function write_info(display [1] : logical)
;;                  return_val [1] : string
;; Arguments:
;;     display -- logical with attributes listing the items to put in the infobox
;; Return value:
;;     infobox -- string to write to plot
;; Description:
;;     Read the attributes in displayinfo and returns them in a formatted string
local infobox, info_id
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering write_info", verbosity, 8)
    infobox = ""
    infobox = infobox + "~F22~DEBUG INFO:~F21~ ~C~"
    display_varatts = getvaratts(display)
    exceptions = (/"description", "description_ycoord"/)
    do info_id = 0, dimsizes(display_varatts) - 1
        info_tag = display_varatts(info_id)
        if (any(info_tag .eq. exceptions)) then
            continue
        end if
        ;; ~C~ = new line
        infobox = infobox + info_tag + ": " + display@$info_tag$ + "~C~"
    end do
    info_output(">>>>>>>> Leaving write_info", verbosity, 8)
    return(infobox)
end
