; Read a mask from file and apply it to data
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "./interface_scripts/info_output.ncl"

;************************************************************
undef("applyMask")
procedure applyMask(data,
                    masks[*] : string)
;; Arguments:
;;    data  --  the data to be masked, must have coordinate variables
;;              named "lat" and "lon"
;;    masks --  list of mask file names to apply
;; Description:
;;    Applies the intersection of the masks in the list "masks" to "data".
;;    The requested
;;    masks are expected to be available in the local subfolder 'masks/'.
local intersection_mask, mymask, a, lsdata, coord_index, lat_index, lon_index
begin
    
    funcname = "applyMask"

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering apply_mask.ncl", verbosity, 2)
    ;; Find lat/lon coordinate index in data array
    lat_index = -1
    lon_index = -1
    do coord_index = 0, dimsizes(dimsizes(data)) - 1
        if (data!coord_index .eq. "lat") then
            lat_index = coord_index
        end if
        if (data!coord_index .eq. "lon") then
            lon_index = coord_index
        end if
    end do
    if (lat_index .eq. -1 .or. lon_index .eq. -1) then
        error_msg("f", "apply_mask.ncl", funcname, \
                  "Could not find either lat or lon index")
    end if

    intersection_mask = conform(data, 1.0, -1)
    intersection_mask@_FillValue = 1e+20
    intersection_mask!lat_index = "lat"
    intersection_mask&lat = data&lat
    intersection_mask!lon_index = "lon"
    intersection_mask&lon = data&lon
    ;; Loop and create the intersection of all masks
    do curr_mask_index = 0, dimsizes(masks) - 1
        curr_mask = masks(curr_mask_index)
        a = addfile("masks/" + curr_mask, "r")
        lsdata = a->LSMASK
        mymask  = landsea_mask(lsdata, intersection_mask&lat, intersection_mask&lon)
        intersection_mask = mask(intersection_mask, mymask .eq. 1, True)
        delete(a)
        delete(lsdata)
    end do
    data = mask(data, intersection_mask .eq. 1, True)
    delete(intersection_mask)
    delete(mymask)
    info_output(">>>>>>>> Leaving apply_mask.ncl", verbosity, 2)
end
