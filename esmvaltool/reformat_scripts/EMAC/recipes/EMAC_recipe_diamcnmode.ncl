;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: diamcnmode
;; RECIPE:   particle median dry diameter in each size mode
;; HISTORY:  20140916-A_righ_ma: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
local modes, nmodes, dims, im, nn, imodes, sigma
begin

    ;; Define all possible modes and corresponding sigmas
    modes = get_aerosol_modes()

    ;; Number and dry diameters
    modes = "dryrad" + modes

    ;; Count number of modes
    nmodes = 0
    do im = 0, dimsizes(modes) - 1
        dd = find_and_read_var(modes(im), date, False)
        if (typeof(dd).ne."string") then
            info_output("Adding mode " + modes(im), verbosity, 2)
            nmodes = nmodes + 1
            dims = dimsizes(dd)
        end if
        delete(dd)
    end do

    xx = new(array_append_record(dims, nmodes, 0), float)
    xx!4 = "mode"
    xx&mode = ispan(1, nmodes, 1)

    ;; Loop over modes
    imodes = 0
    sigma = new(nmodes, float)
    do im = 0, dimsizes(modes) - 1
        dd = find_and_read_var(modes(im), date, False)
        if (typeof(dd).ne."string") then
            dd = 2 * dd  ; radius to diameter
            xx(:, :, :, :, imodes) = dd
            sigma(imodes) = modes@sigma(im)
            imodes = imodes + 1
        end if
        delete(dd)
    end do
    xx&mode@sigma = sigma

    return(xx)

end
