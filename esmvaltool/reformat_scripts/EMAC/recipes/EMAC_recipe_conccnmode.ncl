;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: conccnmode
;; RECIPE:   number concentration of particles in each size mode
;; HISTORY:  20140916-A_righ_ma: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
local modes, nmodes, dims, im, nn, mode, imodes, sigma, pp, tt
begin

    ;; Define all possible modes and corresponding sigmas
    modes = get_aerosol_modes()

    ;; Number and dry diameters
    modes = "N" + modes

    ;; Count number of modes
    nmodes = 0
    do im = 0, dimsizes(modes) - 1
        nn = find_and_read_var(modes(im), date, False)
        if (typeof(nn).ne."string") then
            info_output("Adding mode " + modes(im) + " (" + date + ")", \
                        verbosity, 2)
            nmodes = nmodes + 1
            dims = dimsizes(nn)
        end if
        delete(nn)
    end do

    xx = new(array_append_record(dims, nmodes, 0), float)
    xx!4 = "mode"
    xx&mode = ispan(1, nmodes, 1)

    ;; Loop over modes
    imodes = 0
    sigma = new(nmodes, float)
    do im = 0, dimsizes(modes) - 1
        nn = find_and_read_var(modes(im), date, False)
        if (typeof(nn).ne."string") then
            xx(:, :, :, :, imodes) = nn
            sigma(imodes) = modes@sigma(im)
            imodes = imodes + 1
        end if
        delete(nn)
    end do
    xx&mode@sigma = sigma

    ;; Convert units [1/mol] --> [m-3]
    if (xx@units.eq."1/mol" .or. xx@units.eq."mol-1") then
        pp = find_and_read_var("press", date, True)
        tt = find_and_read_var("tm1", date, True)
        xx = xx * conform(xx, pp, ispan(0, 3, 1))  / Rgas / \
            conform(xx, tt, ispan(0, 3, 1))
        xx@units = "m-3"
    end if

    return(xx)

end
