;;#############################################################################
;; REFORMAT SCRIPT FOR LACE OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    Petzold et al., J. Geophys. Res., 107, D21, 8128,
;;    doi:10.1029/2001JD001043, 2002.
;;
;; Last access
;;    Size distributions fits as published in the source.
;;
;; Download and processing instructions
;;
;; Caveats
;;
;; Modification history
;;    20140918-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "LACE"

    ;; Tier
    TIER = 2

    ;; Input dir
    ;; (not needed)

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Field
    FIELD = "T0M"

    ;; Corresponding standard names
    NAMES = "sizecn"

    ;; Station position
    MINLAT = 51.5
    MAXLAT = 52.7
    MINLON = 13.5
    MAXLON = 14.5

    ;; Period
    PERIOD = "1998/07/31-1998/08/10"

    ;; Region
    REGIONS = (/"Tropopause", "Upper Free Troposphere", \
               "Lower Free Troposphere", "Upper Boundary Layer", \
               "Lower Boundary Layer"/)
    NREGIONS = dimsizes(REGIONS)

    ;; Altitude range [km]
    ALTMIN = (/11.3, 6.1, 4.0, 1.1, 0.3/)
    ALTMAX = (/11.3, 7.3, 4.0, 1.1, 0.3/)

    ;; Types
    CASES = (/"M1", "M2", "M5", "M6", "M6-polluted", "M7", "M7-polluted"/)
    NCASES = dimsizes(CASES)

    ;; Modes
    MODES = (/"Aitken", "Accumulation", "Coarse"/)
    NMODES = dimsizes(MODES)

    ;; Data array
    DATA = new((/NREGIONS, NCASES, NMODES, 3/), float)

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_LACE.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Define data (Table 6 in Petzold et al., 2002)
    ;;                N[cm-3] CMD[um] sigma
    DATA(0, 0, 0, :) = (/105, 0.070, 2.10/)
    DATA(0, 0, 1, :) = (/6, 0.220, 1.60/)
    DATA(0, 0, 2, :) = (/0.02, 1.000, 1.70/)
    DATA(0, 1, 0, :) = (/160, 0.080, 1.70/)
    DATA(0, 1, 1, :) = (/10, 0.200, 1.60/)
    DATA(0, 1, 2, :) = (/0.02, 1.000, 1.70/)
    DATA(0, 2, 0, :) = (/545, 0.032, 2.00/)
    DATA(0, 2, 1, :) = (/1.7, 0.280, 1.35/)
    DATA(0, 2, 2, :) = (/0.08, 0.500, 1.70/)
;;  DATA(0, 3, 0, :) = NA (measurements perturbed by cirrus event)
;;  DATA(0, 3, 1, :) = NA (measurements perturbed by cirrus event)
;;  DATA(0, 3, 2, :) = NA (measurements perturbed by cirrus event)
;;  DATA(0, 4, 0, :) = NA (measurements perturbed by cirrus event)
;;  DATA(0, 4, 1, :) = NA (measurements perturbed by cirrus event)
;;  DATA(0, 4, 2, :) = NA (measurements perturbed by cirrus event)
    DATA(0, 5, 0, :) = (/503, 0.040, 2.00/)
    DATA(0, 5, 1, :) = (/0.8, 0.300, 1.40/)
    DATA(0, 5, 2, :) = (/0.06, 1.000, 1.60/)
;;  DATA(0, 6, 0, :) = NA (measurements perturbed by cirrus event)
;;  DATA(0, 6, 1, :) = NA (measurements perturbed by cirrus event)
;;  DATA(0, 6, 2, :) = NA (measurements perturbed by cirrus event)

    ;;                N[cm-3] CMD[um] sigma
    DATA(1, 0, 0, :) = (/500, 0.047, 2.00/)
    DATA(1, 0, 1, :) = (/20, 0.230, 1.40/)
    DATA(1, 0, 2, :) = (/0.1, 0.900, 1.70/)
    DATA(1, 1, 0, :) = (/480, 0.070, 1.60/)
    DATA(1, 1, 1, :) = (/10, 0.250, 1.40/)
    DATA(1, 1, 2, :) = (/0.3, 0.900, 1.90/)
    DATA(1, 2, 0, :) = (/250, 0.044, 2.00/)
    DATA(1, 2, 1, :) = (/3, 0.280, 1.35/)
    DATA(1, 2, 2, :) = (/2, 0.360, 1.70/)
    DATA(1, 3, 0, :) = (/350, 0.050, 1.80/)
    DATA(1, 3, 1, :) = (/4, 0.250, 1.40/)
    DATA(1, 3, 2, :) = (/0.2, 1.000, 1.60/)
;;  DATA(1, 4, 0, :) = NA (measurements perturbed by cirrus event)
;;  DATA(1, 4, 1, :) = NA (measurements perturbed by cirrus event)
;;  DATA(1, 4, 2, :) = NA (measurements perturbed by cirrus event)
    DATA(1, 5, 0, :) = (/310, 0.034, 2.00/)
    DATA(1, 5, 1, :) = (/0.8, 0.280, 1.40/)
    DATA(1, 5, 2, :) = (/0.16, 0.800, 1.90/)
;;  DATA(1, 6, 0, :) = NA (measurements perturbed by cirrus event)
;;  DATA(1, 6, 1, :) = NA (measurements perturbed by cirrus event)
;;  DATA(1, 6, 2, :) = NA (measurements perturbed by cirrus event)

    ;;                N[cm-3] CMD[um] sigma
    DATA(2, 0, 0, :) = (/330, 0.045, 2.00/)
    DATA(2, 0, 1, :) = (/40, 0.210, 1.40/)
    DATA(2, 0, 2, :) = (/0.2, 0.700, 1.70/)
    DATA(2, 1, 0, :) = (/485, 0.030, 2.00/)
    DATA(2, 1, 1, :) = (/10, 0.210, 1.40/)
    DATA(2, 1, 2, :) = (/0.08, 0.700, 1.80/)
    DATA(2, 2, 0, :) = (/410, 0.057, 2.00/)
    DATA(2, 2, 1, :) = (/250, 0.340, 1.35/)
    DATA(2, 2, 2, :) = (/0.7, 0.900, 1.90/)
    DATA(2, 3, 0, :) = (/270, 0.050, 2.00/)
    DATA(2, 3, 1, :) = (/5, 0.330, 1.35/)
    DATA(2, 3, 2, :) = (/0.1, 1.400, 2.00/)
    DATA(2, 4, 0, :) = (/370, 0.070, 2.00/)
    DATA(2, 4, 1, :) = (/270, 0.300, 1.35/)
    DATA(2, 4, 2, :) = (/1.4, 0.700, 2.00/)
    DATA(2, 5, 0, :) = (/240, 0.042, 2.00/)
    DATA(2, 5, 1, :) = (/6.0, 0.270, 1.35/)
    DATA(2, 5, 2, :) = (/0.15, 1.170, 1.80/)
    DATA(2, 6, 0, :) = (/370, 0.070, 2.00/)
    DATA(2, 6, 1, :) = (/270, 0.320, 1.40/)
    DATA(2, 6, 2, :) = (/0.08, 2.400, 2.00/)

    ;;                N[cm-3] CMD[um] sigma
    DATA(3, 0, 0, :) = (/4400, 0.032, 1.90/)
    DATA(3, 0, 1, :) = (/220, 0.270, 1.33/)
    DATA(3, 0, 2, :) = (/1.5, 2.300, 1.80/)
    DATA(3, 1, 0, :) = (/2000, 0.055, 2.00/)
    DATA(3, 1, 1, :) = (/300, 0.260, 1.26/)
    DATA(3, 1, 2, :) = (/1.0, 1.800, 1.60/)
    DATA(3, 2, 0, :) = (/7000, 0.023, 1.90/)
    DATA(3, 2, 1, :) = (/80, 0.190, 1.35/)
    DATA(3, 2, 2, :) = (/0.5, 2.760, 1.90/)
    DATA(3, 3, 0, :) = (/5000, 0.024, 2.00/)
    DATA(3, 3, 1, :) = (/50, 0.180, 1.35/)
    DATA(3, 3, 2, :) = (/0.8, 1.220, 2.00/)
;;  DATA(3, 4, 0, :) = NA (measurements perturbed by cirrus event)
;;  DATA(3, 4, 1, :) = NA (measurements perturbed by cirrus event)
;;  DATA(3, 4, 2, :) = NA (measurements perturbed by cirrus event)
    DATA(3, 5, 0, :) = (/15800, 0.028, 2.00/)
    DATA(3, 5, 1, :) = (/520, 0.170, 1.35/)
    DATA(3, 5, 2, :) = (/0.8, 2.300, 1.80/)
;;  DATA(3, 6, 0, :) = NA (measurements perturbed by cirrus event)
;;  DATA(3, 6, 1, :) = NA (measurements perturbed by cirrus event)
;;  DATA(3, 6, 2, :) = NA (measurements perturbed by cirrus event)

    ;;                N[cm-3] CMD[um] sigma
    DATA(4, 0, 0, :) = (/5000, 0.030, 2.00/)
    DATA(4, 0, 1, :) = (/250, 0.190, 1.33/)
    DATA(4, 0, 2, :) = (/2.0, 2.300, 1.80/)
    DATA(4, 1, 0, :) = (/4900, 0.035, 2.00/)
    DATA(4, 1, 1, :) = (/450, 0.230, 1.33/)
    DATA(4, 1, 2, :) = (/1.0, 2.200, 1.70/)
    DATA(4, 2, 0, :) = (/7000, 0.023, 1.90/)
    DATA(4, 2, 1, :) = (/80, 0.190, 1.35/)
    DATA(4, 2, 2, :) = (/0.5, 2.760, 1.90/)
;;  DATA(4, 3, 0, :) = NA (measurements perturbed by cirrus event)
;;  DATA(4, 3, 1, :) = NA (measurements perturbed by cirrus event)
;;  DATA(4, 3, 2, :) = NA (measurements perturbed by cirrus event)
;;  DATA(4, 4, 0, :) = NA (measurements perturbed by cirrus event)
;;  DATA(4, 4, 1, :) = NA (measurements perturbed by cirrus event)
;;  DATA(4, 4, 2, :) = NA (measurements perturbed by cirrus event)
    DATA(4, 5, 0, :) = (/18800, 0.027, 2.00/)
    DATA(4, 5, 1, :) = (/550, 0.160, 1.35/)
    DATA(4, 5, 2, :) = (/0.8, 3.200, 1.80/)
;;  DATA(4, 6, 0, :) = NA (measurements perturbed by cirrus event)
;;  DATA(4, 6, 1, :) = NA (measurements perturbed by cirrus event)
;;  DATA(4, 6, 2, :) = NA (measurements perturbed by cirrus event)

    ;; Convert units
    DATA(:, :, :, 0) = DATA(:, :, :, 0) * 1.e6  ; [cm-3] --> [m-3]
    DATA(:, :, :, 1) = DATA(:, :, :, 1) * 1.e-6  ; [um] --> [m]

    ;; Define diameters
    DIAM = set_size_array()
    NDIAM = dimsizes(DIAM)
    delete_VarAtts(DIAM, -1)

    ;; Define output array
    NCOLS = 1
    OUTPUT = new((/NCASES, NDIAM, NCOLS/), float)
    OUTPUT!0 = "case"
    OUTPUT!1 = "diam"
    OUTPUT!2 = "column"
    OUTPUT&case = ispan(1, NCASES, 1)
    OUTPUT&diam = DIAM
    OUTPUT&column = (/"mean"/)
    OUTPUT&case@info          = "Different flights"
    OUTPUT&diam@long_name     = "diameter"
    OUTPUT&diam@units         = "m"
    OUTPUT&diam@standard_name = "particle_dry_diameter"
    OUTPUT = 0.

    ;; Loop over regions
    do rr = 0, NREGIONS - 1

        info_output("Processing region " + REGIONS(rr), verbosity, 1)

        output = OUTPUT

        ;; Loop over cases
        do cc = 0, NCASES - 1

            ;; Calculate and add up size distributions for each mode
            do mm = 0, NMODES - 1
                output(cc, :, 0) = output(cc, :, 0) + \
                    lognormal_dist(DATA(rr, cc, mm, 0), DATA(rr, cc, mm, 1), \
                                   DATA(rr, cc, mm, 2), DIAM)
            end do
        end do

        ;; Replace zeros with missing values
        output = where(output.eq.0, output@_FillValue, output)

        ;; Set variable attributes (CMOR)
        tmp = var_attrib(output, NAMES)
        delete(output)
        output = tmp
        delete(tmp)

        ;; Convert altitude range to level
        tmp = stdatmus_z2tdp(1000 * ALTMIN(rr))
        LEVMAX = tmp(2)
        LEVMAX = 100. * LEVMAX  ; [hPa] --> [Pa]
        tmp = stdatmus_z2tdp(1000 * ALTMAX(rr))
        LEVMIN = tmp(2)
        LEVMIN = 100. * LEVMIN  ; [hPa] --> [Pa]
        delete(tmp)

        ;; Append global attributes
        gAtt = True
        gAtt@comment  = "Size distributions calculated from the given" + \
            " lognormal parameters"
        gAtt@history  = "Created on " + systemfunc("date")
        gAtt@host     = systemfunc("echo $HOST")
        gAtt@user     = systemfunc("echo $USER")
        gAtt@lonmax   = where(MAXLON.lt.0., \
                              MAXLON + 360., \
                              MAXLON)
        gAtt@lonmin   = where(MINLON.lt.0., \
                              MINLON + 360., \
                              MINLON)
        gAtt@latmax   = MAXLAT
        gAtt@latmin   = MINLAT
        gAtt@plevmin  = LEVMIN
        gAtt@plevmax  = LEVMAX
        if (ALTMIN(rr).eq.ALTMAX(rr)) then
            gAtt@info = REGIONS(rr) + " (" + ALTMIN(rr) + " km)"
        else
            gAtt@info = REGIONS(rr) + " (" + ALTMIN(rr) + \
                "-" + ALTMAX(rr) + " km)"
        end if
        gAtt@period   = PERIOD
        gAtt@field    = FIELD
        gAtt@tier     = TIER
        gAtt@source = "Size distributions fits as published in the source."
        gAtt@reference = "Petzold et al., J. Geophys. Res., 107, D21, " + \
            "8128, doi:10.1029/2001JD001043, 2002."
        gAtt@title    = OBSNAME + \
            " particle size distributions reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"

        ;; Outfile
        d1 = str_get_field(str_get_field(PERIOD, 1, "-"), 1, "/") + \
            str_get_field(str_get_field(PERIOD, 1, "-"), 2, "/")
        d2 = str_get_field(str_get_field(PERIOD, 2, "-"), 1, "/") + \
            str_get_field(str_get_field(PERIOD, 2, "-"), 2, "/")
        fout = OUTDIR + "OBS_" + OBSNAME + "_campaign_" + \
            str_sub_str(REGIONS(rr), " ", "") + "_" + FIELD + "_" + \
            NAMES + "_" + d1 + "-" + d2 + ".nc"

        ;; Write variable
        write_nc_profile(fout, NAMES, output, gAtt)

        delete(output)
    end do

    leave_msg(diag_script, "", 1)

end
