;;#############################################################################
;; REFORMAT SCRIPT FOR Melpitz OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    Birmili et al., J. Geophys. Res., 106, D23,
;;    doi:10.1029/2000JD000220, 2001.
;;
;; Last access
;;    Size distributions fits as published in the source.
;;
;; Download and processing instructions
;;
;; Caveats
;;
;; Modification history
;;    20140911-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "Melpitz"

    ;; Tier
    TIER = 2

    ;; Output dir (CMOR-ized data)
    OUTDIR = "/export/pa_data02/ESMVal/obs/Tier" + TIER + "/" + OBSNAME + "/"

    ;; Field
    FIELD = "T0M"

    ;; Corresponding standard names
    NAMES = "sizecn"

    ;; Station position
    MINLAT = 51.5333
    MAXLAT = MINLAT
    MINLON = 12.9333
    MAXLON = MINLON

    ;; Period
    PERIOD = "1996/03/26-1997/08/08"

    ;; Cases
    CASES = (/"NH-maritime", "NH-continental", "Midlatitudes", \
              "SH-maritime", "SH-continental"/)
    NCASES = dimsizes(CASES)

    ;; Modes
    MODES = (/"Nucleation", "Aged nucleation", "Aitken", "Accumulation", \
              "Accumulation 2"/)
    NMODES = dimsizes(MODES)

    ;; Data array
    DATA = new((/NCASES, NMODES, 3/), float)

end

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

begin

    verbosity = 1
    diag_script = "reformat_obs_Melpitz.ncl"

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Define data (Table 2 in Birmili et al., 2001)
    ;;             Dg[nm] N[cm-3] sigma
    DATA(0, 0, :) = (/5.9, 630, 1.33/)
    DATA(0, 1, :) = (/15.8, 2010, 1.61/)
    DATA(0, 2, :) = (/58, 2300, 1.73/)
    DATA(0, 3, :) = (/199, 500, 1.52/)
    DATA(0, 4, :) = (/380, 67, 1.41/)
    DATA(1, 0, :) = (/5.9, 630, 1.33/)
    DATA(1, 1, :) = (/15.8, 2010, 1.61/)
    DATA(1, 2, :) = (/76, 3300, 1.71/)
    DATA(1, 3, :) = (/206, 1200, 1.57/)
    DATA(1, 4, :) = (/350, 330, 1.44/)
    DATA(2, 0, :) = (/5.9, 630, 1.33/)
    DATA(2, 1, :) = (/15.8, 2010, 1.61/)
    DATA(2, 2, :) = (/68, 2600, 1.70/)
    DATA(2, 3, :) = (/203, 720, 1.53/)
    DATA(2, 4, :) = (/320, 360, 1.44/)
    DATA(3, 0, :) = (/5.9, 630, 1.33/)
    DATA(3, 1, :) = (/15.8, 2010, 1.61/)
    DATA(3, 2, :) = (/69, 2400, 1.72/)
    DATA(3, 3, :) = (/218, 510, 1.53/)
    DATA(3, 4, :) = (/300, 300, 1.46/)
    DATA(4, 0, :) = (/5.9, 630, 1.33/)
    DATA(4, 1, :) = (/15.8, 2010, 1.61/)
    DATA(4, 2, :) = (/66, 4900, 1.55/)
    DATA(4, 3, :) = (/187, 1200, 1.55/)
    DATA(4, 4, :) = (/280, 200, 1.44/)

    ;; Convert units
    DATA(:, :, 0) = DATA(:, :, 0) * 1.e-9  ; [nm] --> [m]
    DATA(:, :, 1) = DATA(:, :, 1) * 1.e6  ; [cm-3] --> [m-3]

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
    OUTPUT&column = (/"median"/)
    OUTPUT&case@info          = "Different air masses"
    OUTPUT&diam@long_name     = "diameter"
    OUTPUT&diam@units         = "m"
    OUTPUT&diam@standard_name = "particle_dry_diameter"
    OUTPUT = 0.

    ;; Loop over cases
    do cc = 0, NCASES - 1
        ;; Calculate and add up size distributions for each mode
        do mm = 0, NMODES - 1
            OUTPUT(cc, :, 0) = OUTPUT(cc, :, 0) + \
                lognormal_dist(DATA(cc, mm, 1), DATA(cc, mm, 0), \
                               DATA(cc, mm, 2), DIAM)
        end do
    end do

    ;; Replace zeros with missing values
    OUTPUT = where(OUTPUT.eq.0, OUTPUT@_FillValue, OUTPUT)

    ;; Set variable attributes (CMOR)
    tmp = var_attrib(OUTPUT, NAMES)
    delete(OUTPUT)
    OUTPUT = tmp
    delete(tmp)

    ;; Append global attributes
    gAtt = True
    gAtt@comment       = "Size distributions calculated from the given" + \
        " lognormal parameters"
    gAtt@history       = "Created on " + systemfunc("date")
    gAtt@host          = systemfunc("echo $HOST")
    gAtt@user          = systemfunc("echo $USER")
    gAtt@lonmax        = where(MAXLON.lt.0., \
                               MAXLON + 360., \
                               MAXLON)
    gAtt@lonmin        = where(MINLON.lt.0., \
                               MINLON + 360., \
                               MINLON)
    gAtt@latmax        = MAXLAT
    gAtt@latmin        = MINLAT
    gAtt@info          = ""
    gAtt@period        = PERIOD
    gAtt@field         = FIELD
    gAtt@tier          = TIER
    gAtt@source        = "Size distributions fits as published in the source."
    gAtt@reference     = "Birmili et al., J. Geophys. Res., 106, D23, " + \
        "doi:10.1029/2000JD000220, 2001"
    gAtt@title         = OBSNAME + \
        " particle size distributions reformatted for the ESMValTool"
    gAtt@conventions   = "CF/CMOR"

    ;; Outfile
    fout = OUTDIR + "OBS_" + OBSNAME + "_campaign_1_" + \
        FIELD + "_" + NAMES + ".nc"

    ;; Write variable
    write_nc_profile(fout, NAMES, OUTPUT, gAtt)

    delete(OUTPUT)

end
