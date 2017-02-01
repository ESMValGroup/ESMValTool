;;#############################################################################
;; REFORMAT SCRIPT FOR Putaud OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    Putaud et al., EUR 20411 EN, Joint Research Center, JRC23835, 2002.
;;
;; Last access
;;    Size distributions fits as published in the source.
;;
;; Download and processing instructions
;;
;; Caveats
;;
;; Modification history
;;    20140925-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "Putaud"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    ;; (not needed)

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Field
    FIELD = "T0M"

    ;; Corresponding standard names
    NAMES = "sizecn"

    ;; Region
    STATIONS = (/"Jungfraujoch", "Aspvreten", "Harwell", "Hohenpeissenberg", \
                "Melpitz", "Ispra", "London Bloomsbury", "Leipzig", \
                "Milano Bresso", "London Marylebone", "Copenhagen Jagtvej"/)
    STATIONS@category = (/"free tropop.", "natural", "rural", "rural", \
                         "near-city", "near-city", "urban", "urban", "urban", \
                         "kerbside", "kerbside"/)
    STATIONS@latitude = (/47.55, 58.80, 51.57, 47.80, 51.53, 45.82, 51.52, \
                         51.35, 45.53, 51.52, 55.70/)
    STATIONS@longitude = (/7.98, 17.40, -1.32, 11.02, 12.93, 8.63, -0.13, \
                          12.43, 9.20, -0.15, 12.55/)
    NSTATIONS = dimsizes(STATIONS)

    ;; Types
    CASES = (/"morning", "afternoon", "night"/)
    NCASES = dimsizes(CASES)

    ;; Modes
    MODES = (/"Mode 1", "Mode 2", "Mode 3"/)
    NMODES = dimsizes(MODES)

    ;; Seasons
    SEASONS = (/"Summer", "Winter"/)
    NSEASONS = dimsizes(SEASONS)

    ;; Month
    MONTHS = (/(/6, 7, 8/), (/1, 2, 12/)/)

    ;; Data array
    DATA = new((/NSEASONS, NSTATIONS, NCASES, NMODES, 3/), float)

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_Putaud.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Define data (Annex 3 in Putaud et al., 2002)

    ;; Summer
    ;;                N[cm-3] CMD[um] sigma
;;    DATA(0, 0, 0, 0, :) = NA
    DATA(0, 0, 0, 1, :) = (/343, 0.059, 1.76/)
    DATA(0, 0, 0, 2, :) = (/180, 0.149, 1.61/)
;;    DATA(0, 0, 1, 0, :) = NA
    DATA(0, 0, 1, 1, :) = (/358, 0.053, 1.71/)
    DATA(0, 0, 1, 2, :) = (/286, 0.139, 1.54/)
;;    DATA(0, 0, 2, 0, :) = NA
    DATA(0, 0, 2, 1, :) = (/389, 0.060, 1.75/)
    DATA(0, 0, 2, 2, :) = (/178, 0.154, 1.56/)

    DATA(0, 1, 0, 0, :) = (/71, 0.024, 1.53/)
    DATA(0, 1, 0, 1, :) = (/1534, 0.073, 1.63/)
    DATA(0, 1, 0, 2, :) = (/324, 0.232, 1.39/)
    DATA(0, 1, 1, 0, :) = (/185, 0.026, 1.56/)
    DATA(0, 1, 1, 1, :) = (/1364, 0.085, 1.61/)
    DATA(0, 1, 1, 2, :) = (/276, 0.246, 1.38/)
    DATA(0, 1, 2, 0, :) = (/160, 0.040, 1.56/)
    DATA(0, 1, 2, 1, :) = (/1718, 0.082, 1.70/)
    DATA(0, 1, 2, 2, :) = (/224, 0.245, 1.38/)

    DATA(0, 2, 0, 0, :) = (/993, 0.027, 1.54/)
    DATA(0, 2, 0, 1, :) = (/639, 0.059, 1.38/)
    DATA(0, 2, 0, 2, :) = (/1415, 0.104, 1.86/)
    DATA(0, 2, 1, 0, :) = (/2089, 0.028, 1.53/)
    DATA(0, 2, 1, 1, :) = (/553, 0.060, 1.36/)
    DATA(0, 2, 1, 2, :) = (/1459, 0.114, 1.69/)
    DATA(0, 2, 2, 0, :) = (/1790, 0.047, 1.76/)
    DATA(0, 2, 2, 1, :) = (/553, 0.070, 1.42/)
    DATA(0, 2, 2, 2, :) = (/773, 0.153, 1.61/)

    DATA(0, 3, 0, 0, :) = (/367, 0.012, 1.52/)
    DATA(0, 3, 0, 1, :) = (/1192, 0.036, 2.00/)
    DATA(0, 3, 0, 2, :) = (/1185, 0.105, 1.99/)
    DATA(0, 3, 1, 0, :) = (/468, 0.013, 1.61/)
    DATA(0, 3, 1, 1, :) = (/1629, 0.030, 2.00/)
    DATA(0, 3, 1, 2, :) = (/1241, 0.115, 1.90/)
    DATA(0, 3, 2, 0, :) = (/94, 0.024, 1.74/)
    DATA(0, 3, 2, 1, :) = (/1668, 0.051, 1.94/)
    DATA(0, 3, 2, 2, :) = (/557, 0.161, 1.64/)

    DATA(0, 4, 0, 0, :) = (/3467, 0.013, 2.00/)
    DATA(0, 4, 0, 1, :) = (/1142, 0.056, 1.49/)
    DATA(0, 4, 0, 2, :) = (/1460, 0.117, 1.84/)
    DATA(0, 4, 1, 0, :) = (/2938, 0.013, 1.95/)
    DATA(0, 4, 1, 1, :) = (/3989, 0.032, 2.00/)
    DATA(0, 4, 1, 2, :) = (/1356, 0.123, 1.73/)
    DATA(0, 4, 2, 0, :) = (/1011, 0.014, 1.73/)
    DATA(0, 4, 2, 1, :) = (/2228, 0.058, 1.65/)
    DATA(0, 4, 2, 2, :) = (/1023, 0.158, 1.66/)

    DATA(0, 5, 0, 0, :) = (/632, 0.010, 1.93/)
    DATA(0, 5, 0, 1, :) = (/2575, 0.03, 1.87/)
    DATA(0, 5, 0, 2, :) = (/3609, 0.095, 1.89/)
    DATA(0, 5, 1, 0, :) = (/638, 0.010, 1.69/)
    DATA(0, 5, 1, 1, :) = (/2508, 0.027, 1.61/)
    DATA(0, 5, 1, 2, :) = (/4863, 0.093, 1.90/)
    DATA(0, 5, 2, 0, :) = (/256, 0.018, 1.73/)
    DATA(0, 5, 2, 1, :) = (/2084, 0.047, 1.71/)
    DATA(0, 5, 2, 2, :) = (/2859, 0.120, 1.81/)

;;    DATA(0, 6, 0, 0, :) = NA
    DATA(0, 6, 0, 1, :) = (/3461, 0.025, 1.49/)
    DATA(0, 6, 0, 2, :) = (/7412, 0.059, 2.00/)
;;    DATA(0, 6, 1, 0, :) = NA
    DATA(0, 6, 1, 1, :) = (/4933, 0.027, 1.56/)
    DATA(0, 6, 1, 2, :) = (/5572, 0.067, 2.00/)
;;    DATA(0, 6, 2, 0, :) = NA
    DATA(0, 6, 2, 1, :) = (/3047, 0.033, 1.85/)
    DATA(0, 6, 2, 2, :) = (/3373, 0.073, 2.00/)

    DATA(0, 7, 0, 0, :) = (/12921, 0.016, 1.99/)
    DATA(0, 7, 0, 1, :) = (/3360, 0.050, 1.65/)
    DATA(0, 7, 0, 2, :) = (/2431, 0.113, 1.97/)
    DATA(0, 7, 1, 0, :) = (/7751, 0.011, 1.78/)
    DATA(0, 7, 1, 1, :) = (/5422, 0.041, 2.00/)
    DATA(0, 7, 1, 2, :) = (/599, 0.189, 1.55/)
    DATA(0, 7, 2, 0, :) = (/3711, 0.015, 2.00/)
    DATA(0, 7, 2, 1, :) = (/5118, 0.043, 1.96/)
    DATA(0, 7, 2, 2, :) = (/676, 0.185, 1.62/)

    DATA(0, 8, 0, 0, :) = (/3664, 0.014, 1.37/)
    DATA(0, 8, 0, 1, :) = (/12914, 0.025, 1.56/)
    DATA(0, 8, 0, 2, :) = (/13961, 0.059, 2.00/)
    DATA(0, 8, 1, 0, :) = (/3187, 0.014, 1.40/)
    DATA(0, 8, 1, 1, :) = (/10981, 0.027, 1.58/)
    DATA(0, 8, 1, 2, :) = (/11527, 0.061, 2.00/)
    DATA(0, 8, 2, 0, :) = (/2917, 0.014, 1.40/)
    DATA(0, 8, 2, 1, :) = (/5233, 0.024, 1.68/)
    DATA(0, 8, 2, 2, :) = (/6178, 0.068, 2.00/)

    DATA(0, 9, 0, 0, :) = (/24092, 0.021, 1.77/)
    DATA(0, 9, 0, 1, :) = (/26919, 0.037, 2.00/)
    DATA(0, 9, 0, 2, :) = (/6872, 0.111, 1.74/)
    DATA(0, 9, 1, 0, :) = (/23586, 0.020, 2.00/)
    DATA(0, 9, 1, 1, :) = (/13734, 0.030, 1.57/)
    DATA(0, 9, 1, 2, :) = (/18452, 0.074, 1.97/)
    DATA(0, 9, 2, 0, :) = (/16333, 0.017, 1.84/)
    DATA(0, 9, 2, 1, :) = (/8192, 0.029, 1.79/)
    DATA(0, 9, 2, 2, :) = (/9670, 0.072, 2.00/)

;;    DATA(0, 10, 0, 0, :) = NA
;;    DATA(0, 10, 0, 1, :) = NA
;;    DATA(0, 10, 0, 2, :) = NA
;;    DATA(0, 10, 1, 0, :) = NA
;;    DATA(0, 10, 1, 1, :) = NA
;;    DATA(0, 10, 1, 2, :) = NA
;;    DATA(0, 10, 2, 0, :) = NA
;;    DATA(0, 10, 2, 1, :) = NA
;;    DATA(0, 10, 2, 2, :) = NA

    ;; Winter
    ;;                N[cm-3] CMD[um] sigma
    DATA(1, 0, 0, 0, :) = (/123, 0.033, 1.80/)
    DATA(1, 0, 0, 1, :) = (/40, 0.066, 1.96/)
    DATA(1, 0, 0, 2, :) = (/33, 0.115, 1.48/)
    DATA(1, 0, 1, 0, :) = (/184, 0.031, 1.99/)
    DATA(1, 0, 1, 1, :) = (/49, 0.077, 2.00/)
    DATA(1, 0, 1, 2, :) = (/28, 0.116, 1.45/)
    DATA(1, 0, 2, 0, :) = (/102, 0.033, 1.78/)
    DATA(1, 0, 2, 1, :) = (/75, 0.048, 2.00/)
    DATA(1, 0, 2, 2, :) = (/31, 0.121, 1.49/)

    DATA(1, 1, 0, 0, :) = (/287, 0.024, 1.51/)
    DATA(1, 1, 0, 1, :) = (/716, 0.058, 1.55/)
    DATA(1, 1, 0, 2, :) = (/425, 0.198, 1.52/)
    DATA(1, 1, 1, 0, :) = (/374, 0.022, 1.72/)
    DATA(1, 1, 1, 1, :) = (/711, 0.061, 1.66/)
    DATA(1, 1, 1, 2, :) = (/354, 0.209, 1.44/)
    DATA(1, 1, 2, 0, :) = (/20, 0.011, 1.42/)
    DATA(1, 1, 2, 1, :) = (/1170, 0.049, 1.75/)
    DATA(1, 1, 2, 2, :) = (/410, 0.199, 1.49/)

    DATA(1, 2, 0, 0, :) = (/758, 0.033, 1.42/)
    DATA(1, 2, 0, 1, :) = (/338, 0.061, 1.34/)
    DATA(1, 2, 0, 2, :) = (/871, 0.116, 1.75/)
    DATA(1, 2, 1, 0, :) = (/670, 0.032, 1.39/)
    DATA(1, 2, 1, 1, :) = (/426, 0.058, 1.36/)
    DATA(1, 2, 1, 2, :) = (/854, 0.113, 1.71/)
    DATA(1, 2, 2, 0, :) = (/520, 0.036, 1.37/)
    DATA(1, 2, 2, 1, :) = (/400, 0.063, 1.34/)
    DATA(1, 2, 2, 2, :) = (/746, 0.130, 1.63/)

    DATA(1, 3, 0, 0, :) = (/363, 0.014, 1.64/)
    DATA(1, 3, 0, 1, :) = (/1119, 0.048, 2.00/)
    DATA(1, 3, 0, 2, :) = (/247, 0.180, 1.50/)
    DATA(1, 3, 1, 0, :) = (/1173, 0.010, 1.63/)
    DATA(1, 3, 1, 1, :) = (/1839, 0.035, 2.00/)
    DATA(1, 3, 1, 2, :) = (/633, 0.138, 1.79/)
    DATA(1, 3, 2, 0, :) = (/80, 0.013, 1.40/)
    DATA(1, 3, 2, 1, :) = (/1086, 0.040, 2.00/)
    DATA(1, 3, 2, 2, :) = (/432, 0.140, 1.74/)

    DATA(1, 4, 0, 0, :) = (/3433, 0.013, 1.78/)
    DATA(1, 4, 0, 1, :) = (/1925, 0.059, 1.76/)
    DATA(1, 4, 0, 2, :) = (/507, 0.200, 1.65/)
    DATA(1, 4, 1, 0, :) = (/3784, 0.012, 1.81/)
    DATA(1, 4, 1, 1, :) = (/1762, 0.051, 1.70/)
    DATA(1, 4, 1, 2, :) = (/603, 0.177, 1.62/)
    DATA(1, 4, 2, 0, :) = (/1756, 0.021, 2.00/)
    DATA(1, 4, 2, 1, :) = (/995, 0.067, 1.57/)
    DATA(1, 4, 2, 2, :) = (/743, 0.163, 1.87/)

;;    DATA(1, 5, 0, 0, :) = NA
    DATA(1, 5, 0, 1, :) = (/6889, 0.041, 1.84/)
    DATA(1, 5, 0, 2, :) = (/6067, 0.125, 1.62/)
;;    DATA(1, 5, 1, 0, :) = NA
    DATA(1, 5, 1, 1, :) = (/4731, 0.035, 1.91/)
    DATA(1, 5, 1, 2, :) = (/7416, 0.110, 1.71/)
;;    DATA(1, 5, 2, 0, :) = NA
    DATA(1, 5, 2, 1, :) = (/2908, 0.052, 1.74/)
    DATA(1, 5, 2, 2, :) = (/9541, 0.118, 1.62/)

;;    DATA(1, 6, 0, 0, :) = NA
    DATA(1, 6, 0, 1, :) = (/15390, 0.027, 1.67/)
    DATA(1, 6, 0, 2, :) = (/9736, 0.068, 2.00/)
;;    DATA(1, 6, 1, 0, :) = NA
    DATA(1, 6, 1, 1, :) = (/14511, 0.027, 1.69/)
    DATA(1, 6, 1, 2, :) = (/11314, 0.066, 2.00/)
;;    DATA(1, 6, 2, 0, :) = NA
    DATA(1, 6, 2, 1, :) = (/5633, 0.026, 1.69/)
    DATA(1, 6, 2, 2, :) = (/5477, 0.070, 2.00/)

    DATA(1, 7, 0, 0, :) = (/21319, 0.018, 1.98/)
    DATA(1, 7, 0, 1, :) = (/1911, 0.079, 1.51/)
    DATA(1, 7, 0, 2, :) = (/1724, 0.161, 1.80/)
    DATA(1, 7, 1, 0, :) = (/16231, 0.016, 2.00/)
    DATA(1, 7, 1, 1, :) = (/2474, 0.063, 1.65/)
    DATA(1, 7, 1, 2, :) = (/2085, 0.145, 1.90/)
    DATA(1, 7, 2, 0, :) = (/9775, 0.017, 2.00/)
    DATA(1, 7, 2, 1, :) = (/1869, 0.064, 1.63/)
    DATA(1, 7, 2, 2, :) = (/1580, 0.153, 1.90/)

;;    DATA(1, 8, 0, 0, :) = NA
    DATA(1, 8, 0, 1, :) = (/13577, 0.031, 1.66/)
    DATA(1, 8, 0, 2, :) = (/15806, 0.114, 1.83/)
;;    DATA(1, 8, 1, 0, :) = NA
    DATA(1, 8, 1, 1, :) = (/12013, 0.028, 1.71/)
    DATA(1, 8, 1, 2, :) = (/19267, 0.096, 2.00/)
;;    DATA(1, 8, 2, 0, :) = NA
    DATA(1, 8, 2, 1, :) = (/4398, 0.029, 1.71/)
    DATA(1, 8, 2, 2, :) = (/12060, 0.123, 1.74/)

    DATA(1, 9, 0, 0, :) = (/11442, 0.018, 2.00/)
    DATA(1, 9, 0, 1, :) = (/40878, 0.025, 1.60/)
    DATA(1, 9, 0, 2, :) = (/24044, 0.066, 2.00/)
    DATA(1, 9, 1, 0, :) = (/11478, 0.023, 2.00/)
    DATA(1, 9, 1, 1, :) = (/37570, 0.026, 1.61/)
    DATA(1, 9, 1, 2, :) = (/31623, 0.066, 2.00/)
    DATA(1, 9, 2, 0, :) = (/12069, 0.021, 1.56/)
    DATA(1, 9, 2, 1, :) = (/13212, 0.031, 2.00/)
    DATA(1, 9, 2, 2, :) = (/4531, 0.087, 2.00/)

;;    DATA(1, 10, 0, 0, :) = NA
    DATA(1, 10, 0, 1, :) = (/23020, 0.025, 1.71/)
    DATA(1, 10, 0, 2, :) = (/9274, 0.075, 1.95/)
;;    DATA(1, 10, 1, 0, :) = NA
    DATA(1, 10, 1, 1, :) = (/13342, 0.024, 1.74/)
    DATA(1, 10, 1, 2, :) = (/5830, 0.073, 2.00/)
;;    DATA(1, 10, 2, 0, :) = NA
    DATA(1, 10, 2, 1, :) = (/5798, 0.034, 1.97/)
    DATA(1, 10, 2, 2, :) = (/648, 0.165, 1.64/)

    ;; Convert units
    DATA(:, :, :, :, 0) = DATA(:, :, :, :, 0) * 1.e6  ; [cm-3] --> [m-3]
    DATA(:, :, :, :, 1) = DATA(:, :, :, :, 1) * 1.e-6  ; [um] --> [m]

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
    OUTPUT&case@info          = "Different times of the day"
    OUTPUT&diam@long_name     = "diameter"
    OUTPUT&diam@units         = "m"
    OUTPUT&diam@standard_name = "particle_dry_diameter"
    OUTPUT = 0.

    ;; Loop over seasons
    do ss = 0, NSEASONS - 1

        ;; Loop over regions
        do rr = 0, NSTATIONS - 1

            output = OUTPUT

            ;; Loop over cases
            do cc = 0, NCASES - 1

                ;; Calculate and add up size distributions for each mode
                do mm = 0, NMODES - 1
                    size = lognormal_dist(DATA(ss, rr, cc, mm, 0), \
                                          DATA(ss, rr, cc, mm, 1), \
                                          DATA(ss, rr, cc, mm, 2), DIAM)
                    if (.not.all(ismissing(size))) then
                        output(cc, :, 0) = output(cc, :, 0) + size
                    end if
                    delete(size)
                end do
            end do

            ;; Replace zeros with missing values
            output = where(output.eq.0, output@_FillValue, output)

            ;; Set variable attributes (CMOR)
            tmp = var_attrib(output, NAMES)
            delete(output)
            output = tmp
            delete(tmp)

            ;; Append global attributes
            gAtt = True
            gAtt@comment  = "Size distributions calculated from the given" + \
                " lognormal parameters"
            gAtt@history  = "Created on " + systemfunc("date")
            gAtt@host     = systemfunc("echo $HOST")
            gAtt@user     = systemfunc("echo $USER")
            gAtt@lonmax   = where(STATIONS@longitude(rr).lt.0., \
                                  STATIONS@longitude(rr) + 360., \
                                  STATIONS@longitude(rr))
            gAtt@lonmin   = gAtt@lonmax
            gAtt@latmax   = STATIONS@latitude(rr)
            gAtt@latmin   = gAtt@latmax
            gAtt@info     = STATIONS(rr) + " - " + SEASONS(ss)
            gAtt@month    = MONTHS(ss, :)
            gAtt@field    = FIELD
            gAtt@tier     = TIER
            gAtt@source = "Size distributions fits as published in the source."
            gAtt@reference = "Putaud et al., EUR 20411 EN, Joint Research " + \
                "Center, JRC23835, 2002."
            gAtt@title    = OBSNAME + \
                " particle size distributions reformatted for the ESMValTool"
            gAtt@conventions = "CF/CMOR"

            if (.not.all(ismissing(output))) then

                ;; Outfile
                fout = OUTDIR + "OBS_" + OBSNAME + "_campaign_" + \
                    str_sub_str(STATIONS(rr), " ", "") + SEASONS(ss) + "_" + \
                    FIELD + "_" + NAMES + "_200201-200212.nc"

                ;; Write variable
                write_nc_profile(fout, NAMES, output, gAtt)
                delete(output)

            end if
        end do
    end do

    leave_msg(diag_script, "", 1)

end
