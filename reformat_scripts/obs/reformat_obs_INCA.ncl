;;#############################################################################
;; REFORMAT SCRIPT FOR INCA OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 3: restricted dataset.
;;
;; Source
;;    Minikin et al., Geophys. Res. Lett., 30, 1503,
;;    doi:10.1029/2002GL016458, 2003.
;;
;; Last access
;;    Original data from A. Lauer.
;;
;; Download and processing instructions
;;    Data has been converted with dos2unix and untabified.
;;
;; Caveats
;;    Altitudes are converted to pressure levels assuming standard pressures.
;;    Resulting profiles are slightly different from those published in the
;;    above paper, possibly due to binning.
;;
;; Modification history
;;    20140414-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "INCA"

    ;; Tier
    TIER = 3

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Field
    FIELD = "T0M"

    ;; Selected variables
    VARIABLES = (/"acc", "N14", "N5"/)

    ;; Corresponding standard names
    NAMES = (/"conccnSTPd120", "conccnSTPd14", "conccnSTPd5"/)

    ;; Campaigns (see INDIR/readme.txt)
    CAMPAIGN = (/"Prestwick", "Punta-Arenas"/)
    CODE     = (/"pw", "pa"/)
    MINLAT   = (/54., -58./)
    MAXLAT   = (/61., -50./)
    MINLON   = (/351., 276./)
    MAXLON   = (/4., 291./)
    PERIOD   = (/"2000/09/29-2000/10/13", "2000/03/23-2000/04/13"/)

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_INCA.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Generate data array
    NLEVS = 14
    NCOLS = 6
    OUTPUT = new((/NLEVS, NCOLS/), float)
    OUTPUT!0 = "plev"
    OUTPUT!1 = "column"
    tmp = stdatmus_z2tdp(ispan(500, 13500, 1000))
    OUTPUT&plev = todouble(100. * tmp(2, :))
    OUTPUT&column = (/"N", "mean", "stddev", "25%", "median", "75%"/)
    OUTPUT&plev@long_name     = "pressure"
    OUTPUT&plev@positive      = "down"
    OUTPUT&plev@axis          = "Z"
    OUTPUT&plev@units         = "Pa"
    OUTPUT&plev@standard_name = "air_pressure"
    delete(tmp)

    ;; Loop over campaigns
    do cpID = 0, dimsizes(CAMPAIGN) - 1

        ;; Loop over variables
        do vID = 0, dimsizes(VARIABLES) - 1

            info_output("Processing " + VARIABLES(vID), verbosity, 1)

            infile = INDIR + "prfl" + "_" + CODE(cpID) + "_" + \
                VARIABLES(vID) + ".txt"

            ;; Read header
            head = readAsciiHead(infile, 1)

            ;; Read data
            data = readAsciiTable(infile, 1, "string", 1)

            ;; Remove ","
            data = str_sub_str(data, ",", "")

            ;; Number of lines
            nlines = dimsizes(data(:, 0))

            ;; Copy metadata
            output = OUTPUT

            ;; Read number of observations
            tmp = tofloat(str_get_field(data(:, 0), 7, " "))
            output(0:nlines - 1, 0) = (/tmp/)
            delete(tmp)

            ;; Read mean
            tmp = tofloat(str_get_field(data(:, 0), 5, " "))
            output(0:nlines - 1, 1) = (/tmp/)
            delete(tmp)

            ;; Read standard deviation
            tmp = tofloat(str_get_field(data(:, 0), 6, " "))
            output(0:nlines - 1, 2) = (/tmp/)
            delete(tmp)

            ;; Read 25 percentile
            tmp = tofloat(str_get_field(data(:, 0), 3, " "))
            output(0:nlines - 1, 3) = (/tmp/)
            delete(tmp)

            ;; Read median
            tmp = tofloat(str_get_field(data(:, 0), 2, " "))
            output(0:nlines - 1, 4) = (/tmp/)
            delete(tmp)

            ;; Read 75 percentile
            tmp = tofloat(str_get_field(data(:, 0), 4, " "))
            output(0:nlines - 1, 5) = (/tmp/)
            delete(tmp)
            delete(data)

            ;; Convert units (except # of measurements)
            output(:, 1::) = \
                output(:, 1::) * 1.e6  ; [1/cm3-STP] --> [1/m3-STP]

            ;; Set variable attributes (CMOR)
            tmp = var_attrib(output, NAMES(vID))
            delete(output)
            output = tmp
            delete(tmp)

            ;; Append global attributes
            gAtt = True
            gAtt@history       = "Created on " + systemfunc("date")
            gAtt@host          = systemfunc("echo $HOST")
            gAtt@user          = systemfunc("echo $USER")
            gAtt@lonmax        = where(MAXLON(cpID).lt.0., \
                                       MAXLON(cpID) + 360., \
                                       MAXLON(cpID))
            gAtt@lonmin        = where(MINLON(cpID).lt.0., \
                                       MINLON(cpID) + 360., \
                                       MINLON(cpID))
            gAtt@latmax        = MAXLAT(cpID)
            gAtt@latmin        = MINLAT(cpID)
            gAtt@info          = CAMPAIGN(cpID)
            gAtt@period        = PERIOD(cpID)
            gAtt@field         = FIELD
            gAtt@tier          = TIER
            gAtt@source        = "Original data from Axel Lauer (DLR)"
            gAtt@reference     = "Minikin et al., Geophys. Res. Lett., " + \
                "30, 1503, doi:10.1029/2002GL016458, 2003"
            gAtt@title         = OBSNAME + \
                " vertical profiles reformatted for the ESMValTool"
            gAtt@conventions   = "CF/CMOR"

            ;; Outfile
            d1 = str_get_field(str_get_field(PERIOD(cpID), 1, "-"), 1, "/") + \
                str_get_field(str_get_field(PERIOD(cpID), 1, "-"), 2, "/")
            d2 = str_get_field(str_get_field(PERIOD(cpID), 2, "-"), 1, "/") + \
                str_get_field(str_get_field(PERIOD(cpID), 2, "-"), 2, "/")
            fout = OUTDIR + "OBS_" + OBSNAME + "_campaign_" + \
                CAMPAIGN(cpID) + "_" + FIELD + "_" + NAMES(vID) + \
                "_" + d1 + "-" + d2 + ".nc"

            ;; Write variable
            write_nc_profile(fout, NAMES(vID), output, gAtt)

            delete(output)

        end do
    end do

    leave_msg(diag_script, "", 1)

end
