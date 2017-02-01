;;#############################################################################
;; REFORMAT SCRIPT FOR Texas OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 3: restricted dataset.
;;
;; Source
;;    Schwarz et al., J. Geophys. Res., 111, D16207,
;;    doi:10.1029/2006JD007076, 2006.
;;
;; Last access
;;    Data provided by J. Schwarz.
;;
;; Download and processing instructions
;;    Data has been untabified.
;;
;; Caveats
;;
;; Modification history
;;    20140909-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "Texas"

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

    ;; Standard names
    NAMES = (/"mmrbc", "mmraer"/)

    ;; Corresponding columns
    COLUMNS = (/4, 2/)

    ;; Campaigns
    CAMPAIGN = (/"10", "12"/)
    PERIOD   = (/"2004/11/10-2004/11/10", "2004/11/12-2004/11/12"/)

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_Texas.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Loop over campaigns
    do cpID = 0, dimsizes(CAMPAIGN) - 1

        info_output("Campaign: " + CAMPAIGN(cpID), verbosity, 1)

        ;; Find lat and lon ranges
        infile = INDIR + "200411" + CAMPAIGN(cpID) + "_AVE_Flight_Track.txt"
        data = readAsciiTable(infile, 1, "string", 1)
        lat = tofloat(str_get_field(data(:, 0), 3, " "))
        lon = tofloat(str_get_field(data(:, 0), 4, " "))
        MINLAT = min(lat)
        MINLON = min(lon)
        MAXLAT = max(lat)
        MAXLON = max(lon)
        delete(lat)
        delete(lon)
        delete(data)

        ;; Read data
        infile = INDIR + "200411" + CAMPAIGN(cpID) + "_BC_Data.txt"
        data = readAsciiTable(infile, 1, "string", 1)
        NLEVS = dimsizes(data(:, 0))
        NCOLS = 2
        OUTPUT = new((/NLEVS, NCOLS/), float)
        OUTPUT!0 = "plev"
        OUTPUT!1 = "column"
        OUTPUT&plev = todouble(ispan(1, NLEVS, 1))  ; dummy
        OUTPUT&column = (/"mean", "stddev"/)
        OUTPUT&plev@long_name     = "pressure"
        OUTPUT&plev@positive      = "down"
        OUTPUT&plev@axis          = "Z"
        OUTPUT&plev@units         = "Pa"
        OUTPUT&plev@standard_name = "air_pressure"

        ;; Loop over variables
        do vID = 0, dimsizes(NAMES) - 1

            ;; Copy metadata
            output = OUTPUT

            ;; Read pressure coordinate and convert units [hPa] --> [Pa]
            output&plev = todouble(str_get_field(data(:, 0), 1, " "))
            output&plev = output&plev * 100.

            ;; Read mean and convert units [ng/kg] --> [kg/kg]
            output(:, 0) = \
                tofloat(str_get_field(data(:, 0), COLUMNS(vID), " "))
            output(:, 0) = output(:, 0) * 1.e-12

            ;; Read standard deviation and convert units [ng/kg] --> [kg/kg]
            output(:, 1) = \
                tofloat(str_get_field(data(:, 0), COLUMNS(vID) + 1, " "))
            output(:, 1) = output(:, 1) * 1.e-12

            ;; Set variable attributes (CMOR)
            tmp = var_attrib(output, NAMES(vID))
            delete(output)
            output = tmp
            delete(tmp)

            ;; Append global attributes
            gAtt = True
            gAtt@comment       = "The data includes a scaling factor of " + \
                "1.7 over the entire profile to account for mass at sizes " + \
                "outside the SP2 detection range"
            gAtt@history       = "Created on " + systemfunc("date")
            gAtt@host          = systemfunc("echo $HOST")
            gAtt@user          = systemfunc("echo $USER")
            gAtt@lonmax        = where(MAXLON.lt.0., MAXLON + 360., MAXLON)
            gAtt@lonmin        = where(MINLON.lt.0., MINLON + 360., MINLON)
            gAtt@latmax        = MAXLAT
            gAtt@latmin        = MINLAT
            gAtt@info          = "Nov. " + CAMPAIGN(cpID)
            gAtt@period        = PERIOD(cpID)
            gAtt@field         = FIELD
            gAtt@tier          = TIER
            gAtt@source        = "Data provided by Joshua Schwarz (NOAA)"
            gAtt@reference     = "Schwarz et al., J. Geophys. Res., 111," + \
                " D16207, doi:10.1029/2006JD007076, 2006"
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
        delete(OUTPUT)
        delete(data)

    end do

    leave_msg(diag_script, "", 1)

end
