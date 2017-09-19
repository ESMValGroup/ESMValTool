;;#############################################################################
;; REFORMAT SCRIPT FOR CR-AVE OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 3: restricted dataset.
;;
;; Source
;;    Schwarz et al., J. Geophys. Res., 113, D03203,
;;    doi:10.1029/2007JD009042, 2008.
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
    OBSNAME = "CR-AVE"

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
    NAME = "mmrbc"

    ;; Campaigns
    CAMPAIGN = "1"
    PERIOD   = "2006/02/06-2006/02/09"
    MINLAT   = -1.
    MAXLAT   = 11.
    MINLON   = 275.
    MAXLON   = 281.

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_CR-AVE.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    infile = INDIR + OBSNAME + "_200702_679_MMR.txt"

    ;; Read data
    data = readAsciiTable(infile, 1, "string", 5)
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

    ;; Read pressure coordinate and convert units [hPa] --> [Pa]
    OUTPUT&plev = todouble(str_get_field(data(:, 0), 1, " "))
    OUTPUT&plev = OUTPUT&plev * 100.

    ;; Read mean and convert units [ng/kg] --> [kg/kg]
    OUTPUT(:, 0) = \
        tofloat(str_get_field(data(:, 0), 3, " "))
    OUTPUT(:, 0) = OUTPUT(:, 0) * 1.e-12

    ;; Read standard deviation and convert units [ng/kg] --> [kg/kg]
    OUTPUT(:, 1) = \
        tofloat(str_get_field(data(:, 0), 4, " "))
    OUTPUT(:, 1) = OUTPUT(:, 1) * 1.e-12

    ;; Set variable attributes (CMOR)
    tmp = var_attrib(OUTPUT, NAME)
    delete(OUTPUT)
    OUTPUT = tmp
    delete(tmp)

    ;; Append global attributes
    gAtt = True
    gAtt@comment     = "The data is adjusted by a factor 1.1 to account " + \
        "for the undetected fine-mode BC mass"
    gAtt@history     = "Created on " + systemfunc("date")
    gAtt@host        = systemfunc("echo $HOST")
    gAtt@user        = systemfunc("echo $USER")
    gAtt@lonmax      = where(MAXLON.lt.0., MAXLON + 360., MAXLON)
    gAtt@lonmin      = where(MINLON.lt.0., MINLON + 360., MINLON)
    gAtt@latmax      = MAXLAT
    gAtt@latmin      = MINLAT
    gAtt@info        = ""
    gAtt@period      = PERIOD
    gAtt@field       = FIELD
    gAtt@tier        = TIER
    gAtt@source      = "Data provided by Joshua Schwarz (NOAA)"
    gAtt@reference   = "Schwarz et al., J. Geophys. Res., 113," + \
        " D03203, doi:10.1029/2007JD009042, 2008"
    gAtt@title       = OBSNAME + \
        " vertical profiles reformatted for the ESMValTool"
    gAtt@conventions = "CF/CMOR"

    ;; Outfile
    d1 = str_get_field(str_get_field(PERIOD, 1, "-"), 1, "/") + \
        str_get_field(str_get_field(PERIOD, 1, "-"), 2, "/")
    d2 = str_get_field(str_get_field(PERIOD, 2, "-"), 1, "/") + \
        str_get_field(str_get_field(PERIOD, 2, "-"), 2, "/")
    fout = OUTDIR + "OBS_" + OBSNAME + "_campaign_" + \
        CAMPAIGN + "_" + FIELD + "_" + NAME + "_" + d1 + "-" + d2 + ".nc"

    ;; Write variable
    write_nc_profile(fout, NAME, OUTPUT, gAtt)
    delete(OUTPUT)
    delete(data)

    leave_msg(diag_script, "", 1)

end
