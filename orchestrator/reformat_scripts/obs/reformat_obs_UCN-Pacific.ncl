;;#############################################################################
;; REFORMAT SCRIPT FOR UCN-Pacific OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    Clarke and Kapustin, J. Atmos. Sci., 59, 363-382
;;    doi:10.1175/1520-0469(2002)059<0363:APASPI>2.0.CO;2, 2002.
;;
;; Last access
;;    Data provided by Axel Lauer (extracted from Fig. 9 in the source).
;;
;; Download and processing instructions
;;
;; Caveats
;;
;; Modification history
;;    20140508-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "UCN-Pacific"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Field
    FIELD = "T0M"

    ;; Variable
    VARIABLE = "conccnSTPd3"

    ;; Campaigns (see INDIR/readme.txt)
    CAMPAIGN = (/"70S-20S", "20S-20N", "20N-70N"/)
    CODE     = str_lower(CAMPAIGN)
    MINLAT   = (/-70., -20., 20./)
    MAXLAT   = (/-20., 20., 70./)
    MINLON   = (/130., 130., 130./)
    MAXLON   = (/266., 288., 288./)

    ;; Selected months (apply to all campaigns!)
    MONTH    = (/3, 5, 9, 11/)

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_UCN-Pacific.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Generate data array
    NLEVS = 13
    NCOLS = 2
    OUTPUT = new((/NLEVS, NCOLS/), float)
    OUTPUT!0 = "plev"
    OUTPUT!1 = "column"
    OUTPUT&column = (/"mean", "stddev"/)
    OUTPUT&plev = todouble(ispan(1, NLEVS, 1))  ; dummy
    OUTPUT&plev@long_name     = "pressure"
    OUTPUT&plev@positive      = "down"
    OUTPUT&plev@axis          = "Z"
    OUTPUT&plev@units         = "Pa"
    OUTPUT&plev@standard_name = "air_pressure"

    ;; Loop over campaigns
    do cpID = 0, dimsizes(CAMPAIGN) - 1

        infile = INDIR + "clarke" + "_" + CODE(cpID) + ".asc"
        info_output(infile, verbosity, 1)

        ;; Read data
        data = readAsciiTable(infile, 1, "string", 0)

        ;; Copy metadata
        output = OUTPUT

        ;; Read level coordinate
        output&plev = todouble(str_get_field(data(:, 0), 1, " "))
        output&plev = 100. * output&plev

        ;; Read mean
        output(:, 0) = tofloat(str_get_field(data(:, 0), 2, " "))

        ;; Read stddev
        output(:, 1) = tofloat(str_get_field(data(:, 0), 3, " "))

        ;; Convert units
        output = output * 1.e6  ; [1/cm3-STP] --> [1/m3-STP]

        ;; Set variable attributes (CMOR)
        tmp = var_attrib(output, VARIABLE)
        delete(output)
        output = tmp
        delete(tmp)

        ;; Append global attributes
        gAtt = True
        gAtt@comment       = "Altitude conversion based on " + \
            "standard pressures"
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
        gAtt@month         = MONTH
        gAtt@field         = FIELD
        gAtt@tier          = TIER
        gAtt@source        = "Data as published in the source."
        gAtt@reference     = "Clarke and Kapustin, J. Atmos. Sci., " + \
            "59, 363-382 doi:10.1175/1520-0469(2002)" + \
            "059<0363:APASPI>2.0.CO;2, 2002"
        gAtt@title         = OBSNAME + \
            " vertical profiles reformatted for the ESMValTool"
        gAtt@conventions   = "CF/CMOR"

        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_campaign_" + \
            CAMPAIGN(cpID) + "_" + FIELD + "_" + VARIABLE + "_199005-199903.nc"

        ;; Write variable
        write_nc_profile(fout, VARIABLE, output, gAtt)

        delete(output)

    end do

    leave_msg(diag_script, "", 1)

end
