;;#############################################################################
;; REFORMAT SCRIPT FOR HIPPO OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 3: restricted dataset.
;;
;; Source
;;    Schwarz et al., Geophys. Res. Lett., 40, 5542-5547,
;;    doi:10.1002/2013GL057775, 2013.
;;
;; Last access
;;    Data provided by J. Schwarz.
;;
;; Download and processing instructions
;;    Data has been converted with mac2unix and untabified.
;;    nan values have been manually added to the last row to fill empty spaces.
;;
;; Caveats
;;
;; Modification history
;;    20140910-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "HIPPO"

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

    ;; Period (consider full year, as data represents a climatology)
    PERIOD   = "2009/01/01-2011/12/31"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_HIPPO.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    infile = INDIR + "AllHIPPO_latBand_Observations.txt"

    ;; Read header and set different regions/campaigns
    head = readAsciiHead(infile, 1)
    ncols = str_fields_count(head, " ")
    NCAMPAIGNS  = ncols / 3
    CAMPAIGN = new(NCAMPAIGNS, string)
    delete(ncols)

    ;; Read data
    data = readAsciiTable(infile, 1, "string", 1)
    data = str_sub_str(data, "nan", "-999")
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

    ;; Loop over campaigns
    do cpID = 0, NCAMPAIGNS - 1

        ;; Define campaign region and name
        tmp = str_get_field(head, 1 + 3 * cpID, " ")
        sign = where(str_get_cols(tmp, 10, 10).eq."N", 1., -1.)
        MINLAT = sign * tofloat(str_get_cols(tmp, 11, 12))
        sign = where(str_get_cols(tmp, 14, 14).eq."N", 1., -1.)
        MAXLAT = sign * tofloat(str_get_cols(tmp, 15, 16))
        CAMPAIGN = abs(MINLAT) + str_get_cols(tmp, 10, 10) + \
            abs(MAXLAT) + str_get_cols(tmp, 14, 14)
        delete(tmp)
        delete(sign)
        MINLON = 160.  ; same for all
        MAXLON = 240.  ; same for all

        info_output(CAMPAIGN(cpID), verbosity, 1)

        ;; Copy metadata
        output = OUTPUT

        ;; Read pressure coordinate and convert units [hPa] --> [Pa]
        output&plev = todouble(str_get_field(data(:, 0), 2 + 3 * cpID, " "))
        output&plev = output&plev * 100.

        ;; Read mean and convert units [ng/kg] --> [kg/kg]
        output(:, 0) = \
            tofloat(str_get_field(data(:, 0), 1 + 3 * cpID, " "))
        output(:, 0) = output(:, 0) * 1.e-12

        ;; Read standard deviation and convert units [ng/kg] --> [kg/kg]
        output(:, 1) = \
            tofloat(str_get_field(data(:, 0), 3 + 3 * cpID, " "))
        output(:, 1) = output(:, 1) * 1.e-12

        ;; Set missing values
        output = where(output.le.0., output@_FillValue, output)

        ;; Filter valid values
        validx = ind(output&plev.ge.0.)
        valid = output(validx, :)
        delete(output)
        delete(validx)

        ;; Set variable attributes (CMOR)
        tmp = var_attrib(valid, NAME)
        delete(valid)
        valid = tmp
        delete(tmp)

        ;; Append global attributes
        gAtt = True
        gAtt@history     = "Created on " + systemfunc("date")
        gAtt@host        = systemfunc("echo $HOST")
        gAtt@user        = systemfunc("echo $USER")
        gAtt@lonmax      = where(MAXLON.lt.0., MAXLON + 360., MAXLON)
        gAtt@lonmin      = where(MINLON.lt.0., MINLON + 360., MINLON)
        gAtt@latmax      = MAXLAT
        gAtt@latmin      = MINLAT
        gAtt@info        = CAMPAIGN(cpID)
        gAtt@period      = PERIOD
        gAtt@field       = FIELD
        gAtt@tier        = TIER
        gAtt@source      = "Data provided by Joshua Schwarz (NOAA)"
        gAtt@reference   = "Schwarz et al., Geophys. Res. Lett., 40, " + \
            "5542-5547, doi:10.1002/2013GL057775, 2013"
        gAtt@title       = OBSNAME + \
            " vertical profiles reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"

        ;; Outfile
        d1 = str_get_field(str_get_field(PERIOD, 1, "-"), 1, "/") + \
            str_get_field(str_get_field(PERIOD, 1, "-"), 2, "/")
        d2 = str_get_field(str_get_field(PERIOD, 2, "-"), 1, "/") + \
            str_get_field(str_get_field(PERIOD, 2, "-"), 2, "/")
        fout = OUTDIR + "OBS_" + OBSNAME + "_campaign_" + \
            CAMPAIGN(cpID) + "_" + FIELD + "_" + NAME + \
            "_" + d1 + "-" + d2 + ".nc"

        ;; Write variable
        write_nc_profile(fout, NAME, valid, gAtt)
        delete(valid)

    end do
    delete(data)

    leave_msg(diag_script, "", 1)

end
