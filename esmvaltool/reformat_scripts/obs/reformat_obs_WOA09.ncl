;;#############################################################################
;; REFORMAT SCRIPT FOR WORLD OCEAN ATLAS 09 MONTHLY CLIMATOLOGY DATA
;;#############################################################################
;;
;; Source
;;  ftp://ftp.nodc.noaa.gov/pub/WOA09/NetCDFdata
;;
;; Last access
;;    20150521
;;
;; Download and processing instructions
;;    Download following files with any ftp program:
;;    ftp://ftp.nodc.noaa.gov/pub/WOA09/NetCDFdata/salinity_monthly_1deg.nc
;;    ftp://ftp.nodc.noaa.gov/pub/WOA09/NetCDFdata/temperature_monthly_1deg.nc
;;
;;    This script will rename and modify the original netcdf files so they are
;;    accepted by the ESMValTool as T3M fields.
;     Furthermore it will create TO2Ms files with the top layer.
;;
;; Modification history
;;    20150915-A_laue_ax: added global attributes TIER and REFERENCE
;;    20150521-A_vanu_be: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "WOA09-monthly"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; variable names as in file name, also used for long and standard names in
    ;; surface files
    VARSFN = (/"salinity", "temperature"/)

    ;; Selected variable (standard name)
    VARS = (/"so", "to"/)

    ;; Name in the raw data
    NAME = (/"s_an", "t_an"/)

    ;; Units in the input  data and corresponding conversion factor or addition
    ;; to CMOR standard units
    RAWUNITS = (/"1", "degrees Celsius"/)
    CONVERSION = (/1., 1./)
    ADDITION = (/0., 273.15/)
end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_WOA09.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not. isfilepresent_esmval(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Loop over variables
    do vID = 0, dimsizes(VARS) - 1
        info_output("", verbosity, 1)
        info_output("Processing " + VARS(vID), verbosity, 1)
        FIELD = "T3M"

        ;; open file
        fname = INDIR + VARSFN(vID) + "_monthly_1deg.nc"
        f = addfile(fname, "r")

        ;; Read variable
        output = f->$NAME(vID)$

        ;; Rename vertical dimension
        output!1 = "lev"

        ;; Convert units
        if (output@units.ne.RAWUNITS(vID)) then
            error_msg("f", get_script_name(), "", \
                      "Unexpected units for " + NAME(vID))
        else
            output = output * CONVERSION(vID) + ADDITION(vID)
        end if

        ;; Set variable attributes
        tmp = var_attrib(output, VARS(vID))
        delete(output)
        output = tmp
        delete(tmp)

        ;; Set global attributes
        gAtt = True
        gAtt@history       = "Created on " + systemfunc("date")
        gAtt@host          = systemfunc("echo $HOST")
        gAtt@user          = systemfunc("echo $USER")
        gAtt@period        = "climatology"
        gAtt@field         = FIELD
        gAtt@tier          = TIER
        gAtt@source        = "ftp://ftp.nodc.noaa.gov/pub/WOA09/NetCDFdata"
        gAtt@inputfile     = fname
        gAtt@reference     = \
            "Locarnini, R. A., et al., 2010. World Ocean Atlas 2009, " + \
            "Volume 1: Temperature. S. Levitus, Ed., " + \
            "NOAA Atlas NESDIS 68, U.S. Gov. Printing Office, " + \
            "Washington, D.C., 184 pp.;" + \
            "Antonov, J. I., et al., 2010. World Ocean Atlas 2009, " + \
            "Volume 2: Salinity. S. Levitus Ed., " + \
            "NOAA Atlas NESDIS 69, U.S. Gov. Printing Office, " + \
            "Washington, D.C., 184 pp."
        gAtt@title = OBSNAME + \
            " climatology data for the ESMVal evaluation tool"

        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_insitu_1_" + FIELD + \
            "_" + VARS(vID) + "_000001-000012.nc"

        ;; Write variable
        write_nc(fout, VARS(vID), output, gAtt)
        delete(fout)

        ;; Extract the top layer as a surface field
        FIELD = "TO2Ms"
        sfc = output(:, 0, :, :)

        ;; Set variable attributes
        sfc@long_name = "Sea Surface " + str_capital(VARSFN(vID))
        sfc@standard_name = "sea_surface_" + VARSFN(vID)

        ;; modify global attributes
        gAtt@field         = FIELD

        ;; Write variable
        fout = OUTDIR + "OBS_" + OBSNAME + "_insitu_1_" + FIELD + \
            "_" + VARS(vID) + "s_000001-000012.nc"
        write_nc(fout, VARS(vID) + "s", sfc, gAtt)
        delete(fout)

        ;; clean
        delete(gAtt)
        delete(output)
        delete(sfc)
        delete(f)
    end do

    leave_msg(diag_script, "", 1)

end
