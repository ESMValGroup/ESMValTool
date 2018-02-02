;;#############################################################################
;; REFORMAT SCRIPT FOR Putaud OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    Asmi et al., Atmos. Chem. Phys., 11, 5505, doi:10.5194/acp-11-5505-2011,
;;    2011.
;;
;; Last access
;;    Size distributions fits as published in the source.
;;
;; Download and processing instructions
;;
;; Caveats
;;
;; Modification history
;;    20161019-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "Asmi11"

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
    NAMES = "sizecnSTP"

    ;; Stations (Table 2 in Asmi et al. 2011)
    STATIONS = (/"ASP", "BIR", "PAL", "PLA", "SMR", \
                 "VHL", "BOE", "KPO", "MPZ", \
                 "OBK", "HPB", "WAL", "CBW", \
                 "HWL", "MHD", "FKL", "ISP", \
                 "ZEP", "PDD", "SCH", "ZSF", \
                 "JFJ", "BEO", "CMN"/)
    STATIONS@name = (/"Aspvreten", "Birkenes", "Pallas", "Preila", "SMEAR II",\
                      "Vavihill", "Boesel", "K-Puszta", "Melpitz", \
                      "Kosetice", "Hohenpeissenberg", "Waldhof", \
                      "Cabauw", "Harwell", "Mace Head", "Finokalia", \
                      "JRC-Ispra", "Zeppelin", "Puy de Dome", "Schauinsland", \
                      "Zugspitze", "Jungfraujoch", "BEO Moussala", \
                      "Monte Cimone"/)
    STATIONS@latitude = (/58.80, 58.38, 67.97, 55.92, 61.85, \
                         56.02, 53.00, 46.97, 51.53, \
                         49.58, 47.80, 52.52, 51.30, \
                         51.57, 53.32, 35.33, 45.82, \
                         78.92, 45.77, 47.92, 47.42, \
                         46.53, 42.17, 44.18/)
    STATIONS@longitude = (/17.38, 8.25, 24.12, 21.00, 24.28, \
                          13.15, 57.95, 19.55, 12.90, \
                          15.08, 11.00, 10.77, 4.92, \
                          358.68, 350.12, 25.67, 8.63, \
                          11.90, 2.95, 7.92, 10.98, \
                          7.98, 23.58, 10.68/)
    STATIONS@altitude = (/30., 190., 560., 5., 181., \
                         172., 16., 125., 87., \
                         534., 988., 70., 60., \
                         60., 5., 250., 209., \
                         474., 1465., 1210., 2650., \
                         3580., 2971., 2165/)

    NSTATIONS = dimsizes(STATIONS)

    ;; Modes
    MODES = (/"Mode 1", "Mode 2"/)
    NMODES = dimsizes(MODES)

    ;; Data array
    DATA = new((/NSTATIONS, NMODES, 3/), float)

end

load "./interface_scripts/logging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

begin

    verbosity = config_user_info@verbosity
    diag_script = "reformat_obs_Asmi11.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Define data (Table 6 in Asmi et al. 2011)

    ;; Mode 1                           N[cm-3] CMD[um] sigma
    DATA(ind(STATIONS.eq."ASP"), 0, :) = (/1233, 64, 1.87/)
    DATA(ind(STATIONS.eq."BIR"), 0, :) = (/ 539, 53, 1.63/)
    DATA(ind(STATIONS.eq."PAL"), 0, :) = (/ 221, 48, 1.71/)
    DATA(ind(STATIONS.eq."PLA"), 0, :) = (/3175, 78, 2.04/)
    DATA(ind(STATIONS.eq."SMR"), 0, :) = (/1008, 57, 1.96/)
    DATA(ind(STATIONS.eq."VHL"), 0, :) = (/1756, 56, 1.86/)
    DATA(ind(STATIONS.eq."BOE"), 0, :) = (/ 432, 21, 1.36/)
    DATA(ind(STATIONS.eq."KPO"), 0, :) = (/2788, 54, 2.20/)
    DATA(ind(STATIONS.eq."OBK"), 0, :) = (/2479, 51, 2.26/)
    DATA(ind(STATIONS.eq."HPB"), 0, :) = (/1609, 45, 2.39/)
    DATA(ind(STATIONS.eq."MPZ"), 0, :) = (/3516, 58, 2.49/)
    DATA(ind(STATIONS.eq."WAL"), 0, :) = (/3837, 58, 2.37/)
    DATA(ind(STATIONS.eq."CBW"), 0, :) = (/1721, 35, 1.40/)
    DATA(ind(STATIONS.eq."HWL"), 0, :) = (/2196, 46, 2.06/)
    DATA(ind(STATIONS.eq."MHD"), 0, :) = (/ 303, 35, 1.79/)
    DATA(ind(STATIONS.eq."FKL"), 0, :) = (/1359, 86, 1.94/)
    DATA(ind(STATIONS.eq."ISP"), 0, :) = (/2828, 31, 2.19/)
    DATA(ind(STATIONS.eq."ZEP"), 0, :) = (/  50, 39, 1.65/)
    DATA(ind(STATIONS.eq."PDD"), 0, :) = (/1223, 48, 2.25/)
    DATA(ind(STATIONS.eq."SCH"), 0, :) = (/1709, 51, 2.32/)
    DATA(ind(STATIONS.eq."ZSF"), 0, :) = (/ 729, 45, 2.29/)
    DATA(ind(STATIONS.eq."JFJ"), 0, :) = (/ 221, 47, 1.70/)
    DATA(ind(STATIONS.eq."BEO"), 0, :) = (/ 514, 60, 1.76/)
    DATA(ind(STATIONS.eq."CMN"), 0, :) = (/ 716, 52, 1.98/)

    ;; Mode 2                           N[cm-3] CMD[um] sigma
    DATA(ind(STATIONS.eq."ASP"), 1, :) = (/ 144, 222, 1.36/)
    DATA(ind(STATIONS.eq."BIR"), 1, :) = (/ 142, 170, 1.41/)
    DATA(ind(STATIONS.eq."PAL"), 1, :) = (/  81, 176, 1.45/)
    DATA(ind(STATIONS.eq."PLA"), 1, :) = (/  84, 207, 1.27/)
    DATA(ind(STATIONS.eq."SMR"), 1, :) = (/ 144, 201, 1.36/)
    DATA(ind(STATIONS.eq."VHL"), 1, :) = (/ 253, 195, 1.41/)
    DATA(ind(STATIONS.eq."BOE"), 1, :) = (/4710,  57, 2.44/)
    DATA(ind(STATIONS.eq."KPO"), 1, :) = (/1590, 134, 1.73/)
    DATA(ind(STATIONS.eq."OBK"), 1, :) = (/2019, 131, 1.78/)
    DATA(ind(STATIONS.eq."HPB"), 1, :) = (/ 554, 148, 1.66/)
    DATA(ind(STATIONS.eq."MPZ"), 1, :) = (/ 338, 203, 1.60/)
    DATA(ind(STATIONS.eq."WAL"), 1, :) = (/ 175, 194, 1.38/)
    DATA(ind(STATIONS.eq."CBW"), 1, :) = (/4450,  60, 1.76/)
    DATA(ind(STATIONS.eq."HWL"), 1, :) = (/ 256, 196, 1.55/)
    DATA(ind(STATIONS.eq."MHD"), 1, :) = (/ 152, 174, 1.70/)
    DATA(ind(STATIONS.eq."FKL"), 1, :) = (/ 228, 189, 1.41/)
    DATA(ind(STATIONS.eq."ISP"), 1, :) = (/4344,  91, 1.87/)
    DATA(ind(STATIONS.eq."ZEP"), 1, :) = (/ 109, 150, 1.61/)
    DATA(ind(STATIONS.eq."PDD"), 1, :) = (/ 207, 147, 1.54/)
    DATA(ind(STATIONS.eq."SCH"), 1, :) = (/ 322, 166, 1.48/)
    DATA(ind(STATIONS.eq."ZSF"), 1, :) = (/ 132, 157, 1.49/)
    DATA(ind(STATIONS.eq."JFJ"), 1, :) = (/  84, 136, 1.46/)
    DATA(ind(STATIONS.eq."BEO"), 1, :) = (/ 330, 156, 1.53/)
    DATA(ind(STATIONS.eq."CMN"), 1, :) = (/ 237, 157, 1.50/)

    ;; Convert units
    ;; Warning: Asmi et al. (2011) use STP concentrations and define STP with
    ;; T = 293.15 K, while in the ESMValTool T = 273.15 K. This difference is
    ;; accounted for in the conversion below
    DATA(:, :, 0) = DATA(:, :, 0) * 1.e6  ; [cm-3]--> [m-3]
    DATA(:, :, 0) = DATA(:, :, 0) * 293.15 / 273.15  ; [STP_293] --> [STP_273]
    DATA(:, :, 1) = DATA(:, :, 1) * 1.e-9  ; [nm] --> [m]

    ;; Define diameters
    DIAM = set_size_array()
    NDIAM = dimsizes(DIAM)
    delete_VarAtts(DIAM, -1)

    ;; Define output array
    NCOLS = 1
    OUTPUT = new((/NDIAM, NCOLS/), float)
    OUTPUT!0 = "diam"
    OUTPUT!1 = "column"
    OUTPUT&diam = DIAM
    OUTPUT&column = (/"median"/)
    OUTPUT&diam@long_name     = "diameter"
    OUTPUT&diam@units         = "m"
    OUTPUT&diam@standard_name = "particle_dry_diameter"
    OUTPUT = 0.

    ;; Loop over regions
    do stID = 0, NSTATIONS - 1

        output = OUTPUT

        ;; Calculate and add up size distributions for each mode
        do mm = 0, NMODES - 1
            size = lognormal_dist(DATA(stID, mm, 0), \
                                  DATA(stID, mm, 1), \
                                  DATA(stID, mm, 2), DIAM)
            if (.not.all(ismissing(size))) then
                output(:, 0) = output(:, 0) + size
            end if
            delete(size)
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
        gAtt@lonmax   = where(STATIONS@longitude(stID).lt.0., \
                              STATIONS@longitude(stID) + 360., \
                              STATIONS@longitude(stID))
        gAtt@lonmin   = gAtt@lonmax
        gAtt@latmax   = STATIONS@latitude(stID)
        gAtt@latmin   = gAtt@latmax
        gAtt@info     = STATIONS(stID)
        gAtt@month    = ispan(1, 12, 1)
        gAtt@field    = FIELD
        gAtt@tier     = TIER
        gAtt@source = "Size distributions fits as published in the source."
        gAtt@reference = "Asmi et al., Atmos. Chem. Phys., 11, 5505-5538, " + \
            "doi:10.5194/acp-11-5505-2011, 2011."
        gAtt@title    = OBSNAME + \
            " particle size distributions reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"

        if (.not.all(ismissing(output))) then

            ;; Outfile
            fout = OUTDIR + "OBS_" + OBSNAME + "_ground_" + STATIONS(stID) + \
                "_" + FIELD + "_" + NAMES + "_200801-200912.nc"

            ;; Write variable
            write_nc_profile(fout, NAMES, output, gAtt)
            delete(output)

        end if
    end do

    leave_msg(diag_script, "", 1)

end
