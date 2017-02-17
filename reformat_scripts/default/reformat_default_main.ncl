;;#############################################################################
;; DEFAULT REFORMAT ROUTINES
;; Author: Mattia Righi (DLR, Germany)
;;#############################################################################
;;
;; Description
;;    Checks and reformats model data to the CMOR standard.
;;
;; Caveats
;;
;; Modification history
;;    20150630-a_laue_ax: saving the following details to the output NetCDF
;;                        files (to be added to the ESMValTool log-files):
;;                        - ESMValTtool version
;;                        - list of input files processed
;;                        - name of the "FIXFILE" applied (if any)
;;    20150415-A_righ_ma: extended for ocean variables.
;;    20130528-A_righ_ma: written.
;;
;;#############################################################################

load "./interface_data/ncl.interface"
load "./interface_scripts/messaging.ncl"

scriptname = "reformat_scripts/default/reformat_default_main.ncl"

;;### GLOBAL SETTINGS #########################################################
begin

    ;; Selected variable (standard name)
    VAR  = getenv("ESMValTool_variable")

    ;; Selected field
    FIELD  = getenv("ESMValTool_field")

    ;; Time period (YEAR1-01 to YEAR2-12 will be processed)
    YEAR1 = toint(getenv("ESMValTool_start_year"))
    YEAR2 = toint(getenv("ESMValTool_end_year"))

    ;; Model name
    MODEL = getenv("ESMValTool_model")

    ;; Project
    PROJECT = getenv("ESMValTool_project")

    ;; Ensemble
    ENSEMBLE = getenv("ESMValTool_ensemble")

    ;; Input path
    INPATH = getenv("ESMValTool_indir_path")

    ;; Input file
    INFILE = getenv("ESMValTool_infile_path")

    ;; Output file
    OUTFILE = getenv("ESMValTool_outfile_fullpath")

    ;; Grid type (initialize)
    GRID = ""

    ;; Area file (for irregular grids only)
    AREAFILE = getenv("ESMValTool_areafile_path")

    ;; Mask file for land variables (regular grids)
    LMASKFILE = getenv("ESMValTool_lmaskfile_path")

    ;; Mask file for ocean variables (irregular grids)
    OMASKFILE = getenv("ESMValTool_omaskfile_path")

    ;; Porosity file for mrso variable
    POROFILE = getenv("ESMValTool_porofile_path")

    ;; Fix file (for model-specific fixes)
    FIXFILE = "./reformat_scripts/fixes/" + \
        getenv("ESMValTool_project_basename") + "_" + MODEL + "_fix.ncl"

end

;;### MAIN ####################################################################

load "./reformat_scripts/default/reformat_default_func.ncl"

begin

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, "", 2)

    info_output("INFILE = " + INFILE, verbosity, 1)

    ;; Check input path
    check_inpath()

    ;; Find grid type
    GRID = find_grid_type()
    info_output("GRID TYPE = " + GRID, verbosity, 2)

    ;; Read variable and extract given time range
    filelist = new(1000, string)  ; list of files processed by read_var
    fixfile  = new(1, string)  ; list of files processed by FIXFILE
    VARPROC = read_var(VAR, filelist, fixfile)  ; filelist, fixfile are output

    ;; Check rank/field
    tmp = check_rank(VARPROC)
    delete(VARPROC)
    VARPROC = tmp
    delete(tmp)

    ;; Check fill value
    check_fill(VARPROC)

    ;; Reformat all coordinates
    tmp = reformat_coord(VARPROC)
    delete(VARPROC)
    VARPROC = tmp
    delete(tmp)

    ;; Check units
    check_units(VARPROC, VAR)

    ;; Set variable attributes (CMOR)
    tmp = var_attrib(VARPROC, VAR)
    delete(VARPROC)
    VARPROC = tmp
    delete(tmp)

    ;; Write output
    version  = getenv("0_ESMValTool_version")
    write_output(VARPROC, VAR, filelist, fixfile, version)
    delete(fixfile)
    delete(filelist)

    ;; Add latitude and longitude information (for irregular grids only)
    if (GRID.eq."irregular") then
        add_irregulargrid_info(VARPROC)
    end if

    ;; Add land-mask information
    if (isStrSubset(INFILE, "Lmon")) then
        add_lmask_info(VARPROC)
    end if

    ;; Add ocean-mask information
    if (isStrSubset(INFILE, "Omon")) then
        add_omask_info(VARPROC)
    end if

    ;; Add porosity information
    if (VAR.eq."mrso") then
        add_porosity_info(VARPROC)
    end if

    ;; Add depth information
    if (VAR.eq."mrsos") then
        add_depth_info(VARPROC)
    end if

    ;; Add region information (for ocean variables using the basin coord)
    if (iscoord(VARPROC, "basin")) then
        add_regionbasin_info()
    end if

    delete(VARPROC)
    info_output("OUTFILE = " + OUTFILE, verbosity, 1)

    leave_msg(scriptname, "", 2)

end
