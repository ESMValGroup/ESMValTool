;;#############################################################################
;; REFORMAT ROUTINES FOR THE ESMVal TOOL FOR ECEARTH
;; Author: Bert van Ulft (KNMI, the Netherlands)
;; EMBRACE project
;;#############################################################################
;;
;; Description
;;    Tool to reformat EC-Earth output so it conforms to the tool's
;;    requirements. Written for the ocean part.
;;
;; Modification history
;;    20151029-A_laue_ax: changed dimension of file list from 1,000 to 10,000
;;    20141205-A_vanu_be: written.
;;
;;#############################################################################

;; Load Temporary file with current namelist configuration
load "./interface_scripts/logging.ncl"

scriptname = "reformat_scripts/default/reformat_ECEARTH_main.ncl"

;;### GLOBAL SETTINGS #########################################################
begin

    enter_msg(scriptname, "")

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

    ;; grid description files
    HGRIDFILE = getenv("ESMValTool_hgridfile_path")
    ZGRIDFILE = getenv("ESMValTool_zgridfile_path")
    LSMFILE   = getenv("ESMValTool_lsmfile_path")

    ;; Fix file (for model-specific fixes)
    FIXFILE = "./reformat_scripts/fixes/" + \
              getenv("ESMValTool_project_basename") + "_" + MODEL + "_fix.ncl"

end

;;### MAIN ####################################################################

;; use functions from default reformat scripts whenever possible
load "./reformat_scripts/default/reformat_default_func.ncl"
load "./reformat_scripts/ECEARTH/reformat_ECEARTH_func.ncl"

begin

    ;; Check field type is supported
    if (.not. isStrSubset(FIELD, "T3") .and. .not. \
        (isStrSubset(FIELD, "T") .and. isStrSubset(FIELD, "2") .and. \
         isStrSubset(FIELD, "s"))) then
        error_msg("f", scriptname, "", \
                  "variable " + VAR + " has an unsupported field " + \
                  "type: " + FIELD + ".")
    end if

    ;; Find EC-Earth name, has grid (U, V, T, W) as attribute
    ECEARTH_VAR = find_name(VAR, "ECEARTH")

    ;; add staggered grid type to input filename
    sub_staggergrid(ECEARTH_VAR)

    log_info("  INFILE = " + INFILE)

    ;; Check input path
    check_inpath()

    ;; set grid type
    GRID = "irregular"
    log_debug("GRID TYPE set to " + GRID)

    ;; Read variable and extract given time range
    filelist = new(10000, string)
    fixfile  = new(1, string)
    VARPROC = read_var(ECEARTH_VAR, filelist, fixfile)

    ;; Check rank/field
    tmp = check_rank(VARPROC)
    delete(VARPROC)
    VARPROC = tmp
    delete(tmp)

    ;; Check fill value
    check_fill(VARPROC)

    ;; Filter out landpoints
    mask_land(VARPROC)

    ;; Rename coordinate variables so they can be reformatted with the default
    ;; reformat_coord routine
    rename_time(VARPROC, ECEARTH_VAR)
    rename_lev(VARPROC, ECEARTH_VAR)

    ;; add i and j dims and vars and assign them as coordinate variables
    add_ijcoordvars(VARPROC)

    ;; Reformat all coordinates
    tmp = reformat_coord(VARPROC)
    delete(VARPROC)
    VARPROC = tmp
    delete(tmp)

    ;; Convert (if required) and check units
    convert_units(VARPROC, VAR)
    check_units(VARPROC, VAR)

    ;; Set variable attributes (CMOR)
    tmp = var_attrib(VARPROC, VAR)
    delete(VARPROC)
    VARPROC = tmp
    delete(tmp)

    ;; Write output
    version  = diag_script_info@version
    write_output(VARPROC, VAR, filelist, fixfile, version)

    ;; Add latitude and longitude
    add_ECEARTHgrid_info(VARPROC, ECEARTH_VAR)

    ;; clean up
    delete(VARPROC)

    log_info("  OUTFILE = " + OUTFILE)

    leave_msg(scriptname, "")

end
