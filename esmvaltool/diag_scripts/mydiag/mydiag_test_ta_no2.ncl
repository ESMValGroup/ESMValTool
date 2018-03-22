;##############################################################################
;# SANDBOX FOR IMPLEMENTING NEW DIAGNOSTICS                                   #
;##############################################################################
;# Author: Klaus-Dirk Gottschaldt (DLR, Germany)                              #
;# ESMVal project
;##############################################################################
;#                                                                            #
;# Description:                                                               #
;#    This script is a driver for averaging the input from the first model    #
;#    over time and creating a global contour plot.                           #
;#    Output to NetCDF files is optional (controlled via diag_script_info)    #
;#    Time averaging and plotting is done in separate routines.               #
;#    All information necessary for plotting is added as attributes to the    #
;#    variable containing the data. Attribute names "diag_script", "var"      #
;#    and those starting with "res_" are reserved for that purpose.           #
;#    The data variable (incl. metadata) can be passed to the plot script     #
;#    directly (most efficient), or via a netCDF file (might be useful for    #
;#    non-NCL routines).                                                      #
;#    A workstation ("wks" = file to plot to) will by default be generated    #
;#    in the plotting routine. However, the default will be overridden, if    #
;#    a "graphic" object is passed to the plotting routine. The latter        #
;#    option is intended to be used in special cases, e.g. with               #
;#    res_gsnDraw = False and res_gsnFrame = False.                           #
;#    The script has been tested with input from ./variable_defs/MyVar.ncl    #
;#                                                                            #
;# Required diag_script_info attributes (diagnostics specific):               #
;#    projection: map projection to use                                       #
;#    netCDF_file: full path and name for *.nc output                         #
;#                                                                            #
;# Required variable_info attributes (variable specific):                     #
;#    MyDiag_title: string for annotation                                     #
;#                                                                            #
;# Caveats:                                                                   #
;#                                                                            #
;# Modification history:                                                      #
;#    20151027-A_laue_ax: moved call to 'write_references' to the beginning   #
;#                        of the code                                         #
;#    20150325-A_laue_ax: modified reference tags used for acknowledgements   #
;#                        (projects, observations, etc)                       #
;#    20140312-A_righ_ma: adapted to new time_operations structure.           #
;#    20130618-A_gott_kl: adjusted to r738 of                                 #
;#                     https://svn.dlr.de/ESM-Diagnostic/sources/trunk        #
;#    20121130-A_gott_kl: written.                                            #
;#                                                                            #
;##############################################################################

load "interface_scripts/interface.ncl"

load "./diag_scripts/shared/statistics.ncl"

load "./diag_scripts/shared/plot/style.ncl"
load "./diag_scripts/shared/plot/contour_maps.ncl"


begin
    ;;###########################################################
    ;;# Fetch general parameters, set in namelist_MyDiag.yml,   #
    ;;# passed via environment variables by python code         #
    ;;###########################################################
    enter_msg(diag_script, "")

    ;; 'diag_script' and 'variables' are fetched from the
    ;; above 'ncl.interface' file
    log_info("++++++++++++++++++++++++++++++++++++++++++")
    log_info(diag_script + " (var: " + variables(0) + ")")
    log_info("++++++++++++++++++++++++++++++++++++++++++")

    ;; 'models', 'variables' and 'field_types' are fetched from the
    ;; above 'ncl.interface' file
    dim_MOD = dimsizes(model_info@model)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    ;;###########################################
    ;;# References Information                  #
    ;;###########################################
    write_references(diag_script, \  ; script name
                     "A_gott_kl", \  ; authors
                     "A_eval_ma", \  ; contributors
                     "D_0000",    \  ; diag_references
                     "",          \  ; obs_references
                     (/"P_dlrveu", "P_esmval"/))   ; proj_references

end

begin
    ;;#############################################
    ;;# Get parameters from ./variable_defs/*.ncl #
    ;;# passed via the 'info' attribute           #
    ;;#############################################
    if(isvar("MyParam")) then
        delete(MyParam)
    end if
    ;;if(isatt(variable_info, "MyDiag_title")) then
    ;;    MyParam = variable_info@MyDiag_title
    ;;    log_debug(" MyParam = " + MyParam)
    ;;else
    ;;    error_msg("f", diag_script, "", "missing variable attribute " + \
    ;;              "'MyDiag_title'")
    ;;end if

    ;;###########################################
    ;;# Get data and average time               #
    ;;###########################################
    ;; get data from first model
    imod = 0  ; NCL array indicies start from zero
    log_debug("processing " + model_info@model(imod))

    ;; See ./interface_scripts/data_handling.ncl
    A0 = read_data(imod, var0, field_type0)

    ;; Check dimensions
    dims = getvardims(A0)
    if (dimsizes(dims) .lt. 2) then
        error_msg("f", diag_script, "", dimsizes(dims) + " dimensions, " + \
                  "need 2 or 3")
    end if
    idx = ind(dims .eq. "lat")
    if (ismissing(idx)) then
        error_msg("f", diag_script, "", "no lat dimension")
    end if
    idx = ind(dims .eq. "lon")
    if (ismissing(idx)) then
        error_msg("f", diag_script, "", "no lon dimension")
    end if

    ;; Average over time (function in ./diag_scripts/shared/statistics.ncl)
    data1 = time_operations(A0, -1, -1, "average", "annualclim", True)

    ;;###########################################
    ;;# Open file for plot                      #
    ;;###########################################
    ;; If you really need to do this in the diagnostic sript, please see
    ;; http://www.ncl.ucar.edu/Document/Graphics/Interfaces/gsn_open_wks.shtml
    ;; and ~/diag_scripts/shared/plot/contour_maps.ncl (look for "wks") on how to do
    ;; this. However, wks should be defined in the plotting routine whenever
    ;; possible.

    ;;###########################################
    ;;# Style dependent annotation              #
    ;;###########################################
    ;; Retrieve unique strings describing the data
    ;; Function in ./diag_scripts/shared/plot/style.ncl
    annots = project_style(diag_script_info, "annots")

    ;;###########################################
    ;;# Plot ressources                         #
    ;;###########################################
    ;; Add to data1, as attributes with prefix "res_".
    ;; Provide attributes only that differ from the defaults given at
    ;; http://www.ncl.ucar.edu/Document/Graphics/Resources/list_alpha_res.shtml
    data1@res_gsnMaximize       = True            ; Use full page for the plot
    data1@res_cnFillOn          = True            ; Color plot desired
    data1@res_cnLineLabelsOn    = False           ; Contour lines
    data1@res_cnLinesOn         = False
    ;; Colors
    ;; http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
    data1@res_cnFillPalette     = diag_script_info@colormap  ; Set colormap
    ;; Annotation
    data1@res_tiMainString      = "Test diag"  ; Main title
    data1@res_gsnLeftString     = annots(imod)
    data1@res_gsnRightString    = variable_info@long_name + " in " \
                                                          + variable_info@units
    ;; Specified in ./nml/cfg_MyDiag/cfg_MyDiag.ncl
    data1@res_mpProjection      = diag_script_info@projection
    ;; Set explicit contour levels
    if(isatt(diag_script_info, "explicit_cn_levels")) then
        data1@res_cnLevelSelectionMode = "ExplicitLevels"
        data1@res_cnLevels = diag_script_info@explicit_cn_levels
    end if

    ;;###########################################
    ;;# Other Metadata: diag_script, var        #
    ;;###########################################
    ;; Add to data1, as attributes without prefix
    if(isatt(data1, "diag_script")) then  ; Add to existing entries
        temp = data1@diag_script
        delete(data1@diag_script)
        data1@diag_script = array_append_record(temp, (/diag_script/), 0)
        delete(temp)
    else  ; Add as new attribute
        data1@diag_script = (/diag_script/)
    end if
    data1@var = var0  ; Overwrite existing entry
    if(isatt(variable_info, "long_name")) then
        data1@var_long_name = variable_info@long_name
    end if
    if(isatt(variable_info, "units")) then
        data1@var_units     = variable_info@units
    end if

    ;;###########################################
    ;;# Optional output to netCDF               #
    ;;###########################################
    ;Optional output (controlled by diag_script_info)
    if(isatt(diag_script_info, "ncdf")) then
        data1@ncdf = diag_script_info@ncdf

        ;; Function in ~/interface_scripts/auxiliary.ncl
        ncdf_outfile = ncdf_write(data1, diag_script_info@ncdf)
    end if

    ;;###########################################
    ;;# Create the plot                         #
    ;;###########################################
    ;; There are two options to pass data to the plotting routine: directly
    ;; or via a netCDF file. All information needed to create a graphics file
    ;; are passed with the data. It is possible to provide "dummy_for_wks" and
    ;; "dummy_for_var" as arguments in the call.
    ;; You may invoke a wks here, and/or provide the variable name (see example
    ;; in next section). The routine "contour_map" is a function in
    ;; ~/diag_scripts/shared/plot/aux_plotting.ncl
    ;;
    ;; Option 1: pass data directly
    ;; -> may only work with NCL (i.e. ok here)
    ;; Dummies trigger defaults
     map = contour_map("dummy_for_wks", data1, "dummy_for_var")
    ;; Option 2: pass data via netCDF file
    ;;-> useful for non-NCL routines & metrics
;    map = contour_map("dummy_for_wks", ncdf_outfile,"dummy_for_var")

    ;;######################################################
    ;;# Separate output to a graphics file, Panelling      #
    ;;######################################################
    ;; Plotting routines by default create a graphics output file.
    ;; The default may be overridden by specifying res_gsnDraw=False (and
    ;; possibly res_gsnFrame=False).
    ;; Then file output must be done explicitely here, e.g. with the
    ;; "panelling" routine (function in ./diag_scripts/shared/plot/aux_plotting.ncl).
    ;; The separation between plotting and output might be useful, when
    ;; multiple calls to plotting routine(s) (e.g. in a loop) are needed to
    ;; create a (set of) plot(s).
    ;; This option may (currently) only work within NCL, because the plotting
    ;; routines return only a pointer within the current wks. Thus plots need
    ;; to be collected in an existing wks, like pages in a book. Variable
    ;; "pres" may be used to override panelling defaults.
    ;;   -> NOTE: The "panelling" routine does need a valid wks!!!!!
    ;; To test this option, please disable everything in the "Create the plot"
    ;; section and enable the following lines.
;    data1@res_gsnDraw           = False  ; Do not draw yet
;    data1@res_gsnFrame          = False  ; Don't advance frame.
    ;; function in aux_plotting.ncl
;    wks = get_wks("dummy_for_wks", diag_script, var0)
;    nmap = 3
;    maps = new(nmap,graphic) ; collect individual "map"s in a graphic array
;    do i = 0,nmap-1  ; this stupid loop creates 3 times the same plot
;        maps(i) = contour_map(wks, data1,var0)
;    end do
;    pres                  = True   ; needed to override panelling defaults
;    pres@gsnPanelLabelBar = False  ; no general label bar desired here
    ; number of vertical and horizontal plots per page: nmap, 1
;    outfile = panelling(wks, maps, nmap, 1, pres)
;    log_info(" Wrote " + outfile)

    leave_msg(diag_script, "")
end
