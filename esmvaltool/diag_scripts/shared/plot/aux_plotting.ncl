;;#############################################################################
;; AUXILIARY PLOTTING ROUTINES
;;#############################################################################
;; Please consider using or extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;
;;    procedure create_legend_lines
;;    function output_type
;;    procedure copy_VarAtt_sel
;;    function panelling
;;    function get_plot_dir
;;    function get_outfile_name
;;    function get_wks
;;    procedure add_markers
;;    procedure add_num_markers
;;    procedure add_errorbar
;;    procedure horizontal_whiskers
;;    procedure add_prediction_error
;;
;;#############################################################################


load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/logging.ncl"

;##############################################################################
undef("create_legend_lines")
procedure create_legend_lines(labels:string, \
                              styles, \
                              outfile:string,\
                              opt)
;;
;; Arguments
;;    labels:  labels of the lines
;;    styles: style of the lines in the plot
;;              @colors
;;              @dashes  -> optional (line dashes
;;              @thicks  -> optional (line/marker thickness)
;;              @markers -> optional (marker index)
;;              @sizes   -> optional (marker size)
;;    outfile: outfile directory
;;    opt: "lines" or "markers" legend
;;
;; Description
;;    Draws an extra plot with a legend, specified by labels and styles.
;;
;; Caveats
;;
;; Modification history:
;;    20150511-A_laue_ax: added safe lower limits for panelling plot
;;    20150508-A_righ_ma: added lines/markers option.
;;    20150120-A_gott_kl: remove pre-existing file type suffix
;;    20140305-A_righ_ma: modified to plot always as epsi format.
;;    20140219-A_fran_fr: written.
;;
local funcname, scriptname, region, temp, outfile, n_suff, wks_legend
begin

    funcname = "create_legend_lines"
    scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
    enter_msg(scriptname, funcname)

    ;; Flag if RGB or RGBA is used for colors
    if (dimsizes(dimsizes(styles@colors)).ne.1) then
        RGB = True
    else
        RGB = False
    end if

    ;; Open workstation
    wks_legend = gsn_open_wks("epsi", outfile)

    ;; General resources
    res = True
    res@tmBorderThicknessF = 0.
    res@tmBorderLineColor = "white"
    res@tmXTOn = False
    res@tmXBOn = False
    res@tmYLOn = False
    res@tmYROn = False
    res@trXMinF = 0
    res@trXMaxF = 1
    res@trYMinF = 0
    res@trYMaxF = 1
    res@vpXF = 0.0
    res@vpYF = 1.0

    ;; Resources for polyline and text
    resL = True
    resT = True
    resT@txFontHeightF = 0.013
    resT@txJust = "CenterLeft"

    ;; Set number of lines and labels
    dim_LAB = dimsizes(labels)

    ;; Rows and columns
    deltax = 0.06
    deltay = 0.03
    nrowmax = 32
    ncols = 3
    res@vpWidthF = 0.8
    if (dim_LAB .gt. nrowmax * ncols) then
        ncols = 4
        res@vpWidthF = 1.0
        if (dim_LAB .gt. nrowmax * ncols) then
            error_msg("f", "xy_line.ncl", funcname, "too many models")
        end if
    end if
    nrows = toint(dim_LAB/tofloat(ncols))

    line_xpos = new((/ncols, 2/), float)
    do jj = 0, ncols - 1
        line_xpos(jj, 0) = 0.01 + jj / tofloat(ncols)
    end do
    line_xpos(:, 1) = line_xpos(:, 0) + deltax
    text_pos = line_xpos(:, 1) + 0.01
    res@vpHeightF = 0.04 + deltay * (nrows - 1)

    ;; Draw blank plot
    blank = gsn_blank_plot(wks_legend, res)

    ;; Loop over labels, draw lines and text
    ypos = 0.98
    ii = 0
    do while (ii.lt.dim_LAB)
        ;; Set color
        if (RGB) then
            resL@gsLineColor = styles@colors(ii, :)
            resL@gsMarkerColor = styles@colors(ii, :)
        else
            resL@gsLineColor = styles@colors(ii)
            resL@gsMarkerColor = styles@colors(ii)
        end if
        ;; Set dash pattern
        if (isatt(styles, "dashes")) then
            resL@gsLineDashPattern    = styles@dashes(ii)
        end if
        ;; Set thickness
        if (isatt(styles, "thicks")) then
            resL@gsLineThicknessF    = styles@thicks(ii)
            resL@gsMarkerThicknessF  = styles@thicks(ii)
        end if
        ;; Set marker index
        if (isatt(styles, "markers")) then
            resL@gsMarkerIndex       = styles@markers(ii)
        end if
        ;; Set marker size
        if (isatt(styles, "sizes")) then
            resL@gsMarkerSizeF = styles@sizes(ii)
        end if
        ;; Draw legend
        if (opt.eq."lines") then
            gsn_polyline_ndc(wks_legend, line_xpos(mod(ii, ncols), :), \
                             (/ypos, ypos/), resL)
            gsn_text_ndc(wks_legend, labels(ii), text_pos(mod(ii, ncols)), \
                         ypos, resT)
        else if (opt.eq."markers") then
            gsn_polymarker_ndc(wks_legend, avg(line_xpos(mod(ii, ncols), :)), \
                             ypos, resL)
            gsn_text_ndc(wks_legend, labels(ii), \
                         text_pos(mod(ii, ncols)) - 0.02, ypos, resT)
        else
            error_msg("f", "xy_line.ncl", funcname, "unrecognized option " + \
                      opt)
        end if
        end if
        if (mod(ii, ncols).eq.ncols - 1) then
            ypos = ypos - deltay  ; go to the next row
        end if
        ii = ii + 1
    end do
    draw(blank)
    frame(wks_legend)

    leave_msg(scriptname, funcname)

end

;;#############################################################################
undef("output_type")
function output_type()
;;
;; Arguments
;;
;; Return value
;;    A string with the output file type
;;
;; Description
;;    Provides a default, if file type is not explicitly specified
;;
;; Caveats
;;
;; Modification history
;;    20131028-A_gott_kl: written.
;;
local funcname, scriptname, file_type
begin

    funcname = "output_type"
    scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
    enter_msg(scriptname, funcname)

    file_type = config_user_info@output_file_type
    if(ismissing(file_type)) then
        file_type = "ps"
    end if
    file_type = str_lower(file_type)

    leave_msg(scriptname, funcname)
    return(file_type)

end

;;#############################################################################
undef("copy_VarAtt_sel")
procedure copy_VarAtt_sel(var1, var2, sel: string)
;;
;; Arguments
;;    var1: variable of any type and dimension with source attributes 
;;    var2: variable of any type and dimension, shall receive selected 
;;          attributes from var1  
;;    sel: string (or list of strings) that specify the BEGINNING letters
;;         of the attributes to copy  
;;
;; Return value
;;    var2 gets additional attributes 
;;
;; Description
;;    Intended to copy selected plot ressources for use in a different 
;;    plot routine that may not allow all ressources from var1. 
;;    Written for function legend_lines (legends.ncl).
;;
;; Caveats:
;;
;; Modification history:
;;    20141227-A_gott_kl written.
;;
local funcname, scriptname, atts, i1, i2, n
begin

    funcname = "copy_VarAtt_sel"
    scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
    enter_msg(scriptname, funcname)

    ;; check if var1 has attributes at all 
    atts = getvaratts(var1)
    if(all(atts.eq."missing")) then
         log_debug("warning: no attributes to copy "+\
                     "(function copy_VarAtt_sel)")
         return(var2)  ; return var2 unchanged
    end if

    ;; copy selected attributes from var1
    ;; Existing attributes of the same name will be overwritten in var2.
    do i2 = 0, dimsizes(sel) - 1 
        n = strlen(sel(i2))
        do i1 = 0, dimsizes(atts) - 1
            if(str_get_cols(atts(i1), 0, n - 1) .eq. sel(i2)) then
                if(isatt(var2,atts(i1))) then
                    delete(var2@$atts(i1)$)
                end if
                var2@$atts(i1)$ = var1@$atts(i1)$  
            end if
        end do
    end do

    leave_msg(scriptname, funcname)

end

;;#############################################################################
undef("panelling")
function panelling(wks,
                   plots : graphic,
                   nvert[1] : integer,
                   nhori[1] : integer,
                   pres_in[1] : logical)
;;
;; Arguments
;;    wks: graphics object associated with a file.
;;    plots: graphics variable containing one or more plots.
;;    nvert: Maximum allowed number of plots per page (vertical).
;;    nhori: Maximum allowed number of plots per page (horizontal).
;;    pres_in: if it is a logical variable, attributes are used to override
;;             panelling defaults
;;
;; Return value
;;    A string with the output file name
;;
;; Description
;;    Writes a graphics file
;;
;; Caveats
;;    Might not be fit to be used with non-NCL routines.
;;
;; Modification history
;;    20131112-A_gott_kl: written.
;;
local funcname, scriptname, plots, info, wks, outfile
begin

    funcname = "panelling"
    scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
    enter_msg(scriptname, funcname)

    ;; check if valid wks has been provided
    if (typeof(wks).ne."graphic") then
        error_msg("f", "aux_plotting.ncl", funcname, "panelling of " + \
                  "multiple plots needs valid wks")
    end if

    ;; get outfile name
    if (isatt(wks, "fullname")) then
        outfile = wks@fullname
    else
        outfile = wks@name
        error_msg("w", "aux_plotting.ncl", funcname, "wks@fullname " + \
                  "missing, consider to use get_wks to open wks")
    end if

    ;; panel parameters
    npages = toint(ceil(tofloat(dimsizes(plots)) / tofloat(nvert * nhori)))

    ;; Draw the panel
    pres = True
    pres@gsnMaximize       = True
    ;pres@gsnPanelLabelBar  = True
    pres@gsnPanelYWhiteSpacePercent = 5
    pres@gsnPanelXWhiteSpacePercent = 5
    if (isatt(diag_script_info, "PanelTop")) then
        top = tofloat(diag_script_info@PanelTop)
    else
        top = 0.99  ; default
    end if
    pres@gsnPanelTop = top
    if (npages .gt. 1) then
        pres@gsnPaperOrientation   = "portrait"
        if(typeof(pres_in).eq."logical") then  ; Override defaults
            copy_VarMeta(pres_in, pres)
        end if
        ia = 0
        ie = nvert * nhori - 1
        ie = min((/dimsizes(plots)-1, ie/))  ; a_laue_ax: added for cloud diag
        do ipage = 0, npages - 1
            gsn_panel(wks, plots(ia:ie), (/nvert, nhori/), pres)
            ia = ie + 1
            ie = ia + (nvert * nhori) - 1
            ie = min((/dimsizes(plots)-1, ie/))  ; a_laue_ax
        end do
    else
        pres@gsnPaperOrientation   = "auto"
        if(typeof(pres_in).eq."logical") then  ; Override defaults
            copy_VarMeta(pres_in, pres)
        end if
        nrow = toint(ceil(tofloat(dimsizes(plots)) / tofloat(nhori)))
        ncol = min((/nhori, dimsizes(plots)/))
        gsn_panel(wks, plots, (/nrow, ncol/), pres)
    end if

    leave_msg(scriptname, funcname)
    return(outfile)

end

;;#############################################################################
undef("get_plot_dir")
function get_plot_dir()
;;
;; Arguments
;;
;; Return value
;;    A string root path for plots
;;
;; Description
;;    Provides a default, if plot_dir is not explicitely specified
;;
;; Caveats
;;
;; Modification history
;;    20131104-A_gott_kl: written.
;;
local funcname, scriptname, plot_dir
begin

    funcname = "get_plot_dir"
    scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
    enter_msg(scriptname, funcname)

    plot_dir = config_user_info@plot_dir
    if(ismissing(plot_dir)) then
        plot_dir = "./plots/"
    end if
    if (.not. fileexists(plot_dir)) then
        system("mkdir -p " + plot_dir)
    end if

    leave_msg(scriptname, funcname)
    return(plot_dir)

end

;;#############################################################################
undef("get_outfile_name")
function get_outfile_name(add_specs[1] : string)
;;
;; Arguments
;;    add_specs: string containing specific elements to be added to the file
;;               name if not of type string, only variable name is used.
;;
;; Return value
;;    Complete outfile name incl. path, additional specifications and file type
;;    plot_dir + add_specs + file_type
;;
;; Description
;;    Fetches file_type, plot_dir
;;    Fetches string with additional elements specified within diag script
;;
;; Modification history
;;    20131204-A_senf_da: generalized naming.
;;    20131104-A_gott_kl: written.
;;
local funcname, scriptname, outfile, file_type, plot_dir, diag_script_base, \
    add_specs, output_dir
begin

    funcname = "get_outfile_name"
    scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
    enter_msg(scriptname, funcname)

    ;; Collect parts of outfile name
    file_type = output_type()
    plot_dir = get_plot_dir()

    ;; Assemble final name
    if (typeof(add_specs).ne."string") then
        outfile = plot_dir + var + "." + file_type
    else
        outfile = plot_dir + add_specs + "." + file_type
    end if

    leave_msg(scriptname, funcname)
    return(outfile)
end

;;#############################################################################
undef("get_wks")
function get_wks(wks_in,
                 diag_script[1]: string,
                 add_specs[1]: string)
;;
;; Arguments
;;    wks_in: dummy or graphic object
;;    diag_script: name of diagnostic script
;;    add_specs: own specificactions to be added to file name, i.e. variable
;;               name, etc. - needed for function get_outfile_name
;;
;; Return value
;;    wks: graphic object
;;
;; Description
;;    Provides a default wks, if wks_in is not of type "graphic".
;;    Attribute wks@fullname is used to transfer the output file name, since
;;    wks@name cuts off the path to the file name.
;;
;; Caveats
;;
;; Modification history
;;    20131113-A_gott_kl: written.
;;
local funcname, scriptname, wks_in, diag_script, add_specs, file_type
begin

    funcname = "get_wks"
    scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
    enter_msg(scriptname, funcname)

    ;; Check if a valid wks has been provided, otherwise invoke default
    if(typeof(wks_in).ne."graphic") then  ; Use defaults
        ;; Function in ~/diag_scripts/shared/plot/aux_plotting.ncl
        file_type = output_type()
        copy_VarMeta(wks_in,file_type) 

        ;; Function in ~/diag_scripts/shared/plot/aux_plotting.ncl
        outfile = get_outfile_name(add_specs)
        wks = gsn_open_wks(file_type, outfile)
        wks@fullname = outfile
    else
        wks = wks_in
        copy_VarMeta(wks_in, wks)
    end if

    leave_msg(scriptname, funcname)
    return(wks)

end

;;#############################################################################
undef("add_markers")
procedure add_markers(wks[1] : graphic,
                      plot[1] : graphic,
                      res_in[1] : logical,
                      xpos_in : numeric,
                      ypos_in : numeric)
;;
;; Arguments
;;    wks: valid workstation, e.g. created by get_wks
;;    plot: plot identifier, e.g. created by gsn_*
;;    res_in: plot ressources that may override local function defaults
;;    xpos_in: horizontal marker position(s)  
;;    ypos_in: vertical marker position(s) 
;;
;; Return value
;;    Attaches polyline IDs as attributes to plot.
;;
;; Description:
;;    Adds markers to an existing plot.
;;
;; Caveats:
;;
;; Modification history:
;;    20140224-A_gott_kl: written for use with profile_plev.ncl in Emmons.ncl
;;
local funcname, scriptname, res, str, atts, iatt, color
begin

    funcname = "add_markers"
    scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
    enter_msg(scriptname, funcname)

    ;; set default plot ressources 
    res               = True
    res@gsMarkerSizeF = 0.006
    res@gsMarkerColor = "black" 
    res@gsMarkerIndex = 16  ; solid circle 
    
    ;; override defaults (after cleaning up for possible vector input)
    if(res_in) then
        atts = getvaratts(res_in)
        if(.not.all(ismissing(atts))) then
            do iatt = 0, dimsizes(atts)-1
                if(isatt(res, atts(iatt))) then
                    delete(res@$atts(iatt)$) 
                end if
            end do
        end if
        copy_VarMeta(res_in, res)
    end if

    ;; plot symbol rim (slightly bigger, white polymarker background)
    color = res@gsMarkerColor  ; remember original color
    res@gsMarkerSizeF = res@gsMarkerSizeF + res@gsMarkerSizeF * 0.2
    res@gsMarkerColor(:) = "white"
    str = unique_string("polymarker")
    plot@$str$ = gsn_add_polymarker(wks, plot, xpos_in, ypos_in, res)

    ;; plot symbols 
    res@gsMarkerColor(:) = color  ; restore original color
    res@gsMarkerSizeF = res@gsMarkerSizeF - res@gsMarkerSizeF * 0.2
    str = unique_string("polymarker")
    plot@$str$ = gsn_add_polymarker(wks, plot, xpos_in, ypos_in, res)

    leave_msg(scriptname, funcname)

end

;;#############################################################################
undef("add_num_markers")
procedure add_num_markers(wks[1] : graphic,
                          plot[1] : graphic,
                          res_in[1] : logical,
                          xpos_in : numeric,
                          ypos_in : numeric)
;;
;; Arguments
;;    wks: valid workstation, e.g. created by get_wks
;;    plot: plot identifier, e.g. created by gsn_*
;;    res_in: plot ressources that may override local function defaults
;;    xpos_in: horizontal marker position(s), needs year as a coordinate  
;;    ypos_in: vertical marker position(s) 
;;
;; Return value
;;    Attaches polytext IDs as attributes to plot.
;;
;; Description:
;;    Adds markers to an existing plot, being the last two intagers of a year.
;;    If a horizontal (vertical) coordinate
;;    has only one element, then this position is used for all markers.
;;
;; Caveats:
;;
;; Modification history:
;;    20150914-A_wenz_sa: written. 
;;
local funcname, scriptname, xpos, ypos, res, str, atts, iatt, color
begin

    funcname = "add_num_markers"
    scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
    enter_msg(scriptname, funcname)

    ;; set default plot ressources 
    res               = True
    res@txFontHeightF = 0.025
    res@txFontColor   = "black" 

    ;; override defaults (after cleaning up for possible vector input)
    if(res_in) then
        atts = getvaratts(res_in)
        if(.not.all(ismissing(atts))) then
            do iatt = 0, dimsizes(atts)-1
                if(isatt(res, atts(iatt))) then
                    delete(res@$atts(iatt)$) 
                end if
            end do
        end if
        copy_VarMeta(res_in, res)
    end if

    atext = new(dimsizes(xpos_in),integer)

       do y = 0, dimsizes(xpos_in)-1
           strs     = str_split_by_length(tostring(xpos_in&year(y)), 1)
           strs_y   = toint(strs(0) + strs(1) + "00")
           atext(y) = xpos_in&year(y) - strs_y
           delete([/strs,strs_y/])
      end do

    ;; create vector, if horizontal or vertical position is given as scalar
    if(dimsizes(xpos_in).eq.1) then
        xpos = new(dimsizes(ypos_in),typeof(xpos_in))
        xpos = (/xpos_in/)
    else
        xpos = (/xpos_in/)  
    end if 
    if(dimsizes(ypos_in).eq.1) then
        ypos = new(dimsizes(xpos_in),typeof(ypos_in))
        ypos = (/ypos_in/)
    else
        ypos = (/ypos_in/)
    end if

    ;; plot symbol 
    str = unique_string("polymarker")
    plot@$str$ = gsn_add_text(wks,plot,sprinti("%-0.2hi",atext),\
                                                 xpos,ypos,res)

    leave_msg(scriptname, funcname)

end

;;#############################################################################
undef("add_errorbar")
procedure add_errorbar(wks[1] : graphic,
                       plot[1] : graphic,
                       res_in[1] : logical,
                       xpos_in : numeric,
                       ypos_in : numeric)
;;
;; Arguments
;;    wks: valid workstation, e.g. created by get_wks
;;    plot: plot identifier, e.g. created by gsn_*
;;    res_in: plot ressources that may override local function defaults
;;    xpos_in: 2D array with horizontal marker position(s) and stddev's, 
;;    ypos_in: 2D array with vertical marker position(s) and stddev's
;;
;; Return value
;;    Attaches polytext IDs as attributes to plot.
;;
;; Description:
;;    Adds error bars to an existing plot. If a horizontal (vertical) coordinate
;;    has only one element, then this position is used for all markers.
;;
;; Caveats:
;;
;; Modification history:
;;    20150914-A_wenz_sa: written. 
;;
local funcname, scriptname, xpos, ypos, res, str, atts, iatt, color
begin

    funcname = "add_errorbar"
    scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
    enter_msg(scriptname, funcname)

    ;; set default plot ressources 
    res                  = True
    res@gsLineThicknessF = 2

    colors = new(dimsizes(xpos_in(0,:)),string)
    colors = "black"
     
    ;; override defaults (after cleaning up for possible vector input)
    if(res_in) then
        atts = getvaratts(res_in)
        if(.not.all(ismissing(atts))) then
            do iatt = 0, dimsizes(atts)-1
                if(isatt(res, atts(iatt))) then
                    delete(res@$atts(iatt)$) 
                end if
            end do
        end if
        copy_VarMeta(res_in, res)

        if isatt(res,"gsLineColor") then
            colors = res@gsLineColor
            delete(res@gsLineColor)
        end if
    end if

    ;; create vector, if horizontal or vertical position is given as scalar
    if(dimsizes(xpos_in(0,:)).eq.1) then
        xpos = new(dimsizes(ypos_in),typeof(xpos_in))
        xpos(:,:) = xpos_in
    else
        xpos = xpos_in   
    end if 

    if(dimsizes(ypos_in(0,:)).eq.1) then
        ypos = new(dimsizes(xpos_in),typeof(ypos_in))
        ypos(:,:) = ypos_in
    else
        ypos = ypos_in
    end if

    Xend = 0.01*((max(ypos(0,:))+0.2*(abs(max(ypos(0,:)))))\
                -(min(ypos(0,:))-0.2*(abs(min(ypos(0,:))))))
    Yend = 0.01*((max(xpos(0,:))+0.2*(abs(avg(xpos(0,:)))))\
                -(min(xpos(0,:))-0.2*(abs(avg(xpos(0,:))))))

    do i =0,dimsizes(xpos(0,:))-1
        res@gsLineColor = colors(i)

        if (.not.any(ismissing((/xpos(0,i),ypos(0,i)/))))
            Xbarstr = unique_string("Xbarerror")
            Ybarstr = unique_string("Ybarerror")
            Xendstr = unique_string("Xenderror")
            Yendstr = unique_string("Yenderror")
            plot@$Xbarstr$ = gsn_add_polyline(wks,plot,\
                                         (/xpos(0,i)+xpos(1,i),\
                                           xpos(0,i)-xpos(1,i)/)\
                                        ,(/ypos(0,i),ypos(0,i)/),res)

            plot@$Ybarstr$ = gsn_add_polyline(wks,plot,\
                                         (/xpos(0,i),xpos(0,i)/)\
                                        ,(/ypos(0,i)+ypos(1,i),\
                                           ypos(0,i)-ypos(1,i)/),res)    

            plot@$Xendstr$ = gsn_add_polyline(wks,plot,\
                                         (/xpos(0,i)+xpos(1,i),\
                                           xpos(0,i)+xpos(1,i),\
                                           xpos(0,i)+xpos(1,i),\
                                           xpos(0,i)-xpos(1,i),\
                                           xpos(0,i)-xpos(1,i),\
                                           xpos(0,i)-xpos(1,i)/)\
                                        ,(/ypos(0,i)-Xend,ypos(0,i)+Xend,\
                                           ypos(0,i),ypos(0,i),\
                                           ypos(0,i)-Xend,ypos(0,i)+Xend/),res)

            plot@$Yendstr$ = gsn_add_polyline(wks,plot,\
                                         (/xpos(0,i)-Yend,xpos(0,i)+Yend,\
                                           xpos(0,i),xpos(0,i),\
                                           xpos(0,i)-Yend,xpos(0,i)+Yend/)\
                                         ,(/ypos(0,i)+ypos(1,i),\
                                            ypos(0,i)+ypos(1,i),\
                                            ypos(0,i)+ypos(1,i),\
                                            ypos(0,i)-ypos(1,i),\
                                            ypos(0,i)-ypos(1,i),\
                                            ypos(0,i)-ypos(1,i)/),res)
        end if
    end do

    leave_msg(scriptname, funcname)

end


;;#############################################################################
undef("horizontal_whiskers")
procedure horizontal_whiskers(wks[1] : graphic,
                              plot[1] : graphic,
                              res_in[1] : logical,
                              xmin_in : numeric,
                              xmax_in : numeric,
                              ypos_in: numeric)
;;
;; Arguments
;;    wks: valid workstation, e.g. created by get_wks.
;;    plot: plot identifier, e.g. created by gsn_*.
;;    res_in: plot ressources that may override local function defaults.
;;    xmin_in: vector of whiskers' left ends (same size as xmax & y).
;;    xmax_in: vector of whiskers' right ends (same size as xmin & y).
;;    ypos_in: vector of whiskers' vertical positions (must have same size as
;;             xmax & xmin).
;;
;; Return value
;;    Attaches polyline IDs as attributes to plot.
;;
;; Description
;;    Creates vectors suitable as input for gsn_add_polyline:
;;      x = (/xmin1,xmax1,_FillValue,xmin2,xmax2,_FillValue, .../)
;;      y = (/ypos1,ypos1,_FillValue,ypos2,ypos2,_FillValue, .../)
;;    The separation by _FillValue results in individual whiskers.
;;    No whisker is created where xmin, xmax or ypos is missing.
;;
;; Caveats
;;
;; References
;;    www.ncl.ucar.edu/Document/Graphics/Interfaces/gsn_add_polyline.shtml
;;
;; Modification history
;;    20140224-A_gott_kl: written.
;;
local funcname, scriptname, n, i, j, xmin, xmax, ypos, imiss, ivalid, result, \
    x, y, res, str
begin

    funcname = "horizontal_whiskers"
    scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
    enter_msg(scriptname, funcname)

    ;; check sizes of input vectors
    n = dimsizes(ypos_in)
    if(dimsizes(xmin_in).ne.n .or. dimsizes(xmax_in).ne.n) then
        error_msg("f", "aux_plotting.ncl", funcname, "all input vectors " + \
                  "must have same size")
    end if

    ;; unify _FillValue and type
    ypos = tofloat(ypos_in)  ; #####################
    copy_VarMeta(ypos_in, ypos)
    ypos@_FillValue = default_fillvalue("float")
    ypos = where(ismissing(ypos_in), ypos@_FillValue, ypos)
    xmin = tofloat(xmin_in)  ; #####################
    copy_VarMeta(xmin_in, xmin)
    xmin@_FillValue = default_fillvalue("float")
    xmin = where(ismissing(xmin_in), xmin@_FillValue, xmin)
    xmax = tofloat(xmax_in)  ; #####################
    copy_VarMeta(xmax_in, xmax)
    xmax@_FillValue = default_fillvalue("float")
    xmax = where(ismissing(xmax_in), xmax@_FillValue, xmax)

    ;; find positions where no input element is missing
    imiss = ispan(0, n - 1, 1)
    imiss@_FillValue = default_fillvalue("integer")
    imiss = where(ismissing(ypos), imiss@_FillValue, imiss)
    imiss = where(ismissing(xmin), imiss@_FillValue, imiss)
    imiss = where(ismissing(xmax), imiss@_FillValue, imiss)
    ivalid = ind(.not. ismissing(imiss))

    if(all(ismissing(ivalid))) then
        log_debug(" No valid data to plot whiskers: Skipping")
    else
        ;; create vectors suited for gsn_add_polyline
        x = new(3 * dimsizes(ivalid) - 1, typeof(xmin))
        x@_FillValue = xmin@_FillValue
        x(:) = x@_FillValue
        y = new(3 * dimsizes(ivalid) - 1, typeof(ypos))
        y@_FillValue = xmax@_FillValue
        y(:) = y@_FillValue
        do i = 0, dimsizes(ivalid) - 1
            j = ivalid(i)
            x(i * 3)   = xmin(j)
            x(i * 3 + 1) = xmax(j)
            y(i * 3)   = ypos(j)
            y(i * 3 + 1) = ypos(j)
        end do

        ;; resources
        res = True
        if(res_in) then
            copy_VarMeta(res_in, res)
        end if

        ;; plot whiskers
        ;; (see Example 2 of gsn_add_polyline help page):
        ;; Attaching the return value of gsn_add_polyline
        ;; to the plot variable ensures that this
        ;; reference can live outside this script
        str = unique_string("polyline")
        plot@$str$ = gsn_add_polyline(wks, plot, x, y, res)
    end if

    leave_msg(scriptname, funcname)

end

;;#############################################################################
undef("add_prediction_error")
procedure add_prediction_error(wks[1] : graphic,
                          plot[1] : graphic,
                          res_in[1] : logical,
                          xpos_in : numeric,
                          ypos_in : numeric)
;;
;; Arguments
;;    wks: valid workstation, e.g. created by get_wks
;;    plot: plot identifier, e.g. created by gsn_*
;;    res_in: plot ressources that may override local function defaults
;;    xpos_in: array defining the x-position of the four corners 
;;    ypos_in: vertical marker position(s) 
;;
;; Return value
;;    Attaches polytext IDs as attributes to plot.
;;
;; Description:
;;    Adds markers to an existing plot, being the last two intagers of a year.
;;    If a horizontal (vertical) coordinate
;;    has only one element, then this position is used for all markers.
;;
;; Caveats:
;;
;; Modification history:
;;    20150914-A_wenz_sa: written. 
;;
local funcname, scriptname, xpos, ypos, res, str, atts, iatt, color
begin

    funcname = "add_prediction_error"
    scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
    enter_msg(scriptname, funcname)

    ;; set default plot ressources 
    res               = True
    ;;Create a plot with filled area between both curves
    res@xyLineColor       = "grey85"
    res@gsFillColor       = "grey85"
    res@gsLineDashPattern = 0

    ;; override defaults (after cleaning up for possible vector input)
    if(res_in) then
        atts = getvaratts(res_in)
        if(.not.all(ismissing(atts))) then
            do iatt = 0, dimsizes(atts)-1
                if(isatt(res, atts(iatt))) then
                    delete(res@$atts(iatt)$) 
                end if
            end do
        end if
        copy_VarMeta(res_in, res)
    end if

    ;; plot symbol 
    str = unique_string("polygon")
    plot@$str$ = gsn_add_polygon(wks,plot,xpos_in,ypos_in,res)

    leave_msg(scriptname, funcname)

end

