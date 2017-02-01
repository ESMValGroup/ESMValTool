;;#############################################################################
;; SouthernOcean_transport.ncl
;; Author: Bert van Ulft (KNMI, the Netherlands)
;; EMBRACE project
;;#############################################################################
;;
;; Description
;;    Calculates the volume transport across a section
;;    Calculating and derivation of broken line is based on cdftransport.F90
;;        from the CDFTOOLS package
;;    CDFTOOLS: http://servforge.legi.grenoble-inp.fr/projects/CDFTOOLS
;;
;; Required diag_script_info attributes (diagnostic specific)
;;    styleset    : plot style
;;    lon1        : longitude of start of section
;;    lat1        : latitude of start of section
;;    lon2        : longitude of end of section
;;    lat2        : latitude of end of section
;;    section_name: name used in title
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;
;; Caveats
;;    Only for EC-Earth/NEMO files as only for these files we have all
;;        required grid info on the staggered grids
;;    The grid points nearest to the lat and lon in the cfg file are used.
;;        .i.e. it may deviate a bit if the specified coordinates are not
;;        exactly on a grid point
;;    It is assumed that the u and v are relative to the geographic grid, i.e.
;;        no rotation is performed here. As this is true for EC-Earth on the
;;        Southern hemisphere, abort for positive lats.
;;
;; Modification history
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20141028-A_vanu_be: written.
;;
;;#############################################################################

;; Load additional ncl code
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

;; Load Temporary file with current namelist configuration
load "./interface_data/ncl.interface"

;; Load auxiliary NCL routines
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/style.ncl"

load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./plot_scripts/ncl/xy_line.ncl"

;;#############################################################################
;; Collect information
;;#############################################################################
begin
    nl = str_get_nl()
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(diag_script, "", 4)
    diagsc_base = systemfunc("basename " + diag_script)

    ;; print info, part 1
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script, verbosity, 1)

    ;; check that we have two variables
    dim_VAR = dimsizes(variables)
    if (dim_VAR .ne. 2) then
        error_msg("f", diag_script, "", "must have 2 input " + \
                  "variables (u and v), but got " + dim_VAR + ".")
    end if
    varu = variables(0)     ; u variable
    varv = variables(1)     ; v variable

    write_references(diag_script,           \  ; script name
                     "A_vanu_be",           \  ; authors
                     (/""/),                \  ; contributors
                     "D_0000",              \  ; diag_references
                     "",                    \  ; obs_references
                     (/"P_embrace"/))          ; proj_references

    ;; print info, part 2
    info_output("    transport calculated from: " + varu + " & " + varv, \
                verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ;; Check field types
    field_type0 = field_types(0)
    if (field_type0 .ne. "T3M") then
        error_msg("f", diag_script, "", "can't use field of type " + \
                  field_type0 + " for diag_script " + diag_script + ".")
    end if
    if (.not. all(field_types(1:) .eq. field_type0)) then
        error_msg("f", diag_script, "", "field types of all " + \
                  "variables must be the same for diag_script " + \
                  diag_script + " but they are not: " + nl + \
                  str_join(str_insert(variables + ": " + field_types, \
                           "field_type for var ", 0), nl))
    end if

    ;; Get environment variables
    file_type = getenv("ESMValTool_output_file_type")
    if(ismissing(file_type)) then
        file_type = "PS"
    end if

    ;; Load var_info
    vardeffile = "interface_data/" + varu + "_info.tmp"
    loadscript(vardeffile)
    vardeffile = "interface_data/" + varv + "_info.tmp"
    loadscript(vardeffile)

    ;; output file names
    outfile_xy = "transport_" + varu + "_" + varv

    ;; check all required diag_script_info attributes have been set
    required_attributes = (/"styleset",\
                            "lon1",\
                            "lat1",\
                            "lon2",\
                            "lat2", \
                            "section_name"/)
    exit_if_missing_atts(diag_script_info, required_attributes)

    ;; Check we are on the Southern hemisphere
    if (diag_script_info@lat1 .gt. 0 .or. diag_script_info@lat2 .gt. 0) then
        error_msg("f", diag_script, "", "You should only calculate " + \
                  "transport on the Southern hemisphere.")
    end if
end

;;#############################################################################
;; Prepare for analysis and plotting
;;#############################################################################
begin
    ;; Prepare annotation
    annots = project_style(diag_script_info, "annots")

    ;; open workstation for plotting
    wks_xy = get_wks(0, diag_script, outfile_xy)

    diag_script_info@multi_model_mean = "n"  ; required for xy_line scr.
end

;;#############################################################################
;; Read data
;;#############################################################################
begin
    dim_MOD = dimsizes(models@name)
    do idx_mod = 0, dim_MOD - 1
        datau = read_data(idx_mod, varu, field_type0)
        datav = read_data(idx_mod, varv, field_type0)
    end do
end

;;#############################################################################
;; Exact calculation as in CDFTOOLS, only for ECEARTH project files
;; because dimensions of staggered grid are needed
;;#############################################################################
begin
    ;; declare storage array
    MAXTSTEPS  = 150 * 12
    dvoltrpsum = new((/dim_MOD, MAXTSTEPS/), double)
    dvoltrpsum(:, :) = 0.0
    times = new((/dim_MOD, MAXTSTEPS/), double)
    ntimes = new(dim_MOD, integer)
    ntimes(:) = 0

    do idx_mod = 0, dim_MOD - 1
        if (models@project(idx_mod) .eq. "ECEARTH") then
            ;; collect the required data, including all grid sizes
            extru = extract_area(idx_mod, datau, varu,  varu)
            latu  = extract_area(idx_mod, datau, "lat", varu)
            lonu  = extract_area(idx_mod, datau, "lon", varu)
            e2u   = extract_area(idx_mod, datau, "e2u", varu)
            e3u   = extract_area(idx_mod, datau, "e3u", varu)
            extrv = extract_area(idx_mod, datav, varv,  varv)
            latv  = extract_area(idx_mod, datav, "lat", varv)
            lonv  = extract_area(idx_mod, datav, "lon", varv)
            e1v   = extract_area(idx_mod, datav, "e1v", varv)
            e3v   = extract_area(idx_mod, datav, "e3v", varv)

            ;; check time and lev axes of u and v are the same
            if (any(extru&time .ne. extrv&time) .or. \
                extru&time@units .ne. extrv&time@units) then
                error_msg("f", diag_script, "", "time axes in u " + \
                          "and v files differ for model " + \
                          annots(idx_mod) + ".")
            end if
            if (any(extru&lev .ne. extrv&lev) .or. \
                extru&lev@units .ne. extrv&lev@units) then
                error_msg("f", diag_script, "", "lev axes in u " + \
                          "and v files differ for model " + \
                          annots(idx_mod) + ".")
            end if

            ;; store number of time steps and their values
            ntimes(idx_mod) = dimsizes(extru&time)
            if (ntimes(idx_mod) .gt. MAXTSTEPS) then
                error_msg("f", diag_script, "", "model " + \
                          annots(idx_mod) + " has more time steps " + \
                          "than can be handled. Increase MAXTSTEPS " + \
                          "in " + diag_script + ".")
            end if
            times(idx_mod, :ntimes(idx_mod) - 1) = extru&time
            if (.not. isatt(times, "units")) then
                times@units = extru&time@units
            else if (times@units .ne. extru&time@units) then
                error_msg("f", diag_script, "", "the input files " + \
                          "must have the same time units for all " + \
                          "models.")
            end if
            end if

            ;; find indices of points closest to supplied coordinates
            ;; take x (or i) from u and y (or j) from v
            ix_x = getind_latlon2d(latv, lonu, \
                                   (/diag_script_info@lat1, diag_script_info@lat2/),\
                                   (/diag_script_info@lon1, diag_script_info@lon2/))

            ;; code below taken based on cdftransport (CDFTOOLS)
            ;; Find broken line between P1 (ii0,ij0) and P2 (ii1, ij1)
            ii0   = ix_x(0, 1)
            ij0   = ix_x(0, 0)
            ii1   = ix_x(1, 1)
            ij1   = ix_x(1, 0)
            jpseg = 1000  ; max array size
            yypt  = new((/jpseg, 2/), integer)
            rxx   = new(jpseg, integer)
            ryy   = new(jpseg, integer)

            ;; compute direction of integrations and signs
            ;; The transport across the section is the dot product of
            ;; integral(line){(Mx,My)*dS}
            ;; Mx=integral(u*dz)  My=integral(v*dz)) and dS=(dy,-dx)}

            ;; By defining the direction of the integration as
            ;; positive to the east  or if ii1=ii0
            ;; positive to the north or if ij1=ij0
            if (ii1 .lt. ii0) then
                idirx = -1
            else
                idirx = 1
            end if
            if (ij1 .lt. ij0) then
                idiry = -1
            else
                idiry = 1
            end if

            ;;Then dS=(e2u*idiry,-e1v*idirx)
            ;;This will produce the following sign convention:
            ;;    West-to-east  (dx>0, dy=0) =>
            ;;                  -My*dx (-ve for a northward flow)
            ;;    South-to-north(dy>0, dx=0) =>
            ;;                   Mx*dy (+ve for an eastward flow)
            norm_u =  idiry
            norm_v = -idirx

            ;; Compute equation:  ryj = aj rxi + bj [valid in the
            ;; (i,j) plane]
            if (ii1 - ii0 .ne. 0) then
                aj = (ij1 - ij0) / tofloat(ii1 - ii0)
                bj = ij0 - aj * ii0
            else
                aj = 10000.0  ; flag value
                bj = 0.0
            end if

            ;; Compute equation:  rxi = ai ryj + bi [valid in the
            ;; (i,j) plane]
            if ((ij1 - ij0) .ne. 0) then
                ai = (ii1 - ii0) / tofloat(ij1 - ij0)
                bi = ii0 - ai * ij0
            else
                ai = 10000.0  ; flag value
                bi = 0.0
            end if

            ;; Compute the integer pathway: a succession of F points
            np = 0
            ;; Choose the straight line with the smallest slope
            if (abs(aj) .le. 1) then
                ;; Here, the best line is y(x)
                ;; If ii1 < ii0 swap points [ always describe section
                ;; from left to right ]
                if (ii1 .lt. ii0) then
                    iitmp = ii0
                    ijtmp = ij0
                    ii0   = ii1
                    ij0   = ij1
                    ii1   = iitmp
                    ij1   = ijtmp
                end if

                ;; iist,ijst is the grid offset to pass from F point to
                ;; either U/V point
                if (ij1 .ge. ij0) then     ; line heading NE
                    iist = 1
                    ijst = 1
                else                       ; line heading SE
                    iist = 1
                    ijst = 0
                end if

                ;; compute the nearest ji point on the line crossing at
                ;; ji
                do ji = ii0, ii1
                    if (np .ge. jpseg) then
                        error_msg("f", diag_script, "", "too many " + \
                                  "sections. Increase value of " + \
                                  "jpseg in code or reduce section " + \
                                  "length.")
                    end if
                    ij = toint(decimalPlaces(aj * ji + bj, 0, True))
                    yypt(np, 0) = ji
                    yypt(np, 1) = ij
                    np = np + 1
                end do
            else
                ;; Here, the best line is x(y)
                ;; If ij1 < ij0 swap points [ always describe section
                ;; from bottom to top ]
                if (ij1 .lt. ij0) then
                    iitmp = ii0
                    ijtmp = ij0
                    ii0   = ii1
                    ij0   = ij1
                    ii1   = iitmp
                    ij1   = ijtmp
                end if

                ;; iist,ijst is the grid offset to pass from F point to
                ;; either U/V point
                if (ii1 .ge. ii0) then
                    iist = 1
                    ijst = 1
                else
                    iist = 0
                    ijst = 1
                end if

                ;; compute the nearest ji point on the line crossing at
                ;; jj
                do jj = ij0, ij1
                    if (np .ge. jpseg) then
                        error_msg("f", diag_script, "", "too many " + \
                                  "sections. Increase value of " + \
                                  "jpseg in code or reduce section " + \
                                  "length.")
                    end if
                    ii = toint(decimalPlaces(ai * jj + bi, 0, True))
                    yypt(np, 0) = ii
                    yypt(np, 1) = jj
                    np = np + 1
                end do
            end if

            ;; Look for intermediate points to be added.
            ;; The final positions are saved in rxx,ryy
            rxx(0) = yypt(0, 0)
            ryy(0) = yypt(0, 1)
            nn     = 0

            do jk = 1, np - 1
                ;;.. distance between 2 neighbour points
                rd = sqrt((yypt(jk, 0) - yypt(jk - 1, 0)) ^ 2 + \
                          (yypt(jk, 1) - yypt(jk - 1, 1)) ^ 2)
                ;; intermediate points required if rd > 1
                if (rd .gt. 1) then
                    nn = nn + 1
                    if (nn .ge. jpseg) then
                        error_msg("f", diag_script, "", "too many " + \
                                  "sections. Increase value of " + \
                                  "jpseg in code or reduce section " + \
                                  "length.")
                    end if

                    if (abs(aj) .le. 1) then
                        ;; use y(x)
                        ;; possible intermediate points:
                        ylptmp1 = (/yypt(jk - 1, 0) + 1, \      ; M1 i
                                    yypt(jk - 1, 1) + 0/)       ; M1 j
                        if (aj .lt. 0) then
                            ylptmp2 = (/yypt(jk - 1, 0) + 0, \  ; M2 i
                                        yypt(jk - 1, 0) - 1/)   ; M2 j
                        else
                            ylptmp2 = (/yypt(jk - 1, 0) + 0, \  ; M2 i
                                        yypt(jk - 1, 0) + 1/)   ; M2 j
                        end if

                        ;;  M1 is the candidate point:
                        zxm = ylptmp1(0)
                        zym = ylptmp1(1)
                        za0 = aj
                        zb0 = bj
                        za1 = -1.0 / za0
                        zb1 = zym - za1 * zxm
                        ;; P1 is the projection of M1 on the strait line
                        zxp = -(zb1 - zb0) / (za1 - za0)
                        zyp = za0 * zxp + zb0
                        ;; zd1 is the distance M1P1
                        zd1 = (zxm - zxp) * (zxm - zxp) + \
                              (zym - zyp) * (zym - zyp)

                        ;; M2 is the candidate point:
                        zxm = ylptmp2(0)
                        zym = ylptmp2(1)
                        za1 = -1. / za0
                        zb1 = zym - za1 * zxm
                        ;; P2 is the projection of M2 on the strait line
                        zxp = -(zb1 - zb0) / (za1 - za0)
                        zyp = za0 * zxp + zb0
                        ;; zd2 is the distance M2P2
                        zd2 = (zxm - zxp) * (zxm - zxp) + \
                              (zym - zyp) * (zym - zyp)
                    else
                        ;; use x(y)
                        ;; possible intermediate points:
                        if (ai .lt. 0) then
                            ylptmp1 = (/yypt(jk - 1, 0) - 1, \  ; M1 i
                                        yypt(jk - 1, 0) + 0/)   ; M1 j
                        else
                            ylptmp1 = (/yypt(jk - 1, 0) + 1, \  ; M1 i
                                        yypt(jk - 1, 0) + 0/)   ; M1 j
                        end if
                        ylptmp2 = (/yypt(jk - 1, 0) + 0, \      ; M2 i
                                    yypt(jk - 1, 1) + 1/)       ; M2 j

                        ;; M1 is the candidate point:
                        zxm = ylptmp1(0)
                        zym = ylptmp1(1)
                        za0 = ai
                        zb0 = bi
                        za1 = -1.0 / za0
                        zb1 = zxm - za1 * zym
                        zyp = -(zb1 - zb0) / (za1 - za0)
                        zxp = za0 * zyp + zb0
                        zd1 = (zxm - zxp) * (zxm - zxp) + \
                              (zym - zyp) * (zym - zyp)

                        ;; M2 is the candidate point:
                        zxm = ylptmp2(0)
                        zym = ylptmp2(1)
                        za1 = -1. / za0
                        zb1 = zxm - za1 * zym
                        zyp = -(zb1 - zb0) / (za1 - za0)
                        zxp = za0 * zyp + zb0
                        zd2 = (zxm - zxp) * (zxm - zxp) + \
                              (zym - zyp) * (zym - zyp)
                    end if

                    ;; choose the smallest (zd1, zd2)
                    if (zd2 .le. zd1) then
                        ;; use M2
                        rxx(nn) = ylptmp2(0)
                        ryy(nn) = ylptmp2(1)
                    else
                        ;; use M1
                        rxx(nn) = ylptmp1(0)
                        ryy(nn) = ylptmp1(1)
                    end if
                end if

                nn = nn + 1
                if (nn .ge. jpseg) then
                    error_msg("f", diag_script, "", "too many " + \
                              "sections. Increase value of jpseg " + \
                              "in code or reduce section length.")
                end if
                rxx(nn) = yypt(jk, 0)
                ryy(nn) = yypt(jk, 1)
            end do

            ;; segment jseg is a line between (rxx(jseg),ryy(jseg))
            ;; and (rxx(jseg+1),ryy(jseg+1))
            do jseg = 0, nn - 1
                ii0 = rxx(jseg)
                ij0 = ryy(jseg)
                if (rxx(jseg) .eq. rxx(jseg + 1)) then
                    ;; meridional segment, use U velocity
                    do jtime = 0, ntimes(idx_mod) - 1
                        do jlev = 0, dimsizes(extru&lev) - 1
                            myval = extru(jtime, jlev, ij0 + ijst, ii0)
                            if (.not. ismissing(myval)) then
                                dvoltrpsum(idx_mod, jtime) = \
                                    dvoltrpsum(idx_mod, jtime) +\
                                    (myval * \
                                     e2u(ij0 + ijst, ii0) *\
                                     e3u(jlev, ij0 + ijst, ii0) *\
                                     norm_u)
                            end if
                        end do
                    end do
                else if (ryy(jseg) .eq. ryy(jseg + 1)) then
                    ;; zonal segment, use V velocity
                    do jtime = 0, ntimes(idx_mod) - 1
                        do jlev = 0, dimsizes(extrv&lev) - 1
                            myval = extrv(jtime, jlev, ij0, ii0 + iist)
                            if (.not. ismissing(myval)) then
                                dvoltrpsum(idx_mod, jtime) = \
                                    dvoltrpsum(idx_mod, jtime) +\
                                    (myval * \
                                     e1v(ij0, ii0 + iist) *\
                                     e3v(jlev, ij0, ii0 + iist) *\
                                     norm_v)
                            end if
                        end do
                    end do
                else
                    ;; should never happen !
                    error_msg("f", diag_script, "", "rxx(jseg), " + \
                              "ryy(jseg), rxx(jseg + 1), " + \
                              "ryy(jseg + 1)" + nl + \
                              rxx(jseg) + ", " + \
                              ryy(jseg) + ", " + \
                              rxx(jseg + 1) + ", " + \
                              ryy(jseg + 1))
                end if
                end if
            end do

            ;; clean
            delete([/extru, latu, lonu, e2u, e3u, \
                     extrv, latv, lonv, e1v, e3v, \
                     yypt, rxx, ryy, yypt/])
        else
            info_output("skipped transport calculation for model " + \
                        models@name(idx_mod) + " because not all " + \
                        "sizes are known for project " + \
                        models@project(idx_mod) + ".", verbosity, 1)
        end if
    end do
end

;;#############################################################################
;; Create plot
;;#############################################################################
begin
    ;; go from m3/s to Sverdrup
    dvoltrpsum = dvoltrpsum / 1000000

    ;; define plot settings
    res_xy = True
    res_xy@tiMainString    = "Volume transport for section " + \
                             diag_script_info@section_name
    res_xy@tiXAxisString   = "date (yyyy-mm)"
    res_xy@tiYAxisString   = "volume transport (Sv)"
    res_xy@trXMinF = min(times) - 0.05 * (max(times) - \
                                          min(times))
    res_xy@trXMaxF = max(times) + 0.25 * (max(times) - \
                                          min(times))
    res_xy@tmXTOn          = False
    res_xy@tmXBLabelAngleF = 45
    res_xy@tmXBLabelJust   = "CenterRight"
    res_xy@tmYLFormat      = "f"      ; remove trailing ".0"

    ;; format time labels, there must be an easier way of doing this
    ;; find begin and end date, convert to year, month etc, then create
    ;; arrays with years, months, convert them back to "days since" and
    ;; let time_axis_labels figure out the steps
    mm_time = round((/floor(min(times)), ceil(max(times))/), 3)
    mm_time@units = times@units
    ;; adjust dates
    utc_date = toint(cd_calendar(mm_time, 0))
    utc_date(1, 0) = utc_date(1, 0) + 1
    nmonths = (utc_date(1, 0) - utc_date(0, 0)) * 12 + 1
    years  = new(nmonths, integer)
    months = new(nmonths, integer)
    ones   = new(nmonths, integer)
    zeros  = new(nmonths, integer)
    ones(:)  = 1
    zeros(:) = 0
    iy = utc_date(0, 0)
    im = 1
    do i = 0, nmonths - 1
        years(i) = iy
        months(i) = im
        im = im + 1
        if (im .gt. 12) then
            im = 1
            iy = iy + 1
        end if
    end do
    adj_mm_time = round(cd_inv_calendar(years, months, ones, zeros, \
                                        zeros, zeros, times@units, 0), 3)
    adj_mm_time@units = times@units
    if (isatt(adj_mm_time, "_FillValue")) then
        delete(adj_mm_time@_FillValue)
    end if
    restick                = True
    restick@ttmFormat      = "%Y-%N"
    restick@ttmMajorStride = 12
    restick@ttmMinorStride = 1
    time_axis_labels(adj_mm_time, res_xy, restick)

    ;; add line to plot
    debuginfo = False
    xy_line(wks_xy, dvoltrpsum, times, 0, res_xy, debuginfo)

    ;; clean
    delete([/adj_mm_time, dvoltrpsum, mm_time, months, ntimes, ones, \
             res_xy, restick, times, utc_date, years, zeros/])
end

;;#############################################################################
;; References Information
;;#############################################################################
begin
    leave_msg(diag_script, "", 4)
end
