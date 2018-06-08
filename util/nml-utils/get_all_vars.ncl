;;#############################################################################
;; UTILITY SCRIPT TO GET ALL VARIABLES USED BY THE TOOL
;; Author: Mattia Righi (DLR, Germany)
;;#############################################################################
;;
;; Description
;;    Check all the namelists and creates a list with all variables required by
;;    the tool, together with the corresponding MIP and frequency.
;;
;; Caveats
;;    The scripts only considers the CMIP5 models.
;;
;; Modification history
;;    20170203-A_righ_ma: written.
;;
;;#############################################################################

undef("emptystr")
function emptystr(nn[1]:integer)
local funcname, scriptname, verbosity, ii
begin

    str = " "
    do ii = 1, nn - 1
        str = str + " "
    end do

    return(str)

end

;;#############################################################################

undef("unique_values")
function unique_values(arr)
local loc, idx
begin

    loc = arr
    loc@_FillValue = default_fillvalue(typeof(loc))
    do ii = 0, dimsizes(loc) - 1
        idx = ind(loc(ii).eq.loc)
        if (.not.all(ismissing(idx))) then
            if (dimsizes(idx).gt.1) then
                loc(idx(1::)) = loc@_FillValue
            end if
        end if
        delete(idx)
    end do
    
    idx = ind(.not.ismissing(loc))

    out = loc(idx)
    return(out)

end

;;#############################################################################

undef("isfilepresent_esmval")
function isfilepresent_esmval(file_path[1]: string)
local ASCII_ZERO,   FILEEXISTS_MIN_VERSION,  file_is_present,\
    ncl_version,  temp_char,               temp_int_1,\
    temp_int_2,   temp_str,                use_fileexists
begin
    ASCII_ZERO = 48  ; Used for char to int single digit conversion
    FILEEXISTS_MIN_VERSION = (/6, 2, 1/)  ;; 'fileexists' only available
                                          ;; in 6.2.1+

    ;; Get NCL version, to determine whether to use 'fileexists' or
    ;; 'isFilePresent'

    ncl_version = (/0, 0, 0/)  ; Start with an undefined value
    temp_str = get_ncl_version()
    temp_char = stringtochar(temp_str)
    ncl_version(0) = chartointeger(temp_char(0))\
                      - ASCII_ZERO  ; Convert ascii to 0-9
    ncl_version(1) = chartointeger(temp_char(2))\
                      - ASCII_ZERO  ; Convert ascii to 0-9
    ncl_version(2) = chartointeger(temp_char(4))\
                      - ASCII_ZERO  ; Convert ascii to 0-9

    ;; Determine whether version is sufficient for use of 'fileexists'
    temp_int_1 = ncl_version(0) * 100\
                 + ncl_version(1) * 10\
                 + ncl_version(2)
    temp_int_2 = FILEEXISTS_MIN_VERSION(0) * 100\
                 + FILEEXISTS_MIN_VERSION(1) * 10\
                 + FILEEXISTS_MIN_VERSION(2)
    use_fileexists = temp_int_1.ge.temp_int_2

    ;; Set default return value
    file_is_present = False

    if (use_fileexists) then  ; Only for version >= 6.2.1
        if (fileexists(file_path)) then
            file_is_present = True
        end if
    else
        if (isfilepresent(file_path)) then  ; Only for version < 6.2.1
            file_is_present = True
        end if
    end if
    return file_is_present
end

;;#############################################################################

begin

    nmls = systemfunc("ls ../../nml/*.xml")
    nnmls = dimsizes(nmls)

    outfile = "ESMValTool_vars_CMIP5.dat"
    system("rm -f " + outfile)

    print("")
    print(nnmls + " namelists found")
    print("")

    ;; ### Loop over namelists ################################################
    do nn = 0, nnmls - 1

        print("Reading namelist " + nmls(nn))

        data = readAsciiTable(nmls(nn), 1, "string", 0)
        data := data(:,0)
        ndata = dimsizes(data)

        ;; ### Read the <MODELS> section ######################################
        id1 = 0
        do ii = 0, ndata - 1
            if (isStrSubset(data(ii), "<MODELS>"))
                id1 = ii
                break
            end if
        end do
        id2 = 0
        do ii = 0, ndata - 1
            if (isStrSubset(data(ii), "</MODELS>"))
                id2 = ii
                break
            end if
        end do
        if (id1.eq.0 .or. id2.eq.0) then
            print("  This namelist does not contain a MODEL section, skipping.")
            delete(data)
            continue
        end if

        subdata = data(id1+1:id2-1)
        nsubdata = dimsizes(subdata)
        do ii = 0, nsubdata - 1
            if (isStrSubset(subdata(ii), "<!--")) then
                continue
            end if
            if (isStrSubset(subdata(ii), "<model>")) then
                if (isStrSubset(str_get_field(subdata(ii), 2, " "), \
                    "CMIP5")) then
                    info = str_get_field(subdata(ii), 4, " ") + "-" + \
                        str_get_field(subdata(ii), 5, " ")
                    if (.not.isdefined("freqmip4all")) then
                        freqmip4all = info
                    else
                        freqmip4all := array_append_record(freqmip4all, info, 0)
                    end if
                    delete(info)
                end if
            end if
        end do
        if (isdefined("freqmip4all"))
            freqmip4all := unique_values(freqmip4all)
        end if
        delete(subdata)

        ;; ### Read the <DIAGNOSTICS> section #################################
        EOF = False
        id1 = 0
        id2 = 0
        do while (.not.EOF) 

            do ii = id1+1, ndata - 1
                if (isStrSubset(data(ii), "<diag>"))
                    id1 = ii
                    break
                end if
                if (ii.eq.ndata - 1) then
                    EOF = True
                end if
            end do

            do ii = id2+1, ndata - 1
                if (isStrSubset(data(ii), "</diag>"))
                    id2 = ii
                    break
                end if
                if (ii.eq.ndata - 1) then
                    EOF = True
                end if               
            end do

            if (EOF) then
                continue
            end if

            ;; Search for variable and variable-specific models
            subdata = data(id1+1:id2-1)
            nsubdata = dimsizes(subdata)
            varfound = False

            do ii = 0, nsubdata - 1
                if (isStrSubset(subdata(ii), "variable_def_dir")) then
                    continue
                end if
                if (isStrSubset(subdata(ii), "<variable")) then
                    qq = str_index_of_substr(subdata(ii), "<variable", -1) + strlen("<variable") - 1
                    tmp = str_get_cols(subdata(ii), qq, strlen(subdata(ii)))
                    tmp := str_squeeze(str_get_field(str_get_field(tmp, 2, ">"), 1, "<"))
                    if (varfound) then
                        infov = infov + "," + tmp
                    else
                        infov = tmp
                    end if
                    delete(tmp)
                    delete(qq)
                    if (isStrSubset(subdata(ii), "MIP=")) then
                        qq = str_index_of_substr(subdata(ii), "MIP=", -1) + strlen("MIP=")
                        tmp = str_get_cols(subdata(ii), qq, strlen(subdata(ii)))
                        mip_var_def = str_get_field(tmp, 1, str_get_dq())
                        delete(tmp)
                        delete(qq)
                    end if
                    varfound = True
                end if
                if (isStrSubset(subdata(ii), "<model>")) then
                    if (isStrSubset(str_get_field(subdata(ii), 2, " "), \
                        "CMIP5")) then
                        infom = str_get_field(subdata(ii), 4, " ") + "-" + \
                            str_get_field(subdata(ii), 5, " ")
                        if (.not.isdefined("freqmip")) then
                            freqmip = infom
                        else
                            freqmip := array_append_record(freqmip, infom, 0)
                        end if
                        delete(infom)
                    end if
                end if
            end do

            ;; Variable found, check if derived and add to the list
            if (isdefined("infov")) then
                vars = (/"DUMMY"/)
                nvars = str_fields_count(infov, ",")
                do vv = 0, nvars - 1
                    if (nvars.gt.1) then
                        cvar = str_get_field(infov, vv+1, ",")
                    else
                        cvar = infov
                    end if
                    deffile = "../../variable_defs/" + cvar + ".ncl"
                    if (.not.isfilepresent_esmval(deffile)) then
                        print("  Warning! Variable definition file for " + cvar + \
                            " not found. Assuming this variable is not derived.")
                        derived = "none"
                    else
                        derived = systemfunc("more " + deffile + " | grep 'Requires:' " + \
                            "| awk -F 'Requires:' '{print $2}'")
                            derived = str_squeeze(derived)
                    end if
                    if (all(derived.ne.(/"none", "None"/))) then
                        ndvars = str_fields_count(derived, ",")
                        do dv = 0, ndvars - 1
                            xx = str_squeeze(str_get_field(\
                                str_get_field(derived, dv+1, ","), 1, ":"))
                            vars := array_append_record(vars, xx, 0)
                        end do
                    else
                        vars := array_append_record(vars, cvar, 0)
                    end if
                    delete(cvar)
                end do
                vars := vars(ind(vars.ne."DUMMY"))
                delete(infov)
                do vv = 0, dimsizes(vars) - 1
                    if (isdefined("freqmip")) then
                        if (isdefined("mip_var_def")) then
                            out = \
                                str_sub_str(freqmip, "MIP_VAR_DEF", mip_var_def)
                        else
                            out = freqmip
                        end if
                        do ff = 0, dimsizes(out) - 1
                            freq = str_get_field(out(ff), 1, "-")
                            mip = str_get_field(out(ff), 2, "-")
                            str = vars(vv) + "," + freq + "," + mip
                            system("echo '" + str + "' >> " + outfile)
                        end do
                        delete(out)
                    end if
                    if (isdefined("freqmip4all")) then
                        if (isdefined("mip_var_def")) then
                            out = \
                                str_sub_str(freqmip4all, "MIP_VAR_DEF", mip_var_def)
                        else
                            out = freqmip4all
                        end if
                        do ff = 0, dimsizes(out) - 1
                            freq = str_get_field(out(ff), 1, "-")
                            mip = str_get_field(out(ff), 2, "-")
                            str = vars(vv) + "," + freq + "," + mip
                            system("echo '" + str + "' >> " + outfile)
                        end do
                        delete(out)
                    end if
                end do
                delete(vars)
            end if
            delete(subdata)
            if (isdefined("freqmip")) then
                delete(freqmip)
            end if

        end do
        delete(data)
        if (isdefined("freqmip4all"))
            delete(freqmip4all)
        end if

    end do

    ;; Remove Stderr variables (observations only!)
    system("sed -i '/Stderr/d' "+ outfile)

    ;; Manually remove some non-conformal mips
    system("sed -i '/seawifs/d' "+ outfile)      ; namelist GlobalOcean
    system("sed -i '/takahashi14/d' "+ outfile)  ; namelist GlobalOcean
    system("sed -i '/ETH/d' "+ outfile)       
    system("sed -i '/?mon/d' "+ outfile)
    system("sed -i '/obs/d' "+ outfile)
    system("sed -i '/reanalysis/d' "+ outfile)
  
    ;; Add some header info
    str = "List of CMIP5 variables required by the ESMValTool"
    system("echo " + str + " > tmp.dat")
    str = "Based on git " + systemfunc("git log | head -n 1")
    system("echo " + str + " >> tmp.dat")
    str = systemfunc("git log | head -n 3 | tail -n 1")
    system("echo " + str + " >> tmp.dat")
    system("echo ' ' >>  tmp.dat")

    ;; Remove duplicates and sort alphabetically
    system("more " + outfile + " | sort | uniq >> tmp.dat")
    system("mv tmp.dat " + outfile)

end
