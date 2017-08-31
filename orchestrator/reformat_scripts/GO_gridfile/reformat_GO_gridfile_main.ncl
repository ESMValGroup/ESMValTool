
;*********************************************************************
;This is the main ncl program called from the script file reformat.py
;
;contact  I. Cionni     i.cionni@niwa.co.nz
;
; Changes
; 
; 2009-06-11   CAF  added date range to output filename
; 2009-06-12   CAF  fixed dimension problem with T2Dz
; 2010-09-15   HS/SWW   fixed lev->plev through out file, changed
;                       input file, fil_in
;                       and general tidy-up including comments
; --------------------------------------------------------------------
; attribute.ncl
; 2010-11-30   HS   handle multiple input files and read data in blocks with
;                   limited size (maxDataSize). Data gets fixed if
;                   project_name & model specified file exists or if its field
;                   and rank do not match. Output a file with only selected
;                   time range, without "Define Mode".
; 2011-08-26   HS   bug fix - an error on selecting proper time range.
; 2012-06-07   HS   change the location of fix_file (function fix_data).
; 2013-05-17 op_kg  fix_data(...,"ocean-latlon") moved out of main loop to
;                   avoid crashing when the first file is not in the selected
;                   time period and thus fout (needed by fix_data) is still empty
; 2013-06-01 op_kg  introduced fix_file_default & default_fix
;*********************************************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "interface_scripts/messaging.ncl"
load "interface_scripts/auxiliary.ncl"

;*********************************************************************
; Check that all these functions are in the same directory of attribute.ncl
;*********************************************************************

;;;load "reformat_scripts/original/time_conform.ncl"
load "reformat_scripts/constants.ncl"
load "reformat_scripts/GO/depth_conform.ncl"
load "reformat_scripts/GO/reformat_GO_func.ncl"
;;;load "reformat_scripts/original/coord_function.ncl"
;;load "reformat_scripts/default/reformat_default_func.ncl"
;;;load "reformat_scripts/original/att_function.ncl"
;;load "reformat/ORIGINAL/global_attribute.ncl"
load "interface_scripts/messaging.ncl"



;; Note - loadscript arguments are only executed _after_ a end-statement,
;;        hence the need for the additional begin/end block below
begin
    
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering reformat_GO_main.ncl", verbosity, 2)
    project = getenv("ESMValTool_project")
    model= getenv("ESMValTool_model")
    info_output("In reformatGO_main.ncl", verbosity, 3)
    ;; If ProjectName_ModelName.ncl exists, then call it.
    ;; Otherwise use ProjectName_default.ncl, if it exists.
    fix_file = "./reformat_scripts/fixes/"  + project + "_" + model + ".ncl"
    fix_file_default = "./reformat_scripts/fixes/"  + project + "_default.ncl"
    info_output("fix_file " + fix_file, verbosity, 3)
    info_output("fix_file_default " + fix_file_default, verbosity, 3)
    data_to_fix = isfilepresent_esmval(fix_file)
    default_fix = isfilepresent_esmval(fix_file_default) 
    if data_to_fix then
        info_output(" loading "+fix_file, verbosity, 2)
        loadscript(fix_file)
    else if (default_fix .and. .not. data_to_fix) then
        info_output(" loading "+fix_file_default, verbosity, 2)
        loadscript(fix_file_default)
    end if
    end if

end

begin

    fields_list = (/"TO2Ms: Monthly-mean 2-d ocean or sea ice data (longitude, latitude, time:month)","TO3M: Monthly mean 3d ocean data"/)

    ;;===========================================================
    ;; Get environment variables
    ;;===========================================================
 verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    outfile_path       = getenv("ESMValTool_outfile_fullpath")
    field              = getenv("ESMValTool_field")
    infile_path        = getenv("ESMValTool_infile_path")
    variable           = getenv("ESMValTool_variable")
    max_data_blocksize = stringtofloat(getenv("ESMValTool_max_data_blocksize"))
    start_date         = stringtofloat(getenv("ESMValTool_start_year"))
    end_date           = stringtofloat(getenv("ESMValTool_end_year"))
    references         = getenv("ESMValTool_references")
    contact            = getenv("ESMValTool_contact")
;;    experiment         = getenv("ESMValTool_experiment")
    experiment         = getenv("ESMValTool_case_name")
    history            = getenv("ESMValTool_history")
    comment            = getenv("ESMValTool_comment")
    source             = getenv("ESMValTool_source")
    areafile = getenv("ESMValTool_areafile_path")  ;;defined in main nml (e.g. for project CMIP5_gridfile)
    ;; Find GO name 
    GO_var = find_name(variable)

    info_output("outfile_path " + outfile_path, verbosity, 3)
    info_output("start_date " + start_date, verbosity, 3)
    info_output("variable " + variable, verbosity, 3)
    info_output("experiment " + experiment, verbosity, 3)
    fn_out = outfile_path
    if (isfilepresent_esmval(fn_out)) then
        system("/bin/rm " + fn_out)
    end if

    ;; Check field of input-file
  
    att_field = str_match(fields_list, field)
    if (ismissing(att_field)) then
        print("fatal:" + field \
              + " is wrong field name: please check your setting!")
        status_exit(1)
    end if

    ;;===========================================================
    ;; Input and output files
    ;;===========================================================
    ;; Get input-file(s): i.e. all .nc files in input directory with variable
    ;; and then experiment in hteir name
    info_output("infile_path " + infile_path, verbosity, 3)
    info_output("variable " + variable, verbosity, 3)
    info_output("GO variable " + GO_var, verbosity, 3)
    info_output(infile_path +  GO_var + "*.nc", verbosity, 3)
    fils = systemfunc("ls " + infile_path  + GO_var + "*" + model + "*" + experiment + "*.nc" )  ; File paths
;;    print ( fils )

    if(all(ismissing(fils))) then
        print ("fatal: Files not found. " + infile_path)
        status_exit(1)
    end if

    ;; Create output file
    fout = addfile(fn_out, "c")
    ;; Make time an UNLIMITED dimension
    filedimdef(fout, "time", -1, True)

    first_write = True
    file_pointer = 0

    ;;===========================================================
    ;; Main loop - go through file(s)
    ;;===========================================================
    dim_fils = dimsizes(fils)
    do i = 0, dim_fils - 1
        info_output("file names: " + fils(i), verbosity, 3)
        fin = addfile(fils(i), "r")
        vNames = getfilevarnames(fin)
        info_output("vNames: " + vNames, verbosity, 3)

        if (sum(where(vNames .eq. GO_var, 1, 0)) .eq. 0) then
            print ("fatal: var = " + GO_var + " is not in the file")
            status_exit(1)
        end if

        dNames0 = getfilevardims(fin, GO_var)
        vards = getfilevardimsizes(fin, GO_var)
        rank = dimsizes(vards)
        varType = getfilevartypes(fin, GO_var)

        ;; Calculate variable data size
        ;; datasize = systemfunc("ls -s " + fils(i) + " | cut -d ' ' -f1")
        if (varType .eq. "double") then
            datasize = 8.
        else
            datasize = 4.
        end if

        do j = 0, rank - 1, 1
            datasize = datasize * vards(j)
        end do

        datasize = datasize / 1000000  ; In Mbyte
       
	if (inlist("time", dNames0)) then  
	       time_name = "time"
	else if (inlist("time_counter", dNames0)) then  
	       time_name = "time_counter"
	else if (inlist("TIME", dNames0)) then  
	       time_name = "TIME"
        else
              print("fatal:Couldn't identify time dimension name in reformatGO_main.ncl")
              status_exit(1)
        end if
        end if
        end if
        ;; Set time start/end indices for reading blocks
        time_len = vards(ind(dNames0 .eq. time_name))
        info_output("Reading times: " + time_len, verbosity, 3)
        if (datasize .gt. max_data_blocksize) then
            nblocks = floattoint(ceil(datasize / max_data_blocksize))
            time_block = time_len / nblocks
            if (mod(time_len, nblocks) .ne. 0) then
                time_block = time_block + 1
            end if
            i_st = new(nblocks, "integer")
            i_ed = i_st
            do j = 0, nblocks - 1
                i_st(j) = j * time_block
                i_ed(j) = (j + 1) * time_block - 1
            end do
            i_ed(nblocks - 1) = time_len - 1
        else
            i_st = 0
            i_ed = time_len - 1
            nblocks = 1
        end if
        delete(time_len)
        delete(vards)
    
	;; Query spatial coordinate names
	if (inlist("j", dNames0)) then  
	    lat_name = "j"
	else if (inlist("y", dNames0)) then  
	    lat_name = "y"
	else if (inlist("lat", dNames0)) then  
	    lat_name = "lat"
        else if (inlist("rlat", dNames0)) then
            lat_name = "rlat"
        else if (inlist("LATITUDE", dNames0)) then
             lat_name = "LATITUDE"
        else if (inlist("nav_lat", dNames0)) then
            lat_name = "nav_lat"
        else
            print("fatal:Couldn't identify lat-coordinate name in attribute.ncl")
            status_exit(1)
	end if
        end if
	end if
        end if
	end if
        end if
        lev_name = "missing"
        if ( rank .ge. 3 ) then 
	 if (inlist("depth", dNames0)) then  
	    lev_name = "depth"
	 else if (inlist("DEPTH", dNames0)) then  
	    lev_name = "DEPTH"
	 else if (inlist("deptht", dNames0)) then  
	    lev_name = "deptht"
	 else if (inlist("lev", dNames0)) then  
	    lev_name = "lev"
	 else if (inlist("lon", dNames0)) then  
	    lev_name = "plev" 
         else
            print("Warning:Couldn't identify depth-coordinate name in attribute.ncl")
            
	 end if
	 end if
         end if
         end if
	 end if
        end if
         if (inlist("i", dNames0)) then  
	    lon_name = "i"
  	 else if (inlist("x", dNames0)) then  
	    lon_name = "x"
	 else if (inlist("lon", dNames0)) then  
	    lon_name = "lon"
         else if (inlist("LONGITUDE", dNames0)) then  
	    lon_name = "LONGITUDE"
         else if (inlist("rlon", dNames0)) then
            lon_name = "rlon"
         else if (inlist("nav_lon", dNames0)) then
            lon_name = "nav_lon"
         else
            print("fatal:Couldn't identify lon-coordinate name in attribute.ncl")
            status_exit(1)
	 end if
	 end if
         end if
	 end if 
        end if
         end if
        info_output("nblocks " + nblocks, verbosity, 3)
        info_output("var " + GO_var, verbosity, 3)
   
     
        ;; Read in blocks
        do jj = 0, nblocks - 1

            ;; Get a block of data
            if (rank .eq. 4) then

                varData0 = fin->$GO_var$($time_name$|i_st(jj):i_ed(jj), $lev_name$|:, $lat_name$|:, $lon_name$|:)
            end if
            if (rank .eq. 3) then
                if isStrSubset(field, "s") then
                    varData0 = fin->$GO_var$($time_name$|i_st(jj):i_ed(jj), $lat_name$|:, $lon_name$|:)
                else
                    varData0 = fin->$GO_var$($time_name$|i_st(jj):i_ed(jj), $lev_name$|:, $lat_name$|:)
                end if
            end if
            if (rank .eq. 2) then
                varData0 = fin->$GO_var$($time_name$|i_st(jj):i_ed(jj), $lat_name$|:)
            end if
            if (rank .eq. 1) then
                varData0 = fin->$GO_var$($time_name$|i_st(jj):i_ed(jj))
            end if
              
; Check for attribute _FillValue mis-spelt as _Fillvalue;
;this is a common mistake in PlankTOM output 
            if(isatt(varData0, "_Fillvalue")) then
             fill = varData0@_Fillvalue
             varData0@_FillValue = fill
             delete_VarAtts(varData0,"_Fillvalue" )
            end if
        
       
;            print ("varData0 fill is " + varData0@_FillValue )
;            if(isatt(varData0, "_Fillvalue")) then
;              varData0 = where(varData0.ge.varData0@_Fillvalue, 1.e20, varData0)
;            end if
;            varData0@_FillValue = 1.e20
           
;            print ("and data is " + varData0(1,1,1) )
        
      ;; If project_name_model.ncl exists, then call fix_data.
            if data_to_fix then
                if isdefined("fix_data") then
                    info_output("Calling fix_data", verbosity, 3)
                    temp = fix_data(varData0, "default")
                    delete(varData0)
                    varData0 = temp
                    delete(temp)
                end if
            end if
             

            ;; Select time dimension and check for range

	   
            tim = varData0&$time_name$
                
            time0 = tim
;;;            time0 = time_conform(tim)
            tim_ut = cd_calendar(time0, 0)  ; NB: Only check for year, not months
     
        

            ind_t = ind(tim_ut(:,0) .ge. start_date .and. tim_ut(:,0) .lt. (end_date + 1))
;;            print("fatal: " + tim_ut(:,0) )
;;            print("start: " +  start_date )
;;            print("end: " + end_date )
            delete(tim)
            delete(time0)
            delete(tim_ut)
             
            if (.not.all(ismissing(ind_t)))
                if (rank .eq. 4) then
                    varData1 = varData0(ind_t, :, :, :)
                end if
                if (rank .eq. 3) then
                    varData1 = varData0(ind_t, :, :)
                end if
                if (rank .eq. 2) then
                    varData1 = varData0(ind_t, :)
                end if
                if (rank .eq. 1) then
                    varData1 = varData0(ind_t)
                end if
                delete(varData0)
                dNames = getvardims(varData1)
                info_output("dNames  varData1 is: " + dNames, verbosity, 3)

                ;; double2float, convert2flt
 ;               if (varType .eq. "double") then
;                    temp = tofloat(varData1)
;                    delete(varData1)
;                    varData1 = temp
;                    delete(temp)
;                end if
               
                
                varData = (/varData1/)
 ;               varData1!0=dNames(0)
;                varData1!1=dNames(1)
;                varData1!2=dNames(2)
;                varData1!3=dNames(3) 
                copy_VarAtts(varData1, varData)
               
                ;; dimensions: check name, type and attributes
                dNames = getvardims(varData1)
                info_output("dNames is: " + dNames, verbosity, 3)
 
                dDNames = getvardims(dNames(1))               
;;                mylats=fin->lat(:,:) 
;;                print("mylats is "+ mylats(10,10))
                do k = 0, rank - 1, 1
                    varDimX = varData1&$dNames(k)$
                    ;; Dimensions also need to be float
                    if (.not.isfloat(varDimX)) then
                        temp = tofloat(varDimX)
                        delete(varDimX)
                        varDimX = temp
                        delete(temp)
                    end if
                    if (dNames(k) .eq. time_name) then
;;                        time = time_conform(varDimX)
 ;;                       time = coord_list("time", time)
                        varData!k = "time"
                          varData&time = varData1&$time_name$
;;                        delete(time)
                    end if
                    if (dNames(k) .eq. lev_name) then
              
;;                        varDimX = depth_conform(varDimX)
;;                        depth = coord_list("plev", varDimX)
                        varData!k = "depth"
                          varData&depth = varData1&$lev_name$
;;                        delete(depth)
                    end if
                    if (dNames(k) .eq. lat_name) then
;;                        lat = coord_list("lat", varDimX)
                        varData!k = "lat"
                       varData&lat =  varData1&$lat_name$
;;                        delete(lat)
                    end if
                    if (dNames(k) .eq. lon_name) then
;;                        lon = coord_list("lon", varDimX)
                        varData!k = "lon"
                        varData&lon =  varData1&$lon_name$
;;                        delete(lon)
                    end if
                    delete(varDimX)
                end do
                delete(varData1)
                delete(dNames)
                info_output("Lat size is " + dimsizes(dimsizes(varData&lat)), verbosity, 3)
                info_output("Lat size is " + dimsizes(dimsizes(varData&lat)), verbosity, 3)

                ;; Computes the zonal average of a variable [4D -> 3D data (*2*z)]
                if ((rank .eq. 4) .and. (isStrSubset(field, "2")) \
                    .and. (isStrSubset(field, "z"))) then
                    temp = dim_avg_Wrap(varData)
                    delete(varData)
                    varData = temp
                    delete(temp)
                    rankN = 3
                else
                    rankN = rank
                end if
    ;; Check units
              check_units(varData, variable)
              units = varData@units
              info_output("Units here are " + units, verbosity, 3)
         
                ;;=========================================================
                ;; Write to output file, fout
                ;;=========================================================
                if first_write then

                    ;; Create global attributes of the file
                    fAtt = True                         ; Assign file attributes
                    fAtt@fields = att_field
 ;;                   global_attribute(model, references, contact, references,\
 ;;                                    ensemble, history, comment, source, fAtt)
 ;;                   fileattdef(fout, fAtt)      ; Copy file attributes

                    ;; Set variable attributes
 ;;                   varData = var_list(GO_var, varData)
                    cell_methods = ""
                    if isStrSubset(field, "z") then
                        cell_methods = "lon:mean (zonal mean), "
                    end if
                    if isStrSubset(field, "D") then
                        cell_methods = cell_methods + "time:mean (of day)"
                    end if
                    if isStrSubset(field, "M") then
                        cell_methods = cell_methods + "time:mean (of calendar months)"
                    end if
                    if isStrSubset(field, "I") then
                        cell_methods = cell_methods + "time:istantaneous"
                    end if
                    if isStrSubset(field, "A") then
                        cell_methods = cell_methods + "time:mean (of year)"
                    end if
                    varData@cell_methods = cell_methods
                    varData@units = units
                   
  
                    info_output("fill value now is " + varData@_FillValue, verbosity, 3)
                    fout->$variable$ = varData
 

                    first_write = False

                else
                   info_output("size of varData " + dimsizes(varData) , verbosity, 3)
              
                    end_point = file_pointer + dimsizes(ind_t) - 1

                    fout->time(file_pointer:end_point) = (/varData&time/)
                    info_output("lhs size " + dimsizes( fout->$variable$(file_pointer:end_point, :, :) ), verbosity, 3)

                    if (rankN .eq. 4) then
                        fout->$variable$(file_pointer:end_point, :, :, :) = (/varData/)
                    end if
                    if (rankN .eq. 3) then
                        fout->$variable$(file_pointer:end_point, :, :) = (/varData/)
                    end if
                    if (rankN .eq. 2) then
                        fout->$variable$(file_pointer:end_point, :) = (/varData/)
                    end if
                    if (rankN .eq. 1) then
                        fout->$variable$(file_pointer:end_point) = (/varData/)
                    end if
                end if
                file_pointer = file_pointer + dimsizes(ind_t)

                delete(varData)
            else
                delete(varData0)
            end if

            delete(ind_t)

        end do  ; Loop (jj - blocks)

        delete(dNames0)
        delete(vNames)
        delete(i_st)
        delete(i_ed)

        keep_fin = fin ;keep latest fin to extract grid parameters for "ocean-latlon"
        delete(fin)
    end do  ; Loop (i - files)

    ;; If project_name_model.ncl exists, then call fix_data.
    ;; The switch "ocean-latlon" controls additional operations in fix_data.
    ;; WARNING: Only the latest input file is passed to fix_data for
    ;;   extraction of grid parameters. This implies that all input files 
    ;;   for the current model are on the same grid, but this is not checked
    ;;
    ;; 1st priority: use model specific fix script 
    ;;   (Note: Must include possible default operations!)
    info_output("data_to_fix: " + data_to_fix, verbosity, 3)
    info_output("fix_data: " +  isdefined("fix_data"), verbosity, 3)
    info_output("default_fix: " + default_fix, verbosity, 3)
    info_output("fix_default: " +  isdefined("fix_default"), verbosity, 3)

    if data_to_fix then
    if isdefined("fix_data") then
        info_output("Calling fix_data", verbosity, 3)
        file_handles = new(2, "file")
        file_handles(0) = keep_fin
        file_handles(1) = fout
        asciiwrite("./interface_data/lat_name.txt", lat_name)
        asciiwrite("./interface_data/lon_name.txt", lon_name)
        if ( rankN .eq. 4 ) then 
          asciiwrite("./interface_data/depth_name.txt", lev_name)
        end if
        asciiwrite("./interface_data/fn_out.txt", fn_out)
        temp = fix_data(file_handles, "ocean-latlon" )
        delete(temp)
    end if
    end if

    ;; 2nd priority: use project specific default
    ;;  (e.g. for passing non-gaussian grid parameters to the plot section)
    if (default_fix .and. .not. data_to_fix) then
    if isdefined("fix_default" )
      if (areafile .ne.  "null") then
        file_handles = new(2, "file")
        file_handles(0) = keep_fin
        file_handles(1) = fout
        asciiwrite("./interface_data/lat_name.txt", lat_name)
        asciiwrite("./interface_data/lon_name.txt", lon_name)
        if ( rankN .eq. 4 ) then 
          asciiwrite("./interface_data/depth_name.txt", lev_name)
        end if
        asciiwrite("./interface_data/fn_out.txt", fn_out)
        temp = fix_default(file_handles, "ocean-latlon")
        delete(temp)
      end if
    end if
    end if
   
    delete(keep_fin)
    delete(fils)
    delete(fout)

    ;; If first_write = True, then there is no data within the specified
    ;; date range in the input file(s). So delete output file again and
    ;; print a "fatal"-message so that the run will crash.
    if first_write then
        system("/bin/rm " + fn_out)
        print("fatal: No data was within the specified date range. Please check your setting!")
        status_exit(1)
    end if
    info_output(">>>>>>>> Leaving reformat_GO_main.ncl", verbosity, 2)

end
