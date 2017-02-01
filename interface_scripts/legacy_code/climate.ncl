;
; 2009-03-10  CAF  now search for input variable instead of assuming it's
;                  index of 0
;                  also added hybrid level information to climate files
; 2009-07-17  CAF  write to $variable$ instead of $var$
; 2010-12-16  HS   file formats problems are dealt in earlier stages now.
;                  bug fix! move cell_methods to later part, so get used.
;                  check month of data@time instead of assuming it's correct.
;                  also use dim_avg_n_Wrap and removed unnecessary calls
; 2011-01-31  HS   bug fix - an error when start_year and end_year are same.
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./interface_scripts/info_output.ncl"

;; Temporary files update by Python on the fly 
load "./interface_data/ncl.interface"

begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering climate.ncl", verbosity, 2)

    error = systemfunc("rm -f " + mfile)  ; mfile fetched from ncl.interface 
    error = systemfunc("rm -f " + sfile)  ; sfile fetched from ncl.interface
    error = systemfunc("rm -f " + afile)  ; afile fetched from ncl.interface

    fout_mfile = addfile(mfile, "c")
    fout_sfile = addfile(sfile, "c")
    fout_afile = addfile(afile, "c")

    ; 'infilename' fetched from ncl.interface
    fin = addfile(infilename, "r")

    var = getfilevarnames(fin)

    ; 'base_variable' fetched from ncl.interface
    index = ind(var .eq. base_variable)
    data = fin->$var(index)$
    atts = getvaratts(fin)

    s = dimsizes(data)
    rank = dimsizes(s)
    sm = s
    sm(0) = 12
    ss = s
    ss(0) = 4
    sa = s
    sa(0) = 1
    datam = new(sm, typeof(data))
    datas = new(ss, typeof(data))
    dataa = new(sa, typeof(data))
    fatt_mfile = True
    fatt_sfile = True
    fatt_afile = True

    dimNames = getvardims(data)
    if (dimNames(0) .ne. "time")
        error_msg("f", "climate.ncl", "", "first dimension of " + var(0) + \
                  " in file " + infilename + " must be time")
    end if

    do i = 0, dimsizes(dimNames) - 1, 1
        datam!i = dimNames(i)
        datas!i = dimNames(i)
        dataa!i = dimNames(i)
    end do

    if (sum(where(dimNames .eq. "time", 1, 0)) .gt. 0) then
        datam&time = data&time(0:11)
        datas&time = data&time(0:11:3)
        dataa&time = data&time(0)
    end if

    ; Get month info from data@time
    time = data&time
    utc_date = cd_calendar(time, 0)
    utc_mm = floattointeger(utc_date(:, 1))
    delete(time)
    delete(utc_date)

    if (rank .eq. 1)
        do mn = 0, 11, 1
            ind_mm = ind(utc_mm .eq. (mn + 1))
            if (.not.all(ismissing(ind_mm)))
                if (dimsizes(ind_mm) .eq. 1)
                    datam(mn) = data(ind_mm)
                else
                    datam(mn) = dim_avg_Wrap(data(ind_mm))
                end if
            end if
            delete(ind_mm)
        end do
        datas(0) = (datam(0) + datam(1) + datam(11)) / 3.0
        datas(1) = dim_avg_Wrap(datam(2:4))
        datas(2) = dim_avg_Wrap(datam(5:7))
        datas(3) = dim_avg_Wrap(datam(8:10))
        dataa(0) = dim_avg_Wrap(datam)
    end if
    if (rank .eq. 2)
        do mn = 0, 11, 1
            ind_mm = ind(utc_mm .eq. (mn + 1))
            if (.not.all(ismissing(ind_mm)))
                if (dimsizes(ind_mm) .eq. 1)
                    datam(mn, :) = data(ind_mm, :)
                else
                    datam(mn, :) = dim_avg_n_Wrap(data(ind_mm, :), 0)
                end if
            end if
            delete(ind_mm)
        end do
        datas(0, :) = (datam(0, :) + datam(1, :) + datam(11, :)) / 3.0
        datas(1, :) = dim_avg_n_Wrap(datam(2:4, :), 0)
        datas(2, :) = dim_avg_n_Wrap(datam(5:7, :), 0)
        datas(3, :) = dim_avg_n_Wrap(datam(8:10, :), 0)
        dataa(0, :) = dim_avg_n_Wrap(datam, 0)
    end if
    if (rank .eq. 3)
        do mn = 0, 11, 1
            ind_mm = ind(utc_mm .eq. (mn + 1))
            if (.not.all(ismissing(ind_mm)))
                if (dimsizes(ind_mm) .eq. 1)
                    datam(mn, :, :) = data(ind_mm, :, :)
                else
                    datam(mn, :, :) = dim_avg_n_Wrap(data(ind_mm, :, :), 0)
                end if
            end if
            delete(ind_mm)
        end do
        datas(0, :, :) = (datam(0, :, :) + datam(1, :, :) + datam(11, :, :)) / 3.0
        datas(1, :, :) = dim_avg_n_Wrap(datam(2:4, :, :), 0)
        datas(2, :, :) = dim_avg_n_Wrap(datam(5:7, :, :), 0)
        datas(3, :, :) = dim_avg_n_Wrap(datam(8:10, :, :), 0)
        dataa(0, :, :) = dim_avg_n_Wrap(datam, 0)
    end if
    if (rank .eq. 4)
        do mn = 0, 11, 1
            ind_mm = ind(utc_mm .eq. (mn + 1))
            if (.not.all(ismissing(ind_mm)))
                if (dimsizes(ind_mm) .eq. 1)
                    datam(mn, :, :, :) = data(ind_mm, :, :, :)
                else
                    datam(mn, :, :, :) = dim_avg_n_Wrap(data(ind_mm, :, :, :), 0)
                end if
            end if
            delete(ind_mm)
        end do
        datas(0, :, :, :) = (datam(0, :, :, :) + datam(1, :, :, :) + datam(11, :, :, :)) / 3.0
        datas(1, :, :, :) = dim_avg_n_Wrap(datam(2:4, :, :, :), 0)
        datas(2, :, :, :) = dim_avg_n_Wrap(datam(5:7, :, :, :), 0)
        datas(3, :, :, :) = dim_avg_n_Wrap(datam(8:10, :, :, :), 0)
        dataa(0, :, :, :) = dim_avg_n_Wrap(datam, 0)
    end if

    datam@cell_methods = "time: mean within years time: mean over years"
    datas@cell_methods = "time: mean within years time: mean over years"
    dataa@cell_methods = "time: mean within years time: mean over years"

    datam&time@calendar = "climatology_bounds"
    datas&time@calendar = "climatology_bounds"
    dataa&time@calendar = "climatology_bounds"

    if (all(ismissing(atts))) then
      info_output("climate.ncl: no global attributes in "+infilename, verbosity, 3)
      fatt_mfile@source = infilename
      fatt_sfile@source = infilename
      fatt_afile@source = infilename
    else
      do  k = 0, dimsizes(atts) - 1, 1
          fatt_mfile@$atts(k)$ = fin@$atts(k)$
          fatt_sfile@$atts(k)$ = fin@$atts(k)$
          fatt_afile@$atts(k)$ = fin@$atts(k)$
      end do
    end if
    fileattdef(fout_mfile, fatt_mfile)
    fileattdef(fout_sfile, fatt_sfile)
    fileattdef(fout_afile, fatt_afile)

    fout_mfile->$base_variable$ = datam
    fout_sfile->$base_variable$ = datas
    fout_afile->$base_variable$ = dataa

    if (sum(where(var .eq. "hyam", 1, 0)) .gt. 0) then
        fout_mfile->hyam = fin->hyam
        fout_sfile->hyam = fin->hyam
        fout_afile->hyam = fin->hyam
    end if
    if (sum(where(var .eq. "hybm", 1, 0)) .gt. 0) then
        fout_mfile->hybm = fin->hybm
        fout_sfile->hybm = fin->hybm
        fout_afile->hybm = fin->hybm
    end if
    if (sum(where(var .eq. "hyai", 1, 0)) .gt. 0) then
        fout_mfile->hyai = fin->hyai
        fout_afile->hyai = fin->hyai
        fout_sfile->hyai = fin->hyai
    end if
    if (sum(where(var .eq. "hybi", 1, 0)) .gt. 0) then
        fout_mfile->hybi = fin->hybi
        fout_sfile->hybi = fin->hybi
        fout_afile->hybi = fin->hybi
    end if
    if (sum(where(var .eq. "p0", 1, 0)) .gt. 0) then
        fout_mfile->p0 = fin->p0
        fout_sfile->p0 = fin->p0
        fout_afile->p0 = fin->p0
    end if

    delete(data)
    delete(datam)
    delete(datas)
    delete(dataa)
    delete(utc_mm)
    info_output(">>>>>>>> Leaving climate.ncl", verbosity, 2)
end
