;*********************************************************************
; uniform_coordinates.ncl
;*********************************************************************
; PROCEDURE uniform_coordinates
; fields supported:
; "T2Ms", "T3M", "T2Mz", "T1Ms", "T2Ds", "T2Dz", "T2Is", "T2Iz", "T0As", \
; "TO2Ms", "C3M", "C2Mz", "C2Ms"
;  AG" add "T0M"
;
; FUNCTION convert_doubletofloat
;
;*********************************************************************
undef("uniform_coordinates")
procedure uniform_coordinates(data:numeric, FieldNo:string)
; Check "lat" is monotonically increasing, otherwise flip
; Check "plev" is monotonically decreasing, otherwise flip

local idx, plot_fields, ind_field, latX, plevX
begin

    funcname = "uniform_coordinates"

    plot_fields = (/"T2Ms", "T3M", "T2Mz", "T1Ms", "T2Ds", "T2Dz", \
                    "T2Is", "T2Iz", "T0As", "TO2Ms", "T1Ds", \
                    "C3M", "C2Mz", "C2Ms", "C2Ds", "C1Ds", "T0M"/)
    ind_field = ind(plot_fields .eq. FieldNo)
    if (all(ismissing(ind_field))) then
        error_msg("f", "uniform_coordinates.ncl", funcname, FieldNo + \
                  " not upported")
    end if
    ;; "T3M", "C3M"
    if ((FieldNo .eq. "T3M") .or. (FieldNo .eq. "C3M")) then
        latX = data&lat
        if (latX(0) .gt. latX(1)) then
            data = data(:, :, ::-1, :)
        end if
        plevX = data&plev
        if (plevX(0) .lt. plevX(1)) then
            data = data(:, ::-1, :, :)
        end if
    end if
    ;; "T2Mz", "T2Dz", "C2Mz", "T2Iz"
    if ((FieldNo .eq. "T2Mz") .or. (FieldNo .eq. "T2Dz") .or. \
        (FieldNo .eq. "C2Mz") .or. (FieldNo .eq. "T2Iz")) then
        latX = data&lat
        if (latX(0) .gt. latX(1)) then
            data = data(:, :, ::-1)
        end if
        plevX = data&plev
        if (plevX(0) .lt. plevX(1)) then
            data = data(:, ::-1, :)
        end if
    end if
    ;; "T2Ms", "T2Ds", "C2Ms", "T2Is", "TO2Ms", "C2Ds"
    if ((FieldNo .eq. "T2Ms") .or. (FieldNo .eq. "T2Ds") .or. \
        (FieldNo .eq. "C2Ms") .or. (FieldNo .eq. "T2Is") .or. \
        (FieldNo .eq. "TO2Ms") .or. (FieldNo .eq. "C2Ds")) then
        latX = data&lat
        if (latX(0) .gt. latX(1)) then
            data = data(:, ::-1, :)
        end if
    end if
    ;; "T1Ms", "T0As", "T1Ds", "C1Ds"
    if ((FieldNo .eq. "T1Ms") .or. (FieldNo .eq. "T0As") .or. \
        (FieldNo .eq. "T1Ds") .or. (FieldNo .eq. "C1Ds")) then
        latX = data&lat
        if (latX(0) .gt. latX(1)) then
            data = data(:, ::-1)
        end if
    end if

end

;*********************************************************************
undef("convert_doubletofloat")
function convert_doubletofloat(data:numeric)
local result, data_att, data_flt_att, iatt, vardims, data_dim, \
      data_flt_dim, idim, nvars, dim_type, idx_dbl
begin

    if ((typeof(data) .eq. "double") .or. (typeof(data) .eq. "float")) then
        if (typeof(data) .eq. "double") then
            result = doubletofloat(data)

            ;; VarAtts
            data_att = getvaratts(data)
            data_flt_att = getvaratts(data)
            if (.not. all(ismissing(data_att))) then
                if (all(ismissing(data_flt_att))) then
                    copy_VarAtts(data, result)
                else
                    ;; copy attributes if does not exist
                    do iatt = 0, dimsizes(data_att) - 1
                        if(.not. isatt(result, data_att(iatt))) then
                            result@$data_att(iatt)$ = data@$data_att(iatt)$
                        end if
                    end do
                end if
            end if
            delete(data_att)
            delete(data_flt_att)

            ;; VarDims
            vardims = getvardims(data)
            do idim = 0, dimsizes(vardims) - 1
                data_dim = data&$vardims(idim)$
                data_flt_dim = doubletofloat(data_dim)
                copy_VarAtts(data_dim, data_flt_dim)
                result!idim = vardims(idim)
                result&$vardims(idim)$ = data_flt_dim
                delete(data_dim)
                delete(data_flt_dim)
            end do
            delete(vardims)
        else
            ;; VarDims
            vardims = getvardims(data)
            nvars = dimsizes(vardims)
            dim_type = new(nvars, "string")
            do idim = 0, nvars - 1
                dim_type(idim) = typeof(data&$vardims(idim)$)
            end do
            idx_dbl = ind(dim_type .eq. "double")
            if (all(ismissing(idx_dbl))) then
                ;; no conversion reqiured
                result = data
            else
                result = (/data/)
                copy_VarAtts(data, result)
                do idim = 0, dimsizes(vardims) - 1
                    if (dim_type(idim) .eq. "double") then
                        data_dim = data&$vardims(idim)$
                        data_flt_dim = doubletofloat(data_dim)
                        copy_VarAtts(data_dim, data_flt_dim)
                        delete(data_dim)
                    else
                        data_flt_dim = data&$vardims(idim)$
                    end if
                    result!idim = vardims(idim)
                    result&$vardims(idim)$ = data_flt_dim
                    delete(data_flt_dim)
                end do
            end if
            delete(vardims)
            delete(dim_type)
            delete(idx_dbl)
        end if
    else
        ;; no conversion reqiured
        result = data
    end if

    return(result)

end
