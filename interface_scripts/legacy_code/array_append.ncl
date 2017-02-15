
;
;  Append arrays of different rankings
;
;  2008/11/26  CAF  Created
;
undef("array_start_record")
function array_start_record(var, index)
local dim, rank, result
begin
    dim = dimsizes(var)
    rank = dimsizes(dim)
    if (rank .eq. 1) then
        result = var(index:)
    end if
    if (rank .eq. 2) then
        result = var(index:, :)
    end if
    if (rank .eq. 3) then
        result = var(index:, :, :)
    end if
    if (rank .eq. 4) then
        result = var(index:, :, :, :)
    end if
    return(result)
end

undef("array_mid_record")
function array_mid_record(var, index1, index2)
local dim, rank, result
begin
    dim = dimsizes(var)
    rank = dimsizes(dim)
    if (rank .eq. 1) then
        result = var(index1:index2)
    end if
    if (rank .eq. 2) then
        result = var(index1:index2, :)
    end if
    if (rank .eq. 3) then
        result = var(index1:index2, :, :)
    end if
    if (rank .eq. 4) then
        result = var(index1:index2, :, :, :)
    end if
    return(result)
end

undef("array_end_record")
function array_end_record(all_var, var, index)
local dim, rank, result
begin
    dim = dimsizes(var)
    rank = dimsizes(dim)
    if (rank .eq. 1) then
        result = array_append_record(all_var, var(:index), 0)
    end if
    if (rank .eq. 2) then
        result = array_append_record(all_var, var(:index, :), 0)
    end if
    if (rank .eq. 3) then
        result = array_append_record(all_var, var(:index, :, :), 0)
    end if
    if (rank .eq. 4) then
        result = array_append_record(all_var, var(:index, :, :, :), 0)
    end if
    return(result)
end
