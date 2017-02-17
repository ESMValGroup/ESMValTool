undef("find_name")
function find_name(var:string)
;;
;; Arguments
;;    var: input variable.
;;
;; Return value
;;    The GO name of the selected variable.
;;
;; Description
;;    Finds the corresponding GO name of the input variable var.
;;
;; Caveats
;;    Needs an external input table "names_GO.dat"
;;
;; Reference
;;
;; Modification history; developed from 
;;    20130502-A_RiMa: written refoprmat_EMAC_func.ncl
;;
;;
local funcname, fname, data, std, emac, type, vert, id
begin

    funcname = "find_name"

    fname = "./reformat_scripts/GO/names_GO.dat"

    data = readAsciiTable(fname, 1, "string", \
                          dimsizes(readAsciiHead(fname, "STANDARD")) + 1)

    std  = str_squeeze(str_get_field(data(:, 0), 1, "|"))
    go_var = str_squeeze(str_get_field(data(:, 0), 2, "|"))
;;    type = str_squeeze(str_get_field(data(:, 0), 3, "|"))
;;    vert = str_squeeze(str_get_field(data(:, 0), 4, "|"))

    id = ind(std.eq.var)
   


    if (ismissing(id)) then
        error_msg("f", "reformat_GO_func.ncl", funcname, \
                  "cannot find GO name for " + var)
    else
        out = go_var(id)
;;        out@type = type(id)
;;        out@vert = vert(id)
    end if
    return(out)

end