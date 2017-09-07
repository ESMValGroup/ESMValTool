;
;  Requires: clwvi:T2*s,clivi:T2*s
;

variable_info = True
variable_info@derived = True
variable_info@long_name = "liquid water path"
variable_info@units = "kg m-2"

; Contour levels (in g m-2)
; note that variable will be plotted in kg m-2, so contour level will be divided by 1000
variable_info@map_ref_Levels = ispan(10, 200, 10)
variable_info@map_diff_Levels = ispan(-50, 50, 10)

; misc
variable_info@title = "LWP"

undef("calculate")
function calculate(index [1] : integer,
                   variable [1] : string,
                   field_type [1] : string)
;;                 return_val [1] : logical
;; Arguments:
;;    index    - index to current infile defined ini the
;;               'interface_data/ncl.interface'-file
;;    variable - Current variable as string
;;    field_type  - string with field type classification
;; Return value:
;;    data_new - logical
local tmp, dum, i, verbosity
begin

;   print (">>>>>>>>> lwp.ncl <<<<<<<<<<")
;   print (models@name(index))
;   print (models@project(index))

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))

    project = models@project(index)
    model = models@name(index)

    info_output("<<<<<<<< Entering lwp.ncl", verbosity, 4)
    data_new = True
    tmp1 = read_data(index, "clwvi", "T2Ms")
    dum1 = extract_data(index, tmp1, -1, 0, 0)

    dum = dum1

    ; The variable "clwvi" (condensed water path) is supposed to contain the sum of liquid water ("lwp") and
    ; ice water (variable "clivi"). However, "clwvi" contains only liquid water (lwp) for some CMIP5 models.
    ; For these models, the ice water must not be subtracted from "clwvi" in order to calculate lwp.

    sub_iwp = 1 ; default: lwp = clwvi - clivi

    if ((project.eq."CMIP5").or.(project.eq."CMIP5_ETHZ")) then
       if (model.eq."CESM1-CAM5-1-FV2") then
          sub_iwp = 0
       end if
       if (model.eq."CESM1-CAM5") then
          sub_iwp = 0
       end if
       if (model.eq."CMCC-CESM") then
          sub_iwp = 0
       end if
       if (model.eq."CMCC-CM") then
          sub_iwp = 0
       end if
       if (model.eq."CMCC-CMS") then
          sub_iwp = 0
       end if
       if (model.eq."IPSL-CM5A-MR") then
          sub_iwp = 0
       end if
       if (model.eq."IPSL-CM5A-LR") then
          sub_iwp = 0
       end if
       if (model.eq."IPSL-CM5B-LR") then
          sub_iwp = 0
       end if
       if (model.eq."CCSM4") then
          sub_iwp = 0
       end if
       if (model.eq."IPSL-CM5A-MR") then
          sub_iwp = 0
       end if
       if (model.eq."MIROC-ESM") then
          sub_iwp = 0
       end if
       if (model.eq."MIROC-ESM-CHEM") then
          sub_iwp = 0
       end if
       if (model.eq."MIROC-ESM") then
          sub_iwp = 0
       end if
       if (model.eq."CSIRO-Mk3-6-0") then
          sub_iwp = 0
       end if
       if (model.eq."MPI-ESM-MR") then
          sub_iwp = 0
       end if
       if (model.eq."MPI-ESM-LR") then
          sub_iwp = 0
       end if
       if (model.eq."MPI-ESM-P") then
          sub_iwp = 0
       end if
    end if

    if (sub_iwp.eq.0) then
       print("lwp.ncl INFO:")
       print("   assuming that variable clwvi from " + project + " model " + model + " contains only liquid water")
    end if

    ; observational data (UWisc) contain LWP only

    if (project.eq."OBS") then
       if (model.eq."UWisc") then
          sub_iwp = 0
       end if
    end if

    if (sub_iwp.eq.1) then
       tmp2 = read_data(index, "clivi", "T2Ms")
       dum2 = extract_data(index, tmp2, -1, 0, 0)
       dum = dum1 - dum2  ; calculate lwp: subtract ice water path from total condensed water path (= water + ice)
    end if

    dum@long_name = variable_info@long_name
    dum@units = variable_info@units
    add_data_var(index, data_new, dum, variable)

    info_output(">>>>>>>> Leaving lwp.ncl", verbosity, 4)
    return(data_new)
end
