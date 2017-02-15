;;#############################################################################
;; AUXILIARY ROUTINES FOR AEROSOL DIAGNOSTICS
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    function obs_model_timeseries
;;
;;#############################################################################

load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

;;#############################################################################
undef("obs_model_timeseries")
function obs_model_timeseries(ymin_mod[1]:integer,
                              ymax_mod[1]:integer,
                              data[1]:logical,
                              variable[1]:string)
;;
;; Arguments
;;    ymin_mod: models start date.
;;    ymax_ref: models end data.
;;    data: model data.
;;    variable: variable to read
;;
;; Return value
;;    A 3-D array with colocated model and observations values for the given
;;    variable for the given time range.
;;
;; Description
;;    Creates time-series of co-located model and observations for all
;;    stations.
;;
;; Caveats
;;
;; Reference
;;
;; Modification history
;;    20140331-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, timeunits, stations, nstations, \
    stname, stlat, stlon, stalt, ymin_ref, ymax_ref, ntime_ref, time_ref, yy, \
    mm, date_ref, year_ref, mnames, stID, f, xx, date_obs, year_obs, id1_ref, \
    id1_obs, id2_ref, id2_obs, id1_mod, id2_mod, date_mod, year_mod, id_lon, \
    id_lat
begin

    funcname = "obs_model_timeseries"
    scriptname = "diag_scripts/aux/aerosol/aux_aerosol.ncl"
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Expected time units
    timeunits = "days since 1950-01-01 00:00:00"

    ;; Read list of stations
    stations = systemfunc("ls " + diag_script_info@datapath + \
                          "/OBS_" + diag_script_info@network + \
                          "_ground_*_T0M_" + variable + "_*.nc")

    nstations = dimsizes(stations)
    stname = new(nstations, string)
    stlat  = new(nstations, double)
    stlon  = new(nstations, double)
    stalt  = new(nstations, float)
    info_output(nstations + " stations found", verbosity, 2)

    ;; Find overlap between model range and observations range
    f = addfile(stations(0), "r")
    date = cd_calendar(f->time, 0)
    ymin_obs = min(date(:, 0))
    ymax_obs = max(date(:, 0))
    ymin_ref = toint(max((/ymin_mod, ymin_obs/)))
    ymax_ref = toint(min((/ymax_mod, ymax_obs/)))
    delete(date)

    ;; Return missing if no overlap
    if (ymin_ref.gt.ymax_ref) then
        data_array = new(1, float)
        return(data_array)
    end if

    ;; Define a reference time coordinate for all data
    ntime_ref = (ymax_ref - ymin_ref + 1) * 12
    time_ref = new(ntime_ref, double)
    time_ref@units = timeunits
    do yy = ymin_ref, ymax_ref
        do mm = 1, 12
            time_ref(12 * (yy - ymin_ref) + mm - 1) \
                = cd_inv_calendar(yy, mm, 15, 0, 0, 0, timeunits, 0)
        end do
    end do
    date_ref = cd_calendar(time_ref, 0)
    year_ref = date_ref(:, 0)

    ;; Create data arrays
    data_array = new((/1 + dimsizes(models@name), ntime_ref, nstations/), \
                     float)
    data_array!0 = "model"
    data_array!1 = "time"
    data_array!2 = "station"
    mnames = array_append_record(diag_script_info@network, \
                                 unique_labels_min((/"name", "case_name", \
                                 "ensemble", "experiment"/)), 0)
    data_array&model = mnames
    data_array&time = time_ref

    ;; Loop over stations, extract data
    do stID = 0, nstations - 1

        f = addfile(stations(stID), "r")
        stname(stID) = f@station
        stlat(stID)  = todouble(f@latitude)
        stlon(stID)  = todouble(f@longitude)
        stlon = where(stlon.lt.0., stlon + 360., stlon)
        if (isatt(f, "altitude"))
            stalt(stID) = tofloat(f@altitude)
        end if

        xx = f->$variable$
        date_obs = cd_calendar(xx&time, 0)
        year_obs = date_obs(:, 0)
        delete(date_obs)

        ;; Skip station if no overlap in time interval
        if (max(year_obs) .lt. min(year_ref) .or. \
            min(year_obs) .gt. max(year_ref)) then
            info_output("No data for station " + stname(stID) \
                        + " in the time range " \
                        + min(year_ref) + "-" + max(year_ref), verbosity, 2)
            continue
        end if

        ;; Extract time interval matching reference time coordinate
        if (min(year_obs) .lt. min(year_ref)) then
            id1_ref = 0
            id1_obs = min(ind(year_obs .ge. min(year_ref)))
        else
            id1_ref = min(ind(year_ref .ge. min(year_obs)))
            if (ismissing(id1_ref)) then
                data_array = data_array@_FillValue
            end if
            id1_obs = 0
        end if
        if (year_ref(id1_ref) .ne. year_obs(id1_obs)) then
            error_msg("f", "aux_aerosol", funcname, "time coordinate error")
        end if

        if (max(year_obs) .gt. max(year_ref)) then
            id2_ref = dimsizes(year_ref) - 1
            id2_obs = max(ind(year_obs.le. max(year_ref)))
        else
            id2_ref = max(ind(year_ref .le. max(year_obs)))
            id2_obs = dimsizes(year_obs) - 1
        end if

        if (year_ref(id2_ref) .ne. year_obs(id2_obs)) then
            error_msg("f", "aux_aerosol", funcname, "time coordinate error")
        end if

        ;; Assign to observations array
        data_array(0, id1_ref:id2_ref, stID) = xx(id1_obs:id2_obs)

        delete(id1_ref)
        delete(id2_ref)
        delete(id1_obs)
        delete(id2_obs)
        delete(year_obs)
        delete(xx)

        ;; Loop over models, extract data
        do mID = 0, dimsizes(models@name) - 1

            xx = extract_data(mID, data, - 1, 0, 0)

            if (mID.eq.0) then
                copy_VarAtts(xx, data_array)
            end if

            date_mod = cd_calendar(xx&time, 0)
            year_mod = date_mod(:, 0)
            delete(date_mod)

            ;; Extract  time interval matching reference time coordinate
            if (min(year_mod) .lt. min(year_ref)) then
                id1_ref = 0
                id1_mod = min(ind(year_mod .ge. min(year_ref)))
            else
                id1_ref = min(ind(year_ref .ge. min(year_mod)))
                id1_mod = 0
            end if
            if (year_ref(id1_ref) .ne. year_mod(id1_mod)) then
                error_msg("f", "aux_aerosol", funcname, \
                          "time coordinate error")
            end if
            if (max(year_mod) .gt. max(year_ref)) then
                id2_ref = dimsizes(year_ref) - 1
                id2_mod = max(ind(year_mod .le. max(year_ref)))
            else
                id2_ref = max(ind(year_ref .le. max(year_mod)))
                id2_mod = dimsizes(year_mod) - 1
            end if
            if (year_ref(id2_ref) .ne. year_mod(id2_mod)) then
                error_msg("f", "aux_aerosol", funcname, \
                          "time coordinate error")
            end if
            delete(year_mod)

            ;; Extract longitude
            id_lon = stlon(stID)
            if (all(xx&lon .lt. stlon(stID))) then
                if ((abs(360. + min(xx&lon) - stlon(stID))) .lt. \
                    abs(stlon(stID) - max(xx&lon))) then
                        id_lon = min(xx&lon)
                else
                    id_lon = max(xx&lon)
                end if
            end if
            if (all(xx&lon .gt. stlon(stID))) then
                if (abs(min(xx&lon) - stlon(stID)) .lt. \
                    abs(360. - max(xx&lon) + stlon(stID))) then
                    id_lon = min(xx&lon)
                else
                    id_lon = max(xx&lon)
                end if
            end if

            ;; Extract latitude
            id_lat = stlat(stID)
            if (all(xx&lat .lt. stlat(stID)) .or. \
                all(xx&lat .gt. stlat(stID))) then
                if (all(xx&lat .lt. stlat(stID))) then
                    id_lat = max(xx&lat)
                else
                    id_lat = min(xx&lat)
                end if
            end if

            ;; Assign to model array
            data_array(mID + 1, id1_ref:id2_ref, stID) = \
                (/xx(id1_mod:id2_mod, {id_lat}, {id_lon})/)
            delete(id1_ref)
            delete(id2_ref)
            delete(id1_mod)
            delete(id2_mod)
            delete(id_lat)
            delete(id_lon)
            delete(xx)

        end do  ; models
    end do  ; stations

    ;; Append station info
    data_array&station = stname + " ; " + stlat + " ; "\
                                + stlon + " ; " + stalt

    leave_msg(scriptname, funcname, 4)
    return(data_array)

end
