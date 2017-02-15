;;#############################################################################
;; SPECIFIC ROUTINES FOR DEFINING DIAGNOSTIC DESCRIPTIONS
;;#############################################################################
;; Please consider using or extending existing routines before adding new ones.
;; Check the header of each routine for documentation.      
;;
;; Contents:
;;    function SAMonsoon_daily_precip_propagation_descr
;;    function SAMonsoon_daily_precip_filter_descr
;;    function SAMonsoon_daily_precip_filter_norm_descr
;;    function SAMonsoon_teleconnections_descr
;;    function description_for_evapotranspiration
;;
;;#############################################################################

undef("SAMonsoon_daily_precip_propagation_descr")
function SAMonsoon_daily_precip_propagation_descr(di [1] : logical)
;;
;; Arguments
;;    di: logcial with dynamic diagnostic script information
;;
;; Return value
;;    description: a long text string (with line breaks) describing the
;;                 diagnostic at hand
;;
;; Caveats
;;
;; References
;;
;; Modificiation history
;;    20150702-A_eval_ma: written.
;;
local description
begin
    ;; ~C~ is newline
    description = "~F22~Diag def:~F21~ " + "Lagged correlations of daily data band pass filtered "+\
                  "~C~                 over " + di@min_days + "-" + di@max_days + " days. Implementation is similar to the NCL"+\
                  "~C~                 MJO Clivar script for lagged correlations with a few"+\
                  "~C~                 changes such as using smoothed precip anomalies"+\
                  "~C~                 and changing the season to JJAS. See the original"+\
                  "~C~                 script mjoclivar_9.ncl at, "+\
                  "~C~                        https://www.ncl.ucar.edu/Applications/mjoclivar.shtml."+\
                  "~C~"+\
                  "~C~                 The reference time series is computed from the"+\
                  "~C~                 weighted area average over"+\
                  "~C~                        lat: " + write_tuple_with_dash(di@prop_ref_latrange) + ", lon: " + write_tuple_with_dash(di@prop_ref_lonrange)+\
                  "~C~"+\
                  "~C~                 Latrange for the (lag, lon) plot is: " + write_tuple_with_dash(di@prop_lag_latrange) +\
                  "~C~                 Lonrange for the (lag, lat) plot is: " + write_tuple_with_dash(di@prop_ref_lonrange) 
    return description
end

undef("SAMonsoon_daily_precip_filter_descr")
function SAMonsoon_daily_precip_filter_descr(di [1] : logical)
local description
begin
    ;; ~C~ is newline
    description = "~F22~Diag def:~F21~ " + di@season + " average of standard deviation across"+\
                  "~C~                  indicated years of filtered daily rainfall. Filtering" +\
                  "~C~                  done by the NCL function band_pass_latlon(...). Data is" +\
                  "~C~                  normalized with climatoloy and areas with climatology precip"+\
                  "~C~                  less than 1 mm/day are masked away" +\
                  "~C~~F22~mean def:~F21~ NCL function wgt_areaave(...) computed with "+\
                  "~C~                  cosine weights in latitude."+\
                  "~C~~F22~rmse def:~F21~ NCL function wgt_arearmse(...) computed with "+\
                  "~C~                  cosine weights in latitude."+\
                  "~C~~F22~corr def:~F21~ NCL function pattern_cor(...) computed with "+\
                  "~C~                  cosine weights in latitude. Correlation is between"+\
                  "~C~                  reference and model data."+\
                  "~C~~F22~bob def:~F21~ Area avearge over Bay of Bengal normalized with climatology, "+\
                  "~C~              lat=" + di@bob_latrange(0) + "-" + di@bob_latrange(1)+\
                  ", lon=" + di@bob_lonrange(0) + "-" + di@bob_lonrange(1)+\
                  "~C~~F22~eio def:~F21~ Area avearge over eastern equatorial Indian Ocean normalized"+\
                  "~C~              with climatology, lat=" + di@eio_latrange(0) + "-" + di@eio_latrange(1)+\
                  ", lon=" + di@eio_lonrange(0) + "-" + di@eio_lonrange(1)
    return description
end

undef("SAMonsoon_daily_precip_filter_norm_descr")
function SAMonsoon_daily_precip_filter_norm_descr(di [1] : logical)
local description
begin
    ;; ~C~ is newline
    description = "~F22~Diag def:~F21~ " + di@season + " average of standard deviation across"+\
                  "~C~                  indicated years of filtered daily rainfall. Filtering" +\
                  "~C~                  done by the NCL function band_pass_latlon(...). Data is" +\
                  "~C~                  normalized with climatoloy and areas with climatology precip"+\
                  "~C~                  less than 1 mm/day are masked away" +\
                  "~C~~F22~mean def:~F21~ NCL function wgt_areaave(...) computed with "+\
                  "~C~                  cosine weights in latitude."+\
                  "~C~~F22~rmse def:~F21~ NCL function wgt_arearmse(...) computed with "+\
                  "~C~                  cosine weights in latitude."+\
                  "~C~~F22~corr def:~F21~ NCL function pattern_cor(...) computed with "+\
                  "~C~                  cosine weights in latitude. Correlation is between"+\
                  "~C~                  reference and model data."+\
                  "~C~~F22~bob def:~F21~ Area avearge over Bay of Bengal normalized with climatology, "+\
                  "~C~              lat=" + di@bob_latrange(0) + "-" + di@bob_latrange(1)+\
                  ", lon=" + di@bob_lonrange(0) + "-" + di@bob_lonrange(1)+\
                  "~C~~F22~eio def:~F21~ Area avearge over eastern equatorial Indian Ocean normalized"+\
                  "~C~              with climatology, lat=" + di@eio_latrange(0) + "-" + di@eio_latrange(1)+\
                  ", lon=" + di@eio_lonrange(0) + "-" + di@eio_lonrange(1)
    return description
end

undef("SAMonsoon_teleconnections_descr")
function SAMonsoon_teleconnections_descr(di [1] : logical)
local description
begin
    ;; ~C~ is newline
    description = "~F22~Diag def:~F21~ Map of correlations between Monsoon precip and area averaged SST "+\
                  "~C~             timeseries. Computed with the NCL function 'escorc(...)'"+\
                  "~C~~F22~mean def:~F21~ NCL function wgt_areaave(...) computed with "+\
                  "~C~                  cosine weights in latitude."+\
                  "~C~~F22~rmse def:~F21~ NCL function wgt_arearmse(...) computed with "+\
                  "~C~                  cosine weights in latitude."+\
                  "~C~~F22~corr def:~F21~ NCL function pattern_cor(...) computed with "+\
                  "~C~                  cosine weights in latitude. Correlation is between"+\
                  "~C~                  reference and model data."+\
                  "~C~~F22~sa34 def:~F21~ NCL function pattern_cor(...) computed with "+\
                  "~C~                       cosine weights in latitude. Correlation is between"+\
                  "~C~                       precip in lat: " + di@monsoon_latrange(0) + "-" + di@monsoon_latrange(1) +\
                  ", lon: " + di@monsoon_lonrange(0) + "-" + di@monsoon_lonrange(1) + " and " + " SST in Nino3.4" +\
                  "~C~~F22~en34 def:~F21~ NCL function pattern_cor(...) computed with "+\
                  "~C~                       cosine weights in latitude. Correlation is between"+\
                  "~C~                       precip and SST in Nino3.4" +\
                  "~C~~F22~Diag ref:~F21~ http://link.springer.com/article/10.1007%2Fs00382-012-1607-6"
    return description
end

undef("description_for_evapotranspiration")
function description_for_evapotranspiration(di [1] : logical)
local description
begin
    ;; ~C~ is newline
    description = "~F22~Diag def:~F21~ Evapotranspiration in units 'mm day-1'"+\
                  "~C~                     et = hfls * 24 * 3600 / et_lhv"+\
                  "~C~                 where,"+\
                  "~C~                            hfls = 'Surface Upward Latent Heat Flux' [W m-2]"+\
                  "~C~                          et_lhv = 2.465e6 [J kg -1]"+\
                  "~C~            (*_lhv = latent heat of vaporisation)"+\
                  "~C~~F22~Diag ref:~F21~ http://onlinelibrary.wiley.com/doi/10.1002/2013GL058055/abstract)"
    return description
end
