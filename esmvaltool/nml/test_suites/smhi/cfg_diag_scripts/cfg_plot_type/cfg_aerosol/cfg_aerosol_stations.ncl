; Settings used by aerosol_stations plot_types
plot_type_info = True

; Selected network (IMPROVE, EMEP, EANET, CASTNET, AERONET)
plot_type_info@network  = "IMPROVE"

; Path to station data
plot_type_info@datapath = "/nobackup/rossby16/sm_maeva/Data/EMAC/"

; How to deal with missing values in station data
;   "strong" = skip station if a missing value is present within the period
;   "weak"   = skip station if the series of valid values is non-continuous
;   "none"   = skip station only if all values are missing within the period
;   [value]  = skip station if number of missing values is > [value]%
plot_type_info@skipmode = "none"

; Plot mode (lin or log) for aerosol_stations_scatter
plot_type_info@scattermode = "log"

; Lat-lon limits for aerosol_stations_map
if (plot_type_info@network.eq."IMPROVE".or.\
    plot_type_info@network.eq."CASTNET") then
  plot_type_info@minlon = -170.
  plot_type_info@maxlon =  -50.
  plot_type_info@minlat =   10.
  plot_type_info@maxlat =   70.
end if

if (plot_type_info@network.eq."EMEP") then
  plot_type_info@minlon = -20
  plot_type_info@maxlon =  50.
  plot_type_info@minlat =  30.
  plot_type_info@maxlat =  80.
end if

if (plot_type_info@network.eq."EANET") then
  plot_type_info@minlon =   70.
  plot_type_info@maxlon =  160.
  plot_type_info@minlat =  -10.
  plot_type_info@maxlat =   60.
end if
