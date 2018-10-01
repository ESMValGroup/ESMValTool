#####################################################################
#
#  HyInt configuration file
#
# About: configuration file for ESMValTool HyInt namelist.
#        General configuration for running the HyInt diagnostic 
#        (models and time period are defined in the HyInt namelist).
#        In order to optimize recursive analysis, users can spearately select: 
#        a) grid and region of HyInt pre-processing and dignostic calculation 
#        b) region and years to be plotted
#
#####################################################################

# Replace here global values from namelist_hyint.xml
# write_netcdf=T
# write_plots=T
run_regridding=T
force_regridding=F
run_diagnostic=T
force_diagnostic=F
etccdi_preproc=F
run_timeseries=T
# output_file_type="png"
diag_script_cfg="./esmvaltool/diag_scripts/aux/hyint/hyint_parameters.r"

# Define time dimension name
timedimname="time"

# Pre-processing options
seasons <- c("ALL")   # seasons to be analysed: "ALL", "DJF", ...

rgrid   <- F #"REF" #"r320x160"  # a) FALSE to keep original resolution
                          # b) set desired regridding resolution in cdo format e.g., "r320x160"
                          # c) "REF" to use resolution of reference model

rotLongitude <- "full" # a) "full" to convert input arrays from 0/360 to -180/180 longitude grid
                     # b) "no" to leave input data on its original grid

rlonlatdata <- F #c(-180,180,-90,90) # region where the pre-processing and diagnostic calculation is to be performed. NOTE: lon(-180/180) 
                           # (keep global coverage to allow re-use of same data set for plots over different regions)
                           # Set FALSE to keep geographical coverage of the original dataset
grid_file <- "grid_file"  # suffix for grid file
topography_file <- "topo" # suffix for topography file (needed for filtering ocean/land or land elevation)
topography_highres <- "/home/arnone/work/data/Elevation/GMTED2010_15n030_0125deg.nc"
etccdi_dir <- "/work/users/arnone/data/ETCCDI/historical+rcp85/"


# Diagnostic options
#norm_years=c(1976,2005)  # reference normalization period - set in namelist
external_norm=F#"HIST"       # a) F=use internal data to normalize over the norm_years period
                      # b) list of names of normalization files (one per input data file or one for all)
                      # c) "HIST" to automatically generate the name of the historical experiment associated with the model name 
#c("/home/arnone/work/esmtest/work/HyInt/EC-Earth/1976_2005/ALL/HyInt_EC-Earth_historical_r8i1p1_r320x160_1976_2005_ALL.nc")

external_r95=external_norm    # a) F=use internal data to evaluate r95 threshold over the norm_years period  
                              # b) list of names of files (one per input data file or one for all) 
                              # c) "HIST" to automatically generate the name of the historical experiment associated with the model name

# Plotting options
# Plot_type set in namelist
# plot_type <- 14  # 1) lon/lat maps per individual field/exp/single year, 2) lon/lat maps per individual field exp-ref-diff/single year, 
                # 3) lon/lat maps multi-field/exp-ref-diff/single year,  4) lon/lat maps multifield/exp/multiyear,  
                # 11) timeseries over required individual region/exp, 12) timeseries over multiple regions/exp
                # 13) timeseries with multiple models, 14) summary trend coefficients multiple regions  
                #  15) summary trend coefficients multiple models 
npancol=2 # number of columns for trend/tseries multi-panel figures
npanrow=3 # number of rows for trend/tseries multi-panel figures
ryearplot <- 2006 # c(1997,2002,2003) # years to be plotted for experiments (maps over individual years): 
                  # a) actual years, b) "FIRST" = first year in dataset or c) "ALL"  = all years in dataset. E.g., c(1998,2000,2005)   
rmultiyear_mean <- T # T to plot multiyear mean (this override ryearplot)
ryearplot_ref <- c("EXP") # year to be plotted for reference dataset: options a) "EXP" == same as experiments, b) one year only, e.g. c(1998)    
force_ref <- F # set TRUE to force plotting of reference data as any other experiment

# user defined extra label for figure file name
#label= "test" set in namelist 

map_continents <- -2 # thickness of continents: positive values in white, negative values in gray
map_continents_regions <- F # T to plot also regional boundaries

# colorbar
add_colorbar=F # T to add colorbar
legend_distance=3

# timeseries options
weight_tseries=T  # T to calculate area weighted time averages
trend_years= F#c(2006,2100,1976,2005) # a) F=all; 
                         # b) c(year1,year2) to apply trend calculation and plotting only to a limited time interval (year1<=years<=year2) 
                         # c) c(year1,year2,year3,year4) to apply trend to two separate time intervals (year1<=years<=year2) and (year3<=years<=year4)
removedesert=F      # T to remove (flag as NA) grid points with mean annual pr < 0.5 mm/day (desertic areas, Giorgi et al. 2014)
maskSeaLand=F # T to mask depending on seaLandElevation threshold
seaLandElevation=0 # a) 0 land; b) positive value: land above given elevation;
                   # c) negative value: ocean below given depth. The topography/bathymetry file is generated with cdo from ETOPO data. 
reverse_maskSeaLand=F # T to reject what selected, F to keep what selected
highreselevation=F#500 # a) F: neglect; b) value: threshold of minimum elevation to be overplotted with contour lines of elevation
highreselevation_only=F # T to plot only high resolution elevation contours
oplot_grid=F # T to plot grid points over maps

# timeseries and trend plotting options
lm_trend=T         # T to calculate linear trend
add_trend=T        # T to add linear trend to plot
add_trend_sd=F     # T to add stdev range to timeseries
add_trend_sd_shade=F   # T to add shade of stdev range to timeseries 
add_tseries_lines=F    # T to plot lines of timeseries over points
add_zeroline=T         # T to plot a dashed line at y=0 
trend_years_only=F # T to limit timeseries plotting to trend_years[1:2] time interval
scale100years=T    # T to plot trends as 1/100 years
scalepercent=F     # T to plot trends as percent change (this is not applied to HY-INT)
add_legend=5       # a) F=no legend; b) n>0 list disposed in n column; c) <0 horizontal legend 
xy_legend=c(0.03,0.4) # position of legend in fraction of plotting panel 
tag_legend=c(T,F,F) # 1=model name, 2=model experiment, 3=model ensemble (select one or more)

# region box matrix (predefined following Giorgi et al. 2011,2014): add here further regions and select those needed through iregion
region_names=c("World","World60","Tropics","South-America","Africa","North-America","India","Europe","East-Asia","Australia")
region_codes=c("Globe","GL","TR","SA","AF","NA","IN","EU","EA","AU")
# Select one or more index values to define regions to be used. Default c(1) == global. 
# selregions=c(1:10) # set in namelist 

boxregion=F  #-2 # !=0 plot region boxes over maps with thickness = abs(boxregion) and white (>0) or grey (<0).  This automatically works on global maps only.

regions=matrix(nrow=length(region_names),ncol=4)
# c(lon1,lon2,lat1,lat2) NOTE: lon(-180/180)
regions[1,]=c(-180,180,-90,90) # First row = global
regions[2,]=c(-180,180,-60,60) # GL |lat|<60
regions[3,]=c(-180,180,-30,30)
regions[4,]=c(-90,-30,-60,10)
regions[5,]=c(-20,60,-40,35)
regions[6,]=c(-140,-60,10,60)
regions[7,]=c(60,100,0,35)
regions[8,]=c(-10,30,35,70)
regions[9,]=c(100,150,20,50)
regions[10,]=c(110,160,-40,-10)

# define fields for timeseries calculation and plotting
hyint_list=c("int_norm","dsl_norm","wsl_norm","hyint","int","dsl","wsl","pa_norm","r95_norm")
etccdi_yr_list=c("altcddETCCDI","altcsdiETCCDI","altcwdETCCDI","altwsdiETCCDI","cddETCCDI",
              "csdiETCCDI","cwdETCCDI","dtrETCCDI","fdETCCDI","gslETCCDI","idETCCDI","prcptotETCCDI",
              "r10mmETCCDI","r1mmETCCDI","r20mmETCCDI","r95pETCCDI","r99pETCCDI","rx1dayETCCDI",
              "rx5dayETCCDI","sdiiETCCDI","suETCCDI","tn10pETCCDI","tn90pETCCDI","tnnETCCDI",
              "tnxETCCDI","trETCCDI","tx10pETCCDI","tx90pETCCDI","txnETCCDI","txxETCCDI","wsdiETCCDI")

etccdi_list_import=etccdi_yr_list
field_names=c(hyint_list,etccdi_yr_list)

# Select one or more fields to be plotted with the required order  
# selfields=c(8,4,1,9,3,2) # set in the namelist

# define titles and units
title_unit_m=matrix(nrow=length(field_names),ncol=4)
title_unit_m[1,]=c("SDII","Norm. annual mean INT","Norm. annual mean precipitation intensity","")
title_unit_m[2,]=c("DSL","Norm. annual mean DSL","Norm. annual mean dry spell length","")
title_unit_m[3,]=c("WSL","Norm. annual mean WSL","Norm. annual mean wet spell length","")
title_unit_m[4,]=c("HY-INT","HY-INT","Hydroclimatic intensity","")
title_unit_m[5,]=c("ABS_INT","Annual mean INT","Annual mean precipitation intensity","(mm/day)")
title_unit_m[6,]=c("ABS_DSL","Annual mean DSL","Annual mean dry spell length","(days)")
title_unit_m[7,]=c("ABS_WSL","Annual mean WSL","Annual mean wet spell length","(days)")
title_unit_m[8,]=c("PA","Normalized precipitation area","Norm. precipitation area","")
title_unit_m[9,]=c("R95","Norm. heavy precipitation index","Norm. % of total precip. above 95% percentile of reference distribution","")

# define levels for contour/yrange for abs. values: (minlev,maxlev,minlev_diff,maxlev_diff) and nlev
nlev=24
levels_m=matrix(nrow=length(field_names),ncol=4)

levels_m[1,]=c(0.9,1.1,-1.2,1.2)
levels_m[1,]=c(0.5,1.3,-1.2,1.2)
levels_m[2,]=c(0.9,1.1,-1.2,1.2)
levels_m[2,]=c(0.6,1.4,-1.2,1.2)
levels_m[3,]=c(0.9,1.1,-1.2,1.2)
levels_m[3,]=c(0.7,1.3,-1.2,1.2)
levels_m[4,]=c(0.5,1.5,-1.2,1.2)
levels_m[5,]=c(0,10,-5,5)
levels_m[6,]=c(0,20,-5,5)
levels_m[7,]=c(0,10,-3,3)
levels_m[8,]=c(0.5,1.5,-1.2,1.2)
levels_m[9,]=c(0.5,1.5,-2,2)
levels_m[10,]=c(0,200,-5,5)
levels_m[11,]=c(-5,15,-5,5)
levels_m[12,]=c(0,20,-5,5)
levels_m[13,]=c(0,20,-5,5)
levels_m[14,]=c(0,200,-5,5)
levels_m[15,]=c(-10,30,-5,5)
levels_m[16,]=c(0,80,-5,5)
levels_m[17,]=c(0,15,-4,4)
levels_m[18,]=c(0,200,-10,10)
levels_m[19,]=c(0,400,-10,10)
levels_m[20,]=c(-10,200,-10,10)
levels_m[21,]=c(0,3000,-100,100)
levels_m[22,]=c(0,80,-10,10)
levels_m[23,]=c(0,300,-10,10)
levels_m[24,]=c(0,50,-2,2)
levels_m[25,]=c(0,800,-20,20)
levels_m[26,]=c(0,300,-10,10)
levels_m[27,]=c(0,100,-10,10)
levels_m[28,]=c(0,200,-10,10)
levels_m[29,]=c(0,15,-5,5)
levels_m[30,]=c(0,300,-20,20)
levels_m[31,]=c(-5,25,-5,5)
levels_m[32,]=c(0,300,-5,5)
levels_m[33,]=c(-40,40,-5,5)
levels_m[34,]=c(0,40,-5,5)
levels_m[35,]=c(-20,300,-5,5)
levels_m[36,]=c(-5,25,-2,2)
levels_m[37,]=c(-20,140,-4,4)
levels_m[38,]=c(-30,30,-5,5)
levels_m[39,]=c(0,50,-2,2)
levels_m[40,]=c(-20,320,-2,2)

# define levels for contour/yrange for trends (minlev,maxlev)
ntlev=24
tlevels_m=matrix(nrow=length(field_names),ncol=2)
tlevels_m[1,]=c(-0.05,0.2)*0.01
tlevels_m[2,]=c(-0.1,0.4)*0.01#c(-0.1,0.3)*0.01
tlevels_m[3,]=c(-0.1,0.1)*0.01
tlevels_m[4,]=c(0,0.4)*0.01#c(0,0.5)*0.01
tlevels_m[5,]=c(0,1.5)*0.01
tlevels_m[6,]=c(-1,6)*0.01
tlevels_m[7,]=c(-0.8,0.8)*0.01
tlevels_m[8,]=c(-0.3,0.5)*0.01
tlevels_m[9,]=c(0,0.6)*0.01#c(-0.2,0.6)*0.01
tlevels_m[10,]=c(0,200)*0.01
tlevels_m[11,]=c(0,12)*0.01
tlevels_m[12,]=c(0,20)*0.01
tlevels_m[13,]=c(0,20)*0.01
tlevels_m[14,]=c(0,15)*0.01 # c(0,30)*0.01
tlevels_m[15,]=c(-70,0)*0.01
tlevels_m[16,]=c(-4,4)*0.01#c(-5,5)*0.01
tlevels_m[17,]=c(-1,0)*0.01
tlevels_m[18,]=c(-70,10)*0.01
tlevels_m[19,]=c(-10,90)*0.01
tlevels_m[20,]=c(-60,0)*0.01
tlevels_m[21,]=c(-20,120)*0.01
tlevels_m[22,]=c(0,10)*0.01
tlevels_m[23,]=c(-15,5)*0.01
tlevels_m[24,]=c(0,6)*0.01
tlevels_m[25,]=c(0,100)*0.01
tlevels_m[26,]=c(0,60)*0.01
tlevels_m[27,]=c(0,15)*0.01#c(0,30)*0.01
tlevels_m[28,]=c(0,50)*0.01
tlevels_m[29,]=c(0,15)*0.01
tlevels_m[30,]=c(0,140)*0.01
tlevels_m[31,]=c(-30,0)*0.01
tlevels_m[32,]=c(0,100)*0.01
tlevels_m[33,]=c(0,8)*0.01
tlevels_m[34,]=c(0,8)*0.01
tlevels_m[35,]=c(0,150)*0.01
tlevels_m[36,]=c(-30,0)*0.01
tlevels_m[37,]=c(0,160)*0.01
tlevels_m[38,]=c(2,8)*0.01
tlevels_m[39,]=c(0,8)*0.01
tlevels_m[40,]=c(-100,300)*0.01

autolevels=F # undefine levels if you wish to plot actual range of target values
if (autolevels) {  
 tlevels_m[]=NA
 levels_m[]=NA
}

# Specific settings for PNG output
png_height=720
png_width=960
png_width=1440
png_units="px"
png_pointsize=12
png_bg="white"

# Specific settings for PDF and EPS output (in inches)
pdf_width=24
pdf_height=12

# Specific settings for x11 output (in inches)
x11_width=7
x11_height=8

# color palette to be used
# palette0 is taken from tim.colors of field to avoid library dependencies...
palette0=colorRampPalette(c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF",
        "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF",
        "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF",
        "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF",
        "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF",
        "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F",
        "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50",
        "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00",
        "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00",
        "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000",
        "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000",
        "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000",
        "#AF0000", "#9F0000", "#8F0000", "#800000"))
palette1=colorRampPalette(c("white","orange","darkred"))
palette2=colorRampPalette(c("blue","white","red"))
palette3=colorRampPalette(c("darkblue","blue","dodgerblue","white","orange","red","darkred"))
palette_giorgi2011=colorRampPalette(c("white","khaki1","darkseagreen2","mediumseagreen","lightskyblue1",                
                      "lightskyblue","deepskyblue2","dodgerblue2","dodgerblue3","royalblue4"))
