# This is a config file for CCI data and CMIP5 land cover diagnostics

# key list of data sets from CCI and CMIP, depracted
#could be written automatically? It's a copy of namelist_lc_ESACCI.xml with smaller changes
#translatorlist={ 'bare soil' : [['Bare_Soil'],['bare soil']],
#                 'natural grass' : [['Natural_Grass'],['grass']],
#                 'managed grass and crops' : [['Managed_Grass'],['crop','pasture']],
#                 'shrubs' : [['Shrub_Broadleaf_Evergreen','Shrub_Needleleaf_Evergreen','Shrub_Broadleaf_Deciduous','Shrub_Needleleaf_Deciduous'],['shrub']],
#                 'forest' : [['Tree_Broadleaf_Evergreen','Tree_Needleleaf_Evergreen','Tree_Broadleaf_Deciduous','Tree_Needleleaf_Deciduous'],['tree']]
#} 

#translatorlist={ 'bare soil' : [['Bare_Soil'],['bare soil']],
#                 'grass and crop' : [['Managed_Grass','Natural_Grass'],['crop','pasture','grass']],
#                 'shrub and forest' : [['Tree_Broadleaf_Evergreen','Tree_Needleleaf_Evergreen','Tree_Broadleaf_Deciduous','Tree_Needleleaf_Deciduous','Shrub_Broadleaf_Evergreen','Shrub_Needleleaf_Evergreen','Shrub_Broadleaf_Deciduous','Shrub_Needleleaf_Deciduous'],['tree','shrub']],
#}

translatorlist={ 'bare soil' : [['Bare_Soil'],['bare soil']],
                 'grass and crop' : [['Managed_Grass','Natural_Grass'],['crop','pasture','grass']],
                 'shrub and forest' : [['Tree_Broadleaf_Evergreen','Tree_Needleleaf_Evergreen','Tree_Broadleaf_Deciduous','Tree_Needleleaf_Deciduous','Shrub_Broadleaf_Evergreen','Shrub_Needleleaf_Evergreen','Shrub_Broadleaf_Deciduous','Shrub_Needleleaf_Deciduous'],['tree','shrub']],
}

# general flags for regionalization
regionalization = False
shape = "continents"
shapeNames = 2 #column of the name values
#start_year=2005
#stop_year=2005

# flags for basic diagnostics
globmeants = True
mima_globmeants=[0,100]
cmap_globmeants='YlGn'
mima_ts=[0,100]
mima_mts=[0,100]
portrait = False
globmeandiff = True
mima_globmeandiff=[-100,100]
mima_globmeandiff_r=[-1,1]
trend = False
trend_p = False

# flags for specific diagnostics
single_years=False #TODO rename variable
mima_single_year=[-100,100]
std_factor=4

#plotting parameters
projection={'projection': 'robin', 'lon_0': 180., 'lat_0': 0.}
