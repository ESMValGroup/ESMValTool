from netCDF4 import Dataset

file_ecv = "pr_day_EC-EARTH_historical_r12i1p1_19500101-19741231.nc"

ecv = Dataset(file_ecv, mode='r')
# print ecv

pr = ecv.variables['pr']
lat = ecv.variables['lat']
lon = ecv.variables['lon']
time = ecv.variables['time']

# print lat[:]

# if lat and lon are 1-dimensional: take all parwise combinations
if len(lat.shape)==1 and len(lon.shape)==1:
    coord_points = [(x,y) for x in lat[:] for y in lon[:]] 

# if lat and lon are matrices
if len(lat.shape)==2:
    # check if lat.shape == lon.shape, issue error otherwise
    coord_points = [(lat[x,y],lon[x,y]) for x in xrange(lat.shape[0]) for y in xrange(lon.shape[1])]

from shapely.geometry import MultiPoint
points = MultiPoint(coord_points)

import shapefile
sf = shapefile.Reader("./niger-hype_select.shp")
#first feature of the shapefile
# feature = sf.shapeRecords()[0]
# first = feature.shape.__geo_interface__  
# print first

shapes = sf.shapes()

# print sf.record(1)
# len(shapes)
# shapes[0].shapeType
# shapes[0].bbox
# fields = sf.fields[1:]
records = sf.records()
# len(records)
# len(records[1])
# shapes[1].bbox

fields = sf.fields[1:]
attr = [[records[i][j] for i in xrange(len(records))] for j in xrange(len(fields))]

from shapely.geometry import shape
# shp_geom = shape(first)
# print shp_geom
# LINESTRING (0 0, 25 10, 50 50)
# print type(shp_geom)

# shp_geom = shape(first)

from shapely.geometry.multipolygon import MultiPolygon
mulpol = MultiPolygon([shape(pol) for pol in shapes])
cent=mulpol.centroid

# from shapely.geometry import MultiPoint
cent = MultiPoint([pol.centroid for pol in mulpol])

from shapely.ops import nearest_points
selected_points=[]
for i in cent:
    nearest = nearest_points(i, points)
    selected_points.append(list(nearest[1].coords)[0])

# test=pr[:,selected_points[0][0],selected_points[0][1]]

# create table of results (time x one column for each catchment)
pr_out = []
for point in selected_points:
    pri = pr[:,point[0],point[1]]
    pr_out.append(pri)

# conversion to array for netcdf (nor sure it's needed?)
import numpy as np
test = np.array(pr_out)

# tries copy/creation netcdf
outnetcdf = Dataset("test.nc", "w", format="NETCDF4")

time_out = outnetcdf.createDimension("time", len(time))
catchment_id = outnetcdf.createDimension("catchment_id", len(mulpol))

time_out = outnetcdf.createVariable('time',time.datatype,time.dimensions)
catchment_id = outnetcdf.createVariable('catchment_id','i4')

# pr_out = outnetcdf.createVariable('pr_out','f4',("time"


##ex loop to recycle the variables from the other netcdf file_ecv
# for name, variable in src.variables.iteritems():

    ##take out the variable you don't want
    # if name == 'some_variable': 
      # continue

    # x = dst.createVariable(name, variable.datatype, variable.dimensions)
    # dst.variables[x][:] = src.variables[x][:]
    
    
    
# write strings to netcdf (for attibutes table from shapefile): https://stackoverflow.com/questions/23780324/multicharacter-strings-to-netcdf-with-python-netcdf4