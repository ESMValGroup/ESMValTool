########################################################################
# EXAMPLE FOR MASKING AND TIME AND AREA OPERATIONS
# V.Predoi, University of Reading, May 2017
########################################################################
import iris
from iris.analysis import Aggregator
import iris.plot as iplt
import iris.quickplot as qplt
import numpy as np
import matplotlib.pyplot as plt
from mask_suite import *
from time_area_ops_suite import *
from regrid_Bill import regrid
import cartopy.io.shapereader as shpreader
from cartopy.mpl.patch import geos_to_path
from shapely.geometry import Point
import cartopy.crs as ccrs
import iris.plot as iplt
import sys
import subprocess

print('All output from running this code can be found in example_output directory.')
dirname = 'example_output'
mkd = 'mkdir -p ' + dirname
proc = subprocess.Popen(mkd, stdout=subprocess.PIPE, shell=True)
(out, err) = proc.communicate()

outfile = open('example_output/Example_explained.txt', 'w')
print >> outfile
print >> outfile, "########################################################"
print >> outfile, "# Example how to use new Mask and Time-Area Functions  #"
print >> outfile, "########################################################"
print >> outfile
print >> outfile, "FUNCTION FILES USED: mask_suite, time_area_ops_suite  #"

# let's load some data
print >> outfile, "########################################################"
print >> outfile, "#       1. Get some data and build some cubes           #"
print >> outfile, "########################################################"
print >> outfile, "Let's get some data: IRIS example data E1_north_america.nc  #"
filenc = iris.sample_data_path('E1_north_america.nc')
cube = iris.load_cube(filenc)
box = area_slice(cube, 225., 315., 40., 60.)
# let's apply the counts mask and get two different cubes
# apply the counts mask with different windows and thresholds 
# to simulate the effect of having two different files
print >> outfile, "We make two different cubes out of the data  #"
print >> outfile, "by thresholding on temperature - one on 260 and the other on 270 K's #"
th_box = mask_threshold(box, 1.0)
th_box1 = mask_threshold(box, 260.0)
th_box2 = mask_threshold(box, 270.0)
# here 'cube' is the reference cube so let's regrid
# cube1 and cube2 on cube
print >> outfile, "We assume the cubes have already been regridded to a reference grid  #"
print >> outfile, "If not, use regrid()  from the main regridding suite #"
# let's say we have a region of interest
# and select it out of the cubes
print >> outfile, "########################################################"
print >> outfile, "#        2. Mask out a region of interest              #"
print >> outfile, "########################################################"
print >> outfile, "Let's now select a region of interest - #"
print >> outfile, " - since data shows North America we will select Canada as our area of interest  #"
shpfilename = shpreader.natural_earth(resolution='110m',
                                      category='cultural',
                                      name='admin_0_countries')
Canada_geom = maskgeometry(shpfilename, 'name', 'Canada')
print >> outfile, "But first let's build a box around the rough corners of Canada, save us computation #"
mcube = mask_2d(th_box,Canada_geom)
mcube1 = mask_2d(th_box1,Canada_geom)
mcube2 = mask_2d(th_box2,Canada_geom)
tmcube = mask_threshold(mcube, 1.0)
tmcube1 = mask_threshold(mcube1, 1.0)
tmcube2 = mask_threshold(mcube2, 1.0)
print >> outfile, "Min/Max temperatures: Canada all, Canada thresholded 260K, Canada thresholded 270K:"
print >> outfile, np.min(tmcube.data), np.max(tmcube.data)
print >> outfile, np.min(tmcube1.data), np.max(tmcube1.data)
print >> outfile, np.min(tmcube2.data), np.max(tmcube2.data)
############
cube = iris.load_cube(filenc)
qplt.contourf(cube[0,:,:], cmap='RdYlBu_r',bbox_inches='tight')
plt.title('Original data')
plt.gca().coastlines()
plt.tight_layout()
plt.savefig('example_output/Original_Data.png')
plt.close()
print >> outfile, "Plotted Original_Data at initial time t0 #"
#############
qplt.contourf(tmcube[0,:,:], cmap='RdYlBu_r',bbox_inches='tight')
plt.title('Boxed and selected Canada')
plt.gca().coastlines()
plt.tight_layout()
plt.savefig('example_output/T_t0_Tgt0_Canada.png')
plt.close()
print >> outfile, "Plotted Canada data at time t0 T_t0_Tgt0_Canada.png #"
############
qplt.contourf(tmcube1[0,:,:], cmap='RdYlBu_r',bbox_inches='tight')
plt.title('Canada: area of temp. > 260')
plt.gca().coastlines()
plt.tight_layout()
plt.savefig('example_output/T_t0_Tgt260_Canada.png')
plt.close()
print >> outfile, "Plotted Canada data at time t0 and T>260K T_t0_Tgt260_Canada.png #"
############
qplt.contourf(tmcube2[0,:,:], cmap='RdYlBu_r',bbox_inches='tight')
plt.title('Canada: area of temp. > 270')
plt.gca().coastlines()
plt.tight_layout()
plt.savefig('example_output/T_t0_Tgt270_Canada.png')
plt.close()
print >> outfile, "Plotted Canada data at time t0 and T>270K T_t0_Tgt270_Canada.png #"
############
print >> outfile, "########################################################"
print >> outfile, "# 3. Missing data: dividing the time in dt time blocks #"
print >> outfile, "########################################################"
print >> outfile, "Now let's look at 5-year averages, in this case the time block dt=5 years #"
s = window_counts(tmcube, 1.0, 5, 68)
s1 = window_counts(tmcube1, 260.0, 5, 68)
s2 = window_counts(tmcube2, 270.0, 5, 68)
plt.subplot(311)
plt.hist(s[0],bins=50)
plt.xlim([1,50])
plt.xlabel('Counts per lon-lat grid point')
plt.ylabel('Number of lon-lat grid points')
plt.grid()
plt.title('Number of temperature data points every 5 years\nCanada: no temp. threshold')
plt.subplot(312)
plt.hist(s1[0],bins=50)
plt.xlim([1,50])
plt.xlabel('Counts per lon-lat grid point')
plt.ylabel('Number of lon-lat grid points')
plt.grid()
plt.title('Number of temperature data points every 5 years\nCanada: T>260K')
plt.subplot(313)
plt.hist(s2[0],bins=50,log=True)
plt.xlim([1,50])
plt.xlabel('Counts per lon-lat grid point')
plt.ylabel('Number of lon-lat grid points')
plt.grid()
plt.title('Number of temperature data points every 5 years\nCanada: T>270K')
plt.savefig('example_output/Five_Year_Averages.png')
plt.close()
print >> outfile, "Plotted Canada Five year Averages example_output/Five_Year_Averages.png#"
print >> outfile, "#  Counts per lon-lat grid point - 68 percentiles are: #"
print >> outfile, "--------------------------------------------------------"
print >> outfile, int(s[3])
print >> outfile, int(s1[3])
print >> outfile, int(s2[3])
print >> outfile, "#################################################################"
print >> outfile, "#   4. Missing data: thresholding on a minimum number           #"
print >> outfile, "#      of data points per time block dt per lon-lat grid point  #"
print >> outfile, "#################################################################"
print >> outfile, "Suppose now we want to apply a mask that discards all the grid #"
print >> outfile, "points that have less than 20 temperature data points in each 5 year block #"
c = mask_cube_counts(tmcube, 1.0, 20, 5)[2]
c1 = mask_cube_counts(tmcube1, 1.0, 20, 5)[2]
c2 = mask_cube_counts(tmcube2, 1.0, 20, 5)[2]
###############
qplt.contourf(c[0,:,:], cmap='RdYlBu_r',bbox_inches='tight')
plt.title('Canada: masked >20 points/5yr window/grid point')
plt.gca().coastlines()
plt.tight_layout()
plt.savefig('example_output/T_t0_Tgt0_20pts_Canada.png')
plt.close()
print >> outfile, "Plotted Canada data at time t0 and >20 points/5years/grid point T_t0_Tgt0_20pts_Canada.png #"
qplt.contourf(c1[0,:,:], cmap='RdYlBu_r',bbox_inches='tight')
plt.title('Canada: masked >20 points/5yr window/grid point and T>260K')
plt.gca().coastlines()
plt.tight_layout()
plt.savefig('example_output/T_t0_Tgt260_20pts_Canada.png')
plt.close()
print >> outfile, "Plotted Canada data at time t0 and >20 points/5years/grid point and T>260K T_t0_Tgt260_20pts_Canada.png #"
qplt.contourf(c2[0,:,:], cmap='RdYlBu_r',bbox_inches='tight')
plt.title('Canada: masked >20 points/5yr window/grid point and T>270K')
plt.gca().coastlines()
plt.tight_layout()
plt.savefig('example_output/T_t0_Tgt270_20pts_Canada.png')
plt.close()
print >> outfile, "Plotted Canada data at time t0 and >20 points/5years/grid point and T>270K T_t0_Tgt270_20pts_Canada.png #"
###################
print >> outfile, "########################################################"
print >> outfile, "#  5. Time averaging after applying missing data mask  #"
print >> outfile, "########################################################"
print >> outfile, "Let us plot the time averages of the original and masked cubes:"
cube = iris.load_cube(filenc)
qplt.contourf(time_average(cube), cmap='RdYlBu_r',bbox_inches='tight')
plt.title('Original data: time avg')
plt.gca().coastlines()
plt.tight_layout()
plt.savefig('example_output/Original_Data_TimeAvg.png')
plt.close()
print >> outfile, "Plotted Canada data average time Original_Data_TimeAvg.png #"
qplt.contourf(time_average(c), cmap='RdYlBu_r',bbox_inches='tight')
plt.title('Canada: masked >20 points/5yr window/grid point: time avg')
plt.gca().coastlines()
plt.tight_layout()
plt.savefig('example_output/T_tavg_Tgt0_20pts_Canada.png')
plt.close()
print >> outfile, "Plotted Canada data average time and >20 points/5years/grid point T_tavg_Tgt0_20pts_Canada.png #"
qplt.contourf(time_average(c1), cmap='RdYlBu_r',bbox_inches='tight')
plt.title('Canada: masked >20 points/5yr window/grid point and T>260K: time avg')
plt.gca().coastlines()
plt.tight_layout()
plt.savefig('example_output/T_tavg_Tgt260_20pts_Canada.png')
plt.close()
print >> outfile, "Plotted Canada data average time and >20 points/5years/grid point and T>260K T_tavg_Tgt260_20pts_Canada.png #"
qplt.contourf(mask_threshold(time_average(c2),1.0), cmap='RdYlBu_r',bbox_inches='tight')
plt.title('Canada: masked >20 points/5yr window/grid point and T>270K: time avg')
plt.gca().coastlines()
plt.tight_layout()
plt.savefig('example_output/T_tavg_Tgt270_20pts_Canada.png')
plt.close()
print >> outfile, "Plotted Canada data average time and >20 points/5years/grid point and T>270K T_tavg_Tgt270_20pts_Canada.png #"
print >> outfile, "########################################################"
print >> outfile, "#  6. Proportions after applying missing data mask  #"
print >> outfile, "########################################################"
print >> outfile, "Let us plot the proportions temperatures >270K over latitude"
cp = proportion_greater(c, 'latitude', 270.)
cp1 = proportion_greater(c1, 'latitude', 270.)
cp2 = proportion_greater(mask_threshold(c2,1.0), 'latitude', 270.)
qplt.contourf(proportion_greater(cube,'latitude',270), cmap='RdYlBu_r',bbox_inches='tight')
plt.title('Original data: proportion T>270K over latitude')
plt.grid()
plt.tight_layout()
plt.savefig('example_output/Original_Data_PropTgt270LAT.png')
plt.close()
qplt.contourf(cp, cmap='RdYlBu_r',bbox_inches='tight')
plt.title('Canada: masked >20 points/5yr window/grid point:\nproportion T>270K over latitude')
plt.grid()
plt.tight_layout()
plt.savefig('example_output/T_20pts_PropTgt270LAT_Canada.png')
plt.close()
qplt.contourf(cp1, cmap='RdYlBu_r',bbox_inches='tight')
plt.title('Canada: masked >20 points/5yr window/grid point and T>260K:\nproportion T>270K over latitude')
plt.grid()
plt.tight_layout()
plt.savefig('example_output/T_Tgt260_20pts_PropTgt270LAT_Canada.png')
plt.close()
qplt.contourf(cp2, cmap='RdYlBu_r',bbox_inches='tight')
plt.title('Canada: masked >20 points/5yr window/grid point and T>270K:\nproportion T>270K over latitude')
plt.grid()
plt.tight_layout()
plt.savefig('example_output/T_Tgt270_20pts_PropTgt270LAT_Canada.png')
plt.close()
print >> outfile, "Proportion plots are labelled _PropTgt270LAT_"


