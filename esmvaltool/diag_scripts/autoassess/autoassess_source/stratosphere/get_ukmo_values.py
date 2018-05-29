'''
Routine to get Met Office metrics for QBO
'''

import numpy as np

import iris
import iris.coord_categorisation as icc
from iris.fileformats.pp import STASH

from utility.method_constraint import MethodConstraint
from stratosphere.strat_metrics_1 import mean_and_strength

# for now, using datetime objects for constraints needs to be switched on:
iris.FUTURE.cell_datetime_objects = True

# Callback to normalise data
def stash_callback(cube, field, filename):
    if str(cube.attributes['STASH']) == 'm??s16i203':
        cube.attributes['STASH'] = STASH(model=1, section=16, item=203)
    if cube.coords('pseudo_level'):
        cube.remove_coord('pseudo_level')
    lats = np.linspace(-90, 90, cube.coord('latitude').points.size)
    lons = np.linspace(0, 360, cube.coord('longitude').points.size,
                       endpoint=False)
    cube.coord('latitude').points = lats
    cube.coord('longitude').points = lons
    return

# Set up glob pattern for files to read in
files = '/data/nwp1/maguu/ppassm/ucormean/ppassm_ucormean_y??_m??.pp'

# Create constraints for data to extract from files
stash = iris.AttributeConstraint(STASH='m01s16i203')
proc = MethodConstraint('time', method='mean')
level = iris.Constraint(pressure=100., latitude=0.0)
year = iris.Constraint(time=lambda cell: 1992<=cell.point.year<=2001)

# Load data
constraint = stash & proc & level & year
cube = iris.load_cube(files, constraint, callback=stash_callback)

# Create zonal mean
zmean = cube.collapsed('longitude', iris.analysis.MEAN)
# Aggregate over years to create seasonal cycle
icc.add_month(cube, 'time', name='month')
tzmean = zmean.aggregated_by('month', iris.analysis.MEAN)
# Calculate mean and strength using same code as metrics
(tmean, tstrength) = mean_and_strength(tzmean)

# Write out results
print(tzmean.data)
print(tmean)
print(tstrength)

