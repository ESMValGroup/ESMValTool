"""Module for permafrost metrics."""

import os

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

import iris
import iris.coord_categorisation
import iris.quickplot as qplt
import iris.util as ut

from esmvaltool.diag_scripts.autoassess.loaddata import load_run_ss
# from esmvaltool.diag_scripts.shared._supermeans import get_supermean
from . import permafrost_koven_sites


# main permafrost subroutine
def land_permafrost_top(run):
    """
    Make permafrost metrics.

    Code development Eleanor Burke.

    Arguments:
        run - dictionary containing model run metadata
              (see auto_assess/model_run.py for description)

    Returns:
        metrics - dictionary of metrics names and values
        also produces image files in the current working directory

    """
    metrics = dict()

    # Load the whole monthly soil temperature time-sequence as a single cube.
    # STASH m01s08i225
    period = "monthly"
    soiltemp = load_run_ss(
        run, period, 'soil_temperature')  # has dims(time, depth, lat, long)

    # check soil depths
    expected_soil_depths = [0.05, 0.225, 0.675, 2.0]
    soil_depths = soiltemp.coord('depth')
    if np.array_equal(soil_depths, expected_soil_depths):
        msg = ('Soil has changed levels from usual {}, '
               'following not supported: {}'.format(expected_soil_depths,
                                                    soil_depths))
        raise Exception(msg)

    # load the whole monthly air temperature (STASH m01s03i236)
    airtemp = load_run_ss(run, period, 'air_temperature')

    # get the land fraction mask using whole cube
    landfrac = get_landfr_mask(run)

    # extract northern latitudes
    airtemp = airtemp.extract(iris.Constraint(latitude=lambda cell: cell > 0))
    soiltemp = soiltemp.extract(iris.Constraint(
        latitude=lambda cell: cell > 0))

    # calculate the permafrost area and fraction less than zero
    # permafrost_area returns a dict, which is added to the main metrics dict
    # by metrics.update()
    metrics.update(permafrost_area(soiltemp, airtemp, landfrac, run))

    # calculate the koven temperature metrics
    metrics.update(koven_temp_offsets(soiltemp, airtemp))
    metrics.update(koven_temp_atten(soiltemp, airtemp))

    return metrics


def permafrost_area(soiltemp, airtemp, landfrac, run):
    """Calculate the permafrost area and make a plot."""
    # Define parameters of the test to calculate the existence of permafrost
    thresh_temperature = 273.2
    frozen_months = 24
    prop_months_frozen = 0.5  # frozen for at least half of the simulation

    # make a mask of land fraction over non iced areas and extract northern
    # latitudes
    nonice = get_nonice_mask(run)
    mask = iris.analysis.maths.multiply(nonice, landfrac)
    mask = mask.extract(iris.Constraint(latitude=lambda cell: cell > 0))

    # extract northern high latitudes [and deeepst soil level]
    soiltemp = soiltemp.extract(iris.Constraint(depth=2.0))  # from 1m to 3m

    # Make an aggregator to define the permafrost extent
    # I dont really understand this but it works
    frozen_count = iris.analysis.Aggregator(
        'frozen_count', num_frozen, units_func=lambda units: 1)

    # Calculate the permafrost locations
    pf_periods = soiltemp.collapsed(
        'time',
        frozen_count,
        threshold=thresh_temperature,
        frozen_length=frozen_months)
    tot_time = len(soiltemp.coord('time').points)
    pf_periods = pf_periods / float(tot_time)
    pf_periods.rename('Fraction of months layer 4 (-1m to -3m) soil is frozen')

    # mask out non permafrost points, sea points and ice points
    pf_periods.data = np.ma.masked_less(pf_periods.data, prop_months_frozen)
    # set all non-masked values to 1 for area calculation
    # (may be a better way of doing this but I'm not sure what it is)
    pf_periods = pf_periods / pf_periods
    # mask for land area also
    pf_periods = pf_periods * mask

    # calculate the area of permafrost
    # Generate area-weights array. This method requires bounds on lat/lon
    # coords, add some in sensible locations using the "guess_bounds"
    # method.
    for coord in ['latitude', 'longitude']:
        if not pf_periods.coord(coord).has_bounds():
            pf_periods.coord(coord).guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(pf_periods)
    # calculate the areas not masked in pf_periods
    pf_area = pf_periods.collapsed(
        ['longitude', 'latitude'], iris.analysis.SUM, weights=grid_areas).data

    # what is the area where the temperature is less than 0 degrees C?
    airtemp = airtemp.collapsed('time', iris.analysis.MEAN)
    # if more than 2 dims, select the ground level
    if airtemp.ndim > 2:
        airtemp = airtemp[0]
    airtemp_below_zero = np.where(airtemp.data < 273.2, 1, 0)
    airtemp_area = np.sum(airtemp_below_zero * grid_areas)

    pf_prop = pf_area / airtemp_area
    pf_area = pf_area / 1e12

    # Figure Permafrost extent north america
    plt.figure(figsize=(8, 8))
    ax = plt.axes(
        projection=ccrs.Orthographic(
            central_longitude=-80.0, central_latitude=60.0))
    qplt.pcolormesh(pf_periods)
    ax.gridlines()
    ax.coastlines()
    levels = [thresh_temperature]
    qplt.contour(airtemp, levels, colors='k', linewidths=3)
    plt.title('Permafrost extent & zero degree isotherm ({})'.format(
        run['runid']))
    plt.savefig('pf_extent_north_america_' + run['runid'] + '.png')

    # Figure Permafrost extent asia
    plt.figure(figsize=(8, 8))
    ax = plt.axes(
        projection=ccrs.Orthographic(
            central_longitude=100.0, central_latitude=50.0))
    qplt.pcolormesh(pf_periods)
    ax.gridlines()
    ax.coastlines()
    levels = [thresh_temperature]
    qplt.contour(airtemp, levels, colors='k', linewidths=3)
    plt.title('Permafrost extent & zero degree isotherm ({})'.format(
        run['runid']))
    plt.savefig('pf_extent_asia_' + run['runid'] + '.png')

    # defining metrics for return up to top level
    metrics = {
        'permafrost area': pf_area,
        'fraction area permafrost over zerodeg': pf_prop,
    }

    return metrics


# define the frozen area
def num_frozen(data, threshold, axis, frozen_length):
    """
    Count valid frozen points.

    Function to calculate the number of points in a sequence where the value
    is less than freezing for at least a certain number of timepoints.

    Generalised to operate on multiple time sequences arranged on a specific
    axis of a multidimensional array.
    """
    if axis < 0:
        # just cope with negative axis numbers
        axis += data.ndim

    # Threshold the data to find the 'significant' points.
    data_hits = data < threshold
    # Make an array with data values "windowed" along the time axis.
    hit_windows = ut.rolling_window(data_hits, window=frozen_length, axis=axis)
    # Find the windows "full of True-s" (along the added 'window axis').
    full_windows = np.all(hit_windows, axis=axis + 1)
    # Count points fulfilling the condition (along the time axis).
    frozen_point_counts = np.sum(full_windows, axis=axis, dtype=int)

    return frozen_point_counts


# land fraction
def get_landfr_mask(run):
    """Get the land fraction mask."""
    supermean_data_dir = os.path.join(run['data_root'], run['runid'],
                                      run['_area'] + '_supermeans')
    # m01s03i395
    # TODO: replacing time-varying mask with fixed sftfx
    # cube = get_supermean('land_area_fraction', 'ann', supermean_data_dir)
    # replaced momentarily with:
    name_constraint = iris.Constraint(name='land_area_fraction')
    cubes_path = os.path.join(supermean_data_dir, 'cubeList.nc')
    cubes = iris.load(cubes_path)
    cube = cubes.extract_cube(name_constraint)

    return cube


# land ice mask
def get_nonice_mask(run):
    """
    Get the land points without ice.

    Need to read the soil moisture data from the supermeans
    """
    # TODO: currently set to mrsofc: soil_moisture_content_at_field_capacity
    supermean_data_dir = os.path.join(run['data_root'], run['runid'],
                                      run['_area'] + '_supermeans')

    # m01s08i223
    # TODO: original code
    # cube = get_supermean('moisture_content_of_soil_layer', 'ann',
    #                      supermean_data_dir)
    # replaced with new time-invariant variable
    name_constraint = iris.Constraint(
        name='soil_moisture_content_at_field_capacity')
    cubes_path = os.path.join(supermean_data_dir, 'cubeList.nc')
    cubes = iris.load(cubes_path)
    cube = cubes.extract_cube(name_constraint)

    # TODO: mrsofc does not have depth
    # cube = cube.extract(iris.Constraint(depth=2.0))  # layer from 1m to 3m

    # make it into a mask of ones - extract first layer
    # use masked_values for floating point fuzzy equals
    cube.data = np.ma.masked_values(cube.data, 0.0)
    cube = cube / cube

    return cube


def extract_sites(ex_points, cube):
    """Extract points for the sites given."""
    tempsite = cube.interpolate(ex_points, iris.analysis.Linear())
    tempsite = np.diagonal(tempsite.data)
    tempsite = np.ma.masked_array(tempsite)
    tempsite = np.ma.masked_less(tempsite, 0.0)
    return tempsite


def koven_temp_offsets(soiltemp, airtemp):
    """Define thermal offsets in Koven et al 2013."""
    # read in list of observed lats and lons from Koven paper
    ex_points = permafrost_koven_sites.site_points

    # interpolate to depth required
    # the soil temperatures are for the middle of the layer not the bottom of
    # the layer
    linear = iris.analysis.Linear()
    soiltemp_surf = soiltemp.interpolate([('depth', 0.0)], linear)
    soiltemp_1m = soiltemp.interpolate([('depth', 1.0)], linear)

    # extract points for eachsite
    airtemp_1d = extract_sites(ex_points, airtemp)
    if len(airtemp_1d.shape) > 2:
        airtemp_1d = airtemp_1d[:, :, 0]
    soiltemp_surf_1d = extract_sites(ex_points, soiltemp_surf)
    soiltemp_1m_1d = extract_sites(ex_points, soiltemp_1m)

    # assign metrics
    metrics = {}
    metrics['offset 1m minus surface'] = np.median(soiltemp_1m_1d -
                                                   soiltemp_surf_1d)
    metrics['offset surface minus air'] = np.median(soiltemp_surf_1d -
                                                    airtemp_1d)
    return metrics


def make_monthly_amp(cube):
    """Make monthly climatology."""
    iris.coord_categorisation.add_month(cube, 'time', name='month')
    cube_clim = cube.aggregated_by('month', iris.analysis.MEAN)
    cube_ampl = cube_clim.collapsed('time', iris.analysis.MAX) - \
        cube_clim.collapsed('time', iris.analysis.MIN)
    return cube_ampl


def koven_temp_atten(soiltemp, airtemp):
    """Define thermal attenuation ratios as in Koven et al 2013."""
    # read in list of observed lats and lons from Koven paper
    ex_points = permafrost_koven_sites.site_points

    # make amplitudes
    airtemp_ampl = make_monthly_amp(airtemp)
    soiltemp_ampl = make_monthly_amp(soiltemp)

    # interpolate the log to the correct depth
    soiltemp_log = iris.analysis.maths.log(soiltemp_ampl)
    linear = iris.analysis.Linear()
    soiltemp_log_surf = soiltemp_log.interpolate([('depth', 0.0)], linear)
    soiltemp_ampl_surf = iris.analysis.maths.exp(soiltemp_log_surf)
    soiltemp_log_1m = soiltemp_log.interpolate([('depth', 1.0)], linear)
    soiltemp_ampl_1m = iris.analysis.maths.exp(soiltemp_log_1m)

    # extract points for eachsite
    airtemp_ampl_1d = extract_sites(ex_points, airtemp_ampl)
    if len(airtemp_ampl_1d.shape) > 1:
        airtemp_ampl_1d = airtemp_ampl_1d[:, 0]
    soiltemp_ampl_surf_1d = extract_sites(ex_points, soiltemp_ampl_surf)
    soiltemp_ampl_1m_1d = extract_sites(ex_points, soiltemp_ampl_1m)

    # assign metrics
    metrics = {}
    metrics['attenuation 1m over surface'] = np.median(
        soiltemp_ampl_1m_1d / soiltemp_ampl_surf_1d)
    metrics['attenuation surface over air'] = np.median(
        soiltemp_ampl_surf_1d / airtemp_ampl_1d)

    return metrics
