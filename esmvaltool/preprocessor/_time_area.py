"""
Time operations on cubes

Allows for selecting data subsets using certain time bounds;
constructing seasonal and area averages.
"""
import iris
import iris.coord_categorisation
import numpy as np


# slice cube over a restricted time period
def time_slice(mycube, yr1, mo1, day_1, yr2, mo2, day_2):
    """
    Slice cube on time

    Function that returns a subset of the original cube (slice)
    given two dates of interest date1 and date2
    date1 and date2 should be given in a yr,mo,d (int)format e.g.
    time_slice(cube,2006,2,2,2010,1,1) or
    time_slice(cube,'2006','2','2','2010','1','1');

    Returns a cube
    """
    import datetime
    time_units = mycube.coord('time').units
    if time_units.calendar == '360_day':
        if day_1 > 30:
            day_1 = 30
        if day_2 > 30:
            day_2 = 30
    my_date1 = datetime.datetime(int(yr1), int(mo1), int(day_1))
    my_date2 = datetime.datetime(int(yr2), int(mo2), int(day_2))

    t_1 = time_units.date2num(my_date1)
    t_2 = time_units.date2num(my_date2)
    # TODO replace the block below for when using iris 2.0
    # my_constraint = iris.Constraint(time=lambda t: (
    #     t_1 < time_units.date2num(t.point) < t_2))
    my_constraint = iris.Constraint(time=lambda t: (
        t_1 < t.point < t_2))
    cube_slice = mycube.extract(my_constraint)
    return cube_slice


def extract_season(cube, season):
    """
    Slice cube to get only the data belonging to a specific season

    Parameters
    ----------
    cube: iris.cube.Cube
        Original data
    season: str
        Season to extract. Available: DJF, MAM, JJA, SON
    """
    iris.coord_categorisation.add_season(cube, 'time', name='clim_season')
    season_cube = cube.extract(iris.Constraint(clim_season=season.lower()))
    return season_cube


def extract_month(mycube, month):
    """
    Slice cube to get only the data belonging to a specific month

    Parameters
    ----------
    cube: iris.cube.Cube
        Original data
    month: int
        Month to extract as a number from 1 to 12
    """
    season_cube = mycube.extract(iris.Constraint(month_number=month))
    return season_cube


# get the time average
def time_average(cube):
    """
    Compute time average

    Get the time average over the entire cube. The average is weighted by the
    bounds of the time coordinate.

    Arguments
    ---------
        cube: iris.cube.Cube
            input cube.

    Returns
    -------
    iris.cube.Cube
        time averaged cube.
    """
    time = cube.coord('time')
    time_thickness = time.bounds[..., 1] - time.bounds[..., 0]

    # The weights need to match the dimensionality of the cube.
    slices = [None for i in cube.shape]
    coord_dim = cube.coord_dims('time')[0]
    slices[coord_dim] = slice(None)
    time_thickness = np.abs(time_thickness[tuple(slices)])
    ones = np.ones_like(cube.data)
    time_weights = time_thickness * ones

    return cube.collapsed('time', iris.analysis.MEAN,
                          weights=time_weights)


# get the seasonal mean
def seasonal_mean(cube):
    """
    Function to compute seasonal means with MEAN

    Chunks time in 3-month periods and computes means over them;

    Arguments
    ---------
        cube: iris.cube.Cube
            input cube.

    Returns
    -------
    iris.cube.Cube
        Seasonal mean cube
    """
    iris.coord_categorisation.add_season(cube, 'time', name='clim_season')
    iris.coord_categorisation.add_season_year(
        cube, 'time', name='season_year')
    annual_seasonal_mean = cube.aggregated_by(['clim_season', 'season_year'],
                                              iris.analysis.MEAN)

    def spans_three_months(time):
        """Check for three months"""
        return (time.bound[1] - time.bound[0]) == 2160

    three_months_bound = iris.Constraint(time=spans_three_months)
    return annual_seasonal_mean.extract(three_months_bound)
