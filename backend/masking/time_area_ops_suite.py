########################################################################
 # TIME AND AREA OPERATIONS
 # V.Predoi, University of Reading, May 2017
########################################################################

# slice cube over a restricted time period
import datetime
import iris.unit
def time_slice(mycube,yr1,mo1,d1,yr2,mo2,d2):
    """
    Function that returns a subset of the original cube (slice)
    given two dates of interest date1 and date2
    date1 and date2 should be given in a yr,mo,d (int)format e.g.
    time_slice(cube,2006,2,2,2010,1,1) or time_slice(cube,'2006','2','2','2010','1','1')
    """
    myDate1 = datetime.datetime(int(yr1),int(mo1),int(d1))
    myDate2 = datetime.datetime(int(yr2),int(mo2),int(d2))
    t1 = mycube.coord('time').units.date2num(myDate1)
    t2 = mycube.coord('time').units.date2num(myDate2)
    myConstraint = iris.Constraint(time=lambda t: t.point > t1 and t.point < t2)
    cubeslice = mycube.extract(myConstraint)
    return cubeslice

# get the area average
import iris.analysis.cartography
def area_average(mycube, coord1, coord2):
    """
    Function that determines the area average
    Can be used with coord1 and coord2 (strings,
    usually 'longitude' and 'latitude' but depends on the cube)
    """
    mycube.coord(coord1).guess_bounds()
    mycube.coord(coord2).guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(mycube)
    new_cube = mycube.collapsed([coord1, coord2], iris.analysis.MEAN, weights=grid_areas)
    return new_cube

# get the seasonal mean
import iris.coord_categorisation
def seasonal_mean(mycube):
    iris.coord_categorisation.add_season(mycube, 'time', name='clim_season')
    iris.coord_categorisation.add_season_year(mycube, 'time', name='season_year')
    annual_seasonal_mean = mycube.aggregated_by(['clim_season', 'season_year'], iris.analysis.MEAN)
    spans_three_months = lambda time: (time.bound[1] - time.bound[0]) == 2160
    three_months_bound = iris.Constraint(time=spans_three_months)
    return annual_seasonal_mean.extract(three_months_bound)
