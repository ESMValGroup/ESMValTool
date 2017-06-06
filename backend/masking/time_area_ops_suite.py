########################################################################
 # TIME AND AREA OPERATIONS
 # V.Predoi, University of Reading, May 2017
########################################################################
import iris

# slice cube over a restricted time period
def time_slice(mycube,yr1,mo1,d1,yr2,mo2,d2):
    """
    Function that returns a subset of the original cube (slice)
    given two dates of interest date1 and date2
    date1 and date2 should be given in a yr,mo,d (int)format e.g.
    time_slice(cube,2006,2,2,2010,1,1) or time_slice(cube,'2006','2','2','2010','1','1');
    Returns a cube
    """
    import datetime
    import iris.unit
    myDate1 = datetime.datetime(int(yr1),int(mo1),int(d1))
    myDate2 = datetime.datetime(int(yr2),int(mo2),int(d2))
    t1 = mycube.coord('time').units.date2num(myDate1)
    t2 = mycube.coord('time').units.date2num(myDate2)
    myConstraint = iris.Constraint(time=lambda t: t.point > t1 and t.point < t2)
    cubeslice = mycube.extract(myConstraint)
    return cubeslice

# slice cube over a restricted area (box)
def area_slice(mycube, long1, long2, lat1, lat2):
    """
    Function that subsets a cube on a box (long1,long2,lat1,lat2)
    This function is a restriction of masked_cube_lonlat();
    Returns a cube
    """
    sublon = iris.Constraint(longitude=lambda cell: float(long1) <= cell <= float(long2))
    sublat = iris.Constraint(latitude=lambda cell: float(lat1) <= cell <= float(lat2))
    region_subset = mycube.extract(sublon & sublat)
    return region_subset

# get the time average
def time_average(mycube):
    """
    Function to get the time average over MEAN;
    Returns a cube
    """
    var_mean = mycube.collapsed('time', iris.analysis.MEAN)
    return var_mean

# get the probability a value is greater than a threshold
def proportion_greater(mycube, coord1, threshold):
    """
    Function that returns the probability
    that a cetain variable coord1 (string) is greater than
    a threshold threshold (float or string), across a cube mycube;
    Returns a cube
    """
    thr = float(threshold)
    result = mycube.collapsed(coord1, iris.analysis.PROPORTION,
                              function=lambda values: values > thr)
    return result

# get zonal means
def zonal_means(mycube, coord1, mean_type):
    """
    Function that returns zonal means along a coordinate coord1;
    the type of mean is controlled by mean_type variable (string):
        'mean' -> MEAN
        'stdev' -> STD_DEV
        'variance' -> VARIANCE

    Returns a cube
    """
    if mean_type=='mean':
        result = mycube.collapsed(coord1, iris.analysis.MEAN)
    elif mean_type=='stdev':
        result = mycube.collapsed(coord1, iris.analysis.STD_DEV)
    elif mean_type=='variance':
        result = mycube.collapsed(coord1, iris.analysis.VARIANCE)
    return result

# get the area average
def area_average(mycube, coord1, coord2):
    """
    Function that determines the area average
    Can be used with coord1 and coord2 (strings,
    usually 'longitude' and 'latitude' but depends on the cube);
    Returns a cube
    """
    import iris.analysis.cartography
    mycube.coord(coord1).guess_bounds()
    mycube.coord(coord2).guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(mycube)
    result = mycube.collapsed([coord1, coord2], iris.analysis.MEAN, weights=grid_areas)
    return result

# get the seasonal mean
def seasonal_mean(mycube):
    """
    Function to compute seasonal means with MEAN
    Chunks time in 3-month periods and computes means over them;
    Returns a cube
    """
    import iris.coord_categorisation
    iris.coord_categorisation.add_season(mycube, 'time', name='clim_season')
    iris.coord_categorisation.add_season_year(mycube, 'time', name='season_year')
    annual_seasonal_mean = mycube.aggregated_by(['clim_season', 'season_year'], iris.analysis.MEAN)
    spans_three_months = lambda time: (time.bound[1] - time.bound[0]) == 2160
    three_months_bound = iris.Constraint(time=spans_three_months)
    return annual_seasonal_mean.extract(three_months_bound)

# operate along a trajectory line
def trajectory_cube(mycube, long1, long2, lat1, lat2, plong1, plong2, plat1, plat2,samplecounts):
    """
    Function that subsets a cube on a box (long1,long2,lat1,lat2)
    then creates a trajectory with waypoints (plong1,plong2,plat1, plat2),
    populates it with samplecounts number of points
    and subsets the cube along the trajectory
    """
    from iris.analysis import trajectory
    sublon = iris.Constraint(longitude=lambda cell: float(long1) <= cell <= float(long2))
    sublat = iris.Constraint(latitude=lambda cell: float(lat1) <= cell <= float(lat2))
    wspd_subset = mycube.extract(sublon & sublat)
    pnts = [{'longitude': float(plong1), 'latitude': float(plat1)}, {'longitude': float(plong2), 'latitude': float(plat2)}]
    traj = trajectory.Trajectory(pnts, sample_count=int(samplecounts))
    lon = [d['longitude'] for d in traj.sampled_points]
    lat = [d['latitude'] for d in traj.sampled_points]
    sampled_points = [('longitude', lon),('latitude', lat)]
    section = trajectory.interpolate(wspd_subset, sampled_points)
    lon, lat = wspd_subset.coord('longitude').points, wspd_subset.coord('latitude').points
    return section, lon, lat
