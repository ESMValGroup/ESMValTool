'''
Module to hold functions useful to area assessments
'''

import numpy as np
import numpy.ma as ma

import iris
import iris.analysis.cartography as iac

import iris_updates as newiris


def area_average(cube,
                 weighted=True,
                 mask=None,
                 logicmask=False,
                 coords=None,
                 aggregator=iris.analysis.MEAN,
                 **aggkeys):
    '''
    Routine to calculate weighted horizontal area aggregations

    Routine defaults to longitude and latitude, but can be configured to
    collapse over any coordinate in the cube

    Inputs:

    cube = cube to aggregate

    Keywords:

    weighted = perform area weighted aggregation (default: True)
    mask = cube containing mask data (default: None)
    logicmask = Does mask contain logical data (default: False)
    aggregator = aggregator for collapsed method (default: iris.analysis.MEAN)
    coords = list of coordinates to collapse cube over
             (default: ["latitude", "longitude"])
    "coord" = (coord_min, coord_max) - range of coordinate to collapse over
    **kwargs = any keywords required for the aggregator

    Return:

    aggregated cube

    '''
    if coords is None:
        coords = ['latitude', 'longitude']

    # Make sure that aggregator is an Aggregator instance
    assert isinstance(aggregator, iris.analysis.Aggregator)
    # If doing weighted aggregation make sure that aggregator
    # is a WeightAggregator instance
    if weighted:
        assert isinstance(aggregator, iris.analysis.WeightedAggregator)

    # Extract region specification if available
    intkeys = {}
    for coord in coords:
        if coord in aggkeys:
            intkeys[coord] = aggkeys.pop(coord)

    # Extract region if required
    if intkeys:
        newcube = cube.intersection(ignore_bounds=True, **intkeys)
        # For some reason cube.intersection() promotes dtype of coordinate
        # arrays to float64, whereas cube.extract() doesn't. Need to make
        # sure behaviour is identical.
        for coord in intkeys.keys():
            newcube.coord(coord).points = \
                newcube.coord(coord).points.astype(np.float32, copy=False)
    else:
        newcube = cube.copy()

    # If doing area-weighted aggregation then calculate area weights
    if weighted:
        # Coords need bounding
        for coord in coords:
            if not newcube.coord(coord).has_bounds():
                # Test to make sure latitude bounds do not wrap over pole.
                if coord in ['latitude', 'grid_latitude']:
                    newiris.guess_bounds(
                        newcube.coord(coord), bound_min=-90., bound_max=90.)
                else:
                    newcube.coord(coord).guess_bounds()
        aggkeys['weights'] = iac.area_weights(newcube)

    # Apply mask
    if mask:
        # Extract region of mask to match data
        if intkeys:
            newmask = mask.intersection(ignore_bounds=True, **intkeys)
        else:
            newmask = mask.copy()
        # Apply mask to weights if they exist, else apply to data
        # Do I really need two methods here?
        if 'weights' in aggkeys:
            if logicmask:
                aggkeys['weights'] = ma.array(
                    data=aggkeys['weights'], mask=newmask.data)
            else:
                aggkeys['weights'] *= newmask.data
        else:
            if logicmask:
                newcube.data = ma.array(data=newcube.data, mask=newmask.data)
            else:
                newcube.data *= newmask.data

    return newcube.collapsed(coords, aggregator, **aggkeys)
