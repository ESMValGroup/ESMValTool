import numpy as np

import iris
import iris.coords
import iris.coord_systems
import iris.cube
import iris.fileformats.pp


def _create_1x1():
    cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS) 

    ydata = np.linspace(-89.5, 89.5, 180) 
    lats = iris.coords.DimCoord(ydata,
                                standard_name='latitude',
                                units='degrees_north',
                                coord_system=cs)

    xdata = np.linspace(0.5, 359.5, 360)
    lons = iris.coords.DimCoord(xdata,
                                standard_name='longitude',
                                units='degrees_east',
                                coord_system=cs)

    shape = (ydata.size, xdata.size)
    coords = [(lats, 0), (lons, 1)]
    target_1x1 = iris.cube.Cube(np.zeros(shape),
                                dim_coords_and_dims=coords)


def regrid(source_cube, target_grid, scheme):

    if scheme is None:
        if target_grid is not None:
            raise Exception('target grid must be none if no scheme is given')
        return source_cube

    if target_grid is None:
        raise Exception('target grid must not be none')
    elif target_grid == '1x1':
        # XXX: Ideally, it would be better to create and cache this once,
        # and re-use it multiple times, rather than re-create it everytime.
        target_grid = _create_1x1()
    elif not isinstance(target_cube, iris.cube.Cube):
        emsg = 'Expecting a {!r} instance or "1x1", got {!r}.'
        raise ValueError(emsg.format(iris.cube.Cube.__name__,
                                     type(target_cube)))
        
    schemes = dict(Linear=iris.analysis.Linear(),
                    Nearest=iris.analysis.Nearest(),
                    AreaWeighted=iris.analysis.AreaWeighted()
                   )

    #horizontal regridding
    result_cube = source_cube.regrid(target_grid, schemes[scheme])

    #vertical regridding (if needed)

    print 'source', source_cube
    print 'target',target_grid_cube
    print 'scheme',scheme
    print 'result',result_cube
    print 'done!'

    return result_cube
