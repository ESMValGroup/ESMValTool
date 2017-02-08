import iris

def regrid(source_cube, target_grid_cube, scheme):

    if scheme is None:
        if target_grid_cube is not None:
            raise Exception('target grid must be none if no scheme is given')
        return source_cube

    if target_grid_cube is None:
        raise Exception('target grid must not be none')

    schemes = dict(Linear=iris.analysis.Linear(),
                    Nearest=iris.analysis.Nearest(),
                    AreaWeighted=iris.analysis.AreaWeighted()
                   )

    #horizontal regridding
    result_cube = source_cube.regrid(target_grid_cube, schemes[scheme])

    #vertical regridding (if needed)


    print 'source', source_cube
    print 'target',target_grid_cube
    print 'scheme',scheme
    print 'result',result_cube
    print 'done!'

    return result_cube