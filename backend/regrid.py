import iris

def regrid(source_cube, target_grid_cube, scheme):

    schemes = dict(Linear=iris.analysis.Linear(),
                    Nearest=iris.analysis.Nearest(),
                    AreaWeighted=iris.analysis.AreaWeighted()
                   )

    result_cube = source_cube.regrid(target_grid_cube, schemes[scheme])

    print 'source', source_cube
    print 'target',target_grid_cube
    print 'scheme',scheme
    print 'result',result_cube
    print 'done!'

    return result_cube