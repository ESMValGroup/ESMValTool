"""
*********************************************************************
 dummy_python.py
*********************************************************************
 PYTHON script
 dummy_python.py
 alexander.loew@mpimet.mpg.de, February 2014
*********************************************************************
 This script is a dummy python script template for
 diagnostics implemented in python

 ensure to define a main() routine in your script as this is
 important for the calling routine interface

*********************************************************************
"""

from esmval_lib import ESMValProject

import iris
import os


def regridder(source_cube, target_grid_cube, scheme):

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


def main(project_info):
    print('Hello, here is the regrid routine from the direct python interface!')

    # create instance of a wrapper that allows easy access to data
    E = ESMValProject(project_info)

    # get filenames of preprocessed climatological mean files
    model_filenames = E.get_clim_model_filenames(variable='ta', monthly=True)
    print model_filenames


    model_files = ESMValProject(project_info).get_clim_model_filenames(variable='ta')

    preprocess = project_info['DIAGNOSTICS'].diagnostics[0].preprocess


    print preprocess

    print iris.__version__

    #assumes target grid is one of the models
    target_grid_model = preprocess['regrid_target_grid']
    scheme = preprocess['regrid_scheme']

    reference_model_filename = model_files[target_grid_model]

    target_grid_cube = iris.load_cube(reference_model_filename)
    out_dir = os.path.join(E.get_work_dir(), 'regridded')
    if not os.path.exists(out_dir):
        os.mkdirs(out_dir)

    for model,file in model_files.items():
        #if not equal to target grid dimensions with some tolerance
        if model != target_grid_model:
            source_cube = iris.load_cube(file)
            result_cube = regridder(source_cube, target_grid_cube, scheme)

            print result_cube

            result_file_name = os.path.join(out_dir, os.path.basename(file))

            iris.save(result_cube, result_file_name)






    print('Do something here!')
    print('ENDED SUCESSFULLY!!')
    print('')

