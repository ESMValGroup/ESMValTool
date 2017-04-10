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
import numpy as np
import matplotlib.pyplot as plt
import os

def main(project_info):
    print('Hello, here is the dummy routine from the direct python interface!')

    # create instance of a wrapper that allows easy access to data
    E = ESMValProject(project_info)

    # get filenames of preprocessed climatological mean files
    model_filenames = E.get_clim_model_filenames(variable='ta', monthly=True)
    print model_filenames

    plot_dir = project_info['GLOBAL']['plot_dir']
    if plot_dir[-1] != os.sep:
        plot_dir += os.sep

    wrk_dir = project_info['GLOBAL']['wrk_dir']
    if wrk_dir[-1] != os.sep:
        wrk_dir += os.sep

    f = plt.figure()
    ax = f.add_subplot(111)
    ax.plot(np.random.random(100))
    f.savefig(plot_dir + 'test_image.png')

    o=open(plot_dir + 'python_test_out.tab', 'w')
    o.write('WELCOME TO THE EXTERNAL PYTHON SCRIPT\n')
    o.write('HELLO WORLD')
    o.close()


    # create empty file
    o=open(plot_dir + 'empty.txt', 'w')
    o.close()

    # write something to the logfile
    o = open(wrk_dir + 'refs-acknows_dummy_python_new.log', 'w')
    o.write('Luke Skywalker was here!')
    o.close()


    #~ print 'Project info keys:'
    #~ for k in project_info['GLOBAL'].keys():
        #~ print k

    print('Do something here!')
    print('ENDED SUCESSFULLY!!')
    print('')

