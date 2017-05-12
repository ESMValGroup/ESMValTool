"""
This script is a dummy python script template for
diagnostics implemented in python

ensure to define a main() routine in your script as this is
important for the calling routine interface

The docstrings should follow the convention used by the NUMPY project
https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt

"""

#from esmval_lib import ESMValProject

def main(project_info):
    """
    This is the main routine

    Parameters
    ----------
    project_info : dict
        Dicrtionary with project information
    """
    print('Hello, here is the dummy routine from the direct python interface!')

    # create instance of a wrapper that allows easy access to data
    #E = ESMValProject(project_info)

    # get filenames of preprocessed climatological mean files
    #model_filenames = E.get_clim_model_filenames(variable='ta', monthly=True)
    #print model_filenames

    print('Do something here!')
    print('ENDED SUCESSFULLY!!')
    print('')

