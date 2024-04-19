"""
ESMValTool diagnostic for ESA CCI LST V3 data - Uncertainity Propagation
"""

# Expected keys for OBS LST data in loaded_data:
# ts_day
# ts_night
# lst_unc_loc_sfc_day
# lst_unc_loc_sfc_night
# lst_unc_loc_atm_day
# lst_unc_loc_atm_night
# lst_unc_sys_day
# lst_unc_sys_night
# lst_unc_ran_day
# lst_unc_ran_night

import logging
import iris
import matplotlib.pyplot as plt
import numpy as np

import iris.plot as iplt

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(__name__)


def _get_input_cubes(metadata):
    """Load the data files into cubes.

    Based on the hydrology diagnostic originally.

    Inputs:
    metadata = List of dictionaries made from the preprocessor config

    Outputs:
    inputs = Dictionary of cubes
    ancestors = Dictionary of filename information
    """
    inputs = {}
    ancestors = {}
    for attributes in metadata:

        short_name = attributes['short_name']
        filename = attributes['filename']
        logger.info("Loading variable %s", short_name)
        cube = iris.load_cube(filename)
        cube.attributes.clear()
        inputs[short_name] = cube
        ancestors[short_name] = [filename]

    return inputs, ancestors


def _get_provenance_record(attributes, ancestor_files):
    """Create the provenance record dictionary.

    Inputs:
    attributes = dictionary of ensembles/models used, the region bounds
                 and years of data used.
    ancestor_files = list of data files used by the diagnostic.

    Outputs:
    record = dictionary of provenance records.
    """
    caption = "This needs a new caption.".format(**attributes)

    record = {
        'caption': caption,
        'statistics': ['mean', 'stddev'],
        'domains': ['reg'],
        'plot_types': ['times'],
        'authors': ['king_robert'],
        # 'references': [],
        'ancestors': ancestor_files
    }

    return record


def _diagnostic(config):
    """Perform the control for the ESA CCI LST diagnostic.

    Parameters
    ----------
    config: dict
        the preprocessor nested dictionary holding
        all the needed information.

    Returns
    -------

    """
    # this loading function is based on the hydrology diagnostic
    input_metadata = config['input_data'].values()
    loaded_data = {}
    ancestor_list = []

    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        cubes, ancestors = _get_input_cubes(metadata)
        loaded_data[dataset] = cubes

    # Methodolgy:
    # calcualte for day and night seperately
    # calls to propagation equations for each compoent
    # ts eq 1 eq_arithmetic_mean
    # lst_unc_loc_atm eq 7 eq_weighted_sqrt_mean
    # lst_unc_sys eq 5 = eq 1 eq_arithmetic_mean, no spatial propagation here
    # lst_unc_loc_sfc DEPENDS ON SENSOR
    # see Table 1 and then table 4
    # MODIS Aqua and Terra = GSW eq 6 = eq 1 eq_arithmetic_mean.
    # lst_unc_ran eq 4 eq_propagate_random_with_sampling
    #
    # Then combine to get a total day and night vale
    # total uncert  eq 9 eq_sum_in_quadrature
    # a total all day uncertainty can then be obtained using
    # eq 9 eq_sum_in_quadrature
    # and a mean eq 1 eq_arithmetic_mean to get an 'average' lst for the day

    lst_unc_variables = ['lst_unc_loc_atm', 'lst_unc_sys',
                         'lst_unc_loc_sfc', 'lst_unc_ran']

    # This will be a dictionary of variables and cubes
    # of their propagated values
    propagated_values = {}

    # These define the total number of points in the data
    lat_len = len(loaded_data['ESACCI-LST']['ts_day'].coord('latitude').points)
    lon_len = len(loaded_data['ESACCI-LST']['ts_day'].coord('longitude').points)

    # n_fill and n_use ad dictionaries with keys 'day' and 'night'
    # the item for each key is an array
    # These give the cloud/don't use pixel numbers,
    # and the number of useable pixels for each timestep
    n_fill = {}
    n_use = {}
    n_fill['day'] = np.array(
        [np.sum(item)
         for item in loaded_data['ESACCI-LST']['ts_day'].data.mask])
    n_fill['night'] = np.array(
        [np.sum(item)
         for item in loaded_data['ESACCI-LST']['ts_night'].data.mask])

    n_use['day'] = (lat_len*lon_len) - n_fill['day']
    n_use['night'] = (lat_len*lon_len)

    # This loop call the propagation equations, once for 'day' and 'night'
    for time in ['day', 'night']:

        propagated_values[f'ts_{time}'] = \
        eq_arithmetic_mean(loaded_data['ESACCI-LST'][f'ts_{time}'])

        propagated_values[f'lst_unc_loc_atm_{time}'] = \
        eq_weighted_sqrt_mean(loaded_data['ESACCI-LST'][f'lst_unc_loc_atm_{time}'],
        n_use[f'{time}'])
        # no spatial propagation of the systamatic uncertainity
        propagated_values[f'lst_unc_sys_{time}'] = \
        loaded_data['ESACCI-LST'][f'lst_unc_sys_{time}']

        propagated_values[f'lst_unc_loc_sfc_{time}'] = \
        eq_arithmetic_mean(loaded_data['ESACCI-LST'][f'lst_unc_loc_sfc_{time}'])

        propagated_values[f'lst_unc_ran_{time}'], \
        propagated_values[f'lst_sampling_{time}'] = \
        eq_propagate_random_with_sampling(loaded_data['ESACCI-LST'][f'lst_unc_ran_{time}'],
                                          loaded_data['ESACCI-LST'][f'ts_{time}'],
                                          n_use[f'{time}'], n_fill[f'{time}'])

    # Combines all uncertainty types to get total uncertainty
    # for 'day' and 'night'
    for time in ['day', 'night']:
        time_cubelist = iris.cube.CubeList([propagated_values[f'{variable}_{time}']
                                            for variable in lst_unc_variables])
        propagated_values[f'lst_total_unc_{time}'] = \
                                            eq_sum_in_quadrature(time_cubelist)

    test_plot(propagated_values)


def test_plot(propagated_values):
    """This is a very simple plot to just test the method
    """

    for time in ['day', 'night']:
        # one plot for day and night seperately

        plt.figure(figsize=(12, 10))

        plt.subplot(211)
        iplt.plot(propagated_values[f'ts_{time}'], c='b')

        plt.subplot(212)
        iplt.plot(propagated_values[f'lst_unc_loc_atm_{time}'], c='r',
                  label=f'lst_unc_loc_atm_{time}')
        iplt.plot(propagated_values[f'lst_unc_loc_sfc_{time}'], c='g',
                  label=f'lst_unc_loc_sfc_{time}')
        iplt.plot(propagated_values[f'lst_unc_sys_{time}'], c='b',
                  label=f'lst_unc_sys_{time}')
        iplt.plot(propagated_values[f'lst_unc_ran_{time}'], c='m',
                  label=f'lst_unc_ran_{time}')

        iplt.plot(propagated_values[f'lst_sampling_{time}'], c='c',
                  label=f'lst_sampling_{time}')
        iplt.plot(propagated_values[f'lst_total_unc_{time}'], c='k',
                  label=f'lst_total_unc_{time}')

        plt.legend()

        plt.savefig(f'test_{time}.png')

# These are the propagation equations


def eq_propagate_random_with_sampling(cube_unc_ran, cube_ts, n_fill, n_use):
    """Propagate radom uncertatinty using the sampling uncertainty
    ATBD eq 4

    the sum in quadrature of the arithmetic mean of lst_unc_ran and
    the sampling uncertainty

    Sampling uncertainty is
    n_cloudy * Variance of LST / n_total-1

    see ATBD section 3.13.3

    Inputs:
    cube_unc_ran: The cube with the lst_unc_ran day/night data
    cube_ts:      The lst use for day/night as appropriate
    """

    # total number of pixed
    n_total = n_fill + n_use

    # the mean of the random uncertainty
    unc_ran_mean = eq_arithmetic_mean(cube_unc_ran)

    # calculate the sampling error
    # variance of the lst * n_fill/n_total-1
    lst_variance = cube_ts.collapsed(['latitude', 'longitude'],
                                     iris.analysis.VARIANCE)
    factor = n_fill/(n_total - 1)
    unc_sampling = iris.analysis.maths.multiply(lst_variance,
                                                factor)

    # apply the ATBD equation
    # note the square of random uncertainty is needed
    output = eq_sum_in_quadrature(iris.cube.CubeList([unc_ran_mean**2,
                                                      unc_sampling]))
    output = iris.analysis.maths.exponentiate(output, 0.5)

    return output, unc_sampling


def eq_arithmetic_mean(cube):
    """Arithmetic mean of cube, across latitude and longitude
    ATBD eq 1
    """

    out_cube = cube.collapsed(['latitude', 'longitude'], iris.analysis.MEAN)

    return out_cube


def eq_sum_in_quadrature(cubelist):
    """Sum in quadrature
    ATBD eq 9

    Input:
    cubelist : A cubelist of 1D cubes
    """

    # dont want to in-place replace the input
    newlist = cubelist.copy()
    for cube in newlist:
        iris.analysis.maths.exponentiate(cube, 2, in_place=True)

    cubes_sum = 0
    for cube in newlist:
        cubes_sum += cube

    output = iris.analysis.maths.exponentiate(cube, 0.5,
                                              in_place=False)

    return output


def eq_weighted_sqrt_mean(cube, n_use):
    """Mean with square root of n factor
    ATBD eq 7

    Inputs:
    cube:
    n_use: the number of useable pixels
    NEED TO IMPLIMENT A CHECK ON MASKS BEING THE SAME?
    """
    output = iris.analysis.maths.multiply(cube.collapsed(['latitude',
                                                          'longitude'],
                                                         iris.analysis.MEAN),
                                          1/np.sqrt(n_use))
    return output


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        _diagnostic(config)
