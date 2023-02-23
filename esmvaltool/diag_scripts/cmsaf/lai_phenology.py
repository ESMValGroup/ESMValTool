"""
ESMValTool diagnostic for GLOBMAP LAI 4GST.

xxxxxxxxxxxxx add text here
"""

import datetime
import logging
import numpy as np
from scipy.stats import linregress

import iris
import iris.coord_categorisation as icc
#  import iris.quickplot as qplt
import iris.plot as iplt

import cf_units

import matplotlib.pyplot as plt
#  from matplotlib.dates import DateFormatter
import matplotlib.dates as mdates
from matplotlib.patches import Patch
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

from esmvaltool.diag_scripts.shared import (
    #  ProvenanceLogger,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(__name__)

ARGMAX = iris.analysis.Aggregator("argmax",
                                  np.argmax,
                                  units_func=lambda units: 1)

sin_ifunc = iris.analysis.maths.IFunc(np.sin,
                                       lambda cube: cf_units.Unit('1'))
cos_ifunc = iris.analysis.maths.IFunc(np.cos,
                                       lambda cube: cf_units.Unit('1'))
def atan_units_func(cube1, cube2):
    return cf_units.Unit('1')
atan_ifunc = iris.analysis.maths.IFunc(np.arctan2,
                                        atan_units_func)

def mod_units_func(cube1):
    return cf_units.Unit('1')
def mod_func(X):
    return np.mod(X,12)
mod_ifunc = iris.analysis.maths.IFunc(mod_func, mod_units_func)

def mod_12(x):
    return x%12
# can use same single cube units as mod_ifunc
mod_12_ifunc = iris.analysis.maths.IFunc(mod_12, mod_units_func)

# [regression pattern], id number
grow_season_types = {'EVG': [[0, 0, 0, 0], 0],
                     'TGS': [[1, -1, 1, -1], 1],
                     'SGS-S': [[1, 1, -1, -1], 2],
                     'SGS-D': [[-1, -1, 1, 1], 3],
                     }
# EVG is defined differently hence 0 0 0 0 placeholder

fig_fontsizes = {'title': 30,
                 'labels': 24,
                 'ticks': 20,
                 }

# https://davidmathlogic.com/colorblind 'Wong' colour map
# from https://www.nature.com/articles/nmeth.1618
colour_list = ['#000000', '#E69F00', '#56B4E9', '#009E73',
               '#F0E442', '#0072B2', '#D55E00', '#CC79A7'
               ]

# for color legends of the types + undefined
legend_elements = [Patch(facecolor=colour_list[3 + 0], edgecolor='black',
                         label='EVG', alpha=0.4),
                   Patch(facecolor=colour_list[3 + 1], edgecolor='black',
                         label='TGS', alpha=0.4),
                   Patch(facecolor=colour_list[3 + 2], edgecolor='black',
                         label='SGS-S', alpha=0.4),
                   Patch(facecolor=colour_list[3 + 3], edgecolor='black',
                         label='SGS-D', alpha=0.4),
                   Patch(facecolor=colour_list[3 + 4], edgecolor='black',
                         label='Undef', alpha=0.4),
                   ]

def _get_input_cubes(metadata):
    """Load the data files into cubes.

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

    xxxxxxxxxxxxxxxxxx update this to GLOBMAP LAI stuff
    Inputs:
    attributes = xxxxxxx
    ancestor_files = list of data files used by the diagnostic.

    Outputs:
    record = dictionary of provenance records.
    """
    caption = "xxxxxxx do this"

    record = {
        'caption': caption,
        'statistics': ['mean', 'stddev'],  # this will need changing
        'domains': ['reg'],
        'plot_types': ['times'],
        'authors': ['king_robert'],
        #  'references': [],
        'ancestors': ancestor_files
    }

    return record

def calc_ave_month(cube, collapse_coord = 'year'):
    """
    Calculate the average month, respecting the year boundary
    
    Uses the vector modulo average technique
    """
    # Needs month number 0-11 not 1-12
    cube = iris.analysis.maths.subtract(cube, 1, in_place=False)

    print(cube)

    # Convert each month to a unit vector
    data_cos = iris.analysis.maths.apply_ufunc(np.cos, 2*np.pi*cube/12,
                                               in_place=False)
    data_sin = iris.analysis.maths.apply_ufunc(np.sin, 2*np.pi*cube/12,
                                               in_place=False)

    # Find mean of each component
    data_cos_mean = data_cos.collapsed(collapse_coord,iris.analysis.MEAN)
    data_sin_mean = data_sin.collapsed(collapse_coord,iris.analysis.MEAN)

    # Wanted angle is arctan of Y_mean and X_mean
    # use acrtan2 to get sector correct
    data_atan = atan_ifunc(data_sin_mean, data_cos_mean, new_name='mean vector month')
    

    # convert radians to 0-11
    # modulo 12 and %12 mean output is in 0 <= month < 12
    # ie doesnt inc 12 itself 
    data_atan = data_atan*12/np.pi/2

    data_mod = mod_ifunc(data_atan, new_name='mean month')
    
    mean_month = mod_12_ifunc(data_mod, new_name='mean month')

    return mean_month

def calc_phenology(data):
    """
    write this
    """

    month_of_max = {}
    value_of_max = {}

    for dataset in data.keys():
        if dataset == 'GLOBMAP_LAI':
            variable = 'lairpk'
        else:
            variable = 'lai'

        month_of_max[dataset] = iris.cube.CubeList()
        value_of_max[dataset] = iris.cube.CubeList()
        for year in np.arange(1991,1999): # make the whole range...
            
            this_year = data[dataset][variable].extract(iris.Constraint(year=year))

            # Calculate the month of max, and value of max LAI
            this_max = this_year.collapsed('time', ARGMAX)
            this_max_value = this_year.collapsed('time', iris.analysis.MAX)
            this_month = this_year.coord('month_number').points


            value_of_max[dataset].append(this_max_value)

            wanted = this_month[this_max.data]

            lat_coord = this_max.coord('latitude')
            lon_coord = this_max.coord('longitude')

            max_month = iris.cube.Cube(wanted,
                                       long_name='month of max lai',
                                       units=1,
                                       dim_coords_and_dims = ((lat_coord, 0),
                                                              (lon_coord, 1)),
                                       )
            # add year
            year_coord = iris.coords.AuxCoord(year, long_name='year', units=1)
            max_month.add_aux_coord(year_coord)

            month_of_max[dataset].append(max_month)

        month_of_max[dataset] = month_of_max[dataset].merge_cube()
        value_of_max[dataset] = value_of_max[dataset].merge_cube()

        # temp save loc
        iris.save(month_of_max[dataset],
                  '/home/users/robking/plots/lai_month_max_{dataset}.nc')
        iris.save(value_of_max[dataset],
                  f'/home/users/robking/plots/lai_value_max_{dataset}.nc')

    return month_of_max, value_of_max


def plot_lai_months(month_of_max, value_of_max):
    """
    Write this....
    """
    # first calc mean of 'month' and value
    # then plot map
    # this will be at the data's orignal resolutions

    ave_month = {}
    ave_value = {}
    for dataset in month_of_max.keys():
        # assume month_... and value_... have same structure

        this_mean_value = value_of_max[dataset].collapsed('time', iris.analysis.MEAN)
        ave_value[dataset] = this_mean_value

        plt.figure()
        iplt.pcolormesh(this_mean_value,
                        cmap='tab10',
                        vmin=0,
                        vmax=8
                        )

        plt.suptitle(f'Max Value of LAI - {dataset}',
                     fontsize=fig_fontsizes['title']
                     )

        plt.colorbar()

        plt.savefig(f'lai_max_{dataset}.png')
    
        ############
        XXX = calc_ave_month(month_of_max[dataset])
        ave_month[dataset] = XXX

        plt.figure()
        iplt.pcolormesh(XXX,
                        cmap='hsv',
                        )

        plt.suptitle(f'Max Month of LAI - {dataset}',
                     fontsize=fig_fontsizes['title'],
                     )

        plt.colorbar()

        plt.savefig(f'lai_max_month_{dataset}.png')

        YYY = month_of_max[dataset].collapsed('year', iris.analysis.MEAN)
        plt.figure()
        iplt.pcolormesh(YYY,
                        cmap='hsv',
                        )

        plt.suptitle(f'Max Month of LAI - {dataset}',
                     fontsize=fig_fontsizes['title'],
                     )

        plt.colorbar()

        plt.savefig(f'lai_max_month_old_{dataset}.png')

        
        plt.figure()
        iplt.pcolormesh(XXX-YYY,
                        
                        )

        plt.suptitle(f'Diff in Max Month of LAI methods - {dataset}',
                     fontsize=fig_fontsizes['title'],
                     )

        plt.colorbar()

        plt.savefig(f'lai_max_month_diff_{dataset}.png')

        plt.close('all')


    return ave_month, ave_value

def plot_zonal(ave_month, ave_value):

    fig, ax = plt.subplots(1, 2, figsize=(20, 20))

    count = 0
    for dataset in ave_month.keys():
        count += 1
        colour = colour_list[count%8]
        
        months_working =  ave_value[dataset].copy()
        mask_use = np.ma.masked_less(ave_value[dataset].data, 0.1)
        months_working.data = np.ma.array(ave_value[dataset].data,
                                          mask = mask_use.mask)




        this_month = calc_ave_month(months_working, collapse_coord = 'longitude')
        this_value = ave_value[dataset].collapsed('longitude', iris.analysis.MEAN)
        # note month number plot will still jump about dec - jan -dec 

        plt.sca(ax[0])
        plt.plot(this_month.data,
                 this_month.coord('latitude').points,
                 label=dataset,
                 color=colour,
                 linewidth=3
                 )

        plt.sca(ax[1])
        plt.plot(this_value.data,
                 this_value.coord('latitude').points,
                 label=dataset,
                 color=colour,
                 linewidth=3
                 )


    ax[0].set_title('Ave. Month of Max', fontsize=fig_fontsizes['labels'])
    ax[1].set_title('Ave. Value of Max', fontsize=fig_fontsizes['labels'])

    
    ax[0].set_xticks(np.arange(0,12),
                     ['Jan','Feb','Mar','Apr','May','Jun',
                      'Jul','Aug','Sep','Oct','Nov','Dec'],
                     fontsize=fig_fontsizes['ticks']
                     )
    ax[0].set_xlim((-0.5,12.5))

    ax[1].set_xticks(np.arange(0,9),np.arange(0,9),
                     fontsize=fig_fontsizes['ticks']
                 )
    ax[1].set_xlim((0,8))

    ax[0].set_yticks(np.arange(-60,76,5),np.arange(-60,76,5),
                     fontsize=fig_fontsizes['ticks']
                 )

    ax[1].set_yticks(np.arange(-60,76,5),
                     ['' for item in np.arange(-60,76,5)],
                     fontsize=fig_fontsizes['ticks']
                 )

    ax[0].set_ylim((-60,75))
    ax[1].set_ylim((-60,75))

    ax[0].grid()
    ax[1].grid()

    ax[0].set_xlabel('Month', fontsize=fig_fontsizes['labels'])
    ax[1].set_xlabel('LAI', fontsize=fig_fontsizes['labels'])


    plt.legend()

    plt.savefig(f'lai_max_zonal.png')
    
def _diagnostic(config):
    """Perform the .......

    Parameters
    ----------
    config: dict
        the preprocessor nested dictionary holding
        all the needed information.

    Returns
    -------
    xxxxxxxxxxxxxxxxxxxxxxxx
    """
    # this loading function is based on the hydrology diagnostic
    input_metadata = config['input_data'].values()

    loaded_data = {}
    ancestor_list = []  # what is this actually used for?
    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        cubes, ancestors = _get_input_cubes(metadata)
        loaded_data[dataset] = cubes
        # ancestor_list.append(ancestors['ts'][0])

    print(loaded_data)
    # OBS loaded_data['GLOBMAP_LAI']['lairpk']

    # note LST data needed some time coords changing
    # CMIP data had 360 day calendar, CCI data has 365 day calendar

    for dataset in loaded_data.keys():
        for this_var in loaded_data[dataset]:
            tcoord = loaded_data[dataset][this_var].coord('time')
            if tcoord.units.calendar != 'gregorian':
                loaded_data[dataset][this_var].coord('time').units = \
                cf_units.Unit(tcoord.units.origin,
                calendar='gregorian')

            icc.add_month_number(loaded_data[dataset][this_var], 'time')
            icc.add_year(loaded_data[dataset][this_var], 'time')

            if dataset == 'GLOBMAP_LAI':
                loaded_data[dataset][this_var] /= 100.

    month_of_max, value_of_max = calc_phenology(loaded_data)

    ave_month, ave_value = plot_lai_months(month_of_max, value_of_max)

    plot_zonal(ave_month, ave_value)

    # Provenance

    # Get this information form the data cubes
    # add this back in using LST diagnostic as template


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        _diagnostic(config)
