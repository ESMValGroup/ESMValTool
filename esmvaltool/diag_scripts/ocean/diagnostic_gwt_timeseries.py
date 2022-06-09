"""
GWT Time series diagnostics
=======================

Diagnostic to produce figures of the time development of the GWT field from
cubes. These plost show time on the x-axis and cube value (ie temperature) on
the y-axis.

Two types of plots are produced: individual model timeseries plots and
multi model time series plots. The inidivual plots show the results from a
single cube, even if this is a mutli-model mean made by the _multimodel.py
preproccessor. The multi model time series plots show several models
on the same axes, where each model is represented by a different line colour.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has a time component, no depth component, and no
latitude or longitude coordinates.

An approproate preprocessor for a 3D+time field would be::

  preprocessors:
    prep_timeseries_1:# For Global Volume Averaged
      volume_statistics:
        operator: mean


An approproate preprocessor for a 3D+time field at the surface would be::

    prep_timeseries_2: # For Global surface Averaged
      extract_levels:
        levels:  [0., ]
        scheme: linear_extrap
      area_statistics:
        operator: mean


An approproate preprocessor for a 2D+time field would be::

    prep_timeseries_2: # For Global surface Averaged
      area_statistics:
        operator: mean


This tool is part of the ocean diagnostic tools package in the ESMValTool.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk
"""

import logging
import os
import datetime
import iris
import matplotlib.pyplot as plt
logging.getLogger('matplotlib.font_manager').disabled = True
from matplotlib import gridspec
from matplotlib.colors import to_rgba
import numpy as np
from itertools import product
import cf_units
import glob
import shelve
import csv
import scipy

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvalcore.preprocessor._time import regrid_time
from esmvalcore.preprocessor._multimodel import multi_model_statistics, _multicube_statistics
from esmvalcore.preprocessor._io import _fix_cube_attributes

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


# List of ensemble members to skip.
# data_dict_skips = {
# #'    'CanESM5':['r11i1p2f1', 'r12i1p1f1', 'r12i1p2f1'], # hist data doesn't exist
# #    'CanESM5-CanOE': ['r1i1p2f1', 'r2i1p2f1', 'r3i1p2f1',], # hist data doesn't exist
# #    'ACCESS-ESM1-5': ['r2i1p1f1', 'r3i1p1f1'], #  # hist data doesn't exist
# #    'CESM2-WACCM':['r2i1p1f1', 'r3i1p1f1'], # no hist data.
#     'CESM2-WACCM-FV2': ['r1i1p1f1', 'r2i1p1f1', 'r3i1p1f1'], # no SSP data at all.
#     'CESM2-FV2': ['r1i1p1f1', 'r2i1p1f1', 'r3i1p1f1'], # no SSP data at all.
#
#
#     'MPI-ESM-1-2-HAM': ['r1i1p1f1', 'r2i1p1f1', ], # no hist data in r2i1p1f1, both only run unil 2050-ish.
# #    'MIROC-ES2L': ['r1i1p1f2', 'r2i1p1f2'],  # no hist data.
# #    'MPI-ESM1-2-LR': ['r2i1p1f1', 'r3i1p1f1'], # no hist data.
#     'NorCPM1': ['r1i1p1f1', ], #fgco2 off by 1e-10, no SSP data, weird TLS, weird NBP.
#     'NorESM2-LM': ['r2i1p1f1'], # missing historical run. #years 2015-2020, and sometimes after 2055 too.
#     }
"""
 mod_exp_ens_skips = {
    # Access existed but wasn't included in the original run:
    ('ACCESS-ESM1-5', 'ssp126', 'r1i1p1f1'): True, # hist missing in early run/
    ('ACCESS-ESM1-5', 'ssp126', 'r2i1p1f1'): True, # hist missing in early run/
    ('ACCESS-ESM1-5', 'ssp126', 'r3i1p1f1'): True, # hist missing in early run/

    ('ACCESS-ESM1-5', 'ssp245', 'r1i1p1f1'): True, # hist missing in early run/
    ('ACCESS-ESM1-5', 'ssp245', 'r2i1p1f1'): True, # hist missing in early run/
    ('ACCESS-ESM1-5', 'ssp245', 'r3i1p1f1'): True, # hist missing in early run/

    ('ACCESS-ESM1-5', 'ssp370', 'r1i1p1f1'): True, # hist missing in early run/
    ('ACCESS-ESM1-5', 'ssp370', 'r2i1p1f1'): True, # hist missing in early run/
    ('ACCESS-ESM1-5', 'ssp370', 'r3i1p1f1'): True, # hist missing in early run/

    ('CESM2-WACCM', 'ssp245', 'r2i1p1f1'): True, # Delete later.
    ('CESM2-WACCM', 'ssp245', 'r3i1p1f1'): True, # Delete later.
    ('CESM2-WACCM', 'ssp370', 'r1i1p1f1'): True, # Delete later.
    ('CESM2-WACCM', 'ssp370', 'r2i1p1f1'): True, # Delete later.
    ('CESM2-WACCM', 'ssp370', 'r3i1p1f1'): True, # Delete later.
    ('CNRM-ESM2-1', 'ssp370', '*'): True, #
    ('CNRM-ESM2-1', 'ssp245', '*'): True, #
    ('CNRM-ESM2-1', 'ssp126', '*'): True, #
    ('CanESM5-CanOE', 'ssp126', '*'): True, # Delete later.
    ('CanESM5-CanOE', 'ssp245', '*'): True, # Delete later.
    ('CanESM5-CanOE', 'ssp370', '*'): True, # Delete later.

    ('IPSL-CM6A-LR', 'ssp245', 'r10i1p1f1'):True,
    ('IPSL-CM6A-LR', 'ssp370', 'r10i1p1f1'): True,
    ('IPSL-CM6A-LR', 'ssp245', 'r11i1p1f1'):True,
    ('IPSL-CM6A-LR', 'ssp370', 'r11i1p1f1'): True,

    ('MIROC-ES2L', 'ssp119', 'r1i1p1f2'): True,
    ('MIROC-ES2L', 'ssp126', 'r1i1p1f2'): True,
    ('MIROC-ES2L', 'ssp245', 'r1i1p1f2'): True,
    ('MIROC-ES2L', 'ssp370', 'r1i1p1f2'): True,
    ('MIROC-ES2L', 'ssp585', 'r1i1p1f2'): True,
    ('MIROC-ES2L', 'ssp119', 'r2i1p1f2'): True,

    ('MPI-ESM1-2-LR', 'ssp126', 'r2i1p1f1'): True,
    ('MPI-ESM1-2-LR', 'ssp245', 'r2i1p1f1'): True,
    ('MPI-ESM1-2-LR', 'ssp370', 'r2i1p1f1'): True,
    ('MPI-ESM1-2-LR', 'ssp126', 'r3i1p1f1'): True,
    ('MPI-ESM1-2-LR', 'ssp245', 'r3i1p1f1'): True,
    ('MPI-ESM1-2-LR', 'ssp370', 'r3i1p1f1'): True,





    ('CanESM5', 'ssp119', 'r11i1p2f1'): True, # no hist run
    ('CanESM5', 'ssp119', 'r12i1p2f1'): True, # no hist run
    ('CanESM5', 'ssp126', 'r11i1p2f1'): True, # no hist run
    ('CanESM5', 'ssp126', 'r12i1p2f1'): True, # no hist run
    ('CanESM5', 'ssp245', 'r11i1p2f1'): True, # no hist run
    ('CanESM5', 'ssp245', 'r12i1p2f1'): True, # no hist run
    ('CanESM5', 'ssp370', 'r11i1p2f1'): True, # no hist run
    ('CanESM5', 'ssp370', 'r12i1p2f1'): True, # no hist run
    ('CanESM5', 'ssp585', 'r11i1p2f1'): True, # no hist run
    ('CanESM5', 'ssp585', 'r12i1p2f1'): True, # no hist run

#    ('IPSL-CM6A-LR', 'nbpgt', 'historical-ssp245', 'r10i1p1f1'):True,
#    ('IPSL-CM6A-LR', 'nbpgt', 'historical-ssp245', 'r10i1p1f1'): True,

    ('NorESM2-LM', 'ssp245', 'r2i1p1f1',): True, # no histroical run.

"""
mod_exp_ens_skips = {
#    ('ACCESS-ESM1-5', 'historical-ssp126', 'r2i1p1f1') : True, # No historical run on jasmin
#    ('ACCESS-ESM1-5', 'historical-ssp245', 'r2i1p1f1') : True, # No historical run on jasmin
#    ('ACCESS-ESM1-5', 'historical-ssp245', 'r3i1p1f1') : True, # No historical run on jasmin
#    ('ACCESS-ESM1-5', 'historical-ssp370', 'r2i1p1f1') : True, # No historical run on jasmin
#    ('CESM2-WACCM', 'historical-ssp245', 'r2i1p1f1') : True, # No historical run on jasmin
#    ('CESM2-WACCM', 'historical-ssp245', 'r3i1p1f1') : True, # No historical run on jasmin
#    ('CNRM-ESM2-1', 'historical-ssp370', 'r1i1p1f2') : True, # No historical run on jasmin
#    ('CNRM-ESM2-1', 'historical-ssp370', 'r3i1p1f2') : True, # No historical run on jasmin
#    ('CNRM-ESM2-1', 'historical-ssp370', 'r4i1p1f2') : True, # No historical run on jasmin
#    ('CNRM-ESM2-1', 'historical-ssp370', 'r5i1p1f2') : True, # No historical run on jasmin
#    ('CanESM5-CanOE', 'historical-ssp126', 'r1i1p2f1') : True, # No historical run on jasmin
#    ('CanESM5-CanOE', 'historical-ssp126', 'r2i1p2f1') : True, # No historical run on jasmin
#    ('CanESM5-CanOE', 'historical-ssp245', 'r1i1p2f1') : True, # No historical run on jasmin
#    ('CanESM5-CanOE', 'historical-ssp245', 'r2i1p2f1') : True, # No historical run on jasmin
#    ('CanESM5-CanOE', 'historical-ssp370', 'r1i1p2f1') : True, # No historical run on jasmin
#    ('CanESM5-CanOE', 'historical-ssp370', 'r2i1p2f1') : True, # No historical run on jasmin
#    ('IPSL-CM6A-LR', 'historical-ssp245', 'r10i1p1f1') : True, # No historical run on jasmin
#    ('IPSL-CM6A-LR', 'historical-ssp245', 'r11i1p1f1') : True, # No historical run on jasmin
#    ('IPSL-CM6A-LR', 'historical-ssp370', 'r10i1p1f1') : True, # No historical run on jasmin
#    ('MIROC-ES2L', 'historical-ssp119', 'r1i1p1f2') : True, # No historical run on jasmin
#    ('MIROC-ES2L', 'historical-ssp119', 'r2i1p1f2') : True, # No historical run on jasmin
#    ('MIROC-ES2L', 'historical-ssp245', 'r1i1p1f2') : True, # No historical run on jasmin
#    ('MIROC-ES2L', 'historical-ssp370', 'r1i1p1f2') : True, # No historical run on jasmin
#    ('MPI-ESM1-2-LR', 'historical-ssp126', 'r2i1p1f1') : True, # No historical run on jasmin
#    ('MPI-ESM1-2-LR', 'historical-ssp245', 'r2i1p1f1') : True, # No historical run on jasmin
#    ('MPI-ESM1-2-LR', 'historical-ssp245', 'r3i1p1f1') : True, # No historical run on jasmin
#    ('MPI-ESM1-2-LR', 'historical-ssp370', 'r2i1p1f1') : True, # No historical run on jasmin

    ('IPSL-CM5A2-INCA', 'historical', 'r1i1p1f1',): True, # no SSP runs.
    ('IPSL-CM5A2-INCA', 'piControl', 'r1i1p1f1'): True,

    ('CESM2-WACCM-FV2', 'historical', 'r1i1p1f1',): True, # no SSP runs.
    ('CESM2-FV2', 'historical', 'r1i1p1f1',): True, # no SSP runs.

   # ('CESM2-WACCM', 'ssp370', 'r2i1p1f1'): True, # SSP run ends at 2054.
    ('CESM2-WACCM', 'ssp370', 'r2i1p1f1'): True, # SSP run ends at 2054.
    ('CESM2-WACCM', 'ssp370', 'r3i1p1f1'): True, # SSP run ends at 2054.

    ('MPI-ESM-1-2-HAM', 'historical', 'r1i1p1f1',): True, # no SSP runs.
    ('MPI-ESM-1-2-HAM', 'historical', 'r2i1p1f1',): True, # no SSP runs.
    ('MPI-ESM-1-2-HAM', 'ssp370', 'r1i1p1f1',): True, # SSP run ends at 2054.
    ('MPI-ESM-1-2-HAM', 'ssp370', 'r2i1p1f1',): True, # SSP run ends at 2054.

    ('NorCPM1', 'historical', 'r1i1p1f1',): True, # no SSP runs.

    ('NorESM2-LM', 'ssp370', 'r2i1p1f1',): True, # SSP run ends at 2054.

    ('*', 'ssp534-over', ''): True, # no LUE data available.
    }

def extend_mod_exp_ens_skips(mod_exp_ens_skips):
    new_dict = {}
    for (mod, exp, ens), boole in mod_exp_ens_skips.items():
        new_exp = 'historical-'+exp
        new_dict[(mod, new_exp, ens)] = True
        new_exp2 = exp.replace('historical-', '')
        new_dict[(mod, new_exp2, ens)] = True

        if ens == '*':
            for r, p, f in product(range(1,20), range(1,5), range(1,5)):
                new_ens = ''.join(['r', str(r),'i1p', str(p),'f', str(f)])
                new_dict[(mod, exp, new_ens)] = True
                new_dict[(mod, new_exp, new_ens)] = True
                new_dict[(mod, new_exp2, new_ens)] = True


    mod_exp_ens_skips.update(new_dict)
    return mod_exp_ens_skips
mod_exp_ens_skips = extend_mod_exp_ens_skips(mod_exp_ens_skips)

# For models (like UKESM), where the hist and ssp have different ensemble ids:
# This one exists, because in this one run, I accidentatlly asked for two ensembmes in the historiacal,. whoops.
data_dict_linked_hist= {# model: {ssp:hist'
    'UKESM1-0-LL': {#'r5i1p1f2':'r5i1p1f3','r6i1p1f2':'r6i1p1f3','r7i1p1f2':'r7i1p1f3',
                    #'r5i1p1f3':'r5i1p1f3','r6i1p1f3':'r6i1p1f3','r7i1p1f3':'r7i1p1f3',
                    ('r5i1p1f3', 'r5i1p1f2'): 'r5i1p1f3', #'r5i1p1f3'),
                    ('r6i1p1f3', 'r6i1p1f2'): 'r6i1p1f3', #'r6i1p1f3'),
                    ('r7i1p1f3', 'r7i1p1f2'): 'r7i1p1f3', #'r7i1p1f3'),
                    'r5i1p1f3_r5i1p1f2': 'r5i1p1f3', #'r5i1p1f3'),
                    'r6i1p1f3_r6i1p1f2': 'r6i1p1f3', #'r6i1p1f3'),
                    'r7i1p1f3_r7i1p1f2': 'r7i1p1f3', #'r7i1p1f3'),

    }
}
data_dict_linked_ens = {# model: {ssp:hist'
    'UKESM1-0-LL': {#'r5i1p1f2':'r5i1p1f3','r6i1p1f2':'r6i1p1f3','r7i1p1f2':'r7i1p1f3',
                    #'r5i1p1f3':'r5i1p1f3','r6i1p1f3':'r6i1p1f3','r7i1p1f3':'r7i1p1f3',
                    ('r5i1p1f3', 'r5i1p1f2'): 'r5i1p1f3',
                    ('r6i1p1f3', 'r6i1p1f2'): 'r6i1p1f3',
                    ('r7i1p1f3', 'r7i1p1f2'): 'r7i1p1f3',

                    'r5i1p1f3_r5i1p1f2': 'r5i1p1f3', #'r5i1p1f3'),
                    'r6i1p1f3_r6i1p1f2': 'r6i1p1f3', #'r6i1p1f3'),
                    'r7i1p1f3_r7i1p1f2': 'r7i1p1f3', #'r7i1p1f3'),
                    #('r5i1p1f2', 'r5i1p1f3'): 'r5i1p1f3',
                    #('r6i1p1f2', 'r6i1p1f3'): 'r6i1p1f3',
                    #('r7i1p1f2', 'r7i1p1f3'): 'r7i1p1f3',
                   },
#    'CanESM5': {'r11i1p2f1': 'r11i1p1f1', 'r12i1p2f1': 'r12i1p1f1'},
    }

exp_colours = {'historical':'black',
               'ssp119':'green',
               'ssp126':'dodgerblue',
               'ssp245':'darkorange',
               'ssp370':'purple',
               'ssp585': 'red',
               'ssp434':'magenta',
               'ssp534-over':'orange',}

exp_colours_fill = {'historical':'black',
               'ssp119':'green',
               'ssp126':'blue',
               'ssp245':'orange',
               'ssp370':'purple',
               'ssp585': 'darkred',}

exp_colours_alpha = {'historical': 0.5,
               'ssp119': 0.6,
               'ssp126': 0.3,
               'ssp245': 0.3,
               'ssp370': 0.3,
               'ssp585': 0.6, }


exp_colours_dark = {'historical':'black',
               'ssp119':'darkgreen',
               'ssp126':'royalblue',
               'ssp245':'saddlebrown',
               'ssp370':'indigo',
               'ssp585': 'darkred',
               'ssp434':'magenta',
               'ssp534-over':'orange',}



cmap_ssp_land = {
        #'historical':'black',
        'SSP119': 'darkgreen',
        'SSP126': 'blue',
        'SSP245': 'saddlebrown',
        'SSP370': 'indigo',
        'SSP585': 'darkred',
        }
cmap_ssp_ocean = {
        'SSP119': 'green',
        'SSP126': 'royalblue',
        'SSP245': 'darkorange',
        'SSP370': 'purple',
        'SSP585': 'red',
    }
cmap_ssp_air = {
        'SSP119': 'lightgreen',
        'SSP126': 'dodgerblue',
        'SSP245': 'sandybrown',
        'SSP370': 'orchid',
        'SSP585': 'lightcoral',
    }

def sspify(ssp):
    sspdict = {
        'historical': 'Historical',
        'SSP119': 'SSP1-1.9',
        'SSP126': 'SSP1-2.6',
        'SSP245': 'SSP2-4.5',
        'SSP370': 'SSP3-7.0',
        'SSP585': 'SSP5-8.5',
        'ssp119': 'SSP1-1.9',
        'ssp126': 'SSP1-2.6',
        'ssp245': 'SSP2-4.5',
        'ssp370': 'SSP3-7.0',
        'ssp585': 'SSP5-8.5',
    }
    return sspdict[ssp]


def make_mean_of_dict_list(dict_list, short_name, metric='mean'):
    """
    Takes the mean of a list of cubes (not an iris.cube.CubeList).
    Assumes all the cubes are the same shape.
    Assumes annual 1D data.
    """
    year_counts = {}
    full_times = {}
    # perform checks:
    for i, ddict in enumerate(dict_list):
        years = ddict['time']

        count = len(years)
        if year_counts.get(count, False):
            year_counts[count].append(i)
        else:
            year_counts[count]= [i, ]

        for year in years:
            try:
                full_times[year] += 1
            except:
                full_times[year] = 1

    if len(year_counts.keys())>1:
        print('make_mean_of_dict_list: ERROR: more than one dataset length:', year_counts)
        print('make_mean_of_dict_list: ERROR: list of years:', full_times)
        assert 0

    print('make_mean_of_dict_list: INFO:,', metric,' everything is right size:', year_counts)
    output_datas = {}
    dict_mean=dict_list[0]
    for i, ddict in enumerate(dict_list):
        years = ddict['time']

        for yr, year in enumerate(years):
            d = ddict[short_name][yr]
            if output_datas.get(year, False):
                output_datas[year].append(d)
            else:
                 output_datas[year] = [d, ]
    datas = []
    years = []
    for yr in sorted(output_datas.keys()):
        if metric=='mean':
            print('calculating average', yr, np.sum(output_datas[yr])/float(len(output_datas[yr])),':',output_datas[yr])
            datas.append(np.sum(output_datas[yr])/float(len(output_datas[yr])))
        elif metric=='min':
            print('calculating minimum:', yr, np.min(output_datas[yr]))
            datas.append(np.min(output_datas[yr]))
        elif metric=='max':
            print('calculating maximum:', yr, np.max(output_datas[yr]))
            datas.append(np.max(output_datas[yr]))
        else:
            assert 0
        years.append(yr)

    dict_mean[short_name] =  np.array(datas)
    dict_mean['time'] = np.array(years)
    return dict_mean



def make_mean_of_cube_list_iris(cube_list, metric='mean'):
    assert 0
    #return multi_model_statistics(cube_list, span='full', statistics=['mean',], keep_input_datasets=False)
#   return _multicube_statistics
    print('\n\nmake_mean_of_cube_list', cube_list, metric)
    operation = iris.analysis.MEAN # or MAX or MIN or STD etc
    #eaned_cubes_list = [cube.collapsed(DIM, operation) for cube in orig_cubes]

    #_fix_cube_attributes(cube_list)

    for c,cu in enumerate(cube_list):
        cu.data = np.ma.array(cu.data.astype('float64')) # cast to masked array f64
        try: iris.coord_categorisation.add_year(cu, 'time', name='year')
        except: pass

        print('cube_list[',c,']:', cube_list[c])
        #print(cu.coords())
        print(cu) #beMetadata())
        continue
        for co,coo in enumerate(cu.coords()):
            if c > 0:
                print(coo.standard_name)
                print(coo.points - cube_list[0].coords(coo.standard_name))
        #print('data type: ', type(cube_list[c].data))
        #print('data:', cu.data)

    meaned_cube = iris.cube.CubeList(cube_list).concatenate_cube()

    print(meaned_cube)
    mean_cube = meaned_cube.collapsed('time', operation)#.data

    #ean_cube  = iris.cube.CubeList([cube for cube in cube_list]).merge()
    print(mean_cube)
    return mean_cube



def make_mean_of_cube_list(cube_list, metric='mean'):
    """
    Takes the mean of a list of cubes (not an iris.cube.CubeList).
    Assumes all the cubes are the same shape.
    Assumes annual 1D data.
    """
    # Fix empty times
    full_times = {}
    times = []
    if len(cube_list)==1:
        return cube_list[0]

    year_counts = {}

    # perform checks:
    for i, cube in enumerate(cube_list):
        try: iris.coord_categorisation.add_year(cube, 'time', name='year')
        except: pass
        years = cube.coord('year')

        count = len(years.points)
        if year_counts.get(count, False):
            year_counts[count].append(i)
        else:
            year_counts[count]= [i, ]

        # print('make_mean_of_cube_list: performing checks:',)
        for year in cube.coord('year').points:
            try:
                full_times[year] += 1
            except:
                full_times[year] = 1

    if len(year_counts.keys())>1:
        print('make_mean_of_cube_list: ERROR: more than one dataset length:', year_counts)
        print('make_mean_of_cube_list: ERROR: list of years:', full_times)
        #assert 0

    print('make_mean_of_cube_list: INFO:, everything is right size:', year_counts)
    output_datas = {}
    cube_mean=cube_list[0].copy()
    for i, cube in enumerate(cube_list):
        print('make_mean_of_cube_list:', metric, i, 'of', len(cube_list))
        years = cube.coord('year')
        for yr, year in enumerate(years.points):
            if output_datas.get(year, False):
                output_datas[year].append(cube.data[yr])
            else:
                output_datas[year] = [cube.data[yr], ]
    datas = []
    for yr in sorted(output_datas.keys()):
        if metric=='mean':
            print('calculating', metric, yr, np.sum(output_datas[yr])/float(len(output_datas[yr])),':',output_datas[yr])
            #datas.append(np.sum(output_datas[yr])/float(len(output_datas[yr])))
            datas.append(np.mean(output_datas[yr]))
        if metric=='min':
            print('calculating', metric, yr, np.min(output_datas[yr]))
            datas.append(np.min(output_datas[yr]))
        if metric=='max':
            print('calculating', metric, yr, np.max(output_datas[yr]))
            datas.append(np.max(output_datas[yr]))

    cube_mean.data = np.array(datas)
    # need to align the time range...
    #or i, cube in enumerate(cube_list[1:]):
    #   cube_mean.data+=cube.data
    #ube_mean.data = cube_mean.data/ float(len(cube_list))
    if metric=='mean' and np.array_equal(cube_mean.data, cube_list[0].data):
        assert 0
    if cube_mean.data.max()> np.max([c.data.max() for c in cube_list]):
        print('something is not working here:', metric, cube_mean.data.max(), '>', [c.data.max() for c in cube_list])
        assert 0

    return cube_mean



def timeplot(cube, **kwargs):
    """
    Create a time series plot from the cube.

    Note that this function simple does the plotting, it does not save the
    image or do any of the complex work. This function also takes and of the
    key word arguments accepted by the matplotlib.pyplot.plot function.
    These arguments are typically, color, linewidth, linestyle, etc...

    If there's only one datapoint in the cube, it is plotted as a
    horizontal line.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube

    """
    cubedata = np.ma.array(cube.data)
    if len(cubedata.compressed()) == 1:
        plt.axhline(cubedata.compressed(), **kwargs)
        return

    times = diagtools.cube_time_to_float(cube)
    plt.plot(times, cubedata, **kwargs)

ssp_title_dict = {
    'ssp119': 'SSP1 1.9',
    'ssp125': 'SSP1 2.5',
    'ssp126': 'SSP1 2.6',
    'ssp245': 'SSP2 4.5',
    'ssp370': 'SSP3 7.0',
    'ssp585': 'SSP5 8.5',
    }

def quick_ts_plot(cfg, data_dict, path, short_namei=None, dataseti=None, expi = None):
#    if os.path.exists(path):
#        return
    title = ''
    for t in dataseti, short_namei, expi:
        if t: title = ' '.join([title, t])
    print('quick_ts_plot: INFO: plotting:', title)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    something = 0
    for (dataset, short_name, exp, ensemble), cube in data_dict.items():
        if dataseti and dataset!= dataseti: continue
        if short_namei and short_namei != short_name: continue
        #print((dataset, short_name, exp, ensemble), (dataseti, short_namei, expi))

        if expi and expi != exp: continue

        #print('adding', (dataset, short_name, exp, ensemble))
        times = cube_to_years(cube)

        if isinstance(cube, dict):
            data = cube[short_name]
        else:
            data = cube.data

        label = ''
        for t in (dataset, short_name, exp, ensemble):
            if t in [short_namei, dataseti,expi]: continue
            label = ' '.join([label, t])
        plt.plot(times, data, label=label)
        something+=1
    if not something:
        plt.close()
        return
    plt.title(title)
    plt.legend(fontsize="x-small")

    plt.savefig(path)
    plt.close()


def plot_data_dict(cfg, data_dict):
    """
    Make a simple plot of every model-shortname pair.
    """
    datasets = {}
    short_names = {}
    exps = {}
    image_extention = diagtools.get_image_format(cfg)

    for (dataset, short_name, exp, ensemble), cube in data_dict.items():
        if short_name in ['areacello', 'areacella']: continue
        datasets[dataset] = True
        short_names[short_name] = True
        exps[exp] = True

    for dataseti, short_namei in product(datasets.keys(), short_names.keys()):
        path = diagtools.folder([cfg['plot_dir'], 'data_dict_checks' ])
        path += '_'.join(['data_dict_checks', short_namei, dataseti]) + image_extention
        quick_ts_plot(cfg, data_dict, path, short_namei=short_namei, dataseti=dataseti)
        for expi in exps:
            continue
            path = diagtools.folder([cfg['plot_dir'], 'data_dict_checks' ])
            path += '_'.join(['data_dict_checks', short_namei, dataseti, expi]) + image_extention
            quick_ts_plot(cfg, data_dict, path, short_namei=short_namei, dataseti=dataseti,expi=expi)



def moving_average(cube, window):
    """
    Calculate a moving average.

    The window is a string which is a number and a measuremet of time.
    For instance, the following are acceptable window strings:

    * ``5 days``
    * ``12 years``
    * ``1 month``
    * ``5 yr``

    Also note the the value used is the total width of the window.
    For instance, if the window provided was '10 years', the the moving
    average returned would be the average of all values within 5 years
    of the central value.

    In the case of edge conditions, at the start an end of the data, they
    only include the average of the data available. Ie the first value
    in the moving average of a ``10 year`` window will only include the average
    of the five subsequent years.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube
    window: str
        A description of the window to use for the

    Returns
    ----------
    iris.cube.Cube:
        A cube with the movinage average set as the data points.

    """
    if 'time' not in [coord.name() for coord in cube.coords()]:
        return

    window = window.split()
    window_len = int(window[0]) / 2.
    win_units = str(window[1])

    if win_units not in [
            'years', 'yrs',
            'year', 'yr'
    ]:
        raise ValueError("Moving average window units not recognised: " +
                         "{}".format(win_units))

    times = cube.coord('time').units.num2date(cube.coord('time').points)

    cal_dt = diagtools.guess_calendar_datetime(cube)

    output = []

    times = np.array([
        cal_dt(time_itr.year, time_itr.month, time_itr.day, time_itr.hour,
               time_itr.minute) for time_itr in times
    ])

    for time_itr in times:
        tmin = cal_dt(time_itr.year - window_len, time_itr.month,
                      time_itr.day, time_itr.hour, time_itr.minute)

        tmax = cal_dt(time_itr.year + window_len, time_itr.month,
                      time_itr.day, time_itr.hour, time_itr.minute)

        #if time_itr.year - times.min().year < window_len:

        arr = np.ma.masked_where((times < tmin) + (times > tmax), cube.data)
        output.append(arr.mean())
    cube.data = np.array(output)
    return cube


def calculate_anomaly(cube, anomaly, calc_average=False):
    """
    Calculate the anomaly using a specified time range.

    The anomaly window is a list which includes a starting year and and end
    year to indicate the start and end of the time period in which to calculate
    the anomaly.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube
    anomaly: list
        A start year and end year to calculate an anomaly.
    calc_average: bool
        Flag to return the average of the anomaly period, instead of subtracting
        it from the cube.
    Returns
    ----------
    iris.cube.Cube:
        A cube with the anomaly calculated.
    """
    start_year = int(np.array(anomaly).min())
    end_year = int(np.array(anomaly).max())
    end_day = 31
    time_units = cube.coord('time').units
    if time_units.calendar == '360_day':
        end_day = 30

    start_date = datetime.datetime(int(start_year), 1, 1)
    end_date = datetime.datetime(int(end_year), 12, end_day)

    t_1 = time_units.date2num(start_date)
    t_2 = time_units.date2num(end_date)
    constraint = iris.Constraint(
        time=lambda t: t_1 < time_units.date2num(t.point) < t_2)

    new_cube = cube.extract(constraint)
    if new_cube is None:
        return None
    mean = new_cube.data.mean()
    if calc_average:
       return mean
    cube.data = cube.data - mean
    return cube


def get_threshold_exceedance_date(cube, threshold, last_year=2090.):
    """
    Calculate the threshold exceedance date.

    Assumes that model run ends 2100, and uses 20 year window.
    """
    if isinstance(threshold, float) or isinstance(threshold, int):
        loc = np.where(cube.data > threshold)[0]
        if not len(loc): return None
        times = cube.coord('time').units.num2date(
            cube.coord('time').points)
        time = times[loc[0]]
        if time.year > last_year:
            return None
        return time

    threshold_os_dict = { '1.5os': 1.5,'1.6os': 1.6,'1.7os': 1.7,'1.8os': 1.8,'1.9os': 1.9,'2.os': 2.,'3.os': 3.}
    if threshold not in threshold_os_dict.keys():
        assert 0
    loc = np.where(cube.data > threshold_os_dict[threshold])[0]
    if not len(loc): return None # never reaches threshold
    times = cube.coord('time').units.num2date(
        cube.coord('time').points)
    index = loc[-1]
    if index >  len(times)-10:
       # last time step
       return None
    time = times[index+1] # one time step after last time tas is above threshold.
    if time.year > last_year:
        return None
    return time


def print_exceedance_dates(cfg, exceedance_dates, window = 10, short_name = 'tas',mip='Amon', preprocessor='prep_1', grid = 'gn'):
    """
    prints the exceedance_dates in a format ready to go into a recipe.
    exceednace key: (metadata['exp'], metadata['ensemble'], threshold)
    Like this:
      tas_ssp119_15:
        short_name: tas
        preprocessor: prep_1
        additional_datasets:
          - {dataset: UKESM1-0-LL, project: CMIP6, mip: Amon, exp: [historical, ssp119], ensemble: r1i1p1f2, start_year: 2014, end_year: 2034, grid: gn}
          - {dataset: UKESM1-0-LL, project: CMIP6, mip: Amon, exp: [historical, ssp119], ensemble: r3i1p1f2, start_year: 2013, end_year: 2033, grid: gn}
          - {dataset: UKESM1-0-LL, project: CMIP6, mip: Amon, exp: [historical, ssp119], ensemble: r4i1p1f2, start_year: 2017, end_year: 2037, grid: gn}
    """
    # Define the contents
    exps = set()
    ensembles = set()
    thresholds = set()

    for (exp, ens, thresh) in exceedance_dates.keys():
        if exp == 'historical':
            continue
        # print(exp, exp.find('-')+1, exp[exp.find('-')+1:])
        # exp = exp[exp.find('-'):]
        exps.add(exp)
        ensembles.add(ens)
        thresholds.add(thresh)

    txt='      # Procedurally generated recipe contents:'
    # Add the historical ensembles:
    lines = []
    lines.append('\n') #
    lines.append('      '+ '_'.join([short_name, 'historical'])+':') #  tas_ssp119_15:
    lines.append('        short_name: '+ short_name)
    lines.append('        preprocessor: '+ preprocessor)
    lines.append('        additional_datasets:')
    for ens in sorted(ensembles):
        lines.append('         - {'
                     'dataset: UKESM1-0-LL, '
                     'project: CMIP6, '
                     'mip: ' + mip + ', '
                     'exp: historical, '
                     'ensemble: ' + ens + ', '
                     'start_year: 1850, '
                     'end_year: 1900, '
                     'grid: ' + grid + '}'
                     )
    txt += '\n'.join(lines)

    # For each combination of short_name, threshold:
    for exp, thresh in product(sorted(exps), sorted(thresholds)):
        ssp = exp[exp.find('-')+1:]
        lines = []
        lines.append('\n') #
        lines.append('      '+ '_'.join([short_name, ssp, str(thresh)])+':') #  tas_ssp119_15:
        lines.append('        short_name: '+ short_name)
        lines.append('        preprocessor: '+ preprocessor)
        lines.append('        additional_datasets:')

        # matches = []
        for ens in sorted(ensembles):
            print(exp, thresh, ens)
            try:
                exceedance_date = float(exceedance_dates[(exp, ens, thresh)])
            except:
                continue

            start_year = str(int(exceedance_date - window))
            end_year = str(int(exceedance_date + window))

            # What if end year is after 2100?
            if int(end_year)> 2099:
                continue

            lines.append('         - {'
                         'dataset: UKESM1-0-LL, '
                         'project: CMIP6, '
                         'mip: ' + mip + ', '
                         'exp: [historical, ' + ssp + '], '
                         'ensemble: ' + ens + ', '
                         'start_year: ' + start_year + ', '
                         'end_year: ' + end_year + ', '
                         'grid: ' + grid + '}'
                         )
        if len(lines) == 5:
            continue
        txt += '\n'.join(lines)

    txt += '\n'
    print(txt)
    fn = cfg['work_dir']+'/new_recipe.yml'
    print('Saved to: ', fn)
    out = open(fn, 'w')
    out.write(txt)
    out.close()


def marine_gt(data_dict, short, gt): #, cumul=False):
    """
    Calculate global from the data dictionary.
    """
    areas = {}
    test_data_dict(data_dict)
    for (dataset, short_name, exp, ensemble), cube in data_dict.items():
        if short_name == 'areacello':
            area = cube.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
            areas[dataset] =  area

    tmp_dict = {}
    for (dataset, short_name, exp, ensemble), cube in data_dict.items():
        if short_name != short:
            continue
        if (dataset, gt, exp, ensemble) in data_dict.items():
            continue
        if dataset in ['CMIP6',]: continue

        area = areas.get(dataset, None)
        if not area:
            print(dataset, 'area not found', sorted(areas.keys()))
            assert 0
        #if not isinstance(area, (type(np.array([0])), float, list, )):
        # assume properly masked! (done in preprocessor)
        # print(dataset, 'areacella cube:', area)
        #area = area.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
        #areas[dataset] = area.data[0]

        cubegt = cube.copy()

        time_units = cubegt.coord('time').units
        if time_units.calendar == '360_day':
            sec_peryear = 360*24*60*60.
        else:
            sec_peryear = 365*24*60*60.

        if cube.units == cf_units.Unit('kg m-2 s-1'):
            cubegt.data = cube.data * area.data * 1.E-12 * sec_peryear
            cubegt.units = cf_units.Unit('Pg yr^-1') #-1}$')
        elif cube.units == cf_units.Unit('kg m-2'):
            cubegt.data = cube.data * area.data * 1.E-12
            cubegt.units = cf_units.Unit('Pg')
        elif cube.units == cf_units.Unit('mol m-2 s-1'):
            cubegt.data = cube.data * area.data * 12.0107* 1.E-15 * sec_peryear
            cubegt.units = cf_units.Unit('Pg yr^-1')
        else:
            print('Units not Recognised:', cube.units)
            assert 0
    #    if cumul:
   #         #print(cubegt.data, np.ma.masked_invalid(cubegt.data), np.cumsum(np.ma.masked_invalid(cubegt.data)))
  #          #assert 0
 #           cubegt.data = np.cumsum(np.ma.masked_invalid(cubegt.data))
#            cubegt.units = cf_units.Unit('Pg yr^-1')
        tmp_dict[(dataset, gt, exp, ensemble)] = cubegt
    data_dict.update(tmp_dict)
    return data_dict


def calculate_cumulative(data_dict, short_name, cumul_name, new_units=''):
    """
    Calculate the cumulative sum of the annual data.
    """
    hist_datas = {}
    tmp_dict={}
    for (dataset, short, exp, ensemble), cube in data_dict.items():
        if short_name != short:
            continue
        if exp not in ['historical', ]:
            continue
        print('calculate_cumulative', (dataset, short, exp, ensemble)) #, cube.units)
        hist_cumul_cube = cube.copy()
        times = diagtools.cube_time_to_float(hist_cumul_cube)
        #print(hist_cumul_cube.data)
        #print(np.cumsum(np.ma.masked_invalid(hist_cumul_cube.data)))

        hist_cumul_cube.data = np.cumsum(np.ma.masked_invalid(hist_cumul_cube.data))
        tmp_dict[(dataset, cumul_name, exp, ensemble)] = hist_cumul_cube
        print('adding to hist_datas:', (dataset, ensemble))
        hist_datas[(dataset, ensemble)] = {'time': times, 'data': hist_cumul_cube.data}

    # For models (like UKESM), where the hist and ssp have different ensemble ids:
    for mod, ens in data_dict_linked_hist.items():
        print(mod, ens)
        for ens_ssp, ens_hist in ens.items():
             print(mod, ens_ssp, ens_hist)
             hist_datas[(mod, ens_ssp)] = hist_datas[(mod, ens_hist)]

    #calculate the cumulative value, and add the historical point to it.
    # Fails are ensembles thaqt need to be ignored to
    fails = {}
    for (dataset, short, exp, ensemble), cube in data_dict.items():
        if short_name != short:
            continue
        if exp in ['historical', ]:
            continue
        print('calculate_cumulative', (dataset, short, exp, ensemble), cube.data.shape )
        cumul_cube = cube.copy()
        times = diagtools.cube_time_to_float(cumul_cube)
        #f isinstance(ensemble, str):
        hist_dat = hist_datas.get((dataset,ensemble), False)
        if fails.get((dataset, exp, ensemble), False): continue
        if hist_dat is False:
            print('looking for ', ensemble)
            hist_dat = hist_datas.get((dataset,ensemble), False)
            for ens in ensemble:
                if hist_dat: continue
                #if len(hist_dat): continue
                hist_dat = hist_datas.get((dataset, ens), False)
            #if not hist_dat:
                print('unable to find', (dataset, exp, ensemble))
            #    assert 0
        if hist_dat is False:
            print('failed:', cube)
            print('Unable to find hist_dat:', (dataset, short, exp, ensemble))
            for (dax, shortx, expx, ensemblex) in data_dict.keys():
                if dax != dataset: continue
                if shortx != short_name: continue
                #if expx != 'historical': continue
                print('candidates:', (dax, shortx, expx, ensemblex))
            #ssert 0
            fails[(dataset, exp, ensemble)] = True
            continue

        hist_point = get_threshold_point(hist_dat, np.min(times))
        print('found hist point:', hist_point, hist_dat['data'][hist_point], np.min(times))
        if hist_point is None:
            print('Problem with hist point:', hist_point)
            print('hist_dat:', hist_dat)
            print('times', times)
            assert 0
        hist_cumul = hist_dat['data'][hist_point]

        #hist_point = get_threshold_point(hist_datas[ensemble], np.min(times))
        #hist_cumul = hist_datas[ensemble]['data'][hist_point]
        print('iter 4:', cumul_cube.data.shape, hist_cumul.shape)

        cumul_cube.data = cumul_cube.data
        print(hist_cumul, np.cumsum(np.ma.masked_invalid(cumul_cube.data)), times)
        print((dataset, short, exp, ensemble))
        print('iter 5:', cumul_cube.data.shape, hist_cumul.shape)
        cumul_cube.data = np.cumsum(np.ma.masked_invalid(cumul_cube.data)) + hist_cumul

        if len(new_units) >0:
            cumul_cube.units = cf_units.Unit(new_units)
        tmp_dict[(dataset, cumul_name, exp, ensemble)] = cumul_cube
    if len(fails.keys()):
        print ('Add this to the list at the start of this code:')
        for ind, bol in fails.items():
            print('     ', ind, ': True, # No historical run on jasmin')
        assert 0
    data_dict.update(tmp_dict)
    return data_dict


def fgco2gt(data_dict): return marine_gt(data_dict, short='fgco2', gt='fgco2gt')
def intppgt(data_dict): return marine_gt(data_dict, short='intpp', gt='intppgt')
def epc100gt(data_dict): return marine_gt(data_dict, short='epc100', gt='epc100gt')
def intdicgt(data_dict): return marine_gt(data_dict, short='intdic', gt='intdicgt')
def intpocgt(data_dict): return marine_gt(data_dict, short='intpoc', gt='intpocgt')
def fricgt(data_dict): return marine_gt(data_dict, short='fric', gt='fricgt')
def frocgt(data_dict): return marine_gt(data_dict, short='froc', gt='frocgt')
def frcgt(data_dict): return marine_gt(data_dict, short='frc', gt='frcgt')
def fgco2gt_cumul(data_dict):
    data_dict = marine_gt(data_dict, short='fgco2', gt='fgco2gt')
    return calculate_cumulative(data_dict, short_name='fgco2gt', cumul_name='fgco2gt_cumul', new_units = 'Pg')

def test_data_dict( data_dict):
    for index, cube in data_dict.items():
        if len(index) != 4:
            print('test_data_dict', index, 'fail')
            assert 0


def land_gt(data_dict, short='npp', gt='nppgt'):
    """
    Calculate land_gt from the data dictionary.
    """
    areas = {}
    test_data_dict(data_dict)
    for (dataset, short_name, exp, ensemble), cube in data_dict.items():
        if short_name == 'areacella':
            area = cube.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
            areas[dataset] = area

    tmp_dict={}
    for (dataset, short_name, exp, ensemble), cube in data_dict.items():
        area = areas.get(dataset, None)
        if dataset in ['CMIP6',]: continue
        if not area:
            print(dataset, 'area not found', sorted(areas.keys()))
            print('other keddys:',  (dataset, short_name, exp, ensemble))
            assert 0
        # assume properly masked! (done in preprocessor)
        #print(dataset, 'areacella cube:', area)
        #print(area.data)
        #rint('land_gt:', short_name, exp, ensemble, [short,gt])
        if short_name != short:
        #   print('land_gt:', short_name,'!=', short)
            continue
        if (dataset, gt, exp, ensemble) in data_dict.keys() or (gt, exp, ensemble) in tmp_dict:
            print('land_gt:', gt, 'already calculated')
            continue
        #print('land_gt:', dataset, short_name, exp, ensemble, [short,gt])
        cubegt = cube.copy()

        time_units = cubegt.coord('time').units
        if time_units.calendar == '360_day':
            sec_peryear = 360*24*60*60.
        else:
            sec_peryear = 365*24*60*60.

        cubegt.data = cube.data * area.data * 1.E-12 * sec_peryear

        cubegt.units = cf_units.Unit('Pg yr^-1')

        #print('land_gt:', (dataset, gt, exp, ensemble), cubegt.data.mean())
        tmp_dict[(dataset, gt, exp, ensemble)] = cubegt

    data_dict.update(tmp_dict)
    #if short=='nbp':
    #    print(data_dict[(dataset, gt, exp, ensemble)])
    return data_dict


def rhgt(data_dict):
    """
    Calculate rhgt from the data dictionary.
    """
    return land_gt(data_dict, short='rh', gt='rhgt')

def nppgt(data_dict): return land_gt(data_dict, short='npp', gt='nppgt')
def gppgt(data_dict): return land_gt(data_dict, short='gpp', gt='gppgt')
def nbpgt(data_dict):
    return land_gt(data_dict, short='nbp', gt='nbpgt')

def nbpgt_cumul(data_dict):
    data_dict = land_gt(data_dict, short='nbp', gt='nbpgt')
    return calculate_cumulative(data_dict, short_name='nbpgt', cumul_name='nbpgt_cumul', new_units='Pg')


def frc(data_dict):
    """
    Calculate total flux to sea floor from the data dictionary.
    """
    #data_dict = fric(data_dict)
    exps = {}
    ensembles = {}
    datasets = {}
    for (dataset, short_name, exp, ensemble)  in data_dict.keys():
        exps[exp] = True
        ensembles[ensemble] = True
        datasets[dataset] = True

    for dataset, exp, ensemble in product(datasets, exps, ensembles):
        if (dataset, 'fric', exp, ensemble) not in data_dict: continue
        if (dataset, 'froc', exp, ensemble) not in data_dict: continue
        cube = data_dict[(dataset, 'fric', exp, ensemble)].copy()
        cube2 = data_dict[(dataset, 'froc', exp, ensemble)]
        cube.data = cube.data + cube2.data
        data_dict[(dataset, 'frc', exp, ensemble)] = cube
    return data_dict


def exchange(data_dict, inverse=False):
    """
    Calculate exchange from the data dictionary.
    """
    data_dict = rhgt(data_dict)
    data_dict = nppgt(data_dict)

    exps = {}
    ensembles = {}
    datasets = {}
    for (dataset, short_name, exp, ensemble)  in data_dict.keys():
        exps[exp] = True
        ensembles[ensemble] = True
        datasets[dataset] = True

    for exp, ensemble, dataset in product(exps, ensembles, datasets):
        if (dataset,'nppgt', exp, ensemble) not in data_dict: continue
        if inverse == False:
            cube = data_dict[(dataset,'nppgt', exp, ensemble)].copy()
            cube2 = data_dict[(dataset,'rhgt', exp, ensemble)]
            cube.data = cube.data - cube2.data
            data_dict[(dataset,'exchange', exp, ensemble)] = cube

        if inverse == True:
            cube = data_dict[(dataset,'rhgt', exp, ensemble)].copy()
            cube2 = data_dict[(dataset,'nppgt', exp, ensemble)]
            cube.data = cube.data - cube2.data
            data_dict[(dataset,'inverse_exchange', exp, ensemble)] = cube
    return data_dict


def inverse_exchange(data_dict,):
    """
    reverses calculation of exchange.
    """
    return exchange(data_dict, inverse=True)


def tas_norm(data_dict):
    """
    Calculate tas_norm from the data dictionary.
    """
    exps = {}
    ensembles = {}
    baselines={}
    datasets = {}
    for (dataset, short_name, exp, ensemble), cube  in data_dict.items():

        exps[exp] = True
        ensembles[ensemble] = True
        datasets[dataset] = True
        if short_name != 'tas': continue
        if exp != 'historical': continue
        if data_dict.get((dataset, 'tas_norm', exp, ensemble), False): continue

        baselines[(dataset, short_name, ensemble)] = calculate_anomaly(cube, [1850, 1900], calc_average=True)
        print('baseline tas:',(dataset, short_name, ensemble), baselines[(dataset, short_name, ensemble)])

    for exp, ensemble, dataset in product(exps, ensembles, datasets):
        if not (dataset, 'tas', exp, ensemble) in data_dict.keys(): continue
        if data_dict.get((dataset, 'tas_norm', exp, ensemble), False): continue

        cube = data_dict[(dataset, 'tas', exp, ensemble)].copy()
        if isinstance(ensemble, str):
            bline = baselines.get((dataset, 'tas', ensemble), None)
            if not bline:
                print('baseline not found', dataset, 'tas', ensemble)
                for blindex, bl in baselines.items():
                    if blindex[0] != dataset: continue
                    if blindex[1] != 'tas': continue
                    print('looking for:', (dataset, 'tas', ensemble), 'candidate:', blindex, bl)
                continue
                assert 0
            cube.data = cube.data - bline
        else:
            print('looking for ', dataset, ensemble)
            print(baselines.keys())
            dat = baselines.get((dataset, 'tas', ensemble), [])
            for ens in ensemble:
                if dat: continue
                if len(dat): continue
                print('loading', ens, dat, baselines.get((dataset, 'tas', ensemble), []),  baselines.get((dataset, 'tas', ens), []))
                dat =  baselines.get((dataset, 'tas', ens), [])
            if not dat:
                print('unable to find', (dataset, 'tas', ensemble))
                assert 0
            print(dat)
            cube.data = cube.data - dat
        data_dict[(dataset, 'tas_norm', exp, ensemble)] = cube
    return data_dict


def norm_co2(data_dict, short='nppgt'):
    """
    Weight a value according to the ratio of the forcing co2 for each year
    against the average co2 forcing in 1850-1900.
    """
    print(data_dict.keys())
    print(data_dict[('co2', 'historical', 'r1i1p1f2' )]['time'][:50],
          data_dict[('co2', 'historical', 'r1i1p1f2' )]['co2'][:50])
    baseline = np.mean(data_dict[('co2', 'historical', 'r1i1p1f2' )]['co2'][:50])
    new_data_dict = {}
    for (dataset, short_name, exp, ensemble), cube  in data_dict.items():
        if short_name != short: continue
        #if exp != 'historical': continue
        cube = data_dict[(dataset, short, exp, ensemble)].copy()
        print('norm_co2:',dataset, short_name, exp, ensemble, 'baseline:',baseline)
        out = []
        co2_data= data_dict[(dataset, 'co2', exp, ensemble )]['co2']
        if len(cube.data) != len(co2_data):
            times = cube.coord('time').units.num2date(cube.coord('time').points)
            print('times do not match', (dataset, short_name, exp, ensemble), len(cube.data), '!=', len(co2_data))
            for t1, t2 in zip(times, data_dict[(dataset, 'co2', exp, ensemble )]['time']):
                print(short_name, exp, ensemble, short_name+':', t1, 'co2:', t2)
            assert 0
        for d,co2 in zip(cube.data, data_dict[(dataset, 'co2', exp, ensemble)]['co2']):
            out.append(d*baseline/co2)
        cube.data = np.ma.array(out)
        new_data_dict[(dataset, short+'_norm', exp, ensemble)] = cube
    data_dict.update(new_data_dict)
    return data_dict


def norm_co2_nppgt(data_dict): return norm_co2(data_dict, short='nppgt')
def norm_co2_rhgt(data_dict): return norm_co2(data_dict, short='rhgt')
def norm_co2_exchange(data_dict): return norm_co2(data_dict, short='exchange')
def norm_co2_fgco2gt(data_dict): return norm_co2(data_dict, short='fgco2gt')


def print_data_dict(data_dict):
    for i, index in enumerate( sorted(data_dict.keys())):
        (dataset, short_name, exp, ensemble) = index

        try: print('d_d:',i, index, data_dict[index].data.mean() )
        except: print('d_d:',i, index)
        if 'CMIP6' in index and 'ensemble_mean' in index:
            print('d_d', index, 'data:', data_dict[index])


def load_scenario_carbon(cfg, data_dict):
    data_dict = load_co2_forcing(cfg, data_dict) # co2
    data_dict = calc_atmos_carbon(cfg, data_dict)
    return data_dict


def standardized_ens(ens):
    if isinstance(ens, str): return ens
    return '_'.join(ens)

def standardized_exps(exp):
    experiments = ['historical', 'piControl', 'ssp119', 'ssp126', 'ssp245', 'ssp370',
                   'ssp434', 'ssp585', 'ssp534-over']
    if exp in experiments:
        return exp
#    if isinstance(exp, str):
#        print('exp not recognised:', exp)
#        assert 0

    if isinstance(exp, list):
        exp = tuple(exp)
    if isinstance(exp, tuple):
        exp = '_'.join(exp)

    if exp in experiments:
        return exp

    exp = exp.replace('-', '_')
    if exp in ['historical_ssp534_over', 'historical_ssp585_ssp534_over', 'ssp534_over']:
        return 'ssp534-over'

    for exp1 in experiments:
        if '_'.join(['historical', exp1]) == exp:
            return exp1
    print('exp not recognised:', exp)
    assert 0


def load_timeseries(cfg, short_names):
    """
    Load times series as a dict.

    Dict is :
    data_dict[(dataset, short_name, exp, ensemble) ] = cube
    assume only one model
    """
    data_dict_shelve = diagtools.folder([cfg['work_dir'], 'gwt_timeseries'])+'data_dict.shelve'
    overwrite_shelve=False

    data_dict = {}

    if not overwrite_shelve and glob.glob(data_dict_shelve+'*'):
        print('loading:', data_dict_shelve )
        sh = shelve.open(data_dict_shelve)
        data_dict = sh['data_dict']
        sh.close()
        return data_dict

    transforms = {
        'fgco2gt': ['fgco2', 'areacello'],
        'gppgt': ['gpp', 'areacella'],
        'nppgt': ['npp', 'areacella'],
        'nbpgt': ['nbp', 'areacella'],
        'rhgt': ['rh', 'areacella'],
        'epc100gt': ['epc100', 'areacello'],
        'intppgt': ['intpp', 'areacello'],
        'intdicgt': ['intdic', 'areacello'],
        'intpocgt': ['intpoc', 'areacello'],
        'fricgt': ['fric', 'areacello'],
        'frocgt': ['froc', 'areacello'],
        'frc': ['fric', 'froc', 'areacello'],
        'frcgt': ['frc', ],
        'exchange': ['rh', 'npp', 'areacella'],
        'inverse_exchange': ['rh', 'npp', 'areacella'],
        'tas_norm': ['tas', ],
        'nppgt_norm': ['nppgt', ],
        'rhgt_norm': ['rhgt', ],
        'exchange_norm': ['exchange', ],
        'fgco2gt_norm': ['fgco2gt', ],
        'fgco2gt_cumul': ['fgco2', 'areacello' ],
        'nbpgt_cumul' : ['nbp', 'areacella' ],
        # 'tls': ['nbp', 'nbpgt', 'luegt']
        }

    transforms_functions = {
        'fgco2gt': fgco2gt,
        'fgco2gt_cumul': fgco2gt_cumul,
        'gppgt': gppgt,
        'nppgt': nppgt,
        'nbpgt': nbpgt,
        'nbpgt_cumul': nbpgt_cumul,
        'rhgt': rhgt,
        'epc100gt': epc100gt,
        'intppgt': intppgt,
        'intdicgt': intdicgt,
        'intpocgt': intpocgt,
        'fricgt': fricgt,
        'frocgt': frocgt,
        'frc': frc,
        'frcgt': frcgt,
        'exchange': exchange,
        'inverse_exchange': inverse_exchange,
        'tas_norm': tas_norm,
        'nppgt_norm':norm_co2_nppgt,
        'rhgt_norm':norm_co2_rhgt,
        'exchange_norm':norm_co2_exchange,
        'fgco2gt_norm':norm_co2_fgco2gt,
        }

    short_names_to_load = short_names.copy()

    for sn in short_names:
        if sn in transforms:
            short_names_to_load.extend(transforms[sn])

    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info('load_timeseries:\t%s', metadata_filename)

        metadatas = diagtools.get_input_files(cfg, index=index)
        for fn in sorted(metadatas):
            short_name = metadatas[fn]['short_name']
            exp = standardized_exps(metadatas[fn]['exp'])
            dataset = metadatas[fn]['dataset']
            ensemble = standardized_ens(metadatas[fn]['ensemble'])
            if isinstance(exp, list): exp = tuple(exp)
            if isinstance(ensemble, list): ensemble = tuple(ensemble)

            if data_dict.get((dataset, short_name, exp, ensemble), False): continue
            if isinstance(ensemble ,list): ensemble = tuple(ensemble)
            print(dataset, short_name, exp, ensemble)

            if mod_exp_ens_skips.get((dataset, exp, ensemble), False):
                continue
            if exp == 'ssp534-over': continue

            # if dataset in data_dict_skips.keys():
            #     if ensemble in data_dict_skips[dataset]:
            #         continue

            if short_name not in short_names_to_load:
                continue
            if data_dict.get((dataset, short_name, exp, ensemble), False):
                continue
            cube = iris.load_cube(fn)
            #cube = diagtools.bgc_units(cube, short_name)

            print('load_timeseries:\t%s successfull loaded data:', (dataset, short_name, exp, ensemble), 'mean:', cube.data.mean())

            data_dict[(dataset, short_name, exp, ensemble)] = cube

#    if 'co2' in short_names_to_load:
#        data_dict = load_co2_forcing(cfg, data_dict)

    if 'luegt' in short_names_to_load or 'tls' in short_names_to_load:
        data_dict = load_luegt(cfg, data_dict)
        #if 'tls' in short_names_to_load:
        #    data_dict = calc_tls(cfg, data_dict)
        #    print(data_dict.keys())
        #    assert 0
    #print(short_names_to_load)

#    if set(['emissions', 'cumul_emissions']) & set(short_names_to_load):
#        data_dict = load_emissions_forcing(cfg, data_dict)

    for sn in short_names_to_load:
        if sn in transforms:
            data_dict = transforms_functions[sn](data_dict)

    if 'tls' in short_names_to_load:
        data_dict = calc_tls(cfg, data_dict)
        #rint(data_dict.keys())

#    if 'atmos_carbon' in short_names_to_load:
#        data_dict = calc_atmos_carbon(cfg, data_dict)

    #print_data_dict(data_dict)
    data_dict = calc_model_mean(cfg, short_names_to_load, data_dict)
    save_data_dict(data_dict, data_dict_shelve)
    return data_dict


def calc_model_mean(cfg, short_names_in, data_dict):

    debug = False
    calculate_model_mean = True
    if calculate_model_mean:
        short_names, exps, datasets = {}, {}, {}
        for (dataset, short_name, exp, ensemble) in  data_dict.keys():
            if mod_exp_ens_skips.get((dataset, exp, ensemble), False):
                continue
            if exp == 'ssp534-over': continue

            # if dataset in data_dict_skips.keys():
            #     if ensemble in data_dict_skips[dataset]:
            #         continue
            if short_name in ['areacello', 'areacella']:continue
            if short_name not in short_names_in:continue
            if dataset in ['CMIP6', ]: continue
            short_names[short_name] = True
            exps[exp] = True
            datasets[dataset] = True

        print(calc_model_mean,short_names, exps, datasets)
        for short_name, exp, dataset in product(short_names.keys(), exps.keys(), datasets.keys()):
            if short_name in ['areacello', 'areacella']:
                continue
            if short_name not in short_names_in:continue

            if dataset in ['CMIP6', ]: continue
            if debug:
                if dataset == 'ACCESS-ESM1-5' and exp == 'ssp126' and short_name =='tas': pass
                else: continue
            else:
                if  data_dict.get((dataset, short_name, exp, 'ensemble_mean'), False): # exists already
                    continue

            print('calculate_model_mean',  short_name, exp, dataset)
            cubes = []
            for (dataset_i,short_name_i, exp_i, ensemble_i),cube in  data_dict.items():
                if dataset != dataset_i: continue
                if short_name != short_name_i: continue
                if exp_i != exp: continue
                if ensemble_i == 'ensemble_mean': continue
                if short_name in ['co2', 'luegt', 'tls', 'atmos_carbon']: # 'emissions', 'cumul_emissions',
                     pass
                     #continue
                #rint("calculate_model_mean: including:", dataset_i,short_name_i, exp_i, ensemble_i)
                else:
                    cube = regrid_time(cube, 'yr')
                cubes.append(cube)

            if debug: print('debug', short_name, exp, dataset, len(cubes)) #ubes)
            if not len(cubes):
                continue

            if len(cubes) == 1:
                data_dict[(dataset, short_name, exp, 'ensemble_mean')] = cubes[0]
                data_dict[(dataset, short_name, exp, 'ensemble_min')] = cubes[0]
                data_dict[(dataset, short_name, exp, 'ensemble_max')] = cubes[0]


                continue
            if isinstance(cubes[0], dict):
                data_dict[(dataset, short_name, exp, 'ensemble_mean')] = make_mean_of_dict_list(cubes, short_name)
                data_dict[(dataset, short_name, exp, 'ensemble_min')] = make_mean_of_dict_list(cubes, short_name, metric='min')
                data_dict[(dataset, short_name, exp, 'ensemble_max')] = make_mean_of_dict_list(cubes, short_name, metric='max')
            else:
                data_dict[(dataset, short_name, exp, 'ensemble_mean')] = make_mean_of_cube_list(cubes)
                data_dict[(dataset, short_name, exp, 'ensemble_min')] = make_mean_of_cube_list(cubes, metric='min')
                data_dict[(dataset, short_name, exp, 'ensemble_max')] = make_mean_of_cube_list(cubes, metric='max')

            if debug:
                print('debug:, mean cube data:', data_dict[(dataset, short_name, exp, 'ensemble_mean')].data)
                for c,cu  in enumerate(cubes):
                    print('debug:',c, cu.data)

            # test:
            if len(cubes) > 1:
                for c,cu in enumerate(cubes):
                    if isinstance(cu, dict):continue
                    if np.array_equal(data_dict[(dataset, short_name, exp, 'ensemble_mean')].data, cu.data):
                        print('ERROR: ensemble mean cube is identical to one of the input cubes.')
                        print('mean cube:', data_dict[(dataset, short_name, exp, 'ensemble_mean')].data[:3], '...')
                        print('mean cube:', cu.data[:3], '...')
                        print('source:', short_name, exp, dataset, len(cubes))
                        assert 0

    #print_data_dict(data_dict)
    #save_data_dict(data_dict, data_dict_shelve)
    #assert 0

    calculate_cmip6_mean = True
    if calculate_cmip6_mean:
        # Calculate the mean of multiple datasets ensemble_means
        short_names, exps, datasets = {}, {}, {}
        for (dataset, short_name, exp, ensemble) in  data_dict.keys():
            if short_name in ['areacello', 'areacella']:continue
            if short_name not in short_names_in:continue

            if mod_exp_ens_skips.get((dataset, exp, ensemble), False):
                continue
            if exp == 'ssp534-over': continue

            # if dataset in data_dict_skips.keys():
            #     if ensemble in data_dict_skips[dataset]:
            #         continue

            short_names[short_name] = True
            exps[exp] = True
            datasets[dataset] = True

        for short_name, exp in product(short_names.keys(), exps.keys()):
            cubes = []
            if short_name in ['areacello', 'areacella']:continue
            print("calculate_cmip6_mean: including:", short_name, exp)
            if data_dict.get(('CMIP6', short_name, exp, 'ensemble_mean'), False): continue

            for (dataset_i,short_name_i, exp_i, ensemble_i),cube in  data_dict.items():
                if dataset_i == 'CMIP6': continue
                if short_name != short_name_i: continue
                if exp_i != exp: continue
                if ensemble_i != 'ensemble_mean': continue
                if short_name in ['co2', 'luegt', 'tls', 'atmos_carbon']: #  'emissions', 'cumul_emissions',
                     pass #ue
                else:
                    cube = regrid_time(cube, 'yr')
                cubes.append(cube)
                print("calculate_cmip6_mean: including:", dataset_i,short_name_i, exp_i, ensemble_i)

            print('calculate_cmip6_mean:', short_name, exp, len(cubes))
            if not len(cubes):
                continue
            if len(cubes) == 1:
                data_dict[('CMIP6', short_name, exp, 'ensemble_mean')] = cubes[0]
                continue

            print('calculating:make_mean_of_cube_list', ('CMIP6', short_name, exp, 'ensemble_mean'))
            if isinstance(cubes[0], dict):
                data_dict[('CMIP6', short_name, exp, 'ensemble_mean')] = make_mean_of_dict_list(cubes, short_name)
                data_dict[('CMIP6', short_name, exp, 'ensemble_min')] = make_mean_of_dict_list(cubes, short_name, metric='min')
                data_dict[('CMIP6', short_name, exp, 'ensemble_max')] = make_mean_of_dict_list(cubes, short_name, metric='max')
            else:
                data_dict[('CMIP6', short_name, exp, 'ensemble_mean')] = make_mean_of_cube_list(cubes)
                data_dict[('CMIP6', short_name, exp, 'ensemble_min')] = make_mean_of_cube_list(cubes, metric='min')
                data_dict[('CMIP6', short_name, exp, 'ensemble_max')] = make_mean_of_cube_list(cubes, metric='max')

#    if 'tls' in short_names_in:
#        data_dict = calc_tls(cfg, data_dict)
        #rint(data_dict.keys())

#    if 'atmos_carbon' in short_names_in:
#        data_dict = calc_atmos_carbon(cfg, data_dict)

        #rint(data_dict.keys())

    #print_data_dict(data_dict)
    #assert 0
    return data_dict


def save_data_dict(data_dict, data_dict_shelve):
    print('saving::', data_dict_shelve )
    sh = shelve.open(data_dict_shelve)
    sh['data_dict'] = data_dict
    sh.close()


def load_thresholds(cfg, data_dict, short_names = ['tas', ], thresholds = [1.5, 2., 3., 4., 5.], ):
    """
    Load thresholds  as a dict.

    Dict is :
    data_dict[(dataset, short_name, exp, ensemble) ] = {threshold: year}
    """
    # If the threshold shelve exists already, just use that.
    overwrite_shelve = False
    thresholds_shelve = diagtools.folder([cfg['work_dir'], 'thresholds_dict'])+'thresholds.shelve'
    if not overwrite_shelve and glob.glob(thresholds_shelve+'*'):
        print('opening:', thresholds_shelve)
        sh = shelve.open(thresholds_shelve)
        thresholds_dict = sh['thresholds_dict']
        sh.close()
        return thresholds_dict

    # Calculate baseline tas for 1850-1900.
    thresholds_dict = {}
    baselines = {}
    baseline_cubes = {}

    debug = False # testing for 'ACCESS-ESM1-5', ssp126

    for (dataset, short_name, exp, ensemble), cube in data_dict.items():
        if short_name not in short_names:
             continue

        if exp != 'historical': continue
        if debug:
            if dataset == 'ACCESS-ESM1-5': pass
            else: continue

        baseline_cubes[(dataset, short_name, ensemble)] = cube
        baselines[(dataset, short_name, ensemble)] = calculate_anomaly(cube, [1850, 1900], calc_average=True)
        if debug: print('calculate anomaly:', cube)
        if debug: print('baseline: ', baselines[(dataset, short_name, ensemble)])


    # Do the work here.
    for (dataset, short_name, exp, ensemble), cube in data_dict.items():
        if short_name not in short_names:
             continue
        if debug:
            if dataset == 'ACCESS-ESM1-5' and  exp == 'ssp126': pass
            else: continue
        print('calculating threshold:', (dataset, short_name, exp, ensemble))
        if debug:
            print('debug mode:', (dataset, short_name, exp, ensemble))

        # find the baseline:
        baseline = baselines.get((dataset, 'tas', ensemble), False)

        if baseline is False and dataset in data_dict_linked_ens.keys():
           new_ens = data_dict_linked_ens[dataset].get(ensemble, False)
           baseline = baselines.get((dataset, 'tas', new_ens), False)


#        if baseline is False:
#            baseline = baselines.get((dataset, 'tas', ensemble.replace('f2', 'f3')), False)

        if baseline is False:
             print('No Baseline found', (dataset, short_name, exp, ensemble), 'available baselines are:')
             for bs_index in baselines.keys():
                 if bs_index[0] != dataset: continue
                 print(bs_index, ':', baselines.get(bs_index, 0.))
             assert 0

        for ens in ensemble:
            if baseline: continue
            baseline =  baselines.get((dataset, 'tas', ens), False)

        # calculate the moving average/normalised tmep
        cube2 = moving_average(cube.copy(), '21 years')
        cube2.data = cube2.data - baseline
        if debug:
            print('cube2.data', cube2.data)

        # ERROR: So the exceedence year is always set to the first cube in the series.
        # is that because the mean of multiple cubes is failing, or because the calculation is only done with the first cube?
        # looks like the ensemble mean cube is identical to the first cube.

        # Calculate the exceedance year:
        thresholds_dict[(dataset, short_name, exp, ensemble)] = {}
        for threshold in thresholds:
            time = get_threshold_exceedance_date(cube2, threshold)

            thresholds_dict[(dataset, short_name, exp, ensemble)][threshold] = time
            print("load_thresholds: Found Threshold:",(dataset, short_name, exp, ensemble), threshold, ':', time)
            if time is None: continue
            if float(threshold)>=2. and time.year == 2015:
                thresholds_dict[(dataset, short_name, exp, ensemble)][threshold] = None
                print('Bad time!')
                print('baseline', baseline)
                print('baseline_cube', baseline_cubes[(dataset, 'tas', ensemble)].data - baseline)

                print('data', cube2.data)
                #assert 0

    for i, index in enumerate(thresholds_dict.keys()):
        print('thresholds:', i, index)
        for thresh, year in thresholds_dict[index].items():
            if year:
                print('\t',thresh, year.year)
            else: print('\t',thresh, year)
    #assert 0

    print('Saving:', thresholds_shelve)
    sh = shelve.open(thresholds_shelve)
    sh['thresholds_dict'] = thresholds_dict
    sh.close()
    return thresholds_dict


def cube_to_years(cube):
    """
    Convert from time coordinate into years.

    Takes an iris time coordinate and returns a list of floats.
    Parameters
    ----------
    cube: iris.cube.Cube
        the opened dataset as a cube.

    Returns
    -------
    list
        List of floats showing the time coordinate in decimal time.

    """
    if isinstance(cube, dict):
        return np.array([int(t)+0.5 for t in cube['time']])
    if not cube.coords('year'):
        iris.coord_categorisation.add_year(cube, 'time')
    return cube.coord('year').points


def load_co2_forcing(cfg, data_dict):
    """
    Load annual CO2 data from the auxiliary datasets.

    Unlike the rest of data_dcit, it's isn't loaded as a cube, but rather as a
    dict.
    """
    fold = cfg['auxiliary_data_dir']+'/atmos_co2_forcing/'
    files = glob.glob(fold+'*.dat')
    print(files)
    hist_datas = []
    hist_times = []
    ssp585_datas = []
    ssp585_times = []
    exps={}
    datasets = {}
    ensembles = {'ensemble_mean': True}
    for (dataset, short_name, exp, ensemble)  in data_dict.keys():
        exps[exp] = True
        ensembles[ensemble] = True
        datasets[dataset] = True
    # load the co2 from the file.
    for fn in files:
        #f data_dict.get((dataset, 'co2', key, 'ensemble_mean' ), False): continue
        open_fn = open(fn, 'r')
        key = os.path.basename(fn).replace('_co2.dat', '')
        key = standardized_exps(key)
        times = []
        data = []
        for line in open_fn.readlines()[1:]:
            line = line.split(' ')
            for x in range(len(line)):
                if '' in line: line.remove('')
                if '\n' in line: line.remove('\n')
            t = float(line[0]) + 0.5
            if t > 2100.: continue
            if t < 1850.: assert 0
            times.append(t)
            data.append(float(line[1]))
        if key == 'historical':
            hist_datas = np.array(data).copy()
            hist_times = np.array(times).copy()
        if key == 'ssp585':
            ssp585_datas = np.array(data).copy()
            ssp585_times = np.array(times).copy()
        for ens, dataset  in product(ensembles.keys(), datasets.keys()): # in range(20):
            #ns='r'+str(e)
            #for ens in ['r1', 'r2',  'r3', 'r4', 'r8']:
            data_dict[(dataset, 'co2', key, ens)] = {'time': times, 'co2':data}
            print('load_co2_forcing:\t%s successfull loaded data:', ('co2', key, ens), 'mean:', np.array(data).mean())
        data_dict[(dataset, 'co2', key, 'ensemble_mean' )] = {'time': times, 'co2':data}
        open_fn.close()

    # Check for historical-ssp scenarios pairs.
    tmp_dict = {}
    for (dataset, short_name, exp, ensemble), ssp_cube in data_dict.items():
        continue
        if short_name in ['co2', 'areacella', 'areacello',]:
            continue
        if (dataset, 'co2', 'historical-'+exp, 'historical-ssp585-'+exp, ensemble ) in tmp_dict.keys():
            continue
        if exp == 'historical':
            continue
        ssp_only = exp.replace('historical-', '')
        if ssp_only == 'ssp585-ssp534-over':
            ssp_only = 'ssp534-over'
        #if data_dict.get((dataset, 'co2', ssp_only, ensemble), False): continue

        new_times = []
        new_datas = []
        print((dataset, short_name, exp,(ssp_only), ensemble))
        ssp_times = cube_to_years(ssp_cube)
        min_time = np.array(ssp_times).min()
        if exp == 'historical-ssp585-ssp534-over':
            print(exp, len(new_times), len(new_datas))
            new_times = list(np.ma.masked_where(hist_times<min_time, hist_times).compressed())
            new_datas = list(np.ma.masked_where(hist_times<min_time, hist_datas).compressed())
            print(exp, len(new_times), len(new_datas) )
            new_times.extend(np.ma.masked_where(ssp585_times >= 2040., ssp585_times).compressed())
            new_datas.extend(np.ma.masked_where(ssp585_times >= 2040., ssp585_datas).compressed())
            print(exp, len(new_times), len(new_datas))
            new_times.extend(data_dict[(dataset, 'co2', ssp_only, ensemble)]['time'])
            new_datas.extend(data_dict[(dataset, 'co2', ssp_only, ensemble)]['co2'])
            print(exp, len(new_times), len(new_datas))
        else:
            if min_time > np.array(hist_times).max():
                print(short_name, exp, ensemble, 'no overlap', ('ssp:', min_time, '>', 'hist max:', np.array(hist_times).max()))
                # no overlap
                new_times = data_dict[(dataset, 'co2', ssp_only, ensemble)]['time']
                new_datas = data_dict[(dataset, 'co2', ssp_only, ensemble)]['co2']
            else:
                # Some overlap
                print(short_name, exp, ensemble,'some overlap', (min_time, '<=', np.array(hist_times).max()))
                new_times = list(np.ma.masked_where(hist_times<min_time, hist_times).compressed())
                new_datas = list(np.ma.masked_where(hist_times<min_time, hist_datas).compressed())
                new_times.extend(data_dict[(dataset, 'co2', ssp_only, ensemble)]['time'])
                new_datas.extend(data_dict[(dataset, 'co2', ssp_only, ensemble)]['co2'])

        if len(new_times) != len(ssp_times):
            print('New times do not match old times:',dataset,short_name, exp, ensemble, len(new_times), '!=', len(ssp_times),'\nnew:',new_times, '\nssp:',ssp_times)
            assert 0
        print('co2', exp, ensemble, len(new_times), len(new_datas))
        tmp_dict[(dataset, 'co2', exp, ensemble )] ={'time': new_times, 'co2':new_datas}

    data_dict.update(tmp_dict)
    # make sure the ensemble mean is set for all co2.
    for (dataset, short_name, exp, ensemble), ssp_cube in data_dict.items():
        if short_name not in ['co2', ]: continue
        tmp_dict[(dataset, short_name, exp, 'ensemble_mean')] = ssp_cube
    data_dict.update(tmp_dict)


    # Save the co2 image:
    path = diagtools.folder(cfg['plot_dir'])
    image_extention = diagtools.get_image_format(cfg)
    path += 'co2_forcing' + image_extention
    if not os.path.exists(path):
#       exp_colours = {'historical':'black',
#                      'ssp119':'green',
#                      'ssp126':'dodgerblue',
#                      'ssp245':'blue',
#                      'ssp370':'purple',
#                      'ssp434':'magenta',
#                      'ssp585': 'red',
#                      'ssp534-over':'orange'}
        for key in exp_colours.keys():
            for (z,a,b,c), dat in data_dict.items():
                if a != 'co2': continue
                if b != key: continue
                plt.plot(dat['time'],
                     dat['co2'],
                     c=exp_colours[key],
                     label=key)
                break

        plt.legend()
        plt.savefig(path)
        plt.close()
    image_extention = diagtools.get_image_format(cfg)
    path += 'co2_forcing_hists' + image_extention
    if not os.path.exists(path):
#        exp_colours = {'historical':'black',
#                       'historical-ssp119':'green',
#                       'historical-ssp126':'dodgerblue',
#                       'historical-ssp245':'blue',
#                       'historical-ssp370':'purple',
#                       'historical-ssp434':'magenta',
#                       'historical-ssp585': 'red',
#                       'historical-ssp585-ssp534-over':'orange'}
        for key in exp_colours.keys():
            for (z,a,b,c), dat in data_dict.items():
                if a != 'co2': continue
                if b != key: continue
                plt.plot(dat['time'],
                     dat['co2'],
                     c=exp_colours[key],
                     label=key)
                break
        plt.legend()
        plt.savefig(path)
        plt.close()
    return data_dict


def load_emissions_forcing_wrong(cfg, data_dict):
    """
    Load annual CO2 data from the auxiliary datasets.

    Unlike the rest of data_dcit, it's isn't loaded as a cube, but rather as a
    dict.
    """
    assert 0
    # fold = cfg['auxiliary_data_dir']+'/emissions/'
    # files = glob.glob(fold+'*.txt')
    # #print(files)
    # #hist_datas = []
    # #hist_times = []
    # #ssp585_datas = []
    # #ssp585_times = []
    # exps = {}
    # ensembles = {'ensemble_mean': True, 'ensemble_min':True, 'ensemble_max':True}
    # datasets = {'CMIP6':True}
    #
    # for (dataset, short_name, exp, ensemble)  in data_dict.keys():
    #     if mod_exp_ens_skips.get((dataset, exp, ensemble), False):
    #         continue
    #     # if dataset in data_dict_skips.keys():
    #     #     if ensemble in data_dict_skips[dataset]:
    #     #         continue
    #     exps[exp] = True
    #     ensembles[ensemble] = True
    #     datasets[dataset] = True
    #
    # # load the co2 from the file.
    # for fn in files:
    #     times = []
    #     data = []
    #
    #     open_fn = open(fn, 'r')
    #     scenario = os.path.basename(fn)
    #     scenario = scenario.replace('UKESM1_', '')
    #     scenario = scenario.replace('.txt', '')
    #     scenario = scenario.replace('historical_', 'historical-')
    #
    #     for line in open_fn.readlines()[2:]:
    #         line = [x.replace('\n', '') for x in line.split(' ')]
    #         t = float(line[0]) + 0.5
    #         #if t > 2100.: continue
    #         times.append(t)
    #         data.append(float(line[3]))
    #         # print (fn, line)
    #     # for t,d in zip(times,data): print(scenario, t,d)
    #     if scenario.find('historical')==-1:
    #         # no need to double up.
    #         continue
    #
    #     # historical-ssp:
    #     times=np.array(times)
    #     cumsumdata = np.cumsum(data)
    #     for dataset, ensemble in product(datasets,ensembles):
    #         data_dict[(dataset, 'emissions', scenario, ensemble)] = {'time': times, 'emissions':data}
    #         #if scenario.find('ssp119')> -1:
    #         #    print('load_emissions_forcing:', scenario, len(times), len(data))
    #         data_dict[(dataset, 'cumul_emissions', scenario, ensemble)] = {'time': times, 'cumul_emissions':cumsumdata}
    #         if dataset == 'CMIP6': print('load_emissions_forcing', dataset, 'cumul_emissions', scenario, ensemble)
    #
    #     # historical only:
    #     histdata = np.ma.masked_where(times > 2015., data).compressed()
    #     histcumsumdata = np.ma.masked_where(times > 2015., cumsumdata).compressed()
    #     histtimes = np.ma.masked_where(times > 2015., times).compressed()
    #     for dataset, ensemble in product(datasets,ensembles):
    #         data_dict[(dataset, 'emissions', 'historical', ensemble)] = {'time': histtimes, 'emissions':histdata}
    #         data_dict[(dataset, 'cumul_emissions', 'historical', ensemble)] = {'time': histtimes, 'cumul_emissions':histcumsumdata}
    #
    #     # SSP only:
    #     scenario = scenario.replace('historical-', '')
    #     times=np.array(times)
    #     data = np.ma.masked_where(times <2015., data).compressed()
    #     cumsumdata = np.ma.masked_where(times <2015., cumsumdata).compressed()
    #     times = np.ma.masked_where(times <2015., times).compressed()
    #     for dataset, ensemble in product(datasets,ensembles):
    #         data_dict[(dataset, 'emissions', scenario, ensemble)] = {'time': times, 'emissions':data}
    #         data_dict[(dataset, 'cumul_emissions', scenario, ensemble)] = {'time': times, 'cumul_emissions':cumsumdata}
    #         if dataset == 'CMIP6': print('load_emissions_forcing', dataset, 'cumul_emissions', scenario, ensemble)
    #
    # return data_dict


def calc_tls(cfg, data_dict):
    """
    Load True  Land Sink by adding Land Use Emissions from file and nbp
    Net biome production.
    """
    exps = {} #'ssp119':True, 'ssp126':True, 'ssp245':True, 'ssp370':True, 'ssp585':True, 'ssp534-over':True, 'historical':True}
    ensembles = {'ensemble_mean':True, 'ensemble_min' :True, 'ensemble_max' :True}

    datasets = {}
    #for (dataset, short_name, exp, ensemble)  in data_dict.keys():
    #    exps[exp] = True
    #    ensembles[ensemble] = True
    #    datasets[dataset] = True
    #print(exps, ensembles,datasets)
    tmp_data_dict = {}
    short = 'tls'
    for (dataset, short, exp, ensemble), cube in data_dict.items():
        print('calc_tls:',(dataset, short, exp, ensemble))
        if short not in ['nbpgt_cumul', ]: # 'nbpgt_cumul']:
             continue
        times = diagtools.cube_time_to_float(cube)
        data = cube.data.copy()
        nbp_dict = {int(t):d for t,d in zip(times, data)}

        luegt_dict = data_dict.get((dataset, 'luegt', exp, ensemble), False)
        if not luegt_dict:
            continue
        for t, d in zip(luegt_dict['time'], luegt_dict['luegt']):
            if not nbp_dict.get(t, nbp_dict.get(int(t), False)):
                print('error:', (dataset, short, exp, ensemble), t,'and', int(t), 'not in', nbp_dict.keys())
                continue
            nbp_dict[int(t)] += d

        new_times = np.array(sorted(nbp_dict.keys()))
        new_data = np.array([nbp_dict[t] for t in new_times])

        tmp_data_dict[(dataset, 'tls', exp, ensemble)] = {'time':new_times+0.5, 'tls': new_data }

    #print(tmp_data_dict.keys(), tmp_data_dict)
    if not len(tmp_data_dict): assert 0
    data_dict.update(tmp_data_dict)
    return data_dict



def calc_atmos_carbon(cfg, data_dict):
    """
    Load remaining atmospheric carbon.
    The old way is wrong:
    cumul_emissions - nbpgt_cumul - fgco2gt_cumul

    This way uses the 1ppm = 2.13 Pg C, from
    https://cdiac.ess-dive.lbl.gov/pns/convert.html
    """
    tmp_data_dict = {}
    new_short = 'atmos_carbon'
    for (dataset, short, exp, ensemble), tmp_data in data_dict.items():
        if short not in ['co2', ]:
            continue
        # print('calc_atmos_carbon:',(dataset, short, exp, ensemble), tmp_data)

        tmp_data = zip_time(tmp_data, short)
        tmp_times, tmp_dat = unzip_time(tmp_data)
        if not len(tmp_dat):
            print('calc_atmos_carbon: ERROR:',(dataset, short, exp, ensemble), 'no data:', tmp_dat)
            assert 0

        tmp_times = np.array(tmp_times)+0.5
        tmp_dat = (np.array(tmp_dat) - 283.15316772 )* 2.13 # Convert ppm into Pg C relative to 1850.

        tmp_data_dict[(dataset, new_short, exp, ensemble)] = {'time': tmp_times, new_short: tmp_dat}
        if exp[:3] == 'ssp':
            tmp_data_dict[(dataset, new_short, 'historical-'+exp, ensemble)] = {'time': tmp_times, new_short: tmp_dat}


    data_dict.update(tmp_data_dict)
    #print(tmp_data_dict)
    return data_dict


def calc_emissions(cfg, data_dict):
    """
    Using the other values, we calculate emissions.
    # emmissions
    """
    return data_dict




def calc_atmos_carbon_old(cfg, data_dict):
    """
    Load remaining atmospheric carbon.
    cumul_emissions - nbpgt_cumul - fgco2gt_cumul
    """
    print('This method is wrong because the cumulative emission data is UKESM only')
    assert 0
    # tmp_data_dict = {}
    # new_short = 'atmos_carbon'
    # for (dataset, short, exp, ensemble), tmp_data in data_dict.items():
    #     print('calc_atmos_carbon:',(dataset, short, exp, ensemble))
    #     if short not in ['cumul_emissions', ]:
    #         continue
    #
    #     # tmp_data = {time:times, 'cumul_emissions': dat}
    #     tmp_data = zip_time(tmp_data, short)
    #     tmp_data = {int(t):d for t, d in tmp_data.items()}
    #     if not len(tmp_data):
    #         print('calc_atmos_carbon: ERROR:',(dataset, short, exp, ensemble), 'no data:', tmp_data)
    #         assert 0
    #
    #     # subtrack nbp & ocean carbon flux.
    #     found = 0
    #     for cube_key in  ['fgco2gt_cumul', 'nbpgt_cumul']:
    #         tmp_cube_data = data_dict.get((dataset, cube_key, exp, ensemble), False)
    #         if not tmp_cube_data:
    #             for ens in ensemble:
    #                 if tmp_cube_data: continue
    #                 tmp_cube_data = data_dict.get((dataset, cube_key, exp, ens), False)
    #             if not tmp_cube_data:
    #                  print('Fail to find:', (dataset, cube_key, exp, ensemble))
    #                  # print(data_dict.keys())
    #                  continue
    #             if not len(tmp_cube_data):
    #                  print('Fail to find:', (dataset, cube_key, exp, ensemble))
    #                  continue
    #
    #         tmp_cube_data = {'time': [int(t) for t in diagtools.cube_time_to_float(tmp_cube_data)],
    #                          cube_key: tmp_cube_data.data}
    #         tmp_cube_data = zip_time(tmp_cube_data, cube_key)
    #         found+=1
    #         for t,d in tmp_cube_data.items():
    #             #t = int(t)
    #             if t not in tmp_data: #t(t,False):
    #                 print('calc_atmos_carbon: ERROR: unable to find year', t, 'in',(dataset, cube_key, exp, ensemble))
    #                 print('output:',tmp_data)
    #                 print('input:', cube_key, tmp_cube_data)
    #                 assert 0
    #
    #             tmp_data[t] = tmp_data[t] - d
    #     if found !=2: continue #assert 0
    #     tmp_times, tmp_dat = unzip_time(tmp_data)
    #     tmp_times = tmp_times+0.5
    #
    #     tmp_data_dict[(dataset, new_short, exp, ensemble)] = {'time': tmp_times, new_short: tmp_dat}
    #
    # data_dict.update(tmp_data_dict)
    # print(tmp_data_dict)
    # return data_dict


def load_luegt(cfg, data_dict):
    """
    Load Land Use Emissions from file.
    """
    # load from file
    tmp_data_dict = {}
    short = 'luegt'
    # open files
    # for each column, create a dict
    # add year: value to the column.
    # edit file header so that it closely matches other data.
    # set one for each dataset

    exps = {'ssp119':True, 'ssp126':True, 'ssp245':True, 'ssp370':True, 'ssp585':True, 'historical':True}
    ensembles = {'ensemble_mean':True, 'ensemble_min' :True, 'ensemble_max' :True}
    datasets = {}
    for (dataset, short_name, exp, ensemble)  in data_dict.keys():
        exps[exp] = True
        ensembles[ensemble] = True
        datasets[dataset] = True
    print(exps, ensembles, datasets)

    #for exp, ensemble,dataset in product(exps.keys(), ensembles.keys(),datasets.keys()):
    #    if data_dict.get((dataset, short, exp, ensemble), False: continue

    # assert 0
    # This data was added to the esmvaltool/diagnostics/ocean/aux_data directory.
    aux_fn = cfg['auxiliary_data_dir']+'/land_usage/landusage_ssp.txt'
    data = {}
    years = []
    header = {}
    print('load_luegt', ensembles, exps)
    with open(aux_fn) as csvDataFile:
        csvReader = csv.reader(csvDataFile, delimiter=';')
        for row_number, row in enumerate(csvReader):
            print('load_luegt: ', row_number, row)
            if row_number == 0:
                for d, da in enumerate(row):
                    print('load_luegt: header: ', d, da)
                    da = da.replace(' ', '')
                    header[d] = da
                    data[da] = []
                print('load_luegt: header:', header)
                continue

            years.append(float(row[0])+0.5)
            for d, da in enumerate(row):

                #if d == 0: continue
                if da == ' ': continue
                data[header[d]].append(float(da))

    for exp, ensemble,dataset in product(exps.keys(), ensembles.keys(),datasets.keys()):
        if data_dict.get((dataset, short, exp, ensemble), False): continue
        print(exp, ensemble)
        da =  data.get(exp, [])
        if not da:
            da =  data.get(exp.replace('historical-',''), [])
        if not da:
            da =  data.get('-'.join(['historical',exp]), [])
        if not da:
            print('problem', exp, 'and', exp.replace('historical-',''), 'not int', data.keys())
            #assert 0
            continue
        if len(da) != len(years):
            print("data and time don't match:", len(da),  len(years))
            y2 = years[:len(da)]
        else:
            y2 = years
        tmp_data_dict[(dataset, short, exp, ensemble)] = {'time': y2, short: da}

    #print(tmp_data_dict.keys())
    #assert 0
    data_dict.update(tmp_data_dict)
    #assert 0
    return data_dict


def get_threshold_point(cube, year):
    """
    get the location of the year provided.
    """
    if isinstance(cube, dict):
        if np.min(np.abs(np.array(cube['time']) - year)) > 5.:
             print('Problem with thresholds:', np.min(cube['time']), '-', np.max(cube['time']), year)
             return None
        arg_min = np.argmin(np.abs(np.array(cube['time']) - year))
        #print(arg_min, year, np.array(cube['time']))
    else:
        try:
            time_units = cube.coord('time').units
        except:
            print('get_threshold_point:', cube, year)
            assert 0
        date = datetime.datetime(int(year), 6, 1)
        t_1 = time_units.date2num(date)
        years = np.array(diagtools.cube_time_to_float(cube))
        if np.min(np.abs(years - year)) > 5.:
             print('Problem with thresholds:', np.min(years), '-', np.max(years), year)
             return None
        arg_min = np.argmin( np.abs(cube.coord('time').points - t_1))
    #print('get_threshold_point',t_1, date, arg_min)
    return arg_min


def get_long_name(name):
    """
    Get a title friendly longname.
    """
    longnames = {
        'tas' : 'Temperature, K',
        'tas_norm' : 'Normalised Temperature, '+r'$\degree$' + 'C',
        'co2' : 'Atmospheric CO'+r'$_{2}$',
        # 'emissions' : 'Anthropogenic emissions',
        # 'cumul_emissions': 'Cumulative Anthropogenic emissions',
        'luegt': 'Land Use Emissions',
        'tls': 'True Land Sink',
        'rh': 'Heterotrophic respiration',
        'intpp' : 'Marine Primary Production',
        'intdic' : 'Dissolved Inorganic Carbon',
        'intpoc' : 'Particulate Organic Carbon',
        'epc100' : 'POC flux at 100m',
        'npp'   : 'Net Primary Production on Land',
        'fgco2' : 'Air sea Flux of CO'+r'$_{2}$',
        'frc':  'Carbon Flux at sea floor',
        'fric':  'Inorganic Carbon Flux at sea floor',
        'froc':  'Organic Carbon Flux at sea floor',
        'nbp': 'Net Biome Production',
        'nbpgt_cumul': 'Cumulative Global Total Net Biome Production, Pg',
        'fgco2gt_cumul': 'Cumulative Global Total Air sea Flux of CO'+r'$_{2}$'+', Pg',
        'atmos_carbon': 'Remnant Atmospheric Anthropogenic CO'+r'$_{2}$'+', Pg',
        # units:
        'Pg yr^-1': 'Pg yr'+r'$^{-1}$',

    }
    long_name = ''
    if name.find('gt_norm') > -1:
        long_name += 'Normalised Global Total '
        name = name[name.find('gt_norm')]
    elif name[-2:] == 'gt':
        long_name += 'Global Total '
        name = name[:-2]

    return long_name + longnames.get(name, name)


def load_ensemble(data_dict, short_name, exp, ensemble):
    dat = data_dict.get(('CMIP6', short_name, exp, ensemble), False)
    for ens in ensemble:
        if dat: continue
        dat = data_dict.get(('CMIP6', short_name, exp, ens), False)
    return dat


# def make_ts_envellope_figure(
#         cfg, data_dict, thresholds_dict,
#         x='time',
#         y='npp',
#         plot_dataset='CMIP6',
#         fill = 'ensemble_means',
#         markers='thresholds',
#         fig=None,
#         ax=None,
#     ):
#     """
#     make a time series envellope figure.
#     x axis and y axis are determined by the short_names provuided in x and y
#     vars.
#     Markers are placed at certain points when the tas goes above thresholds.
#
#     fill = 'ensemble_means': plots only ensemble means
#     fill = 'all_ensembles:' plots all ensembles, nut not ensemble means.
#
#     plot_dataset = 'CMIP6': All datasets
#     plot_dataset = 'model': Only an ensemble of that model.
#
#     Parameters
#     ----------
#     cfg: dict
#         the opened global config dictionairy, passed by ESMValTool.
#     """
#     if x !='time':return
#     draw_line=True
#     exps = {}
#     ensembles = {}
#     datasets = {}
#     for (dataset, short_name, exp, ensemble)  in data_dict.keys():
#          exps[exp] = True
#          ensembles[ensemble] = True
#          datasets[dataset] = True
#
#     if plot_dataset != 'CMIP6':
#         if isinstance(plot_dataset, str):
#             datasets = {plot_dataset:True}
#         else:
#             datasets = {pd:True for pd in plot_dataset}
#             plot_dataset = '_'.join(plot_dataset)
#
#     # set path:
#     image_extention = diagtools.get_image_format(cfg)
#     path = diagtools.folder([cfg['plot_dir'], 'envellope_ts_figure'])
#     path += '_'.join([x, y, markers, fill, plot_dataset]) + image_extention
#     #if os.path.exists(path): return
#
#     exp_colours = {'historical':'black',
#                    'ssp119':'green',
#                    'ssp126':'dodgerblue',
#                    'ssp245':'blue',
#                    'ssp370':'purple',
#                    'ssp434':'magenta',
#                    'ssp585': 'red',
#                    'ssp534-over':'orange',
#                    'historical-ssp119':'green',
#                    'historical-ssp126':'dodgerblue',
#                    'historical-ssp245':'blue',
#                    'historical-ssp370':'purple',
#                    'historical-ssp434':'magenta',
#                    'historical-ssp585': 'red',
#                    'historical-ssp585-ssp534-over':'orange'}
#
#     marker_styles = {1.5: '*', 2.:'o', 3.:'D', 4.:'s', 5.:'X'}
#
#     #if ensemble_mean: ensembles = ['ensemble_mean', ]
#     exps = sorted(exps.keys())
#     exps.reverse()
#     print(exps, ensembles)
#     print(data_dict.keys())
#
#     if fig == None:
#         fig = plt.figure()
#         ax = fig.add_subplot(111)
#         make_figure_here = True
#     else:
#         make_figure_here = False
#         plt.sca(ax)
#
#     x_label,y_label = [], []
#     #number_of_lines=0
#
#     label_dicts = {
#         'tas': ' '.join(['Temperature, K',]), # ''.join([r'$\degree$', 'C'])]),
#         'tas_norm': ' '.join(['Normalised Temperature,', ''.join([r'$\degree$', 'C'])]),
#         'co2': ' '.join(['Atmospheric co2, ppm']),
#         'emissions': ' '.join(['Anthropogenic emissions, Pg/yr']),
#         'cumul_emissions': ' '.join(['Cumulative Anthropogenic emissions, Pg']),
#         'luegt': ' '.join(['Land use emissions, Pg']),
#         'tls': ' '.join(['True Land Sink, Pg']),
#         'atmos_carbon': 'Remnant Anthropogenic CO2, Pg',
#     }
#
#     # Only drawy line and markers for CMIP6 ensemble_mean.
#     #dataset_1 = plot_dataset
#     #ensemble_1 = 'ensemble_mean'
#
#     envellopes = {exp: {} for exp in exps}
#     for exp_1, ensemble_1,dataset_1 in product(exps, ensembles,datasets):
#         x_data, y_data = [], []
#         x_times, y_times = [], []
#         for (dataset, short_name, exp, ensemble), cube in data_dict.items():
#             if short_name not in [x,y]: continue
#             if exp != exp_1: continue
#             if ensemble != ensemble_1: continue
#             if dataset != dataset_1: continue
#
#             if fill == 'ensemble_means' and ensemble_1 != 'ensemble_mean':
#                 #  plots only ensemble means.
#                 continue
#             #if fill == 'all_ensembles' and ensemble_1 == 'ensemble_mean':
#             #    #  plots all ensembles, nut not ensemble means.
#             #    continue
#
#             print('found:', (dataset, short_name, exp, ensemble))
#             if x == 'time' and short_name == y:
#                 x_label = 'Year'
#                 if isinstance(cube, iris.cube.Cube):
#                     x_data = diagtools.cube_time_to_float(cube)
#                     x_times = x_data.copy()
#                     #print('setting x time to ',short_name, exp, ensemble)
#                 else:
#                     x_data = cube['time']
#                     x_times = x_data.copy()
#                     #print('setting x time to ',short_name, exp, ensemble)
#             elif x == short_name and x in label_dicts.keys():
#                 x_data = cube[x].copy()
#                 x_times = cube['time'].copy()
#                 #print('setting x time to ',short_name, exp, ensemble)
#                 x_label = label_dicts[x]
#
#             elif short_name == x:
#                 x_data = np.array(cube.data.copy())
#                 x_times = diagtools.cube_time_to_float(cube)
#                 #print('setting x axis to ',short_name, exp, ensemble, np.min(x_data), np.max(x_data))
#                 x_label = ' '.join([get_long_name(x), str(cube.units)])
#
#             if y == 'time':
#                #print('what kind of crazy person plots time on y axis?')
#                 assert 0
#             elif y == short_name and y in label_dicts.keys():
#                 #['co2', 'emissions', 'cumul_emissions', 'luegt']:
#                 y_data = cube[y].copy()
#                 y_times = cube['time'].copy()
#                 y_label = label_dicts[y]
#                 #print('setting y time to ',short_name, exp, ensemble)
#             elif short_name == y:
#                 #print(short_name, 'is a cube for y ts plot..')
#                 y_data = cube.data.copy()
#                 y_times = diagtools.cube_time_to_float(cube)
#                 #print('setting y time to ',short_name, exp, ensemble, y_data)
#                 y_label = ' '.join([get_long_name(y), str(cube.units)])
#
#         if 0 in [len(x_data), len(y_data), len(x_times), len(y_times)]:
#             #print('no data found',(exp_1, ensemble_1,dataset_1),  x,'vs',y, 'x:', len(x_data), 'y:',len(y_data))
#             continue
#
#         if len(x_data) != len(x_times) or len(y_data) != len(y_times):
#             print('x:', len(x_data), len(x_times), 'y:', len(y_data), len(y_times))
#             assert 0
#
#         # calculate the datasety for the fills
#         for xd, yd in zip(x_data, y_data):
#             if 'x' == 'time': xd = int(xd)
#             if xd < 1851.:
#                 print(plot_dataset, fill, (exp_1, ensemble_1, dataset_1, short_name), xd, yd)
#             if plot_dataset=='CMIP6':
#                 if dataset_1 == 'CMIP6':
#                     # include everything except 'CMIP6 ensemble means.'
#                     continue
#             else:
#                 if dataset_1 != plot_dataset:
#                     continue
#             if fill == 'ensemble_means' and ensemble_1 != 'ensemble_mean':
#                 #  plots only ensemble means.
#                 continue
#             if fill == 'all_ensembles' and ensemble_1 == 'ensemble_mean':
#                 #  plots all ensembles, nut not ensemble means.
#                 continue
#
#             if envellopes[exp_1].get(xd, False):
#                 envellopes[exp_1][xd].append(yd)
#             else:
#                 envellopes[exp_1][xd] = [yd, ]
#             print('made it!',exp_1, ensemble_1, fill, xd,yd)
#
#         # Only draw line for plot_dataset
#         if draw_line and dataset_1 == plot_dataset:
#             x_times = np.ma.array(x_times)
#             y_times = np.ma.array(y_times)
#             #number_of_lines+=1
#             if plot_dataset=='CMIP6':
#                 lw=1.3
#             else:
#                 lw=0.5
#             plt.plot(x_data,
#                  y_data,
#                  lw=lw,
#                  color=exp_colours[exp_1])
#
#         # Only plot markers for plot_dataset ensemble_mean
#         if markers == 'thresholds' and dataset_1 == plot_dataset and ensemble_1 == 'ensemble_mean':
#             try: threshold_times = thresholds_dict[(dataset_1, 'tas', exp_1, ensemble_1)]
#             except:
#                threshold_times = {}
#             ms = 8
#             for threshold, time in threshold_times.items():
#                 if not time:
#                     continue
#                 x_point = get_threshold_point({'time':x_times}, time.year)
#                 y_point = get_threshold_point({'time':y_times}, time.year)
#
#                 plt.plot(x_data[x_point],
#                          y_data[y_point],
#                          marker_styles[threshold],
#                          markersize = ms,
#                          fillstyle='none',
#                          color=exp_colours[exp_1])
#                 plt.plot(x_data[x_point],
#                          y_data[y_point],
#                          'o',
#                          markersize = 2,
#                          #fillstyle='none',
#                          color=exp_colours[exp_1])
#
#     #add_env = True
#     for exp, ranges in envellopes.items():
#          times, mins, maxs = [], [], []
#          # This won't work if the ranges aren't the same!
#          for t in sorted(ranges.keys()):
#              times.append(t)
#              mins.append(np.min(ranges[t]))
#              maxs.append(np.max(ranges[t]))
#          print('enveloppe', exp, times, mins, maxs)
#          if not len(times):
#              continue
#          ax.fill_between(
#              times,
#              mins,
#              maxs,
#              lw=0,
#              alpha=0.5,
#              color=exp_colours[exp])
#
#     #print(x,y,fill, plot_dataset, envellopes)
#
#     if not envellopes:
#         print('No lines plotted')
#         plt.close()
#         return
#
#     exp_colours_leg = {'historical':'black',
#                    'ssp119':'green',
#                    'ssp126':'dodgerblue',
#                    'ssp245':'blue',
#                    'ssp370':'purple',
#                    #'ssp434':'magenta',
#                    'ssp585': 'red',}
#                    #'ssp534-over':'orange'}
#
#     plot_details = {}
#     for exp,color in sorted(exp_colours_leg.items()):
#         plot_details[exp] = {
#                     'c': color,
#                     'ls': '-',
#                     'lw': 2.,
#                     'label': exp
#                 }
#     for thres,ms in sorted(marker_styles.items()):
#         plot_details[str(thres)] = {
#                     'c': 'black',
#                     'marker': ms,
#                     'fillstyle':'none',
#                     'label': '>' + str(thres)+u'\u00B0C'
#                 }
#
#     diagtools.add_legend_outside_right(
#                 plot_details, plt.gca(), column_width=0.175)
#
#     # set labels:
#     plt.xlabel(x_label)
#     plt.ylabel(y_label)
#
#     # set title:
#     if x == 'time':
#         title = get_long_name(y)
#     else:
#         title = ' '.join([get_long_name(x), 'by', get_long_name(y)])
#
#     if plot_dataset != 'all_models':
#         title += ' '+plot_dataset
#
#     plt.title(title)
#
#     print('saving figure:', path)
#     if make_figure_here:
#         plt.savefig(path)
#         plt.close()
#     else:
#         return fig, ax





def make_ts_figure(cfg, data_dict, thresholds_dict, x='time', y='npp',
    markers='thresholds',
    draw_line=True,
    do_moving_average=False,
    plot_dataset='all_models',
    ensemble_mean = False, # does nothing now.
    plot_styles = ['ensemble_mean', ], #['CMIP6_range', 'CMIP6_mean'],
    fig=None,
    ax=None,
    do_legend=True,
    plot_thresholds = [1.5, 2., 3., 4., 5.,],
    skip_historical_ssp=False,
    experiments = ['historical', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'],
    #short_time_range = False,
    ):
    """
    make a 2D figure.
    x axis and y axis are determined by the short_names provuided in x and y
    vars.
    Markers are placed at certain points when the tas goes above thresholds.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    plot_dataset:
        either 'all_models'
        or an individual model name.
    plot_styles: ambition:
        ['ensemble_mean', ] default behaviour before
        CMIP6_mean: Only CMIP6 multi model mean
        CMIP6_range: range bewteen each model mean shown
        CMIP6_full_range: full range between individual model ensemble menmbers
        all_models_means: Each individual models mean is plotted.
        all_models_range: Each individual models range is plotted
    """
    exps = {}
    ensembles = {}
    datasets = {}
    for (dataset, short_name, exp, ensemble)  in data_dict.keys():
         if not exp in experiments: continue
         exps[exp] = True
         ensembles[ensemble] = True
         datasets[dataset] = True
         print(dataset, short_name, exp, ensemble)

    if plot_dataset != 'all_models':
        if isinstance(plot_dataset, str):
            datasets = {plot_dataset:True}
        else:
            datasets = {pd:True for pd in plot_dataset}
            plot_dataset = '_'.join(plot_dataset)

    # set path:
    image_extention = diagtools.get_image_format(cfg)
    path = diagtools.folder([cfg['plot_dir'], 'ts_figure',plot_dataset])
    print(path, plot_styles)
    ensemble_mean_txt = '_'.join(plot_styles)
    path += '_'.join([x, y, markers, ensemble_mean_txt, plot_dataset]) + image_extention

    if do_moving_average:
        path = path.replace(image_extention, '_21ma'+image_extention)
    #if os.path.exists(path): return

#    exp_colours = {'historical':'black',
#                   'ssp119':'green',
#                   'ssp126':'dodgerblue',
#                   'ssp245':'blue',
#                   'ssp370':'purple',
#                   'ssp434':'magenta',
#                   'ssp585': 'red',
#                   'ssp534-over':'orange',
#                   'historical-ssp119':'green',
#                   'historical-ssp126':'dodgerblue',
#                   'historical-ssp245':'blue',
#                   'historical-ssp370':'purple',
#                   'historical-ssp434':'magenta',
#                   'historical-ssp585': 'red',
#                   'historical-ssp585-ssp534-over':'orange'}

    marker_styles = {1.5: '*', 2.:'o', 3.:'D', 4.:'s', 5.:'X'}

    if plot_styles == ['ensemble_mean', ]:
        ensembles = ['ensemble_mean', ]

    exps = sorted(exps.keys())
    exps.reverse()
    print(exps, ensembles)
    print(data_dict.keys())

    if fig == None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        make_figure_here = True
    else:
        make_figure_here = False
        plt.sca(ax)

    x_label,y_label = [], []
    print('\n\n\n\n\nStaring plot:',
        'x:', x,
        'y:', y,
        'markers:', markers,
        'draw_line:', draw_line,
        'do_moving_average:', do_moving_average,
        'plot_styles:', plot_styles)

    number_of_lines=0

    label_dicts = {
#        'tas': ' '.join(['Temperature, K',]), # ''.join([r'$\degree$', 'C'])]),
#        'tas_norm': ' '.join(['Normalised Temperature,', ''.join([r'$\degree$', 'C'])]),
        'co2': ' '.join(['Atmospheric co2, ppm']),
        # 'emissions': ' '.join(['Anthropogenic emissions, Pg/yr']),
        # 'cumul_emissions': ' '.join(['Cumulative Anthropogenic emissions, Pg']),
        'luegt': ' '.join(['Land use emissions, Pg']),
        'tls': ' '.join(['True Land Sink, Pg']),
        'atmos_carbon': 'Remnant Anthropogenic CO2, Pg',
    }

    for exp_1, ensemble_1,dataset_1, plot_style in product(exps, ensembles, datasets, plot_styles):
        x_data, y_data = [], []
        x_times, y_times = [], []

            #plot_styles: ambition:
            #    ['ensemble_mean', ] default behaviour before
            #    CMIP6_mean: Only CMIP6 multi model mean
            #    CMIP6_range: range bewteen each model mean shown
            #    CMIP6_full_range: full range between individual model ensemble menmbers
            #    all_models_means: Each individual models mean is plotted.
            #    all_models_range: Each individual models range is plotted
            #    all_ensembles: Every single ensemble member is shown.

        if plot_style == 'all_ensembles':
            if ensemble_1 in ['ensemble_mean', 'ensemble_min', 'ensemble_max']: continue
            if dataset_1 == 'CMIP6': continue
            #print(plot_style, exp_1, ensemble_1,dataset_1, plot_style)
            #assert 0

        if plot_style == 'all_models_means':
            if ensemble_1 != 'ensemble_mean': continue
            if dataset_1 == 'CMIP6': continue

        if plot_style == 'CMIP6_mean':
            if ensemble_1 != 'ensemble_mean': continue
            if dataset_1 != 'CMIP6': continue

        if skip_historical_ssp:
            if exp_1.find('historical-')>-1: continue
            if exp_1.find('historical_')>-1: continue

        if plot_style in ['CMIP6_range', 'all_models_range', ]:
            if plot_style in ['CMIP6_range', ] and dataset_1 != 'CMIP6': continue
            if plot_style in ['all_models_range', ] and dataset_1 == 'CMIP6': continue
            # just so it only gets plotted once.
            if ensemble_1 != 'ensemble_min': continue

            if x != 'time':
               print('plotting a range with time on the x axes is not possible.')
               assert 0

            if (dataset_1, y, exp_1, 'ensemble_min') not in data_dict.keys(): continue
            mins_data = data_dict[(dataset_1, y, exp_1, 'ensemble_min')]
            maxs_data = data_dict[(dataset_1, y, exp_1, 'ensemble_max')]
            x_label = 'Year'
            if isinstance(mins_data, iris.cube.Cube):
                x_times = diagtools.cube_time_to_float(mins_data)
                y_data_mins = np.array(mins_data.data.copy())
                y_data_maxs = np.array(maxs_data.data.copy())
            else:
                x_times = mins_data['time']
                y_data_mins = mins_data[y].copy()
                y_data_maxs = maxs_data[y].copy()

            x_times = np.array(x_times)
            if exp_1 in ['historical', ]:
                x_times = np.ma.masked_where(x_times>2015., x_times)
            else:
                x_times = np.ma.masked_where(x_times<2015., x_times)

            y_data_mins = np.ma.masked_where(x_times.mask, y_data_mins)
            y_data_maxs = np.ma.masked_where(x_times.mask, y_data_maxs)



            if y in label_dicts.keys():
                y_label = label_dicts[y]
            else:
                y_label = ' '.join([get_long_name(y) ])#, get_long_name(str(mins_data.units))])

            if np.array_equal(y_data_mins, y_data_maxs) or np.abs(np.mean(y_data_mins - y_data_maxs)/np.mean(y_data_maxs)) < 1e-6:
                plt.plot(x_times, y_data_mins, c=exp_colours[exp_1], alpha=0.4, lw=2.5)
            else:
                plt.fill_between(x_times, y_data_mins, y_data_maxs, fc = exp_colours_fill[exp_1], alpha=exp_colours_alpha[exp_1])
                plt.plot(x_times, y_data_mins, c=exp_colours[exp_1], lw=0.4)
                plt.plot(x_times, y_data_maxs, c=exp_colours[exp_1], lw=0.4)
            continue

        for (dataset, short_name, exp, ensemble), cube in data_dict.items():
            if short_name not in [x,y]: continue
            if exp != exp_1: continue
            if ensemble != ensemble_1: continue
            if dataset != dataset_1: continue

            print('Everything matches', plot_style, (dataset, short_name, exp, ensemble),'vs', [x,y], (exp_1, ensemble_1))
            print('make_ts_figure: found', plot_style, dataset, short_name, exp, ensemble, x,y)
            if x == 'time' and short_name == y:
                x_label = 'Year'
                if isinstance(cube, iris.cube.Cube):
                    x_data = diagtools.cube_time_to_float(cube)
                    x_times = x_data.copy()
                    print('setting x time to ',short_name, exp, ensemble)
                else:
                    x_data = cube['time']
                    x_times = x_data.copy()
                    print('setting x time to ',short_name, exp, ensemble)
            elif x == short_name and x in label_dicts.keys():
                x_data = cube[x].copy()
                x_times = cube['time'].copy()
                print('setting x time to ',short_name, exp, ensemble)
                x_label = label_dicts[x]

            elif short_name == x:
                x_data = np.array(cube.data.copy())
                x_times = diagtools.cube_time_to_float(cube)
                print('setting x axis to ',short_name, exp, ensemble, np.min(x_data), np.max(x_data))
                x_label = ' '.join([get_long_name(x), ] ) # get_long_name(str(cube.units))])

            if y == 'time':
                print('what kind of crazy person plots time on y axis?')
                assert 0
            elif y == short_name and y in label_dicts.keys():
                #['co2', 'emissions', 'cumul_emissions', 'luegt']:
                y_data = cube[y].copy()
                y_times = cube['time'].copy()
                y_label = label_dicts[y]
                print('setting y time to ',short_name, exp, ensemble)
            elif short_name == y:
                print(short_name, 'is a cube for y ts plot..')
                y_data = cube.data.copy()
                y_times = diagtools.cube_time_to_float(cube)
                print('setting y time to ',short_name, exp, ensemble, y_data)
                y_label = ' '.join([get_long_name(y), ] ) # get_long_name(str(cube.units))])

            print('make_ts_figure: loaded x data', short_name, exp, ensemble, x, np.mean(x_data))
            print('make_ts_figure: loaded y data', short_name, exp, ensemble, y, np.mean(y_data))

        if 0 in [len(x_data), len(y_data), len(x_times), len(y_times)]:
            print('no data found',(exp_1, ensemble_1,dataset_1),  x,'vs',y, 'x:', len(x_data), 'y:',len(y_data))
            continue

        if len(x_data) != len(x_times) or len(y_data) != len(y_times):
            print('x:', len(x_data), len(x_times), 'y:', len(y_data), len(y_times))
            assert 0

        label = ' '.join([exp_1, ]) #ensemble_1])
        # masks fromn the year 2005 of hist data, so that they can line up properly.

        hist_ssp_sync = False #
        if draw_line:
            x_times = np.ma.array(x_times)
            y_times = np.ma.array(y_times)
            number_of_lines+=1

            #plot_styles: ambition:
            #    ['ensemble_mean', ] default behaviour before
            #    CMIP6_mean: Only CMIP6 multi model mean
            #    CMIP6_range: range bewteen each model mean shown
            #    CMIP6_full_range: full range between individual model ensemble menmbers
            #    all_models_means: Each individual models mean is plotted.
            #    all_models_range: Each individual models range is plotted
            #    all_ensembles: Every single ensemble member is shown.

            if plot_styles == ['ensemble_mean', ]:
                lw = 1.3
            elif plot_style=='CMIP6_mean':
                lw = 1.7
            elif plot_style=='all_models_means':
                lw = 1.0
            else:
                lw = 0.5


            if exp_1 in ['historical', ]:
                x_times = np.ma.masked_where(x_times>2015., x_times)
                y_times = np.ma.masked_where(y_times>2015., y_times)
            else:
                x_times = np.ma.masked_where(x_times<2015., x_times)
                y_times = np.ma.masked_where(y_times<2015., y_times)

            x_data = np.ma.masked_where(x_times.mask, x_data)
            y_data = np.ma.masked_where(y_times.mask, y_data)
            plt.plot(x_data, y_data,
                     lw=lw,
                     color=exp_colours[exp_1])


#                if hist_ssp_sync:
#                    histx_t = np.ma.masked_where(x_times > 2005., x_times)
#                    histy_t = np.ma.masked_where(y_times > 2005., y_times)
#                    histx_d = np.ma.masked_where(histx_t.mask, x_data).compressed()
#                    histy_d = np.ma.masked_where(histy_t.mask, y_data).compressed()
#             plt.plot(histx_d, histy_d,
#                      lw=lw,
#                     color=exp_colours[exp_1])
#
#                else:
#             plt.plot(x_data, y_data,
#                        lw=lw,
#                        color=exp_colours[exp_1])
#
#            else:
#                if hist_ssp_sync:
#                    tdatcx = np.ma.masked_where((2004. > x_times) + (x_times > 2015.), x_times).compressed()
#                    tdatcy = np.ma.masked_where((2004. > y_times) + (y_times > 2015.), y_times).compressed()
#                    datcx = np.ma.masked_where((2004. > x_times) + (x_times > 2015.), x_data).compressed()
#                    datcy = np.ma.masked_where((2004. > y_times) + (y_times > 2015.), y_data).compressed()
#                    if len(tdatcx) == len(tdatcy):
#                        plt.plot(
#                            datcx, # np.ma.masked_where((2004 > x_times) + (x_times > 2015), x_data).compressed(),
#                            datcy, # np.ma.masked_where((2004 > y_times) + (y_times > 2015), y_data).compressed(),
#                            lw=lw,
#                            color=exp_colours['historical'])
#
#                    xdatc = np.ma.masked_where((x_times < 2015.) + (x_times > 2100.), x_data).compressed()
#                    ydatc = np.ma.masked_where((y_times < 2015.) + (y_times > 2100. ), y_data).compressed()
#                    xtdatc = np.ma.masked_where((x_times < 2015.) + (x_times > 2100.), x_times).compressed()
#                    ytdatc = np.ma.masked_where((y_times < 2015.) + (y_times > 2100. ), y_times).compressed()
#
#                    if len(xtdatc) == len(ytdatc):
#                        plt.plot(xdatc, # np.ma.masked_where(x_times < 2015., x_data).compressed(),
#                             ydatc, # np.ma.masked_where(y_times < 2015., y_data).compressed(),
#                             lw=lw,
#                             color=exp_colours[exp_1])
#                else:
#                    if len(x_data) != len(y_data):
#                        print('WARNING: x!=y:', len(x_data), '!=', len(y_data), 'x:', x, 'y:',y)
#                        print(x, 'x_times:', x_times)
#                        print(y, 'y_times:', y_times)
#
#                    plt.plot(x_data,
#                         y_data,
#                         lw=lw,
#                         color=exp_colours[exp_1])

        # plot_style == ''

        if markers == 'thresholds':
            try: threshold_times = thresholds_dict[(dataset_1, 'tas', exp_1, ensemble_1)]
            except:
               threshold_times = {}
            ms = 6
            for threshold, time in threshold_times.items():
                if not time:
                    continue
                if threshold not in plot_thresholds:
                    continue
                x_point = get_threshold_point({'time':x_times}, time.year)
                #_point = get_threshold_point(cube, time.year)
                y_point = get_threshold_point({'time':y_times}, time.year)

                print('thresholds:', dataset_1, exp_1, ensemble_1,x,y, threshold, time, x_point, y_point, len(x_data),len(y_data))
                #assert 0
                plt.plot(x_data[x_point],
                         y_data[y_point],
                         marker_styles[threshold],
                         markersize = ms,
                         fillstyle='full',
                         zorder=10,
                         color=to_rgba(exp_colours[exp_1], alpha=0.5),
                         markeredgecolor=exp_colours_dark[exp_1]
                          )
                #plt.plot(x_data[x_point],
                #         y_data[y_point],
                #         'o',
                #         markersize = 2,
                #         #fillstyle='none',
                #         color=exp_colours[exp_1])

#      if plot_style == '' ## fills:



    if not number_of_lines:
        print('No lines plotted')
        plt.close()
        return fig, ax
        #assert 0

    if do_legend:
        exp_colours_leg = {'historical':'black',
                       'ssp119':'green',
                       'ssp126':'dodgerblue',
                       'ssp245':'blue',
                       'ssp370':'purple',
                       'ssp434':'magenta',
                       'ssp585': 'red',
                       'ssp534-over':'orange'}
        plot_details = {}
        for exp,color in sorted(exp_colours_leg.items()):
            plot_details[exp] = {
                        'c': color,
                        'ls': '-',
                        'lw': 2.,
                        'label': exp
                    }
        for thres,ms in sorted(marker_styles.items()):
            plot_details[str(thres)] = {
                        'c': 'black',
                        'marker': ms,
                        'fillstyle':'none',
                        'label': '>' + str(thres)+u'\u00B0C'
                    }
        diagtools.add_legend_outside_right(
                plot_details, plt.gca(), column_width=0.175)

    # set labels:
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    # set title:
    if x == 'time':
        title = '' #get_long_name(y)
    else:
        title = ' '.join([get_long_name(x), 'by', get_long_name(y)])

    if make_figure_here and plot_dataset != 'all_models':
        title += ' '+plot_dataset

    plt.title(title)

    if make_figure_here:
        print('saving figure:', path)
        plt.savefig(path)
        plt.close()
    else:
        return fig, ax


def prepare_percentages_data( cfg,
    data_dict,
    thresholds_dict,
    threshold = '2.0',
    land_carbon = 'tls',
    #ensemble_key = 'all',
    ):
    #print("I think the problem is somewhere in here.")
    #assert 0


    #print(thresholds_dict)
    #print(119, thresholds_dict.get(('CMIP6', 'tas', 'ssp119', 'ensemble_mean'), None))
    #print(126, thresholds_dict[('CMIP6', 'tas', 'ssp126', 'ensemble_mean')])
    #print(245, thresholds_dict[('CMIP6', 'tas', 'ssp245', 'ensemble_mean')])
    #print(370, thresholds_dict[('CMIP6', 'tas', 'ssp370', 'ensemble_mean')])
    #print(585, thresholds_dict[('CMIP6', 'tas', 'ssp585', 'ensemble_mean')])

    data_dict_shelve = diagtools.folder([cfg['work_dir'], 'percentages_dicts'])
    data_dict_shelve+='_'.join(['allocations', threshold])+'.shelve'
    overwrite_shelve=True

    if not overwrite_shelve and glob.glob(data_dict_shelve+'*'):
        print('loading:', data_dict_shelve)
        sh = shelve.open(data_dict_shelve)
        remnants = sh['remnants']
        landcs = sh['landcs']
        fgco2gts = sh['fgco2gts']
        sh.close()

        # (t_dataset, t_exp, t_ens, threshold)
        print(remnants.keys())
        for k in sorted(remnants.keys()):
          if 'CMIP6' in k and 'ensemble_mean' in k:
              print(k, remnants[k])

        if threshold == '2.0':
           print(585, remnants[('CMIP6', 'SSP585', 'ensemble_mean', '2.0')])
           print(126, remnants[('CMIP6', 'SSP126', 'ensemble_mean', '2.0')])

        return remnants, landcs, fgco2gts

    #emissions = {}
    #threshold_years = {}
    remnants = {}
    landcs = {}
    fgco2gts = {}

    for (t_dataset, t_short, t_exp, t_ens), threshold_times in thresholds_dict.items():
        if t_short != 'tas': continue
        t_exp = standardized_exps(t_exp)
        t_ens = standardized_ens(t_ens)
        if t_exp == 'ssp534-over': continue
        print((t_dataset, t_short, t_exp, t_ens), threshold_times)
        #if t_dataset != 'CMIP6': continue
        #if t_ens != 'ensemble_mean': continue
        #if ensemble_key == 'all': pass
        #if ensemble_key = 'ensemble_mean' and t_ens != 'ensemble_mean': continue

        print("prepare_percentages_data", t_short, t_exp, t_ens)
#         cumul_emissions = data_dict.get((t_dataset, 'cumul_emissions', t_exp, t_ens), None) #dict
#         if cumul_emissions is None:  # Because not all sceanrios have emissions data.
#            print('couldnt find cumul_emissions:', (t_dataset, 'cumul_emissions', t_exp, t_ens))
#            for (t_dataset1, t_short1, t_exp1, t_ens1), datad in data_dict.items():
#                 if t_dataset1 != t_dataset: continue
#                 if t_short1 != 'cumul_emissions': continue
#                 print('candidate:', (t_dataset1, t_short1, t_exp1, t_ens1))
# #           for (t_dataset1, t_short1, t_exp1, t_ens1), datad in data_dict.items():
# #                if t_dataset1 != t_dataset: continue
# #                #if t_short1 != 'cumul_emissions': continue
# #                print('candidate 1:', (t_dataset1, t_short1, t_exp1, t_ens1))
#
#            assert 0

        atmos_carbon = data_dict.get((t_dataset, 'atmos_carbon',  t_exp, t_ens), None) #dict
        if atmos_carbon is None:  # Because not all sceanrios have emissions data.
           print('couldnt find atmos carbon:', t_dataset, 'atmos_carbon',  t_exp, t_ens)
           for (t_dataset1, t_short1, t_exp1, t_ens1), threshold_times in thresholds_dict.items():
                if t_dataset1 != t_dataset: continue
                if t_short1 != 'atmos_carbon': continue
                print('candidate:', (t_dataset, 'atmos_carbon',  t_exp, t_ens))
           assert 0

        # if land_carbon == 'nbpgt':
        #     landc_cumul = data_dict[(t_dataset, 'nbpgt_cumul', t_exp, t_ens)] # cube
        if land_carbon == 'tls':
            landc_cumul = data_dict.get((t_dataset, 'tls', t_exp, t_ens), None)# dict
        if landc_cumul is None:
           print('couldnt find land carbon:', t_dataset, 'tls:',  t_exp, t_ens)
           for (t_dataset1, t_short1, t_exp1, t_ens1), threshold_times in data_dict.items():
                if t_dataset1 != t_dataset: continue
                if t_short1 != 'tls': continue
                print('candidate:', (t_dataset, 'tls',  t_exp, t_ens))
           assert 0


        fgco2gt_cumul = data_dict[(t_dataset, 'fgco2gt_cumul', t_exp, t_ens)] # cube
        print('prepare_percentages_data: vworking here right now', threshold_times)
        for thresh, time in threshold_times.items():
            # print("prepare_percentages_data", t_short, t_exp, t_ens, thresh, threshold, time)
            if float(threshold) != float(thresh):
                #print(threshold, '!=', thresh)
                continue
            if not time:
                print('time', 'isn t there' , time)
                continue
            print("prepare_percentages_data", t_short, t_exp, t_ens, thresh, threshold, 'time:',time)

            fl_threshold = float(threshold)
            if fl_threshold > 1850.: # Threshold is a specific point in time.
                # e_xpoint = get_threshold_point(cumul_emissions, fl_threshold)
                a_xpoint = get_threshold_point(atmos_carbon, fl_threshold)
                n_xpoint = get_threshold_point(landc_cumul, fl_threshold)
                f_xpoint = get_threshold_point(fgco2gt_cumul, fl_threshold)
                if None in [a_xpoint, n_xpoint, f_xpoint]: continue

            print("prepare_percentages_data",threshold, time)
            # e_xpoint = get_threshold_point(cumul_emissions, time.year)
            a_xpoint = get_threshold_point(atmos_carbon, time.year)
            n_xpoint = get_threshold_point(landc_cumul, time.year)
            f_xpoint = get_threshold_point(fgco2gt_cumul, time.year)

            # emission = cumul_emissions['cumul_emissions'][e_xpoint]
            remnant = atmos_carbon['atmos_carbon'][a_xpoint]
            fgco2gt = fgco2gt_cumul.data[f_xpoint]

            if isinstance(landc_cumul, dict):
                landc = landc_cumul[land_carbon][n_xpoint]
            else:
                landc = landc_cumul.data[n_xpoint]
            t_exp = t_exp.replace('historical-', '').upper()

#            if t_dataset.find('UKESM')>-1 and  threshold == '4.0':
#                if landc+fgco2gt+ remnant<600.:
#                    print('ERROR:', t_short, t_exp, t_ens, thresh, threshold, time, (landc,fgco2gt, remnant))
#                    assert 0
            #if t_ens
            unique_key = (t_dataset, t_exp, t_ens, threshold)
            # emissions[unique_key] = emission
            remnants[unique_key] = remnant
            fgco2gts[unique_key] = fgco2gt
            landcs[unique_key] = landc

    print('prepare_percentages_data: saving:', data_dict_shelve)
    sh = shelve.open(data_dict_shelve)
    sh['remnants'] = remnants
    sh['landcs'] = landcs
    sh['fgco2gts'] = fgco2gts
    sh.close()

    return remnants, landcs, fgco2gts

def load_ecs_data():
    """
     Load zelinka data.
    """
    ECS_data = {}
    #ERF_data = {}
    with open('zelinka_ecs.csv', 'r') as file:
        reader = csv.reader(file)
        # MODEL, ECS, ERF
        for each_row in reader:
            print(each_row)
            if not len(each_row): continue
            if each_row[0] == '#': continue
            if each_row[0][0] == '#': continue

            if each_row[0] == 'MODEL': continue
            ECS_data[each_row[0]] = float(each_row[1])
            #ERF_data[each_row[0]] = float(each_row[2])
    return ECS_data

latexendline = lambda a: ''.join([a, '\\\\', '\n'])
latexhline = lambda a: ''.join([a, '\hline', '\n'])
latexrow = lambda a: ' & '.join(a)

def make_count_and_sensitivity_table(cfg, data_dict, thresholds_dict ):
    """
    Need two tables:
        One for the
            Model, historical, ssp.... ERT

        One showing the models that hit Each threshold:
            Model, ssp1, ssp
            2, 3, 4
            ukesm, 2,3,4 ; 2,3, , x
    This is in latex format.
    """
    short_name1='tas'
    experiments = ['historical', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
    table_data = {}
    datasets = {}

    for (dataset, short_name, exp, ensemble),cube  in data_dict.items():
        if exp not in experiments: continue
        if short_name!= short_name1: continue
        if ensemble in ['ensemble_mean','ensemble_min', 'ensemble_max' ]: continue
        if not cube.data.mean(): continue
        if mod_exp_ens_skips.get((dataset, exp, ensemble), False):
            continue
        if exp == 'ssp534-over': continue

        try: table_data[(dataset, exp)] +=1
        except: table_data[(dataset, exp)] =1
        datasets[dataset] = True

    # Load zelinka data.
    ECS_data = load_ecs_data()
#    #ERF_data = {}
#    with open('zelinka_ecs.csv', 'r') as file:
#        reader = csv.reader(file)
#        # MODEL, ECS, ERF
#        for each_row in reader:
#            print(each_row)
#            if not len(each_row): continue
#            if each_row[0] == '#': continue
#            if each_row[0][0] == '#': continue
#
#            if each_row[0] == 'MODEL': continue
#            ECS_data[each_row[0]] = float(each_row[1])
#            #ERF_data[each_row[0]] = float(each_row[2])

    header_list = ['Model', ]

    header_list += [sspify(exp) for exp in experiments]
    header_list += ['ECS', ]#'ERF']
    print(header_list)
    header = latexrow(header_list)
    header = latexendline(header)

    txt = '    \\begin{tabular}{|l|cccccc|c|}\n'
    txt = latexhline(txt)
    txt = ''.join([txt, header, '\hline \n'])
    for dataset in sorted(datasets.keys()):
        row = [dataset, ]
        for exp in experiments:
            row.append(str(table_data.get((dataset, exp), ' ')))
        ecs_val = str(ECS_data.get(dataset, False))
        if not ecs_val:
            print('Cant find ECS:', dataset, ecs_val)
            assert 0
        row.append(str(ECS_data.get(dataset, '--')))
        #row.append(str(ERF_data.get(dataset, '--')))
        row_str = latexrow(row)
        row_str= latexendline(row_str)
        txt = ''.join([txt, row_str])
        print(row_str)
    txt = latexhline(txt)

    # now Calculating totals:
    ens_total_row = ['Total number of Ensembles', ]
    total_row = ['Total number of Models', ]
    for exp in experiments:
        model_count = 0
        ens_count = 0
        for dataset in sorted(datasets.keys()):
            ens = table_data.get((dataset, exp), 0)
            if ens:
                model_count += 1
                ens_count += ens
        ens_total_row.append(str(ens_count))
        total_row.append(str(model_count))
    ens_total_row.append(' \\\\ \n')
    total_row.append(' \\\\ \n')
    ens_total_row = latexrow(ens_total_row)
    total_row = latexrow(total_row)
    txt = ''.join([txt, ens_total_row, total_row,])

    txt = latexhline(txt)

    # Now calculated the weighted ECS, ERF:
    ecs_mean_row = ['Weighted ECS',]
    for exp in experiments:
        weighted_ecs = []
        for dataset in sorted(datasets.keys()):
            a = table_data.get((dataset, exp), 0.)
            if not a: continue
            print('dataset:', dataset, 'in', ECS_data.keys(), (dataset in ECS_data))
            weighted_ecs.append(ECS_data.get(dataset, '--'))
            #weight_erf.append(ERF_data.get(dataset, '--'))

        print(exp, 'weighted_ecs', weighted_ecs)
        print(np.mean(weighted_ecs))
        ecs_mean_row.append(str(round(np.mean(weighted_ecs),2)))
    ecs_mean_row.append(' \\\\ \n')
    ecs_mean_row = latexrow( ecs_mean_row)
    txt = ''.join([txt, ecs_mean_row])
    txt = latexhline(txt)
    txt = ''.join([txt, '% Table generated by make_count_and_sensitivity_table function'])

        #print(exp, 'weight_erf', weight_erf, np.mean(weight_erf))


    print(txt)
    csvpath = diagtools.folder([cfg['plot_dir'], 'Table_1' ])
    csvpath += '_'.join(['latex_table1']) + '.txt'

    csv_file = open(csvpath,'w')
    csv_file.write(txt)
    csv_file.close()



def make_ensemble_barchart_pane(
    cfg,
    data_dict,
    thresholds_dict,
    threshold = '2.0',
    year0 = None,
    do_legend=True,
    land_carbon = 'nbpgt',
    #atmos='atmos_carbon',
    ensemble_key = 'ensemble_mean',
    plot_style = 'percentages',
    group_by = 'group_by_model',
    sorting='alphabetical',
    fig=None,
    ax=None):
    """
    Make a bar chart (of my favourite pies)

    Ensemble plot
    """
    if fig == None:
        fig, ax = plt.subplots()
        make_figure_here = True
    else:
        make_figure_here = False
        plt.sca(ax)

    # single pane, single threshold
    remnants, landcs, fgco2gts = prepare_percentages_data(cfg,
        data_dict,
        thresholds_dict,
        threshold = threshold,
        land_carbon = 'tls')

    datasets, exps, ensembles, thresholds = {}, {}, {}, {}
    # only_one_ensemble = {}
    ecs_keys = {}
    for unique_key, remnant in remnants.items():
        (t_dataset, t_exp, t_ens, t_threshold) = unique_key
        datasets[t_dataset] = True
        exps[t_exp] = True
        ensembles[t_ens] = True
        thresholds[t_threshold] = True
        ecs_keys[(t_dataset, t_exp)] = True

    exps = sorted(list(exps.keys()))
    thresholds = sorted(list(thresholds.keys()))

    ensembles = sorted(list(ensembles.keys()))
    ensembles.insert(0,  ensembles.pop(ensembles.index('ensemble_mean')))

    datasets = sorted(list(datasets.keys()))
    datasets.insert(0,  datasets.pop(datasets.index('CMIP6')))

    # Build a list of keys.
    unique_key_order = []
    # order is dataset [CMIP then others], ssps, ensemble means, then ensemble members
    # Order by model:
    if group_by == 'ecs':

        ECS_data = load_ecs_data()
        ECS_data_ssp = {}
        for exp in exps:
            ECS_data_ssp[exp] = {'CMIP6': []}
            for dat in datasets:
                if dat=='CMIP6': continue
                ecs = ECS_data.get(dat, False)
                if not ecs: continue
                ECS_data_ssp[exp]['CMIP6'].append(ecs)
                ECS_data_ssp[exp][dat] = ecs
            ECS_data_ssp[exp]['CMIP6'] = np.mean(ECS_data_ssp[exp]['CMIP6'])

        print('ECS_data_ssp:', ECS_data_ssp)
        for thr in thresholds:
            for exp in exps:
                new_dataset_order = sorted((value, key) for (key,value) in ECS_data_ssp[exp].items())
                new_dataset_order=[k for v, k in new_dataset_order]

#               new_dataset_order = {d:ecs for (d, e), ecs in ECS_data_ssp.items() if e == exp}
                print(thr, exp, 'new_dataset_order:', new_dataset_order)
 #              new_dataset_order = sorted(new_dataset_order.items(), key=lambda x:x[1])
                #print(ECS_data_ssp)
                #assert 0
                #dataset_order = sortedECS_data_ssp.items()
                #or dset,ecs in new_dataset_order.items():

                for dset in new_dataset_order:
                    print(dset,ecs)
                    for ens in ensembles:
                        if ensemble_key == 'ensemble_mean' and ens != 'ensemble_mean':
                            continue
                        if ensemble_key != 'ensemble_mean' and ens == 'ensemble_mean' and dset!= 'CMIP6': continue
                        if dset.find('UKESM')>-1: print(dset, exp, ens, thr)
                        unique_key = (dset, exp, ens, thr)
                        if unique_key in unique_key_order: continue
                        if unique_key not in remnants: continue
                        print('Key added', unique_key)
                        unique_key_order.append(unique_key)

    if group_by == 'group_by_model':
        for thr in thresholds:
            for dset in datasets:
               for exp in exps:
                    for ens in ensembles:
                        if ensemble_key == 'ensemble_mean' and ens != 'ensemble_mean':
                            continue
                        if ensemble_key != 'ensemble_mean' and ens == 'ensemble_mean' and dset!= 'CMIP6': continue
                        if dset.find('UKESM')>-1: print(dset, exp, ens, thr)
                        unique_key = (dset, exp, ens, thr)
                        if unique_key in unique_key_order: continue
                        if unique_key not in remnants: continue
                        print('Key added', unique_key)
                        unique_key_order.append(unique_key)

    if group_by == 'group_by_ssp':
        for thr in thresholds:
            for exp in exps:
                for dset in datasets:
                    for ens in ensembles:
                        if ensemble_key == 'ensemble_mean' and ens != 'ensemble_mean':
                            continue
                        if ensemble_key != 'ensemble_mean' and ens == 'ensemble_mean' and dset!= 'CMIP6': continue
                        if dset.find('UKESM')>-1: print(dset, exp, ens, thr)
                        unique_key = (dset, exp, ens, thr)
                        if unique_key in unique_key_order: continue
                        if unique_key not in remnants: continue
                        print('Key added', unique_key)
                        unique_key_order.append(unique_key)

    #totals={}
    labels, land, air, ocean = [], [], [], []
    xvalues = []
    widths = []
    emissions_bottoms = []

    previous_datasets = {}
#    if ensemble_key == 'ensemble_mean':
#        gap = 0.3
#    else: gap = 0.6
    quit = False

    if ensemble_key == 'ensemble_mean':
        print(group_by, unique_key_order)
        #assert 0

    make_table_here = True
    if make_table_here:
        csvpath = diagtools.folder([cfg['plot_dir'], 'ensemble_barcharts_csv' ])
        csvpath += '_'.join(['ensemble_barchart'+str(threshold), land_carbon, ensemble_key]) + '.txt'
        out_txt = 'Count, model, exp, ensemble, threshold, atmos, land, ocean, total, \n'
        for i, unique_key in enumerate(unique_key_order):

            if ensemble_key == 'ensemble_mean' and t_ens != 'ensemble_mean': continue
            if t_threshold != threshold: continue
            remnant = remnants[unique_key]
            landc = landcs[unique_key]
            oceanc = fgco2gts[unique_key]
            total = remnant + landc + oceanc

            (t_dataset, t_exp, t_ens, t_threshold) = unique_key
            line1 = ', '.join([str(i), t_dataset, t_exp, t_ens, t_threshold])
            line2 = ', '.join([str(v) for v in [remnant, landc, oceanc, total, '\n']])
            out_txt  += line1 +', '+line2

        csv_file = open(csvpath,'w')
        csv_file.write(out_txt)
        csv_file.close()

        csvpath = diagtools.folder([cfg['plot_dir'], 'ensemble_counts_csv' ])
        csvpath += '_'.join(['ensemble_barchart_count'+str(threshold), land_carbon, ensemble_key]) + '.txt'
        out_txt = '# Count for thrshold: '+str(threshold)+'\n'
        out_txt += '# land_carbon:' + land_carbon+'\n'
        out_txt += '# ensemble_key:' + ensemble_key+'\n'
        #out_txt += 'Model, scenario, count,\n'
        counts = {}
        models = {}
        ssps = {}
        for i, unique_key in enumerate(unique_key_order):

            if ensemble_key == 'ensemble_mean' and t_ens != 'ensemble_mean': continue
            if t_threshold != threshold: continue
            (t_dataset, t_exp, t_ens, t_threshold) = unique_key
            models[t_dataset] = True
            ssps[t_exp] = True
            count_key = (t_dataset, t_exp)
            if counts.get(count_key, False):
                counts[count_key]+=1
            else:
                counts[count_key] = 1
        out_txt  = 'Model, ' + ', '.join(sorted(ssps.keys()))+'\n'
        for moded_csv in sorted(models.keys()):
           line1 = moded_csv+', '
           for exp_csv in sorted(ssps.keys()):
               if (moded_csv, exp_csv) in counts:

                   line1 += str(counts[(moded_csv, exp_csv)])+', '
               else:
                   line1 += '0, '
           out_txt += line1 +'\n'
       #for (t_dataset, t_exp) in sorted(counts.keys()):
       #     line1 = ','.join([t_dataset, t_exp, str(counts[(t_dataset, t_exp)])])
       #     out_txt  += line1 +',\n'

        csv_file = open(csvpath,'w')
        csv_file.write(out_txt)
        csv_file.close()


    for i, unique_key in enumerate(unique_key_order):
    #for unique_key, remnant in sorted(remnants.items()):
        (t_dataset, t_exp, t_ens, t_threshold) = unique_key

        if ensemble_key == 'ensemble_mean' and t_ens != 'ensemble_mean': continue

        #if 'CMIP6' in unique_key: continue

        if t_threshold != threshold: continue
        #if unique_key not in remnants: continue

        remnant = remnants[unique_key]
        landc = landcs[unique_key]
        oceanc = fgco2gts[unique_key]
        total = remnant + landc + oceanc

        adding_gaps = False
        dataset_blank = (i>0 and group_by == 'group_by_model' and t_dataset not in labels[-1])
        exp_blank =  (i>0 and group_by == 'group_by_ssp' and t_exp not in labels[-1])

        if adding_gaps and (dataset_blank or exp_blank):
            xvalues.append((xvalues[-1]+0.4))
            widths.append(1.)
            land.append(0.)
            ocean.append(0.)
            air.append(0.)
            emissions_bottoms.append(0.)
            labels.append([''.join(['.' for k in range(i)]), ])


#        if i == 0:
#            print(i, unique_key)
#            assert 0
        # create bar label.
        #if ensemble_key == 'ensemble_mean':
        #    label_keys = [t_dataset, t_exp]
        #else:
        label_keys = [t_dataset, t_exp, t_ens]
 #       if 'CMIP6' in label_keys and 'SSP119' in t_exp and

        labels.append(label_keys)

        xvalues.append(i+0.5)
        widths.append(1.)

#        if i == 0:
#            xvalues.append(0.)
#            widths.append(1.)
#        else:
#            xvalues.append(xvalues[-1]+1.)
#            widths.append(1.)

        if plot_style == 'percentages':
            land.append(100. * landc/total)
            ocean.append(100. * oceanc/total)
            air.append(100. * remnant/total)
            emissions_bottoms.append(100* (landc +oceanc)/total)
        else:
            land.append(landc)
            ocean.append(oceanc)
            air.append(remnant)
            emissions_bottoms.append(landc+oceanc)
            if t_dataset.find('UKESM')>-1 and  threshold == '4.0':
                if landc+oceanc + remnant<600.:
                    print('ERROR:', label_keys, landc, oceanc, remnant, landc+oceanc + remnant )
                    quit = True
                    assert 0

        # Adds a blank line between models.
        #if t_dataset not in previous_dats:
        #    land.append(0.)
        #    ocean.append(0.)
        #    air.append(0.)
        #    emissions_bottoms.append(0.)
        #    experiments.append('')
        #previous_dats[t_dataset] = True

        #emissions_bottoms[unique_key]
    #totals = [a + f + b for a,f,b in zip(remnant, fgco2gts, landcs)]

    # if atmos=='atmos_carbon':
    #     emissions_diff = remnant
    # else:
    #     emissions_diff = [e - f - b for e,f,b in zip(emissions, fgco2gts, landcs )]
    #emissions_diff = remnant
    # emissions_bottoms = [f + b for f,b in zip(fgco2gts, landcs )]
    # totals = [a + f + b for a,f,b in zip(remnant, fgco2gts, landcs)]
    #
    # for i, exp  in enumerate(experiments):
    #     print("Final values:",threshold, i, exp, totals[i], '=',emissions_diff[i],('or', remnant[i]), '+', landcs[i], '+', fgco2gts[i], '(a, l, o)')
    # Add bars:
    #if quit: assert 0
    #for x,w,l in zip(xvalues, widths, labels): print(x,w,l)

    label_strs = []
    linewidths = []
    edge_colours = []
    colours_land, colours_ocean, colours_air = [],[],[]

    ssp_land = {
        #'historical':'black',
        'SSP119': 'darkgreen',
        'SSP126': 'blue',
        'SSP245': 'saddlebrown',
        'SSP370': 'indigo',
        'SSP585': 'darkred',
        }
    ssp_ocean = {
        'SSP119': 'green',
        'SSP126': 'royalblue',
        'SSP245': 'darkorange',
        'SSP370': 'purple',
        'SSP585': 'red',
    }
    ssp_air = {
        'SSP119': 'lightgreen',
        'SSP126': 'dodgerblue',
        'SSP245': 'sandybrown',
        'SSP370': 'orchid',
        'SSP585': 'lightcoral',
    }
    full_alpha = (0.,0.,0.,0.)
    for i, lablist in enumerate(labels):

        if lablist[0][0] == '.':
            # empty bar
            label_strs.append(' ')
            linewidths.append(0.)
            edge_colours.append(full_alpha)
            colours_land.append('white')
            colours_ocean.append('white')
            colours_air.append('white')
            continue

        (t_dataset, t_exp, t_ens) = lablist

        colours_land.append(ssp_land[t_exp])
        colours_ocean.append(ssp_ocean[t_exp])
        colours_air.append(ssp_air[t_exp])

        if t_dataset == 'CMIP6':
            edge_colours.append('k')
            linewidths.append(0.)
        else:
            edge_colours.append(full_alpha)
            linewidths.append(0.)

        if ensemble_key == 'ensemble_mean':
            if group_by == 'group_by_model':
                label_strs.append(' '.join([t_dataset, t_exp]))

            if group_by == 'group_by_ssp':
                 if t_dataset == 'CMIP6':
                     label_strs.append('Multi-model mean')
                 else:
                     label_strs.append(' '.join([t_dataset, ]))
            if group_by == 'ecs':
                 ecs= str(round(ECS_data_ssp[t_exp][t_dataset],2))
                 if t_dataset == 'CMIP6':
                     label_strs.append(' '.join(['Multi-model mean', ])) # ecs
                 else:
                     label_strs.append(' '.join([t_dataset, ])) # ecs

        else:
            if group_by == 'group_by_model':
               if t_dataset  not in labels[i-1]:
                   label_strs.append(t_dataset)
               else: label_strs.append(' ')
            if group_by == 'group_by_ssp':
               if t_exp not in labels[i-1]:
                   label_strs.append(t_exp)
               else: label_strs.append(' ')

            if group_by == 'ecs':
               if t_exp not in labels[i-1]:
                   label_strs.append(t_exp)
               else: label_strs.append(' ')

    #labels = [' '.join(label) for label in labels]
    if len(label_strs)== len(land) == len(widths) == len(xvalues) == len(colours_land): pass
    else:
        print('ERROR:',  len(labels),len(land),len(widths), len(xvalues), len(label_strs), len(colours_land))
        assert 0


    ax.bar(xvalues, land, width=widths, label='Land', color=colours_land, tick_label = label_strs, edgecolor=edge_colours, linewidth=linewidths)
    ax.bar(xvalues, ocean, width=widths, bottom = land,  label='Ocean', color=colours_ocean, edgecolor=edge_colours, linewidth=linewidths)
    ax.bar(xvalues, air, width=widths, bottom = emissions_bottoms,  label='Atmos', color=colours_air, edgecolor=edge_colours, linewidth=linewidths)

    # ax.bar(xvalues, land, width=widths, label='Land', color='green', tick_label = label_strs)
    # ax.bar(xvalues, ocean, width=widths, bottom = land,  label='Ocean', color='dodgerblue')
    # ax.bar(xvalues, air, width=widths, bottom = emissions_bottoms,  label='Atmos', color='grey')
    #ax.set_xlabel('Scenarios')
    plt.xticks(rotation=90, fontsize='xx-small')
    plt.xlim([-0.1, np.sum(widths)+0.1])

    # # Add bars:
    # horizontal = False
    # if horizontal:
    #     ax.barh(experiments, land, label='Land', color='mediumseagreen')
    #     ax.barh(experiments, ocean, left = land,  label='Ocean', color='dodgerblue')
    #     ax.barh(experiments, air, left = emissions_bottoms,  label='Atmos', color='silver')
    #     #ax.set_ylabel('Scenarios')
    #     #lt.xticks(rotation=90)
    # else:
    #     ax.bar(experiments, land, width=1.0, label='Land', color='green')
    #     ax.bar(experiments, ocean, width=1.0, bottom = land,  label='Ocean', color='dodgerblue')
    #     ax.bar(experiments, air, width=1.0, bottom = emissions_bottoms,  label='Atmos', color='grey')
    #     #ax.set_xlabel('Scenarios')
    #     plt.xticks(rotation=90)
    # Add percentages:
    # add_pc_text = True
    # if add_pc_text:
    #     def pc(a,b): return str("{0:.1f}%".format(100.*a/b))
    #     for e, exp in enumerate(experiments):
    #         print(e, exp, landcs[e], fgco2gts[e], emissions_diff[e])
    #         t = totals[e]
    #         ax.text(landcs[e]/2, e, pc(landcs[e], t),
    #                 color='#002200', #'darkgreen',
    #                 fontsize=8 , # fontweight='bold',
    #                 verticalalignment='center',horizontalalignment='center')
    #
    #         ax.text(landcs[e]+fgco2gts[e]/2., e, pc(fgco2gts[e], t),
    #                 color='darkblue', fontsize=8 , # fontweight='bold',
    #                 verticalalignment='center',horizontalalignment='center')
    #
    #         ax.text(emissions_bottoms[e]+emissions_diff[e]/2., e, pc(emissions_diff[e],t ),
    #                 color='black', fontsize=8 , # fontweight='bold',
    #                 verticalalignment='center',horizontalalignment='center')

    # Add mean year:
    # add_mean_year = False
    # if add_mean_year:
    #     for e, exp in enumerate(experiments):
    #         t = totals[e]
    #         yr = threshold_years[e]
    #         ax.text(t *1.05, e, str(int(yr)),
    #                 color='black', fontsize=8 , # fontweight='bold',
    #                 verticalalignment='center',horizontalalignment='center')

    if float(threshold) > 1850.:
        ax.set_title(str(threshold))
    else:
        ax.set_title(str(threshold)+r'$\degree$'+' warming')

    if plot_style == 'percentages':
        ax.set_ylim([0., 100.,])
        ax.set_ylabel('Fractional Carbon Allocation, %')
    else:
        ax.set_ylabel('Total Carbon Allocation, Pg')


    if do_legend:
        ax.legend()

    if make_figure_here:
        image_extention = diagtools.get_image_format(cfg)
        path = diagtools.folder([cfg['plot_dir'], 'single_barcharts' ])
        path += '_'.join(['ensemble_barchart'+str(threshold), land_carbon, atmos, group_by]) + image_extention
        plt.savefig(path)
        plt.close()
    else:
        return fig, ax


def make_ensemble_barchart(
        cfg,
        data_dict,
        thresholds_dict,
        plot_style='percentages',
        ensemble_key = 'ensemble_mean',
        group_by = 'group_by_model',
        thresholds = ['4.0', '3.0', '2.0'],
#        sorting='alphabetical',
    ):
    """
    Make a barchat for the whole ensemble
    Vertical bars

    """
    fig = plt.figure()
    fig.set_size_inches(12 , 6)
    gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.230,height_ratios=[8.  , 1])
    ax_2=  fig.add_subplot(gs[0, 0])
    ax_3 =  fig.add_subplot(gs[0, 1])
    ax_4 =  fig.add_subplot(gs[0, 2])
    ax_leg = fig.add_subplot(gs[0, 3])

    axes = [ax_4, ax_3, ax_2]
    for ax, threshold in zip(axes, thresholds):
        make_ensemble_barchart_pane(cfg, data_dict, thresholds_dict,threshold = threshold,fig=fig, ax=ax, do_legend=False,
            plot_style= plot_style,
            ensemble_key = ensemble_key,
            group_by = group_by,
#            sorting=sorting,
            )
        plt.xticks(rotation=90)
        ax.tick_params(axis = 'x', labelsize = 'x-small')

        if ax in [ax_3, ax_4]:
            ax.set_yticks([])
            ax.set_yticklabels([])
            ax.set_ylabel('')
            ax.spines['left'].set_visible(False)

    xranges = []
    ranges = []
    for ax in [ax_2, ax_3, ax_4]:
        plt.sca(ax)
        ranges.append(ax.get_ylim())
        xranges.append(np.max(ax.get_xlim()) - np.min(ax.get_xlim()))

    for ax in [ax_2, ax_3, ax_4]:
        plt.sca(ax)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
#        if ax in [ax_3, ax_2]:
#            ax.spines['left'].set_visible(False)

        ax.set_ylim([np.min(ranges), np.max(ranges)])


    # Set x axis size to make all panes the same aspect ratio (maybe)
    #gs.update(width_ratios=(xr[1] for xr in xranges))

    print('New width ratios:', xranges)
    xranges.append(np.sum(xranges)/15.) # for legend.
    gs.set_width_ratios(xranges)
    fig.subplots_adjust(wspace=0.02 )

    # Adding legend.
    # from: https://stackoverflow.com/questions/31908982/python-matplotlib-multi-color-legend-entry
    """
    cmap_ssp_land = {
        #'historical':'black',
        'SSP119': 'darkgreen',
        'SSP126': 'blue',
        'SSP245': 'saddlebrown',
        'SSP370': 'indigo',
        'SSP585': 'darkred',
        }
    cmap_ssp_ocean = {
        'SSP119': 'green',
        'SSP126': 'royalblue',
        'SSP245': 'darkorange',
        'SSP370': 'purple',
        'SSP585': 'red',
    }
    cmap_ssp_air = {
        'SSP119': 'lightgreen',
        'SSP126': 'dodgerblue',
        'SSP245': 'sandybrown',
        'SSP370': 'orchid',
        'SSP585': 'lightcoral',
    }"""
    plt.sca(ax_leg)

    ms = 10
    mland, = ax.plot([], [], c='black', marker='s', markersize=ms,
              linestyle='none', markeredgecolor="black")
    mocean, = ax.plot([], [], c='grey', marker='s', markersize=ms,
              linestyle='none', markeredgecolor="black")
    mair, = ax.plot([], [], c='silver', marker='s', markersize=ms,
              linestyle='none', markeredgecolor="black")

    dummies = {}
    for ssp in cmap_ssp_air.keys():
        dummies[(ssp, 'land')], = ax.plot([], [], c=cmap_ssp_land[ssp], marker='s', markersize=ms,
              fillstyle='bottom',linestyle='none', markeredgecolor=cmap_ssp_ocean[ssp])
        #dummies[(ssp, 'ocean')], = ax.plot([], [], c=cmap_ssp_ocean[ssp], marker='s', markersize=20,
        #      fillstyle='cent', linestyle='none', markeredgecolor="black")
        dummies[(ssp, 'air')], = ax.plot([], [], c=cmap_ssp_air[ssp], marker='s', markersize=ms,
              fillstyle='top', linestyle='none', markeredgecolor=cmap_ssp_ocean[ssp])

    keys, labels = [], []
    for k, label in zip([mair , mocean, mland],['Atmosphere', 'Ocean', 'Land']):
        keys.append(k)
        labels.append(label)

    for ssp in cmap_ssp_air.keys():
        #keys.append((dummies[(ssp, 'land')], dummies[(ssp, 'ocean')], dummies[(ssp, 'sea')]))
        keys.append((dummies[(ssp, 'land')], dummies[(ssp, 'air')]))

        labels.append(sspify(ssp))

    legd = ax_leg.legend(keys, labels,
        bbox_to_anchor=(2.6, 0.5), #
        numpoints=1, labelspacing=1.2,
        loc='center right', ) #tsize=16)

    legd.draw_frame(False)
    legd.get_frame().set_alpha(0.)
    ax_leg.get_xaxis().set_visible(False)
    ax_leg.get_yaxis().set_visible(False)
    plt.axis('off')

    image_extention = diagtools.get_image_format(cfg)
    path = diagtools.folder([cfg['plot_dir'], 'ensemble_barcharts_new'])
    path += '_'.join(['ensemble_barcharts', plot_style, ensemble_key, group_by]) + image_extention
    print('Save image:', path)
    plt.savefig(path)
    plt.close()
    #assert 0


def make_bar_chart(cfg, data_dict, thresholds_dict, threshold = '2.0',
    year0 = None,
    do_legend=True,
    land_carbon = 'nbpgt',
    atmos='atmos_carbon',
    exp_order = ['historical', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'],
    plot_dataset='CMIP6',
    fig=None,
     ax=None):
    """
    Make a bar chart (of my favourite pies)
    """
    #emissions = []
    threshold_years = []
    remnant = []
    landcs = []
    fgco2gts = []
    experiments = []

    for exp1 in exp_order:
      for (t_dataset, t_short, t_exp, t_ens), threshold_times in thresholds_dict.items():
        # print((t_dataset, t_short, t_exp, t_ens), threshold_times)
        if t_dataset != plot_dataset: continue
        if t_short != 'tas': continue
        if t_ens != 'ensemble_mean': continue
        if t_exp != exp1: continue

        print("make_bar_chart", t_short, t_exp, t_ens)
        atmos_carbon = data_dict.get((t_dataset, atmos,  t_exp, t_ens), None) #dict

        if atmos_carbon is None:  # Because not all sceanrios have emissions data.
           print('Did not find atmoshs_carbon:', (t_dataset, atmos,  t_exp, t_ens))
           assert 0

        if land_carbon == 'nbpgt':
            landc_cumul = data_dict[(t_dataset, 'nbpgt_cumul', t_exp, t_ens)] # cube
        elif land_carbon == 'tls':
            landc_cumul = data_dict[(t_dataset, 'tls', t_exp, t_ens)] # dict

        fgco2gt_cumul = data_dict[(t_dataset, 'fgco2gt_cumul', t_exp, t_ens)] # cube
        fl_threshold = float(threshold)
        if fl_threshold > 1850.:

            # e_xpoint = get_threshold_point(cumul_emissions, fl_threshold)
            a_xpoint = get_threshold_point(atmos_carbon, fl_threshold)
            n_xpoint = get_threshold_point(landc_cumul, fl_threshold)
            f_xpoint = get_threshold_point(fgco2gt_cumul, fl_threshold)
            if None in [a_xpoint, n_xpoint, f_xpoint]: continue

            # emissions.append(cumul_emissions['cumul_emissions'][e_xpoint])
            remnant.append(atmos_carbon['atmos_carbon'][a_xpoint])
            if isinstance(landc_cumul, dict):
                landcs.append(landc_cumul[land_carbon][n_xpoint])
            else:
                landcs.append(landc_cumul.data[n_xpoint])
            fgco2gts.append(fgco2gt_cumul.data[f_xpoint])

            t_exp = t_exp.replace('historical-', '').upper()
            experiments.append(t_exp)
            continue
        print('threshold_times: vworking here right now', threshold_times)
        for thresh, time in threshold_times.items():
            print("make_bar_chart", t_short, t_exp, t_ens, thresh, threshold, time)
            if float(threshold) != float(thresh):
                print(threshold, '!=', thresh)
                continue
            if not time:
                print('time', 'isn t there' , time)
                continue
            #if threshold != thresh: continue
            print("make_bar_chart",threshold, time)

            # e_xpoint = get_threshold_point(cumul_emissions, time.year)
            a_xpoint = get_threshold_point(atmos_carbon, time.year)

            n_xpoint = get_threshold_point(landc_cumul, time.year)
            f_xpoint = get_threshold_point(fgco2gt_cumul, time.year)

            print("make_bar_chart",thresh, time, 'atmos_carbon', a_xpoint, 'land:', n_xpoint, 'ocean', f_xpoint)
            if year0:
                assert 0
                #e_baseline = get_threshold_point(cumul_emissions, year0)
                #n_baseline = get_threshold_point(landc_cumul, year0)
                #f_baseline = get_threshold_point(fgco2gt_cumul, year0)

               # emissions.append(cumul_emissions['cumul_emissions'][e_xpoint] - cumul_emissions['cumul_emissions'][e_baseline])
               # if isinstance(landc_cumul, dict):
               #     landcs.append(landc_cumul[land_carbon][n_xpoint] - landc_cumul[land_carbon][n_baseline])
               #     fgco2gts.append(fgco2gt_cumul.data[f_xpoint] - landc_cumul[land_carbon][f_baseline])
               # else:
               #     landcs.append(landc_cumul.data[n_xpoint] - landc_cumul.data[n_baseline])
               #     fgco2gts.append(fgco2gt_cumul.data[f_xpoint] - landc_cumul.data[f_baseline])
            else:
                # print("make_bar_chart",thresh, time,'cumul_emissions', cumul_emissions['cumul_emissions'][e_xpoint], cumul_emissions['time'][e_xpoint])
                print("make_bar_chart",thresh, time,'land:          ', landc_cumul[land_carbon][n_xpoint], landc_cumul['time'][n_xpoint])
                print("make_bar_chart",thresh, time,'ocean:         ', fgco2gt_cumul.data[f_xpoint], fgco2gt_cumul.coord('time').points[f_xpoint])

                # emissions.append(cumul_emissions['cumul_emissions'][e_xpoint])
                remnant.append(atmos_carbon['atmos_carbon'][a_xpoint])
                threshold_years.append(time.year)
                if isinstance(landc_cumul, dict):
                    landcs.append(landc_cumul[land_carbon][n_xpoint])
                else:
                    landcs.append(landc_cumul.data[n_xpoint])
                fgco2gts.append(fgco2gt_cumul.data[f_xpoint])

            t_exp = t_exp.replace('historical-', '').upper()
            experiments.append(t_exp)

    if fig == None:
        fig, ax = plt.subplots()
        make_figure_here = True
    else:
        make_figure_here = False

    if not len(experiments):
        print("make_bar_chart", landcs, fgco2gts, remnant,  experiments)
        print("make_bar_chart",thresholds_dict.keys())
        print("make_bar_chart",'looking for:', threshold)
        return fig, ax


    #experiments = [exp.replace('historical-', '').upper() for exp in experiments]
    if atmos=='atmos_carbon':
        emissions_diff = remnant
    else:
        assert 0
        # emissions_diff = [e - f - b for e,f,b in zip(emissions, fgco2gts, landcs )]
    emissions_bottoms = [f + b for f,b in zip(fgco2gts, landcs )]

    totals = [a + f + b for a,f,b in zip(remnant, fgco2gts, landcs)]

    for i, exp  in enumerate(experiments):
        print("Final values:",threshold, i, exp, totals[i], '=',emissions_diff[i],('or', remnant[i]), '+', landcs[i], '+', fgco2gts[i], '(a, l, o)')
    #assert 0
    # Add bars:
    horizontal = True
    if horizontal:
        ax.barh(experiments, landcs, label='Land', color='mediumseagreen')
        ax.barh(experiments, fgco2gts, left = landcs,  label='Ocean', color='dodgerblue')
        ax.barh(experiments, emissions_diff, left = emissions_bottoms,  label='Atmos', color='silver')
#        ax.set_ylabel(str(threshold)+r'$\degree$'+' warming')

    else:
        ax.barh(experiments, landcs, label='Land', color='green')
        ax.barh(experiments, fgco2gts, bottom = landcs,  label='Ocean', color='dodgerblue')
        ax.barh(experiments, emissions_diff, bottom = emissions_bottoms,  label='Atmos', color='grey')
#        ax.set_xlabel('Scenarios')

    # Add percentages:
    add_pc_text = True
    if add_pc_text:
        def pc(a,b): return str("{0:.1f}%".format(100.*a/b))
        for e, exp in enumerate(experiments):
            print(e, exp, landcs[e], fgco2gts[e], emissions_diff[e])
            t = totals[e]
            ax.text(landcs[e]/2, e, pc(landcs[e], t),
                    color='#002200', #'darkgreen',
                    fontsize=8 , # fontweight='bold',
                    verticalalignment='center',horizontalalignment='center')

            ax.text(landcs[e]+fgco2gts[e]/2., e, pc(fgco2gts[e], t),
                    color='darkblue', fontsize=8 , # fontweight='bold',
                    verticalalignment='center',horizontalalignment='center')

            ax.text(emissions_bottoms[e]+emissions_diff[e]/2., e, pc(emissions_diff[e],t ),
                    color='black', fontsize=8 , # fontweight='bold',
                    verticalalignment='center',horizontalalignment='center')

    # Add mean year:
    add_mean_year = False
    if add_mean_year:
        for e, exp in enumerate(experiments):
            t = totals[e]
            yr = threshold_years[e]
            ax.text(t *1.05, e, str(int(yr)),
                    color='black', fontsize=8 , # fontweight='bold',
                    verticalalignment='center',horizontalalignment='center')

    if float(threshold) > 1850.:
#        ax.set_title('Carbon Allocation at '+str(threshold))
        ax.set_ylabel(str(threshold)) #+r'$\degree$'+' warming')

    else:
        ax.set_ylabel(str(threshold)+r'$\degree$'+'C GWT')
#        ax.set_title('Carbon Allocation at '+str(threshold)+r'$\degree$'+' warming')
    ax.set_xlabel('Anthropogenic Carbon Allocation, Pg')

    if do_legend:
        ax.legend(loc='lower right')

    if make_figure_here:
        image_extention = diagtools.get_image_format(cfg)
        path = diagtools.folder([cfg['plot_dir'], 'barcharts' ])
        path += '_'.join(['barchart'+str(threshold), land_carbon, atmos]) + image_extention
        plt.savefig(path)
        print('Saved:', path)
        plt.close()
    else:
        return fig, ax




def zip_time(dict, key):
    """
    zip data_dict style dict into a time:value style dict.
    """
    out_dict = {t:d for t,d in zip(dict['time'], dict[key]) if t < 2100.}
    return out_dict


def unzip_time(dict_t):
    """
    zip data_dict style dict into a time:value style dict.
    """
    times = np.array(sorted(list(dict_t.keys())))
    data = np.array([dict_t[t] for t in times])

    return times, data


def align_times(list_of_pairs):
    """
    try to find the range of the set of data by interpolating over the shared range.
    """
    min_times = []
    max_times = []
    for times, dats in list_of_pairs:
        min_times.append(np.min(times))
        max_times.append(np.max(times))
    ranges = [np.max(min_times), np.min(max_times)] # overlap
    output = []
    out_times = np.arange(ranges[0], ranges[1], 0.1 )
    for times, dats in list_of_pairs:
        interp = scipy.interpolate.interp1d(times, dats)
        out_data = interp(out_times)
        output.append([out_times, out_data])
    return output





def make_cumulative_timeseries(cfg, data_dict,
      thresholds_dict,
      ssp='historical-ssp126',
      dataset = 'CMIP6',
      ensemble = 'ensemble_mean',
      plot_type = 'pc', # 'pc', 'area_over_zero are the two in the megaplot.
      plot_thresholds = [1.5, 2.,3.,4.,5.],
      do_leg= True,
      fig = None, gs = None, ax= None,
):
    """
    Make a plot showing the time series of carbon allocation.
    """
    if None in [fig, ax]:
        fig = plt.figure()
        fig.set_size_inches(12, 6)
        gs = gridspec.GridSpec(1, 1,figure=fig )# width_ratios=[1,1], wspace=0.5, hspace=0.5)
        ax =  fig.add_subplot(gs[0, 0])
        save = True
    else:
        save = False
    plt.sca(ax)

   # load data.
    cube_keys = ['fgco2gt_cumul', 'nbpgt_cumul']
    colours = {#'cumul_emissions': 'silver',
        'atmos_carbon': 'silver',
        'fgco2gt_cumul':'dodgerblue',
        'nbpgt_cumul':'orange',
        'tls':'mediumseagreen',
        'luegt': 'purple'}

    #colours = {'cumul_emissions': 'grey', 'fgco2gt_cumul':'blue', 'nbpgt_cumul':'orange', 'tls':'green', 'luegt':'red'}
    if ensemble == 'ensemble_mean':
        ensembles = ['ensemble_mean',]
    else: assert 0
    if isinstance(dataset, str):
         datasets = [dataset, ]
    else: assert 0

    if ssp[:3] == 'ssp':
         exps = ['historical', ssp, '-'.join(['historical', ssp])]
    elif ssp == 'historical': exps = ['historical', ]

    data = {k:{} for k in colours.keys()}
    for ssp_it, ensemble, key, dset in product(exps, ensembles, colours.keys(),datasets):
        print('load data', (dset, key, ssp_it, ensemble))
        tmp_data = data_dict.get((dset, key, ssp_it, ensemble), False)
        if not tmp_data:
            print('Did not find:', (dset, key, ssp_it, ensemble))
            continue
        print('load data: found:', (dset, key, ssp_it, ensemble))
        if key in cube_keys:
            tmp_data = {'time': diagtools.cube_time_to_float(tmp_data),
                             key: tmp_data.data}

        tmp_data = zip_time(tmp_data, key)
        data[key] = data[key] | tmp_data # combine two dicts (python 3.9 and later)
        print('load data', (key, ssp_it, ensemble),':', len(data[key].keys()))

    print(dataset, ssp, 'loaded data:',data.keys())
    # add atmos_stock:
    tmp_times, tmp_dat = unzip_time(data['atmos_carbon'])
    data['atmos_carbon'] = zip_time({'time':tmp_times, 'atmos_carbon':tmp_dat}, 'atmos_carbon')

    # found=0
    # for key in ['nbpgt_cumul', 'fgco2gt_cumul']:
    #     print('adding atmospheric stock', key)
    #     key_times, key_dat = unzip_time(data[key])
    #     print(dataset, ssp, tmp_dat, key_dat, key_times, tmp_times)
    #     if len(tmp_times) != len(key_times):
    #         print('error:', len(tmp_times), '!=', len(key_times))
    #         print('carbon_times:', tmp_times[0], tmp_times[-1])
    #         print(key, 'times', key_times[0], key_times[-1])
    #         assert 0
    #         continue
    #     tmp_dat = tmp_dat - key_dat
    #     found+=1
    # if found==0:return fig, ax

    #data['atmos_carbon'] = zip_time({'time':tmp_times, 'atmos_carbon':tmp_dat}, 'atmos_carbon')
    #colours['atmos_carbon'] = 'purple'

    thresholds = {}
    for ssp_it, ensemble, dset in product(exps, ensembles, datasets):
        dicts = thresholds_dict.get((dset, 'tas', ssp_it, ensemble), False)
        if not dicts:continue
        thresholds= thresholds|dicts

    # plot simple time series:
#    if plot_type == 'simple_ts':
#        for key, dat in data.items():
#            times, dat_list = unzip_time(dat)
#            #if key == 'tls':
#            #    times, dat = times[:-1], dat[:-1]
#            plt.plot(times,
#                dat_list,
#                lw=2,
#                color=colours[key],
#                label = key)
#        plt.axhline(y=0., c='k', ls='--')
#        if do_leg: plt.legend()
#
#    # plot simple time series:
#    if plot_type == 'sink_source':
#        for key, dat in data.items():
#            times, dat_list = unzip_time(dat)
#            if key in ['fgco2gt_cumul', 'nbpgt_cumul', 'tls', 'luegt']:
#                dat_list = -1*dat_list
#            plt.plot(times,
#                dat_list,
#                lw=2,
#                color=colours[key],
#                label = key)
#        plt.axhline(y=0., c='k', ls='--')
#        if do_leg: plt.legend()

    lat, lad = unzip_time(data['tls'])
    ont, ond = unzip_time(data['fgco2gt_cumul'])
    #        lut, lud = unzip_time(data['luegt'])
    nbt, nbd = unzip_time(data['nbpgt_cumul'])
    att, atd = unzip_time(data['atmos_carbon'])
    [[lat, lad], [ont, ond ], [nbt, nbd], [att, atd]] = align_times([[lat, lad], [ont, ond ], [nbt, nbd], [att, atd]])
    emt = lat[:]
    emd = lad + ond + atd
    aligned_interpolated_data_shelve = diagtools.folder([cfg['work_dir'], 'aligned_interpolated_data'])+'aligned_interpolated_data.shelve'
    sh = shelve.open(aligned_interpolated_data_shelve)
    sh['lat'] = lat
    sh['lad'] = lad
    sh['ont'] = ont
    sh['ond'] = ond
    sh['nbt'] = nbt
    sh['nbd'] = nbd
    sh['att'] = att
    sh['atd'] = atd
    sh['emt'] = emt
    sh['emd'] = emd
    sh.close()


    # plot simple time series:
    if plot_type in ['area', 'area_over_zero']:
        #colours = {'cumul_emissions': 'grey', 'fgco2gt_cumul':'blue', 'nbpgt_cumul':'orange', 'tls':'green'}
        print(data.keys())
        # emt, emd = unzip_time(data['cumul_emissions'])
        lat, lad = unzip_time(data['tls'])
        ont, ond = unzip_time(data['fgco2gt_cumul'])
#        lut, lud = unzip_time(data['luegt'])
        nbt, nbd = unzip_time(data['nbpgt_cumul'])
        att, atd = unzip_time(data['atmos_carbon'])

        [[lat, lad], [ont, ond ], [nbt, nbd], [att, atd]] = align_times([[lat, lad], [ont, ond ], [nbt, nbd], [att, atd]])


        if plot_type in ['area',] :
            ond = -1.* ond
            lad = -1.*(lad + ond)
        if plot_type in ['area_over_zero',] :
            atd_t = atd + lad + ond # air land sea.
            lad_t = lad + ond
            ond_t = ond
#            lad = lad+ond
 #           atd = atd + lad# air land sea.

            # print('emissions times:', emt[:5], emt[-5:])
#            print('LUE times:', lut[:5], lut[-5:])
            #print('emissions data:', emd[:5], emd[-5:])
#            print('LUE data :', lud[:5], lud[-5:])
            # print('sizes:', len(emt), len(emd), len(lut), len(lud))
#            if np.max(lut) > np.max(att):
#                lut = lut[:-1] # Remove 2015.5 (not in histor period)
#                lud = lud[:-1]

#            plt.plot(emt, emd+lud, 'k-', lw=1.3,label='Emissions+LUE')
#            plt.plot(lat, lad, 'k--', lw=1.3, label='LUE')
        # water:
        # plt.plot(
        #     ont,
        #     ond,
        #     lw=2,
        #     color=colours['fgco2gt_cumul'],
        #     label = )
        ax.fill_between(# zero and blue line
            ont,
            np.zeros_like(ond_t),
            ond_t, # zero and blue line
            lw=0,
            label='Ocean',
            color=colours['fgco2gt_cumul'])

        # land:
        # plt.plot(
        #     lat,
        #     lad,
        #     lw=2,
        #     color=colours['tls'],
        #     label = 'True Land Sink')
        ax.fill_between( # blue line and land line
            lat,
            ond_t, # ocean line
            lad_t, # land + ocean
            #lw=2,
            color=colours['tls'],
            label = 'Land')

        # atmos
        # plt.plot(
        #     att,
        #     atd,
        #     lw=2,
        #     color=colours['cumul_emissions'],
        #     label = 'Atmosphere')
        ax.fill_between( # atmos anthro stock
            lat,
            lad_t, # land + ocean
            atd_t, # land ocean and air
            #lw=2,
            color=colours['atmos_carbon'],
            label = 'Atmosphere')
        #ax.set_xlim([2010., 2100])
        ax.set_ylabel('Cumulative carbon, Pg')
        if plot_type in ['area',]:
            plt.axhline(y=0., c='k', ls='--')

    # plot simple time series:
    if plot_type in ['pc', 'triple'] :
        colours = {'atmos_carbon': 'silver', 'fgco2gt_cumul':'dodgerblue', 'nbpgt_cumul':'orange', 'tls':'mediumseagreen'}
        # emt, emd = unzip_time(data['cumul_emissions'])
        lat, lad = unzip_time(data['tls'])
        ont, ond = unzip_time(data['fgco2gt_cumul'])
        #lut, lud = unzip_time(data['luegt'])
        nbt, nbd = unzip_time(data['nbpgt_cumul'])
        att, atd = unzip_time(data['atmos_carbon'])
        [[lat, lad], [ont, ond ], [nbt, nbd], [att, atd]] = align_times([[lat, lad], [ont, ond ], [nbt, nbd], [att, atd]])

        total = ond + lad + atd
        water_land_line = (ond/total)*100.
        land_air_line = 100. *(ond + lad)/total

        water_land_line = np.ma.masked_invalid(water_land_line)
        land_air_line = np.ma.masked_invalid(land_air_line)
        # water:
        plt.plot(
            [], [],
            lw=6,
            color=colours['fgco2gt_cumul'],
            label = 'Air sea Flux of CO2')
        ax.fill_between(# zero and blue line
            att,
            0.* water_land_line, # zero and blue line
            water_land_line,
            color=colours['fgco2gt_cumul'])

        # land:
        plt.plot([],[],
            lw=6,
            color=colours['tls'],
            label = 'True Land Sink')
        ax.fill_between( # blue line and land line
            att,
            water_land_line,
            land_air_line,
            lw=2,
            color=colours['tls'])

        # atmos
        plt.plot([], [],
            lw=6,
            color=colours['atmos_carbon'],
            label = 'Atmosphere')
        ax.fill_between( # atmos anthro stock
            att,
            land_air_line,
            land_air_line*0. +100.,
            lw=2,
            color=colours['atmos_carbon'])

        #plt.axhline(y=0., c='k', ls='--')
        ax.set_ylabel('Percentage')
        #x.set_xlim([2010., 2100])
        ax.set_ylim([0., 100.])

    if ssp in ['historical', ]:
        if plot_type in ['pc', ]:
            plt.plot([1959.,1980., 2000., 2012.], [56., 56., 56., 56.,], c='k', ls=':', label = 'Raupach 2014' )
            plt.plot([1959.,1980., 2000., 2012.], [25., 25., 25., 25.,], c='navy', ls='-.', label = 'Watson 2020')

        else:
            plt.plot([], [], c='k', ls=':', label = 'Raupach 2014' )
            plt.plot([], [], c='navy', ls='-.', label = 'Watson 2020'  )

    if do_leg:plt.legend(fontsize='small')

    threshold_stlye='new'
    if threshold_stlye=='old':
        print(thresholds)
        for thres, dt in thresholds.items():
            if dt is None: continue
            print('adding threshold lineL', thres, dt)
            if thres not in plot_thresholds: continue

            plt.axvline(x=float(dt.year)+0.5, c='k', ls=':' )
            x = float(dt.year) +0.02 *(np.max(ax.get_xlim()) - np.min(ax.get_xlim()))
            y = np.max(ax.get_ylim())- 0.11 * (np.max(ax.get_ylim()) - np.min(ax.get_ylim()))
            plt.text(x, y, 'GWT: ' +str(thres), ha='right', va='top', rotation=90) #fontsize=8, fontweight='bold',rotation=90)


   # threshold_colours = {2.0: 'darkorange', 3.0:'red', 4.0:'purple',}
    threshold_colours = {1.5: 'black', 2.0: 'darkblue', 3.0:'darkred', 4.0:'purple',}

    if threshold_stlye=='new':
        print(thresholds)
        for thres, dt in thresholds.items():
            if dt is None: continue
            if thres not in plot_thresholds: continue
            print('adding new threshold lineL', thres, dt)
            index = np.argmin(np.abs(dt.year- att))
            #print(dt.year, emt, index, emt[index], atd[index], emd[index],lud[index])
            x=float(dt.year)+0.5
            if plot_type == 'pc':
                ymax=1.
            elif plot_type == 'area_over_zero':
                 ymax = atd_t[index]/2000.
            plt.axvline(
                x=x,
                ymin=0.,
                ymax = ymax,
                c=threshold_colours[thres],
                alpha=0.7,
                lw=1.7,
                ls='-',)

            if plot_type == 'pc':
                if thres == 1.5:
                    txt = ''.join('1.5', r'$\degree$', 'C - ', str(int(dt.year)))
                else: txt = str(int(thres)) + r'$\degree$'+'C - '+str(int(dt.year))
                plt.text(x+2.4, 5, txt,
                    c = threshold_colours[thres],
                    ha='left', # Left puts txt on right of line. Right puts txt on left of line.
                    va='bottom', # top puts txt below x axes.
                    rotation=90)
                   # fontsize=8, fontweight='bold',rotation=90)

    plt.title(ssp_title_dict.get(ssp, None))

    # save plot.
    if save:
        image_extention = diagtools.get_image_format(cfg)
        path = diagtools.folder([cfg['plot_dir'], 'allocation_timeseries'])
        path += '_'.join(['allocation_timeseries', plot_type, ssp]) + image_extention
        print('saving figure:', path)
        plt.savefig(path)
        plt.close()
    else: return fig, ax



def make_cumulative_timeseries_megaplot(cfg, data_dict,
       thresholds_dict,
       ssps= ['historical', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'],
       plot_types = ['pair', 'area_over_zero'],
       ensemble = 'ensemble_mean',
       dataset='CMIP6',):
    """
    Single plot massive thing.
    """
    fig = plt.figure()
    fig.set_size_inches(12, 8)
    gs = gridspec.GridSpec(5, 4, figure=fig,
          height_ratios=[1,1,0.25,1,1], hspace=0.1,
          width_ratios=[1,1,1,0.4], wspace=0.130 )# width_ratios=[1,1], wspace=0.5, hspace=0.5)

    ssp_points = {
        'historical':[0, 0], # row, column
        'ssp119':[0, 1],
        'ssp126':[0, 2],
        'ssp245':[3, 0],
        'ssp370':[3, 1],
        'ssp585':[3, 2],
    }
    ax_leg = fig.add_subplot(gs[:, -1])

    for ssp in ssps:
        ax1 =  fig.add_subplot(gs[ssp_points[ssp][0], ssp_points[ssp][1]])
        ax2 =  fig.add_subplot(gs[ssp_points[ssp][0]+1,ssp_points[ssp][1]])

        fig, ax1 = make_cumulative_timeseries(cfg, data_dict,
            thresholds_dict,
            ssp=ssp,
            ensemble = ensemble,
            dataset = dataset,
            plot_type = 'area_over_zero',
            plot_thresholds = [2., 3., 4.],
            fig = fig,
            ax= ax1,
            do_leg=False,
        )
        ax1.set_ylim([0., 2000.])
        ax1.set_xticks([])
        ax1.set_title('')
        ax1.text(0.01, 1.01, sspify(ssp) , horizontalalignment='left',
            verticalalignment='bottom', transform=ax1.transAxes)

        fig, ax2 = make_cumulative_timeseries(cfg, data_dict,
            thresholds_dict,
            ssp=ssp,
            ensemble = ensemble,
            dataset = dataset,
            plot_type = 'pc',
            plot_thresholds = [2., 3., 4.],
            fig = fig,
            ax= ax2,
            do_leg=False,
        )
        ax2.set_title('')
        #ax1.grid(axis='y')
        #ax2.grid(axis='y')

        if ssp_points[ssp][1]>0:
            ax1.set_ylabel('')
            ax2.set_ylabel('')
            ax1.set_yticks([])
            ax2.set_yticks([])

#       if ssp_points[ssp][0]>0:
#           ax1.set_xlabel('')
#           ax2.set_xlabel('')

        if ssp == 'historical':
            ax1.set_xlim([1860., 2015.])
            ax2.set_xlim([1860., 2015.])
        else:
            ax1.set_xlim([2015., 2100.])
            ax2.set_xlim([2015., 2100.])

    plt.sca(ax_leg)
    colours = { #'cumul_emissions': 'silver',
        'atmos_carbon': 'silver',
        'fgco2gt_cumul':'dodgerblue',
        'nbpgt_cumul':'orange',
        'tls':'mediumseagreen',
        'luegt': 'purple'}

    plt.plot([],[], c='silver', lw=8, ls='-', label = 'Atmosphere')
    plt.plot([],[], c='mediumseagreen', lw=8, ls='-', label = 'Land')
    plt.plot([],[], c='dodgerblue', lw=8, ls='-', label = 'Ocean')

    #threshold_colours = {2.0: 'darkorange', 3.0:'red', 4.0:'purple',}
    threshold_colours = {2.0: 'darkblue', 3.0:'darkred', 4.0:'purple',}
    #plt.plot([],[], c=threshold_colours[2.0], ls='-', lw=1.7, label='2'+r'$\degree$'+ ' GWT', alpha=0.7,)
    #plt.plot([],[], c=threshold_colours[3.0], ls='-', lw=1.7, label='3'+r'$\degree$'+ ' GWT', alpha=0.7,)
    #plt.plot([],[], c=threshold_colours[4.0], ls='-', lw=1.7, label='4'+r'$\degree$'+ ' GWT', alpha=0.7,)

    #plt.plot([],[], 'k-', lw=1.3,label='Emissions+LUE')
    #plt.plot([],[], 'k--', lw=1.3, label='LUE')
    plt.plot([], [], c='k', ls=':', label = 'Raupach (2014)' )
    plt.plot([], [], c='navy', ls='-.', label = 'Watson (2020)'  )

    legd = ax_leg.legend( #keys, labels,
        bbox_to_anchor=(1.7, 0.5),
        numpoints=1, labelspacing=1.2,
        loc='center right', ) #tsize=16)

    legd.draw_frame(False)
    legd.get_frame().set_alpha(0.)
    ax_leg.get_xaxis().set_visible(False)
    ax_leg.get_yaxis().set_visible(False)
    plt.axis('off')

    if dataset == 'CMIP6':
        plt.suptitle('Anthropogenic Carbon Allocation Timeseries')
    else:
        plt.suptitle(dataset+' Anthropogenic Carbon Allocation Timeseries')

    image_extention = diagtools.get_image_format(cfg)
    path = diagtools.folder([cfg['plot_dir'], 'allocation_timeseries_megaplot'])
    path += '_'.join(['allocation_timeseries', 'megaplot', dataset]) + image_extention
    print('saving figure:', path)
    plt.savefig(path)
    plt.close()




def make_cumulative_timeseries_pair(cfg, data_dict,
      thresholds_dict,
      ssp='ssp126',
      ensemble = 'ensemble_mean',
    ):

    fig = plt.figure()
    fig.set_size_inches(6, 6)
    gs = gridspec.GridSpec(2, 1,figure=fig )# width_ratios=[1,1], wspace=0.5, hspace=0.5)
    ax1 =  fig.add_subplot(gs[0, 0])
    ax2 =  fig.add_subplot(gs[1, 0])

    fig, ax1 = make_cumulative_timeseries(cfg, data_dict,
        thresholds_dict,
        ssp=ssp,
        ensemble = ensemble,
        dataset = 'CMIP6',
        plot_type = 'area_over_zero',
        fig = fig, ax= ax1,
    )
    fig, ax2 = make_cumulative_timeseries(cfg, data_dict,
        thresholds_dict,
        ssp=ssp,
        ensemble = ensemble,
        dataset = 'CMIP6',
        plot_type = 'pc',
        fig = fig, ax= ax2,
        do_leg=False,
    )
    if ssp == 'historical':
        ax1.set_xlim([1850., 2015.])
        ax2.set_xlim([1850., 2015.])
#        ax1.tick_params(labeltop=False, labelright=True)
#        ax2.tick_params(labeltop=False, labelright=True)

    else:
        ax1.set_xlim([2015., 2100.])
        ax2.set_xlim([2015., 2100.])
#        ax1.tick_params(labeltop=False, labelright=True)
#        ax2.tick_params(labeltop=False, labelright=True)
    ax1.grid(axis='y')
    ax2.grid(axis='y')

    image_extention = diagtools.get_image_format(cfg)
    path = diagtools.folder([cfg['plot_dir'], 'allocation_timeseries'])
    path += '_'.join(['allocation_timeseries', 'pair', ssp]) + image_extention
    print('saving figure:', path)
    plt.savefig(path)
    plt.close()


def do_emission_scatterplot(cfg, data_dict, thresholds_dict, x_axis = 'emissions', fig = None, ax = None):

    save_single_plot= False
    if fig is None or ax is None:
        fig = plt.figure()
        fig.set_size_inches(7 , 6)
        ax = fig.add_subplot(111)
        save_single_plot = True
    else:
        plt.sca(ax)
        #gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.5)

    Model_colours = {}
    marker_styles = {1.5: None, 2.: 's', 3.: 'o', 4:'^', 5.:None}

    aligned_interpolated_data_shelve = diagtools.folder([cfg['work_dir'], 'aligned_interpolated_data'])+'aligned_interpolated_data.shelve'
    if not len(glob.glob(aligned_interpolated_data_shelve+'*')):
        print('ERROR:, missing data, run make_cumulative_timeseries first!')
        assert 0
    sh = shelve.open(aligned_interpolated_data_shelve)
    lat = sh['lat']
    lad = sh['lad']
    ont = sh['ont']
    ond = sh['ond']
    nbt = sh['nbt']
    nbd = sh['nbd']
    att = sh['att']
    atd = sh['atd']
    emt = sh['emt']
    emd = sh['emd']
    sh.close()

    if x_axis == 'emissions':
        y_time = emt
        y_data = emd
    .
    if x_axis == 'land':
        y_time = lat
        y_data = lad

    if x_axis == 'ocean':
        y_time = ont
        y_data = ond

    if x_axis == 'atmosphere':
        y_time = att
        y_data = atd

    for (dataset, short_name, exp, ensemble), thresholds in thresholds_dict.items():
        if short_name != 'tas': continue
        if ensemble != 'ensemble_mean': continue
        colour = Model_colours[dataset]

        for threshold, time in thresholds.items():
            if time is None: continue # no GWT.
            ms = marker_styles.get(threshold, None)
            if not ms: continue
            x = time.year + 0.5
            index - np.argmin(np.abs(y_time - x))
            y = y_data[index]
            #label =
            plt.scatter(x,y, c = colour, maker = ms)


    if not save_single_plot: return fig, ax

    image_extention = diagtools.get_image_format(cfg)
    path = diagtools.folder([cfg['plot_dir'], 'emissions_scatter',])
    path += '_'.join(['emssions_scatter', x_axis]) + image_extention
    plt.savefig(path)
    plt.close()



def make_cumulative_vs_threshold(cfg, data_dict,
    thresholds_dict,
    land_carbon = 'nbpgt',
    LHS_panes = [{'x':'atmos_carbon', 'y':'tas_norm'}, ],
    thresholds = ['4.0', '3.0', '2.0'],
    plot_dataset='CMIP6',
):
    """
    Make a specific kind of figure.
    land_carbon can be nbpgt or tls
    """
    #standalone
    #for threshold in thresholds:
    #    make_bar_chart(cfg, data_dict, thresholds_dict, threshold = threshold, land_carbon = land_carbon, fig=None, ax=None)

    fig = plt.figure()
    if len(LHS_panes) ==0:
        fig.set_size_inches(7 , 6)
        gs = gridspec.GridSpec(3, 1, figure=fig, hspace=0.5)
    else:
        fig.set_size_inches(14, 6)
        gs = gridspec.GridSpec(3, 2,figure=fig, width_ratios=[1,1], wspace=0.5, hspace=0.5)

    # Get a big line graph on LHS.
    if len(LHS_panes) ==1:
        ax_ts =  fig.add_subplot(gs[:, 0])
        fig, ax = make_ts_figure(cfg, data_dict, thresholds_dict,
            x=LHS_panes[0]['x'],
            y=LHS_panes[0]['y'],
            markers='thresholds',
            draw_line=True,
            do_moving_average=False,
            ensemble_mean = True,
            fig=fig,
            ax = ax_ts,)

    if len(LHS_panes) ==3:
        for i, LHS_pane in enumerate(LHS_panes):
            ax_ts =  fig.add_subplot(gs[i, 0])
            fig, ax = make_ts_figure(cfg, data_dict, thresholds_dict,
                x=LHS_pane['x'],
                y=LHS_pane['y'],
                markers='thresholds',
                draw_line=True,
                do_moving_average=False,
                ensemble_mean = True,
                fig=fig,
                ax = ax_ts,)
            if LHS_pane['x'] == 'time':
                ax_ts.set_xlim([2000., 2100.])


    if len(LHS_panes) ==0:
        ax_4 =  fig.add_subplot(gs[0, 0])
        ax_3 =  fig.add_subplot(gs[1, 0])
        ax_2 =  fig.add_subplot(gs[2, 0])
    else:
        ax_4 =  fig.add_subplot(gs[0, 1])
        ax_3 =  fig.add_subplot(gs[1, 1])
        ax_2 =  fig.add_subplot(gs[2, 1])

    axes = [ax_4, ax_3, ax_2]
    # make_bar_chart(cfg, data_dict, thresholds_dict, threshold = '1.5', fig=None, ax=None)
    legends= {ax_4:False, ax_3:False, ax_2:True}
    for ax, threshold in zip(axes, thresholds):
        make_bar_chart(cfg, data_dict, thresholds_dict,
                       threshold = threshold, land_carbon = land_carbon,fig=fig, ax=ax, do_legend=legends[ax],
                       plot_dataset = plot_dataset)

    ranges = []
    height_ratios = []

    for ax in [ax_4, ax_3, ax_2]:
        plt.sca(ax)
        ranges.append(ax.get_xlim())
        height_ratios.append(np.max(ax.get_ylim()) - np.min(ax.get_ylim()))

    for ax in [ax_4, ax_3, ax_2]:
        plt.sca(ax)
        ax.set_xlim([np.min(ranges), np.max(ranges)])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if ax in [ax_3, ax_2]:
            ax.set_xlabel('')
            ax.spines['left'].set_visible(False)
        if ax == ax_2:
            ax.legend(loc='upper right')

    if len(LHS_panes) ==0:
        # Set y axis size to make all panes the same aspect ratio (maybe)
        print('New height ratios:', height_ratios)
        gs.set_height_ratios(height_ratios)
        fig.subplots_adjust(wspace=0.05)
    # plt.suptitle('Carbon Allocation')

    if plot_dataset != 'all_models':
        plt.suptitle(' '.join([plot_dataset, 'Carbon Allocation', ]))
    image_extention = diagtools.get_image_format(cfg)
    path = diagtools.folder([cfg['plot_dir'], 'emissions_figures_new', plot_dataset])
    path += '_'.join(['emssions_figure', land_carbon, str(len(LHS_panes)), '_'.join(thresholds), plot_dataset]) + image_extention
    plt.savefig(path)
    plt.close()
    #assert 0
    # 4 degree threshold data.
    # for each scenario at 4 degrees, we want:
    # total emissions at time_4
    # nbp and


def timeseries_megapane(cfg, data_dict, thresholds_dict, key,
    plot_styles = ['CMIP6_range', 'CMIP6_mean'],
    experiments = ['historical', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'],
    #experiments = ['historical', 'historical-ssp119', 'historical-ssp126', 'historical-ssp245', 'historical-ssp370', 'ssp585'],

    fig = None,
    ax = None,
    ):
    """
    Single time series pane for the megaplot.
    """
    make_ts_figure(cfg, data_dict, thresholds_dict, x='time', y=key,
        markers='thresholds',
        draw_line=True,
        do_moving_average=False,
        plot_dataset='all_models',
        ensemble_mean = True,
        fig=fig,
        ax=ax,
        do_legend = False,
        plot_thresholds = [2., 3., 4.,],
        plot_styles=plot_styles,
        skip_historical_ssp = True,
        experiments = experiments,
        #short_time_range = False,
        )
    ax.set_xlim([1850., 2100.])
    return fig, ax


def timeseries_megaplot(cfg, data_dict, thresholds_dict,
        panes = ['tas_norm', 'atmos_carbon', 'fgco2gt_cumul', 'nbpgt_cumul', ],
        plot_styles = ['CMIP6_range', 'CMIP6_mean'],
        experiments = ['historical', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'],
        ):
    """
    4 pane time series plot which shows the:
    temperature                      -  Anthropogenic emissions
    Global total air sea flux of co2 -  Net biome production?

    or something like that.
    """
    fig = plt.figure()
    fig.set_size_inches(12, 8)

    if isinstance(plot_styles, str):
        plot_styles = [plot_styles, ]

    axes = {}
    if len(panes) == 1:
        gs = gridspec.GridSpec(1, 2,figure=fig, width_ratios=[1, 0.1], wspace=0.352, hspace=0.352)
        axes[panes[0]] = fig.add_subplot(gs[0, 0]) # row, column

    if len(panes) == 4:
        gs = gridspec.GridSpec(2, 3,figure=fig, width_ratios=[1,1, 0.1], wspace=0.352, hspace=0.352)
        axes[panes[0]] = fig.add_subplot(gs[0, 0]) # row, column
        axes[panes[1]] = fig.add_subplot(gs[0, 1])
        axes[panes[2]] = fig.add_subplot(gs[1, 0])
        axes[panes[3]] = fig.add_subplot(gs[1, 1])

    if len(panes) in [5, 6]:
        gs = gridspec.GridSpec(2, 4,figure=fig, width_ratios=[1, 1, 1, 0.21], wspace=0.2, hspace=0.6)
        axes[panes[0]] = fig.add_subplot(gs[0, 0]) # row, column
        axes[panes[1]] = fig.add_subplot(gs[0, 1])
        axes[panes[2]] = fig.add_subplot(gs[0, 2])
        axes[panes[3]] = fig.add_subplot(gs[1, 0])
        axes[panes[4]] = fig.add_subplot(gs[1, 1])
        if len(panes) in [6, ]:
            axes[panes[5]] = fig.add_subplot(gs[1, 2])

    # add legend column ion RHS:
    axes['legend'] = fig.add_subplot(gs[:, -1])

    for key, ax in axes.items():
        if key == 'legend': continue
        fig, ax = timeseries_megapane(cfg, data_dict, thresholds_dict, key,
            plot_styles=plot_styles,
            experiments=experiments,
            fig= fig, ax = ax)

    # l;egend pane:
    plt.sca(axes['legend'])
    keys, labels = [], []
#    exp_colours = {'historical':'black',
#                   'ssp119':'green',
#                   'ssp126':'dodgerblue',
#                   'ssp245':'blue',
#                   'ssp370':'purple',
#                   'ssp434':'magenta',
#                   'ssp585': 'red',
#                   'ssp534-over':'orange',}
                   # 'historical-ssp119':'green',
                   # 'historical-ssp126':'dodgerblue',
                   # 'historical-ssp245':'blue',
                   # 'historical-ssp370':'purple',
                   # 'historical-ssp434':'magenta',
                   # 'historical-ssp585': 'red',
                   # 'historical-ssp585-ssp534-over':'orange'}

    if 'CMIP6_mean' in plot_styles:
        plt.plot([],[], ls='-', c='k', lw=1.9, label = 'CMIP6 mean')

    if 'CMIP6_range' in plot_styles:
        plt.plot([],[], ls='-', c='k', alpha =0.5, lw=8., label = 'CMIP6 range')

    if 'all_models_means' in plot_styles:
        plt.plot([],[], ls='-', c='k', lw=1.0, label = 'Model mean')

    if 'all_models_range' in plot_styles:
        plt.plot([],[], marker='s', ms=6, markeredgecolor='black', color=(1,1,1,0.5), label = 'Model range')

    if 'all_ensembles' in plot_styles:
        plt.plot([],[], ls='-', c='k', lw=0.5, label = 'Ensemble')

    for exp in experiments:
        plt.plot([],[], ls='-', c=exp_colours[exp], lw=4., label = sspify(exp))
        #keys.append(plt.plot([],[], ls='-', c=exp_colours[exp], lw=4.), label = )
        #labels.append(sspify(exp))

    marker_styles = {1.5: '*', 2.:'o', 3.:'D', 4.:'s', 5.:'X'}
    for gwt in [2., 3., 4.]:
        lab = ''.join([str(int(gwt)), r'$\degree$', 'C'])
        plt.scatter([], [], marker=marker_styles[gwt], c='k', label=lab)
        #keys.append(plt.scatter([], [], marker=marker_styles[gwt], c='k',))
        #labels.append(''.join([str(int(gwt)), r'$\degree$', 'C']))


    legd = axes['legend'].legend(#keys, labels,
        bbox_to_anchor=(2.5, 0.5),
        #numpoints=1, labelspacing=1.2,
        loc='center right', ) #tsize=16)

    legd.draw_frame(False)
    legd.get_frame().set_alpha(0.)
    axes['legend'].get_xaxis().set_visible(False)
    axes['legend'].get_yaxis().set_visible(False)
    plt.axis('off')

    image_extention = diagtools.get_image_format(cfg)
    if len(panes)==1:
        path = diagtools.folder([cfg['plot_dir'], 'timeseries_megaplots', panes[0]])
    else:
        path = diagtools.folder([cfg['plot_dir'], 'timeseries_megaplots'])

    path += '_'.join(['timeseries_megaplot', '_'.join(panes), '_'.join(plot_styles), '_'.join(experiments) ])
    path = ''.join([path, image_extention])
    print('Save image:', path)
    plt.savefig(path)
    plt.close()


def fix_indices(data_dict):
    new_dict = {}
    for (dset, key, ssp_it, ensemble), d in data_dict.items():
        ssp_it = ssp_it.replace('historical-', '')
    data_dict.update(new_dict)
    return data_dict


def main(cfg):
    """
    Load the config file and some metadata, then make plots.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

    """

    #jobtype = 'debug'
    jobtype = 'cumulative_plot'

    if jobtype == 'cumulative_plot':
        short_names = ['tas', 'tas_norm',
                       'co2', #'emissions', #'cumul_emissions',
                       'nbp', 'nbpgt', 'nbpgt_cumul',
                       #'gpp', 'gppgt',
                       #'intpp',  'intppgt',
                       'fgco2','fgco2gt', 'fgco2gt_cumul',
                       'luegt', #  land-use emissions gt
                       'tls', #true land sink = nbp + land-use emissions
                       'atmos_carbon', # remant antrho carbon in atmosphere
                       ]
        short_names_x = ['time', ] #'atmos_carbon', 'tas_norm', 'cumul_emissions', 'fgco2gt_cumul','nbpgt_cumul', ]

        short_names_y = short_names.copy()

#    if jobtype == 'bulk':
#        short_names = ['tas', 'tas_norm', 'co2', 'emissions', #'cumul_emissions'
#                       'nbp', 'nbpgt', 'gpp', 'gppgt',
#                       'intpp', 'fgco2', 'intppgt','fgco2gt',
#                       'fgco2gt_cumul',
#                       'nbpgt_cumul'
#                       'luegt', #  land-use emissions gt
#                       'tls', #true land sink = nbp + land-use emissions
#                       'atmos_carbon', # remant antrho carbon in atmosphere
#                       ]
#        short_names_x = ['time', 'co2', 'tas_norm', 'fgco2gt_cumul','nbpgt_cumul', ]
#        short_names_y = short_names.copy()

    if jobtype == 'debug':
        short_names = [
                       #'emissions', 'cumul_emissions',
                       'co2', 'atmos_carbon',
                       'tas',
#                      'nbp', #'nbpgt', 'nbpgt_cumul',
#                      #'gpp', 'gppgt',
#                      'fgco2',#'fgco2gt', 'fgco2gt_cumul',
#                      'bgc', #'bgcgt',
#                      #'luegt', #  land-use emissions gt
#                       'tls', #true land sink = nbp + land-use emissions
#                      'atmos_carbon', # remant antrho carbon in atmosphere
                       ]
        short_names_x = ['time', ]#ssions',]#'time', 'emissions', 'cumul_emissions',]#'cumul_emissions', 'gppgt'] #'time', ]#'co2', 'emissions', 'tas_norm', 'fgco2gt', 'nbpgt']
        short_names_y = short_names.copy()


#    if jobtype == 'full':
#        short_names = ['tas', 'tas_norm', 'co2',
#                       'npp', 'nppgt', 'rhgt', 'exchange',
#                       'nppgt_norm','rhgt_norm','exchange_norm','fgco2gt_norm', 'intppgt_norm',
#                       'intpp', 'fgco2', 'epc100', 'intdic', 'intpoc', 'fric', 'froc',
#                       'intppgt','fgco2gt', 'epc100gt', 'intdicgt', 'intpocgt', 'fricgt', 'frocgt',
#                       ]
#        short_names_x = ['time', 'co2', 'tas', 'tas_norm', 'fgco2gt',
#                        'intpp', 'epc100', 'intdic', 'intpoc', 'fric', 'froc', 'nppgt', 'fgco2gt', 'rhgt', 'exchange']
#        short_names_y = ['nppgt', 'nppgt_norm','rhgt_norm','exchange_norm','fgco2gt_norm', 'co2',
#                         'intpp', 'fgco2', 'epc100', 'intdic', 'intpoc', 'fric', 'froc', 'fgco2gt', 'intppgt','epc100gt', 'intdicgt', 'intpocgt', 'fricgt', 'frocgt',]
#
    pairs = []

    for do_ma in [True, ]:#False]:
        data_dict = load_timeseries(cfg, short_names)
        data_dict = load_scenario_carbon(cfg, data_dict)
        # data_dict = calc_emissions(cfg, data_dict)

        #data_dict = calc_model_mean(cfg, short_names, data_dict)

        #data_dict = fix_indices(data_dict)

        thresholds_dict = load_thresholds(cfg, data_dict)
        datasets = {}
        for (dataset, short_name, exp, ensemble),cube  in data_dict.items():
            datasets[dataset] = True

        #plot_data_dict(cfg, data_dict)
        #ake_bar_chart(cfg, data_dict, thresholds_dict, threshold = '4.0', land_carbon = 'tls')
        #ake_bar_chart(cfg, data_dict, thresholds_dict, threshold = '3.0', land_carbon = 'tls')
        #ake_bar_chart(cfg, data_dict, thresholds_dict, threshold = '2.0', land_carbon = 'tls')


        do_emission_scatter = True
        if do_emission_scatter:
            do_emission_scatterplot(cfg, data_dict, thresholds_dict)
        #assert 0

        if jobtype in ['cumulative_plot', 'debug']:
            #make_cumulative_vs_threshold(cfg, data_dict, thresholds_dict, land_carbon = 'tls')
            #make_cumulative_vs_threshold(cfg, data_dict, thresholds_dict, land_carbon = 'nbpgt')
            #make_cumulative_timeseries(cfg, data_dict, thresholds_dict, ssp='historical-ssp585',)
            #make_cumulative_timeseries(cfg, data_dict, thresholds_dict, ssp='historical',)

            do_count_and_sensitivity_table = False
            if do_count_and_sensitivity_table:
                make_count_and_sensitivity_table(cfg, data_dict, thresholds_dict)

            do_timeseries_megaplot = False
            if do_timeseries_megaplot:
                # master
                timeseries_megaplot(cfg, data_dict, thresholds_dict,plot_styles=['CMIP6_range', 'CMIP6_mean'],
                        panes = ['tas_norm', 'atmos_carbon', 'fgco2gt_cumul', 'nbpgt_cumul', ],) # defaults

                for plot_styles in [
                       'CMIP6_range', 'all_models_range', 'all_models_means',
                       'all_ensembles',  'CMIP6_mean'
                       ]:
                    for pane in ['tas_norm','atmos_carbon','tls', 'fgco2gt_cumul',]:
                        timeseries_megaplot(cfg, data_dict, thresholds_dict,plot_styles=plot_styles,
                            panes = [pane, ])

#                timeseries_megaplot(cfg, data_dict, thresholds_dict,plot_styles=['CMIP6_range', 'CMIP6_mean', 'all_models_means', 'all_ensembles'],
#                        panes = ['atmos_carbon', ],) # defaults
#                timeseries_megaplot(cfg, data_dict, thresholds_dict,plot_styles=['CMIP6_range', 'CMIP6_mean', 'all_models_means', 'all_ensembles'],
#                        panes = ['luegt', ],) # defaults
#
#                timeseries_megaplot(cfg, data_dict, thresholds_dict,plot_styles=['CMIP6_range', 'CMIP6_mean', 'all_models_means', 'all_ensembles'],
#                        panes = ['atmos_carbon', 'nbpgt_cumul', 'luegt', 'tas', 'nbpgt','tls'],) # defaults
#
#                for plot_styles in [['CMIP6_range', 'CMIP6_mean'],
#                       'CMIP6_range', ['all_models_range', 'all_models_means'],
#                       'all_ensembles',  'CMIP6_mean',
#                       ['all_ensembles','all_models_range',],
#                       ]:
#            #plot_styles: ambition:
#            #    ['ensemble_mean', ] default behaviour before
#            #    CMIP6_mean: Only CMIP6 multi model mean
#            #    CMIP6_range: range bewteen each model mean shown
#            #    CMIP6_full_range: full range between individual model ensemble menmbers
#            #    all_models_means: Each individual models mean is plotted.
#            #    all_models_range: Each individual models range is plotted
#            #    all_ensembles: Every single ensemble member is shown.
#                    continue
#                    # panes = ['tas_norm', 'atmos_carbon', 'fgco2gt_cumul', 'nbpgt_cumul', ],
#                    timeseries_megaplot(cfg, data_dict, thresholds_dict,plot_styles=plot_styles,
#                        panes = ['tas_norm', 'atmos_carbon', 'fgco2gt_cumul', 'nbpgt_cumul', ],) # defaults
#                    timeseries_megaplot(cfg, data_dict, thresholds_dict,plot_styles=plot_styles,
#                        panes = ['tas', 'co2', 'tls', 'fgco2gt', 'luegt', 'nbpgt'])
#                    timeseries_megaplot(cfg, data_dict, thresholds_dict,plot_styles=plot_styles,
#                        panes = ['tas', 'atmos_carbon','tls', 'fgco2', 'lue', 'nbp'])

            do_cumulative_plot = False
            if do_cumulative_plot:

                plot_styles = ['percentages', 'values']
                ens_styles = ['ensemble_mean', 'all_ens']
                group_bys = ['ecs', ] #'group_by_ssp', 'ecs'] # 'group_by_model'
                for plot_style, ens, group_by in product(plot_styles, ens_styles, group_bys):
                    make_ensemble_barchart(cfg, data_dict, thresholds_dict, plot_style=plot_style, ensemble_key=ens, group_by=group_by)

            do_horizontal_plot = False
            if do_horizontal_plot:
                # Horizontal bar charts with allocartions:
                for plotdataset in sorted(datasets.keys()):
                    make_cumulative_vs_threshold(cfg, data_dict, thresholds_dict, land_carbon = 'tls', LHS_panes = {}, plot_dataset=plotdataset)
                    make_cumulative_vs_threshold(cfg, data_dict, thresholds_dict, land_carbon = 'tls', LHS_panes = {}, thresholds=['2075', '2050', '2025'], plot_dataset=plotdataset)


            do_cumulative_ts_megaplot = False
            # Massive plot that has like 12 panes.
            if do_cumulative_ts_megaplot:
                make_cumulative_timeseries_megaplot(cfg, data_dict,
                                       thresholds_dict,
                                       ssps= ['historical', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'],
                                       plot_types = ['pair', 'area_over_zero'],
                                       ensemble = 'ensemble_mean')
                datasets = {}
                for (dataset, short_name, exp, ensemble),cube  in data_dict.items():
                    datasets[dataset] = True
                for dataset in datasets.keys():
                    make_cumulative_timeseries_megaplot(cfg, data_dict,
                             thresholds_dict,
                             ssps= ['historical', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585'],
                             plot_types = ['pair', 'area_over_zero'],
                             ensemble = 'ensemble_mean',
                             dataset=dataset)

            do_make_cumulative_timeseries_pair = False
#           # This is all done in a single plot with make_cumulative_timeseries_megaplot
            if do_make_cumulative_timeseries_pair:
                plot_types = ['pair', ]
                ssps = ['ssp119',] #'historical', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
                for pt, exp in product(plot_types, ssps):
                      make_cumulative_timeseries_pair(cfg, data_dict,
                          thresholds_dict,
                          ssp=exp,
                          ensemble = 'ensemble_mean')
#                      continue
#                   make_cumulative_timeseries(cfg, data_dict, thresholds_dict, ssp=exp, plot_type = pt)
#               #continue



        return
        for x in short_names_x:
            for y in short_names_y:
                make_ts_figure(cfg, data_dict, thresholds_dict, x=x, y=y,
                           markers='thresholds', do_moving_average=False,
                           plot_dataset='all_models',
                           ensemble_mean=True )

        for x in short_names_x:
            for y in short_names_y:
                for plot_dataset in datasets:
                    #if plot_dataset.find('UKESM')==-1: continue
                    if x == y:
                        continue
                    print('main:', do_ma, x, y)
                    make_ts_figure(cfg, data_dict, thresholds_dict, x=x, y=y,
                               markers='thresholds', do_moving_average=False,
                               plot_dataset=plot_dataset,
                               ensemble_mean=False)

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
