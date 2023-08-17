"""
Mpas's plotting tool. PAS
============================

Basically the same as time series, just with a few additions.

Diagnostic to produce figures of the time development of a field from
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
import glob
import itertools
import sys

import iris
import iris.quickplot as qplt

import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from math import log10, floor
import numpy as np

import cartopy
import cartopy.crs as ccrs
from shelve import open as shopen

#from sigfig import round


from esmvalcore.preprocessor._time import extract_time, annual_statistics, regrid_time

from esmvalcore.preprocessor._regrid import extract_levels, regrid


from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))



# Wrong values:
#central_longitude = -14.25 #W
#central_latitude = -7.56

# right value:
central_longitude = -14.416667 #W #-160.+3.5
central_latitude = -7.933333 # S

#pH log:
pH_log=False


ipcc_colours={
    'historical': 'blue',
    'hist': 'blue',
    'ssp126': 'green',
    'ssp245': 'gold',
    'ssp370': 'orange',
    'ssp585': 'red',}
ipcc_colours_dark={
    'historical': 'darkblue',
    'hist': 'blue', #'darkblue',
    'ssp126':'green', # 'darkgreen',
    'ssp245': '#EEBC1D',
    'ssp370': 'darkorange',
    'ssp585': 'red', } #'darkred',}

# official
#colours2_rgb = [(23, 60, 102), (0, 173, 207), (247, 148, 32), (231, 29, 37), (149, 27, 30)]
#colours2_rgb = [(a/256., b/256. ,c/256.) for (a, b, c) in colours2_rgb]
#f = 0.75
#coloursscatter_rgb = [(f*a, f*b, f*c) for (a, b, c) in colours2_rgb] # darker version

ipcc_colours={
#    'historical': (23, 60, 102), # a dark navy
#    'hist': (23, 60, 102), # a dark navy
    'historical': (35, 87, 145), # new darkish blue
    'hist': (35, 87, 145), # new darkish blue
    'ssp126': (0, 173, 207), # light sky
    'ssp245': (247, 148, 32), # orange
    'ssp370': (231, 29, 37), # red
    'ssp585': (149, 27, 30), } # dark maroon red
ipcc_colours = {k: (v1/256., v2/256.,v3/256.) for k, (v1,v2,v3) in ipcc_colours.items()}
color_factor = 0.75
ipcc_colours_dark = {k: (color_factor*v1, color_factor*v2, color_factor* v3) for k, (v1,v2,v3) in ipcc_colours.items()}



long_name_dict = {
    'thetao': 'Temperature',
    'tos': 'Surface Temperature',
    'tob': 'Seafloor Temperature',
    'sos': 'Surface Salinity',
    'uo': 'Zonal Velocity',
    'vo': 'Meridional Velocity',
    'ph': 'Surface pH',
    'chl': 'Surface chlorophyll',
    'zos': 'Sea Surface Height',
    'no3': 'Dissolved Nitrate',
    'o2': 'Dissolved Oxygen',
    'intpp': 'Integrated Primary production'}


suptitles = {
    'tos': 'Temperature, '+r'$\degree$'+'C',
	'sal': 'Salinity',
	'chl': 'Chlorophyll concentration, mg m'+r'$^{-3}$',
	'ph': 'pH',
	'intpp': 'Integrated Primary Production, mol m'+r'$^{-2}$'+' d'+r'$^{-1}$',
	'po4': 'Phosphate Concentration, mmol m'+r'$^{-3}$',
	'no3': 'Nitrate Concentration, mmol m'+r'$^{-3}$',
	'mld': 'Mixed Layer Depth, m',
	'o2': 'Dissolved Oxygen Concentration at 500m, mmol m'r'$^{-3}$',
}
short_suptitles = {
    'tos': 'Temperature, '+r'$\degree$'+'C',
    'sal': 'Salinity',
    'chl': 'Chl., mg m'+r'$^{-3}$',
    'ph': 'pH',
    'intpp': 'Int PP, mol m'+r'$^{-2}$'+' d'+r'$^{-1}$',
    'po4': 'PO'+r'$_{4}$'+', mmol m'+r'$^{-3}$',
    'no3': 'NO'+r'$_{3}$'+', mmol m'+r'$^{-3}$',
    'mld': 'MLD, m',
    'o2': 'O'+r'$_{2}$'+', mmol m'+r'$^{-3}$',
}
anomaly_titles = {
    'tos': 'Temperature Anomaly '+r'$\degree$'+'C',
    'sal': 'Salinity Anomaly',
    'chl': 'Chl. Anomaly, mg m'+r'$^{-3}$',
    'ph': 'pH Anomaly',
    'intpp': 'IntPP Anomaly, mol m'+r'$^{-2}$'+' d'+r'$^{-1}$',
    'po4': 'PO'+r'$_{4}$'+' Anomaly, mmol m'+r'$^{-3}$',
    'no3': 'NO'+r'$_{3}$'+' Anomaly, mmol m'+r'$^{-3}$',
    'mld': 'MLD Anomaly, m',
    'o2': 'O'+r'$_{2}$'+' Anomaly, mmol m'+r'$^{-3}$',
}

for ndict2 in [suptitles, short_suptitles, anomaly_titles]:
    ndict2['thetao'] = ndict2['tos']
    ndict2['mlotst'] = ndict2['mld']
    ndict2['so'] = ndict2['sal']
    ndict2['sos'] = ndict2['sal']
    ndict2['pH'] = ndict2['ph']

#odels_to_skip = ['GISS-E2-1-G', ] #SST is strangely high and there are very few SSP ensembles.
models_to_skip = {'all': ['GISS-E2-1-G', 'MPI-ESM-1-2-HAM',],
    'chl': ['MPI-ESM1-2-LR', ],
    'o2': ['MPI-ESM-1-2-HAM', ],
    'ph': ['MPI-ESM-1-2-HAM', ],
    'sal': ['MPI-ESM-1-2-HAM', 'FIO-ESM-2-0',],
}

# used to
hard_wired_obs = {
    ('so', 'timeseries', 'min'): {1980:35.8812, 2010:35.8812}, # time range assumed from https://www.ncei.noaa.gov/sites/default/files/2020-04/woa18_vol2.pdf page 1
    ('so', 'timeseries', 'max'): {1980:36.50431, 2010:36.50431},
#    ('so', 'clim', 'min'): {1980:35.8812, 2010:35.8812},
#    ('so', 'clim', 'max'): {1980:36.50431, 2010:36.50431},

    ('ph', 'timeseries', 'min'): {1973:8.04368, 2013:8.04368}, # doi:10.5194/essd-8-325-2016
    ('ph', 'timeseries', 'max'): {1973:8.070894, 2013:8.070894},
    ('ph', 'clim', 'min'): {1973:8.04368, 2013:8.04368},
    ('ph', 'clim', 'max'): {1973:8.070894, 2013:8.070894},

    ('mld', 'timeseries', 'min'): {2005:34.53, 2017:34.53}, # 34.53295 56.90122
    ('mld', 'timeseries', 'max'): {2005:56.90   , 2017:56.90   },
#    ('mld', 'clim', 'min'): {1961:51.72, 2008:51.72},
#    ('mld', 'clim', 'max'): {1961:8.070894, 2008:8.070894},

#   ('o2', 'timeseries', 'min'): {1961:93.3, 2017:93.3}, # 750m
#   ('o2', 'timeseries', 'max'): {1961:123.1, 2017:123.1}, # 750m
#   ('o2', 'timeseries', 'min'): {1961:107.3, 2017:107.3}, # 500m
#   ('o2', 'timeseries', 'max'): {1961:123.5, 2017:123.5}, # 500m

   ('o2', 'timeseries', 'min'): {1961:89.1 , 2017:89.1 }, # 500m
   ('o2', 'timeseries', 'max'): {1961:113.1, 2017:113.1}, # 500m

# [110.69357735770089, 108.27589634486607, 117.27530343191964, 119.49379185267857, 116.9248046875, 107.2939453125, 121.88762555803571, 107.13483537946429, 115.68745640345982, 123.4818115234375, 109.77983747209821, 113.68933977399554]

#    ('o2', 'clim', 'min'): {1961:93.3, 2017:93.3}, # clim exists seprately.
#    ('o2', 'clim', 'max'): {1961:123.1, 2017:123.1},

    ('no3', 'timeseries', 'min') : {1971:0.039, 2008:0.039},
    ('no3', 'timeseries', 'max') : {1971:0.411, 2008:0.411},

    ('po4', 'timeseries', 'min') : {1971:0.126, 2008:0.126},
    ('po4', 'timeseries', 'max') : {1971:0.204, 2008:0.204},


    }


#for key in ['sos', 'sal','psu']:
for key, a,b in itertools.product(['sos', 'sal','psu'], ['timeseries', ],['min', 'max']): #'clim'
    hard_wired_obs[key, a, b] = hard_wired_obs['so', a, b]

for a,b in itertools.product(['timeseries', ],['min', 'max']): #'clim' #:
    hard_wired_obs['mlotst', a, b] = hard_wired_obs['mld', a, b]



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
        'all': 'All SSPs',
        'obs': 'Observations',
    }
    return sspdict[ssp]


def get_shelve_path(field, pane='timeseries'):
    # from standalone calc ai obs.
    shelve_path = diagtools.folder(['aux', 'obs_shelves'])
    shelve_path += '_'.join([field, pane]) +'.shelve'
    return shelve_path


def fix_chl(cube):
    """
    Several datasets seem to have really strange CHL values,
    so I'm fixing them
    """
    print('pre fix_chl:',cube.data.min(), cube.data.max())
    # cmax = cube.data.max()
    for i in range(3):
        if cube.data.max()<30.:
            cube.data *= 1000.

    print('post fix_chl:',cube.data.min(), cube.data.max())
    return cube


def round_sig(x, sig=2):
    """
    round to x significant figures.
    from https://stackoverflow.com/a/3413529
    """
    return round(x, sig-int(floor(log10(abs(x))))-1)


def timeplot(cube, **kwargs):
    """
    Create a time series plot from the cube.

    Note that this function simple does the plotting, it does not save the
    image or do any of the complex work. This function also takes and of the
    key word arguments accepted by the matplotlib.plt.plot function.
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
    if window  in ['annual', ]:
        window = '1 year'
    window = window.split()
    window_len = int(window[0]) / 2.
    win_units = str(window[1])

    if win_units not in [
            'days', 'day', 'dy', 'months', 'month', 'mn', 'years', 'yrs',
            'year', 'yr'
    ]:
        raise ValueError("Moving average window units not recognised: " +
                         "{}".format(win_units))

    times = cube.coord('time').units.num2date(cube.coord('time').points)

    datetime = diagtools.guess_calendar_datetime(cube)

    output = []

    times = np.array([
        datetime(time_itr.year, time_itr.month, time_itr.day, time_itr.hour,
                 time_itr.minute) for time_itr in times
    ])

    for time_itr in times:
        if win_units in ['years', 'yrs', 'year', 'yr']:
            tmin = datetime(time_itr.year - window_len, time_itr.month,
                            time_itr.day, time_itr.hour, time_itr.minute)
            tmax = datetime(time_itr.year + window_len, time_itr.month,
                            time_itr.day, time_itr.hour, time_itr.minute)

        elif win_units in ['months', 'month', 'mn']:
            tmin = datetime(time_itr.year, time_itr.month - window_len,
                            time_itr.day, time_itr.hour, time_itr.minute)
            tmax = datetime(time_itr.year, time_itr.month + window_len,
                            time_itr.day, time_itr.hour, time_itr.minute)

        elif win_units in ['days', 'day', 'dy']:
            tmin = datetime(time_itr.year, time_itr.month,
                            time_itr.day - window_len, time_itr.hour,
                            time_itr.minute)
            tmax = datetime(time_itr.year, time_itr.month,
                            time_itr.day + window_len, time_itr.hour,
                            time_itr.minute)
        else: assert 0
        arr = np.ma.masked_where((times < tmin) + (times > tmax), cube.data)
        output.append(arr.mean())
    cube.data = np.array(output)
    return cube



def write_csv_ts(cfg, metadatas, data_values, single_model, short_name):
    """
    Write a CSV file for the time series
    """
    out_shelve_dir = diagtools.folder([cfg['work_dir'], 'timeseries_csv', ])


    for (variable_group, short_name, dataset, scenario, ensemble), data in sorted(data_values.items()):
        # if 'all' not in plotting: continue
        out_csv = out_shelve_dir+ '_'.join(['time_series',single_model, variable_group, short_name, dataset, scenario, ensemble])+'.csv'
        if os.path.exists(out_csv):continue

        txt = ' '.join([ '#', short_name, dataset, scenario, ensemble, '\n'])
        txt += 'time, value\n'
        times = [t for t in sorted(data.keys())]
        values = [data[t] for t in times]
        for t,v in zip(times,values):
            txt = ''.join([txt, str(t),', ', str(v),'\n'])

        print('write_csv_ts: save:', out_csv)
        fn = open(out_csv, 'w')
        fn.write(txt)
        fn.close()


        # plt.plot(times, values, ls='-', c=color, lw=0.4) #, label=dataset)

    #assert 0


def find_hist_cube(hist_cubes, dataset, short_name, scenario, ensemble, return_missing=False):
    """
    Tool to figure out which historical cube it is.
    """
    # If the exact match exists, use that:
    if (dataset, short_name, 'historical', ensemble) in hist_cubes:
        return hist_cubes[(dataset, short_name, 'historical', ensemble)]

    # Figure out which historical ensemble members exist:
    sets = []
    for d, s, sce, e in hist_cubes.keys():
        if sce != 'historical': continue
        #print('find_hist_cube: hist_cubes:', len(hist_cubes.keys()))
        if d == dataset:
            #print(d, s, sce, e, ('looking for ', ensemble))
            sets.append((d, s, sce, e))

    if dataset == 'UKESM1-0-LL':
        return hist_cubes[(dataset, short_name, 'historical', ensemble.replace('f2', 'f3'))]
    else:
        print('could not find:', (dataset, short_name, scenario, ensemble), 'but did find ', len(sets), 'ensemble members:', sets)

    if return_missing:
        return (dataset, short_name, 'historical', ensemble)
    assert 0



def make_hist_anomaly_figure(cfg, metadatas, short_name, anomaly_table, future_range, style='std', single_model_marker=True):
    """
    A figure for the historical anomaly in figure 2.
    """
    backgroundcolor='white'
    print('hist, anom', short_name, anomaly_table)

    hist_anom_data = {}
    hist_anom_data['temp'] = ['Sea Surface Temperature',  r'$\degree$ C',]
    hist_anom_data['tos'] = hist_anom_data['temp']
    hist_anom_data['sal'] = ['Salinity', 'PSU', ]
    hist_anom_data['sos'] = hist_anom_data['sal']
    hist_anom_data['mld'] = ['Mixed Layer Depth', 'm', ]
    hist_anom_data['mlotst'] = hist_anom_data['mld']
    hist_anom_data['ph'] = ['pH', '', ]

    hist_anom_data['o2'] = ['Oxygen (500m)', r'mmol m$^{-3}$', ]
    hist_anom_data['nit'] = ['Nitrate', r'$\mu$mol m$^{-3}$', ]
    hist_anom_data['no3'] =  hist_anom_data['nit']
    hist_anom_data['pho'] = ['Phosphate', r'$\mu$mol m$^{-3}$',]
    hist_anom_data['po4'] = hist_anom_data['pho']
    hist_anom_data['chl']  = ['Chlorophyll', r'mg m$^{-3}$', ]
    hist_anom_data['intpp'] = ['Int. Production', r'mmol m$^{-2}$ d$^{-1}$',]


    fig = plt.figure()
    axes = fig.add_subplot(1, 1, 1)
    fig.set_facecolor(backgroundcolor)
    axes.set_facecolor(backgroundcolor)

    fig.set_size_inches(3.0, 3.5)


    print(hist_anom_data[short_name])
    if hist_anom_data[short_name][1] in ['', None]:
            plt.title(hist_anom_data[short_name][0]) # ' Anomaly']))
    else:
        plt.title(''.join([hist_anom_data[short_name][0],', ',  hist_anom_data[short_name][1]]))

    xes = [0.75, 1.5, 2., 2.5, 3.,  ]
    a = []
    b = []
    c = []
    c05 = []
    c95 = []

    # hist_val = 0.
    # hist_stf = 0.

    ssps = ['Historical', 'SSP1-2.6', 'SSP2-4.5', 'SSP3-7.0', 'SSP5-8.5']
    #columns = {'Historical':2, 'SSP1-2.6':3, 'SSP2-4.5':4, 'SSP3-7.0':5, 'SSP5-8.5':6}
    # colours = ['blue', 'green' ,'gold', 'orange', 'red']
    # colours2 = ['navy', 'darkgreen' ,'goldenrod', 'darkorange', 'darkred']

    colours2_rgb = [(23, 60, 102), (0, 173, 207), (247, 148, 32), (231, 29, 37), (149, 27, 30)]
    colours = [(a/256., b/256. ,c/256.) for (a, b, c) in colours2_rgb]
    f = 0.75/256.
    colours2 = [(f*a, f*b, f*c) for (a, b, c) in colours2_rgb]


    tmp = ' '.join(['#',short_name, str(int(future_range[0])),'-', str(int(future_range[1])), '\n'])
    tmp += 'SSP, min, mean, median, max \n'
    for x, ssp in zip(xes, ssps):
        for new_ssp, (mean, std, dat) in anomaly_table.items():
            if sspify(new_ssp) != ssp:
                continue

            if ssp.lower() == 'historical': # no historical here
            #    hist_val = v
            #    hist_stf = s
                continue

            #print(short_name, ssp, v,s)
            a.append(x)
            b.append(mean)
            c.append(std)
            c05.append(np.percentile(dat, 5))
            c95.append(np.percentile(dat, 95))

            if single_model_marker:
                plt.scatter([x for d in dat], dat, color='k', marker='d', zorder=3,)
            tmp+= ','.join([new_ssp, str(np.min(dat)), str(np.mean(dat)), str(np.median(dat)), str(np.max(dat)), ])+'\n'

    fn = open(''.join([short_name, '-',  str(int(future_range[0])),'-', str(int(future_range[1])),'.csv']), 'w')
    fn.write(tmp)
    fn.close()
    print('a:', a,'\nb:', b, '\nc:', c)
    #    #axes.axhline(hist_val, xmin=0.1, xmax=0.9, c=colours2[0])#', zorder=2.)
    axes.axhline(0., linewidth=2.5, xmin=0.1, xmax=0.9, c=colours2[0], zorder=-2.)

    plt.xlim([min(a)-0.6, max(a)+0.6])

    plt.xticks([], [])
    axes.set_xticks([])
    axes.set_xticks([], minor=True)
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['bottom'].set_visible(False)
    plt.subplots_adjust(left=0.20)

    plt.scatter(a, b, color=colours2[1:], marker='_', s=700,zorder=1,) # marker

    if style == 'std':
        plt.errorbar(a, b, yerr=c, ecolor=colours[1:], capsize=0., elinewidth=22.,  fmt="None", zorder=-1) # error bars only (fmt=None)

        impath = diagtools.folder(''.join([cfg['plot_dir'], '/mean_anomaly_plots']))
        impath += ''.join(['hist_anom_', short_name] )

    if style == '5-95':
        plt.errorbar(a, b, yerr=[np.array(b) - np.array(c05), np.array(c95) - np.array(b) ], ecolor=colours[1:], capsize=0., elinewidth=15.,  fmt="None", zorder=-1) # error bars only (fmt=None)

        impath = diagtools.folder(''.join([cfg['plot_dir'], '/mean_anomaly_plots_pc5_pc95']))
        impath += ''.join(['hist_anom_pc_', short_name] )

    impath += '-'.join([str(int(t)) for t in future_range])
    if single_model_marker:
        impath+='_mark'

    print('Saving:', impath)
    plt.savefig(impath+diagtools.get_image_format(cfg), dpi=300., transparent=True)
    plt.savefig(impath+'_white'+diagtools.get_image_format(cfg), dpi=300.) #ransparent=True)

    plt.savefig(impath+'.svg', transparent=True)
    plt.close()


def multi_model_time_series(
        cfg,
        metadatas,
        ts_dict = {},
        moving_average_str='',
        colour_scheme = 'IPCC',
        hist_time_range = None,
        ssp_time_range = None,
        plotting = [ 'means',  '5-95'], #'medians', 'all_models', 'range',
        fig = None,
        ax = None,
        save = False,
        single_model = 'all',
        ):
    """
    Make a time series plot showing several preprocesssed datasets.

    This tool loads several cubes from the files, checks that the units are
    sensible BGC units, checks for layers, adjusts the titles accordingly,
    determines the ultimate file name and format, then saves the image.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.
    single_model:
        set to all, plots all models, otherwise it's just an single model results

    """
    overwrite = True
    impath = diagtools.folder(cfg['plot_dir']+'/individual_panes_timeseries')
    impath += '_'.join(['multi_model_ts', moving_average_str] )
    impath += '_'+'_'.join(plotting)
    if single_model == 'all': pass
    else: impath += '_'+single_model
    #if ukesm == 'not': impath += '_noUKESM'
    impath += diagtools.get_image_format(cfg)
    if not overwrite and os.path.exists(impath): return fig, ax


    # Determine short name & model list:
    variable_groups = {}
    datasets = {}
    short_names = {}
    ensembles = {}
    scenarios = {}
    for variable_group, filenames  in ts_dict.items():
        for fn in sorted(filenames):
            if metadatas[fn]['mip'] in ['Ofx', 'fx']: continue
            dataset = metadatas[fn]['dataset']
            short_name = metadatas[fn]['short_name']
            if dataset in models_to_skip['all']: continue
            if dataset in models_to_skip.get(short_name, {}): continue

            if single_model == 'all': pass
            elif dataset != single_model: continue
            variable_groups[variable_group] = True
            datasets[dataset] = True
            short_names[metadatas[fn]['short_name']] = True
            ensembles[metadatas[fn]['ensemble']] = True
            scenarios[metadatas[fn]['exp']] = True

    if len(short_names.keys()) == 0:
        # no data in time series (ie, tos doesn't exist, but thetao does?)
        if save:
            return
        else:
            return fig, ax

    if len(short_names.keys()) != 1:
        print(short_names.keys())
        assert 0
    short_name = list(short_names.keys())[0]
    if short_name == 'thetao':
        short_name = 'tos'
    if short_name == 'so':
        short_name = 'sos'

    # create a csv of the model contents:
    # short_name, model, scenario, ensemble member
    make_csvs = True
    if make_csvs and single_model=='all':
        # load the netcdfs and populate the shelve dicts
        lines = []
        mod_scen_counts = {}
        total_counts = {}
        models_per_sceanrio = {}
        ensembles_per_sceanrio = {}

        for variable_group, filenames  in ts_dict.items():
            for fn in sorted(filenames):
                if metadatas[fn]['mip'] in ['Ofx', 'fx']: continue
                dataset = metadatas[fn]['dataset']
                short_name = metadatas[fn]['short_name']
                if dataset in models_to_skip['all']: continue
                if dataset in models_to_skip.get(short_name, {}): continue

                short_name = metadatas[fn]['short_name']
                scenario = metadatas[fn]['exp']
                ensemble = metadatas[fn]['ensemble']
                lines.append(', '.join([short_name, dataset, scenario, ensemble,'\n']))

                sds = (short_name, dataset, scenario )
                if sds in mod_scen_counts: mod_scen_counts[sds] +=1
                else: mod_scen_counts[sds] = 1

                ss = (short_name, scenario )
                if ss in total_counts: total_counts[ss] +=1
                else: total_counts[ss] = 1

                if len(models_per_sceanrio.get(ss, [])):
                    models_per_sceanrio[ss][dataset] = True
                else: models_per_sceanrio[ss] = {dataset: True}

                # if ensembles_per_sceanrio.get(ss, 0):
                #     ensembles_per_sceanrio[ss] +=1
                # else: ensembles_per_sceanrio[ss] = 1

        model_table = 'short_name, model, scenario, ensemble member, \n'
        for line in sorted(lines):
            model_table = ''.join([model_table, line])
        model_path  = diagtools.folder(cfg['work_dir']+'/model_table')
        model_path += short_name+'_full.csv'
        mp_fn = open(model_path, 'w')
        mp_fn.write(model_table)
        mp_fn.close()
        print('output:', model_table)
        print('saved model path:', model_path)

        # csv for the individual model scenario count.
        mod_scen_str = 'short_name, model, scenario, count, \n'
        for (short_name, dataset, scenario ) in sorted(mod_scen_counts.keys()):
            count = mod_scen_counts[(short_name, dataset, scenario)]
            line = ', '.join([short_name, dataset, scenario, str(count), '\n'])
            mod_scen_str = ''.join([mod_scen_str, line])
        model_path  = diagtools.folder(cfg['work_dir']+'/model_table')
        model_path += short_name+'_model_scenario_count.csv'
        mp_fn = open(model_path, 'w')
        mp_fn.write(mod_scen_str)
        mp_fn.close()

        mod_scen_str = 'short_name, scenario, model count, ensemble count, \n'
        for (short_name, scenario ) in sorted(models_per_sceanrio.keys()):
            modelcount = models_per_sceanrio[(short_name, scenario)]
            ensemble_count = total_counts[(short_name, scenario)]

            line = ', '.join([short_name, scenario, str(count), str(ensemble_count), '\n'])
            mod_scen_str = ''.join([mod_scen_str, line])
        model_path  = diagtools.folder(cfg['work_dir']+'/model_table')
        model_path += short_name+'_scenarios_counts.csv'
        mp_fn = open(model_path, 'w')
        mp_fn.write(mod_scen_str)
        mp_fn.close()

    ####
    # Load the data for each layer as a separate cube
    out_shelve = diagtools.folder([cfg['work_dir'], 'timeseries_shelves', ])
    if single_model == 'all':
        out_shelve += 'time_series_'+short_name+'.shelve'
    else:
        out_shelve += 'time_series_'+short_name+'_'+single_model+'.shelve'

    model_cubes = {}
    if len(glob.glob(out_shelve+'*')):
       print('loading from shelve:', out_shelve)
       sh = shopen(out_shelve)
       #model_cubes = sh['model_cubes']
       data_values= sh['data_values']
       model_cubes_paths = sh['model_cubes_paths']
       sh.close()
    else:
       model_cubes_paths = {}
       data_values = {}

    changes = 0
    variable_groups = {}

    # load the netcdfs and populate the shelve dicts
    for variable_group, filenames  in ts_dict.items():
        for fn in sorted(filenames):
            if metadatas[fn]['mip'] in ['Ofx', 'fx']: continue
            dataset = metadatas[fn]['dataset']
            short_name = metadatas[fn]['short_name']
            if dataset in models_to_skip['all']: continue
            if dataset in models_to_skip.get(short_name, {}): continue

            if single_model == 'all': pass
            elif dataset != single_model: continue

            short_name = metadatas[fn]['short_name']
            scenario = metadatas[fn]['exp']
            ensemble = metadatas[fn]['ensemble']
            nc_index = (variable_group, short_name, dataset, scenario, ensemble)

            # if already loaded, then skip.
            if model_cubes_paths.get(fn, False): continue
            if data_values.get(nc_index, False): continue

            #print('loading', nc_index, 'fn:', fn)
            # load netcdf as a cube
            cube = iris.load_cube(fn)
            print('\npre:',short_name, dataset, scenario, ensemble, cube.units, cube.data.min(), cube.data.max())

            cube = diagtools.bgc_units(cube, short_name)
            if short_name == 'chl':
                 cube = fix_chl(cube)
            print('post',short_name, dataset, scenario, ensemble, cube.units, cube.data.min(), cube.data.max())
            #model_cubes = add_dict_list(model_cubes, variable_group, cube)

            # Take a moving average, if needed.
            if not moving_average_str: pass
            elif moving_average_str == 'annual':
                cube = annual_statistics(cube.copy(), operator='mean')
                cube = regrid_time(cube, 'yr')
            else:
                cube = moving_average(cube, moving_average_str)

            # load times:
            times = diagtools.cube_time_to_float(cube)

            data = {}
            for t, d in zip(times, cube.data):
                if moving_average_str == 'annual':
                    t = int(t) + 0.5
                data[t] = d
                #data_values = add_dict_list(data_values, t, d)
                #model_data_groups[dataset] = add_dict_list(model_data_groups[dataset], t, d)
            data_values[nc_index] = data
            model_cubes_paths = add_dict_list(model_cubes_paths, fn, True)
            changes+=1

    if changes:
       print('adding ', changes, 'new files to', out_shelve)
       sh = shopen(out_shelve)
       sh['data_values'] = data_values
       sh['model_cubes_paths'] = model_cubes_paths
       sh.close()

    make_csvs_data = False
    if make_csvs_data :
        write_csv_ts(cfg,metadatas, data_values, single_model, short_name)

    if make_csvs and single_model=='all':
        # want a table that includes:
        # field: historical mean 2000-2010 +/1 std, each scenario +/- std
        model_values = {}
        ensembles_dict = {}
        model_ensembles_dict={}
        for (variable_group, short_name, dataset, scenario, ensemble), data in sorted(data_values.items()):
            times = np.array([t for t in sorted(data.keys())])
            values = np.array([data[t] for t in times])
            if scenario == 'historical':
                times = np.ma.masked_outside(times, 2000., 2010.)
            else:
                times = np.ma.masked_outside(times, 2040., 2050.)
            value = np.ma.masked_where(times.mask, values).mean()
            key = (short_name, dataset, scenario)
            model_values = add_dict_list(model_values, key, value)

            model_ensembles_dict = add_dict_list(model_ensembles_dict, key, ensemble)
            ensembles_dict = add_dict_list(ensembles_dict, (short_name, scenario), ensemble)

        scenario_values = {}
        for (short_name, dataset, scenario), means in model_values.items():
            model_mean =  np.mean(means)
            scenario_values = add_dict_list(scenario_values, (short_name, scenario), model_mean)

        header = ['field', ]
        line = [short_name, ]
        for (short_name, scenario) in sorted(scenario_values.keys()):
            header.append(scenario)
            model_means = np.array(scenario_values[(short_name, scenario)])
            try: line.append(' '.join([str(round(model_means.mean(),sigfigs=4)),'\pm', str(round(model_means.std(), sigfigs=3))]))
            except: line.append(' '.join([str(round(model_means.mean(),4)),'\pm', str(round(model_means.std(), 3))]))

        header.append('\n')
        line.append('\n')
        header = ', '.join(header)
        line1 = ', '.join(line  )
        line2 = ' & '.join(line).replace('\pm', '$\pm$').replace('\n', '\\\\ \n')
        csv_str = ''.join([header, line1, line2])
        csv_path  = diagtools.folder(cfg['work_dir']+'/model_table')
        csv_path += 'diff_table_'+short_name+'.csv'
        mp_fn = open(csv_path, 'w')
        mp_fn.write(csv_str)
        mp_fn.close()
        print('csv_str:',csv_str)
        print('saved to path:', csv_path)

        #anomaly table:
        # want a table that includes:
        # field: historical mean anomaly against 1850-1900 for years 2000-2010 +/1 std, each scenario +/- std
        model_values = {}
        ensembles_dict = {}
        model_ensembles_dict={}
        anomalgy_basis = {}
        for (variable_group, short_name, dataset, scenario, ensemble), data in sorted(data_values.items()):
            if scenario.lower().find('historical')==-1: continue
            times = np.array([t for t in sorted(data.keys())])
            values = np.array([data[t] for t in times])
            times = np.ma.masked_outside(times, 1850., 1900.)
            value = np.ma.masked_where(times.mask, values).mean()
            anomalgy_basis[ (short_name, dataset, ensemble[:3])] = value
            print('anomalgy_basis:', (short_name, dataset, ensemble), ':', value)

        for (variable_group, short_name, dataset, scenario, ensemble), data in sorted(data_values.items()):
            times = np.array([t for t in sorted(data.keys())])
            values = np.array([data[t] for t in times])
            if scenario == 'historical':
                times = np.ma.masked_outside(times, 2000., 2010.)
            else:
                times = np.ma.masked_outside(times, 2040., 2050.)

            if (short_name, dataset, ensemble[:3]) in  anomalgy_basis:
                anom =  anomalgy_basis[(short_name, dataset, ensemble[:3])]
            elif (short_name, dataset, 'r1i') in  anomalgy_basis:
                anom =  anomalgy_basis[(short_name, dataset, 'r1i')]
            else:
                print('unalbe to include:', (short_name, dataset, ensemble[:3]) )
                continue
            value = np.ma.masked_where(times.mask, values).mean() - anom
            key = (short_name, dataset, scenario)
            model_values = add_dict_list(model_values, key, value)

            model_ensembles_dict = add_dict_list(model_ensembles_dict, key, ensemble)
            ensembles_dict = add_dict_list(ensembles_dict, (short_name, scenario), ensemble)
            #print('mean:', variable_group, short_name, dataset, scenario, ensemble, value)
        #assert 0

        scenario_values = {}
        for (short_name, dataset, scenario), means in model_values.items():
            model_mean =  np.mean(means)
            scenario_values = add_dict_list(scenario_values, (short_name, scenario), model_mean)

        header = ['field', ]
        line = [short_name, ]
        for (short_name, scenario) in sorted(scenario_values.keys()):
            header.append(scenario)
            model_means = np.array(scenario_values[(short_name, scenario)])
            try: line.append(' '.join([str(round(model_means.mean(),sigfigs=4)),'\pm', str(round(model_means.std(), sigfigs=3))]))
            except:line.append(' '.join([str(round(model_means.mean(),4)),'\pm', str(round(model_means.std(), 3))]))

        header.append('\n')
        line.append('\n')
        header = ', '.join(header)
        line   = ', '.join(line  )

        csv_str = ''.join([header, line])
        csv_path  = diagtools.folder(cfg['work_dir']+'/model_table')
        csv_path += 'anomaly_table_'+short_name+'.csv'
        mp_fn = open(csv_path, 'w')
        mp_fn.write(csv_str)
        mp_fn.close()
        print('csv_str:',csv_str)
        print('saved to path:', csv_path)


        # historiocal anomaly table:
        # want a table that includes:
        # field: historical mean anomaly against 1850-1900 for years 2000-2010 +/1 std, each scenario +/- std
        model_values = {}
        ensembles_dict = {}
        model_ensembles_dict={}
        anomalgy_basis = {}
        # Calculate the anomaly against region (2000-2010 hist)
        #future_range=[2090., 2100.]
        future_range=[2040., 2050.]

        for (variable_group, short_name, dataset, scenario, ensemble), data in sorted(data_values.items()):
            if scenario.lower().find('historical')==-1: continue
            times = np.array([t for t in sorted(data.keys())])
            values = np.array([data[t] for t in times])
            times = np.ma.masked_outside(times, 2000., 2010.)
            value = np.ma.masked_where(times.mask, values).mean()
            anomalgy_basis[ (short_name, dataset, ensemble[:3])] = value
            print('anomalgy_basis:', (short_name, dataset, ensemble), ':', value)

        for (variable_group, short_name, dataset, scenario, ensemble), data in sorted(data_values.items()):
            times = np.array([t for t in sorted(data.keys())])
            values = np.array([data[t] for t in times])
            if scenario == 'historical':
                times = np.ma.masked_outside(times, 2000., 2010.)
            else:
                times = np.ma.masked_outside(times, future_range[0], future_range[1])

            if (short_name, dataset, ensemble[:3]) in  anomalgy_basis:
                anom =  anomalgy_basis[(short_name, dataset, ensemble[:3])]
            elif (short_name, dataset, 'r1i') in  anomalgy_basis:
                anom =  anomalgy_basis[(short_name, dataset, 'r1i')]
            else:
                print('unalbe to include:', (short_name, dataset, ensemble[:3]) )
                continue
            # calculate the single model ensemble diff against the historical anomaly.
            value = np.ma.masked_where(times.mask, values).mean() - anom
            key = (short_name, dataset, scenario)
            model_values = add_dict_list(model_values, key, value)

            model_ensembles_dict = add_dict_list(model_ensembles_dict, key, ensemble)
            ensembles_dict = add_dict_list(ensembles_dict, (short_name, scenario), ensemble)

        # calculate the single model means:
        scenario_values = {}
        for (short_name, dataset, scenario), means in model_values.items():
            model_mean =  np.mean(means)
            scenario_values = add_dict_list(scenario_values, (short_name, scenario), model_mean)

        #  anomaly_table: this is the data used for figure 2 (summary)
        # generate the output line: # This
        header = ['field', ]
        line = [short_name, ]
        line1 = [short_name, ]
        line2 = [short_name, ]
        line3 = [short_name, ]
        anomaly_table = {}
        for (short_name, scenario) in sorted(scenario_values.keys()):
            header.append(scenario)
            model_means = np.array(scenario_values[(short_name, scenario)])
            try:
                value1 = str(round(model_means.mean(),sigfigs=4))
                value2 = str(round(model_means.std(), sigfigs=3))
                value3 = str(round(model_means.min(), sigfigs=4))
                value4 = str(round(model_means.max(), sigfigs=4))
            except:
                value1 = str(round(model_means.mean(), 4))
                value2 = str(round(model_means.std(), 3))
                value3 = str(round(model_means.min(), 4))
                value4 = str(round(model_means.max(), 4))


            anomaly_table[scenario] = (model_means.mean(), model_means.std(), model_means)

            line.append(''.join([' ', value1, ' $\pm$ ', value2,]))
            line1.append(''.join(['"', value1, ' $\pm$ ', value2, '"',]))
            line2.append(''.join(['"', value1, ' + ', value4, ' - ', value3, '"',]))
            l3 = [str(round(v, 3)) for v in model_means]
            l3.insert(0,value1)
            line3.append(' '.join(l3))
            #try: line.append(' '.join([str(round(model_means.mean(),sigfigs=4)),'\pm', str(round(model_means.std(), sigfigs=3))]))
            #except:line.append(' '.join([str(round(model_means.mean(),4)),'\pm', str(round(model_means.std(), 3))]))

        header.append('\n')
        line.append('\n')
        line1.append('# summary hist anom\n')
        line2.append('# min max ranges \n')
        line3.append('# single job means \n')
        header = ', '.join(header)
        line   = ' & '.join(line  ) # latex
        line1 = ', '.join(line1) # python
        line2 = ', '.join(line2) # python
        line3 = ', '.join(line3) # python

        csv_str = ''.join([header, line, line1, line2, line3])
        csv_path  = diagtools.folder(cfg['work_dir']+'/model_table')
        csv_path += 'historical_anomaly_table_'+short_name+'.csv'
        mp_fn = open(csv_path, 'w')
        mp_fn.write(csv_str)
        mp_fn.close()
        print('csv_str:',csv_str)
        print('saved to path:', csv_path)

        make_hist_anomaly_figure(cfg, metadatas, short_name, anomaly_table, future_range, style='std')
        for single_model_marker in [True, False]:
            make_hist_anomaly_figure(cfg, metadatas, short_name, anomaly_table, future_range, style='5-95', single_model_marker=single_model_marker)

        # number of ensembles
        # Number of models, number of ensembles, (ensembles per model)
        header = ['field', ]
        line = [short_name, ]
        line2 = [' ', ]
        # model_ensembles_dict = add_dict_list(model_values, key, ensemble)
        # ensembles_dict = add_dict_list(ensembles_dict, key, ensemble)

        for (short_name, scenario) in sorted(ensembles_dict.keys()):
            model_ensembles = {}
            for (sn1, d1, sc1) in sorted(model_ensembles_dict.keys()):
                if short_name != sn1: continue
                if sc1 != scenario: continue
                model_ensembles[d1] = model_ensembles_dict[(sn1, d1, sc1)]

            num_models = str(len(model_ensembles.keys()))
            num_ensembles = str(len(ensembles_dict[(short_name, scenario)]))
            num_mod_ensembles = [str(len(model_ensembles[mod])) for mod in sorted(model_ensembles.keys())]
            header.append(scenario)
            line.append(''.join([num_models, ' (',num_ensembles,')']))
            line2.append('-'.join(num_mod_ensembles))

        header.append('\n')
        line.append('\n')
        line2.append('\n')

        header = ', '.join(header)
        line1   = ', '.join(line)
        line1ltc   = ' & '.join(line).replace('\n', '\\\\ \n')

        line2  = ', '.join(line2)

        csv_str = ''.join([header, line1, line1ltc, line2])
        csv_path  = diagtools.folder(cfg['work_dir']+'/model_table')
        csv_path += 'ensemble_count_table_'+short_name+'.csv'
        mp_fn = open(csv_path, 'w')
        mp_fn.write(csv_str)
        mp_fn.close()
        print('csv_str:',csv_str)
        print('saved to path:', csv_path)

    # Make the figure
    if fig is None or ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        save = True
    else:
        plt.sca(ax)

    if single_model == 'all':
        if short_name in ['tos', 'thetao', 'sal', 'sos', 'so',  'chl', 'no3', 'nit', 'po4', 'pho', 'ph', 'pH', ]:
               title = 'Surface Time Series'
        elif short_name in ['o2',]:
               title = 'Time Series (500m)'
        else:  title = 'Time Series'

    else:
        title = ' '.join([single_model,  'Time series'])

    z_units = ''
    plot_details = {}
    if colour_scheme in ['viridis', 'jet']:
        cmap = plt.cm.get_cmap(colour_scheme)

    # Plot every single ensemble member
    # plottings:
    #     all: every single ensemble member shown as thin line.
    #     model_means: the mean of every model ensmble is shown as a thicker dotted line
    #     model_median: the median of every model ensmble is shown as a thicker dotted line
    #     model_range: the range of every model ensmble is shown as a coloured band (or line if only one)
    #     model_5_95: the range of every model ensmble is shown as a coloured band (or line if only one)
    #     Global_mean: the mean of every model_means is shown as a thick solid line
    #     Global_range: the mean of every model_means is shown as a thick solid line

    for (variable_group, short_name, dataset, scenario, ensemble), data in sorted(data_values.items()):
        if 'all' not in plotting: continue
        if colour_scheme in ['viridis', 'jet']:
            if len(metadata) > 1:
                color = cmap(index / (len(metadata) - 1.))
            else:
                color = 'blue'

        if colour_scheme in 'IPCC':
            color = ipcc_colours[scenario]

        times = [t for t in sorted(data.keys())]
        values = [data[t] for t in times]

        plt.plot(times, values, ls='-', c=color, lw=0.4) #, label=dataset)

    datas_ds = {}
    mm_data = {}

    # single model lines:
    for dataset_x,scenario_x in itertools.product(sorted(datasets.keys()), sorted(scenarios.keys())):
        if not len(set(plotting) & set(
            ('model_means', 'model_median', 'model_range', 'model_5_95', 'Global_mean', 'Global_range', ))):
            continue # none needed.

        if colour_scheme in 'IPCC':
            color = ipcc_colours[scenario_x]
            fill_color = ipcc_colours_dark[scenario_x]
        dat_scen_data = {}
        for (variable_group, short_name, dataset, scenario, ensemble), data in sorted(data_values.items()):
            if dataset_x!= dataset: continue
            if scenario_x != scenario: continue
            for t,d in data.items():
                dat_scen_data = add_dict_list(dat_scen_data, t, d)

        if not len(dat_scen_data):
            print('nothinhg found for ', dataset_x, scenario_x)
            continue

        # calculate model mean.
        if len(set(plotting) & set(('model_means', 'Global_mean', 'Global_range'))):
            model_mean = {}
            for t, ds in dat_scen_data.items():
                model_mean[t] = np.mean(ds)
            # model_mean = {t:d for t,d in zip(times, model_mean)}
            datas_ds[(dataset_x, scenario_x)] = model_mean
            print('calculate model mean',(dataset_x, scenario_x))

        if 'model_means' in plotting:
            model_mean = datas_ds.get((dataset_x, scenario_x), {})
            times = [t for t in sorted(model_mean.keys())]
            model_mean = [model_mean[t] for t in times]
            plt.plot(times, model_mean, ls=':', c=color, lw=2.) #, label=dataset)
            mm_data[(variable_group, short_name, dataset, scenario, ensemble)] = {t:d for t,d in zip(times, model_mean)}

        if 'model_medians' in plotting:
            times = [t for t in sorted(dat_scen_data.keys())]
            model_median = []
            for t, ds in dat_scen_data.items():
                model_median.append(np.median(ds))
            plt.plot(times, model_median, ls=':', c=color, lw=2.) #, label=dataset)

        if 'model_range' in plotting: # for each model plots the entire range at each point in time.
            times = [t for t in sorted(dat_scen_data.keys())]
            model_mins = []
            model_maxs = []
            for t, ds in dat_scen_data.items():
                model_mins.append(np.min(ds))
                model_maxs.append(np.max(ds))

            if len(datasets.keys()) > 15: alpha = 0.1
            elif len(datasets.keys()) > 10: alpha = 0.2
            elif len(datasets.keys()) > 5: alpha = 0.3
            else: alpha = 0.35

            if model_mins == model_maxs:
                plt.plot(times, model_mins, ls='-', c=color, lw=2., alpha=alpha)
                mm_data[(variable_group, short_name, dataset_x, scenario_x, 'ensemble_min_max')] = {t:d for t,d in zip(times, model_mins)}
            else:
                plt.fill_between(times, model_mins, model_maxs, color=fill_color, ec=None, alpha=alpha)
                mm_data[(variable_group, short_name, dataset_x, scenario_x, 'ensemble_min')] = {t:d for t,d in zip(times, model_mins)}
                mm_data[(variable_group, short_name, dataset_x, scenario_x, 'ensemble_max')] = {t:d for t,d in zip(times, model_maxs)}

        if 'model_5_95' in plotting:
            assert 0
            # times = [t for t in sorted(dat_scen_data.keys())]
            # model_mins = []
            # model_maxs = []
            # for t, ds in dat_scen_data.items():
            #     model_mins.append(np.min(ds))
            #     model_maxs.append(np.max(ds))
            # if model_mins == model_maxs:
            #     plt.plot(times, model_mins, ls=':', c=color, lw=2.)
            # else:
            #     plt.fill_between(times, model_mins, model_maxs, color = fill_color, ec=None, alpha=0.3)

    # global model lines:
    for scenario_x in sorted(scenarios):
        if not len(set(plotting) & set(('Global_mean', 'Global_range', ))):
            continue # none needed.

        if colour_scheme in 'IPCC':
            color = ipcc_colours[scenario_x]
            fill_color = ipcc_colours_dark[scenario_x]

        global_model_means = {}
        for dataset_x in datasets:
            model_mean = datas_ds.get((dataset_x, scenario_x), {})
            for t,d in model_mean.items():
                global_model_means = add_dict_list(global_model_means, t, d)

        if not len(dat_scen_data):
            print('nothinhg found for ', dataset_x, scenario_x)

        if 'Global_mean' in plotting:
            times = [t for t in sorted(global_model_means.keys())]
            global_mean = [np.mean(global_model_means[t]) for t in times]
            color = ipcc_colours[scenario_x]

            plt.plot(times, global_mean, ls='-', c='black', lw=2.5) #, label=dataset)
            plt.plot(times, global_mean, ls='-', c=color, lw=2.) #, label=dataset)
            mm_data[(variable_group, short_name, 'Multimodel', scenario_x, 'mean')] = {t:d for t,d in zip(times, global_mean)}

        if 'Global_range' in plotting:
            #print(plotting)

            times = [t for t in sorted(global_model_means.keys())]
            model_mins = []
            model_maxs = []
            for t, ds in global_model_means.items():
                #print('global_range',scenario_x, t, ds)
                model_mins.append(np.min(ds))
                model_maxs.append(np.max(ds))

            if model_mins == model_maxs:
                plt.plot(times, model_mins, ls='-', c=color, lw=3.,alpha=0.3,)
                mm_data[(variable_group, short_name, 'Multimodel', scenario_x, 'minmax')] = {t:d for t,d in zip(times, model_mins)}

            else:
                plt.fill_between(times, model_mins, model_maxs, color= fill_color, ec=None, alpha=0.3)
                mm_data[(variable_group, short_name, 'Multimodel', scenario_x, 'min')] = {t:d for t,d in zip(times, model_mins)}
                mm_data[(variable_group, short_name, 'Multimodel', scenario_x, 'max')] = {t:d for t,d in zip(times, model_maxs)}

    ax.set_xlim(1846., 2104.)
    make_csvs_data = True
    if make_csvs_data :
        write_csv_ts(cfg,metadatas, mm_data, single_model, short_name)
        write_csv_ts(cfg,metadatas, data_values, single_model, short_name)

    plot_obs = True
    if plot_obs:
        shpath = get_shelve_path(short_name, 'ts')
        sh = shopen(shpath)
        # times, data = sh['times'], sh['data']
        if hard_wired_obs.get((short_name, 'timeseries', 'min'), False):
            ctimes = [t for t in hard_wired_obs[(short_name, 'timeseries', 'min')]]
            mins =  [hard_wired_obs[(short_name, 'timeseries', 'min')][t] for t in ctimes]
            maxs =  [hard_wired_obs[(short_name, 'timeseries', 'max')][t] for t in ctimes]
            plt.fill_between(ctimes, mins, maxs, fc=(0,0,0,0), ec='k',lw=2. )

        if 'annual_times' in sh.keys():
            annual_times, annual_data = sh['annual_times'], sh['annual_data']
            plt.plot(annual_times, annual_data, ls='-', c='k', lw=1.5)
        # clim = sh['clim']
        sh.close()

    y_lims = ax.get_ylim()
    if hist_time_range:
        plt.axvspan(2000., 2010., color= 'k', alpha=0.3,zorder=0) #abel = 'Historical period')
#        plt.fill_betweenx(y_lims, hist_time_range[0], hist_time_range[1], color= 'k', alpha=0.4, label = 'Historical period')
    if ssp_time_range:
        plt.axvspan(2040., 2050., color= 'purple', alpha=0.3,zorder=0) # label = 'SSP period')
#        plt.fill_betweenx(y_lims, ssp_time_range[0], ssp_time_range[1], color= 'purple', alpha=0.4, label = 'SSP period')

#    if short_name in ['o2', 'no3', 'intpp', 'mlotst',]: #'ph'
#        legd = plt.legend(loc='upper left')
#    else:
#        legd = plt.legend(loc='lower left')
#    legd.draw_frame(False)
#    legd.get_frame().set_alpha(0.)

    # Add title, legend to plots
    plt.title(title)
    plt.xlabel('Year') #, fontsize='small')
    plt.ylabel(short_suptitles[short_name])# , fontsize='small')

    # Saving files:
    if save:
        # Resize and add legend outside thew axes.
        fig.set_size_inches(9., 6.)
        diagtools.add_legend_outside_right(
            plot_details, plt.gca(), column_width=0.15)

        logger.info('Saving plots to %s', impath)
        plt.savefig(impath)
        plt.close()
    else:
        return fig, ax


def save_clim_csv(cfg, keys, times, data):
    """
    Save clim data as a csv.
    """
    out_shelve_dir = diagtools.folder([cfg['work_dir'], 'clim_csv', ])
    out_csv = out_shelve_dir + '_'.join(keys)+'.csv'
    if os.path.exists(out_csv):return

    txt = '# '+ ' '.join(keys) + '\n'
    if 'historical' in keys:
        txt += '# time range: 2000-2010\n'
    else:
        txt += '# time range: 2040-2050\n'
    txt += 'month, value\n'

    for t,v in zip(times, data):
        txt = ''.join([txt, str(t),', ', str(v),'\n'])

    print('save_clim_csv: save:', out_csv)
    fn = open(out_csv, 'w')
    fn.write(txt)
    fn.close()


def multi_model_clim_figure(
        cfg,
        metadatas,
        ts_dict = {},
        figure_style = 'plot_all_years',
        hist_time_range = [1990., 2000.],
        ssp_time_range = [2040., 2050.],
        plotting = [ 'means',  '5-95', ], #'all_models'] #'medians', , 'range',
        fig = None,
        ax = None,
        save=False,
        single_model='all'
    ):
    """
    produce a monthly climatology figure.
    """
    if fig is None or ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        save=True
    plt.sca(ax)

    ####
    short_names = {}
    models = {}
    for variable_group, filenames  in ts_dict.items():
        for fn in sorted(filenames):
            if metadatas[fn]['mip'] in ['Ofx', 'fx']: continue
            dataset = metadatas[fn]['dataset']
            short_name = metadatas[fn]['short_name']
            if dataset in models_to_skip['all']: continue
            if dataset in models_to_skip.get(short_name, {}): continue

            if single_model == 'all': pass
            elif metadatas[fn]['dataset'] != single_model: continue
            model = metadatas[fn]['dataset']
            models[model] = True
            short_names[metadatas[fn]['short_name']] = True

    if len(short_names.keys()) == 0:
        # no data in time series (ie, tos doesn't exist, but thetao does?)
        if save:
            return
        else:
            return fig, ax

    if len(short_names) > 1:#
        assert 0
    short_name = list(short_names.keys())[0]

    ####
    # Load the data for each layer as a separate cube
    out_shelve = diagtools.folder([cfg['work_dir'], 'clim_shelves', ])
    if single_model == 'all':
        out_shelve += 'clim_'+short_name+'.shelve'
    else:
        out_shelve += 'clim_'+short_name+'_'+single_model+'.shelve'

    changes = 0

    if len(glob.glob(out_shelve+'*')):
       print('loading from shelve:', out_shelve)
       sh = shopen(out_shelve)
       #model_cubes = sh['model_cubes'] # these aren't cubes as you can't shelve a cube.

       model_cubes= sh['model_cubes']
       omov_cubes= sh['omov_cubes']
       model_cubes_paths = sh['model_cubes_paths']
       omov_cubes_paths = sh['omov_cubes_paths']
       units = sh['units']
       sh.close()
    else:
        model_cubes = {}
        omov_cubes = {}
        omov_cubes_paths = {}
        model_cubes_paths = {}
        units = ''

    for variable_group, filenames  in ts_dict.items():
        for fn in sorted(filenames):
            print(variable_group, fn)
            #assert 0
            if metadatas[fn]['mip'] in ['Ofx', 'fx']: continue
            dataset = metadatas[fn]['dataset']
            short_name = metadatas[fn]['short_name']
            if dataset in models_to_skip['all']: continue
            if dataset in models_to_skip.get(short_name, {}): continue

            #if ukesm == 'only' and dataset != 'UKESM1-0-LL': continue
            if single_model == 'all': pass
            elif metadatas[fn]['dataset'] != single_model: continue

            model = metadatas[fn]['dataset']
            scenario = metadatas[fn]['exp']
            short_name = metadatas[fn]['short_name']

            if fn in model_cubes_paths.get(variable_group,[]) and fn in omov_cubes_paths.get((variable_group, model), []):
                print('already loaded', fn)
                continue

            print('loading', fn)
            cube = iris.load_cube(fn)

            cube = diagtools.bgc_units(cube, metadatas[fn]['short_name'])
            if short_name == 'chl':
                 cube = fix_chl(cube)

            units = cube.units

            if not cube.coords('year'):
                iris.coord_categorisation.add_year(cube, 'time')

            if not cube.coords('month_number'):
                iris.coord_categorisation.add_month_number(cube, 'time', name='month_number')

            if scenario == 'historical':
                cube = extract_time(cube, hist_time_range[0], 1, 1, hist_time_range[1], 12, 31)
            else:
                cube = extract_time(cube, ssp_time_range[0], 1, 1, ssp_time_range[1], 12, 31)

            # Only interested in the mean over the time range provided.
            cube = cube.aggregated_by(['month_number', ], iris.analysis.MEAN)
            #cube_min = cube.copy().aggregated_by(['month_number', ], iris.analysis.MIN)
            #cube_max = cube.copy().aggregated_by(['month_number', ], iris.analysis.MAX)
            cube = {'month_number': cube.coord('month_number').points, 'data': cube.data}

            model_cubes = add_dict_list(model_cubes, variable_group, cube)
            model_cubes_paths = add_dict_list(model_cubes_paths, variable_group, fn)

            omov_cubes = add_dict_list(omov_cubes, (variable_group, model ), cube)
            omov_cubes_paths= add_dict_list(omov_cubes_paths, (variable_group, model ), fn)
            changes+=1

    if changes > 0:
       print('Saving new shelve:', out_shelve)
       sh = shopen(out_shelve)
       sh['model_cubes'] = model_cubes
       sh['omov_cubes'] = omov_cubes
       sh['model_cubes_paths'] = model_cubes_paths
       sh['omov_cubes_paths'] = omov_cubes_paths
       sh['units'] = units
       sh.close()

    if {'OMOC_modelmeans','OneModelOneVote','OMOC_modelranges'}.intersection(set(plotting)):

        omoc_means = {}
        variable_groups = {}
        for (variable_group, model), cubes in omov_cubes.items():
            variable_groups[variable_group] = True
            data_values = {} # one of these for each model and scenario.
            for i, cube in enumerate(cubes):
                print('calculating:', i, (variable_group, model), i, 'of', len(cubes))
                fn = omov_cubes_paths[(variable_group, model)][i]
                if metadatas[fn]['mip'] in ['Ofx', 'fx']: continue
                scenario = metadatas[fn]['exp']

                months =  cube['month_number']
                for t, d in zip(months, cube['data']):
                    data_values = add_dict_list(data_values, t, d)
            omoc_means[(variable_group, model, scenario)] = {t: np.mean(data_values[t]) for t in months}
            color = ipcc_colours[scenario]

            if 'OMOC_modelmeans' in plotting:
                dat = omoc_means[(variable_group, model, scenario)]
                times = sorted(dat.keys())
                mean = [np.mean(dat[t]) for t in times]
                color = ipcc_colours[scenario]
                plt.plot(times, mean, ls='-', c=color, lw=1.)
                save_clim_csv(cfg, ['clim', variable_group, model, scenario, 'OMOC_modelmeans'], times, mean)

            if 'OMOC_modelranges' in plotting:
                dat = omoc_means[(variable_group, model, scenario)]
                times = sorted(dat.keys())
                mins = [np.min(dat[t]) for t in times]
                maxs = [np.max(dat[t]) for t in times]
                if mins == maxs:
                    plt.plot(times, mean, ls='-', c=color, lw=1.5)
                    save_clim_csv(cfg, ['clim', variable_group, model, scenario, 'OMOC_modelranges'], times, mean)

                else:
                    plt.fill_between(times, mins, maxs, color=color, alpha=0.15)
                    save_clim_csv(cfg, ['clim', variable_group, model, scenario, 'OMOC_modelmin'], times, mins)
                    save_clim_csv(cfg, ['clim', variable_group, model, scenario, 'OMOC_modelmax'], times, maxs)

        for variable_group_master in variable_groups.keys():
            omoc_mean = {}

            for (variable_group, model, scenario),model_mean in omoc_means.items():
                if variable_group!= variable_group_master: continue
                scen = scenario
                for t, d in model_mean.items():
                    omoc_mean = add_dict_list(omoc_mean, t, d)

            times = sorted(omoc_mean.keys())
            mean = [np.mean(omoc_mean[t]) for t in times]
            color = ipcc_colours[scen]
            if 'OneModelOneVote' in plotting:
                plt.plot(times, mean, ls='-', c=color, lw=2.)
                save_clim_csv(cfg, ['clim', variable_group_master, 'multimodel', scen, 'OneModelOneVote'], times, mean)

    #assert 0

    #labels = []
    for variable_group, cubes in model_cubes.items():
        data_values = {}
        scenario = ''
        for i, cube in enumerate(cubes):
            fn = model_cubes_paths[variable_group][i]
            if metadatas[fn]['mip'] in ['Ofx', 'fx']: continue
            scenario = metadatas[fn]['exp']

            dataset = metadatas[fn]['dataset']
            short_name = metadatas[fn]['short_name']
            if dataset in models_to_skip['all']: continue
            if dataset in models_to_skip.get(short_name, {}): continue

            months =  cube['month_number']
            for t, d in zip(months, cube['data']):
                data_values = add_dict_list(data_values, t, d)
                #years_range = sorted({yr:True for yr in cube.coord('year').points}.keys())
            color = ipcc_colours[scenario]
            label =  scenario

            if 'all_models' in plotting:
                plt.plot(months, cube['data'], ls='-', c=color, lw=0.6)

        times = sorted(data_values.keys())
        if 'means' in plotting:
            mean = [np.mean(data_values[t]) for t in times]
            plt.plot(times, mean, ls='-', c=color, lw=2.)

        if 'medians' in plotting:
            mean = [np.median(data_values[t]) for t in times]
            plt.plot(times, mean, ls='-', c=color, lw=2.)
#
        if 'range' in plotting:
            times = sorted(data_values.keys())
            mins = [np.min(data_values[t]) for t in times]
            maxs = [np.max(data_values[t]) for t in times]
            plt.fill_between(times, mins, maxs, color=color, alpha=0.15)

        if '5-95' in plotting:
            times = sorted(data_values.keys())
            mins = [np.percentile(data_values[t], 5.) for t in times]
            maxs = [np.percentile(data_values[t], 95.) for t in times]
            plt.fill_between(times, mins, maxs, color=color, alpha=0.15)

    time_str = '_'.join(['-'.join([str(t) for t in hist_time_range]), 'vs',
                         '-'.join([str(t) for t in ssp_time_range])])

    if save and single_model != 'all':
        plt.title('Climatology - '+single_model)
    else:
        if short_name in ['tos', 'thetao', 'sal', 'sos', 'so',  'chl', 'no3', 'nit', 'po4', 'pho', 'ph', 'pH', ]:
               title = 'Surface Climatology'
        elif short_name in ['o2',]:
               title = 'Climatology (500m)'
        else:  title = 'Climatology'

        plt.title(title)

    plot_obs = True
    if plot_obs and short_name not in ['o2', ]:
        if hard_wired_obs.get((short_name, 'clim', 'min'), False):
            ctimes = [t for t in hard_wired_obs[(short_name, 'clim', 'min')]]
            mins =  [hard_wired_obs[(short_name, 'clim', 'min')][t] for t in ctimes]
            maxs =  [hard_wired_obs[(short_name, 'clim', 'max')][t] for t in ctimes]
            plt.fill_between([times[0], times[-1]], mins, maxs, edgecolor='k', linewidth=2., facecolor=(0.,0.,0.,0.,)) #alpha=0.15)

        elif short_name in ['so' ,'sos', 'sal', 'psu', 'mld','mlotst', 'po4','pho','nit', 'no3']:
           if  short_name in ['so' ,'sos', 'sal', 'psu']:
               ob_csv = './aux/obs_csv/sos_clim.csv'
           if short_name in ['mld', 'mlotst']:
               ob_csv = './aux/obs_csv/mld_clim.csv'
           if short_name in ['po4','pho']:
               ob_csv = './aux/obs_csv/po4_clim.csv'
           if short_name in ['nit', 'no3']:
               ob_csv = './aux/obs_csv/no3_clim.csv'

           times = []
           clim = []
           with open(ob_csv) as file:
               for line in file:
                   print(line)
                   if line[0] == '#': continue
                   line = line.split(',')
                   if len(line) !=2: continue
                   times.append(float(line[0]))
                   clim.append(float(line[1].replace('\n','')))
           plt.plot(times, clim, ls='-', c='k', lw=2.2, zorder=3)
           print(times, clim)
        else:

            shpath = get_shelve_path(short_name, 'ts')
            sh = shopen(shpath)
            if 'clim' in sh.keys():
                clim = sh['clim']
                plt.plot(times, clim, ls='-', c='k', lw=1.5)
            sh.close()
    # Ticks!
    ax.set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,])
    ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D']) #, fontsize='small')
    ax.set_xlim([0.8,12.2])

    plt.xlabel('Month', fontsize='small')
    ax.tick_params(axis='both', labelsize='small')

    plt.ylabel(short_suptitles[short_name], fontsize='small')

    if save:
        plt.legend()
        # save and close.
        path = diagtools.folder(cfg['plot_dir']+'/multi_model_clim')
        path += '_'.join(['multi_model_clim', time_str])
        path += '_' +'_'.join(plotting)
        if single_model == 'all': pass
        else: path+= '_'+single_model
        path += diagtools.get_image_format(cfg)

        logger.info('Saving plots to %s', path)
        plt.savefig(path)
        plt.close()
    else:
        return fig, ax





###################################

def extract_depth_range(data, depths, drange='surface', threshold=-1000.):
    """
    Extract either the top 1000m or everything below there.

    drange can be surface or depths.

    Assuming everything is 1D at this stage!
    """
    depths = np.array(depths)
    data = np.ma.array(data)
    data = np.ma.masked_invalid(data)
    if drange == 'surface':
        data = np.ma.masked_where(data.mask + (depths < threshold), data)
    elif drange == 'depths':
        data = np.ma.masked_where(data.mask + (depths > threshold), data)
    return data


def plot_z_line(depths, data, ax0, ax1, ls='-', c='blue', lw=2., label = '', zorder=1):

    for ax, drange in zip([ax0, ax1], ['surface', 'depths']):
        data2 = extract_depth_range(data, depths, drange=drange)
        if not len(data2.compressed()): continue
        ax.plot(data2, depths,
            lw=lw,
            ls=ls,
            c=c,
            label= label, zorder= zorder)
    return ax0, ax1


def plot_z_area(depths, data_min, data_max, ax0, ax1, color='blue', alpha=0.15, label = '', zorder=1):
    for ax, drange in zip([ax0, ax1], ['surface', 'depths']):
        data_min2 = extract_depth_range(data_min, depths, drange=drange)
        data_max2 = extract_depth_range(data_max, depths, drange=drange)
        ax.fill_betweenx(depths, data_min2, data_max2,
            color=color,
            alpha = alpha,
            label= label,
            zorder=zorder)
    return ax0, ax1

def save_profile_csv(cfg, keys, depths, data):
    """
    Save profile data as a csv.
    """
    out_shelve_dir = diagtools.folder([cfg['work_dir'], 'profile_csv', ])
    out_csv = out_shelve_dir + '_'.join(keys)+'.csv'
    #f os.path.exists(out_csv):return

    txt = '# '+ ' '.join(keys) + '\n'
    if 'historical' in keys:
        txt += '# time range: 2000-2010\n'
    else:
        txt += '# time range: 2040-2050\n'
    txt += 'depth, value\n'

    for z, v in zip(depths, data):
        txt = ''.join([txt, str(z),', ', str(v),'\n'])

    print('save_profile_csv: save:', out_csv)
    fn = open(out_csv, 'w')
    fn.write(txt)
    fn.close()



def make_multi_model_profiles_plotpair(
        cfg,
        metadatas,
        profile_fns = {},
        #short_name,
        obs_metadata={},
        obs_filename='',
        hist_time_range = [1990., 2000.],
        ssp_time_range = [2040., 2050.],
        figure_style = 'difference',
        plotting = [ 'means_split', '5-95_split', ], #'all_models'] #'medians', , 'range',
        fig = None,
        ax = None,
        save = False,
        draw_legend=True,
        single_model = 'all',

    ):
    """
    Make a profile plot for an individual model.

    The optional observational dataset can also be added.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.
    filename: str
        The preprocessed model file.
    obs_metadata: dict
        The metadata dictionairy for the observational dataset.
    obs_filename: str
        The preprocessed observational dataset file.

    """
    if not profile_fns:
        return fig, ax

    # Load cube and set up units
    if fig is None or ax is None:
        save = True
        fig = plt.figure()

        gs = matplotlib.gridspec.GridSpec(ncols=1, nrows=1) # master
        gs0 =gs[0,0].subgridspec(ncols=2, nrows=2, height_ratios=[2, 1], hspace=0.)
        ax0=fig.add_subplot(gs0[0,0]) # surface
        ax1=fig.add_subplot(gs0[1,0]) # depths

        #gs1 =gs[0,1].subgridspec(ncols=1, nrows=2, height_ratios=[2, 1], hspace=0.)
        ax2=fig.add_subplot(gs0[0,1]) # surface
        ax3=fig.add_subplot(gs0[1,1]) # depths

    else:
        # ax is actually a grudsoec,
        gs0 = ax.subgridspec(ncols=2, nrows=2, height_ratios=[2, 1], hspace=0., wspace=0.08)
        ax0=fig.add_subplot(gs0[0,0]) # surface
        ax1=fig.add_subplot(gs0[1,0]) # depths
        ax2=fig.add_subplot(gs0[0,1]) # surface
        ax3=fig.add_subplot(gs0[1,1]) # depths

    ax0 = add_letter(ax0, 'c')
    ax2 = add_letter(ax2, 'd')

    models = {}
    model_cubes = {}
    model_cubes_paths = {}
    model_hist_cubes = {}
    debug_txt = ''
    for variable_group, filenames in profile_fns.items():
        for i, fn in enumerate(filenames):
            dataset = metadatas[fn]['dataset']
            short_name = metadatas[fn]['short_name']
            exp = metadatas[fn]['exp']
            ensemble = metadatas[fn]['ensemble']
            if dataset in models_to_skip['all']: continue
            if dataset in models_to_skip.get(short_name, {}): continue
            models[dataset] = True
            if single_model == 'all':
                pass
            elif dataset != single_model:
                continue

            cube = iris.load_cube(fn)
            try: cube.coord('depth')
            except: continue
            short_name = metadatas[fn]['short_name']

            cube = diagtools.bgc_units(cube, metadatas[fn]['short_name'])
            if short_name == 'chl':
                cube = fix_chl(cube)
                if metadatas[fn]['dataset'] == 'CESM2' and single_model != 'CESM2':
                     continue
            scenario = metadatas[fn]['exp']
            debug_txt = ', '.join([debug_txt, 'a', dataset, short_name, exp, ensemble, str(cube[:,9].data.mean()), '\n'])
            print(dataset, short_name, exp, ensemble, cube[:,9].data.mean())
            if scenario == 'historical':
                cube = extract_time(cube, hist_time_range[0], 1, 1, hist_time_range[1], 12, 31)
            else:
                cube = extract_time(cube, ssp_time_range[0], 1, 1, ssp_time_range[1], 12, 31)
            print(dataset, short_name, exp, ensemble, cube[9].data)
            debug_txt = ', '.join([debug_txt, 'b', dataset, short_name, exp, ensemble, str(cube[:,9].data.mean()), '\n'])

            # Only interested in the mean over the time range provided for each model.
            cube = cube.collapsed('time', iris.analysis.MEAN)
            print(dataset, short_name, exp, ensemble, cube[9].data)
            debug_txt = ', '.join([debug_txt, 'c', dataset, short_name, exp, ensemble, str(cube[9].data), '\n'])

            print(cube.shape)
            # ensure everyone uses the same levels:
            cube = extract_levels(cube,
                scheme='linear',
                levels =  [0.5, 1.0, 5.0, 10.0, 50.0, 100.0, 150., 200.0, 250., 300.0, 350., 400.0, 450., 500.0,
                           600.0, 650., 700.0, 750., 800.0, 850., 900.0, 950., 999.0,
                           1001., 1250., 1500.0, 1750., 2000.0, 2250., 2500.0, 2750., 3000.0, 3250., 3500.0, 3750.,
                           4000.0, 4250., 4500.0, 4750., 5000.0]
                )

            print(cube.shape)
            #assert 0
            print(dataset, short_name, exp, ensemble, cube[13])
            debug_txt = ', '.join([debug_txt, 'd', dataset, short_name, exp, ensemble, str(cube[13].data), '\n\n'])

            model_cubes = add_dict_list(model_cubes, variable_group, cube)
            model_cubes_paths = add_dict_list(model_cubes_paths, variable_group, fn)
            if exp == 'historical':
                model_hist_cubes[(dataset, short_name, exp, ensemble)] = cube

    if not len(model_cubes):
        return fig, ax

    # calculate file count:
    header = ['field', ]
    line = [short_name, '(profile)', ]

    data_dict = {}
    for variable_group, cubes in model_cubes.items():
#        data_values = {}
        scenario = ''
        color = ''
        single_model_means = {}
        single_model_diffs = {}




    for i, cube in enumerate(cubes):
        fn = model_cubes_paths[variable_group][i]
        metadata = metadatas[fn]
        dataset = metadatas[fn]['dataset']
        short_name = metadatas[fn]['short_name']
        ensemble = metadatas[fn]['ensemble']

        if dataset in models_to_skip['all']: continue
        if dataset in models_to_skip.get(short_name, {}): continue

        if dataset not in single_model_means.keys():
            single_model_means[dataset] = {}
            single_model_diffs[dataset] = {}

        if single_model == 'all': pass
        elif metadata['dataset'] != single_model: continue

        scenario = metadata['exp']
        color = ipcc_colours[scenario]

        depths = -1.* cube.coord('depth').points

        # Try to fin historical cube
        hist_cube = find_hist_cube(model_hist_cubes, dataset, short_name, scenario, ensemble)

        for z, d in zip(depths, cube.data):
            single_model_means[dataset] = add_dict_list(single_model_means[dataset], z, d)
        if 'all_models'  in plotting:
             ax0, ax1 = plot_z_line(depths, cube.data, ax0, ax1, ls='-', c=color, lw=0.8, label = scenario)

        for z, d, hist in zip(depths, cube.data, hist_cube.data):
            single_model_diffs[dataset] = add_dict_list(single_model_diffs[dataset], z, d - hist)
        if 'all_models'  in plotting:
            ax2, ax3 = plot_z_line(depths, cube.data - hist_cube.data, ax2, ax3, ls='-', c=color, lw=0.8, label = scenario)
        # ./output_csv/output_csv_*/model_table/ensemble_count_table_*.csv
            mmm = {}

        mmm_diff = {}
        header.append('&')
        header.append(scenario)
        num_ens = 0
        # ensemble mean:
        for dataset, data_values in single_model_means.items():

            depths = sorted(data_values.keys())
            model_mean = [np.ma.masked_invalid(data_values[z]).mean() for z in depths]
            print(depths, model_mean, data_values)
            num_ens += len(data_values[depths[0]])
            for z, dep  in enumerate(depths):
                print(variable_group, scenario, dataset, z, dep, len(data_values[dep]), data_values[dep],  model_mean[z], 'mean:', model_mean[z])

            for z,d in zip(depths, model_mean):
                 mmm = add_dict_list(mmm, z, d)
                 if z == 500.:
                     debug_txt = ', '.join([debug_txt,'single_model_means:', dataset, variable_group, 'mean', z, d, '\n'])

                 print('mmm',variable_group, dataset, z, d)

        # anomaly:
        for dataset, diff_values in single_model_diffs.items():
            print(dataset, diff_values)
            ddepths = sorted(diff_values.keys())
            model_anom_mean = [np.ma.masked_invalid(diff_values[z]).mean() for z in ddepths]
 #          print(depths, model_anom_mean, diff_values)
#            for z, dep  in enumerate(depths):
#                print(variable_group, scenario, dataset, z, dep, len(diff_values[dep]), diff_values[dep],  model_anom_mean[z], 'mean:', model_anom_mean[z])

            for z,d in zip(ddepths, model_anom_mean):
                 mmm_diff = add_dict_list(mmm_diff, z, d)
                 if z == 500.:
                     debug_txt = ', '.join([debug_txt,'single_model_means:', dataset, variable_group, 'mean', z, d, '\n'])

                 print('mmm',variable_group, dataset, z, d)

        scen_str = ''.join([str(int(len(single_model_means.keys()))), ' (',str(int(num_ens)), ') &'])
        line.append(scen_str)

        depths = sorted(mmm.keys())
        mean = np.ma.array([np.mean(mmm[z]) for z in depths])
        mean = np.ma.masked_invalid(mean)

        ddepths = sorted(mmm_diff.keys())
        mean_diff = np.ma.array([np.mean(mmm_diff[z]) for z in ddepths])
        mean_diff = np.ma.masked_invalid(mean_diff)
        print(variable_group, scenario, ddepths, mean_diff, mmm_diff)
        for z, dep  in enumerate(depths):
            print(variable_group, scenario, z, dep, len(mmm[dep]), mmm[dep], 'mmm:', mean[z],  mean_diff[z])

        print(depths, mean, mean_diff)
        debug_txt = ', '.join([debug_txt,'\n\nmulti_model_means:', variable_group, scenario, ', '.join([str(d) for d in mean]), '\n'])

        for test_val in [0., None, np.nan, np.inf, False]:
             if test_val in mean:
                 print('Fail:', test_val, 'in mean:', mean)
                 assert 0
             if test_val in depths:
                 print('Fail:', test_val, 'in depths:', depths)
                 assert 0


        data_dict[scenario] = {'mean':mean, 'depths': depths, 'scenario': scenario, 'variable_group': variable_group, 'anomaly':mean_diff, 'ddepths':ddepths}

        if 'means' in plotting:
            ax0, ax1 = plot_z_line(depths, mean, ax0, ax1, ls='-', c=color, lw=2., label = scenario)
            save_profile_csv(cfg, ['profile', variable_group,single_model, model, scenario, 'means'], depths, mean)

        if scenario in ['historical', 'hist']:
            zorder = 5
        else: zorder=1

        if 'means_split' in plotting:
            ax0, ax1 = plot_z_line(depths, mean, ax0, ax1, ls='-', c=color, lw=2., label = scenario, zorder=zorder)
            save_profile_csv(cfg, ['profile', variable_group, single_model, scenario, 'means_split'], depths, mean)

        if 'medians' in plotting:
            assert 0
            # medians = [np.median(data_values[z]) for z in depths]
            # data_dict[scenario]['medians'] = medians
            # ax0, ax1 = plot_z_line(depths, medians, ax0, ax1, ls='-', c=color, lw=2., label = scenario)
            # save_profile_csv(cfg, ['profile', variable_group, single_model, scenario, 'medians'], depths, mean)

        if 'range' in plotting:
            assert 0
            # mins = [np.min(data_values[z]) for z in depths]
            # maxs = [np.max(data_values[z]) for z in depths]
            # data_dict[scenario]['mins'] = mins
            # data_dict[scenario]['maxs'] = maxs
            # ax0, ax1 =  plot_z_area(depths, mins, maxs,ax0, ax1, color= ipcc_colours[scenario], alpha=0.15)
            # save_profile_csv(cfg, ['profile', variable_group, single_model, scenario, 'mins'], depths, mins)
            # save_profile_csv(cfg, ['profile', variable_group, single_model, scenario, 'maxs'], depths, maxs)

        if '5-95' in plotting:
            assert 0
            # mins = [np.percentile(data_values[z], 5.) for z in depths]
            # maxs = [np.percentile(data_values[z], 95.) for z in depths]
            # data_dict[scenario]['5pc'] = mins
            # data_dict[scenario]['95pc'] = maxs
            # ax0, ax1 =  plot_z_area(depths, mins, maxs, ax0, ax1, color= ipcc_colours[scenario], alpha=0.15)

        if '5-95_split' in plotting:
            assert 0
            # mins = [np.percentile(data_values[z], 5.) for z in depths]
            # maxs = [np.percentile(data_values[z], 95.) for z in depths]
            # data_dict[scenario]['5pc'] = mins
            # data_dict[scenario]['95pc'] = maxs
            # if scenario in ['historical', 'hist']:
            #     ax0, ax1 =  plot_z_area(depths, mins, maxs, ax0, ax1, color= ipcc_colours[scenario], alpha=0.15)

    print(debug_txt)


        #./output_csv/output_csv_*/model_table/ensemble_count_table_*.csv
    csv_str = ''.join([''.join(header), '\\\\\n', ''.join(line), '\\\\\n'])
    print(csv_str)
    csv_path  = diagtools.folder(cfg['work_dir']+'/model_table')
    csv_path += 'ensemble_count_table_'+short_name+'_profile.csv'
    mp_fn = open(csv_path, 'w')
    mp_fn.write(csv_str)
    mp_fn.close()
    print('csv_str:',csv_str)
    print('saved to path:', csv_path)

    #assert 0
    # plot rhs
    #try:
    hist_data = np.array(data_dict['historical']['mean'])

    for scenario, ddict in data_dict.items():
        if not len(hist_data):
            continue
        color = ipcc_colours[scenario]
        if 'means_split' in plotting: # and scenario not in ['historical', 'hist']:

            print(ddict['ddepths'],  ddict['anomaly'])
            ax2, ax3 = plot_z_line(
                ddict['ddepths'],
                ddict['anomaly'],
                ax2, ax3, ls='-', c=color, lw=2., label = scenario)

        if '5-95_split' in plotting:
            assert 0
            # mins = np.array(ddict['5pc']) - hist_data
            # maxs = np.array(ddict['95pc']) - hist_data
            # ax2, ax3 =  plot_z_area(ddict['depths'], mins, maxs, ax2, ax3, color=color, alpha=0.15)

    # Add observational data.
    plot_obs = True
    if single_model=='all':
        plot_obs_diff = False
    else:
        plot_obs_diff = True
    obs_filename = 'aux/obs_ncs/'+short_name+'_profile.nc'
    if plot_obs and os.path.exists(obs_filename):
        obs_cube = iris.load_cube(obs_filename)
        obs_depths = -1.* obs_cube.coord('depth').points
        obs_key = 'Observations'
        ax0, ax1 = plot_z_line(obs_depths, obs_cube.data, ax0, ax1, ls='-', c='k', lw=2.5, label = obs_key)

        if plot_obs_diff:
            # interpolated historical data in obs_depths
            hist_interp = np.interp(obs_depths, depths, hist_data) # calculate histo
            ax2, ax3 = plot_z_line(obs_depths, obs_cube.data - hist_interp, ax2, ax3, ls='-', c='k', lw=2.5, label = obs_key)

    # set x axis limits:
    xlims = np.array([ax0.get_xlim(), ax1.get_xlim()])
    ax0.set_xlim([xlims.min(), xlims.max()])
    ax1.set_xlim([xlims.min(), xlims.max()])

    xlims = np.array([ax2.get_xlim(), ax3.get_xlim()])
    ax2.set_xlim([xlims.min(), xlims.max()])
    ax3.set_xlim([xlims.min(), xlims.max()])

    ylims = np.array([ax0.get_ylim(), ax1.get_ylim(), ax2.get_ylim(), ax3.get_ylim()])
    ax0.set_ylim([-1000., ylims.max()])
    ax1.set_ylim([ylims.min(), -1001.])
    ax2.set_ylim([-1000., ylims.max()])
    ax3.set_ylim([ylims.min(), -1001.])

    # hide between pane axes and ticks:
    ax0.spines['bottom'].set_visible(False)
    ax0.xaxis.set_ticks_position('none')
    ax0.xaxis.set_ticklabels([])
    ax1.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.xaxis.set_ticks_position('none')
    ax2.xaxis.set_ticklabels([])
    ax3.spines['top'].set_visible(False)

    # Hide RHS y axis ticks:
    #ax2.yaxis.set_ticks_position('none')
    ax2.yaxis.set_ticklabels([])
    #ax3.yaxis.set_ticks_position('none')
    ax3.yaxis.set_ticklabels([])

    # draw a line between figures
    ax0.axhline(-999., ls='--', lw=1.5, c='black')
    ax1.axhline(-1001., ls='--', lw=1.5, c='black')
    ax2.axhline(-999., ls='--', lw=1.5, c='black')
    ax3.axhline(-1001., ls='--', lw=1.5, c='black')

    if pH_log and short_name in ['ph', 'pH']:
#        #ax0.set_xscale('log')
#        #ax1.set_xscale('log')
#        ax0.semilogx()
#        ax1.semilogx()
        plt.sca(ax0)
        plt.xscale("log")
        plt.sca(ax1)
        plt.xscale("log")

    ax0.set_title('Depth Profile', fontsize=9) #'small')
    ax2.set_title('Anomaly Profile', fontsize= 9) #'small')

    #title_time_str = 'Difference against Historical'
    #title_time_str = ' '.join(['-'.join([str(int(t)) for t in ssp_time_range]), 'against',
    #                     '-'.join([str(int(t)) for t in hist_time_range])])
    #title_time_str = 'Mid-Century Anomaly'
    title_time_str = anomaly_titles[short_name]

    ax1.set_xlabel(short_suptitles[short_name], fontsize='small')
    ax3.set_xlabel(title_time_str, fontsize='small')

    ax0.tick_params(axis='y', labelsize='small')
    ax1.tick_params(axis='y', labelsize='small')

    ax0.set_ylabel('Depth', loc='bottom', fontsize='small')


    if not save:
        return fig, ax
    # Saving files:

    time_str = '_'.join(['-'.join([str(t) for t in hist_time_range]), 'vs',
                         '-'.join([str(t) for t in ssp_time_range])])
    plt.suptitle(' '.join([single_model, time_str]))
    # Determine image filename:
    path = diagtools.folder(cfg['plot_dir']+'/profiles_pair')
    path += '_'.join(['multi_model_profile', time_str])
    path += '_'+ '_'.join(plotting)
    if single_model == 'all': pass
    else: path += '_'+single_model
    path += diagtools.get_image_format(cfg)

    logger.info('Saving plots to %s', path)
    plt.savefig(path)
    plt.close()



#####################################
# Map sections:
def single_map_figure(cfg, cube, variable_group, exp='', model='', ensemble='', time_range=[], region = 'midatlantic'):
    """
    Make a figure of a single cube.
    """

    time_str = '-'.join([str(t) for t in time_range])
    path = diagtools.folder([cfg['plot_dir'],'Maps', variable_group])
    if ensemble in ['', 'AllEnsembles', 'All',]:
        path = diagtools.folder([path, 'Means'])
    else:
        path = diagtools.folder([path, model])
    path += '_'.join(['maps',model, exp, ensemble, region, time_str])
    path += diagtools.get_image_format(cfg)

    if os.path.exists(path):
        return

    fig = plt.figure()
    fig.set_size_inches(10,6)
    gs = matplotlib.gridspec.GridSpec(ncols=1, nrows=1)

    if region == 'global':
        #central_longitude = -14.25 #W #-160.+3.5
        proj = ccrs.Robinson(central_longitude=central_longitude)

    if region in ['midatlantic', 'tightmidatlantic', 'verytightMPA']:
        proj=cartopy.crs.PlateCarree()

    ax0 = fig.add_subplot(gs[0,0], projection=proj)

    qplot = iris.plot.contourf(
        cube,
        #nspaces[sbp_style],
        linewidth=0,
        extend='both' ,#,
        )

    plt.colorbar()

    if region in ['midatlantic', 'tightmidatlantic', 'verytightMPA',]:

        lat_bnd = 20.
        lon_bnd = 30.
        if region == 'tightmidatlantic':
            lat_bnd = 15.
            lon_bnd = 25.
        if region == 'verytightMPA':
            lat_bnd = 7.5
            lon_bnd = 12.5


        #central_longitude = -14.25 #W #-160.+3.5
        #central_latitude = -7.56
        ax0.set_extent([central_longitude-lon_bnd,
                       central_longitude+lon_bnd,
                       central_latitude-lat_bnd,
                       central_latitude+lat_bnd, ])
        # Compute the required radius in projection native coordinates:
        #r_ortho = compute_radius(proj, 3., proj=proj, lat = central_latitude, lon=central_longitude,)
        ax0.add_patch(mpatches.Circle(xy=[central_longitude, central_latitude], radius=2.88, ec='black', fc=(1.,1.,1.,0.),transform=proj, zorder=30))

    try:
        plt.gca().coastlines()
    except AttributeError:
        logger.warning('Not able to add coastlines')

    print(exp, cube.data.min(),cube.data.max())
    title = ' '.join([model, exp, ensemble, str(cube.units) ]) # long_names.get(sbp_style, sbp_style,)])
    plt.title(title)

    if region == 'midatlantic':
        plt.text(0.95, 0.9, exp,  ha='right', va='center', transform=ax0.transAxes, color=ipcc_colours[exp], fontweight='bold')
    if region in ['tightmidatlantic', 'vertightMPA']:
        plt.text(0.95, 1.01, exp,  ha='right', va='bottom', transform=ax0.transAxes, color=ipcc_colours[exp], fontweight='bold')


    logger.info('Saving plots to %s', path)
    plt.savefig(path)
    plt.close()



def prep_cube_map(fn, metadata, time_range, region):
    """
    Load and prepare a netcdf to become a map cube.
    """
    cube = iris.load_cube( fn)

    cube = diagtools.bgc_units(cube, metadata['short_name'])
    if metadata['short_name'] == 'chl':
         cube = fix_chl(cube)

    if 'time' in [c.name for c in cube.coords()]:
        print('prep_cube_map: extract time: ', time_range)

        cube = extract_time(cube, time_range[0], 1, 1, time_range[1], 12, 31)
        cube = cube.collapsed('time', iris.analysis.MEAN)

    print('map plot', region, cube)
    cube = regrid_intersect(cube, region=region)
    if metadata['short_name'] == 'chl':
         cube = fix_chl(cube)
    return cube


def multi_model_map_figure(
        cfg,
        metadatas,
        maps_fns = {},
        figure_style = 'hist_and_ssp',
        hist_time_range = [2000., 2010.],
        ssp_time_range = [2040., 2050.],
        region='midatlantic',
        obs_metadata={},
        obs_filename='',
        colour_scheme = 'IPCC',
        fig = None,
        ax = None,
        save = False,
        single_model = 'all',
    ):
    """
    produce a monthly climatology figure.
    """
    #central_longitude = -14.25 #W #-160.+3.5
    #central_latitude = -7.56
    if region == 'global':
        #central_longitude = -14.25 #W #-160.+3.5
        proj = ccrs.Robinson(central_longitude=central_longitude)

    if region in ['midatlantic', 'tightmidatlantic', 'verytightMPA',]:
        proj=cartopy.crs.PlateCarree()

    if fig is None or ax is None:
        fig = plt.figure()
        fig.set_size_inches(10,6)
        gs = matplotlib.gridspec.GridSpec(ncols=1, nrows=1)
        ax = gs[0, 0]
        save = True
        fig.set_size_inches(11,5)

    # IPCC colour schemes downloaded from:
    # https://github.com/IPCC-WG1/colormaps/tree/master
    if colour_scheme in ['standard', 'IPCC']:
        seq_cmap = 'viridis'
        div_cmap ='BrBG'
    elif colour_scheme in ['temp', 'prec', 'wind', 'chem', 'cryo', 'slev', 'misc1', 'misc2', 'misc3', 'IPCC']:
        if colour_scheme in ['misc1', 'misc2', 'misc3']:
            key = colour_scheme[:4]
            num = colour_scheme[-1]
            rgb_data_seq = np.loadtxt('colormaps/continuous_colormaps_rgb_0-1/'+key+'_seq_'+num+'.txt')
            seq_cmap = mcolors.LinearSegmentedColormap.from_list('colormap', rgb_data_seq)
            rgb_data_div = np.loadtxt('colormaps/continuous_colormaps_rgb_0-1/misc_div.txt')
            div_cmap = mcolors.LinearSegmentedColormap.from_list('colormap', rgb_data_div)

        else:
            key = colour_scheme #'temp' # prec #wind #misc #slev  #chem
            rgb_data_seq = np.loadtxt('colormaps/continuous_colormaps_rgb_0-1/'+key+'_seq.txt')
            seq_cmap = mcolors.LinearSegmentedColormap.from_list('colormap', rgb_data_seq)
            rgb_data_div = np.loadtxt('colormaps/continuous_colormaps_rgb_0-1/'+key+'_div.txt')
            div_cmap = mcolors.LinearSegmentedColormap.from_list('colormap', rgb_data_div)

    elif colour_scheme in ['ph',]:
            # No divergent cmap !
            seq_cmap='viridis_r'
#           rgb_data_seq = np.loadtxt('colormaps/continuous_colormaps_rgb_0-1/temp_seq.txt')
#           seq_cmap = mcolors.LinearSegmentedColormap.from_list('colormap', rgb_data_seq)

            rgb_data_div = np.loadtxt('colormaps/continuous_colormaps_rgb_0-1/slev_seq.txt')
            div_cmap = mcolors.LinearSegmentedColormap.from_list('colormap', rgb_data_div)
    elif colour_scheme in ['temp_1',]:
            # No divergent cmap !
            rgb_data_seq = np.loadtxt('colormaps/continuous_colormaps_rgb_0-1/temp_seq.txt')
            seq_cmap = mcolors.LinearSegmentedColormap.from_list('colormap', rgb_data_seq[::-1])
            rgb_data_seq = np.loadtxt('colormaps/continuous_colormaps_rgb_0-1/misc_seq_3.txt')
            div_cmap = mcolors.LinearSegmentedColormap.from_list('colormap', rgb_data_seq)

            #rgb_data_div = np.loadtxt('colormaps/continuous_colormaps_rgb_0-1/slev_seq.txt')
            #seq_cmap = 'magma' #'inferno'
            #div_cmap = 'afmhot'
            #div_cmap = 'afmhot_r' #Reds' #mcolors.LinearSegmentedColormap.from_list('colormap', rgb_data_div)

    if figure_style=='five_means':
        subplots_nums = {231: 'historical', 232: 'ssp126', 233: 'ssp245', 235: 'ssp370', 236: 'ssp585'}
        subplot_style = {231: 'hist', 232:'mean', 233: 'mean', 235: 'mean', 236: 'mean'}
        cmaps = {231: seq_cmap, 232:seq_cmap, 233: seq_cmap, 234: seq_cmap, 235: seq_cmap, 236: seq_cmap}

    elif figure_style=='hist_and_ssp':
        # This is the final style used in the publication.
        subplots_nums = {231: 'historical', 232: 'ssp126', 233: 'ssp245', 235: 'ssp370', 236: 'ssp585'}
        subplot_style = {231: 'hist', 232:'diff', 233: 'diff', 235: 'diff', 236: 'diff'}
        cmaps = {231: seq_cmap, 232:div_cmap, 233: div_cmap, 234: seq_cmap, 235: div_cmap, 236: div_cmap}

    else:
        assert 0

    gs1 =ax.subgridspec(ncols=3, nrows=2, )
    subplots = {}
    subplots['historical'] = fig.add_subplot(gs1[0,0], projection=proj)

    plot_obs = True
    if plot_obs:
        subplots['obs'] = fig.add_subplot(gs1[1,0], projection=proj)

    subplots['ssp126'] = fig.add_subplot(gs1[0,1], projection=proj)
    subplots['ssp245'] = fig.add_subplot(gs1[1,1], projection=proj)
    subplots['ssp370'] = fig.add_subplot(gs1[0,2], projection=proj)
    subplots['ssp585'] = fig.add_subplot(gs1[1,2], projection=proj)

    subplots['historical'] = add_letter(subplots['historical'], 'e')
    if plot_obs:
        subplots['obs'] = add_letter(subplots['obs'], 'f' )
    subplots['ssp126'] = add_letter(subplots['ssp126'], 'g')
    subplots['ssp245'] = add_letter(subplots['ssp245'], 'h')
    subplots['ssp370'] = add_letter(subplots['ssp370'], 'i')
    subplots['ssp585'] = add_letter(subplots['ssp585'], 'j')

    model_cubes = {}
    model_cubes_paths = {}
    model_cubes_anom = {}
    model_cubes_paths_anom = {}

    variablegroup_model_cubes = {}
    exps = {}
    models = {}
    short_names={}
    hist_fns = {}
    # Run through all input files and generate list of details.
    for variable_group, filenames in maps_fns.items():
        for i, fn in enumerate(filenames):
            dataset = metadatas[fn]['dataset']
            short_name = metadatas[fn]['short_name']
            scenario = metadatas[filenames[0]]['exp']
            exps[scenario] = variable_group
            ensemble = metadatas[filenames[0]]['ensemble']
            if dataset in models_to_skip['all']:
                continue
            if dataset in models_to_skip.get(short_name, {}):
                continue
            if single_model == 'all':
                pass
            elif dataset != single_model:
                continue
            models[dataset] = True
            short_names[metadatas[fn]['short_name']] = True
            hist_fns[(dataset, short_name, scenario,  ensemble)] = fn
    if len(short_names.keys()) == 0:
        if save:
            return
        else:
            return fig, ax

    if len(short_names) > 1:#
        assert 0
    short_name = list(short_names.keys())[0]

    # Load and do basic map manipulation data
    for variable_group, filenames in maps_fns.items():
        print('variable_group:', variable_group, len(filenames))
        work_dir = diagtools.folder([cfg['work_dir'], 'var_group_means',single_model, ])
        vg_path = work_dir+'_'.join([variable_group, single_model, region ])+'.nc'
        if os.path.exists(vg_path):
            print('path exists:', vg_path)
            variablegroup_model_cubes[variable_group] = iris.load_cube(vg_path)
            scenario= metadatas[filenames[0]]['exp']
            short_name =  metadatas[filenames[0]]['short_name']
            print('loaded:', variablegroup_model_cubes.keys(), variable_group)
            continue

        for model_itr in models:
            print('variable_group:', variable_group, model_itr, 'looking for ', single_model)

            if model_itr in models_to_skip['all']:
                 continue
            if single_model == 'all':
                 pass
            elif model_itr != single_model:
                 continue

            work_dir = diagtools.folder([cfg['work_dir'], 'model_group_means','single_model_'+single_model, model_itr])
            model_path = work_dir+'_'.join([variable_group, model_itr, region])+'.nc'
            scenario = metadatas[filenames[0]]['exp']
            short_name =  metadatas[filenames[0]]['short_name']

            model_cubes_paths = add_dict_list(model_cubes_paths, (variable_group, model_itr), filenames[0])
            #assert 0
            if scenario == 'historical':
                time_range = hist_time_range
            else:
                time_range = ssp_time_range

            if os.path.exists(model_path):
                print('path exists:', model_path)
                model_cubes[(variable_group, model_itr)] = iris.load_cube(model_path)
                model_cubes[(variable_group, model_itr)] = regrid_intersect(model_cubes[(variable_group, model_itr)], region)
                continue
            file_count = 0
            for i, fn in enumerate(filenames):
#                if fn.find('CMIP6_UKESM1-0-LL_Omon_ssp370_r9i1p1f2_tos_2049-2050.nc') > -1:
#                    # this file doesn't work for some reason.
#                continue
                model = metadatas[fn]['dataset']
                ensemble = metadatas[fn]['ensemble']
                if model != model_itr: continue
                if model_itr in models_to_skip['all']: continue
                short_name = metadatas[fn]['short_name']
                if model in models_to_skip.get(short_name, {}): continue


#                if single_model == 'all': pass
#                elif model != single_model: continue
                print('loading:', i, variable_group ,scenario, fn)

                cube = prep_cube_map(fn, metadatas[filenames[0]], time_range, region)
                cube = regrid_intersect(cube, region)
                if short_name == 'chl':
                    cube = fix_chl(cube)
                    if cube.data.max()>100000.:
                        print('WHAT!? Thats too much ChlorophYLL!', cube.data.max(), cube.units)
                        #cube.data = cube.data/1000.
                        assert 0
                    cube = fix_chl(cube)
                model_cubes = add_dict_list(model_cubes, (variable_group, model), cube)

                single_map_figure(cfg, cube, variable_group,
                    exp=scenario,
                    model=model,
                    ensemble=ensemble,
                    time_range=time_range,
                    region = region)

            if not model_cubes.get((variable_group, model_itr), False):
                continue
            # make the model mean cube here:
            model_cubes[(variable_group, model_itr)] = diagtools.make_mean_of_cube_list_notime(model_cubes[(variable_group, model_itr)])
            single_map_figure(
                cfg,
                model_cubes[(variable_group, model_itr)],
                variable_group,
                exp=scenario,
                model=model_itr,
                ensemble='AllEnsembles',
                time_range=time_range,
                region = region)
            iris.save(model_cubes[(variable_group, model_itr)], model_path)

            # add model mean cube to dict.
            variablegroup_model_cubes = add_dict_list(
                variablegroup_model_cubes,
                variable_group,
                model_cubes[(variable_group, model_itr)])

        print('making make_mean_of_cube_list_notime:', variable_group)
        print('variablegroup_model_cubes:',variablegroup_model_cubes)
        cube_list_group = variablegroup_model_cubes.get(variable_group, False)
        if not cube_list_group:
            print('didnt find variable_group in variablegroup_model_cubes')
            print('variable_group:', variable_group)
            print('variablegroup_model_cubes', variablegroup_model_cubes)
            print(single_model)
            #assert 0
            continue
        variablegroup_model_cubes[variable_group] = diagtools.make_mean_of_cube_list_notime(variablegroup_model_cubes[variable_group])

        # make the model mean cube here:
        single_map_figure(
            cfg,
            variablegroup_model_cubes[variable_group],
            variable_group,
            exp=scenario,
            model='AllModels',
            ensemble='AllEnsembles',
            time_range=time_range,
            region = region)

        print('saving:', vg_path)
        iris.save(variablegroup_model_cubes[variable_group], vg_path)

    # Calculate anomaly ensemble menmbers.
    variablegroup_model_cubes_anom = {}
    for variable_group, filenames in maps_fns.items():
        print('variable_group:', variable_group, len(filenames))
        work_dir = diagtools.folder([cfg['work_dir'], 'var_group_anomalies', single_model, ])
        vg_path = work_dir+'_'.join([variable_group, single_model, region ])+'_anom.nc'
        if os.path.exists(vg_path):
            print('path exists:', vg_path)
            variablegroup_model_cubes_anom[variable_group] = iris.load_cube(vg_path)
            scenario= metadatas[filenames[0]]['exp']
            short_name =  metadatas[filenames[0]]['short_name']
            print('loaded:', variablegroup_model_cubes_anom.keys(), variable_group)
            continue

        for model_itr in models:
            print('variable_group:', variable_group, model_itr, 'looking for ', single_model, 'anomaly')
            if model_itr in models_to_skip['all']:
                 continue
            if single_model == 'all':
                 pass
            elif model_itr != single_model:
                 continue

            work_dir = diagtools.folder([cfg['work_dir'], 'model_group_anoms','single_model_'+single_model, model_itr])
            model_path = work_dir+'_'.join([variable_group, model_itr, region])+'_anom.nc'
            scenario = metadatas[filenames[0]]['exp']
            short_name =  metadatas[filenames[0]]['short_name']
            if scenario == 'historical':
                continue

            if scenario == 'historical':
                time_range = hist_time_range
            else:
                time_range = ssp_time_range

            if os.path.exists(model_path):
                print('path exists:', model_path)
                model_cubes_anom[(variable_group, model_itr)] = iris.load_cube(model_path)
                model_cubes_anom[(variable_group, model_itr)] = regrid_intersect(model_cubes_anom[(variable_group, model_itr)], region)
                continue
            file_count = 0
            for i, fn in enumerate(filenames):
                model = metadatas[fn]['dataset']
                short_name = metadatas[fn]['short_name']
                ensemble = metadatas[fn]['ensemble']
                model_cubes_paths_anom = add_dict_list(model_cubes_paths_anom, (variable_group, model_itr), fn)

                if metadatas[fn]['exp'] != scenario:
                    assert 0
                if model != model_itr:
                    continue
                if model_itr in models_to_skip['all']:
                    continue
                if model in models_to_skip.get(short_name, {}):
                    continue

                print('loading:', i, variable_group ,scenario, fn)

                cube = prep_cube_map(fn, metadatas[filenames[0]], time_range, region)
                cube = regrid_intersect(cube, region)

                hist_fn = find_hist_cube(hist_fns, model, short_name, scenario, ensemble)
                hist_cube = prep_cube_map(hist_fn, metadatas[fn], hist_time_range, region)
                hist_cube = regrid_intersect(hist_cube, region)

                if short_name == 'chl':
                    cube = fix_chl(cube)
                    hist_cube = fix_chl(hist_cube)
                    if cube.data.max()>100000. or hist_cube.data.max()>100000.:
                        print('WHAT!? Thats too much ChlorophYLL!', cube.data.max(), cube.units, hist_cube.data.max(),
                              hist_cube.units, )
                        #cube.data = cube.data/1000.
                        assert 0
                    cube = fix_chl(cube)
                    hist_cube = fix_chl(hist_cube)
                cube = cube - hist_cube
                if cube.data.min() == cube.data.max():
                    print('Problem!:', variable_group, model_itr,cube.data.min() , scenario,model_cubes_paths_anom[(variable_group, model_itr)])
                    assert 0

                model_cubes_anom = add_dict_list(model_cubes_anom, (variable_group, model), cube)
                # single_map_figure(cfg, cube, variable_group,
                #     exp=scenario,
                #     model=model,
                #     ensemble=ensemble,
                #     time_range=time_range,
                #     region = region)

            if not model_cubes_anom.get((variable_group, model_itr), False):
                continue
            # make the model mean cube here:
            model_cubes_anom[(variable_group, model_itr)] = diagtools.make_mean_of_cube_list_notime(model_cubes_anom[(variable_group, model_itr)])
            #single_map_figure(
            #    cfg,
            #    model_cubes_anom[(variable_group, model_itr)],
            #    variable_group,
            #    exp=scenario,
            #    model=model_itr,
            #    ensemble='AllEnsembles',
            #    time_range=time_range,
            #    region = region)
            iris.save(model_cubes_anom[(variable_group, model_itr)], model_path)

            # add model mean cube to dict.
            variablegroup_model_cubes_anom = add_dict_list(
                variablegroup_model_cubes_anom,
                variable_group,
                model_cubes_anom[(variable_group, model_itr)])

        print('making make_mean_of_cube_list_notime:', variable_group)
        print('variablegroup_model_cubes_anom:',variablegroup_model_cubes_anom)
        cube_list_group = variablegroup_model_cubes_anom.get(variable_group, False)
        if not cube_list_group:
            print('didnt find variable_group in variablegroup_model_cubes_anom')
            print('variable_group:', variable_group)
            print('variablegroup_model_cubes_anom', variablegroup_model_cubes_anom)
            print(single_model)
            #assert 0
            continue
        variablegroup_model_cubes_anom[variable_group] = diagtools.make_mean_of_cube_list_notime(variablegroup_model_cubes_anom[variable_group])

        # make the model mean cube here:
        # single_map_figure(
        #     cfg,
        #     variablegroup_model_cubes_anom[variable_group],
        #     variable_group,
        #     exp=scenario,
        #     model='AllModels',
        #     ensemble='AllEnsembles',
        #     time_range=time_range,
        #     region = region)

        print('saving:', vg_path)
        iris.save(variablegroup_model_cubes_anom[variable_group], vg_path)


    # calculate diffs, and range.
    nspaces = {}
#    try:
    hist_variable_group = exps['historical']
    hist_cube = variablegroup_model_cubes[hist_variable_group]
    hist_cube_exists = True
#    except:
#        hist_cube_exists = False
#        hist_variable_group = ''

    diff_cubes = {}

    # Calculate the diff range.
    style_range = {'hist':[], 'mean':[], 'diff':[], }
    obs_filename = 'aux/obs_ncs/'+short_name+'_map.nc'
    if plot_obs and os.path.exists(obs_filename):
        obs_cube = iris.load_cube(obs_filename)
        print('obs map:', obs_cube.data.min(), obs_cube.data.max())
        obs_cube = regrid_intersect(obs_cube, region)
        print('obs map:', obs_cube.data.min(), obs_cube.data.max())
        #method = 'min_max' # hist and obs go between min and max.
        #method = '5pc-95pc' # hist and obs plotted between 5-95 percentiles.
        method = '1pc-99pc' # hist and obs plotted between 5-95 percentiles.
        #print('hist:', hist_cube.data.min(), '->', hist_cube.data.max())
        if method == 'min_max':
            if hist_cube_exists:
                ranges = [hist_cube.data.min(), hist_cube.data.max()]
                ranges.extend([obs_cube.data.min(), obs_cube.data.max()])
            else:
                ranges = [obs_cube.data.min(), obs_cube.data.max()]
        if hist_cube_exists:
            cube_range_list = [ hist_cube.data, obs_cube.data]
        else:
            cube_range_list = [obs_cube.data, ]

        if method == '5pc-95pc':
            ranges = []
            for dat, pc in itertools.product(cube_range_list, [5, 95]):
                ranges.append(np.percentile(dat.compressed(), pc))

        if method == '1pc-99pc':
            ranges = []
            for dat, pc in itertools.product(cube_range_list, [1, 99]):
                ranges.append(np.percentile(dat.compressed(), pc))

        print('hist plot range', method, ranges)

        style_range['hist'].extend([np.min(ranges), np.max(ranges)])

    else:
        style_range['hist'].extend([hist_cube.data.min(), hist_cube.data.max()])
    style_range['historical'] =  style_range['hist']

    # Calculate the diff cubes. (Anomaly is now calculated elsewhere.
    for variable_group, cube in variablegroup_model_cubes_anom.items():
        if variable_group == hist_variable_group:
             continue
        # if hist_cube_exists:
        #     cube = cube - hist_cube
        # else: pass # cube
        print(variable_group, cube.data)
        diff_cubes[variable_group] = cube
        style_range['diff'].extend([cube.data.min(), cube.data.max()])

    extends = {}
    # Create the lin space for maps.
    for style, srange in style_range.items():
        if not len(srange): continue
        print(style, srange)
        style_range[style] = [round_sig(np.array(srange).min(),2), round_sig(np.array(srange).max(),2)]
        #nbins = (style_range[style][1] - style_range[style][0])
        #print(style_range[style], nbins)
        #assert 0


        # Symetric around zero:
        if style in ['diff', ]:
            new_min  = np.min(style_range[style])
            new_max  = np.max(style_range[style])
            abs_max = np.abs(style_range[style]).max()
            if new_min * new_max < 0.: # ie either side of zero!
                nspaces[style] = np.linspace(-abs_max, abs_max, 21)
            else:
                nspaces[style] = np.linspace(new_min, new_max, 21)
            #elif new_min>0. and new_max >0.:
            #    nspaces[style] = np.linspace(0., abs_max, 21)
            #elif new_min<0. and new_max <0.:
            #    nspaces[style] = np.linspace(-1.* abs_max, 0., 21)
            #else: assert 0

        else:
            nspaces[style] = np.linspace(style_range[style][0], style_range[style][1], 11)

    shared_cmap = {'hist':[], 'ssp':[]}
    shaped_ims = {'hist':[], 'ssp':[]}
    for sbp, exp in subplots_nums.items():
        ax0 = subplots[exp]
        plt.sca(ax0)
        sbp_style = subplot_style[sbp]
        if figure_style=='hist_and_ssp':
            print(sbp, exp, sbp_style)

            if exp in ['historical', 'hist']: #d hist_cube_exists:
                cube = hist_cube
            else:
                variable_group = exps[exp]
                if variable_group not in diff_cubes: continue
                cube = diff_cubes[variable_group]
        else:
            assert 0
        print('plotting', exp, sbp, single_model)
        qplot = iris.plot.contourf(
            cube,
            nspaces[sbp_style],
            linewidth=0,
            cmap=cmaps[sbp],
            extend='both', # was both
            zmin=style_range[sbp_style][0],
            zmax=style_range[sbp_style][1],
            )
        if sbp_style == 'hist':
            shared_cmap['hist'].append(ax0)
            shaped_ims['hist'].append(qplot)

        if sbp_style == 'diff':
            shared_cmap['ssp'].append(ax0)
            shaped_ims['ssp'].append(qplot)

#        if region in ['midatlantic', 'tightmidatlantic']:
#            lat_bnd = 20.
#            lon_bnd = 30.
#            if region == 'tightmidatlantic':
#                lat_bnd = 15.
#                lon_bnd = 25.
#            ax0.set_extent([central_longitude-lon_bnd,
#                           central_longitude+lon_bnd,
#                           central_latitude-lat_bnd,
#                           central_latitude+lat_bnd, ])

        # Compute the required radius in projection native coordinates:
        #r_ortho = compute_radius(proj, 3., proj=proj, lat = central_latitude, lon=central_longitude,)
        #ax0.add_patch(mpatches.Circle(xy=[central_longitude, central_latitude], radius=2.88, color='black', alpha=0.3, transform=proj, zorder=30))

        # Add circle:
        rad = 2.88
        ax0.add_patch(mpatches.Circle(xy=[central_longitude, central_latitude], radius=rad, ec='black', fc=(1.,1.,1.,0.), transform=proj, zorder=30))
        # add square:
        square_lon_cc = [central_longitude-rad, central_longitude+rad, central_longitude+rad, central_longitude-rad, central_longitude-rad]
        square_lat_cc = [central_latitude-rad, central_latitude-rad, central_latitude+rad, central_latitude+rad, central_latitude-rad ]
        ax0.plot(square_lon_cc, square_lat_cc,
             #[-17.25, -11.25, -11.25, -17.25, -17.25 ], [-10.56, -10.56, -4.56, -4.56, -10.56],
             color='black', linewidth=1,
             transform=ccrs.Geodetic(),
             )

        gl = ax0.gridlines(draw_labels=False, linewidth=0.5, alpha=0.4, color='k',linestyle='--')

        try:
            plt.gca().add_feature(cartopy.feature.LAND, zorder=2, facecolor='#D3D3D3')
            plt.gca().coastlines(zorder=3)

        except AttributeError:
            logger.warning('Not able to add coastlines')

        # Add title to plot
        #long_names = {
        #   'diff':'difference',
        #   'hist':'mean',
        #   'historial':'mean',
        #}
        print(exp, cube.data.min(),cube.data.max())

        title = ' '.join([sspify(exp),]) # long_names.get(sbp_style, sbp_style,)])
        if single_model !='all':
            title = ' '.join([title, single_model])
        if region in ['midatlantic', ]:
            plt.text(0.95, 0.9, title,  ha='right', va='center',
                     transform=ax0.transAxes,
                     color=ipcc_colours[exp],fontweight='bold')
            #plt.text(0.99, 0.94, title,  ha='right', va='center', transform=ax0.transAxes,color=ipcc_colours[exp],fontweight='bold')

        elif region in ['tightmidatlantic', 'verytightMPA',]:
            plt.text(0.95, 1.01, title,  ha='right', va='bottom',
                     transform=ax0.transAxes,
                     color=ipcc_colours[exp],fontweight='bold')

        else:
            plt.title(title)

    obs_filename = 'aux/obs_ncs/'+short_name+'_map.nc'
    print('maps: obs file:', obs_filename)
    if plot_obs and os.path.exists(obs_filename):
        obs_cube = iris.load_cube(obs_filename)
        if short_name == 'chl':
            obs_cube.var_name = 'chl'
            obs_cube.data = obs_cube.data*1000.
#       obs_cube = diagtools.bgc_units(obs_cube, short_name)

        # obs_cube = obs_cube.collapsed('time', iris.analysis.MEAN)
        obs_cube = regrid_intersect(obs_cube, region=region)

        obs_key = 'Observations' #obs_metadata['dataset']
        ax0 = subplots['obs']
        plt.sca(ax0)
        sbp_style = 'hist'

        qplot = iris.plot.contourf(
            obs_cube,
            nspaces[sbp_style],
            linewidth=0,
            cmap=cmaps[234],
            extend='both' , # was both
            zmin=style_range[sbp_style][0],
            zmax=style_range[sbp_style][1],
            )
#        if region in ['midatlantic', 'tightmidatlantic']:
#            lat_bnd = 20.
#            lon_bnd = 30.
#            if region == 'tightmidatlantic':
#                lat_bnd = 15.
#                lon_bnd = 25.
#
#            ax0.set_extent([central_longitude-lon_bnd,
#                           central_longitude+lon_bnd,
#                           central_latitude-lat_bnd,
#                           central_latitude+lat_bnd, ])

        print('obs:', obs_cube.data.min(),obs_cube.data.max())
        #plt.gca().coastlines()
        plt.gca().add_feature(cartopy.feature.LAND, zorder=2, facecolor='#D3D3D3')
        plt.gca().coastlines(zorder=3)

        # Compute the required radius in projection native coordinates:
        #r_ortho = compute_radius(proj, 3., proj=proj, lat = central_latitude, lon=central_longitude,)
        rad=2.88
        ax0.add_patch(mpatches.Circle(xy=[central_longitude, central_latitude],
             radius=rad, ec='black', fc=(1.,1.,1.,0.), transform=proj, zorder=30))

        #tart_longitude= -17.25 #14.25W
        #nd_longitude = -11.25 #
        #tart_latitude= -10.56 # -7.56
        #nd_latitude= -4.56 #
        square_lon_cc = [central_longitude-rad, central_longitude+rad, central_longitude+rad, central_longitude-rad, central_longitude-rad]
        square_lat_cc = [central_latitude-rad, central_latitude-rad, central_latitude+rad, central_latitude+rad, central_latitude-rad ]

        ax0.plot(
             square_lon_cc, square_lat_cc,
             #[17.25, -11.25, -11.25, -17.25, -17.25 ], [-10.56, -10.56, -4.56, -4.56, -10.56],
             color='black', linewidth=1,
             transform=ccrs.Geodetic(),
             )

        gl = ax0.gridlines(draw_labels=False, linewidth=0.5, alpha=0.4, color='k',linestyle='--')

        shared_cmap['hist'].append(ax0)
        shaped_ims['hist'].append(qplot)

        title = ' '.join(['Observations',]) # long_names.get(sbp_style, sbp_style,)])
        if region in ['midatlantic', ]:
            plt.text(0.99, 0.92, title,  ha='right', va='center', transform=ax0.transAxes,color='black',fontweight='bold',fontsize='small')
        elif region in ['tightmidatlantic', 'verytightMPA']:
            plt.text(0.95, 1.01, title,  ha='right', va='bottom', transform=ax0.transAxes,color='black',fontweight='bold',) #fontsize='small')
        else:
            plt.title(title)
    else:
        # no hist plot.
        ax0 = subplots['obs']
        plt.sca(ax0)
        plt.axis('off')

    if len(shaped_ims['hist']):
        fig.colorbar(shaped_ims['hist'][0], ax=shared_cmap['hist'], label=short_suptitles[short_name], ) #, label = 'Historical')

    if len(shaped_ims['ssp']):
        title_time_str = anomaly_titles[short_name]

        fig.colorbar(shaped_ims['ssp'][0], ax=shared_cmap['ssp'], label=title_time_str)

    # save and close.
    time_str = '_'.join(['-'.join([str(t) for t in hist_time_range]), 'vs',
                         '-'.join([str(t) for t in ssp_time_range])])

    if save:
        path = diagtools.folder(cfg['plot_dir']+'/Maps_6panes')
        path += '_'.join(['maps', figure_style, region, time_str, single_model, colour_scheme])
        path += diagtools.get_image_format(cfg)
        logger.info('Saving plots to %s', path)
        plt.savefig(path)
        plt.close()
    else:
        return fig, ax



#####################################
# old Map sections:

def regrid_intersect(cube, region='global'):
    """
    Regird the map
    """
    #central_longitude = -14.25 #W #-160.+3.5
    #central_latitude = -7.56
    cube = regrid(cube, '1x1', 'linear')
    #cube = regrid_to_1x1(cube)
    if region=='global':
        cube = cube.intersection(longitude=(central_longitude-180., central_longitude+180.))
        return cube
    if region=='midatlantic':
        lat_bnd = 20.
        lon_bnd = 30.
    if region=='tightmidatlantic':
        lat_bnd = 15.
        lon_bnd = 25.
    if region=='verytightMPA':
        lat_bnd = 7.5
        lon_bnd = 12.5
    cube = cube.intersection(longitude=(central_longitude-lon_bnd, central_longitude+lon_bnd),
                             latitude=(central_latitude-lat_bnd, central_latitude+lat_bnd), )


    return cube


def compute_radius(ortho, radius_degrees, proj= ccrs.PlateCarree(), lat=0, lon=0):
    """
    catlculate the correct radius:
    from:
    https://stackoverflow.com/questions/52105543/drawing-circles-with-cartopy-in-orthographic-projection
    """
    phi1 = lat + radius_degrees if lat <= 0 else lat - radius_degrees
    _, y1 = ortho.transform_point(lon, phi1, proj)
    return abs(y1)



#####################################
def add_legend(ax):
    """
    Add a legend in a subplot.
    """

    #rows = 25
    #ncols = int(legend_size / nrows) + 1
    #ax1.set_position([
    #    box.x0, box.y0, box.width * (1. - column_width * ncols), box.height
    #])
    # Add emply plots to dummy axis.
    plt.sca(ax)

    for exp in sorted(ipcc_colours):
        if exp=='hist': continue
        plt.plot([], [], c=ipcc_colours[exp], lw=2.5, ls='-', label=sspify(exp))

    plt.plot([], [], c='k', lw=2.5, ls='-', label='Observations')

    plt.fill_betweenx([], [], [], color= 'k', alpha=0.3, label = '2000-2010')
    plt.fill_betweenx([], [], [], color= 'purple', alpha=0.3, label = '2040-2050')

    # plt.plot([], [], c='k', alpha = 0.25, lw=7.5, ls='-', label='5-95 pc')

    legd = ax.legend(
        loc='center left',
        ncol=1,
        prop={'size': 10},
        bbox_to_anchor=(-2., 0.5))
    legd.draw_frame(False)
    legd.get_frame().set_alpha(0.)
    plt.axis('off')
    plt.xticks([])
    plt.yticks([])

    return ax


def do_gridspec(cfg, ):
    """
    Create a grid spec for the master plot.
    """
    #fig = None
    fig = plt.figure()
    fig.set_size_inches(12., 9.)

    gs = matplotlib.gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[5., 5.]) #hspace=0.5,
    gs0 =gs[0,0].subgridspec(ncols=5, nrows=2,
        width_ratios=[2, 2, 1., 1., 0.25],
        height_ratios=[3, 4.6],
        hspace=0.5, wspace=0.55)

    subplots = {}
    subplots['timeseries'] = fig.add_subplot(gs0[0:2,0:2])
    subplots['climatology'] = fig.add_subplot(gs0[0, 2:4])
    #subplots['clim_diff'] = fig.add_subplot(gs0[0, 3])
    subplots['profile'] = gs0[1, 2:4]
    add_profile_label = False
    if add_profile_label:
        fig.add_subplot(subplots['profile'], frameon=False)
        # hide tick and tick label of the big axes
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        plt.grid(False)
        plt.ylabel('Profile')

    subplots['legend'] = fig.add_subplot(gs0[:, -1])
    subplots['legend'] = add_legend(subplots['legend'])
    #subplots['prof_diff'] = fig.add_subplot(gs0[1, 3])
    #
    # maps:
    subplots['maps'] = gs[1,0]
    # gs1 =gs[1,0].subgridspec(ncols=3, nrows=2, )
    # subplots['map_hist'] = fig.add_subplot(gs1[0,0])
    # subplots['map_obs'] = fig.add_subplot(gs1[1,0])
    # subplots['map_ssp125'] = fig.add_subplot(gs1[0,1])
    # subplots['map_ssp245'] = fig.add_subplot(gs1[1,1])
    # subplots['map_ssp379'] = fig.add_subplot(gs1[0,2])
    # subplots['map_ssp585'] = fig.add_subplot(gs1[1,2])
    #plt.show()
    #plt.savefig(cfg['plot_dir']+'/tmp.png')
    #assert 0
    return fig, subplots


def add_dict_list(dic, key, item):
   try: dic[key].append(item)
   except: dic[key]= [item, ]
   return dic


def add_letter(ax, letter, par=True ):
    """
    Adds  Letter to top left above the pane.
    """
    if par:#parenthesis.
        letter = ''.join([letter, ')'])
    ax.text(0., 1.0, letter,
       horizontalalignment='left',
       verticalalignment='bottom',
       transform=ax.transAxes)
    return ax


def main(cfg):
    """
    Load the config file and some metadata, then pass them the plot making
    tools.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

    """
    # Here is plan the parts of the plot:
    # Plot A:
    # top left: time series
    # top right: climatology
    # middle: profile:
    # right: maps

    # Strategy is to get all the individual plots made.
    # individual plots of the ensemble mean
    # Then stitch them together into a big single plot.
    metadatas = diagtools.get_input_files(cfg, )

    time_series_fns = {}
    profile_fns = {}
    maps_fns = {}
    models = {'all': True}
    model_scenarios = {}
    variable_groups = {}

    # This should be the only time that it iterates over metadata.items
    for fn, metadata in metadatas.items():
        variable_group = metadata['variable_group']
        dataset = metadata['dataset']
        short_name = metadata['short_name']
        if dataset in models_to_skip['all']: continue
        if dataset in models_to_skip.get(short_name, {}): continue
        if dataset in models_to_skip:
            continue

        scenario = metadata['exp']
        if not model_scenarios.get(dataset, False):
            model_scenarios[dataset] = {scenario: True}
        else:
            model_scenarios[dataset][scenario] = True



    models_without_futures = []
    for model, scenarios in model_scenarios.items():
        if len(scenarios.keys()) < 2:
            print(model, 'does not have enough scenarios:', scenarios)
            models_without_futures.append(model)
        else:
            print(model, 'is fine: ', scenarios)


    model_hist_cubes= {'ts':{}, 'profile':{}, 'map':{}}

    for fn, metadata in metadatas.items():
        #print(os.path.basename(fn),':',metadata['variable_group'])
        variable_group = metadata['variable_group']

        dataset = metadata['dataset']
        short_name = metadata['short_name']
        exp = metadata['exp']
        ensemble = metadata['ensemble']
        if dataset in models_to_skip['all']: continue
        if dataset in models_to_skip.get(short_name, {}): continue

        if dataset in models_to_skip:
            continue

        if dataset in models_without_futures:
            continue

        models[dataset] = True
        variable_groups[variable_group] = True
        if variable_group.find('_ts_')>-1:
            time_series_fns = add_dict_list(time_series_fns, variable_group, fn)
            model_hist_cubes['ts'][(dataset, short_name, exp, ensemble)] = fn
        if variable_group.find('_profile_')>-1:
            profile_fns = add_dict_list(profile_fns, variable_group, fn)
            model_hist_cubes['profile'][(dataset, short_name, exp, ensemble)] = fn

        if variable_group.find('_map_')>-1:
            model_hist_cubes['map'][(dataset, short_name, exp, ensemble)] = fn
            maps_fns = add_dict_list(maps_fns, variable_group, fn)


    # Find missing historical datasets:
    missing_cubes = {'ts':{}, 'profile':{}, 'map':{}}
    missing = 0
    for files_dict, models_hists in zip([time_series_fns, profile_fns, maps_fns], 
        [model_hist_cubes['ts'], model_hist_cubes['profile'], model_hist_cubes['map']]):
        for outer_vg, files in files_dict.items():
            for fn in files:
                dataset = metadatas[fn]['dataset']
                short_name = metadatas[fn]['short_name']
                ensemble = metadatas[fn]['ensemble']
                scenario = metadatas[fn]['exp']
                variable_group = metadatas[fn]['variable_group']
                if variable_group.find('_ts_')>-1:
                    vg_type = 'ts'
                if variable_group.find('_profile_')>-1:
                    vg_type = 'profile'
                if variable_group.find('_map_')>-1:
                    vg_type = 'map'

                if scenario == 'historical':
                    continue
                hist_cube = find_hist_cube(models_hists, dataset, short_name, scenario, ensemble, return_missing=True)
                if isinstance(hist_cube, tuple): 
                    missing_cubes[vg_type][hist_cube] = True

    for vg_type, missing_tuples in missing_cubes.items():
        if len(missing_tuples):
            print('The following historical jobs are missing:', vg_type)
            for miss_tuple in missing_tuples.keys():
                print(vg_type, 'dataset: ',miss_tuple[0], 'exp: historical, ensemble:',miss_tuple[3])
                missing+=1
    if missing >0:
        assert 0
    #jobs:
    #Make UKESM-only version of everything (or single model?)
    #prepare data for export.
    #fix chlorophyll map.
    #add text about Being right for the wrong reason
    #maybe make some bias correction methods?
    #    subtract the time series by the mean of the observatioinal; time range


    #for ddict in [time_series_fns, profile_fns, maps_fns]:
    #    for kkey, ffiles in ddict.items():
    #        print(kkey, ffiles)
    #        for fn in ffiles:
    #            print('-----\n', fn)
    #            cube = iris.load_cube(fn)
    #            print(cube.data.min(), cube.data.max())
    #assert 0


    # Individual plots - standalone
    do_standalone = True
    if do_standalone:

        # plottings:
        #     all: every single ensemble member shown as thin line.
        #     model_means: the mean of every model ensmble is shown as a thicker dotted line
        #     model_median: the median of every model ensmble is shown as a thicker dotted line
        #     model_range: the range of every model ensmble is shown as a coloured band (or line if only one)
        #     model_5_95: the range of every model ensmble is shown as a coloured band (or line if only one)
        #     Global_mean: the mean of every model_means is shown as a thick solid line
        #     Global_range: the mean of every model_means is shown as a thick solid line
        plottings = [['all', 'model_means', 'Global_mean' ],
                     ['model_range',],
                     ['Global_range', ],
                     ['model_range', 'Global_mean',],
                     ['Global_mean', 'Global_range'],
                     ['model_means', 'Global_range'],
                     ] #'medians', 'all_models', 'range',
        for plotting in plottings:
            continue
            for single_model in models.keys():
                continue
                multi_model_time_series(
                    cfg,
                    metadatas,
                    ts_dict = time_series_fns,
                    moving_average_str='annual',
                    hist_time_range = [2000., 2010.],
                    ssp_time_range = [2040., 2050.],
                    plotting = plotting,
                    single_model=single_model,
                )

            for single_model in ['all', ]: #'not', 'only', 'all']:
                multi_model_time_series(
                    cfg,
                    metadatas,
                    ts_dict = time_series_fns,
                    moving_average_str='annual',
                    hist_time_range = [2000., 2010.],
                    ssp_time_range = [2040., 2050.],
                    plotting = plotting,
                    single_model=single_model,
                )

        # Climatology plot
        plottings =  [['OneModelOneVote',],['means', 'OneModelOneVote',], ['OMOC_modelmeans','OneModelOneVote','OMOC_modelranges'],] # ['all_models',],[ 'means',  '5-95'], ]#  ['means',],  ['5-95',], ['all_models', ]]
        for plotting in plottings:
            continue
            multi_model_clim_figure(
                cfg,
                metadatas,
                time_series_fns,
                hist_time_range = [2000., 2010.],
                ssp_time_range = [2040., 2050.],
                plotting=plotting,
                single_model='all',
            )
        # maps:
        for single_model in ['all', ]: #sorted(models.keys()): # ['all', 'only']:
            #continue
            for reg in ['tightmidatlantic', ]: #'midatlantic', 'verytightMPA']:
              for cmap in ['slev', ]:#'temp', 'prec', 'wind', 'chem', 'cryo', 'standard', 'slev', 'misc1', 'misc2', 'misc3']:
                multi_model_map_figure(
                    cfg,
                    metadatas,
                    maps_fns = maps_fns,
                    figure_style = 'hist_and_ssp',
                    hist_time_range = [2000., 2010.],
                    ssp_time_range = [2040., 2050.],
                    colour_scheme = cmap,
                    region=reg,
                    single_model=single_model )
            #assert 0

        # Profile pair
        plottings =  [['all_models',], ['means_split',], ]#['5-95_split',], ['means_split', '5-95_split', ], ['range',] ]
        for plotting in plottings:
            #continue
            for single_model in sorted(models.keys()): # ['all', 'only']:
                make_multi_model_profiles_plotpair(
                    cfg,
                    metadatas,
                    profile_fns = profile_fns,
                    #short_name,
                    #obs_metadata={},
                    #obs_filename='',
                    hist_time_range = [2000., 2010.],
                    ssp_time_range = [2040., 2050.],
                    plotting = plotting,
                    single_model = single_model,
                    #figure_style = 'difference',
                    fig = None,
                    ax = None,
                    save = False,
                    draw_legend=True
                )

    #print(time_series_fns)
    #print(profile_fns)
    #print(maps_fns)

    do_whole_plot = True
    for single_model in ['all', ]: # sorted(models.keys()): # ['all', 'only']:
    #for single_model in sorted(models.keys()): # ['all', 'only']:

        if not  do_whole_plot: continue
        if single_model in models_to_skip['all']: continue
        if single_model in models_without_futures: continue
        #if single_model != 'all':continue

        fig, subplots = do_gridspec(cfg, )
        hist_time_range = [2000., 2010.] #[1990., 2015.]
        ssp_time_range = [2040., 2050.]

        # a
        fig, subplots['timeseries'] = multi_model_time_series(
            cfg,
            metadatas,
            ts_dict = time_series_fns,
            moving_average_str='annual',
            #colour_scheme = 'viridis',
            hist_time_range = hist_time_range,
            ssp_time_range = ssp_time_range,
            plotting= ['model_range', 'Global_mean',], #['Global_mean', 'Global_range'],#['means', '5-95',],
            fig = fig,
            ax = subplots['timeseries'],
            single_model = single_model
            )
        subplots['timeseries'] = add_letter(subplots['timeseries'], 'a')

        # b
        fig, subplots['climatology'] =  multi_model_clim_figure(
            cfg,
            metadatas,
            time_series_fns,
            hist_time_range = hist_time_range,
            ssp_time_range = ssp_time_range,
            fig = fig,
            ax =  subplots['climatology'],
            plotting=['OneModelOneVote', ], #'means',],
            single_model = single_model,
        )
        subplots['climatology'] = add_letter(subplots['climatology'], 'b')

        # c, d
        fig, subplots['profile'] = make_multi_model_profiles_plotpair(
            cfg,
            metadatas,
            profile_fns = profile_fns,
            #short_name,
            #obs_metadata={},
            #obs_filename='',
            hist_time_range = hist_time_range,
            ssp_time_range = ssp_time_range,
            plotting=['means_split',],
            #figure_style = 'difference',
            fig = fig,
            ax = subplots['profile'],
            single_model = single_model,
            )

        if pH_log and short_name in ['ph', 'pH']:
            #subplots['timeseries'].set_yscale('log')
            #subplots['climatology'].set_yscale('log')
            #subplots['timeseries'].semilogy()
            #subplots['climatology'].semilogy()
            plt.sca(subplots['timeseries'])
            plt.yscale("log")
            plt.sca(subplots['climatology'])
            plt.yscale("log")

            #subplots['profile'].set_xscale('log')

        # cmaps:
        # 'slev',
        # 'temp', 'prec', 'wind', 'chem', 'cryo', 'standard', 'slev', 'misc1', 'misc2', 'misc3']:
        # e, f, g, h, i, j
        cmaps = {
            'tos': 'temp_1', 'thetao': 'temp_1', # temperature
            'sos': 'temp', 'so': 'temp',
            'mld': 'slev', 'mlotst': 'slev',
            'ph': 'ph', # only one way!
            'o2' : 'chem',
            'chl':'prec',
            'intpp': 'prec', #'misc1',
            'no3': 'chem',
            'po4': 'chem',
        }
        fig, subplots['maps'] = multi_model_map_figure(
            cfg,
            metadatas,
            maps_fns = maps_fns,
            figure_style = 'hist_and_ssp',
            hist_time_range = hist_time_range,
            ssp_time_range = ssp_time_range,
            region='tightmidatlantic', #midatlantic',
            fig = fig,
            colour_scheme = cmaps.get(short_name, 'temp'),
            ax =  subplots['maps'],
            single_model = single_model,
            )

        joining = '\n'
        if 'tos_ts_hist' in time_series_fns.keys():
            short_name = 'tos'
        elif 'chl_ts_hist' in time_series_fns.keys():
            short_name = 'chl'
        elif 'ph_ts_hist' in time_series_fns.keys():
            short_name = 'pH'
            joining = ' -'
        elif 'intpp_ts_hist' in time_series_fns.keys():
            short_name = 'intpp'
        elif 'po4_ts_hist' in time_series_fns.keys():
            short_name = 'po4'
        elif 'no3_ts_hist' in time_series_fns.keys():
            short_name = 'no3'
        elif 'mld_ts_hist' in time_series_fns.keys() or 'mlotst_ts_hist' in time_series_fns.keys():
            short_name = 'mld'
        elif 'o2_ts_hist' in time_series_fns.keys():
            short_name = 'o2'
        elif 'so_ts_hist' in time_series_fns.keys() or 'sos_ts_hist' in time_series_fns.keys():
            short_name = 'sal'
        else:
            print('suptitle not found:', time_series_fns.keys())
            assert 0

        suptitle = suptitles[short_name] # dict at start.

        suptitle += ' '.join([
                            joining,
                            'Recent Past', '('+ '-'.join([str(int(t)) for t in hist_time_range]) +')',
                            'vs Mid-Century Future', '('+'-'.join([str(int(t)) for t in ssp_time_range])+')' ])
        if single_model != 'all':
            suptitle = single_model+' '+suptitle
#        plt.suptitle(suptitle)

        # save and close.
        path = diagtools.folder(cfg['plot_dir']+'/whole_plot_modelrange')
        path += '_'.join(['multi_model_whole_plot', short_name, single_model])
        # if full_plot == 'ukesm':
        #     path += '_ukesm'
        path += diagtools.get_image_format(cfg)

        logger.info('Saving plots to %s', path)
        plt.savefig(path)
#        if single_model == 'all':
        plt.savefig(path.replace('.png', '_300dpi.png'), dpi=300, bbox_inches='tight')
        plt.savefig(path.replace('.png', '_300dpi_trans.png'), dpi=300, transparent=True)
        plt.savefig(path.replace('.png', '.pdf'), bbox_inches='tight')
        plt.savefig(path.replace('.png', '.svg'))
        plt.close()

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
