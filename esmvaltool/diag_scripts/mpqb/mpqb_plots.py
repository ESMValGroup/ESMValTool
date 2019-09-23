
import iris
import matplotlib.pyplot as plt
import copy

cmap_dict = {
    'sm' : {
        'sequential' : 'YlOrBr',
        'diverging'  : 'RdYlBu',
         },
    'lai' : {
        'sequential' : 'YlOrBr',
        'diverging'  : 'BrBG',
        }
}


metrics_plot_dictionary = {
    'pearsonr' : {
        'title' : 'pearsonr',
        'vmin' : -1.,
        'vmax' : 1.,
        'cmap' : 'diverging', 
    },
    'rmsd' : {
        'title' : 'rmsd',
        'vmin' : 0.,
        'vmax' : 1.,
        'cmap' : 'sequential',
    },
    'absdiff' : {
        'title' : 'absdiff',
        'vmin' : -.5,
        'vmax' : .5,
        'cmap' : 'diverging',
    },
    'reldiff' : {
        'title' : 'reldiff',
        'cmap' : 'diverging',
        'vmin' : -150,
        'vmax' : 150,
    },
    'theilsenmk' : {
        'title' : 'theilsen',
        'cmap' : 'diverging',
        'vmin' : -0.05,
        'vmax' : 0.05,
    },
    'timemean' : {
        'title' : 'timemean',
        'cmap' : 'sequential',
        'vmin' : 0.0,
        'vmax' : 0.5,
    }
}

def get_plot_config(ecv_name):
    plot_config = copy.deepcopy(metrics_plot_dictionary)
    for metricname in plot_config:
        plot_config[metricname]['cmap'] = cmap_dict[ecv_name][plot_config[metricname]['cmap']]
    return plot_config

def mpqb_mapplot(cube,filename,**plotkwargs):
    plottitle = plotkwargs.pop('title')
    fig = plt.figure()
    ax = fig.add_subplot(projection=iris.plot.default_projection(cube))
    iris.quickplot.pcolormesh(cube,**plotkwargs)
    plt.gca().coastlines()
    plt.title(plottitle)
    fig.savefig(filename)
    plt.close(fig)
