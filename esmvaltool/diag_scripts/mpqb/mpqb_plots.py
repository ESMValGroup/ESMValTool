
import iris
import matplotlib.pyplot as plt

metrics_plot_dictionary = {
    'pearsonr' : {
        'title' : 'pearsonr',
        'vmin' : -1.,
        'vmax' : 1.,
        'cmap' : 'RdYlBu_r', # diverging
    },
    'rmsd' : {
        'title' : 'rmsd',
        'vmin' : 0.,
        'vmax' : 1.,
        'cmap' : 'YlOrBr', # sequential
    },
    'absdiff' : {
        'title' : 'absdiff',
        'vmin' : -.5,
        'vmax' : .5,
        'cmap' : 'RdYlBu_r', # diverging
    },
    'reldiff' : {
        'title' : 'reldiff',
        'cmap' : 'RdYlBu_r', # diverging
        'vmin' : -150,
        'vmax' : 150,
    },
    'theilsen' : {
        'title' : 'theilsen',
        'cmap' : 'RdYlBu', # diverging, blue -> wettening, red drying
        'vmin' : -0.001,
        'vmax' : 0.001,
    },
    'timemean' : {
        'title' : 'timemean',
        'cmap' : 'YlGnBu', # sequential 
        'vmin' : 0.0,
        'vmax' : 0.5,
    }
}



def mpqb_mapplot(cube,filename,**plotkwargs):
    plottitle = plotkwargs.pop('title')
    fig = plt.figure()
    ax = fig.add_subplot(projection=iris.plot.default_projection(cube))
    iris.quickplot.pcolormesh(cube,**plotkwargs)
    plt.gca().coastlines()
    plt.title(plottitle)
    fig.savefig(filename)
    plt.close(fig)
