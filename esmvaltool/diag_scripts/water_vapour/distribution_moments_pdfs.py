import logging
import cf_units
import iris
from iris.util import equalise_attributes, unify_time_units
import iris.plot as iplt
import numpy as np
import matplotlib.pyplot as plt
import esmvalcore.preprocessor as eprep
import os
from fastkde import fastKDE

from esmvaltool.diag_scripts.water_vapour import distribution_moments_maps as dist_mom

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, select_metadata, group_metadata, get_diagnostic_filename, save_data, ProvenanceLogger
import esmvaltool.diag_scripts.shared.plot as eplot
from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import ProvenanceLogger

# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))

def plot_pdfs(data_dic, dataset, cfg):

    st_file = eplot.get_path_to_mpl_style(cfg.get('mpl_style'))
    plt.style.use(st_file)

    border_dic = {'mean': (0, 80), 'std': (0,25), 
                  'skew': (-2.5, 2.5), 'p95': (0,80)}
    seasons_dic = {'annual': 'k', 'winter': (68/255, 117/255, 191/255),
                                  'spring': (29/255, 149/255, 131/255),
                                  'summer': (214/255, 121/255, 96/255),
                                    'fall': (244/255, 182/255, 0) } 

    seass = list(data_dic.keys()); n_seas = len(seass)
    n_stat = len(data_dic[seass[0]].keys())
    n_cols = int(np.ceil(n_stat/2))

    fig_pdf, ax_pdf = plt.subplots(ncols=n_cols, nrows=2)
    fig_pdf.set_dpi(300) ; fig_pdf.set_size_inches(w=5.5*n_cols,h=7.5)
    ax_pdf = ax_pdf.flatten()

    for n,stat in enumerate(data_dic[seass[0]].keys()):
        for seas in seass:
            pdf, bins = fastKDE.pdf(data_dic[seas][stat].data.ravel())
            ax_pdf[n].plot(bins, pdf, c=seasons_dic[seas])
        ax_pdf[n].set_xlim(border_dic[stat][0], border_dic[stat][1])
        ylims = ax_pdf[n].get_ylim()
        ax_pdf[n].set_ylim(0, ylims[1])
        ax_pdf[n].set_xlabel('TCWV kg/m$^2$')
        ax_pdf[n].set_title(stat)
        if n%n_cols==0:
            ax_pdf[n].set_ylabel('pdf')

    fig_pdf.suptitle(dataset+' TCWV moments distributions in '+cfg['region']+' ('+cfg['time_range']+')')
    for n, seas in enumerate(seass):
        x = 1/(n_seas+1) - 0.05
        fig_pdf.text(x+n*x, 0.01, seas, c=seasons_dic[seas], fontsize='large')
    plt.tight_layout()
    fig_pdf.subplots_adjust(bottom=0.15)
    fig_pdf.savefig(os.path.join(cfg['plot_dir'], 'pdf_'+dataset+'_'+ cfg.get('region')+ diagtools.get_image_format(cfg)))
    plt.close(fig_pdf)

    return


def main(cfg):

    input_data = cfg['input_data']; dist_mom.reform_inp_data(input_data)

    datasets = group_metadata(input_data.values(), 'dataset', sort=True)

    for dataset in datasets.keys():
        dtst_dic = {} 
        groups = group_metadata(datasets[dataset], 'variable_group')
        for group in groups.keys():
            filepaths = list(group_metadata(groups[group], 'filename').keys())
            n_real = len(filepaths)
            if n_real == 1: 
                mod_cb = iris.load_cube(filepaths[0])
                mod_weight = np.ones(mod_cb.shape)
                dtst_dic[group] = dist_mom.calculate_single_cube_stats(mod_cb)
            else:
                mod_cblst = iris.load(filepaths)
                equalise_attributes(mod_cblst)
                [cb.add_aux_coord(iris.coords.AuxCoord(n, long_name='n_order',
                        var_name='n_order')) for n, cb in enumerate(mod_cblst)]
                mod_cb = mod_cblst.merge_cube()
                mod_weight = np.ones(mod_cb.shape)/n_real
                dtst_dic[group] = dist_mom.calculate_cube_stats(mod_cb, mod_weight)
        plot_pdfs(dtst_dic, dataset, cfg)

    logger.info('Success')


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)