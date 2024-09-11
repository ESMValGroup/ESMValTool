"""Diagnostic script to plot extent and differences.

based on code from Anton Steketee's COSIMA recipes notebook
https://cosima-recipes.readthedocs.io/en/latest/Examples/Sea_Ice_Area_Concentration_Volume_with_Obs.html
"""

import logging
import calendar
from cartopy.crs import SouthPolarStereo

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np

from esmvalcore.preprocessor import extract_month
from esmvaltool.diag_scripts.shared import (run_diagnostic, save_figure,
                                            group_metadata, save_data)

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def map_fig_diff(mod_dict, obs_si, months):
    """Create map with model dictionary: labels, cubes"""

    # fig set up, width for 2 models, check len mod_dict
    figure = plt.figure(figsize=(9,len(months)*3.5))
    j=0  # to iterate through positions on figure

    for mon in months: #[2,9]
        i=1
        for mod_label, mod_si in mod_dict.items(): #iterate over model subset
            
            mod_cube = extract_month(mod_si, mon)
            obs_cube = extract_month(obs_si, mon)
            
            out = mod_cube.copy()
            diff = mod_cube.data - obs_cube.data
            out.data = diff
        
            ax = plt.subplot(len(months), 3, i+j*3, projection=crs.SouthPolarStereo(true_scale_latitude=-70))
            
            diffmap = iplt.contourf(out, levels=np.arange(-90,91,20), cmap='RdBu') #?
            
            iplt.contour(obs_cube, levels=[15], colors=['yellow'])
            iplt.contour(mod_cube, levels=[15], linewidths=1.0, colors=['black'])
            
            plt.title(calendar.month_abbr[mon]+' '+mod_label)

            i+=1
        j+=1

    line_cdr = mlines.Line2D([], [], color='yellow', label="Observed Extent")
    line_mod = mlines.Line2D([], [], color='black', label="Modelled Extent")
    
    plt.legend(handles=[line_cdr,line_mod], loc='center left', bbox_to_anchor=(1.2,0.5))
    cax = plt.axes([0.7,0.55,0.04,0.3])
    _ = plt.colorbar(diffmap, cax=cax, label='Difference in \nSea Ice Concentration')
    
    plt.subplots_adjust(left=0.05, bottom=0.05, 
                    right=0.95, top=0.95,
                    wspace=0.05, hspace=0.05)

    return figure


def main(cfg):
    """Compute."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()
    
    groups = group_metadata(input_data, 'variable_group', sort='project')
    for group_name in groups:
        mod_dict={}
        ancestor_filels = []
        logger.info("Processing variable %s", group_name)
        for attributes in groups[group_name]:
            if attributes['project'].startswith('OBS'):
                logger.info("Processing OBS dataset %s", attributes['dataset'])
                obs_si = iris.load_cube(attributes['filename'])
            else:
                logger.info("Processing dataset %s", attributes['dataset'])
                mod_dict[attributes['dataset']] = iris.load_cube(attributes['filename'])
            ancestor_filels.append(attributes['filename'])

        logger.info("creating map differences")
        mapfig = map_fig_diff(mod_dict, obs_si, cfg['months'])

        provenance_record = get_provenance_record(ancestor_filels)
        save_figure(group_name, provenance_record, cfg, figure=mapfig)


def get_provenance_record(ancestor_files):
    """Build provenance dictionary."""
    record = {
        'ancestors': ancestor_files,
        'authors': ['chun_felicity', 'steketee_anton'],
        'caption': '',
        'domains': ['shpolar'],
        'plot_types': ['polar'],
        'references': [],
        'statistics': ['diff'],
    }
    return record


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
