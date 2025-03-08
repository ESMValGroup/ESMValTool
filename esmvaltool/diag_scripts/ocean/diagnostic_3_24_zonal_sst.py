"""
Diagnostics for fig 3.24 in Chapter 3 of IPCC AR6 WGI
========================


Author: Lee de Mora (PML)
        ledm@pml.ac.uk
Revised and corrected (15.01.2021): Elizaveta Malinina (CCCma)
                        elizaveta.malinina-rieger@canada.ca

"""
import logging
import os
import sys

import iris
import iris.plot
import matplotlib.pyplot as plt
import numpy as np
import esmvalcore.preprocessor as eprep

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import (run_diagnostic, 
                                            group_metadata, 
                                            select_metadata, 
                                            save_figure)
import esmvaltool.diag_scripts.shared.plot as eplot

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

class Data4Analyis:

    def __init__(self, name: str, group: list, cfg: dict, 
                mask_type: str | bool, mask_meta: list[dict]| None):

        self.name = name
        self.bias_type = 'bias' if cfg.get('bias') else None
        self.mask = cfg.get('mask').get('flag') if cfg.get('mask') else False
        self.mask_type = mask_type
        self.determine_reference(group)
        self.obtain_data(group, mask_meta)
        stats = cfg.get('statistics')
        if not(stats): 
            raise ValueError("statistics dictionary should be provided in the recipe. "
                             "The keywords 'best_guess' and 'borders' should be provided.")
        self.calculate_statistics(stats)


    def determine_reference(self, group: list):

        if self.bias_type or self.mask:
            reference = list(group_metadata(group, 'reference_dataset').keys())
            if len(reference) > 1: 
                raise ValueError('More then one reference dataset was given')
            else: 
                self.reference = reference[0] 
            if reference:
                ref_metadata = select_metadata(group, dataset=self.reference)[0]
                self.ref_cube = iris.load_cube(ref_metadata['filename'])
                group.remove(ref_metadata)
            else: 
                raise ValueError('No reference dataset has been provided')


    def obtain_data(self, group: list, mask_meta: list[dict]| None):

        files = list(group_metadata(group, 'filename', sort=True).keys())
        data_cblst = iris.load(files)
        self.data = data_cblst

        if self.mask:
            self.mask_data(mask_meta)

        if self.bias_type:
            self.data = eprep.bias(self.data, reference=self.ref_cube)

        return
    
    def mask_data(self, mask_meta: list[dict]| None):

        if self.mask_type == 'simple':
            mask = self.ref_cube.data.mask
        elif self.mask_type == 'resolved':
            mask_files = list(group_metadata(mask_meta, 'filename', sort=True).keys())
            if len(mask_files) > 1:
                raise ValueError("More than one dataset for the resolved mask "
                                 "has been provided. Only one is supported.")
            mask_cb = iris.load_cube(mask_files[0])
            data_coord = self.data[0].dim_coords
            if len(data_coord)>1:
                raise ValueError("The data cubes have more than one "
                            "cordinates. Only flat cubes are supported.")
            data_coord = data_coord[0]
            dim = [c.name() for c in mask_cb.dim_coords].index(data_coord.name())
            mask = list()
            for i in range(mask_cb.shape[dim]):
                mask.append(mask_cb.data[(slice(None),) * dim + (i,)].count() <= 0.5 * len(mask_cb.data[(slice(None),) * dim + (i,)]))
            # add masking of a ref cube 
        else:
            raise ValueError(f"Mask type {self.mask_type} is not supported. "
                             "Only 'simple' and 'resolved' are supported.")


        for n_cb in range(len(self.data)):
            self.data[n_cb].data.mask = self.data[n_cb].data.mask | mask

        return

    def calculate_statistics(self, stats: dict):

        if len(self.data)>1:
            bg_dic = eprep.multi_model_statistics(
                                        self.data, span='full', 
                                        statistics=[stats['best_guess']],
                                        ignore_scalar_coords=True)
            self.best_guess = bg_dic[list(bg_dic.keys())[0]]
            bord_dic = eprep.multi_model_statistics(
                                        self.data, span='full', 
                                        statistics=stats['borders'], 
                                        ignore_scalar_coords=True)
            self.border1 = bord_dic[list(bord_dic.keys())[0]]
            self.border2 = bord_dic[list(bord_dic.keys())[1]]
        else:
            self.best_guess = self.data[0]
            self.border1 = None
            self.border2 = None 



def create_provenance(caption: str):
    """Creates provenance dictionary."""

    provenance_dic = {
        'authors': ['malinina_elizaveta', 'demora_lee'],
        'caption': caption,
        'references': ['eyring21ipcc']
    }

    return provenance_dic    


def plot_bias_plot(data_list : list[Data4Analyis], cfg: dict):

    # add caption 
    prov_dic = create_provenance('')

    st_file = eplot.get_path_to_mpl_style(cfg.get('mpl_style'))
    plt.style.use(st_file)
    
    fig = plt.figure(figsize=(6, 2.5))
    
    for data in data_list:
        if data.best_guess.dim_coords[0].name() == 'longitude':
            # to make the Pacific and Atlantic continuous
            data.best_guess = data.best_guess.intersection(longitude=(20., 380.))
            data.border1 = data.border1.intersection(longitude=(20., 380.))
            data.border2 = data.border2.intersection(longitude=(20., 380.))
        iris.plot.plot(data.best_guess)
        iris.plot.fill_between(data.best_guess.dim_coords[0], data.border1, data.border2, alpha = 0.2, linewidth=0)
    
    # add ocean borders for the equatorial 
    # add horizontal line for bias
    # add reading style colors
    # add plotting ref cube for the non bias

    plt.tight_layout()
    fig_path = os.path.join(cfg['plot_dir'],'sst_bias')
    save_figure(fig_path, prov_dic, cfg, fig, close=True)

    return

def main(cfg : dict):

    input_data = cfg['input_data']

    groups = group_metadata(input_data.values(), 'variable_group', sort=True)

    mask_type = cfg.get('mask').get('type') if cfg.get('mask') else False
    mask_group = cfg.get('mask').get('group') if mask_type == 'resolved' else None
    mask_meta = groups.pop(mask_group) if mask_group else None

    data_list = []
    for group in groups.keys(): 
        group_data = Data4Analyis(name=group, group=groups[group], cfg=cfg, 
                                  mask_type=mask_type, mask_meta=mask_meta)
        data_list.append(group_data)

    plot_bias_plot(data_list, cfg)


    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)