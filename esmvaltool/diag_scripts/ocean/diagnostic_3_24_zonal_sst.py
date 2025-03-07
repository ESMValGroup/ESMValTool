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
import matplotlib.pyplot as plt
import numpy as np
import esmvalcore.preprocessor as eprep

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic,  group_metadata, select_metadata
import esmvaltool.diag_scripts.shared.plot as eplot

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

class Data4Analyis:

    def __init__(self, name: str, group: list, cfg: dict):

        self.name = name
        self.bias_type = 'bias' if cfg.get('bias') else None
        self.mask = cfg.get('mask').get('flag') if cfg.get('mask') else False
        self.mask_type = cfg.get('mask').get('type') if cfg.get('mask') else False
        self.determine_reference(group)
        self.obtain_data(group)
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


    def obtain_data(self, group: list):

        files = list(group_metadata(group, 'filename', sort=True).keys())
        data_cblst = iris.load(files)
        self.data = data_cblst

        if self.mask:
            self.mask_data()

        if self.bias_type:
            self.data = eprep.bias(self.data, reference=self.ref_cube)

        return
    
    def mask_data(self):

        mask = self.ref_cube.data.mask

        for n_cb in range(len(self.data)):
            self.data[n_cb].data.mask = self.data[n_cb].data.mask & mask

        # add complex mask

        return

    def calculate_statistics(self, stats: dict):

        if len(self.data)>1:
            bg_dic = eprep.multi_model_statistics(self.data, span='full', statistics=[stats['best_guess']], ignore_scalar_coords=True)
            self.best_guess = bg_dic[list(bg_dic.keys())[0]]
            bord_dic = eprep.multi_model_statistics(self.data, span='full', statistics=stats['borders'], ignore_scalar_coords=True)
            self.border1 = bord_dic[list(bord_dic.keys())[0]]
            self.border2 = bord_dic[list(bord_dic.keys())[1]]
        else:
            self.best_guess = self.data[0]
            self.border1 = None
            self.border2 = None 


def plot_bias_plot(data_list:list, cfg: dict):

    # if the longitude, redefine for the oceans 

    return

def main(cfg : dict):

    input_data = cfg['input_data']

    groups = group_metadata(input_data.values(), 'variable_group', sort=True)

    data_list = []
    for group in groups.keys(): 
        group_data = Data4Analyis(name=group, group=groups[group], cfg=cfg)
        data_list.append(group_data)

    #  plot 


    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)