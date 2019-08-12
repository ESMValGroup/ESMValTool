import os
import logging

import matplotlib.pyplot as plt

import iris
import iris.cube
import iris.analysis
import iris.util
import iris.quickplot as qp

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata

logger = logging.getLogger(os.path.basename(__file__))


class JetLatitude(object):
    def __init__(self, config):
        self.cfg = config
        self.filenames = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.output_name = None
        self.target_grid = self.cfg.get('target_grid')
        self.grid_cube = None
        self.sftlf = None

    def compute(self):
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        for alias in data:
            ua = iris.load_cube(data[alias][0]['filename'])
            logger.debug(ua)
            qp.pcolor(ua)
            plt.savefig(os.path.join(
                self.cfg[n.PLOT_DIR], '{}.png'.format(alias)))


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        JetLatitude(config).compute()


if __name__ == '__main__':
    main()
