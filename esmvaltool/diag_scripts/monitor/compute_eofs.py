"""Diagnostic to compute and plot the first EOF of an arbitrary input.

It is part of the monitoring recipe. Imports class MonitorBase and
class PlotMap from mapgenerator in order to create plots, build
paths and record provenance.
"""
import logging

import iris
import matplotlib.pyplot as plt
from eofs.iris import Eof
from mapgenerator.plotting.plotmap import PlotMap

import esmvaltool.diag_scripts.shared
from esmvaltool.diag_scripts.monitor.monitor_base import MonitorBase
from esmvaltool.diag_scripts.shared import group_metadata

logger = logging.getLogger(__name__)


class Eofs(MonitorBase):
    """Diagnostic to compute EOFs and plot them.

    It is also an example on how to derive the monitor class to use its
    plotting capabilities in diagnostics that can not be done using only
    the preprocessor.
    """

    def compute(self):
        """Compute the diagnostic."""
        for module in ['matplotlib', 'fiona']:
            module_logger = logging.getLogger(module)
            module_logger.setLevel(logging.WARNING)

        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        # Loop over datasets
        for alias in data:
            # Loop over variables
            variables = group_metadata(data[alias], 'variable_group')
            for var_name, var_info in variables.items():
                logger.info('Plotting variable %s', var_name)
                var_info = var_info[0]
                # Load variable
                cube = iris.load_cube(var_info['filename'])
                # Initialise solver
                solver = Eof(cube, weights='coslat')
                # Get variable options as defined in monitor_config.yml
                variable_options = self._get_variable_options(
                    var_info['variable_group'], '')
                # Initialise PlotMap class
                plot_map = PlotMap(loglevel='INFO')
                # Compute EOF
                eof = solver.eofs(neofs=1)[0, ...]
                # Set metadata
                eof.long_name = var_info.get('eof_name', eof.long_name)
                eof.standard_name = None
                # Plot EOF map using plot_cube from PlotMap
                plot_map.plot_cube(eof, save=False, **variable_options)
                # Get filename for the EOF plot
                filename = self.get_plot_path('eof', var_info, 'png')
                # Save figure
                plt.savefig(filename,
                            bbox_inches='tight',
                            pad_inches=.2,
                            dpi=plot_map.dpi)
                plt.close(plt.gcf())
                # Record provenance for EOF plot
                self.record_plot_provenance(filename, var_info, 'eof')

                # Compute PC
                pcomp = solver.pcs(npcs=1, pcscaling=1)[:, 0]
                # Set metadata
                pcomp.long_name = var_info.get('pc_name', pcomp.long_name)
                pcomp.standard_name = None
                # Get filename for the PC plot
                filename = self.get_plot_path('pc', var_info)
                # Plot PC timeseries using plot_cube from MonitorBase
                self.plot_cube(pcomp, filename)
                # Record provenance for the PC plot
                self.record_plot_provenance(filename, var_info, 'pc')


def main():
    """Run EOFs diagnostic."""
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        Eofs(config).compute()


if __name__ == "__main__":
    main()
