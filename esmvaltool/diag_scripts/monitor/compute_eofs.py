"""Diagnostic to compute and plot the first EOF of an arbitrary input.

Description
-----------
This diagnostic can be used to compute and show Empirical Orthogonal Functions
(EOFs) and Principal Components (PCs) of arbitrary input. It creates a map plot
of the first EOF and the associated PC time series.

Configuration options in recipe
-------------------------------
cartopy_data_dir: str, optional (default: None)
    Path to cartopy data dir. Defaults to None. See
    https://scitools.org.uk/cartopy/docs/latest/.
config_file: str, optional
    Path to the monitor configuration file. Defaults to ``monitor_config.yml``
    in the same folder as the diagnostic script. More information on the
    monitor configuration file can be found :ref:`here <monitor_config_file>`.
plot_filename: str, optional
    Filename pattern for the plots.
    Defaults to ``{plot_type}_{real_name}_{dataset}_{mip}_{exp}_{ensemble}``.
    All tags (i.e., the entries in curly brackets, e.g., ``{dataset}``, are
    replaced with the corresponding tags).
plot_folder: str, optional
    Path to the folder to store figures. Defaults to
    ``{plot_dir}/../../{dataset}/{exp}/{modeling_realm}/{real_name}``.  All
    tags (i.e., the entries in curly brackets, e.g., ``{dataset}``, are
    replaced with the corresponding tags).  ``{plot_dir}`` is replaced with the
    default ESMValTool plot directory (i.e.,
    ``output_dir/plots/diagnostic_name/script_name/``, see
    :ref:`esmvalcore:user configuration file`).
rasterize_maps: bool, optional (default: True)
    If ``True``, use `rasterization
    <https://matplotlib.org/stable/gallery/misc/rasterization_demo.html>`_ for
    map plots to produce smaller files. This is only relevant for vector
    graphics (e.g., ``output_file_type=pdf,svg,ps``).

.. hint::

   Extra arguments given to the recipe are ignored, so it is safe to use yaml
   anchors to share the configuration of common arguments with other monitor
   diagnostic script.

"""
import logging
from copy import deepcopy

import iris
import matplotlib.pyplot as plt
from eofs.iris import Eof
from mapgenerator.plotting.plotmap import PlotMap

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.monitor.monitor_base import MonitorBase
from esmvaltool.diag_scripts.shared import group_metadata

logger = logging.getLogger(__name__)


class Eofs(MonitorBase):
    """Diagnostic to compute EOFs and plot them.

    It is also an example on how to derive the monitor class to use its
    plotting capabilities in diagnostics that can not be done using only
    the preprocessor.
    """

    def __init__(self, config):
        """Initialize class member."""
        super().__init__(config)

        # Get default settings
        self.cfg = deepcopy(self.cfg)
        self.cfg.setdefault('rasterize_maps', True)

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
                # Use rasterization if desired
                # Note: plt.gca() is the colorbar here, use plt.gcf().axes to
                # access the correct axes
                if self.cfg['rasterize_maps']:
                    self._set_rasterized(plt.gcf().axes[0])
                # Get filename for the EOF plot
                filename = self.get_plot_path('eof', var_info)
                # Save figure
                plt.savefig(filename,
                            bbox_inches='tight',
                            pad_inches=.2,
                            dpi=plot_map.dpi)
                plt.close(plt.gcf())
                # Record provenance for EOF plot
                caption = (f"{eof.long_name} of dataset {var_info[n.DATASET]} "
                           f"(project {var_info[n.PROJECT]}).")
                self.record_plot_provenance(filename, var_info, 'eof',
                                            caption=caption)

                # Compute PC
                pcomp = solver.pcs(npcs=1, pcscaling=1)[:, 0]
                # Set metadata
                pcomp.long_name = var_info.get('pc_name', pcomp.long_name)
                pcomp.standard_name = None
                # Get filename for the PC plot
                filename = self.get_plot_path('pc', var_info, add_ext=False)
                # Plot PC timeseries using plot_cube from MonitorBase
                self.plot_cube(pcomp, filename)
                # Record provenance for the PC plot
                caption = (f"{pcomp.long_name} of dataset "
                           f"{var_info[n.DATASET]} "
                           f"(project {var_info[n.PROJECT]}).")
                self.record_plot_provenance(
                    self._add_file_extension(filename),
                    var_info,
                    'pc',
                    caption=caption,
                )


def main():
    """Run EOFs diagnostic."""
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        Eofs(config).compute()


if __name__ == "__main__":
    main()
