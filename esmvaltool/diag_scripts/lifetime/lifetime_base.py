"""Base class for lifetime diagnostics."""

import logging
import os
import re

import matplotlib.pyplot as plt
import yaml

from esmvaltool.diag_scripts.shared import ProvenanceLogger, names

logger = logging.getLogger(__name__)


def _replace_tags(paths, variable):
    """Replace tags in the config-developer's file with actual values."""
    if isinstance(paths, str):
        paths = set((paths.strip('/'), ))
    else:
        paths = set(path.strip('/') for path in paths)
    tlist = set()
    for path in paths:
        tlist = tlist.union(re.findall(r'{([^}]*)}', path))
    if 'sub_experiment' in variable:
        new_paths = []
        for path in paths:
            new_paths.extend(
                (re.sub(r'(\b{ensemble}\b)', r'{sub_experiment}-\1', path),
                 re.sub(r'({ensemble})', r'{sub_experiment}-\1', path)))
            tlist.add('sub_experiment')
        paths = new_paths

    for tag in tlist:
        original_tag = tag
        tag, _, _ = _get_caps_options(tag)

        if tag == 'latestversion':  # handled separately later
            continue
        if tag in variable:
            replacewith = variable[tag]
        else:
            raise ValueError(f"Dataset key '{tag}' must be specified for "
                             f"{variable}, check your recipe entry")
        paths = _replace_tag(paths, original_tag, replacewith)
    return paths


def _replace_tag(paths, tag, replacewith):
    """Replace tag by replacewith in paths."""
    _, lower, upper = _get_caps_options(tag)
    result = []
    if isinstance(replacewith, (list, tuple)):
        for item in replacewith:
            result.extend(_replace_tag(paths, tag, item))
    else:
        text = _apply_caps(str(replacewith), lower, upper)
        result.extend(p.replace('{' + tag + '}', text) for p in paths)
    return list(set(result))


def _get_caps_options(tag):
    lower = False
    upper = False
    if tag.endswith('.lower'):
        lower = True
        tag = tag[0:-6]
    elif tag.endswith('.upper'):
        upper = True
        tag = tag[0:-6]
    return tag, lower, upper


def _apply_caps(original, lower, upper):
    if lower:
        return original.lower()
    if upper:
        return original.upper()
    return original


class LifetimeBase():
    """Base class for lifetime diagnostic.

    It contains the common methods for path creation, provenance
    recording, option parsing and to create some common plots.

    """

    def __init__(self, config):
        self.cfg = config
        plot_folder = config.get(
            'plot_folder',
            '{plot_dir}/../../{dataset}/{exp}/{modeling_realm}/{real_name}',
        )
        plot_folder = plot_folder.replace('{plot_dir}',
                                          self.cfg[names.PLOT_DIR])
        self.plot_folder = os.path.abspath(plot_folder)
        self.plot_filename = config.get(
            'plot_filename',
            '{plot_type}_{real_name}_{dataset}_{mip}_{exp}_{ensemble}')
        self.plots = config.get('plots', {})
        default_config = os.path.join(os.path.dirname(__file__),
                                      "lifetime_config.yml")
        # cartopy_data_dir = config.get('cartopy_data_dir', None)
        # if cartopy_data_dir:
        #     cartopy.config['data_dir'] = cartopy_data_dir
        with open(config.get('config_file', default_config),
                  encoding='utf-8') as config_file:
            self.config = yaml.safe_load(config_file)

    def _add_file_extension(self, filename):
        """Add extension to plot filename."""
        return f"{filename}.{self.cfg['output_file_type']}"

    def get_plot_path(self, plot_type, var_info, add_ext=True):
        """Get plot full path from variable info.

        Parameters
        ----------
        plot_type: str
            Name of the plot
        var_info: dict
            Variable information from ESMValTool
        add_ext: bool, optional (default: True)
            Add filename extension from configuration file.

        """
        return os.path.join(
            self.get_plot_folder(var_info),
            self.get_plot_name(plot_type, var_info, add_ext=add_ext),
        )

    def get_plot_folder(self, var_info):
        """Get plot storage folder from variable info.

        Parameters
        ----------
        var_info: dict
            Variable information from ESMValTool

        """
        info = {
            'real_name': self._real_name(var_info['variable_group']),
            **var_info
        }
        folder = os.path.expandvars(
            os.path.expanduser(
                list(_replace_tags(self.plot_folder, info))[0]
            )
        )
        if self.plot_folder.startswith('/'):
            folder = '/' + folder
        if not os.path.isdir(folder):
            os.makedirs(folder, exist_ok=True)
        return folder

    def get_plot_name(self, plot_type, var_info, add_ext=True):
        """Get plot filename from variable info.

        Parameters
        ----------
        plot_type: str
            Name of the plot
        var_info: dict
            Variable information from ESMValTool
        add_ext: bool, optional (default: True)
            Add filename extension from configuration file.

        """
        info = {
            "plot_type": plot_type,
            'real_name': self._real_name(var_info['variable_group']),
            **var_info
        }
        file_name = list(_replace_tags(self.plot_filename, info))[0]
        if add_ext:
            file_name = self._add_file_extension(file_name)
        return file_name

    @staticmethod
    def _set_rasterized(axes=None):
        """Rasterize all artists and collection of axes if desired."""
        if axes is None:
            axes = plt.gca()
        if not isinstance(axes, list):
            axes = [axes]
        for single_axes in axes:
            for artist in single_axes.artists:
                artist.set_rasterized(True)
            for collection in single_axes.collections:
                collection.set_rasterized(True)

    @staticmethod
    def _real_name(variable_group):
        for subfix in ('Ymean', 'Ysum', 'mean', 'sum'):
            if variable_group.endswith(subfix):
                variable_group = variable_group.replace(subfix, '')
        return variable_group
