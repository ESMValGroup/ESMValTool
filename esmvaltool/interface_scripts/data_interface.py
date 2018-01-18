"""
Create file(s) for the interface between esmvaltool and diagnostic scripts.

Completely rewritten tool to be able to deal with
the new yaml parser and simplified interface_scripts
toolbox. Author: Valeriu Predoi, University of Reading,
Initial version: August 2017
contact: valeriu.predoi@ncas.ac.uk
"""
import logging
import os

from ..version import __version__
from .data_finder import get_output_file

logger = logging.getLogger(__name__)


def write_ncl_settings(settings, filename, mode='wt'):
    """Write settings to NCL file."""
    logger.debug("Writing NCL configuration file %s", filename)

    def _format(value):
        """Format string or list as NCL"""
        if value is None or isinstance(value, str):
            txt = '"{}"'.format(value)
        elif isinstance(value, (list, tuple)):
            txt = '(/{}/)'.format(', '.join(_format(v) for v in value))
        else:
            txt = str(value)
        return txt

    def _format_dict(name, dictionary):
        """Format dict as NCL"""
        lines = ['{} = True'.format(name)]
        for key, value in sorted(dictionary.items()):
            lines.append('{}@{} = {}'.format(name, key, _format(value)))
        txt = '\n'.join(lines)
        return txt

    def _header(name):
        """Delete any existing NCL variable known as `name`."""
        return ('if (isvar("{name}")) then\n'
                '    delete({name})\n'
                'end if\n'.format(name=name))

    lines = []
    for key, value in sorted(settings.items()):
        txt = _header(name=key)
        if isinstance(value, dict):
            txt += _format_dict(name=key, dictionary=value)
        else:
            txt += '{} = {}'.format(key, _format(value))
        lines.append(txt)
    with open(filename, mode) as file:
        file.write('\n\n'.join(lines))
        file.write('\n')


def get_legacy_ncl_interface(variables, settings, namelist_file, script):
    """Get a dictionary with the contents of the former ncl.interface file."""
    esmvaltool_root = os.path.dirname(os.path.dirname(__file__))
    project_root = os.path.dirname(esmvaltool_root)

    # perfmetrics_main.ncl expects diag_script to be a relative path
    diag_path = os.path.join(esmvaltool_root, 'diag_scripts')
    if script.startswith(diag_path):
        script = os.path.relpath(script, diag_path)

    ncl_interface = {
        'diag_script':
        script,
        'variable_def_dir':
        os.path.join(esmvaltool_root, 'variable_defs'),
        'in_refs':
        [os.path.join(project_root, 'doc', 'MASTER_authors-refs-acknow.txt')],
        'out_refs':
        [os.path.join(settings['run_dir'], 'references-acknowledgements.txt')],
        'yml': [namelist_file],
        'yml_name': [os.path.basename(namelist_file)],
        'output_file_type': [settings['output_file_type']],
        'plot_dir': [settings['plot_dir']],
        'work_dir': [settings['work_dir']],
        'regridding_dir': [],
        'write_netcdf': [settings['write_netcdf']],
        'read_from_vault': [],
        'cwd': [os.getcwd()],
        'force_processing': [],
        'show_debuginfo': [],
        'show_diag_description': [],
        # climate.ncl variables
        'infilename': [],
        'mfile': [],
        'sfile': [],
        'afile': [],
        'base_variable': [],
        'max_data_filesize': [settings['max_data_filesize']],
        'fx_keys': [],
        'fx_values': [],
        'str_vault_sep':
        "-",
        # Data structures to hold information on the
        # reference/acknowledgement output file
        'ref_auth': [],
        'ref_contr': [],
        'ref_diag': [],
        'ref_obs': [],
        'ref_proj': [],
        'ref_script': [],
    }

    # Add info that requires a variable
    if variables:
        # select random variable
        tmp, single_variable = variables.popitem()
        variables[tmp] = single_variable

        # ordered variables
        variable_keys = sorted(variables)

        infiles = [
            get_output_file_template(v, settings['preproc_dir'])
            for v in single_variable
        ]

        # ref_model (only works if the same for all variables in diagnostic)
        ref_model = [single_variable[0].get('reference_model')]
        if 'alternative_model' in single_variable[0]:
            ref_model.append(single_variable[0]['alternative_model'])

        ncl_interface.update({
            'dictkeys': {
                'dictkeys': [get_dict_key(v) for v in single_variable],
            },
            'figfiles_suffix':
            [get_figure_file_names(v) for v in single_variable],
            'infile_paths': [os.path.dirname(p) for p in infiles],
            'infiles': [os.path.basename(p) for p in infiles],
            'fullpaths': [p for p in infiles],
            'var_attr_mip': [v.get('mip') for v in single_variable],
            'var_attr_exp': [v.get('exp') for v in single_variable],
            'var_attr_ref':
            ref_model,
            'var_attr_exclude':
            [v.get('exclude', "False") for v in single_variable],
            'model_attr_skip': [v.get('skip') for v in single_variable],
            'variables':
            [variables[k][0]['short_name'] for k in variable_keys],
            'derived_var': [
                variables[k][0]['short_name'] for k in variable_keys
            ],
            'field_types': [variables[k][0]['field'] for k in variable_keys],
            'derived_field_type': [
                variables[k][0]['field'] for k in variable_keys
            ],
            'models': {
                'project': [v['project'] for v in single_variable],
                'name': [v['model'] for v in single_variable],
                'mip': [v.get('mip', 'mip') for v in single_variable],
                'experiment': [v.get('exp') for v in single_variable],
                'ensemble':
                [v.get('ensemble', 'ensemble') for v in single_variable],
                'start_year': [v['start_year'] for v in single_variable],
                'end_year': [v['end_year'] for v in single_variable],
                'freq': [v.get('freq') for v in single_variable],
                'dir': [v.get('dir') for v in single_variable],
                'level': [v.get('level') for v in single_variable],
                'case_name': [v.get('case_name') for v in single_variable],
            },
            'model_attr_id': [v.get('model') for v in single_variable],
        })

    # remove items with no value
    ncl_interface = {k: v for k, v in ncl_interface.items() if v}

    # remove sub-dicts with lists containing no values
    for key, value in ncl_interface.items():
        if isinstance(value, dict):
            ncl_interface[key] = {k: v for k, v in value.items() if any(v)}

    return ncl_interface


def write_legacy_ncl_interface(variables, settings, namelist_file, script):
    """Write legacy ncl interface files."""
    # get legacy ncl interface dictionary
    ncl_interface = get_legacy_ncl_interface(variables, settings,
                                             namelist_file, script)
    # add namelist script settings
    ncl_interface['diag_script_info'] = settings

    # write ncl.interface
    run_dir = settings['run_dir']
    if not os.path.isdir(run_dir):
        os.makedirs(run_dir)
    interface_file_tmp = os.path.join(run_dir, 'interface.ncl')
    write_ncl_settings(ncl_interface, interface_file_tmp)
    interface_file = os.path.join(run_dir, 'ncl.interface')
    os.rename(interface_file_tmp, interface_file)
    logger.info("with configuration file %s", interface_file)

    # variable info files
    for name, variable in variables.items():
        info_file_tmp = os.path.join(run_dir, name + '_info.ncl')
        common_items = {
            k: v
            for k, v in variable[0].items()
            if all(v == w.get(k) for w in variable)
        }
        variable_info = {'variable_info': common_items}
        write_ncl_settings(variable_info, info_file_tmp, mode='at')
        info_file = os.path.splitext(info_file_tmp)[0] + '.tmp'
        logger.info("and configuration file %s", info_file)
        os.rename(info_file_tmp, info_file)


def get_legacy_ncl_env(settings, namelist_basename):
    """Get legacy ncl environmental variables."""
    project_root = os.sep.join(__file__.split(os.sep)[:-3])
    prefix = 'ESMValTool_'

    env = dict(os.environ)
    for key in ('work_dir', 'plot_dir', 'output_file_type', 'write_plots',
                'write_netcdf'):
        env[prefix + key] = settings[key]
    env[prefix + 'interface_data'] = settings['run_dir']
    env['0_ESMValTool_version'] = __version__
    env[prefix
        + 'verbosity'] = 100 if settings['log_level'].lower() == 'debug' else 1
    env[prefix + 'in_refs'] = os.path.join(project_root, 'doc',
                                           'MASTER_authors-refs-acknow.txt')
    env[prefix + 'out_refs'] = os.path.join(settings['run_dir'],
                                            'references-acknowledgements.txt')
    env[prefix + 'yml_name'] = namelist_basename

    return env


def get_output_file_template(variable, preproc_dir):
    """Get output filename for backward compatible ncl interface."""
    variable = dict(variable)
    variable['short_name'] = '${VARIABLE}'
    variable['field'] = '${FIELD}'
    return get_output_file(variable, preproc_dir)


def get_figure_file_names(variable):
    """Get plot filename."""
    #     return "_".join([
    #         model['project'],
    #         model['name'],
    #         model['mip'],
    #         model['exp'],
    #         model['ensemble'],
    #         str(model['start_year']) + "-" + str(model['end_year']),
    #     ])
    return "_".join([
        variable['project'],
        variable['model'],
        str(variable['start_year']) + "-" + str(variable['end_year']),
    ])


def get_dict_key(model):
    """Return a unique key based on the model entries provided.

    This function creates and returns a key used in a number of NCL
    scripts to refer to a specific dataset, see e.g., the variable
    'cn' in 'interface_scripts/read_data.ncl'
    """
    # allow for different projects

    # CMIP5
    if model['project'] == 'CMIP5':
        dict_key = "_".join([
            model['project'],
            model['model'],
            model['mip'],
            model['exp'],
            model['ensemble'],
            str(model['start_year']),
            str(model['end_year']),
        ])
    else:
        dict_key = "_".join([
            model['project'],
            model['model'],
            str(model['start_year']),
            str(model['end_year']),
        ])

    return dict_key
