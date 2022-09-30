import importlib
import inspect
import os

import esmvaltool.cmorizers.data.downloaders.datasets as ddt
import esmvaltool.cmorizers.data.formatters.datasets as fdt


def test_formatters_have_required_interface():
    formatters_folder = os.path.dirname(fdt.__file__)
    arg_names = ('in_dir', 'out_dir', 'cfg', 'cfg_user', 'start_date',
                 'end_date')
    unused_arg_names = ('_', '__', '___')

    error = False

    for formatter in os.listdir(formatters_folder):
        if not formatter.endswith('.py') or formatter == '__init__.py':
            continue
        module = formatter[:-3]
        member = importlib.import_module(
            f".{module}",
            package="esmvaltool.cmorizers.data.formatters.datasets")
        spec = inspect.getfullargspec(member.__getattribute__('cmorization'))
        try:
            assert len(spec.args) == len(arg_names)
            for x, arg in enumerate(spec.args):
                assert arg == arg_names[x] or arg in unused_arg_names
        except AssertionError:
            print(f'Bad args in {os.path.join(formatters_folder, formatter)}: '
                  f'{spec.args}')
            error = True
    assert not error


def test_downloaders_have_required_interface():
    formatters_folder = os.path.dirname(ddt.__file__)
    arg_names = ('config', 'dataset', 'dataset_info', 'start_date', 'end_date',
                 'overwrite')
    unused_arg_names = ('_', '__', '___')

    error = False

    for formatter in os.listdir(formatters_folder):
        if not formatter.endswith('.py') or formatter == '__init__.py':
            continue
        module = formatter[:-3]
        member = importlib.import_module(
            f".{module}",
            package="esmvaltool.cmorizers.data.downloaders.datasets")
        spec = inspect.getfullargspec(
            member.__getattribute__('download_dataset'))
        try:
            assert len(spec.args) == len(arg_names)
            for x, arg in enumerate(spec.args):
                assert arg == arg_names[x] or arg in unused_arg_names
        except AssertionError:
            print(f'Bad args in {os.path.join(formatters_folder, formatter)}: '
                  f'{spec.args}')
            error = True

    assert not error
