import importlib
import inspect
import os

import esmvaltool.cmorizers.data.formatters.datasets as fdt


def test_formatters_have_required_interface():
    formatters_folder = os.path.dirname(fdt.__file__)
    unused_arg_names = ('_', '__', '___')
    arg_names = ('in_dir', 'out_dir', 'cfg', 'cfg_user', 'start_date',
                 'end_date')

    for formatter in os.listdir(formatters_folder):
        if not formatter.endswith('.py') or formatter == '__init__.py':
            continue
        module = formatter[:-3]
        member = importlib.import_module(
            f".{module}",
            package="esmvaltool.cmorizers.data.formatters.datasets")
        spec = inspect.getfullargspec(member.__getattribute__('cmorization'))
        try:
            assert len(spec.args) == 6
            for x, arg in enumerate(spec.args):
                assert arg == arg_names[x] or arg in unused_arg_names
        except AssertionError:
            print(f'Bad args in {os.path.join(formatters_folder, formatter)}: '
                  f'{spec.args}')
