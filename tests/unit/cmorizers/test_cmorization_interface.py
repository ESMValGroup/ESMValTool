import importlib
import inspect
import os

import esmvaltool.cmorizers.data.downloaders.datasets as ddt
import esmvaltool.cmorizers.data.formatters.datasets as fdt


def test_formatters_have_required_interface():
    formatters_folder = os.path.dirname(fdt.__file__)
    arg_names = (
        "in_dir",
        "out_dir",
        "cfg",
        "cfg_user",
        "start_date",
        "end_date",
    )
    unused_arg_names = ("_", "__", "___")

    error = False

    # aeronet.py and noaa_gml_surface_flask*.py need pys2index
    do_not_run_formatters = [
        '__init__.py',
        'aeronet.py',
        'noaa_gml_surface_flask.py',
        'noaa_gml_surface_flask_ch4.py',
        'noaa_gml_surface_flask_co2.py',
        'noaa_gml_surface_flask_n2o.py',
    ]
    all_formatters = os.listdir(formatters_folder)
    to_run_formatters = [
        f for f in all_formatters if f not in do_not_run_formatters
    ]

    for formatter in to_run_formatters:
        if not formatter.endswith(".py"):
            continue
        module = formatter[:-3]
        member = importlib.import_module(
            f".{module}",
            package="esmvaltool.cmorizers.data.formatters.datasets",
        )
        spec = inspect.getfullargspec(member.__getattribute__("cmorization"))
        try:
            assert len(spec.args) == len(arg_names)
            for x, arg in enumerate(spec.args):
                assert arg == arg_names[x] or arg in unused_arg_names
        except AssertionError:
            print(
                f"Bad args in {os.path.join(formatters_folder, formatter)}: "
                f"{spec.args}",
            )
            print(f"Expected {arg_names}.")
            error = True
    assert not error


def test_downloaders_have_required_interface():
    formatters_folder = os.path.dirname(ddt.__file__)
    arg_names = (
        "original_data_dir",
        "dataset",
        "dataset_info",
        "start_date",
        "end_date",
        "overwrite",
    )
    unused_arg_names = ("_", "__", "___")

    error = False

    for formatter in os.listdir(formatters_folder):
        if not formatter.endswith(".py") or formatter == "__init__.py":
            continue
        module = formatter[:-3]
        member = importlib.import_module(
            f".{module}",
            package="esmvaltool.cmorizers.data.downloaders.datasets",
        )
        spec = inspect.getfullargspec(
            member.__getattribute__("download_dataset"),
        )
        try:
            assert len(spec.args) == len(arg_names)
            for x, arg in enumerate(spec.args):
                assert arg == arg_names[x] or arg in unused_arg_names
        except AssertionError:
            print(
                f"Bad args in {os.path.join(formatters_folder, formatter)}: "
                f"{spec.args}",
            )
            error = True

    assert not error
