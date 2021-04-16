import logging
import sys
from pathlib import Path

import pytest
import yaml

from esmvaltool.diag_scripts import shared


def test_get_plot_filename():

    cfg = {
        'plot_dir': '/some/path',
        'output_file_type': 'png',
    }
    filename = shared.get_plot_filename('test', cfg)
    assert filename == '/some/path/test.png'


def test_get_diagnostic_filename():

    cfg = {
        'work_dir': '/some/path',
    }
    filename = shared.get_diagnostic_filename('test', cfg)
    assert filename == '/some/path/test.nc'


def test_get_diagnostic_filename_ext():

    cfg = {
        'work_dir': '/some/path',
    }
    filename = shared.get_diagnostic_filename('test', cfg, extension='csv')
    assert filename == '/some/path/test.csv'


def test_provenance_logger(tmp_path):

    record = {'attribute1': 'xyz'}
    with shared.ProvenanceLogger({'run_dir': str(tmp_path)}) as prov:
        prov.log('output.nc', record)

    provenance = yaml.safe_load(
        (tmp_path / 'diagnostic_provenance.yml').read_bytes())

    assert provenance == {'output.nc': record}


def test_provenance_logger_twice(tmp_path):

    record1 = {'attribute1': 'xyz'}
    with shared.ProvenanceLogger({'run_dir': str(tmp_path)}) as prov:
        prov.log('output1.nc', record1)

    record2 = {'attribute2': 'xyz'}
    with shared.ProvenanceLogger({'run_dir': str(tmp_path)}) as prov:
        prov.log('output2.nc', record2)

    provenance = yaml.safe_load(
        (tmp_path / 'diagnostic_provenance.yml').read_bytes())

    assert provenance == {'output1.nc': record1, 'output2.nc': record2}


def test_provenance_logger_duplicate_raises(tmp_path):

    record = {'attribute1': 'xyz'}
    with shared.ProvenanceLogger({'run_dir': str(tmp_path)}) as prov:
        prov.log('output.nc', record)
        with pytest.raises(KeyError):
            prov.log('output.nc', record)


def test_select_metadata():

    metadata = [
        {
            'short_name': 'pr',
            'filename': 'test_pr.nc',
        },
        {
            'short_name': 'ta',
            'filename': 'test_ta.nc',
        },
    ]

    result = shared.select_metadata(metadata, short_name='ta')

    assert result == [{'short_name': 'ta', 'filename': 'test_ta.nc'}]


def test_group_metadata():

    metadata = [
        {
            'short_name': 'pr',
            'filename': 'test_pr.nc',
        },
        {
            'short_name': 'ta',
            'filename': 'test_ta.nc',
        },
    ]

    result = shared.group_metadata(metadata, 'short_name')

    assert result == {
        'ta': [
            {
                'short_name': 'ta',
                'filename': 'test_ta.nc'
            },
        ],
        'pr': [
            {
                'short_name': 'pr',
                'filename': 'test_pr.nc'
            },
        ],
    }


def test_group_metadata_sorted():

    metadata = [
        {
            'short_name': 'ta',
            'dataset': 'dataset2',
        },
        {
            'short_name': 'ta',
            'dataset': 'dataset1',
        },
    ]

    result = shared.group_metadata(metadata, 'short_name', sort='dataset')

    assert result == {
        'ta': [
            {
                'short_name': 'ta',
                'dataset': 'dataset1'
            },
            {
                'short_name': 'ta',
                'dataset': 'dataset2'
            },
        ],
    }


def test_group_metadata_sorted_true():

    metadata = [
        {
            'short_name': 'ta',
        },
        {
            'short_name': 'pr',
        },
    ]

    result = shared.group_metadata(metadata, 'short_name', sort=True)

    assert result == {
        'pr': [
            {
                'short_name': 'pr',
            },
        ],
        'ta': [
            {
                'short_name': 'ta',
            },
        ],
    }


def test_sorted_metadata():

    metadata = [
        {
            'short_name': 'ta',
            'dataset': 'dataset2',
        },
        {
            'short_name': 'pr',
            'dataset': 'dataset2',
            'random_attribute': 1,
        },
        {
            'short_name': 'ta',
            'dataset': 'dataset1',
        },
    ]

    result = shared.sorted_metadata(metadata, sort=['short_name', 'dataset'])

    assert result == [
        {
            'short_name': 'pr',
            'dataset': 'dataset2',
            'random_attribute': 1,
        },
        {
            'short_name': 'ta',
            'dataset': 'dataset1'
        },
        {
            'short_name': 'ta',
            'dataset': 'dataset2'
        },
    ]


@pytest.mark.parametrize('as_iris', [True, False])
def test_extract_variables(as_iris):

    cfg = {
        'input_data': {
            'file1.nc': {
                'short_name': 'ta',
                'standard_name': 'air_temperature',
                'long_name': 'Air Temperature',
                'units': 'K',
            },
            'file2.nc': {
                'short_name': 'ta',
                'standard_name': 'air_temperature',
                'long_name': 'Air Temperature',
            },
            'file3.nc': {
                'short_name': 'pr',
                'standard_name': 'precipitation_flux',
                'long_name': 'Precipitation',
                'extra_attribute': 1,
            },
            'file4.nc': {
                'short_name': 'toz',
                'standard_name': '',
                'long_name': 'Total Ozone Column',
            },
        }
    }

    if as_iris:
        expected = {
            'ta': {
                'var_name': 'ta',
                'standard_name': 'air_temperature',
                'long_name': 'Air Temperature',
                'units': 'K',
            },
            'pr': {
                'var_name': 'pr',
                'standard_name': 'precipitation_flux',
                'long_name': 'Precipitation',
            },
            'toz': {
                'var_name': 'toz',
                'standard_name': None,
                'long_name': 'Total Ozone Column',
            },
        }
    else:
        expected = {
            'ta': {
                'short_name': 'ta',
                'standard_name': 'air_temperature',
                'long_name': 'Air Temperature',
                'units': 'K',
            },
            'pr': {
                'short_name': 'pr',
                'standard_name': 'precipitation_flux',
                'long_name': 'Precipitation',
            },
            'toz': {
                'short_name': 'toz',
                'standard_name': '',
                'long_name': 'Total Ozone Column',
            },
        }

    result = shared.extract_variables(cfg, as_iris)

    assert result == expected


def test_variables_available():

    cfg = {
        'input_data': {
            'file1.nc': {
                'short_name': 'ta'
            },
        }
    }
    assert shared.variables_available(cfg, ['ta']) is True
    assert shared.variables_available(cfg, ['pr']) is False


def test_get_input_data_files(tmp_path):

    metadata1 = {'file1.nc': {'short_name': 'ta', 'dataset': 'dataset1'}}
    metadata_dir1 = tmp_path / 'preproc' / 'ta'
    metadata_dir1.mkdir(parents=True)
    metadata_file1 = metadata_dir1 / 'metadata.yml'
    metadata_file1.write_text(yaml.safe_dump(metadata1))

    metadata2 = {'file2.nc': {'short_name': 'tas', 'dataset': 'dataset1'}}
    metadata_dir2 = tmp_path / 'work_dir'
    metadata_dir2.mkdir()
    metadata_file2 = metadata_dir2 / 'tas_metadata.yml'
    metadata_file2.write_text(yaml.safe_dump(metadata2))

    cfg = {'input_files': [str(metadata_file1), str(metadata_dir2)]}
    input_data = shared._base._get_input_data_files(cfg)

    assert input_data == {
        'file1.nc': metadata1['file1.nc'],
        'file2.nc': metadata2['file2.nc'],
    }


def create_settings(path):

    settings = {
        'log_level': 'debug',
        'work_dir': str(path / 'work_dir'),
        'plot_dir': str(path / 'plot_dir'),
        'run_dir': str(path / 'run_dir'),
        'script': 'diagnostic.py',
        'input_files': [],
        'example_setting': 1,
    }

    return settings


def write_settings(settings):

    run_dir = Path(settings['run_dir'])
    run_dir.mkdir()

    settings_file = run_dir / 'settings.yml'
    settings_file.write_text(yaml.safe_dump(settings))

    return str(settings_file)


def test_run_diagnostic(tmp_path, monkeypatch):

    settings = create_settings(tmp_path)
    settings_file = write_settings(settings)

    monkeypatch.setattr(sys, 'argv', ['', settings_file])

    # Create files created by ESMValCore
    for filename in ('log.txt', 'profile.bin', 'resource_usage.txt'):
        file = Path(settings['run_dir']) / filename
        file.touch()

    with shared.run_diagnostic() as cfg:
        assert 'example_setting' in cfg


@pytest.mark.parametrize('flag', ['-l', '--log-level'])
def test_run_diagnostic_log_level(tmp_path, monkeypatch, flag):
    """Test if setting the log level from the command line works."""
    settings = create_settings(tmp_path)
    settings_file = write_settings(settings)

    monkeypatch.setattr(sys, 'argv', ['', flag, 'error', settings_file])

    with shared.run_diagnostic():
        assert shared._base.logger.getEffectiveLevel() == logging.ERROR


def create_run_content(settings):
    """Create some files to make it look like the diagnostic ran."""
    for dir_name in 'work_dir', 'plot_dir':
        dir_path = Path(settings[dir_name])
        dir_path.mkdir()
        (dir_path / 'example_output.txt').touch()

    for filename in ('log.txt', 'profile.bin', 'diagnostic_provenance.yml',
                     'resource_usage.txt', 'tmp.nc'):
        file = Path(settings['run_dir']) / filename
        file.touch()


def test_rerun_diagnostic_raises(tmp_path, monkeypatch):
    """Test if re-running the diagnostic script fails when output exists."""
    settings = create_settings(tmp_path)
    settings_file = write_settings(settings)

    create_run_content(settings)

    monkeypatch.setattr(sys, 'argv', ['', settings_file])

    with pytest.raises(FileExistsError):
        with shared.run_diagnostic():
            pass


@pytest.mark.parametrize('flag', ['-i', '--ignore', '-f', '--force'])
def test_rerun_diagnostic_flag(tmp_path, monkeypatch, flag):
    """Test if re-running the diagnostic script works."""
    exist = flag in {'-i', '--ignore'}

    settings = create_settings(tmp_path)
    settings_file = write_settings(settings)

    create_run_content(settings)

    monkeypatch.setattr(sys, 'argv', ['', flag, settings_file])

    with shared.run_diagnostic():
        assert not (Path(settings['run_dir']) /
                    'diagnostic_provenance.yml').exists()
        for file in (
                Path(settings['run_dir']) / 'tmp.nc',
                Path(settings['work_dir']) / 'example_output.txt',
                Path(settings['plot_dir']) / 'example_output.txt',
        ):
            assert file.exists() == exist
