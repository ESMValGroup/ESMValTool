"""Provenance checks for Seaborn plots and plotted data."""

import pandas as pd
import pytest

from esmvaltool.diag_scripts import seaborn_diag


@pytest.mark.parametrize("override", [False, True])
def test_plot_provenance_ancestors(monkeypatch, tmp_path, override):
    """Use configured ancestors, retaining the existing default otherwise."""
    recorded = []

    class Recorder:
        def __init__(self, _cfg):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            pass

        def log(self, path, record):
            recorded.append((path, record["ancestors"]))

    monkeypatch.setattr(seaborn_diag, "ProvenanceLogger", Recorder)
    monkeypatch.setattr(seaborn_diag, "_save_nc_data", lambda *_: None)
    monkeypatch.setattr(seaborn_diag.plt, "savefig", lambda *_a, **_k: None)
    monkeypatch.setattr(seaborn_diag.plt, "close", lambda: None)

    cfg = {
        "seaborn_kwargs": {},
        "seaborn_func": "barplot",
        "plot_object_methods": {},
        "suptitle": None,
        "legend_title": None,
        "write_netcdf": True,
        "plot_filename": "test",
        "savefig_kwargs": {},
        "input_data": {"source": {"filename": "plot_input.nc"}},
        "plot_dir": str(tmp_path),
        "work_dir": str(tmp_path),
        "output_file_type": "png",
    }
    if override:
        cfg["provenance_ancestors"] = ["original_preprocessed.nc"]
    seaborn_diag._create_plot(lambda **_kwargs: object(), pd.DataFrame(), cfg)

    expected = ["original_preprocessed.nc" if override else "plot_input.nc"]
    assert len(recorded) == 2
    assert all(ancestors == expected for _, ancestors in recorded)
