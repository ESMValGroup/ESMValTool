"""Tests for the impact/bias_and_change.py diagnostic."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from esmvaltool.diag_scripts.impact.bias_and_change import (
    DEFAULT_ALIAS_FACETS,
    METRIC_LABELS,
    VARIABLE_LABELS,
    build_vegalite_spec,
    make_wide,
    save_vegalite_specs,
)

CORDEX_ALIAS_FACETS = {
    "dataset": "model",
    "rcm_version": "version_realization",
    "driver": "driver",
    "ensemble": "member",
}


def make_tidy(aliases: list[str]) -> pd.DataFrame:
    """Create a tidy dataframe as produced by the diagnostic."""
    records = [
        {
            "dataset": alias,
            "variable": variable,
            METRIC_LABELS["bias"]: 1.0 + i,
            METRIC_LABELS["change"]: 2.0 + i,
        }
        for i, alias in enumerate(aliases)
        for variable in VARIABLE_LABELS.values()
    ]
    return pd.DataFrame.from_records(records).set_index(
        ["dataset", "variable"],
    )


@pytest.fixture
def cmip_df() -> pd.DataFrame:
    facets = {
        "CMIP5_ACCESS1-0_r1i1p1": {
            "project": "CMIP5",
            "model": "ACCESS1-0",
            "member": "r1i1p1",
        },
        "CMIP6_MIROC6_r1i1p1f1": {
            "project": "CMIP6",
            "model": "MIROC6",
            "member": "r1i1p1f1",
        },
        "CMIP6_MIROC6_r2i1p1f1": {
            "project": "CMIP6",
            "model": "MIROC6",
            "member": "r2i1p1f1",
        },
    }
    return make_wide(make_tidy(list(facets)), facets)


@pytest.fixture
def cordex_df() -> pd.DataFrame:
    facets = {
        f"ICON-CLM-202407-1-1_v1-r1_{driver}_r1i1p1f1": {
            "model": "ICON-CLM-202407-1-1",
            "version_realization": "v1-r1",
            "driver": driver,
            "member": "r1i1p1f1",
            # Added by load_data even if not in alias_facets
            "project": "CORDEX-CMIP6",
        }
        for driver in ("MIROC6", "MPI-ESM1-2-HR")
    }
    return make_wide(make_tidy(list(facets)), facets)


def test_make_wide(cmip_df):
    assert cmip_df.index.name == "dataset"
    assert list(cmip_df.columns) == [
        "tas_bias",
        "pr_bias",
        "tas_change",
        "pr_change",
        "project",
        "model",
        "member",
    ]
    # Precipitation is converted from kg/m2/s to mm/day
    assert cmip_df["pr_bias"].iloc[0] == 86400.0
    assert cmip_df["tas_bias"].iloc[0] == 1.0


def test_build_vegalite_spec_cmip(cmip_df):
    project_df = cmip_df[cmip_df["project"] == "CMIP6"]
    notes = ["A note."]
    spec = build_vegalite_spec(
        project_df,
        "CMIP6",
        DEFAULT_ALIAS_FACETS,
        notes,
        "recipe_impact_20260101_120000",
    )

    assert spec["usermeta"] == {
        "project": "CMIP6",
        "notes": notes,
        "recipe_output": "recipe_impact_20260101_120000",
    }
    assert [p["name"] for p in spec["params"]] == ["brush", "panzoom", "query"]
    query = spec["params"][2]
    assert query["select"]["fields"] == ["dataset"]

    values = spec["data"]["values"]
    assert [v["dataset"] for v in values] == [
        "CMIP6_MIROC6_r1i1p1f1",
        "CMIP6_MIROC6_r2i1p1f1",
    ]
    assert values[0]["model"] == "MIROC6"

    tas, pr = spec["hconcat"]
    assert tas["encoding"]["x"]["field"] == "tas_bias"
    assert tas["encoding"]["y"]["field"] == "tas_change"
    assert pr["encoding"]["x"]["field"] == "pr_bias"
    assert tas["encoding"]["fill"]["condition"]["field"] == "model"
    tooltip_fields = [t["field"] for t in tas["encoding"]["tooltip"]]
    assert tooltip_fields[:4] == ["dataset", "project", "model", "member"]

    # The specification can be serialized to JSON
    json.dumps(spec)


def test_build_vegalite_spec_cordex(cordex_df):
    spec = build_vegalite_spec(
        cordex_df,
        "CORDEX-CMIP6",
        CORDEX_ALIAS_FACETS,
        [],
        "recipe_impact_cordex-cmip6_20260101_120000",
    )
    values = spec["data"]["values"]
    # Runs of the same model are distinguished by their dataset identifier
    assert len({v["dataset"] for v in values}) == 2
    assert {v["model"] for v in values} == {"ICON-CLM-202407-1-1"}
    tooltip_fields = [
        t["field"] for t in spec["hconcat"][0]["encoding"]["tooltip"]
    ]
    assert "driver" in tooltip_fields
    assert "project" in tooltip_fields


def test_build_vegalite_spec_axis_titles(cmip_df):
    project_df = cmip_df[cmip_df["project"] == "CMIP6"]
    spec = build_vegalite_spec(
        project_df,
        "CMIP6",
        DEFAULT_ALIAS_FACETS,
        [],
        "recipe_impact_20260101_120000",
    )
    tas = spec["hconcat"][0]
    assert tas["encoding"]["x"]["title"] == METRIC_LABELS["bias"]
    assert tas["encoding"]["y"]["title"] == METRIC_LABELS["change"]

    spec = build_vegalite_spec(
        project_df,
        "CMIP6",
        DEFAULT_ALIAS_FACETS,
        [],
        "recipe_impact_20260101_120000",
        axis_titles={"bias": "Bias with respect to ERA5"},
    )
    for view in spec["hconcat"]:
        assert view["encoding"]["x"]["title"] == "Bias with respect to ERA5"
        assert view["encoding"]["y"]["title"] == METRIC_LABELS["change"]


def test_build_vegalite_spec_nan(cmip_df):
    cmip_df.loc["CMIP5_ACCESS1-0_r1i1p1", "pr_change"] = np.nan
    project_df = cmip_df[cmip_df["project"] == "CMIP5"]
    spec = build_vegalite_spec(
        project_df,
        "CMIP5",
        DEFAULT_ALIAS_FACETS,
        [],
        "recipe_impact_20260101_120000",
    )
    assert spec["data"]["values"][0]["pr_change"] is None
    # NaN is not valid JSON
    json.dumps(spec, allow_nan=False)


def test_save_vegalite_specs(mocker, tmp_path, cmip_df):
    def get_diagnostic_filename(basename, _cfg, extension):
        return str(tmp_path / f"{basename}.{extension}")

    module = "esmvaltool.diag_scripts.impact.bias_and_change"
    mocker.patch(
        f"{module}.get_diagnostic_filename",
        side_effect=get_diagnostic_filename,
    )
    log_provenance = mocker.patch(f"{module}.log_provenance")

    cfg = {
        "work_dir": str(
            tmp_path
            / "recipe_impact_20260101_120000"
            / "work"
            / "bias_and_change"
            / "visualize",
        ),
    }
    save_vegalite_specs(cmip_df, DEFAULT_ALIAS_FACETS, ["note"], [], cfg)

    assert sorted(p.name for p in tmp_path.glob("*.json")) == [
        "vegalite_spec_CMIP5.json",
        "vegalite_spec_CMIP6.json",
    ]
    spec = json.loads((tmp_path / "vegalite_spec_CMIP6.json").read_text())
    assert spec["usermeta"] == {
        "project": "CMIP6",
        "notes": ["note"],
        "recipe_output": "recipe_impact_20260101_120000",
    }
    assert len(spec["data"]["values"]) == 2
    assert log_provenance.call_count == 2
