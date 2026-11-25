"""Tests for selecting ENSO teleconnection masks."""

import pytest

from esmvaltool.diag_scripts.enso_metrics import (
    teleconnections_metrics as diag,
)


@pytest.fixture(params=["CMIP6", "CMIP7"])
def model(request):
    return {
        "project": request.param,
        "dataset": "model-a",
        "exp": "historical",
        "ensemble": "r1i1p1f1",
        "variable_group": "tos_enso",
    }


def test_select_model_mask(model):
    mask = {**model, "variable_group": "enso_mask", "filename": "mask.nc"}
    other_inputs = [
        {
            **mask,
            "project": "CMIP7" if model["project"] == "CMIP6" else "CMIP6",
        },
        {**mask, "dataset": "model-b"},
        {**mask, "exp": "piControl"},
        {**mask, "ensemble": "r2i1p1f1"},
        {**mask, "variable_group": "tos_enso"},
    ]

    assert diag.select_model_mask([*other_inputs, mask], model) is mask


@pytest.mark.parametrize("count", [0, 2])
def test_select_model_mask_requires_one_match(model, count):
    mask = {**model, "variable_group": "enso_mask"}
    with pytest.raises(
        ValueError,
        match=f"Expected one ENSO mask for model-a; found {count}",
    ):
        diag.select_model_mask([mask] * count, model)


def test_main_uses_each_model_mask_and_records_provenance(
    model,
    mocker,
    tmp_path,
):
    variables = [
        "tos_enso",
        "pr_global",
        "tos_global",
        "ts_enso_events",
        "enso_mask",
    ]
    input_data = {}
    for dataset in ["model-b", "model-a", "reference"]:
        for variable in variables:
            filename = f"{dataset}_{variable}.nc"
            input_data[filename] = {
                **model,
                "project": "OBS"
                if dataset == "reference"
                else model["project"],
                "dataset": dataset,
                "variable_group": variable,
                "filename": filename,
            }
    mocker.patch.object(diag.iris, "load_cube", side_effect=lambda path: path)
    compute = mocker.patch.object(
        diag,
        "compute_telecon_metrics",
        return_value=({"DJF": 1.0}, {"DJF": None}, {"DJF": None}),
    )
    mocker.patch.object(
        diag,
        "get_diagnostic_filename",
        return_value=str(tmp_path / "metrics.csv"),
    )
    save = mocker.patch.object(diag, "save_figure")

    diag.main({"input_data": input_data})

    assert compute.call_count == 4
    for call in compute.call_args_list:
        input_pair, _, mask, _ = call.args
        dataset = next(
            name for name in input_pair if name.startswith("model-")
        )
        assert mask == f"{dataset}_enso_mask.nc"
    assert save.call_count == 8
    for call in save.call_args_list:
        basename, provenance, _ = call.args
        dataset = basename.split("_")[0]
        assert f"{dataset}_enso_mask.nc" in provenance["ancestors"]
