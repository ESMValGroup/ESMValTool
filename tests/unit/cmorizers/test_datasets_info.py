import datetime
from pathlib import Path

import pytest
import yamale
import yaml

import esmvaltool
import esmvaltool.cmorizers.data.formatters.datasets
from esmvaltool.cmorizers.data.cmorizer import DATASETS_FILE
from esmvaltool.cmorizers.data.utilities import read_cmor_config


def test_only_datasets_are_present():
    recipe = yamale.make_data(DATASETS_FILE)
    schema = yamale.make_schema(
        DATASETS_FILE.parent / "datasets_schema.yml",
    )
    yamale.validate(schema, recipe)


def test_latest_version_format():
    cfg = yaml.safe_load(DATASETS_FILE.read_text(encoding="utf-8"))
    for dataset_info in cfg["datasets"].values():
        datetime.datetime.strptime(
            str(dataset_info["last_access"]),
            "%Y-%m-%d",
        )


FORMATTERS_DIR = Path(
    esmvaltool.cmorizers.data.formatters.datasets.__file__,
).parent


@pytest.mark.parametrize(
    "dataset",
    [
        dataset
        for dataset in yaml.safe_load(
            Path(DATASETS_FILE).read_text(encoding="utf-8"),
        )["datasets"]
        # Only Python CMORizers use the configuration file.
        if (
            FORMATTERS_DIR / f"{dataset.lower().replace('-', '_')}.py"
        ).is_file()
    ],
)
def test_datasets_have_valid_config(dataset: str) -> None:
    """Test that all datasets have a valid type."""
    # This list should be synchronized with the valid values mentioned in
    # doc/sphinx/source/input.rst are used.
    valid_types = {
        "sat",
        "ground",
        "clim",
        "reanaly",
        "campaign",
    }
    valid_projects = {
        "OBS",
        "OBS6",
    }

    cfg = read_cmor_config(dataset)
    assert "attributes" in cfg
    attributes = cfg["attributes"]
    assert "type" in attributes
    assert attributes["type"] in valid_types
    assert "project_id" in attributes
    assert attributes["project_id"] in valid_projects


def test_datasets_are_added_to_test_recipe():
    with open(DATASETS_FILE) as file:
        cfg = yaml.safe_load(file)

    recipes_folder = Path(esmvaltool.__file__).parent / "recipes"
    recipe_path = recipes_folder / "examples" / "recipe_check_obs.yml"
    recipe = yaml.safe_load(recipe_path.read_text(encoding="utf-8"))

    tested_datasets = set()
    for diagnostic in recipe.get("diagnostics", {}).values():
        for dataset in diagnostic.get("additional_datasets", {}):
            tested_datasets.add(dataset["dataset"])
        for variable in diagnostic.get("variables", {}).values():
            if variable is None:
                continue
            for dataset in variable.get("additional_datasets", {}):
                tested_datasets.add(dataset["dataset"])

    info_datasets = set(cfg["datasets"].keys())

    if tested_datasets.symmetric_difference(info_datasets):
        for dataset in tested_datasets - info_datasets:
            print(f"Dataset {dataset} missing from datasets.yml")
        for dataset in info_datasets - tested_datasets:
            print(f"Dataset {dataset} missing from recipe_check_obs.yml")
        assert False
