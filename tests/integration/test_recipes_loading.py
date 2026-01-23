"""Test recipes are well formed."""

from dataclasses import dataclass
from pathlib import Path

import esmvalcore
import esmvalcore._recipe.check
import esmvalcore._recipe.recipe
import esmvalcore.cmor.check
import esmvalcore.dataset
import esmvalcore.typing
import pytest
import pytest_mock
import yaml
from esmvalcore.config import CFG, Session, _config

try:
    # Since ESValCore v2.14.0
    from esmvalcore.io import local

    USE_DATA_SOURCES = True
except ImportError:
    # Prior to ESMValCore v2.14.0
    from esmvalcore import local

    USE_DATA_SOURCES = False

import esmvaltool


@pytest.fixture
def session(mocker: pytest_mock.MockerFixture, tmp_path: Path) -> Session:
    session = CFG.start_session("test")
    session["auxiliary_data_dir"] = str(tmp_path / "auxiliary_data_dir")
    session["check_level"] = esmvalcore.cmor.check.CheckLevels["DEFAULT"]
    if USE_DATA_SOURCES:
        for project in session["projects"]:
            session["projects"][project]["data"] = {
                "mock_data_source": {
                    "type": "tests.integration.test_recipes_loading.MockDataSource",
                    "priority": 1,
                    "rootpath": tmp_path / project,
                },
            }
    else:
        # Mock input file finding.
        mocker.patch.dict(
            CFG,
            drs={},
            search_esgf="never",
            rootpath={"default": str(tmp_path)},
        )
        mocker.patch.object(
            local,
            "glob",
            autospec=True,
            side_effect=lambda *_, **__: [
                "test_0001-1849.nc",
                "test_1850-9999.nc",
            ],
        )
        # The patched_datafinder fixture does not return the correct input
        # directory structure, so make sure it is set to flat for every project
        for project in _config.CFG:
            mocker.patch.dict(_config.CFG[project]["input_dir"], default="/")

    return session


@dataclass
class MockDataSource:
    name: str
    project: str
    priority: int
    rootpath: Path
    debug_info: str = "Mock data source for testing."

    def find_data(
        self,
        **facets: esmvalcore.typing.FacetValue,  # noqa: ARG002
    ) -> local.LocalFile:
        file1 = local.LocalFile(self.rootpath / "test_0001-1849.nc")
        file1.facets["timerange"] = "0001/1849"
        file1.facets["version"] = "v20230101"
        file2 = local.LocalFile(self.rootpath / "test_1850-9999.nc")
        file2.facets["timerange"] = "1850/9999"
        file2.facets["version"] = "v20230101"
        return [file1, file2]


def _get_recipes():
    recipes_path = Path(esmvaltool.__file__).absolute().parent / "recipes"
    recipes = sorted(recipes_path.glob("**/recipe*.yml"))
    ids = tuple(str(p.relative_to(recipes_path)) for p in recipes)
    return recipes, ids


RECIPES, IDS = _get_recipes()


@pytest.mark.parametrize("recipe_file", RECIPES, ids=IDS)
def test_recipe_valid(recipe_file, session, mocker):
    """Check that recipe files are valid ESMValTool recipes."""
    # Do not remove unexpanded supplementaries. These cannot be expanded
    # because the mocked file finding above does not produce facets.
    mocker.patch.object(
        esmvalcore.dataset.Dataset,
        "_remove_unexpanded_supplementaries",
        autospec=True,
        spec_set=True,
    )

    # Mock vertical levels
    mocker.patch.object(
        esmvalcore._recipe.recipe,
        "get_reference_levels",
        autospec=True,
        spec_set=True,
        side_effect=lambda *_, **__: [1, 2],
    )

    # Mock valid NCL version
    mocker.patch.object(
        esmvalcore._recipe.check,
        "ncl_version",
        autospec=True,
        spec_set=True,
    )

    # Mock interpreters installed
    def which(executable):
        if executable in ("ncl", "python", "Rscript"):
            path = "/path/to/" + executable
        else:
            path = None
        return path

    mocker.patch.object(
        esmvalcore._task,
        "which",
        autospec=True,
        side_effect=which,
    )

    # Create a shapefile for extract_shape preprocessor if needed
    recipe = yaml.safe_load(recipe_file.read_text())
    for preproc in recipe.get("preprocessors", {}).values():
        extract_shape = preproc.get("extract_shape")
        if extract_shape and "shapefile" in extract_shape:
            filename = (
                Path(session["auxiliary_data_dir"])
                / extract_shape["shapefile"]
            )
            filename.parent.mkdir(parents=True, exist_ok=True)
            filename.touch()

    esmvalcore._recipe.recipe.read_recipe_file(recipe_file, session)
