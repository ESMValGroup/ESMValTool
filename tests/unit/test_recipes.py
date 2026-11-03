"""Recipe tests."""

import textwrap
from functools import cache
from pathlib import Path

import pytest
import yaml

import esmvaltool
from tests.integration.test_recipes_loading import IDS, RECIPES

ESMVALTOOL_ROOT = Path(esmvaltool.__file__).absolute().parent
REFERENCES_PATH = ESMVALTOOL_ROOT / "references"

CONFIG_REFERENCES_PATH = ESMVALTOOL_ROOT / "config-references.yml"
AUTHORS = yaml.safe_load(CONFIG_REFERENCES_PATH.read_text())["authors"]


@cache
def load_recipe(recipe_file: Path) -> dict:
    """Load a recipe only once and use the fast C-based loader."""
    return yaml.load(
        recipe_file.read_text(encoding="utf-8"),
        Loader=yaml.CSafeLoader,
    )


@pytest.mark.parametrize("recipe_file", RECIPES, ids=IDS)
def test_reference_tags(recipe_file):
    """Check bibtex file is added to REFERENCES_PATH."""
    recipe = load_recipe(recipe_file)
    tags = recipe.get("documentation", {}).get("references", [])
    msg = textwrap.dedent("""
        The tag '{tag}' is mentioned in recipe '{recipe}'.
        However, its reference file '{tag}.bibtex' is not available in {path}.
        Please check instructions on how to add references at
        https://docs.esmvaltool.org/en/latest/community/diagnostic.html#adding-references
    """)
    for tag in tags:
        bibtex_file = REFERENCES_PATH / f"{tag}.bibtex"
        assert bibtex_file.is_file(), msg.format(
            tag=tag,
            recipe=recipe_file,
            path=REFERENCES_PATH,
        )


@pytest.mark.parametrize("recipe_file", RECIPES, ids=IDS)
def test_maintainers(recipe_file):
    """Check recipe maintainers."""
    recipe = load_recipe(recipe_file)

    # Make sure that 'documentation' and 'maintainer' entries are present
    msg = "'documentation' missing in recipe"
    assert "documentation" in recipe, msg
    msg = "'maintainer' missing in 'documentation' of recipe"
    assert "maintainer" in recipe["documentation"], msg
    maintainers = recipe["documentation"]["maintainer"]

    # Make sure that maintainer entry is list and non-empty
    msg = f"'maintainer' entry needs to be a list, got {type(maintainers)}"
    assert isinstance(maintainers, list), msg
    msg = "'maintainer' needs to contain at least one element"
    assert maintainers, msg

    # Make sure that 'unmaintained' is single entry if used
    if "unmaintained" in maintainers:
        msg = (
            f"If 'unmaintained' is given as maintainer, this has to be "
            f"the sole entry, got {maintainers}"
        )
        assert len(maintainers) == 1, msg

    # Check that maintainers are valid
    invalid_maintainers = [m for m in maintainers if m not in AUTHORS]
    msg = (
        f"Got invalid maintainers: {invalid_maintainers}. Valid entries "
        f"are authors from {CONFIG_REFERENCES_PATH}."
    )
    assert not invalid_maintainers, msg
