"""Test the esmvaltool.utils.testing.regression.compare module."""
import os

from esmvaltool.utils.testing.regression.compare import (
    get_recipe_name_from_file
)


def test_get_recipe_name_from_file(tmp_path):
    """Test the string extractor for recipe name."""
    path_to_file = tmp_path / "recipe_python_20220317_162441" / "run"
    os.makedirs(path_to_file)
    obtained = get_recipe_name_from_file(path_to_file)
    print("Obtained: %s", obtained)
    print("Expected: name_from_file0/recipe_python")
    assert obtained == "name_from_file0/recipe_python"
