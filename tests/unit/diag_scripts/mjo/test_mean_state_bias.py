# Copyright (C) 2026 ESMValTool development team
"""Tests for dataset-independent mean-state bias labels."""

from esmvaltool.diag_scripts.mjo import mean_state_bias


class FakeCube:
    """Minimal cube for testing label propagation."""

    def __init__(self):
        self.long_name = None
        self.var_name = None

    def collapsed(self, *_args):
        return self

    def __sub__(self, _other):
        return FakeCube()

    def rename(self, name):
        self.long_name = name


def test_bias_names_follow_recipe_aliases(monkeypatch):
    """Changing input aliases updates the saved bias descriptions."""
    aliases = {
        "obs_pr": "Reference-A",
        "obs_ua": "Reference-A",
        "obs_va": "Reference-A",
        "historical_pr": "Model-B_historical",
        "control_pr": "Model-C_piControl",
    }
    input_data = {
        name: {
            "variable_group": name,
            "filename": f"{name}.nc",
            "dataset": "unused-default",
            "alias": alias,
        }
        for name, alias in aliases.items()
    }
    monkeypatch.setattr(mean_state_bias.iris, "load_cube", lambda _: FakeCube())
    monkeypatch.setattr(mean_state_bias.iris.util, "squeeze", lambda cube: cube)
    monkeypatch.setattr(
        mean_state_bias,
        "regrid",
        lambda *_args, **_kwargs: FakeCube(),
    )

    cubes, input_aliases, _ = mean_state_bias._load_inputs(
        {"input_data": input_data},
    )
    _, historical_bias, control_bias = mean_state_bias.calculate_biases(
        cubes,
        input_aliases,
    )

    assert input_aliases == aliases
    assert historical_bias.long_name == (
        "Model-B_historical minus Reference-A precipitation"
    )
    assert control_bias.long_name == (
        "Model-C_piControl minus Reference-A precipitation"
    )
