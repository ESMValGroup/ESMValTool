from __future__ import annotations

import xarray as xr

from rf_shap.cube.tabular import SampleTable, table_to_cube


def reconstruct_prediction(pred, table: SampleTable, name: str = "prediction") -> xr.DataArray:
    template = table.extras.get("template")
    if template is None:
        raise ValueError("SampleTable is missing a cube template; cannot restore prediction")
    return table_to_cube(pred, table.index, template, name=name)
