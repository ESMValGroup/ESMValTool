"""Calculate AR6 regional trends without expanding every region at once.

The standard decomposed ``extract_shape`` operation creates one full-sized
masked cube per region. For high-resolution monthly data, that intermediate
can exceed memory before the following area mean reduces it. This diagnostic
performs the same preprocessor operations one region at a time.
"""

import logging
from pathlib import Path

import iris
from esmvalcore import preprocessor as pp
from iris.cube import CubeList

from esmvaltool.diag_scripts import seaborn_diag
from esmvaltool.diag_scripts.shared import ProvenanceLogger, run_diagnostic

logger = logging.getLogger(__name__)


def _regional_trend(cube, regions, units):
    """Apply the recipe's original preprocessor sequence by region."""
    means = CubeList()
    for region in regions:
        logger.info("Computing area mean for AR6 region %s", region)
        selected = pp.extract_shape(
            cube,
            "ar6",
            decomposed=True,
            ids={"Acronym": [region]},
        )
        means.append(pp.area_statistics(selected, operator="mean"))

    regional_means = means.merge_cube()
    trend = pp.linear_trend(regional_means)
    return pp.convert_units(trend, units)


def main(cfg):
    """Calculate trend input files and pass them to the existing plotter."""
    input_data = {}
    original_files = []
    for index, metadata in enumerate(cfg["input_data"].values()):
        source = metadata["filename"]
        original_files.append(source)
        logger.info("Calculating regional trends for %s", source)
        cube = iris.load_cube(source)
        trend = _regional_trend(cube, cfg["regions"], cfg["trend_units"])
        output = Path(cfg["work_dir"]) / f"regional_trend_{index}.nc"
        iris.save(trend, output)
        input_data[str(output)] = {**metadata, "filename": str(output)}
        provenance = {
            "ancestors": [source],
            "authors": ["schlund_manuel", "weigel_katja"],
            "caption": "Linear trend of area means for 46 IPCC AR6 land regions",
        }
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(str(output), provenance)

    cfg["input_data"] = input_data
    cfg["provenance_ancestors"] = original_files
    seaborn_diag.main(cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
