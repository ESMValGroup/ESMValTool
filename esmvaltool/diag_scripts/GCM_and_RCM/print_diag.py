# quick diagnostic for printing output from pre-processor
from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
    group_metadata,
)

import os
import logging

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    # The config object is a dict of all the metadata from the pre-processor
    logger.debug(cfg)

    projects = group_metadata(cfg["input_data"].values(), "project")

    for k, p in projects.items():
        m_list = set()
        for ds in p:
            if k == "CORDEX":
                ds_str = f"{ds['driver']} - {ds['dataset']}"
            else:
                ds_str = ds["dataset"]
            m_list.add(ds_str)
        print(f"{k} - {len(m_list)} models:")
        print(m_list)


if __name__ == "__main__":
    with run_diagnostic() as cfg:
        main(cfg)
