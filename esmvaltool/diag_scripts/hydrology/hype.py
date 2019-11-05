"""PCR-GLOBWB diagnostic."""
import logging
from pathlib import Path

import dask.array as da
import iris

import numpy
# daskframe is possible if
# pip install 'fsspec>=0.3.3'
import pandas


from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata, run_diagnostic)

logger = logging.getLogger(Path(__file__).name)


def get_provenance_record(attributes):
    """Create a provenance record."""

    ancestor_file = attributes['filename']

    record = {
        'caption': "Forcings for the Hype hydrological model.",
        'domains': ['global'],
        'authors': [
            'pelupessy_inti',
        ],
        'projects': [
            'ewatercycle',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': [ancestor_file],
    }
    return record


def get_output_stem(attributes):
    
    short_to_stem=dict(tas="Tobs",
                       tasmin="TMINobs",
                       tasmax="TMAXobs",
                       pr="Pobs")
    
    shortname=attributes["short_name"]
    if shortname in short_to_stem:
        stem = short_to_stem[shortname]
    else:
        stem = Path(attributes['filename']).stem + '_hype'
        
    return stem


def get_data_times_and_ids(attributes):
    input_file = attributes['filename']

    cube = iris.load_cube(input_file)

    data = numpy.array(cube.core_data()).T

    print("units:", attributes["units"])

    # ad-hoc fix of precipitation:
    if attributes["short_name"] == "pr":
        if attributes["units"]=="kg m-2 s-1":
            data*=86400
        elif not attributes["units"]=="mm day-1":
            raise Exception("possible units error")
              

    # Round times to integer number of days
    time_coord = cube.coord('time')
    time_coord.points = da.floor(time_coord.core_points())
    time_coord.bounds = None

    times = [str(x.point.date()) for x in time_coord.cells()]
    ids = cube.coord('shape_id').core_points()

    return data,times,ids


def main(cfg):
    """Process data for use as input to the HYPE hydrological model."""
    input_data = cfg['input_data'].values()
    grouped_input_data = group_metadata(input_data,
                                        'standard_name',
                                        sort='dataset')

    for standard_name in grouped_input_data:
        logger.info("Processing variable %s", standard_name)
        for attributes in grouped_input_data[standard_name]:
            logger.info("Processing dataset %s", attributes['dataset'])

            output_file = get_diagnostic_filename(get_output_stem(attributes),
                                                  cfg, 'txt')

            data, times, ids = get_data_times_and_ids(attributes)

            frame = pandas.DataFrame(data, index=times, columns=ids)

            frame.to_csv(output_file, sep=' ', index_label="DATE",
                         float_format='%.3f')

            # Store provenance
            provenance_record = get_provenance_record(attributes)
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(output_file, provenance_record)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
