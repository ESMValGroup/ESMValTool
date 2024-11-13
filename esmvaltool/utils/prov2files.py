"""Print out the input files used to generate a result."""
import argparse

from prov.model import ProvDerivation, ProvDocument


def prov2files(filename):
    """Figure out what file was generated from which source files.

    Parameters
    ----------
    filename: str
        Name of the file containing the provenance.

    Returns
    -------
    (str, list[str])
        A tuple, the first entry is the name of the result
        and the second entry a list of files used to compute
        that result.
    """
    provenance = ProvDocument.deserialize(filename, format='xml')

    source_files = set()
    generated_files = set()
    for rec in provenance.get_records(ProvDerivation):
        # Find all derivation relations
        generated, used = rec.args[:2]
        source_files.add(used.localpart)
        generated_files.add(generated.localpart)

    # Filter out intermediate files
    intermediate_files = source_files & generated_files
    source_files = source_files - intermediate_files
    result_files = generated_files - intermediate_files

    if not len(result_files) == 1:
        # If this changes, need to rewrite this function so it
        # builds a provenance graph.
        raise ValueError("Invalid provenance file encountered,"
                         " ESMValTool provenance describes one result only.")
    return result_files.pop(), sorted(source_files)


def main():
    """Print out a list of files."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        'provenance_files',
        nargs='+',
        type=str,
        help='Path to one or more files containing provenance.')
    args = parser.parse_args()

    for filename in args.provenance_files:
        if not filename.endswith('_provenance.xml'):
            print("Skipping", filename,
                  "does it contain ESMValTool provenance?")
            continue
        result, files = prov2files(filename)
        print(f"{result} was derived from:")
        print('\n'.join(files))
        print('')


if __name__ == '__main__':
    main()
