import copy
import os
import subprocess
import sys
from itertools import product

import yaml


def expand(filename):
    yield filename

    options_file = os.path.join(os.path.dirname(__file__), 'options.yml')
    options_key = os.path.basename(filename)
    with open(options_file) as file:
        options = yaml.safe_load(file)[options_key]

    with open(filename) as file:
        recipe = yaml.safe_load(file)

    for key in options:
        for diag_name, diagnostic in recipe['diagnostics'].items():
            for script_name, script in diagnostic['scripts'].items():
                if key in script:
                    for value in options[key]:
                        if script[key] != value:
                            outrecipe = copy.deepcopy(recipe)
                            outrecipe['diagnostics'][diag_name]['scripts'][
                                script_name][key] = value
                            outfile = os.path.splitext(
                                os.path.basename(filename))[0]
                            outfile = outfile + '_' + str(key) + '_' + str(
                                value).replace(os.sep, '-') + '.yml'
                            print("Creating", outfile)
                            with open(outfile, 'w') as file:
                                yaml.safe_dump(outrecipe, file)
                            yield outfile


def matrix_expand(filename):
    options_file = os.path.join(os.path.dirname(__file__), 'options.yml')
    options_key = os.path.basename(filename)
    with open(options_file) as file:
        options = yaml.safe_load(file)[options_key]

    with open(filename) as file:
        recipe = yaml.safe_load(file)

    for diag_name, diagnostic in recipe['diagnostics'].items():
        for script_name, script in diagnostic['scripts'].items():
            outrecipe = copy.deepcopy(recipe)
            keys = list(options)
            for values in product(*[options[k] for k in keys]):
                outfile = os.path.splitext(os.path.basename(filename))[0]
                for i, key in enumerate(keys):
                    value = values[i]
                    if key in script:
                        outrecipe['diagnostics'][diag_name]['scripts'][
                            script_name][key] = value
                        outfile = outfile + '_' + str(key) + '_' + str(
                            value).replace(os.sep, '-')
                outfile += '.yml'
                print("Creating", outfile)
                with open(outfile, 'w') as file:
                    yaml.safe_dump(outrecipe, file)
                yield outfile


def submit(recipe):
    bsub_job_template = os.path.join(
        os.path.dirname(__file__), 'bsub.template')
    with open(bsub_job_template) as file:
        content = file.read()
    bsub_jobfile = os.path.splitext(os.path.basename(recipe))[0] + '.bsub'

    with open(bsub_jobfile, 'w') as file:
        file.write(content.format(recipe=recipe))

    print("Submitting", bsub_jobfile)
    with open(bsub_jobfile) as stdin:
        subprocess.run('bsub', stdin=stdin, check=True)


def main():
    for filename in sys.argv[1:]:
        for recipe in matrix_expand(filename):
            submit(recipe)


if __name__ == '__main__':
    main()
