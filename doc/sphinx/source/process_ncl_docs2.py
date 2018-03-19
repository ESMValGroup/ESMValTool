"""
This script is part of the ESMValTool distribution.  It's been added as part of
the incorporation of the Sphinx documentation generator.  Sphinx was originally
developed for documenting Python code, and one of its features is that it is
able - using the so-called autodoc extension - to extract documentation strings
from Python source files and use them in the documentation it generates.

The autodoc feature apparently does not exist for NCL source files (such as
those which are used in ESMValTool), but it has been mimicked
(or - more-or-less - reverse-engineered) here via this script, which walks
through a subset of the ESMValTool NCL scripts, extracts function names,
argument lists and descriptions (from the comments immediately following the
function definition), and assembles them in a subdirectory of
doc/sphinx/source. These output files are in the so-called reStructuredText
format (see, e.g., http://docutils.sourceforge.net/rst.html), which is the
markup language used by Sphinx; running make in doc/sphinx builds the
ESMValTool documentation from them, as noted above.

Created on July 14, 2015

@author: jeremy.walton@metoffice.gov.uk
"""

import os
import glob
import re
import string
import collections


def make_param_details(params):
    """
    Create a list of parameter names and types from the params string.
    :param params:
    :return:
    """

    # We'll store the parameter names and their types in a dictionary.
    # Note that it has to be an ordered dictionary, because later on we want
    # to pull the entries out in the same order
    # that we added them.
    param_details = collections.OrderedDict()
    for param in params:

        # Extract the type if it's specified,
        # otherwise default to integer (say).
        if ':' in param:
            [pname, ptype] = param.split(':')
        else:
            pname = param
            ptype = 'integer'

        # If the parameter is an array,
        # we only want its name in the description.
        pname = pname.split('[')[0]
        pname = pname.strip()

        # Tie the name and the type of the parameter together.
        param_details[pname] = ptype

    return param_details


def process_params(params, inp, oup):
    """
    Extract the parameter names and types from the params string, pull their
    descriptions out from the input file and reformat the lot in the output.
    """
    # Get the names and types.
    param_details = make_param_details(params)

    # We assume we're at the line before the first parameter description.
    # Bump it, then check to see if we're really at the right location and
    # issue a warning if not.
    line = next(inp)
    param_keys = list(param_details.keys())
    if param_keys[0] not in line:
        print("Warning - parameter " + param_keys[0] +
              " not found in this line:\n" + line)

    # Want a blank line just before parameter descriptions.
    oup.write('\n')

    # Loop over all parameters in the argument list.
    for i, pname in enumerate(param_keys):

        # Now assemble the description from the line(s).
        if pname in line:

            # Get the text in the line which follows the first occurrence
            # (reading from the left) of the parameter name, then strip
            # trailing spaces (including the CR).
            pdesc = line.split(pname, 1)[1]
            pdesc = pdesc.rstrip()

            # The description could continue on the following lines, which
            # need to be concatenated together.  For all except the last
            # parameter, the end of the description is signaled by the name of
            # the next parameter.  For the last (or maybe the  only) parameter,
            #  it's signaled by a blank line.
            line = next(inp)
            if i < len(param_keys)-1:
                pnext = param_keys[i + 1]
                if pnext not in line:
                    # Do the concatenation, stripping whitespace
                    # (including the CR) as we go.
                    while pnext not in line:
                        pdesc += " " + line.replace(';;', '  ', 1).strip()
                        line = next(inp)
            else:
                while not line.replace(';;', '  ', 1).isspace():
                    pdesc += " " + line.replace(';;', '  ', 1).strip()
                    line = next(inp)

            # Ensure the description starts with a colon.
            if pdesc[0] != ':':
                pdesc = ':' + pdesc

            # Write out the complete description of this parameter.
            oup.write('   :param ' + param_details[pname] + ' '
                      + pname + pdesc + '\n')

    # Want a blank line just after parameter descriptions.
    oup.write('\n')


def find_argument(inp):
    """
    Find the start of the Arguments list.
    """

    line = next(inp)
    count = 1
    while 'Arguments' not in line:
        line = next(inp)

        # We assume we're going to find this within two lines of the original
        # location of the input
        # - stop looking if we don't.
        count += 1
        if count > 2:
            return False

    return True


def parse_file(in_filename, out_filename):
    """
    Processes an ncl file and produces an rst file as output, which contains
    documentation of the ncl functions in a form suitable for input to
    the Sphinx documentation generator.
    :param in_filename:
    :param out_filename:
    :return:
    """

    # Open the files.
    try:
        inp = open(in_filename, "r")
    except IOError:
        print("Couldn't open", in_filename)
        return

    try:
        oup = open(out_filename, "w")
    except IOError:
        print("Couldn't open", out_filename)
        return

    # We assume the file name has the form /path/to/foo.ncl, and the
    # module name is foo.  Pull it out, and write it to the output file
    # as the title.
    mod_name = os.path.splitext(os.path.basename(in_filename))[0]

    oup.write(':mod:' + '`' + mod_name + '`' + '\n')
    oup.write("=" * (7+len(mod_name)) + '\n')

    for line in inp:

        # Is this the start of a function?
        if re.match('^function', line) or re.match('^procedure', line):

            # The function could have parameters on the following lines.
            # Concatenate them up until the closing bracket, stripping
            # whitespace (including the CR) as we go.
            fname = line.rstrip()
            while ')' not in fname:
                line = next(inp)
                fname += " " + line.strip()

            # Some ncl files have backslashes in the function declaration to
            # indicate continuation to the next line (even though this isn't
            # necessary in ncl).  These will mess up our processing of
            # the argument list, and don't look good in the doc. so we pull
            # them out here.
            fname = fname.replace('\\', '')

            # Write the line out from the word 'function' onwards, and suitably
            # decorated for rst. Need the CR at the end, as we've been pulling
            # that off throughout the assembly of this line.
            oup.write('.. function:: ' + fname[len('function')+1:] + '\n')

            # Now extract the list of parameters from the function declaration.
            # First, pullout the text between the brackets, then split that
            # into individual parameter names.
            plist = fname.split('(')[1].split(')')[0]
            params = plist.split(',')

            # Position the input just after the line containing 'Arguments'.
            if not find_argument(inp):
                print("Warning - argument list not found for " + fname)
            else:

                # Here's where we check whether this function has any
                # parameters.  If it doesn't, then we don't need to
                # process any.
                if len(plist) > 0:
                    # Read the parameter descriptions and reformat them
                    # before writing them out.
                    process_params(params, inp, oup)

                # We assume the first batch of comments immediately following
                # the function arepart of the documentation.
                line = next(inp)
                while re.match('^;;', line):

                    # Write out this line, replacing the comments with spaces.
                    oup.write(line.replace(';;', '  ', 1))
                    line = next(inp)

    # Close the files.
    inp.close()
    oup.close()


def create_doc_files_from_ncl():
    # Do some rudimentary checking of where this script is being run from,
    # because we're going to be using relative paths below to find the
    # directories containing the input & output.
    file_path = os.path.dirname(os.path.realpath(__file__))
    esmval_root_folder = os.path.abspath(os.path.join(file_path, '..', '..',
                                                      '..'))

    # List the directories containing input files, then loop over them.
    ncl_folders = {'diag_scripts': 'esmvaltool/diag_scripts/lib/ncl',
                   'plot_scripts': 'esmvaltool/plot_scripts/ncl'}
    for ncl_folder in ncl_folders:
        in_dir = os.path.join(esmval_root_folder, ncl_folders[ncl_folder])
        # Form the output directory name from the input directory name
        # (NB we assume the latter are all named ../../../foo/bar, where foo
        # is the useful part of the name.
        out_dir = os.path.join(esmval_root_folder, "doc/sphinx/source/",
                               ncl_folder)
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)

        # Find all the ncl files in the input directory, and loop over them.
        in_files = glob.glob(os.path.join(in_dir, '*.ncl'))
        index_file = open(os.path.join(out_dir, 'index.rst'), 'w')
        write_index_header(index_file, ncl_folder)

        for nclFile in in_files:
            print("Processing " + nclFile)
            rst_filename = os.path.basename(nclFile).replace('.ncl', '.rst')
            rst_file = os.path.join(out_dir, rst_filename)
            parse_file(nclFile, rst_file)
            index_file.write('   ')
            index_file.write(os.path.basename(nclFile).replace('.ncl', ''))
            index_file.write('\n')


def write_index_header(index_file, ncl_folder):
    index_file.write(ncl_folder.upper())
    index_file.write('\n')
    index_file.write('-' * len(ncl_folder))
    index_file.write('\n')
    index_file.write('\n')
    index_file.write('.. toctree::\n   :maxdepth: 2\n\n')


if __name__ == '__main__':
    create_doc_files_from_ncl()
