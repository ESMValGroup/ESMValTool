"""Utility script for inspecting and converting ncl color tables."""
import argparse
import glob
import logging
import os
import subprocess
import tempfile

import matplotlib
matplotlib.use("Agg")  # noqa
import matplotlib.pyplot as plt
import numpy as np
import yaml
from jinja2 import Template

from esmvaltool.diag_scripts.shared.plot import __file__ as plot_path

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
CONSOLE_HANDLER = logging.StreamHandler()
CONSOLE_HANDLER.setLevel(logging.INFO)
CONSOLE_HANDLER.setFormatter(
    logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(CONSOLE_HANDLER)

PATH_TO_COLORTABLES = os.path.join(os.path.dirname(plot_path), "rgb")

NCL_SCRIPT = """
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin
  print("Autogenerated NCL script is running.")
  wks = gsn_open_wks("pdf","{{ outdir }}/available_colormaps_for_ncl")
  {% for n in list_of_snippets %}{{n}} {% endfor %}
end
"""

COLOR_SNIPPET = """
  cmap := read_colormap_file("{{ path }}")
  opt = True
  opt@Frame = False
  draw_color_palette(wks, cmap, opt)
  gsn_text_ndc(wks, "{{ name }}", 0.5, 0.9, True)
  frame(wks)
"""


def load_ncl_color_map(name, colorpath):
    """Load ncl color map to a list that is returned."""
    def _format(content):
        out = []
        for item in content.split("\n"):
            item = item.strip()
            if item and not ('ncolors' in item or item.startswith('#')
                             or item.startswith(';')):
                out.append([int(elem) / 256
                            for elem in item.split()[0:3]] + [1])
        return out

    filename = "{0}/{1}.rgb".format(colorpath, name)
    if not os.path.exists(filename):
        raise ValueError("Path {0} does not exist.".format(filename))
    with open(filename, 'r') as ncl_color_map:
        return _format(ncl_color_map.read())


def get_color_map(name, colorpath):
    """Convert colormap from ncl to python.

    Parameters
    ----------
    name: str
        Name of ncl color map
    Returns
    -------
    matplotlib.colors.ListedColorMap object
    """
    colors = load_ncl_color_map(name, colorpath)
    logger.debug("RGB values for '%s':\n%s", name, yaml.dump(colors))
    return matplotlib.colors.ListedColormap(colors, name=name, N=None)


def list_ncl_color_maps(colorpath):
    """Get list of all available ncl color maps."""
    from os import walk

    def _format(name):
        return os.path.splitext(os.path.basename(name))[0]

    out = []
    for (_, _, filenames) in walk(colorpath):
        out.extend([
            _format(filename) for filename in filenames
            if 'rgb' in filename.split('.')
        ])
    return out


def plot_example_for_colormap(name, colorpath, outdir='./'):
    """Create plots of given color map using python."""
    logger.info("Plotting example for '%s'", name)
    fig = plt.figure(1)
    axis = fig.add_axes([0.1, 0.3, 0.5, 0.5])
    np.random.seed(12345678)
    data = np.random.randn(30, 30)
    psm = axis.pcolormesh(
        data,
        cmap=get_color_map(name, colorpath),
        rasterized=True,
        vmin=-4,
        vmax=4)
    fig.colorbar(psm, ax=axis)
    plt.savefig(os.path.join(outdir, "{0}.png".format(name)))
    plt.close()


def main_plot_python_cm(colorpath, outpath):
    """Execute functions for python plots."""
    for name in list_ncl_color_maps(colorpath):
        plot_example_for_colormap(name, colorpath, outdir=outpath)


def main_plot_ncl_cm(colorpath, outpath):
    """Execute functions for ncl plots."""
    t_color_snippet = Template(COLOR_SNIPPET)
    template = Template(NCL_SCRIPT)
    list_of_snippets = []
    for path in glob.glob(colorpath + "/*rgb"):
        _, tail = os.path.split(path)
        list_of_snippets.append(t_color_snippet.render(path=path, name=tail))
    with tempfile.NamedTemporaryFile(mode='w', suffix='ncl') as fname:
        fname.write(
            template.render(
                list_of_snippets=sorted(list_of_snippets), outdir=outpath))
        subprocess.check_call(["ncl", fname.name])


def get_args():
    """Define the commandline arguments."""
    parser = argparse.ArgumentParser(description="""
        Utility module for inspecting and converting
        ncl color tables.""")
    parser.add_argument(
        '-c',
        '--colorpath',
        metavar='COLOR_TABLE_DIR',
        default=PATH_TO_COLORTABLES,
        help="Set directory to <COLOR_TABLE_DIR>. Default is {0}".format(
            PATH_TO_COLORTABLES))
    parser.add_argument(
        '-n',
        '--ncl',
        dest='n',
        action='store_true',
        help="""
            Create report of all ncl color maps in
            <COLOR_TABLE_DIR> using ncl.""")
    parser.add_argument(
        '-p',
        '--python',
        dest='p',
        action='store_true',
        help="""Create example plots of all ncl color maps
        in <COLOR_TABLE_DIR> using python.""")
    parser.add_argument(
        '-o',
        '--outpath',
        metavar='OUTDIR',
        default='./',
        help="Set out directory to <OUTDIR>. Default is the current directory")
    return parser.parse_args()


def main(args):
    """Call functions according to command line arguments."""
    colorpath = args.colorpath
    outpath = args.outpath

    if not os.path.isdir(colorpath):
        logger.warning("Path '%s' is invalid", colorpath)
        raise OSError

    if not os.path.isdir(outpath) and not os.path.exists(outpath):
        logger.info("Creating directory '%s'", outpath)
        os.mkdir(outpath)

    if args.n:
        logger.info("Creating report with ncl")
        main_plot_ncl_cm(colorpath, outpath)
    if args.p:
        logger.info("Creating report with python")
        main_plot_python_cm(colorpath, outpath)


def run():
    """Run the program."""
    args = get_args()
    main(args)


if __name__ == '__main__':
    run()
