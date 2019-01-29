"""Utility script for inspecting ncl color tables."""
import argparse
import glob
import logging
import os

from jinja2 import Template

from esmvaltool.diag_scripts.shared.plot import __file__ as plot_path

logger = logging.getLogger(__name__)

PATH_TO_COLORTABLES = os.path.join(os.path.dirname(plot_path), "rgb")

TEMPDIR = "/tmp"

ncl_script = """
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin
  wks = gsn_open_wks("pdf","{{ outdir }}/test")
  {% for n in list_of_snippets %}{{n}} {% endfor %}
end
"""

color_snippet = """
  cmap := read_colormap_file("{{ path }}")
  opt = True
  opt@Frame = False
  draw_color_palette(wks, cmap, opt)
  gsn_text_ndc(wks, "{{ name }}", 0.5, 0.9, True)
  frame(wks)
"""


def main(args):
    """Main."""
    t_color_snippet = Template(color_snippet)
    t = Template(ncl_script)

    if not os.path.isdir(args.colorpath):
        logger.warning("Path '%s' is invalid" % args.colorpath)
        raise OSError
    list_of_snippets = []
    for path in glob.glob(args.colorpath + "/ipcc-ar6*rgb"):
        head, tail = os.path.split(path)
        list_of_snippets.append(t_color_snippet.render(path=path, name=tail))
    filename = os.path.join(TEMPDIR, "show_color_tables.ncl")
    with open(filename, "w") as f:
        f.write(
            t.render(
                list_of_snippets=sorted(list_of_snippets),
                outdir=args.outpath))
    if os.system("ncl " + filename) != 0:
        logger.warning("External command failed")


def get_args():
    """Define the commandline arguments."""
    parser = argparse.ArgumentParser(
        description="Utility script for inspecting ncl color tables.")
    parser.add_argument(
        '-c',
        '--colorpath',
        metavar='COLOR_TABLE_DIR',
        default=PATH_TO_COLORTABLES,
        help="Set directory to <COLOR_TABLE_DIR>.")
    parser.add_argument(
        '-o',
        '--outpath',
        metavar='OUTDIR',
        default=TEMPDIR,
        help="Set out directory to <OUTDIR>.")

    return parser.parse_args()


def run():
    main(get_args())


if __name__ == '__main__':
    run()
