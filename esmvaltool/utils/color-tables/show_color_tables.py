import glob
import os

from jinja2 import Template

PATH_TO_COLORTABLES = "../../../esmvaltool/diag_scripts/shared/plot/rgb"

ncl_script = """
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin
  wks = gsn_open_wks("pdf","test")
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
t_color_snippet = Template(color_snippet)
t = Template(ncl_script)

list_of_snippets = []
for path in glob.glob(PATH_TO_COLORTABLES + "/*rgb"):
    head, tail = os.path.split(path)
    list_of_snippets.append(t_color_snippet.render(path=path, name=tail))
with open("show_color_tables_new.ncl", "w") as f:
    f.write(t.render(list_of_snippets=list_of_snippets))
