# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 13:02:38 2016

@author: bmueller
"""

import subprocess
import sys
from optparse import OptionParser
import datetime
import os
import xml.sax
from flask import Flask, render_template

main_ESMValTool = "../"

# path to ESMVal python toolbox library
sys.path.append(main_ESMValTool+"/interface_scripts")
sys.path.append(main_ESMValTool+"/diag_scripts/lib/python")
sys.path.append(main_ESMValTool+"/diag_scripts")

# path to ESMVal REPORT python toolbox library
sys.path.append("./lib")
sys.path.append("./utils")

import xml_parsers
import projects
from html_writer import HTML_writer
from report_cleanup import cleanup_report

# additional info
version = "0.1"
authors = ["B. Mueller (LMU)", "b.mueller@iggf.geo.uni.muenchen.de"]

# Check command arguments.
usage = "%prog nml/namelist-file.xml"
description = """Reporting service for the ESMValTool \
                    - Earth System Model Evaluation Tool."""

parser = OptionParser(usage=usage, description=description)
options, args = parser.parse_args()  # args equals namelist
if len(args) == 0:
    parser.print_help()
    sys.exit(0)

# Get command arguments.
input_xml_full_path = args[0]

# Parse input namelist into project_info-dictionary.
Project = xml_parsers.namelistHandler()
parser = xml.sax.make_parser()
parser.setContentHandler(Project)

# Current working directory adjustments depending on namelist
CWD = os.getcwd()
os.chdir(main_ESMValTool)
if input_xml_full_path.split(os.sep)[0] == \
    (main_ESMValTool[:-1] if main_ESMValTool[-1] == os.sep
     else main_ESMValTool):
    input_xml_full_path = input_xml_full_path[1:]
else:
    input_xml_full_path = CWD + input_xml_full_path[2:]

parser.parse(input_xml_full_path)
os.chdir(CWD)

# Project_info is a dictionary with all info from the namelist.
project_info = Project.project_info

# Additional entries to 'project_info'. The 'project_info' construct
# is one way by which Python passes on information to the NCL-routines.
project_info['RUNTIME'] = {}

# Input xml path/file
project_info['RUNTIME']['xml'] = input_xml_full_path
input_xml_file = os.path.basename(input_xml_full_path)
project_info['RUNTIME']['xml_name'] = input_xml_file

# Master references-acknowledgements file (hard coded)
in_refs = os.path.join(os.getcwd(), 'doc/MASTER_authors-refs-acknow.txt')
project_info['RUNTIME']['in_refs'] = in_refs

# Current working directory
project_info['RUNTIME']['cwd'] = CWD

# cleanup directories
do_print = False
template_dir = './templates/'
image_dir = './static/images/'
cleanup_report(template_dir, image_dir, do_print)

# try to process the data
try:
    timestamp = datetime.datetime.utcnow()
    ESMValTool_log = subprocess.check_output("cd " + main_ESMValTool +
                                             " && python " + "main.py " +
                                             project_info['RUNTIME']['xml'],
                                             shell=True)

    case = "pre"
except Exception, TOOL_error:
    ESMValTool_log = str(TOOL_error.output)
    case = "post"

# run server
app = Flask(__name__)
host_add = "127.0.0.1"
host_port = 5000
full_host = host_add + ":" + str(host_port)

# initialize HTML production
D = HTML_writer()

if case == "pre":

    # Home
    left_box = '<p>This is the output from the ESMValTool run. \
                Please check if any errors occur!</p>'
    D.home_html(left_box,
                middle_box=ESMValTool_log,
                right_box=[project_info['RUNTIME']['xml_name']])
    # TODO Errors optional

    k0 = 1

    Key_list = list()
    var_list = list()

    for currDiag in project_info['DIAGNOSTICS']:

        requested_vars = currDiag.get_variables_list()
        requested_vars = list(
            diag_o.__dict__['var'] for diag_o in requested_vars
            )
        var_list.extend(requested_vars)

        # Update currDiag-specific models
        project_info['MODELS'] = projects.remove_diag_specific_models(
            project_info['MODELS'])
        diag_specific_models = currDiag.get_diag_models()
        projects.add_model(project_info, diag_specific_models)
        variables = currDiag.get_variables()
        field_types = currDiag.get_field_types()

        project_info['RUNTIME']['currDiag'] = currDiag
        for derived_var, derived_field in zip(variables, field_types):
            project_info['RUNTIME']['derived_var'] = derived_var
            project_info['RUNTIME']['derived_field_type'] = derived_field

        thisfile = main_ESMValTool + \
            ".".join(
                project_info['RUNTIME']['currDiag'].
                diag_script_cfg.split(".")[1:]
                )

#        with open(thisfile) as f:
#            cfg = f.readlines()

        cfg = thisfile

        Key = str(k0).zfill(3)
        k0 += 1

        DKEY = "Auto_Diag_"+Key
        # Diagnostic
        D.diag_html(DKEY, project_info['GLOBAL']['plot_dir'], cfg, [DKEY],
                    full_host,
                    time=(timestamp -
                          datetime.datetime.utcfromtimestamp(0)).
                    total_seconds()
                    )

        Key_list.append({'key': DKEY, 'name': DKEY})

    Key_list = {"diagnostic": Key_list}
    Key_list["home"] = None

    # Index
    D.index_html(list_of_vars=list(set(var_list)),
                 list_of_keys=Key_list,
                 version=version,
                 list_of_authors=authors)

elif case == "post":

    # Home
    home_text = "THIS IS A TAG-NAMELIST GENERATED REPORT!" + \
                "THERE IS NO DIRECT LOG FROM THE ESMVaLTool!"
    title_alt = 'ESMValTool'
    ESMValTool_Image = '<figure><a ' + \
        'title=' + title_alt + '>' + \
        '<img style="background-color:white;margin-left:auto;' + \
        'margin-right:auto; display:block;" ' + \
        'alt=' + title_alt + \
        ' src="' + '../static/images/ESMValTool-logo.png' + \
        '" width="2000"' + '></a></figure>'
    D.home_html(left_box=home_text,
                middle_box=ESMValTool_Image,
                right_box=[project_info['RUNTIME']['xml_name']])

    Key_list = list()

    for t in project_info['TAGS'].tags:
        DKEY = "tag_"+"_".join(t.split(","))
        Key_list.append({'key': DKEY, 'name': t.split(",")})
        # Tag-sites
        D.tag_html(DKEY, project_info['FOLDERS'].folders,
                   t.split(","), full_host)

    Key_list = {"diagnostic": Key_list}
    Key_list["home"] = None

    # Index
    D.index_html(list_of_vars=project_info['TAGS'].tags,
                 list_of_keys=Key_list,
                 version=version,
                 list_of_authors=authors)


@app.route('/')
def home():
    return render_template('home.html')


@app.route('/<tagname>')
def diagnostic(tagname):
    return render_template(tagname+'.html')

if __name__ == '__main__':
    app.debug = False
    app.run(host=host_add, port=host_port)
