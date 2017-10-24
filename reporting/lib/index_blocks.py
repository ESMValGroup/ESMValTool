#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 10:41:38 2017

@author: bmueller
"""

import re


def index_head(inputs):
    """ this block gives the head of the html index """
#    input_desc=[inputs[0] if len(inputs)==1 else ", ".join(inputs)][0]
    # 'Diagnostics results for ' + input_desc
    Headline = 'Results from the ESMValTool'
    lines = [
        '<!DOCTYPE html>',
        '<html>',
        '<head>',
        '<title>' + Headline + '</title>',
        '<meta http-equiv="Content-Type" content="application/xhtml+xml;' +
        ' charset=utf-8" />',
        '<meta name="description" content="' + Headline + '" />',
        '<meta name="keywords" content="' + Headline + '" />',
        '<meta name="robots" content="index, follow" />',
        '<link rel="stylesheet" type="text/css"' +
        ' href="{{url_for("static",filename="css/styles.css")}}"' +
        ' media="screen" />',
        '<link rel="stylesheet" type="text/css"' +
        ' href="{{url_for("static",filename="css/styles_print.css")}}"' +
        ' media="print" />',
        '</head>',
        '<body>',
        '<div id="header">',
        '<h1>' + Headline + '</h1>',
        '<ul>'
    ]
    return lines


def index_tabs(inputs, actual_tab=None):
    """ this block gives the head of the html index and the tocs """
    lines = []

    active_class = 'class="active"'

    def tabline(hook, info=None):
        """ returns ref for tag """
        if info is None:
            res = ['<li><a href="{{url_for("' + hook + '")}}"' + ' ' +
                   (active_class if actual_tab == hook else '') + '>' +
                   hook.upper() + '</a></li>']
        elif isinstance(info, list):
            res = []
            for i in info:
                if isinstance(i['name'], (str, unicode)):
                    res.extend(['<li><a href="{{url_for("' +
                                hook + '",tagname="' + i['key'] + '")}}"' +
                                ' ' + (active_class if actual_tab == i['key']
                                       else '') + '>' +
                                i['name'].upper() + '</a></li>'])
                elif isinstance(i['name'], list):
                    res.extend(['<li><a href="{{url_for("' +
                                hook + '",tagname="' + i['key'] + '")}}"' +
                                ' ' + (active_class if actual_tab == i['key']
                                       else '') + '>' +
                                "/".join(i['name']).upper() + '</a></li>'])
                else:
                    assert False, "this info name is not expected!"
        else:
            assert False, "this info is not expected!"

        return res

    for i in inputs:
        lines.extend(tabline(i, inputs[i]))

    return lines


def index_foot(version, inputs):
    """ this block gives the head of the html index """
    inputs = ['<a href="mailto:' + string + '">' + string +
              '</a>' if re.search('@', string)
              else string for string in inputs]
    input_desc = [inputs[0] if len(inputs) == 1 else ", ".join(inputs)][0]
    lines = [
        '<form>',
        '<li><a href="javascript:window.print()"> print </a></li>',
        '</form>',
        '</ul>',
        '</div>',
        '<div class="colmask threecol">',
        '<div class="colmid">',
        '<div class="colleft">',
        '<div class="col1">',
        '<!-- Column 1 (middle) start -->',
        '<div class="container">',
        '{% block content1 %}',
        '{% endblock content1 %}',
        '</div>',
        '<!-- Column 1 end -->',
        '</div>',
        '<div class="col2">',
        '<!-- Column 2 (left) start -->',
        '{% block content2 %}',
        '{% endblock content2 %}',
        '<!-- Column 2 end -->',
        '</div>',
        '<div class="col3">',
        '<!-- Column 3 (right) start -->',
        '{% block content3 %}',
        '{% endblock content3 %}',
        '<!-- Column 3 end -->',
        '</div>',
        '</div>',
        '</div>',
        '</div>',
        '<div id="footer">',
        '<p>ESMValTool Reporting service v.' + version +
        ' | Contact(s): ' + input_desc + '</p>',
        '</div>',
        '</body>',
        '</html>'
    ]
    return lines
