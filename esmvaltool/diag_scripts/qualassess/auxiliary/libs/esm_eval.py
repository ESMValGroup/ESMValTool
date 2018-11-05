# -*- coding: utf-8 -*-
"""
Created on Fri Feb 09 09:29:45 2018

@author: hass_bg
"""

import csv
from .reporting import do_report


def create_esmeval_dict(esm_eval_csv, ecv_name, work_dir):
    """
    - esm_eval_csv    is a csv file containing the information about the ESM
                      evaluation assessment
    - ecv_name        is a string that defines which ECV is looked at in the
                      single assessment report
    """

    with open(esm_eval_csv, 'rb') as csvfile:
        s = csv.reader(csvfile, delimiter=",")
        my_list = list(s)

    # check which ECV is supposed to be displayed and where the
    result = [l for l in my_list for j in range(len(l)) if l[j] == ecv_name]

    # create a dictionary that will be handed over to the 'do_report' routine
    esm_eval_dict = {}
    if len(result) >= 1:
        #esm_eval_dict[my_list[0][0]] = result[0][0]
        for num_entries in range(len(result)):
            for num_keys in range(len(my_list[0])):
                esm_eval_dict[my_list[0][num_keys]
                              ] = result[num_entries][num_keys]

    do_report(esm_eval_dict, "ESM Evaluation", work_dir)

#create_esmeval_dict("D:\project\ESM_Eval_list.csv", 'ozone')
