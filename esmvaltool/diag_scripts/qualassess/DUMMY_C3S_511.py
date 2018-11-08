#!/usr/bin/env python

"""
;#############################################################################
; DUMMY Diagnostics
; Author: Benjamin Mueller (LMU Munich, GER)
; C3S_511 project
;#############################################################################
; Description
;    Produces various general diagnostic plots and statistics for the
;    reference data sets of the QA4ECV Project
;
; Required diag_script_info attributes (diagnostics specific)
;    none
;
; Optional diag_script_info attributes (diagnostic specific)
;    none
;
; Required variable_info attributes (variable specific)
;    none
;
; Optional variable_info attributes (variable specific)
;    none
;
; Caveats
;
; Modification history
;    20180619-A_muel_bn: port to version 2
;    20161128-A_laue_ax: added call to write_references
;    20160818-A_muel_bn: Routines written.
;
;#############################################################################
"""
# Basic Python packages
import logging
import os

from esmvaltool.diag_scripts.shared import run_diagnostic
from auxiliary.c3s_511_basic import Basic_Diagnostic_SP

logger = logging.getLogger(os.path.basename(__file__))

def main(cfg):
    logger.info('>>>>>>>> DUMMY_C3S_511.py is running! <<<<<<<<<<<<')

    for filename, attributes in cfg['input_data'].items():
            logger.info("Processing variable %s from model %s",
                        attributes['standard_name'], attributes['dataset'])
            logger.debug("Preparing diagnostic")
            Diag = Basic_Diagnostic_SP()
            Diag.set_info(cfg=cfg)
            logger.debug("Loading %s", filename)
            Diag.read_data()
            logger.debug("Running computation")
            Diag.run_diagnostic()

    logger.info('>>>>>>>> ENDED SUCESSFULLY!! <<<<<<<<<<<<')

if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
