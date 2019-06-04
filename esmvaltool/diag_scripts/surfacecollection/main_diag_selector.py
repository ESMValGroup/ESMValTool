#!/usr/bin/env python

"""
;#############################################################################
; diagnostic selector
; Author: Benjamin Mueller (LMU Munich, GER)
; Multiple projects
;#############################################################################
; Description
;    Produces various diagnostics based on model data and a reference data set for:
;    - soil moisture
;
; Required diag_script_info attributes (diagnostics specific)
;    request: list of names of requested sub-diagnostics
;              (minimal: None; returns all available request names)
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
;    20190603-A_muel_bn: port to version 2
;
;#############################################################################
"""
# Basic Python packages
import logging
import os

from esmvaltool.diag_scripts.shared import run_diagnostic
#from auxiliary.c3s_511_basic import Basic_Diagnostic_SP

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    logger.info('>>>>>>>> diagnostic selector is running! <<<<<<<<<<<<')

    for filename, attributes in cfg['input_data'].items():
        logger.info("Processing variable %s from data set %s",
                    attributes['short_name'], attributes['dataset'])
        #logger.info("Preparing diagnostic")
        #Diag = Basic_Diagnostic_SP()
        #Diag.set_info(cfg=cfg)
        #logger.info("Loading %s", filename)
        #Diag.read_data()
        #logger.info("Running computation")
        #Diag.run_diagnostic()

    logger.info('>>>>>>>> ENDED SUCCESSFULLY!! <<<<<<<<<<<<')


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
