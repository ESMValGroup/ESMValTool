#!/usr/bin/env python

"""
;#############################################################################
; albedo Diagnostics
; Author: Benjamin Mueller (LMU Munich, GER)
; C3S_511 project
;#############################################################################
; Description
;    Produces various general diagnostic plots and statistics for albedo
;    data sets of the CDS
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
;    20190401-A_muel_bn: Dummy ported to work with albedo
;
;#############################################################################
"""
# Basic Python packages
import logging
import os

from esmvaltool.diag_scripts.shared import run_diagnostic
from auxiliary.c3s_511_albedo import albedo_Diagnostic_SP

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    logger.info('>>>>>>>> albedo_C3S_511_SPQB.py is running! <<<<<<<<<<<<')
    


    for filename, attributes in cfg['input_data'].items():
        if not attributes['short_name'] in ["bdalb","bhalb"]:
            logger.error("This is the wrong variable: " + 
                         attributes['short_name'])
            assert False, "not assessing the correct ECVs bdalb or bhalb"
        logger.info("Processing variable %s from data set %s",
                    attributes['short_name'], attributes['dataset'])
        logger.info("Preparing diagnostic")
        Diag = albedo_Diagnostic_SP()
        Diag.set_info(cfg=cfg)
        logger.info("Loading %s", filename)
        Diag.read_data()
        logger.info("Running computation")
        Diag.run_diagnostic()

    logger.info('>>>>>>>> ENDED SUCCESSFULLY!! <<<<<<<<<<<<')


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
