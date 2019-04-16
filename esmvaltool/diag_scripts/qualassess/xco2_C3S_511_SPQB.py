#!/usr/bin/env python

"""
;#############################################################################
; xco2 Diagnostics
; Author: Birgit Hassler (DLR Oberpfaffenhofen, GER)
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
;    20190416-A_hass_bg: adjusted for xco2
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
from auxiliary.c3s_511_xco2 import xco2_Diagnostic_SP

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    logger.info('>>>>>>>> xco2_C3S_511_SPQB.py is running! <<<<<<<<<<<<')
#
#    if len(cfg['input_data'].items()) == 0:
#        logger.info("Preparing diagnostic")
#        Diag = ta_Diagnostic_SP()
#        logger.info("Running computation")
#        Diag.run_diagnostic(cfg=cfg)

    for filename, attributes in cfg['input_data'].items():
        logger.info("Processing variable %s from data set %s",
                    attributes['short_name'], attributes['dataset'])
        logger.info("Preparing diagnostic")
        Diag = xco2_Diagnostic_SP()
        Diag.set_info(cfg=cfg)
        logger.info("Loading %s", filename)
        Diag.read_data()
        logger.info("Running computation")
        Diag.run_diagnostic()

    logger.info('>>>>>>>> ENDED SUCCESSFULLY!! <<<<<<<<<<<<')


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)