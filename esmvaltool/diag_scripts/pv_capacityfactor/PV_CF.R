# #############################################################################
# diagnostic.R
# Authors:       Irene Cionni (ENEA, Italy)
# #############################################################################
# Description
# This script calculates the capacity factor used in pv_capacity_factor.R
#
# Required
#
# Optional
#
# Caveats
#
# Modification history
#    20210401-cionni_irene: written for v2.0
#    20210621-weigel_katja: formattiong updates
#    20210714-weigel_katja: removed unused parts
#
# ############################################################################library(ggplot2)
library(plyr)
rsds2cf <- function(rsds, t2m) {
  cf <- 0.9 * (1.0 - 0.00042 * (t2m - (273.15 + 25.))) * rsds / 1000
  return(cf)
}