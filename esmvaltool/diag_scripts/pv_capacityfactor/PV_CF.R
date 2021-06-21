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
#
# ############################################################################library(ggplot2)
library(plyr)
rsds2CF <- function(rsds, t2m) {
# power <- pc$fun(wind)
  CF <- 0.9 * (1.0 - 0.00042 * (t2m - (273.15 + 25.))) * rsds / 1000 # nolint
}
WPD <- function(wind, ro) {
  return(0.5 * ro * wind^3)
}
bump <- function(x) {
  f <- function(y) {
    exp(-1 / y^2)
  }
  return(f(x) / (f(x) + f(1 - x)))
}
