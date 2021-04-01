library(ggplot2)
library(plyr)
rsds2CF <- function(rsds, t2m) {
# power <- pc$fun(wind)
  CF <- 0.9*(1-0.00042*(t2m-(273.15+25.)))*rsds/1000 # nolint
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
