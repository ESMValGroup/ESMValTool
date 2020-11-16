library(ggplot2)
library(plyr)

read_pc <- function(file) {
  pc <- list()
  pc$points <- rbind(c(0, 0), read.delim(file, comment.char = "#"))
  pc$fun <- approxfun(
    pc$points$WindSpeed,
    pc$points$Power,
    method = "linear",
    yleft = NA,
    yright = 0
  )
  attr <- strsplit(
    trimws(system(
      paste(
        "perl -e 'open FH,\"",
        file,
        "\";while(<FH>){@parts= /^# (.+): (.+) /;print \"@parts \";}'", # nolint
        sep = ""
      ),
      intern = TRUE
    )),
    "\\s+"
  )
  attr <- matrix(unlist(attr), ncol = 2, byrow = T)
  pc$attr <- as.list(attr[, 2])
  names(pc$attr) <- attr[, 1]
  pc$attr$Filename <- file # nolint
  pc$attr$RatedPower <- as.numeric(pc$attr$RatedPower) # nolint
  return(pc)
}
read_xml_pc <- function(file) {
  xml <- xmlTreeParse(file, useInternalNodes = TRUE) # nolint
  xml_data <- xmlToList(xml) # nolint
  pc <- list()
  pcs <- xml_data$wind_turbine_properties$power_curves
  for (i in seq_along(pcs)) {
    if (pcs[[i]]$air_density == 1.225) {
      pc$points <- ldply(
        pcs[[i]]$power_curve_table, # nolint
        data.frame
      )[, c(2, 3)]
      colnames(pc$points) <- c("WindSpeed", "Power") # nolint
      pc$points <- transform(
        pc$points,
        WindSpeed = as.numeric(as.character(WindSpeed)),
        # nolint
        Power = as.numeric(as.character(Power))
      )
      pc$points <- rbind(c(0, 0), pc$points)
      break
    }
  }
  pc$fun <- approxfun(
    pc$points$WindSpeed,
    # nolint
    pc$points$Power,
    # nolint
    method = "linear",
    yleft = NA,
    yright = 0,
    ties = "ordered"
  )
  pc$attr$Diameter <-
    xml_data$wind_turbine_properties$rotor_diameter # nolint
  pc$attr$CutIn <- NA # nolint
  pc$attr$CutOut <- NA # nolint
  pc$attr$ReCutIn <- NA # nolint
  pc$attr$RatedSpeed <- NA # nolint
  pc$attr$RatedPower <-
    xml_data$wind_turbine_properties$rated_power # nolint
  pc$attr$IECClass <- NA # nolint
  pc$attr$Control <- NA # nolint
  pc$attr$Density <- 1.225 # nolint
  pc$attr$Name <- file # nolint
  pc$attr$Filename <- file # nolint
  pc$attr$RatedPower <- as.numeric(pc$attr$RatedPower) # nolint
  return(pc)
}
plot_pc <- function(pc) {
  plot <- ggplot(pc$points, aes(x = WindSpeed, y = Power)) + # nolint
    geom_point() +
    stat_function(fun = pc$fun) +
    xlim(0, 35)
  return(plot)
}
plot_pc_list <- function(list_pcs) {
  list_funs <- lapply(list_pcs, function(x) {
    function(y) {
      x$fun(y) / x$attr$RatedPower
    } # nolint
  })
  names <- lapply(list_pcs, function(x)
    x$attr$Name) # nolint
  plot <- ggplot(NULL, aes(x = x, colour = Turbine)) # nolint
  for (i in seq_along(list_pcs)) {
    plot <- plot + stat_function(
      data = data.frame(
        x = 0:30,
        Turbine = factor(names[[i]])
      ),
      fun = list_funs[[i]]
    )
  }
  plot <-
    plot + xlab("Wind speed (m/s)") + ylab("Capacity Factor (%)") +
    ggtitle("Selected power curves")
  return(plot)
}
get_list_turbines <- function() {
  files <- list.files()
  turb_list <- list()
  for (i in seq(files)) {
    file <- files[i]
    turb_list[[i]] <- read_xml_pc(file)
  }
  names(turb_list) <- files
  return(turb_list)
}
wind2power <- function(wind, pc) {
  power <- pc$fun(wind)
}
wind2CF <- function(wind, pc) {
  power <- pc$fun(wind)
  CF <- power / pc$attr$RatedPower # nolint
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
