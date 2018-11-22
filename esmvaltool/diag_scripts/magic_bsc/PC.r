######################################3
# Power Curve routines
# llledo@bsc.es - Apr 2016
######################################3
library(ggplot2)
library(XML) # nolint
library(plyr)

#=======================
# Read a powercurve file
# Create the approximation function
# TODO: Document file format
#=======================
read_pc <- function(file)
{	#----------------
  # Create empty list
  #----------------
  pc <- list()

  #----------------
  # Add pc points
  #----------------
  pc$points <- rbind(c(0,0), read.delim(file,comment.char="#"))

  #----------------
  # Create an approximating function
  #----------------
  pc$fun <- approxfun(pc$points$WindSpeed,pc$points$Power,method="linear",yleft=NA,yright=0)

  #----------------
  # Read other turbine characteristics (using perl's grep)
  #----------------
  attr <- strsplit(trimws(system(paste(
    "perl -e 'open FH,\"", file, "\";while(<FH>){@parts= /^# (.+): (.+) /;print \"@parts \";}'",
    sep=""),
    intern=T)),
    "\\s+")
  attr <- matrix(unlist(attr),ncol=2,byrow=T)
  pc$attr <- as.list(attr[,2])
  names(pc$attr) <- attr[,1]
  pc$attr$Filename <- file

  # TODO: check attributes!
  pc$attr$RatedPower <- as.numeric(pc$attr$RatedPower)

  return(pc)
}

#=======================
# Read a powercurve file in windographer xml format
# Create the approximation function
#=======================
read_xml_pc <- function(file)
{	#----------------
  # Read xml and convert to nested list
  #----------------
  xml <- xmlTreeParse(file,useInternalNodes=T)
  xml_data <- xmlToList(xml)

  #----------------
  # Create empty list
  #----------------
  pc <- list()

  #----------------
  # Loop all pcs at different densities
  # Add pc points for standard density
  #----------------
  pcs <- xml_data$wind_turbine_properties$power_curves
  for (i in 1:length(pcs))
  {	# get pc for standard air density
    if(pcs[[i]]$air_density == 1.225)
    {	pc$points <- ldply(pcs[[i]]$power_curve_table,data.frame)[,c(2,3)]
    colnames(pc$points) <- c("WindSpeed","Power")
    pc$points <- transform(
      pc$points,
      WindSpeed=as.numeric(as.character(WindSpeed)),
      Power=as.numeric(as.character(Power))
    )
    pc$points <- rbind(c(0,0),pc$points)
    break
    }
  }
  #----------------
  # Create an approximating function
  #----------------
  pc$fun <- approxfun(
    pc$points$WindSpeed,
    pc$points$Power,
    method="linear",
    yleft=NA,
    yright=0,
    ties="ordered"
  )
  #----------------
  # Read other turbine characteristics
  #----------------
  pc$attr$Diameter <- xml_data$wind_turbine_properties$rotor_diameter
  pc$attr$CutIn <- NA
  pc$attr$CutOut <- NA
  pc$attr$ReCutIn <- NA
  pc$attr$RatedSpeed <- NA
  pc$attr$RatedPower <- xml_data$wind_turbine_properties$rated_power
  pc$attr$IECClass <- NA
  pc$attr$Control <- NA
  pc$attr$Density <- 1.225
  pc$attr$Name <- file
  pc$attr$Filename <- file
  # TODO: check attributes!
  pc$attr$RatedPower <- as.numeric(pc$attr$RatedPower)
  return(pc)
}
#=======================
# Plot a pc, using points
# and approximating function
#=======================
plot_pc <- function(pc)
{	plot <- ggplot(pc$points,aes(x=WindSpeed,y=Power)) +
  #geom_line() +
  geom_point() +
  stat_function(fun = pc$fun) +
  xlim(0,35)
return(plot)
}

#=======================
# Plot a list of pcs, using
# approximating function
#=======================
plot_pc_list <- function(list_pcs){
  list_funs <- lapply(list_pcs, function(x) { function(y) {x$fun(y)/x$attr$RatedPower} })
  names <- lapply(list_pcs, function(x) x$attr$Name )
  plot <- ggplot(NULL, aes(x=x, colour = Turbine))
  for(i in 1:length(list_pcs)){
    plot <- plot + stat_function(data=data.frame(x=0:30,Turbine=factor(names[[i]])),fun = list_funs[[i]])
  }
  plot <- plot + xlab("Wind speed (m/s)") + ylab("Capacity Factor (%)") + ggtitle("Selected power curves")
  return(plot)
}

#=======================
# Read windographer pcs
#=======================
get_list_turbines <- function(){
  files <- list.files()
  turb_list <- list()
  for (i in seq(files)){
    file <- files[i]
    turb_list[[i]]<- read_xml_pc(file)
  }
  names(turb_list) <- files
  return(turb_list)
}

#=======================
# plotly or png a selection of turbines
#=======================
# library(plotly)
# ggp <- plot_pc_list(selected)
# saveWidget(ggplotly(ggp), "selected_5_turbines.html")
# ggsave("5turbines.png",ggp,width=12,height=6)


#=======================
# Cross wind with power curve, and obtain power.
#=======================
wind2power <- function(wind,pc){
  power <- pc$fun(wind)
}

#=======================
# Cross wind with power curve, and obtain CF.
#=======================
wind2CF <- function(wind,pc){
  power <- pc$fun(wind)
  CF <- power / pc$attr$RatedPower
}

#=======================
# Wind power density
#=======================
WPD <- function(wind,ro){
  return(0.5 * ro * wind ^ 3)
}
# proves
# hist(density=T,breaks=30,wind2power(rweibull(1000000,shape=1.8,scale=10),pc))
#=======================
# Theoretical bump function
#=======================
bump <- function(x){
  f <- function(y) {
    exp(-1/y^2)
  }
  return(f(x) / (f(x) + f(1 - x)))
}
