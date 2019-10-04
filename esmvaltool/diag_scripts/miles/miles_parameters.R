##########################################################
#-------------Plot configurations------------------------#
##########################################################

# plot all the warnings
options(warn = 1)

# Specific settings for PNG output
png_width <- 900
png_height <- 900
png_units <- "px"
png_pointsize <- 12
png_bg <- "white"

# Specific settings for PDF and EPS output (in inches)
pdf_width <- 12
pdf_height <- 12

# aspect ratio
af <- 1

# Type of projection ("no" for standard plotting")
# All projection from mapproj package should be supported
# but error may arise for non-polar plots
# DEFAULT IS POLAR PLOT
map_projection <- "azequalarea" # Alternative: "azequidistant"

# Number of panels per figure (rows and column): default for polar plots
panels <- c(3, 1)

# Latitudinal range for plots
lat_lim <- c(25, 90)

# if not regular projection (i.e. if using polar)
if (map_projection != "no") {
  af <- round(sqrt(3), 2)
  pdf_height <- pdf_height / af
  pdf_width <- 3 * pdf_width / af
  png_height <- png_height / af
  png_width <- 3 * png_width / af
  panels <- rev(panels)
}

# Custom parameteres for plots
plotpar <-
  list(
    mfrow = panels,
    cex.main = 2.5,
    cex.axis = 1.5,
    cex.lab = 1.5,
    mar = c(5, 5, 5, 7),
    oma = c(1, 1, 3, 2)
  )

# imagescale3 color bar details
imgscl_colorbar <- 1.4
imgscl_label <- 1.5
imgscl_line <- 3

# color palette to be used
# palette0 is taken from tim.colors of field to avoid library dependencies...
palette0 <- colorRampPalette(
  c(
    "#00008F",
    "#00009F",
    "#0000AF",
    "#0000BF",
    "#0000CF",
    "#0000DF",
    "#0000EF",
    "#0000FF",
    "#0010FF",
    "#0020FF",
    "#0030FF",
    "#0040FF",
    "#0050FF",
    "#0060FF",
    "#0070FF",
    "#0080FF",
    "#008FFF",
    "#009FFF",
    "#00AFFF",
    "#00BFFF",
    "#00CFFF",
    "#00DFFF",
    "#00EFFF",
    "#00FFFF",
    "#10FFEF",
    "#20FFDF",
    "#30FFCF",
    "#40FFBF",
    "#50FFAF",
    "#60FF9F",
    "#70FF8F",
    "#80FF80",
    "#8FFF70",
    "#9FFF60",
    "#AFFF50",
    "#BFFF40",
    "#CFFF30",
    "#DFFF20",
    "#EFFF10",
    "#FFFF00",
    "#FFEF00",
    "#FFDF00",
    "#FFCF00",
    "#FFBF00",
    "#FFAF00",
    "#FF9F00",
    "#FF8F00",
    "#FF8000",
    "#FF7000",
    "#FF6000",
    "#FF5000",
    "#FF4000",
    "#FF3000",
    "#FF2000",
    "#FF1000",
    "#FF0000",
    "#EF0000",
    "#DF0000",
    "#CF0000",
    "#BF0000",
    "#AF0000",
    "#9F0000",
    "#8F0000",
    "#800000"
  )
)
palette1 <- colorRampPalette(c("white", "orange", "darkred"))
palette2 <- colorRampPalette(c("blue", "white", "red"))
palette3 <- colorRampPalette(c(
  "darkblue",
  "blue",
  "dodgerblue",
  "white",
  "orange",
  "red",
  "darkred"
))

# additional color palette used for extradiagnostics histogram
KOL <- c(
  "black",
  "darkgreen",
  "blue",
  "darkorange",
  "red",
  "violet",
  "grey50",
  "black"
)
