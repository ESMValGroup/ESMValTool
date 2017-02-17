begin.ref.year <- 1970
end.ref.year <- 2000
timescale <- 3  ## Valid values are 3, 6 and 12 months
seasons <- c("ann", "djf", "mam", "jja", "son")
spi_colorbar_max <- 0.75
my.colors <- colorRampPalette(c("brown","orange","white","lightblue","blue"))

## Specific settings for PNG output
png_width <- 1600
png_height <- 960
png_units <- "px"
png_pointsize <- 12
png_bg <- "white"
