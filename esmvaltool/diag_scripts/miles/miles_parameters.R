## Analysis options
# seasons <- c("DJF") # set in namelist

# select which EOFs you want to compute
# "NAO": the 4 first  EOFs of North Atlantic, i.e. North Atlantic Oscillation as EOF1
# "AO" : the 4 first EOFs of Northern Hemispiere, i.e. Arctic Oscillation as EOF1 
# teles <- c("NAO") # set in namelist

# select how many clusters for k-means over the North Atlantic
# NB: only 4 clusters supported so far.  
nclusters=4

# Specific settings for PNG output
png_width=960
png_height=960
png_units="px"
png_pointsize=12
png_bg="white"

# Specific settings for PDF and EPS output (in inches)
pdf_width=12
pdf_height=12

#color palette to be used
#palette0 is taken from tim.colors of field to avoid library dependencies...
palette0=colorRampPalette(c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF",
        "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF",
        "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF",
        "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF",
        "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF",
        "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F",
        "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50",
        "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00",
        "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00",
        "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000",
        "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000",
        "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000",
        "#AF0000", "#9F0000", "#8F0000", "#800000"))
palette1=colorRampPalette(c("white","orange","darkred"))
palette2=colorRampPalette(c("blue","white","red"))
palette3=colorRampPalette(c("darkblue","blue","dodgerblue","white","orange","red","darkred"))

