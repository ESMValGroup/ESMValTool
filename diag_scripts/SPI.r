# #############################################################################
# SPI.r
# Author: Boris Orlowsky (ETH, Switzerland)
# Contributor: Martin Evaldsson (SMHI, Sweden)
# #############################################################################
# Description
#    For each month, the precipitation over the preceding TIMESCALE months, x,
#    is summed. Then a two-parameter Gamma distribution of cumulative
#    probability, Gamma_\alpha,\beta , is fitted to the strictly positive
#    TIMESCALE month sums, such that the probability of a non-zero
#    precipitation sum being below a certain value x corresponds to
#    Gamma_\alpha,\beta (x). We estimate shape parameter \alpha and scale
#    parameter \beta with a maximum likelihood approach. If the estimation
#    does not converge, \alpha and \beta are approximated using empirical
#    relations (Bordi et al., 2001). Accounting for TIMESCALE month periods
#    of no precipitation, occurring at a frequency q, the total cumulative
#    probability distribution of a precipitation sum below x, H(x), becomes
#    H(x) = q + (1 - q)*Gamma_\alpha,\beta (x). In the last step, a
#    precipitation sum x is assigned to its corresponding SPI value by
#    computing the quantile q_N(0,1) of the standard normal distribution at
#    probability H(x). The SPI of a precipitation sum x, thus, corresponds
#    to the quantile of the standard normal distribution which is assigned
#    by preserving the probability of the original precipitation sum, H(x).
#
# Required diag_script_info attributes (diagnostics specific)
#
# Optional diag_script_info attributes (diagnostic specific)
#
# Required variable_info attributes (variable specific)
#
# Optional variable_info attributes (variable specific)
#
# Caveats
#    No output of processed files to log-file yet!
#
# Modification history
#    20151113-A_laue_ax: added header
#
# ############################################################################

source('diag_scripts/aux/SPI/SPI_auxiliary.r')
source('interface_data/r.interface')
source(diag_script_cfg)
source('diag_scripts/lib/R/info_output.r')

## Do not print warnings
options(warn=-1)

var0 <- variables[1]
field_type0 <- field_types[1]

info_output(paste0("<<<<<<<< Entering ", diag_script), verbosity, 4)
info_output("+++++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
info_output(paste0("plot - ", diag_script, " (var: ", variables[1], ")"), verbosity, 1)
info_output("+++++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

library(tools)
diag_base = file_path_sans_ext(diag_script)

## Create working dirs if they do not exist
dir.create(plot_dir, showWarnings = FALSE)
dir.create(file.path(plot_dir, diag_base), showWarnings = FALSE)
dir.create(climo_dir, showWarnings = FALSE)

##
## Run it all
##
for (model_idx in c(1:length(models_name))) {
    fullpath_filename <- interface_get_fullpath(var0, field_type0, model_idx)
    filename <- interface_get_infile(var0, field_type0, model_idx)
    spi.info <- generate.spi(fullpath_filename=fullpath_filename,
                             filename=filename,
                             start_year=models_start_year[model_idx],
                             end_year=models_end_year[model_idx],
                             timescale=timescale,
                             begin.ref.year=begin.ref.year,
                             end.ref.year=end.ref.year, 
                             model_idx)
    
    evaluate.spi(spi.info)
    
    ##
    ## A simple map of the trend
    ##
    requested_data_type <- "_trends-n-pvals_wrt_"
    in.filename <- processed_spi_info(timescale, requested_data_type,
                                      begin.ref.year, end.ref.year,
                                      filename)

    suppressPackageStartupMessages(require("ncdf4"))
    eval.nc <- nc_open(file.path(climo_dir, in.filename))
    lon <- ncvar_get(eval.nc, "Lon")
    lat <- ncvar_get(eval.nc, "Lat")
    trends <- ncvar_get(eval.nc, "Decadal_Trend")
    p.values <- ncvar_get(eval.nc, "p_Value")
    
    ## Check if lats are ascending
    if (lat[1] > lat[length(lat)]) {
        lat <- rev(lat)
        trends <- trends[, dim(trends)[2]:1, ]
        p.values <- p.values[, dim(p.values)[2]:1, ]
    }
    
    ## Same range for all seasons
    max <- max(abs(trends), na.rm = TRUE)
    
    ## Choose the world
    {
        if (min(lon) < -175) 
            my.world <- "world"
        else my.world <- "world2"
    }
    
    suppressPackageStartupMessages(require("fields"))
    suppressPackageStartupMessages(require("maps"))
    
    ## Output figure filename
    n_diag_script <- nchar(diag_script)
    diag_script_base <- substr(diag_script, 0, n_diag_script - 2)
    description <- "_figure_"
    aux_info <- processed_spi_info(timescale, description,
                                   begin.ref.year, end.ref.year,
                                   "")
    figure_filename <- interface_get_figure_filename(diag_script_base,
                                                     var0,
                                                     field_type0,
                                                     aux_info,
                                                     model_idx)
    figure_filename <- paste(figure_filename, output_file_type, sep=".")
    figure_filename <- file.path(plot_dir, diag_base, figure_filename)
    
    ## Chose output format for figure
    if (tolower(output_file_type) == "png") {
        png(filename = figure_filename,
            width <- png_width,
            height <- png_height,
            units <- png_units,
            pointsize <- png_pointsize,
            bg <- png_bg) 
    } else if (tolower(output_file_type) == "pdf") {
        pdf(file <- figure_filename)
    } else if (tolower(output_file_type) == "eps") {
        setEPS()
        postscript(figure_filename)
    }
    par(mfrow = c(2, 3))
    for (i in 1:5) {
        ## Plot trends
        image(lon, lat, trends[, , i], main = paste(seasons[i], 
              "SPI trends [1/decade] (", models_name[model_idx], ")"),
               xlab = "lon", ylab = "lat", 
               zlim = c(-spi_colorbar_max, spi_colorbar_max),
               col = my.colors(11))
    
        ## Add significance level as contour
        contour(lon, lat, p.values[, , i], levels = 0.1, add = TRUE)
    
        ## Brave new world
        map(my.world, add = TRUE, int = FALSE)
    }
    
    ## Add colorbar
    plot.new()
    image.plot(legend.only = TRUE, zlim = c(-spi_colorbar_max, spi_colorbar_max), 
               col = my.colors(11), 
        smallplot <- c(0.1, 0.2, 0.3, 0.7))
    
    invisible(dev.off())
}
info_output(paste0(">>>>>>>> Leaving ", diag_script), verbosity, 4)
