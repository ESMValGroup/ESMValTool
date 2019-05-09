# These are default values, loaded before the recipe is read
regrid_dataset <- NA
base_range <- NA
analysis_range <- NA
climdex_parallel <- 4
mip_name <- "cmip"
ts_col_list <- c("dodgerblue2", "darkgreen",
                 "firebrick2", "darkorchid", "aquamarine3")
ts_png_width <- 640
ts_png_height <- 480
ts_png_units <- "px"
ts_png_pointsize <- 12
ts_png_bg <- "white"
ts_lty_list <- c(1, 4, 2, 3, 5)
ts_lwd_list <- c(2, 2, 2, 2, 2)
ts_data <- TRUE
normalize <- FALSE
timeseries_idx <- c(
  "tn10pETCCDI_yr", "tn90pETCCDI_yr",
  "tx10pETCCDI_yr", "tx90pETCCDI_yr"
)
gleckler_idx <- c(
  "tn10pETCCDI_yr", "tn90pETCCDI_yr",
  "tx10pETCCDI_yr", "tx90pETCCDI_yr"
)

ts_plt <- TRUE
glc_plt <- TRUE
glc_arr <- FALSE
gl_mar_par <- c(10, 4, 3, 14)
gl_png_res <- 480
gl_png_units <- "px"
gl_png_pointsize <- 12
gl_png_bg <- "white"
gl_rmsespacer <- 0.01
gl_scaling_factor <- 0.9
gl_text_scaling_factor <- 1.0
gl_xscale_spacer_rmse <- 0.05
gl_xscale_spacer_rmsestd <- 0.25
gl_symb_scaling_factor <- 1.0
gl_symb_xshift <- 0.2
gl_symb_yshift <- 0.275
gl_text_symb_scaling_factor <- 0.5
