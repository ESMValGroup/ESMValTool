# #############################################################################
# make_timeseries_plot.r
#
# Author: Marit Sandstad (CICERO, Norway)
#       : Christian W. Mohr (CICERO, Norway)
#
# #############################################################################
# Description
#    Code to plot a timeseries plot for a set of climdex indices
#
# Modification history
#    20180816-A_cwmohr: adding input procedure and plotting for/of observation data
#    20180725-A_cwmohr: modification of time cropping
#    20180618-A_cwmohr: alpha levels for polygon plotting, second y-axis, 
#    20180131-A_laue_ax: clean-up of code, adaptation to ESMValTool standards,
#                        added tagging, bugfixes: time axis, cdo, filenames
#    20170920-A_maritsandstad: Creation
#
# #############################################################################

source('diag_scripts/aux/ExtremeEvents/common_climdex_preprocessing_for_plots.r')


##
##
## Method to call all preprocessing and loop through
## all models and indices then call plotting script
## to produce time series plots for a list of indices
## @param path is the path to where the original indices 
## are stored
## @param idx_list lists the indices to be considered in
## this run. Defaults are the indices from the IPCC
## report.
##

## TESTING FUNCTION AREA ##
#path = out_dir
#idx_list =  timeseries_idx
#model_list = models_name 
#plot_dir = paste(plot_dir, "/", diag_base, sep = "")
###########################


timeseries_main <- function(path = "../work/ExtremeEvents", idx_list = c("tn10pETCCDI_yr", "tn90pETCCDI_yr", "tx10pETCCDI_yr", "tx90pETCCDI_yr"), 
                            model_list = c("IPSL-CM5A-LR", "MIROC-ESM"), obs_list = c("MERRA2"), plot_dir = "./plot", normalize=FALSE){

    ## For file structure and files
    tsGrid <-  paste(path, "/tsGridDef", sep = "")
    time_cropped <- paste(path, "/timeCropped", sep = "")
    landmask <- paste(path, "/landSeaMask.nc", sep = "")
    regridded <- paste(path, "/regridded", sep = "")
    land <- paste(path, "/Land", sep = "")

    # Initial nc-file time crop, regrid, land and plot purge
    cmd <- paste("rm -f ", time_cropped, "/*.nc", " ", regridded, "/*.nc", " ", land, "/*.nc", sep = "")
    print(cmd)
    system(cmd)
    
    # Initial grid and landmask creation reset
    gridAndLandmask = TRUE
    
    ## Loop over the indices to produce a plot for each index
    idx = idx_list[1]
    for (idx in idx_list){
    
        ## Combine the list of models and observations
        modelObs_list <- unique(c(model_list, obs_list))
        
        ## Find the model files 
        modelAndObs <-  basename(Sys.glob(file.path(path, paste(idx, "*.nc", sep = ""))))
        
        
      if(chk.ts_data){  
        ## Time crop
        #returnvalue <- setTimeForFilesEqual(path = path, idx = idx, model_list = model_list, time_cropped = time_cropped)
        returnvalue <- setTimeForFilesEqual(path = path, idx = idx, model_list = modelObs_list, time_cropped = time_cropped) # This is a temporary solution
        
        max_start <- returnvalue[1]
        min_end <- returnvalue[2]

        ## If there is no overlap in the files the index
        ## should be skipped
        if (max_start >= min_end){
            print(paste("No time overlap in files for index", idx))
            break()
        }
        
        ## Find the new model files after time cropping
        modelAndObs <-  basename(Sys.glob(file.path(time_cropped, paste(idx, "*.nc", sep = ""))))
        
        ## !New Grid and landseamask for each idx 
        ## !(or just the first idx set) should be
        ## !produced here
        if (gridAndLandmask){
          createGrid(path = path, loc = tsGrid)
          createLandSeaMask(regrid = tsGrid, loc = path, landmask = landmask)
          gridAndLandmask = FALSE
        }
        
        ## Loop over each file so it can be regridded
        ## and landseaMasked  
        for (m in modelAndObs){
          print(paste(time_cropped, "/", m, sep = ""))
          regridAndLandSeaMask(idx_raw =paste(time_cropped, "/", m, sep = ""), regrid = tsGrid,
                               landmask = landmask, regridded = regridded, land = land, loc = path)
        }
        
        
#               
#         ## THIS SECTION WILL BE REACTIVATED LATER ONCE OBSERVATION DATA IS AVAILABLE
#
#         ## Find the new model files after time cropping
#         models <-  basename(Sys.glob(file.path(time_cropped, paste(idx, "*.nc", sep = ""))))
#
#         ## Find the start year (to be used in plotting)
#         start_yr <- strtoi(substr(models[1], nchar(models[1]) - 11, nchar(models[1]) - 8))
#         ## !New Grid and landseamask for each idx 
#         ## !(or just the first idx set) should be
#         ## !produced here
#         if (gridAndLandmask){
#            createGrid(path = path, loc = tsGrid)
#            createLandSeaMask(regrid = tsGrid, loc = path, landmask = landmask)
#            gridAndLandmask = FALSE
#         }
# 
#         ## Loop over each file so it can be regridded
#         ## and landseaMasked  
#         for (m in models){
#             print(paste(time_cropped, "/", m, sep = ""))
#             regridAndLandSeaMask(idx_raw =paste(time_cropped, "/", m, sep = ""), regrid = tsGrid,
#                                  landmask = landmask, regridded = regridded, land = land, loc = path)
#         }
#         
#         
#         ## Find the observation files
#         ## This section will be later activated once there are real observation.
#         modelAndObs <-  basename(Sys.glob(file.path(path, paste(idx, "*.nc", sep = ""))))
#         o=obs_list
#         for (o in obs_list){
#           obs <- modelAndObs[grep(o, modelAndObs)]
#           print(paste(path, "/", obs, sep = ""))
#           regridAndLandSeaMask(idx_raw = paste(path, "/", obs, sep = ""), regrid = tsGrid,
#                                landmask = landmask, regridded = regridded, land = land, loc = path)
#         }

        ## Then do the preprocessing
        #time_series_preprocessing(path = time_cropped, land = land, idx = idx, plot_dir = plot_dir)
        
        time_series_preprocessing(land = land, idx = idx, model_list = model_list, obs_list = obs_list, plot_dir = plot_dir, normalize=normalize)
      }

        ## Find the start year (to be used in plotting)
        start_yr <- strtoi(substr(modelAndObs[1], nchar(modelAndObs[1]) - 11, nchar(modelAndObs[1]) - 8))

        ## Produce plot for this index
        timeseries_plot(plot_dir = plot_dir, idx = idx, obs_list = obs_list, start_yr = start_yr, normalize=normalize)
    }
}



##
## Method that preprocesses idx-files for a single index
## in order to get data to plot time series plot
## for this index
## @param path is the path to index file location
## @param idx is the index to be processed.
##
time_series_preprocessing <- function(land = "./Land", idx = 'tnnETCCDI_yr', 
                                      model_list = model_list, obs_list = obs_list, normalize = FALSE, plot_dir = './plot'){

  ## List of indices which are never normalized:
  pidx <- c('tn10pETCCDI_yr', 'tx10pETCCDI_yr', 'tn90pETCCDI_yr', 'tx90pETCCDI_yr', "csdiETCCDI_yr", "wsdiETCCDI_yr",
            'tn10pETCCDI_mon', 'tx10pETCCDI_mon', 'tn90pETCCDI_mon', 'tx90pETCCDI_mon', "csdiETCCDI_mon", "wsdiETCCDI_mon")
  
  
  # Getting a list of all the files for the index
  modelsAndObs <- basename(Sys.glob(file.path(land, paste(idx, "*.nc", sep = ""))))  
  modelsAndObsSplitList <- strsplit(modelsAndObs, split = "_")
  modelsAndObsSplit <- unlist(lapply(modelsAndObsSplitList, function(x){x[3]}))
  
  # Extracting only the model files
  models <- modelsAndObs[which(modelsAndObsSplit %in% model_list)]
  print("These are the models:")
  print(models)
  
  
  ## Extracting only the observation files
  obs_order <- which(modelsAndObsSplit %in% obs_list)
  obs <- modelsAndObs[obs_order]
  print("These are the observations:")
  print(obs)
    
  #### NORMALIZE VALUES ####
  if(normalize){
    #File string to be filled with file names that can be used 
    # For the aggregated statistics
    # (ensmean, enspctl)   
    file_string_models = ""
    m=models[1]
    for (m in models){
      print(m)      
      if(idx %in% pidx){
        
        # Fieldmean results
        cmd <- paste('cdo -O fldmean', ' ', land, '/', m,  ' ',
                     land, '/', 'fldm_', m , sep = "")
        print(cmd)
        system(cmd)
        
        ## add the preprocessed file to the filestring        
        file_string_models <- paste(file_string_models, land,
                                    '/fldm_', m, " ", sep = "")
        
      }else{
        # Subtracting timemeans from land files:
        cmd <- paste('cdo -O sub', ' ', land, '/', m, ' ',  land, '/tm_', m, ' ',
                     land, '/', 'norm_', m,  sep = "")
        print(cmd)
        system(cmd)
        
        # Detrended results:
        cmd <- paste('cdo -O detrend', ' ', land, '/norm_', m,  ' ', land, '/',
                     'detrend_', m,  sep = "")
        print(cmd)
        system(cmd)
        
        # Timstd of detrend
        cmd <- paste('cdo -O timstd', ' ', land, '/detrend_', m,  ' ', land, '/',
                     'detrend_std_', m,  sep = "")
        print(cmd)
        system(cmd)
        
        # Divide normalized by timstded detrend
        cmd <- paste('cdo -O div', ' ', land, '/norm_', m,  ' ', land, '/', 'detrend_std_',
                     m , ' ', land, '/', 'detrend_standard_', m, sep = "")
        print(cmd)
        system(cmd)
        
        # Fieldmean results
        cmd <- paste('cdo -O fldmean', ' ', land, '/detrend_standard_', m,  ' ',
                     land, '/', 'detrend_std_fldm_', m , sep = "")
        print(cmd)
        system(cmd)
        
        ## add the preprocessed file to the filestring        
        file_string_models <- paste(file_string_models, land,
                                    '/detrend_std_fldm_', m, " ", sep = "")
      }
      
    }
    #Find model ensemble mean
    cmd <- paste("cdo -O ensmean ", file_string_models, plot_dir, "/", idx,
                 "_ensmean_for_timeseries.nc", sep = "")
    print(cmd)
    system(cmd)
    
    #Find ensemble 25th percentile
    cmd <- paste("cdo -O enspctl,25 ", file_string_models, plot_dir, "/", idx,
                 "_25enspctl_for_timeseries.nc", sep = "")
    print(cmd)
    system(cmd)
    
    #Find ensemble 75th percentile
    cmd <- paste("cdo -O enspctl,75 ", file_string_models, plot_dir, "/", idx,
                 "_75enspctl_for_timeseries.nc", sep = "")
    print(cmd)
    system(cmd)
    

    n <- 0
    for (o in obs){
      print(o)
      
      if(idx %in% pidx){
        # Fieldmean results
        cmd <- paste('cdo -O fldmean', ' ', land, '/', o,  ' ',
                     land, '/', 'fldm_', o , sep = "")
        print(cmd)
        system(cmd)
        
        # Copy obs file to plot
        n  <- n+1
        cmd <- paste("cp ", land, "/fldm_", o, " ", plot_dir, "/", idx, "_",
                     modelsAndObsSplit[obs_order[n]], "_for_timeseries.nc", sep = "")
        print(cmd)
        system(cmd)
        
      }else{
        # Subtracting timemeans from land files:
        cmd <- paste('cdo -O sub', ' ', land, '/', o, ' ',  land, '/tm_', o, ' ',
                     land, '/', 'norm_', o,  sep = "")
        print(cmd)
        system(cmd)
        
        # Detrended results:
        cmd <- paste('cdo -O detrend', ' ', land, '/norm_', o,  ' ', land, '/',
                     'detrend_', o,  sep = "")
        print(cmd)
        system(cmd)
        
        # Timstd of detrend
        cmd <- paste('cdo -O timstd', ' ', land, '/detrend_', o,  ' ', land, '/',
                     'detrend_std_', o,  sep = "")
        print(cmd)
        system(cmd)
        
        # Divide normalized by timstded detrend
        cmd <- paste('cdo -O div', ' ', land, '/norm_', o,  ' ', land, '/', 'detrend_std_',
                     o , ' ', land, '/', 'detrend_standard_', o, sep = "")
        print(cmd)
        system(cmd)
        
        # Fieldmean results
        cmd <- paste('cdo -O fldmean', ' ', land, '/detrend_standard_', o,  ' ',
                     land, '/', 'detrend_std_fldm_', o , sep = "")
        print(cmd)
        system(cmd)
        
        # Copy obs file to plot
        n  <- n+1
        cmd <- paste("cp ", land, "/detrend_std_fldm_", o, " ", plot_dir, "/", idx, "_",
                     modelsAndObsSplit[obs_order[n]], "_for_timeseries.nc", sep = "")
        print(cmd)
        system(cmd)
      }
    }
  }
  
  #### ABOSOLUTE VALUES ####
  ## Non-normalized values fieldmeans
  if(!normalize){
    file_string_models = ""
    m=models[1]
    for (m in models){
      print(m)
      # Fieldmean results
      cmd <- paste('cdo -O fldmean', ' ', land, '/', m,  ' ',
                   land, '/', 'fldm_', m , sep = "")
      print(cmd)
      system(cmd)
      
      ## add the preprocessed file to the filestring        
      file_string_models <- paste(file_string_models, land,
                                  '/fldm_', m, " ", sep = "")
    }
    #Find model ensemble mean
    cmd <- paste("cdo -O ensmean ", file_string_models, plot_dir, "/", idx,
                 "_ensmean_for_timeseries.nc", sep = "")
    print(cmd)
    system(cmd)
    
    #Find ensemble 25th percentile
    cmd <- paste("cdo -O enspctl,25 ", file_string_models, plot_dir, "/", idx,
                 "_25enspctl_for_timeseries.nc", sep = "")
    print(cmd)
    system(cmd)
    
    #Find ensemble 75th percentile
    cmd <- paste("cdo -O enspctl,75 ", file_string_models, plot_dir, "/", idx,
                 "_75enspctl_for_timeseries.nc", sep = "")
    print(cmd)
    system(cmd)
    
    ## Extracting only the observation files
    obs_order <- which(modelsAndObsSplit %in% obs_list)
    obs <- modelsAndObs[obs_order]
    
    print("These are the observations:")
    print(obs)
    n <- 0
    for (o in obs){
      print(o)
      # Fieldmean results
      cmd <- paste('cdo -O fldmean', ' ', land, '/', o,  ' ',
                   land, '/', 'fldm_', o , sep = "")
      print(cmd)
      system(cmd)
        
      # Copy obs file to plot
      n  <- n+1
      cmd <- paste("cp ", land, "/fldm_", o, " ", plot_dir, "/", idx, "_",
                   modelsAndObsSplit[obs_order[n]], "_for_timeseries.nc", sep = "")
      print(cmd)
      system(cmd)
    }
  }
}




##
##
## Method to plot the time series plot
## of single idx for already preprocessed data
## yearly data is assumed
## @param path - path to directory containing ensemble mean
## and percentile data.
## @param idx name of index to be processed
## @param start_yr start year for data to be used to convert
## values from days after start year format to
## year.
##

timeseries_plot <- function(plot_dir = "./plot", idx = "tn10pETCCDI_yr", obs_list, start_yr = 2006, normalize = FALSE){
    ## Loading ncdf4 library
    library(ncdf4)

    ## Loading scales library (required for transpareny in plots)
    library(scales)
    
    ## Loading in climdex dataframe
    source('nml/cfg_ExtremeEvents/cfg_climdex.r')
    
    
    # Drawing parameters
    #col_list <- c("dodgerblue2", "darkgreen", "firebrick2")
    #lty_list <- c(1, 4, 2)
    #lwd_list <- c(2, 2, 2)
    leg_names <- c(CMIP_name, obs_list)
    
    ## Reading the netcdf data files into R
    ## First ensemble mean file
    ensm <- nc_open(paste(plot_dir, "/", idx, "_ensmean_for_timeseries.nc", sep = ""))
    #ensm <- nc_open("/net/pdo/div/pdo/extreme/chriswm/plots/ExtremeEvents/tn10pETCCDI_yr_ensmean_for_timeseries.nc")  
   
    ## Then 25th percentile file
    enspctl25 <- nc_open(paste(plot_dir, "/", idx, "_25enspctl_for_timeseries.nc", sep = ""))
    ## Finally 75th percentile file    
    enspctl75 <- nc_open(paste(plot_dir, "/", idx, "_75enspctl_for_timeseries.nc", sep = ""))

    ## Reading in time variable and converting to years:
    ts  <- nc.get.time.series(ensm)  # from ncdf4.helpers
    time_conv <- format(ts, "%Y")    # extract years

    ## Stripping off the _yr tail to the index name
    #idx_name <- strsplit(idx, "_")[[1]][1]
    idx_no <- which(idx_df$idxETCCDI_time == idx)
    idx_name <- paste(idx_df$idxETCCDI[idx_no], "ETCCDI", sep="")
    
    ## Reading in the y-variables to be plotted
    ## First the ensemble mean
    idx_ensm <- ncvar_get(ensm, idx_name)
    ## Then the 25th percentile
    idx_ens25 <- ncvar_get(enspctl25, idx_name)
    ## Finally the 75th percentile
    idx_ens75 <- ncvar_get(enspctl75, idx_name)
    
    
    ## Maximum and minimum x and y values
    max.x <- max(time_conv)
    min.x <- min(time_conv)
    max.y <- max(idx_ensm, idx_ens25, idx_ens75)
    min.y <- min(idx_ensm, idx_ens25, idx_ens75)

    
    ## Reading in the observations and plotting via a loop
    obsdata_list <- list()
    n=0
    for(o in obs_list){
      n=n+1
      nc_obs <- nc_open(paste(plot_dir, "/", idx, "_", o, "_for_timeseries.nc", sep = ""))
      ts_obs  <- nc.get.time.series(nc_obs)  # from ncdf4.helpers
      time_conv_obs <- format(ts_obs, "%Y")    # extract years
      idx_obs <- ncvar_get(nc_obs, idx_name)
      nc_close(nc_obs)
      obsdata_list[[n]] <- list(o, as.numeric(time_conv_obs), idx_obs)
      max.x <- max(max.x, time_conv_obs)
      min.x <- min(min.x, time_conv_obs)
      max.y <- max(max.y, idx_obs)
      min.y <- min(min.y, idx_obs)
      
      if(n>length(ts_col_list)){
        print("Error: There are more observations, than available color plotting parameters.")
        print("Update cfg_ExtermeEvents.r file.")
        dev.off()
        break()
      }
      #lines(time_conv_obs, idx_obs, col = ts_col_list[n], lty= ts_lty_list[n], lwd=ts_lwd_list[n]) # plot observation
    }
    
    
    ## Setting the x- and y-range limits for plotting
    xrng <- as.numeric(c(min.x, max.x))
    yrng <- c(min.y, max.y)
    print(xrng)
    print(yrng)
    #diff.yrng <- abs(diff(yrng))
    #yrng[1] <- yrng[1] - diff.yrng*0.5
    #yrng[2] <- yrng[2] + diff.yrng*0.5
    
    ## Making name string for the plot
    plotname <- paste(plot_dir, "/", idx,"_", length(obs_list),"-obs_ensmean_timeseriesplot", sep = "")

    ## Setting device to write the plot to
    figure_filename <- paste(plotname, output_file_type, sep = ".")

    ## Chose output format for figure
    if (tolower(output_file_type) == "png") {
        png(filename = figure_filename,
        width = ts_png_width,
        height = ts_png_height,
        units = ts_png_units,
        pointsize = ts_png_pointsize,
        bg = ts_png_bg) 
    } else if (tolower(output_file_type) == "pdf") {
        pdf(file <- figure_filename)
    } else if (tolower(output_file_type) == "eps") {
        setEPS()
        postscript(figure_filename)
    }
    
    n=1
    ## Parameters for plot
    par(mfrow = c(1,1), mar = c(4.5, 4.5, 2, 3))
    ## Plotting first the ensemblemean
    plot(time_conv, idx_ensm, type = "l", col = ts_col_list[n], lty= ts_lty_list[n], xlim = xrng, ylim=yrng, lwd = ts_lwd_list[n], ann = FALSE, xaxs = "i", yaxt="n")
    ## Then making a transparent polygon between the 25th and 75 percentile
    polygon(c(time_conv,rev(time_conv)), c(idx_ens75, rev(idx_ens25)), col = alpha(ts_col_list[n],0.1),
            border = NA)
    
    ## Plotting observations and plotting via a loop
    n=0
    for(o in obs_list){
      n=n+1
      print(obsdata_list[[n]])
      lines(obsdata_list[[n]][[2]], obsdata_list[[n]][[3]], col = ts_col_list[n+1], lty= ts_lty_list[n+1], lwd=ts_lwd_list[n+1]) # plot observation
    }
    
    ## Produce a legend
    legend("top", legend = leg_names, col = ts_col_list, lty = ts_lty_list, lwd=ts_lwd_list, bty="n", ncol = 3)
    
    
    ## Produce a first y-axis
    axis(side = 2, at = pretty(yrng, 5))
    axis(side = 2, at = pretty(yrng, 5))
    pretty(yrng, 10)
    
    
    axis(side = 2, at = pretty(yrng, 5))
    
    ## Produce a second y-axis
    axis(side = 4, at = pretty(yrng, 5))
    
    ## Producing a title from info in netcdf file
    #title(main = ensm$var[[which(names(ensm$var) == idx_name, arr.ind = TRUE)]]$longname, font.main = 2)
    title(main = idx_df$name[idx_no], font.main = 2)
    
    ## Choosing x-label
    title(xlab = "Year")
  
    ## Chosing y-label from idx_ylab list
    title(ylab = idx_ylab[idx_no])
    #title(ylab = paste("Ensemble mean", ensm$var[[which(names(ensm$var) == idx_name, arr.ind = TRUE)]]$units))
    ## Resetting plotting device to default
    dev.off()

    # TO DO: add corresponding metadata to plot:
    #        - list of variables
    #        - list of input files
    tags = union(tags, c("PT_times", "ST_extreme", "DM_global"))
    caption = paste("Extreme events (", idx, ")", sep = "")
    id = paste("id", diag_base, idx, sep = "_")
    contrib_authors <- c("A_sand_ma", "A_broe_bj", "A_laue_ax")
    var <- ncatt_get(ensm, 0, "variable")
#    all_global_atts <- ncatt_get(ensm, 0)
    ESMValMD(figure_filename, tags, caption, id, var, models_name, fullpath_filenames,
             diag_script, contrib_authors)

## Close Ensemble files
    nc_close(ensm)
    nc_close(enspctl25)
    nc_close(enspctl75)

}


