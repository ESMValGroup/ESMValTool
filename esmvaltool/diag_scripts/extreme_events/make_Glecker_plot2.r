# #########################################################################################################
# make_Glecker_plot2.r
#
# Author: Christian W. Mohr (CICERO, Norway)
#         Marit Sandstad (CICERO, Norway)
#
#
# #########################################################################################################
# Description:
# Code to plot Glecker polygon diagram to compare climdex index
#	performance between models and reanalysis.	
#
# Modification history
#
#    20180601-A_cwmohr: re-creation (complete new script incorparating segments from "make_timeseries_plot.r & make_Glecker_plot.r")
#
# #########################################################################################################

source('diag_scripts/aux/ExtremeEvents/common_climdex_preprocessing_for_plots.r')


##
## Method that preprocesses idx-file for a single file
## in order to get the data to plot in a Glecker diagram
## specific to this index
## @param path is the path to the index file
## @param idx is the index to be prepocessed
##
#Glecker_preprocessing_main <- function(path = '../../climdex/ensMeanMed/rcp26', idx_list = c('sdii', 'r95p', 'rxday', 'cdd', 'tr', 'fd', 'txn', 'txx', 'tnn', 'tnx')){
#  for (idx in idx_list){
#    model_paths <- list.files(path, pattern = paste(idx, "*.nc", sep = ""), full.names = TRUE)
#    models <- list.files(path, pattern = paste(idx, "*.nc", sep = ""))  	 
#    for (model_idx in c(1:length(models))){
#      testRegridAndLandSeaMask(model_paths[model_idx])
#      
#    }
#  }
#}




gleckler_main <- function(path = "./", idx_list = c("tn10pETCCDI_yr", "tn90pETCCDI_yr", "tx10pETCCDI_yr", "tx90pETCCDI_yr", 'tnnETCCDI_yr', 'tnxETCCDI_yr',
                                                    'txnETCCDI_yr', 'txxETCCDI_yr'),
                          model_list = c("EC-EARTH", "HadGEM2-ES", "MPI-ESM-LR", "IPSL-CM5A-LR", "MIROC-ESM"),
                          obs_list = c("CanESM2"), plot_dir = "../plot/ExtremeEvents/", promptInput=promptInput){
  
  
  #### CLIMDEX PREPROCESSING ####
  
  ## For file structure and files
  tsGrid <-  paste(path, "/tsGridDef", sep = "")
  time_cropped <- paste(path, "/timeCropped", sep = "")
  landmask <- paste(path, "/landSeaMask.nc", sep = "")
  regridded <- paste(path, "/regridded", sep = "")
  land <- paste(path, "/Land", sep = "")
  
  
  nmodel = length(model_list) # number of models
  nidx = length(idx_list) # number of indices
  nobs = length(obs_list) # number of observations
  
  
  
  if(promptInput=="y"){
    ## Initial nc-file time crop, regrid, land and plot purge
    cmd <- paste("rm -f ", time_cropped, "/*.nc", " ", regridded, "/*.nc", " ", land, "/*.nc", sep = "")
    print(cmd)
    system(cmd)
    
    ## Initial grid and landmask creation reset
    gridAndLandmask = TRUE
    
    ## Combine model and observation list
    modelAndObs_list <- unique(c(model_list, obs_list))
    
    ## Loop over the indices to produce a plot for each index
    #idx = idx_list[1]
    for (idx in idx_list){
      
      ## Time crop
      returnvalue <- setTimeForFilesEqual(path = path, idx = idx, model_list = modelAndObs_list, time_cropped = time_cropped)
      
      max_start <- returnvalue[1]
      min_end <- returnvalue[2]
      
      ## If there is no overlap in the files the index
      ## should be skipped
      if (max_start >= min_end){
        print(paste("No time overlap in files for index", idx))
        break()
      }
      
      ## Find the new model and observation names (after time cropping)
      modelsAndObs <-  basename(Sys.glob(file.path(time_cropped, paste(idx, "*.nc", sep = ""))))
      split_modelsAndObs <- strsplit(modelsAndObs, split = "_")
      modelsAndObs_index <- unlist(lapply(split_modelsAndObs, function(x){x[3]}))
      
      ## new models
      models <- modelsAndObs[which(modelsAndObs_index %in% model_list)]
      
      ## new observations
      obs <- modelsAndObs[which(modelsAndObs_index %in% obs_list)]
      
      ## Find the start year (to be used in plotting)
      start_yr <- strtoi(substr(models[1], nchar(models[1]) - 11, nchar(models[1]) - 8))
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
      for (mo in modelsAndObs){
        print(paste(time_cropped, "/", mo, sep = ""))
        regridAndLandSeaMask(idx_raw =paste(time_cropped, "/", mo, sep = ""), regrid = tsGrid,
                             landmask = landmask, regridded = regridded, land = land, loc = path)
      }
      
    }
    
    #### Gleckler Array Processing ####
    RMSErelarr <- gleckler_array(path = land, idx_list=idx_list, model_list=model_list, obs_list=obs_list)
    
    ## Save Array 
    saveRDS(object = RMSErelarr, file = paste0(plot_dir, "/Gleclker-Array_", nidx, "-idx_", nmodel,"-models_", nobs, "-obs",  ".RDS"))
    saveRDS(object = returnvalue, file = paste0(plot_dir, "/Gleclker-years.RDS"))
  }
  
  
  
  #### Gleckler Plotting ####
  RMSErelarr <- readRDS(file = paste0(plot_dir, "/Gleclker-Array_", nidx, "-idx_", nmodel,"-models_", nobs, "-obs",  ".RDS"))
  year_range <- readRDS(file = paste0(plot_dir, "/Gleclker-years.RDS"))
  gleckler_plotting(arr = RMSErelarr, idx_list=idx_list, model_list=model_list, obs_list=obs_list, plot_dir = plot_dir, syear=year_range[1], eyear=year_range[2])
  
}






#### Computing the RMSEs ####


gleckler_array <- function(path = land, idx_list=gleckler_idx, model_list=model_list, obs_list=obs_list){
  
  library(ncdf4)
  

  ## Produce an array to hold all the model and reanalysis means
  
  
  ## Input data for testing the plotting routine
  nidx <- length(idx_list) # number of indices
  nmodel <- length(model_list) # number of models
  nobs <- length(obs_list) # number of reanalyses
  
  
  ## Check point for reanalysis data
  if(nobs == 0){
    print("No reanalysis datasets provided")
    break()
  }
  
    
  ## Function to calculate area mean
  area.mean <- function(x, lat){
    #x <- tm_model_idx
    #x <- matrix(1:9, nrow = 3, ncol = 3)
    
    #x[2,2] <- NA
    
    nlon <- dim(x)[1]
    nlat <- dim(x)[2]
      
    meanlat <- apply(x, 2, function(x){mean(x, na.rm = TRUE)})
    
    #model_lat_b <- seq(90, -90, length.out = nlat+1)
    #lat <- model_lat_b[-(nlat+1)]+diff(model_lat_b)/2
    #lat <- model_lat
    fi <- lat*3.14159/180
    
    
    wgt.prod <- meanlat*cos(fi)
    
    
    nan.check <- is.nan(wgt.prod) # At some latitudes there is no land and therfore no data. 
                                  # The mean of missing data is not a number, and hench results in NaNs. 
                                  # These NaN must be removed in order to calculate the correct area mean.
    
    gl <- sum(wgt.prod[!nan.check])
    sumcos <- sum(cos(fi)[!nan.check])
    ar.m <- gl/sumcos
    return(ar.m)
  }
    
  
  ## Function to calculate the RMSE between the model and observed climatology (RMSExy)
  ## Equation 1, from Sillmann et. al 2013
  RMSE <- function(model=tm_model_idx, obs=tm_obs_idx, lat = model_lat){
    RMSE <- sqrt(area.mean((model-obs)^2, lat))
    return(RMSE)
  }
  
  RMSEarr <- array(NA, dim=c(nidx + 1, nmodel + 3, nobs)) # Array for the RMSE spaces in the array are created so that the RSMEall, ENSmean, ENSmedian and CMIP RMSE can be created
  RMSErelarr <- RMSEarr
  ensmodel_list <- list()

  i=2
  m=1
  o=1
  lat_collect <- TRUE
  
  for(i in seq_along(idx_list)){
    for(m in seq_along(model_list)){
      ## Read in model annual climatology
      tm_model <- nc_open(Sys.glob(file.path(path, paste("tm_",idx_list[i], "_", model_list[m],"*.nc", sep = ""))))
      idxs <- unlist(strsplit(idx_list[i], split = "_"))[1]
      tm_model_idx <- ncvar_get(tm_model, idxs)
      
      if(lat_collect){
        model_lat <- ncvar_get(tm_model, "lat") # extract latitudes for area mean calculations
        lat_collect <- FALSE
      }
      
      nc_close(tm_model)
      ensmodel_list[[m]] <- tm_model_idx
    }
    ## Create a new array for adding the time mean model matices
    ensarr <- array(NA, dim=c(nrow(tm_model_idx), ncol(tm_model_idx), length(ensmodel_list)+2))
    
    ## Copy each matrix from the multimodel list to the array "ensarr". 
    ## Notice the "+2" on the 3rd dimention. This is so later the model ensemble mean and median matrices can be added to the array.
    for(n in seq_along(ensmodel_list)){
      ensarr[,,n+2] <- ensmodel_list[[n]]
    }
    
        
    ## Calculate the ensemble mean and median of all the model time mean matrices
    ensMean <- apply(ensarr, c(1,2), function(x){mean(na.omit(x))})
    ensMedian <- apply(ensarr, c(1,2), function(x){median(na.omit(x))})
    
    ## Place the ensemble model mean and medians into the first two matrices (3-dimention) of the array "ensarr" 
    ensarr[,,1] <- ensMean
    ensarr[,,2] <- ensMedian
    
    j=1
    ## Calculate the RMSE for all the models and the ensemble mean and median
    for(j in 1:dim(ensarr)[3]){
      ## Read in reannalysis annual climatology
      for(o in seq_along(obs_list)){
        tm_obs <- nc_open(Sys.glob(file.path(path, paste("tm_",idx_list[i], "_", obs_list[o],"*.nc", sep = ""))))
        tm_obs_idx <- ncvar_get(tm_obs, idxs)
        nc_close(tm_obs)
        RMSEarr[i+1,j, o] <- RMSE(model = ensarr[,,j], obs = tm_obs_idx, lat = model_lat) # Calculate each RMSE and place value in RMSE-array
        
        ## Calculate the model standard deviation. Later used for calculating the RMSEmedian,std.
        ## Denominator in equation 3, from Sillmann et. al 2013
        RMSEarr[i+1,ncol(RMSEarr), o] <- sqrt(area.mean((tm_obs_idx - area.mean(tm_obs_idx, lat = model_lat))^2, lat = model_lat))
      }
    }
  }
  
  ## Calculate the RMSE median for the models
  tmpRMSEarr <- RMSEarr[,-c(1,2, ncol(RMSEarr)),]
  if(length(dim(tmpRMSEarr)) == 3){
    RMSEmed <- apply(tmpRMSEarr, c(1,3), function(x){median(x, na.rm = TRUE)})
  }else{
    RMSEmed <- apply(tmpRMSEarr, 1, function(x){median(x, na.rm = TRUE)})
  }
  
  
  
  ## Function to calculate the relatvie RMSE (RMSE'xy) between the model and observed climatology
  ## Equation 2, from Sillmann et. al 2013
  RMSErel <- function(RMSE, RMSEmed){
    RMSErel <- (RMSE - RMSEmed)/RMSEmed
    return(RMSErel)
  }
  
  ## Calculating the relative RMSE (RMSE'xy)
  m=1
  for(m in 1:(ncol(RMSEarr)-1)){
    RMSErelarr[,m,] <- RMSErel(RMSE = RMSEarr[,m,], RMSEmed = RMSEmed)
  }
  
  ## Calculating the RMSE median,std. Equation 3, from Sillmann et. al 2013
  RMSErelarr[,ncol(RMSErelarr),] <- RMSEmed/RMSEarr[,ncol(RMSEarr),]
  
  
  
  ## Calculating the RSME mean
  tmpRMSEarr <- RMSErelarr[,-ncol(RMSErelarr),]
  if(length(dim(tmpRMSEarr))==3){
    RMSErelarr[1,-ncol(RMSErelarr),] <- apply(tmpRMSEarr, c(2,3), function(x){mean(x, na.rm=TRUE)})
  }else{
    RMSErelarr[1,-ncol(RMSErelarr),] <- apply(tmpRMSEarr, c(2), function(x){mean(x, na.rm=TRUE)})
  }
  print(RMSErelarr)
  return(RMSErelarr)
}


#### Plotting Routine ####
gleckler_plotting <- function(arr = RMSErelarr, idx_list, model_list, obs_list, plot_dir = "../plots/ExtremeEvents/", syear=max_start, eyear=min_end){
  
  
  nidx <- length(idx_list) # number of indices
  nmodel <- length(model_list) # number of models
  nobs <- length(obs_list) # number of reanalyses
  

  ## Numbers for color scale
  sclseq <- seq(-0.55, 0.55, 0.1)
  
  ## Colour scale
  library(RColorBrewer)
  glc <- brewer.pal(length(sclseq)-2, "RdYlBu")
  glc <- c("#662506", glc, "#3f007d")
  glc <- rev(glc)
  
  ## Numbers for black & white scale
  sclseq_bw <- seq(0.05, 1.15, 0.1);sclseq_bw
  glbw <- gray(seq(0,1, length.out = length(sclseq_bw)))
  glbw <- rev(glbw)
  
  
  ## Determinng what shapes should be plotted, based on number of observations
  if(nobs==1){
    # One reanalysis references
    x1 <- c(0, 1, 1, 0)
    y1 <- c(0, 0, 1, 1)
    xs <- list(x1)
    ys <- list(y1)
    
    # text coordinates
    xtx <- 0.50
    ytx <- -0.25
    rotx <- 0            # text rotation in degrees
  }
  
  if(nobs==2){
    # Two reanalysis references
    x1 <- c(0, 1, 1) # lower triangle
    y1 <- c(0, 0, 1) # lower triangle
    x2 <- c(0, 1, 0) # upper triangle
    y2 <- c(0, 1, 1) # upper triangle
    
    xs <- list(x1, x2)
    ys <- list(y1, y2)
    
    # text coordinates
    xtx <- c(0.75, 0.25)
    ytx <- c(-0.25, 1.25)
    rotx <- c(0, 0)            # text rotation in degrees
  }
  
  if(nobs==3){
    # Three reanalysis references
    x1 <- c(0, 0.5, 0.5, 0)         # bottom left
    y1 <- c(0, 0, 0.5, 1)           # bottom left
    x2 <- c(0.5, 1, 1, 0.5)         # bottom right
    y2 <- c(0, 0, 1, 0.5)           # bottom right
    x3 <- c(0, 0, 0.5, 1, 1)        # top
    y3 <- c(1, 0.75, 0.5, 0.75, 1)  # top
    
    xs <- list(x1, x2, x3)
    ys <- list(y1, y2, y3)
    
    # text coordinates
    xtx <- c(-0.25, 1.25, 0.5)
    ytx <- c(0.25, 0.25, 1.25)
    rotx <- c(90, 90, 0)            # text rotation in degrees
  }
  
  if(nobs==4){
    # Four reanalysis references
    x1 <- c(0, 0.5, 1) # bottom triangle
    y1 <- c(0, 0.5, 0) # bottom triangle
    x2 <- c(0, 0.5, 0) # left triangle
    y2 <- c(0, 0.5, 1) # left triangle
    x3 <- c(0, 0.5, 1) # top triangle
    y3 <- c(1, 0.5, 1) # top triangle
    x4 <- c(1, 0.5, 1) # right triangle
    y4 <- c(1, 0.5, 0) # right triangle
    
    xs <- list(x1, x2, x3, x4)
    ys <- list(y1, y2, y3, y4)
    
    # text coordinates
    xtx <- c(0.5, -0.25, 0.5, 1.25)
    ytx <- c(-0.25, 0.5, 1.25, 0.5)
    rotx <- c(0, 90, 0, 90)            # text rotation in degrees
  }
  
  
  if(!(nobs %in% c(1,2,3,4))){
    if(nobs == 0){
      print("No reanalysis dataset provided")
      break()
    }else{
      print("Too many reanalysis datasets provided. Please choose between 1 and 4 datasets")
      break()
    }
  }
  
  #plot_dir <- "../plots/ExtremeEvents/"
  print("--- Creating Gleckler plot ---")
  #res=2000
  #mar.par=c(7, 4, 3, 10)
  img.adj <- gl_mar.par*0.05
  #pt.fct <- 0.02
  #pts <- res*pt.fct
  width.fct <- ((nmodel+3)/(nidx+1))+sum(img.adj[c(2,4)])
  height.fct <- 1+sum(img.adj[c(1,3)])
    
  figure_filename <- paste(plot_dir,"/Gleckler_", CMIP_name, "_", nmodel,"-models_", nidx, "-idx_", nobs, "-obs_", syear, "-", eyear, ".",output_file_type, sep="")
  
  ## Chose output format for figure
  if (tolower(output_file_type) == "png") {
    png(filename = figure_filename, 
        width = gl_png_res*(width.fct/height.fct), 
        height = gl_png_res,
        units = gl_png_units, 
        pointsize = gl_png_pointsize, 
        bg = gl_png_bg)
  } else if (tolower(output_file_type) == "pdf") {
    pdf(file <- figure_filename)
  } else if (tolower(output_file_type) == "eps") {
    setEPS()
    postscript(figure_filename)
  }
  
    
  #jpeg(filename = imagefilename, width = gl_png_res*(width.fct/height.fct), height = gl_png_res, units = gl_png_units, pointsize = gl_png_pointsize, bg = gl_png_bg)
  #jpeg(filename = imagefilename, width = res*(width.fct/height.fct), height = res, units = "px", pointsize = pts)
  par(mfrow=c(1,1), mar=gl_mar.par, xpd=FALSE, oma=rep(0, 4))
  #par(mfrow=c(1,1), mar=mar.par, xpd=FALSE, oma=rep(0, 4))
  #gl_RMSEspacer <- 0.01
  plot(x = c(0,1+gl_RMSEspacer), y=c(0,1), type="n", ann=FALSE, xaxs = "i", yaxs = "i", bty="n", xaxt = "n", yaxt = "n")
  #plot(x = c(0,1+spacer), y=c(0,1), type="n", ann=TRUE, xaxs = "i", yaxs = "i", bty="n") #for testing purpose

  ## Array dimentions
  xn = ncol(arr)
  yn = nrow(arr)
  
  ## Testing array plotting
  xi = 1  #model
  yj = 2  #index
  zk = 1  #obs
  
  
  ## Plotting RMSE of models, ensemble mean and median and RSMEall
  for(xi in 1:(xn-1)){
    for(yj in 1:yn){
      for(zk in 1:nobs){
        polygon(x = (xs[[zk]]/xn)+((xi-1)/xn), y = (ys[[zk]]/yn)+((yn-yj)/yn), col=glc[which.min(abs(sclseq - arr[yj,xi,zk]))])
      }
    }
  }
    
  ## Plotting RMSE median standard diviation
  for(yj in 2:yn){
    for(zk in 1:nobs){
      polygon(x = (xs[[zk]]/xn)+((xn-1)/xn)+gl_RMSEspacer, y = (ys[[zk]]/yn)+((yn-yj)/yn), col=glbw[which.min(abs(sclseq_bw - arr[yj,xn,zk]))])
    }
  }
    
  
  ## Produce the borders for the Glecker plot
  par(xpd=TRUE)
  rect(xleft=0, ybottom=0, xright=(1-1/xn), ytop=(1-1/yn), density = NULL, angle = 45,
       col = NA, border = 1, lty = par("lty"), lwd = 4)
  rect(xleft=0, ybottom=(1-1/yn), xright=(1-1/xn), ytop=1, density = NULL, angle = 45,
       col = NA, border = 1, lty = par("lty"), lwd = 4)
  

  ## Scale for Gleckler plot
  gleckler_scale <- function(sclseq, glc, scaling_factor, text.scaling_factor, xscale_spacer){
    par(xpd=TRUE)
    ## Square legend
    sqrxs <- c(0, 1, 1, 0)
    sqrys <- c(0, 0, 1, 1)
    
    ## up-triangle legend
    utrixs <- c(0, 1, 0.5)
    utriys <- c(0, 0, 1)
    
    ## down-triangle legend
    dtrixs <- c(0.5, 1, 0)
    dtriys <- c(0, 1, 1)
    
    ## Legend number shifter
    seq_shift <- mean(diff(sclseq)/2) # Shifts the legend numbers so that they represent the border values
    
    ## y-scale spacer
    yscale_spacer <- (1-scaling_factor)/2

    #txt_spacer <- 10/100
    #plot(x = c(0,1+spacer), y=c(0,1), type="n", ann=TRUE, xaxs = "i", yaxs = "i", bty="n")
    exlen <- length(glc)
    for(a in 1:exlen){
      if(a == 1){
        xtmp <- scaling_factor*(dtrixs/exlen)+1+xscale_spacer
        ytmp <- scaling_factor*(dtriys/exlen+(a-1)/exlen)+yscale_spacer
        polygon(x = xtmp, y = ytmp, col=glc[a])
        text(x = max(xtmp), y = max(ytmp), round(sclseq[a]+seq_shift,1), cex=text.scaling_factor, pos = 4)
      }else if(a == exlen){
        xtmp <- scaling_factor*(utrixs/exlen)+1+xscale_spacer
        ytmp <- scaling_factor*(utriys/exlen+(a-1)/exlen)+yscale_spacer
        polygon(x = xtmp, y = ytmp, col=glc[a])
      }else{
        xtmp <- scaling_factor*(sqrxs/exlen)+1+xscale_spacer
        ytmp <- scaling_factor*(sqrys/exlen+(a-1)/exlen)+yscale_spacer
        polygon(x = xtmp, y = ytmp, col=glc[a])
        text(x = max(xtmp), y = max(ytmp), round(sclseq[a]+seq_shift,1), cex=text.scaling_factor, pos = 4)
      }
    }
  }

  ## Plot scales
  gleckler_scale(sclseq, glc, scaling_factor=gl_scaling_factor, text.scaling_factor = gl_text.scaling_factor, xscale_spacer=gl_xscale_spacer_RMSE)
  
  gleckler_scale(sclseq_bw, glbw, scaling_factor=gl_scaling_factor, text.scaling_factor = gl_text.scaling_factor, xscale_spacer=gl_xscale_spacer_RMSEstd)
  
  
  ## Plotting symbol legend
  exlen <- length(glc)
  xsym1 <- gl_scaling_factor*(0.5/exlen)+1+gl_xscale_spacer_RMSE
  exlen <- length(glbw)
  xsym2 <- gl_scaling_factor*(0.5/exlen)+1+gl_xscale_spacer_RMSEstd
  x.max_adj <- max(gl_symb.scaling_factor*(xs[[zk]]/xn))
  x.min_adj <- min(gl_symb.scaling_factor*(xs[[zk]]/xn))
  xmidadj <- (x.max_adj - x.min_adj)/2
  
  gl_symb.xshift <- (xsym1 + xsym2)/2 - xmidadj
  

  for(zk in 1:nobs){
    xsym <- gl_symb.scaling_factor*(xs[[zk]]/xn) + gl_symb.xshift
    ysym <- gl_symb.scaling_factor*(ys[[zk]]/yn)-gl_symb.yshift
    print(paste("xs:", xsym))
    print(paste("ys:", ysym))
    polygon(x = xsym, y = ysym, col="white", border=1)
    
    xtxsym <- gl_symb.scaling_factor*(xtx[[zk]]/xn)+gl_symb.xshift
    ytxsym <- gl_symb.scaling_factor*(ytx[[zk]]/yn)-gl_symb.yshift
    
    text(x = xtxsym, y = ytxsym, labels = obs_list[zk], adj = 0.5, cex = gl_text.symb.scaling_factor, srt=rotx[zk])
  }


  
  
  
  
  
  ## Label adjusting parameters
  axlabsize <- 0.8
  lineadj <- -0.5
  
  ## Add model labels
  col_names <- c("ENSMEAN", "ENSMEDIAN", model_list)
  xtcks1 <- seq((0.5/xn), ((xn-1)/xn), by = (1/xn))
  #axis(side = 1, at = xtcks1, labels = paste("model-", letters[1:length(xtcks1)], sep=""), las=2, cex.axis=axlabsize, tick = FALSE, line = lineadj)
  axis(side = 1, at = xtcks1, labels = col_names, las=2, cex.axis=axlabsize, tick = FALSE, line = lineadj)
  
  xtcks2 <- ((xn-1)/xn)+gl_RMSEspacer+(0.5/xn)
  axis(side = 1, at = xtcks2, labels = expression("RMSE"["std"]), las=2, cex.axis=axlabsize, tick = FALSE, line = lineadj)
  
  source("nml/cfg_ExtremeEvents/cfg_climdex.r")
  ## Add index labels
  row_names <- vector(mode = "character", length = length(idx_list))
  for(i in seq_along(idx_list)){
    row_names[i] <- idx_df$idxETCCDI[which(idx_df$idxETCCDI_time %in% idx_list[i])]
  }
  row_names <- rev(c(expression("RSME"["all"]), row_names))
  ytcks1 <- seq((1/yn)*0.5, 1, by = (1/yn))
  axis(side = 2, at = ytcks1, labels = row_names, las=2, cex.axis=axlabsize, tick = FALSE, line = lineadj)
    
  mtext(text = paste(CMIP_name," global land ", syear, "-", eyear, sep=""), side = 3, line = 1, font = 2, cex=1.1)
  
  dev.off() 
}











