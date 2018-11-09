library(yaml)
library(ncdf4)
library(SPEI)
library(RColorBrewer)

getnc <- function(yml, m, lat=FALSE) {
  id <- nc_open(yml[m][[1]]$filename, readunlim=FALSE)
  if(lat){
    v <- ncvar_get(id, 'lat')
  }else{
    v <- ncvar_get(id, yml[m][[1]]$short_name)
  }
  nc_close(id)
  return(v)
}
ncwrite <- function(yml, m, data, wdir){
  fnam <- strsplit(yml[m][[1]]$filename,"/")[[1]]
  pcs <- strsplit(fnam[length(fnam)],"_")[[1]]
  pcs[which(pcs==yml[m][[1]]$short_name)] = "spi"
  onam <- paste0(wdir,"/",paste(pcs,collapse="_"))
  system(paste0("cp ", yml[m][[1]]$filename, " ",onam))
  idw <- nc_open(onam, write=TRUE)
    ncvar_put(idw, yml[m][[1]]$short_name, data)
  nc_close(idw)
}

args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
wdir <- params$work_dir
dir.create(wdir, recursive = TRUE)
pdir <- params$plot_dir
dir.create(pdir, recursive = TRUE)
var1_input <- read_yaml(params$input_files[1])
nmods <- length(names(var1_input))

histbrks <- c(-99999, -2, -1.5, -1, 1, 1.5, 2, 99999)
histnams <- c("Extremely dry", "Moderately dry", "Dry",
              "Neutral",
              "Wet", "Moderately wet", "Extremely wet")

print(names(var1_input))
for(n in 1:nmods){
   if(var1_input[n][[1]]$cmor_table == "OBS"){
      lat <- getnc(var1_input, n, lat=TRUE)
      if(max(lat)>90){
        print(paste0("Latitude must be [-90,90]: min=",
              min(lat), " max=", max(lat)))
        stop("Aborting!")
      }
      ref <- getnc(var1_input, n, lat=FALSE)
      refmsk <- apply(ref,c(1,2),FUN=mean)
      refmsk[refmsk > 10000] = NA
      refmsk[!is.na(refmsk)]=1
      nref <- n
   }
}

histarr <- array(NA, c(nmods, length(histnams)))
for(mod in 1:nmods){
   v1 <- getnc(var1_input, mod)
   print(var1_input[mod][[1]]$cmor_table)
   d <- dim(v1)
   v1_spi <- v1*NA
   for(i in 1:d[1]){
    tmp <- v1[i,,]
    v1_spi[i,,] <- t(spi(t(tmp), 1, na.rm=TRUE,
                     distribution='PearsonIII')$fitted)
   }
   v1_spi[is.infinite(v1_spi)] = NA
   v1_spi[v1_spi > 10000] = NA
   print(summary(as.numeric(v1_spi)))
   ncwrite(var1_input, mod, v1_spi, wdir)
   for(t in 1:d[3]){
      tmp <- v1_spi[,,t]
      tmp[is.na(refmsk)] = NA
      v1_spi[,,t] <- tmp
   }#t
   v1_spi[is.infinite(v1_spi)] = NA
   v1_spi[v1_spi > 10000] = NA
   # Weight against latitude
   h <- c(1:length(histnams))*0
   for(j in 1:d[2]){
     h <- h + hist(v1_spi[j,,], breaks=histbrks,
                   plot=FALSE)$counts * cos(lat[j]*pi/180.)
   }#j
   #print(summary(as.numeric(v1_spi)))
   #print(histbrks)
   #h <- hist(v1_spi, breaks=histbrks,
   #                      plot=FALSE)$density
   histarr[mod,] <- h/sum(h)
}#mod
save(histarr, file=paste0(params$work_dir,
                          "/histarr.rsav"))

bhistarr <- array(NA,c(nmods-1,7))
marr <- c(1:nmods)[c(1:nmods) != nref]
cnt <- 1
for(m in marr){
  bhistarr[cnt,] <- histarr[m,]-histarr[nref,]
  cnt <- cnt + 1
}
parr <- c(nref,marr)

mnam <- c(1:nmods)*NA
for(m in 1:nmods){mnam[m] <- var1_input[m][[1]]$dataset}

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                           rownames(qual_col_pals)))
cols=c("black",sample(col_vector, nmods-1))

png(paste0(params$plot_dir,"/histplot.png"),
    width=1000, height=500)
 par(mfrow=c(2,1),oma=c(3,3,3,13),mar=c(2,1,1,1))
 barplot(histarr[parr,], beside=1, names.arg=histnams,
         col=cols, xaxs="i")
 box()
 mtext("Probability", side=2, line=2.1)
 barplot(bhistarr, beside=1, names.arg=histnams,
         col=cols[2:nmods], xaxs="i")
 box()
 mtext("Absolute difference", side=2, line=2.1)
 mtext("Standardized precipitation index", outer=TRUE, cex=2, font=2)
 par(fig=c(0.8,.95,0.1,0.9),new=T,oma=c(1,1,1,1)*0,mar=c(0,0,0,0))
 legend("topright", mnam[parr], fill=cols)
dev.off()
