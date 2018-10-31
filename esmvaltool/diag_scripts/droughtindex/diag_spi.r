library(yaml)
library(ncdf4)
library(SPEI)

getnc <- function(yml, m) {
  id <- nc_open(yml[m][[1]]$filename, readunlim=FALSE)
  v <- ncvar_get(id, yml[m][[1]]$short_name)
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
var1_input <- read_yaml(params$input_files[1])
nmods <- length(names(var1_input))

histbrks <- c(-99999, -2, -1.5, -1, 1, 1.5, 2, 99999)
histnams <- c("Extremely dry", "Moderately dry", "Dry",
              "Neutral",
              "Wet", "Moderately wet", "Extremely wet")

print(names(var1_input))
for(n in 1:nmods){
   if(var1_input[n][[1]]$cmor_table == "OBS"){
      ref <- getnc(var1_input, n)
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
   print(summary(as.numeric(v1_spi)))
   v1_spi[is.infinite(v1_spi)] = NA
   print(summary(as.numeric(v1_spi)))
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
   # Should weight against latitude!
   print(summary(as.numeric(v1_spi)))
   print(histbrks)
   histarr[mod,] <- hist(v1_spi, breaks=histbrks,
                         plot=FALSE)$density
}#mod

save(histarr, file='/home/a001850/tmp/histarr.rsav')

cols <- c("black","red") #,"green")

png("/home/a001850/tmp/tst.png", width=1000, height=500)
 barplot(histarr, beside=1, names.arg=histnams,
         col=cols)
 mtext("Probability", side=2, line=2.1)

 #plot(histarr[nref,], pch=8,
 #     font=2, cex=2, axes=FALSE,
 #     xlab="", ylab="")
 #axis(1,labels=histnams, at=1:7)
 #axis(2)
 #box()
 #if(nmods>1){
 #   pch=1
 #   for(n in 1:nmods){
 #     if(n != nref){
 #       points(histarr[n,], pch=pch, col="red")
 #     }
 #   }
 #}
dev.off()
