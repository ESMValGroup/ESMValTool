# CCMVal ozone report48#
# John Scinocca and David Stephenson, October 17 2009
# John.scinocca@ec.gc.ca, D.B.Stephenson@exeter.ac.uk
#
# This code was developed using R version 2.6.2 (2008-02-08)
# It requires the extra R library, mgcv
#
#############################################################


####  control switches/variables
args = commandArgs(trailingOnly = TRUE)

output_dir=args[1]
bayear=as.integer(args[2])
nmodel=as.integer(args[3])
yearstart=as.integer(args[4])
yearend=as.integer(args[5])
ntime=as.integer(args[6])
fil=args[7]
npn=8+as.integer(nmodel)-1
name=args[8:npn]
npn1=8+as.integer(nmodel)+as.integer(nmodel)-1
npn=npn+1
ensemble=args[npn:npn1]
J <-nmodel
tdim=ntime
# number of years in each model time series
year.sub <- array(0,c(tdim,nmodel))
# arrays to hold IMT estimates and confidence intervals for each model
mod.fits <- array(0,c(tdim,nmodel))
err.fits <- array(0,c(tdim,nmodel,2))

# maximum number of ensembles
nem <- 10

# set up some arrarys to hold input data
year <- array(0,c(ntime,nem,nmodel))
data <- array(0,c(ntime,nem,nmodel))
ensemble <- array(1,c(nmodel))
tlen <- array(1,c(nem,nmodel))
#name <- array("",c(nmodel))
# Put all the model years in a big vector t and all the model
# response values in a big vector y

J <- nmodel
t <- NULL
y <- NULL
n <- NULL
for (nn in 1:nmodel){
  for (mm in 1:ensemble[nn]) {
    t <- c(t,year[1:tlen[mm,nn],mm,nn])
    y <- c(y,data[1:tlen[mm,nn],mm,nn])
    n <- c(n,tlen[mm,nn])
  }
}

# store model names in nom and number of ensembles for each
# model in K
nom <- NULL
K <- NULL
for (mm in 1:nmodel){
  nom <- c(nom,name[mm])
  K <- c(K,ensemble[mm])
}
# store minimum start year and maximum end year for each model
enstart <- array(yearend,c(nmodel))
enend  <- array(0,c(nmodel))
# np - maximum time series length
# models - long vector with model number at each data point
#          in time series
# runs - long vector with ensemble number at each data point
#          in time series
# il - total number of time series
# modelp,runp - like models, runs except each time series has
#               full length np
np<-ntime
models<-NULL
runs<-NULL
il<-0	# index indicating which of the 36 time series
modelp<-NULL
runp<-NULL
for (nn in 1:J)
{
for (mm in 1:K[nn])
{
il<-il+1
models<-c(models,rep(nn,n[il]))
runs<-c(runs,rep(mm,n[il]))
modelp<-c(modelp,rep(nn,np))
runp<-c(runp,rep(mm,np))
}}

pathin=output_dir
flin=paste(pathin,"RAW","_",fil,".ascii",sep="")
fleout=paste(pathin,"TSAM","_",fil,"_",bayear,sep="")

dat=read.table(flin)
names(dat)=c("ozone","year","model","run")
dat$model<-as.factor(dat$model)
dat$run<-as.factor(dat$run)
models<-dat[,3]
y<-dat[,1]
t<-dat[,2]
runs<-dat[,4]
###################################################
# Load up GAM library - note that it should be version
# mgcv 1.4-1 or later (safest to use latest version from
# www.r-project.org.

library(mgcv)

################################################################
# perform first gam call to obtain initial smooth IMT estimate of
# ensemble time series for each model

mod<-gam(ozone~model+s(year)+s(year,by=model),data=dat)
out<-predict(mod,dat,se.fit=T)


# calculate quantities for baseline adjustment in bayear
# modave - value of IMT estimate for each ensemble
# maa - average IMT estimate at bayear over all models
# ya - baseline adjusted time series data to bayear
modave <- array(0,c(J))
nel <- NULL
for (j in 1:J){
  nel=0
  for (k in 1:K[j]){
    
    tt<-dat$year[(dat$model==j) & (dat$run==k)]
    yy<-out$fit[(dat$model==j) & (dat$run==k)]
#   average IMT estimate over each model's ensemble members.
#   It turns out that the call to gam above calculates one
#   IMT estimate for each ensemble and outputs the same IMT estimate
#   for each member of an ensemble.  So this average over ensemble is
#   a bit redundant
    nn=0
    for (iy in tt){
      nn=nn+1
       if( iy == bayear ) {
        nel=nel+1
        modave[j] <- modave[j]+yy[nn]

      }
    }
    if ( nel == 0 ) {
      #print(k)
      #print(j)
      #print(tt)
    }
  }
  modave[j] <- modave[j]/nel
}

maa <- 0
for (j in 1:J){
  maa <- maa+modave[j]
}

maa <- maa/J
#print(maa)

cm <- modave[models]-maa


ya <- y-cm
gamout1 <- cbind(maa)
#print(fleout)
write.table(gamout1, file = paste(fleout,"_baseline",".ascii",sep=""),quote = FALSE,row.names = FALSE,col.names = FALSE)

# create new data frame of year, model, ensemble number
# to help with plotting
new1<-data.frame(year=rep(yearstart:yearend),modelp,runp)
names(new1)=c("year","model","run")
new1$model<-new1$model
new1$run<-as.factor(new1$run)
#print(new1)
############################################################

# plots raw data with initial IMT estimates 
out<-predict(mod,new1,se.fit=T)
# Put baseline adjusted time-series data in a data frame
dat<-data.frame(cbind(ya,t,models,runs))
names(dat)=c("ozone","year","model","run")
dat$model<-as.factor(dat$model)
dat$run<-as.factor(dat$run)
#print(dat)
##############################################################
# Fit GAM models to the data
#

# Fit one trend curve to all models: Y_jk(t)=f(t)+e
# Null model: y=c+g(t)+e
modnull<-gam(ozone~s(year),data=dat)
# Fit different trends to different models: Y_jk(t)=f(t)+g_j(t)+e
mod<-gam(ozone~model+s(year)+s(year,by=model),data=dat)
#####################################################
#####################################################
# Weighted fit that minimises Sum w*e2
# Note: Only needed if the residuals for different models have
# very different variances.
# Not sure if mean(w) need to be 1?
e<-residuals(mod)
s<-sd(e)
for (j in 1:J)
{
#w[j]<-1/var(e[dat$model==j])
w<-1/var(e[dat$model==j])
}
w<-J*w/sum(w)  # scaled so that mean(w)=1
m<-as.numeric(as.vector(dat$model))
ww<-w[m]
#modw<-gam(ozone~s(year)+s(year,by=model),data=dat,weights=ww)
#modw<-gam(ozone~s(year)+s(year,by=model),data=dat,weights=ww)
#####################################################

# Model summary

# Which is the best model? i.e. smallest AIC
anova(modnull,mod)
AIC(modnull,mod)
# plot out baseline adjusted time series data
new<-data.frame(year=rep(yearstart:yearend),modelp,runp)
names(new)=c("year","model","run")
new$model<-as.factor(new$model)
new$run<-as.factor(new$run)
# plot out baseline adjusted time series data and final IMT estimates
out<-predict(mod,new,se.fit=T)
##############################################################
# Model checking: Residual diagnostics
enull<-residuals(modnull)
e<-residuals(mod)
###########################################################
# Extract mean and standard error of trends for each model
# newy - IMT estimate
# sy - standard error for each model
# da - a vector of 1 where model data is present
# wy - quadratic taper weighting (eq. 9.21 and 9.22 of Appendix)
# ty - vector of years from 1950-2100
newy<-matrix(nrow=tdim,ncol=J)
sy<-matrix(nrow=tdim,ncol=J)
da<-matrix(0,nrow=tdim,ncol=J)
wy<-matrix(0,nrow=tdim,ncol=J)
out<-predict(mod,new,se.fit=T)
ty<-seq(min(new$year),max(new$year))
for (j in 1:J)
{
lmod<-as.vector((as.vector(new$model)==j)&(as.vector(new$run)==1))
newy[,j]<-out$fit[lmod]
sy[,j]<-out$se.fit[lmod]
yrs<-dat$year[dat$model==j]
da[yrs-yearstart+1,j]<-1
z<-(-1+2*(ty-min(yrs))/(max(yrs)-min(yrs)))
wy[,j]<-(1-(abs(z))^2)*da[,j]
}
# set IMT estimate and standard errors to NA where da=0 (i.e. where
# no data exists in original time series.
newy[da==0]<-NA
sy[da==0]<-NA
# squared standard error
sy2<-sy^2
########################################################################
# script below follows Apendix Section 9A.4
# all questions should be directed to
# David Stephenson, D.B.Stephenson@exeter.ac.uk
## Assume the model trends h_j(t) are normally distributed
# about the TRUE trend h(t) with variance la+sy_j2 where
# sy_j is the prediction error of model j trend estimate and
# la is a between model discrepancy parameter that is found
# by tuning the scaled residuals(h_j-h)/(la+sy_j2) to
# have unit variance.
#
# To do this we do the following steps:
# 1. Estimate h(t) using prediction error weighted sum of h_j
# 2. Find la iteratively so that scaled residuals have unit variance
# 3. Re-calculate h(t) using weights proportional to v_j/(la+sy_j2)
# where v_j are prior weights (metrics) specified for each model.
#
#

# Set up some weights for first estimate of h(t)
wt<-1/sy2                       # set up an (n x J) matrix for model weights
w<-apply(wt,1,sum,na.rm=T)       # sum of weights for each time
w[w==0]<-NA
wi<-1/w
wi[is.na(wi)]<-0
wj<-wi*wt		# normalize the weights so they sum to one over all models
wj[w==0,]<-NA


# Make a weighted sum of the trends
ny<-apply(wj*newy,1,sum,na.rm=T)
ny[is.na(w)]<-NA
ey2<-apply(wj*wj*sy2,1,sum,na.rm=T)
ey<-sqrt(ey2)
ry<-newy-rep(ny,J)


# Use the residuals to estimate the BMV
# Choose BMV so that standardised residuals have unit variance

vf <- function(x){
s2<-x+sy2;
s2[s2<0]<-0;
rys<-ry/sqrt(s2);
v<-var(as.vector(rys),na.rm=T);
(v-1)^2}
out<-nlm(vf, 0.9*var(as.vector(ry),na.rm=T))
la2<-out$estimate
la<-sqrt(la2)

sy2x<-la2+sy2		# Use la2 estimate to define a more realistic new variance

# Check if the scaled residuals are normally distributed
# with mean zero and variance of one.
v<-rep(1,J)            # all models equally weighted 


v<-v/sum(v)     # normalize v to have unit sum
selectedmodels<-seq(1,J)[v!=0]  # models with non-zero weight are selected
norm<-length(selectedmodels)

# Use the fancy weights with priors v for the different models
wt<-sy2x
wtx<-t(v*t(wy))
wtx<-matrix(as.vector(wtx)/as.vector(wt),ncol=J)
wtx[wy==0]<-NA
w<-apply(wtx,1,sum,na.rm=T)      # sum of weights for each time
w[w==0]<-NA
wi<-1/w
wi[is.na(wi)]<-0
wj<-wi*wtx
wj[w==0,]<-NA

# Make a new weighted sum of the trends
# and standard error on the weighted sum
ny<-apply(wj*newy,1,sum,na.rm=T)
ey2<-apply(wj*wj*sy2x,1,sum,na.rm=T)
ey<-sqrt(ey2)
ny[is.na(w)]<-NA
ey[is.na(w)]<-NA
ey2[is.na(w)]<-NA
s<-sd(residuals(mod))
t1<-min(ty[!is.na(ny)])
t2<-max(ty[!is.na(ny)])

#################################################################3

selectedmodels<-c((1:J)) # all models

# IMT estimates with prediction intervals for each model
par(mfrow=c(2,1),mar=c(4,4.4,1,1),cex=0.7,cex.axis=1.5,cex.lab=1.5)
for (j in selectedmodels)
{
nny<-newy[,j]
ney<-sy[,j]
cintu<-as.vector(na.omit(nny+1.96*ney))
cintl<-as.vector(na.omit(nny-1.96*ney))
t1<-min(ty[!is.na(nny)])
t2<-max(ty[!is.na(nny)])
xx <-c(t1:t2,rev(t1:t2))
yy <- c(cintu, rev(cintl))
err.fits[,j,1] <- as.vector((nny-1.96*ney))
err.fits[,j,2] <- as.vector((nny+1.96*ney))
mod.fits[,j] <- nny

err.fits[is.na(err.fits)] <- 0
mod.fits[is.na(mod.fits)] <- 0

# output final TSAM MMT estimate from 1961-2098
#tlo=2005
#thi=2048

tlo=2000
thi=2050

cry <- ty
for ( i in 1:ntime ){
  cry[i]=ifelse(cry[i] > thi , NA , cry[i] )
  cry[i]=ifelse(cry[i] < tlo , NA , cry[i] )
}

cny <- ny
cny[is.na(cry)] <- NA



# output MMT estimates and IMT estimates
s1 <- c(as.vector(ny-1.96*ey))
s2 <- c(as.vector(ny+1.96*ey))

s1[is.na(s1)] <- 0
s2[is.na(s2)] <- 0
ny[is.na(ny)] <- 0
# output MMT estimates and IMT estimates
s1 <- c(as.vector(ny-1.96*ey))
s2 <- c(as.vector(ny+1.96*ey))


s1[is.na(s1)] <- 0
s2[is.na(s2)] <- 0
ny[is.na(ny)] <- 0
ss1 <- c(as.vector(ny-1.96*sqrt(s^2+ey2)))
ss2 <- c(as.vector(ny+1.96*sqrt(s^2+ey2)))

ss1[is.na(ss1)] <- 0
ss2[is.na(ss2)] <- 0


gamout <- cbind(c(yearstart:yearend),ny,s1,s2,ss1,ss2)
#print(gamout)
#gamout <- cbind(c(yearstart:yearend),ny,s1,s2)
write.table(gamout, file = paste(fleout,"_MMT_estimate.ascii",sep=""),quote = FALSE,row.names = FALSE,col.names = FALSE)

for (j in selectedmodels){
  gamout <- cbind(c(yearstart:yearend),mod.fits[,j],err.fits[,j,1],err.fits[,j,2])
  write.table(gamout, file = paste(fleout,"_IMT_estimate_",name[j],".ascii",sep=""),quote = FALSE,row.names = FALSE,col.names = FALSE)
}

}


