#check for present library paths
RLIBPATH=.libPaths()
#check if we can write in the present R libaries paths
if (any(file.access(RLIBPATH,2)==0))
        {
        #if possible, use the standard one for following instalation
        RLIBLOC=RLIBPATH[which(file.access(RLIBPATH,2)==0)[1]]
        } else {
        #if not possible, create a local library in the home directory
        RLIBLOC=Sys.getenv("R_LIBS_USER")
        dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
        }

#preference is for web-based installation
web=1

if (web==1)
{
# web based installation
#packages=c("spam","maps","fields","ncdf4")
packages=c("maps","ncdf4","PCICt")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]

if (length(new.packages)==0) {print("All packages are there, no need to install anything")}
if (length(new.packages)!=0) {print(paste("Installing",length(new.packages),"R packages... follow on-screen instruction"))}

for (pack in new.packages)
{
print(paste("Installing... ",pack))
install.packages(pack,repos="http://cran.irsn.fr",type="source")
}
}


if (web==0)
{
#manual installation
#packages=c("spam_1.0-1.tar.gz","maps_2.3-9.tar.gz","fields_7.1.tar.gz","ncdf4_1.13.tar.gz")
packages=c("maps_2.3-9.tar.gz","ncdf4_1.13.tar.gz")
PROGDIR<-Sys.getenv(c("PROGDIR"))
R_PACKDIR=paste0(PROGDIR,"/R_packages")
for (pack in new.packages)
{
packageinstall=paste(R_PACKDIR,"/",pack,sep="")
install.packages(packageinstall,repos=NULL,type="source")
}
}

