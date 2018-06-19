#check for present library paths
RLIBPATH=.libPaths()

#check if we can write in the present R libaries paths
if (any(file.access(RLIBPATH,2)==0)) {
    #if possible, use the standard one for following instalation
    RLIBLOC=RLIBPATH[which(file.access(RLIBPATH,2)==0)[1]]
} else {
    #if not possible, create a local library in the home directory
    RLIBLOC=Sys.getenv("R_LIBS_USER")
    dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
}

cat("\nINFO: Installing packages to --> ", RLIBLOC, "\n\n")

# list of dependencies
pkgMirror <- 'https://cloud.r-project.org'
# chooseCRANmirror(ind=1)
dependencies <-read.table("requirements.r", header=FALSE)

packageList <- dependencies[!(dependencies %in% installed.packages()[,"Package"])]

if (length(packageList)==0) {
    print("All packages are there, no need to install anything")
} else {
    print(paste("Installing ",length(packageList),"R packages... follow on-screen instruction."))
}

for (pack in packageList) {
    print(paste("Installing package --> ", pack))
    install.packages(pack, repos=pkgMirror, type="source")
}

