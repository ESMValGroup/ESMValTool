#check for present library paths
RLIBPATH = .libPaths()

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

# define the R mirror to download packages
pkgMirror <- 'https://cloud.r-project.org'

# get the script path
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
print(script.dirname)

# read the dependencies
dependencies <-read.table(paste(script.dirname,"r_requirements.txt", sep="/"), header=FALSE) # TODO: find a solution for script directory

packageList <- dependencies[!(dependencies %in% installed.packages()[,"Package"])]

if (length(packageList)==0) {
    print("All packages are installed!")
} else {
    print(paste("These packages will be installed: ", length(packageList)))
}

for (pack in packageList) {
    print(paste("Installing package --> ", pack))
    # install.packages(pack, repos = pkgMirror, type="source")
    install.packages(pack, repos = pkgMirror)
}
