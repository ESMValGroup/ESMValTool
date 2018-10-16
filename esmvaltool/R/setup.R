#check for present library paths
RLIBPATH <- .libPaths()

#check if we can write in the present R libaries paths
if (any(file.access(RLIBPATH, 2) == 0)) {
    #if possible, use the standard one for following instalation
    RLIBLOC <- RLIBPATH[which(file.access(RLIBPATH, 2) == 0)[1]]
} else {
    #if not possible, create a local library in the home directory
    RLIBLOC <- Sys.getenv("R_LIBS_USER")
    dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE,
               recursive = TRUE)
}

cat("\nINFO: Installing packages to --> ", RLIBLOC, "\n\n")

# define the R mirror to download packages
pkg_mirror <- "https://cloud.r-project.org"

# get the script path
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- sub(file_arg_name, "",
                   initial_options[grep(file_arg_name, initial_options)])
script_dirname <- dirname(script_name)
print(script_dirname)

# read the dependencies
# TODO: find a solution for script directory
dependencies <- read.table(paste(script_dirname, "r_requirements.txt",
                                 sep = "/"), header = TRUE)
package_list <- as.list(as.vector(dependencies[["package"]]))
package_list <- package_list[!(package_list %in%
                               installed.packages()[, "Package"])]

if (length(package_list) == 0) {
    print("All packages are installed!")
} else {
    print(paste("These packages will be installed: ", length(package_list)))
    print(package_list)
}

for (pack in package_list) {
    print(paste("Installing package --> ", pack))
    install.packages(pack, repos = pkg_mirror, dependencies = TRUE)
}

print("List of installed packages:")
installed.packages()
