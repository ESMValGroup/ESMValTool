##############################################################################
# GENERAL ROUTINES FOR WRITING META DATA
##############################################################################
# Please consider using of extending existing routines before adding new ones.
# Check the header of each routine for documentation.
#
# Contents:
#    function ESMValMD
#
##############################################################################

source('interface_data/r.interface')
source('diag_scripts/lib/R/info_output.r') 

#############################################################################
ESMValMD <- function(filename,
                  tags,
                  caption,
                  id,
                  varnames,
                  modelnames,
                  infiles,
                  diag_name,
                  contrib_authors) {
#
# Arguments
#    filename: file name of the figure (including path)
#    tags: list of tags
#    caption: figure caption
#    id: id string
#    varnames: list of variable names
#    modelnames: list of model names
#    infiles: list of input files (= climo files)
#    diag_name: name of diagnostic
#    contib_authors: list of contributing authors
#
# Description
#    Creates meta data string and calls Python function to write meta data to
#    figure files.
#
# Caveats
#
# References
#
# Modification history:
#    20170922-A_bock_ls: written.
#

    funcname <- "ESMValMD"
    scriptname <- "diag_scripts/lib/R/meta_data.r"
    info_output(paste0("<<<<<<<< Entering ", funcname, " (", scriptname, ")"), verbosity, 6)

    m_tags = paste("M", modelnames, sep="_")
    v_tags = paste("V", varnames, sep="_")
    tags_plus = union(tags, m_tags)
    tags_plus = union(tags_plus, v_tags)
    rm(m_tags)
    rm(v_tags)

    n = length(tags_plus)
    str = ""
    str[1] = filename
    str[2] = "both"  
    str[3] = tags_plus[1]
    if (n > 1) {
        for (i in 2:n) {
            str[3] = paste(str[3], tags_plus[i], sep=",")
        }
    }
    str[4] = caption
    str[5] = id

    n = length(infiles)
    str[6] = infiles[1]
    if (n > 1) {
        for (i in 2:n) {
            str[6] = paste(str[6], infiles[i], sep=",")
        }
    }
    str[7] = diag_name

    n = length(contrib_authors)
    str[8] = contrib_authors[1]
    if (n > 1) {
        for (i in 2:n) {
            str[8] = paste(str[8], contrib_authors[i], sep=",")
        }
    }

    ascii_file = paste(filename, "list.txt", sep="_")
    write(str, ascii_file, sep=("\n"))
    rm(str)

    system(paste("python diag_scripts/lib/python/running_MD_for_ncl_with_file.py '",ascii_file, "'", sep=""))
#    system("rm '" + ascii_file + "'")

    info_output(paste0(">>>>>>>> Leaving ", funcname, " (", scriptname, ")"), verbosity, 6)
}
