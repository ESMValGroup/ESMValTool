###############################################################################
# CREATES A TARBALL OF THE CURRENT TOOL VERSION
###############################################################################
#
#   all files under version control are included
#   .svn files are excluded
#   to_be_checked/ directories are excluded
#   links are escluded
#
# csh util/create-tar.csh
#
# This creates a file ESMValTool_vX.Y.tar.
#
# Written by Mattia Righi (DLR, Germany)
#
###############################################################################

# Check current path
set pwd = `pwd | sed 's|/.*/||'`
if (`echo $pwd | grep util` != "") then
  echo "  ERROR: this script must be executed from the root path"
  echo "  csh util/cleanup/create-tar.csh"
  exit 1
endif

# Set tar file name
set ver = `more main.py | grep "version =" | awk -F '"' '{print $2}'`
set tar = "ESMValTool_v"${ver}.tar
if (-e $tar) rm -f $tar

# Set list of files under version control
set tmp = `svn ls --recursive`

# Remove "to_be_checked" paths from the list
set svn = ()
foreach t ($tmp)
  if (`echo $t | grep to_be_checked` == "" && ! -d $t && -e $t && ! -l $t) then
    set svn = ($svn $t)
  endif
end

# Set command
set exe = "tar cvf $tar $svn --exclude-vcs --transform 's,^,/ESMValTool_v${ver}/,'"

# Run tar
eval $exe
