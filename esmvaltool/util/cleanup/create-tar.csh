###############################################################################
# CREATES A TARBALL OF THE CURRENT TOOL VERSION
###############################################################################
#
# csh util/create-tar.csh
#
# This creates a file ESMValTool_vX.Y.Z.tar
#    all files under version control are included
#    .svn files are excluded
#    to_be_checked/ directories are excluded
#    links are escluded
#    the generic config_private.xml file is included
#    all the others config_private_* files are excluded
#    the include statements in all nmls is set to the generic config_private.xml
#    a generic ncl.interface file is added
#
# Written by Mattia Righi (DLR, Germany)
#
###############################################################################

# Check current path
echo -n "Checking current path... "
set pwd = `pwd | sed 's|/.*/||'`
if (`echo $pwd | grep util` != "") then
  echo "  ERROR: this script must be executed from the root path"
  echo "  csh util/cleanup/create-tar.csh"
  exit 1
endif
echo "OK"

# Check that all changes have been committed
echo -n "Checking that there no uncommited changes..."
set nn = `svn st | wc -l`
if ($nn != 0) then
  echo "ERROR: there are uncommitted changes"
  exit 1
endif
echo "OK"

# Add a generic ncl.interface file
echo "Adding a generic ncl.interface file to interface_data"
cp util/cleanup/ncl.interface interface_data

# Clean-up ~ files
echo "Removing ~ files"
find -name "*~" -exec rm -f "{}" \;

# Temporary replace include statement in nmls with generic config file
set nmls = `find . -type f -name "*namelist*.xml"`
foreach nml ($nmls)
  sed -i 's/config_private_DLR-PA2.xml/config_private.xml/g' $nml
end

# Set tar file name
echo -n "Setting tar file name... "
set ver = `more main.py | grep "version =" | awk -F '"' '{print $2}'`
set tar = "ESMValTool_v"${ver}.tar
if (-e $tar) rm -f $tar
echo $tar

# Set list of files under version control
set tmp = `svn ls --recursive`

# Create list of files
echo "Creating list of files"
set svn = ()
foreach t ($tmp)
  if (`echo $t | grep to_be_checked` == "" && \
      `echo $t | grep _DLR-PA2` == "" && \
      ! -d $t && -e $t && ! -l $t) then
    set svn = ($svn $t)
  endif
end
set svn = ($svn interface_data/ncl.interface)

# Set command
set exe = "tar cf $tar $svn --exclude-vcs --transform 's,^,/ESMValTool_v${ver}/,'"

# Run tar
echo "Running tar..."
eval $exe

# Revert changes to include statements in nmls
foreach nml ($nmls)
  svn revert $nml > /dev/null
end

if ($status == 0) echo "Done!"
