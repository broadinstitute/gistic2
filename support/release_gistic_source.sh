#!/bin/bash

################################
# GISTIC source release script #
################################
#
# This is the script to release GISTIC 2.0 source code, which
# is delivered in a directory tarball containing executable, source, 
# documentation, and a working example. To run it, type (e.g.):
#
#      release_gistic_source 2_0_23
#
# The argument is version text that will be appended to the tarball name.
#
#

if [ -z $1 ]
then
    echo "USAGE: $0 <release_tag>"
    echo "where <release_tag> is version text, e.g. 2_0_23" 
    exit
fi
TAG=$1

#-- define subversion branch work directory used for release --#
#
SVNROOT="/home/unix/$USER/CancerGenomeAnalysis"

#-- define master release directory --#
#
# Parent directory for the release image. It also contains files not 
# under version control that will be copied to the release image.
#
RELDIR="/xchip/gistic/GISTIC2.0"
# this directory will hold the release image for the compressed archive
RELEASE_IMAGE="$RELDIR/release_image"

#-- ensure that the release image directory does not already exist --#
if [ -d $RELEASE_IMAGE ]
then
    echo "release image directory [$RELEASE_IMAGE] already exists - aborting!"
    exit
fi
# create release image directory
mkdir $RELEASE_IMAGE

#-- get MCRInstaller.bin from Matlab release directory --
# (note that this will change with the target MATLAB version)
mkdir $RELEASE_IMAGE/MCR_Installer
cp /broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2014a/toolbox/compiler/deploy/glnxa64/MCRInstaller.zip $RELEASE_IMAGE/MCR_Installer 
#cp /broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2012b/toolbox/compiler/deploy/glnxa64/MCRInstaller.zip $RELEASE_IMAGE/MCR_Installer

# location of version control-based files for the release other than source
SVNRELDIR="$SVNROOT/trunk/gistic2core/support"
SVNSRCDIR="$SVNROOT/trunk/gistic2core/source"

#-- copy release files under version control to the release directory
# (https://svnrepos/CancerGenomeAnaysis/branches/gistic-2.0-genepattern/matlab/snp/gistic2/release_files/)

cp $SVNRELDIR/INSTALL.txt $RELEASE_IMAGE/
cp $SVNRELDIR/README.txt $RELEASE_IMAGE/
cp $SVNRELDIR/LICENSE.txt $RELEASE_IMAGE/
cp $SVNRELDIR/gistic2 $RELEASE_IMAGE/
cp $SVNRELDIR/run_gistic_example $RELEASE_IMAGE/
cp $SVNRELDIR/GISTICDocumentation_standalone.htm $RELEASE_IMAGE/
mkdir $RELEASE_IMAGE/GISTICDocumentation_standalone_files
cp $SVNRELDIR/GISTICDocumentation_standalone_files/* $RELEASE_IMAGE/GISTICDocumentation_standalone_files

#-- copy working example files to the release image --
mkdir $RELEASE_IMAGE/examplefiles
cp -p $SVNRELDIR/examplefiles/*.txt $RELEASE_IMAGE/examplefiles/
#-- copy refgene files to the release image --
mkdir $RELEASE_IMAGE/refgenefiles
cp -p $SVNRELDIR/refgenefiles/*.mat $RELEASE_IMAGE/refgenefiles/

#-- make empty directories --#
mkdir $RELEASE_IMAGE/example_results
mkdir $RELEASE_IMAGE/MATLAB_Compiler_Runtime

#-- copy source code --
rm -rf $RELEASE_IMAGE/source
cp -LR $SVNSRCDIR $RELEASE_IMAGE/source

#-- compile gp_gistic2_from_seg from source --#
cd $RELEASE_IMAGE/source
source /broad/software/scripts/useuse
reuse -q .matlab-2014a
mcc -v -m -w enable -I $RELEASE_IMAGE/source gp_gistic2_from_seg
mv gp_gistic2_from_seg ..

# clean up files generated by the compiler that we don't need
rm readme.txt
rm run_gp_gistic2_from_seg.sh
rm mccExcludedFiles.log

#Version updates (other than gistic_version.m):
# - name of tar file in INSTALL.txt
# - date and version number in GISTICDocumentation_standalone.htm

cd $RELEASE_IMAGE
# fix permissions
chmod -R 777 *
# compress image
tar czf GISTIC_$TAG.tar.gz .
chmod -R 666 GISTIC_$TAG.tar.gz