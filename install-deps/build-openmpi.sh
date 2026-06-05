#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/openmpi-4.1.6
# Location where dependency sources will be downloaded
DEP_SOURCES=$GKYLSOFT/dep_src/

mkdir -p $DEP_SOURCES
cd $DEP_SOURCES

if [ "$DOWNLOAD_PKGS" = "yes" ]
then
    echo "Downloading OpenMPI .."
    # delete old checkout and builds
    rm -rf openmpi-4.1.6.tar* openmpi-4.1.6
    curl -L https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.6.tar.gz > openmpi-4.1.6.tar.gz
fi

if [ "$BUILD_PKGS" = "yes" ]
then
    echo "Building OpenMPI .."
    gunzip -c openmpi-4.1.6.tar.gz | tar xf -
    cd openmpi-4.1.6

    ./configure --prefix=$PREFIX --enable-mpi-fortran=yes CC=$CC CXX=$CXX
    make -j6 install 

    # soft-link 
    ln -sfn $PREFIX $GKYLSOFT/openmpi
fi
