#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/openmpi-5.0.7
# Location where dependency sources will be downloaded
DEP_SOURCES=$GKYLSOFT/dep_src/

mkdir -p $DEP_SOURCES
cd $DEP_SOURCES

if [ "$DOWNLOAD_PKGS" = "yes" ]
then
    echo "Downloading OpenMPI .."
    # delete old checkout and builds
    rm -rf openmpi-5.0.7.tar* openmpi-5.0.7
    curl -L https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.7.tar.gz > openmpi-5.0.7.tar.gz
fi

if [ "$BUILD_PKGS" = "yes" ]
then
    echo "Building OpenMPI .."
    gunzip -c openmpi-5.0.7.tar.gz | tar xf -
    cd openmpi-5.0.7

    ./configure --prefix=$PREFIX --enable-mpi-fortran=none CC=$CC CXX=$CXX
    make -j6 install 

    # soft-link 
    ln -sfn $PREFIX $GKYLSOFT/openmpi
fi
