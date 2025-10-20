#!/bin/bash

source ./build-opts.sh

# Edit to suite your system
PREFIX=$GKYLSOFT/superlu-7.0.0
# Location where dependency sources will be downloaded
DEP_SOURCES=$GKYLSOFT/dep_src/

mkdir -p $DEP_SOURCES
cd $DEP_SOURCES

if [ "$DOWNLOAD_PKGS" = "yes" ]
then
    echo "Downloading SuperLU .."
    # delete old checkout and builds
    rm -rf superlu-*
    curl -L https://github.com/xiaoyeli/superlu/archive/refs/tags/v7.0.0.tar.gz > superlu-7.0.0.tar.gz
fi

if [ "$BUILD_PKGS" = "yes" ]
then
    echo "Building SuperLU .."
    gunzip -f superlu-7.0.0.tar.gz
    tar xvf superlu-7.0.0.tar

    cd superlu-7.0.0
    mkdir build
    cd build

    # configure build
    cmake .. -DCMAKE_C_FLAGS="-g -O3 -fPIC" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_INSTALL_LIBDIR=lib -Denable_tests=NO -Denable_internal_blaslib=NO -DXSDK_ENABLE_Fortran=NO

    # build and install
    make -j 32 VERBOSE=1
    make install

    # soft-link 
    ln -sfn $PREFIX $GKYLSOFT/superlu
fi
